/*=================================================================
              ppsl   PseudoPotential SuperLattice calculation
  =================================================================

   This program performs the PseudoPotential SuperLattice as a perturbation
   calculation on a user-defined cell, the atomic species of which are defined
   in the file atoms.xyz (XYZ format file).
 
   Note this code is written for clarity of understanding and not
   solely computational speed.  

   Input files:
		ank?.r		bulk eigenvectors a_nk(G)
		Ek?.r		bulk eigenvalues E_nk, for each k
		atoms.xyz	atomic species and positions of the
				unperturbed lattice
		atomsp.xyz	atomic species and positions of the 
				perturbed lattice
		G.r		reciprocal lattice vectors
		k.r		electron wave vectors (k)

   Output files:
   		Exi.r		superlattice eigenvalues E_xi


   Paul Harrison, October 1998

   Electric field capability added 27th October 1998.

   Modifications, March 1999

   Modifications, January 2000
	Use of LAPACK diagonalisation routine

								*/

#include <stdbool.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "struct.h"
#include "maths.h"
#include "qclsim-constants.h"

#include "ppff.h"

typedef struct
{
 char	type[12];
 vector	r;
}atom;

int main(int argc,char *argv[])
{
complex double	V();		/* potential component of Hdash			*/
complex double	VF();		/* field potential component of Hdash		*/
complex double	*read_ank();	/* reads in all the bulk eigenvectors ank(G)	*/
double	*read_Enk();	/* reads in all the bulk eigenvalues Enk	*/
atom	*read_atoms();	/* read in atomic positions/species		*/
void	clean_Hdash();	/* removes round-up errors => makes H' Hermitian*/
void	zheev_();	/* Matrix diagonalization routine (LAPACK)	*/
void	write_VF();	/* writes the FT of the electric field potential*/
int	read_ank0();	/* deduces the number of bands in calculation	*/
vector	*read_kxi();	/* reads in kxi points				*/
vector	*read_rlv();	/* function to read reciprocal lattice vectors	*/

complex double	*Ank;		/* coefficients of eigenvectors			*/
complex double	*ank;		/* coefficients of bulk eigenvectors		*/
complex double	*Hdash;		/* components of H' (see notes)			*/
complex double	*WORK;		/* LAPACK: workspace				*/
double	A0;		/* Lattice constant				*/
double	*Enk;		/* bulk energy eigenvalues			*/
double	*Exi;		/* energy eigenvalues				*/
double	F;		/* the electric field strength			*/
double	m_per_au;	/* unit conversion factor, m/a.u.		*/
double	q;		/* the carrier charge, e=-e_0, h=+e_0		*/
double	*RWORK;		/* LAPACK: workspace				*/
int	INFO;		/* LAPACK: information integer			*/
int	N;		/* number of reciprocal lattice vectors		*/
int	Nn;		/* number of bands in calculation		*/
int	Nkxi;		/* number of k points in calculation		*/
int	n_atoms;	/* number of atoms in (large) cell		*/
int     n_min;          /* lowest output band                           */
int     n_max;          /* highest output band                          */
int	i;		/* loop index for matrix rows			*/
int	iE;		/* loop index for energy eigenvalues		*/
int	iG;		/* index over G					*/
int	iGdash;		/* index over G'				*/
int	ikxi;		/* index for kxi 				*/
int	ikxidash;	/* index for kxidash 				*/
int	in;		/* index over the bulk band n			*/
int	indash;		/* index over the bulk band n'			*/
int	j;		/* loop index for matrix columns		*/
int	LWORK;		/* LAPACK: workspace				*/
int	OH;		/* the order of Hdash				*/
char	filename[12];	/* character string for Energy output filename	*/
char	JOBZ='V';	/* LAPACK: compute eigenvectors and eigenvalues	*/
char	p;		/* the particle					*/
char	UPLO='L';	/* LAPACK: diagonalise `L'ower triangle 	*/
FILE	*FExi;		/* pointer to output energy file `Exi.r'	*/
atom	*atoms;		/* the type and position of the atoms		*/
atom	*atomsp;	/* the type and position of the perturbed atoms	*/
vector	*G;		/* reciprocal lattice vectors			*/
vector	*kxi;		/* set of electron wave vectors			*/
vector	g;		/* G'-G+kxi'-kxi				*/
bool	o;		/* if set, output the Fourier Transform VF(g)	*/

/* default values	*/

A0=5.65e-10;
F=0.0;
q=e;
n_min=0;
n_max=3;
o=false;
p='h';
m_per_au=4*pi*eps0*gsl_pow_2(hBar/e)/me;

while((argc>1)&&(argv[1][0]=='-'))
{
 switch(argv[1][1])
 {
  case 'A':
	   A0=atof(argv[2])*1e-10;
           break;
  case 'f':
	   F=atof(argv[2])*1e+5;	/* convert kV/cm->V/m	*/
           break;
  case 'n':
           n_min=atoi(argv[2])-1;	/* Note -1=>top VB=4, CB=5 */
           break;
  case 'm':
           n_max=atoi(argv[2])-1;	/* Note -1=>top VB=4, CB=5 */
           break;
  case 'o':
           o=true;
 	   argv--;
	   argc++;
           break;
  case 'p':
           p=*argv[2];
           switch(p)
           {
            case 'e': q=-e;break;
            case 'h': q=+e;break;
            default:  printf("Usage:  ppsl [-p particle (e or \033[1mh\033[0m)]\n");
                      exit(0);
           }
           break;

  default :
	   printf("Usage:  ppsl [-A lattice constant (\033[1m5.65\033[0mA)][-f (\033[1m0\033[0mkV/cm)]\n");
	   printf("             [-n # lowest band \033[1m1\033[0m][-m highest band \033[1m4\033[0m], output eigenvalues\n");
	   printf("             [-o output field FT][-p particle (e or \033[1mh\033[0m)]\n");
	   exit(0);
 }
 argv++;
 argv++;
 argc--;
 argc--;
}
strcpy(filename,"atoms.xyz");
atoms=read_atoms(&n_atoms,filename);	/* read in atomic basis	*/
strcpy(filename,"atomsp.xyz");
atomsp=read_atoms(&n_atoms,filename);	/* read in perturbation	*/

G=read_rlv(A0,&N);	/* read in reciprocal lattice vectors	*/
Nn=read_ank0(N);	/* reads a single ank.r file just to 
			   deduce the number of bands Nn	*/

kxi=read_kxi(A0,&Nkxi);	/* read in set of kxi points		*/

ank=read_ank(N,Nn,Nkxi);/* read in bulk eigenvectors		*/
Enk=read_Enk(Nn,Nkxi);	/* read in bulk eigenvalues		*/

/* Output the  Fourier Transform of the electric field	*/

if(o) write_VF(A0,F,q,atoms,n_atoms);

/* Allocate memory */

LWORK=2*Nn*Nkxi;		/* LAPACK workspace definition		*/

WORK=(complex double *)calloc(LWORK,sizeof(complex double));
 if(WORK==0){fprintf(stderr,"Cannot allocate memory!\n");exit(0);}

RWORK=(double *)calloc(3*Nn*Nkxi-2,sizeof(double));
 if(RWORK==0){fprintf(stderr,"Cannot allocate memory!\n");exit(0);}

/* end of LAPACK definitions	*/

Hdash=(complex double *)calloc(Nn*Nkxi*Nn*Nkxi,sizeof(complex double));
 if(Hdash==0){fprintf(stderr,"Cannot allocate memory!\n");exit(0);}

Exi=(double *)calloc(Nn*Nkxi,sizeof(double));
 if(Exi==0){fprintf(stderr,"Cannot allocate memory!\n");exit(0);}

Ank=(complex double *)calloc(Nn*Nkxi*Nn*Nkxi,sizeof(complex double));
 if(Ank==0){fprintf(stderr,"Cannot allocate memory!\n");exit(0);}


/* Create H' matrix elements	*/

for(i=0;i<Nn*Nkxi;i++)	/* index down rows, recall order of Hdash is Nn*Nkxi */
{
 for(j=0;j<=i;j++)	/* index across cols, including diagonal elements */
 {
  /* Deduce the indices of of the bulk eigenvectors for each H' element	*/

  ikxidash=i/Nn;
  ikxi=j/Nn;

  indash=i%Nn;
  in=j%Nn;

  /* Initialise the matrix element before the sum and add energy
     eigenvalues as specified by delta functions 	*/

  if((indash==in)&&(ikxidash==ikxi))
   Hdash[i*Nn*Nkxi+j] = Enk[ikxi*Nn+in];
  else
   Hdash[i*Nn*Nkxi+j] = 0;

  for(iGdash=0;iGdash<N;iGdash++)	/* sum over G'	*/
  {
   for(iG=0;iG<N;iG++)			/* sum over G	*/
   {
    /* Calculate appropriate g vector	*/

    g.x=((G+iGdash)->x)-((G+iG)->x)+((kxi+ikxidash)->x)-((kxi+ikxi)->x);
    g.y=((G+iGdash)->y)-((G+iG)->y)+((kxi+ikxidash)->y)-((kxi+ikxi)->y);
    g.z=((G+iGdash)->z)-((G+iG)->z)+((kxi+ikxidash)->z)-((kxi+ikxi)->z);

    /* Add on potential term	*/

    *(Hdash+i*Nn*Nkxi+j)=
                      *(Hdash+i*Nn*Nkxi+j)
		      +
		      (
                                 ( conj(*(ank+ikxidash*N*Nn+iGdash*Nn+indash))

				  *

                                  *(ank+ikxi*N*Nn+iG*Nn+in)
				  )
                                 *
				 (
				 V(A0,m_per_au,atoms,atomsp,n_atoms,g) +
				 VF(A0,F,q,atoms,n_atoms,g)
				 )
		      );
   } /* end iG */
  } /* end iGdash */


 } /* end j */
} /* end i */

/* Clean up matrix H'	*/

clean_Hdash(Hdash,Nn*Nkxi);

for(i=0;i<Nn*Nkxi;i++)		/* LAPACK routine zheev() writes vectors */
 for(j=0;j<=i;j++)		/* over original matrix, so copy matrix  */
  *(Ank+i*Nn*Nkxi+j)=*(Hdash+i*Nn*Nkxi+j);/* to eigenvector space	 */

/* call LAPACK diagonalisation routine	*/

OH=Nn*Nkxi;			/* the order of Hdash	*/
ctranspose(Ank,Nn*Nkxi);	/* because of FORTRAN LAPACK	*/
zheev_(&JOBZ,&UPLO,&OH,Ank,&OH,Exi,WORK,&LWORK,RWORK,&INFO);

/* Output eigenvalues in a separate file for each k point */

FExi=fopen("Exi.r","w");
for(iE=n_min;iE<=n_max;iE++)fprintf(FExi,"%10.6f\n",*(Exi+iE)/e);
fclose(FExi);


free(G);
free(kxi);
free(Hdash);
free(Exi);
free(Ank);
free(Enk);
free(atoms);
free(atomsp);

return EXIT_SUCCESS;
}/* end main */




void
clean_Hdash(Hdash,O)

/* This function accounts for any computational errors in the diagonal of
   H', it ensures that the imaginary components of these elements are zero	
 									*/

complex	double *Hdash;	/* The matrix to be `cleaned'	*/
int	O;	/* the order of the matrix	*/

{
 int	i;	/* index over diagonal elements	*/

 for(i=0;i<O;i++)
 {
  if(cimag(Hdash[i*O+i])/creal(Hdash[i*O+i])<(1e-10))
  Hdash[i*O+i] = creal(Hdash[i*O+i]);	   
 }

}



atom
*read_atoms(n_atoms,filename)

/* This function reads the atomic species and positions defined in XYZ
   format from the `filename', either `atoms.xyz', or `atomsp.xyz'	*/

int	*n_atoms;
char	filename[];
{
 int    ia=0;
 FILE 	*Fatoms;        /* file pointer to wavefunction file	*/
 atom	*atoms;		/* atomic definitions			*/

 if((Fatoms=fopen(filename,"r"))==0)
 {
  fprintf(stderr,"Error: Cannot open input file '%s'!\n",filename);
  exit(0);
 }

 /* Read in the first line and hence the number of atoms	*/

 fscanf(Fatoms,"%i",n_atoms);
 
 /* Allocate memory for atom definitions	*/

 atoms=(atom *)calloc(*n_atoms,sizeof(atom));
 if(atoms==0)
 {
  fprintf(stderr,"Cannot allocate memory!\n");
  exit(0);
 }
 
 while((fscanf(Fatoms,"%s %lf %lf %lf",atoms[ia].type,
        &(atoms+ia)->r.x,&(atoms+ia)->r.y,&(atoms+ia)->r.z))!=EOF)
 {
  /* Convert atomic positions from Angstrom into S.I. units	*/

  (atoms+ia)->r.x*=1e-10;(atoms+ia)->r.y*=1e-10;(atoms+ia)->r.z*=1e-10;
  ia++;
 }
 fclose(Fatoms);

 return(atoms);
}



int
read_ank0(N)

/* This function reads the first of the ank.r files, i.e., ank0.r, just
   just with the purpose of deducing the number of bands included in the
   calculation
   									*/

int	N;		/* The number of terms in each eigenvector	*/

{
 int	n;		/* counter for the number of elements in file	*/
 int	Nn;		/* number of bands in calculation		*/
 FILE	*Fank;		/* file pointer to eigenvectors file		*/

 if((Fank=fopen("ank0.r","r"))==0)
{
 fprintf(stderr,"Error: Cannot open input file 'ank0.r'!\n");
 exit(0);
}

/* Deduce number of complexes in file and hence number of bands	*/

n=0;
while(fscanf(Fank,"%*lf %*lf")!=EOF)
 n++;

/* The number of bands Nn is therefore the total number of elements divided
   by the number of terms in each eigenvector				*/

Nn=n/N;

return(Nn);

}


complex double
*read_ank(N,Nn,Nkxi)

/* This function reads the eigenvectors (a_nk(G)) from the file a_nk.r
   created by the code pplb.c into the array a_nk[N][Nn]		*/

int     N;		/* The number of terms in each eigenvector	*/
int     Nn;		/* The number of bands in file			*/
int	Nkxi;		/* number of k points in calculation		*/
{
 int	in;		/* index across bands				*/
 int	iG;		/* index across G vectors			*/
 int	ikxi;		/* index across kxi				*/
 char	filename[9];	/* eigenfunction output filename                */
 complex double       *ank;
 FILE   *Fank;		/* file pointer to eigenvectors file		*/

/* Allocate memory for eigenvectors	*/

ank=(complex double *)calloc(N*Nn*Nkxi,sizeof(complex double));
if(ank==0){fprintf(stderr,"Cannot allocate memory!\n");exit(0);}

/* Finally read eigenvectors into structure	*/

for(ikxi=0;ikxi<Nkxi;ikxi++)
{
 /* Open the bulk eigenvector file at each of the bulk k (kxi) points	*/

 sprintf(filename,"ank%i.r",ikxi);
 if((Fank=fopen(filename,"r"))==0)
 {fprintf(stderr,"Error: Cannot open input file 'ank%i.r'!\n",ikxi);exit(0);}

 /* Read in data	*/

 for(iG=0;iG<N;iG++)
  for(in=0;in<Nn;in++)
  {
   double temp_re=0.0;
   double temp_im=0.0;
   fscanf(Fank,"%lf %lf", &temp_re, &temp_im);
   ank[ikxi*N*Nn+iG*Nn+in] = temp_re + I * temp_im;
  }

 fclose(Fank);
}

return(ank);
}



double
*read_Enk(Nn,Nkxi)

/* This function reads the eigenvalues (E_nk) from the files Ek?.r
   created by the code pplb.c 						*/

int     Nn;		/* The number of bands in file			*/
int	Nkxi;		/* number of k points in calculation		*/
{
 int	in;		/* index across bands				*/
 int	ikxi;		/* index across kxi				*/
 char	filename[9];	/* eigenfunction output filename                */
 double	*Enk;
 FILE   *FEnk;		/* file pointer to eigenvectors file		*/

/* Allocate memory for eigenvalues	*/

Enk=(double *)calloc(Nn*Nkxi,sizeof(double));
if(Enk==0){fprintf(stderr,"Cannot allocate memory!\n");exit(0);}

/* Read eigenvalues into structure	*/

for(ikxi=0;ikxi<Nkxi;ikxi++)
{
 /* Open the bulk eigenvalues file at each of the bulk k (kxi) points	*/

 sprintf(filename,"Ek%i.r",ikxi);
 if((FEnk=fopen(filename,"r"))==0)
 {fprintf(stderr,"Error: Cannot open input file 'Ek%i.r'!\n",ikxi);exit(0);}

 /* Read in data	*/

 for(in=0;in<Nn;in++)
 {
  fscanf(FEnk,"%lf",(Enk+ikxi*Nn+in));
  *(Enk+ikxi*Nn+in)*=e;		/* convert from eV->S.I.	*/
 }

 fclose(FEnk);
}

return(Enk);
}


vector
*read_rlv(A0,N)

/* This function reads the reciprocal lattice vectors (defined in
   the file G.r) then converts into SI units */

double	A0;
int     *N;
{
 int    i=0;
 vector *G;
 FILE   *FG;           /* file pointer to wavefunction file */

 if((FG=fopen("G.r","r"))==0)
 {
  fprintf(stderr,"Error: Cannot open input file 'G.r'!\n");
  exit(0);
 }

 *N=0;
 while(fscanf(FG,"%*lf %*lf %*lf")!=EOF)
  (*N)++;
 rewind(FG);

 G=(vector *)calloc(*N,sizeof(vector));
 if (G==0)  {
  fprintf(stderr,"Cannot allocate memory!\n");
  exit(0);
 }

 while(fscanf(FG,"%lf %lf %lf",&(G+i)->x,&(G+i)->y,&(G+i)->z)!=EOF)
  {
   (G+i)->x*=(2*pi/A0);(G+i)->y*=(2*pi/A0);(G+i)->z*=(2*pi/A0);
   i++;
  }

 fclose(FG);
 
 return(G);
}


vector
*read_kxi(A0,Nkxi)

/* This function reads the set of electron momenta vectors kxi (defined
   in the file k.r---as needed to produce the bulk band structure) into the 
   allocated memory space accessed via the pointer kxi, then converting them
   from units of (2*pi/A0) into S.I. 				*/

double	A0;
int     *Nkxi;
{
 int    i=0;
 vector *kxi;
 FILE   *Fkxi;           /* file pointer to wavefunction file */

 if((Fkxi=fopen("k.r","r"))==0)
 {
  fprintf(stderr,"Error: Cannot open input file 'k.r'!\n");
  exit(0);
 }

 *Nkxi=0;
 while(fscanf(Fkxi,"%*lf %*lf %*lf")!=EOF)
  (*Nkxi)++;
 rewind(Fkxi);

 kxi=(vector *)calloc(*Nkxi,sizeof(vector));
 if (kxi==0)  {
  fprintf(stderr,"Cannot allocate memory!\n");
  exit(0);
 }

 while(fscanf(Fkxi,"%lf %lf %lf",&(kxi+i)->x,&(kxi+i)->y,&(kxi+i)->z)!=EOF)
  {
   (kxi+i)->x*=(2*pi/A0);(kxi+i)->y*=(2*pi/A0);(kxi+i)->z*=(2*pi/A0);
   i++;
  }

 fclose(Fkxi);
 
 return(kxi);
}



complex double
V(A0,m_per_au,atoms,atomsp,n_atoms,g)

double A0;	/* Lattice constant 		*/
double m_per_au;/* conversion factor from SI to a.u. */
atom   *atoms;	/* atomic definitions		*/
atom   *atomsp;	/* atomic definitions		*/
int    n_atoms; /* number of atoms in structure */
vector g;       /* the vector, g=G'-G+kxi'-kxi	*/
{
 extern double Vf();   /* Fourier transform of potential	*/
 complex double v;     /* potential					*/
 double	vf;	/* storage for returned data from Vf()		*/
 double	vfdash;	/* storage for returned data from Vf()		*/
 int	ia;	/* index across atoms				*/
 vector t;      /* general vertor representing atom within cell	*/

 v=0;

 for(ia=0;ia<n_atoms;ia++)
 {
  t.x=(atoms+ia)->r.x;
  t.y=(atoms+ia)->r.y;
  t.z=(atoms+ia)->r.z;
  vf=Vf(A0,m_per_au,g.x*g.x+g.y*g.y+g.z*g.z,(atoms+ia)->type);
  vfdash=Vf(A0,m_per_au,g.x*g.x+g.y*g.y+g.z*g.z,(atomsp+ia)->type);
  v+=cos(g.x*t.x+g.y*t.y+g.z*t.z)*(vfdash-vf) - I * sin(g.x*t.x+g.y*t.y+g.z*t.z)*(vfdash-vf);
 }

 /* These last divisions represent Omega_c/Omega_sl included here for
    convenience	*/
 v/=(double)(n_atoms/2);

 return(v);
}



complex double
VF(A0,F,q,atoms,n_atoms,g)

double A0;	/* Lattice constant 		*/
double F;	/* the electric field strength	*/
double q;	/* the carrier charge		*/
atom   *atoms;	/* atomic definitions		*/
int    n_atoms; /* number of atoms in structure */
vector g;       /* the vector, g=G'-G+kxi'-kxi	*/
{
 complex double v;	/* potential					*/
 complex double v1;	/* first term in potential			*/
 complex double v2;	/* second term in potential			*/
 double	z0;	/* the midpoint of the unit cell		*/
 double	zero;	/* an effective `0'				*/
 int	n_z;	/* number of lattice constants in each period	*/

 /* Calculate the number of lattice constants in each superlattice period */

 n_z=n_atoms/4;

 /* The smallest reciprocal lattice vector is therefore pi/(n_z*A0), thus 
    anything less than this must be zero				*/
 
 zero=pi/(n_z*A0)/10;

 v=0;

 /* Immitates the behaviour of the product of the two delta-functions	*/ 

 if((fabs(g.x)>zero)||(fabs(g.y)>zero)) return(v);
 else
 {
  z0=(((atoms+n_atoms-1)->r.z)+(atoms->r.z))/2;

  /* Evaluate the first and second terms in the integral	*/

  if(fabs(g.z)>zero)
  {
   v1=1/gsl_pow_2(g.z) + I*(((atoms+n_atoms-1)->r.z)-z0)/g.z;
   v1*=cexp(-g.z*((atoms+n_atoms-1)->r.z));

   v2=1/gsl_pow_2(g.z) + I*((atoms->r.z)-z0)/g.z;
   v2*=cexp(-g.z*(atoms->r.z));
  }
  else
  {
   v1=gsl_pow_2((atoms+n_atoms-1)->r.z)/2-((atoms+n_atoms-1)->r.z)*z0;
   v2=gsl_pow_2(atoms->r.z)/2-(atoms->r.z)*z0;
  }


  /* Combine the two */
  v=v1-v2;

  /* Multiply through by the scaling factor	*/
  v*=-q*F/((double)(n_z)*A0);

  return(v);
 }
}


void
write_VF(A0,F,q,atoms,n_atoms)

double A0;	/* Lattice constant 		*/
double F;	/* the electric field strength	*/
double q;	/* the carrier charge		*/
atom   *atoms;	/* atomic definitions		*/
int    n_atoms; /* number of atoms in structure */

{
 complex double	VF();	
 vector g;		/* the vector, g=G'-G+kxi'-kxi	*/
 FILE	*FVFg;		/* pointer to output file	*/
 int	i;		/* index			*/

 FVFg=fopen("VFg.r","w");

 for(i=-500;i<500;i++)
 {
  g.x=0;g.y=0;g.z=((float)i*4/500)*pi/A0;
  fprintf(FVFg,"%le %le\n",g.z/(pi/A0),cabs(VF(A0,F,q,atoms,n_atoms,g)));
 }

 fclose(FVFg);

}
