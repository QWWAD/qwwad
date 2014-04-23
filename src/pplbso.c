/*=================================================================
     pplbso   Pseudo-Potential Large Basis with Spin-Orbit coupling
  =================================================================

   This program adds spin-orbit coupling to the original
   PseudoPotential Large Basis calculation.
 
   Note this code is written for clarity of understanding and not
   solely computational speed.  

   Input files:
		atoms.xyz	atomic species and positions
		G.r		reciprocal lattice vectors
		k.r		electron wave vectors (k)

   Output files:
		ank.r		eigenvectors	
		Ek?.r		eigenenergies for each k


   Paul Harrison, April 2000                                
 
								*/
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "struct.h"
#include "maths.h"
#include "qclsim-constants.h"

#include "ppff.h"	/* the PseudoPotential Form Factors	*/
#include "ppsop.h"	/* the Spin-Orbit Parameters		*/

typedef struct
{
 char	type[12];
 vector	r;
}atom;

int main(int argc,char *argv[])
{
complex	V();		/* potential component of H_GG			*/
complex	Vso();		/* spin-orbit component of H_GG			*/
atom	*read_atoms();	/* read in atomic positions/species		*/
void	zheev_();	/* Matrix diagonalization routine (LAPACK)	*/
void	write_ank();	/* writes eigenvectors to file			*/
vector	*read_rlv();	/* function to read reciprocal lattice vectors	*/

complex	double *ank;		/* coefficients of eigenvectors			*/
complex	double *Vc;		/* matrix containing the crystal potential	*/
complex	double T_GG;		/* kinetic energy component of H_GG		*/
complex	double *WORK;		/* LAPACK: workspace				*/
double	A0;		/* Lattice constant				*/
double	*E;		/* energy eigenvalues				*/
double	m_per_au;	/* unit conversion factor, m/a.u.		*/
double	*RWORK;		/* LAPACK: workspace				*/
int	N;		/* number of reciprocal lattice vectors		*/
int	Ns;		/* order of H_GG with spin (2*N)		*/
int     n_min;          /* lowest output band				*/
int     n_max;          /* highest output band				*/
int	n_atoms;	/* number of atoms in (large) cell		*/
int	INFO;		/* LAPACK: information integer			*/
int	i;		/* loop index for matrix rows			*/
int	iE;		/* loop index for energy eigenvalues		*/
int	ik;		/* loop index for k-vectors			*/
int	j;		/* loop index for matrix columns		*/
int	LWORK;		/* LAPACK: workspace				*/
char	filenameE[9];	/* character string for Energy output filename	*/
char	JOBZ='V';	/* LAPACK: compute eigenvectors and eigenvalues	*/
char	UPLO='L';	/* LAPACK: diagonalise `L'ower triangle 	*/
FILE	*Fk;		/* pointer to k.r file				*/
FILE	*FEk;		/* pointer to Ek.r file				*/
atom	*atoms;		/* the type and position of the atoms		*/
bool	ev;		/* flag, if set output eigenvectors 		*/
vector	*G;		/* reciprocal lattice vectors			*/
vector	k;		/* electron wave vector				*/
vector	q;		/* G'-G						*/

/* default values	*/

A0=5.65e-10;
ev=false;
n_min=0;
n_max=-1;
m_per_au=4*pi*eps0*gsl_pow_2(hBar/e)/me;

while((argc>1)&&(argv[1][0]=='-'))
{
 switch(argv[1][1])
 {
  case 'A':
	   A0=atof(argv[2])*1e-10;
           break;
  case 'n':
           n_min=atoi(argv[2])-1;         /* Note -1=>top VB=4, CB=5 */
           break;
  case 'm':
           n_max=atoi(argv[2])-1;         /* Note -1=>top VB=4, CB=5 */
           break;
  case 'w':
           ev=true;
	   argv--;
	   argc++;
           break;
  default :
	   printf("Usage:  pplbso [-A lattice constant (\033[1m5.65\033[0mA)]\n");
	   printf("               [-n # lowest band \033[1m1\033[0m][-m highest band \033[1m4\033[0m], output states\n");
	   printf("               [-w output eigenvectors (wavefunctions) in range n->m]\n");
	   exit(0);
 }
 argv++;
 argv++;
 argc--;
 argc--;
}

if((Fk=fopen("k.r","r"))==0)
 {fprintf(stderr,"Error: Cannot open input file 'k.r'!\n");exit(0);}

atoms=read_atoms(&n_atoms);		/* read in atomic basis	*/

G=read_rlv(A0,&N);	/* read in reciprocal lattice vectors	*/

Ns=2*N;			/* Matrix order with spin		*/

/* Allocate memory: SO doubles the order of the H_GG matrix	*/

LWORK=2*Ns;		/* LAPACK workspace definition		*/

WORK=(complex *)calloc(LWORK,sizeof(complex));
 if(WORK==0){fprintf(stderr,"Cannot allocate memory!\n");exit(0);}

RWORK=(double *)calloc(3*Ns-2,sizeof(double));
 if(RWORK==0){fprintf(stderr,"Cannot allocate memory!\n");exit(0);}

/* end of LAPACK definitions	*/

Vc=(complex *)calloc(Ns*Ns,sizeof(complex));
 if(Vc==0){fprintf(stderr,"Cannot allocate memory!\n");exit(0);}

E=(double *)calloc(Ns,sizeof(double));
 if(E==0){fprintf(stderr,"Cannot allocate memory!\n");exit(0);}

ank=(complex *)calloc(Ns*Ns,sizeof(complex));
 if(ank==0){fprintf(stderr,"Cannot allocate memory!\n");exit(0);}

/* Create crystal potential Vc	*/

for(i=0;i<N;i++)        /* index down rows */
 for(j=0;j<=i;j++)       /* index across cols, note creates leading diag. */
 {
  q.x=((G+i)->x)-((G+j)->x);
  q.y=((G+i)->y)-((G+j)->y);
  q.z=((G+i)->z)-((G+j)->z);

  *(Vc+i*Ns+j)=V(A0,m_per_au,atoms,n_atoms,q);    

  /* Now copy to lower triangle of bottom right quarter 
     (block 4, s'=-1/2, s=1/2)						*/ 

  *(Vc+(i+N)*Ns+N+j)=*(Vc+i*Ns+j);		
 }

/* The crystal potential is the same for all k-points so keep in matrix Vc.
   Use the fact that the LAPACK zheev() overwrites the Hamiltonian matrix 
   and start constructing H_GG in the eigenvector ank space		*/

/* Add k-dependent elements to matrix H_GG' */

ik=0;	/* Initialise k-loop index	*/
while((fscanf(Fk,"%lf %lf %lf",&k.x,&k.y,&k.z))!=EOF)
{
 k.x*=(2*pi/A0);k.y*=(2*pi/A0);k.z*=(2*pi/A0);

 for(i=0;i<Ns;i++) 		/* LAPACK routine zheev() writes vectors*/
  for(j=0;j<=i;j++)		/* over original matrix, so copy crystal*/
   *(ank+i*Ns+j)=*(Vc+i*Ns+j);	/* potential to eigenvector space	*/

 /* Now add kinetic energy components to leading diagonal		*/

 for(i=0;i<N;i++)        
 {
  T_GG=hBar*(gsl_pow_2(((G+i)->x)+k.x)+gsl_pow_2(((G+i)->y)+k.y)+gsl_pow_2(((G+i)->z)+k.z))*hBar/(2*me);
  ank[i*Ns+i] += T_GG;
  ank[(i+N)*Ns+N+i] = ank[i*Ns+i];	 	/* bottom right quarter	*/
 }

 /* Add spin-orbit components to lower triangle	*/

 for(i=0;i<Ns;i++) 		
  for(j=0;j<=i;j++)	
  {
   ank[i*Ns+j] = Vso(A0,m_per_au,atoms,G,n_atoms,k,i,j,N) + *(ank+i*Ns+j);
  }

 /* call LAPACK diagonalisation routine	*/

 ctranspose(ank,Ns);		/* because of FORTRAN LAPACK	*/
 zheev_(&JOBZ,&UPLO,&Ns,ank,&Ns,E,WORK,&LWORK,RWORK,&INFO);

 /* Output eigenvalues in a separate file for each k point */

 sprintf(filenameE,"Ek%i.r",ik);
 FEk=fopen(filenameE,"w");
 for(iE=n_min;iE<=n_max;iE++)fprintf(FEk,"%10.6f\n",*(E+iE)/e);
 fclose(FEk);

 /* Output eigenvectors */

 if(ev){
	ctranspose(ank,Ns);		/* because of FORTRAN LAPACK	*/
	write_ank(ank,ik,N,Ns,n_min,n_max);
        }

 ik++;	/* increment loop counter	*/
 
}/* end while*/
fclose(Fk);

free(G);
free(Vc);
free(E);
free(ank);
free(atoms);

free(WORK);free(RWORK);		/* The LAPACK definitions	*/

return EXIT_SUCCESS;
}/* end main */




atom
*read_atoms(n_atoms)

/* This function reads the atomic species (defined in the file as.r)
   into memory (addressed by the pointer as) and returns the start
   address of this block of memory and the number of lines	   */

int	*n_atoms;
{
 int    ia=0;
 FILE 	*Fatoms;        /* file pointer to wavefunction file       */
 atom	*atoms;		/* atomic definitions			*/

 if((Fatoms=fopen("atoms.xyz","r"))==0)
 {
  fprintf(stderr,"Error: Cannot open input file 'atoms.xyz'!\n");
  exit(0);
 }

 /* Read in the first line and hence the number of atoms	*/

 int n_read = fscanf(Fatoms,"%i",n_atoms);
 
 /* Allocate memory for atom definitions	*/
 if (n_read != 1)
 {
     fprintf(stderr, "Could not read number of atoms");
     exit(EXIT_FAILURE);
 }

 atoms=(atom *)calloc(*n_atoms,sizeof(atom));
 if(atoms==0)
 {
  fprintf(stderr,"Cannot allocate memory!\n");
  exit(0);
 }
 
 while((fscanf(Fatoms,"%s %lf %lf %lf",&(atoms+ia)->type,
        &(atoms+ia)->r.x,&(atoms+ia)->r.y,&(atoms+ia)->r.z))!=EOF)
 {
  /* Convert atomic positions from Angstrom into S.I. units	*/

  (atoms+ia)->r.x*=1e-10;(atoms+ia)->r.y*=1e-10;(atoms+ia)->r.z*=1e-10;
  ia++;
 }
 fclose(Fatoms);

 return(atoms);
}



vector
*read_rlv(A0,N)

/* This function reads the reciprocal lattice vectors (defined in
   the file G.r) into the array G[] and then converts into SI units */

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
 while(fscanf(FG,"%*f %*f %*f")!=EOF)
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



complex 
V(A0,m_per_au,atoms,n_atoms,q)

double A0;	/* Lattice constant */
double m_per_au;/* conversion factor from SI to a.u. */
atom   *atoms;	/* atomic definitions	*/
int    n_atoms; /* number of atoms in structure */
vector q;       /* a reciprocal lattice vector, G'-G */
{
 extern double Vf();   /* Fourier transform of potential	*/
 complex double v;     /* potential					*/
 double	vf;	/* storage for returned data from Vf()		*/
 int	ia;	/* index across atoms				*/
 vector t;      /* general vertor representing atom within cell	*/

 v = 0;

 for(ia=0;ia<n_atoms;ia++)
 {
  t.x=(atoms+ia)->r.x;
  t.y=(atoms+ia)->r.y;
  t.z=(atoms+ia)->r.z;
  vf=Vf(A0,m_per_au,q.x*q.x+q.y*q.y+q.z*q.z,(atoms+ia)->type);
  v += cos(q.x*t.x+q.y*t.y+q.z*t.z)*vf - I * sin(q.x*t.x+q.y*t.y+q.z*t.z)*vf;
 }

 v/=(double)(n_atoms/2);

 return(v);
}



complex 
Vso(A0,m_per_au,atoms,G,n_atoms,k,i,j,N)

double A0;	/* Lattice constant 				*/
double m_per_au;/* conversion factor from SI to a.u. 		*/
atom   *atoms;	/* atomic definitions				*/
vector *G;	/* the reciprocal lattice vectors		*/
int    n_atoms; /* number of atoms in structure 		*/
vector k;       /* the electron momentum			*/
int    i;	/* index					*/
int    j;	/* index					*/
int    N;	/* the number of reciprocal lattice vectors	*/
{
 double lambda();   /* the spin orbit interaction		*/
 complex double Lambda;/* complex lambda()---for ease of multiplication*/
 complex double v;     /* potential					*/
 complex double vso;	/* intermediate value of spin-orbit interaction	*/
 int	ia;	/* index across atoms				*/
 vector A;	/* the vector product (G'+k)x(G+k)		*/
 vector q;	/* a reciprocal lattice vector, G'-G 		*/
 vector t;      /* general vector representing atom within cell	*/

 vso = 0;

 if((i<N)&&(j<N))	/* Block 1	*/
 {
  A=vvprod(vadd(*(G+i),k),vadd(*(G+j),k));
  A.x/=gsl_pow_2(2*pi/A0);A.y/=gsl_pow_2(2*pi/A0);A.z/=gsl_pow_2(2*pi/A0);
  vso = A.z;
  q=vsub(*(G+i),*(G+j));
 }
 if((i<N)&&(j>=N))	/* Block 2	*/
 {
  A=vvprod(vadd(*(G+i),k),vadd(*(G+j-N),k));
  A.x/=gsl_pow_2(2*pi/A0);A.y/=gsl_pow_2(2*pi/A0);A.z/=gsl_pow_2(2*pi/A0);
  vso = A.x - I * A.y;
  q=vsub(*(G+i),*(G+j-N));
 }
 if((i>=N)&&(j<N))	/* Block 3	*/
 {
  A=vvprod(vadd(*(G+i-N),k),vadd(*(G+j),k));
  A.x/=gsl_pow_2(2*pi/A0);A.y/=gsl_pow_2(2*pi/A0);A.z/=gsl_pow_2(2*pi/A0);
  vso = A.x + I * A.y;
  q=vsub(*(G+i-N),*(G+j));
 }
 if((i>=N)&&(j>=N))	/* Block 4	*/
 {
  A=vvprod(vadd(*(G+i-N),k),vadd(*(G+j-N),k));
  A.x/=gsl_pow_2(2*pi/A0);A.y/=gsl_pow_2(2*pi/A0);A.z/=gsl_pow_2(2*pi/A0);
  vso = -A.z;
  q=vsub(*(G+i-N),*(G+j-N));
 }	
 vso *= I;			/* -i(G'+k)x(G+k).sigma	*/

 v=0;				/* Initialise for sum	*/
 for(ia=0;ia<n_atoms;ia++)
 {
  t.x=(atoms+ia)->r.x;
  t.y=(atoms+ia)->r.y;
  t.z=(atoms+ia)->r.z;
  Lambda = lambda((atoms+ia)->type);
  v += Lambda * cexp(-q.x*t.x-q.y*t.y-q.z*t.z);
 }

 v *= vso;
 v/=(double)(n_atoms);

 return(v);
}



void
write_ank(ank,ik,N,Ns,n_min,n_max)

/* This function writes the eigenvectors (a_nk(G)) to the files ank.r */

complex	double *ank;
int	ik;		/* k point identifier				*/
int	N;
int	Ns;
int     n_min;          /* lowest output band				*/
int     n_max;          /* highest output band				*/
{
 int	iG;		/* index over G vectors				*/
 int	in;		/* index over bands				*/
 char	filename[9];	/* eigenfunction output filename		*/
 FILE 	*Fank;		/* file pointer to eigenvectors file		*/

sprintf(filename,"ank%i.r",ik);
Fank=fopen(filename,"w");
 
for(iG=0;iG<N;iG++)
{
 for(in=n_min;in<=n_max;in++)
  fprintf(Fank,"%20.16le %20.16le ", creal(ank[iG*Ns+in]), cimag(ank[iG*Ns+in]));
 fprintf(Fank,"\n");
}


fclose(Fank);
}

