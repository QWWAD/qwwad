/*=================================================================
              pplb   Pseudo-Potential Large Basis 
  =================================================================

   This program performs a PseudoPotential Large Basis calculation
   on a user-defined cell, the atomic species of which are defined
   in the file atoms.xyz (XYZ format file).
 
   Note this code is written for clarity of understanding and not
   solely computational speed.  

   Input files:
		atoms.xyz	atomic species and positions
		G.r		reciprocal lattice vectors
		k.r		electron wave vectors (k)

   Output files:
		ank.r		eigenvectors	
		Ek?.r		eigenenergies for each k


   Paul Harrison, November 1994                                
 
   Major modifications, 1996
	Command line arguments

   Major Modifications, July 1998
	Move to XYZ format 
	Reading in atomic positions from a general basis

   Tweak, July 1999
   	Include form factor file directly
								*/

#include <Nag/nag.h>
#include <stdio.h>
#include <Nag/nag_stdlib.h>
#include <Nag/nagf02.h>
#include <math.h>
#include <stdlib.h>
#include "struct.h"
#include "maths.h"
#include "const.h"
#include "bools.h"

#include "ppff.c"

typedef struct
{
 char	type[12];
 vector	r;
}atom;

main(int argc,char *argv[])
{
Complex	V();		/* potential component of H_GG			*/
atom	*read_atoms();	/* read in atomic positions/species		*/
void	f02axc();	/* Matrix diagonalization routine (NAG)		*/
void	write_ank();	/* writes eigenvectors to file			*/
vector	*read_rlv();	/* function to read reciprocal lattice vectors	*/

Complex	*ank;		/* coefficients of eigenvectors			*/
Complex	*H_GG;		/* components of H_G'G (see notes)		*/
Complex	T_GG;		/* kinetic energy component of H_GG		*/
Complex	*V_GG;		/* potential energy of H_GG, diagonal elements	*/
			/* stored here for efficiency			*/
double	A0;		/* Lattice constant				*/
double	*E;		/* energy eigenvalues				*/
double	m_per_au;	/* unit conversion factor, m/a.u.		*/
int	N;		/* number of reciprocal lattice vectors		*/
int     n_min;          /* lowest output band				*/
int     n_max;          /* highest output band				*/
int	n_atoms;	/* number of atoms in (large) cell		*/
int	i;		/* loop index for matrix rows			*/
int	iE;		/* loop index for energy eigenvalues		*/
int	ik;		/* loop index for k-vectors			*/
int	j;		/* loop index for matrix columns		*/
char	filenameE[9];	/* character string for Energy output filename	*/
FILE	*Fk;		/* pointer to k.r file				*/
FILE	*FEk;		/* pointer to k.r file				*/
atom	*atoms;		/* the type and position of the atoms		*/
boolean	ev;		/* flag, if set output eigenvectors 		*/
vector	*G;		/* reciprocal lattice vectors			*/
vector	k;		/* electron wave vector				*/
vector	q;		/* G'-G						*/

/* default values	*/

A0=5.65e-10;
ev=false;
n_min=0;
n_max=3;
m_per_au=4*pi*epsilon_0*sqr(hbar/e_0)/m0;

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
	   printf("Usage:  pplb [-A lattice constant (\033[1m5.65\033[0mA)]\n");
	   printf("             [-n # lowest band \033[1m1\033[0m][-m highest band \033[1m4\033[0m], output states\n");
	   printf("             [-w output eigenvectors (wavefunctions) in range n->m]\n");
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

/* Allocate memory */

H_GG=(Complex *)calloc(N*N,sizeof(Complex));
 if(H_GG==0){fprintf(stderr,"Cannot allocate memory!\n");exit(0);}

V_GG=(Complex *)calloc(N,sizeof(Complex));
 if(V_GG==0){fprintf(stderr,"Cannot allocate memory!\n");exit(0);}

E=(double *)calloc(N,sizeof(double));
 if(E==0){fprintf(stderr,"Cannot allocate memory!\n");exit(0);}

ank=(Complex *)calloc(N*N,sizeof(Complex));
 if(ank==0){fprintf(stderr,"Cannot allocate memory!\n");exit(0);}

for(i=0;i<N;i++)        /* index down rows */
{
 for(j=0;j<i;j++)       /* index across cols, creates off diagonal elements */
 {
  q.x=((G+i)->x)-((G+j)->x);
  q.y=((G+i)->y)-((G+j)->y);
  q.z=((G+i)->z)-((G+j)->z);

  *(H_GG+i*N+j)=V(A0,m_per_au,atoms,n_atoms,q);    
 }
q.x=0;q.y=0;q.z=0;      /* creates diagonal potentials */

*(V_GG+i)=V(A0,m_per_au,atoms,n_atoms,q);
}

/* Add diagonal elements to matrix H_GG' */

ik=0;	/* Initialise k-loop index	*/
while((fscanf(Fk,"%lf %lf %lf",&k.x,&k.y,&k.z))!=EOF)
{
 k.x*=(2*pi/A0);k.y*=(2*pi/A0);k.z*=(2*pi/A0);
 for(i=0;i<N;i++)        /* add kinetic energy to diagonal elements */
 {
  T_GG.re=hbar*(sqr(((G+i)->x)+k.x)+sqr(((G+i)->y)+k.y)+sqr(((G+i)->z)+k.z))*hbar/(2*m0);
  T_GG.im=0;
  *(H_GG+i*N+i)=*(V_GG+i);
  ((H_GG+i*N+i)->re)+=T_GG.re;
 }

 /* Call NAG C routine to diagonalise matrix	*/

 f02axc(N,H_GG,N,E,ank,N,NAGERR_DEFAULT);

 /* Output eigenvalues in a separate file for each k point */

 sprintf(filenameE,"Ek%i.r",ik);
 FEk=fopen(filenameE,"w");
 for(iE=n_min;iE<=n_max;iE++)fprintf(FEk,"%10.6f\n",*(E+iE)/e_0);
 fclose(FEk);

 /* Output eigenvectors */

 if(ev)write_ank(ank,ik,N,n_min,n_max);

 ik++;	/* increment loop counter	*/
 
}/* end while*/
fclose(Fk);

free(G);
free(H_GG);
free(V_GG);
free(E);
free(ank);
free(atoms);

}/* end main */




atom
*read_atoms(n_atoms)

/* This function reads the atomic species (defined in the file as.r)
   into memory (addressed by the pointer as) and returns the start
   address of this block of memory and the number of lines	   */

int	*n_atoms;
{
 int    *as;
 int    ia=0;
 FILE 	*Fatoms;        /* file pointer to wavefunction file       */
 atom	*atoms;		/* atomic definitions			*/

 if((Fatoms=fopen("atoms.xyz","r"))==0)
 {
  fprintf(stderr,"Error: Cannot open input file 'atoms.xyz'!\n");
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



Complex 
V(A0,m_per_au,atoms,n_atoms,q)

double A0;	/* Lattice constant */
double m_per_au;/* conversion factor from SI to a.u. */
atom   *atoms;	/* atomic definitions	*/
int    n_atoms; /* number of atoms in structure */
vector q;       /* a reciprocal lattice vector, G'-G */
{
 extern double Vf();   /* Fourier transform of potential	*/
 Complex v;     /* potential					*/
 double	vf;	/* storage for returned data from Vf()		*/
 int	ia;	/* index across atoms				*/
 vector t;      /* general vertor representing atom within cell	*/

 v.re=0;v.im=0;

 for(ia=0;ia<n_atoms;ia++)
 {
  t.x=(atoms+ia)->r.x;
  t.y=(atoms+ia)->r.y;
  t.z=(atoms+ia)->r.z;
  vf=Vf(A0,m_per_au,q.x*q.x+q.y*q.y+q.z*q.z,(atoms+ia)->type);
  v.re+=cos(q.x*t.x+q.y*t.y+q.z*t.z)*vf;
  v.im-=sin(q.x*t.x+q.y*t.y+q.z*t.z)*vf;
 }

 v.re/=(double)(n_atoms/2);
 v.im/=(double)(n_atoms/2);

 return(v);
}



void
write_ank(ank,ik,N,n_min,n_max)

/* This function writes the eigenvectors (a_nk(G)) to the files ank.r */

Complex	*ank;
int	ik;		/* k point identifier				*/
int	N;
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
  fprintf(Fank,"%20.16le %20.16le ",(ank+iG*N+in)->re,(ank+iG*N+in)->im);
 fprintf(Fank,"\n");
}


fclose(Fank);
}

