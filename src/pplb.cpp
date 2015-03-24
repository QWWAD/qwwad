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
  
   Modifications, December 1999
   	Use of LAPACK diagonalisation routine
								*/

#if HAVE_CONFIG_H
# include <config.h>
#endif

#include <complex>
#include <valarray>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <gsl/gsl_math.h>

#include <armadillo>

#include "struct.h"
#include "maths.h"
#include "qclsim-constants.h"
#include "qclsim-linalg.h"

#include "ppff.h"

typedef struct
{
 char	type[12];
 vector	r;
}atom;

atom   * read_atoms(size_t *n_atoms);
vector * read_rlv(double A0, size_t *N);

std::complex<double> V(double     A0,
                       double     m_per_au,
                       atom      *atoms,
                       size_t     n_atoms,
                       arma::vec &q);
void
write_ank(arma::cx_mat &ank,
          int           ik,
          int           N,
          int           n_min,
          int           n_max);

int main(int argc,char *argv[])
{
double	A0;		/* Lattice constant				*/
double	m_per_au;	/* unit conversion factor, m/a.u.		*/
size_t	N;		/* number of reciprocal lattice vectors		*/
int     n_min;          /* lowest output band				*/
int     n_max;          /* highest output band				*/
int	iE;		/* loop index for energy eigenvalues		*/
int	ik;		/* loop index for k-vectors			*/
char	filenameE[9];	/* character string for Energy output filename	*/
FILE	*Fk;		/* pointer to k.r file				*/
FILE	*FEk;		/* pointer to Ek.r file				*/
atom	*atoms;		/* the type and position of the atoms		*/
bool	ev;		/* flag, if set output eigenvectors 		*/
vector	*G;		/* reciprocal lattice vectors			*/
vector	k;		/* electron wave vector				*/

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

size_t	n_atoms;	/* number of atoms in (large) cell		*/
atoms=read_atoms(&n_atoms);		/* read in atomic basis	*/

G=read_rlv(A0,&N);	/* read in reciprocal lattice vectors	*/

// components of H_G'G (see notes)
arma::cx_mat H_GG(N,N);

// potential energy of H_GG, diagonal elements
std::valarray< std::complex<double> > V_GG(N);

for(unsigned int i=0;i<N;i++)        /* index down rows */
{
    for(unsigned int j=0;j<N;j++)       /* index across cols, creates off diagonal elements */
    {
        arma::vec q(3);		/* G'-G						*/
        q(0) = G[i].x - G[j].x;
        q(1) = G[i].y - G[j].y;
        q(2) = G[i].z - G[j].z;

        H_GG(i,j) = V(A0,m_per_au,atoms,n_atoms,q);    
    }

    V_GG[i] = H_GG(i,i); // Record potentials from diagonal
}

/* Add diagonal elements to matrix H_GG' */

ik=0;	/* Initialise k-loop index	*/
while((fscanf(Fk,"%lf %lf %lf",&k.x,&k.y,&k.z))!=EOF)
{
 k.x*=(2*pi/A0);k.y*=(2*pi/A0);k.z*=(2*pi/A0);
 for(unsigned int i=0;i<N;i++)        /* add kinetic energy to diagonal elements */
 {
     // kinetic energy component of H_GG [QWWAD3, 15.91]
     std::complex<double> T_GG=hBar*hBar/(2*me) * (gsl_pow_2(G[i].x+k.x)+gsl_pow_2(G[i].y+k.y)+gsl_pow_2(G[i].z+k.z));
     H_GG(i,i) = T_GG + V_GG[i];
 }

 // Find the eigenvalues & eigenvectors of the Hamiltonian matrix
 arma::vec E(N); // Energy eigenvalues
 arma::cx_mat ank(N,N); // coefficients of eigenvectors
 arma::eig_sym(E, ank, H_GG);

 /* Output eigenvalues in a separate file for each k point */

 sprintf(filenameE,"Ek%i.r",ik);
 FEk=fopen(filenameE,"w");
 for(iE=n_min;iE<=n_max;iE++)fprintf(FEk,"%10.6f\n",E(iE)/e);
 fclose(FEk);

 /* Output eigenvectors */

 if(ev){
	write_ank(ank,ik,N,n_min,n_max);
        }

 ik++;	/* increment loop counter	*/
}/* end while*/
fclose(Fk);

free(G);
free(atoms);

return EXIT_SUCCESS;
}/* end main */

/* This function reads the atomic species (defined in the file as.r)
   into memory (addressed by the pointer as) and returns the start
   address of this block of memory and the number of lines	   */
atom *read_atoms(size_t *n_atoms)
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

 int n_read = fscanf(Fatoms,"%lu",n_atoms);
 
 /* Allocate memory for atom definitions	*/
 if (n_read == 1)
   atoms=(atom *)calloc(*n_atoms,sizeof(atom));
 else
 {
   fprintf(stderr, "Could not read number of atoms!\n");
   exit(EXIT_FAILURE);
 }

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

/* This function reads the reciprocal lattice vectors (defined in
   the file G.r) into the array G[] and then converts into SI units */
vector * read_rlv(double A0, size_t *N)
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

/**
 * \param A0	   Lattice constant
 * \param m_per_au conversion factor from SI to a.u.
 * \param atoms    atomic definitions
 * \param n_atoms  number of atoms in structure
 * \param q        a reciprocal lattice vector, G'-G
 */
std::complex<double> V(double     A0,
                       double     m_per_au,
                       atom      *atoms,
                       size_t     n_atoms,
                       arma::vec &q)
{
    std::complex<double> v = 0.0;     /* potential					*/

    for(unsigned int ia=0; ia<n_atoms; ++ia)
    {
        arma::vec t(3);      /* general vector representing atom within cell	*/
        t(0) = atoms[ia].r.x;
        t(1) = atoms[ia].r.y;
        t(2) = atoms[ia].r.z;
        const double q_dot_q = dot(q,q);
        const double q_dot_t = dot(q,t);
        const double vf = Vf(A0,m_per_au,q_dot_q,atoms[ia].type);
        v += exp(std::complex<double>(0.0,-q_dot_t)) * vf; // Add contribution to potential from this atom [QWWAD3, 15.92]
    }

    v *= 2.0/n_atoms;

    return v;
}

/** This function writes the eigenvectors (a_nk(G)) to the files ank.r
 * \param double	ank
 * \param ik            k point identifier
 * \param N
 * \param n_min         lowest output band
 * \param n_max         highest output band
 */
void
write_ank(arma::cx_mat &ank,
          int           ik,
          int           N,
          int           n_min,
          int           n_max)
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
  fprintf(Fank,"%20.16le %20.16le ",ank(iG,in).real(), ank[iG*N+in].imag());
 fprintf(Fank,"\n");
}


fclose(Fank);
}
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :
