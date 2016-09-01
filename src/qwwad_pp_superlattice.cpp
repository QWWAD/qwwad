/**
 * \file   qwwad_pp_superlattice.cpp
 * \brief  PseudoPotential SuperLattice calculation
 * \author Paul Harrison  <p.harrison@shu.ac.uk>
 * \author Alex Valavanis <a.valavanis@leeds.ac.uk>
 *
 * \details
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
*/

#include <cstdio>
#include <cmath>
#include <cstdlib>
#include <complex>
#include "struct.h"
#include "maths.h"
#include "qwwad/constants.h"
#include "qwwad/file-io.h"
#include "qwwad/linear-algebra.h"

#include "ppff.h"

using namespace QWWAD;
using namespace constants;

static std::complex<double> i1(0,1);

static std::complex<double>
V(double           A0,
  double           m_per_au,
  std::vector<atom> const &atoms,
  std::vector<atom> const &atomsp,
  arma::vec const &g);

static std::complex<double>
VF(double           A0,
   double           F,
   double           q,
   std::vector<atom> const &atoms,
   arma::vec const &g);

static std::valarray<std::complex<double>> read_ank(int N,
                                                    int Nn,
                                                    int Nkxi);

static double * read_Enk(int Nn,
                         int Nkxi);

static int read_ank0(int N);

static std::vector<arma::vec> read_kxi(double  A0);

static void write_VF(double  A0,
                     double  F,
                     double  q,
                     std::vector<atom> const &atoms);

static void clean_Hdash(std::complex<double> *Hdash,
                        int                   order);

int main(int argc,char *argv[])
{
std::complex<double>	*WORK;		/* LAPACK: workspace				*/
double	A0;		/* Lattice constant				*/
double	*Enk;		/* bulk energy eigenvalues			*/
double	*Exi;		/* energy eigenvalues				*/
double	F;		/* the electric field strength			*/
double	m_per_au;	/* unit conversion factor, m/a.u.		*/
double	q;		/* the carrier charge, e=-e_0, h=+e_0		*/
double	*RWORK;		/* LAPACK: workspace				*/
int	INFO;		/* LAPACK: information integer			*/
int	Nn;		/* number of bands in calculation		*/
int	Nkxi;		/* number of k points in calculation		*/
int     n_min;          /* lowest output band                           */
int     n_max;          /* highest output band                          */
int	i;		/* loop index for matrix rows			*/
int	iE;		/* loop index for energy eigenvalues		*/
unsigned int	iG;		/* index over G					*/
unsigned int	iGdash;		/* index over G'				*/
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
auto const atoms=read_atoms(filename);	/* read in atomic basis	*/
strcpy(filename,"atomsp.xyz");
auto const atomsp=read_atoms(filename);	/* read in perturbation	*/

auto G = read_rlv(A0); // read in reciprocal lattice vectors
auto N = G.size(); // number of reciprocal lattice vectors
Nn=read_ank0(N);	/* reads a single ank.r file just to 
			   deduce the number of bands Nn	*/

auto kxi = read_kxi(A0); // read in set of kxi points
Nkxi = kxi.size();

auto ank=read_ank(N,Nn,Nkxi);/* read in bulk eigenvectors		*/
Enk=read_Enk(Nn,Nkxi);	/* read in bulk eigenvalues		*/

/* Output the  Fourier Transform of the electric field	*/

if(o) write_VF(A0,F,q,atoms);

/* Allocate memory */

LWORK=2*Nn*Nkxi;		/* LAPACK workspace definition		*/

WORK=(std::complex<double> *)calloc(LWORK,sizeof(std::complex<double>));
 if(WORK==0){fprintf(stderr,"Cannot allocate memory!\n");exit(0);}

RWORK=(double *)calloc(3*Nn*Nkxi-2,sizeof(double));
 if(RWORK==0){fprintf(stderr,"Cannot allocate memory!\n");exit(0);}

/* end of LAPACK definitions	*/
auto Hdash=(std::complex<double> *)calloc(Nn*Nkxi*Nn*Nkxi,sizeof(std::complex<double>));
 if(Hdash==0){fprintf(stderr,"Cannot allocate memory!\n");exit(0);}

Exi=(double *)calloc(Nn*Nkxi,sizeof(double));
 if(Exi==0){fprintf(stderr,"Cannot allocate memory!\n");exit(0);}

// coefficients of eigenvectors
auto Ank= std::valarray<std::complex<double>>(Nn*Nkxi*Nn*Nkxi);

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
       // Calculate appropriate g vector [QWWAD4, 16.38]
       auto g = G[iGdash] - G[iG] + kxi[ikxidash] - kxi[ikxi];

       /* Add on potential term	*/

       Hdash[i*Nn*Nkxi+j] +=
           (
            ( conj(   ank[ikxidash*N*Nn+iGdash*Nn+indash]   )

              *

              ank[ikxi*N*Nn+iG*Nn+in]
            )
            *
            (
             V(A0,m_per_au,atoms,atomsp,g) +
             VF(A0,F,q,atoms,g)
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
  Ank[i*Nn*Nkxi+j]=Hdash[i*Nn*Nkxi+j]; // to eigenvector space

/* call LAPACK diagonalisation routine	*/

OH=Nn*Nkxi;			/* the order of Hdash	*/
ctranspose(Ank); //,Nn*Nkxi);	/* because of FORTRAN LAPACK	*/
zheev_(&JOBZ,&UPLO,&OH,&Ank[0],&OH,Exi,WORK,&LWORK,RWORK,&INFO);

/* Output eigenvalues in a separate file for each k point */

FExi=fopen("Exi.r","w");
for(iE=n_min;iE<=n_max;iE++)fprintf(FExi,"%10.6f\n",*(Exi+iE)/e);
fclose(FExi);

free(Hdash);
free(Exi);
free(Enk);

return EXIT_SUCCESS;
}/* end main */




/**
 * \brief removes round-up errors => makes H' Hermitian
 *
 * \param Hdash The matrix to be `cleaned'
 * \param N 	the order of the matrix
 *
 * \details This function accounts for any computational errors in the diagonal of
 *          H', it ensures that the imaginary components of these elements are zero	
 */
static void clean_Hdash(std::complex<double> *Hdash,
                        int                   N)
{
 for(int i=0; i<N; i++)
 {
  if(Hdash[i*N+i].imag() / Hdash[i*N+i].real() < 1e-10){
      Hdash[i*N+i] = Hdash[i*N+i].real();
  }
 }

}

/**
 * \brief Reads the first of the ank.r files, i.e., ank0.r, just
 *        with the purpose of deducing the number of bands included in the
 *        calculation
 *
 * \param N The number of terms in each eigenvector
 */
static int read_ank0(int N)
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
while(fscanf(Fank,"%*f %*f")!=EOF)
 n++;

/* The number of bands Nn is therefore the total number of elements divided
   by the number of terms in each eigenvector				*/

Nn=n/N;

return(Nn);

}

/**
 * \brief Reads the eigenvectors (a_nk(G)) from the file a_nk.r
 *        created by the code pplb.c into the array a_nk[N][Nn]
 *
 * \param N    The number of terms in each eigenvector
 * \param Nn   The number of bands in file
 * \param Nkxi number of k points in calculation
 */
static std::valarray<std::complex<double>> read_ank(int N,
                                                    int Nn,
                                                    int Nkxi)
{
 char	filename[9];	/* eigenfunction output filename                */
 FILE   *Fank;		/* file pointer to eigenvectors file		*/

 /* Allocate memory for eigenvectors	*/
 std::valarray<std::complex<double> > ank(N*Nn*Nkxi);

 /* Finally read eigenvectors into structure	*/
 for(int ikxi=0;ikxi<Nkxi;ikxi++)
 {
     /* Open the bulk eigenvector file at each of the bulk k (kxi) points	*/
     sprintf(filename,"ank%i.r",ikxi);
     if((Fank=fopen(filename,"r"))==0)
     {fprintf(stderr,"Error: Cannot open input file 'ank%i.r'!\n",ikxi);exit(0);}

     /* Read in data	*/

     // Loop over all G vectors
     for(int iG=0;iG<N;iG++)
     {
         // Loop over all bands
         for(int in=0;in<Nn;in++)
         {
             double temp_re=0.0;
             double temp_im=0.0;
             int n_read = fscanf(Fank,"%lf %lf", &temp_re, &temp_im);
             if (n_read == 2)
             {
                 ank[ikxi*N*Nn+iG*Nn+in] = std::complex<double>(temp_re, temp_im);
             }
             else
             {
                 fprintf(stderr, "Could not read number.\n");
                 exit(EXIT_FAILURE);
             }
         }
     }

     fclose(Fank);
 }

return(ank);
}


/**
 * \brief Reads the eigenvalues (E_nk) from the files Ek?.r
 *        created by the code pplb.c
 *
 * \param Nn   The number of bands in file
 * \param Nkxi number of k points in calculation
 */
static double * read_Enk(int Nn,
                         int Nkxi)
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
  int n_read = fscanf(FEnk,"%lf",(Enk+ikxi*Nn+in));
  if (n_read == 1)
      Enk[ikxi*Nn+in] *= e;		/* convert from eV->S.I.	*/
  else
  {
      fprintf(stderr, "Could not read number\n");
      exit(EXIT_FAILURE);
  }
 }

 fclose(FEnk);
}

return(Enk);
}

/**
 * \brief Reads the set of electron momenta vectors kxi (defined
 *        in the file k.r---as needed to produce the bulk band structure) into the 
 *        allocated memory space accessed via the pointer kxi, then converting them
 *        from units of (2*pi/A0) into S.I.
 *
 * \param[in] A0 lattice constant
 *
 * \return array of momenta
 */
static std::vector<arma::vec> read_kxi(double  A0)
{
    std::valarray<double> x;
    std::valarray<double> y;
    std::valarray<double> z;
    read_table("k.r", x, y, z);

    std::vector<arma::vec> k;
    size_t N = x.size();
    for(unsigned int i = 0; i < N; ++i)
    {
        arma::vec temp(3);
        temp(0) = x[i];
        temp(1) = y[i];
        temp(2) = z[i];
        temp *= 2*pi/A0;
        k.push_back(temp);
    }

    return k;
}

/**
 * \brief Potential component of Hdash
 *
 * \param A0	   Lattice constant
 * \param m_per_au conversion factor from SI to a.u.
 * \param atoms    atomic definitions
 * \param atomsp   atomic definitions
 * \param g        the vector, g=G'-G+kxi'-kxi
 */
static std::complex<double>
V(double           A0,
  double           m_per_au,
  std::vector<atom> const &atoms,
  std::vector<atom> const &atomsp,
  arma::vec const &g)
{
 std::complex<double> v;     /* potential					*/
 double	vf;	/* storage for returned data from Vf()		*/
 double	vfdash;	/* storage for returned data from Vf()		*/
 v=0;

 // Sum over atoms [QWWAD4, 16.50]
 auto const n_atoms = atoms.size();
 for(unsigned int ia=0;ia<n_atoms;ia++)
 {
     // general vector representing atom within cell
     auto const t = atoms[ia].r;
     auto const g_dot_g = arma::dot(g,g);
     auto const g_dot_t = arma::dot(g,t);

     vf=Vf(A0, m_per_au,   g_dot_g, atoms[ia].type);
     vfdash=Vf(A0,m_per_au,g_dot_g, atomsp[ia].type);
     v += exp(-i1*g_dot_t)*(vfdash-vf);
 }

 /* These last divisions represent Omega_c/Omega_sl included here for
    convenience	*/
 v/=(double)(n_atoms/2);

 return(v);
}

/**
 * \brief field potential component of Hdash
 *
 * \param A0      Lattice constant
 * \param F       the electric field strength
 * \param q       the carrier charge
 * \param atoms   atomic definitions
 * \param g       the vector, g=G'-G+kxi'-kxi
 */
static std::complex<double>
VF(double           A0,
   double           F,
   double           q,
   std::vector<atom> const &atoms,
   arma::vec const &g)
{
    // Calculate the number of lattice constants in each superlattice period
    int const n_z=atoms.size()/4;

    /* The smallest reciprocal lattice vector is therefore pi/(n_z*A0), thus 
       anything less than this must be zero				*/
    double const zero=pi/(n_z*A0)/10;

    std::complex<double> v=0.0; // potential

    // Imitates the behaviour of the product of the two delta-functions
    auto const g_x = g(0);
    auto const g_y = g(1);

    if(fabs(g_x) > zero ||
       fabs(g_y) > zero)
    {
        v = 0.0;
    }
    else
    {
        // Get midpoint of unit cell
        const auto z_first = atoms.front().r(2); // z-position of first atom
        const auto z_last  = atoms.back().r(2); // z-position of last atom
        double const z0=(z_last + z_first)/2;

        std::complex<double> v1;	/* first term in potential			*/
        std::complex<double> v2;	/* second term in potential			*/

        const auto g_z = g(2); // z-component of wave vector

        if(fabs(g_z)>zero)
        {
            // Evaluate the first and second terms in the integral [QWWAD4, 16.66]
            v1 = std::complex<double>(1.0/(g_z*g_z), (z_last-z0)/g_z);
            v1 *= exp(-i1 * g_z * z_last);

            v2 = std::complex<double>(1.0/(g_z*g_z), (z_first-z0)/g_z);
            v2 *= exp(-i1 * g_z * z_first);
        }
        else
        {
            // For case where g_z = 0, avoid singularity using QWWAD4, 16.70
            v1 = z_last *z_last /2 - z_last  * z0;
            v2 = z_first*z_first/2 - z_first * z0;
        }

        // Combine the two
        v=v1-v2;

        // Multiply through by the scaling factor
        v*=-q*F/((double)(n_z)*A0);
    }

    return v;
}

/**
 * \brief writes the FT of the electric field potential
 *
 * \param A0      Lattice constant
 * \param F       the electric field strength
 * \param q       the carrier charge
 * \param atoms   atomic definitions
 */
static void write_VF(double  A0,
                     double  F,
                     double  q,
                     std::vector<atom> const &atoms)
{
 FILE	*FVFg;		/* pointer to output file	*/
 int	i;		/* index			*/

 FVFg=fopen("VFg.r","w");

 for(i=-500;i<500;i++)
 {
     arma::vec g(3); // The vector, g=G'-G+kxi'-kxi
     g(0)=0;
     g(1)=0;
     g(2)=((float)i*4/500)*pi/A0;

     fprintf(FVFg,"%le %le\n",g(2)/(pi/A0),abs(VF(A0,F,q,atoms,g)));
 }

 fclose(FVFg);

}
