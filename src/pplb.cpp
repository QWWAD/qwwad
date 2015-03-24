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
#include "qwwad-options.h"

#include "ppff.h"

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

Options configure_options(int argc, char* argv[])
{
    Options opt;

    std::string doc("Large-basis pseudopotential calculation for user-defined cell");

    opt.add_numeric_option("latticeconst,A", 5.65, "Lattice constant [angstrom]");
    opt.add_size_option   ("nmin,n",            4, "Lowest output band index (VB = 4, CB = 5)");
    opt.add_size_option   ("nmax,m",            5, "Highest output band index (VB = 4, CB = 5)");
    opt.add_switch        ("printev,w",            "Print eigenvectors to file");

    opt.add_prog_specific_options_and_parse(argc, argv, doc);

    return opt;
}

int main(int argc,char *argv[])
{
    const Options opt = configure_options(argc, argv);

    const double A0    = opt.get_numeric_option("latticeconst") * 1e-10; // Lattice constant [m]
    const int    n_min = opt.get_size_option("nmin"); // Lowest output band
    const int    n_max = opt.get_size_option("nmax"); // Highest output band
    const bool   ev    = opt.get_switch("printev");   // Print eigenvectors?

// Read desired wave vector points from file
std::valarray<double> kx;
std::valarray<double> ky;
std::valarray<double> kz;
read_table("k.r", kx, ky, kz);
std::vector<arma::vec> k;
size_t nk = k.size();

for(unsigned int ik = 0; ik < nk; ++ik)
{
    arma::vec _k(3);
    _k(0) = kx[ik];
    _k(1) = ky[ik];
    _k(2) = kz[ik];
    _k *= 2.0*pi/A0;
    k.push_back(_k);
}

size_t	n_atoms;	/* number of atoms in (large) cell		*/
std::string filename("atoms.xyz");
atom *atoms = read_atoms(&n_atoms, filename.c_str());		/* read in atomic basis	*/

std::vector<arma::vec> G=read_rlv(A0);	/* read in reciprocal lattice vectors	*/
size_t	N = G.size();		/* number of reciprocal lattice vectors		*/

// components of H_G'G (see notes)
arma::cx_mat H_GG(N,N);

// potential energy of H_GG, diagonal elements
std::valarray< std::complex<double> > V_GG(N);

double m_per_au=4*pi*eps0*gsl_pow_2(hBar/e)/me; // Unit conversion factor, m/a.u
for(unsigned int i=0;i<N;i++)        /* index down rows */
{
    for(unsigned int j=0;j<N;j++)       /* index across cols, creates off diagonal elements */
    {
        arma::vec q = G[i] - G[j];		/* G'-G						*/
        H_GG(i,j) = V(A0,m_per_au,atoms,n_atoms,q);    
    }

    V_GG[i] = H_GG(i,i); // Record potentials from diagonal
}

/* Add diagonal elements to matrix H_GG' */
for(unsigned int ik = 0; ik < nk; ++ik)
{
 for(unsigned int i=0;i<N;i++)        /* add kinetic energy to diagonal elements */
 {
     // kinetic energy component of H_GG [QWWAD3, 15.91]
     arma::vec G_plus_k = G[i] + k[ik];
     const double G_plus_k_sq = dot(G_plus_k, G_plus_k);
     std::complex<double> T_GG=hBar*hBar/(2*me) * G_plus_k_sq;
     H_GG(i,i) = T_GG + V_GG[i];
 }

 // Find the eigenvalues & eigenvectors of the Hamiltonian matrix
 arma::vec E(N); // Energy eigenvalues
 arma::cx_mat ank(N,N); // coefficients of eigenvectors
 arma::eig_sym(E, ank, H_GG);

 /* Output eigenvalues in a separate file for each k point */
 char	filenameE[9];	/* character string for Energy output filename	*/
 sprintf(filenameE,"Ek%i.r",ik);
 FILE *FEk=fopen(filenameE,"w");
 for(int iE=n_min; iE<=n_max; iE++) fprintf(FEk,"%10.6f\n",E(iE)/e);
 fclose(FEk);

 /* Output eigenvectors */

 if(ev){
	write_ank(ank,ik,N,n_min,n_max);
        }

 ik++;	/* increment loop counter	*/
}/* end while*/

free(atoms);

return EXIT_SUCCESS;
}/* end main */


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
        const double q_dot_q = dot(q,q);
        const double q_dot_t = dot(q, atoms[ia].r);
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
     fprintf(Fank,"%20.16le %20.16le ",ank(iG,in).real(), ank(iG,in).imag());
 fprintf(Fank,"\n");
}

fclose(Fank);
}
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :
