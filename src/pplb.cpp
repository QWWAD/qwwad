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
#include "qwwad/constants.h"
#include "qwwad/linear-algebra.h"
#include "qwwad/options.h"

#include "pplb-functions.h"
#include "ppff.h"

Options configure_options(int argc, char* argv[])
{
    Options opt;

    std::string doc("Large-basis pseudopotential calculation for user-defined cell");

    opt.add_option<double>("latticeconst,A", 5.65, "Lattice constant [angstrom]");
    opt.add_option<size_t>("nmin,n",            4, "Lowest output band index (VB = 4, CB = 5)");
    opt.add_option<size_t>("nmax,m",            5, "Highest output band index (VB = 4, CB = 5)");
    opt.add_option<bool>  ("printev,w",            "Print eigenvectors to file");

    opt.add_prog_specific_options_and_parse(argc, argv, doc);

    return opt;
}

int main(int argc,char *argv[])
{
    const auto opt = configure_options(argc, argv);

    const auto A0    = opt.get_option<double>("latticeconst") * 1e-10; // Lattice constant [m]
    const auto n_min = opt.get_option<size_t>("nmin")-1;               // Lowest output band
    const auto n_max = opt.get_option<size_t>("nmax")-1;               // Highest output band
    const auto ev    = opt.get_option<bool>  ("printev");              // Print eigenvectors?

    // Read desired wave vector points from file
    std::valarray<double> kx;
    std::valarray<double> ky;
    std::valarray<double> kz;
    read_table("k.r", kx, ky, kz);
    size_t nk = kx.size(); // Number of wave vector samples to compute

    // Copy wave vector components in each direction into a list of wave vectors 
    std::vector<arma::vec> k(nk, arma::vec(3));

    for(unsigned int ik = 0; ik < nk; ++ik)
    {
        k[ik](0) = kx[ik];
        k[ik](1) = ky[ik];
        k[ik](2) = kz[ik];
        k[ik] *= 2.0*pi/A0;
    }

    size_t	n_atoms;	/* number of atoms in (large) cell		*/
    std::string filename("atoms.xyz");
    atom *atoms = read_atoms(&n_atoms, filename.c_str());		/* read in atomic basis	*/

    const auto G = read_rlv(A0); // read in reciprocal lattice vectors
    const auto N = G.size(); // number of reciprocal lattice vectors

    const auto m_per_au = 4.0*pi*eps0*hBar*hBar/(e*e*me); // Unit conversion factor, m/a.u

    // Compute crystal potential matrix. Note that this is independent of wave-vector
    // so we only need to do this once.
    arma::cx_mat V_GG(N,N);

    for(unsigned int i=0;i<N;i++)        /* index down rows */
    {
        // Fill in the upper triangle of the matrix
        for(unsigned int j=i;j<N;j++)
        {
            const auto q = G[i] - G[j];
            V_GG(i,j) = V(A0,m_per_au,atoms,n_atoms,q);

            // Fill in the lower triangle by taking the Hermitian transpose of the elements
            V_GG(j,i) = conj(V_GG(i,j));
        }
    }

    /* Add diagonal elements to matrix H_GG' */
    for(unsigned int ik = 0; ik < nk; ++ik)
    {
        if(opt.get_verbose())
            std::cout << "Calculating energy at k = " << std::endl
                << k[ik] << " (" << ik + 1 << "/" << nk << ")" << std::endl;

        // Construct the complete Hamiltonian matrix now, using crystal potential and
        // kinetic energy on the diagonals
        auto H_GG = V_GG;

        for(unsigned int i=0;i<N;i++)
        {
            // kinetic energy component of H_GG [QWWAD3, 15.77]
            arma::vec G_plus_k = G[i] + k[ik];
            const double G_plus_k_sq = dot(G_plus_k, G_plus_k);
            std::complex<double> T_GG=hBar*hBar/(2*me) * G_plus_k_sq;
            H_GG(i,i) += T_GG;
        }

        // Find the eigenvalues & eigenvectors of the Hamiltonian matrix
        arma::vec E(N); // Energy eigenvalues
        arma::cx_mat ank(N,N); // coefficients of eigenvectors
        arma::eig_sym(E, ank, H_GG);

        /* Output eigenvalues in a separate file for each k point */
        char	filenameE[9];	/* character string for Energy output filename	*/
        sprintf(filenameE,"Ek%i.r",ik);
        FILE *FEk=fopen(filenameE,"w");

        for(auto iE=n_min; iE<=n_max; iE++)
            fprintf(FEk,"%10.6f\n",E(iE)/e);

        fclose(FEk);

        /* Output eigenvectors */

        if(ev){
            write_ank(ank,ik,N,n_min,n_max);
        }
    }/* end while*/

    free(atoms);

    return EXIT_SUCCESS;
}/* end main */

// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :
