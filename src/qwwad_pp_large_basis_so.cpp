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

#if HAVE_CONFIG_H
# include <config.h>
#endif

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <complex>
#include "struct.h"
#include "maths.h"
#include "qwwad/constants.h"
#include "qwwad/file-io.h"
#include "qwwad/options.h"

#include "ppff.h"	/* the PseudoPotential Form Factors	*/
#include "pplb-functions.h"
#include "ppsop.h"	/* the Spin-Orbit Parameters		*/

using namespace QWWAD;
using namespace constants;

std::complex<double>
Vso(const std::vector<atom>      &atoms,
    const std::vector<arma::vec> &G,
    arma::vec const              &k,
    const unsigned int            i,
    const unsigned int            j,
    const size_t                  N);

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

    std::string filename("atoms.xyz");
    const auto atoms = read_atoms(filename.c_str()); // read in atomic basis

    const auto G  = read_rlv(A0); // read in reciprocal lattice vectors
    const auto N  = G.size(); // number of reciprocal lattice vectors
    const auto Ns = 2*N;      // order of H_GG with spin (2*N)

    const auto m_per_au = 4.0*pi*eps0*hBar*hBar/(e*e*me); // Unit conversion factor, m/a.u

    // Compute crystal potential matrix. Note that this is independent of wave-vector
    // so we only need to do this once.
    arma::cx_mat V_GG(Ns,Ns);

    for(unsigned int i=0;i<N;i++) // index down rows within Block 1 and 2 of matrix
    {
        // Create upper triangle of Block 1
        for(unsigned int j=i; j<N; j++)
        {
            const auto q = G[i] - G[j];
            V_GG(i,j) = V(A0,m_per_au,atoms,q);

            // Copy elements to upper triangle of all other blocks
            V_GG(i+N, j) = V_GG(i, j+N) = V_GG(i+N, j+N) = V_GG(i,j);

            // Since matrix is Hermitian, fill in the lower triangles in all 4 blocks
            V_GG(j, i) = V_GG(j+N, i) = V_GG(j, i+N) = V_GG(j+N, i+N) = conj(V_GG(i,j));
        }
    }

    /* Add k-dependent elements to matrix H_GG' */
    for(unsigned int ik = 0; ik < nk; ++ik)
    {
        if(opt.get_verbose())
            std::cout << "Calculating energy at k = " << std::endl
                << k[ik] << " (" << ik + 1 << "/" << nk << ")" << std::endl;

        auto H_GG = V_GG; // Complete Hamiltonian matrix

        for(unsigned int i=0;i<N;i++)        /* add kinetic energy to diagonal elements */
        {
            // kinetic energy component of H_GG [QWWAD3, 15.77]
            arma::vec G_plus_k = G[i] + k[ik];
            const double G_plus_k_sq = dot(G_plus_k, G_plus_k);
            std::complex<double> T_GG=hBar*hBar/(2*me) * G_plus_k_sq;
            H_GG(i, i) += T_GG; // Block 1
            H_GG(i+N, i+N) += T_GG; // Block 4
        }

        /* Add spin-orbit components to lower triangle	*/
        for(unsigned int i=0;i<Ns;i++)
        {
            for(unsigned int j=0;j<=i;j++)	
            {
                H_GG(i,j) += Vso(atoms,G,k[ik],i,j,N);
            }
        }

        // Find the eigenvalues & eigenvectors of the Hamiltonian matrix
        arma::vec E(N); // Energy eigenvalues
        arma::cx_mat ank(Ns,Ns); // coefficients of eigenvectors
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
    }

    return EXIT_SUCCESS;
}/* end main */

/**
 * \brief Spin-orbit component of H_GG
 *
 * \param[in] atoms   atomic definitions
 * \param[in] G	      reciprocal lattice vectors
 * \param[in] k       the electron momentum
 * \param[in] i	      index
 * \param[in] j       index
 * \param[in] N       number of reciprocal lattice vectors
 */
std::complex<double>
Vso(const std::vector<atom>      &atoms,
    const std::vector<arma::vec> &G,
    arma::vec const              &k,
    const unsigned int            i,
    const unsigned int            j,
    const size_t                  N)
{
    std::complex<double> vso = 0; // intermediate value of spin-orbit interaction
    arma::vec q;	/* a reciprocal lattice vector, G'-G 		*/

    if((i<N)&&(j<N))	/* Block 1	*/
    {   
        const arma::vec A = cross(G[i] + k, G[j] + k);
        vso = A(2);
        q = G[i] - G[j];
    }
    if((i<N)&&(j>=N))	/* Block 2	*/
    {
        const arma::vec A = cross(G[i] + k, G[j-N] + k);
        vso = std::complex<double>(A(0), -A(1));
        q = G[i] - G[j-N];
    }
    if((i>=N)&&(j<N))	/* Block 3	*/
    {
        const arma::vec A = cross(G[i-N] + k, G[j] + k);
        vso = std::complex<double>(A(0), A(1));
        q = G[i-N] - G[j];
    }
    if((i>=N)&&(j>=N))	/* Block 4	*/
    {
        const arma::vec A = cross(G[i-N] + k, G[j-N] + k);
        vso = -A(2);
        q = G[i-N] - G[j-N];
    }	
    vso *= std::complex<double>(0,-1); /* -i(G'+k)x(G+k).sigma	*/

    // Potential term
    std::complex<double> v=0;				/* Initialise for sum	*/

    for(auto const atom : atoms)
    {
        const double q_dot_t = dot(q, atom.r);
        const auto Lambda = lambda(atom.type);
        v += Lambda * exp(std::complex<double>(0.0,-q_dot_t)); // [QWWAD3, 15.81]
    }

    v *= vso;
    v/=atoms.size();

    return v;
}
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :
