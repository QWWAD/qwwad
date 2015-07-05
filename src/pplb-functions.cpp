#include "pplb-functions.h"

/** Writes the eigenvectors (a_nk(G)) to the files ank.r
 * \param[in] ank    Eigenvector coefficients
 * \param[in] ik     k point identifier
 * \param[in] N      Number of eigenvectors
 * \param[in] n_min  lowest output band
 * \param[in] n_max  highest output band
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
        for(in=n_min;in<=n_max && in < N;in++)
            fprintf(Fank,"%20.16le %20.16le ",ank(iG,in).real(), ank(iG,in).imag());
        fprintf(Fank,"\n");
    }

    fclose(Fank);
}

/**
 * \brief Get potential component of H_GG
 *
 * \param[in] A0       Lattice constant
 * \param[in] m_per_au conversion factor from SI to a.u.
 * \param[in] atoms    atomic definitions
 * \param[in] n_atoms  number of atoms in structure
 * \param[in] q        a reciprocal lattice vector, G'-G
 */
std::complex<double> V(double           A0,
                       double           m_per_au,
                       atom            *atoms,
                       size_t           n_atoms,
                       arma::vec const &q)
{
    std::complex<double> v = 0.0; // potential
    const double q_dot_q = dot(q,q);

    // Loop over all atoms in the set and add contribution from each
    for(unsigned int ia=0; ia<n_atoms; ++ia)
    {
        const double q_dot_t = dot(q, atoms[ia].r);
        const double vf = Vf(A0,m_per_au,q_dot_q,atoms[ia].type);
        v += exp(std::complex<double>(0.0,-q_dot_t)) * vf; // [QWWAD3, 15.76]
    }

    v *= 2.0/n_atoms;

    return v;
}
