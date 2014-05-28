/**
 * \file    efshoot.cpp
 * \brief   Solve Schroedinger's equation using shooting method
 * \author  Paul Harrison  <p.harrison@shu.ac.uk>
 * \author  Alex Valavanis <a.valavanis@leeds.ac.uk>
 *
 * \details This program uses a shooting technique to calculate the
   uncorrelated one particle energies of any user supplied
   potential.  The potential is read from the file v.r

   Paul Harrison, July 1992                  

   The program has been updated to be run entirely from the command 
   line and stripped down to calculate the energies only.  It now also
   includes support for non-parabolicity.

   Paul Harrison, December 1996
 */

#include <iostream>
#include <cstdlib>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>
#include "qclsim-constants.h"
#include "qclsim-linalg.h"
#include "qwwad-options.h"

using namespace Leeds;
using namespace constants;

struct shoot_params
{
    std::valarray<double>       &z;       ///< Spatial sampling points [m]
    const std::valarray<double> &V;       ///< Potential profile [J]
    const std::valarray<double> &m0;      ///< Band-edge effective mass at each point [kg]
    const std::valarray<double> &alpha;   ///< Nonparabolicity parameter at each point [1/J]
    const bool                   np_flag; ///< True if nonparabolicity is to be used
};

double psi_at_inf(double  E,
                  void   *params); 

double shoot_wavefunction(std::valarray<double>       &psi,
                          const double                 E,
                          const std::valarray<double> &z,
                          const std::valarray<double> &V,
                          const std::valarray<double> &m0,
                          const std::valarray<double> &alpha,
                          const bool                   np_flag);

/**
 * Handler for command-line options
 */
class EFShootOptions : public Options
{
    public:
        EFShootOptions(int argc, char* argv[])
        {
            try
            {
                program_specific_options->add_options()
                    ("nonparabolic,a", po::bool_switch()->default_value(false),
                     "Include nonparabolicity effects.  If selected, the nonparabolicity parameter at each point is read "
                     "from alpha.r.")

                    ("particle,p", po::value<char>()->default_value('e'),
                     "Particle to be used: 'e', 'h' or 'l'")

                    ("states,s", po::value<size_t>()->default_value(1),
                     "Number of states to find")

                    ("dE,d", po::value<double>()->default_value(1e-3),
                     "Minimum separation (in energy) between states [meV]")
                    ;

                std::string doc("Find the eigenstates of an arbitrary 1D potential using a "
                                "shooting method.  The potential profile is read from "
                                "v.r, and the band-edge effective mass (at each point) from "
                                "m.r. "
                                "The energies are written to the file \"E*.r\", and the "
                                "wavefunctions are written to \"wf_*i.r\" where the '*' "
                                "is replaced by the particle ID in each case and the "
                                "'i' is replaced by the number of the state");

                add_prog_specific_options_and_parse(argc, argv, doc);	
            }
            catch(std::exception &e)
            {
                std::cerr << e.what() << std::endl;
                exit(EXIT_FAILURE);
            }
        }

        /// \returns the particle ID
        char get_particle() const {return vm["particle"].as<char>();}

        /// \returns the number of states to find
        size_t get_n_states() const {return vm["states"].as<size_t>();}

        /// \returns the minimum energy spacing between states [J]
        double get_dE() const {return vm["dE"].as<double>()*1e-3*e;}

        /// \returns true if nonparabolicity effects are to be included
        bool nonparabolic() const {return vm["nonparabolic"].as<bool>();}
};

int main(int argc,char *argv[])
{
    const EFShootOptions opt(argc, argv);

    const char   p       = opt.get_particle();
    const bool   np_flag = opt.nonparabolic();
    const double delta_E = opt.get_dE();
    const size_t nst     = opt.get_n_states();

    std::valarray<double> z;
    std::valarray<double> V;
    read_table_xy("v.r", z, V);

    std::valarray<double> z_tmp;
    std::valarray<double> m;
    read_table_xy("m.r", z, m);

    const size_t nz = V.size();

    // Read nonparabolicity data if needed
    std::valarray<double> alpha(nz);

    if(np_flag)
        read_table_xy("alpha.r", z, alpha);

    double Elo=V.min() + delta_E;    // first energy estimate

    shoot_params params = {z, V, m, alpha, np_flag};
    gsl_function f;
    f.function = &psi_at_inf;
    f.params   = &params;
    gsl_root_fsolver *solver = gsl_root_fsolver_alloc(gsl_root_fsolver_brent);

    std::vector<State> solutions;

    for(unsigned int ist=0; ist < nst; ++ist)  
    {
        // Shift the lower estimate up past the last state we found
        if(ist > 0)
        {
            const double E_last = solutions[ist-1].get_E()*e/1000;
            Elo = E_last + delta_E;
        }

        // Value for y=f(x) at bottom of search range
        const double y1 = GSL_FN_EVAL(&f,Elo);

        // Find the range in which the solution lies by incrementing the
        // upper limit of the search range until the function changes sign.
        // Since this coarse search can require many iterations, we keep the
        // lower limit fixed to minimise the amount of computation at each step.
        // The Brent algorithm is extremely fast, so it really doesn't matter that
        // the range we find here is large.
        //
        // Note the end stop to prevent infinite loop in absence of solution
        //
        // TODO: Make the cut-off configurable
        double y2 = y1;
        double Ehi = Elo;
        do
        {
            Ehi += delta_E;
            y2=GSL_FN_EVAL(&f, Ehi);
        }while(y1*y2>0);

        double E = (Elo + Ehi)/2;
        gsl_root_fsolver_set(solver, &f, Elo, Ehi);
        int status = 0;

        // Improve the estimate of the solution using the Brent algorithm
        // until we hit a desired level of precision
        do
        {
            status = gsl_root_fsolver_iterate(solver);
            E   = gsl_root_fsolver_root(solver);
            Elo = gsl_root_fsolver_x_lower(solver);
            Ehi = gsl_root_fsolver_x_upper(solver);
            status = gsl_root_test_interval(Elo, Ehi, 1e-12*e, 0);
        }while(status == GSL_CONTINUE);

        std::valarray<double> psi(z.size());
        const double psi_inf = shoot_wavefunction(psi, E, z, V, m, alpha, np_flag);

        // Check that wavefunction is tightly bound
        // TODO: Implement a better check
        if(gsl_fcmp(fabs(psi_inf), 0, 1) == 1)
            std::cerr << "Warning: Wavefunction is not tightly bound" << std::endl;

        E *= 1e3/e; // Rescale to meV
        solutions.push_back(State(E,psi));
    }

    // Output energies to file
    char energy_filename[9];
    sprintf(energy_filename,"E%c.r",p);

    char wf_prefix[9];
    sprintf(wf_prefix,"wf_%c",p);
    State::write_to_file(energy_filename,
                         wf_prefix,
                         ".r",
                         solutions,
                         z,
                         true);

    return EXIT_SUCCESS;
}

/**
 * \brief Find the wavefunction just beyond the right-hand side of the system
 *
 * \details This function returns the value of the wavefunction (psi)
 *          at +infinity for a given value of the energy.  The solution
 *          to the energy occurs for psi(+infinity)=0.
 *
 * \param[in] E      Energy [J]
 * \param[in] params System parameters (type ShootParams)
 *
 * \returns The wavefunction amplitude immediately to the right of the structure
 */
double psi_at_inf(double  E,
                  void   *params)
{
    const shoot_params *p = reinterpret_cast<shoot_params *>(params);
    std::valarray<double> psi(p->z.size());

    const double psi_inf = shoot_wavefunction(psi, E, p->z, p->V, p->m0, p->alpha, p->np_flag);
    return psi_inf;
}

/**
 * \brief Computes wavefunction iteratively from left to right of structure
 *
 * \details The value of the wavefunction is taken to be zero at the point
 *          immediately to the left of the potential profile. Subsequent
 *          values are computed using QWWAD3, Eq. 3.53.
 *
 * \param[out] wf      Array to which wavefunction will be written [m^{-1/2}]
 * \param[in]  E       Energy at which to compute wavefunction
 * \param[in]  z       Spatial positions [m]
 * \param[in]  V       Potential profile [J]
 * \param[in]  m0      Band-edge effective mass [kg]
 * \param[in]  alpha   Nonparabolicity parameter [1/J]
 * \param[in]  np_flag True if nonparabolicity effects are to be considered
 *
 * \returns The wavefunction amplitude at the point immediately to the right of the structure
 */
double shoot_wavefunction(std::valarray<double>       &wf,
                          const double                 E,
                          const std::valarray<double> &z,
                          const std::valarray<double> &V,
                          const std::valarray<double> &m0,
                          const std::valarray<double> &alpha,
                          const bool                   np_flag)
{
    const size_t nz = z.size();
    wf.resize(nz);
    const double dz = z[1] - z[0];

    // Effective mass (including nonparabolicity)
    std::valarray<double> m(m0);

    // Recalculate effective mass if non-parabolicity is specified
    if(np_flag)
        m = m0*(1.0+alpha*(E-V));

    // boundary conditions (psi[-1] = psi[n] = 0)
    wf[0]   = 1.0;
    double wf_next = 1.0;

    for(unsigned int i=0; i < nz; i++) // last potential not used
    {
        double wf_prev = 0;

        // Compute m(z + dz/2)
        double m_prev = 0.0;
        double m_next = 0.0;

        if(i != 0)
        {
            wf_prev = wf[i-1];
            m_prev = (m[i] + m[i-1])/2.0;
        }
        else
        {
            m_prev = m[i];
        }

        if(i != nz - 1)
            m_next = (m[i] + m[i+1])/2.0;
        else
            m_next = m[i];

        wf_next = (2*m_next*dz*dz/hBar/hBar*(V[i]-E)+
                1.0 + m_next/m_prev)*wf[i]
                - wf_prev * m_next/m_prev;
        wf_prev += 0;

        if(i != nz-1) wf[i+1] = wf_next;
    } 

    // Normalise the wavefunction
    const std::valarray<double> pd = wf*wf;
    const double norm = trapz(pd,dz);
    wf/= sqrt(norm);

    return wf_next/sqrt(norm);
}
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :
