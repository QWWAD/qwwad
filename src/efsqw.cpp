/**
 * \file    efsqw.cpp
 * \brief   Calculate the 1-particle energies and wavefunctions in a single quantum well
 * \author  Paul Harrison  <p.harrison@shu.ac.uk>
 * \author  Alex Valavanis <a.valavanis@leeds.ac.uk>
 *
 * \details One particle energies and wavefunctions in a single
 *          quantum well.  Three basic theoretical approaches are
 *          contained within this sourcecode
 *
   (i) Constant mass

   \f[
   -\frac{\hbar^2}{2m^*} \frac{\mathrm{d}^2}{\mathrm{d}z^2}\psi + V(z) \psi = E \psi,   -\frac{l_w}{2} < z < +\frac{l_w}{2},  V=0
   \f]

   (ii) Different masses in well and barrier

   \f[
   -\frac{\hbar^2}{2m^*} \frac{\mathrm{d}^2}{\mathrm{d}z^2}\psi + V(z) \psi = E \psi,   -l_w < z < +l_w,  V=0
   \f]

   with the additional constraint of the boundary conditions

   \f[
   \psi \mathrm{and} \frac{1}{m}\frac{d\psi}{dz}, \mathrm{continuous}
   \f]
 
   this represents the Hamiltonian \f$P_z\frac{1}{m] P_z + V\f$

   (iii) Different masses in well and barrier

   \f[
   -\frac{\hbar^2}{2m^*} \frac{\mathrm{d}^2}{\mathrm{d}z^2}\psi + V(z) \psi = E \psi,   -l_w < z < +l_w,  V=0
   \f]

   with the additional constraint of the boundary conditions

   \f[
   \psi \mathrm{and} \frac{d\psi}{dz}, \mathrm{continuous}
   \f]

   this represents the Hamiltonian \f$\frac{1}{\sqrt{m}}P_zP_z\frac{1}{\sqrt{m}] + V\f$

   The code is based around point (ii). Point (i) is 
   implemented by the user selecting m_b=m_w.  Point
   (iii) is implemented by allowing different m_b and m_w
   in the evaluation of k and K, but m_b is artificially
   forced to equal m_w for the boundary conditions.

   The system is solved by expressing the
   standard condition as a function f(x)=0 and
   implementing a Newton-Raphson iterative method
   where the independent variable x is the energy.
   The first estimate of the energy is found by 
   following the function along the x axis until it
   changes sign then using a midpoint rule.
*/

#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <valarray>
#include "qclsim-constants.h"
#include "qwwad-schroedinger-finite-well.h"
#include "qwwad-options.h"

using namespace Leeds;
using namespace constants;

#define   N  1001      /* number of wavefunction sampling
                          points, must be odd to provide an
                          even number of strips for numerical 
                          integration.                            */

/**
 * Handler for command-line options
 */
class EFSQWOptions : public Options
{
    public:
        EFSQWOptions(int argc, char* argv[])
        {
            try
            {
                program_specific_options->add_options()
                    ("well-width,a", po::value<double>()->default_value(100),
                     "Width of quantum well [angstrom].")
                    
                    ("barrier-width,b", po::value<double>()->default_value(200),
                     "Width of barrier [angstrom]. Note that this is only used "
                     "for the purposes of outputting the data. The calculation here "
                     "assumes that the barriers are infinitely thick.  As such, the "
                     "wavefunctions do not decay to precisely zero at the boundaries.")

                    ("alt-KE,k", po::bool_switch()->default_value(false),
                     "Use alternative kinetic energy operator (1/m)PP")

                    ("well-mass,m", po::value<double>()->default_value(0.067),
                     "Effective mass in well (relative to that of a free electron)")

                    ("barrier-mass,n", po::value<double>()->default_value(0.067),
                     "Effective mass in barrier (relative to that of a free electron)")

                    ("output-equations", po::bool_switch()->default_value(false),
                     "Output the matching equations for the system. The left-hand "
                     "side of the equation is output to 'lhs.r' and each branch of "
                     "the right-hand side is output to a set of 'rhs_i.r' files, "
                     "where 'i' is the index of the state that lies on that branch. "
                     "An rhs file is output for all the bound states in the system and "
                     "one additional branch (with no real solution)")

                    ("particle,p", po::value<char>()->default_value('e'),
                     "Particle to be used: 'e', 'h' or 'l'")

                    ("states,s", po::value<size_t>()->default_value(1),
                     "Number of states to find")

                    ("potential", po::value<double>()->default_value(100),
                     "Barrier potential [meV]")
                    ;

                std::string doc("Find the eigenstates of a single finite quantum well. "
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

        /// \returns the width of the quantum well [m]
        double get_well_width() const {return vm["well-width"].as<double>()*1e-10;}
        
        /// \returns the width of the barriers [m]
        double get_barrier_width() const {return vm["barrier-width"].as<double>()*1e-10;}

        /// \returns true if normal kinetic energy operator is to be used
        bool get_T_flag() const {return !vm["alt-KE"].as<bool>();}

        /// \returns true if we are required to output the matching equations for the system
        bool output_equations() const {return vm["output-equations"].as<bool>();}

        /// \returns the effective mass in the quantum well [kg]
        double get_well_mass() const {return vm["well-mass"].as<double>()*me;}

        /// \returns the effective mass in the quantum well [kg]
        double get_barrier_mass() const {return vm["barrier-mass"].as<double>()*me;}

        /// \returns the particle ID
        char get_particle() const {return vm["particle"].as<char>();}

        /// \returns the number of spatial points
        size_t get_n_states() const {return vm["states"].as<size_t>();}

        /// \returns the effective mass in the quantum well [J]
        double get_potential() const {return vm["potential"].as<double>()*1e-3*e;}
};

int main(int argc,char *argv[])
{
    EFSQWOptions opt(argc, argv);
    const double a      = opt.get_well_width();
    const double b      = opt.get_barrier_width();
    const double m_w    = opt.get_well_mass();
    const double m_b    = opt.get_barrier_mass();
    const char   p      = opt.get_particle();
    const double V      = opt.get_potential();   
    const size_t state  = opt.get_n_states();
    const bool   T_flag = opt.get_T_flag();

    SchroedingerSolverFiniteWell se(a, b, V, m_w, m_b, N, T_flag,state);

    if(opt.output_equations())
    {
        const size_t nst    = se.get_n_bound();
        const double v_max  = (nst+1)*pi/2;

        const size_t nv = (nst+1)*1000;
        const double dv = v_max/(nv-1);
        std::valarray<double> v(nv);
        std::valarray<double> lhs(nv);

        // Output lhs data in a single contiguous chunk
        for (unsigned int iv = 0; iv < nv; ++iv)
        {
            v[iv]   = iv*dv;
            lhs[iv] = se.get_lhs(v[iv]);
        }

        Leeds::write_table_xy("lhs.r", v, lhs);

        // Output a separate file for each rhs branch
        for (unsigned int ibranch = 0; ibranch < nst+1; ++ibranch)
        {
            const size_t nv_branch = 1000;

            std::valarray<double> v_branch(nv_branch);
            std::valarray<double> rhs(nv_branch);

            // Set the spacing between points so that we don't quite reach the
            // asymptote at the "end" of the branch
            const double dv_branch = (pi/2.0*0.999999)/(nv_branch-1);

            for (unsigned int iv_branch = 0; iv_branch < nv_branch; ++iv_branch)
            {
                v_branch[iv_branch] = pi/2.0 * ibranch  + iv_branch*dv_branch;
                rhs[iv_branch]      = se.get_rhs(v_branch[iv_branch]);
            }

            // Set filename
            std::ostringstream oss;
            oss << "rhs_" << ibranch+1 << ".r";
            Leeds::write_table_xy(oss.str().c_str(), v_branch, rhs);
        }
    }

    // Dump to file
    char energy_filename[9];
    sprintf(energy_filename,"E%c.r",p);

    char wf_prefix[9];
    sprintf(wf_prefix,"wf_%c",p);
    State::write_to_file(energy_filename,
                         wf_prefix,
                         ".r",
                         se.get_solutions(true),
                         se.get_z(),
                         true);

    return EXIT_SUCCESS;
}
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :
