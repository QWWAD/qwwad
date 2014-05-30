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
#include <gsl/gsl_math.h>
#include "qclsim-constants.h"
#include "qclsim-linalg.h"
#include "qwwad-options.h"
#include "qwwad-schroedinger-shooting.h"

using namespace Leeds;
using namespace constants;

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

                    ("try-energy", po::value<double>(),
                     "Calculate a trial wavefunction at a given energy [meV] and write to "
                     "file.")
                    ;

                std::string doc("Find the eigenstates of an arbitrary 1D potential using a "
                                "shooting method.  The potential profile is read from "
                                "v.r, and the band-edge effective mass (at each point) from "
                                "m.r. "
                                "The energies are written to the file \"E*.r\", and the "
                                "wavefunctions are written to \"wf_*i.r\" where the '*' "
                                "is replaced by the particle ID in each case and the "
                                "'i' is replaced by the number of the state.\n"
                                "\n"
                                "Alternatively, the wavefunction can be calculated at a specified "
                                "trial energy and dumpted to file \"wf_*E.r\".");

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

        /// \returns true if we want to just calculate a trial wavefunction at a given energy
        bool try_energy() const {return (vm.count("try-energy") == 1);}

        /// \returns the energy at which to calculate a trial wavefunction [J]
        double get_trial_energy() const {return vm["try-energy"].as<double>() * 1e-3*e;}
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

    SchroedingerSolverShooting se(m, alpha, V, z, delta_E, nst);

    if (opt.try_energy())
    {
        const double E_trial = opt.get_trial_energy();
        const std::valarray<double> psi = se.trial_wavefunction(E_trial);

        char wf_filename[9];
        sprintf(wf_filename, "wf_%cE.r",p);
        write_table_xy(wf_filename, z, psi);
    }
    else
    {
        std::vector<State> solutions = se.get_solutions(true);

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
    }

    return EXIT_SUCCESS;
}
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :
