/**
 * \file   efiw.cpp Calculates eigenstates of an infinite square well
 *
 * \author Paul Harrison  <p.harrison@shu.ac.uk>
 * \author Alex Valavanis <a.valavanis@leeds.ac.uk>
 */

#include <iostream>
#include <cstdlib>
#include <cmath>
#include <valarray>
#include <gsl/gsl_math.h>
#include "qclsim-constants.h"
#include "qclsim-linalg.h"
#include "qclsim-schroedinger.h"
#include "qwwad-options.h"

using namespace Leeds;
using namespace constants;

/**
 * Handler for command-line options
 */
class EFIWOptions : public Options
{
    public:
        EFIWOptions(int argc, char* argv[])
        {
            try
            {
                program_specific_options->add_options()
                    ("width,L", po::value<double>()->default_value(100), 
                     "Width of quantum well [angstrom].")

                    ("mass,m", po::value<double>()->default_value(0.067),
                     "Effective mass (relative to that of a free electron)")

                    ("npoint,N", po::value<size_t>()->default_value(100),
                     "Number of spatial points for output file")

                    ("particle,p", po::value<char>()->default_value('e'),
                     "Particle to be used: 'e', 'h' or 'l'")

                    ("states,s", po::value<size_t>()->default_value(1),
                     "Number of states to find")
                    
                    ("band-edge", po::value<double>()->default_value(0),
                     "Band edge potential [eV]")

                    ("alpha", po::value<double>()->default_value(0),
                     "Nonparabolicity parameter [eV^{-1}]")
                    ;

                std::string doc("Find the eigenstates of an infinite quantum well. "
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

        /**
         * \returns the width of the quantum well [m]
         */
        double get_width() const {return vm["width"].as<double>()*1e-10;}

        /**
         * \returns the effective mass [kg]
         */
        double get_mass() const {return vm["mass"].as<double>()*me;}

        /**
         * \returns the particle ID
         */
        char get_particle() const {return vm["particle"].as<char>();}

        /**
         * \returns the number of spatial points
         */
        size_t get_n_points() const {return vm["npoint"].as<size_t>();}

        /**
         * \returns the number of spatial points
         */
        size_t get_n_states() const {return vm["states"].as<size_t>();}

        /// \returns the band edge [J]
        double get_band_edge() const {return vm["band-edge"].as<double>()*e;}

        /// \returns the nonparabolicity parameter [J^{-1}]
        double get_alpha() const {return vm["alpha"].as<double>()/e;}
};

int main(int argc, char *argv[])
{
    EFIWOptions opt(argc, argv);

    const double L     = opt.get_width();    // well width [m]
    const char   p     = opt.get_particle(); // particle (e, h or l)
    const double m     = opt.get_mass();     // effective mass [kg]
    const size_t N     = opt.get_n_points(); // number of spatial steps
    const size_t s     = opt.get_n_states(); // number of states
    const double alpha = opt.get_alpha();
    const double V     = opt.get_band_edge();

    SchroedingerSolverInfWell se(m, L, N, alpha, V, s);

    std::vector<State> solutions = se.get_solutions(true);

    // Dump to file
    char energy_filename[9];
    sprintf(energy_filename,"E%c.r",p);

    char wf_prefix[9];
    sprintf(wf_prefix,"wf_%c",p);
    State::write_to_file(energy_filename,
                         wf_prefix,
                         ".r",
                         solutions,
                         se.get_z(),
                         true);

    return EXIT_SUCCESS;
}
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :
