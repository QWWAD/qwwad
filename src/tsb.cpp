/**
 * \file   tsb.cpp
 * \brief  Calculate transmission coefficient for single barrier
 * \author Paul Harrison  <p.harrison@shu.ac.uk>
 * \author Alex Valavanis <a.valavanis@leeds.ac.uk>
 */

#include <iostream>
#include <cstdlib>
#include <valarray>
#include <gsl/gsl_math.h>
#include "qclsim-constants.h"
#include "qclsim-fileio.h"
#include "qwwad-options.h"

using namespace Leeds;
using namespace constants;

/**
 * Handler for command-line options
 */
class TSBOptions : public Options
{
    public:
        TSBOptions(int argc, char* argv[])
        {
            try
            {
                program_specific_options->add_options()
                    ("barrier-width,L", po::value<double>()->default_value(100),
                     "Width of barrier [angstrom].")

                    ("mass,m", po::value<double>()->default_value(0.067),
                     "Effective mass (relative to that of a free electron). This is "
                     "assumed to be constant throughout the system")

                    ("potential", po::value<double>()->default_value(100),
                     "Barrier potential [meV]")

                    ("energy-step,d", po::value<double>()->default_value(0.1),
                     "Energy step [meV]")
                    ;

                std::string doc("Find the transmission coefficient through a single tunnelling barrier. "
                                "The values are written as a function of energy to the file T.r.");

                add_prog_specific_options_and_parse(argc, argv, doc);	
            }
            catch(std::exception &e)
            {
                std::cerr << e.what() << std::endl;
                exit(EXIT_FAILURE);
            }
        }

        /// \returns the step in energy required for the output file [J]
        double get_energy_step() const {return vm["energy-step"].as<double>()*1e-3*e;}

        /// \returns the width of the barriers [m]
        double get_barrier_width() const {return vm["barrier-width"].as<double>()*1e-10;}

        /// \returns the effective mass [kg]
        double get_mass() const {return vm["mass"].as<double>()*me;}

        /// \returns the effective mass in the quantum well [J]
        double get_potential() const {return vm["potential"].as<double>()*1e-3*e;}
};

int main(int argc,char *argv[])
{
    const TSBOptions opt(argc, argv);

    const double dE = opt.get_energy_step();   // [J]
    const double L  = opt.get_barrier_width(); // [m]
    const double m  = opt.get_mass();          // [kg]
    const double V  = opt.get_potential();     // [J]

    const double E_cutoff = V * 10; // Cut-off energy for plot
    const size_t nE = floor(E_cutoff/dE); // Number of points in plot

    std::valarray<double> E(nE); // Array of energies
    std::valarray<double> T(nE); // Array of transmission coefficients

    for(unsigned int iE = 0; iE < nE; ++iE)
    {
        E[iE] = iE*dE; // Find energy
        const double k=sqrt(2*m*E[iE])/hBar; // Wave-vector in `well' material

        if(gsl_fcmp(E[iE], V, dE/1000) == -1) // If E < V
        {
            const double K = sqrt(2*m*(V-E[iE]))/hBar; // Decay constant in barrier
            T[iE] = 1/(1+gsl_pow_2((k*k+K*K)/(2*k*K) * sinh(K*L))); // [QWWAD4, 2.199]
        }
        else // if E > V
        {
            const double kdash = sqrt(2*m*(E[iE]-V))/hBar; // Wave-vector above barrier
            T[iE] = 1/(1+gsl_pow_2((k*k-kdash*kdash)/(2*k*kdash) * sin(kdash*L))); // [QWWAD4, 2.201]
        }
    }

    // Rescale to meV for output
    E/=(1e-3*e);
    write_table_xy("T.r", E, T);

    return EXIT_SUCCESS;
}
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :
