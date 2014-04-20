/**
 * \file   tdb.cpp
 * \brief  Calculate transmission coefficient for double barrier structure
 * \author Paul Harrison  <p.harrison@shu.ac.uk>
 * \author Alex Valavanis <a.valavanis@leeds.ac.uk>
 *
 * \details calculates the transmission coefficient T(E)
 *          for a double barrier structure as a function of the 
 *          incident energy E of a particle.  The barriers can be of 
 *          unequal width and the particle can have a different mass in
 *          the barrier material to the `well', the barrier heights are
 *          however equal.
 */

#include <cstdlib>
#include <cmath>
#include <valarray>
#include "double-barrier.h"
#include "qclsim-constants.h"
#include "qclsim-fileio.h"
#include "qwwad-options.h"

using namespace Leeds;
using namespace constants;

/**
 * Handler for command-line options
 */
class TDBOptions : public Options
{
    public:
        TDBOptions(int argc, char* argv[])
        {
            try
            {
                program_specific_options->add_options()
                    ("left-barrier-width,a", po::value<double>()->default_value(100),
                     "Width of left barrier [angstrom].")

                    ("well-width,b", po::value<double>()->default_value(100),
                     "Width of well [angstrom].")

                    ("right-barrier-width", po::value<double>()->default_value(100),
                     "Width of right barrier [angstrom].")

                    ("well-mass,m", po::value<double>()->default_value(0.067),
                     "Effective mass in well (relative to that of a free electron).")

                    ("barrier-mass,n", po::value<double>()->default_value(0.067),
                     "Effective mass in both barriers (relative to that of a free electron).")

                    ("potential", po::value<double>()->default_value(100),
                     "Barrier potential [meV]")

                    ("energy-step,d", po::value<double>()->default_value(0.01),
                     "Energy step [meV]")
                    ;

                std::string doc("Find the transmission coefficient through a double tunnelling barrier. "
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

        /// \returns the width of the left barrier [m]
        double get_left_barrier_width() const {return vm["left-barrier-width"].as<double>()*1e-10;}

        /// \returns the width of the well [m]
        double get_well_width() const {return vm["well-width"].as<double>()*1e-10;}

        /// \returns the width of the right barrier [m]
        double get_right_barrier_width() const {return vm["right-barrier-width"].as<double>()*1e-10;}

        /// \returns the effective mass in the well [kg]
        double get_well_mass() const {return vm["well-mass"].as<double>()*me;}

        /// \returns the effective mass in the barriers [kg]
        double get_barrier_mass() const {return vm["barrier-mass"].as<double>()*me;}

        /// \returns the effective mass in the quantum well [J]
        double get_potential() const {return vm["potential"].as<double>()*1e-3*e;}
};

int main(int argc,char *argv[])
{
    const TDBOptions opt(argc, argv);

    const double dE  = opt.get_energy_step();         // [J]
    const double L1  = opt.get_left_barrier_width();  // [m]
    const double L2  = opt.get_well_width();          // [m]
    const double L3  = opt.get_right_barrier_width(); // [m]
    const double m_w = opt.get_well_mass();           // [kg]
    const double m_b = opt.get_barrier_mass();        // [kg]
    const double V   = opt.get_potential();           // [J]

    const size_t nE = floor(V/dE); // Number of points in plot

    std::valarray<double> E(nE); // Array of energies
    std::valarray<double> T(nE); // Array of transmission coefficients

    // Loop over energy
    for(unsigned int iE = 0; iE < nE; ++iE)
    {
        E[iE] = iE*dE;
        T[iE] = get_transmission_coefficient(E[iE], m_w, m_b, V, L1, L2, L3);
    }

    // Rescale to meV for output
    E/=(1e-3*e);
    write_table_xy("T.r", E, T);

    return EXIT_SUCCESS;
}
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :
