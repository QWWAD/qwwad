/**
 * \file   ivdb.cpp
 * \brief  Calculate I-V for double barrier structure
 * \author Paul Harrison  <p.harrison@shu.ac.uk>
 * \author Alex Valavanis <a.valavanis@leeds.ac.uk>
 *
 * \details This code reads in the transmission coefficient T(E) for a 
 *          double barrier and uses a very simple model to calculate 
 *          an I-V curve.
 */

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <valarray>
#include <gsl/gsl_math.h>
#include "dos-functions.h"
#include "double-barrier.h"
#include "qclsim-constants.h"
#include "qclsim-fermi.h"
#include "qclsim-fileio.h"
#include "qwwad-options.h"

using namespace Leeds;
using namespace constants;

/**
 * Handler for command-line options
 */
class IVDBOptions : public Options
{
    public:
        IVDBOptions(int argc, char* argv[])
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
                     "Effective mass in well (relative to that of a free electron). This is used "
                     "to compute the distribution of carriers at the input.")

                    ("barrier-mass,n", po::value<double>()->default_value(0.067),
                     "Effective mass in both barriers (relative to that of a free electron).")

                    ("temperature,T", po::value<double>()->default_value(300),
                     "Temperature of carrier distribution [K]")

                    ("potential", po::value<double>()->default_value(100),
                     "Barrier potential [meV]")

                    ("energy-step,d", po::value<double>()->default_value(0.01),
                     "Energy step [meV]")
                    ;

                std::string doc("Find the current through a double barrier structure as a "
                                "function of voltage.  The values are written in (V, I) format "
                                "to the file IV.r. The transmission coefficient is also written "
                                "to the file T.r as a function of energy.");

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

        /// \returns the temperature of the carrier distribution [K]
        double get_temperature() const {return vm["temperature"].as<double>();}
};

int main(int argc,char *argv[])
{
    IVDBOptions opt(argc, argv);
    const double dE  = opt.get_energy_step();         // [J]
    const double L1  = opt.get_left_barrier_width();  // [m]
    const double L2  = opt.get_well_width();          // [m]
    const double L3  = opt.get_right_barrier_width(); // [m]
    const double m_w = opt.get_well_mass();           // [kg]
    const double m_b = opt.get_barrier_mass();        // [kg]
    const double T   = opt.get_temperature();         // [K]
    const double Vb  = opt.get_potential();           // [J]

    const size_t nE = floor(Vb/dE); // Number of points in table of energies
    const double Ef=2*1e-3*e;      // Just set fixed Fermi energy to represent some fixed density

    std::valarray<double> Tx(nE); // Transmission coefficient
    std::valarray<double> E(nE);  // Energy [J]

    // Compute transmission coefficient at each energy
    for(unsigned int iE = 0; iE < nE; ++iE)
    {
        E[iE] = iE*dE;
        Tx[iE] = get_transmission_coefficient(E[iE], m_w, m_b, Vb, L1, L2, L3);
    }

    write_table_xy("T.r", E, Tx);

    const size_t nF = 100; // Number of field points
    std::valarray<double> V(nF); // Voltage drop across structure [V]
    std::valarray<double> current(nF);

    // Loop over field
    for(unsigned int iF=0; iF<nF; ++iF)
    {
        const double F = iF*1e+5;            // Electric field [Vm^{-1}]
        const double DeltaE=e*F*(L1+0.5*L2); // Field induced shift in band energy
        V[iF] = F*(L1+L2+L3);

        current[iF] = 0; // Initialise current for integration
        for(unsigned int iE=0; iE<nE; ++iE)
        {
            if(E[iE] > DeltaE)	/* only add contribution to current if E>DeltaE	*/
            {
                const double rho   = calculate_dos_3D(m_w, E[iE] - DeltaE); // Bulk density of states
                const double _f_FD = f_FD(Ef + DeltaE, E[iE], T); // Fermi function

                current[iF] += Tx[iE]*_f_FD*rho*dE;
            }
        }
    }

    write_table_xy("IV.r", V, current);

    return EXIT_SUCCESS;
}
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :
