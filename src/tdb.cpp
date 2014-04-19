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
#include <armadillo>
#include "qclsim-constants.h"
#include "qclsim-fileio.h"
#include "qwwad-options.h"

using namespace arma;
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

                    ("right-barrier-width,c", po::value<double>()->default_value(100),
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

    // Calculate interfaces
    const double I2=L1;
    const double I3=L1+L2;
    const double I4=L1+L2+L3;

    // Loop over energy
    for(unsigned int iE = 0; iE < nE; ++iE)
    {
        E[iE] = iE*dE; // Find energy
        const double k=sqrt(2*m_w*E[iE])/hBar;
        const double K=sqrt(2*m_b*(V-E[iE]))/hBar;

        // Define transfer matrices
        cx_mat M1(2,2);
        M1(0,0) = 1;
        M1(0,1) = 1;
        M1(1,0) = cx_double(0.0,  k/m_w);
        M1(1,1) = cx_double(0.0, -k/m_w);

        cx_mat M2(2,2);
        M2(0,0) = 1;
        M2(0,1) = 1;
        M2(1,0) = +K/m_b;
        M2(1,1) = -K/m_b;

        cx_mat M3(2,2);
        M3(0,0) = exp(+K*I2);
        M3(0,1) = exp(-K*I2);
        M3(1,0) =  K*exp(+K*I2)/m_b;
        M3(1,1) = -K*exp(-K*I2)/m_b;

        cx_mat M4(2,2);
        M4(0,0) = cx_double(cos(k*I2),         +sin(k*I2));
        M4(0,1) = cx_double(cos(k*I2),         -sin(k*I2));
        M4(1,0) = cx_double(-k*sin(+k*I2)/m_w, +k*cos(k*I2)/m_w);
        M4(1,1) = cx_double(+k*sin(-k*I2)/m_w, -k*cos(k*I2)/m_w);

        cx_mat M5(2,2);
        M5(0,0) = cx_double(cos(k*I3),         +sin(k*I3));
        M5(0,1) = cx_double(cos(k*I3),         -sin(k*I3));
        M5(1,0) = cx_double(-k*sin(+k*I3)/m_w, +k*cos(k*I3)/m_w);
        M5(1,1) = cx_double(+k*sin(-k*I3)/m_w, -k*cos(k*I3)/m_w);

        cx_mat M6(2,2);
        M6(0,0) = exp(+K*I3);
        M6(0,1) = exp(-K*I3);
        M6(1,0) =  K*exp(+K*I3)/m_b;
        M6(1,1) = -K*exp(-K*I3)/m_b;

        cx_mat M7(2,2);
        M7(0,0) = exp(+K*I4);
        M7(0,1) = exp(-K*I4);
        M7(1,0) =  K*exp(+K*I4)/m_b;
        M7(1,1) = -K*exp(-K*I4)/m_b;

        cx_mat M8(2,2);
        M8(0,0) = cx_double(cos(k*I4),         +sin(k*I4));
        M8(0,1) = cx_double(cos(k*I4),         -sin(k*I4));
        M8(1,0) = cx_double(-k*sin(+k*I4)/m_w, +k*cos(k*I4)/m_w);
        M8(1,1) = cx_double(+k*sin(-k*I4)/m_w, -k*cos(k*I4)/m_w);

        // Little hack to stop nonsense output when E = 0
        if (iE > 0)
        {
            cx_mat M = inv(M1) * M2 * inv(M3) * M4 * inv(M5) * M6 * inv(M7) * M8;
            T[iE] = 1/(norm(M(1,1))); // Transission coeff
        }
    }

    // Rescale to meV for output
    E/=(1e-3*e);
    write_table_xy("T.r", E, T);

    return EXIT_SUCCESS;
}
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :
