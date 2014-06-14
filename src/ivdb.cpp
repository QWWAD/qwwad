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
#include <iostream>
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
 * Configure command-line options for the program
 */
Options configure_options(int argc, char* argv[])
{
    Options opt;
    opt.add_numeric_option("left-barrier-width,l",    100, "Width of left barrier [angstrom].");
    opt.add_numeric_option("well-width,b",            100, "Width of well [angstrom].");
    opt.add_numeric_option("right-barrier-width,r",   100, "Width of right barrier [angstrom].");
    opt.add_numeric_option("well-mass,m",           0.067, "Effective mass in well (relative to that of a free electron). This is used to compute the distribution of carriers at the input.");
    opt.add_numeric_option("barrier-mass,n",        0.067, "Effective mass in barrier (relative to that of a free electron).");
    opt.add_numeric_option("potential",               100, "Barrier potential [meV]");
    opt.add_numeric_option("temperature",             300, "Temperature of carrier distribution [K]");
    opt.add_numeric_option("energy-step,d",          0.01, "Energy step [meV]");

    std::string summary("Find the current through a double barrier structure as a function of bias voltage.");
    std::string details("The following output text file is created:\n"
                        "  'IV.r'   \tCurrent--voltage characteristics:\n"
                        "           \tCOLUMN 1: voltage [V].\n"
                        "           \tCOLUMN 2: current [a.u.].\n"
                        "  'T.r'    \tTransmission coefficient as a function of energy:\n"
                        "           \tCOLUMN 1: energy of incident carrier [meV].\n"
                        "           \tCOLUMN 2: transmission coefficient.\n"
                        "\n"
                        "Examples:\n"
                        "   Compute current--voltage characteristics through a pair of 100-meV, 200-angstrom barriers separated by a 50-angstrom well at 50 K temperature:\n\n"
                        "   ivdb --left-barrier-width 200 --right-barrier-width 200 --potential 100 --well-width 50 --temperature 50\n"
                        "\n"
                        "   Compute current--voltage characteristics through a pair of 1-eV, 100-angstrom barrier, with effective mass = 0.1 m0, using resolution of 1 meV:\n\n"
                        "   ivdb --potential 1000 --left-barrier-width 100 --right-barrier-width 100 --barrier-mass 0.1 --energy-step 1\n");

    opt.add_prog_specific_options_and_parse(argc, argv, summary, details);

    return opt;
}

int main(int argc,char *argv[])
{
    const Options opt = configure_options(argc, argv);

    const double dE  = opt.get_numeric_option("energy-step") * 1e-3 * e;      // [J]
    const double L1  = opt.get_numeric_option("left-barrier-width") * 1e-10;  // [m]
    const double L2  = opt.get_numeric_option("well-width") * 1e-10;          // [m]
    const double L3  = opt.get_numeric_option("right-barrier-width") * 1e-10; // [m]
    const double m_w = opt.get_numeric_option("well-mass") * me;              // [kg]
    const double m_b = opt.get_numeric_option("barrier-mass") * me;           // [kg]
    const double Vb  = opt.get_numeric_option("potential") * e / 1000;        // [J]
    const double T   = opt.get_numeric_option("temperature");                 // [K]

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
