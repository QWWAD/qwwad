/**
 * \file   qwwad_tx_double_barrier_iv.cpp
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
#include "qwwad/dos-functions.h"
#include "qwwad/double-barrier.h"
#include "qwwad/constants.h"
#include "qwwad/fermi.h"
#include "qwwad/file-io.h"
#include "qwwad/options.h"

using namespace QWWAD;
using namespace constants;

/**
 * Configure command-line options for the program
 */
Options configure_options(int argc, char** argv)
{
    Options opt;

    std::string summary("Find the current through a double barrier structure as a function of bias voltage.");

    opt.add_option<double>("leftbarrierwidth,l",    100, "Width of left barrier [angstrom].");
    opt.add_option<double>("wellwidth,b",           100, "Width of well [angstrom].");
    opt.add_option<double>("rightbarrierwidth,r",   100, "Width of right barrier [angstrom].");
    opt.add_option<double>("wellmass,m",          0.067, "Effective mass in well (relative to that "
                                                         "of a free electron). This is used to "
                                                         "compute the distribution of carriers at the input.");
    opt.add_option<double>("barriermass,n",       0.067, "Effective mass in barrier (relative to "
                                                         "that of a free electron).");
    opt.add_option<double>("barrierpotential",      100, "Barrier potential [meV]");
    opt.add_option<double>("Te",                    300, "Temperature of carrier distribution [K]");
    opt.add_option<double>("dE,d",                 0.01, "Energy step [meV]");

    opt.add_prog_specific_options_and_parse(argc, argv, summary);

    return opt;
}

int main(int argc,char *argv[])
{
    const auto opt = configure_options(argc, argv);

    const auto dE  = opt.get_option<double>("dE") * 1e-3 * e;               // [J]
    const auto L1  = opt.get_option<double>("leftbarrierwidth") * 1e-10;    // [m]
    const auto L2  = opt.get_option<double>("wellwidth") * 1e-10;           // [m]
    const auto L3  = opt.get_option<double>("rightbarrierwidth") * 1e-10;   // [m]
    const auto m_w = opt.get_option<double>("wellmass") * me;               // [kg]
    const auto m_b = opt.get_option<double>("barriermass") * me;            // [kg]
    const auto Vb  = opt.get_option<double>("barrierpotential") * e / 1000; // [J]
    const auto Te  = opt.get_option<double>("Te");                          // [K]

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

    write_table("T.r", E, Tx);

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
                const double _f_FD = f_FD(Ef + DeltaE, E[iE], Te); // Fermi function

                current[iF] += Tx[iE]*_f_FD*rho*dE;
            }
        }
    }

    write_table("IV.r", V, current);

    return EXIT_SUCCESS;
}
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :
