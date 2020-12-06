#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <valarray>
#include <gsl/gsl_math.h>
#include <gsl/gsl_deriv.h>
#include <gsl/gsl_sf_debye.h>
#include "qwwad/constants.h"
#include "qwwad/debye.h"
#include "qwwad/options.h"
#include "qwwad/file-io.h"

using namespace QWWAD;
using namespace constants;

/**
 * Configure command-line options for the program
 */
Options configure_options(int argc, char** argv)
{
    Options opt;

    std::string summary("Tabulate the specific heat capacity of a material at a range of temperatures.");

    opt.add_option<double>     ("debye",            360, "Debye temperature [K]");
    opt.add_option<double>     ("molarmass,m", 0.144645, "Molar mass [kg/mol]");
    opt.add_option<size_t>     ("natoms",             2, "Number of atoms in molecular unit");
    opt.add_option<double>     ("Tmin",               1, "Minimum temperature to compute [K]");
    opt.add_option<double>     ("Tmax",             300, "Maximum temperature to compute [K]");
    opt.add_option<double>     ("Tstep",              1, "Step in temperature [K]");
    opt.add_option<std::string>("filename",       "c.r", "File in which to save specific heat capacity data");
    opt.add_option<bool>       ("approx"               , "Use quick low/high temperature appriximation");

    opt.add_prog_specific_options_and_parse(argc, argv, summary);

    return opt;
}

int main(int argc,char *argv[])
{
    const auto opt = configure_options(argc, argv);

    const auto T_D    = opt.get_option<double>("debye");     // Debye temperature for material [K]
    const auto Tmin   = opt.get_option<double>("Tmin");      // Minimum temperature for loop [K]
    const auto Tmax   = opt.get_option<double>("Tmax");      // Maximum temperature for loop [K]
    const auto dT     = opt.get_option<double>("Tstep");     // Temperature step for loop [K]
    const auto M      = opt.get_option<double>("molarmass"); // Molar mass [kg/mol]
    const auto natoms = opt.get_option<size_t>("natoms");    // Number of atoms in molecular unit
    const auto approx = opt.get_option<bool>  ("approx");

    const auto nT = 1 + (Tmax-Tmin)/dT;

    std::valarray<double> T(nT);  // Array of temperatures [K]
    std::valarray<double> cp(nT); // Array of spec. heat. capacity [J/(kg K)]

    DebyeModel dm(T_D, M, natoms);

    // Loop over temperature
    for(unsigned int iT = 0; iT < nT; ++iT)
    {
        T[iT]  = Tmin + iT*dT;

        if(approx)
            cp[iT] = dm.get_cp_approx(T[iT]);
        else
            cp[iT] = dm.get_cp(T[iT]);
    }

    write_table(opt.get_option<std::string>("filename").c_str(), T, cp);

    return EXIT_SUCCESS;
}

// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :
