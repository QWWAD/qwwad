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
#include "qclsim-constants.h"
#include "qwwad-options.h"
#include "qclsim-fileio.h"
#include "qwwad-debye.h"

using namespace Leeds;
using namespace constants;

/**
 * Configure command-line options for the program
 */
Options configure_options(int argc, char* argv[])
{
    Options opt;

    std::string summary("Tabulate the specific heat capacity of a material at a range of temperatures.");

    opt.add_numeric_option("debye",            360, "Debye temperature [K]");
    opt.add_numeric_option("molarmass,m", 0.144645, "Molar mass [kg/mol]");
    opt.add_size_option   ("natoms",             2, "Number of atoms in molecular unit");
    opt.add_numeric_option("Tmin",               1, "Minimum temperature to compute [K]");
    opt.add_numeric_option("Tmax",             300, "Maximum temperature to compute [K]");
    opt.add_numeric_option("Tstep",              1, "Step in temperature [K]");
    opt.add_string_option ("filename",       "c.r", "File in which to save specific heat capacity data");

    opt.add_prog_specific_options_and_parse(argc, argv, summary);

    return opt;
}

int main(int argc,char *argv[])
{
    const Options opt = configure_options(argc, argv);

    const double T_D    = opt.get_numeric_option("debye");     // Debye temperature for material [K]
    const double Tmin   = opt.get_numeric_option("Tmin");      // Minimum temperature for loop [K]
    const double Tmax   = opt.get_numeric_option("Tmax");      // Maximum temperature for loop [K]
    const double dT     = opt.get_numeric_option("Tstep");     // Temperature step for loop [K]
    const double M      = opt.get_numeric_option("molarmass"); // Molar mass [kg/mol]
    const size_t natoms = opt.get_size_option("natoms");       // Number of atoms in molecular unit

    const size_t nT = 1 + (Tmax-Tmin)/dT;

    std::valarray<double> T(nT);  // Array of temperatures [K]
    std::valarray<double> cp(nT); // Array of spec. heat. capacity [J/(kg K)]

    DebyeModel dm(T_D, M, natoms);

    // Loop over temperature
    for(unsigned int iT = 0; iT < nT; ++iT)
    {
        T[iT]  = Tmin + iT*dT;
        cp[iT] = dm.get_cp(T[iT]);
//        cp[iT] = dm.get_cp_low_T(T[iT]);
    }

    write_table(opt.get_string_option("filename").c_str(), T, cp);

    return EXIT_SUCCESS;
}

// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :
