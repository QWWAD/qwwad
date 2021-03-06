/**
 * \file   critical-thickness.cpp
 * \brief  Calculate critical thickness for a thin layer in a cubic crystal
 * \author Alex Valavanis <a.valavanis@leeds.ac.uk>
 */

#include <cstdlib>
#include <cmath>
#include <iostream>
#include <gsl/gsl_sf_lambert.h>
#include "qwwad/constants.h"
#include "qwwad/file-io.h"
#include "qwwad/options.h"

using namespace QWWAD;
using namespace constants;

/**
 * Configure command-line options for the program
 */
auto configure_options(int argc, char** argv) -> Options
{
    Options opt;

    std::string summary("Find the critical thickness of a thin film.");

    opt.add_option<double>("substrate", 0.0,     "Substrate alloy fraction [0 to 1]");
    opt.add_option<double>("C110",      165.773, "Elastic constant C11 for alloy=0");
    opt.add_option<double>("C111",      128.528, "Elastic constant C11 for alloy=1");
    opt.add_option<double>("C120",       63.924, "Elastic constant C12 for alloy=0");
    opt.add_option<double>("C121",       48.260, "Elastic constant C12 for alloy=1");
    opt.add_option<double>("a0",          5.431, "Lattice constant for alloy=0");
    opt.add_option<double>("a1",          5.633, "Lattice constant for alloy=1");

    opt.add_prog_specific_options_and_parse(argc, argv, summary);

    return opt;
}

auto main(int argc,char *argv[]) -> int
{
    const auto opt = configure_options(argc, argv);

    const auto xs    = opt.get_option<double>("substrate");
    const auto C11_0 = opt.get_option<double>("C110");
    const auto C11_1 = opt.get_option<double>("C111");
    const auto C12_0 = opt.get_option<double>("C120");
    const auto C12_1 = opt.get_option<double>("C121");
    const auto a_0   = opt.get_option<double>("a0");
    const auto a_1   = opt.get_option<double>("a1");

    const size_t nx = 100;

    std::vector<double> x(nx);  // Array of energies
    std::vector<double> hc(nx); // Array of critical thickness

    const double a_subst = a_0*(1.0-xs) + a_1*xs;

    // Loop over alloy fraction
    for(unsigned int ix = 0; ix < nx; ++ix)
    {
        x[ix] = ix/100.0;
        const double C11 = C11_0*(1-x[ix]) + C11_1*x[ix];
        const double C12 = C12_0*(1-x[ix]) + C12_1*x[ix];
        const double a = a_0*(1-x[ix]) + a_1*x[ix];
        const double b = a/sqrt(2); // Magnitude of Burgers vector
        const double f = fabs(a - a_subst)/a_subst;
        const double nu = 2*C12/C11;

        const double A = b/(2*pi*f) * (1-0.25*nu)/((1+nu)*0.5);
        const double B = exp(1)/b;

        hc[ix] = -A * gsl_sf_lambert_Wm1(-1/(A*B));
    }

    std::string outfile("hc.r");

    try {
        write_table(outfile, x, hc);
    } catch (std::runtime_error &e) {
        std::cerr << "Could not write to " << outfile << std::endl;
        std::cerr << e.what();
    }

    return EXIT_SUCCESS;
}
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :
