/**
 * \file  qwwad_ef_parabolic_well.cpp
 * \brief Generate parabolic alloy profile
 *
 * \author Paul Harrison  <p.harrison@shu.ac.uk>
 * \author Alex Valavanis <a.valavanis@leeds.ac.uk>
 *
 * \details This program produces the function of alloy concentration
 *          x with displacement z for a parabolic quantum well,
 *          the output is in a format suitable for conversion into
 *          electron and hole potentials using efxv.      
 *
 *  ----------b----------+ a +----------b---------- x_max
 *                       |   |
 *                       \___/  x_min
 */

#include <cstdio>
#include <cstdlib>
#include <gsl/gsl_math.h>
#include "qwwad/file-io.h"
#include "qwwad/options.h"

using namespace QWWAD;

/**
 * Configure command-line options for the program
 */
auto configure_options(int argc, char** argv) -> Options
{
    Options opt;

    std::string doc("Generate a parabolic alloy profile surrounded by thick barriers.");

    constexpr double W_DEFAULT = 100.0;
    constexpr size_t N_DEFAULT = 301;
    constexpr double X_DEFAULT = 0;
    constexpr double Y_DEFAULT = 0.1;

    opt.add_option<double>("wellwidth,a",    W_DEFAULT, "Width at top of quantum well [angstrom].");
    opt.add_option<double>("barrierwidth,b", W_DEFAULT, "Width of barriers [angstrom].");
    opt.add_option<size_t>("nz,N",           N_DEFAULT, "Number of spatial points for output file.");
    opt.add_option<double>("xmin,x",         X_DEFAULT,   "Minimum alloy fraction.");
    opt.add_option<double>("xmax,y",         Y_DEFAULT, "Maximum alloy fraction.");

    opt.add_prog_specific_options_and_parse(argc, argv, doc);

    return opt;
};

auto main(int argc,char *argv[]) -> int
{
    const auto opt = configure_options(argc, argv);

    const auto a     = opt.get_option<double>("wellwidth") * 1e-10;    // [m]
    const auto b     = opt.get_option<double>("barrierwidth") * 1e-10; // [m]
    const auto nz    = opt.get_option<size_t>("nz");      // Number of points for output file
    const auto x_min = opt.get_option<double>("xmin"); // Minimum alloy fraction
    const auto x_max = opt.get_option<double>("xmax"); // Maximum alloy fraction

    const auto dz = (a+2*b)/(nz-1); // Find width of each spatial interval [m]

    std::vector<double> z(nz); // array of spatial points [m]
    std::vector<double> x(nz); // alloy concentration at each point

    // Loop through spatial points and compute alloy fractions
    for(unsigned int iz = 0; iz < nz; ++iz)
    {
        z[iz] = iz*dz;

        if(gsl_fcmp(z[iz], b, dz/10) == -1 || gsl_fcmp(z[iz], b+a, dz/10) == 1) { // Barriers
            x[iz] = x_max;
	} else {
            x[iz] = x_min+gsl_pow_2(z[iz]-(b+a/2))*(x_max-x_min)/gsl_pow_2(a/2);
	}
    }

    write_table("x.r", z, x);

    return EXIT_SUCCESS;
}
