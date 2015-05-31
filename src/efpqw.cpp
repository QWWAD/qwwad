/**
 * \file  efpqw.cpp
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
Options configure_options(int argc, char* argv[])
{
    Options opt;

    std::string doc("Generate a parabolic alloy profile surrounded by thick barriers.");

    opt.add_option<double>("well-width,a",    100, "Width at top of quantum well [angstrom].");
    opt.add_option<double>("barrier-width,b", 100, "Width of barriers [angstrom].");
    opt.add_option<size_t>("nz,N",            301, "Number of spatial points for output file.");
    opt.add_option<double>("xmin,x",          0,   "Minimum alloy fraction.");
    opt.add_option<double>("xmax,y",          0.1, "Maximum alloy fraction.");

    opt.add_prog_specific_options_and_parse(argc, argv, doc);

    return opt;
};

int main(int argc,char *argv[])
{
    const auto opt = configure_options(argc, argv);

    const auto a     = opt.get_option<double>("well-width") * 1e-10;    // [m]
    const auto b     = opt.get_option<double>("barrier-width") * 1e-10; // [m]
    const auto nz    = opt.get_option<size_t>("nz");      // Number of points for output file
    const auto x_min = opt.get_option<double>("xmin"); // Minimum alloy fraction
    const auto x_max = opt.get_option<double>("xmax"); // Maximum alloy fraction

    const auto dz = (a+2*b)/(nz-1); // Find width of each spatial interval [m]

    std::valarray<double> z(nz); // array of spatial points [m]
    std::valarray<double> x(nz); // alloy concentration at each point

    // Loop through spatial points and compute alloy fractions
    for(unsigned int iz = 0; iz < nz; ++iz)
    {
        z[iz] = iz*dz;

        if(gsl_fcmp(z[iz], b, dz/10) == -1 || gsl_fcmp(z[iz], b+a, dz/10) == 1) // Barriers
            x[iz] = x_max;
        else
            x[iz] = x_min+gsl_pow_2(z[iz]-(b+a/2))*(x_max-x_min)/gsl_pow_2(a/2);
    }

    write_table("x.r", z, x);

    return EXIT_SUCCESS;
}
