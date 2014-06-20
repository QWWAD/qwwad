/*=========================================================
      efpqw 
  =========================================================*/

/* This program produces the function of alloy concentration
   x with displacement z for a parabolic quantum well,
   the output is in a format suitable for conversion into
   electron and hole potentials using efxv.      

   ----------b----------+ a +----------b---------- x_max
                        |   |
                        \___/  x_min

   Paul Harrison,  October 1995                            
   
   Major modifications, 15th May 1998                      */

#include <cstdio>
#include <cstdlib>
#include <gsl/gsl_math.h>
#include "qclsim-fileio.h"
#include "qwwad-options.h"

using namespace Leeds;

/**
 * Configure command-line options for the program
 */
Options configure_options(int argc, char* argv[])
{
    Options opt;

    opt.add_numeric_option("well-width,a",    100, "Width at top of quantum well [angstrom].");
    opt.add_numeric_option("barrier-width,b", 100, "Width of barriers [angstrom].");
    opt.add_size_option   ("nz,N",            301, "Number of spatial points for output file.");
    opt.add_numeric_option("xmin,x",          0,   "Minimum alloy fraction.");
    opt.add_numeric_option("xmax,y",          0.1, "Maximum alloy fraction.");

    std::string doc("Generate a parabolic alloy profile surrounded by thick barriers.");

    std::string details("The following output text files are created:\n"
                        "  'x.r'   \tAlloy at each point:\n"
                        "          \tCOLUMN 1: spatial position [m].\n"
                        "          \tCOLUMN 2: alloy fraction.\n"
                        "\n"
                        "Examples:\n"
                        "   Generate 100-angstrom-wide parabolic-graded alloy with values ranging from 0 to 0.4, surrounded by 200-angstrom barriers:\n\n"
                        "   efpqw --well-width 100 --xmin 0 --xmax 0.4 --barrier-width 200\n"
                        "\n"
                        "   Generate 100-angstrom-wide parabolic-graded alloy with values ranging from 0.1 to 0.3, surrounded by 500-angstrom barriers using 500 points in output file:\n\n"
                        "   efpqw --well-width 100 --xmin 0.1 --xmax 0.3 --barrier-width 500 --nz 500");

    opt.add_prog_specific_options_and_parse(argc, argv, doc, details);

    return opt;
};

int main(int argc,char *argv[])
{
    Options opt = configure_options(argc, argv);

    const double a     = opt.get_numeric_option("well-width") * 1e-10;    // [m]
    const double b     = opt.get_numeric_option("barrier-width") * 1e-10; // [m]
    const size_t nz    = opt.get_size_option("nz");      // Number of points for output file
    const double x_min = opt.get_numeric_option("xmin"); // Minimum alloy fraction
    const double x_max = opt.get_numeric_option("xmax"); // Maximum alloy fraction

    const double dz = (a+2*b)/(nz-1); // Find width of each spatial interval [m]

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

    write_table_xy("x.r", z, x);

    return EXIT_SUCCESS;
}/* end main */



