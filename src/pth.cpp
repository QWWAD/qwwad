/**
 * \file   pth.cpp
 * \brief  Calculate the eigenstates of a Poschl-Teller hole
 * \author Paul Harrison  <p.harrison@shu.ac.uk>
 * \author Alex Valavanis <a.valavanis@leeds.ac.uk>
 */

#include <cstdio>
#include <cstdlib>
#include <valarray>

#include <gsl/gsl_math.h>
#include "qwwad-options.h"
#include "qclsim-constants.h"
#include "qclsim-fileio.h"

using namespace Leeds;
using namespace constants;

/**
 * Configure command-line options for the program
 */
Options configure_options(int argc, char* argv[])
{
    Options opt;

    std::string doc("Generate a Poeschl--Teller potential profile and finds eigenstate energies analytically.");

    opt.add_numeric_option("alpha,a",    0.1,   "Width parameter [1/angstrom].");
    opt.add_numeric_option("lambda,l",   2.0,   "Depth parameter.");
    opt.add_numeric_option("length,L",   300,   "Length of potential profile [angstrom].");
    opt.add_numeric_option("mass,m",     0.067, "Effective mass (relative to free electron).");
    opt.add_size_option   ("nz,N",       301,   "Number of spatial points for output file.");
    opt.add_char_option   ("particle,p", 'e',   "ID of particle to be used: 'e', 'h' or 'l', for electrons, heavy holes or light holes respectively.");

    opt.add_prog_specific_options_and_parse(argc, argv, doc);

    return opt;
};

int main(int argc,char *argv[])
{
    const Options opt = configure_options(argc, argv);

    const double alpha  = opt.get_numeric_option("alpha") * 1e10;   // Width parameter [1/m]
    const double lambda = opt.get_numeric_option("lambda");         // Depth parameter
    const double L      = opt.get_numeric_option("length") * 1e-10; // Length of potential [m]
    const double m      = opt.get_numeric_option("mass") * me;      // effective mass [kg]
    const size_t nz     = opt.get_size_option("nz"); // Number of spatial points for output file
    const char   p      = opt.get_char_option("particle"); // Particle ID

    if(opt.get_verbose())
    {
        std::cout << "alpha  = " << alpha  << " m^{-1}" << std::endl;
        std::cout << "lambda = " << lambda << std::endl;
        std::cout << "mass   = " << m      << " kg" << std::endl;
    }

    std::valarray<double> z(nz); // Spatial samples [m]
    std::valarray<double> V(nz); // Potential [J]

    // Generate potential profile
    const double dz = L / (nz -1); // Size of each spatial cell [m]

    for(unsigned int iz = 0; iz < nz; ++iz)
    {
        // Construct array of positions centered at zero
        z[iz] = iz*dz - L/2.0;
        V[iz] = -gsl_pow_2(hBar*alpha/cosh(alpha*z[iz]))*lambda*(lambda-1)/(2*m);
    }

    write_table_xy("v.r", z, V); // Write potential profile to file

    // Number of bound states
    size_t nst = ceil(lambda-1); // principal quantum number

    std::valarray<double> E(nst); // Energy of each state

    // Compute energy of state analytically
    for(unsigned int ist = 0; ist < nst; ++ist)
        E[ist] = -gsl_pow_2(hBar*alpha*(lambda-1-(float)ist))/(2*m);

    E *= 1000/e; // Convert to meV

    // Write energy to file
    char filename[9];
    sprintf(filename,"E%c.r",p);
    write_table_x(filename, E, true);

    return EXIT_SUCCESS;
}
