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
#include <gsl/gsl_sf_hyperg.h>
#include "qwwad-options.h"
#include "qclsim-constants.h"
#include "qclsim-fileio.h"
#include "qclsim-linalg.h"

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

    const std::valarray<double> sinh_alpha_z = sinh(alpha*z);
    const std::valarray<double> _x = -pow(sinh_alpha_z,2.0);

    write_table_xy("v.r", z, V); // Write potential profile to file

    // Number of bound states
    size_t nst = ceil(lambda-1); // principal quantum number

    std::vector<State> _solutions;

    // Compute solutions analytically
    for(unsigned int ist = 0; ist < nst; ++ist)
    {
        // Energy is found using [QWWAD4, 3.60]
        const double kappa = alpha * (lambda-1-ist);
        const double E = -gsl_pow_2(hBar * kappa) / (2.0*m);

        // Wavefunction is taken from "Practical Quantum Mechanics", Flugge (1970).
        // The solution is a hypergeometric function, whose arguments and scaling
        // factor depend on whether the state has odd or even parity.
        
        double arg1, arg2, arg3;    // Arguments for hypergeometric function
        std::valarray<double> fact; // Prefactor for hypergeometric function

        // Flugge, 39.24
        const double a = 0.5 * (lambda - kappa/alpha);
        const double b = 0.5 * (lambda + kappa/alpha);

        if(ist % 2) // Odd-parity states
        {
            // From Flugge, 39.10b
            arg1 = a+0.5;
            arg2 = b+0.5;
            arg3 = 1.5;
            fact = pow(cosh(alpha*z),lambda) * sinh_alpha_z;
        }
        else // Even-parity states
        {
            // From Flugge, 39.10a
            arg1 = a;
            arg2 = b;
            arg3 = 0.5;
            fact = pow(cosh(alpha*z),lambda);
        }

        std::valarray<double> psi(nz); // Wavefunction amplitude at each point [m^{-0.5}]

        for(unsigned int iz = 0; iz < nz; ++iz)
        {
            if(abs(_x[iz]) < 1)
            {
                psi[iz] = fact[iz] *
                          gsl_sf_hyperg_2F1(arg1, arg2, arg3, _x[iz]);
            }
            // If the argument is too large, we need to apply a linear
            // transformation such that |_x| < 1
            else if(gsl_fcmp(_x[iz]/(_x[iz]-1), 1, 0.0025) == -1)
            {
                psi[iz] = fact[iz] *
                          pow(1-_x[iz],-arg2) *
                          gsl_sf_hyperg_2F1(arg2, arg3-arg1, arg3, _x[iz]/(_x[iz]-1));
            }
            // In case we're *very* close to _x = 1, GSL can't cope, so we
            // need to simplify things further, and just pass a large number
            // as the argument.
            // This seems to be OK, but might need a little investigation
            else
            {
                psi[iz] = fact[iz] *
                          pow(1-_x[iz],-arg2) *
                          gsl_sf_hyperg_2F1(arg2, arg3-arg1, arg3, 0.99);
            }
        }

        _solutions.push_back(State(E * 1000/e, psi));
        _solutions.back().normalise(z);
    }

    // Write energy to file
    char energy_filename[9];
    sprintf(energy_filename,"E%c.r",p);

    char wf_prefix[9];
    sprintf(wf_prefix,"wf_%c",p);

    State::write_to_file(energy_filename,
                         wf_prefix,
                         ".r",
                         _solutions,
                         z,
                         true);

//    write_table_x(filename, E, true);

    return EXIT_SUCCESS;
}
