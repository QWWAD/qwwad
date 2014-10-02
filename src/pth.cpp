/**
 * \file   pth.cpp
 * \brief  Calculate the eigenstates of a Poschl-Teller hole
 * \author Paul Harrison  <p.harrison@shu.ac.uk>
 * \author Alex Valavanis <a.valavanis@leeds.ac.uk>
 */

#include <cstdlib>

#include "qwwad-options.h"
#include "qclsim-constants.h"
#include "qclsim-fileio.h"
#include "qwwad-schroedinger-poeschl-teller.h"

using namespace Leeds;
using namespace constants;

/**
 * Configure command-line options for the program
 */
Options configure_options(int argc, char* argv[])
{
    Options opt;

    std::string doc("Generate a Poeschl--Teller potential profile and finds eigenstates analytically.");

    opt.add_numeric_option("alpha,a",    0.1,   "Width parameter [1/angstrom].");
    opt.add_numeric_option("lambda,l",   2.0,   "Depth parameter.");
    opt.add_numeric_option("length,L",   300,   "Length of potential profile [angstrom].");
    opt.add_numeric_option("mass,m",     0.067, "Effective mass (relative to free electron).");
    opt.add_size_option   ("nz,N",       301,   "Number of spatial points for output file.");
    opt.add_size_option   ("nst-max",     0,    "Maximum number of subbands to find.  The default (0) means "
                                                "that all states will be found up to maximum confining potential.");
    opt.add_char_option   ("particle,p", 'e',   "ID of particle to be used: 'e', 'h' or 'l', for electrons, "
                                                "heavy holes or light holes respectively.");

    opt.add_prog_specific_options_and_parse(argc, argv, doc);

    return opt;
};

int main(int argc,char *argv[])
{
    const Options opt = configure_options(argc, argv);

    const double alpha   = opt.get_numeric_option("alpha") * 1e10;   // Width parameter [1/m]
    const double lambda  = opt.get_numeric_option("lambda");         // Depth parameter
    const double L       = opt.get_numeric_option("length") * 1e-10; // Length of potential [m]
    const double m       = opt.get_numeric_option("mass") * me;      // effective mass [kg]
    const size_t nz      = opt.get_size_option("nz"); // Number of spatial points for output file
    const size_t nst_max = opt.get_size_option("nst-max");
    const char   p       = opt.get_char_option("particle"); // Particle ID

    if(opt.get_verbose())
    {
        std::cout << "alpha  = " << alpha  << " m^{-1}" << std::endl;
        std::cout << "lambda = " << lambda << std::endl;
        std::cout << "mass   = " << m      << " kg" << std::endl;
    }

    SchroedingerSolverPoeschlTeller se(alpha, lambda, L, m, nz, nst_max);

    // Dump to file
    char energy_filename[9];
    sprintf(energy_filename,"E%c.r",p);

    char wf_prefix[9];
    sprintf(wf_prefix,"wf_%c",p);

    State::write_to_file(energy_filename,
                         wf_prefix,
                         ".r",
                         se.get_solutions(true),
                         se.get_z(),
                         true);

    Leeds::write_table_xy("v.r", se.get_z(), se.get_V());
    
    return EXIT_SUCCESS;
}
