/**
 * \file   qwwad_ef_poeschl_teller.cpp
 * \brief  Calculate the eigenstates of a Poschl-Teller hole
 * \author Paul Harrison  <p.harrison@shu.ac.uk>
 * \author Alex Valavanis <a.valavanis@leeds.ac.uk>
 */

#include <cstdlib>
#include <iostream>

#include "qwwad/options.h"
#include "qwwad/constants.h"
#include "qwwad/file-io.h"
#include "qwwad/schroedinger-solver-poeschl-teller.h"

using namespace QWWAD;
using namespace constants;

/**
 * Configure command-line options for the program
 */
Options configure_options(int argc, char** argv)
{
    Options opt;

    std::string doc("Generate a Poeschl--Teller potential profile and finds eigenstates analytically.");

    opt.add_option<double>("widthparameter,a", 0.1,   "Width parameter [1/angstrom].");
    opt.add_option<double>("depthparameter,l", 2.0,   "Depth parameter.");
    opt.add_option<double>("length,L",         300,   "Length of potential profile [angstrom].");
    opt.add_option<double>("mass,m",           0.067, "Effective mass (relative to free electron).");
    opt.add_option<size_t>("nz,N",             301,   "Number of spatial points for output file.");
    opt.add_option<size_t>("nstmax",            0,    "Maximum number of subbands to find.  The default (0) means "
                                                      "that all states will be found up to maximum confining potential.");
    opt.add_option<char>  ("particle,p",       'e',   "ID of particle to be used: 'e', 'h' or 'l', for electrons, "
                                                      "heavy holes or light holes respectively.");

    opt.add_prog_specific_options_and_parse(argc, argv, doc);

    return opt;
};

int main(int argc,char *argv[])
{
    const auto opt = configure_options(argc, argv);

    const auto alpha   = opt.get_option<double>("widthparameter") * 1e10;   // Width parameter [1/m]
    const auto lambda  = opt.get_option<double>("depthparameter");         // Depth parameter
    const auto L       = opt.get_option<double>("length") * 1e-10; // Length of potential [m]
    const auto m       = opt.get_option<double>("mass") * me;      // effective mass [kg]
    const auto nz      = opt.get_option<size_t>("nz");             // Number of spatial points for output file
    const auto nst_max = opt.get_option<size_t>("nstmax");
    const auto p       = opt.get_option<char>("particle"); // Particle ID

    if(opt.get_verbose())
    {
        std::cout << "alpha  = " << alpha  << " m^{-1}" << std::endl;
        std::cout << "lambda = " << lambda << std::endl;
        std::cout << "mass   = " << m      << " kg" << std::endl;
    }

    SchroedingerSolverPoeschlTeller se(alpha, lambda, L, m, nz, nst_max);

    // Dump to file
    std::ostringstream energy_filename;
    energy_filename << "E" << p << ".r";

    std::ostringstream wf_prefix;
    wf_prefix << "wf_" << p;

    Eigenstate::write_to_file(energy_filename.str(),
                              wf_prefix.str(),
                              ".r",
                              se.get_solutions(true),
                              true);

    write_table("v.r", se.get_z(), se.get_V());
    
    return EXIT_SUCCESS;
}
