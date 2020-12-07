/**
 * \file   efiw.cpp Calculates eigenstates of an infinite square well
 *
 * \author Paul Harrison  <p.harrison@shu.ac.uk>
 * \author Alex Valavanis <a.valavanis@leeds.ac.uk>
 */

#include <iostream>
#include <cstdlib>
#include <cmath>
#include <valarray>
#include <gsl/gsl_math.h>
#include "qwwad/constants.h"
#include "qwwad/schroedinger-solver-infinite-well.h"
#include "qwwad/options.h"

using namespace QWWAD;
using namespace constants;

/**
 * \brief Configure command-line options for the program
 */
auto configure_options(int argc, char** argv) -> Options
{
    Options opt;

    std::string doc("Find the eigenstates of an infinite quantum well.");

    opt.add_option<double>("wellwidth,L",   100, "Width of quantum well [angstrom].");
    opt.add_option<double>("mass,m",      0.067, "Effective mass (relative to free electron).");
    opt.add_option<size_t>("nz,N",          100, "Number of spatial points for output file.");
    opt.add_option<size_t>("nst,s",           1, "Number of states to find.");
    opt.add_option<char>  ("particle,p",    'e', "ID of particle to be used: 'e', 'h' or 'l', for "
                                                 "electrons, heavy holes or light holes respectively.");
    opt.add_option<double>("vcb",          0.00, "Band-edge potential [eV]");
    opt.add_option<double>("alpha",        0.00, "Non-parabolicity parameter [eV^{-1}]");
    opt.add_option<double>("Emin",               "Lower cut-off energy for solutions [meV]");
    opt.add_option<double>("Emax",               "Upper cut-off energy for solutions [meV]");
    opt.add_option<double>("barrierwidth", 0.00, "Width of barriers [angstrom]");

    opt.add_prog_specific_options_and_parse(argc, argv, doc);

    return opt;
};

auto main(int argc, char *argv[]) -> int
{
    const auto opt = configure_options(argc, argv);

    const auto L     = opt.get_option<double>("wellwidth") * 1e-10;    // well width [m]
    const auto Lb    = opt.get_option<double>("barrierwidth") * 1e-10; // barrier width [m]
    const auto p     = opt.get_option<char>  ("particle");             // particle ID (e, h or l)
    const auto m     = opt.get_option<double>("mass") * me;            // effective mass [kg]
    const auto N     = opt.get_option<size_t>("nz");                   // number of spatial steps
    const auto s     = opt.get_option<size_t>("nst");                  // number of states
    const auto alpha = opt.get_option<double>("alpha") / e;            // Non-parabolicity [1/J]
    const auto V     = opt.get_option<double>("vcb") * e;              // band_edge potential [J]

    SchroedingerSolverInfWell se(m, L, N, alpha, V, s);
    se.set_padding_width(Lb);

    // Set cut-off energies if desired
    if(opt.get_argument_known("Emin"))
    {
        se.set_E_min(opt.get_option<double>("Emin") * e/1000);
    }
    if(opt.get_argument_known("Emax"))
    {
        se.set_E_max(opt.get_option<double>("Emax") * e/1000);
    }

    const auto solutions = se.get_solutions(true);

    // Dump to file
    std::ostringstream energy_filename;
    energy_filename << "E" << p << ".r";

    std::ostringstream wf_prefix;
    wf_prefix << "wf_" << p;
    Eigenstate::write_to_file(energy_filename.str(),
                              wf_prefix.str(),
                              ".r",
                              solutions,
                              true);

    return EXIT_SUCCESS;
}
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :
