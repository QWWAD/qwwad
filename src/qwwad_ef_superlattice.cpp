/**
 * \file    qwwad_ef_superlattice.cpp
 * \brief   Calculate the energy levels for a Kronig-Penney superlattice
 * \author  Alex Valavanis <a.valavanis@leeds.ac.uk>
 */

#include <iostream>
#include <cstdlib>
#include <valarray>
#include "qwwad/options.h"
#include "qwwad/constants.h"
#include "qwwad/file-io.h"
#include "qwwad/schroedinger-solver-kronig-penney.h"

using namespace QWWAD;
using namespace constants;

/**
 * Configure command-line options for the program
 */
Options configure_options(int argc, char* argv[])
{
    Options opt;

    std::string summary("Find the eigenstates of an infinite Kronig-Penney superlattice.");

    opt.add_option<double>("wavevector,k",      0,   "Wave-vector expressed relative to pi/L, where L is period length");
    opt.add_option<double>("wellwidth,a",     100,   "Width of quantum well [angstrom].");
    opt.add_option<double>("barrierwidth,b",  100,   "Width of barrier [angstrom].");
    opt.add_option<double>("wellmass,m",      0.067, "Effective mass in well (relative to that of a free electron)");
    opt.add_option<double>("barriermass,n",   0.067, "Effective mass in barrier (relative to that of a free electron)");
    opt.add_option<bool>  ("outputpotential",        "Output the potential profile for the system to v.r");
    opt.add_option<char>  ("particle,p",       'e',  "ID of particle to be used: 'e', 'h' or 'l', for "
                                                     "electrons, heavy holes or light holes respectively.");
    opt.add_option<size_t>("nz,N",             1000, "Number of spatial points for output file.");
    opt.add_option<size_t>("nst,s",              1,  "Number of states to find");
    opt.add_option<double>("barrierpotential", 100,  "Barrier potential [meV]");
    opt.add_option<double>("Ecutoff",                "Cut-off energy for solutions [meV]");
    
    opt.add_prog_specific_options_and_parse(argc, argv, summary);

    return opt;
}

int main(int argc,char *argv[])
{
    const auto opt = configure_options(argc, argv);

    const auto a    = opt.get_option<double>("wellwidth") * 1e-10;
    const auto b    = opt.get_option<double>("barrierwidth") * 1e-10;
    const auto m_w  = opt.get_option<double>("wellmass") * me;
    const auto m_b  = opt.get_option<double>("barriermass") * me;
    const auto p    = opt.get_option<char>  ("particle");         // particle ID (e, h or l)
    const auto V    = opt.get_option<double>("barrierpotential") * e / 1000;
    const auto nst  = opt.get_option<size_t>("nst");
    const auto k    = opt.get_option<double>("wavevector") * pi/(a+b);   // [1/m]
    const auto N    = opt.get_option<size_t>("nz");               // number of spatial steps

    SchroedingerSolverKronigPenney se(a, b, V, m_w, m_b, k, N, 4, nst);

    // Set cut-off energy if desired
    if(opt.get_argument_known("Ecutoff"))
        se.set_E_cutoff(opt.get_option<double>("Ecutoff") * e/1000);

    // Dump to file
    char energy_filename[9];
    sprintf(energy_filename,"E%c.r",p);

    char wf_prefix[9];
    sprintf(wf_prefix,"wf_%c",p);
    Eigenstate::write_to_file(energy_filename,
                              wf_prefix,
                              ".r",
                              se.get_solutions(true),
                              true);

    // Write potential profile to file if wanted
    if(opt.get_option<bool>("outputpotential"))
        write_table("v.r", se.get_z(), se.get_V());

    return EXIT_SUCCESS;
}
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :
