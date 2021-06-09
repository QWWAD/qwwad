/**
 * \file    qwwad_ef_square_well.cpp
 * \brief   Calculate the 1-particle energies and wavefunctions in a single quantum well
 * \author  Paul Harrison  <p.harrison@shu.ac.uk>
 * \author  Alex Valavanis <a.valavanis@leeds.ac.uk>
 */

#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <valarray>
#include "qwwad/constants.h"
#include "qwwad/schroedinger-solver-finite-well.h"
#include "qwwad/options.h"
#include "qwwad/file-io.h"

using namespace QWWAD;
using namespace constants;

/**
 * Configure command-line options for the program
 */
auto configure_options(int argc, char** argv) -> Options
{
    Options opt;

    std::string summary("Find the eigenstates of a single finite quantum well. ");

    opt.add_option<double>("wellwidth,a",     100,   "Width of quantum well [angstrom].");
    opt.add_option<double>("barrierwidth,b",  200,   "Width of barrier [angstrom]. Note that this is only used "
                                                     "for the purposes of outputting the data. The calculation here "
                                                     "assumes that the barriers are infinitely thick.  As such, the "
                                                     "wavefunctions do not decay to precisely zero at the boundaries.");
    opt.add_option<double>("wellmass,m",      0.067, "Effective mass in well (relative to that of a free electron)");
    opt.add_option<double>("barriermass,n",   0.067, "Effective mass in barrier (relative to that of a free electron)");
    opt.add_option<bool>  ("outputequations",        "Output the matching equations for the system. The left-hand "
                                                     "side of the equation is output to 'lhs.r' and each branch of "
                                                     "the right-hand side is output to a set of 'rhs_i.r' files, "
                                                     "where 'i' is the index of the state that lies on that branch. "
                                                     "An rhs file is output for all the bound states in the system and "
                                                     "one additional branch (with no real solution)");
    opt.add_option<bool>  ("outputpotential",        "Output the potential profile for the system to v.r");
    opt.add_option<char>  ("particle,p",       'e',  "ID of particle to be used: 'e', 'h' or 'l', for electrons, "
                                                     "heavy holes or light holes respectively.");
    opt.add_option<size_t>("nz,N",             1000, "Number of spatial points for output file.");
    opt.add_option<size_t>("nst,s",              1,  "Number of states to find");
    opt.add_option<double>("barrierpotential",  100, "Barrier potential [meV]");
    opt.add_option<double>("Emin",                   "Lower cut-off energy for solutions [meV]");
    opt.add_option<double>("Emax",                   "Upper cut-off energy for solutions [meV]");

    opt.add_prog_specific_options_and_parse(argc, argv, summary);

    return opt;
}

auto main(int argc,char *argv[]) -> int
{
    const auto opt   = configure_options(argc, argv);
    const auto a     = opt.get_option<double>("wellwidth") * 1e-10;
    const auto b     = opt.get_option<double>("barrierwidth") * 1e-10;
    const auto m_w   = opt.get_option<double>("wellmass") * me;
    const auto m_b   = opt.get_option<double>("barriermass") * me;
    const auto p     = opt.get_option<char>  ("particle");         // particle ID (e, h or l)
    const auto V     = opt.get_option<double>("barrierpotential") * e / 1000;
    const auto state = opt.get_option<size_t>("nst");
    const auto N     = opt.get_option<size_t>("nz");               // number of spatial steps

    SchroedingerSolverFiniteWell se(a, b, V, m_w, m_b, N, state);

    // Set cut-off energies if desired
    if(opt.get_argument_known("Emin")) {
        se.set_E_min(opt.get_option<double>("Emin") * e/1000);
    }

    if(opt.get_argument_known("Emax")) {
        se.set_E_max(opt.get_option<double>("Emax") * e/1000);
    }

    if(opt.get_option<bool>("outputequations")) {
        const auto nst    = se.get_n_bound();
        const auto v_max  = (nst+1)*pi/2;

        const auto nv = (nst+1)*1000;
        const auto dv = v_max/(nv-1);
        std::valarray<double> v(nv);
        std::valarray<double> lhs(nv);

        // Output lhs data in a single contiguous chunk
        for (unsigned int iv = 0; iv < nv; ++iv) {
            v[iv]   = iv*dv;
            lhs[iv] = se.get_lhs(v[iv]);
        }

        write_table("lhs.r", v, lhs);

        // Output a separate file for each rhs branch
        for (unsigned int ibranch = 0; ibranch < nst+1; ++ibranch) {
            const size_t nv_branch = 1000;

            std::valarray<double> v_branch(nv_branch);
            std::valarray<double> rhs(nv_branch);

            // Set the spacing between points so that we don't quite reach the
            // asymptote at the "end" of the branch
            const auto dv_branch = (pi/2.0*0.999999)/(nv_branch-1);

            for (unsigned int iv_branch = 0; iv_branch < nv_branch; ++iv_branch) {
                v_branch[iv_branch] = pi/2.0 * ibranch  + iv_branch*dv_branch;
                rhs[iv_branch]      = SchroedingerSolverFiniteWell::get_rhs(v_branch[iv_branch]);
            }

            // Set filename
            std::ostringstream oss;
            oss << "rhs_" << ibranch+1 << ".r";
            write_table(oss.str().c_str(), v_branch, rhs);
        }
    }

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

    // Write potential profile to file if wanted
    if(opt.get_option<bool>("outputpotential")) {
        write_table("v.r", se.get_z(), se.get_V());
    }

    return EXIT_SUCCESS;
}
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :
