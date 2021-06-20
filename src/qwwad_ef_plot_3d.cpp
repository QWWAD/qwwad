/**
 * \file   qwwad_ef_plot.cpp
 * \author Alex Valavanis <a.valavanis@leeds.ac.uk>
 * \brief  Generates a plot file for wavefunction data
 */

#if HAVE_CONFIG_H
# include "config.h"
#endif

#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <iostream>
#include <gsl/gsl_math.h>

#include "qwwad/constants.h"
#include "qwwad/eigenstate.h"
#include "qwwad/file-io.h"
#include "qwwad/wf_options.h"

using namespace QWWAD;
using namespace constants;

/**
 * Configure command-line options for the program
 */
auto configure_options(int argc, char** argv) -> WfOptions
{
    WfOptions opt;

    std::string summary("Translate wavefunction data into a prettier plottable form.");

    opt.add_option<std::string>("totalpotentialfile", "v.r",   "Name of file containing the total confining potential.");
    opt.add_option<std::string>("plotfile",           "vwf.r", "Name of file to which plottable data will be written.");
    opt.add_option<size_t>     ("nstmax",                10,   "Maximum number of states to plot.");
    opt.add_option<std::string>("style",                "pd",  "Style of plot: 'pd' = probability density, 'wf' = wave functions.");
    opt.add_option<bool>       ("scalebynstates",              "Scale the wavefunctions by the number of states");

    opt.add_prog_specific_options_and_parse(argc, argv, summary);

    return opt;
}

auto main(int argc, char* argv[]) -> int
{
    const auto opt = configure_options(argc, argv);

    // Read the eigenstates
    const auto states = Eigenstate::read_from_file(opt.get_energy_filename(),
                                                   opt.get_wf_prefix(),
                                                   opt.get_wf_ext(),
                                                   1000.0/e,
                                                   true);

    if(opt.get_verbose()) {
        std::cout << "Read data for " << states.size() << " wave functions" << std::endl;
    }

    // Read the potential profile
    arma::vec V;
    arma::vec z;
    const auto totalpotentialfile = opt.get_option<std::string>("totalpotentialfile");
    read_table(totalpotentialfile.c_str(), z, V);

    const auto nz_st = states[0].get_position_samples().size();

    if(V.size() != nz_st) {
        std::ostringstream oss;
        oss << "Different number of data points in " << totalpotentialfile
            << " (" << V.size() << " lines) and " << opt.get_wf_filename(1)
            << " (" << nz_st << " lines).  Do these correspond to the same eigenvalue problem?";
        throw std::length_error(oss.str());
    }

    // Find minimum and maximum energies for plot
    auto E_min = states.front().get_energy();

    if(V.min() < E_min) {
        E_min = V.min();
    }
    
    auto E_max = states.back().get_energy();

    if(V.max() > E_max) {
        E_max = V.max();
    }

    E_min -= 0.1*(E_max - E_min);
    E_max += 0.1*(E_max - E_min);

    const auto nE = 500; // Number of energy samples to plot
    const auto E_range = E_max - E_min;
    const auto dE = E_range/nE;

    arma::mat plotdata = arma::zeros(nE, nz_st);

    // Add probability densities to the plot
    for (auto const &state : states) {
        const auto E = state.get_energy();
        const auto index_E = int((E-E_min)/dE);
        auto PD = state.get_PD();
        PD /= PD.max();

        plotdata.row(index_E) += PD.t();
    }

    for (unsigned iz = 0; iz != z.size(); ++iz) {
        const auto E = V(iz);
        const auto index_E = int((E-E_min)/dE);
        plotdata(index_E, iz) = 1.5;
    }

    std::ofstream mapfile("vwf.xyz");
    mapfile << plotdata;

    // Generate MATLAB script to plot datafile
    std::ofstream gpfile("plotfile.m");
    gpfile << "figure;\n"
           << "imgdata = load('vwf.xyz');\n"
           << "xscale = linspace(" << z.min()*1e10 << "," << z.max()*1e10 << "," << z.size() << ");\n"
           << "yscale = linspace(" << E_min*1000/e << "," << E_max*1000/e << "," << nE       << ");\n"
           << "colormap hot;\n"
           << "surf(xscale, yscale, imgdata, 'EdgeColor', 'None', 'FaceColor', 'interp');\n"
           << "view(2);\n"
           << "xlabel('Position [Angstrom]');\n"
           << "ylabel('Energy [meV]');\n"
           << "\n"
           << "Vdata = load('v.r');\n"
           << "z = Vdata(:,1) * 1e10;\n"
           << "V = Vdata(:,2) * 1000/1.6e-19;\n"
           << "hold on\n"
           << "plot(z,V, 'LineWidth', 2)\n";

    return EXIT_SUCCESS;
}
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :
