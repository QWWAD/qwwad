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
    opt.add_option<std::string>("style",                "pd",  "Style of plot: 'pd' = probability density, "
                                                               "'real' = real part of wave functions, "
                                                               "'imag' = imaginary part of wave functions.");
    opt.add_option<bool>       ("scalebynstates",              "Scale the wavefunctions by the number of states");
    opt.add_option<bool>       ("showbaseline",                "Show an extra baseline at the energy of the state");

    opt.add_prog_specific_options_and_parse(argc, argv, summary);

    return opt;
}

/** 
 * \brief Find scaling factor to fit wavefunctions nicely on plot 
 *
 * \param[in] states The set of states in the system
 * \param[in] V      The potential profile
 * \param[in] style  The type of plot to generate (PD, real, imag)
 * \param[in] scalebynstates Scale by the number of states if true
 *
 * \returns The factor by which to scale probability densities in plot
 *
 * \details The scaling factor, \f$\sigma\f$ is calculated as
 *            \f$ \sigma = \frac{V_{max} - V_{min}}{N|\psi|^2_{max}} \f$,
 *          i.e., the range of potentials in the plot, divided by the
 *          number of states and the maximum probability density.
 */
static auto scaling_factor(const std::vector<Eigenstate> &states,
                           const arma::vec               &V,
                           const std::string             &style,
                           bool                           scalebynstates) -> double
{
    double scale = V.max() - V.min();

    // If the potential range is zero, use the state energy range instead
    if (scale <= 0) {
        scale = states[states.size()-1].get_energy() - states[0].get_energy();
    }

    if (style == "pd") {
        // Scale by the maximum probability density
        scale /= (Eigenstate::psi_squared_max(states) * 5);
    } else if (style == "real" || style == "imag") {
        // Scale by the square root of this
        scale /= (sqrt(Eigenstate::psi_squared_max(states)) * 10);
    }

    // Scale by number of states if desired
    if(scalebynstates) {
        scale /= states.size();
    }

    return scale;
}

/**
 * \brief Outputs scaled plot of probability densities in
 * 	  quantum well system
 * \todo  Add option to to control how much/if any of the wavefunction tales to cut off
 */
static void output_plot(const WfOptions               &opt,
                        const std::vector<Eigenstate> &states,
                        const arma::vec               &V,
                        const arma::vec               &z)
{
    const auto scalebynstates = opt.get_option<bool>("scalebynstates");
    const auto showbaseline = opt.get_option<bool>("showbaseline");
    const auto nst_max = opt.get_option<size_t>     ("nstmax");
    const auto style   = opt.get_option<std::string>("style");
    const auto plot_file = opt.get_option<std::string>("plotfile");

    const auto dz    = z[1] - z[0];
    const auto nz    = V.size();

    const auto scale = scaling_factor(states, V, style, scalebynstates);

    // Open plot file
    std::ofstream plot_stream(plot_file);

    if(!plot_stream) {
        std::ostringstream oss;
        oss << "Cannot create plot output file " << plot_file << std::endl;
        exit(EXIT_FAILURE);
    }

    // Output conduction band profile
    for(unsigned int iz=0; iz < nz; iz++) {
        plot_stream << z[iz]*1e10 << "\t" << V[iz]/(1e-3*e) << std::endl;
    }

    unsigned int nst_plotted=0; // Counter to limit number of plotted states

    // Output the probability densities
    for(const auto st : states) {
        const auto PD  = st.get_PD(); // Probability density at each point

        // Real and imaginary parts of wavefunction at each point
        const arma::vec psi_real = real(st.get_wavefunction_samples());
        const arma::vec psi_imag = imag(st.get_wavefunction_samples());

        arma::vec plot_data; // The function to be plotted (either PD, real or imag part)

        if(style == "pd") {
            plot_data = PD;
        } else if (style == "real") {
            plot_data = psi_real;
        } else if (style == "imag") {
            plot_data = psi_imag;
        } else {
            std::ostringstream oss;
            oss << "Unrecognised plot style: " << style << std::endl;
            throw std::runtime_error(oss.str());
        }

        if(nst_plotted < nst_max) {
            plot_stream << std::endl; // Separate each PD plot by a blank line
            const auto E = st.get_energy();

            double P_left = 0.0; // probability of electron being found on left of a point

            if (showbaseline) {
                for(unsigned int iz = 0; iz < nz; iz++) {
                    plot_stream << z[iz]*1e10 << "\t" << E/(1e-3*e) << std::endl;
                }

                plot_stream << std::endl;
            }

            // Plot scaled probability density
            for(unsigned int iz = 0; iz < nz; iz++) {
                P_left += PD[iz]*dz;

                // Only plot the bits of the wavefunction with appreciable
                // amplitude
                // TODO: Make this configurable
                if(P_left>0.0001 && P_left<0.9999) {
                    plot_stream << z[iz]*1e10 << "\t" << (plot_data[iz]*scale + E)/(1e-3*e) << std::endl;
                }
            }

            ++nst_plotted;
        }
    }
}

auto main(int argc, char* argv[]) -> int
{
    const auto opt = configure_options(argc, argv);

    const auto states = Eigenstate::read_from_file(opt.get_energy_filename(),
                                                   opt.get_wf_prefix(),
                                                   opt.get_wf_ext(),
                                                   1000.0/e,
                                                   true);

    if(opt.get_verbose()) {
        std::cout << "Read data for " << states.size() << " wave functions" << std::endl;
    }

    arma::vec V;
    arma::vec z;
    const auto totalpotentialfile = opt.get_option<std::string>("totalpotentialfile");
    read_table(totalpotentialfile.c_str(), z, V);

    const auto nz_st = states[0].get_position_samples().size();

    if(V.size() != nz_st)
    {
        std::ostringstream oss;
        oss << "Different number of data points in " << totalpotentialfile
            << " (" << V.size() << " lines) and " << opt.get_wf_filename(1)
            << " (" << nz_st << " lines).  Do these correspond to the same eigenvalue problem?";
        throw std::length_error(oss.str());
    }

    output_plot(opt, states, V, z);

    return EXIT_SUCCESS;
}
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :
