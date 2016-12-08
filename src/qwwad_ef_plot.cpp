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
WfOptions configure_options(int argc, char* argv[])
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

/** 
 * \brief Find scaling factor to fit wavefunctions nicely on plot 
 *
 * \param[in] states The set of states in the system
 * \param[in] V      The potential profile
 *
 * \returns The factor by which to scale probability densities in plot
 *
 * \details The scaling factor, \f$\sigma\f$ is calculated as
 *            \f$ \sigma = \frac{V_{max} - V_{min}}{N|\psi|^2_{max}} \f$,
 *          i.e., the range of potentials in the plot, divided by the
 *          number of states and the maximum probability density.
 */
static double scaling_factor(const std::vector<Eigenstate> &states,
                             const arma::vec               &V,
                             const bool                     scalebynstates)
{
    double scale = V.max() - V.min();

    // If the potential range is zero, use the state energy range instead
    if (scale <= 0)
    {
        scale = states[states.size()-1].get_energy() - states[0].get_energy();
    }

    // Scale by the maximum probability density
    scale /= (Eigenstate::psi_squared_max(states) * 5);

    // Scale by number of states if desired
    if(scalebynstates)
    {
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
    const auto dz    = z[1] - z[0];
    const auto scalebynstates = opt.get_option<bool>("scalebynstates");
    const auto scale = scaling_factor(states, V, scalebynstates);
    const auto nz    = V.size();
    const auto plot_file = opt.get_option<std::string>("plotfile");

    // Open plot file
    FILE* plot_stream = fopen(plot_file.c_str(), "w");
    if(plot_stream==NULL)
    {
        std::ostringstream oss;
        oss << "Cannot create plot output file " << plot_file << std::endl;
        exit(EXIT_FAILURE);
    }

    // Output conduction band profile
    for(unsigned int iz=0; iz < nz; iz++)
        fprintf(plot_stream, "%e\t%e\n", z[iz]*1e10, V[iz]/(1e-3*e));

    unsigned int nst_plotted=0; // Counter to limit number of plotted states

    const auto nst_max = opt.get_option<size_t>     ("nstmax");
    const auto style   = opt.get_option<std::string>("style");

    // Output the probability densities
    for(const auto st : states)
    {
        if(nst_plotted < nst_max)
        {
            fprintf(plot_stream, "\n"); // Separate each PD plot by a blank line
            const auto PD  = st.get_PD(); // Probability density at each point
            const auto psi = st.get_wavefunction_samples(); // Wavefunction

            double P_left = 0.0; // probability of electron being found on left of a point

            // Plot scaled probability density
            for(unsigned int iz = 0; iz < nz; iz++)
            {
                P_left += PD[iz]*dz;

                // Only plot the bits of the wavefunction with appreciable
                // amplitude
                // TODO: Make this configurable
                if(P_left>0.0001 && P_left<0.9999)
                {
                    const auto E = st.get_energy();

                    if (style == "wf")
                    {
                        double scale_wf = (V.max()-V.min())/((psi.max() - psi.min()) * states.size() * 2);
                        fprintf(plot_stream, "%e\t%e\n", z[iz]*1e10,
                                (psi[iz]*scale_wf + E)/(1e-3*e));
                    }
                    else
                        fprintf(plot_stream, "%e\t%e\n", z[iz]*1e10,
                                (PD[iz]*scale + E)/(1e-3*e));
                }
            }

            ++nst_plotted;
        }
    }

    fclose(plot_stream);
}

int main(int argc, char* argv[])
{
    const auto opt = configure_options(argc, argv);

    const auto states = Eigenstate::read_from_file(opt.get_energy_filename(),
                                                   opt.get_wf_prefix(),
                                                   opt.get_wf_ext(),
                                                   1000.0/e,
                                                   true);

    if(opt.get_verbose())
        std::cout << "Read data for " << states.size() << " wave functions" << std::endl;

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
