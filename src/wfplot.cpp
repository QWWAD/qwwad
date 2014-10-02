/**
 * \file   wfplot.cpp
 * \author Alex Valavanis <a.valavanis@leeds.ac.uk>
 * \date   2012-01-30
 * \brief  Generates a plot file for wavefunction data
 */

#if HAVE_CONFIG_H
# include "config.h"
#endif

#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <gsl/gsl_math.h>

#include "qclsim-constants.h"
#include "qclsim-linalg.h"
#include "wf_options.h"

using namespace Leeds;
using namespace Leeds::constants;

/**
 * Configure command-line options for the program
 */
WfOptions configure_options(int argc, char* argv[])
{
    WfOptions opt(WF_OPTION_MODE_IN);

    std::string summary("Translate wavefunction data into a prettier plottable form.");

    opt.add_string_option("plot-file", "vwf.r", "Name of file to which plottable data will be written.");
    opt.add_size_option  ("nst-max",   10,      "Maximum number of states to plot.");
    opt.add_switch       ("plot-wf",            "Plot wavefunctions rather than probability density.");

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
static double scaling_factor(const std::vector<State>    &states,
                             const std::valarray<double> &V)
{
    double Vrange = V.max() - V.min();

    // If the potential range is zero, use the state energy range instead
    if (Vrange <= 0)
        Vrange = states[states.size()-1].get_E() - states[0].get_E();

    return Vrange/(State::psi_squared_max(states)*states.size());
}

/**
 * \brief Outputs scaled plot of probability densities in
 * 	  quantum well system
 * \todo  Add option to to control how much/if any of the wavefunction tales to cut off
 */
static void output_plot(const WfOptions             &opt,
                        const std::vector<State>    &states,
                        const std::valarray<double> &V,
                        const std::valarray<double> &z)
{
    const double dz    = z[1] - z[0];
    const double scale = scaling_factor(states, V);
    const size_t nz    = V.size();
    std::string plot_file = opt.get_string_option("plot-file");

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

    const size_t nst_max = opt.get_size_option("nst-max");
    const bool   plot_wf = opt.get_switch("plot-wf");

    // Output the probability densities
    for(std::vector<State>::const_iterator ist = states.begin(); ist != states.end(); ++ist)
    {
        if(nst_plotted < nst_max)
        {
            fprintf(plot_stream, "\n"); // Separate each PD plot by a blank line
            const std::valarray<double> PD  = ist->psi_squared(); // Probability density at each point
            const std::valarray<double> psi = ist->psi_array(); // Wavefunction

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
                    if (plot_wf)
                    {
                        double scale_wf = (V.max()-V.min())/((psi.max() - psi.min()) * states.size() * 2);
                        fprintf(plot_stream, "%e\t%e\n", z[iz]*1e10,
                                (psi[iz]*scale_wf + ist->get_E())/(1e-3*e));
                    }
                    else
                        fprintf(plot_stream, "%e\t%e\n", z[iz]*1e10,
                                (PD[iz]*scale + ist->get_E())/(1e-3*e));
                }
            }

            ++nst_plotted;
        }
    }

    fclose(plot_stream);
}

int main(int argc, char* argv[])
{
    const WfOptions opt = configure_options(argc, argv);

    std::vector<State> states = State::read_from_file(opt.get_energy_input_path(),
                                                      opt.get_wf_input_prefix(),
                                                      opt.get_wf_input_ext(),
                                                      1000.0/e,
                                                      true);

    std::valarray<double> V;
    std::valarray<double> z;    
    read_table_xy(opt.get_potential_input_path().c_str(), z, V);

    if(V.size() != states[0].size())
    {
        std::ostringstream oss;
        oss << "Different number of data points in " << opt.get_potential_input_path()
            << " (" << V.size() << " lines) and " << opt.get_wf_input_path(1)
            << " (" << states[0].size() << " lines).  Do these correspond to the same eigenvalue problem?";
        throw std::length_error(oss.str());
    }

    output_plot(opt, states, V, z);

    return EXIT_SUCCESS;
}
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :
