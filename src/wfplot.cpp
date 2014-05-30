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
 * \brief User-settings for wavefunction plotter
 */
class WfPlotOptions : public WfOptions
{
    public:
        WfPlotOptions(int argc, char* argv[]);
        std::string get_plot_file() const {return vm["plot-file"].as<std::string>();}
        inline unsigned int get_upper() const {return vm["upper-state"].as<unsigned int>();}
        bool plot_wf() const {return vm["plot-wf"].as<bool>();}

        void print() const {}
};

/**
 * \brief Constructor: Define and parse all user options
 *
 * \param[in] argc Number of command-line arguments
 * \param[in] argv Array of command-line arguments
 */
WfPlotOptions::WfPlotOptions(int argc, char* argv[])
{
    program_specific_options->add_options()
        ("plot-file",
         po::value<std::string>()->default_value("Vwf.dat"),
         "Name of file to which pretty plot of probability densities "
         "will be written.")
        
        ("upper-state,u", po::value<unsigned int>()->default_value(10),
         "Maximum number of states to plot")

        ("plot-wf", po::bool_switch()->default_value(false),
         "Plot the wavefunctions rather than the probability density in each state")
        ;

    static char doc[] = 
        "Generates plot of wavefunction data.";

    add_prog_specific_options_and_parse(argc, argv, doc);

    if(get_verbose()) print();		
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
static void output_plot(const WfPlotOptions         &opt,
                        const std::vector<State>    &states,
                        const std::valarray<double> &V,
                        const std::valarray<double> &z)
{
    const double dz    = z[1] - z[0];
    const double scale = scaling_factor(states, V);
    const size_t nz    = V.size();

    // Open plot file
    FILE* plot_stream = fopen(opt.get_plot_file().c_str(), "w");
    if(plot_stream==NULL)
    {
        std::ostringstream oss;
        oss << "Cannot create plot output file " << opt.get_plot_file().c_str() << std::endl;
        exit(EXIT_FAILURE);
    }

    // Output conduction band profile
    for(unsigned int iz=0; iz < nz; iz++)
        fprintf(plot_stream, "%e\t%e\n", z[iz]*1e10, V[iz]/(1e-3*e));

    unsigned int nst_plotted=0; // Counter to limit number of plotted states

    // Output the probability densities
    for(std::vector<State>::const_iterator ist = states.begin(); ist != states.end(); ++ist)
    {
        if(nst_plotted < opt.get_upper())
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
                    if (opt.plot_wf())
                    {
                        double scale_wf = (V.max()-V.min())/((psi.max() - psi.min()) * states.size());
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
    const WfPlotOptions opt(argc, argv);

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
