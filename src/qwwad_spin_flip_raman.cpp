/**
 * \file   qwwad_spin_flip_raman.cpp
 * \brief  Calculate spin-flip Raman spectrum
 * \author Paul Harrison  <p.harrison@shu.ac.uk>
 * \author Alex Valavanis <a.valavanis@leeds.ac.uk>
 */

#include <cstdio>
#include <cstdlib>
#include <strings.h>
#include <cmath>
#include <gsl/gsl_randist.h>
#include "struct.h"
#include "maths.h"
#include "qwwad/constants.h"
#include "qwwad/file-io.h"
#include "qwwad/options.h"

using namespace QWWAD;
using namespace constants;

/**
 * \brief Configure command-line options for the program
 */
Options configure_options(int argc, char* argv[])
{
    Options opt;

    std::string doc("Find the spin-flip Raman spectrum.");

    opt.add_option<double>     ("linewidth,l",             1, "Linewidth of all transitions [cm^{-1}].");
    opt.add_option<double>     ("wavenumbermin,s",         0, "Minimum wavenumber for spectrum [cm^{-1}].");
    opt.add_option<double>     ("wavenumbermax,u",       100, "Maximum wavenumber for spectrum [cm^{-1}].");
    opt.add_option<double>     ("wavenumberstep,t",        1, "Step in wavenumber for spectrum [cm^{-1}].");
    opt.add_option<std::string>("spinflipfile",     "e_sf.r", "Table of spin-flip energies vs donor position.");
    opt.add_option<std::string>("spectrumfile",        "I.r", "Filename for output spectrum.");

    opt.add_prog_specific_options_and_parse(argc, argv, doc);

    return opt;
};

int main(int argc, char *argv[])
{
    const auto opt = configure_options(argc, argv);

    // Spectral parameters [1/cm]
    const auto E_min     = opt.get_option<double>("wavenumbermin");
    const auto E_max     = opt.get_option<double>("wavenumbermax");
    const auto E_step    = opt.get_option<double>("wavenumberstep");
    const auto linewidth = opt.get_option<double>("linewidth");

    // Read spin-flip energies as a function of donor position
    std::valarray<double> r_d;  // donor positions
    std::valarray<double> E_sf; // spin-flip energy for each r_d
    const auto spinflipfile = opt.get_option<std::string>("spinflipfile");
    read_table(spinflipfile, r_d, E_sf);

    E_sf *= 1e-3*e/(h*c)/100.0; // Convert from meV to cm^{-1}

    const auto N_rd=r_d.size(); // Number of donor positions

    // Standard deviation from linewidths (FWHM)
    const auto sigma=linewidth/(2*sqrt(2*log(2)));

    std::vector<double> E_plot;
    std::vector<double> intensity_plot; // intensity of Raman signal at Ei

    for(auto E=E_min; E < E_max; E += E_step) // Spectral energy
    {
        auto intensity=0.0;

        for(unsigned int i_i=0;i_i<N_rd;i_i++)
            intensity += gsl_ran_gaussian_pdf(E-E_sf[i_i], sigma);

        E_plot.push_back(E);
        intensity_plot.push_back(intensity);
    }

    const auto spectrumfile = opt.get_option<std::string>("spectrumfile");
    write_table(spectrumfile, E_plot, intensity_plot);

    return EXIT_SUCCESS;
}
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :
