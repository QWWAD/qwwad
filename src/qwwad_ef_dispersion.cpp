/**
 * \file    qwwad_ef_dispersion_relation.cpp
 * \brief   Prints out the dispersion relations for each subband
 * \author  Alex Valavanis  <a.valavanis@leeds.ac.uk>
 * \author  Jonathan Cooper <jdc.tas@gmail.com>
 */

#include <iostream>
#include <valarray>

#include "qwwad/wf_options.h"
#include "qwwad/subband.h"
#include "qwwad/file-io.h"
#include "qwwad/constants.h"
using namespace QWWAD;
using namespace constants;

/**
 * Configure command-line options for the program
 */
auto configure_options(int argc, char** argv) noexcept -> WfOptions
{
    WfOptions opt;

    std::string summary = "Compute the dispersion relation for a set of subbands.";

    const double nkbt_def = 5.0;   // Default cut-off [multiple of kT]
    const double Te_def   = 100.0; // Default carrier temperature [K]
    const size_t nk_def   = 100;   // Default number of k points
    const double me_GaAs  = 0.067; // Default effective mass (rel. to free electron)

    opt.add_option<double>     ("nkbt",         nkbt_def, "Maximum energy to print out for the subband [multiple of kT].");
    opt.add_option<double>     ("Te",           Te_def,   "Carrier temperature [K].");
    opt.add_option<size_t>     ("nk",           nk_def,   "Number of k-space points to print out");
    opt.add_option<std::string>("disp-prefix",    "dr_e", "Filename prefix to which dispersion curves will be written");
    opt.add_option<std::string>("disp-ext",       ".r",   "Filename extension to which dispersion curves will be written");
    opt.add_option<bool>       ("relative",               "Output dispersion relative to subband minima. If not specified, "
                                                          "the dispersion is given relative to the band edge.");
    opt.add_option<double>     ("mass",         me_GaAs,  "In-plane effective mass (relative to free electron).");
    opt.add_option<double>     ("alpha",          0.0,    "In-plane non-parabolicity parameter [1/eV].");
    opt.add_option<double>     ("vcb",            0.0,    "Conduction band edge [eV].");

    opt.add_prog_specific_options_and_parse(argc, argv, summary);

    return opt;
};

auto main (int argc, char* argv[]) -> int
{
    const auto opt = configure_options(argc, argv);

    const auto nk   = opt.get_option<size_t>("nk");
    const auto nkbt = opt.get_option<double>("nkbt");
    const auto Te   = opt.get_option<double>("Te");

    std::vector<Subband> subbands;

    try {
        subbands = Subband::read_from_file(opt.get_energy_filename(),
                                           opt.get_wf_prefix(),
                                           opt.get_wf_ext(),
                                           opt.get_option<double>("mass") * me,
                                           opt.get_option<double>("alpha") / e,
                                           opt.get_option<double>("vcb") * e);
    } catch (std::runtime_error &e) {
        std::cerr << "Could not read subbands from file " << opt.get_energy_filename() << std::endl;
        exit(EXIT_FAILURE);
    }

    const double J_to_meV = 1e-3*e;

    // Loop over subbands
    unsigned int ist = 1;
    for(auto const sb : subbands)
    {
        std::valarray<double> k(nk);
        std::valarray<double> Ek(nk);

        // Calculate maximum wavevector
        double k_max = 0.0;

        try {
            k_max = sb.get_k_at_Ek(nkbt*kB*Te);
        } catch (std::domain_error &e) {
            std::cerr << e.what();
            exit(EXIT_FAILURE);
        }

        // Calculate wavevector spacing
        const auto dk = k_max/nk;

        // Loop over wavevectors and find corresponding energies
        for(unsigned int ik=0; ik<nk; ik++)
        {
            k[ik]  = ik*dk;

            try {
                Ek[ik] = sb.get_Ek_at_k(k[ik]);
            } catch (std::domain_error &e) {
                std::cerr << e.what();
                exit(EXIT_FAILURE);
            }
        }

        // If absolute energies are required then offset energies by the energy of the subband
        // minima
        if(!opt.get_option<bool>("relative"))
            Ek += sb.get_E_min();

        Ek /= J_to_meV; // Rescale to meV

        // If verbose option selected output some information about subband
        if(opt.get_verbose())
        {
            std::cout << "Subband " << ist << " at " << sb.get_E_min()/J_to_meV << " meV." << std::endl
                      << "D.o.s effective mass: " << sb.get_effective_mass() << std::endl
                      << "D.o.s nonparabolicity parameter: " << sb.get_alpha() << std::endl
                      << "Wavevector at " << nkbt << "*kB*Te: " << k_max << std::endl
                      << std::endl;
        }

        // Construct filename and output
        std::stringstream filename;
        filename << opt.get_option<std::string>("disp-prefix") << ist << opt.get_option<std::string>("disp-ext");

        try {
            write_table(filename.str().c_str(), k, Ek);
        } catch (std::runtime_error &e) {
            std::cerr << e.what();
            exit(EXIT_FAILURE);
        }

        // Increment state counter
        ist++;
    }

    return EXIT_SUCCESS;
}
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :
