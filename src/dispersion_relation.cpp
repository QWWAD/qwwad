/**
 * \file    dispersion_relation.cpp
 * \brief   Prints out the dispersion relations for each subband
 * \author  Alex Valavanis  <a.valavanis@leeds.ac.uk>
 * \author  Jonathan Cooper <jdc.tas@gmail.com>
 */

#include <iostream>
#include <valarray>

#include "wf_options.h"
#include "qclsim-subband.h"

using namespace Leeds;

#include "qclsim-constants.h"
using namespace constants;

/**
 * Configure command-line options for the program
 */
WfOptions configure_options(int argc, char* argv[])
{
    WfOptions opt;

    opt.add_numeric_option("nkbt",           5.0,        "Maximum energy to print out for the subband [multiple of kT].");
    opt.add_numeric_option("Te",             100,        "Electron temperature [K].");
    opt.add_size_option   ("nk",             100,        "Number of k-space points to print out");
    opt.add_string_option ("disp-prefix",    "dr_e",     "Filename prefix to which dispersion curves will be written");
    opt.add_string_option ("disp-ext",       ".r",       "Filename extension to which dispersion curves will be written");
    opt.add_switch        ("nonparabolic",               "Use non-parabolic dispersion relation.");
    opt.add_switch        ("relative",                   "Output dispersion relative to subband minima. If not specified, "
                                                         "the dispersion is given relative to the band edge.");
    opt.add_string_option ("mass-file",      "m_perp.r", "Filename from which to read in-plane mass profile.");
    opt.add_string_option ("alpha-file",     "alpha.r",  "Filename from which to read in-plane non-parabolicity profile.");
    opt.add_string_option ("potential-file", "v.r",      "Filename from which to read conduction band profile.");

    std::string summary = "Compute the dispersion relation for a set of subbands.";

    std::string details("Input files:\n"
                        "  'E*.r'   \tEnergy of each state:\n"
                        "           \tCOLUMN 1: state index.\n"
                        "           \tCOLUMN 2: energy [meV].\n"
                        "  'wf_*i.r'\tWave function amplitude at each position\n"
                        "           \tCOLUMN 1: position [m].\n"
                        "           \tCOLUMN 2: wave function amplitude [m^{-1/2}].\n"
                        "\n"
                        "\tIn each case, the '*' is replaced by the particle ID and the 'i' is replaced by the number of the state.\n"
                        "\n"
                        "Examples:\n"
                        "   Compute the ground state in a 150-angstrom well with effective mass = 0.1 m0:\n\n"
                        "   efiw --width 150 --mass 0.1\n"
                        "\n"
                        "   Compute the first three heavy-hole states in a 200-angstrom well, using effective mass = 0.62 m0:\n\n"
                        "   efiw --width 200 --mass 0.62 --particle h");

    opt.add_prog_specific_options_and_parse(argc, argv, summary, details);

    return opt;
};

int main (int argc, char* argv[])
{
    WfOptions opt = configure_options(argc, argv);

    const size_t nk   = opt.get_size_option("nk");
    const double nkbt = opt.get_numeric_option("nkbt");
    const double Te   = opt.get_numeric_option("Te");

    std::vector<Leeds::Subband> subbands;

    if(!opt.get_switch("nonparabolic"))
    {
        subbands = Leeds::Subband::read_from_file(opt.get_energy_input_path(),
                                                  opt.get_wf_input_prefix(),
                                                  opt.get_wf_input_ext(),
                                                  opt.get_string_option("mass-file"));
    }
    else
    {
        subbands = Leeds::Subband::read_from_file(opt.get_energy_input_path(),
                                                  opt.get_wf_input_prefix(),
                                                  opt.get_wf_input_ext(),
                                                  opt.get_string_option("mass-file"),
                                                  opt.get_string_option("alpha-file"),
                                                  opt.get_string_option("potential-file"));
    }

    // Loop over subbands
    unsigned int ist = 1;
    for(std::vector<Leeds::Subband>::iterator subband = subbands.begin(); subband < subbands.end(); ++subband)
    {
        std::valarray<double> k(nk);
        std::valarray<double> Ek(nk);

        // Calculate maximum wavevector
        double k_max = subband->k(nkbt*kB*Te);

        // Calculate wavevector spacing
        double dk = k_max/nk;

        // Loop over wavevectors and find corresponding energies
        for(unsigned int ik=0; ik<nk; ik++)
        {
            // Calculate wavevector
            k[ik] = ik*dk;
            // Calculate energy
            Ek[ik] = subband->Ek(k[ik]);
        }

        // If absolute energies are required then offset energies by the energy of the subband
        // minima
        if(!opt.get_switch("relative"))
            Ek += subband->get_E();

        Ek *= 1/(1e-3*e);

        // If verbose option selected output some information about subband
        if(opt.get_verbose())
        {
            std::cout << "Subband " << ist << " at " << subband->get_E()/(1e-3*e) << "eV." << std::endl;
            std::cout << "D.o.s effective mass: " << subband->get_md_0() << std::endl;
            std::cout << "D.o.s nonparabolicity parameter: " << subband->get_alphad() << std::endl;
            std::cout << "Wavevector at " << nkbt << "*kB*Te: " << k_max << std::endl;
            std::cout << std::endl;
        }

        // Construct filename and output
        std::stringstream filename;
        filename << opt.get_string_option("disp-prefix") << ist << opt.get_string_option("disp-ext");
        Leeds::write_table_xy(filename.str().c_str(), k, Ek);

        // Increment state counter
        ist++;
    }

    return EXIT_SUCCESS;
}
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :
