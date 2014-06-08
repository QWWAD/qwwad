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
    opt.add_numeric_option("mass",           0.067,      "In-plane effective mass (relative to free electron).");
    opt.add_numeric_option("alpha",          0.0,        "In-plane non-parabolicity parameter [1/eV].");
    opt.add_numeric_option("vcb",            0.0,        "Conduction band edge [eV].");

    std::string summary = "Compute the dispersion relation for a set of subbands.";

    std::string details("Default input files:\n"
                        "  'E*.r'    \tEnergy of each state:\n"
                        "            \tCOLUMN 1: state index.\n"
                        "            \tCOLUMN 2: energy [meV].\n"
                        "  'wf_*i.r' \tWave function amplitude at each position\n"
                        "            \tCOLUMN 1: position [m].\n"
                        "            \tCOLUMN 2: wave function amplitude [m^{-1/2}].\n"
                        "\n"
                        "Default output files:\n"
                        "  'dr_*i.r' \tDispersion relation for subband i\n"
                        "            \tCOLUMN 1: In-plane wave-vector [1/m]\n"
                        "            \tCOLUMN 2: Energy [meV]\n"
                        "\n"
                        "\tIn each case, the '*' is replaced by the particle ID and the 'i' is replaced by the number of the state.\n"
                        "\n"
                        "Examples:\n"
                        "   Compute the dispersion relation up to 10 kT with an electron temperature of 10 K:\n\n"
                        "   dispersion_relation --nkbt 10 --Te 10\n"
                        "\n"
                        "   Compute the non-parabolic dispersion relation using 1000 data points:\n\n"
                        "   dispersion_relation --nk 1000 --nonparabolic");

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
                                                  opt.get_numeric_option("mass") * me);
    }
    else
    {
        subbands = Leeds::Subband::read_from_file(opt.get_energy_input_path(),
                                                  opt.get_wf_input_prefix(),
                                                  opt.get_wf_input_ext(),
                                                  opt.get_numeric_option("mass") * me,
                                                  opt.get_numeric_option("alpha") / e,
                                                  opt.get_numeric_option("vcb") * e);
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

        Ek /= (1e-3*e); // Rescale to meV

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
