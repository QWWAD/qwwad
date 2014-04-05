/**
 * \file    dispersion_relation.cpp
 * \brief   Prints out the dispersion relations for each subband
 * \author  Jonathan Cooper <el06jdc@leeds.ac.uk>
 * \date    2013-02-14
 */

#include <iostream>
#include <valarray>

#include "wf_options.h"
#include "qclsim-subband.h"

using namespace Leeds;

#include "qclsim-constants.h"
using namespace constants;

class DispersionOptions : public WfOptions
{
    public:
        void print() const{};
        inline bool use_relative_energy_scale() const {return vm["relative"].as<bool>();}
        inline double get_nkbt()                const {return vm["nkbt"].as<double>();}
        inline double get_Te()                  const {return vm["Te"].as<double>();}
        inline unsigned int get_nk()            const {return vm["nk"].as<unsigned int>();}
        inline std::string get_disp_prefix()    const {return vm["disp-prefix"].as<std::string>();}
        inline std::string get_disp_ext()       const {return vm["disp-ext"].as<std::string>();}

        // TODO All this should be in some generic program options class and not defined in every
        // program that uses them!
	std::string get_m_d_filename() const {return vm["massd-file"].as<std::string>();}
        std::string get_alphad_filename() const {return vm["alphad-file"].as<std::string>();}
        std::string get_potential_filename() const {return vm["potential-file"].as<std::string>();}
        std::string get_populations_filename() const {return vm["population-file"].as<std::string>();}		
        std::string get_fermi_energy_filename() const {return vm["fermienergy-file"].as<std::string>();}		
        bool get_parabolic() const {return vm["parabolic"].as<bool>();}

        DispersionOptions(int argc, char* argv[])
        {
            program_specific_options->add_options()
                ("nkbt", po::value<double>()->default_value(5.0),
                 "Set the maximum energy to print out for the subband")

                ("Te", po::value<double>()->default_value(100),
                 "Set the electron temperature.")

                ("nk", po::value<unsigned int>()->default_value(100),
                 "Set the number of k-space points to print out")

                ("disp-prefix", po::value<std::string>()->default_value("dr_e"),
                 "Set filename prefix to which the dispersion curves for all states will be written")

                ("disp-ext", po::value<std::string>()->default_value(".dat"),
                 "Set filename ext to which the dispersion curves for all states will be written")

                ("parabolic", po::bool_switch()->default_value(false),
                 "Use parabolic dispersion relation.  If this flag is not set, we use nonparabolic dispersion.")

                ("relative", po::bool_switch(),
                 "Output dispersions relative to the subband minima. If this is not specified, "
                 "the dispersion curve is given relative to the band edge.")

                ("fermienergy-file",
                 po::value<std::string>()->default_value("Ef.dat"),
                 "Set filename from which fermi energies will be read from.")

                ("population-file",
                 po::value<std::string>()->default_value("populations.dat"),
                 "Set filename from which subband populations will be read from.")

                ("massd-file",
                 po::value<std::string>()->default_value("massd.dat"),
                 "Set filename from which to read the d.o.s. mass for one period of the strucutre")

                ("alphad-file",
                 po::value<std::string>()->default_value("alphad.dat"),
                 "Set filename from which to read the in-plane nonparabolicity parameter.")

                ("potential-file",
                 po::value<std::string>()->default_value("Vtotal.dat"),
                 "Set filename from which to read the conduction band profile.")
            ;

            std::string doc = "Prints the dispersion relations for all states.";

            add_prog_specific_options_and_parse(argc, argv, doc);
        }
};

int main (int argc, char* argv[])
{
    DispersionOptions opt(argc, argv);

    std::vector<Leeds::Subband> subbands;

    if(opt.get_parabolic())
    {
        subbands = Leeds::Subband::read_from_file(opt.get_energy_input_path(),
                                                  opt.get_wf_input_prefix(),
                                                  opt.get_wf_input_ext(),
                                                  opt.get_populations_filename(),
                                                  opt.get_fermi_energy_filename(),
                                                  opt.get_m_d_filename());
    }
    else
    {
        subbands = Leeds::Subband::read_from_file(opt.get_energy_input_path(),
                                                  opt.get_wf_input_prefix(),
                                                  opt.get_wf_input_ext(),
                                                  opt.get_populations_filename(),
                                                  opt.get_fermi_energy_filename(),
                                                  opt.get_m_d_filename(),
                                                  opt.get_alphad_filename(),
                                                  opt.get_potential_filename());
    }

    // Loop over subbands
    unsigned int ist = 1;
    for(std::vector<Leeds::Subband>::iterator subband = subbands.begin(); subband < subbands.end(); ++subband)
    {
        std::valarray<double> k(opt.get_nk());
        std::valarray<double> Ek(opt.get_nk());

        // Calculate maximum wavevector
        double k_max = subband->k(opt.get_nkbt()*kB*opt.get_Te());
        // Calculate wavevector spacing
        double dk = k_max/opt.get_nk();

        // Loop over wavevectors and find corresponding energies
        for(unsigned int ik=0; ik<opt.get_nk(); ik++)
        {
            // Calculate wavevector
            k[ik] = ik*dk;
            // Calculate energy
            Ek[ik] = subband->Ek(k[ik]);
        }

        // If absolute energies are required then offset energies by the energy of the subband
        // minima
        if(!opt.use_relative_energy_scale()){
            Ek += subband->get_E();
        }

        Ek *= 1/(1e-3*e);

        // If verbose option selected output some information about subband
        if(opt.get_verbose())
        {
            std::cout << "Subband " << ist << " at " << subband->get_E()/(1e-3*e) << "eV." << std::endl;
            std::cout << "D.o.s effective mass: " << subband->get_md_0() << std::endl;
            std::cout << "D.o.s nonparabolicity parameter: " << subband->get_alphad() << std::endl;
            std::cout << "Wavevector at " << opt.get_nkbt() << "*kB*Te: " << k_max << std::endl;
            std::cout << std::endl;
        }

        // Construct filename and output
        std::stringstream filename;
        filename << opt.get_disp_prefix() << ist << opt.get_disp_ext();
        Leeds::write_table_xy(filename.str().c_str(), k, Ek);

        // Increment state counter
        ist++;
    }

    return EXIT_SUCCESS;
}
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :
