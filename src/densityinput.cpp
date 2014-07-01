/**
 * \file   densityinput.cpp
 * \brief  Generate approximate subband populations
 * \author Andrew Grier <el09a2g@leeds.ac.uk>
 * \author Alex Valavanis <a.valavanis@leeds.ac.uk>
 */

#if HAVE_CONFIG_H
# include "config.h"
#endif

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <string>
#include <fstream>
#include <iostream>
#include "qclsim-constants.h"
#include "qclsim-linalg.h"
#include "qclsim-fermi.h"
#include "wf_options.h"

using namespace Leeds;
using namespace Leeds::constants;

/// Type of carrier energy distribution
enum DistributionType {
    DIST_GROUND, ///< All carriers are in ground state
    DIST_EVEN, ///< Carriers spread evenly through structure
    DIST_FERMI ///< Carriers follow a thermal distribution
};

class DensityinputOptions : public WfOptions 
{
    private:
        DistributionType distType; ///< The type of carrier distribution to use
    public:
        DensityinputOptions(int argc, char* argv[]);
        void print() const {}

        // Program options
        inline double get_temp() const {return vm["temperature"].as<double>();}
        inline double get_Ef() const {return vm["Ef"].as<double>();}
        inline size_t get_nval() const {return vm["nval"].as<size_t>();}
        DistributionType get_dist_type() const {return distType;}

        bool use_fixed_Ef() const {return vm["fixed-Ef"].as<bool>();}

        // Output file
        std::string get_population_filename() const {return vm["population-file"].as<std::string>();}		
};

DensityinputOptions::DensityinputOptions(int argc, char* argv[]) :
    distType(DIST_EVEN)
{
    add_string_option ("doping-file", "d.r", "Filename from which doping profile is read [m^{-3}]");
    add_numeric_option("mass",        0.067, "In-plane effective mass (relative to free electron)");

    program_specific_options->add_options()
        // Output files
        ("population-file",
         po::value<std::string>()->default_value("N.r"),
         "Set filename to which subband populations are written [m^{-2}]")

        // Program options
        ("temperature,T", po::value<double>()->default_value(100),
         "Temperature of electron distribution [K]")

        ("nval",
         po::value<size_t>()->default_value(1),
         "Split population between a number of equivalent valleys")

        ("type", po::value<std::string>()->default_value("even"),
         "Type of carrier distribution across states")
        
         ("fixed-Ef,e", po::bool_switch()->default_value(false),
         "Define a fixed fermi energy (eV) relative to potential at top of stack")

        ("Ef,f", po::value<double>()->default_value(-0.5),
         "Fermi energy offset [eV]")
        ;

    std::string doc = "Generate an approximate carrier energy distribution (even, thermal or all in ground-state)";
    add_prog_specific_options_and_parse(argc, argv, doc);

    // Parse the distribution type parameter
    if(vm.count("type"))
    {
        const char* arg = vm["type"].as<std::string>().c_str();

        if(!strcmp(arg, "fermi") ||
                !strcmp(arg, "Fermi") ||
                !strcmp(arg, "thermal") ||
                !strcmp(arg, "Thermal"))
            distType = DIST_FERMI;
        else if(!strcmp(arg, "ground") ||
                !strcmp(arg, "Ground") ||
                !strcmp(arg, "cold") ||
                !strcmp(arg, "Cold"))
            distType = DIST_GROUND;
        else if(!strcmp(arg, "even") ||
                !strcmp(arg, "Even") ||
                !strcmp(arg, "equal") ||
                !strcmp(arg, "Equal"))
            distType = DIST_EVEN;
        else
        {
            std::cerr << "Unrecognised distribution type: " << arg << std::endl;
            exit(EXIT_FAILURE);
        }
    }

    if(get_verbose()) print();
}

int main(int argc, char *argv[]){
    DensityinputOptions opt(argc, argv);
    const double _md = opt.get_numeric_option("mass") * me; // Density-of-states mass [kg]

    std::vector<State> states = State::read_from_file(opt.get_energy_input_path(),
                                                      opt.get_wf_input_prefix(),
                                                      opt.get_wf_input_ext(),
                                                      1000.0/e,
                                                      true);

    std::valarray<double> V;   ///< Potential profile
    std::valarray<double> z;   ///< Spatial points [m]
    std::valarray<double> d;   ///< Volume doping profile [m^{-3}] 
    read_table_xy(opt.get_potential_input_path().c_str(), z, V);
    read_table_xy(opt.get_string_option("doping-file").c_str(), z, d);

    // Integrate doping profile to find sheet doping [m^{-2}]
    const double dz = z[1] - z[0];
    const double nz = z.size();
    const double n2D = trapz(d,dz);

    const size_t nst = states.size();
    std::valarray<double> pop(nst); // Population of each subband

    // Generate distribution depending on the user-specified type
    switch(opt.get_dist_type())
    {
        // Split sheet density evenly over subbands
        case DIST_EVEN:
            for(unsigned int i=0; i < nst; i++)
                pop[i] = n2D/(nst*opt.get_nval());
            break;

        // Thermal distribution of subband populations
        case DIST_FERMI:
            {
                // Fermi energy for entire system
                double Ef = 0.0;
                if(opt.use_fixed_Ef()==0)
                    Ef = find_fermi_global(states, _md, n2D, opt.get_temp());
                else
                    Ef = V[nz-1]+opt.get_Ef()*e;

                if(opt.get_verbose())
                    printf("Fermi energy = %g J\n",Ef);

                for(unsigned int i=0; i<nst; i++)
                {
                    double temp = find_pop(states[i].get_E(), Ef, _md, opt.get_temp()) / opt.get_nval();
                    if(temp>1)
                        pop[i] = temp;
                    else
                        pop[i] = 1.0;
                }
            }
            break;

        case DIST_GROUND:
            pop[0] = n2D/opt.get_nval();

            // For all excited states, use a population of 1 electron/cm^2
            // i.e. a ridiculously low value.  This prevents the
            // physically-impossible zero population states from breaking
            // later calculations of Fermi energy etc
            for(unsigned int i=1;i<nst;i++)
                pop[i] = 1.0;
            break;

        default:
            std::cerr << "Unknown distribution type" << std::endl;
            exit(EXIT_FAILURE);
    }

    write_table_x(opt.get_population_filename().c_str(), pop);

    return EXIT_SUCCESS;
}
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :
