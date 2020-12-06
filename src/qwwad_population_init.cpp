/**
 * \file   qwwad_population_input.cpp
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
#include "qwwad/constants.h"
#include "qwwad/maths-helpers.h"
#include "qwwad/file-io.h"
#include "qwwad/fermi.h"
#include "qwwad/wf_options.h"

using namespace QWWAD;
using namespace constants;

/// Type of carrier energy distribution
enum DistributionType {
    DIST_GROUND, ///< All carriers are in ground state
    DIST_EVEN, ///< Carriers spread evenly through structure
    DIST_FERMI ///< Carriers follow a thermal distribution
};

class DensityinputOptions : public Options 
{
    private:
        DistributionType distType = DIST_EVEN; ///< The type of carrier distribution to use
    public:
        DensityinputOptions(int argc, char** argv);

        // Program options
        DistributionType get_dist_type() const {return distType;}
};

DensityinputOptions::DensityinputOptions(int argc, char** argv)
{
    std::string doc("Generate an approximate carrier energy distribution for a set of subbands");

    add_option<std::string>("dopingfile",      "d.r",  "Filename from which doping profile is read [m^{-3}]");
    add_option<std::string>("energyfile",      "Ee.r", "Filename from which subband energies are read [meV]");
    add_option<double>     ("mass",            0.067,  "In-plane effective mass (relative to free electron)");
    add_option<double>     ("Te",                100,  "Temperature of carrier distribution [K]");
    add_option<std::string>("populationfile",  "N.r",  "Set filename to which subband populations are written [m^{-2}]");
    add_option<size_t>     ("nval",                1,  "Split population between a number of equivalent valleys");
    add_option<std::string>("type",           "even",  "Type of carrier distribution across states. Permitted "
                                                       "options are: fermi, ground or even");

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
}

int main(int argc, char *argv[])
{
    DensityinputOptions opt(argc, argv);
    const auto nval = opt.get_option<size_t>("nval"); // Number of equivalent valleys
    
    // Integrate doping profile to find sheet doping [m^{-2}]

    arma::vec z;   ///< Spatial points [m]
    arma::vec d;   ///< Volume doping profile [m^{-3}]
    read_table(opt.get_option<std::string>("dopingfile").c_str(), z, d);

    const double dz  = z[1] - z[0];  // Spatial step [m]
    const double n2D = trapz(d,dz); // Sheet doping [m^{-2}]

    arma::uvec _inx; // State indices
    arma::vec E;          // Energies of subband minima [J]
    read_table(opt.get_option<std::string>("energyfile").c_str(), _inx, E);
    E *= e/1000; // Rescale to J

    const size_t nst = E.size();    // Number of subbands
    arma::vec pop(nst); // Population of each subband

    // Generate distribution depending on the user-specified type
    switch(opt.get_dist_type())
    {
        // Split sheet density evenly over subbands
        case DIST_EVEN:
            for(unsigned int i=0; i < nst; i++)
                pop[i] = n2D/(nst*nval);
            break;

        // Thermal distribution of subband populations
        case DIST_FERMI:
            {
                const auto _md = opt.get_option<double>("mass") * me; // Density-of-states mass [kg]
                const auto T   = opt.get_option<double>("temperature");

                // Fermi energy for entire system [J]
                double Ef = find_fermi_global(E, _md, n2D, T);

                if(opt.get_verbose())
                    std::cout << "Fermi energy = " << Ef << " J (" << Ef *1000/e << " meV)." << std::endl;

                for(unsigned int i=0; i<nst; i++)
                    pop[i] = find_pop(E[i], Ef, _md, T) / nval;
            }
            break;

        case DIST_GROUND:
            pop[0] = n2D/nval;

            // For all excited states, use a population of 1 electron/cm^2
            // i.e. a ridiculously low value.  This prevents the
            // physically-impossible zero population states from breaking
            // later calculations of Fermi energy etc
            for(unsigned int i=1;i<nst;i++)
                pop[i] = 1.0;
            break;
    }

    write_table(opt.get_option<std::string>("populationfile").c_str(), pop);

    return EXIT_SUCCESS;
}
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :
