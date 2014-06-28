/** \file   chargedensity.cpp
 *  \brief  Calculates charge density profile for a 2D heterostructure
 *  \author Alex Valavanis <a.valavanis@leeds.ac.uk>
 *  \date   2012-02-08
 */

#if HAVE_CONFIG_H
# include "config.h"
#endif

#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <fstream>
#include <gsl/gsl_math.h>

#include "qclsim-fileio.h"
#include "qclsim-constants.h"
#include "qclsim-linalg.h"
#include "wf_options.h"

using namespace Leeds;
using namespace Leeds::constants;

class ChargeDensityOptions : public WfOptions
{
    public:
        ChargeDensityOptions(int argc, char* argv[]);
        void print() const {};

        size_t get_nper() const {
            return vm["nper"].as<size_t>();}

        std::string get_population_filename() const {
            return vm["population-file"].as<std::string>();}

        std::string get_doping_filename() const {
            return vm["doping-file"].as<std::string>();}

        std::string get_alloy_1per_filename() const {
            return vm["alloy-1per-file"].as<std::string>();}
        
        std::string get_degeneracy_filename() const {
            return vm["degeneracy-file"].as<std::string>();}

        std::string get_charge_filename() const {
            return vm["charge-file"].as<std::string>();}

        std::string get_edensity_filename() const {
            return vm["edensity-file"].as<std::string>();}
};

/** Holds standard-input data */
class ChargeDensityData
{
    public:
        ChargeDensityData(const ChargeDensityOptions& opt);
        std::vector<Leeds::State> states;  ///< Wf and energy data
        std::valarray<double> z;   ///< Spatial profile
        std::valarray<double> z_1per;   ///< Spatial profile
        std::valarray<double> x;   ///< alloy profile
        std::valarray<double> pop; ///< Subband population [m^{-2}]
        std::valarray<double> n3D; ///< Spatial profile of doping in structure [Cm^{-3}]
        std::valarray<unsigned int> nval; ///< Degeneracy of subbands
        size_t get_nz_1per() const {return nz_1per;};
    private:
        size_t nz_1per; ///< Number of spatial points in a single period
};

ChargeDensityOptions::ChargeDensityOptions(int argc, char* argv[])
{
    program_specific_options->add_options()
        ("nper,p", po::value<size_t>()->default_value(1),
         "Number of periods the wavefunctions cross")

        ("doping-file", 
         po::value<std::string>()->default_value("doping-profile.dat"),
         "Filename from which to read volume doping profile [m^{-3}]")

        ("alloy-1per-file", 
         po::value<std::string>()->default_value("alloy-profile_1per.dat"),
         "Filename from which to read alloy profile")
        
        ("population-file",
         po::value<std::string>()->default_value("populations.dat"),
         "Filename from which to read subband populations [m^{-2}]")

        ("degeneracy-file",
         po::value<std::string>()->default_value("degeneracy.dat"),
         "Filename from which to read subband degeneracies")

        ("charge-file",
         po::value<std::string>()->default_value("charge-density.dat"),
         "Filename to which charge density profile will be written")

        ("edensity-file",
         po::value<std::string>()->default_value("edensity.dat"),
         "Filename to which electron density profile will be written")
        ;

    std::string doc = 
        "Find charge density in a 2D heterostructure.  "
        "Input data should be provided for ONE period of the doping "
        "profile, and for all wavefunctions localised in that period. "
        "If the wavefunctions have tails that extend over multiple periods "
        "of a heterostructure, then use the --nper flag to specify how many. "
        "A table of charge density and electron density [m^-3] at each position [m] "
        "is given as the output.";

    add_prog_specific_options_and_parse(argc, argv, doc);

    if(vm["nper"].as<size_t>() < 1)
        throw std::domain_error("Number of periods must be one or more.");

    if(get_verbose()) print();
}

/** Checks that number of valleys is valid */
static void check_nval(unsigned int* pnval)
{
    if(*pnval != 1 and
       *pnval != 2 and
       *pnval != 3 and
       *pnval != 4 and
       *pnval != 6)
    {
        std::ostringstream oss;
        oss << "Invalid number of equivalent valleys (" << *pnval << ") detected."
            " It must be one of {1,2,3,4,6}";
        throw std::domain_error(oss.str());
    }
}

ChargeDensityData::ChargeDensityData(const ChargeDensityOptions& opt) :
    states(State::read_from_file(opt.get_energy_input_path(),
                                 opt.get_wf_input_prefix(),
                                 opt.get_wf_input_ext())),
    z(states[0].psi_array().size()),
    z_1per(states[0].psi_array().size()/opt.get_nper()),
    x(states[0].psi_array().size()/opt.get_nper()),
    pop(states.size()),
    n3D(states[0].psi_array().size()),
    nval(states.size()),
    nz_1per(z.size()/opt.get_nper())
{
    std::valarray<double> V(states[0].psi_array().size());
    read_table_xy(opt.get_potential_input_path().c_str(), z, V);

    std::valarray<double> z_n3D;

    // Read data for each subband
    read_table_x(opt.get_population_filename().c_str(), pop);
    read_table_x(opt.get_degeneracy_filename().c_str(), nval);
    read_table_xy(opt.get_doping_filename().c_str(), z_n3D, n3D);
    read_table_xy(opt.get_alloy_1per_filename().c_str(), z_1per, x);

    // Check that all input files have same size
    const size_t nst = pop.size();

    if(nval.size() != nst)
    {
        std::ostringstream oss;
        oss << "Different lengths of data were read from " << opt.get_population_filename()
            << " (" << nst << " lines) and " << opt.get_degeneracy_filename()
            << " (" << nval.size() << " lines).";
        throw std::length_error(oss.str());
    }

    // Loop through subbands and check that data is sensible
    for(unsigned int ist=0; ist < nst; ist++)
    {
        check_positive(&pop[ist]);
        check_nval(&nval[ist]);
    }

    // Check that doping profile looks sensible
    for(unsigned int iz=0; iz<nz_1per; iz++)
        check_not_negative(&n3D[iz]);
}

/**
 * \brief     main function for program
 *
 * \param[in] argc The number of command-line arguments
 * \param[in] argv Array of command-line arguments
 * 
 * \returns   Exit code of 0 signifies successful completion
 */
int main(int argc, char* argv[])
{
    const ChargeDensityOptions opt(argc, argv);
    const ChargeDensityData data(opt);
    const size_t nper    = opt.get_nper();     // Number of periods over which wavefunction spreads
    const size_t nz_1per = data.get_nz_1per(); // Number of points in a single period
    const size_t nst     = data.states.size(); // Number of states localised in one period

    // Find electron density at each point in a single period
    // by summing the "tails" of wavefunctions in each period
    std::valarray<double> edensity_1per(0.0, nz_1per);

    for(unsigned int iper = 0; iper < nper; iper++)
    {
        for(unsigned int ist = 0; ist < nst; ist++)
        {
            std::valarray<double> PD     = data.states[ist].psi_squared();
            std::valarray<double> PD_per = PD[std::slice(iper*nz_1per, nz_1per, 1)];
            edensity_1per += e * data.pop[ist] * data.nval[ist] * PD_per;
        }
    }

    // Charge density is obtained by subtracting electron density from doping density
    const std::valarray<double> rho_1per = e * data.n3D[std::slice(0, nz_1per, 1)] - edensity_1per;

    // Output position, charge density and electron density for a single period [Cm^-3]
    write_table_xy(opt.get_charge_filename().c_str(), data.z_1per, rho_1per);
    write_table_xy(opt.get_edensity_filename().c_str(), data.z_1per, edensity_1per);

    return EXIT_SUCCESS;
}

// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :
