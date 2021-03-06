/** 
 * \file   qwwad_charge_density.cpp
 * \brief  Calculates charge density profile for a 2D heterostructure
 * \author Alex Valavanis <a.valavanis@leeds.ac.uk>
 */

#if HAVE_CONFIG_H
# include "config.h"
#endif

#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <fstream>
#include <gsl/gsl_math.h>

#include "qwwad/data-checker.h"
#include "qwwad/file-io.h"
#include "qwwad/constants.h"
#include "qwwad/eigenstate.h"
#include "qwwad/wf_options.h"

using namespace QWWAD;
using namespace constants;

class ChargeDensityOptions : public WfOptions
{
public:
    ChargeDensityOptions(int argc, char** argv);
};

ChargeDensityOptions::ChargeDensityOptions(int argc, char** argv)
{
    std::string doc("Find charge density in a 2D heterostructure.");

    add_option<size_t>     ("nper,p",                    1, "Number of periods crossed by the wavefunction. This is "
                                                            "only relevant for periodic heterostructures.");
    add_option<std::string>("dopingfile",            "d.r", "File from which to read volume doping profile [m^{-3}]");
    add_option<std::string>("populationfile",        "N.r", "File from which to read subband populations [m^{-2}]");
    add_option<std::string>("degeneracyfile",               "File from which to read subband degeneracies. If not "
                                                            "specified, all states are taken to be non-degenerate");
    add_option<std::string>("chargefile",           "cd.r", "File to which charge density profile will be written");
    add_option<std::string>("carrierdensityfile", "dens.r", "File to which electron density profile will be written");
    add_option<bool>       ("ptype",                        "Dopants are to be treated as acceptors, and wavefunctions "
                                                            "treated as hole states");

    add_prog_specific_options_and_parse(argc, argv, doc);

    if(get_option<size_t>("nper") < 1) {
        throw std::domain_error("Number of periods must be one or more.");
    }
}

/** Holds standard-input data */
class ChargeDensityData
{
    private:
        size_t _nper;   ///< Number of periods crossed by wavefunction
        std::vector<Eigenstate> states_;  ///< Wf and energy data
        arma::vec  pop_;  ///< Subband population [m^{-2}]
        arma::uvec nval_; ///< Degeneracy of subbands

    public:
        ChargeDensityData(const ChargeDensityOptions& opt);

        [[nodiscard]] inline auto get_states() const -> decltype(states_) {return states_;}
        [[nodiscard]] inline auto get_pop()    const -> decltype(pop_)    {return pop_;}
        [[nodiscard]] inline auto get_nval()   const -> decltype(nval_)   {return nval_;}
};

ChargeDensityData::ChargeDensityData(const ChargeDensityOptions& opt) :
    _nper(opt.get_option<size_t>("nper")),
    states_(Eigenstate::read_from_file(opt.get_energy_filename(),
                                       opt.get_wf_prefix(),
                                       opt.get_wf_ext())),
    pop_(arma::zeros(states_.size())),
    nval_(arma::ones<arma::uvec>(states_.size()))
{
    // Read population of each subband
    const auto population_file = opt.get_option<std::string>("populationfile");
    read_table(population_file.c_str(), pop_);
    const size_t nst = pop_.size();

    // Check that populations are all positive
    DataChecker::check_positive(pop_);
    
    // Read state degeneracy if specified
    if(opt.get_argument_known("degeneracyfile"))
    {
        const auto degeneracy_file = opt.get_option<std::string>("degeneracyfile");
        read_table(degeneracy_file.c_str(), nval_);

        // Check that all input files have same size
        if(nval_.size() != nst)
        {
            std::ostringstream oss;
            oss << "Different lengths of data were read from " << population_file
                << " (" << nst << " lines) and " << degeneracy_file
                << " (" << nval_.size() << " lines).";
            throw std::length_error(oss.str());
        }
    }
}

auto main(int argc, char* argv[]) -> int
{
    const ChargeDensityOptions opt(argc, argv);

    const auto nper = opt.get_option<size_t>("nper"); // Number of periods over which wavefunction spreads

    const ChargeDensityData data(opt);
    const auto states = data.get_states();
    const auto pop    = data.get_pop();
    const auto nval   = data.get_nval();
 
    const auto nst  = states.size();          // Number of states localised in one period

    // Read doping profile (for entire multi-period structure)
    arma::vec z; // Spatial location [m]
    arma::vec d; // Spatial profile of doping in structure [Cm^{-3}]
    read_table(opt.get_option<std::string>("dopingfile").c_str(), z, d);

    // Compute the spatial points for a single period
    const size_t nz_1per = z.size() / nper; // Number of points in a single period

    // Get doping for one period by slicing it from total profile
    const arma::vec z_1per = z.subvec(0, nz_1per-1);
    const auto d_1per = d.subvec(0, nz_1per-1);

    // Find carrier density at each point in a single period
    // by summing the "tails" of wavefunctions in each period.
    // This implements the summation in [QWWAD4, 3.108]
    // [m^{-3}]
    arma::vec carrier_density_1per = arma::zeros<arma::vec>(nz_1per);

    for(unsigned int iper = 0; iper < nper; iper++)
    {
        for(unsigned int ist = 0; ist < nst; ist++)
        {
            // Find probability density function for carrier over the entire structure
            const auto PD     = states[ist].get_PD();

            // Grab the part of the PDF that lies in this period
            const arma::vec PD_per = PD.subvec(iper*nz_1per, (iper+1)*nz_1per-1);

            // Add this into the total carrier density profile
            carrier_density_1per += pop[ist] * static_cast<double>(nval[ist]) * PD_per;
        }
    }

    // Charge density is obtained by subtracting carrier density from doping density
    // [QWWAD4, 3.108]. Note q = -e by default (for electrons). [C m^{-3}]
    arma::vec rho_1per = e*(d - carrier_density_1per);

    // Invert charge profile if it's a p-type system
    if (opt.get_option<bool>("ptype")) {
        rho_1per *= -1;
    }

    // Output position, charge density [Cm^{-3}] and carrier density [m^{-3}] for a single period

    try {
        write_table(opt.get_option<std::string>("chargefile").c_str(), z_1per, rho_1per);
    } catch (std::runtime_error &e) {
        std::cerr << "Error writing file" << std::endl;
        std::cerr << e.what() << std::endl;
    }

    try {
        write_table(opt.get_option<std::string>("carrierdensityfile").c_str(), z_1per, carrier_density_1per);
    } catch (std::runtime_error &e) {
        std::cerr << "Error writing file" << std::endl;
        std::cerr << e.what() << std::endl;
    }

    return EXIT_SUCCESS;
}

// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :
