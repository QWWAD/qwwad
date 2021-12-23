/**
 * \file  qwwad_fermi_distribution.cpp
 * \brief Calculate subband populations
 *
 * \details This program generates the Fermi-Dirac distribution function
 *          for an individual subband given its population and the lattice
 *          temperature.
 *
 *          Input files:
 *            Ep.r  Subband energies file, p=e,h,l
 *
 *          Output files:
 *            FDX.r F-D distribution for subband X
 *
 * \author Paul Harrison  <p.harrison@shu.ac.uk>
 * \author Alex Valavanis <a.valavanis@leeds.ac.uk>
 */

#include <iostream>
#include <sstream>
#include <cstdio>
#include <cstdlib>
#include <gsl/gsl_math.h>
#include "qwwad/constants.h"
#include "qwwad/fermi.h"
#include "qwwad/file-io.h"
#include "qwwad/options.h"

using namespace QWWAD;
using namespace constants;

// Conversion factor for 1e10 cm^{-2} to m^{-2}
constexpr double pop_unit_conv = 10000*1e10;

static auto calc_dist(double       Emin,
                      double       Ef,
                      double       m,
                      double       T,
                      size_t       nE,
                      unsigned int s,
                      double       alpha,
                      double       V) -> double;

/**
 * Handler for command-line options
 */
class SBPOptions : public Options
{
    public:
        SBPOptions(int argc, char** argv)
        {
            std::string summary("Find the Fermi-Dirac distribution functions for a set of subbands.");

            const double m_GaAs = 0.067;
            const double Te_def = 300;
            const double nE_def = 1000;

            add_option<bool>  ("fd,f",                        "Output Fermi-Dirac distribution.");
            add_option<double>("mass,m",              m_GaAs, "Effective mass (relative to free electron)");
            add_option<double>("vcb",                 0.00,   "Band-edge potential [eV]");
            add_option<double>("alpha",               0.00,   "Non-parabolicity parameter [eV^{-1}]");
            add_option<char>  ("particle,p",          'e',    "ID of particle to be used: 'e', 'h' or 'l', for "
                    "electrons, heavy holes or light holes respectively.");
            add_option<double>("Te",                  Te_def, "Carrier temperature [K].");
            add_option<size_t>("nenergy,n",           static_cast<size_t>(nE_def), "Number of energy samples to print out");
            add_option<double>("global-population,N", 0.0,    "Use equilibrium population for the entire system "
                    "instead of reading subband "
                    "populations from file [x1e10 cm^{-2}]");

            add_prog_specific_options_and_parse(argc, argv, summary);
        }

        /// \returns the global population [m^{-2}]
        [[nodiscard]] auto get_global_pop() const -> double
        {
            return get_option<double>("global-population") * pop_unit_conv;
        }

        /// \returns true if the system is in thermal equilibrium
        [[nodiscard]] auto equilibrium() const -> bool
        {
            const double epsilon = 1e-6; // Small number for comparison
            return get_argument_known("global-population") && (gsl_fcmp(get_global_pop(),0,epsilon) != 0);
        }
};

auto main(int argc,char *argv[]) -> int
{
    SBPOptions opt(argc, argv);

    const auto FD_flag = opt.get_option<bool>  ("fd");
    const auto m       = opt.get_option<double>("mass") * me;
    const auto alpha   = opt.get_option<double>("alpha") / e;     // Non-parabolicity [1/J]
    const auto V       = opt.get_option<double>("vcb") * e;       // band_edge potential [J]
    const auto p       = opt.get_option<char>  ("particle");        // particle ID (e, h or l)
    const auto T       = opt.get_option<double>("Te");
    const auto nE      = opt.get_option<size_t>("nenergy");

    // Read energies from file
    std::ostringstream Efile;
    Efile << "E" << p << ".r";
    arma::uvec idx;
    arma::vec  E;
    read_table(Efile.str().c_str(), idx, E);
    E*=e/1000; // Rescale to J

    const auto nst = E.size();

    arma::vec Ef(nst); // Fermi energies for each subband [J]
    arma::vec N(nst); // Population of each subband [m^{-2}]

    if(opt.equilibrium())
    {
        const auto N_total   = opt.get_global_pop();
        const auto Ef_global = find_fermi_global(E, m, N_total, T, alpha, V);

        for(unsigned int i=0; i<nst; ++i)
        {
            Ef[i] = Ef_global;
            N[i]  = find_pop(E[i], Ef[i], m, T, alpha, V);
        }

        N /= pop_unit_conv; // Rescale to 1e10 cm^{-2}
        write_table("N-out.r", N, true, 17);
    }
    else
    {
        // reads subband populations file
        read_table("N.r", idx, N);

        if(N.size() != nst) {
            std::cerr << "Populations file, N.r contains data for " << N.size() << " states but " << Efile.str() << " has " << nst << std::endl;
            exit(EXIT_FAILURE);
        }

        for(unsigned int i=0; i<nst; ++i) { // i=0 => ground state
            Ef[i] = find_fermi(E[i],m,N[i],T,alpha,V);
        }
    }

    if(FD_flag) {
        for(unsigned int i=0; i<nst; ++i) {
            const auto subband_pop = calc_dist(E[i],Ef[i],m,T,nE,i, alpha,V);

            if(opt.get_verbose()) {
                std::cout << "Ne = " << subband_pop/pop_unit_conv << std::endl;
            }
        }
    }

    Ef *= 1000.0/e; // Rescale to meV

    try {
        write_table("Ef.r", Ef, true, 17);
    } catch (std::runtime_error &e) {
        std::cerr << "Error writing file" << std::endl;
        std::cerr << e.what() << std::endl;
    }

    return EXIT_SUCCESS;
}

/**
 * \brief calculates the probability of occupation of the subband energies
 *
 * \param[in] Emin  subband minima [J]
 * \param[in] Ef    Fermi energy [J]
 * \param[in] m     effective mass at band edge [kg]
 * \param[in] T     temperature [K]
 * \param[in] nE    number of energies to output FD
 * \param[in] s     number of subband
 * \param[in] alpha nonparabolicity [1/J]
 * \param[in] V     band edge [J]
 *
 * \returns Total population of the subband [m^{-2}]
 */
static auto calc_dist(const double       Emin,
                        const double       Ef,
                        const double       m,
                        const double       T,
                        const size_t       nE,
                        const unsigned int s,
                        const double       alpha,
                        const double       V) -> double
{
    // output filename for FD distribs
    std::ostringstream filename_stream;
    filename_stream << "FD" << s+1 << ".r";

    auto Emax=Ef+10*kB*T; // Cut-off energy for plot [J]

    if(Emax<Emin) {
        Emax=Emin+10*kB*T;
    }

    arma::vec E(nE); // Array of energies for plot
    arma::vec f(nE); // Occupation probabilities

    const auto dE=(Emax-Emin)/static_cast<double>(nE-1); // Energy increment for integration
    for(unsigned int i=0; i<nE; i++) {
        E[i] = Emin + i*dE;
        f[i] = f_FD(Ef, E[i], T);
    }

    E/=(1e-3*e); // Convert to meV for output

    write_table(filename_stream.str(), E, f);
    return find_pop(Emin, Ef, m, T,alpha,V);
}
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :
