/**
 * \file   srelo.cpp
 * \brief  Scattering Rate Electron-LO phonon
 * \author Paul Harrison  <p.harrison@shu.ac.uk>
 * \author Alex Valavanis <a.valavanis@leeds.ac.uk>
 */

#include <complex>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <iostream>
#include "qwwad/constants.h"
#include "qwwad/scattering-calculator-LO.h"
#include "qwwad/file-io.h"
#include "qwwad/subband.h"
#include "qwwad/options.h"

using namespace QWWAD;
using namespace constants;

void ff_output(const arma::vec &Kz,
               const arma::vec &Gifsqr,
               unsigned int     i,
               unsigned int     f);

static Options configure_options(int argc, char* argv[])
{
    Options opt;

    std::string doc("Find the polar LO-phonon scattering rate.");

    opt.add_option<bool>  ("outputff,a",            "Output form-factors to file.");
    opt.add_option<bool>  ("noblocking,b",          "Disable final-state blocking.");
    opt.add_option<bool>  ("noscreening,S",         "Disable screening.");
    opt.add_option<double>("latticeconst,A",  5.65, "Lattice constant in growth direction [angstrom]");
    opt.add_option<double>("ELO,E",          36.0,  "Energy of LO phonon [meV]");
    opt.add_option<double>("epss,e",         13.18, "Static dielectric constant");
    opt.add_option<double>("epsinf,f",       10.89, "High-frequency dielectric constant");
    opt.add_option<double>("mass,m",         0.067, "Band-edge effective mass (relative to free electron)");
    opt.add_option<char>  ("particle,p",       'e', "ID of particle to be used: 'e', 'h' or 'l', for "
                                                    "electrons, heavy holes or light holes respectively.");
    opt.add_option<double>("Te",               300, "Carrier temperature [K].");
    opt.add_option<double>("Tl",               300, "Lattice temperature [K].");
    opt.add_option<size_t>("nki",              101, "Number of initial wave-vector samples.");
    opt.add_option<size_t>("nKz",              101, "Number of phonon wave-vector samples.");

    opt.add_prog_specific_options_and_parse(argc, argv, doc);

    return opt;
}

int main(int argc,char *argv[])
{
    const auto opt = configure_options(argc, argv);

    const auto ff_flag     =  opt.get_option<bool>  ("outputff");             // True if formfactors are wanted
    const auto A0          =  opt.get_option<double>("latticeconst") * 1e-10; // Lattice constant [m]
    const auto Ephonon     =  opt.get_option<double>("ELO") * e/1000;         // Phonon energy [J]
    const auto epsilon_s   =  opt.get_option<double>("epss")   * eps0;        // Static permittivity [F/m]
    const auto epsilon_inf =  opt.get_option<double>("epsinf") * eps0;        // High-frequency permittivity [F/m]
    const auto m           =  opt.get_option<double>("mass")*me;              // Band-edge effective mass [kg]
    const auto p           =  opt.get_option<char>  ("particle");             // Particle ID
    const auto Te          =  opt.get_option<double>("Te");                   // Carrier temperature [K]
    const auto Tl          =  opt.get_option<double>("Tl");                   // Lattice temperature [K]
    const auto b_flag      = !opt.get_option<bool>  ("noblocking");           // Include final-state blocking by default
    const auto S_flag      = !opt.get_option<bool>  ("noscreening");          // Include screening by default
    const auto nki         =  opt.get_option<size_t>("nki");                  // number of ki calculations
    const auto nKz         =  opt.get_option<size_t>("nKz");                  // number of Kz calculations

    std::ostringstream E_filename; // Energy filename string
    E_filename << "E" << p << ".r";
    std::ostringstream wf_prefix;  // Wavefunction filename prefix
    wf_prefix << "wf_" << p;

    // Read data for all subbands from file
    auto subbands = Subband::read_from_file(E_filename.str(),
                                            wf_prefix.str(),
                                            ".r",
                                            m);

    // Read and set carrier distributions within each subband
    arma::vec  Ef;      // Fermi energies [J]
    arma::uvec indices; // Subband indices (garbage)
    read_table("Ef.r", indices, Ef);
    Ef *= e/1000.0; // Rescale to J

    for(unsigned int isb = 0; isb < subbands.size(); ++isb)
        subbands[isb].set_distribution_from_Ef_Te(Ef[isb], Te);

    // Initialise scattering calculators and set parameters
    ScatteringCalculatorLO em_calculator(subbands, A0, Ephonon, epsilon_s, epsilon_inf, m, Te, Tl, true);
    ScatteringCalculatorLO ab_calculator(subbands, A0, Ephonon, epsilon_s, epsilon_inf, m, Te, Tl, false);
    em_calculator.enable_screening(S_flag);
    ab_calculator.enable_screening(S_flag);
    em_calculator.enable_blocking(b_flag);
    ab_calculator.enable_blocking(b_flag);
    em_calculator.set_phonon_samples(nKz);
    ab_calculator.set_phonon_samples(nKz);
    em_calculator.set_ki_samples(nki);
    ab_calculator.set_ki_samples(nki);

    // Read list of wanted transitions
    arma::uvec i_indices;
    arma::uvec f_indices;

    read_table("rrp.r", i_indices, f_indices);
    const size_t ntx = i_indices.size();

    arma::vec Wabar(ntx);
    arma::vec Webar(ntx);

    // Loop over all desired transitions
    for(unsigned int itx = 0; itx < i_indices.size(); ++itx)
    {
        // Get subband indices.  Note that the -1 is needed because the
        // input file indexes subbands from 1 upward
        unsigned int i = i_indices[itx] - 1;
        unsigned int f = f_indices[itx] - 1;

        // Output form-factors if desired
        if(ff_flag)
        {
            const auto Kz     = em_calculator.get_Kz_table();
            const auto Gifsqr = em_calculator.get_ff_table(i,f);
            ff_output(Kz, Gifsqr, i,f);
        }

        const auto tx_em  = em_calculator.get_transition(i, f);
        const auto tx_ab  = ab_calculator.get_transition(i, f);
        const auto Weif   = tx_em.get_rate_table(); // Emission scattering rate at this wave-vector [1/s]
        const auto Waif   = tx_ab.get_rate_table(); // Absorption scattering rate at this wave-vector [1/s]
        auto Ei_em  = tx_em.get_Ei_total_table();  // Initial TOTAL energies [J]
        auto Ei_ab  = tx_ab.get_Ei_total_table();  // Initial TOTAL energies [J]
        Ei_em *= 1000.0/e; // Rescale to meV
        Ei_ab *= 1000.0/e; // Rescale to meV

        // output scattering rates versus TOTAL carrier energy
        char	filename_em[9];
        sprintf(filename_em, "LOe%i%i.r",i,f);	/* emission	*/
        char	filename_ab[9];
        sprintf(filename_ab,"LOa%i%i.r",i,f);	/* absorption	*/
        write_table(filename_em, Ei_em, Weif);
        write_table(filename_ab, Ei_ab, Waif);

        // Average rates over entire subband
        Wabar[itx] = tx_ab.get_average_rate();
        Webar[itx] = tx_em.get_average_rate();
    } /* end while over states */

    write_table("LOa-if.r", i_indices, f_indices, Wabar);
    write_table("LOe-if.r", i_indices, f_indices, Webar);

    return EXIT_SUCCESS;
}

/**
 * \brief outputs the formfactors into files
 */
void ff_output(const arma::vec &Kz,
               const arma::vec &Gifsqr,
               unsigned int     i,
               unsigned int     f)
{
 char	filename[9];	/* output filename				*/
 sprintf(filename,"G%i%i.r",i,f);	
 write_table(filename, Kz, Gifsqr);
}
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :
