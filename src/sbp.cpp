/**
 * \file  sbp.cpp
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
#include <cstdio>
#include <cstdlib>
#include <gsl/gsl_math.h>
#include "qclsim-constants.h"
#include "qclsim-fermi.h"
#include "qclsim-fileio.h"
#include "qwwad-fileio.h"
#include "qwwad-options.h"

using namespace Leeds;
using namespace constants;

static double calc_dist(const double       Emin,
                        const double       Ef,
                        const double       m,
                        const double       T,
                        const size_t       nE,
                        const unsigned int s,
                        const double       alpha,
                        const double       V);

/**
 * Handler for command-line options
 */
class SBPOptions : public Options
{
    public:
        SBPOptions(int argc, char* argv[])
        {
            try
            {
                std::string summary("Find the Fermi-Dirac distribution functions for a set of subbands.");

                add_switch        ("fd,f",                      "Output Fermi-Dirac distribution.");
                add_numeric_option("mass,m",             0.067, "Effective mass (relative to free electron)");
                add_numeric_option("vcb",                 0.00, "Band-edge potential [eV]");
                add_numeric_option("alpha",               0.00, "Non-parabolicity parameter [eV^{-1}]");
                add_char_option   ("particle,p",          'e',  "ID of particle to be used: 'e', 'h' or 'l', for "
                                                                "electrons, heavy holes or light holes respectively.");
                add_numeric_option("Te",                  300,  "Carrier temperature [K].");
                add_size_option   ("nenergy,n",           1000, "Number of energy samples to print out");
                add_numeric_option("global-population,N", 0.0,  "Use equilibrium population for the entire system "
                                                                "instead of reading subband "
                                                                "populations from file [x1e10 cm^{-2}]");

                add_prog_specific_options_and_parse(argc, argv, summary);
            }
            catch(std::exception &e)
            {
                std::cerr << e.what() << std::endl;
                exit(EXIT_FAILURE);
            }
        }

        /// \returns the global population [m^{-2}]
        double get_global_pop() const {return get_numeric_option("global-population") * 10000 *1e10;}

        /// \returns true if the system is in thermal equilibrium
        bool equilibrium() const {return (vm.count("global-population") > 0 and gsl_fcmp(get_global_pop(),0,1e-6));}
};

int main(int argc,char *argv[])
{
    SBPOptions opt(argc, argv);

    const bool   FD_flag = opt.get_switch("fd");
    const double m       = opt.get_numeric_option("mass") * me;
    const double alpha   = opt.get_numeric_option("alpha") / e;     // Non-parabolicity [1/J]
    const double V       = opt.get_numeric_option("vcb") * e;       // band_edge potential [J]
    const char   p       = opt.get_char_option("particle");         // particle ID (e, h or l)
    const double T       = opt.get_numeric_option("Te");
    const size_t nE      = opt.get_size_option("nenergy");

    std::valarray<double> E = read_E(p); // Reads subband energy file [J]
    const size_t n = E.size();

    std::valarray<double> Ef(n); // Fermi energies for each subband [J]
    std::valarray<double> N(n); // Population of each subband [m^{-3}]

    if(opt.equilibrium())
    {
        const double N_total = opt.get_global_pop();
        const double Ef_global = find_fermi_global(E, m, N_total, T, alpha, V);

        for(unsigned int i=0; i<n; ++i)
        {
            Ef[i] = Ef_global;
            N[i]  = find_pop(E[i], Ef[i], m, T, alpha, V);
        }

        N /= 1e14; // Rescale to 1e10 cm^{-2}
        Leeds::write_table_x("N-out.r", N, true, 17);
    }
    else
    {
        // reads subband populations file
        std::valarray<double> N(n);
        read_table("N.r", N);

        for(unsigned int i=0; i<n; ++i) // i=0 => ground state
            Ef[i]=find_fermi(E[i],m,N[i],T,alpha,V);
    }

    for(unsigned int i=0; i<n; ++i)
    {
        if(FD_flag)
        {
            const double N = calc_dist(E[i],Ef[i],m,T,nE,i, alpha,V);

            if(opt.get_verbose())
                printf("Ne=%20.17le\n", N/1e+14);
        }
    }
    Ef *= 1000.0/e; // Rescale to meV
    Leeds::write_table_x("Ef.r", Ef, true, 17);

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
 */
static double calc_dist(const double       Emin,
                        const double       Ef,
                        const double       m,
                        const double       T,
                        const size_t       nE,
                        const unsigned int s,
                        const double       alpha,
                        const double       V)
{
    char   filename[9]; // output filename for FD distribs
    sprintf(filename,"FD%i.r",s+1);

    double Emax=Ef+10*kB*T; // Cut-off energy for plot [J]

    if(Emax<Emin) Emax=Emin+10*kB*T;

    std::valarray<double> E(nE); // Array of energies for plot
    std::valarray<double> f(nE); // Occupation probabilities

    const double dE=(Emax-Emin)/(nE-1); // Energy increment for integration
    for(unsigned int i=0; i<nE; i++)
    {
        E[i] = Emin + i*dE;
        f[i] = f_FD(Ef, E[i], T);
    }

    E/=(1e-3*e); // Convert to meV for output

    write_table_xy(filename, E, f);
    return find_pop(Emin, Ef, m, T,alpha,V);
}
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :
