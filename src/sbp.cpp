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

void     calc_dist(double Emin, double Ef, double m, double T, int nE, int s, const double alpha, const double V);

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
                program_specific_options->add_options()
                    ("fd,f", po::bool_switch()->default_value(false),
                     "Output Fermi-Dirac distribution.")

                    ("mass,m", po::value<double>()->default_value(0.067),
                     "Effective mass at the band edge (relative to that of a free electron)")

                    ("band-edge", po::value<double>()->default_value(0),
                     "Band edge potential [eV]")

                    ("alpha", po::value<double>()->default_value(0),
                     "Nonparabolicity parameter [eV^{-1}]")

                    ("nenergy,n", po::value<size_t>()->default_value(1000),
                     "Number of energy samples for output file")

                    ("particle,p", po::value<char>()->default_value('e'),
                     "Particle to be used: 'e', 'h' or 'l'")

                    ("temperature,T", po::value<double>()->default_value(300),
                     "Temperature of carrier distribution [K]")
                    ;

                std::string doc("Find the Fermi-Dirac distribution function for an individual "
                        "subband given its population and temperature."
                        "The subband minimum is read from the file \"E*.r\", "
                        "and the distribution is written to \"FDi.r\" where the '*' "
                        "is replaced by the particle ID in each case and the "
                        "'i' is replaced by the number of the state");

                add_prog_specific_options_and_parse(argc, argv, doc);	
            }
            catch(std::exception &e)
            {
                std::cerr << e.what() << std::endl;
                exit(EXIT_FAILURE);
            }
        }

        /// \returns the effective mass [kg]
        bool get_fd_flag() const {return vm["fd"].as<bool>();}

        /// \returns the effective mass [kg]
        double get_mass() const {return vm["mass"].as<double>()*me;}

        /// \returns the band edge [J]
        double get_band_edge() const {return vm["band-edge"].as<double>()*e;}

        /// \returns the nonparabolicity parameter [J^{-1}]
        double get_alpha() const {return vm["alpha"].as<double>()/e;}

        /// \returns the particle ID
        char get_particle() const {return vm["particle"].as<char>();}

        /// \returns the number of energy samples
        size_t get_n_energy() const {return vm["nenergy"].as<size_t>();}

        /// \returns the temperature of the carrier distributeion [K]
        double get_temperature() const {return vm["temperature"].as<double>();}
};

int main(int argc,char *argv[])
{
    SBPOptions opt(argc, argv);

    FILE	*FEf;			/* file pointer to Fermi Energy file	*/

    const bool   FD_flag = opt.get_fd_flag();
    const double m       = opt.get_mass();
    const double alpha   = opt.get_alpha();
    const double V       = opt.get_band_edge();
    const size_t nE      = opt.get_n_energy();
    const char   p       = opt.get_particle();
    const double T       = opt.get_temperature();

    std::valarray<double> E = read_E(p); // reads subband energy file
    const size_t n = E.size();
    std::valarray<double> N=read_populations(n); // reads subband populations file

    if((FEf=fopen("Ef.r","w"))==0)
    {
        fprintf(stderr,"Error: Cannot open output file 'Ef.r'!\n");
        exit(EXIT_FAILURE);
    }

    for(unsigned int s=0; s<n; s++) // s=0 => ground state
    {
        const double Ef=find_fermi(E[s],m,N[s],T,alpha,V);

        fprintf(FEf,"%i %20.17le\n",s+1,Ef/(1e-3*e));

        if(FD_flag) calc_dist(E[s],Ef,m,T,nE,s, alpha,V);
    }

    fclose(FEf);

    return EXIT_SUCCESS;
}

/**
 * \brief calculates the probability of occupation of the subband energies
 *
 * \param[in] Emin subband minima
 * \param[in] Ef   Fermi energy
 * \param[in] m    effective mass
 * \param[in] T    temperature
 * \param[in] nE   number of energies to output FD
 * \param[in] s    number of subband
 */
void calc_dist(double Emin, double Ef, double m, double T, int nE, int s, const double alpha, const double V)
{
    char   filename[9]; // output filename for FD distribs
    sprintf(filename,"FD%i.r",s+1);

    double vmax=Vmax(); // Maximum potential
    double Emax=Ef+10*kB*T; // Cut-off energy for plot [J]

    if(Emax<Emin) Emax=Emin+10*kB*T;
    if(Emax>vmax) Emax=vmax;

    std::valarray<double> E(nE); // Array of energies for plot
    std::valarray<double> f(nE); // Occupation probabilities

    const double dE=(Emax-Emin)/(nE-1); // Energy increment for integration
    for(int i=0; i<nE; i++)
    {
        E[i] = Emin + i*dE;
        f[i] = f_FD(Ef, E[i], T);
    }

    E/=(1e-3*e); // Convert to meV for output

    write_table_xy(filename, E, f);
    printf("Ne=%20.17le\n",find_pop(Emin, Ef, m, T,alpha,V)/1e+14);
}
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :
