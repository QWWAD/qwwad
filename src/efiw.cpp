/**
 * \file   efiw.cpp Calculates eigenstates of an infinite square well
 *
 * \author Paul Harrison  <p.harrison@shu.ac.uk>
 * \author Alex Valavanis <a.valavanis@leeds.ac.uk>
 *
 * \details The well width is passed via the command line.
 */

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <gsl/gsl_math.h>
#include "qwwad-options.h"
#include "struct.h"
#include "maths.h"
#include "const.h"

/**
 * Handler for command-line options
 */
class EFIWOptions : public Options
{
    public:
        EFIWOptions(int argc, char* argv[])
        {
            try
            {
                program_specific_options->add_options()
                    ("width,L", po::value<double>()->default_value(100), 
                     "Width of quantum well [angstrom].")

                    ("mass,m", po::value<double>()->default_value(0.067),
                     "Effective mass (relative to that of a free electron)")

                    ("npoint,N", po::value<size_t>()->default_value(100),
                     "Number of spatial points for output file")

                    ("particle,p", po::value<char>()->default_value('e'),
                     "Particle to be used: 'e', 'h' or 'l'")

                    ("states,s", po::value<size_t>()->default_value(1),
                     "Number of states to find")
                    ;

                std::string doc("Find the eigenstates of an infinite quantum well. "
                                "The energies are written to the file \"E*.r\", and the "
                                "wavefunctions are written to \"wf_*i.r\" where the '*' "
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

        /**
         * \returns the width of the quantum well [m]
         */
        double get_width() const {return vm["width"].as<double>()*1e-10;}

        /**
         * \returns the effective mass [kg]
         */
        double get_mass() const {return vm["mass"].as<double>()*m0;}

        /**
         * \returns the particle ID
         */
        char get_particle() const {return vm["particle"].as<char>();}

        /**
         * \returns the number of spatial points
         */
        size_t get_n_points() const {return vm["npoint"].as<size_t>();}

        /**
         * \returns the number of spatial points
         */
        size_t get_n_states() const {return vm["states"].as<size_t>();}
};

int main(int argc, char *argv[])
{
    EFIWOptions opt(argc, argv);

    const double L = opt.get_width();    // well width [m]
    const char   p = opt.get_particle(); // particle (e, h or l)
    const double m = opt.get_mass();     // effective mass [kg]
    const size_t N = opt.get_n_points(); // number of spatial steps
    const size_t s = opt.get_n_states(); // number of states
    char	filename[9];	/* wavefunction filename	*/

    sprintf(filename,"E%c.r",p);
    FILE *FE=fopen(filename,"w");

    // Loop over all required states
    for(unsigned int is=1; is<=s; is++)
    {
        // Energy of state [J] (QWWAD3, 2.13)
        const double E=gsl_pow_2(pi*hbar*is/L)/(2*m);
        fprintf(FE,"%i %24.17le\n",is,E/(1e-3*e_0));

        sprintf(filename,"wf_%c%i.r",p,is);
        FILE *Fwf=fopen(filename,"w");

        // Loop over spatial locations and find wavefunction
        // amplitude at each point (QWWAD3, 2.15)
        for(unsigned int i=0;i<N;i++)
        {
            const double z=(float)i*L/(float)(N-1); // Position [m]
            const double psi=sqrt(2/L)*sin(is*pi*z/L); // Wavefunction [m^{-0.5}]
            fprintf(Fwf,"%20.17le %20.17le\n",z,psi);
        }

        fclose(Fwf);
    }

    fclose(FE);

    return EXIT_SUCCESS;
}
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :
