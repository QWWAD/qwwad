/**
 * \file    efkpsl.cpp
 * \brief   Calculate the energy levels for a Kronig-Penney superlattice
 * \author  Paul Harrison  <p.harrison@shu.ac.uk>
 * \author  Alex Valavanis <a.valavanis@leeds.ac.uk>
 *
 * \details This program calculates the energy levels for a user 
 *          supplied wave vector along the growth (z-) axis of an 
 *          infinite Kronig-Penney superlattice.
 *          The code is based around:

 Different masses in well and barrier

 -hbar^2 d^2 Y + V Y = E Y,   0<z<a, V=0, a<z<b, V=V
 ------  ----                 
 2m*   dz^2                

 with the additional constraint of the boundary conditions

 Y and 1 dY , continuous
 - --
 m dz

 this represents the Hamiltonian P 1 P + V
 -
 m

 The system is solved by expressing the
 standard condition as a function f(x)=0 and
 implementing a Newton-Raphson iterative method
 where the independent variable x is the energy.
 The first estimate of the energy is found by 
 following the function along the x axis until it
 changes sign then using a midpoint rule.
 */

#include <iostream>
#include <cstdlib>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>
#include <valarray>
#include "qwwad-options.h"
#include "qclsim-constants.h"
#include "qclsim-fileio.h"

using namespace Leeds;
using namespace constants;

/// Parameters needed for computing superlattice eigenstates
struct kpsl_params
{
    double a;           ///< Well width [m]
    double b;           ///< Barrier width [m]
    double m_w;         ///< Well mass [kg]
    double m_b;         ///< Barrier mass [kg]
    double V;           ///< Barrier potential [J]
    double k;           ///< Wave vector [1/m]
};

static double f(double  energy,
                void   *params);

/**
 * Handler for command-line options
 */
class EFKPSLOptions : public Options
{
    public:
        EFKPSLOptions(int argc, char* argv[])
        {
            try
            {
                program_specific_options->add_options()
                    ("well-width,a", po::value<double>()->default_value(100),
                     "Width of quantum well [angstrom].")
                    
                    ("barrier-width,b", po::value<double>()->default_value(100),
                     "Width of barrier [angstrom].")

                    ("wave-vector,k", po::value<double>()->default_value(0),
                     "Wave-vector expressed relative to pi/L, where L is "
                     "the period length")

                    ("well-mass,m", po::value<double>()->default_value(0.067),
                     "Effective mass in well (relative to that of a free electron)")

                    ("barrier-mass,n", po::value<double>()->default_value(0.067),
                     "Effective mass in barrier (relative to that of a free electron)")

//                    ("output-equations", po::bool_switch()->default_value(false),
//                     "Output the matching equations for the system")

                    ("particle,p", po::value<char>()->default_value('e'),
                     "Particle to be used: 'e', 'h' or 'l'")

                    ("states,s", po::value<size_t>()->default_value(1),
                     "Number of states to find")

                    ("potential", po::value<double>()->default_value(100),
                     "Barrier potential [meV]")
                    ;

                std::string doc("Find the eigenstate of an infinite Kronig-Penney superlattice "
                                "that corresponds to the specified wave-vector. "
                                "The energies are written to the file \"E*.r\", "
//                                "and the "
//                                "wavefunctions are written to \"wf_*i.r\" "
                                "where the '*' "
                                "is replaced by the particle ID in each case "
//                                "and the "
//                                "'i' is replaced by the number of the state"
                        );

                add_prog_specific_options_and_parse(argc, argv, doc);	
            }
            catch(std::exception &e)
            {
                std::cerr << e.what() << std::endl;
                exit(EXIT_FAILURE);
            }
        }

        /// \returns the width of the quantum well [m]
        double get_well_width() const {return vm["well-width"].as<double>()*1e-10;}
        
        /// \returns the width of the barriers [m]
        double get_barrier_width() const {return vm["barrier-width"].as<double>()*1e-10;}
        
        /// \returns the wave vector [1/m]
        double get_wave_vector() const {
            const double k = vm["wave-vector"].as<double>();
            const double a = get_well_width();
            const double b = get_barrier_width();
            const double L = a + b; // Total period length [m]

            return k * pi/L;
        }

//        /// \returns true if we are required to output the matching equations for the system
//        bool output_equations() const {return vm["output-equations"].as<bool>();}

        /// \returns the effective mass in the quantum well [kg]
        double get_well_mass() const {return vm["well-mass"].as<double>()*me;}

        /// \returns the effective mass in the quantum well [kg]
        double get_barrier_mass() const {return vm["barrier-mass"].as<double>()*me;}

        /// \returns the particle ID
        char get_particle() const {return vm["particle"].as<char>();}

        /// \returns the number of spatial points
        size_t get_n_states() const {return vm["states"].as<size_t>();}

        /// \returns the effective mass in the quantum well [J]
        double get_potential() const {return vm["potential"].as<double>()*1e-3*e;}
};

int main(int argc,char *argv[])
{
    const EFKPSLOptions opt(argc, argv);

    const double a   = opt.get_well_width();    // [m]
    const double b   = opt.get_barrier_width(); // [m]
    const double k   = opt.get_wave_vector();   // [1/m]
    const double m_w = opt.get_well_mass();     // [kg]
    const double m_b = opt.get_barrier_mass();  // [kg]
    const double V   = opt.get_potential();     // [J]
    const size_t nst = opt.get_n_states();
    const char   p   = opt.get_particle();

    kpsl_params params = {a, b, m_w, m_b, V, k};
    gsl_function F;
    F.function = &f;
    F.params   = &params;
    gsl_root_fsolver *solver = gsl_root_fsolver_alloc(gsl_root_fsolver_brent);

    const double dx=1e-3*e; // arbitrarily small energy increment---0.1meV

    double Elo=dx;    // first energy estimate
    std::valarray<double> E(nst);

    for(unsigned int ist=0; ist<nst; ++ist)
    {
        // Shift the lower estimate up past the last state we found
        if(ist > 0) Elo = E[ist-1] + dx;

        // Value for y=f(x) at bottom of search range
        const double y1 = GSL_FN_EVAL(&F,Elo);

        // Find the range in which the solution lies by incrementing the
        // upper limit of the search range until the function changes sign.
        // Since this coarse search can require many iterations, we keep the
        // lower limit fixed to minimise the amount of computation at each step.
        // The Brent algorithm is extremely fast, so it really doesn't matter that
        // the range we find here is large.
        //
        // Note the end stop to prevent infinite loop in absence of solution
        //
        // TODO: Make the cut-off configurable
        double y2 = y1;
        double Ehi = Elo;
        do
        {
            Ehi+=dx;
            y2=GSL_FN_EVAL(&F, Ehi);
        }while((y1*y2>0)&&(Ehi<100*V));

        E[ist] = (Elo + Ehi)/2;
        gsl_root_fsolver_set(solver, &F, Elo, Ehi);
        int status = 0;

        // Improve the estimate of the solution using the Brent algorithm
        // until we hit a desired level of precision
        do
        {
            status = gsl_root_fsolver_iterate(solver);
            E[ist] = gsl_root_fsolver_root(solver);
            Elo = gsl_root_fsolver_x_lower(solver);
            Ehi = gsl_root_fsolver_x_upper(solver);
            status = gsl_root_test_interval(Elo, Ehi, 1e-9*e, 0);
        }while(status == GSL_CONTINUE);
    }

    // Output energies to file
    char filename[9];
    sprintf(filename,"E%c.r",p);
    E *= 1e3/e; // Rescale to meV
    Leeds::write_table_x(filename, E, true);

    return EXIT_SUCCESS;
}

/**
 * \brief This function is the standard fw result = 0
 *
 * \param[in] energy Local energy
 * \param[in] params Superlattice parameters (cast to kpsl_params)
 */
static double f(double  energy,
                void   *params)
{
    const kpsl_params *p = reinterpret_cast<kpsl_params *>(params);
    const double a   = p->a;
    const double b   = p->b;
    const double m_b = p->m_b;
    const double m_w = p->m_w;
    const double V   = p->V;
    const double k   = p->k;

    double F; // value of function

    // Wave vector inside well [QWWAD3, 2.155]
    const double k_w=sqrt(2*m_w/hBar*energy/hBar);

    if(energy<V)
    {
        // If we're below the barrier, then the wavefunction decays with
        // [QWWAD3, 2.167]
        const double K=sqrt(2*m_b/hBar*(V-energy)/hBar);

        // Matching equation [QWWAD3, 2.169]
        F=cos(k_w*a)*cosh(K*b)-sin(k_w*a)*sinh(K*b)*
            (gsl_pow_2(m_b*k_w)-gsl_pow_2(m_w*K))/(2*m_w*m_b*k_w*K)-cos(k*(a+b));
    }
    else
    {
        // If we're above the barrier, then the wavefunction propagates
        // with [QWWAD3, 2.155]
        const double k_b=sqrt(2*m_b/hBar*(energy-V)/hBar);

        // Matching equation [QWWAD3, 2.166]
        F=cos(k_w*a)*cos(k_b*b)-sin(k_w*a)*sin(k_b*b)*
            (gsl_pow_2(m_b*k_w)+gsl_pow_2(m_w*k_b))/(2*m_w*m_b*k_w*k_b)-cos(k*(a+b));
    }

    return F;
}     
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :
