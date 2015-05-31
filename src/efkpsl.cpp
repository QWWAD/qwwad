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
#include "qwwad/options.h"
#include "qwwad/constants.h"
#include "qwwad/file-io.h"

using namespace QWWAD;
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
 * Configure command-line options for the program
 */
Options configure_options(int argc, char* argv[])
{
    Options opt;

    std::string summary("Find the eigenstates of an infinite Kronig-Penney superlattice.");

    opt.add_option<double>("wave-vector,k",     0,   "Wave-vector expressed relative to pi/L, where L is period length");
    opt.add_option<double>("well-width,a",    100,   "Width of quantum well [angstrom].");
    opt.add_option<double>("barrier-width,b", 100,   "Width of barrier [angstrom].");
    opt.add_option<double>("well-mass,m",     0.067, "Effective mass in well (relative to that of a free electron)");
    opt.add_option<double>("barrier-mass,n",  0.067, "Effective mass in barrier (relative to that of a free electron)");
    opt.add_option<char>  ("particle,p",       'e',  "ID of particle to be used: 'e', 'h' or 'l', for "
                                                     "electrons, heavy holes or light holes respectively.");
    opt.add_option<size_t>("nst,s",              1,  "Number of states to find");
    opt.add_option<double>("potential",        100,  "Barrier potential [meV]");
    
    opt.add_prog_specific_options_and_parse(argc, argv, summary);

    return opt;
}

int main(int argc,char *argv[])
{
    const auto opt = configure_options(argc, argv);

    const auto a    = opt.get_option<double>("well-width") * 1e-10;
    const auto b    = opt.get_option<double>("barrier-width") * 1e-10;
    const auto m_w  = opt.get_option<double>("well-mass") * me;
    const auto m_b  = opt.get_option<double>("barrier-mass") * me;
    const auto p    = opt.get_option<char>("particle");         // particle ID (e, h or l)
    const auto V    = opt.get_option<double>("potential") * e / 1000;
    const auto nst  = opt.get_option<size_t>("nst");
    const auto k    = opt.get_option<double>("wave-vector") * pi/(a+b);   // [1/m]

    kpsl_params params = {a, b, m_w, m_b, V, k};
    gsl_function F;
    F.function = &f;
    F.params   = &params;
    auto solver = gsl_root_fsolver_alloc(gsl_root_fsolver_brent);

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
    write_table(filename, E, true);

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
