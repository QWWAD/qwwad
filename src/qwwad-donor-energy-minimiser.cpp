/**
 * \file   qwwad-donor-energy-minimiser.cpp
 *
 * \author Paul Harrison  <p.harrison@shu.ac.uk>
 * \author Alex Valavanis <a.valavanis@leeds.ac.uk>
 *
 * \brief  Minimisation of donor state energy
 */

#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_min.h>
#include "qclsim-constants.h"
#include "qwwad-donor-energy-minimiser.h"
#include "qwwad-schroedinger-donor-2D.h"

namespace Leeds
{
using namespace constants;

DonorEnergyMinimiser::DonorEnergyMinimiser(SchroedingerSolverDonor *se,
                                           const double             lambda_start,
                                           const double             lambda_step,
                                           const double             lambda_stop) :
    _se(se),
    _lambda_start(lambda_start),
    _lambda_step(lambda_step),
    _lambda_stop(lambda_stop)
{}

/**
 * \brief Find the energy of a carrier using a given Bohr radius
 */
double DonorEnergyMinimiser::find_E_at_lambda(double  lambda,
                                              void   *params)
{
    SchroedingerSolverDonor *se = reinterpret_cast<SchroedingerSolverDonor *>(params);
    se->set_lambda(lambda); // Recalculate at given Bohr radius
    std::vector<State> solutions = se->get_solutions();

    return solutions[0].get_E();
}

/**
 * Find the minimum carrier energy, and corresponding Bohr radius using a fast search algorithm
 */
void DonorEnergyMinimiser::find_E_min_fast()
{
    if(_lambda_stop < 0)
        throw std::domain_error("Upper limit on Bohr radius must be set to a positive value using --lambdastop");

    double __lambda_start = _lambda_start;

    // Set up the numerical solver using GSL
    gsl_function f;
    f.function = &find_E_at_lambda;
    f.params   = _se;

    // First perform a very coarse search for a suitable estimate of a starting point
    const double Elo = GSL_FN_EVAL(&f, __lambda_start);
    const double Ehi = GSL_FN_EVAL(&f, _lambda_stop);

    double E0 = Elo + Ehi; // Set initial estimate as being higher than Elo and Ehi
    _se->set_lambda(__lambda_start);

    const double dlambda = (_lambda_stop - __lambda_start)/4; // Separation between endpoints [m]

    // Search for a suitable lambda value until we find which quadrant the mimimum lies in
    do
    {
        if(_se->get_lambda() >= _lambda_stop)
            throw std::domain_error("Can't find a minimum in this range of Bohr radii");

        _se->set_lambda(_se->get_lambda() + dlambda); // Increment the Bohr radius
        E0 = _se->get_solutions()[0].get_E();
    }
    while((E0 > Elo) || (E0 > Ehi));
    __lambda_start = _se->get_lambda() - dlambda;

    gsl_min_fminimizer *s = gsl_min_fminimizer_alloc(gsl_min_fminimizer_brent);
    gsl_min_fminimizer_set(s, &f, _se->get_lambda(), __lambda_start, _lambda_stop);

    size_t max_iter = 100; // Maximum number of iterations before giving up
    int status = 0;        // Error flag for GSL
    unsigned int iter=0;   // The number of iterations attempted so far

    // Variational calculation (search over lambda)
    do
    {
        ++iter;
        status  = gsl_min_fminimizer_iterate(s);
        const double lambda_lo = gsl_min_fminimizer_x_lower(s);
        const double lambda_hi = gsl_min_fminimizer_x_upper(s);
        status  = gsl_min_test_interval(lambda_lo, lambda_hi, 0.1e-10, 0.0);
        printf("r_d %le lambda %le energy %le meV\n", _se->get_r_d(), _se->get_lambda(), _se->get_solutions()[0].get_E()/(1e-3*e));
    }while((status == GSL_CONTINUE) && (iter < max_iter));

    gsl_min_fminimizer_free(s);
}

/**
 * Find the minimum carrier energy using a linear search
 */
void DonorEnergyMinimiser::find_E_min_linear()
{
    double lambda=_lambda_start; // Initial Bohr radius value [m]
    double E0 = 1e6*e;           // Set minimum energy of single donor to enormous energy [J]

    // Variational calculation (search over lambda)
    do
    {
        const double E = find_E_at_lambda(lambda, _se);
        printf("r_d %le lambda %le energy %le meV\n", _se->get_r_d(),lambda,E/(1e-3*e));

        if (E > E0) break; // Stop looping if we've passed the minimum

        lambda+=_lambda_step; // increments Bohr radius
    }while((_lambda_stop < 0) || lambda < _lambda_stop);

    // We've overshot the minimum, so go back to the previous radius
    _se->set_lambda(lambda-_lambda_step);
}

void DonorEnergyMinimiser::minimise(MinimisationMethod method)
{
    switch(method)
    {
        case MINIMISE_LINEAR:
            find_E_min_linear();
            break;
        case MINIMISE_FAST:
            find_E_min_fast();
            break;
    }
}
} // Namespace Leeds
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :
