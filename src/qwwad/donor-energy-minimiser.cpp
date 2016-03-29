/**
 * \file   donor-energy-minimiser.cpp
 *
 * \author Paul Harrison  <p.harrison@shu.ac.uk>
 * \author Alex Valavanis <a.valavanis@leeds.ac.uk>
 *
 * \brief  Minimisation of donor state energy
 */

#include "donor-energy-minimiser.h"

#include <stdexcept>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_min.h>
#include <gsl/gsl_multimin.h>
#include "constants.h"
#include "schroedinger-solver-donor-2D.h"
#include "schroedinger-solver-donor-variable.h"

namespace QWWAD
{
using namespace constants;

DonorEnergyMinimiser::DonorEnergyMinimiser(SchroedingerSolverDonor *se,
                                           const double             lambda_start,
                                           const double             lambda_step,
                                           const double             lambda_stop) :
    _se(se),
    _lambda_start(lambda_start),
    _lambda_step(lambda_step),
    _lambda_stop(lambda_stop),
    _zeta_start(0.0),
    _zeta_step(0.0),
    _zeta_stop(0.0),
    _lambda_history(std::vector<double>(0)),
    _zeta_history(std::vector<double>(0)),
    _E_history(std::vector<double>(0))
{}

/**
 * \brief Find the energy of a carrier using a given Bohr radius
 */
double DonorEnergyMinimiser::find_E_at_lambda(double  lambda,
                                              void   *params)
{
    auto se = reinterpret_cast<SchroedingerSolverDonor *>(params);
    se->set_lambda(lambda); // Recalculate at given Bohr radius
    const auto solutions = se->get_solutions();

    return solutions[0].get_energy();
}

/**
 * \brief Find the energy of a carrier using a given Bohr radius and symmetry
 */
double DonorEnergyMinimiser::find_E_at_lambda_zeta(const gsl_vector *lambda_zeta,
                                                   void             *params)
{
    auto se = reinterpret_cast<SchroedingerSolverDonorVariable *>(params);
    const double lambda = gsl_vector_get(lambda_zeta, 0);
    const double zeta   = gsl_vector_get(lambda_zeta, 1);

    se->set_lambda_zeta(lambda, zeta); // Set form-factor
    const auto solutions = se->get_solutions();

    return solutions[0].get_energy();
}

/**
 * /brief Find the minimum carrier energy, and corresponding Bohr radius using a fast search algorithm
 */
void DonorEnergyMinimiser::find_E_min_fast()
{
    size_t max_iter = 100; // Maximum number of iterations before giving up
    int status = 0;        // Error flag for GSL
    unsigned int iter=0;   // The number of iterations attempted so far

    // See if we're using a variable symmetry form-factor
    SchroedingerSolverDonorVariable *se_variable = dynamic_cast<SchroedingerSolverDonorVariable *>(_se);
    if(se_variable != NULL) // If it's variable symmetry, try a range of symmetry parameters
    {
        gsl_multimin_function f;
        f.f = &find_E_at_lambda_zeta; // Function to minimise
        f.n = 2; // Number of parameters (lambda, zeta)
        f.params = se_variable;

        gsl_vector *lambda_zeta = gsl_vector_alloc(2); // Parameters to pass to the function
        gsl_vector *step_size   = gsl_vector_alloc(2); // Step size for parameters
        gsl_vector_set(lambda_zeta, 0, _lambda_start);
        gsl_vector_set(lambda_zeta, 1, _zeta_start);
        gsl_vector_set(step_size,   0, _lambda_step);
        gsl_vector_set(step_size,   1, _zeta_step);

        gsl_multimin_fminimizer *s = gsl_multimin_fminimizer_alloc(gsl_multimin_fminimizer_nmsimplex2, 2);
        gsl_multimin_fminimizer_set(s, &f, lambda_zeta, step_size);

        do
        {
            ++iter;
            status  = gsl_multimin_fminimizer_iterate(s);
            status  = gsl_multimin_test_size(s->size, 1e-5); // Second number is effectively the abs. tolerance in symmetry parameter

            // Save search history
            _lambda_history.push_back(_se->get_lambda());
            _zeta_history.push_back(se_variable->get_zeta());
            _E_history.push_back(_se->get_solutions()[0].get_energy());
        }while((status == GSL_CONTINUE) && (iter < max_iter));

        gsl_multimin_fminimizer_free(s);
        gsl_vector_free(lambda_zeta);
    }
    else
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
            E0 = _se->get_solutions()[0].get_energy();
        }
        while((E0 > Elo) || (E0 > Ehi));
        __lambda_start = _se->get_lambda() - dlambda;

        gsl_min_fminimizer *s = gsl_min_fminimizer_alloc(gsl_min_fminimizer_brent);
        gsl_min_fminimizer_set(s, &f, _se->get_lambda(), __lambda_start, _lambda_stop);

        // Variational calculation (search over lambda)
        do
        {
            ++iter;
            status  = gsl_min_fminimizer_iterate(s);
            const double lambda_lo = gsl_min_fminimizer_x_lower(s);
            const double lambda_hi = gsl_min_fminimizer_x_upper(s);
            status  = gsl_min_test_interval(lambda_lo, lambda_hi, 0.1e-10, 0.0);

            // Save search history
            _lambda_history.push_back(_se->get_lambda());
            _zeta_history.push_back(0.0);
            _E_history.push_back(_se->get_solutions()[0].get_energy());
        }while((status == GSL_CONTINUE) && (iter < max_iter));

        gsl_min_fminimizer_free(s);
    }
}

/**
 * Find the minimum carrier energy using a linear search
 */
void DonorEnergyMinimiser::find_E_min_linear()
{
    double lambda=_lambda_start; // Initial Bohr radius value [m]
    double E_min = 1e6*e;        // Set minimum energy of single donor to enormous energy [J]
    bool   E_min_passed = false; // True if we've overshot the minimum
    const bool lambda_stop_auto = (_lambda_stop < 0); // Stop looping automatically if the lambda_stop value is negative
    double lambda_min = lambda;  // The Bohr radius at the minimum energy point [m]

    // Variational calculation (search over lambda)
    do
    {
        double E = 10e6*e; // Energy of current solution [J]

        // See if we're using a variable symmetry form-factor
        SchroedingerSolverDonorVariable *se_variable = dynamic_cast<SchroedingerSolverDonorVariable *>(_se);
        if(se_variable != NULL) // If it's variable symmetry, try a range of symmetry parameters
        {
            double zeta = _zeta_start; // Initial symmetry parameter value
            double E_min_zeta = E_min; // Record the best estimate of the binding energy so far [J]
            bool   E_min_zeta_passed = false; // True if we've overshot the minimum
            const bool zeta_stop_auto = (_zeta_stop < 0); // Stop looping automatically if the zeta_stop value is negative
            double zeta_min = zeta; // The symmetry parameter at the minimum energy point

            gsl_vector *lambda_zeta = gsl_vector_alloc(2);
            gsl_vector_set(lambda_zeta, 0, lambda);

            do
            {
                gsl_vector_set(lambda_zeta, 1, zeta);
                E = find_E_at_lambda_zeta(lambda_zeta, se_variable);

                // Save search history
                _lambda_history.push_back(_se->get_lambda());
                _zeta_history.push_back(se_variable->get_zeta());
                _E_history.push_back(_se->get_solutions()[0].get_energy());

                if (E > E_min_zeta)
                    E_min_zeta_passed = true; // Stop looping if we've passed the minimum
                else
                {
                    // Otherwise, record the new minimum energy and corresponding zeta value
                    E_min_zeta = E;
                    zeta_min   = zeta;
                }

                zeta+=_zeta_step; // increments symmetry parameter
            }while((zeta_stop_auto && !E_min_zeta_passed) // Carry on looping if we haven't found the minimum yet (in auto mode)
                    ||
                    (!zeta_stop_auto && (zeta < _zeta_stop)) // or the symmetry parameter is lower than the stop point, and we're not in auto mode
                  );

            se_variable->set_lambda_zeta(lambda, zeta_min);
            gsl_vector_free(lambda_zeta);
        }
        else // If it's a fixed-symmetry solution, just use this Bohr radius
        {
            E = find_E_at_lambda(lambda, _se);

            // Save search history
            _lambda_history.push_back(_se->get_lambda());
            _zeta_history.push_back(0.0);
            _E_history.push_back(_se->get_solutions()[0].get_energy());
        }

        if (E > E_min)
            E_min_passed = true; // Stop looping if we've passed the minimum
        else
        {
            // Otherwise, record the new minimum energy and corresponding lambda value
            E_min      = E;
            lambda_min = lambda;
        }

        lambda+=_lambda_step; // increments Bohr radius
    }while((lambda_stop_auto && !E_min_passed) // Carry on looping if we haven't found the minimum yet (in auto mode)
            ||
           (!lambda_stop_auto && (lambda < _lambda_stop)) // or the Bohr radius is lower than the stop point, and we're not in auto mode
           );

    // Set the Bohr radius to the value that gives minimum energy
    _se->set_lambda(lambda_min);
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
} // namespace
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :
