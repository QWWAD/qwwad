#include "donor-energy-minimiser-fast.h"

#include <stdexcept>

#include <gsl/gsl_multimin.h>

#include "schroedinger-solver-donor-variable.h"

namespace QWWAD
{

/**
 * /brief Find the minimum carrier energy, and corresponding Bohr radius
 */
void DonorEnergyMinimiserFast::minimise()
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
} // namespace
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :
