#include "donor-energy-minimiser-linear.h"

#include "schroedinger-solver-donor-variable.h"
#include "constants.h"

namespace QWWAD
{
using namespace constants;

/**
 * Find the minimum carrier energy
 */
void DonorEnergyMinimiserLinear::minimise()
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
} // namespace
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :
