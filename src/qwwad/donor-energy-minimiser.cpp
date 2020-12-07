/**
 * \file   donor-energy-minimiser.cpp
 *
 * \author Paul Harrison  <p.harrison@shu.ac.uk>
 * \author Alex Valavanis <a.valavanis@leeds.ac.uk>
 *
 * \brief  Minimisation of donor state energy
 */

#include "donor-energy-minimiser.h"

#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_min.h>
#include "schroedinger-solver-donor-2D.h"
#include "schroedinger-solver-donor-variable.h"

namespace QWWAD
{
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
auto DonorEnergyMinimiser::find_E_at_lambda(double  lambda,
                                              void   *params) -> double
{
    auto se = reinterpret_cast<SchroedingerSolverDonor *>(params);
    se->set_lambda(lambda); // Recalculate at given Bohr radius
    const auto solutions = se->get_solutions();

    return solutions[0].get_energy();
}

/**
 * \brief Find the energy of a carrier using a given Bohr radius and symmetry
 */
auto DonorEnergyMinimiser::find_E_at_lambda_zeta(const gsl_vector *lambda_zeta,
                                                   void             *params) -> double
{
    auto se = reinterpret_cast<SchroedingerSolverDonorVariable *>(params);
    const double lambda = gsl_vector_get(lambda_zeta, 0);
    const double zeta   = gsl_vector_get(lambda_zeta, 1);

    se->set_lambda_zeta(lambda, zeta); // Set form-factor
    const auto solutions = se->get_solutions();

    return solutions[0].get_energy();
}
} // namespace
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :
