/**
 * \file   qwwad-debye.cpp
 *
 * \author Alex Valavanis <a.valavanis@leeds.ac.uk>
 *
 * \brief  A Debye model of specific heat capacity
 */

#include "qwwad-debye.h"

#include <cmath>

#include <gsl/gsl_deriv.h>
#include <gsl/gsl_sf_debye.h>
#include "qwwad/constants.h"

namespace QWWAD
{
using namespace constants;

DebyeModel::DebyeModel(const double T_D,
                       const double M,
                       const size_t natoms)
    : T_D(T_D),
      M(M),
      natoms(natoms)
{}

/**
 * \brief Find the internal energy at a given temperature
 */
double DebyeModel::get_internal_energy(const double T) const
{
    return 3.0 * Na * kB * T * gsl_sf_debye_3(T_D/T) * natoms/M;
}

// Find the molar specific heat capacity by differentiating
// the internal energy of the system with respect to
// temperature.  Note that this is about 50 times faster
// than directly calculating c, because we are able to
// make use of the super-fast gsl_sf_debye_3 function in
// calculating the specific heat capacity.
double DebyeModel::get_cp(const double T) const
{
    gsl_function f;
    f.function = &find_U;
    f.params   = const_cast<DebyeModel *>(this);

    double cp = 0.0;
    double abserr = 0.0;
    gsl_deriv_forward(&f, T, 1, &cp, &abserr);
    return cp;
}

/**
 * \brief Get specific heat, using low-temperature approximation
 */
double DebyeModel::get_cp_low_T(const double T) const
{
    const double pi_sq = pi*pi;
    return 12*pi_sq*pi_sq*Na*kB*T*T*T/(T_D*T_D*T_D*5)*natoms/M;
}

/**
 * \brief Get specific heat, using high-temperature approximation
 */
double DebyeModel::get_cp_high_T() const
{
    return 3*Na*kB*natoms/M;
}

/**
 * \brief Get quick approximation to specific heat
 *
 * \details Use low or high temperature approximation depending on temperature.
 *          The cross-over point between the two models occurs at
 *          \f[
 *            T_0 = T_{\text{D}} \left(\frac{5}{4\pi^4}\right)^{1/3}
 *          \f]
 *          Note that around this transition temperature, this approximate value
 *          can significantly overestimate the specific heat capacity.
 */
double DebyeModel::get_cp_approx(const double T) const
{
    const double pi_sq = pi*pi;
    double cp = 0.0;

    // Temperature at which high and low-temperature models meet
    const double T_match = T_D * std::cbrt(1.25/(pi_sq*pi_sq));

    if (T > T_match)
        cp = get_cp_high_T();
    else
        cp = get_cp_low_T(T);

    return cp;
}

/**
 * \brief A wrapper for compatibility with GSL
 */
double DebyeModel::find_U(double T, void *params)
{
    DebyeModel *dm = reinterpret_cast<DebyeModel *>(params);
    return dm->get_internal_energy(T);
}
} // namespace
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :
