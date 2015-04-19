/**
 * \file   qwwad-debye.cpp
 *
 * \author Alex Valavanis <a.valavanis@leeds.ac.uk>
 *
 * \brief  A Debye model of specific heat capacity
 */

#include "qwwad-debye.h"

#include <gsl/gsl_math.h>
#include <gsl/gsl_deriv.h>
#include <gsl/gsl_sf_debye.h>
#include "qclsim-constants.h"

namespace Leeds {
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
double DebyeModel::get_internal_energy(const double T)
{
    return 3.0 * Na * kB * T * gsl_sf_debye_3(T_D/T) * natoms/M;
}

// Find the molar specific heat capacity by differentiating
// the internal energy of the system with respect to
// temperature.  Note that this is about 50 times faster
// than directly calculating c, because we are able to
// make use of the super-fast gsl_sf_debye_3 function in
// calculating the specific heat capacity.
double DebyeModel::get_cp(const double T)
{
    gsl_function f;
    f.function = &find_U;
    f.params   = this;

    double cp = 0.0;
    double abserr = 0.0;
    gsl_deriv_forward(&f, T, 1, &cp, &abserr);
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
} // namespace Leeds
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :
