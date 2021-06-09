/**
 * \file   schroedinger-solver-donor-2D.cpp
 * \author Paul Harrison  <p.harrison@shu.ac.uk>
 * \author Alex Valavanis <a.valavanis@leeds.ac.uk>
 * \brief  Schrodinger solver functions for donor with 2D hydrogenic wavefunction
 */

#include "schroedinger-solver-donor-2D.h"

#include <cmath>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_min.h>
#include <gsl/gsl_roots.h>
#include "constants.h"

namespace QWWAD
{
using namespace constants;

SchroedingerSolverDonor2D::SchroedingerSolverDonor2D(const double     m,
                                                     const arma::vec &V,
                                                     const arma::vec &z,
                                                     const double     eps,
                                                     const double     r_d,
                                                     const double     lambda,
                                                     const double     dE) :
    SchroedingerSolverDonor(m, V, z, eps, r_d, lambda, dE)
{}

/**
 * \brief Computes the binding energy integral \f$I_1\f$ for a 2D trial wavefunction
 *
 * \returns Binding energy integral [m^2]
 *
 * \details See Eq. 5.39, QWWAD3. The integral is solved analytically as
 *          \f[
 *            I_1 = 2\pi\frac{\lambda^2}{4}
 *          \f]
 */
auto SchroedingerSolverDonor2D::I_1(const double /* z_dash */) const -> double
{
    return 2*pi*gsl_pow_2(_lambda)/4;
}

/**
 * \brief Computes the binding energy integral \f$I_2\f$ for a 2D trial wavefunction
 *
 * \returns Binding energy integral [m]
 *
 * \details See Eq. 5.42, QWWAD3. The integral evaluates to zero
 */
auto SchroedingerSolverDonor2D::I_2(const double /* z_dash */) const -> double
{
    return 0;
}

/**
 * \brief Computes the binding energy integral \f$I_3\f$ for a 2D trial wavefunction
 *
 * \returns Binding energy integral [dimensionless]
 *
 * \details See Eq. 5.52, QWWAD3. The integral is solved analytically as
 *          \f[
 *            I_3 = 2\pi \left(-\frac{1}{4}\right)
 *          \f]
 */
auto SchroedingerSolverDonor2D::I_3(const double /* z_dash */) const -> double
{
    return 2*pi*(-0.25);
}

/**
 * \brief Parameters (except w) that appear in Eq. 5.59, QWWAD3
 */
struct I_4_integrand_params
{
    double z_dash_abs; ///< Absolute value of displacement from donor to carrier [m]
    double lambda;     ///< Bohr radius [m]
};

/**
 * \brief Computes the integrand in Eq. 5.59, QWWAD3
 */
auto I_4_integrand(double w, void *params) -> double
{
    auto * p = reinterpret_cast<I_4_integrand_params *>(params);
    return 2*pi*p->z_dash_abs * exp(-p->z_dash_abs * (1/w-w)/p->lambda) *
        (1-w*w)/(2*w*w);
}

/**
 * \brief Computes the binding energy integral \f$I_4\f$ for a 2D trial wavefunction
 *
 * \param[in] z_dash displacement between electron and donor in z-direction [m]
 *
 * \returns Binding energy integral [m]
 *
 * \details See Eq. 5.59, QWWAD3. The integral is given by
 *          \f[
 *              I_4=2\pi\int_0^1 \exp{\left[-\frac{\vert z^\prime\vert\left(\frac{1}{w}-w\right) }{\lambda}\right]}\;\;\vert z^\prime\vert \frac{1-w^2}{2w^2}\;\;\text{d}w
 *          \f]
 */
auto SchroedingerSolverDonor2D::I_4(const double z_dash) const -> double
{
    const double z_dash_abs = fabs(z_dash); // Magnitude of displacement [m]

    // Set up function to be integrated [QWWAD3, 5.59]
    gsl_function f;
    f.function = &I_4_integrand;
    I_4_integrand_params p = {z_dash_abs, _lambda};
    f.params   = &p;

    // Configure integration algorithm
    const double limit_lo   = 0.0;  // Lower limit for integral
    const double limit_hi   = 1.0;  // Upper limit for integral
    double I4               = 0.0;  // Result of integral
    const double abserr_max = 100;  // Maximum permitted absolute error
    const double relerr_max = 1e-7; // Maximum permitted relative error
    double abserr           = 0.0;  // Absolute error (after calculation)
    size_t neval            = 0;    // Number of times integrand was evaluated

    // Perform the integral
    gsl_integration_qng(&f, limit_lo, limit_hi, abserr_max, relerr_max, &I4, &abserr, &neval);

    return I4;
}

auto
SchroedingerSolverDonor2D::calculate_psi_from_chi() -> std::vector<Eigenstate>
{
    std::vector<Eigenstate> solutions;

    for (auto & st : _solutions_chi) {
        solutions.push_back(st);
    }

    return solutions;
}

} // namespace QWWAD
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :
