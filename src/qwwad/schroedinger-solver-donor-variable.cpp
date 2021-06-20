/**
 * \file   schroedinger-solver-donor-variable.cpp
 * \author Paul Harrison  <p.harrison@shu.ac.uk>
 * \author Alex Valavanis <a.valavanis@leeds.ac.uk>
 * \brief  Schrodinger solver functions for donor with variable symmetrical hydrogenic wavefunction
 */

#include "schroedinger-solver-donor-variable.h"

#include <cmath>
#include <gsl/gsl_integration.h>
#include "constants.h"

namespace QWWAD
{
using namespace constants;
SchroedingerSolverDonorVariable::SchroedingerSolverDonorVariable(const double    m,
                                                                 const arma::vec &V,
                                                                 const arma::vec &z,
                                                                 const double     eps,
                                                                 const double     r_d,
                                                                 const double     lambda,
                                                                 const double     zeta,
                                                                 const double     dE) :
    SchroedingerSolverDonor(m, V, z, eps, r_d, lambda, dE),
    _zeta(zeta)
{}

/**
 * \brief Computes the binding energy integral \f$I_1\f$ for a 3D trial wavefunction
 *
 * \returns Binding energy integral [m^2]
 *
 * \details See Eq. 5.100, QWWAD3. The integral is solved analytically as
 *          \f[
 *             I_1=2\pi\left(
 *             \frac{\zeta\vert z^\prime\vert\lambda}{2}+\frac{\lambda^2}{4}\right)
 *             \mbox{e}^{-\frac{2\zeta\vert z^\prime\vert}{\lambda}}
 *          \f]
 */
auto SchroedingerSolverDonorVariable::I_1(const double z_dash) const -> double
{
    return 2*pi*(_zeta*fabs(z_dash)*_lambda/2+ _lambda*_lambda/4)*
           exp(-2*_zeta*fabs(z_dash)/_lambda);
}

/**
 * \brief Computes the binding energy integral \f$I_2\f$ for a 3D trial wavefunction
 *
 * \returns Binding energy integral [m]
 *
 * \details See Eq. 5.106, QWWAD3. The integral evaluates to
 *          \f[
 *             I_2=2\pi\left(-\frac{\zeta^2 z^\prime}{2}\right)
 *             \mbox{e}^{-\frac{2\zeta\vert z^\prime\vert}{\lambda}}
 *          \f]
 */
auto SchroedingerSolverDonorVariable::I_2(const double z_dash) const -> double
{
    return 2*pi*(-_zeta*_zeta*z_dash/2)*exp(-2*_zeta*fabs(z_dash)/_lambda);
}

struct integral_params
{
    double lambda;
    double zeta;
    double z_dash_abs;
};

/* Eq. 5.116, QWWAD3 */
auto I_33_integrand(double w, void *params) -> double
{
    auto *p = reinterpret_cast<integral_params *>(params);
    double result = exp(-p->zeta*p->z_dash_abs*(1/w+w)/p->lambda)*(1-w*w)/gsl_pow_2(1+w*w);
    return result;
}

/* Eq. 5.117, QWWAD3 */
auto I_34_integrand(double w, void *params) -> double
{
    auto *p = reinterpret_cast<integral_params *>(params);
    double result = exp(-p->zeta*p->z_dash_abs*(1/w+w)/p->lambda)*(1-w*w)/(w*(1+w*w));
    return result;
}

/**
 * \brief Computes the binding energy integral \f$I_3\f$ for a 3D trial wavefunction
 *
 * \returns Binding energy integral [dimensionless]
 *
 * \details See Eq. 5.112, QWWAD3.
 */
auto SchroedingerSolverDonorVariable::I_3(const double z_dash) const -> double
{
    const double z_dash_abs = fabs(z_dash);

    /* Eq. 5.113, QWWAD3 */
    const double I_31=(-1-_zeta*_zeta)/2*exp(-2*_zeta*fabs(z_dash)/_lambda);
    const double I_32=(_zeta*fabs(z_dash)/(2*_lambda)+0.25)*exp(-2*_zeta*fabs(z_dash)/_lambda);

    integral_params p = {_lambda, _zeta, z_dash_abs};

    gsl_function f_33;
    f_33.function = &I_33_integrand;
    f_33.params = &p;

    gsl_function f_34;
    f_34.function = &I_34_integrand;
    f_34.params = &p;

    // Configure integration algorithm
    const double limit_lo   = 0.0;  // Lower limit for integral
    const double limit_hi   = 1.0;  // Upper limit for integral
    double I_33             = 0.0;  // Result of integral
    double I_34             = 0.0;  // Result of integral
    const double abserr_max = 100;  // Maximum permitted absolute error
    const double relerr_max = 1e-7; // Maximum permitted relative error
    double abserr           = 0.0;  // Absolute error (after calculation)
    size_t neval            = 0;    // Number of times integrand was evaluated

    // Perform the integral
    gsl_integration_qng(&f_33, limit_lo, limit_hi, abserr_max, relerr_max, &I_33, &abserr, &neval);
    gsl_integration_qng(&f_34, limit_lo, limit_hi, abserr_max, relerr_max, &I_34, &abserr, &neval);

    I_33*=2*(gsl_pow_3(_zeta)-_zeta)*fabs(z_dash)/_lambda;
    I_34*=(gsl_pow_4(_zeta)-gsl_pow_2(_zeta))*gsl_pow_2(z_dash/_lambda);

    return 2*pi*(I_31+I_32+I_33+I_34);
}

/* Eq. 5.118, QWWAD3 */
static auto I_4_integrand(double w, void *params) -> double
{
    auto *p = reinterpret_cast<integral_params *>(params);
    double result = exp(-2*p->z_dash_abs*sqrt(gsl_pow_2((1-w*w)/(2*w))+gsl_pow_2(p->zeta))/p->lambda)
            *p->z_dash_abs*(1-w*w)/(2*w*w);
    return result;
}

/**
 * \brief Computes the binding energy integral \f$I_4\f$ for a 3D trial wavefunction
 *
 * \param[in] z_dash displacement between electron and donor in z-direction [m]
 *
 * \returns Binding energy integral [m]
 *
 * \details See Eq. 5.92, QWWAD3. The integral is given by
 *          \f[
 *             I_{4}=2\pi\int_0^1
 *             \exp{\left[-\frac{2\vert z^\prime\vert\sqrt{
 *             \left(\frac{1-w^2}{2w}\right)^2+\zeta^2}}{\lambda}\right]}
 *             \vert z^\prime\vert \frac{1-w^2}{2w^2}\;\;\text{d}w
 *          \f]
 */
auto SchroedingerSolverDonorVariable::I_4(const double z_dash) const -> double
{
    const double z_dash_abs = fabs(z_dash);
    integral_params p = {_lambda, _zeta, z_dash_abs};

    gsl_function f;
    f.function = &I_4_integrand;
    f.params = &p;

    // Configure integration algorithm
    const double limit_lo   = 0.0;  // Lower limit for integral
    const double limit_hi   = 1.0;  // Upper limit for integral
    double I_4              = 0.0;  // Result of integral
    const double abserr_max = 1e-5; // Maximum permitted absolute error
    const double relerr_max = 1e-7; // Maximum permitted relative error
    double abserr           = 0.0;  // Absolute error (after calculation)
    size_t neval            = 0;    // Number of times integrand was evaluated

    // Perform the integral
    gsl_integration_qng(&f, limit_lo, limit_hi, abserr_max, relerr_max, &I_4, &abserr, &neval);

    return 2*pi*I_4;
}

auto
SchroedingerSolverDonorVariable::calculate_psi_from_chi() -> std::vector<Eigenstate>
{
    std::vector<Eigenstate> solutions;

    const auto z = get_z();

    for (unsigned int ist = 0; ist < _solutions_chi.size(); ++ist)
    {
        const auto chi = _solutions_chi[0].get_wavefunction_samples();
        const double E = _solutions_chi[0].get_energy();

        auto const psi = chi*exp(-_zeta*abs(z - _r_d)/_lambda);
        solutions.emplace_back(E,z,psi);
    }

    return solutions;
}

} // namespace
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :
