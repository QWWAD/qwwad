/**
 * \file   qwwad-schroedinger-donor-variable.cpp
 * \author Paul Harrison  <p.harrison@shu.ac.uk>
 * \author Alex Valavanis <a.valavanis@leeds.ac.uk>
 * \brief  Schrodinger solver functions for donor with variable symmetrical hydrogenic wavefunction
 */

#include "qwwad-schroedinger-donor-variable.h"

#include <cmath>
#include "qclsim-constants.h"

namespace Leeds
{
using namespace constants;

SchroedingerSolverDonorVariable::SchroedingerSolverDonorVariable(const double                 m,
                                                                 const std::valarray<double> &V,
                                                                 const std::valarray<double> &z,
                                                                 const double                 eps,
                                                                 const double                 r_d,
                                                                 const double                 lambda,
                                                                 const double                 zeta,
                                                                 const double                 dE,
                                                                 const unsigned int           nst_max) :
    SchroedingerSolverDonor(m, V, z, eps, r_d, lambda, dE, nst_max),
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
double SchroedingerSolverDonorVariable::I_1(const double z_dash) const
{
    return 2*pi*(_zeta*fabs(z_dash)*_lambda/2+gsl_pow_2(_lambda)/4)*
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
double SchroedingerSolverDonorVariable::I_2(const double z_dash) const
{
    return 2*pi*(-gsl_pow_2(_zeta)*z_dash/2)*exp(-2*_zeta*fabs(z_dash)/_lambda);
}

/**
 * \brief Computes the binding energy integral \f$I_3\f$ for a 3D trial wavefunction
 *
 * \returns Binding energy integral [dimensionless]
 *
 * \details See Eq. 5.112, QWWAD3.
 */
double SchroedingerSolverDonorVariable::I_3(const double z_dash) const
{
    const size_t N_w=100;

    /* Eq. 5.113, QWWAD3 */
    const double I_31=(-1-gsl_pow_2(_zeta))/2*exp(-2*_zeta*fabs(z_dash)/_lambda);
    const double I_32=(_zeta*fabs(z_dash)/(2*_lambda)+0.25)*exp(-2*_zeta*fabs(z_dash)/_lambda);

    /* perform integrations over `w' for I_33 and I_34, area simply given by
       sum of (height of centre of strip * strip width) */
    const double delta_w=(1.0-0.0)/(float)N_w;
    double I_33=0.0;
    double I_34=0.0;
    for(double w=delta_w/2;w<1;w+=delta_w)
    {
        /* Eq. 5.116, QWWAD3 */
        I_33+=exp(-_zeta*fabs(z_dash)*(1/w+w)/_lambda)*(1-w*w)/gsl_pow_2(1+w*w)*delta_w;

        /* Eq. 5.117, QWWAD3 */
        I_34+=exp(-_zeta*fabs(z_dash)*(1/w+w)/_lambda)*(1-gsl_pow_2(w))/(w*(1+w*w))*delta_w;
    }

    I_33*=2*(gsl_pow_3(_zeta)-_zeta)*fabs(z_dash)/_lambda;
    I_34*=(gsl_pow_4(_zeta)-gsl_pow_2(_zeta))*gsl_pow_2(z_dash/_lambda);

    return 2*pi*(I_31+I_32+I_33+I_34);
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
double SchroedingerSolverDonorVariable::I_4(const double z_dash) const
{
    const size_t N_w=100;
    double I_40=0.0;
    const double delta_w=(1.0-0.0)/(float)N_w;

    /* Eq. 5.118, QWWAD3 */
    for (double w=delta_w/2;w<1;w+=delta_w)
    {
        I_40+=exp(-2*fabs(z_dash)*sqrt(gsl_pow_2((1-w*w)/(2*w))+gsl_pow_2(_zeta))/_lambda)
            *fabs(z_dash)*(1-w*w)/(2*w*w)*delta_w;
    }

    return 2*pi*I_40;
}
} // namespace Leeds
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :
