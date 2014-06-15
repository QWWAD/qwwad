/**
 *  \file     qwwad-schroedinger-finite-well.cpp
 *  \author   Alex Valavanis <a.valavanis@leeds.ac.uk>
 *  \brief    Implementatation of Schrodinger solver functions for finite well
 */

#include "qwwad-schroedinger-finite-well.h"

#include <cmath>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>
#include "qclsim-constants.h"

namespace Leeds
{
using namespace constants;

SchroedingerSolverFiniteWell::SchroedingerSolverFiniteWell(const double l_w,
                                                           const double l_b,
                                                           const double V,
                                                           const double m_w,
                                                           const double m_b,
                                                           const size_t nz,
                                                           const bool   alt_KE,
                                                           const unsigned int nst_max) :
    SchroedingerSolver(std::valarray<double>(nz),
                       std::valarray<double>(nz),
                       nst_max),
    _l_w(l_w),
    _V(V),
    _m_w(m_w),
    _m_b(m_b),
    _m_B(alt_KE?m_b:m_w)
{
    for (unsigned int i_z=0; i_z<nz; i_z++)
        _z[i_z]=i_z*(l_w+2*l_b)/(nz-1)-(l_b+l_w/2);
}

/// Parameters needed for computing square-well eigenstates
struct sqw_params
{
    double a;           ///< Well width [m]
    double m_B;         ///< Boundary mass [kg]
    double m_b;         ///< Barrier mass [kg]
    double m_w;         ///< Well mass [kg]
    double V;           ///< Barrier potential [J]
};

/**
 * \brief Finds the right-hand side of the matching equation for a finite well
 *
 * \param[in] v          Normalised wave-vector for well region
 *
 * \returns The right-hand side of the matching equation
 */
static double SchroedingerSolverFiniteWell_rhs(const double v)
{
    double result = 0;

    // Figure out which branch of the rhs function we're on
    // and whether the wave-function has odd or even symmetry
    const unsigned int branch = floor(v / (pi/2)) + 1;

    const bool odd_parity = (branch % 2 == 0);

    if(odd_parity)
        result = -v*cot(v);
    else
        result = v*tan(v);

    return result;
}

double SchroedingerSolverFiniteWell::get_rhs(const double v) const
{
    return SchroedingerSolverFiniteWell_rhs(v);
}

/**
 * \brief Finds the left-hand side of the matching equation for a finite well
 *
 * \param[in] v    Normalised wave-vector for the well region
 * \param[in] a    Well width [m]
 * \param[in] m_w  Effective mass in well region [kg]
 * \param[in] m_b  Effective mass in barrier region [kg]
 * \param[in] V    Barrier potential [J]
 *
 * \returns the left-hand side of the matching equation
 */
static double SchroedingerSolverFiniteWell_lhs(const double v,
                                               const double a,
                                               const double m_w,
                                               const double m_b,
                                               const double V)
{
    const double k = 2.0*v/a;
    const double E = hBar*hBar*k*k/(2.0*m_w);
    const double u_0_sq = a*a*m_w/(2.0*hBar*hBar)*(m_w*V/m_b + E*(1.0-m_w/m_b));

    double result = 0.0;

    if(gsl_fcmp(v*v, u_0_sq, 1e-12) == -1)
        result = sqrt(u_0_sq - v*v);

    return result;
}

double SchroedingerSolverFiniteWell::get_lhs(const double v) const
{
    return SchroedingerSolverFiniteWell_lhs(v, _l_w, _m_w, _m_B, _V);
}

/**
 * \brief Computes the characteristic equation for a finite well
 *
 * \details This can be used to find the eigenstates for the system,
 *          since the characteristic equation is zero for an eigenstate.
 *          For even parity states, it is given by
 *          \f[
 *            f(E) = k\tan\left(\frac{kl_w}{2}\right) - \kappa
 *          \f]
 *          For odd states, it is
 *          \f[
 *            f(E) = k\cot\left(\frac{kl_w}{2}\right) + \kappa
 *          \f]
 *
 * \param[in] energy local energy
 * \param[in] params square-well parameters
 *
 * \returns The value of the characteristic equation for the well.
 *          This is zero when the energy equals the an eigenvalue.
 */
double SchroedingerSolverFiniteWell_f(double  v,
                                      void   *params)
{
    const sqw_params *p = reinterpret_cast<sqw_params *>(params);
    const double a   = p->a;
    const double m_B = p->m_B;
    const double m_w = p->m_w;
    const double V   = p->V;

    const double lhs = SchroedingerSolverFiniteWell_lhs(v, a, m_w, m_B, V);
    const double rhs = SchroedingerSolverFiniteWell_rhs(v);
    return lhs - rhs;
}

/**
 * \brief calculates the uncorrelated one particle wavefunction
 *
 * \param[in] E           local energy
 * \param[in] odd_parity  true for odd states, false for even
 */
std::valarray<double> SchroedingerSolverFiniteWell::wavef(const double E,
                                                          const bool   odd_parity)
{
    // Define k and K
    const double k=sqrt(2*_m_w/hBar*E/hBar); // wave vector in the well
    const double K=sqrt(2*_m_b/hBar*(_V-E)/hBar); // decay constant in barrier

    const size_t N = _z.size();
    std::valarray<double> psi(N); // wavefunction
    const double dz = _z[1] - _z[0];
    const double epsilon = dz/1000;

    // Compute the normalisation constants
    double A = 0; // Amplitude in well
    double B = 0; // Amplitude in barriers

    if(odd_parity)
    {
        A = 1.0/sqrt(
                _l_w/2.0 - sin(k*_l_w)/(2.0*k) + gsl_pow_2(sin(k*_l_w/2.0))/K
                );
        B = A * exp(K*_l_w/2.0) * sin(k*_l_w/2.0);
    }
    else
    {
        A = 1.0/sqrt(
                _l_w/2.0 + sin(k*_l_w)/(2.0*k) + gsl_pow_2(cos(k*_l_w/2.0))/K
                );
        B = A * exp(K*_l_w/2.0) * cos(k*_l_w/2.0);
    }

    for (unsigned int i_z=0;i_z<N;i_z++)
    {
        // Left barrier
        if (gsl_fcmp(_z[i_z], -_l_w/2, epsilon)==-1)
        {
            if(odd_parity)
                psi[i_z]=-B*exp(-K*fabs(_z[i_z]));
            else
                psi[i_z]=B*exp(-K*fabs(_z[i_z]));
        }
        // Right barrier
        else if (gsl_fcmp(_z[i_z], _l_w/2, epsilon)>=0)
            psi[i_z]=B*exp(-K*_z[i_z]);

        // Find wavefunction within well region
        else if (gsl_fcmp(_z[i_z], -_l_w/2, epsilon) >= 0 && (_z[i_z]<(_l_w/2)))
        {
            if(odd_parity)
                psi[i_z]=A*sin(k*_z[i_z]);
            else
                psi[i_z]=A*cos(k*_z[i_z]);
        }
    }

    return psi;
}

double SchroedingerSolverFiniteWell::get_u0_max() const
{
    return sqrt(_V*_l_w*_l_w*_m_w/(2.0*hBar*hBar));
}

/**
 * \brief Find number of bound states in the well
 */
size_t SchroedingerSolverFiniteWell::get_n_bound() const
{
    // Calculate number of bound states in well
    const double u_0_max = get_u0_max();
    const double nst = ceil(u_0_max/(pi/2.0));
    return nst;
}

void SchroedingerSolverFiniteWell::calculate()
{
    const double u_0_max = get_u0_max();
    const size_t nst = get_n_bound();

    for (unsigned int ist=0; ist < _nst_max && ist < nst; ++ist)
    {
        // deduce parity: false if even parity
        const bool parity_flag = (ist%2 == 1);

        sqw_params params = {_l_w, _m_B, _m_b, _m_w, _V};
        gsl_function F;
        F.function = &SchroedingerSolverFiniteWell_f;
        F.params   = &params;
        gsl_root_fsolver *solver = gsl_root_fsolver_alloc(gsl_root_fsolver_brent);

        // Normally, the root needs to lie within each pi/2 cell so we
        // set the limits for the root accordingly.
        double vlo = ist * pi/2.0;
        double vhi = (ist+0.999999999) * pi/2.0;

        // If this is the highest state in the well, then we need to
        // reduce the range so that the energy doesn't go over the
        // top of the well.
        if (ist == nst - 1)
           vhi = u_0_max;

        double v = 0.5 * (vlo+vhi); // Initial estimate of solution
        gsl_root_fsolver_set(solver, &F, vlo, vhi);
        int status = 0;

        // Improve the estimate of solution using the Brent algorithm
        // until we hit a desired level of precision
        do
        {
            status = gsl_root_fsolver_iterate(solver);
            v = gsl_root_fsolver_root(solver);
            vlo = gsl_root_fsolver_x_lower(solver);
            vhi = gsl_root_fsolver_x_upper(solver);
            status = gsl_root_test_interval(vlo, vhi, 0.001, 0);
        }while(status == GSL_CONTINUE);

        const double k = 2.0*v/_l_w;
        const double E = hBar*hBar*k*k/(2.0*_m_w);

        // Stop if we've exceeded the cut-off energy
        if(_E_cutoff_set && gsl_fcmp(E, _E_cutoff, e*1e-12) == 1)
            break;

        std::valarray<double> psi = wavef(E,parity_flag);
        _solutions.push_back(State(E, psi));
        gsl_root_fsolver_free(solver);
    }
}
} // namespace Leeds
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :
