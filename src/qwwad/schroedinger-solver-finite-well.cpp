/**
 *  \file   schroedinger-solver-finite-well.cpp
 *  \author Alex Valavanis <a.valavanis@leeds.ac.uk>
 *  \brief  Implementatation of Schrodinger solver functions for finite well
 */

#include "schroedinger-solver-finite-well.h"

#include <cmath>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>
#include "constants.h"
#include "maths-helpers.h"

namespace QWWAD
{
using namespace constants;
SchroedingerSolverFiniteWell::SchroedingerSolverFiniteWell(const double l_w,
                                                           const double l_b,
                                                           const double V0,
                                                           const double m_w,
                                                           const double m_b,
                                                           const size_t nz,
                                                           const unsigned int nst_max) :
    _l_w(l_w),
    _V0(V0),
    _m_w(m_w),
    _m_b(m_b)
{
    arma::vec V(nz);
    arma::vec z(nz);
    set_nst_max(nst_max);

    // Generate potential profile
    for (unsigned int iz=0; iz<nz; iz++) {
        z[iz]=iz*(l_w+2*l_b)/(nz-1)-(l_b+l_w/2);

        if(z[iz] < -l_w/2 || z[iz] >= l_w/2) {
            V[iz] = _V0;
        }
    }

    set_V(V);
    set_z(z);
}

/**
 * \brief Finds the right-hand side of the matching equation for a finite well
 *
 * \param[in] v          Normalised wave-vector for well region
 *
 * \returns The right-hand side of the matching equation.
 *
 * \details This is independent of the geometry/properties of the well -
 *          i.e., purely a function of wave-vector.  As such, the method
 *          is defined statically, and can be used without instantiating
 *          a SchroedingerSolverFiniteWell object.
 */
auto SchroedingerSolverFiniteWell::get_rhs(double v) -> double
{
    double result = 0;

    // Figure out which branch of the rhs function we're on
    // and whether the wave-function has odd or even symmetry
    const unsigned int branch = floor(v / (pi/2)) + 1;

    const bool odd_parity = (branch % 2 == 0);

    if(odd_parity) {
        result = -v*cot(v);
    } else {
        result = v*tan(v);
    }

    return result;
}

/**
 * \brief Finds the left-hand side of the matching equation
 *
 * \param[in] v    Normalised wave-vector for the well region
 *
 * \returns the left-hand side of the matching equation
 */
auto SchroedingerSolverFiniteWell::get_lhs(const double v) const -> double
{
    const double k = 2.0*v/_l_w;
    const double E = hBar*hBar*k*k/(2.0*_m_w);
    const double u_0_sq = _l_w*_l_w*_m_w/(2.0*hBar*hBar)*(_m_w*_V0/_m_b + E*(1.0-_m_w/_m_b));

    double result = 0.0;

    if(gsl_fcmp(v*v, u_0_sq, 1e-12) == -1) {
        result = sqrt(u_0_sq - v*v);
    }

    return result;
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
 * \param[in] v      Normalised wave-vector
 * \param[in] params A SchroedingerSolverFiniteWell object for which to solve
 *
 * \returns The value of the characteristic equation for the well.
 *          This is zero when the energy equals the an eigenvalue.
 */
auto SchroedingerSolverFiniteWell::test_matching(double  v,
                                                 void   *params) -> double
{
    const SchroedingerSolverFiniteWell *se = reinterpret_cast<SchroedingerSolverFiniteWell *>(params);

    const double lhs = se->get_lhs(v);
    const double rhs = se->get_rhs(v);
    return lhs - rhs;
}

/**
 * \brief calculates the uncorrelated one particle wavefunction
 *
 * \param[in] E           local energy
 * \param[in] odd_parity  true for odd states, false for even
 */
auto SchroedingerSolverFiniteWell::get_wavefunction(const double E,
                                                    const bool   odd_parity) const -> arma::cx_vec
{
    const auto z = get_z();

    // Define k and K
    const double k=sqrt(2*_m_w/hBar*E/hBar); // wave vector in the well
    const double K=sqrt(2*_m_b/hBar*(_V0-E)/hBar); // decay constant in barrier

    const size_t N = z.size();
    arma::cx_vec psi(N); // wavefunction
    const double dz = z[1] - z[0];
    const double epsilon = dz/1000;

    // Compute the normalisation constants
    double A = 0; // Amplitude in well
    double B = 0; // Amplitude in barriers

    if(odd_parity) {
        A = 1.0/sqrt(
                _l_w/2.0 - sin(k*_l_w)/(2.0*k) + gsl_pow_2(sin(k*_l_w/2.0))/K
                );
        B = A * exp(K*_l_w/2.0) * sin(k*_l_w/2.0);
    } else {
        A = 1.0/sqrt(
                _l_w/2.0 + sin(k*_l_w)/(2.0*k) + gsl_pow_2(cos(k*_l_w/2.0))/K
                );
        B = A * exp(K*_l_w/2.0) * cos(k*_l_w/2.0);
    }

    for (unsigned int i_z=0;i_z<N;i_z++) {
        // Left barrier
        if (gsl_fcmp(z[i_z], -_l_w/2, epsilon)==-1) {
            if(odd_parity) {
                psi[i_z]=-B*exp(-K*fabs(z[i_z]));
            } else {
                psi[i_z]=B*exp(-K*fabs(z[i_z]));
            }
        } else if (gsl_fcmp(z[i_z], _l_w/2, epsilon)>=0) {
        // Right barrier
            psi[i_z]=B*exp(-K*z[i_z]);
        } else if (gsl_fcmp(z[i_z], -_l_w/2, epsilon) >= 0 && (z[i_z]<(_l_w/2))) {
        // Find wavefunction within well region
            if(odd_parity) {
                psi[i_z]=A*sin(k*z[i_z]);
            } else {
                psi[i_z]=A*cos(k*z[i_z]);
            }
        }
    }

    return psi;
}

auto SchroedingerSolverFiniteWell::get_u0_max() const -> double
{
    return sqrt(_V0*_l_w*_l_w*_m_w/(2.0*hBar*hBar));
}

/**
 * \brief Find number of bound states in the well
 */
auto SchroedingerSolverFiniteWell::get_n_bound() const -> size_t
{
    // Calculate number of bound states in well
    const double u_0_max = get_u0_max();
    const double nst = ceil(u_0_max/(pi/2.0));
    return nst;
}

auto
SchroedingerSolverFiniteWell::calculate() -> std::vector<Eigenstate>
{
    std::vector<Eigenstate> solutions;
    const double u_0_max = get_u0_max();
    const size_t nst = get_n_bound();
    const auto nst_max = get_nst_max();
    const auto z = get_z();

    for (unsigned int ist=0; ist < nst_max && ist < nst; ++ist) {
        // deduce parity: false if even parity
        const bool parity_flag = (ist%2 == 1);

        gsl_function F;
        F.function = &test_matching;
        F.params   = this;
        gsl_root_fsolver *solver = gsl_root_fsolver_alloc(gsl_root_fsolver_brent);

        // Normally, the root needs to lie within each pi/2 cell so we
        // set the limits for the root accordingly. Note the tiny increments
        // here so that we avoid the asymptote.
        double vlo = (ist+0.000000001) * pi/2.0;
        double vhi = (ist+0.999999999) * pi/2.0;

        // If this is the highest state in the well, then we need to
        // reduce the range so that the energy doesn't go over the
        // top of the well.
        if (ist == nst - 1) {
           vhi = u_0_max;
        }

        double v; // Best estimate of solution
        gsl_root_fsolver_set(solver, &F, vlo, vhi);
        int status = 0;

        // Improve the estimate of solution using the Brent algorithm
        // until we hit a desired level of precision
        do {
            status = gsl_root_fsolver_iterate(solver);

            if(status != 0) {
                std::cerr << "GSL error in SchroedingerSolverFiniteWell: " << std::endl
                          << "   Singularity in range (" << vlo << "," << vhi << ")" << std::endl;
            }

            v = gsl_root_fsolver_root(solver);
            vlo = gsl_root_fsolver_x_lower(solver);
            vhi = gsl_root_fsolver_x_upper(solver);
            status = gsl_root_test_interval(vlo, vhi, 0.001, 0);
        }while(status == GSL_CONTINUE);

        const double k = 2.0*v/_l_w;
        const double E = hBar*hBar*k*k/(2.0*_m_w);

        // Stop if we've exceeded the cut-off energy
        if(energy_above_range(E)) {
            break;
        }

        const arma::cx_vec psi = get_wavefunction(E,parity_flag);

        // Don't store the solution if it's below the minimum energy 
        if(!energy_below_range(E)) {
            solutions.emplace_back(E, z, psi);
        }

        gsl_root_fsolver_free(solver);
    }

    return solutions;
}
} // namespace
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :
