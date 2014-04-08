/**
 *  \file     qclsim-schroedinger.cpp
 *  \author   Jonathan Cooper 
 *  \author   Alex Valavanis <a.valavanis@leeds.ac.uk>
 *  \brief    Implementatation of Schrodinger solver functions for two-dimensional systems
 */

#include "qclsim-schroedinger.h"

#include <cmath>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>
#include <boost/multi_array.hpp>
#include "qclsim-constants.h"
typedef boost::multi_array<double, 2> matrix;

namespace Leeds
{
using namespace constants;

/**
 * Create tridiagonal Hamiltonian
 * \param[in] nst_max Maximum number of states to find
 *
 * \details If nst_max=0 (the default), all states will be found
 *          that lie within the range of the input potential profile
 */
SchroedingerSolverTridiag::SchroedingerSolverTridiag(const std::valarray<double>& me,
                                                     const std::valarray<double>& V,
                                                     const std::valarray<double>& z,
                                                     const unsigned int           nst_max) :
    SchroedingerSolver(V,z,nst_max),
    diag(std::valarray<double>(z.size())),
    sub(std::valarray<double>(z.size()))
{
    const size_t nz = z.size();
    const double dz = z[1] - z[0];

    for(unsigned int i=0; i<nz; i++){
        double m_minus;
        double m_plus;

        // Calculate mass midpoints for +1/2 and -1/2 avoiding outside addressing
        if(i==0 || i==nz-1){
            m_minus = m_plus = me[i];
        }
        else{
            m_minus = (me[i] + me[i-1])/2;
            m_plus = (me[i+1] + me[i])/2;
        }

        // Calculate a points
        if(i!=nz-1) sub[i] = -gsl_pow_2(hBar/dz)/(2*m_plus);

        // Calculate b points
        diag[i] = 0.5*gsl_pow_2(hBar/dz)*(m_plus+m_minus)/(m_plus*m_minus) + V[i];
    }
}

/**
 * Find solution to eigenvalue problem
 */
void SchroedingerSolverTridiag::calculate()
{
    _solutions = eigen_tridiag(&diag[0], &sub[0], _V.min(), _V.max(), _V.size(), _nst_max);
}

/**
 * Build matrix 'A' from general eigenproblem
 * \param[in] nst_max Maximum number of states to find
 *
 * \details If nst_max=0 (the default), all states will be found
 *          that lie within the range of the input potential profile
 */
SchroedingerSolverFull::SchroedingerSolverFull(const std::valarray<double>& me,
                                               const std::valarray<double>& alpha,
                                               const std::valarray<double>& V,
                                               const std::valarray<double>& z,
                                               const unsigned int           nst_max) :
    SchroedingerSolver(V,z,nst_max),
    A(std::valarray<double>(9*z.size()*z.size()))
{
    const size_t nz = z.size();
    const double dz = z[1] - z[0];
    std::valarray<double> B(3*nz*nz);

    for(unsigned int i=0; i < nz; i++){
        double m_minus;
        double m_plus;
        double alpha_minus;
        double alpha_plus;
        double V_minus;
        double V_plus;

        // Calculate mass midpoints for +1/2 and -1/2 avoiding outside addressing if neccessary
        if (i==0 or i==nz-1)
        {
            m_minus = m_plus = me[i];
            alpha_minus = alpha_plus = alpha[i];
            V_minus = V_plus = V[i];
        }
        else{
            m_minus = (me[i] + me[i-1])/2;
            m_plus = (me[i+1] + me[i])/2;
            alpha_minus = (alpha[i] + alpha[i-1])/2;
            alpha_plus = (alpha[i+1] + alpha[i])/2;
            V_minus = (V[i] + V[i-1])/2;
            V_plus = (V[i+1] + V[i])/2;
        }

        // Calculate a points
        if(i!=0)
            A[(2*nz)+i+((i-1)*3*nz)] = -0.5*gsl_pow_2(hBar/dz)*\
                                       (1-alpha_plus*V_plus)/(m_minus*alpha_plus*alpha_minus);

        // Calculate b points
        A[(2*nz)+i+(i*3*nz)] = 0.5*gsl_pow_2(hBar/dz)/(alpha_plus*alpha_minus)*
            ((1-alpha_minus*V_minus)/m_plus + (1-alpha_plus*V_plus)/m_minus)
            + V[i]*(1 - alpha_minus*V_minus - alpha_plus*V_plus +\
                    alpha_plus*alpha_minus*V_plus*V_minus)/(alpha_plus*alpha_minus);

        // Calculate c points
        if(i!=nz-1)
            A[(2*nz)+i+((i+1)*3*nz)] = -0.5*gsl_pow_2(hBar/dz)*\
                                       (1-alpha_minus*V_minus)/(m_plus*alpha_plus*alpha_minus);

        // Calculate d points
        if(i!=0)
            A[(3*nz*nz)+(2*nz)+i+((i-1)*3*nz)] = -0.5*gsl_pow_2(hBar/dz)/(m_minus*alpha_minus);

        // Calculate e points
        A[(3*nz*nz)+(2*nz)+i+(i*3*nz)] = 0.5*gsl_pow_2(hBar/dz)*
            (1/(m_plus*alpha_plus) + 1/(m_minus*alpha_minus))\
            - (1 - alpha_minus*(V_minus+V[i]) - alpha_plus*(V_plus+V[i])\
                    + alpha_plus*alpha_minus*(V_plus*V_minus+V[i]*V_plus+V[i]*V_minus))/(alpha_plus*alpha_minus);

        // Calculate f points
        if(i!=nz-1)
            A[(3*nz*nz)+(2*nz)+i+((i+1)*3*nz)] = -0.5*gsl_pow_2(hBar/dz)/(m_plus*alpha_plus);

        // Calcualte g points
        A[(6*nz*nz)+(2*nz)+i+(i*3*nz)] = -1/alpha_plus - 1/alpha_minus + V_plus + V[i]+V_minus;

        // Insert identity matrices
        A[(3*nz*nz)+i+(i*3*nz)] = 1;
        A[(6*nz*nz)+nz+i+(i*3*nz)] = 1;
    }
}

/**
 * Find solution to eigenvalue problem
 */
void SchroedingerSolverFull::calculate()
{
    // Find solutions, including all the unwanted "padding" in the eigenvector
    // that comes from the cubic EVP.  See J. Cooper et al., APL 2010
    const std::vector<State> solutions_tmp = eigen_general(&A[0], _V.min(), _V.max(), 3*_z.size(), _nst_max);

    // Now chop off the padding from the eigenvector
    const size_t nst = solutions_tmp.size();
    const size_t nz  = _z.size();

    _solutions.resize(nst, State(nz));
    
    for(unsigned int ist = 0; ist < nst; ist++)
    {
        const double E = solutions_tmp[ist].get_E();
        const std::valarray<double> psi = solutions_tmp[ist].psi_array()[std::slice(0, nz, 1)];

        _solutions[ist] = State(E, psi);
    }
}

/**
 * Get the solutions to the Schroedinger equation.
 *
 * \details The solutions are computed on the first call to this function, but
 *          subsequent calls just recall the values and are hence much faster.
 */
std::vector<State> SchroedingerSolver::get_solutions(const bool convert_to_meV)
{
    // Only calculate if we haven't done so yet
    if(_solutions.empty())
    {
        calculate();

        // Normalise wavefunctions
        for(std::vector<State>::iterator ist = _solutions.begin(); ist != _solutions.end(); ++ist)
            ist->normalise(_z);
    }

    if(convert_to_meV)
    {
        std::vector<State> sol_meV;

        for(std::vector<State>::iterator sol_J = _solutions.begin(); sol_J != _solutions.end(); ++sol_J)
            sol_meV.push_back(State(sol_J->get_E()*1000/e, sol_J->psi_array()));

        return sol_meV;
    }
    else
        return _solutions;
}


/**
 * Assign arrays of potential and position
 *
 * \param[in] nst_max Maximum number of states to find
 *
 * \details If nst_max=0 (the default), all states will be found
 *          that lie within the range of the input potential profile
 */
SchroedingerSolver::SchroedingerSolver(const std::valarray<double> &V,
                                       const std::valarray<double> &z,
                                       const unsigned int           nst_max) :
    _V(V),
    _z(z),
    _nst_max(nst_max),
    _solutions()
{}

/**
 * Create discretised Hamiltonian for system
 *
 * \param[in] me      Effective mass
 * \param[in] alpha   Band nonparabolicity
 * \param[in] V       Confining potential
 * \param[in] z       Spatial coordinates
 * \param[in] nst_max Maximum number of states to find
 *
 * \details If nst_max=0 (the default), all states will be found
 *          that lie within the range of the input potential profile
 */
SchroedingerSolverTaylor::SchroedingerSolverTaylor(const std::valarray<double> &me,
                                       const std::valarray<double> &alpha,
                                       const std::valarray<double> &V,
                                       const std::valarray<double> &z,
                                       const unsigned int           nst_max) :
    SchroedingerSolver(V,z,nst_max),
    AB(std::valarray<double>(2*z.size())),
    BB(std::valarray<double>(2*z.size()))
{
    const size_t nz = z.size();
    const double dz = z[1] - z[0];

    for(unsigned int i=0; i<nz; i++){
        double m_minus;
        double m_plus;
        double alpha_minus;
        double alpha_plus;
        double V_minus;
        double V_plus;

        // Calculate mass midpoints for +1/2 and -1/2 avoiding outside addressing
        if(i==0 or i==nz-1)
        {
            m_minus     = m_plus     = me[i];
            alpha_minus = alpha_plus = alpha[i];
            V_minus     = V_plus     = V[i];
        }
        else
        {
            m_minus     = (me[i] + me[i-1])/2;
            m_plus      = (me[i+1] + me[i])/2;
            alpha_minus = (alpha[i] + alpha[i-1])/2;
            alpha_plus  = (alpha[i+1] + alpha[i])/2;
            V_minus     = (V[i] + V[i-1])/2;
            V_plus      = (V[i+1] + V[i])/2;
        }

        if(i!=nz-1)
        {
            // Calculate a points
            AB[2*(i+1)] = -0.5*gsl_pow_2(hBar/dz)*(1+alpha_plus*V_plus)/m_plus;

            // Calculate d points
            BB[2*(i+1)] = -0.5*gsl_pow_2(hBar/dz)*alpha_plus/m_plus;
        }

        // Calculate b points
        AB[1+(2*i)] = 0.5*gsl_pow_2(hBar/dz)*((1.0+alpha_plus*V_plus)/m_plus + (1.0+alpha_minus*V_minus)/m_minus) + V[i];

        // Calculate e points
        BB[1+(2*i)] = 0.5*gsl_pow_2(hBar/dz)*(alpha_plus/m_plus + alpha_minus/m_minus) + 1;
    }
}

/**
 * Find solutions to Schroedinger's equation for this Hamiltonian
 */
void SchroedingerSolverTaylor::calculate()
{
    _solutions = eigen_banded(&AB[0], &BB[0], _V.min(), _V.max(), _V.size(), _nst_max);
}

SchroedingerSolverInfWell::SchroedingerSolverInfWell(const double       me,
                                                     const double       L,
                                                     const size_t       nz,
                                                     const double       alpha,
                                                     const double       V,
                                                     const unsigned int nst_max) :
    SchroedingerSolver(std::valarray<double>(nz),
                       std::valarray<double>(nz),
                       nst_max),
    _me(me),
    _L(L),
    _alpha(alpha),
    _V(V)
{
    const double dz = L/(nz-1); // Spatial step

    for(unsigned int iz = 0; iz < nz; ++iz)
        _z[iz] = iz*dz;
}

void SchroedingerSolverInfWell::calculate()
{
    double nz = _z.size();
    _solutions.resize(_nst_max, State(nz));

    // Loop over all required states
    for(unsigned int is=1; is<=_nst_max; is++)
    {
        double E = 0;

        // Energy of state [J] (QWWAD3, 2.13)
        if(gsl_fcmp(_alpha, 0, 1e-6) == 1)
            E = 1.0/(2.0*_alpha*_me*_L) * (
                    sqrt(_me*(
                            2.0*_alpha*gsl_pow_2(hBar*pi*is) + _me*_L*_L
                            )
                        )
                    -_me*_L
                    ) + _V;
        else
            E = gsl_pow_2(pi*hBar*is/_L)/(2*_me) + _V;

        std::valarray<double> psi(nz); // Wavefunction amplitude at each point [m^{-0.5}]

        // Loop over spatial locations and find wavefunction
        // amplitude at each point (QWWAD3, 2.15)
        for(unsigned int i=0;i<nz;i++)
            psi=sqrt(2/_L)*sin(is*pi*_z/_L); // Wavefunction [m^{-0.5}]

        _solutions[is-1] = State(E, psi);
    }
}

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
    _l_b(l_b),
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
    bool   parity_flag; ///< True for odd states
};

/**
 * \brief Finds the right-hand side of the matching equation for a finite well
 *
 * \param[in] v          Normalised wave-vector for well region
 * \param[in] odd_parity True if the required state has odd parity
 *
 * \returns The right-hand side of the matching equation
 */
double SchroedingerSolverFiniteWell_rhs(const double v,
                                        const bool   odd_parity)
{
    double result = 0;

    if(odd_parity)
        result = -v*cot(v);
    else
        result = v*tan(v);

    return result;
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
double SchroedingerSolverFiniteWell_lhs(const double v,
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
    const bool   parity_flag = p->parity_flag;

    const double lhs = SchroedingerSolverFiniteWell_lhs(v, a, m_w, m_B, V);
    const double rhs = SchroedingerSolverFiniteWell_rhs(v, parity_flag);
    return lhs - rhs;
}

/**
 * \brief calculates the uncorrelated one particle wavefunctions for the electron and hole and writes to an external file.
 *
 * \param[in] E           local energy
 * \param[in] i_state     state index
 * \param[in] odd_parity  true for odd states, false for even
 */
std::valarray<double> SchroedingerSolverFiniteWell::wavef(const double E,
                                                          const int    i_state,
                                                          const bool   odd_parity)
{
    // Define k and K
    const double k=sqrt(2*_m_w/hBar*E/hBar); // wave vector in the well
    const double K=sqrt(2*_m_b/hBar*(_V-E)/hBar); // decay constant in barrier

    // Determine whether the decay into the barriers is going to underflow
    // out calculations.  If it is, then just treat it as an infinitely sharp decay
    bool sharp_decay=false;
    if (gsl_fcmp(K*_l_b, 700, 1e-6) == 1)
        sharp_decay=true;

    const size_t N = _z.size();
    std::valarray<double> psi(N); // wavefunction
    const double dz = _z[1] - _z[0];
    const double epsilon = dz/1000;

    // If the wavefunction decays sharply into the barriers, assume that the 
    // entire probability density lies in the well and normalise to 1.
    //
    // Otherwise, we let the amplitude = 1 in the barriers and find the
    // amplitude in the well and the normalisation constant for the system.
    double A = 0;
    double norm_int=1.0;  // integral over all space of psi*psi
    if(sharp_decay)
    {
        if(odd_parity)
            A = sqrt(2.0*k/(_l_w*k - sin(_l_w*k)));
        else
            A = sqrt(2.0*k/(_l_w*k + sin(_l_w*k)));
    }
    else
    {
        if(odd_parity)
        {
            A=exp(-K*_l_w/2)/sin(k*_l_w/2);

            norm_int=gsl_pow_2(A)*(_l_w/2-sin(k*_l_w)/(2*k))
                - exp(-K*_l_w)*gsl_expm1(-2*K*_l_b)/K;
        }
        else
        {
            A=exp(-K*_l_w/2)/cos(k*_l_w/2);

            norm_int=gsl_pow_2(A)*(_l_w/2+sin(k*_l_w)/(2*k))
                - exp(-K*_l_w)*gsl_expm1(-2*K*_l_b)/K;
        }
    }

    for (unsigned int i_z=0;i_z<N;i_z++)
    {
        // Fill in the barrier decays only if the decay constant is
        // small enough to compute them accurately.  Otherwise, just
        // leave the barrier wavefunction as zero
        if(!sharp_decay)
        {
            // Left barrier
            if (gsl_fcmp(_z[i_z], -_l_w/2, epsilon)==-1)
            {
                if(odd_parity)
                    psi[i_z]=-exp(-K*fabs(_z[i_z]));
                else
                    psi[i_z]=exp(-K*fabs(_z[i_z]));
            }
            // Right barrier
            else if (gsl_fcmp(_z[i_z], _l_w/2, epsilon)>=0)
                psi[i_z]=exp(-K*_z[i_z]);
        }

        // Find wavefunction within well region
        if (gsl_fcmp(_z[i_z], -_l_w/2, epsilon) >= 0 && (_z[i_z]<(_l_w/2)))
        {
            if(odd_parity)
                psi[i_z]=A*sin(k*_z[i_z]);
            else
                psi[i_z]=A*cos(k*_z[i_z]);
        }
    }

    // normalise wavefunction
    psi/=sqrt(norm_int);
    return psi;
}

void SchroedingerSolverFiniteWell::calculate()
{
    // Calculate number of bound states in well
    double u_0_max = sqrt(_V*_l_w*_l_w*_m_w/(2.0*hBar*hBar));
    const double nst = ceil(u_0_max/(pi/2.0));

    for (unsigned int ist=0; ist < _nst_max && ist < nst; ++ist)
    {
        // deduce parity: false if even parity
        const bool parity_flag = (ist%2 == 1);

        sqw_params params = {_l_w, _m_B, _m_b, _m_w, _V, parity_flag};
        gsl_function F;
        F.function = &SchroedingerSolverFiniteWell_f;
        F.params   = &params;
        gsl_root_fsolver *solver = gsl_root_fsolver_alloc(gsl_root_fsolver_brent);

        // Normally, the root needs to lie within each pi/2 cell so we
        // set the limits for the root accordingly.
        double vlo = ist * pi/2.0;
        double vhi = (ist+1) * pi/2.0;

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
        std::valarray<double> psi = wavef(E,ist+1,parity_flag);
        _solutions.push_back(State(E, psi));
        gsl_root_fsolver_free(solver);
    }
}
} // namespace Leeds
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :
