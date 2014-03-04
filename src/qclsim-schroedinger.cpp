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
                                                     const unsigned int nst_max) :
    SchroedingerSolver(std::valarray<double>(nz),
                       std::valarray<double>(nz),
                       nst_max),
    _me(me),
    _L(L)
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
        // Energy of state [J] (QWWAD3, 2.13)
        double E=gsl_pow_2(pi*hBar*is/_L)/(2*_me);

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
double SchroedingerSolverFiniteWell_f(double  energy,
                                      void   *params)
{
    const sqw_params *p = reinterpret_cast<sqw_params *>(params);
    const double a   = p->a;
    const double m_B = p->m_B;
    const double m_b = p->m_b;
    const double m_w = p->m_w;
    const double V   = p->V;
    const bool   parity_flag = p->parity_flag;

    const double k=sqrt(2*m_w*energy)/hBar; // electron wave vector
    const double K=sqrt(2*m_b*(V-energy))/hBar; // wavefunction decay constant

    double result = 0.0;

    if(parity_flag)
        result = k*cot(k*a/2)/m_w+K/m_B;
    else
        result = k*tan(k*a/2)/m_w-K/m_B;

    return result;
}

/**
 * \brief calculates the uncorrelated one particle wavefunctions for the electron and hole and writes to an external file.
 *
 * \param[in] E           local energy
 * \param[in] i_state     state index
 * \param[in] parity_flag true for odd states, false for even
 */
std::valarray<double> SchroedingerSolverFiniteWell::wavef(const double E,
                                                          const int    i_state,
                                                          const bool   parity_flag)
{
    const size_t N = _z.size();
    double A;         /* In the well the wavefunction psi=Acoskz */
    const double B=1; /* and in the barrier  psi=Bexp(-Kz)       */
    double norm_int;  /* integral over all space of psi*psi      */

    // Define k and K
    const double k=sqrt(2*_m_w/hBar*E/hBar); // wave vector in the well
    const double K=sqrt(2*_m_b/hBar*(_V-E)/hBar); // decay constant in barrier
    
    std::valarray<double> psi(N); // wavefunction

    if(parity_flag) // odd parity wavefunction
    {
        A=B*exp(-K*_l_w/2)/sin(k*_l_w/2);

        for (unsigned int i_z=0;i_z<N;i_z++) // calculate wavefunctions
        {
            if (_z[i_z]<(-_l_w/2)) // Left barrier
            {
                psi[i_z]=-B*exp(-K*fabs(_z[i_z]));
            }
            if ((_z[i_z]>=(-_l_w/2))&&(_z[i_z]<(_l_w/2))) // Well
            {
                psi[i_z]=A*sin(k*_z[i_z]);
            }
            if (_z[i_z]>=(_l_w/2)) // Right barrier
            {
                psi[i_z]=B*exp(-K*_z[i_z]);
            }
        }

        // normalisation integral for odd parity type I
        norm_int=gsl_pow_2(A)*(_l_w/2-sin(k*_l_w)/(2*k))-
            gsl_pow_2(B)*exp(-K*_l_w)*(exp(-2*K*_l_b)-1)/K;
    }
    else // even parity wavefunction
    {
        A=B*exp(-K*_l_w/2)/cos(k*_l_w/2);

        for (unsigned int i_z=0;i_z<N;i_z++) // calculate wavefunctions
        {
            if (_z[i_z]<(-_l_w/2)) // Left barrier
            {
                psi[i_z]=B*exp(-K*fabs(_z[i_z]));
            }
            if ((_z[i_z]>=(-_l_w/2))&&(_z[i_z]<(_l_w/2))) // Well
            {
                psi[i_z]=A*cos(k*_z[i_z]);
            }
            if (_z[i_z]>=(_l_w/2)) // Right barrier
            {
                psi[i_z]=B*exp(-K*_z[i_z]);
            }
        }

        // normalisation integral for even parity type I
        norm_int=gsl_pow_2(A)*(_l_w/2+sin(k*_l_w)/(2*k))+
            gsl_pow_2(B)*exp(-K*_l_w)*(1-exp(-2*K*_l_b))/K;
    }

    // normalise wavefunction
    psi/=sqrt(norm_int);
    return psi;
}

void SchroedingerSolverFiniteWell::calculate()
{
    const double dx = 1e-4*e; // arbitrarily small energy increment [0.1meV]
    double x = dx; // first energy estimate
    for(unsigned int i_state=1; i_state<=_nst_max;i_state++)
    {
        // deduce parity: false if even parity
        const bool parity_flag = (i_state%2 != 1);

        sqw_params params = {_l_w, _m_B, _m_b, _m_w, _V, parity_flag};
        gsl_function F;
        F.function = &SchroedingerSolverFiniteWell_f;
        F.params   = &params;

        gsl_root_fsolver *solver = gsl_root_fsolver_alloc(gsl_root_fsolver_brent);
        int status;
        bool root_found = false;
        double x_end = x;
        bool V_limit_hit = false;

        // Look for a solution in the range (x,x+dx)
        do
        {
            // Set a coarse bracket for search
            x = x_end;
            x_end+=dx;
            double y1=SchroedingerSolverFiniteWell_f(x,     &params);
            double y2=SchroedingerSolverFiniteWell_f(x_end, &params);
            root_found = (y1*y2 < 0);

            // Check that we haven't hit the top of the well
            // (i.e., either y1 or y2 are not-a-number)
            V_limit_hit = std::isunordered(y1,y2);

            if(root_found)
            {
                double E = 0.5 * (x+x_end); // Initial estimate of energy
                double x_lo = x;
                double x_hi = x_end;
                gsl_root_fsolver_set(solver, &F, x_lo, x_hi);

                // Improve the estimate of energy using the Brent algorithm
                // until we hit a desired level of precision
                do
                {
                    status = gsl_root_fsolver_iterate(solver);
                    E = gsl_root_fsolver_root(solver);
                    double x_lo = gsl_root_fsolver_x_lower(solver);
                    double x_hi = gsl_root_fsolver_x_upper(solver);
                    status = gsl_root_test_interval(x_lo, x_hi, 1e-15*e, 0);
                }while(status == GSL_CONTINUE);

                std::valarray<double> psi = wavef(E,i_state,parity_flag);
                _solutions.push_back(State(E, psi));
            }
        }while(!(root_found or V_limit_hit));
        x+=dx;

        gsl_root_fsolver_free(solver);
    }
}
} // namespace Leeds
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :
