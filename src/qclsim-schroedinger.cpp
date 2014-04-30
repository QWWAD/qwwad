/**
 *  \file     qclsim-schroedinger.cpp
 *  \author   Jonathan Cooper 
 *  \author   Alex Valavanis <a.valavanis@leeds.ac.uk>
 *  \brief    Implementatation of Schrodinger solver functions for two-dimensional systems
 */

#include "qclsim-schroedinger.h"

#include <cmath>
#include <gsl/gsl_math.h>
#include "qclsim-constants.h"

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

} // namespace Leeds
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :
