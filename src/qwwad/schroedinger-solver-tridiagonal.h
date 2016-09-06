/**
 * \file   schroedinger-solver-tridiagonal.h
 * \author Jonathan Cooper <jdc.tas@gmail.com>
 * \author Alex Valavanis <a.valavanis@leeds.ac.uk>
 *
 * \brief  Declarations for Schroedinger solver using tridiagonal matrix
 */

#ifndef QWWAD_SCHROEDINGER_SOLVER_TRIDIAGONAL_H
#define QWWAD_SCHROEDINGER_SOLVER_TRIDIAGONAL_H

#include "schroedinger-solver.h"

namespace QWWAD
{
/**
 * Solver for Schroedinger's equation using a tridiagonal Hamiltonian matrix
 */
class SchroedingerSolverTridiag : public SchroedingerSolver
{
private:
    arma::vec _m;   ///< Effective mass at each point
    arma::vec diag; ///< Diagonal elements of matrix
    arma::vec sub;  ///< Sub-diagonal elements of matrix
public:
    SchroedingerSolverTridiag(const decltype(_m) &me,
                              const decltype(_V) &V,
                              const decltype(_z) &z,
                              const unsigned int  nst_max=0);

    std::string get_name() {return "tridiagonal";}
private:
    void calculate();
};
}
#endif
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :
