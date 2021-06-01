/**
 * \file   schroedinger-solver-taylor.h
 * \author Jonathan Cooper <jdc.tas@gmail.com
 * \author Alex Valavanis <a.valavanis@leeds.ac.uk>
 * \brief  Declarations for Schroedinger solver using Taylor expansion
 */

#ifndef QWWAD_SCHROEDINGER_SOLVER_TAYLOR_H
#define QWWAD_SCHROEDINGER_SOLVER_TAYLOR_H

#include "schroedinger-solver.h"

namespace QWWAD
{
/**
 * Schroedinger equation solver (using Taylor approximation)
 */
class SchroedingerSolverTaylor : public SchroedingerSolver
{
private:
    arma::vec _m;     ///< Effective mass at each position
    arma::vec _alpha; ///< Nonparabolicity parameter at each position
    arma::vec AB; ///< Upper triangle of Hamiltonian matrix
    arma::vec BB; ///< Lower triangle of Hamiltonian matrix

public:
    SchroedingerSolverTaylor(const decltype(_m)     &me,
                             const decltype(_alpha) &alpha,
                             const arma::vec        &V,
                             const arma::vec        &z,
                             unsigned int            nst_max=0);
    
    auto get_name() -> std::string override {return "Taylor";}
private:
    auto calculate() -> std::vector<Eigenstate> override;
};
} // namespace
#endif
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :
