/**
 * \file   schroedinger-solver-full.h
 * \author Jonathan Cooper <jdc.tas@gmail.com>
 * \author Alex Valavanis <a.valavanis@leeds.ac.uk>
 * \brief  Declarations for Schroedinger solver using full nonparabolic Hamiltonian
 */

#ifndef QWWAD_SCHROEDINGER_SOLVER_FULL_H
#define QWWAD_SCHROEDINGER_SOLVER_FULL_H

#include "schroedinger-solver.h"

namespace QWWAD
{
/**
 * Schroedinger solver that uses a full generalised matrix
 */
class SchroedingerSolverFull : public SchroedingerSolver
{
private:
    arma::vec _m;     ///< Effective mass at each point
    arma::vec _alpha; ///< Non-parabolicity parameter at each point
    arma::mat _A;     ///< Hamiltonian matrix

public:
    SchroedingerSolverFull(decltype(_m)        m,
                           decltype(_alpha)    alpha,
                           const arma::vec    &V,
                           const arma::vec    &z,
                           const unsigned int  nst_max=0);

    auto get_name() -> std::string override {return "full";}

private:
    auto calculate() -> std::vector<Eigenstate> override;
};
} // namespace
#endif
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :
