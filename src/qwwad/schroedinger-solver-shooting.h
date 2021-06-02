/**
 * \file   schroedinger-solver-shooting.h
 * \author Alex Valavanis <a.valavanis@leeds.ac.uk>
 * \brief  Declarations for Schroedinger solver using shooting method
 */

#ifndef QWWAD_SCHROEDINGER_SOLVER_SHOOTING_H
#define QWWAD_SCHROEDINGER_SOLVER_SHOOTING_H

#include "schroedinger-solver.h"

namespace QWWAD
{
/**
 * Schroedinger solver that uses a shooting method
 */
class SchroedingerSolverShooting : public SchroedingerSolver
{
private:
    arma::vec _me;    ///< Band-edge effective mass [kg]
    arma::vec _alpha; ///< Nonparabolicity parameter [J^{-1}]
    double    _dE;    ///< Minimum energy separation between states [J]

public:
    SchroedingerSolverShooting(decltype(_me)     me,
                               decltype(_alpha)  alpha,
                               const arma::vec  &V,
                               const arma::vec  &z,
                               double            dE,
                               unsigned int      nst_max=0);

    auto get_name() -> std::string override {return "shooting";}

    auto get_solutions_chi(bool convert_to_meV=false) -> std::vector<Eigenstate>;

    static auto psi_at_inf(double  E,
                           void   *params) -> double;

    auto shoot_wavefunction(arma::vec &wf,
                            double     E) const -> double;

private:
    auto calculate() -> std::vector<Eigenstate> override;
};
} // namespace
#endif
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :
