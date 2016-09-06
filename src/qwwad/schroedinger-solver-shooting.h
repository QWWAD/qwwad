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
    double                _dE;    ///< Minimum energy separation between states [J]

public:
    SchroedingerSolverShooting(const decltype(_me)    &me,
                               const decltype(_alpha) &alpha,
                               const decltype(_V)     &V,
                               const decltype(_z)     &z,
                               const double            dE,
                               const unsigned int      nst_max=0);

    std::string get_name() {return "shooting";}

    std::vector<Eigenstate> get_solutions_chi(const bool convert_to_meV=false);

    static double psi_at_inf(double  E,
                             void   *params);

    double shoot_wavefunction(arma::vec    &wf,
                              const double  E) const;

private:
    void calculate();
};
} // namespace
#endif
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :
