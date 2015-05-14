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
public:
    SchroedingerSolverShooting(const std::valarray<double>& me,
                               const std::valarray<double>& alpha,
                               const std::valarray<double>& V,
                               const std::valarray<double>& z,
                               const double                 dE,
                               const unsigned int           nst_max=0);

    std::string get_name() {return "shooting";}

    std::vector<State> get_solutions_chi(const bool convert_to_meV=false);

    static double psi_at_inf(double  E,
                             void   *params);

    double shoot_wavefunction(std::valarray<double> &wf,
                              const double           E) const;

private:
    std::valarray<double> _me;    ///< Band-edge effective mass [kg]
    std::valarray<double> _alpha; ///< Nonparabolicity parameter [J^{-1}]
    double                _dE;    ///< Minimum energy separation between states [J]
    void calculate();
};
} // namespace
#endif
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :
