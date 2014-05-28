/**
 * \file   qwwad-schroedinger-shooting.h
 * \author Alex Valavanis <a.valavanis@leeds.ac.uk>
 * \brief  Declarations for Schroedinger solver using shooting method
 */

#ifndef QWWAD_SCHROEDINGER_SHOOTING_H
#define QWWAD_SCHROEDINGER_SHOOTING_H

#include "qwwad-schroedinger.h"

namespace Leeds {
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

private:
    std::valarray<double> _me;    ///< Band-edge effective mass [kg]
    std::valarray<double> _alpha; ///< Nonparabolicity parameter [J^{-1}]
    double                _dE;    ///< Minimum energy separation between states [J]
    void calculate();
};
} // namespace Leeds
#endif
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :
