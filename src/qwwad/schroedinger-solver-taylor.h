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
public:
    SchroedingerSolverTaylor(const std::valarray<double> &me,
                             const std::valarray<double> &alpha,
                             const std::valarray<double> &V,
                             const std::valarray<double> &z,
                             const unsigned int           nst_max=0);
    
    std::string get_name() {return "Taylor";}
private:
    void calculate();
    std::valarray<double> AB; ///< Upper triangle of Hamiltonian matrix
    std::valarray<double> BB; ///< Lower triangle of Hamiltonian matrix
};
} // namespace
#endif
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :
