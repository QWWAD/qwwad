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
public:
    SchroedingerSolverFull(const std::valarray<double>& me,
                           const std::valarray<double>& alpha,
                           const std::valarray<double>& V,
                           const std::valarray<double>& z,
                           const unsigned int           nst_max=0);

    std::string get_name() {return "full";}

private:
    std::valarray<double> A;
    void calculate();
};
} // namespace
#endif
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :
