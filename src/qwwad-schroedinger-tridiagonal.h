/**
 * \file   qwwad-schroedinger-tridiagonal.h
 * \author Jonathan Cooper <jdc.tas@gmail.com>
 * \author Alex Valavanis <a.valavanis@leeds.ac.uk>
 *
 * \brief  Declarations for Schroedinger solver using tridiagonal matrix
 */

#ifndef QWWAD_SCHROEDINGER_TRIDIAGONAL_H
#define QWWAD_SCHROEDINGER_TRIDIAGONAL_H

#include "qclsim-schroedinger.h"

namespace Leeds
{
/**
 * Solver for Schroedinger's equation using a tridiagonal Hamiltonian matrix
 */
class SchroedingerSolverTridiag : public SchroedingerSolver
{
public:
    SchroedingerSolverTridiag(const std::valarray<double>& me,
                              const std::valarray<double>& V,
                              const std::valarray<double>& z,
                              const unsigned int           nst_max=0);

    std::string get_name() {return "tridiagonal";}
private:
    void calculate();
    std::valarray<double> diag; ///< Diagonal elements of matrix
    std::valarray<double> sub; ///< Sub-diagonal elements of matrix
};
}
#endif
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :
