/**
 * \file   qclsim_poisson_solver.h
 * \brief  Poisson solver declarations
 * \author Jonathan Cooper <el06jdc@leeds.ac.uk>
 * \author Alex Valavanis <a.valavanis@leeds.ac.uk>
 * \date   2013-04-25
 */

#ifndef QCLSIM_POISSON_SOLVER_H
#define QCLSIM_POISSON_SOLVER_H

#if HAVE_CONFIG_H
# include "config.h"
#endif //HAVE_CONFIG_H

#include <valarray>

namespace Leeds {

/**
 * Boundary condition to use for the Poisson equation solver
 */
enum PoissonBoundaryType
{
    /// Zero potential at both ends of the system
    DIRICHLET,

    /// Zero potential at the start of the system
    /// Electric-field identical at each end
    MIXED,

    /// Electric field is zero at both ends of the system
    ZERO_FIELD
};

class Poisson
{
public:
    Poisson(const std::valarray<double>& eps, const double dx, PoissonBoundaryType bt=DIRICHLET);
    
    std::valarray<double> solve(std::valarray<double> phi);
    std::valarray<double> solve(std::valarray<double> phi, double V_drop);
    std::valarray<double> solve_laplace(double V_drop);

private:
    void factorise_dirichlet(const std::valarray<double>& eps, const double dx);
    void factorise_mixed(const std::valarray<double>& eps, const double dx);
    void factorise_zerofield(const std::valarray<double>& eps, const double dx);

    int n; ///< Number of spatial samples
    
    std::valarray<double> diag; ///< Diagonal of Poisson matrix
    std::valarray<double> sub_diag; ///< Sub-diagonal of Poisson matrix
    double corner_point; ///< Value of corner point in matrix resulting form mixed boundary
                         //   conditions

    PoissonBoundaryType boundary_type; ///< Boundary condition type for Poisson solver
};

} // End namespace Leeds
#endif //QCLSIM_POISSON_SOLVER_H
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :
