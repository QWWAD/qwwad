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

namespace QWWAD
{
/**
 * Boundary condition to use for the Poisson equation solver
 */
enum PoissonBoundaryType
{
    /** Zero potential at both ends of the system */
    DIRICHLET,

    /** Zero potential at the start of the system
     *  Electric-field identical at each end */
    MIXED,

    /** Electric field is zero at both ends of the system */
    ZERO_FIELD
};

class Poisson
{
public:
    Poisson(const std::valarray<double>& eps, const double dx, PoissonBoundaryType bt=DIRICHLET);
    
    std::valarray<double> solve(std::valarray<double> phi);
    std::valarray<double> solve(std::valarray<double> phi, double V_drop);
    std::valarray<double> solve_laplace(const double V_drop);

private:
    void factorise_dirichlet();
    void factorise_mixed();
    void factorise_zerofield();
    void compute_half_index_permittivity();

    std::valarray<double> _eps;       ///< Permittivity at each point [F/m]
    std::valarray<double> _eps_minus; ///< Permittivity half a point to left [F/m]
    std::valarray<double> _eps_plus;  ///< Permittivity half a point to right [F/m]

    double _dx;    ///< Spatial step size [m]
    double _L;     ///< Total length of structure [m]
        
    std::valarray<double> diag;     ///< Diagonal of Poisson matrix
    std::valarray<double> sub_diag; ///< Sub-diagonal of Poisson matrix
    double corner_point;            ///< Corner point in matrix resulting from mixed boundary conditions

    PoissonBoundaryType boundary_type; ///< Boundary condition type for Poisson solver
};
} // namespace
#endif //QCLSIM_POISSON_SOLVER_H
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :
