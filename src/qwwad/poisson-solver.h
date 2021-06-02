/**
 * \file   poisson-solver.h
 * \brief  Poisson solver declarations
 * \author Jonathan Cooper <el06jdc@leeds.ac.uk>
 * \author Alex Valavanis <a.valavanis@leeds.ac.uk>
 * \date   2013-04-25
 */

#ifndef QWWAD_POISSON_SOLVER_H
#define QWWAD_POISSON_SOLVER_H

#if HAVE_CONFIG_H
# include "config.h"
#endif //HAVE_CONFIG_H

#include <armadillo>

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

class PoissonSolver
{
private:
    arma::vec _eps;       ///< Permittivity at each point [F/m]

public:
    PoissonSolver(const decltype(_eps) &eps,
                  double                dx,
                  PoissonBoundaryType   bt=DIRICHLET);
    
    [[nodiscard]] auto solve(const arma::vec &rho) const -> arma::vec;
    [[nodiscard]] auto solve(const arma::vec &rho,
                             double           V_drop) const -> arma::vec;
    [[nodiscard]] auto solve_laplace(double V_drop) const -> arma::vec;

private:
    void factorise_dirichlet();
    void factorise_mixed();
    void factorise_zerofield();
    void compute_half_index_permittivity();

    arma::vec _eps_minus; ///< Permittivity half a point to left [F/m]
    arma::vec _eps_plus;  ///< Permittivity half a point to right [F/m]

    double _dx;    ///< Spatial step size [m]
    double _L;     ///< Total length of structure [m]
        
    arma::vec _diag;     ///< Diagonal of Poisson matrix
    arma::vec _sub_diag; ///< Sub-diagonal of Poisson matrix

    double _corner_point; ///< Corner point in matrix resulting from mixed boundary conditions

    arma::vec _D_diag; ///< Diagonal of factorisation matrix, D
    arma::vec _L_sub;  ///< Subdiagonal of factorisation matrix, L

    PoissonBoundaryType _boundary_type; ///< Boundary condition type for Poisson solver
};
} // namespace
#endif //QWWAD_POISSON_SOLVER_H
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :
