/**
 * \file   schroedinger-solver-kronig-penney.h
 * \author Alex Valavanis <a.valavanis@leeds.ac.uk>
 * \author Paul Harrison  <p.harrison@shu.ac.uk>
 * \brief  Schroedinger solver for a Kronig-Penney superlattice
 */

#ifndef QWWAD_SCHROEDINGER_SOLVER_KRONIG_PENNEY_H
#define QWWAD_SCHROEDINGER_SOLVER_KRONIG_PENNEY_H

#include "schroedinger-solver.h"

#include <armadillo>

namespace QWWAD
{
/**
 * Schroedinger solver for a Kronig-Penney superlattice
 */
class SchroedingerSolverKronigPenney : public SchroedingerSolver
{
public:
    SchroedingerSolverKronigPenney(double l_w,
                                   double l_b,
                                   double V,
                                   double m_w,
                                   double m_b,
                                   double k,
                                   size_t nz,
                                   size_t nper = 10,
                                   unsigned int nst_max = 1);

    auto get_name() -> std::string override {return "kronig-penney";}

    static auto test_matching(double  energy,
                              void   *params) -> double;

    [[nodiscard]] auto get_lhs(double E) const -> double;
    [[nodiscard]] auto get_rhs() const -> double;

private:
    double _l_w; ///< Width of well [m]
    double _l_b; ///< Width of well [m]
    double _V0;  ///< Well depth [J]
    double _m_w; ///< Effective mass in well [kg]
    double _m_b; ///< Effective mass in barriers [kg]
    double _k;   ///< Wave vector [1/m]
    size_t _nz;  ///< Number of spatial points in 1 period

    auto calculate() -> std::vector<Eigenstate> override;

    [[nodiscard]] auto get_wavefunction(double E) const -> arma::cx_vec;
};
} // namespace
#endif
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :
