/**
 * \file   schroedinger-solver-finite-well.h
 * \author Alex Valavanis <a.valavanis@leeds.ac.uk>
 * \brief  Declarations for Schroedinger solver for a single finite well
 */

#ifndef QWWAD_SCHROEDINGER_SOLVER_FINITE_WELL_H
#define QWWAD_SCHROEDINGER_SOLVER_FINITE_WELL_H

#include "schroedinger-solver.h"

namespace QWWAD
{
/**
 * Schroedinger solver for a finite square well with infinitely thick barriers
 */
class SchroedingerSolverFiniteWell : public SchroedingerSolver
{
public:
    SchroedingerSolverFiniteWell(double l_w,
                                 double l_b,
                                 double V,
                                 double m_w,
                                 double m_b,
                                 size_t nz,
                                 unsigned int nst_max = 0);

    auto get_name() -> std::string override {return "finite-square-well";}

    [[nodiscard]] auto get_u0_max() const -> double;
    [[nodiscard]] auto get_n_bound() const -> size_t;

    static auto test_matching(double v,
                              void   *params) -> double;

    [[nodiscard]] auto get_lhs(double v) const -> double;
    static auto get_rhs(double v) -> double;
private:
    double _l_w; ///< Width of well [m]
    double _V0;  ///< Well depth [J]
    double _m_w; ///< Effective mass in well [kg]
    double _m_b; ///< Effective mass in barriers [kg]

    auto calculate() -> std::vector<Eigenstate> override;
    
    [[nodiscard]] auto get_wavefunction(double E,
                                        bool   odd_parity) const -> arma::cx_vec;
};
} // namespace
#endif
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :
