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
    SchroedingerSolverFiniteWell(const double l_w,
                                 const double l_b,
                                 const double V,
                                 const double m_w,
                                 const double m_b,
                                 const size_t nz,
                                 const unsigned int nst_max = 0);

    auto get_name() -> std::string override {return "finite-square-well";}

    [[nodiscard]] auto get_u0_max() const -> double;
    [[nodiscard]] auto get_n_bound() const -> size_t;

    static auto test_matching(double v,
                                void   *params) -> double;

    [[nodiscard]] auto get_lhs(const double v) const -> double;
    [[nodiscard]] auto get_rhs(const double v) const -> double;
private:
    double _l_w; ///< Width of well [m]
    double _V0;  ///< Well depth [J]
    double _m_w; ///< Effective mass in well [kg]
    double _m_b; ///< Effective mass in barriers [kg]

    void calculate() override;
    
    [[nodiscard]] auto get_wavefunction(const double E,
                                             const bool   parity_flag) const -> arma::vec;
};
} // namespace
#endif
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :
