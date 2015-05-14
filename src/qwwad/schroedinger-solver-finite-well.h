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

    std::string get_name() {return "finite-square-well";}

    double get_u0_max() const;
    size_t get_n_bound() const;

    static double test_matching(double v,
                                void   *params);

    double get_lhs(const double v) const;
    double get_rhs(const double v) const;
private:
    double _l_w; ///< Width of well [m]
    double _V0;  ///< Well depth [J]
    double _m_w; ///< Effective mass in well [kg]
    double _m_b; ///< Effective mass in barriers [kg]

    void calculate();
    
    std::valarray<double> get_wavefunction(const double E,
                                           const bool   parity_flag) const;
};
} // namespace
#endif
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :
