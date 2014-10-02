/**
 * \file   qwwad-schroedinger-finite-well.h
 * \author Alex Valavanis <a.valavanis@leeds.ac.uk>
 * \brief  Declarations for Schroedinger solver for a single finite well
 */

#ifndef QWWAD_SCHROEDINGER_FINITE_WELL_H
#define QWWAD_SCHROEDINGER_FINITE_WELL_H

#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

#include "qwwad-schroedinger.h"

namespace Leeds {
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
                                 const bool   alt_KE  = false,
                                 const unsigned int nst_max = 0);

    std::string get_name() {return "finite-square-well";}

    double get_u0_max() const;
    size_t get_n_bound() const;
    double get_lhs(const double v) const;
    double get_rhs(const double v) const;
private:
    double _l_w; ///< Width of well [m]
    double _V0;  ///< Well depth [J]
    double _m_w; ///< Effective mass in well [kg]
    double _m_b; ///< Effective mass in barriers [kg]
    double _m_B; ///< Effective mass for use in boundary conditions [kg]

    void calculate();
    
    std::valarray<double> wavef(const double E,
                                const bool   parity_flag);
};
} // namespace Leeds
#endif
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :
