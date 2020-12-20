/**
 * \file   schroedinger-solver-infinite-well.h
 * \author Alex Valavanis <a.valavanis@leeds.ac.uk>
 * \brief  Declarations for infinite-well solver
 */

#ifndef QWWAD_SCHROEDINGER_SOLVER_INFINITE_WELL_H
#define QWWAD_SCHROEDINGER_SOLVER_INFINITE_WELL_H

#include "schroedinger-solver.h"

namespace QWWAD
{
/**
 * Schroedinger solver for an infinite square well
 */
class SchroedingerSolverInfWell : public SchroedingerSolver
{
public:
    SchroedingerSolverInfWell(const double       me,
                              const double       L,
                              const size_t       nz,
                              const double       alpha   = 0,
                              const double       V0      = 0,
                              const unsigned int nst_max = 1);

    void set_padding_width(const double Lb);

    auto get_name() -> std::string override {return "infinite-square-well";}

private:
    double _me;       ///< Effective mass at band-edge (constant) [kg]
    double L_;        ///< Width of quantum well [m]
    double _alpha;    ///< Nonparabolicity [1/J]
    double V0_;       ///< Potential at bottom of well [J]
    size_t _nz;       ///< Number of spatial points
    double Lb_ = 0.0; ///< Width of padding regions [m]

    void make_z_array();
    auto calculate() -> std::vector<Eigenstate> override;
};
} // namespace
#endif
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :
