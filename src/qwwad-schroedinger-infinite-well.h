/**
 * \file   qwwwad-schroedinger-infinite-well.h
 * \author Alex Valavanis <a.valavanis@leeds.ac.uk>
 * \brief  Declarations for infinite-well solver
 */

#ifndef QWWAD_SCHROEDINGER_INFINITE_WELL_H
#define QWWAD_SCHROEDINGER_INFINITE_WELL_H

#include "qwwad-schroedinger.h"

namespace Leeds {
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
                              const double       V       = 0,
                              const unsigned int nst_max = 0);

    std::string get_name() {return "infinite-square-well";}

private:
    double _me;    ///< Effective mass at band-edge (constant) [kg]
    double _L;     ///< Width of quantum well [m]
    double _alpha; ///< Nonparabolicity [1/J]
    double _V;     ///< Band edge [J]

    void calculate();
};
} // namespace Leeds
#endif // QWWAD_SCHROEDINGER_INFINITE_WELL_H
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :
