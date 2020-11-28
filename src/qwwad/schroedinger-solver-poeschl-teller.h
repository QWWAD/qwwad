/**
 * \file   schroedinger-solver-poeschl-teller.h
 * \author Alex Valavanis <a.valavanis@leeds.ac.uk>
 * \brief  Schroedinger solver for a Poeschl-Teller potential hole
 */

#ifndef QWWAD_SCHROEDINGER_SOLVER_POESCHL_TELLER_H
#define QWWAD_SCHROEDINGER_SOLVER_POESCHL_TELLER_H

#include "schroedinger-solver.h"

namespace QWWAD
{
/**
 * Schroedinger solver for a Poeschl-Teller potential hole
 */
class SchroedingerSolverPoeschlTeller : public SchroedingerSolver
{
public:
    SchroedingerSolverPoeschlTeller(const double alpha,
                                    const double lambda,
                                    const double length,
                                    const double mass,
                                    const size_t nz,
                                    const unsigned int nst_max = 0);

    std::string get_name() override {return "poeschl-teller-potential-hole";}
    size_t get_n_bound() const;
private:
    double _alpha;  ///< Width parameter [1/m]
    double _lambda; ///< Depth parameter
    double _length; ///< Length of potential profile [m]
    double _mass;   ///< Effective mass in well [kg]

    void calculate() override;
};
} // namespace
#endif
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :
