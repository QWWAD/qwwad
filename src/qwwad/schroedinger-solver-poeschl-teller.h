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
    SchroedingerSolverPoeschlTeller(double alpha,
                                    double lambda,
                                    double length,
                                    double mass,
                                    size_t nz,
                                    unsigned int nst_max = 0);

    auto get_name() -> std::string override {return "poeschl-teller-potential-hole";}
    [[nodiscard]] auto get_n_bound() const -> size_t;
private:
    double _alpha;  ///< Width parameter [1/m]
    double _lambda; ///< Depth parameter
    double _length; ///< Length of potential profile [m]
    double _mass;   ///< Effective mass in well [kg]

    auto calculate() -> std::vector<Eigenstate> override;
};
} // namespace
#endif
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :
