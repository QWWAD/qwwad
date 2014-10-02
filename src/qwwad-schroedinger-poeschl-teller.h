/**
 * \file   qwwad-schroedinger-poeschl-teller.h
 * \author Alex Valavanis <a.valavanis@leeds.ac.uk>
 * \brief  Schroedinger solver for a Poeschl-Teller potential hole
 */

#ifndef QWWAD_SCHROEDINGER_POESCHL_TELLER_H
#define QWWAD_SCHROEDINGER_POESCHL_TELLER_H

#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

#include "qwwad-schroedinger.h"

namespace Leeds {
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

    std::string get_name() {return "poeschl-teller-potential-hole";}
    size_t get_n_bound() const;
private:
    double _alpha;  ///< Width parameter [1/m]
    double _lambda; ///< Depth parameter
    double _length; ///< Length of potential profile [m]
    double _mass;   ///< Effective mass in well [kg]

    void calculate();
};
} // namespace Leeds
#endif
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :
