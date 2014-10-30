/**
 * \file   qwwad-schroedinger-donor-2D.h
 * \author Alex Valavanis <a.valavanis@leeds.ac.uk>
 * \brief  Declarations for Schroedinger solver for a donor in a 1D potential
 */

#ifndef QWWAD_SCHROEDINGER_DONOR_2D_H
#define QWWAD_SCHROEDINGER_DONOR_2D_H

#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

#include "qwwad-schroedinger.h"

namespace Leeds {
/**
 * Schroedinger solver for a donor in a 1D potential with a 2D hydrogenic wavefunction
 */
class SchroedingerSolverDonor2D : public SchroedingerSolver
{
public:
    SchroedingerSolverDonor2D(const double                 m,
                              const std::valarray<double> &V,
                              const std::valarray<double> &z,
                              const double                 eps,
                              const double                 r_d,
                              const double                 lambda,
                              const double                 dE,
                              const unsigned int           nst_max = 0);

    std::string get_name() {return "donor-2D";}

    std::vector<State> get_solutions_chi(const bool convert_to_meV=false);

    static double psi_at_inf(double  E,
                             void   *params);

    double shoot_wavefunction(const double E,
                              std::valarray<double> &psi,
                              std::valarray<double> &chi) const;

private:
    double _me;     ///< Effective mass at band-edge [kg]
    double _eps;    ///< Permittivity [F/m]
    double _r_d;    ///< Position of donor, relative to lhs [m]
    double _lambda; ///< Donor Bohr radius [m]
    double _dE;     ///< Minimum energy separation between states [J]

    ///< Set of solutions to the Schroedinger equation excluding hydrogenic component
    std::vector<State> _solutions_chi;

    void calculate();
    double I_1() const;
    double I_2() const;
    double I_3() const;
    double I_4(const double z_dash) const;
};
} // namespace Leeds
#endif
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :
