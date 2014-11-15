/**
 * \file   qwwad-schroedinger-donor.h
 * \author Alex Valavanis <a.valavanis@leeds.ac.uk>
 * \brief  Declarations for Schroedinger solver for a hydrogenic donor in a potential
 */

#ifndef QWWAD_SCHROEDINGER_DONOR_H
#define QWWAD_SCHROEDINGER_DONOR_H

#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

#include "qwwad-schroedinger.h"

namespace Leeds {
/**
 * \brief Schroedinger solver for the ground state of a hydrogenic donor in an external potential
 */
class SchroedingerSolverDonor : public SchroedingerSolver
{
public:
    SchroedingerSolverDonor(const double                 m,
                            const std::valarray<double> &V,
                            const std::valarray<double> &z,
                            const double                 eps,
                            const double                 r_d,
                            const double                 lambda,
                            const double                 dE);

    virtual std::string get_name() = 0;

    std::vector<State> get_solutions_chi(const bool convert_to_meV=false);

    static double chi_at_inf(double  E,
                             void   *params);

    double shoot_wavefunction(const double E,
                              std::valarray<double> &chi) const;

    void   set_lambda(const double lambda) {_lambda = lambda; _solutions.clear(); calculate();}
    double get_lambda() const {return _lambda;}
    double get_r_d   () const {return _r_d;}

private:
    double _me;     ///< Effective mass at band-edge [kg]
    double _eps;    ///< Permittivity [F/m]

protected:
    double _r_d;    ///< Position of donor, relative to lhs [m]
    double _lambda; ///< Donor Bohr radius [m]

private:
    double _dE;     ///< Minimum energy separation between states [J]

protected:
    ///< Set of solutions to the Schroedinger equation excluding hydrogenic component
    std::vector<State> _solutions_chi;

    void calculate();
    virtual void   calculate_psi_from_chi() = 0;
    virtual double I_1(const double z_dash) const = 0;
    virtual double I_2(const double z_dash) const = 0;
    virtual double I_3(const double z_dash) const = 0;
    virtual double I_4(const double z_dash) const = 0;
};
} // namespace Leeds
#endif
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :
