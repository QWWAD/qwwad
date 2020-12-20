/**
 * \file   schroedinger-solver-donor.h
 * \author Alex Valavanis <a.valavanis@leeds.ac.uk>
 * \brief  Declarations for Schroedinger solver for a hydrogenic donor in a potential
 */

#ifndef QWWAD_SCHROEDINGER_SOLVER_DONOR_H
#define QWWAD_SCHROEDINGER_SOLVER_DONOR_H

#include "schroedinger-solver.h"

namespace QWWAD {
/**
 * \brief Schroedinger solver for the ground state of a hydrogenic donor in an external potential
 */
class SchroedingerSolverDonor : public SchroedingerSolver
{
public:
    SchroedingerSolverDonor(const double     m,
                            const arma::vec &V,
                            const arma::vec &z,
                            const double     eps,
                            const double     r_d,
                            const double     lambda,
                            const double     dE);

    auto get_name() -> std::string override = 0;

    auto get_solutions_chi(const bool convert_to_meV=false) -> std::vector<Eigenstate>;

    static auto chi_at_inf(double  E,
                             void   *params) -> double;

    auto shoot_wavefunction(const double  E,
                              arma::vec    &chi) const -> double;

    void set_lambda(const double lambda) {_lambda = lambda; refresh_solutions();}
    [[nodiscard]] auto get_lambda() const -> double {return _lambda;}
    [[nodiscard]] auto get_r_d   () const -> double {return _r_d;}

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
    std::vector<Eigenstate> _solutions_chi;

    auto calculate() -> std::vector<Eigenstate> override;
    virtual auto calculate_psi_from_chi() -> std::vector<Eigenstate> = 0;
    [[nodiscard]] virtual auto I_1(const double z_dash) const -> double = 0;
    [[nodiscard]] virtual auto I_2(const double z_dash) const -> double = 0;
    [[nodiscard]] virtual auto I_3(const double z_dash) const -> double = 0;
    [[nodiscard]] virtual auto I_4(const double z_dash) const -> double = 0;
};
} // namespace QWWAD
#endif
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :
