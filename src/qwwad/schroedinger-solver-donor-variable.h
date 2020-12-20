/**
 * \file   schroedinger-solver-donor-variable.h
 * \author Paul Harrison  <p.harrison@leeds.ac.uk>
 * \author Alex Valavanis <a.valavanis@leeds.ac.uk>
 * \brief  Schroedinger solver for a donor in a 1D potential with variable-symmetry hydrogenic wavefunction
 */

#ifndef QWWAD_SCHROEDINGER_SOLVER_DONOR_VARIABLE_H
#define QWWAD_SCHROEDINGER_SOLVER_DONOR_VARIABLE_H

#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

#include "schroedinger-solver-donor.h"

namespace QWWAD
{
/**
 * Schroedinger solver for a donor in a 1D potential with a variable symmetry hydrogenic wavefunction
 */
class SchroedingerSolverDonorVariable : public SchroedingerSolverDonor
{
public:
    SchroedingerSolverDonorVariable(const double     m,
                                    const arma::vec &V,
                                    const arma::vec &z,
                                    const double     eps,
                                    const double     r_d,
                                    const double     lambda,
                                    const double     zeta,
                                    const double     dE);

    auto get_name() -> std::string override {return "donor-variable";}
    void   set_zeta       (const double zeta) {_zeta = zeta; refresh_solutions();}
    void   set_lambda_zeta(const double lambda, const double zeta) {_lambda = lambda; _zeta = zeta; refresh_solutions();}
    [[nodiscard]] auto get_zeta() const -> double {return _zeta;}

private:
    auto calculate_psi_from_chi() -> std::vector<Eigenstate> override;

    [[nodiscard]] auto I_1(const double z_dash) const -> double override;
    [[nodiscard]] auto I_2(const double z_dash) const -> double override;
    [[nodiscard]] auto I_3(const double z_dash) const -> double override;
    [[nodiscard]] auto I_4(const double z_dash) const -> double override;

    double _zeta; ///< Symmetry parameter
};
} // namespace
#endif
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :
