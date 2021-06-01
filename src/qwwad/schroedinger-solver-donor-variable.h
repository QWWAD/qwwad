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
    SchroedingerSolverDonorVariable(double           m,
                                    const arma::vec &V,
                                    const arma::vec &z,
                                    double           eps,
                                    double           r_d,
                                    double           lambda,
                                    double           zeta,
                                    double           dE);

    auto get_name() -> std::string override {return "donor-variable";}
    void   set_zeta       (const double zeta) {_zeta = zeta; refresh_solutions();}
    void   set_lambda_zeta(const double lambda, const double zeta) {_lambda = lambda; _zeta = zeta; refresh_solutions();}
    [[nodiscard]] auto get_zeta() const -> double {return _zeta;}

private:
    auto calculate_psi_from_chi() -> std::vector<Eigenstate> override;

    [[nodiscard]] auto I_1(double z_dash) const -> double override;
    [[nodiscard]] auto I_2(double z_dash) const -> double override;
    [[nodiscard]] auto I_3(double z_dash) const -> double override;
    [[nodiscard]] auto I_4(double z_dash) const -> double override;

    double _zeta; ///< Symmetry parameter
};
} // namespace
#endif
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :
