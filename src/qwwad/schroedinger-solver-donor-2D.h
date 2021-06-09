/**
 * \file   schroedinger-solver-donor-2D.h
 * \author Alex Valavanis <a.valavanis@leeds.ac.uk>
 * \brief  Declarations for Schroedinger solver for a donor in a 1D potential
 */

#ifndef QWWAD_SCHROEDINGER_SOLVER_DONOR_2D_H
#define QWWAD_SCHROEDINGER_SOLVER_DONOR_2D_H

#include "schroedinger-solver-donor.h"

namespace QWWAD {
/**
 * Schroedinger solver for a donor in a 1D potential with a 2D hydrogenic wavefunction
 */
class SchroedingerSolverDonor2D : public SchroedingerSolverDonor
{
public:
    SchroedingerSolverDonor2D(double           m,
                              const arma::vec &V,
                              const arma::vec &z,
                              double           eps,
                              double           r_d,
                              double           lambda,
                              double           dE);

    auto get_name() -> std::string override {return "donor-2D";}

private:
    auto calculate_psi_from_chi() -> std::vector<Eigenstate> override;

    [[nodiscard]] auto I_1(double z_dash) const -> double override;
    [[nodiscard]] auto I_2(double z_dash) const -> double override;
    [[nodiscard]] auto I_3(double z_dash) const -> double override;
    [[nodiscard]] auto I_4(double z_dash) const -> double override;
};
} // namespace QWWAD
#endif
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :
