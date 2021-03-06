/**
 * \file   schroedinger-solver-donor-3D.h
 * \author Paul Harrison  <p.harrison@leeds.ac.uk>
 * \author Alex Valavanis <a.valavanis@leeds.ac.uk>
 * \brief  Schroedinger solver for a donor in a 1D potential with 3D hydrogenic wavefunction
 */

#ifndef QWWAD_SCHROEDINGER_SOLVER_DONOR_3D_H
#define QWWAD_SCHROEDINGER_SOLVER_DONOR_3D_H

#include "schroedinger-solver-donor.h"
#include <iostream>
namespace QWWAD
{
/**
 * Schroedinger solver for a donor in a 1D potential with a 3D hydrogenic wavefunction
 */
class SchroedingerSolverDonor3D : public SchroedingerSolverDonor
{
public:
    SchroedingerSolverDonor3D(double           m,
                              const arma::vec &V,
                              const arma::vec &z,
                              double           eps,
                              double           r_d,
                              double           lambda,
                              double           dE);

    auto get_name() -> std::string override {return "donor-3D";}

private:
    auto calculate_psi_from_chi() -> std::vector<Eigenstate> override;

    [[nodiscard]] auto I_1(double z_dash) const -> double override;
    [[nodiscard]] auto I_2(double z_dash) const -> double override;
    [[nodiscard]] auto I_3(double z_dash) const -> double override;
    [[nodiscard]] auto I_4(double z_dash) const -> double override;
};
} // namespace
#endif
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :
