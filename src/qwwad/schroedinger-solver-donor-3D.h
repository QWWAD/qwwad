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
    SchroedingerSolverDonor3D(const double        m,
                              const decltype(_V) &V,
                              const decltype(_z) &z,
                              const double        eps,
                              const double        r_d,
                              const double        lambda,
                              const double        dE);

    std::string get_name() override {return "donor-3D";}

private:
    void calculate_psi_from_chi() override
    {
        _solutions.clear();

        for (auto ist : _solutions_chi)
        {
            const auto E   = ist.get_energy();
            const auto chi = ist.get_wavefunction_samples();
            const auto psi = exp(-abs(_z - _r_d)/_lambda) * chi;

            const auto psi_state = Eigenstate(E,_z,psi);

            _solutions.push_back(psi_state);
        }
    }

    double I_1(const double z_dash) const override;
    double I_2(const double z_dash) const override;
    double I_3(const double z_dash) const override;
    double I_4(const double z_dash) const override;
};
} // namespace
#endif
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :
