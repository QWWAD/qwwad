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
    SchroedingerSolverDonor2D(const double        m,
                              const decltype(_V) &V,
                              const decltype(_z) &z,
                              const double        eps,
                              const double        r_d,
                              const double        lambda,
                              const double        dE);

    std::string get_name() override {return "donor-2D";}

private:
    void calculate_psi_from_chi() override {
        _solutions.clear();

        for (auto & st : _solutions_chi) {
            _solutions.push_back(st);
        }
    }

    [[nodiscard]] double I_1(const double z_dash) const override;
    [[nodiscard]] double I_2(const double z_dash) const override;
    [[nodiscard]] double I_3(const double z_dash) const override;
    [[nodiscard]] double I_4(const double z_dash) const override;
};
} // namespace QWWAD
#endif
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :
