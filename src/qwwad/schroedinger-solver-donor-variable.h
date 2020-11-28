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
    SchroedingerSolverDonorVariable(const double        m,
                                    const decltype(_V) &V,
                                    const decltype(_z) &z,
                                    const double        eps,
                                    const double        r_d,
                                    const double        lambda,
                                    const double        zeta,
                                    const double        dE);

    std::string get_name() override {return "donor-variable";}
    void   set_zeta       (const double zeta) {_zeta = zeta; _solutions.clear(); calculate();}
    void   set_lambda_zeta(const double lambda, const double zeta) {_lambda = lambda; _zeta = zeta; _solutions.clear(); calculate();}
    double get_zeta() const {return _zeta;}

private:
    void calculate_psi_from_chi() override
    {
        _solutions.clear();

        for (unsigned int ist = 0; ist < _solutions_chi.size(); ++ist)
        {
            const auto chi = _solutions_chi[0].get_wavefunction_samples();
            const double E = _solutions_chi[0].get_energy();

            auto const psi = chi*exp(-_zeta*abs(_z - _r_d)/_lambda);
            _solutions.push_back(Eigenstate(E,_z,psi));
        }
    }

    double I_1(const double z_dash) const override;
    double I_2(const double z_dash) const override;
    double I_3(const double z_dash) const override;
    double I_4(const double z_dash) const override;

    double _zeta; ///< Symmetry parameter
};
} // namespace
#endif
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :
