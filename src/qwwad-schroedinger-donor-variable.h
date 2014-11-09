/**
 * \file   qwwad-schroedinger-donor-variable.h
 * \author Paul Harrison  <p.harrison@leeds.ac.uk>
 * \author Alex Valavanis <a.valavanis@leeds.ac.uk>
 * \brief  Schroedinger solver for a donor in a 1D potential with variable-symmetry hydrogenic wavefunction
 */

#ifndef QWWAD_SCHROEDINGER_DONOR_VARIABLE_H
#define QWWAD_SCHROEDINGER_DONOR_VARIABLE_H

#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

#include "qwwad-schroedinger-donor.h"

namespace Leeds {
/**
 * Schroedinger solver for a donor in a 1D potential with a variable symmetry hydrogenic wavefunction
 */
class SchroedingerSolverDonorVariable : public SchroedingerSolverDonor
{
public:
    SchroedingerSolverDonorVariable(const double                 m,
                                    const std::valarray<double> &V,
                                    const std::valarray<double> &z,
                                    const double                 eps,
                                    const double                 r_d,
                                    const double                 lambda,
                                    const double                 zeta,
                                    const double                 dE,
                                    const unsigned int           nst_max = 0);

    std::string get_name() {return "donor-variable";}

private:
    void calculate_psi_from_chi()
    {
        _solutions.clear();

        for (unsigned int ist = 0; ist < _solutions_chi.size(); ++ist)
        {
            const std::valarray<double> chi = _solutions_chi[0].psi_array();
            const double E = _solutions_chi[0].get_E();

            std::valarray<double> psi = chi*exp(-_zeta*abs(_z - _r_d)/_lambda);
            _solutions.push_back(State(E,psi));
        }
    }

    double I_1(const double z_dash) const;
    double I_2(const double z_dash) const;
    double I_3(const double z_dash) const;
    double I_4(const double z_dash) const;

    double _zeta; ///< Symmetry parameter
};
} // namespace Leeds
#endif
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :
