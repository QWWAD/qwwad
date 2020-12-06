/**
 * \file   qwwad-donor-energy-minimiser.h
 *
 * \author Paul Harrison  <p.harrison@shu.ac.uk>
 * \author Alex Valavanis <a.valavanis@leeds.ac.uk>
 *
 * \brief  Minimisation of donor state energy
 */

#ifndef QWWAD_DONOR_ENERGY_MINIMISER_H
#define QWWAD_DONOR_ENERGY_MINIMISER_H

#include <vector>
#include <gsl/gsl_vector.h>

namespace QWWAD
{
class SchroedingerSolverDonor;

/**
 * \brief A tool for minimising the energy of a donor state
 */
class DonorEnergyMinimiser
{
public:
    DonorEnergyMinimiser(SchroedingerSolverDonor *se,
                         const double             lambda_start,
                         const double             lambda_step,
                         const double             lambda_stop);

    virtual ~DonorEnergyMinimiser() = default;

    virtual void minimise() = 0;
    void set_zeta_params(const double zeta_start,
                         const double zeta_step,
                         const double zeta_stop)
    {
        _zeta_start = zeta_start;
        _zeta_step  = zeta_step;
        _zeta_stop  = zeta_stop;
    };

protected:
    SchroedingerSolverDonor *_se; ///< The solver to be minimised

    double _lambda_start;
    double _lambda_step;
    double _lambda_stop;

    double _zeta_start;
    double _zeta_step;
    double _zeta_stop;

    // Logs of the search history
    std::vector<double> _lambda_history;
    std::vector<double> _zeta_history;
    std::vector<double> _E_history;

    static double find_E_at_lambda(double lambda, void *params);
    static double find_E_at_lambda_zeta(const gsl_vector *lambda_zeta,
                                        void             *params);
public:
    [[nodiscard]] decltype(_lambda_history) get_lambda_history() const {return _lambda_history;}
    [[nodiscard]] decltype(_zeta_history)   get_zeta_history()   const {return _zeta_history;}
    [[nodiscard]] decltype(_E_history)      get_E_history()      const {return _E_history;}
};
} // namespace
#endif // QWWAD_DONOR_ENERGY_MINIMISER_H
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :
