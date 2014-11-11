/**
 * \file   qwwad-donor-energy-minimiser.h
 *
 * \author Paul Harrison  <p.harrison@shu.ac.uk>
 * \author Alex Valavanis <a.valavanis@leeds.ac.uk>
 *
 * \brief  Minimisation of donor state energy
 */

#include <gsl/gsl_vector.h>

namespace Leeds {
class SchroedingerSolverDonor;

enum MinimisationMethod
{
    MINIMISE_FAST,
    MINIMISE_LINEAR
};

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

    void minimise(MinimisationMethod method=MINIMISE_FAST);

private:
    SchroedingerSolverDonor *_se; ///< The solver to be minimised

    double _lambda_start;
    double _lambda_step;
    double _lambda_stop;

    void find_E_min_fast();
    void find_E_min_linear();
    static double find_E_at_lambda(double lambda, void *params);
    static double find_E_at_lambda_zeta(const gsl_vector *lambda_zeta,
                                        void             *params);
};
} // namespace Leeds
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :
