#ifndef QWWAD_DONOR_ENERGY_MINIMISER_FAST_H
#define QWWAD_DONOR_ENERGY_MINIMISER_FAST_H

#include "donor-energy-minimiser.h"

namespace QWWAD
{
class DonorEnergyMinimiserFast : public DonorEnergyMinimiser
{
public:
    DonorEnergyMinimiserFast(std::shared_ptr<SchroedingerSolverDonor> &se,
                             const double                              lambda_start,
                             const double                              lambda_step,
                             const double                              lambda_stop) :
        DonorEnergyMinimiser(se, lambda_start, lambda_step, lambda_stop)
    {};

private:
    void minimise() override;
};
}
#endif //QWWAD_DONOR_ENERGY_MINIMISER_FAST_H
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :
