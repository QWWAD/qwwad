/**
 * \file   qwwad-schroedinger-donor-2D.h
 * \author Alex Valavanis <a.valavanis@leeds.ac.uk>
 * \brief  Declarations for Schroedinger solver for a donor in a 1D potential
 */

#ifndef QWWAD_SCHROEDINGER_DONOR_2D_H
#define QWWAD_SCHROEDINGER_DONOR_2D_H

#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

#include "qwwad-schroedinger-donor.h"

namespace Leeds {
/**
 * Schroedinger solver for a donor in a 1D potential with a 2D hydrogenic wavefunction
 */
class SchroedingerSolverDonor2D : public SchroedingerSolverDonor
{
public:
    SchroedingerSolverDonor2D(const double                 m,
                              const std::valarray<double> &V,
                              const std::valarray<double> &z,
                              const double                 eps,
                              const double                 r_d,
                              const double                 lambda,
                              const double                 dE,
                              const unsigned int           nst_max = 0);

    std::string get_name() {return "donor-2D";}

private:
    double I_1(const double z_dash) const;
    double I_2(const double z_dash) const;
    double I_3(const double z_dash) const;
    double I_4(const double z_dash) const;
};
} // namespace Leeds
#endif
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :
