/**
 * \file   qclsim-fermi.h
 * \brief  Functions for finding the Fermi energy and population
 * \author Alex Valavanis <a.valavanis@leeds.ac.uk>
 */

#ifndef QCLSIM_FERMI_H
#define QCLSIM_FERMI_H

#include <valarray>

namespace Leeds {
double f_FD(const double E_F, const double Ek, const double Te);

double f_FD_ionised(const double E_F, const double Ed, const double Te);

double find_pop(const double Esb,
                const double E_F,
                const double md,
                const double Te,
                const double alpha=0,
                const double V=0);

double find_fermi(const double Esb,
                  const double m,
                  const double N,
                  const double Te,
                  const double alpha=0,
                  const double V=0);

double find_fermi_global(const double                 m,
                         const double                 N,
                         const double                 Te,
                         const std::valarray<double>& E);
} // namespace Leeds
#endif
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :
