/**
 * \file   fermi.h
 * \brief  Functions for finding the Fermi energy and population
 * \author Alex Valavanis <a.valavanis@leeds.ac.uk>
 * \date   2013-04-24
 */

#ifndef FERMI_H
#define FERMI_H

#if HAVE_CONFIG_H
# include "config.h"
#endif

#include <valarray>

class MaterialSpecification;

namespace Leeds {
double f_FD(const double E_F, const double Ek, const double Te);

double f_FD_ionised(const double E_F, const double Ed, const double Te);

double find_pop(const double md, const double E_F, const double Te);

double find_fermi(const double Esb, const double m, const double N, const double Te);

#if 0
double find_fermi(const MaterialSpecification* mat, const double N, 
        const double Te);

double find_fermi_global(const MaterialSpecification* mat, 
        const double N, const double Te, const std::valarray<double>& E);
#endif
} // namespace Leeds
#endif
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :
