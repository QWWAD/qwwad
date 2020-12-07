/**
 * \file   fermi.h
 * \brief  Functions for finding the Fermi energy and population
 * \author Alex Valavanis <a.valavanis@leeds.ac.uk>
 */

#ifndef QWWAD_FERMI_H
#define QWWAD_FERMI_H

#include <vector>
#include "eigenstate.h"

namespace QWWAD
{
auto f_FD(const double E_F, const double Ek, const double Te) -> double;

auto f_FD_ionised(const double E_F, const double Ed, const double Te) -> double;

auto find_pop(const double Esb,
              const double E_F,
              const double md,
              const double Te,
              const double alpha=0,
              const double V=0) -> double;

auto find_fermi(const double Esb,
                const double m0,
                const double N,
                const double Te,
                const double alpha=0,
                const double V=0) -> double;

auto find_fermi_global(const arma::vec &Esb,
                       const double     m0,
                       const double     N,
                       const double     Te,
                       const double     alpha=0,
                       const double     V=0) -> double;

auto find_fermi_global(const std::vector<Eigenstate> &states,
                       const double                   m0,
                       const double                   N,
                       const double                   Te,
                       const double                   alpha=0,
                       const double                   V=0) -> double;
} // namespace
#endif
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :
