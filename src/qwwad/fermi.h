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
auto f_FD(double E_F,
          double Ek,
          double Te) -> double;

auto f_FD_ionised(double E_F,
                  double Ed,
                  double Te) -> double;

auto find_pop(double Esb,
              double E_F,
              double m0,
              double Te,
              double alpha=0,
              double V=0) -> double;

auto find_fermi(double Esb,
                double m0,
                double N,
                double Te,
                double alpha=0,
                double V=0) -> double;

auto find_fermi_global(const arma::vec &Esb,
                       double     m0,
                       double     N,
                       double     Te,
                       double     alpha=0,
                       double     V=0) -> double;

auto find_fermi_global(const std::vector<Eigenstate> &states,
                       double                   m0,
                       double                   N,
                       double                   Te,
                       double                   alpha=0,
                       double                   V=0) -> double;
} // namespace
#endif
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :
