/**
 * Shared functions for large-basis pseudo-potential calculations
 */

#ifndef QWWAD_PPLB_FUNCTIONS_H
#define QWWAD_PPLB_FUNCTIONS_H

#include <complex>
#include <vector>

#include <armadillo>

#include "ppff.h"

void
write_ank(arma::cx_mat &ank,
          int           ik,
          int           N,
          int           n_min,
          int           n_max);

std::complex<double> V(double                   A0,
                       double                   m_per_au,
                       std::vector<atom> const &atoms,
                       arma::vec const         &q);
#endif
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :
