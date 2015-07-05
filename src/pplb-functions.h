/**
 * Shared functions for large-basis pseudo-potential calculations
 */

#ifndef QWWAD_PPLB_FUNCTIONS_H
#define QWWAD_PPLB_FUNCTIONS_H
#include <complex>
#include "ppff.h"

void
write_ank(arma::cx_mat &ank,
          int           ik,
          int           N,
          int           n_min,
          int           n_max);

std::complex<double> V(double           A0,
                       double           m_per_au,
                       atom            *atoms,
                       size_t           n_atoms,
                       arma::vec const &q);
#endif
