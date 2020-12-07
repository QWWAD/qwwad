/*=================================================================
       ppff     PseudoPotential Form Factors
  =================================================================

   This function returns the atomic formfactor (Fourier Transform 
   of the atomic potential) for a variety of atomic species.

   Paul Harrison, July 1998					 */

#ifndef PPFF_H
#define PPFF_H

#include <armadillo>
#include <vector>

typedef struct
{
 char	   *type;
 arma::vec  r;
}atom;

auto Vf(const double  A0,
          const double  m_per_au,
          double        q_sqr,
          const char   *type) -> double;

auto read_atoms(const char * filename) -> std::vector<atom>;

auto read_rlv(double A0) -> std::vector<arma::vec>;
#endif
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :
