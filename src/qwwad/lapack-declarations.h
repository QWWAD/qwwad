/**
 * \file   lapack-declarations.h
 * \brief  External LAPACK function declarations
 * \author Alex Valavanis <a.valavanis@leeds.ac.uk>
 */

#ifndef QWWAD_LAPACK_DECLARATIONS_H
#define QWWAD_LAPACK_DECLARATIONS_H

#include <complex>

extern "C" {

/**
 * Tridiagonal matrix multiplication: \f$B := \alpha A X + \beta B\f$
 */
void dlagtm_(const char   *TRANS,
             int          *N,
             const int    *NRHS,
             double       *ALPHA,
             double const *DL,
             double const *D,
             double const *DU,
             double       *X,
             int          *LDX,
             double       *BETA,
             double const *B,
             int          *LDB);
} // extern
#endif //QWWAD_LAPACK_DECLARATIONS_H
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :
