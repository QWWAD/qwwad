/**
 * \file   lapack-declarations.h
 * \brief  External LAPACK function declarations
 * \author Alex Valavanis <a.valavanis@leeds.ac.uk>
 */

#ifndef QWWAD_LAPACK_DECLARATIONS_H
#define QWWAD_LAPACK_DECLARATIONS_H

#include <complex>

extern "C" {

/** Solve general eigenproblem */
void dgeev_(const char *JOBVL,
            const char *JOBVR,
            const int  *N,
            double      A[],
            const int  *LDA,
            double      WR[],
            double      WI[],
            double      VL[],
            const int  *LDVL,
            double      VR[],
            const int  *LDVR,
            double      WORK[],
            int        *LWORK,
            int        *INFO);

/**
 * Factorise general matrix
 */
void dgetrf_(const int* M, const int* N, double A[], const int* LDA,
             int IPIV[], int* INFO);

/**
 * Solve general matrix
 */
void dgetrs_(const char* TRANS, const int* N, const int* NRHS, double A[],
             const int* LDA, int IPIV[], double B[], const int* LDB, int* INFO);

/**
 * Solve a tridiagonal inversion problem: \f$ x := A^{-1} x\f$
 */
void dgtsv_(const int    *N,
            const int    *NRHS,
            const double *DL,
            const double *D, 
            const double *DU,
            double       *B,
            const int    *LDB,
            int          *INFO);

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

/** Determine double-precision machine parameters */
auto dlamch_(const char* returnValue) -> double;

/**
 * Factorise real symmetric positive definate tridagonal matrix
 */
void dpttrf_(const int* N, double D[], double E[], int* INFO);

/**
 * Solve real symmetric positive definite tridagonal matrix
 */
void dpttrs_(const  int   *N,
             const  int   *NRHS,
             double const *D,
             double const *E,
             double       *B,
             const int    *LDB,
             int          *INFO);


/** Solve symmetric-definite banded eigenproblem*/
void dsbgvx_(const char* JOBZ,
             const char* RANGE,
             const char* UPLO,
             const int* N, const int* KA, const int* KB, double AB[], const int* LDAB,
             double BB[], const int* LDBB, double* Q, const int* LDQ, const double* VL,
             const double* VU, const int* IL, const int* IU, const double* ABSTOL, int* M,
             double W[], double Z[], const int* LDZ, double WORK[], int IWORK[],
             int IFAIL[], int* INFO);

/** Solve real symmetric tridiagonal eigenproblem */
void dstevx_(const char* JOBZ, const char* RANGE, const int* N, double D[],
             double E[], const double* VL, const double* VU,
             const int* IL, const int* IU, const double* ABSTOL, int* M, 
             double W[], double Z[], const int* LDZ, double WORK[], int IWORK[], 
             int IFAIL[], int* INFO);

/**
 * \brief Solve eigenvalue problem for complex matrix
 */
void zheev_(const char *JOBZ, const char *UPLO, const int* N, std::complex<double> ank[],
            const int *LDA, double E[], std::complex<double> WORK[], int *LWORK, double RWORK[], int *INFO);

} // extern
#endif //QWWAD_LAPACK_DECLARATIONS_H
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :
