/**
 * \file   linear-algebra.h
 * \brief  Linear algebra utilities
 * \author Alex Valavanis <a.valavanis@leeds.ac.uk>
 */

#ifndef QWWAD_LINEAR_ALGEBRA_H
#define QWWAD_LINEAR_ALGEBRA_H

#if HAVE_CONFIG_H
# include "config.h"
#endif //HAVE_CONFIG_H

#include <complex>
#include <valarray>
#include <vector>
#include <sstream>
#include <stdexcept>

// Declaration of external LAPACK functions
extern "C" {
// If we don't have the official LAPACK C bindings, we need to declare
// the necessary external functions
//
// TODO: Kill all this once all supported platforms have LAPACK >= 3.4
#if !HAVE_LAPACKE
/** Determine double-precision machine parameters */
double dlamch_(const char* returnValue);

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
 * Solve a tridiagonal inversion problem: \f$ x := A^{-1} x\f$
 */
void dgtsv_(const int* N, const int* NRHS, double DL[], double D[], 
            double DU[], double B[], int* LDB, int* INFO);

/**
 * Factorise real symmetric positive definate tridagonal matrix
 */
void dpttrf_(const int* N, double D[], double E[], int* INFO);

/**
 * Solve real symmetric positive definite tridagonal matrix
 */
void dpttrs_(const int* N, const int* NRHS, double D[], double E[],
             double B[], const int* LDB, int* INFO);

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
 * \brief Solve eigenvalue problem for complex matrix
 */
void zheev_(const char *JOBZ, const char *UPLO, const int* N, std::complex<double> ank[],
            const int *LDA, double E[], std::complex<double> WORK[], int *LWORK, double RWORK[], int *INFO);
#endif // !HAVE_LAPACKE
/**
 * Tridiagonal matrix multiplication: \f$B := \alpha A X + \beta B\f$
 */
void dlagtm_(const char* TRANS, int* N, const int* NRHS, double* ALPHA,
             double DL[], double D[], double DU[], double X[], int* LDX,
             double* BETA, double B[], int* LDB);
} // extern

namespace QWWAD
{
/**
 * \brief A 3D vector
 */
template <class T>
class vector3D {
public:
    T x; ///< The 'x' component
    T y; ///< The 'y' component
    T z; ///< The 'z' component

    vector3D() :
        x(0), y(0), z(0)
    {}

    vector3D(T x, T y, T z) :
        x(x),
        y(y),
        z(z)
    {}

    vector3D& operator+=(const vector3D &rhs)
    {
        x+=rhs.x;
        y+=rhs.y;
        z+=rhs.z;

        return *this;
    }

    vector3D& operator*=(const T &rhs)
    {
        x*=rhs;
        y*=rhs;
        z*=rhs;

        return *this;
    }
};

template <class T>
inline vector3D<T> operator+(vector3D<T> lhs, const vector3D<T> &rhs)
{
    lhs += rhs;
    return lhs;
}

template <class T>
inline vector3D<T> operator*(vector3D<T> lhs, const T &rhs)
{
    lhs *= rhs;
    return lhs;
}

/**
 * \brief A solution to an eigenvalue problem
 */
template <class T>
class EVP_solution {
protected:
    T _E; ///< Eigenvalue
    std::valarray<T> _psi; ///< Eigenvector

public:
    /**
     * Initialise an EVP solution by copying known eigenvalue/
     * eigenvector pair
     */
    EVP_solution(T E, std::valarray<T> psi)
        : _E(E), _psi(psi)
    {}

    EVP_solution(const size_t n) : 
        _E(0),
        _psi(0.0,n)
    {}

    /**
     * Return the value of element i of an eigenvector
     */
    T psi(const unsigned int i) const
    {
        if(i >= _psi.size())
        {
            std::ostringstream oss;
            oss << "Requested element " << i << " in eigenvector of size " << _psi.size() << ".";
            throw std::domain_error(oss.str());
        }

        return _psi[i];
    }

    /**
     * Return the entire eigenvector as an array
     */
    std::valarray<T> psi_array() const
    {
        return _psi;
    }

    size_t size() const
    {
        return _psi.size();
    }

    /** Return the eigenvalue */
    T get_E() const
    {
        return _E;
    }

    /** 
     * Scale the eigenvector by a constant factor
     *
     * \param[in] a A constant scaling factor
     */
    void scale(const T a)
    {
        _psi*=a;
    }

    /**
     * Get the squared value of the eigenvector at a given point
     *
     * \param[in] i The index of the point for which to find eigenvector squared
     */
    T psi_squared(const unsigned int i) const
    {
        if(i >= _psi.size())
        {
            std::ostringstream oss;
            oss << "Requested element " << i << " in an eigenvector of size " << _psi.size();
            throw std::domain_error(oss.str());
        }

        return psi(i)*psi(i);
    }

    /**
     * Get an array of squared eigenvector values
     */
    std::valarray<T> psi_squared() const
    {
        return _psi*_psi;
    }

    /**
     * Find the largest square-magnitude of any element in a set of eigenvectors
     */
    static T psi_squared_max(const std::vector<EVP_solution<T>> &EVP)
    {
        double PDmax = 0.0;

        // Loop through all EVP solutions
        for(auto st: EVP)
        {
            // If this solution has highest square-magnitude so far, store its value
            const auto PD = st.psi_squared();
            PDmax = GSL_MAX_DBL(PDmax, PD.max());
        }

        return PDmax;
    }
};

std::vector< EVP_solution<double> >
eigen_general(double       A[],
              const double VL,
              const double VU,
              int          N,
              unsigned int n_max=0);

std::vector< EVP_solution<double> >
eigen_banded(double       AB[],
             double       BB[],
             const double VL,
             const double VU,
             int          n,
             unsigned int n_max = 0);

std::vector< EVP_solution<double> >
eigen_tridiag(double       D[],
              double       E[],
              const double VL,
              const double VU,
              int          n,
              unsigned int n_max = 0);

std::valarray<double>
solve_cyclic_matrix(std::valarray<double> A_sub,
                    std::valarray<double> A_diag,
                    double                cyclic,
                    std::valarray<double> b);

void matrixProduct(double*      pB,
                   double*      pA,
                   const size_t N);
} // namespace
#endif //QWWAD_LINEAR_ALGEBRA_H
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :
