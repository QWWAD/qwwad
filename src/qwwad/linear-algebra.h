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
#include <utility>

#include <vector>
#include <sstream>
#include <stdexcept>

#include <armadillo>


namespace QWWAD
{
/**
 * \brief A solution to an eigenvalue problem
 */
template <class T>
class EVP_solution {
protected:
    T _E;              ///< Eigenvalue
    arma::Col<T> _psi; ///< Eigenvector

public:
    /**
     * Initialise an EVP solution by copying known eigenvalue/
     * eigenvector pair
     */
    EVP_solution(const T        E,
                 decltype(_psi) psi)
        : _E(E), _psi(std::move(psi))
    {}

    EVP_solution(const size_t n) : 
        _E(0),
        _psi(arma::Col<T>(n))
    {}

    /**
     * Return the value of element i of an eigenvector
     */
    [[nodiscard]] T psi(const unsigned int i) const
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
    [[nodiscard]] inline decltype(_psi) psi_array() const
    {
        return _psi;
    }

    [[nodiscard]] size_t size() const
    {
        return _psi.size();
    }

    /** Return the eigenvalue */
    [[nodiscard]] T get_E() const
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
    [[nodiscard]] T psi_squared(const unsigned int i) const
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
    [[nodiscard]] decltype(_psi) psi_squared() const
    {
        return square(_psi);
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
eigen_general(arma::mat    &A,
              const double VL,
              const double VU,
              unsigned int n_max=0);

std::vector< EVP_solution<double> >
eigen_banded(double       *AB,
             double       *BB,
             const double  VL,
             const double  VU,
             int           n,
             unsigned int  n_max = 0);

std::vector< EVP_solution<double> >
eigen_tridiag(arma::vec   &D,
              arma::vec   &E,
              const double VL,
              const double VU,
              unsigned int n_max = 0);

arma::vec
multiply_vec_tridiag(arma::vec const &M_sub,
                     arma::vec const &M_diag,
                     arma::vec const &M_super,
                     arma::vec const &x,
                     arma::vec const &c);

arma::vec
solve_tridiag(arma::vec const &A_sub,
              arma::vec const &A_diag,
              arma::vec const &A_super,
              arma::vec const &x);

arma::vec
solve_tridiag_LDL_T(arma::vec const &D,
                    arma::vec const &L,
                    arma::vec const &x);

void
factorise_tridiag_LDL_T(arma::vec const &A_diag,
                        arma::vec const &A_subdiag,
                        arma::vec       &D,
                        arma::vec       &L);

arma::vec
solve_cyclic_matrix(arma::vec A_sub,
                    arma::vec A_diag,
                    double    cyclic,
                    arma::vec  b);

void matrixProduct(double*      pB,
                   double*      pA,
                   const size_t N);
} // namespace
#endif //QWWAD_LINEAR_ALGEBRA_H
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :
