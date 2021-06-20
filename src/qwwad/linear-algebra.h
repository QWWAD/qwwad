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
    T E_;              ///< Eigenvalue
    arma::Col<T> _psi; ///< Eigenvector

public:
    /**
     * Initialise an EVP solution by copying known eigenvalue/
     * eigenvector pair
     */
    EVP_solution(const T        E,
                 decltype(_psi) psi)
        : E_(E), _psi(std::move(psi))
    {}

    EVP_solution(const size_t n) : 
        E_(0),
        _psi(arma::Col<T>(n))
    {}

    /**
     * Return the value of element i of an eigenvector
     */
    [[nodiscard]] auto psi(const unsigned int i) const -> T
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
    [[nodiscard]] inline auto psi_array() const -> decltype(_psi)
    {
        return _psi;
    }

    [[nodiscard]] auto size() const -> size_t
    {
        return _psi.size();
    }

    /** Return the eigenvalue */
    [[nodiscard]] auto get_E() const -> T
    {
        return E_;
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
    [[nodiscard]] auto psi_squared(const unsigned int i) const -> T
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
    [[nodiscard]] auto psi_squared() const -> decltype(_psi)
    {
        return square(_psi);
    }

    /**
     * Find the largest square-magnitude of any element in a set of eigenvectors
     */
    static auto psi_squared_max(const std::vector<EVP_solution<T>> &EVP) -> T
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

auto eigen_general(arma::mat    &A,
                   double        VL,
                   double        VU,
                   unsigned int  n_max=0) -> std::vector<EVP_solution<double>>;

auto eigen_banded(double       *AB,
                  double       *BB,
                  double        VL,
                  double        VU,
                  int           n,
                  unsigned int  n_max = 0) -> std::vector<EVP_solution<double>>;

auto eigen_tridiag(arma::vec    &diag,
                   arma::vec    &subdiag,
                   double        VL,
                   double        VU,
                   unsigned int  n_max = 0) -> std::vector<EVP_solution<double>>;

auto multiply_vec_tridiag(arma::vec const &M_sub,
                          arma::vec const &M_diag,
                          arma::vec const &M_super,
                          arma::vec const &x,
                          arma::vec const &c) -> arma::vec;

auto solve_tridiag(arma::vec &A_sub,
                   arma::vec &A_diag,
                   arma::vec &A_super,
                   arma::vec const &b) -> arma::vec;

auto solve_tridiag_LDL_T(arma::vec const &D,
                         arma::vec const &L,
                         arma::vec const &b) -> arma::vec;

void
factorise_tridiag_LDL_T(arma::vec const &A_diag,
                        arma::vec const &A_subdiag,
                        arma::vec       &D,
                        arma::vec       &L);

auto solve_cyclic_matrix(arma::vec A_sub,
                         arma::vec A_diag,
                         double    cyclic,
                         arma::vec  b) -> arma::vec;

void matrixProduct(double* pB,
                   double* pA,
                   size_t  N);
} // namespace
#endif //QWWAD_LINEAR_ALGEBRA_H
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :
