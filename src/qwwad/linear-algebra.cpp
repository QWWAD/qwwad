/**
 * \file   linear-algebra.cpp
 *
 * \brief  Linear algebra utility functions
 *
 * \author Jonathan Cooper <jdc-tas@gmail.com>
 * \author Alex Valavanis <a.valavanis@leeds.ac.uk>
 */

#if HAVE_CONFIG_H
# include "config.h"
#endif //HAVE_CONFIG_H

#include "linear-algebra.h"
#include "lapack-declarations.h"

#include <cstdlib>

#include "maths-helpers.h"
#include <gsl/gsl_math.h>
#include <gsl/gsl_sort.h>

#include <lapack.h>

namespace QWWAD
{
/**
 * \brief Find real solutions to eigenvalue problem from LAPACK
 *
 * \param[in]  A  Input matrix
 * \param[in]  VL Lower limit for eigenvalue search
 * \param[in]  VR Upper limit for eigenvalue search
 * \param[in]  N  Order of matrix
 * \param[in]  n_max Max number of eigenvalues to find
 *
 * \details    Creates standard inputs for Armadillo eigensolver before
 *             executing to return results.  If n_max=0, then all
 *             eigenvalues in the range [VL,VU] will be found.
 */
auto
eigen_general(arma::mat    &A,
              double        VL,
              double        VU,
              unsigned int  n_max) -> std::vector< EVP_solution<double> >
{
    const int N = sqrt(A.size());

    arma::cx_vec eigval;
    arma::cx_mat V_right; // Right eigenvalues

    auto complete = arma::eig_gen(eigval, V_right, A, "balance");

    if (!complete) {
        throw std::runtime_error("Could not solve "
                "eigenproblem. Check all input parameters!");
    }

    // Create buffer for output data (Should never have more than N 'real' solutions since
    // 2*N solutions correspond to psi*E and psi*E^2!)
    std::vector< EVP_solution<double> > solutions(N, EVP_solution<double>(N) );

    // Real and imaginary parts of the computed eigenvalues
    arma::vec WR = real(eigval);
    arma::vec WI = imag(eigval);

    unsigned int nst = 0; // Number of real solutions in the desired range

    // Loop through all possible solutions in solver output
    // and find the real solutions in the desired range
    for(int i=0; i<N; i++) {
        // Only allow real eigenvalues above the lower search limit
        if(WR[i]!=0 and WI[i]==0 and WR[i] > VL) {
            // If we specify a range of eigenvalues, filter the
            // solutions by that range, otherwise keep all of them
            if((n_max > 0) or (WR[i] < VU)){
                arma::vec const psi = real(V_right.col(i));
                solutions[nst] = EVP_solution<double>(WR[i], psi);

                nst++; // Register solution found

                if(nst > (unsigned int)N) {
                    throw std::runtime_error("More solutions found than nz!");
                }
            }
        }
    }

    // Shrink output array to match number of solutions
    solutions.resize(nst, EVP_solution<double>(N));

    // Create temporary storage for sorting the eigenvalues
    arma::vec E_tmp(nst);

    for(unsigned int ist=0; ist < nst; ist++) {
        E_tmp[ist] = solutions[ist].get_E();
    }

    // Indices of the eigenvalues after sorting
    arma::uvec sorted_E_indices = sort_index(E_tmp);

    // Copy the solutions into the desired output array
    std::vector<EVP_solution<double>> solutions_sorted;

    for(auto idx : sorted_E_indices)
    {
        solutions_sorted.push_back(solutions[idx]);

        // Stop if we reach the maximum permitted number of states
        if(n_max > 0 and solutions_sorted.size() == n_max) {
            break;
        }
    }

    return solutions_sorted;
}

/**
 * \brief Find solution to symmetric-definite banded eigenvalue problem A*x=lambda*B*x
 *
 * \param[in]  AB    Upper triangle of symmetric matrix A
 * \param[in]  BB    Upper triangle of symmetric matrix B
 * \param[in]  VL    Lowest value for eigenvalue search
 * \param[in]  VU    Highest value for eigenvalue search
 * \param[in]  n     Order of matrix
 * \param[in]  n_max Max number of eigenvalues to find
 *
 * \details    Creates standard inputs for dsbgvx func. before
 *             executing to return results.  If n_max=0, then all
 *             eigenvalues in the range [VL,VU] will be found.
 */
auto
eigen_banded(double       *AB,
             double       *BB,
             double        VL,
             double        VU,
             int           n,
             unsigned int  n_max) -> std::vector< EVP_solution<double> >
{
    // Workspace to normalise eigenproblem
    arma::vec Q(n*n);

    // LAPACK workspace
    arma::Col<int> ifail(n); // Failure bits for LAPACK
    arma::vec  W(n);    // Temporary storage for eigenvalues
    arma::vec  Z(n*n);  // Temp. storage for eigenvectors
    int                   M;       // Number of solutions found

    // Specify range of solutions by value, unless n_max is given
    char range = (n_max==0) ? 'V' : 'I';

    // Find error tolerance
    char   retval = 'S';   // Return value for LAPACK
    double abstol = 2.0 * dlamch_(&retval
#ifdef LAPACK_FORTRAN_STRLEN_END
            ,1
#endif
            );

    // LAPACK workspace
    arma::vec  work(7*(size_t)n);
    arma::Col<int> iwork(5*(size_t)n);

    // Run LAPACK function to solve eigenproblem
    char jobz  = 'V';   // Task descriptor for LAPACK
    char uplo  = 'U';   // Specifiy if lower or upper triangles are stored
    int  KA    = 1;     // Number of superdiagonals in A
    int  KB    = 1;     // Number of superdiagonals in B
    int  LD    = 2;     // Leading dimension of both AB and BB
    int  IL    = 1;     // Index of first solution to find
    int  IU    = n_max; // Index of last solution to find
    int  info  = 0;     // Output code from LAPACK

    dsbgvx_(&jobz, &range, &uplo, &n, &KA, &KB, AB, &LD, BB, &LD, &Q[0], &n, &VL,
            &VU, &IL, &IU, &abstol, &M, &W[0], &Z[0], &n, &work[0], &iwork(0), &ifail[0], &info
#ifdef LAPACK_FORTRAN_STRLEN_END
            ,1, 1, 1
#endif
            );

    if(info!=0) {
        throw std::runtime_error("Could not solve "
                "eigenproblem. Check all input parameters!");
    }

    // Extract solutions from LAPACK output
    std::vector< EVP_solution<double> > solutions(M, EVP_solution<double>(n) );

    for(int i = 0; i < M; i++){
        const std::vector<double> psi(Z.begin() + n*i,
                                      Z.begin() + n*(i+1));
        solutions[i] = EVP_solution<double>(W[i], psi);
    }

    return solutions;
}

/**
 * \brief Find solution to eigenvalue problem from LAPACK
 *
 * \param[in]  diag    Array holding all diagonal elements of matrix
 * \param[in]  subdiag Array holding all sub-diag. elements of matrix
 * \param[in]  VL      Lowest value for eigenvalue search
 * \param[in]  VU      Highest value for eigenvalue search
 * \param[in]  n_max   Max number of eigenvalues to find
 *
 * \details    Creates standard inputs for dstevx func. before
 *             executing to return results.  If n_max=0, then all
 *             eigenvalues in the range [VL,VU] will be found.
 */
auto
eigen_tridiag(arma::vec    &diag,
              arma::vec    &subdiag,
              double        VL,
              double        VU,
              unsigned int  n_max) -> std::vector< EVP_solution<double> >
{
    const int N    = diag.size();
    const int Nsub = subdiag.size();

    if (Nsub != N-1)
    {
        std::ostringstream oss;
        oss << "Size mismatch for tridiagonal elements: "
            << "(subdiagonal = " << Nsub << "; "
            << "diagonal = " << N << ")";

        throw std::runtime_error(oss.str());
    }

    arma::Col<int> ifail = arma::zeros<arma::Col<int>>(N); // Failure bits for LAPACK
    arma::vec W = arma::zeros(N);     // Temporary storage for eigenvalues
    arma::mat Z = arma::zeros(N,N);   // Temp. storage for eigenvectors
    int M = 0; // Number of solutions found

    // Specify range of solutions by value, unless n_max is given
    char range = (n_max==0) ? 'V' : 'I';

    // If we're checking by range by value, make sure that the upper and lower
    // bounds make sense
    if(n_max == 0 && gsl_fcmp(VL, VU, VL*1e-6) != -1)
    {
        std::ostringstream oss;
        oss << "Range of eigenvalue search is invalid. Lower limit: " << VL << " is greater than upper limit: " << VU;
        throw std::domain_error(oss.str());
    }

    int  info = 0; // Output code from LAPACK
    char jobz='V'; // Task descriptor for LAPACK
    int  IL=1;
    int  IU=n_max;
    arma::vec  work = arma::zeros(5*N); // LAPACK workspace
    arma::Col<int> iwork = arma::zeros<arma::Col<int>>(5*N);

    // Find error tolerance
    char retval='S'; // Return value for LAPACK
    double abstol = 2.0 * dlamch_(&retval
#ifdef LAPACK_FORTRAN_STRLEN_END
            ,1
#endif
            );

    // Run LAPACK function to solve eigenproblem
    dstevx_(&jobz,
            &range,
            &N,
            diag.memptr(),
            subdiag.memptr(),
            &VL, &VU,
            &IL, &IU,
            &abstol,
            &M,
            W.memptr(),
            Z.memptr(),
            &N,
            work.memptr(),
            iwork.memptr(),
            ifail.memptr(),
            &info
#ifdef LAPACK_FORTRAN_STRLEN_END
            , 1, 1
#endif
            );

    if(info!=0)
    {
        std::ostringstream oss;
        oss << "Could not solve eigenvalue problem. LAPACK error code: "
            << info;
        throw std::runtime_error(oss.str());
    }

    // Extract solutions from LAPACK output
    std::vector<EVP_solution<double>> solutions(M, EVP_solution<double>(N));

    for(int i = 0; i < M; i++){
        solutions[i] = EVP_solution<double>(W(i), Z.col(i));
    }

    return solutions;
}

/**
 * \brief Solves a matrix of the cyclic form, generated from the cyclic form of the Poisson solver
 *
 * \param[in]  A_sub  Array holding all sub-diagonal elements of matrix.
 * \param[in]  A_diag Array holding all diagonal elements of matrix
 * \param[in]  cyclic Value of the matrix in the bottom corner which is non-zero due to
 *                    cyclic boundaries.
 * \param[in]  b      Array holding the RHS of the equation to solve (A*X = b)
 *
 * \details    An optimised algorithm for solving the matrix resulting from applying cyclic
 *             boundaries to a discretised PDE problem. For more details of this optimised
 *             algorithm (probably) see Jonathan Cooper's thesis.
 */
auto solve_cyclic_matrix(arma::vec A_sub,
                              arma::vec A_diag,
                              double cyclic,
                              arma::vec b) -> arma::vec
{
    unsigned int ni = A_diag.size();
    arma::vec z(ni);
    arma::vec A_super(A_sub);
    //A_sub = A_super;
    //Forward sweep
    // Initial elements (don't need to set A_diag and b!)
    //A_diag[0] = A_diag[0];
    //b_new[0] = b[0];
    z[0] = 1;
    for(unsigned int i=1; i<ni-1; i++){
        A_diag[i] = A_diag[i] - (A_sub[i-1]*A_super[i-1])/A_diag[i-1];
        b[i] = b[i]-(A_sub[i-1]*b[i-1])/A_diag[i-1];

        // This can probably go in A_sub
        z[i] = -z[i-1]*A_super[i-1]/A_diag[i-1];
    }

    // Last F_dash element
    A_diag[ni-1] = A_diag[ni-1] - (A_sub[ni-2] + cyclic*z[ni-2])*A_super[ni-2]/A_diag[ni-2];

    // Last L_dash element
    double sum = 0.0;
    for(unsigned int i=0; i<ni-2; i++){
        sum += -cyclic*z[i]*b[i]/A_diag[i];
    }

    b[ni-1] = b[ni-1] + sum - (A_sub[ni-2] + cyclic*z[ni-2])*b[ni-2]/A_diag[ni-2];

    // Again no need to initialise first element of 
    // b[ni-1] = b[ni-1]
    for(int i=ni-2; i>-1; i--) {
        b[i] = b[i]-(A_super[i]*b[i+1])/A_diag[i+1];
    }

    return b/A_diag;
}

/**
 * \brief Perform matrix multiplication: y = Mx + c
 *
 * \details x and c are vectors, M is a tridiagonal matrix
 *
 * \param[in] M_sub   Subdiagonal of matrix M
 * \param[in] M_diag  Diagonal of matrix M
 * \param[in] M_super Superdiagonal of matrix M
 * \param[in] x       Vector x
 * \param[in] c       Vector c
 *
 * \return Vector y
 */
auto
multiply_vec_tridiag(arma::vec const &M_sub,
                     arma::vec const &M_diag,
                     arma::vec const &M_super,
                     arma::vec const &x,
                     arma::vec const &c) -> arma::vec
{
    int N = M_diag.size(); // Order of matrix

    double scale = 1;   // Don't scale the c or x vectors
    char   TRANS = 'N'; // Don't transpose the M matrix
    int    NRHS  = 1;   // Only solve for one column-vector

    // Temp copy of x that will be overwritten by LAPACK
    arma::vec x_tmp = x;

    // Perform matrix multiplication using LAPACK
    // result:-> result + RHS*Told
    dlagtm_(&TRANS,
            &N,
            &NRHS,
            &scale,
            M_sub.memptr(),
            M_diag.memptr(),
            M_super.memptr(),
            x_tmp.memptr(),
            &N,
            &scale,
            c.memptr(),
            &N);

    // The x vector has now been overwritten by the resulting
    // y vector, so we can return this
    arma::vec y = x_tmp;
    return y;
}

/**
 * \brief Solve a linear equation Ax = b
 *
 * \param[in] A_sub   The subdiagonal of matrix A
 * \param[in] A_diag  The diagonal of matrix A
 * \param[in] A_super The superdiagonal of matrix A
 * \param[in] b       The right-hand-side vector b
 *
 * \details where b is a vector and A is a tridiagonal matrix
 *
 * \return The vector x
 */
auto
solve_tridiag(arma::vec &A_sub,
              arma::vec &A_diag,
              arma::vec &A_super,
              arma::vec const &b) -> arma::vec
{
    int N    = A_diag.size();
    int NRHS = 1; // Solve for 1 RHS vector only

    arma::vec x_tmp = b;

    int INFO=0;
    dgtsv_(&N,
           &NRHS,
           A_sub.memptr(),
           A_diag.memptr(),
           A_super.memptr(),
           x_tmp.memptr(),
           &N,
           &INFO);

    if(INFO != 0)
    {
        std::ostringstream oss;
        oss << "Cannot solve matrix equation. (LAPACK error code: " << INFO << ")";
        throw std::runtime_error(oss.str());
    }

    return x_tmp;
}

/**
 * \brief Solve a linear equation Ax = b using the L*D*L**T factorisation of A
 *
 * \param[in] D The diagonal of the factorisation matrix D
 * \param[in] L The subdiagonal of the factorisation matrix L
 * \param[in] b The right-hand-side vector b
 *
 * \details where b is a vector and A is a positive definite symmetrical tridiagonal matrix
 *
 * \return The vector x
 */
auto
solve_tridiag_LDL_T(arma::vec const &D,
                    arma::vec const &L,
                    arma::vec const &b) -> arma::vec
{
    int N    = D.size();
    int NRHS = 1; // Solve for 1 RHS vector only

    arma::vec x_tmp = b;

    int INFO=0;
    dpttrs_(&N,
            &NRHS,
            D.memptr(),
            L.memptr(),
            x_tmp.memptr(),
            &N,
            &INFO);

    if(INFO != 0)
    {
        std::ostringstream oss;
        oss << "Cannot solve matrix equation. (LAPACK error code: " << INFO << ")";
        throw std::runtime_error(oss.str());
    }

    return x_tmp;
}

/**
 * \brief L*D*L**T factorisation of a positive definite tridiagonal matrix, A
 *
 * \param[in]  A_diag Diagonal of the matrix A
 * \param[in]  A_sub  Subdiagonal of the matrix A
 * \param[out] D      Diagonal of the factor matrix D
 * \param[out] L      Subdiagonal of the factor matrix L
 */
void
factorise_tridiag_LDL_T(arma::vec const &A_diag,
                        arma::vec const &A_sub,
                        arma::vec       &D,
                        arma::vec       &L)
{
    int info = 0; // Return value for LAPACK
    const int N = A_diag.size(); // Order of the matrix

    // Temporary copies of the input vector.
    // LAPACK overwrites these
    D = A_diag;
    L = A_sub;

    dpttrf_(&N, &D[0], &L[0], &info);

    if(info != 0)
    {
        std::ostringstream oss;
        oss << "Cannot factorise matrix. (LAPACK error code: " << info << ")";
        throw std::runtime_error(oss.str());
    }
}
} // namespace
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :
