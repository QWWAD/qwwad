/**
 * \file   qclsim-linalg.cpp
 *
 * \brief  Linear algebra utility functions for QCLsim
 *
 * \author Jonathan Cooper <el06jdc@leeds.ac.uk>
 * \author Alex Valavanis <a.valavanis@leeds.ac.uk>
 */

#if HAVE_CONFIG_H
# include "config.h"
#endif //HAVE_CONFIG_H

#if HAVE_LAPACKE
# include <lapacke.h>
#endif

#include "qclsim-linalg.h"

#include <cstdlib>

#include "qclsim-maths.h"
#include <gsl/gsl_sort.h>

namespace Leeds
{

/**
 * Find real solutions to eigenvalue problem from LAPACK
 *
 * \param[in]  A  Input matrix
 * \param[in]  VL Lower limit for eigenvalue search
 * \param[in]  VR Upper limit for eigenvalue search
 * \param[in]  N  Order of matrix
 * \param[in]  n_max Max number of eigenvalues to find
 *
 * \details    Creates standard inputs for dgeev func. before
 *             executing to return results.  If n_max=0, then all
 *             eigenvalues in the range [VL,VU] will be found.
 */
std::vector< EVP_solution<double> >
eigen_general(double       A[],
              const double VL,
              const double VU,
              int          N,
              unsigned int n_max)
{
    // Real and imaginary parts of the computed eigenvalues
    std::valarray<double> WR(N);
    std::valarray<double> WI(N);

    // Computed left and right eigenvectors
    std::valarray<double> V_left(N*N);
    std::valarray<double> V_right(N*N);

    // Run LAPACK function to solve eigenproblem
#if HAVE_LAPACKE
    int info = LAPACKE_dgeev(LAPACK_COL_MAJOR,
            'N', // Don't compute left eigenvectors
            'V', // DO compute right eigenvectors
            N, A, N, &WR[0], &WI[0], &V_left[0], N, &V_right[0], N);

#else
    int  info  = 0;   // Output code from LAPACK
    char jobvl = 'N'; // Specify range of solutions by value
    char jobvr = 'V';
    int  lwork = 4*N;
    std::valarray<double> work(4*N); // LAPACK workspace

    dgeev_(&jobvl, &jobvr, &N, A, &N, &WR[0], &WI[0], &V_left[0], &N, &V_right[0], &N,
            &work[0], &lwork, &info);
#endif

    if(info!=0)
        throw std::runtime_error("Could not solve "
                "eigenproblem. Check all input parameters!");

    // Create buffer for output data (Should never have more than N 'real' solutions since
    // 2*N solutions correspond to psi*E and psi*E^2!)
    std::vector< EVP_solution<double> > solutions(N, EVP_solution<double>(N) );

    unsigned int nst = 0; // Number of real solutions in the desired range

    // Loop through all possible solutions in LAPACK output
    // and find the real solutions in the desired range
    for(int i=0; i<N; i++)
    {
        // Only allow real eigenvalues above the lower search limit
        if(WR[i]!=0 and WI[i]==0 and WR[i] > VL)
        {
            // If we specify a range of eigenvalues, filter the
            // solutions by that range, otherwise keep all of them
            if((n_max > 0) or (WR[i] < VU)){
                const std::valarray<double> psi = V_right[std::slice(i*N, N, 1)];
                solutions[nst] = EVP_solution<double>(WR[i], psi);

                nst++; // Register solution found

                if(nst > (unsigned int)N)
                    throw std::runtime_error("More solutions found than nz!");
            }
        }
    }

    // Shrink output array to match number of solutions
    solutions.resize(nst, EVP_solution<double>(N));

    // Create temporary storage for sorting the eigenvalues
    std::valarray<double> E_tmp(nst);

    for(unsigned int ist=0; ist < nst; ist++)
        E_tmp[ist] = solutions[ist].get_E();

    // Indices of the eigenvalues after sorting
    std::valarray<size_t> sorted_E_indices(nst);

    // Find the correct order of the eigenvalues
    gsl_sort_index(&sorted_E_indices[0], &E_tmp[0], 1, nst);

    // Copy the solutions into the desired output array
    std::vector<State> solutions_sorted;

    for(unsigned int ist=0; ist < nst; ist++)
    {
        solutions_sorted.push_back(solutions[sorted_E_indices[ist]]);

        // Stop if we reach the maximum permitted number of states
        if(n_max > 0 and solutions_sorted.size() == n_max)
            break;
    }

    return solutions_sorted;
}

/**
 * Find solution to symmetric-definite banded eigenvalue problem A*x=lambda*B*x
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
std::vector< EVP_solution<double> >
eigen_banded(double       AB[],
             double       BB[],
             double       VL,
             double       VU,
             int          n,
             unsigned int n_max)
{
    // Workspace to normalise eigenproblem
    std::valarray<double> Q(n*n);

    // LAPACK workspace
    std::valarray<int>    ifail(n);// Failure bits for LAPACK
    std::valarray<double> W(n);    // Temporary storage for eigenvalues
    std::valarray<double> Z(n*n);  // Temp. storage for eigenvectors
    int                   M;       // Number of solutions found

    // Specify range of solutions by value, unless n_max is given
    char range = (n_max==0) ? 'V' : 'I';

    // TODO: Once all supported platforms can provide LAPACK > 3.4, kill this conditional code!
#if HAVE_LAPACKE
    // Run LAPACK function to solve eigenproblem
    //
    //   A*x = (lambda)*B*x,
    // 
    // where A and B are symmetric banded matrices
    int info = LAPACKE_dsbgvx(LAPACK_COL_MAJOR, // Use column-major ordering
                              'V',     // Compute eigenvalues and eigenvectors
                              range,
                              'U',     // Store upper triangle of matrix
                              n,       // Order of the A and B matrices
                              1,       // Number of superdiagonals in A
                              1,       // Number of superdiagonals in B
                              AB,      // Upper triangle of A
                              2,       // Leading dimension of A
                              BB,      // Upper triangle of B
                              2,       // Leading dimension of B
                              &Q[0],   // Workspace for normalising eigenproblem
                              n,       // Leading dimension of Q
                              VL, VU,
                              1, n_max,   // Index range for eigenvalue search
                              2.0 * LAPACKE_dlamch('S'), // Error tolerance (2*machine precision)
                              &M,
                              &W[0],
                              &Z[0],
                              n,
                              &ifail[0]);
#else
    // Find error tolerance
    char   retval = 'S';   // Return value for LAPACK
    double abstol = 2.0 * dlamch_(&retval); // Error tolerance

    // LAPACK workspace
    std::valarray<double> work(7*(size_t)n);
    std::valarray<int>    iwork(5*(size_t)n);

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
            &VU, &IL, &IU, &abstol, &M, &W[0], &Z[0], &n, &work[0], &iwork[0], &ifail[0], &info);
#endif // HAVE_LAPACKE

    if(info!=0)
        throw std::runtime_error("Could not solve "
                "eigenproblem. Check all input parameters!");

    // Extract solutions from LAPACK output
    std::vector< EVP_solution<double> > solutions(M, EVP_solution<double>(n) );

    for(int i = 0; i < M; i++){
        const std::valarray<double> psi = Z[std::slice(n*i, n, 1)];
        solutions[i] = EVP_solution<double>(W[i], psi);
    }

    return solutions;
}

/**
 * Find solution to eigenvalue problem from LAPACK
 *
 * \param[in]  D     Array holding all diagonal elements of matrix
 * \param[in]  E     Array holding all sub-diag. elements of matrix
 * \param[in]  VL    Lowest value for eigenvalue search
 * \param[in]  VU    Highest value for eigenvalue search
 * \param[in]  n     Order of matrix
 * \param[in]  n_max Max number of eigenvalues to find
 *
 * \details    Creates standard inputs for dstevx func. before
 *             executing to return results.  If n_max=0, then all
 *             eigenvalues in the range [VL,VU] will be found.
 */
std::vector< EVP_solution<double> >
eigen_tridiag(double       D[],
              double       E[],
              double       VL,
              double       VU,
              int          n,
              unsigned int n_max)
{
    std::valarray<int>    ifail(n); // Failure bits for LAPACK
    std::valarray<double> W(n);     // Temporary storage for eigenvalues
    std::valarray<double> Z(n*n);   // Temp. storage for eigenvectors
    int M; // Number of solutions found

    // Specify range of solutions by value, unless n_max is given
    char range = (n_max==0) ? 'V' : 'I';

    // TODO: Once all supported platforms can provide LAPACK > 3.4, kill this conditional code!
#if HAVE_LAPACKE
    // Run LAPACK function to solve eigenproblem
    int info = LAPACKE_dstevx(LAPACK_COL_MAJOR,
            'V',   // Find eigenvectors and eigenvalues
            range,
            n,     // Order of matrix
            D, E,
            VL, VU,
            1, n_max,   // Index range for eigenvalue search
            2.0 * LAPACKE_dlamch('S'), // Error tolerance (2*machine_precision)
            &M, &W[0], &Z[0], n, &ifail[0]);
#else
    int  info = 0; // Output code from LAPACK
    char jobz='V'; // Task descriptor for LAPACK
    int  IL=1;
    int  IU=n_max;
    std::valarray<double> work(5*n); // LAPACK workspace
    std::valarray<int>    iwork(5*n);

    // Find error tolerance
    char retval='S'; // Return value for LAPACK
    double abstol = 2.0 * dlamch_(&retval); // Error tolerance

    // Run LAPACK function to solve eigenproblem
    dstevx_(&jobz, &range, &n, D, E, &VL, &VU, &IL, &IU, &abstol, &M, &W[0],
            &Z[0], &n, &work[0], &iwork[0], &ifail[0], &info);
#endif

    if(info!=0)
        throw std::runtime_error("Could not solve "
                "eigenproblem. Check all input parameters!");

    // Extract solutions from LAPACK output
    std::vector< EVP_solution<double> > solutions(M, EVP_solution<double>(n) );

    for(int i = 0; i < M; i++){
        const std::valarray<double> psi = Z[std::slice(n*i, n, 1)];
        solutions[i] = EVP_solution<double>(W[i], psi);		
    }

    return solutions;
}

/**
 * Solves a matrix of the cyclic form, generated from the cyclic form of the Poisson solver
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
std::valarray<double> solve_cyclic_matrix(std::valarray<double> A_sub,
                                          std::valarray<double> A_diag,
                                          double cyclic,
                                          std::valarray<double> b)
{
        unsigned int ni = A_diag.size();
	std::valarray<double> z(ni);
        std::valarray<double> A_super(A_sub);
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
	for(int i=ni-2; i>-1; i--)
		b[i] = b[i]-(A_super[i]*b[i+1])/A_diag[i+1];

        return b/A_diag;

}

/**
 * \brief Find the expectation position for a given state
 *
 * \param[in] i State
 * \param[in] z Spatial coordinates [m]
 *
 * \return Expectation position [m]
 */
double z_av(const State& i, const std::valarray<double>& z)
{
    const double dz = z[1] - z[0];
    const std::valarray<double> dz_av = pow(i.psi_array(), 2.0) * z;

    return trapz(dz_av, dz);
}

/** 
 * \brief Find dipole matrix element between a pair of eigenvectors
 *
 * \param[in] i Initial state
 * \param[in] j Final state
 * \param[in] z Spatial coordinates [m]
 *
 * \return Dipole matrix element [m]
 *
 * \todo FIXME: When |i> == |j>, this should just return the expectation
 *       position.  At the moment, we have no way of figuring this out
 *       and the "pivoting" code below generates the wrong value.  It's
 *       not too important however, because the matrix element for an
 *       intrasubband transition is never really used!
 */
double mij(const State& i, const State& j, const std::valarray<double>& z)
{
    const double dz = z[1] - z[0];

    /* Because we have a nonparabolic effective mass, the Schroedinger solutions
     * are NOT part of an orthonormal set. As such, we need to do something to
     * ensure spatial invariance when we calculate dipole matrix element.  We
     * therefore define a pivot point halfway between the expectation positions
     * of the electron in each state.
     */
    const double z0 = 0.5 * (z_av(i,z) + z_av(j,z));

    const std::valarray<double> dmij = i.psi_array() * (z - z0) * j.psi_array();

    return trapz(dmij, dz);
}

} // namespace Leeds
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :
