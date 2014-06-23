/**
 * \file   qclsim_poisson_solver.cpp
 *
 * \brief  Poisson solver class for QCLsim
 *
 * \author Jonathan Cooper <el06jdc@leeds.ac.uk>
 * \author Alex Valavanis  <a.valavanis@leeds.ac.uk>
 */

#include <iostream>
#include "qclsim_poisson_solver.h"

#if HAVE_LAPACKE
# include <lapacke.h>
#endif

#include "qclsim-linalg.h"

#include <error.h>
#include <stdexcept>
#include <iostream>

namespace Leeds
{
/**
 * Create a Poisson solver
 *
 * \param[in] eps Permittivity at each point
 * \param[in] dx  Spatial step
 * \param[in] bt  Poisson boundary condition type
 */
Poisson::Poisson(const std::valarray<double> &eps,
                 const double                 dx,
                 PoissonBoundaryType          bt) :
    n(eps.size()),
    diag(std::valarray<double>(n)),
    sub_diag(std::valarray<double>(n-1)),
    corner_point(0.0),
    boundary_type(bt)
{
    switch(boundary_type)
    {
        case DIRICHLET:
            factorise_dirichlet(eps, dx);
            break;
        case MIXED:
            factorise_mixed(eps, dx);
            break;
        case ZERO_FIELD:
            factorise_zerofield(eps, dx);
            break;
        default:
            throw std::runtime_error("Unknown boundary condition type used to instansiate Poisson solver.");
   }
}

void Poisson::factorise_dirichlet(const std::valarray<double>& eps, const double dx)
{
    // Half indexed eps points
    double eps_i_minus;
    double eps_i_plus;

    // Fill matrix to solve
    unsigned int ni = eps.size();
    for(unsigned int i=0; i < ni; i++)
    {
        // Set eps_plus/minus_half values
        if(i==0)
        {
            eps_i_minus = (eps[i]+eps[ni-1])/2;
            eps_i_plus = (eps[i+1]+eps[i])/2;
        }
        else if(i==ni-1)
        {
            eps_i_minus = (eps[i]+eps[i-1])/2;
            eps_i_plus = (eps[0]+eps[i])/2;
        }
        else
        {
            eps_i_minus = (eps[i]+eps[i-1])/2;
            eps_i_plus = (eps[i+1]+eps[i])/2;
        }

        // Diagonal elemnts
        diag[i] = (1/(dx*dx))*(eps_i_plus + eps_i_minus);

        // Sub-diagonal elements
        if(i>0)
            sub_diag[i-1] = -(1/(dx*dx))*eps_i_minus;
    }

    // Factorise matrix
    int info = 0;
#if HAVE_LAPACKE
    info = LAPACKE_dpttrf(n, &diag[0], &sub_diag[0]);
#else
    dpttrf_(&n, &diag[0], &sub_diag[0], &info);
#endif

    if(info != 0)
        error(EXIT_FAILURE, 0, "Cannot factorise Poisson equation. (Lapack error code: %i)", info);	
} // End Poisson::factorise_hardwall()

void Poisson::factorise_mixed(const std::valarray<double>& eps, const double dx)
{
    // Half indexed eps points
    double eps_i_minus;
    double eps_i_plus;

    // Fill matrix to solve
    unsigned int ni = eps.size();
    for(unsigned int i=0; i < ni; i++)
    {
        // Set eps_plus/minus_half values
        if(i==0)
        {
                eps_i_minus = (eps[i]+eps[ni-1])/2;
                eps_i_plus = (eps[i+1]+eps[i])/2;
        }
        else if(i==ni-1)
        {
                eps_i_minus = (eps[i]+eps[i-1])/2;
                eps_i_plus = (eps[0]+eps[i])/2;
        }
        else
        {
                eps_i_minus = (eps[i]+eps[i-1])/2;
                eps_i_plus = (eps[i+1]+eps[i])/2;
        }

        // Diagonal elements
        if(i<ni-1)
            diag[i] = (1/(dx*dx))*(eps_i_plus + eps_i_minus);
        else
        {
            diag[i] = (1/(dx*dx))*eps_i_minus;
            corner_point = (1/(dx*dx))*eps_i_plus;
        }

        // Sub-diagonal elements
        if(i>0)
            sub_diag[i-1] = -(1/(dx*dx))*eps_i_minus;
    }
}

void Poisson::factorise_zerofield(const std::valarray<double>& eps, const double dx)
{
    // Half indexed eps points
    double eps_i_minus;
    double eps_i_plus;

    // Fill matrix to solve
    unsigned int ni = eps.size();
    for(unsigned int i=0; i < ni; i++)
    {
        // Set eps_plus/minus_half values
        if(i==0)
        {
                eps_i_minus = (eps[i]+eps[ni-1])/2;
                eps_i_plus = (eps[i+1]+eps[i])/2;
        }
        else if(i==ni-1)
        {
                eps_i_minus = (eps[i]+eps[i-1])/2;
                eps_i_plus = (eps[0]+eps[i])/2;
        }
        else
        {
                eps_i_minus = (eps[i]+eps[i-1])/2;
                eps_i_plus = (eps[i+1]+eps[i])/2;
        }

        // Diagonal elements
        if(i==0)
            diag[i] = (1/(dx*dx))*eps_i_plus;
        else if(i==ni-1)
            diag[i] = (1/(dx*dx))*eps_i_minus;
        else
            diag[i] = (1/(dx*dx))*(eps_i_plus + eps_i_minus);

        // Sub-diagonal elements
        if(i>0)
            sub_diag[i-1] = -(1/(dx*dx))*eps_i_plus;
    }
}

std::valarray<double> Poisson::solve(std::valarray<double> rho)
{
    std::valarray<double> phi(rho.size());

    switch(boundary_type)
    {
        case DIRICHLET:
            {
                int nrhs = 1;
                int info = 0;
#if HAVE_LAPACKE
                info = LAPACKE_dpttrs(LAPACK_COL_MAJOR, n, nrhs, &diag[0], &sub_diag[0], &rho[0], n);
#else
                dpttrs_(&n, &nrhs, &diag[0], &sub_diag[0], &rho[0], &n, &info);
#endif
                if(info != 0)
                    error(EXIT_FAILURE, 0, "Cannot solve Poisson equation. (Lapack error code: %i)", info);

                phi = rho;
            }
            break;
        case MIXED:
        case ZERO_FIELD:
            phi = solve_cyclic_matrix(sub_diag, diag, corner_point, rho);
            
            if(phi.size() != static_cast<size_t>(n))
                    error(EXIT_FAILURE, 0, "Cannot solve Poisson equation.");	
            break;
    }

    return phi; 
}

std::valarray<double> Poisson::solve(std::valarray<double> phi, double V_drop)
{
    switch(boundary_type)
    {
        case DIRICHLET:
            {
                int nrhs = 1;
                int info = 0;

                // Add voltage drop
                phi[n-1] += V_drop*diag[n-1];
#if HAVE_LAPACKE
                info = LAPACKE_dpttrs(LAPACK_COL_MAJOR, n, nrhs, &diag[0], &sub_diag[0], &phi[0], n);
#else
                dpttrs_(&n, &nrhs, &diag[0], &sub_diag[0], &phi[0], &n, &info);
#endif
                if(info != 0)
                    error(EXIT_FAILURE, 0, "Cannot solve Poisson equation. (Lapack error code: %i)", info);
            }
            break;
        case MIXED:
            error(EXIT_FAILURE, 0, "Cannot apply bias directly when solving the Poisson equation with "
                    "mixed bounradries. Instead solve cyclic problem without bias, then solve Laplace "
                    "equation and sum the result.");
        default:
            error(EXIT_FAILURE, 0, "Unregonised boundary type for Poisson solver.");
    }

    return phi; 
}


std::valarray<double> Poisson::solve_laplace(double V_drop)
{
    std::valarray<double> phi(0.0, n);
    switch(boundary_type)
    {
        case DIRICHLET:
            {
                int nrhs = 1;
                int info = 0;

                phi[n-1] += V_drop*diag[n-1];
#if HAVE_LAPACKE
                info = LAPACKE_dpttrs(LAPACK_COL_MAJOR, n, nrhs, &diag[0], &sub_diag[0], &phi[0], n);
#else
                dpttrs_(&n, &nrhs, &diag[0], &sub_diag[0], &phi[0], &n, &info);
#endif
                if(info != 0)
                    error(EXIT_FAILURE, 0, "Cannot solve Poisson equation. (Lapack error code: %i)", info);
            }
            break;
        case MIXED:
            error(EXIT_FAILURE, 0, "Cannot solve the Lapalace equation with mixed boundary conditions.");
        default:
            error(EXIT_FAILURE, 0, "Unregonised boundary type for Poisson solver.");
    }

    return phi; 
}

} // End namespace Leeds
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :
