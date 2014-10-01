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

#include <stdexcept>
#include <iostream>

namespace Leeds
{
/**
 * Create a Poisson solver
 *
 * \param[in] eps Permittivity at each point
 * \param[in] dx  Spatial step [m]
 * \param[in] bt  Poisson boundary condition type
 */
Poisson::Poisson(const std::valarray<double> &eps,
                 const double                 dx,
                 PoissonBoundaryType          bt) :
    _eps(eps),
    _eps_minus(eps), // Set the half-index permittivities
    _eps_plus(eps),  // to a default for now
    _dx(dx),
    _L(_eps.size()*_dx),
    diag(std::valarray<double>(_eps.size())),
    sub_diag(std::valarray<double>(_eps.size()-1)),
    corner_point(0.0),
    boundary_type(bt)
{
    compute_half_index_permittivity();

    switch(boundary_type)
    {
        case DIRICHLET:
            factorise_dirichlet();
            break;
        case MIXED:
            factorise_mixed();
            break;
        case ZERO_FIELD:
            factorise_zerofield();
            break;
   }
}

/**
 * \brief Find the permittivity at half-index points
 */
void Poisson::compute_half_index_permittivity()
{
    const size_t ni = _eps.size();

    if (_eps_minus.size() != ni && _eps_plus.size() != ni)
        throw std::runtime_error("Permittivity array not initialized");

    for(unsigned int i=0; i < ni; i++)
    {
        // Assume that start and end values of permittivity are identical
        // (i.e., a periodic structure). This might not be ideal for structures
        // e.g., contained within a cladding material.
        if(i==0)
            _eps_minus[i] = (_eps[i]   + _eps[ni-1]) / 2.0;
        else
            _eps_minus[i] = (_eps[i]   + _eps[i-1] ) / 2.0;

        if(i==ni-1)
            _eps_plus[i]  = (_eps[0]   + _eps[i]   ) / 2.0;
        else
            _eps_plus[i]  = (_eps[i+1] + _eps[i]   ) / 2.0;
    }
}

void Poisson::factorise_dirichlet()
{
    // Fill matrix to solve
    const size_t ni = _eps.size();

    for(unsigned int i=0; i < ni; i++)
    {
        // Diagonal elements b_i [QWWAD4, 3.80]
        diag[i] = (_eps_plus[i] + _eps_minus[i]) / (_dx * _dx);

        // Sub-diagonal elements a_(i+1], c_i [QWWAD4, 3.80]
        if(i>0)
            sub_diag[i-1] = -_eps_minus[i] / (_dx * _dx);
    }

    // Factorise matrix
    int info = 0;
#if HAVE_LAPACKE
    info = LAPACKE_dpttrf(ni, &diag[0], &sub_diag[0]);
#else
    const int N = ni;
    dpttrf_(&N, &diag[0], &sub_diag[0], &info);
#endif

    if(info != 0)
    {
        std::ostringstream oss;
        oss << "Cannot factorise Poisson equation. (Lapack error code: " << info << ")";
        throw std::runtime_error(oss.str());
    }
} // End Poisson::factorise_hardwall()

void Poisson::factorise_mixed()
{
    // Fill matrix to solve
    unsigned int ni = _eps.size();
    for(unsigned int i=0; i < ni; i++)
    {
        // Diagonal elements
        if(i<ni-1)
            diag[i] = (1/(_dx*_dx))*(_eps_plus[i] + _eps_minus[i]);
        else
        {
            diag[i] = (1/(_dx*_dx))*_eps_minus[i];
            corner_point = (1/(_dx*_dx))*_eps_plus[i];
        }

        // Sub-diagonal elements
        if(i>0)
            sub_diag[i-1] = -(1/(_dx*_dx))*_eps_minus[i];
    }
}

void Poisson::factorise_zerofield()
{
    // Fill matrix to solve
    unsigned int ni = _eps.size();
    for(unsigned int i=0; i < ni; i++)
    {
        // Diagonal elements
        if(i==0)
            diag[i] = (1/(_dx*_dx))*_eps_plus[i];
        else if(i==ni-1)
            diag[i] = (1/(_dx*_dx))*_eps_minus[i];
        else
            diag[i] = (1/(_dx*_dx))*(_eps_plus[i] + _eps_minus[i]);

        // Sub-diagonal elements
        if(i>0)
            sub_diag[i-1] = -(1/(_dx*_dx))*_eps_plus[i];
    }
}

std::valarray<double> Poisson::solve(std::valarray<double> rho)
{
    std::valarray<double> phi(rho.size());
    const size_t n = _eps.size();

    switch(boundary_type)
    {
        case DIRICHLET:
            {
                int nrhs = 1;
                int info = 0;
#if HAVE_LAPACKE
                info = LAPACKE_dpttrs(LAPACK_COL_MAJOR, n, nrhs, &diag[0], &sub_diag[0], &rho[0], n);
#else
                const int _N = n;
                dpttrs_(&_N, &nrhs, &diag[0], &sub_diag[0], &rho[0], &n, &info);
#endif
                if(info != 0)
                {
                    std::ostringstream oss;
                    oss << "Cannot solve Poisson equation. (Lapack error code: " << info << ")";
                    throw std::runtime_error(oss.str());
                }

                phi = rho;
            }
            break;
        case MIXED:
        case ZERO_FIELD:
            phi = solve_cyclic_matrix(sub_diag, diag, corner_point, rho);
            
            if(phi.size() != static_cast<size_t>(n))
                    throw std::runtime_error("Cannot solve Poisson equation.");
            break;
    }

    return phi; 
}

std::valarray<double> Poisson::solve(std::valarray<double> phi, double V_drop)
{
    const size_t n = _eps.size();

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
                const int _N = n;
                dpttrs_(&_N, &nrhs, &diag[0], &sub_diag[0], &phi[0], &n, &info);
#endif
                if(info != 0)
                {
                    std::ostringstream oss;
                    oss << "Cannot solve Poisson equation. (Lapack error code: " << info << ")";
                    throw std::runtime_error(oss.str());
                }
            }
            break;
        case MIXED:
            throw std::runtime_error("Cannot apply bias directly when solving the Poisson equation with "
                                     "mixed bounradries. Instead solve cyclic problem without bias, then solve Laplace "
                                     "equation and sum the result.");
        default:
            throw std::runtime_error("Unrecognised boundary type for Poisson solver.");
    }

    return phi; 
}

std::valarray<double> Poisson::solve_laplace(double V_drop)
{
    const size_t n = _eps.size();
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
                const int _N = n;
                dpttrs_(&_N, &nrhs, &diag[0], &sub_diag[0], &phi[0], &n, &info);
#endif
                if(info != 0)
                {
                    std::ostringstream oss;
                    oss << "Cannot solve Poisson equation. (Lapack error code: " << info << ")";
                    throw std::runtime_error(oss.str());
                }
            }
            break;
        case MIXED:
            throw std::runtime_error("Cannot solve the Laplace equation with mixed boundary conditions.");
        default:
            throw std::runtime_error("Unrecognised boundary type for Poisson solver.");
    }

    return phi; 
}

} // End namespace Leeds
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :
