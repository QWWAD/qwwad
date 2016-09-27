/**
 * \file   qclsim_poisson_solver.cpp
 *
 * \brief  Poisson solver class
 *
 * \author Jonathan Cooper <jdc-tas@gmail.com>
 * \author Alex Valavanis  <a.valavanis@leeds.ac.uk>
 */

#include <iostream>
#include "qclsim_poisson_solver.h"

#include "qwwad/lapack-declarations.h"
#include "qwwad/linear-algebra.h"

#include <stdexcept>
#include <iostream>

namespace QWWAD
{
/**
 * Create a Poisson solver
 *
 * \param[in] eps Permittivity at each point
 * \param[in] dx  Spatial step [m]
 * \param[in] bt  Poisson boundary condition type
 */
Poisson::Poisson(const decltype(_eps) &eps,
                 const double          dx,
                 PoissonBoundaryType   bt) :
    _eps(eps),
    _eps_minus(eps), // Set the half-index permittivities
    _eps_plus(eps),  // to a default for now
    _dx(dx), // Size of cells in mesh
    _L(_eps.size() * _dx), // Samples are at CENTRE of each cell so total length of structure is nx dx
    _diag(arma::vec(_eps.size())),
    _sub_diag(arma::vec(_eps.size()-1)),
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
        _diag[i] = (_eps_plus[i] + _eps_minus[i]) / (_dx * _dx);

        // Sub-diagonal elements a_(i+1), c_i [QWWAD4, 3.80]
        if(i>0)
            _sub_diag[i-1] = -_eps_minus[i] / (_dx * _dx);
    }

    // Factorise matrix
    int info = 0;
    const int N = ni;
    dpttrf_(&N, &_diag[0], &_sub_diag[0], &info);

    if(info != 0)
    {
        std::ostringstream oss;
        oss << "Cannot factorise Poisson equation. (Lapack error code: " << info << ")";
        throw std::runtime_error(oss.str());
    }
}

void Poisson::factorise_mixed()
{
    // Fill matrix to solve
    unsigned int ni = _eps.size();
    for(unsigned int i=0; i < ni; i++)
    {
        // Diagonal elements
        if(i<ni-1)
            _diag[i] = (1/(_dx*_dx))*(_eps_plus[i] + _eps_minus[i]);
        else
        {
            _diag[i] = (1/(_dx*_dx))*_eps_minus[i];
            corner_point = (1/(_dx*_dx))*_eps_plus[i];
        }

        // Sub-diagonal elements
        if(i>0)
            _sub_diag[i-1] = -(1/(_dx*_dx))*_eps_minus[i];
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
            _diag[i] = (1/(_dx*_dx))*_eps_plus[i];
        else if(i==ni-1)
            _diag[i] = (1/(_dx*_dx))*_eps_minus[i];
        else
            _diag[i] = (1/(_dx*_dx))*(_eps_plus[i] + _eps_minus[i]);

        // Sub-diagonal elements
        if(i>0)
            _sub_diag[i-1] = -(1/(_dx*_dx))*_eps_plus[i];
    }
}

/**
 * \brief Solves the Poisson equation for a given charge-density with no potential drop
 *
 * \param[in] rho    The charge density profile [C m^{-3}]
 *
 * \return The potential profile [J]
 */
arma::vec Poisson::solve(const arma::vec &rho) const
{
    const auto n = _eps.size();

    if (rho.size() != n)
        throw std::runtime_error("Permittivity and charge density arrays have different sizes");

    auto rhs = rho; // Set right-hand-side to the charge-density

    // Create temporary copies of the diagonal and subdiagonal arrays since this
    // function promises not to change any member variables
    auto diag_tmp     = _diag;
    auto sub_diag_tmp = _sub_diag;

    auto phi = rhs; // Array in which to output the potential [J]

    switch(boundary_type)
    {
        case DIRICHLET:
            {
                int nrhs = 1;
                int info = 0;
                const int _N = n;
                dpttrs_(&_N, &nrhs, &diag_tmp[0], &sub_diag_tmp[0], &rhs[0], &_N, &info);
                if(info != 0)
                {
                    std::ostringstream oss;
                    oss << "Cannot solve Poisson equation. (Lapack error code: " << info << ")";
                    throw std::runtime_error(oss.str());
                }

                phi = rhs;
            }
            break;
        case MIXED:
        case ZERO_FIELD:
            phi = solve_cyclic_matrix(_sub_diag,
                                      _diag, corner_point, rho);
            break;
    }

    return phi; 
}

/**
 * \brief Solves the Poisson equation for a given charge-density and potential drop
 *
 * \param[in] rho    The charge density profile [C m^{-3}]
 * \param[in] V_drop The total potential drop across the structure [J]
 *
 * \return The potential profile [J]
 */
arma::vec Poisson::solve(const arma::vec &rho,
                         const double     V_drop) const
{
    const auto n = _eps.size();

    if (rho.size() != n)
        throw std::runtime_error("Permittivity and charge density arrays have different sizes");

    auto rhs = rho; // Set right-hand-side to the charge-density

    // Create temporary copies of the diagonal and subdiagonal arrays since this
    // function promises not to change any member variables
    auto diag_tmp     = _diag;
    auto sub_diag_tmp = _sub_diag;

    switch(boundary_type)
    {
        case DIRICHLET:
            {
                int nrhs = 1; // number of right-hand-side columns
                int info = 0; // return code for LAPACK

                // We want to fix the potential just BEFORE the structure to 0
                //   i.e., phi[-1] = 0
                // If the voltage drop across the structure is V_drop, then the LAST
                // point in the system takes this value
                //   i.e., phi[n-1] = V_drop = F * length
                // so the first point AFTER the structure has the potential
                //   phi[n] = F * (length + dz) = V_drop + F dz = V_drop (nz + 1) / nz

                const auto next_potential = V_drop * n / (n+1);

                // The boundary condition is then set according to QWWAD4, 3.110.
                rhs[n-1] += _diag[n-1] * next_potential;
                const int _N = n;
                dpttrs_(&_N, &nrhs, &diag_tmp[0], &sub_diag_tmp[0], &rhs[0], &_N, &info);
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
                                     "mixed boundaries. Instead solve cyclic problem without bias, then solve Laplace "
                                     "equation and sum the result.");
        default:
            throw std::runtime_error("Unrecognised boundary type for Poisson solver.");
    }

    // LAPACK outputs the potential in the right-hand-side vector on exit
    const auto phi = rhs;

    return phi;
}

/**
 * \brief Solve the Laplace equation (i.e., the Poisson equation with no charge)
 *
 * \param[in] V_drop the total potential drop across the system [J]
 *
 * \return The potential profile [J]
 */
arma::vec Poisson::solve_laplace(const double V_drop) const
{
    const size_t n = _eps.size();

    // Create an empty charge profile and solve the Poisson equation
    auto rho = arma::zeros<arma::vec>(n);
    const auto phi = solve(rho, V_drop);

    return phi; 
}

} // namespace
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :
