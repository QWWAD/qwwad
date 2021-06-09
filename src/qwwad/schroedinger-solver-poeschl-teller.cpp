/**
 *  \file     schroedinger-solver-poeschl-teller.cpp
 *  \author   Alex Valavanis <a.valavanis@leeds.ac.uk>
 *  \brief    Schrodinger solver functions for Poeschl-Teller potential hole
 */

#include "schroedinger-solver-poeschl-teller.h"

#include <cmath>
#include <gsl/gsl_sf_hyperg.h>
#include "constants.h"

namespace QWWAD
{
using namespace constants;
/**
 * \brief Construct a Poeschl-Teller potential solver
 *
 * \param[in] alpha   Width parameter [1/m]
 * \param[in] lambda  Depth parameter
 * \param[in] length  Length of potential profile [m]
 * \param[in] mass    Effective mass [kg]
 * \param[in] nz      Number of spatial points
 * \param[in] nst_max Maximum number of states to find
 */
SchroedingerSolverPoeschlTeller::SchroedingerSolverPoeschlTeller(const double alpha,
                                                                 const double lambda,
                                                                 const double length,
                                                                 const double mass,
                                                                 const size_t nz,
                                                                 const unsigned int nst_max) :
    _alpha(alpha),
    _lambda(lambda),
    _length(length),
    _mass(mass)
{
    arma::vec V(nz);
    arma::vec z(nz);
    set_nst_max(nst_max);

    // Compute potential profile [QWWAD4, 3.59]
    const double dz = _length / (nz - 1); // Size of each spatial cell [m]

    for (unsigned int iz=0; iz<nz; iz++)
    {
        z[iz] = iz*dz - _length/2.0;
        V[iz] = -gsl_pow_2(hBar*_alpha/cosh(_alpha*z[iz]))*_lambda*(_lambda-1)/(2*mass);
    }

    set_z(z);
    set_V(V);
}

/**
 * \brief Find number of bound states in the well
 */
auto SchroedingerSolverPoeschlTeller::get_n_bound() const -> size_t
{
    return ceil(_lambda-1);
}

auto
SchroedingerSolverPoeschlTeller::calculate() -> std::vector<Eigenstate>
{
    std::vector<Eigenstate> solutions;
    const auto z = get_z();
    const size_t nst = get_n_bound();
    const size_t nz  = z.size();

    const arma::vec sinh_alpha_z = sinh(_alpha*z);
    const arma::vec _x = -arma::square(sinh_alpha_z);

    const auto nst_max = get_nst_max();

    for (unsigned int ist=0; (nst_max == 0 || ist < nst_max) && ist < nst; ++ist) {
        // Energy is found using [QWWAD4, 3.60]
        const double kappa = _alpha * (_lambda-1-ist);
        const double E = -gsl_pow_2(hBar * kappa) / (2.0*_mass);

        // Wavefunction is taken from "Practical Quantum Mechanics", Flugge (1970).
        // The solution is a hypergeometric function, whose arguments and scaling
        // factor depend on whether the state has odd or even parity.
        
        // Arguments for hypergeometric function
        double arg1;
        double arg2;
        double arg3;
        arma::vec fact; // Prefactor for hypergeometric function

        // Flugge, 39.24
        const double a = 0.5 * (_lambda - kappa/_alpha);
        const double b = 0.5 * (_lambda + kappa/_alpha);

        if(ist % 2 != 0) { // Odd-parity states
            // From Flugge, 39.10b
            arg1 = a+0.5;
            arg2 = b+0.5;
            arg3 = 1.5;
            fact = pow(cosh(_alpha*z),_lambda) % sinh_alpha_z;
        } else { // Even-parity states
            // From Flugge, 39.10a
            arg1 = a;
            arg2 = b;
            arg3 = 0.5;
            fact = pow(cosh(_alpha*z),_lambda);
        }

        // Wavefunction amplitude at each point [m^{-0.5}]
        arma::cx_vec psi;

        for(unsigned int iz = 0; iz < nz; ++iz) {
            if(std::abs(_x(iz)) < 1) {
                psi(iz) = fact(iz) *
                          gsl_sf_hyperg_2F1(arg1, arg2, arg3, _x(iz));
            } else if(gsl_fcmp(_x[iz]/(_x[iz]-1), 1, 0.0025) == -1) {
            // If the argument is too large, we need to apply a linear
            // transformation such that |_x| < 1
                psi[iz] = fact[iz] *
                          pow(1-_x[iz],-arg2) *
                          gsl_sf_hyperg_2F1(arg2, arg3-arg1, arg3, _x[iz]/(_x[iz]-1));
            } else {
            // In case we're *very* close to _x = 1, GSL can't cope, so we
            // need to simplify things further, and just pass a large number
            // as the argument.
            // This seems to be OK, but might need a little investigation
                psi[iz] = fact[iz] *
                          pow(1-_x[iz],-arg2) *
                          gsl_sf_hyperg_2F1(arg2, arg3-arg1, arg3, 0.99);
            }
        }

        solutions.emplace_back(E, z, psi);
    }

    return solutions;
}
} // namespace
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :
