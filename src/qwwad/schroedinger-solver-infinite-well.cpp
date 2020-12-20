/**
 * \file schoedinger-solver-infinite-well.cpp
 *
 * \brief Class to calculate states in infinite well
 *
 * \author Paul Harrison  <p.harrison@shu.ac.uk>
 * \author Alex Valavanis <a.valavanis@leeds.ac.uk>
 */

#include "schroedinger-solver-infinite-well.h"
#include <stdexcept>
#include <gsl/gsl_math.h>
#include "constants.h"

namespace QWWAD
{
using namespace constants;
SchroedingerSolverInfWell::SchroedingerSolverInfWell(const double       me,
                                                     const double       L,
                                                     const size_t       nz,
                                                     const double       alpha,
                                                     const double       V0,
                                                     const unsigned int nst_max) :
    _me(me),
    L_(L),
    _alpha(alpha),
    V0_(V0),
    _nz(nz)
{
    set_nst_max(nst_max);
    make_z_array();
}

/**
 * \brief Fills the arrays of spatial points and potential
 *
 * \details Accounts for the psi=0 barrier regions if wanted
 */
void SchroedingerSolverInfWell::make_z_array()
{
    arma::vec z(_nz);
    arma::vec V(_nz);

    const double L_total = L_ + 2.0*Lb_;  // Well + 2 * barrier
    const double dz = L_total/(_nz-1); // Spatial step

    for(unsigned int iz = 0; iz < _nz; ++iz) {
        z[iz] = iz*dz;
        V[iz] = V0_;
    }

    set_z(z);
    set_V(V);
}

auto
SchroedingerSolverInfWell::calculate() -> std::vector<Eigenstate>
{
    const auto z = get_z();
    auto nz = z.size();
    auto nst_max = get_nst_max();

    std::vector<Eigenstate> solutions;

    // Loop over all required states
    for(unsigned int is=1; is <= nst_max; is++)
    {
        double E = 0;

        // Energy of state [J] (QWWAD3, 2.13)
        if(gsl_fcmp(_alpha, 0, 1e-6) == 1) {
            E = 1.0/(2.0*_alpha*_me*L_) * (
                    sqrt(_me*(
                            2.0*_alpha*gsl_pow_2(hBar*pi*is) + _me*L_*L_
                            )
                        )
                    -_me*L_
                    ) + V0_;
        } else {
            E = gsl_pow_2(pi*hBar*is/L_)/(2*_me) + V0_;
        }

        // Stop if we've exceeded the cut-off energy
        if(energy_above_range(E)) {
            break;
        }

        arma::vec psi = arma::zeros(nz); // Wavefunction amplitude at each point [m^{-0.5}]

        // Loop over spatial locations and find wavefunction
        // amplitude at each point (QWWAD3, 2.15)
        for(unsigned int iz=0; iz<nz; ++iz)
        {
            if(z[iz] > Lb_ && z[iz] < Lb_ + L_)
                psi[iz]=sqrt(2/L_)*sin(is*pi*(z[iz]-Lb_)/L_); // Wavefunction [m^{-0.5}]
        }

        // Don't store the solution if it's below the minimum energy
        if(!energy_below_range(E)) {
            solutions.emplace_back(E, z, psi);
        }
    }

    return solutions;
}
    
/**
 * \brief   Adds a pair of "barriers" to the outside of the structure
 *
 * \details The wavefunction is set to zero over this region
 */
void SchroedingerSolverInfWell::set_padding_width(const double Lb)
{
    if(Lb_ < 0)
        throw std::domain_error("Padding width must be positive");

    Lb_ = Lb;

    make_z_array(); // Regenerate spatial points
    calculate();    // Recalculate solutions
}
} // namespace
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :
