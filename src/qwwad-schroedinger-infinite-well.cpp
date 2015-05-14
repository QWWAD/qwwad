/**
 * \file qwwad-schoedinger-infinite-well.h Class to calculate states in infinite well
 *
 * \author Paul Harrison  <p.harrison@shu.ac.uk>
 * \author Alex Valavanis <a.valavanis@leeds.ac.uk>
 */

#include "qwwad-schroedinger-infinite-well.h"
#include <gsl/gsl_math.h>
#include "qwwad/constants.h"

using namespace QWWAD;
using namespace constants;

namespace Leeds {
SchroedingerSolverInfWell::SchroedingerSolverInfWell(const double       me,
                                                     const double       L,
                                                     const size_t       nz,
                                                     const double       alpha,
                                                     const double       V,
                                                     const unsigned int nst_max) :
    SchroedingerSolver(std::valarray<double>(nz),
                       std::valarray<double>(nz),
                       nst_max),
    _me(me),
    _L(L),
    _alpha(alpha),
    _V(V),
    _nz(nz),
    _Lb(0)
{
    make_z_array();
}

/**
 * \brief Creates an array of spatial points
 *
 * \details Accounts for the psi=0 barrier regions if wanted
 */
void SchroedingerSolverInfWell::make_z_array()
{
    const double L_total = _L + 2*_Lb;  // Well + 2 * barrier
    const double dz = L_total/(_nz-1); // Spatial step

    for(unsigned int iz = 0; iz < _nz; ++iz)
        _z[iz] = iz*dz;
}

void SchroedingerSolverInfWell::calculate()
{
    double nz = _z.size();

    // Loop over all required states
    for(unsigned int is=1; is<=_nst_max; is++)
    {
        double E = 0;

        // Energy of state [J] (QWWAD3, 2.13)
        if(gsl_fcmp(_alpha, 0, 1e-6) == 1)
            E = 1.0/(2.0*_alpha*_me*_L) * (
                    sqrt(_me*(
                            2.0*_alpha*gsl_pow_2(hBar*pi*is) + _me*_L*_L
                            )
                        )
                    -_me*_L
                    ) + _V;
        else
            E = gsl_pow_2(pi*hBar*is/_L)/(2*_me) + _V;

        // Stop if we've exceeded the cut-off energy
        if(_E_cutoff_set && gsl_fcmp(E, _E_cutoff, e*1e-12) == 1)
            break;

        std::valarray<double> psi(nz); // Wavefunction amplitude at each point [m^{-0.5}]

        // Loop over spatial locations and find wavefunction
        // amplitude at each point (QWWAD3, 2.15)
        for(unsigned int iz=0; iz<nz; ++iz)
        {
            if(_z[iz] > _Lb && _z[iz] < _Lb + _L)
                psi[iz]=sqrt(2/_L)*sin(is*pi*(_z[iz]-_Lb)/_L); // Wavefunction [m^{-0.5}]
        }

        _solutions.push_back(State(E, psi));
    }
}
    
/**
 * \brief   Adds a pair of "barriers" to the outside of the structure
 *
 * \details The wavefunction is set to zero over this region
 */
void SchroedingerSolverInfWell::set_padding_width(const double Lb)
{
    if(_Lb < 0)
        throw std::domain_error("Padding width must be positive");

    _Lb = Lb;

    make_z_array(); // Regenerate spatial points
    calculate();    // Recalculate solutions
}
} // namespace Leeds
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :
