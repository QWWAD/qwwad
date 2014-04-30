/**
 * \file qwwad-schoedinger-infinite-well.h Class to calculate states in infinite well
 *
 * \author Paul Harrison  <p.harrison@shu.ac.uk>
 * \author Alex Valavanis <a.valavanis@leeds.ac.uk>
 */

#include "qwwad-schroedinger-infinite-well.h"
#include <gsl/gsl_math.h>
#include "qclsim-constants.h"

namespace Leeds {
using namespace constants;

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
    _V(V)
{
    const double dz = L/(nz-1); // Spatial step

    for(unsigned int iz = 0; iz < nz; ++iz)
        _z[iz] = iz*dz;
}

void SchroedingerSolverInfWell::calculate()
{
    double nz = _z.size();
    _solutions.resize(_nst_max, State(nz));

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

        std::valarray<double> psi(nz); // Wavefunction amplitude at each point [m^{-0.5}]

        // Loop over spatial locations and find wavefunction
        // amplitude at each point (QWWAD3, 2.15)
        for(unsigned int i=0;i<nz;i++)
            psi=sqrt(2/_L)*sin(is*pi*_z/_L); // Wavefunction [m^{-0.5}]

        _solutions[is-1] = State(E, psi);
    }
}
} // namespace Leeds
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :
