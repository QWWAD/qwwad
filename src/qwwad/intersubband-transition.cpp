#include "intersubband-transition.h"
#include "constants.h"
#include "maths-helpers.h"

namespace QWWAD {
using namespace constants;

/**
 * \brief Initialise an intersubband transition
 *
 * \param[in] isb Initial subband
 * \param[in] fsb Final subband
 * \param[in] ki  Initial wave-vector samples [1/m]
 * \param[in] Wif Scattering rate at each wave-vector [1/s]
 */
IntersubbandTransition::IntersubbandTransition(const decltype(_isb) isb,
                                               const decltype(_fsb) fsb,
                                               const decltype(_ki)  ki,
                                               const decltype(_Wif) Wif) :
    _isb(isb),
    _fsb(fsb),
    _ki(ki),
    _Wif(Wif)
{
    const auto nki = _ki.size();

    // Tabulate the initial kinetic energy for each initial state in table
    _Eki.resize(nki);

    for(unsigned int iki = 0; iki < nki; ++iki)
        _Eki[iki] = isb.get_Ek_at_k(_ki[iki]);
}

/**
 * \brief Return the average scattering rate
 */
double IntersubbandTransition::get_average_rate() const
{
    const auto nki = _ki.size();
    const auto dki = _ki[1] - _ki[0];

    std::valarray<double> Wbar_integrand_ki(nki);

    for(unsigned int iki=0; iki<nki; ++iki)
    {
        const auto ki   = _ki[iki];
        const auto f_FD = _isb.get_occupation_at_k(ki);
        Wbar_integrand_ki[iki] = _Wif[iki]*ki*f_FD;
    } // End loop over ki

    const auto N = _isb.get_total_population();
    const auto Wif_avg = integral(Wbar_integrand_ki, dki)/(pi*N);
    return Wif_avg;
}
} // namespace
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :
