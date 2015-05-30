#include "intersubband-transition.h"
#include "constants.h"

namespace QWWAD {
using namespace constants;

/**
 * \brief Initialise an intersubband transition
 *
 * \param[in] isb Initial subband
 * \param[in] fsb Final subband
 * \param[in] ki  Initial wave-vector samples [1/m]
 * \param[in] kf  Final wave-vector samples [1/m]
 * \param[in] Wif Scattering rate at each wave-vector [1/s]
 */
IntersubbandTransition::IntersubbandTransition(const decltype(_isb) isb,
                                               const decltype(_fsb) fsb,
                                               const decltype(_ki)  ki,
             //                                  const decltype(_kf)  kf,
                                               const decltype(_Wif) Wif) :
    _isb(isb),
    _fsb(fsb),
    _ki(ki),
//    _kf(kf),
    _Wif(Wif)
{
    const auto nki = _ki.size();

    // Tabulate the initial kinetic energy for each initial state in table
    _Eki.resize(nki);

    for(unsigned int iki = 0; iki < nki; ++iki)
        _Eki[iki] = isb.Ek(_ki[iki]);
}

/**
 * \brief Return the average scattering rate
 */
double IntersubbandTransition::get_average_rate(const double Te) const
{
    const auto nki = _ki.size();
    const auto dki = _ki[1] - _ki[0];

    std::valarray<double> Wbar_integrand_ki(nki);

    for(unsigned int iki=0; iki<nki; ++iki)
    {
        const double ki = _ki[iki];

        Wbar_integrand_ki[iki] = _Wif[iki]*ki*_isb.f_FD_k(ki, Te);
    } // End loop over ki

    const auto N = _isb.get_pop();
    const auto Wif_avg = integral(Wbar_integrand_ki, dki)/(pi*N);
    return Wif_avg;
}
} // namespace
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :
