#include <complex>
#include <utility>

#include "scattering-calculator-LO.h"
#include "constants.h"
#include "maths-helpers.h"

namespace QWWAD {
using namespace constants;

/**
 * \brief Initialise an LO-phonon scattering calculation for a 2D system
 *
 * \param[in] subbands    The energy subbands in the system
 * \param[in] A0          Lattice constant for the crystal [m]
 * \param[in] Ephonon     Phonon energy [J]
 * \param[in] epss        Static dielectric constant [F/m]
 * \param[in] epsinf      High-frequency dielectric constant [F/m]
 * \param[in] m           Effective mass [kg]
 * \param[in] Te          Electron temperature [K]
 * \param[in] Tl          Lattice temperature [K]
 * \param[in] is_emission True if this is an emission process
 */
ScatteringCalculatorLO::ScatteringCalculatorLO(decltype(_subbands)    subbands,
                                               decltype(_A0)          A0,
                                               decltype(_Ephonon)     Ephonon,
                                               decltype(_epss)        epss,
                                               decltype(_epsinf)      epsinf,
                                               decltype(_m)           m,
                                               decltype(_Te)          Te,
                                               decltype(_Tl)          Tl,
                                               decltype(_is_emission) is_emission) :
    _subbands(std::move(subbands)),
    _A0(A0),
    _Ephonon(Ephonon),
    _epss(epss),
    _epsinf(epsinf),
    _m(m),
    _Te(Te),
    _Tl(Tl),
    _is_emission(is_emission),
    _enable_screening(true),
    _enable_blocking(true),
    _nki(101),
    _omega_0(_Ephonon/hBar),
    _N0(1.0/(exp(_Ephonon/(kB*_Tl))-1.0)),
    _prefactor(pi*e*e*_omega_0/_epss*(_epss/_epsinf-1)*
                (_N0 + (_is_emission?1:0))*
                2.0 * _m/(hBar*hBar)*2/(8*pi*pi*pi))
{
    set_phonon_samples(1001);
    calculate_screening_length();
}

/**
 * \brief Find the minimum initial kinetic energy that would allow scattering
 */
auto ScatteringCalculatorLO::get_Eki_min(const unsigned int i,
                                           const unsigned int f) const -> double
{
    const auto isb = _subbands[i];
    const auto fsb = _subbands[f];

    // Subband minima
    const auto Ei = isb.get_E_min();
    const auto Ef = fsb.get_E_min();

    // Subband separation
    const auto Eif = Ei - Ef;

    double Eki_min = 0.0; // Minimum initial kinetic energy

    if(_is_emission && Eif < _Ephonon)
    {
        Eki_min = _Ephonon - Eif;
    }
    else if(!_is_emission && -Eif > _Ephonon)
    {
        Eki_min = -Eif - _Ephonon;
    }

    return Eki_min;
}

/**
 * \brief Find the minimum initial wave-vector that would allow scattering
 */
auto ScatteringCalculatorLO::get_ki_min(const unsigned int i,
                                          const unsigned int f) const -> double
{
    const auto Eki_min = get_Eki_min(i,f);
    const auto isb = _subbands[i];
    const auto ki_min = isb.get_k_at_Ek(Eki_min);
    return ki_min;
}

/**
 * \brief Find a sensible cut-off value for the initial wave vector
 */
auto ScatteringCalculatorLO::get_ki_cutoff(const unsigned int i,
                                             const unsigned int f) const -> double
{
    const auto isb = _subbands[i];

    const auto Eki_min = get_Eki_min(i,f);

    auto Eki_max = Eki_min + 5.0 * kB * _Te;

    const auto Ei_F = isb.get_Ef(); // Fermi energy [J]
    const auto Ei   = isb.get_E_min(); // Subband edge [J]

    // If the Fermi energy is above the subband minimum, then add that on for good measure!
    if(Ei < Ei_F)
        Eki_max += Ei_F;

    return isb.get_k_at_Ek(Eki_max);
}

/**
 * \brief Sets the number of samples of the phonon wave vector to use
 *
 * \param[in] nKz Number of samples
 *
 * \details If the number is different from the currently-used value,
 *          the array of samples is recalculated accordingly.
 *          The map of form-factors will also be cleared and will be
 *          automatically regenerated the next time a scattering rate
 *          is needed.
 */
void ScatteringCalculatorLO::set_phonon_samples(const size_t nKz)
{
    if(nKz != _Kz.size())
    {
        _Kz.resize(nKz);
        _dKz = 2.0/(_A0*nKz);

        for(unsigned int iKz = 0; iKz < nKz; ++iKz)
            _Kz[iKz] = iKz * _dKz;

        ff_table.clear();
    }
}

/**
 * \brief Find the total scattering rate at a given initial wave-vector
 *
 * \param[in] isb The initial subband index
 * \param[in] fsb The final subband index
 * \param[in] ki  The initial wave vector
 *
 * \details The total scattering rate is computed for all possible
 *          transitions to any permitted state in the final subband
 */
auto ScatteringCalculatorLO::get_rate_ki(const unsigned int i,
                                           const unsigned int f,
                                           const double       ki) -> double
{
    const auto ki_min = get_ki_min(i,f);

    double Wif_ki = 0.0;

    if(ki >= ki_min)
    {
        const auto nKz = _Kz.size();
        arma::vec Wif_integrand_dKz(nKz); // Integrand for scattering rate

        const auto isb = _subbands[i];
        const auto fsb = _subbands[f];
        const auto Ei  = isb.get_E_min();
        const auto Ef  = fsb.get_E_min();

        auto Delta = Ef - Ei;

        if(_is_emission)
            Delta += _Ephonon;
        else
            Delta -= _Ephonon;

        auto idx = std::make_pair(i,f);

        if(ff_table.count(idx) == 0)
            make_ff_table(i,f);

        auto Gifsqr = ff_table[idx];

        // Integral over phonon wavevector Kz
        for(unsigned int iKz=0; iKz < nKz; ++iKz)
        {
            auto Kz_2 = _Kz[iKz] * _Kz[iKz];

            // Apply screening if wanted
            if(_enable_screening && iKz != 0)
                Kz_2 *= (1.0 + 2*_lambda_s_sq/Kz_2 + _lambda_s_sq*_lambda_s_sq/(Kz_2*Kz_2));

            const auto Kz_4 = Kz_2 * Kz_2;

            Wif_integrand_dKz[iKz] = Gifsqr[iKz] /
                                     sqrt(Kz_4 + 2.0*Kz_2 * (2.0*ki*ki - 2.0*_m*Delta/(hBar*hBar))+
                                          4.0*_m*_m*Delta*Delta/(hBar*hBar*hBar*hBar));
        } // end integral over Kz

        Wif_ki = _prefactor*pi*integral(Wif_integrand_dKz,_dKz);

        if(_enable_blocking)
        {
            // Initial and final kinetic energy
            const auto Eki = isb.get_Ek_at_k(ki);
            const auto Ekf = Eki - Delta;

            if(Ekf >= 0)
            {
                const auto kf   = fsb.get_k_at_Ek(Ekf);
                const auto f_FD = fsb.get_occupation_at_k(kf);
                Wif_ki *= (1.0 - f_FD);
            }
        }
    }

    return Wif_ki;
}

/**
 * \brief Returns the entire scattering table for an intersubband transition
 *
 * \param[in] i Initial subband index
 * \param[in] f Final subband index
 */
auto ScatteringCalculatorLO::get_transition(const unsigned int i,
                                                              const unsigned int f) -> IntersubbandTransition
{
    // Get the minimum and cut-off initial wave-vectors for the transition
    const auto kimin  = get_ki_min(i, f);
    const auto kimax  = get_ki_cutoff(i, f);
    const auto dki    = (kimax - kimin)/((_nki-1)); // Step length for integration [1/m]

    arma::vec ki(_nki);  // Initial wave vectors [1/m]
    arma::vec Wif(_nki); // Scattering rate at each wave-vector [1/s]

    for (unsigned int iki = 0; iki < _nki; ++iki)
    {
        ki[iki]  = kimin + dki * iki;
        Wif[iki] = get_rate_ki(i, f, ki[iki]);
    }

    const auto isb = _subbands[i];
    const auto fsb = _subbands[f];

    IntersubbandTransition tx(isb, fsb, ki, Wif);

    return tx;
}

/**
 * \brief Compute the squared screening length [QWWAD 3, 10.157]
 *
 * \details If screening is disabled, then the value = 0.
 */
void ScatteringCalculatorLO::calculate_screening_length()
{
    _lambda_s_sq = 0.0;

    if(_enable_screening)
    {
        // Sum over all subbands
        for(const auto jsb : _subbands)
        {
            const auto Ej   = jsb.get_E_min();
            const auto f_FD = jsb.get_occupation_at_E_total(Ej);
            _lambda_s_sq += sqrt(2.0*_m*Ej) * _m * f_FD;
        }

        _lambda_s_sq *= e*e/(pi*pi*hBar*hBar*hBar*_epss);
    }
}

/**
 * \brief Computes the formfactor at a range of phonon wave-vectors
 */
void ScatteringCalculatorLO::make_ff_table(const unsigned int i,
                                           const unsigned int f)
{
    const auto isb = _subbands[i];
    const auto fsb = _subbands[f];

    const auto _nKz = _Kz.size();
    const auto idx = std::make_pair(i,f);
    ff_table[idx].resize(_nKz);

    for(unsigned int iKz=0;iKz < _nKz;iKz++)
    {
        ff_table[idx][iKz] = Gsqr(_Kz[iKz], isb, fsb); // Squared form-factor
    }
}

/**
 * \brief calculates the overlap integral squared between the two states
 */
auto ScatteringCalculatorLO::Gsqr(const double   Kz,
                                    const Subband &isb,
                                    const Subband &fsb) -> double
{
    const auto     z = isb.z_array();
    const auto    dz = z[1] - z[0];
    const auto    nz = z.size();
    const auto psi_i = isb.psi_array();
    const auto psi_f = fsb.psi_array();

    std::complex<double> I(0,1); // Imaginary unit

    // Find form-factor integral
    arma::cx_vec G_integrand_dz(nz);

    for(unsigned int iz=0; iz<nz; ++iz)
        G_integrand_dz[iz] = exp(Kz*z[iz]*I) * psi_i[iz] * psi_f[iz];

    auto G = integral(G_integrand_dz, dz);

    return norm(G);
}

auto ScatteringCalculatorLO::get_ff_table(const unsigned int i,
                                               const unsigned int f) const -> arma::vec
{
    const auto idx = std::make_pair(i,f);
    const auto G = ff_table.at(idx);
    return G;
}
} // namespace
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :
