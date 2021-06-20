/**
 * \file   scattering-calculator-LO.h
 * \author Alex Valavanis <a.valavanis@leeds.ac.uk>
 * \brief  Calculator for scattering rates for electron-LO phonon interactions
 */

#ifndef QWWAD_SCATTERING_CALCULATOR_LO
#define QWWAD_SCATTERING_CALCULATOR_LO

#include <map>
#include <utility>
#include "subband.h"
#include "intersubband-transition.h"

namespace QWWAD {
/**
 * \brief A calculator for electron-phonon scattering rates
 */
class ScatteringCalculatorLO {
private:
    std::vector<Subband> _subbands; ///< The energy subbands in the system

    // Physical properties
    double _A0;      ///< Lattice constant [m]
    double _Ephonon; ///< Phonon energy [J]
    double _epss;    ///< Static dielectric constant [F/m]
    double _epsinf;  ///< High-frequency dielectric constant [F/m]
    double _m;       ///< Effective mass [kg]
    double _Te;      ///< Electron temperature [K]
    double _Tl;      ///< Lattice temperature [K]

    bool _is_emission;      ///< True if this is an emission process
    bool _enable_screening; ///< Allow screening
    bool _enable_blocking;  ///< Allow final-state blocking

    // Precision parameters
    size_t _nki;     ///< Number of initial wave-vector samples

    // Derived properties
    decltype(_A0)      _dKz;         ///< Step size in phonon wave vector [1/m]
    decltype(_Ephonon) _omega_0;     ///< Phonon angular frequency [rad/s]
    decltype(_Ephonon) _N0;          ///< Bose-Einstein factor
    decltype(_Ephonon) _prefactor;   ///< Pre-factor for rates
    decltype(_A0)      _lambda_s_sq; ///< Squared screening length [m^2]

    using map_key = std::pair<unsigned int, unsigned int>;
    arma::vec _Kz; ///< Wave vector samples [1/m]

    /**
     * \brief Table of form factors
     *
     * \details The key refers to the initial and final subband indices
     *          The map contains a table of \f$G_{if}^2(Kz)\f$
     */
    std::map<map_key, arma::vec> ff_table;

    void calculate_screening_length();

public:
    ScatteringCalculatorLO(decltype(_subbands)    subbands,
                           decltype(_A0)          A0,
                           decltype(_Ephonon)     Ephonon,
                           decltype(_epss)        epss,
                           decltype(_epsinf)      epsinf,
                           decltype(_m)           m,
                           decltype(_Te)          Te,
                           decltype(_Tl)          Tl,
                           decltype(_is_emission) is_emission);

   [[nodiscard]] auto get_Eki_min (unsigned int isb,
                                   unsigned int fsb) const -> double;

   [[nodiscard]] auto get_ki_min (unsigned int isb,
                                  unsigned int fsb) const -> double;

   [[nodiscard]] auto get_ki_cutoff(unsigned int isb,
                                    unsigned int fsb) const -> double;

   auto get_rate_ki(unsigned int isb,
                    unsigned int fsb,
                    double       ki) -> double;

   auto get_transition(unsigned int isb,
                       unsigned int fsb) -> IntersubbandTransition;

   [[nodiscard]] inline auto get_screening_length() const {return _lambda_s_sq;}

   inline void set_ki_samples(const decltype(_nki) nki) {_nki = nki;}
   void set_phonon_samples(size_t nKz);

   [[nodiscard]] inline auto get_dKz() const {return _dKz;}

   inline void enable_screening(const bool enabled) {_enable_screening = enabled;}
   inline void enable_blocking (const bool enabled) {_enable_blocking  = enabled;}

   [[nodiscard]] inline auto get_prefactor() const {return _prefactor;}

   void make_ff_table(unsigned int i,
                      unsigned int f);

   static auto Gsqr(double         Kz,
                    const Subband &isb,
                    const Subband &fsb) -> double;

   [[nodiscard]] auto get_ff_table(unsigned int i,
                                   unsigned int f) const -> arma::vec;
   [[nodiscard]] inline auto get_Kz_table() const {return _Kz;}
};
} // namespace
#endif
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :
