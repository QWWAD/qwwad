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

    typedef std::pair<unsigned int, unsigned int> map_key;
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

   [[nodiscard]] double get_Eki_min (const unsigned int isb,
                                     const unsigned int fsb) const;

   [[nodiscard]] double get_ki_min (const unsigned int isb,
                                    const unsigned int fsb) const;

   [[nodiscard]] double get_ki_cutoff(const unsigned int isb,
                                      const unsigned int fsb) const;

   double get_rate_ki(const unsigned int isb,
                      const unsigned int fsb,
                      const double       ki);

   IntersubbandTransition get_transition(const unsigned int isb,
                                         const unsigned int fsb);

   [[nodiscard]] inline decltype(_lambda_s_sq) get_screening_length() const {return _lambda_s_sq;}

   inline void set_ki_samples(const decltype(_nki) nki) {_nki = nki;}
   void set_phonon_samples(const size_t nKz);

   inline decltype(_dKz) get_dKz() {return _dKz;}

   inline void enable_screening(const bool enabled) {_enable_screening = enabled;}
   inline void enable_blocking (const bool enabled) {_enable_blocking  = enabled;}

   [[nodiscard]] inline decltype(_prefactor) get_prefactor() const {return _prefactor;}

   void make_ff_table(const unsigned int i,
                      const unsigned int f);

   double Gsqr(const double   Kz,
               const Subband &isb,
               const Subband &fsb);

   [[nodiscard]] arma::vec get_ff_table(const unsigned int i, const unsigned int f) const;
   [[nodiscard]] inline decltype(_Kz)  get_Kz_table() const {return _Kz;}
};
} // namespace
#endif
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :
