/**
 * \file   intersubband-transition.h
 * \author Alex Valavanis <a.valavanis@leeds.ac.uk>
 * \brief  Table of scattering rates between a pair of subbands
 */

#ifndef QWWAD_INTERSUBBAND_TRANSITION
#define QWWAD_INTERSUBBAND_TRANSITION

#include "../qclsim-subband.h"

namespace QWWAD {
/**
 * \brief A generalised table of scattering rates between a pair of subbands
 */
class IntersubbandTransition {
private:
    Subband _isb; ///< The initial subband
    Subband _fsb; ///< The final subband

    std::valarray<double> _ki;  ///< Array of initial wave-vectors     [1/m]
    //std::valarray<double> _kf;     ///< Array of corresponding energy-conserving final wave-vectors [1/m]

    /**
     * \brief Array of scattering rates [1/s]
     *
     * \details Each element in the array gives the total scattering rate from the
     *          initial state to ALL permitted states in the final subband
     */
    std::valarray<double> _Wif;

    // Derived properties
    std::valarray<double> _Eki; ///< Array of initial kinetic energies [1/m]

public:
    IntersubbandTransition(const decltype(_isb)    isb,
                           const decltype(_fsb)    fsb,
                           const decltype(_ki)     ki,
      //                     const decltype(_kf)     kf,
                           const decltype(_Wif) Wif_ki);

    inline decltype(_ki)  get_ki_table        () const {return _ki;}
    inline decltype(_Eki) get_Eki_table       () const {return _Eki;}
    inline decltype(_Eki) get_Ei_total_table  () const {return _Eki + _isb.get_E();}
    inline decltype(_Wif) get_rate_table () const {return _Wif;}
    inline double get_ki_by_index        (unsigned int ik) const {return _ki[ik];}
//    inline double get_kf_at_ki_by_index  (unsigned int ik) const {return _kf[ik];}
    inline double get_rate_at_ki_by_index(unsigned int ik) const {return _Wif[ik];}

    double get_average_rate(const double Te) const;
};
} // namespace
#endif
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :
