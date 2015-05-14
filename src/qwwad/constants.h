/**
 * \file    constants.h
 * \brief   Frequently-used physical and mathematical constants
 * \author  Alex Valavanis <a.valavanis@leeds.ac.uk>
 * \details This is mostly just a wrapper for stuff in the GNU Scientific
 *          library.  All units are SI
 */

#ifndef PHYSICAL_CONSTANTS_H
#define PHYSICAL_CONSTANTS_H

#if HAVE_CONFIG_H
# include "config.h"
#endif

#include <gsl/gsl_const_mksa.h>
#include <gsl/gsl_const_num.h>
#include <gsl/gsl_math.h>

#if __cplusplus
namespace QWWAD{
namespace constants{
#endif

/** Electron Rest Mass [kg] */
const double me = GSL_CONST_MKSA_MASS_ELECTRON;

/** Electronic charge [C] */
const double e = GSL_CONST_MKSA_ELECTRON_CHARGE;

/** Reduced Planck's constant [Js] */
const double hBar = GSL_CONST_MKSA_PLANCKS_CONSTANT_HBAR;

/** Planck's constant [Js] */
const double h = GSL_CONST_MKSA_PLANCKS_CONSTANT_H;

/** Boltzmann's constant [J/K] */
const double kB = GSL_CONST_MKSA_BOLTZMANN;

/** \f$\pi\f$ */
const double pi = M_PI;

/** Permittivity of free space [F/m] */
const double eps0 = GSL_CONST_MKSA_VACUUM_PERMITTIVITY;

/** Speed of light [m/s] */
const double c = GSL_CONST_MKSA_SPEED_OF_LIGHT;

/** Avogadro's number */
const double Na = GSL_CONST_NUM_AVOGADRO;

/** Rydberg inverse wavelength [1/m] */
/** TODO: Set this using external lib */
/** const double Rinf = GSL_CONST_MKSA_RYDBERG / (h*c); */
const double Rinf = 1.097373156854e7;

/** Bohr magneton */
const double mu_b = GSL_CONST_MKSA_BOHR_MAGNETON;
#if __cplusplus
}// namespace constants
}// namespace Leeds
#endif
#endif

/* vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :
 */
