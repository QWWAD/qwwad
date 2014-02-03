#ifndef DOX_H
#define DOX_H

/**
 * Calculates the diffusion coefficient at a given depth in a structure
 *
 * \param[in] D0 A constant diffusion coefficient [m^2/s]
 * \param[in] x  Diffusant concentration [a.u.]
 * \param[in] z  Depth into structure [m]
 * \param[in] t  Time [s]
 *
 * \returns diffusion coefficient [m^2/s]
 *
 * \details In its most general form this function could express
 *          the dependence of the diffusion coefficient D on all variables
 *          possible, i.e. D itself, the concentration x, the position z
 *          and the time t.  However in this example D depends only on
 *          the position.
 */
double D_of_x(const double D0,
              const double x,
              const double z,
              const double t)
{
 /* Suppress compiler warnings about unused variables */
 (void)D0;
 (void)x;
 (void)t;

 const double sigma = 600;  /* Standard deviation [angstrom] */
 const double z0    = 1800; /* Location of maximum vacancy concentration [angstrom] */
 const double k     = 10;   /* Proportionality constant [angstrom] */

 /* Find diffusion coefficient [angstrom^2/s] at this depth using Eq. 4.16, QWWAD3 */
 const double diff_coeff = k*exp(-sqr((z/1e-10-z0)/sigma)/2);

 return diff_coeff * 1e-20; /* return in m^2/s */
}
#endif
