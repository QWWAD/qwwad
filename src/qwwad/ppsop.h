/*=================================================================
       ppsop     PseudoPotential Spin-Orbit Parameters
  =================================================================

   This function returns the spin-orbit parameter (lambda)
   for a variety of atomic species.

   Paul Harrison, April 2000					 */

#ifndef PPSOP_H
#define PPSOP_H
auto lambda(const char *type) -> double;
#endif //PPSOP_H
