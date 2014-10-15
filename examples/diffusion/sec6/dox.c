double D_of_x(D,x,z,t)

/* In its most general form this user supplied function expresses
   the dependence of the diffusion coefficient D on all variables
   possible, i.e. D itself, the concentration x, the position z
   and the time t.  However in this example D depends only on
   the concentration x.                                           */

double D;  /* diffusion coefficient      */
double x;  /* concentration of diffusant */
double z;  /* position                   */
double t;  /* time                       */

{
 /* Remember D must be in SI units, a typical constant is 1 Angstrom^2/s*/

 return(1e-20*x*x);
}

