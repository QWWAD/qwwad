#include "ppff.h"

#include <cstring>
#include <complex>
#include <armadillo>
#include <gsl/gsl_math.h>
#include "constants.h"
#include "file-io.h"
#include "maths-helpers.h"

using namespace QWWAD;
using namespace constants;

/**
 * This function returns the atomic form factor Vf(q) for the appropriate
 * atomic species
 *
 * \param[in] A0       Lattice constant
 * \param[in] m_per_au number of metres per a.u. of length
 * \param[in] q_sqr    modulus squared of q
 * \param[in]type      atomic species
 */
auto Vf(const double  A0,
        const double  m_per_au,
        double        q_sqr,
        const char   *type) -> double
{
 double A0_au;	/* lattice constant in atomic units	*/
 double a1;
 double a2;
 double a3;
 double a4;
 double a5;
 double a6;
 double Omega;	/* atomic volume, i.e. fcc cube/4	*/
 double Va;	/* the atomic potential			*/
 double	x;	/* an alloy concentration		*/

 /* The Freidel et al. potentials for Si-Ge heterostructures, 
   Phys. Rev. B39 p7974 (1989)					*/

 if(strcmp(type,"SI") == 0)
 {
  q_sqr *= A0*A0/(4*pi*pi);		/* convert q from SI into units of (2*pi/A0) */
  A0_au  = A0/m_per_au;		/* convert A0 from S.I.-->a.u.	*/
  q_sqr *= 4*pi*pi/(A0_au*A0_au);	/* convert q into atomic units  */
  Omega=A0_au*A0_au*A0_au/4;

  a1=106.0686;
  a2=2.2278;
  a3=0.6060;
  a4=-1.9720;
  a5=5.0;
  a6=0.3;

  Va=a1*(q_sqr-a2)*(tanh((a5-q_sqr)/a6)+1)/(2.0*(exp(a3*(q_sqr-a4))+1))/Omega;
 
  return(2*Va*h*c*Rinf);	/*  first factor a.u.--> Rydberg */
				/* second factor Rydberg --> SI  */
 }

 if(strcmp(type,"GE") == 0)
 {
  q_sqr *= A0*A0/(4*pi*pi);		/* convert q from SI into units of (2*pi/A0) */
  A0_au  = A0/m_per_au;		/* convert A0 from S.I.-->a.u.	*/
  q_sqr *= 4*pi*pi/(A0_au*A0_au);	/* convert q into atomic units  */
  Omega=A0_au*A0_au*A0_au/4;

  a1=54.4512;
  a2=2.3592;
  a3=0.7400;
  a4=-0.3800;
  a5=5.0;
  a6=0.3;

  Va=a1*(q_sqr-a2)*(tanh((a5-q_sqr)/a6)+1)/(2.0*(exp(a3*(q_sqr-a4))+1))/Omega;

  return(2*Va*h*c*Rinf);	/* first factor a.u.--> Rydberg	*/
				/* second factor Rydberg --> SI	*/
 }

 /* The continuous potentials of Mader and Zunger, Phys. Rev. B.50
    17393 (1994), for GaAs/AlAs heterostructures		*/

 if(strcmp(type,"GAASmz") == 0)
 {
  q_sqr/=gsl_pow_2(2*pi/A0);		/* convert q from SI into units of (2*pi/A0) */
  A0_au=A0/m_per_au;		/* convert A0 from S.I.-->a.u.	*/
  q_sqr*=gsl_pow_2(2*pi/A0_au);	/* convert q into atomic units  */
  Omega=A0_au*A0_au*A0_au/4;

  Va=(-1.24498*exp(-1.52748*gsl_pow_2(sqrt(q_sqr)-0))
      +0.0366517*exp(-0.959082*gsl_pow_2(sqrt(q_sqr)-2.09782))
      +0.0464357*exp(-0.574047*gsl_pow_2(sqrt(q_sqr)-2.01935))
      -0.0133385*exp(-11.2708*gsl_pow_2(sqrt(q_sqr)-2.93581)))*131.4/Omega;

  return(Va*h*c*Rinf);		/* factor converts Rydberg --> SI */
 }

 if(strcmp(type,"ASGAmz") == 0)
 {
  q_sqr/=gsl_pow_2(2*pi/A0);		/* convert q from SI into units of (2*pi/A0) */
  A0_au=A0/m_per_au;		/* convert A0 from S.I.-->a.u.	*/
  q_sqr*=gsl_pow_2(2*pi/A0_au);	/* convert q into atomic units  */
  Omega=A0_au*A0_au*A0_au/4;

  Va=(-1.0582*exp(-0.959327*gsl_pow_2(sqrt(q_sqr)-0))
      -0.00217627*exp(-6.53145*gsl_pow_2(sqrt(q_sqr)-2.46808))
      -0.0434312*exp(-2.94679*gsl_pow_2(sqrt(q_sqr)-0.851644))
      +0.10569*exp(-0.820922*gsl_pow_2(sqrt(q_sqr)-1.22436)))*145.2/Omega;

  return(Va*h*c*Rinf);		/* factor converts Rydberg --> SI */
 }

 if(strcmp(type,"ALASmz") == 0)
 {
  q_sqr/=gsl_pow_2(2*pi/A0);		/* convert q from SI into units of (2*pi/A0) */
  A0_au=A0/m_per_au;		/* convert A0 from S.I.-->a.u.	*/
  q_sqr*=gsl_pow_2(2*pi/A0_au);	/* convert q into atomic units  */
  Omega=A0_au*A0_au*A0_au/4;

  Va=(-1.32712*exp(-1.59819*gsl_pow_2(sqrt(q_sqr)-0))
      +0.158114*exp(-2.10827*gsl_pow_2(sqrt(q_sqr)-1.77453))
      +0.0601648*exp(-0.527745*gsl_pow_2(sqrt(q_sqr)-2.59550))
      +0.0168167*exp(-11.2708*gsl_pow_2(sqrt(q_sqr)-2.93581)))*111.3
     *(1+0.02*exp(-10*q_sqr))/Omega;

  return(Va*h*c*Rinf);		/* factor converts Rydberg --> SI */
 }

 if(strcmp(type,"ASALmz") == 0)
 {
  q_sqr/=gsl_pow_2(2*pi/A0);		/* convert q from SI into units of (2*pi/A0) */
  A0_au=A0/m_per_au;		/* convert A0 from S.I.-->a.u.	*/
  q_sqr*=gsl_pow_2(2*pi/A0_au);	/* convert q into atomic units  */
  Omega=A0_au*A0_au*A0_au/4;

  Va=(-1.10411*exp(-0.972439*gsl_pow_2(sqrt(q_sqr)-0))
      +0.0174946*exp(-6.53147*gsl_pow_2(sqrt(q_sqr)-2.46793))
      -0.00368081*exp(-5.50601*gsl_pow_2(sqrt(q_sqr)-1.22845))
      +0.0921512*exp(-1.18638*gsl_pow_2(sqrt(q_sqr)-1.35897)))*145.2/Omega;

  return(Va*h*c*Rinf);		/* factor converts Rydberg --> SI */
 }

 /* Alloy definitions for Al(x)Ga(1-x)As using the Mader and Zunger
    potentials, deduces the alloy concentration x from the 7th and 8th 
    characters in the atom type string					*/ 

 if(strncmp(type,"ALGAAS",6) == 0)
 {
  x=(10*(double)(*(type+6)-48)+(double)(*(type+7)-48))/100;
  
  return(x*Vf(A0,m_per_au,q_sqr,"ALASmz")+(1-x)*Vf(A0,m_per_au,q_sqr,"GAASmz"));
 }

 if(strncmp(type,"ASALGA",6) == 0)
 {
  x=(10*(double)(*(type+6)-48)+(double)(*(type+7)-48))/100;
  
  return(x*Vf(A0,m_per_au,q_sqr,"ASALmz")+(1-x)*Vf(A0,m_per_au,q_sqr,"ASGAmz"));
 }



 /* The continuous potentials of Williamson et al. Phys. Rev. B 
    (submitted May 5) (1998)						
  
    Note have to divide the Williamson potentials by a factor of
    2 due to a difference in the definition of the problem	*/

 if(strcmp(type,"GAASw") == 0)
 {
  q_sqr/=gsl_pow_2(2*pi/A0);		/* convert q from SI into units of (2*pi/A0) */
  A0_au=A0/m_per_au;		/* convert A0 from S.I.-->a.u.	*/
  q_sqr*=gsl_pow_2(2*pi/A0_au);	/* convert q into atomic units  */
  Omega=A0_au*A0_au*A0_au/4;

  Va=139478*(q_sqr-2.316)/(3810.60*exp(0.283*q_sqr)-1)/Omega;
  Va/=2;

  return(2*Va*h*c*Rinf);	/* first factor a.u.--> Rydberg	*/
				/* second factor Rydberg --> SI */
 }

 if(strcmp(type,"ASGAw") == 0)
 {
  q_sqr/=gsl_pow_2(2*pi/A0);		/* convert q from SI into units of (2*pi/A0) */
  A0_au=A0/m_per_au;		/* convert A0 from S.I.-->a.u.	*/
  q_sqr*=gsl_pow_2(2*pi/A0_au);	/* convert q into atomic units  */
  Omega=A0_au*A0_au*A0_au/4;

  Va=13.825*(q_sqr-2.878)/(1.169*exp(0.281*q_sqr)-1)/Omega;
  Va/=2;

  return(2*Va*h*c*Rinf);	/* first factor a.u.--> Rydberg */
				/* second factor Rydberg --> SI */
 }
 
  if(strcmp(type,"ALASw") == 0)
 {
  q_sqr/=gsl_pow_2(2*pi/A0);		/* convert q from SI into units of (2*pi/A0) */
  A0_au=A0/m_per_au;		/* convert A0 from S.I.-->a.u.	*/
  q_sqr*=gsl_pow_2(2*pi/A0_au);	/* convert q into atomic units  */
  Omega=A0_au*A0_au*A0_au/4;

  Va=319.275*(q_sqr-2.292)/(13.625*exp(0.315*q_sqr)-1)/Omega;
  Va/=2;

  return(2*Va*h*c*Rinf);	/* first factor a.u.--> Rydberg */
				/* second factor Rydberg --> SI */
 }

 if(strcmp(type,"ASALw") == 0)
 {
  q_sqr/=gsl_pow_2(2*pi/A0);		/* convert q from SI into units of (2*pi/A0) */
  A0_au=A0/m_per_au;		/* convert A0 from S.I.-->a.u.	*/
  q_sqr*=gsl_pow_2(2*pi/A0_au);	/* convert q into atomic units  */
  Omega=A0_au*A0_au*A0_au/4;

  Va=21.993*(q_sqr-2.520)/(1.205*exp(0.336*q_sqr)-1)/Omega;
  Va/=2;

  return(2*Va*h*c*Rinf);	/* first factor a.u.--> Rydberg */
				/* second factor Rydberg --> SI */
 }

 if(strcmp(type,"INASw") == 0)
 {
  q_sqr/=gsl_pow_2(2*pi/A0);		/* convert q from SI into units of (2*pi/A0) */
  A0_au=A0/m_per_au;		/* convert A0 from S.I.-->a.u.	*/
  q_sqr*=gsl_pow_2(2*pi/A0_au);	/* convert q into atomic units  */
  Omega=A0_au*A0_au*A0_au/4;

  Va=107.755*(q_sqr-1.915)/(3.460*exp(0.414*q_sqr)-1)/Omega;
  Va/=2;

  return(2*Va*h*c*Rinf);	/* first factor a.u.--> Rydberg */
				/* second factor Rydberg --> SI */
 }

 if(strcmp(type,"ASINw") == 0)
 {
  q_sqr/=gsl_pow_2(2*pi/A0);		/* convert q from SI into units of (2*pi/A0) */
  A0_au=A0/m_per_au;		/* convert A0 from S.I.-->a.u.	*/
  q_sqr*=gsl_pow_2(2*pi/A0_au);	/* convert q into atomic units  */
  Omega=A0_au*A0_au*A0_au/4;

  Va=49.614*(q_sqr-2.737)/(1.523*exp(0.574*q_sqr)-1)/Omega;
  Va/=2;

  return(2*Va*h*c*Rinf);	/* first factor a.u.--> Rydberg */
				/* second factor Rydberg --> SI */
 }


 /* The Cohen and Bergstresser potentials (in eV) for bulk fcc crystals,
    Phys. Rev. 141 p789 (1966)						*/

 if(strcmp(type,"SIcb") == 0)
 {
  q_sqr/=gsl_pow_2(2*pi/A0);		/* convert q from SI into units of (2*pi/A0) */
  Va=-1.43*(Theta(q_sqr-2.9)-Theta(q_sqr-3.1))		/* Vf(sqrt(3))	*/
     +0.27*(Theta(q_sqr-7.9)-Theta(q_sqr-8.1))		/* Vf(sqrt(8))	*/
     +0.54*(Theta(q_sqr-10.9)-Theta(q_sqr-11.1));	/* Vf(sqrt(11))	*/

  return(Va*e);	/* Convert eV --> SI	*/
 }

 if(strcmp(type,"GEcb") == 0)
 {
  q_sqr/=gsl_pow_2(2*pi/A0);		/* convert q from SI into units of (2*pi/A0) */
  Va=-1.57*(Theta(q_sqr-2.9)-Theta(q_sqr-3.1))		/* Vf(sqrt(3))	*/
     +0.07*(Theta(q_sqr-7.9)-Theta(q_sqr-8.1))		/* Vf(sqrt(8))	*/
     +0.41*(Theta(q_sqr-10.9)-Theta(q_sqr-11.1));	/* Vf(sqrt(11))	*/

  return(Va*e);	/* Convert eV --> SI	*/
 }

 /* These atomic potentials are derived from the Symmetric and
    Antisymmetric potentials of Cohen and Bergstresser.
    Note some of these potentials have definitions for both q=sqrt(4) and
    q=sqrt(8)---these are necessary in the `large basis' method and have
    been deduced by linear interpolation(!) from the existing data.	*/

 if(strcmp(type,"GAAScb") == 0)	/* Ga in GaAs	*/
 {
  q_sqr/=gsl_pow_2(2*pi/A0);		/* convert q from SI into units of (2*pi/A0) */
  Va=-2.04*(Theta(q_sqr-2.9)-Theta(q_sqr-3.1))		/* Vf(sqrt(3))	*/
     -1.51*(Theta(q_sqr-3.9)-Theta(q_sqr-4.1))		/* Vf(sqrt(4))	*/
     -0.10*(Theta(q_sqr-7.9)-Theta(q_sqr-8.1))		/* Vf(sqrt(8))	*/
     +0.34*(Theta(q_sqr-10.9)-Theta(q_sqr-11.1));	/* Vf(sqrt(11))	*/

  return(Va*e);	/* Convert eV --> SI	*/
 }

 if(strcmp(type,"ASGAcb") == 0)	/* As in GaAs	*/
 {
  q_sqr/=gsl_pow_2(2*pi/A0);		/* convert q from SI into units of (2*pi/A0) */
  Va=-1.09*(Theta(q_sqr-2.9)-Theta(q_sqr-3.1))		/* Vf(sqrt(3))	*/
     -0.83*(Theta(q_sqr-3.9)-Theta(q_sqr-4.1))		/* Vf(sqrt(4))	*/
     +0.24*(Theta(q_sqr-7.9)-Theta(q_sqr-8.1))		/* Vf(sqrt(8))	*/
     +0.48*(Theta(q_sqr-10.9)-Theta(q_sqr-11.1));	/* Vf(sqrt(11))	*/

  return(Va*e);	/* Convert eV --> SI	*/
 }

 if(strcmp(type,"INAScb") == 0)	/* In in InAs	*/
 {
  q_sqr/=gsl_pow_2(2*pi/A0);		/* convert q from SI into units of (2*pi/A0) */
  Va=-2.04*(Theta(q_sqr-2.9)-Theta(q_sqr-3.1))		/* Vf(sqrt(3))	*/
     -1.47*(Theta(q_sqr-3.9)-Theta(q_sqr-4.1))		/* Vf(sqrt(4))	*/
     -0.26*(Theta(q_sqr-7.9)-Theta(q_sqr-8.1))		/* Vf(sqrt(8))	*/
     +0.14*(Theta(q_sqr-10.9)-Theta(q_sqr-11.1));	/* Vf(sqrt(11))	*/

  return(Va*e);	/* Convert eV --> SI	*/
 }

 if(strcmp(type,"ASINcb") == 0)	/* As in InAs	*/
 {
  q_sqr/=gsl_pow_2(2*pi/A0);		/* convert q from SI into units of (2*pi/A0) */
  Va=-0.95*(Theta(q_sqr-2.9)-Theta(q_sqr-3.1))		/* Vf(sqrt(3))	*/
     -0.79*(Theta(q_sqr-3.9)-Theta(q_sqr-4.1))		/* Vf(sqrt(4))	*/
     +0.26*(Theta(q_sqr-7.9)-Theta(q_sqr-8.1))		/* Vf(sqrt(8))	*/
     +0.55*(Theta(q_sqr-10.9)-Theta(q_sqr-11.1));	/* Vf(sqrt(11))	*/

  return(Va*e);	/* Convert eV --> SI	*/
 }

 if(strcmp(type,"GAPcb") == 0)	/* Ga in GaP	*/
 {
  q_sqr/=gsl_pow_2(2*pi/A0);		/* convert q from SI into units of (2*pi/A0) */
  Va=-2.31*(Theta(q_sqr-2.9)-Theta(q_sqr-3.1))		/* Vf(sqrt(3))	*/
     -1.56*(Theta(q_sqr-3.9)-Theta(q_sqr-4.1))		/* Vf(sqrt(4))	*/
     -0.06*(Theta(q_sqr-7.9)-Theta(q_sqr-8.1))		/* Vf(sqrt(8))	*/
     +0.34*(Theta(q_sqr-10.9)-Theta(q_sqr-11.1));	/* Vf(sqrt(11))	*/

  return(Va*e);	/* Convert eV --> SI	*/
 }

 if(strcmp(type,"PGAcb") == 0)	/* P in GaP	*/
 {
  q_sqr/=gsl_pow_2(2*pi/A0);		/* convert q from SI into units of (2*pi/A0) */
  Va=-0.68*(Theta(q_sqr-2.9)-Theta(q_sqr-3.1))		/* Vf(sqrt(3))	*/
     -0.61*(Theta(q_sqr-3.9)-Theta(q_sqr-4.1))		/* Vf(sqrt(4))	*/
     +0.47*(Theta(q_sqr-7.9)-Theta(q_sqr-8.1))		/* Vf(sqrt(8))	*/
     +0.61*(Theta(q_sqr-10.9)-Theta(q_sqr-11.1));	/* Vf(sqrt(11))	*/

  return(Va*e);	/* Convert eV --> SI	*/
 }

 if(strcmp(type,"INPcb") == 0)	/* In in InP	*/
 {
  q_sqr/=gsl_pow_2(2*pi/A0);		/* convert q from SI into units of (2*pi/A0) */
  Va=-2.04*(Theta(q_sqr-2.9)-Theta(q_sqr-3.1))		/* Vf(sqrt(3))	*/
     -1.51*(Theta(q_sqr-3.9)-Theta(q_sqr-4.1))		/* Vf(sqrt(4))	*/
     -0.10*(Theta(q_sqr-7.9)-Theta(q_sqr-8.1))		/* Vf(sqrt(8))	*/
     +0.34*(Theta(q_sqr-10.9)-Theta(q_sqr-11.1));	/* Vf(sqrt(11))	*/

  return(Va*e);	/* Convert eV --> SI	*/
 }

 if(strcmp(type,"PINcb") == 0)	/* P in InP	*/
 {
  q_sqr/=gsl_pow_2(2*pi/A0);		/* convert q from SI into units of (2*pi/A0) */
  Va=-1.09*(Theta(q_sqr-2.9)-Theta(q_sqr-3.1))		/* Vf(sqrt(3))	*/
     -0.83*(Theta(q_sqr-3.9)-Theta(q_sqr-4.1))		/* Vf(sqrt(4))	*/
     +0.24*(Theta(q_sqr-7.9)-Theta(q_sqr-8.1))		/* Vf(sqrt(8))	*/
     +0.48*(Theta(q_sqr-10.9)-Theta(q_sqr-11.1));	/* Vf(sqrt(11))	*/

  return(Va*e);	/* Convert eV --> SI	*/
 }

 if(strcmp(type,"CDTEcb") == 0)	/* Cd in CdTe	*/
 {
  q_sqr/=gsl_pow_2(2*pi/A0);		/* convert q from SI into units of (2*pi/A0) */
  Va=-0.175*(Theta(q_sqr-2.9)-Theta(q_sqr-3.1))		/* Vf(sqrt(3))	*/
     -0.12*(Theta(q_sqr-3.9)-Theta(q_sqr-4.1))		/* Vf(sqrt(4))	*/
     -0.03*(Theta(q_sqr-7.9)-Theta(q_sqr-8.1))		/* Vf(sqrt(8))	*/
     +0.00*(Theta(q_sqr-10.9)-Theta(q_sqr-11.1));	/* Vf(sqrt(11))	*/

  return(Va*h*c*Rinf);	/* Convert Rydberg --> SI */
 }

 if(strcmp(type,"TECDcb") == 0)	/* Te in CdTe	*/
 {
  q_sqr/=gsl_pow_2(2*pi/A0);		/* convert q from SI into units of (2*pi/A0) */
  Va=-0.025*(Theta(q_sqr-2.9)-Theta(q_sqr-3.1))		/* Vf(sqrt(3))	*/
     -0.030*(Theta(q_sqr-3.9)-Theta(q_sqr-4.1))		/* Vf(sqrt(4))	*/
     +0.030*(Theta(q_sqr-7.9)-Theta(q_sqr-8.1))		/* Vf(sqrt(8))	*/
     +0.040*(Theta(q_sqr-10.9)-Theta(q_sqr-11.1));	/* Vf(sqrt(11))	*/

  return(Va*h*c*Rinf);	/* Convert Rydberg --> SI */
 }

 /* Chelikowsky and Cohen, Phys. Rev. B14 p556 (1976)	*/

 if(strcmp(type,"GAAScc") == 0)	/* Ga in GaAs	*/
 {
  q_sqr/=gsl_pow_2(2*pi/A0);		/* convert q from SI into units of (2*pi/A0) */
  Va=-0.135*(Theta(q_sqr-2.9)-Theta(q_sqr-3.1))		/* Vf(sqrt(3))	*/
     -0.098*(Theta(q_sqr-3.9)-Theta(q_sqr-4.1))		/* Vf(sqrt(4))	*/
     +0.000*(Theta(q_sqr-7.9)-Theta(q_sqr-8.1))		/* Vf(sqrt(8))	*/
     +0.033*(Theta(q_sqr-10.9)-Theta(q_sqr-11.1));	/* Vf(sqrt(11))	*/

  return(Va*h*c*Rinf);	/* Convert Rydberg --> SI	*/
 }

 if(strcmp(type,"ASGAcc") == 0)	/* As in GaAs	*/
 {
  q_sqr/=gsl_pow_2(2*pi/A0);		/* convert q from SI into units of (2*pi/A0) */
  Va=-0.080*(Theta(q_sqr-2.9)-Theta(q_sqr-3.1))		/* Vf(sqrt(3))	*/
     -0.060*(Theta(q_sqr-3.9)-Theta(q_sqr-4.1))		/* Vf(sqrt(4))	*/
     +0.014*(Theta(q_sqr-7.9)-Theta(q_sqr-8.1))		/* Vf(sqrt(8))	*/
     +0.034*(Theta(q_sqr-10.9)-Theta(q_sqr-11.1));	/* Vf(sqrt(11))	*/

  return(Va*h*c*Rinf);	/* Convert Rydberg --> SI	*/
 }

 /* Cd(1-x)Mn(x)Te potentials from Fei Long et al., J. Appl. Phys. 79,
  * p6939 (1996)							*/

 if(strcmp(type,"CDTE") == 0)	/* Cd in CdTe	*/
 {
  q_sqr/=gsl_pow_2(2*pi/A0);		/* convert q from SI into units of (2*pi/A0) */
  Va=-0.175*(Theta(q_sqr-2.9)-Theta(q_sqr-3.1))		/* Vf(sqrt(3))	*/
     -0.123*(Theta(q_sqr-3.9)-Theta(q_sqr-4.1))		/* Vf(sqrt(4))	*/
     -0.0279*(Theta(q_sqr-7.9)-Theta(q_sqr-8.1))	/* Vf(sqrt(8))	*/
     +0.0213*(Theta(q_sqr-10.9)-Theta(q_sqr-11.1));	/* Vf(sqrt(11))	*/

  return(Va*h*c*Rinf);	/* Convert Rydberg --> SI */
 }

 if(strcmp(type,"TECD") == 0)	/* Te in CdTe	*/
 {
  q_sqr/=gsl_pow_2(2*pi/A0);		/* convert q from SI into units of (2*pi/A0) */
  Va=-0.055*(Theta(q_sqr-2.9)-Theta(q_sqr-3.1))		/* Vf(sqrt(3))	*/
     -0.0518*(Theta(q_sqr-3.9)-Theta(q_sqr-4.1))	/* Vf(sqrt(4))	*/
     +0.0229*(Theta(q_sqr-7.9)-Theta(q_sqr-8.1))	/* Vf(sqrt(8))	*/
     +0.0598*(Theta(q_sqr-10.9)-Theta(q_sqr-11.1));	/* Vf(sqrt(11))	*/

  return(Va*h*c*Rinf);	/* Convert Rydberg --> SI */
 }


 printf("Error atom type '%s' undefined!\n",type);
 exit(EXIT_FAILURE);
}
 
/**
 * \brief Reads the atomic species into memory
 */
auto read_atoms(const char * filename) -> std::vector<atom>
{
 auto *Fatoms=fopen(filename,"r"); // file pointer to wavefunction file
 if(Fatoms == nullptr) {
     std::ostringstream oss;
     oss << "Cannot open input file " << filename;
     throw std::runtime_error(oss.str());
 }

 /* Read in the first line and hence the number of atoms	*/
 int n_atoms=0;
 size_t n_read = fscanf(Fatoms,"%d", &n_atoms);
 
 /* Allocate memory for atom definitions	*/
 std::vector<atom> atoms;
 if (n_read == 1) {
   atoms.resize(n_atoms);
 }
 else {
     throw std::runtime_error("Could not read number of atoms");
 }

 double rx;
 double ry;
 double rz;
 int ia=0;
 while((fscanf(Fatoms,"%s %lf %lf %lf",atoms[ia].type,
        &rx,&ry,&rz))!=EOF)
 {
     arma::vec r(3);
     r(0) = rx;
     r(1) = ry;
     r(2) = rz;

     /* Convert atomic positions from Angstrom into S.I. units	*/
     r *= 1e-10;

     atoms[ia].r = r;
     ia++;
 }
 fclose(Fatoms);

 return atoms;
}   

/* This function reads the reciprocal lattice vectors (defined in
   the file G.r) into the array G[] and then converts into SI units */
auto
read_rlv(double A0) -> std::vector<arma::vec>
{
    arma::vec Gx;
    arma::vec Gy;
    arma::vec Gz;
    read_table("G.r", Gx, Gy, Gz);

    std::vector<arma::vec> G;
    size_t N = Gx.size();
    for(unsigned int iG = 0; iG < N; ++iG)
    {
        arma::vec _G(3);
        _G(0) = Gx[iG];
        _G(1) = Gy[iG];
        _G(2) = Gz[iG];
        _G *= 2.0*pi/A0;
        G.push_back(_G);
    }

    return G;
}
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :
