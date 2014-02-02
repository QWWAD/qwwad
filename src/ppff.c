/*=================================================================
       ppff     PseudoPotential Form Factors
  =================================================================

   This function returns the atomic formfactor (Fourier Transform 
   of the atomic potential) for a variety of atomic species.

   Paul Harrison, July 1998					 */


double
Vf(A0,m_per_au,q_sqr,type)

/* This function returns the atomic form factor Vf(q) for the appropriate
   atomic species */

double	A0;		/* Lattice constant     		*/
double	m_per_au;	/* number of metres per a.u. of length	*/
double	q_sqr;		/* modulus squared of q 		*/
char	type[];		/* atomic species			*/
{
 double	sqrt();
 double A0_au;	/* lattice constant in atomic units	*/
 double a1,a2,a3,a4,a5,a6;
 double Omega;	/* atomic volume, i.e. fcc cube/4	*/
 double Va;	/* the atomic potential			*/
 double	x;	/* an alloy concentration		*/






 /* The Freidel et al. potentials for Si-Ge heterostructures, 
   Phys. Rev. B39 p7974 (1989)					*/

 if(!strcmp(type,"SI"))
 {
  q_sqr/=sqr(2*pi/A0);		/* convert q from SI into units of (2*pi/A0) */
  A0_au=A0/m_per_au;		/* convert A0 from S.I.-->a.u.	*/
  q_sqr*=sqr(2*pi/A0_au);	/* convert q into atomic units  */
  Omega=A0_au*A0_au*A0_au/4;

  a1=106.0686;
  a2=2.2278;
  a3=0.6060;
  a4=-1.9720;
  a5=5.0;
  a6=0.3;

  Va=a1*(q_sqr-a2)*(tanh((a5-q_sqr)/a6)+1)/(2.0*(exp(a3*(q_sqr-a4))+1))/Omega;
 
  return(2*Va*h*c0*Rinf);	/*  first factor a.u.--> Rydberg */
				/* second factor Rydberg --> SI  */
 }

 if(!strcmp(type,"GE"))
 {
  q_sqr/=sqr(2*pi/A0);		/* convert q from SI into units of (2*pi/A0) */
  A0_au=A0/m_per_au;		/* convert A0 from S.I.-->a.u.	*/
  q_sqr*=sqr(2*pi/A0_au);	/* convert q into atomic units  */
  Omega=A0_au*A0_au*A0_au/4;

  a1=54.4512;
  a2=2.3592;
  a3=0.7400;
  a4=-0.3800;
  a5=5.0;
  a6=0.3;

  Va=a1*(q_sqr-a2)*(tanh((a5-q_sqr)/a6)+1)/(2.0*(exp(a3*(q_sqr-a4))+1))/Omega;

  return(2*Va*h*c0*Rinf);	/* first factor a.u.--> Rydberg	*/
				/* second factor Rydberg --> SI	*/
 }

 /* The continuous potentials of Mader and Zunger, Phys. Rev. B.50
    17393 (1994), for GaAs/AlAs heterostructures		*/

 if(!strcmp(type,"GAASmz"))
 {
  q_sqr/=sqr(2*pi/A0);		/* convert q from SI into units of (2*pi/A0) */
  A0_au=A0/m_per_au;		/* convert A0 from S.I.-->a.u.	*/
  q_sqr*=sqr(2*pi/A0_au);	/* convert q into atomic units  */
  Omega=A0_au*A0_au*A0_au/4;

  Va=(-1.24498*exp(-1.52748*sqr(sqrt(q_sqr)-0))
      +0.0366517*exp(-0.959082*sqr(sqrt(q_sqr)-2.09782))
      +0.0464357*exp(-0.574047*sqr(sqrt(q_sqr)-2.01935))
      -0.0133385*exp(-11.2708*sqr(sqrt(q_sqr)-2.93581)))*131.4/Omega;

  return(Va*h*c0*Rinf);		/* factor converts Rydberg --> SI */
 }

 if(!strcmp(type,"ASGAmz"))
 {
  q_sqr/=sqr(2*pi/A0);		/* convert q from SI into units of (2*pi/A0) */
  A0_au=A0/m_per_au;		/* convert A0 from S.I.-->a.u.	*/
  q_sqr*=sqr(2*pi/A0_au);	/* convert q into atomic units  */
  Omega=A0_au*A0_au*A0_au/4;

  Va=(-1.0582*exp(-0.959327*sqr(sqrt(q_sqr)-0))
      -0.00217627*exp(-6.53145*sqr(sqrt(q_sqr)-2.46808))
      -0.0434312*exp(-2.94679*sqr(sqrt(q_sqr)-0.851644))
      +0.10569*exp(-0.820922*sqr(sqrt(q_sqr)-1.22436)))*145.2/Omega;

  return(Va*h*c0*Rinf);		/* factor converts Rydberg --> SI */
 }

 if(!strcmp(type,"ALASmz"))
 {
  q_sqr/=sqr(2*pi/A0);		/* convert q from SI into units of (2*pi/A0) */
  A0_au=A0/m_per_au;		/* convert A0 from S.I.-->a.u.	*/
  q_sqr*=sqr(2*pi/A0_au);	/* convert q into atomic units  */
  Omega=A0_au*A0_au*A0_au/4;

  Va=(-1.32712*exp(-1.59819*sqr(sqrt(q_sqr)-0))
      +0.158114*exp(-2.10827*sqr(sqrt(q_sqr)-1.77453))
      +0.0601648*exp(-0.527745*sqr(sqrt(q_sqr)-2.59550))
      +0.0168167*exp(-11.2708*sqr(sqrt(q_sqr)-2.93581)))*111.3
     *(1+0.02*exp(-10*q_sqr))/Omega;

  return(Va*h*c0*Rinf);		/* factor converts Rydberg --> SI */
 }

 if(!strcmp(type,"ASALmz"))
 {
  q_sqr/=sqr(2*pi/A0);		/* convert q from SI into units of (2*pi/A0) */
  A0_au=A0/m_per_au;		/* convert A0 from S.I.-->a.u.	*/
  q_sqr*=sqr(2*pi/A0_au);	/* convert q into atomic units  */
  Omega=A0_au*A0_au*A0_au/4;

  Va=(-1.10411*exp(-0.972439*sqr(sqrt(q_sqr)-0))
      +0.0174946*exp(-6.53147*sqr(sqrt(q_sqr)-2.46793))
      -0.00368081*exp(-5.50601*sqr(sqrt(q_sqr)-1.22845))
      +0.0921512*exp(-1.18638*sqr(sqrt(q_sqr)-1.35897)))*145.2/Omega;

  return(Va*h*c0*Rinf);		/* factor converts Rydberg --> SI */
 }

 /* Alloy definitions for Al(x)Ga(1-x)As using the Mader and Zunger
    potentials, deduces the alloy concentration x from the 7th and 8th 
    characters in the atom type string					*/ 

 if(!strncmp(type,"ALGAAS",6))
 {
  x=(10*(double)(*(type+6)-48)+(double)(*(type+7)-48))/100;
  
  return(x*Vf(A0,m_per_au,q_sqr,"ALASmz")+(1-x)*Vf(A0,m_per_au,q_sqr,"GAASmz"));
 }

 if(!strncmp(type,"ASALGA",6))
 {
  x=(10*(double)(*(type+6)-48)+(double)(*(type+7)-48))/100;
  
  return(x*Vf(A0,m_per_au,q_sqr,"ASALmz")+(1-x)*Vf(A0,m_per_au,q_sqr,"ASGAmz"));
 }



 /* The continuous potentials of Williamson et al. Phys. Rev. B 
    (submitted May 5) (1998)						
  
    Note have to divide the Williamson potentials by a factor of
    2 due to a difference in the definition of the problem	*/

 if(!strcmp(type,"GAASw"))
 {
  q_sqr/=sqr(2*pi/A0);		/* convert q from SI into units of (2*pi/A0) */
  A0_au=A0/m_per_au;		/* convert A0 from S.I.-->a.u.	*/
  q_sqr*=sqr(2*pi/A0_au);	/* convert q into atomic units  */
  Omega=A0_au*A0_au*A0_au/4;

  Va=139478*(q_sqr-2.316)/(3810.60*exp(0.283*q_sqr)-1)/Omega;
  Va/=2;

  return(2*Va*h*c0*Rinf);	/* first factor a.u.--> Rydberg	*/
				/* second factor Rydberg --> SI */
 }

 if(!strcmp(type,"ASGAw"))
 {
  q_sqr/=sqr(2*pi/A0);		/* convert q from SI into units of (2*pi/A0) */
  A0_au=A0/m_per_au;		/* convert A0 from S.I.-->a.u.	*/
  q_sqr*=sqr(2*pi/A0_au);	/* convert q into atomic units  */
  Omega=A0_au*A0_au*A0_au/4;

  Va=13.825*(q_sqr-2.878)/(1.169*exp(0.281*q_sqr)-1)/Omega;
  Va/=2;

  return(2*Va*h*c0*Rinf);	/* first factor a.u.--> Rydberg */
				/* second factor Rydberg --> SI */
 }
 
  if(!strcmp(type,"ALASw"))
 {
  q_sqr/=sqr(2*pi/A0);		/* convert q from SI into units of (2*pi/A0) */
  A0_au=A0/m_per_au;		/* convert A0 from S.I.-->a.u.	*/
  q_sqr*=sqr(2*pi/A0_au);	/* convert q into atomic units  */
  Omega=A0_au*A0_au*A0_au/4;

  Va=319.275*(q_sqr-2.292)/(13.625*exp(0.315*q_sqr)-1)/Omega;
  Va/=2;

  return(2*Va*h*c0*Rinf);	/* first factor a.u.--> Rydberg */
				/* second factor Rydberg --> SI */
 }

 if(!strcmp(type,"ASALw"))
 {
  q_sqr/=sqr(2*pi/A0);		/* convert q from SI into units of (2*pi/A0) */
  A0_au=A0/m_per_au;		/* convert A0 from S.I.-->a.u.	*/
  q_sqr*=sqr(2*pi/A0_au);	/* convert q into atomic units  */
  Omega=A0_au*A0_au*A0_au/4;

  Va=21.993*(q_sqr-2.520)/(1.205*exp(0.336*q_sqr)-1)/Omega;
  Va/=2;

  return(2*Va*h*c0*Rinf);	/* first factor a.u.--> Rydberg */
				/* second factor Rydberg --> SI */
 }

 if(!strcmp(type,"INASw"))
 {
  q_sqr/=sqr(2*pi/A0);		/* convert q from SI into units of (2*pi/A0) */
  A0_au=A0/m_per_au;		/* convert A0 from S.I.-->a.u.	*/
  q_sqr*=sqr(2*pi/A0_au);	/* convert q into atomic units  */
  Omega=A0_au*A0_au*A0_au/4;

  Va=107.755*(q_sqr-1.915)/(3.460*exp(0.414*q_sqr)-1)/Omega;
  Va/=2;

  return(2*Va*h*c0*Rinf);	/* first factor a.u.--> Rydberg */
				/* second factor Rydberg --> SI */
 }

 if(!strcmp(type,"ASINw"))
 {
  q_sqr/=sqr(2*pi/A0);		/* convert q from SI into units of (2*pi/A0) */
  A0_au=A0/m_per_au;		/* convert A0 from S.I.-->a.u.	*/
  q_sqr*=sqr(2*pi/A0_au);	/* convert q into atomic units  */
  Omega=A0_au*A0_au*A0_au/4;

  Va=49.614*(q_sqr-2.737)/(1.523*exp(0.574*q_sqr)-1)/Omega;
  Va/=2;

  return(2*Va*h*c0*Rinf);	/* first factor a.u.--> Rydberg */
				/* second factor Rydberg --> SI */
 }


 /* The Cohen and Bergstresser potentials (in eV) for bulk fcc crystals,
    Phys. Rev. 141 p789 (1966)						*/

 if(!strcmp(type,"SIcb"))
 {
  q_sqr/=sqr(2*pi/A0);		/* convert q from SI into units of (2*pi/A0) */
  Va=-1.43*(Theta(q_sqr-2.9)-Theta(q_sqr-3.1))		/* Vf(sqrt(3))	*/
     +0.27*(Theta(q_sqr-7.9)-Theta(q_sqr-8.1))		/* Vf(sqrt(8))	*/
     +0.54*(Theta(q_sqr-10.9)-Theta(q_sqr-11.1));	/* Vf(sqrt(11))	*/

  return(Va*e_0);	/* Convert eV --> SI	*/
 }

 if(!strcmp(type,"GEcb"))
 {
  q_sqr/=sqr(2*pi/A0);		/* convert q from SI into units of (2*pi/A0) */
  Va=-1.57*(Theta(q_sqr-2.9)-Theta(q_sqr-3.1))		/* Vf(sqrt(3))	*/
     +0.07*(Theta(q_sqr-7.9)-Theta(q_sqr-8.1))		/* Vf(sqrt(8))	*/
     +0.41*(Theta(q_sqr-10.9)-Theta(q_sqr-11.1));	/* Vf(sqrt(11))	*/

  return(Va*e_0);	/* Convert eV --> SI	*/
 }

 /* These atomic potentials are derived from the Symmetric and
    Antisymmetric potentials of Cohen and Bergstresser.
    Note some of these potentials have definitions for both q=sqrt(4) and
    q=sqrt(8)---these are necessary in the `large basis' method and have
    been deduced by linear interpolation(!) from the existing data.	*/

 if(!strcmp(type,"GAAScb"))	/* Ga in GaAs	*/
 {
  q_sqr/=sqr(2*pi/A0);		/* convert q from SI into units of (2*pi/A0) */
  Va=-2.04*(Theta(q_sqr-2.9)-Theta(q_sqr-3.1))		/* Vf(sqrt(3))	*/
     -1.51*(Theta(q_sqr-3.9)-Theta(q_sqr-4.1))		/* Vf(sqrt(4))	*/
     -0.10*(Theta(q_sqr-7.9)-Theta(q_sqr-8.1))		/* Vf(sqrt(8))	*/
     +0.34*(Theta(q_sqr-10.9)-Theta(q_sqr-11.1));	/* Vf(sqrt(11))	*/

  return(Va*e_0);	/* Convert eV --> SI	*/
 }

 if(!strcmp(type,"ASGAcb"))	/* As in GaAs	*/
 {
  q_sqr/=sqr(2*pi/A0);		/* convert q from SI into units of (2*pi/A0) */
  Va=-1.09*(Theta(q_sqr-2.9)-Theta(q_sqr-3.1))		/* Vf(sqrt(3))	*/
     -0.83*(Theta(q_sqr-3.9)-Theta(q_sqr-4.1))		/* Vf(sqrt(4))	*/
     +0.24*(Theta(q_sqr-7.9)-Theta(q_sqr-8.1))		/* Vf(sqrt(8))	*/
     +0.48*(Theta(q_sqr-10.9)-Theta(q_sqr-11.1));	/* Vf(sqrt(11))	*/

  return(Va*e_0);	/* Convert eV --> SI	*/
 }

 if(!strcmp(type,"INAScb"))	/* In in InAs	*/
 {
  q_sqr/=sqr(2*pi/A0);		/* convert q from SI into units of (2*pi/A0) */
  Va=-2.04*(Theta(q_sqr-2.9)-Theta(q_sqr-3.1))		/* Vf(sqrt(3))	*/
     -1.47*(Theta(q_sqr-3.9)-Theta(q_sqr-4.1))		/* Vf(sqrt(4))	*/
     -0.26*(Theta(q_sqr-7.9)-Theta(q_sqr-8.1))		/* Vf(sqrt(8))	*/
     +0.14*(Theta(q_sqr-10.9)-Theta(q_sqr-11.1));	/* Vf(sqrt(11))	*/

  return(Va*e_0);	/* Convert eV --> SI	*/
 }

 if(!strcmp(type,"ASINcb"))	/* As in InAs	*/
 {
  q_sqr/=sqr(2*pi/A0);		/* convert q from SI into units of (2*pi/A0) */
  Va=-0.95*(Theta(q_sqr-2.9)-Theta(q_sqr-3.1))		/* Vf(sqrt(3))	*/
     -0.79*(Theta(q_sqr-3.9)-Theta(q_sqr-4.1))		/* Vf(sqrt(4))	*/
     +0.26*(Theta(q_sqr-7.9)-Theta(q_sqr-8.1))		/* Vf(sqrt(8))	*/
     +0.55*(Theta(q_sqr-10.9)-Theta(q_sqr-11.1));	/* Vf(sqrt(11))	*/

  return(Va*e_0);	/* Convert eV --> SI	*/
 }

 if(!strcmp(type,"GAPcb"))	/* Ga in GaP	*/
 {
  q_sqr/=sqr(2*pi/A0);		/* convert q from SI into units of (2*pi/A0) */
  Va=-2.31*(Theta(q_sqr-2.9)-Theta(q_sqr-3.1))		/* Vf(sqrt(3))	*/
     -1.56*(Theta(q_sqr-3.9)-Theta(q_sqr-4.1))		/* Vf(sqrt(4))	*/
     -0.06*(Theta(q_sqr-7.9)-Theta(q_sqr-8.1))		/* Vf(sqrt(8))	*/
     +0.34*(Theta(q_sqr-10.9)-Theta(q_sqr-11.1));	/* Vf(sqrt(11))	*/

  return(Va*e_0);	/* Convert eV --> SI	*/
 }

 if(!strcmp(type,"PGAcb"))	/* P in GaP	*/
 {
  q_sqr/=sqr(2*pi/A0);		/* convert q from SI into units of (2*pi/A0) */
  Va=-0.68*(Theta(q_sqr-2.9)-Theta(q_sqr-3.1))		/* Vf(sqrt(3))	*/
     -0.61*(Theta(q_sqr-3.9)-Theta(q_sqr-4.1))		/* Vf(sqrt(4))	*/
     +0.47*(Theta(q_sqr-7.9)-Theta(q_sqr-8.1))		/* Vf(sqrt(8))	*/
     +0.61*(Theta(q_sqr-10.9)-Theta(q_sqr-11.1));	/* Vf(sqrt(11))	*/

  return(Va*e_0);	/* Convert eV --> SI	*/
 }

 if(!strcmp(type,"INPcb"))	/* In in InP	*/
 {
  q_sqr/=sqr(2*pi/A0);		/* convert q from SI into units of (2*pi/A0) */
  Va=-2.04*(Theta(q_sqr-2.9)-Theta(q_sqr-3.1))		/* Vf(sqrt(3))	*/
     -1.51*(Theta(q_sqr-3.9)-Theta(q_sqr-4.1))		/* Vf(sqrt(4))	*/
     -0.10*(Theta(q_sqr-7.9)-Theta(q_sqr-8.1))		/* Vf(sqrt(8))	*/
     +0.34*(Theta(q_sqr-10.9)-Theta(q_sqr-11.1));	/* Vf(sqrt(11))	*/

  return(Va*e_0);	/* Convert eV --> SI	*/
 }

 if(!strcmp(type,"PINcb"))	/* P in InP	*/
 {
  q_sqr/=sqr(2*pi/A0);		/* convert q from SI into units of (2*pi/A0) */
  Va=-1.09*(Theta(q_sqr-2.9)-Theta(q_sqr-3.1))		/* Vf(sqrt(3))	*/
     -0.83*(Theta(q_sqr-3.9)-Theta(q_sqr-4.1))		/* Vf(sqrt(4))	*/
     +0.24*(Theta(q_sqr-7.9)-Theta(q_sqr-8.1))		/* Vf(sqrt(8))	*/
     +0.48*(Theta(q_sqr-10.9)-Theta(q_sqr-11.1));	/* Vf(sqrt(11))	*/

  return(Va*e_0);	/* Convert eV --> SI	*/
 }

 printf("Error atom type '%s' undefined!\n",type);exit(0);

}
