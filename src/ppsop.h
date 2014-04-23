/*=================================================================
       ppsop     PseudoPotential Spin-Orbit Parameters
  =================================================================

   This function returns the spin-orbit parameter (lambda)
   for a variety of atomic species.

   Paul Harrison, April 2000					 */


double
lambda(type)

/* This function returns the spin-orbit parameter lambda for the 
   appropriate atomic species */

char	type[];		/* atomic species			*/
{

 if(!strcmp(type,"SI"))  return(0.000106*h*c*Rinf);
 if(!strcmp(type,"GE"))  return(0.00058*h*c*Rinf);
 if(!strcmp(type,"GAASmz"))return(0.000402*h*c*Rinf);
 if(!strcmp(type,"ASGAmz"))return(0.000402*1.38*h*c*Rinf);
 if(!strcmp(type,"CDTE"))return(0.055*(0.0250-0.0090)*h*c*Rinf);
 if(!strcmp(type,"TECD"))return(0.055*(0.0250+0.0090)*h*c*Rinf);
 if(!strcmp(type,"CDTEcb"))return(0.343*(0.002-0.0002)*h*c*Rinf);
 if(!strcmp(type,"TECDcb"))return(0.343*(0.002+0.0002)*h*c*Rinf);

 printf("Error atom type '%s' undefined in spin-orbit parameter set!\n",type);
 exit(0);

}
