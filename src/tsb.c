/*=========================================================
          tsb Transmission coefficient Single Barrier
  =========================================================*/

/* 

		

   Paul Harrison, April 1998				*/

#include <stdio.h>
#include <stdlib.h>
#include <strings.h>
#include <math.h>
#include "struct.h"
#include "maths.h"
#include "const.h"


main(int argc,char *argv[])
{
double	dE;		/* energy step				*/
double	E;		/* energy				*/
double	k;		/* wave vector in `well' material	*/
double	kdash;		/* wave vector in barrier material	*/
double	K;		/* decay constant in barrier material	*/
double	L;		/* well width				*/
double	m;		/* effective mass 			*/
double	T;		/* transmission coefficient		*/
double	V;		/* barrier height			*/
FILE	*FT;		/* pointer to output file `T.r'		*/


/* default values, appropriate to GaAs-GaAlAs */

L=100e-10;          
m=0.067*m0;
V=100*1e-3*e_0;   

/* Computational values	*/

dE=1e-4*e_0;        /* arbitrarily small energy---0.1meV   */

 
while((argc>1)&&(argv[1][0]=='-'))
{
 switch(argv[1][1])
 {
  case 'd':
	   dE=atof(argv[2])*1e-3*e_0;
	   break;
  case 'L':
	   L=atof(argv[2])*1e-10;
	   break;
  case 'm':
	   m=atof(argv[2])*m0;
	   break;
  case 'V':
	   V=atof(argv[2])*1e-3*e_0;
	   break;
  default:
	   printf("Usage:  tsb [-d energy step (\033[1m0.1\033[0mmeV)][-L barrier width (\033[1m100\033[0mA)]\n");
	   printf("            [-m effective mass (\033[1m0.067\033[0mm0)][-V barrier height (\033[1m100\033[0mmeV)]\n");
	   exit(0);
 }
 argv++;
 argv++;
 argc--;
 argc--;
}

/* Initialise energy	*/

E=0; 

/* Open output file	*/

FT=fopen("T.r","w");

do      /* loop increments energy */
{
 k=sqrt(2*m*E)/hbar;
 K=sqrt(2*m*(V-E))/hbar;
 T=1/(1+sqr((k*k+K*K)/(2*k*K))*sqr(sinh(K*L)));
 fprintf(FT,"%20.17le %20.17le\n",E/(1e-3*e_0),T);
 E+=dE;
}while(E<V);

do      /* loop increments energy */
{
 k=sqrt(2*m*E)/hbar;
 kdash=sqrt(2*m*(E-V))/hbar;
 T=1/(1+sqr((k*k-kdash*kdash)/(2*k*kdash))*sqr(sin(kdash*L)));
 fprintf(FT,"%20.17le %20.17le\n",E/(1e-3*e_0),T);
 E+=dE;
}while(E<10*V);

fclose(FT);

}        /* end main */
