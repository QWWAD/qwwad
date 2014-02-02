/*=========================================================
       sfr    Splin-Flip Raman 
  =========================================================*/

/*
   This programme calculates a hypothetical spin-flip Raman
   spectrum assuming a uniform distribution of dopants at 
   positions r_d and assigning a Gaussian intensity curve, of
   FWHM=linewidth, to each one.  The total signal is therefore
   a simple sum of Gaussians.

   Note the note about the form of the input file `e_sf.r' below.
   Note also the shoddy use of oversized arrays, sorry.

   Breaks the cardinal sin of dealing with energy and linewidths
   in meV and cm^-1 (!).   This is Raman spectroscopy speak!

		Input files

		e_sf.r		Spin-flip energy, best produced
				by splining a data from a donor
				calculation with say 20 donor
				positions, to several hundred
				points with a simple graphics
				programme

		Output files

		I.r		Raman signal intensity versus energy


  Paul Harrison.

  Substantially modified February 1998.

                                                                 */

#include <stdio.h>
#include <stdlib.h>
#include <strings.h>
#include <math.h>
#include "struct.h"
#include "maths.h"
#include "bools.h"
#include "const.h"
#include "io.h"

#define N 1000

main(int argc,char *argv[])
{
double	E;		/* spectral energy			*/
double	E_sf[N];	/* spin-flip energy for each r_d	*/
double	E_min;		/* lower limit of spectrum		*/
double	E_max;		/* upper limit of spectrum		*/
double	E_step;		/* energy increment along spectrum	*/
double	I;		/* intensity of Raman signal at E	*/
double	linewidth;	/* linewidth of signal from each r_d	*/
double	r_d[N];		/* donor positions (r_d)		*/
double	sigma;		/* standard deviation of Gaussians	*/
int	i_i;		/* index over r_d for I sum		*/
int	i_sf;		/* index over r_d for spin-flip energies*/
int	N_rd;		/* number of r_d points in e_sf.r	*/
FILE	*fE_sf;		/* file pointer to input file e_sf.r	*/
FILE	*fI;		/* file pointer to output file I.r	*/


/* default values */

E_min=0;
E_max=100;
E_step=1.0;
linewidth=1.0;

while((argc>1)&&(argv[1][0]=='-'))
{
 switch(argv[1][1])
 {
  case 'l':
           linewidth=atof(argv[2]);
           break;
  case 's':
           E_min=atof(argv[2]);
           break;
  case 't':
           E_step=atof(argv[2]);
           break;
  case 'u':
           E_max=atof(argv[2]);
           break;
  default:
	   printf("Usage:  sfr [-l linewidth (\033[1m1\033[0mcm^-1)]\n");
	   printf("            [-s lower limit of Raman spectra (\033[1m0\033[0mcm^-1)]\n");
	   printf("            [-t increment (\033[1m1\033[0mcm^-1)][-u upper limit (\033[1m100\033[0mcm^-1)]\n");
	   exit(0);
 }
 argv++;
 argv++;
 argc--;
 argc--;
}



/* Read spin-flip energies in from spline file
   in units of wavenumbers                              */

if((fE_sf=fopen("e_sf.r","r"))==0)
{fprintf(stderr,"Error: Cannot open input file 'e_sf.r'!\n");
 fprintf(stderr,"Generate file by saving spline from `e.r+' - `e.r-'\n");
 exit(0);
}

i_sf=0;
while(fscanf(fE_sf,"%lf %lf",&r_d[i_sf],&E_sf[i_sf])!=EOF)
{
 E_sf[i_sf]*=1e-3*e_0/(h*c0)/100;	/* convert from meV to cm^-1	*/
 i_sf++;
}

fclose(fE_sf);

N_rd=i_sf;				/* Assign final index to count	*/

/* Generate standard deviation from linewidths (FWHM)	*/

 sigma=linewidth/(2*sqrt(2*log(2)));
  
/* Calculate intensity versus energy of spectra, using a
     sum of Gaussians.                                     */

fI=fopen("I.r","w");

E=E_min;
do
{
 I=0;
 for(i_i=0;i_i<N_rd;i_i++)
 {
  I+=1/(sigma*sqrt(2*pi))*exp(-0.5*sqr((E-E_sf[i_i])/sigma));
 }
 fprintf(fI,"%le %le\n",E,I); /* E in cm^-1, I in arb. units */
 E+=E_step;
}while(E<E_max);

fclose(fI);



}        /* end main */

