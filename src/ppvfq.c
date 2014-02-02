/*=================================================================
              ppvfq   Pseudo-Potential form factor (Vf) of Q 
  =================================================================

   This program generates the pseudopotential as a function of 
   the reciprocal lattice vector q=G'-G, and stores in a file for
   viewing.
 
   Input files:

   Output files:
   			Vfq`TYPE'.r	Vf(q) in eV versus q in (2pi/A0)

   Paul Harrison, July 1998
								*/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "struct.h"
#include "maths.h"
#include "const.h"

#include "ppff.c"

main(int argc,char *argv[])
{
extern double	Vf();	/* the form factor function			*/
double	A0;		/* Lattice constant				*/
double	G_max;		/* maximum reciprocal lattice vector		*/
double	m_per_au;	/* unit conversion factor, m/a.u.		*/
double	q;		/* q=G'-G					*/
double	v;		/* the potential				*/
int	N;		/* number of points				*/
int	iq;		/* loop index for matrix rows			*/
char	filename[9];	/* character string for output filename		*/
char	type[12];	/* atom type, conforms to standard in ppff.c	*/
FILE	*FVfq;		/* pointer to output file			*/

/* default values	*/

A0=5.65e-10;
sprintf(type,"SI");
G_max=4;
N=100;

/* computational default	*/

m_per_au=4*pi*epsilon_0*sqr(hbar/e_0)/m0;

while((argc>1)&&(argv[1][0]=='-'))
{
 switch(argv[1][1])
 {
  case 'A':
	   A0=atof(argv[2])*1e-10;
           break;
  case 'a':
	   sprintf(type,argv[2]);
	   break;  
  case 'g':
           G_max=atoi(argv[2]);
           break;
  case 'N':
	   N=atoi(argv[2]);
	   break;
  default :
	   printf("Usage:  ppvfq [-A lattice constant (\033[1m5.65\033[0mA)][-a atom type \033[1mSI\033[0m]\n");
	   printf("              [-g maximum value of |G| (\033[1m4\033[0m*2*pi/A0)]\n");
	   printf("              [-N number of q points \033[1m100\033[0m]\n");
	   exit(0);
 }
 argv++;
 argv++;
 argc--;
 argc--;
}

/* Define and open file for output	*/

sprintf(filename,"Vfq%s.r",type);
FVfq=fopen(filename,"w");

for(iq=0;iq<N;iq++)  
{
 q=(float)iq*G_max*2*pi/(A0*(float)N);

 v=Vf(A0,m_per_au,q*q,type);	/* calculate form factor---a function of q^2	*/

 fprintf(FVfq,"%f %le\n",q/(2*pi/A0),v/e_0);	/* write to file	*/
}

/* Close output file	*/

fclose(FVfq);
 
}/* end main */
