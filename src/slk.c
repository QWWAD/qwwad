/*=================================================================
       slk     SuperLattice k vectors
  =================================================================

   This program generates the set of k-vectors `kxi' necessary
   for the perturbative approach to the pseudopotential calculation
   of a superlattice.

   The point within the superlattice `minizone' is given in units
   of (pi/(n_z*A0), all information for which is read in from the
   command line

   Paul Harrison, October 1998		                         */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <strings.h>
#include "struct.h"

main(int argc,char *argv[])
{
float	kxi;		/* the set of electron k-points `kxi'		*/
float	g;		/* superlattice reciprocal lattice vector	*/
float	k;		/* the electron wave vector within the SL BZ	*/
int	i_kxi;		/* index over kxi				*/
int	n_z;		/* number of lattice points along z-axis of cell*/
FILE	*Fk;		/* pointer to output file			*/

/* default values	*/

k=0;
n_z=1;

while((argc>1)&&(argv[1][0]=='-'))
{
 switch(argv[1][1])
 {
  case 'k':
	   k=atof(argv[2]);
	   break;
  case 'z':
           n_z=atoi(argv[2]);
           break;
  default :
	   printf("Usage:  slk  [-k minizone k-point (\033[1m0\033[0m 2*pi/(n_zA0))][-z # cells \033[1m1\033[0m]\n");
	   exit(0);
 }
 argv++;
 argv++;
 argc--;
 argc--;
}

/* Can only now convert kbar (the k-point within the minizone) from units
   of (pi/(n_z*A0) to bulk wavevector units, i.e., units of (2*pi/A0)	*/

k/=((float)2*n_z);

/* Calculate primitive reciprocal lattice vector of the superlattice	*/

g=1/((float)n_z);

/* Open output file for writing	*/

Fk=fopen("k.r","w");

for(i_kxi=0;i_kxi<2*n_z;i_kxi++)
{
 kxi=k+(float)(-n_z+i_kxi)*g;
 fprintf(Fk,"0.0 0.0 %f\n",kxi);
}

fclose(Fk);

}/* end main */

