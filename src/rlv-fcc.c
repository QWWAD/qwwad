/*======================================================================
        rlv-fcc  Reciprocal Lattice Vectors-Face Centred Cubic
  ======================================================================

This program calculates the reciprocal lattice vectors for a
face-centred cubic crystal.  The maximum magnitude of the vectors
begin read in as the only command line argument.

The reciprocal lattice vector is assumed to be constructed from linear
combinations of the reciprocal lattice basis vectors, i.e.

  G=2*pi/A0*(beta_1*b_1+beta_2*b_2+beta_3*b_3)
                    ---        ---        ---

then

  G=2*pi/A0*(gamma_1*i+gamma_2*j+gamma_3*k)
                      -         -         -

where i, j and k are the cartesian basis vectors.  
      -  -     -

Paul Harrison, March 1995                                                */
 
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "struct.h"
#include "maths.h"

main(int argc, char *argv[])
{
int beta_1;          /* coefficients of the reciprocal lattice vectors */
int beta_2;          /*  b_1, b_2 and b_3 constituting the general     */
int beta_3;          /*  reciprocal lattice vector G                   */
float gamma_1;
float gamma_2;
float gamma_3;
double G_max;           /* maximum value of |G|, e.g. [400]=>|G|=4        */
FILE *FG;

/* default values	*/

G_max=4;

while((argc>1)&&(argv[1][0]=='-'))
{
 switch(argv[1][1])
 {
  case 'g':
	   G_max=atof(argv[2]);
           break;
  default :
           printf("Usage:  rlv-fcc [-g maximum value of |G|]\n");
           break;
 }
 argv++;
 argv++;
 argc--;
 argc--;
}





FG=fopen("G.r","w");

for(beta_1=-G_max;beta_1<=G_max;beta_1++)
 for(beta_2=-G_max;beta_2<=G_max;beta_2++)
  for(beta_3=-G_max;beta_3<=G_max;beta_3++)
  {
   gamma_1=(float)(beta_1-beta_2+beta_3);
   gamma_2=(float)(beta_1+beta_2-beta_3);
   gamma_3=(float)(-beta_1+beta_2+beta_3);
   if(Nint(gamma_1*gamma_1+gamma_2*gamma_2+gamma_3*gamma_3)<=G_max*G_max)
    fprintf(FG,"%f %f %f\n",gamma_1,gamma_2,gamma_3);
  }
fclose(FG);
}
