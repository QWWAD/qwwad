/*======================================================================
       rlv-ss     Reciprocal Lattice Vectors for a Single Spiral
  ======================================================================

This program calculates the reciprocal lattice vectors for a
single spiral of a face-centred cubic crystal.  The maximum magnitude of
the vectors begin read in as the only command line argument.

The reciprocal lattice vector is assumed to be constructed from linear
combinations of the reciprocal lattice basis vectors, i.e.

  G=2*pi/A0*(beta1*b_1+beta2*b_2+beta3*b_3)
                   ---       ---       ---

Paul Harrison, July 1998                                                */
 
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "struct.h"
#include "maths.h"

main(int argc, char *argv[])
{
int	beta1;	/* coefficients of the reciprocal lattice vectors	*/
int	beta2;	/*  b_1, b_2 and b_3 constituting the general		*/
int	beta3;	/*  reciprocal lattice vector G				*/
int	beta1max;	/* maximum values of above			*/
int	beta2max;	/* maximum values of above			*/
int	beta3max;	/* maximum values of above			*/
double	omega;	/* storage for scalar triple products			*/
int	n_z;	/* number fcc cells along z-axis of large basis		*/
double	Gmax;	/* maximum value of |G|, e.g. [400]=>|G|=4		*/
vector	a1;	/* lattice vector		*/
vector	a2;	/* lattice vector		*/
vector	a3;	/* lattice vector		*/
vector	b1;	/* reciprocal lattice vector	*/
vector	b2;	/* reciprocal lattice vector	*/
vector	b3;	/* reciprocal lattice vector	*/
vector	c;	/* general vector		*/
FILE	*FG;	/* pointer to output file	*/

/* default values	*/

n_z=1;
Gmax=4;

while((argc>1)&&(argv[1][0]=='-'))
{
 switch(argv[1][1])
 {
  case 'g':
	   Gmax=atof(argv[2]);
           break;
  case 'z':
           n_z=atoi(argv[2]);
           break;
  default :
           printf("Usage:  rlv-ss [-g maximum value of |G| (\033[1m4\033[0m*2*pi/A0)]\n");
           printf("               [-z # fcc cells along z-axis \033[1m1\033[0m]\n");

           break;
 }
 argv++;
 argv++;
 argc--;
 argc--;
}

/* Define the Bravais lattice vectors	*/

a1.x=0.5;a1.y=0.5;a1.z=0;
a2.x=0.5;a2.y=-0.5;a2.z=0;
a3.x=0;a3.y=0;a3.z=n_z;

FG=fopen("G.r","w");

/* Calculate the three reciprocal lattice vectors	*/

c=vvprod(a2,a3);
omega=vsprod(a1,c);
b1.x=c.x/omega;
b1.y=c.y/omega;
b1.z=c.z/omega;

c=vvprod(a3,a1);
omega=vsprod(a2,c);
b2.x=c.x/omega;
b2.y=c.y/omega;
b2.z=c.z/omega;

c=vvprod(a1,a2);
omega=vsprod(a3,c);
b3.x=c.x/omega;
b3.y=c.y/omega;
b3.z=c.z/omega;

/* Calculate the maximum value of each coefficient	*/

beta1max=Nint(Gmax/vmod(b1));
beta2max=Nint(Gmax/vmod(b2));
beta3max=Nint(Gmax/vmod(b3));

for(beta1=-beta1max;beta1<=beta1max;beta1++)
 for(beta2=-beta2max;beta2<=beta2max;beta2++)
  for(beta3=-beta3max;beta3<=beta3max;beta3++)
  {
   c=vadd(vadd(vmult(b1,(double)beta1),
               vmult(b2,(double)beta2)),
               vmult(b3,(double)beta3));
   if(vmod(c)<=Gmax)
    fprintf(FG,"%f %f %f\n",c.x,c.y,c.z);
  }
fclose(FG);
}
