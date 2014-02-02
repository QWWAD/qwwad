/*======================================================================
       rlv-sc Reciprocal Lattice Vectors of a Simple Cube
  ======================================================================

This program calculates the reciprocal lattice vectors for a
simple cubic crystal, with a basis consisting of n_x*n_y*n_z fcc
unit cells.  The maximum magnitude of the vectors begin read in as 
a command line argument.

The reciprocal lattice vector is assumed to be constructed from linear
combinations of the reciprocal lattice basis vectors, i.e.

  G=beta_1*b_1+beta_2*b_2+beta_3*b_3
           ---        ---        ---

where b_1= 2*pi  i   b_2= 2*pi  j , etc.
      --- ------ - ' --- ------ -
	  n_x*A0         n_y*A0

where i, j and k are the cartesian basis vectors.  
      -  -     -

Paul Harrison, June 1996                                

Paul Harrison, Modifications 1998	               */
 
#include <stdio.h>

main(int argc, char *argv[])
{
double	beta_1;		/* coefficients of the reciprocal lattice vectors */
double	beta_2;		/*  b_1, b_2 and b_3 constituting the general     */
double	beta_3;		/*  reciprocal lattice vector G                   */
int	n_x;		/* number fcc cells along x-axis of large basis   */
int	n_y;		/* number fcc cells along y-axis of large basis   */
int	n_z;		/* number fcc cells along z-axis of large basis   */
int	G_max;		/* maximum value of |G|, e.g. [400]=>|G|=4        */
FILE	*FG;

/* default values	*/

G_max=4;
n_x=1;n_y=1;n_z=1;

while((argc>1)&&(argv[1][0]=='-'))
{
 switch(argv[1][1])
 {
  case 'g':
	   G_max=atoi(argv[2]);;
	   break;
  case 'x':
	   n_x=atoi(argv[2]);
           break;
  case 'y':
	   n_y=atoi(argv[2]);
           break;
  case 'z':
	   n_z=atoi(argv[2]);
           break;
  default :
	   printf("Usage: rlv-sc [-g maximum value of |G| \033[1m4\033[0m]\n");
	   printf("              [-x # fcc cells along x-axis \033[1m1\033[0m][-y # \033[1m1\033[0m][-z # \033[1m1\033[0m]\n");
           break;
 }
 argv++;
 argv++;
 argc--;
 argc--;
}

FG=fopen("G.r","w");

for(beta_1=-G_max;beta_1<=G_max;beta_1+=1/(double)n_x)
 for(beta_2=-G_max;beta_2<=G_max;beta_2+=1/(double)n_y)
  for(beta_3=-G_max;beta_3<=G_max;beta_3+=1/(double)n_z)
  {
   if(beta_1*beta_1+beta_2*beta_2+beta_3*beta_3<=(double)(G_max*G_max))
    fprintf(FG,"%20.17lf %20.17lf %20.17lf\n",beta_1,beta_2,beta_3);
  }
fclose(FG);
}
