/*===================================================================
	efiwire Envelope Function Infinite Wire
===================================================================*/

/* This program calculates the eigenfunctions and eigenenergies of
   an infinitely deep rectangular cross-section quantum wire. The 
   relevant parameters are passed via command line arguments.

   Paul Harrison, November 1998				 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "struct.h"
#include "maths.h"
#include "const.h"

main(argc, argv)
int	argc;
char	*argv[2];
{
double	E;	/* Energy				*/
double	Ly;	/* wire width along y-axis		*/
double	Lz;	/* wire width along z-axis		*/
double	m;	/* electron efFcdctive mass		*/
double	psi_y;	/* wave function along y-direction	*/
double	psi_z;	/* wave function along z-direction	*/
double	y;	/* y coordinate				*/
double	z;	/* z coordinate				*/
int   	iy;	/* index				*/
int   	iz;	/* index				*/
int	in_y;	/* index over states			*/
int	in_z;	/* index over states			*/
int   	N;	/* total number of steps		*/
int	s;	/* number of states			*/
char	filename[9];	/* wavefunction filename	*/
char	p;	/* particle (e, h or l)			*/
FILE	*FE;	/* filepointer to energy output file	*/
FILE	*Fcd;	/* filepointer to charge density file	*/

/* default values */

Ly=100e-10;
Lz=100e-10;
N=100;
m=0.067*m0;		/* GaAs electron value	*/
p='e';
s=1;

while((argc>1)&&(argv[1][0]=='-'))
{
 switch(argv[1][1])
 {
  case 'm':
           m=atof(argv[2])*m0;
           break;
  case 'N':
           N=atoi(argv[2]);
           break;
  case 'p':
           p=*argv[2];
           switch(p)
           {
            case 'e': break;
            case 'h': break;
            case 'l': break;
            default:  printf("Usage:  efiw [-p particle (e, h, or l)]\n");
                      exit(0);
           }
           break;

  case 's':
	   s=atoi(argv[2]);
	   break;
  case 'y':
           Ly=atof(argv[2])*1e-10;
           break;
  case 'z':
           Lz=atof(argv[2])*1e-10;
           break;
  default :
	   printf("Usage:  efiwire [-m mass (\033[1m0.067\033[0mm0)][-N number of points \033[1m100\033[0m]\n");
           printf("                [-p particle (\033[1me\033[0m, h, or l)][-s # states \033[1m1\033[0m]\n");
           printf("                [-y width (\033[1m100\033[0mA)][-z width \033[1m100\033[0mA]\n");
	   exit(0);
 }
 argv++;
 argv++;
 argc--;
 argc--;
}

sprintf(filename,"E%c.r",p);
FE=fopen(filename,"w");

for(in_y=1;in_y<=s;in_y++)
 for(in_z=1;in_z<=s;in_z++)
 {
  E=sqr(pi*hbar)/(2*m)*(sqr(in_y/Ly)+sqr(in_z/Lz));
  fprintf(FE,"%i%i %24.17le\n",in_y,in_z,E/(1e-3*e_0));

  sprintf(filename,"cd%i%i.r",in_y,in_z);
  Fcd=fopen(filename,"w");

  for(iy=0;iy<N;iy++)
  {
   y=(float)iy*Ly/(float)(N-1);
   psi_y=sqrt(2/Ly)*sin(in_y*pi*y/Ly);
   for(iz=0;iz<N;iz++)
   {
    z=(float)iz*Lz/(float)(N-1);
    psi_z=sqrt(2/Lz)*sin(in_z*pi*z/Lz);
    fprintf(Fcd,"%20.17le\n",sqr(psi_y*psi_z));
   }
  }
  fclose(Fcd);

 }/* end in_z */

fclose(FE);

}
