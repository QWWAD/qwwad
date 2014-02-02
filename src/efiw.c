/*===================================================================
			iw Infinite Well
===================================================================*/

/* This program calculates the eigenfunctions and eigenenergies of
   an infinite square well. The well width is passed via the
   command line.
*/

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
double	L;	/* well width				*/
double	m;	/* electron efFwfctive mass		*/
double	psi;	/* magnitude of wave function at z	*/
double	z;	/* z coordinate				*/
int   	i;	/* index				*/
int	is;	/* index over states			*/
int   	N;	/* total number of steps		*/
int	s;	/* number of states			*/
char	filename[9];	/* wavefunction filename	*/
char	p;	/* particle (e, h or l)			*/
FILE	*FE;	/* filepointer to E.r			*/
FILE	*Fwf;	/* filepointer to wf_e1.r		*/

/* default values */

L=100e-10;
N=100;
m=0.067*m0;		/* GaAs electron value	*/
p='e';
s=1;

while((argc>1)&&(argv[1][0]=='-'))
{
 switch(argv[1][1])
 {
  case 'L':
           L=atof(argv[2])*1e-10;
           break;
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
  default :
	   printf("Usage:  efiw [-L well width (\033[1m100\033[0mA)][-m mass (\033[1m0.067\033[0mm0)]\n");
	   printf("             [-N total number of points \033[1m100\033[0m]\n");
           printf("             [-p particle (\033[1me\033[0m, h, or l)]\n");
	   printf("             [-s # states \033[1m1\033[0m]\n");
	   exit(0);
 }
 argv++;
 argv++;
 argc--;
 argc--;
}

sprintf(filename,"E%c.r",p);
FE=fopen(filename,"w");

for(is=1;is<=s;is++)
{
 E=sqr(pi*hbar*is/L)/(2*m);
 fprintf(FE,"%i %24.17le\n",is,E/(1e-3*e_0));

 sprintf(filename,"wf_%c%i.r",p,is);
 Fwf=fopen(filename,"w");

 for(i=0;i<N;i++)
 {
  z=(float)i*L/(float)(N-1);
  psi=sqrt(2/L)*sin(is*pi*z/L);
  fprintf(Fwf,"%20.17le %20.17le\n",z,psi);
 }

 fclose(Fwf);

}/* end in */

fclose(FE);

}
