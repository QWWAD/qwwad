/*=========================================================
                           pth
  =========================================================*/

/* this programme calculates the confined energy levels of a 
   Poschl Teller potential hole and writes the potential to a 
   file (v.r) suitable for solution with the shooting method.

   Paul Harrison, May 1998				*/

#include <stdio.h>
#include <stdlib.h>
#include <strings.h>
#include <math.h>
#include "struct.h"
#include "maths.h"
#include "bools.h"
#include "const.h"

main(int argc,char *argv[])
{
double	alpha;		/* with parameter			*/
double	E;		/* energy				*/
double	lambda;		/* depth parameter			*/
double	L;		/* length of potential			*/
double	m;		/* effective mass			*/
double	V;		/* potential				*/
double	z;		/* displacement along growth axis	*/
int	n;		/* principal quantum number		*/
int	N;		/* number of z points per Angstrom	*/
char	p;		/* particle (e, h, or l)		*/
char    filename[9];    /* output filename                      */

FILE	*FE;		/* pointer to energy file 		*/
FILE	*FV;		/* pointer to potential file 		*/

/* default values, appropriate to GaAs-GaAlAs */

alpha=0.1*1e+10;	/* Convert A^-1->m^-1	*/
lambda=2.0;
L=300e-10;
m=0.067*m0;
N=1;
p='e';

while((argc>1)&&(argv[1][0]=='-'))
{
 switch(argv[1][1])
 {
  case 'a':
	   alpha=atof(argv[2])*1e+10;
	   break;
  case 'l':
	   lambda=atof(argv[2]);
	   break;
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
            default:  printf("Usage:  pth [-p particle (\033[1me\033[0m, h, or l)]\n");
                      exit(0);
           }
           break;
  default:
	   printf("Usage:  pth [-a width parameter alpha (\033[1m0.1\033[0m])/A\n");
           printf("            [-l depth parameter lambda \033[1m2.0\033[0m]\n");
           printf("            [-L extent of z- domain (\033[1m300\033[0mA)]\n");
	   printf("            [-m well mass (\033[1m0.067\033[0mm0)]\n");
	   printf("            [-N points/Angstrom (\033[1m1\033[0m/A)]\n");
           printf("            [-p particle (\033[1me\033[0m, h, or l)]\n");
	   exit(0);
 }
 argv++;
 argv++;
 argc--;
 argc--;
}

/* Write potential to file	*/

FV=fopen("v.r","w");

z=-L/2;		/* initialise z position	*/
do      
{
 V=-sqr(hbar*alpha)*lambda*(lambda-1)/(2*m*sqr(cosh(alpha*z)));
 fprintf(FV,"%20.17le %20.17le\n",z,V);
 z+=1/((float)N)*1e-10;
}while(z<=(L/2));

fclose(FV);

/* Write energy levels to file `Ee.r'	*/

sprintf(filename,"E%c.r",p);	/* define output filename	*/
FE=fopen(filename,"w");

n=0;
while((lambda-1-n)>=0)	
{
 E=-sqr(hbar*alpha)*sqr(lambda-1-(float)n)/(2*m);
 
 /* Write data to file, note n->n+1 to conform with other standards	*/

 fprintf(FE,"%i %20.17le\n",n+1,E/(1e-3*e_0));	

 n++;
}

fclose(FE);

}        /* end main */

