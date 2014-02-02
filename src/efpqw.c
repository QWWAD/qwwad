/*=========================================================
      efpqw 
  =========================================================*/

/* This program produces the function of alloy concentration
   x with displacement z for a parabolic quantum well,
   the output is in a format suitable for conversion into
   electron and hole potentials using efxv.      

   ----------b----------+ a +----------b---------- x_max
                        |   |
                        \___/  x_min

   Paul Harrison,  October 1995                            
   
   Major modifications, 15th May 1998                      */

#include <stdio.h>
#include <stdlib.h>
#include <strings.h>
#include <math.h>
#include "struct.h"
#include "maths.h"
#include "bools.h"
#include "const.h"
#include "io.h"


main(int argc,char *argv[])
{
double  a;             /* length a                          */
double  b;             /* end cap width                     */
double  N;             /* number of points per unit length  */
double  x;             /* alloy concentration               */
double  x_min;         /* minimum x value                   */
double  x_max;         /* maximum x value                   */
double  y=0;           /* quaternary alloy concentration
                          note not used, compatibility only */
double  z;             /* displacement                      */ 
FILE   *Fx;            /* file pointer to x versus z data   */
boolean	output_flag;   /* if set, write data to screen      */


/* default values */

a=100e-10;          
b=100e-10;
N=1e+10;
output_flag=false;
x_min=0.0;
x_max=0.100;

while((argc>1)&&(argv[1][0]=='-'))
{
 switch(argv[1][1])
 {
  case 'a':
	   a=atof(argv[2])*1e-10;
	   break;
  case 'b':
	   b=atof(argv[2])*1e-10;
	   break;
  case 'N':
	   N=atof(argv[2])*1e+10;
	   break;
  case 'o':
           output_flag=true;
           argv--;
           argc++;
           break;
  case 'x':
	   x_min=atof(argv[2]);
	   break;
  case 'y':
	   x_max=atof(argv[2]);
	   break;
  default:
	   printf("Usage:  efpqw [-a width at top of well (\033[1m100\033[0mA)][-b barrier width (\033[1m100\033[0mA)]\n");
	   printf("              [-N number of points per Angstrom \033[1m1\033[0m]\n");
	   printf("              [-o output data to screen \033[1mfalse\033[0m]\n");
	   printf("              [-x minimum alloy concentration x \033[1m0.0\033[0m]\n");
	   printf("              [-y maximum alloy concentration x \033[1m0.1\033[0m]\n");
	   exit(0);
 }
 argv++;
 argv++;
 argc--;
 argc--;
}


Fx=fopen("x.r","w");
z=0;
while(z<(2*b+a))
{
 if((z<b)||(z>b+a))
 {
  x=x_max;
 }
 else
 {
  x=x_min+sqr(z-(b+a/2))*(x_max-x_min)/sqr(a/2);
 }
 fprintf(Fx,"%20.17le %le %le\n",z,x,y);
 z+=1/N;               /* z incremented by distance between points */
}/* end while */

fclose(Fx);


}/* end main */



