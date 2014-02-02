/*=========================================================
             effv Envelope Function Field to v
  =========================================================*/

/* This program adds an electric field to the potential in 
   the file v.r.
 
   Paul Harrison, February 1997				 */


#include <stdio.h>
#include <stdlib.h>
#include <strings.h>
#include <math.h>
#include <malloc.h>
#include "struct.h"
#include "maths.h"
#include "bools.h"
#include "const.h"
#include "io.h"

main(int argc,char *argv[])
{
double	F;		/* Electric field			*/
double	q;		/* carrier charge +/- e_0		*/
double	V;		/* CB and VB potentials			*/
double	z;		/* growth direction			*/
double	z0;		/* z-axis centre of structure		*/
int	n;		/* return value of fopen		*/
int	i;		/* index				*/
char	p;		/* particle (e, h, or l)             */
FILE	*Fv;		/* pointer to input potential file	*/
FILE	*FvF;		/* pointer to v.r with field		*/
data11	*V0;		/* pointer to potential data		*/

/* Define global defaults */

F=0.0;
q=-e_0;
p='e';

while((argc>1)&&(argv[1][0]=='-'))
{
 switch(argv[1][1])
 {
  case 'f':
	   F=atof(argv[2])*1e+5;           /* Convert kVcm^-1-->Vm^-1 */
	   break;
  case 'p':
           p=*argv[2];
           switch(p)
           {
            case 'e': q=-e_0;break;
            case 'h': q=e_0;break;
            case 'l': q=e_0;break;
            default:  printf("Usage:  effv [-p particle (\033[1me\033[0m, h, or l)]\n");
                      exit(0);
           }
           break;
  default: 
	   printf("Usage: effv [-f (\033[1m0\033[0mkV/cm)][-p particle (\033[1me\033[0m, h, or l)]\n");
	   exit(0);
 }
 argv++;
 argv++;
 argc--;
 argc--;
}

/* Open zero field potential file */

if((Fv=fopen("v0.r","r"))==0)
 {
  if((system("cp v.r v0.r"))==0) Fv=fopen("v0.r","r");
  else 
   {printf("Error: Cannot open file 'v0.r' or find file 'v.r'!\n");exit(0);}
 }

/* Count number of lines in potential file */

n=0;
while(fscanf(Fv,"%*le %*le")!=EOF)
 n++;
rewind(Fv);

/* Allocate memory */

V0=(data11 *)calloc(n,sizeof(data11));
 if(V0==0){fprintf(stderr,"Cannot allocate memory!\n");exit(0);}

/* Read data into structure */

i=0;
while(fscanf(Fv,"%le %le",&((V0+i)->a),&((V0+i)->b))!=EOF)
 i++;

fclose(Fv);

/* Calculate midpoint z0, of growth (z-) axis */

z0=(((V0+n-1)->a)-(V0->a))/2;

/* Add field in the form qF(z-z0) */

for(i=0;i<n;i++)
 ((V0+i)->b)=((V0+i)->b)+q*F*(((V0+i)->a)-z0);

/* Write data to file 'v.r' */

if((FvF=fopen("v.r","w"))==0)
 printf("Error: Cannot open file 'v.r' to output data!\n");

for(i=0;i<n;i++)
 fprintf(FvF,"%20.17le %20.17le\n",(V0+i)->a,(V0+i)->b);

/* Close files */

fclose(FvF);


}        /* end main */


