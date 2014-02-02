/*=========================================================
        ivdb I-V for Double Barrier
  =========================================================*/

/* This code reads in the transmission coefficient T(E) for a 
   double barrier and uses a very simple model to calculate 
   an I-V curve.

   Paul Harrison, May 1998		*/

#include <stdio.h>
#include <stdlib.h>
#include <strings.h>
#include <math.h>
#include "struct.h"
#include "maths.h"
#include "const.h"

main(int argc,char *argv[])
{
data11	*read_TofE();	/* reads T(E) v. E data into memory	*/

double	DeltaE;		/* change in energy of band minima	*/
double	dE;		/* energy integration step length	*/
double	E;		/* energy				*/
double	Ef;		/* Fermi energy				*/
double  F;              /* Electric field                       */
double	f_FD;		/* Fermi Dirac distribution function	*/
double	I;		/* current				*/
double  L1;             /* left barrier width                   */
double  L2;             /* central well width                   */
double	L3;		/* right hand barrier width		*/
double	m;		/* effective mass			*/
double	rho;		/* bulk density of states		*/
double	T;		/* temperature				*/
double	V;		/* barrier height			*/
int	iE;		/* index over E				*/
int	iF;		/* index over F				*/
int	n;		/* number of lines in `T.r' file	*/
data11	*TofE;		/* transmission coefficient		*/

FILE	*FIV;		/* pointer to output file `IV.r'	*/


/* default values, appropriate to GaAs-GaAlAs */

Ef=2*1e-3*e_0;		/* just set Ef=2meV to represent some fixed density */
F=0.0;
L1=100e-10;
L2=100e-10;
L3=100e-10;
m=0.067*m0;
T=300;

/* Computational values	*/

 
while((argc>1)&&(argv[1][0]=='-'))
{
 switch(argv[1][1])
 {
  case 'a':
           L1=atof(argv[2])*1e-10;
           break;
  case 'b':
           L2=atof(argv[2])*1e-10;
           break;
  case 'c':
           L3=atof(argv[2])*1e-10;
           break;
  case 'm':
	   m=atof(argv[2])*m0;
	   break;
  case 'T':
	   T=atof(argv[2]);
	   break;
  default:
	   printf("Usage:  ivdb [-a left hand barrier width (\033[1m100\033[0mA)][-b well width (\033[1m100\033[0mA)]\n");
	   printf("             [-c right hand barrier width (\033[1m100\033[0mA)]\n");
	   printf("             [-m effective mass (\033[1m0.067\033[0mm0)]\n");
	   printf("             [-T temperature (\033[1m300\033[0mK)]\n");
	   exit(0);
 }
 argv++;
 argv++;
 argc--;
 argc--;
}

TofE=read_TofE(&n);

FIV=fopen("IV.r","w");		/* Open output file for writing	*/

dE=((TofE+1)->a)-(TofE->a);	/* Calculate step length, assume constant */

for(iF=0;iF<100;iF++)		/* Loop for different fields	*/
{
 F=(double)iF*1e+5;		/* Convert kVcm^-1-->Vm^-1	*/
 DeltaE=e_0*F*(L1+0.5*L2);	/* Calculate DeltaE	*/
 V=F*(L1+L2+L3);		/* Calculate voltage	*/

 I=0;				/* Initialise current for integration	*/
 for(iE=0;iE<n;iE++)
 {
  E=(TofE+iE)->a;
  if(E>DeltaE)	/* only add contribution to current if E>DeltaE	*/
  {
   rho=cub(sqrt(2*m)/hbar)*sqrt(E-DeltaE)/(2*sqr(pi));
   f_FD=1/(exp((E-(Ef+DeltaE))/(kb*T))+1);
 
   I+=((TofE+iE)->b)*f_FD*rho*dE;
  /* if(iF==0)printf("%20.17le %20.17le\n",E,rho*f_FD); output f(E)rho(E) */
  }
 }

 fprintf(FIV,"%20.17le %20.17le\n",V,I);	/* output data */

} /* end loop over iE	*/

fclose(FIV);

free(TofE);

}        /* end main */




data11
*read_TofE(n)

/* This function reads the transmission coefficient versus E data
   into memory and returns the start address of this block of 
   memory and the number of lines	   */

int	*n;

{
 int	i;
 data11	*TofE;		/* pointer to start of T(E) data	*/
 FILE 	*FTofE;		/* file pointer to wavefunction file	*/

 
 /* Count number of lines	*/

 if((FTofE=fopen("T.r","r"))==0)
 {fprintf(stderr,"Error: Cannot open input file 'T.r'!\n");exit(0);}

 *n=0;	
 while(fscanf(FTofE,"%*le %*le")!=EOF)
  (*n)++;
 rewind(FTofE);

 /* Allocate memory for T versus E data */

 TofE=(data11 *)calloc(*n,sizeof(data11));
 if(TofE==0){fprintf(stderr,"Cannot allocate memory!\n");exit(0);}

 /* Read in data	*/

  i=0;
  while(fscanf(FTofE,"%le %le",&((TofE+i)->a),&((TofE+i)->b))!=EOF)
  {
   ((TofE+i)->a)*=1e-3*e_0;	/* convert meV->J	*/
   i++;
  }

 fclose(FTofE);					/* Close each file	*/

 return(TofE);
}


