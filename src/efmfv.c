/*=========================================================
         effv Envelope Function Magnetic Field to v
  =========================================================*/

/* This program adds the Zeeman splitting term to the 
   potential in the file `v.r'.

	Input files

		x.r
		v0.r
	or	v.r

	Output files

		v0.r
		v.r

   Paul Harrison, February 1998				 */

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
double	B_J();		/* Brillouin function			*/
double	Seff();		/* effective spin as function of x	*/
double	Teff();		/* effective temperature 		*/

double	A;
double	B;
double	DeltaV;		/* change in potential due to MF	*/
double	J;		/* magnetic ion spin			*/
double	MF;		/* Magnetic field			*/
double	N0alpha;	/* Magnetic parameter			*/
double	N0beta;		/* Magnetic parameter			*/
double	Sz;		/* effective spin			*/
double	T;		/* sample temperature			*/
double	T0;		/* effective temperature		*/
double	V;		/* CB and VB potentials			*/
double	x;		/* the alloy concentration		*/
double	z;		/* growth direction			*/
int	n;		/* return value of fopen		*/
int	i;		/* index				*/
char	p;		/* particle (e, h, or l)		*/
char	s;		/* electron spin state, up or down, +/- */
FILE	*Fv;		/* pointer to input potential file	*/
FILE	*FvF;		/* pointer to v.r with field		*/
FILE	*Fx;		/* pointer to x.r 			*/
data11	*V0;		/* pointer to potential data		*/
data11	*X0;		/* pointer to alloy data		*/

/* Define global defaults */

MF=0.0;
p='e';
s='+';

T=1.8;
N0alpha=0.220*e_0;
N0beta=0.880*e_0;
J=2.5;

while((argc>1)&&(argv[1][0]=='-'))
{
 switch(argv[1][1])
 {
  case 'B':
	   MF=atof(argv[2]);            /* read magnetic field in Tesla */
	   break;
  case 'p':
           p=*argv[2];
           switch(p)
           {
            case 'e': break;
            case 'h': break;
            case 'l': break;
            default:  printf("Usage:  efmfv [-p particle (\033[1me\033[0m, h, or l)]\n");
                      exit(0);
           }
           break;
  case 's':
	   s=*argv[2];
	   break;
  case 'T':
	   T=atof(argv[2]);
	   break;
  default: 
	   printf("Usage: efmfv [-B magnetic field (\033[1m0\033[0mTesla)][-p particle (\033[1me\033[0m, h, or l)]\n");
	   printf("             [-s spin-state (\033[1m+\033[0m/-)][-T temperature (\033[1m1.8\033[0mK)]\n");
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

/* Open structure file `x.r'	*/

if((Fx=fopen("x.r","r"))==NULL)
 {printf("Error: Cannot open file 'x.r'!\n");exit(0);}

/* Count number of lines in potential file */

n=0;
while(fscanf(Fv,"%*le %*le")!=EOF)
 n++;
rewind(Fv);

/* Allocate memory for potential and alloy concentration structures */

V0=(data11 *)calloc(n,sizeof(data11));
 if(V0==0){fprintf(stderr,"Cannot allocate memory!\n");exit(0);}

X0=(data11 *)calloc(n,sizeof(data11));
 if(X0==0){fprintf(stderr,"Cannot allocate memory!\n");exit(0);}

/* Read data into structure */

i=0;
while(fscanf(Fv,"%le %le",&((V0+i)->a),&((V0+i)->b))!=EOF)
{
 /* Note line below assumes magnetic ion in `x' column of `s.r'	*/

 fscanf(Fx,"%le %le %*le",&((X0+i)->a),&((X0+i)->b));
 i++;
}

fclose(Fv);
fclose(Fx);





/* Add field in the form +/- 3A, +/- 3B or +/- B depending on spin */

for(i=0;i<n;i++)
{
 x=(X0+i)->b;
 T0=Teff(x);
 Sz=Seff(N0alpha,N0beta,x);

 A=x*N0alpha*Sz*B_J(J,MF,T,T0)/6;
 B=x*N0beta*Sz*B_J(J,MF,T,T0)/6;

/* calculate change in potential DeltaV due to field	*/

switch(p)
{
 case 'e': 
	  switch(s)
	  {
	   case '+':DeltaV=3*A;break;
	   case '-':DeltaV=-3*A;break;
	  }
	  break;
 case 'h':
	  switch(s)
	  {
	   case '+':DeltaV=3*B;break;
	   case '-':DeltaV=-3*B;break;
	  }
	  break;
 case 'l': 
	  switch(s)
	  {
	   case '+':DeltaV=B;break;
	   case '-':DeltaV=B;break;
	  }
	  break;
}

 ((V0+i)->b)+=DeltaV;
}

/* Write data to file 'v.r' */

if((FvF=fopen("v.r","w"))==NULL)
 printf("Error: Cannot open file 'v.r' to output data!\n");

for(i=0;i<n;i++)
 fprintf(FvF,"%20.17le %20.17le\n",(V0+i)->a,(V0+i)->b);

/* Close files */

fclose(FvF);


}        /* end main */





double
Seff(N0alpha,N0beta,x)

/* This function returns the effective spin <Sz(x)>

   Data taken from 

   J.A.Gaj, C.Bodin-Deshayes, P.Peyla, J.Cibert, G.Feuillet,
   Y.Merle d'Aubigne, R.Romestain and A.Wasiela, Proc. 21st
   Int. Conf. Phys. Semiconductors 1936 (1992)               

   Note Gaj gives	Delta E(saturation)=6A+6B
  
   But Delta E(saturation)=x(N0alpha+N0beta)<Sz>

   Hence <Sz> follows from Gaj's expression	 */

double N0alpha;
double N0beta;
double x;
{
 double A = 2488.0*e_0*1e-3;
 double B = -57880.0*e_0*1e-3;
 double C = 152.7;
 double D = -20760.0*e_0*1e-3;
 double E = 8.083;
 double s;

 s=(A+B*sqr(x)/(1+C*sqr(x))+D*x/(1+E*x))/(N0alpha+N0beta);

 return(s);
}



double
Teff(x)

/* This function returns the effective temperature T0(x)

   (J.A.Gaj, C.Bodin-Deshayes, P.Peyla, J.Cibert, G.Feuillet,
   Y.Merle d'Aubigne, R.Romestain and A.Wasiela, Proc. 21st
   Int. Conf. Phys. Semiconductors 1936 (1992))               */

double x;
{
 double F = 40.7;
 double G = 4.6;

 return(F*x/(1+G*x));
}



double
B_J(J,MF,T,T0)

double J;
double MF;
double T;
double T0;
{
 double	y;
 
 y=((2*J*mu_b*MF)/(kb*(T+T0)));

 return(((2*J+1)/(2*J))*coth(((2*J+1)*y)/(2*J))-(1/(2*J))*coth((y/(2*J))));
}



