/*=========================================================
        tdb Transmission coefficient Double Barrier
  =========================================================*/

/* This code calculates the transmission coefficient T(E)
   for a double barrier structure as a function of the 
   incident energy E of a particle.  The barriers can be of 
   unequal width and the particle can have a different mass in
   the barrier material to the `well', the barrier heights are
   however equal.

   Paul Harrison, 30th April/1st May 1998		*/

#include <stdio.h>
#include <stdlib.h>
#include <strings.h>
#include <math.h>
#include "struct.h"
#include "maths.h"
#include "const.h"

main(int argc,char *argv[])
{
cmat2x2	invcmat2x2();
cmat2x2	cmat2x2mult();

double	dE;		/* energy step				*/
double	E;		/* energy				*/
double	I1,I2,I3,I4;	/* the interfaces			*/
double	k;		/* wave vector in `well' material	*/
double	kdash;		/* wave vector in barrier material	*/
double	K;		/* decay constant in barrier material	*/
double	L1;		/* left barrier width			*/
double	L2;		/* central well width			*/
double	L3;		/* right barrier width			*/
double	m_w;		/* effective mass in well		*/
double	m_b;		/* effective mass in barrier		*/
double	T;		/* transmission coefficient		*/
double	V;		/* barrier height			*/

cmat2x2 M;		/* The transfer matrix			*/
cmat2x2	M1;		/* The individual interface matrices	*/
cmat2x2	M2;
cmat2x2	M3;
cmat2x2	M4;
cmat2x2	M5;
cmat2x2	M6;
cmat2x2	M7;
cmat2x2	M8;

FILE	*FT;		/* pointer to output file `T.r'		*/


/* default values, appropriate to GaAs-GaAlAs */

L1=100e-10;          
L2=100e-10;          
L3=100e-10;          
m_w=0.067*m0;
m_b=0.067*m0;
V=100*1e-3*e_0;   

/* Computational values	*/

dE=0.01*1e-3*e_0;        /* arbitrarily small energy---0.01meV   */

 
while((argc>1)&&(argv[1][0]=='-'))
{
 switch(argv[1][1])
 {
  case 'd':
	   dE=atof(argv[2])*1e-3*e_0;
	   break;
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
	   m_w=atof(argv[2])*m0;
	   break;
  case 'n':
	   m_b=atof(argv[2])*m0;
	   break;
  case 'V':
	   V=atof(argv[2])*1e-3*e_0;
	   break;
  default:
	   printf("Usage:  tdb [-a left hand barrier width (\033[1m100\033[0mA)][-b well width (\033[1m100\033[0mA)]\n");
	   printf("            [-c right hand barrier width (\033[1m100\033[0mA)][-d energy step (\033[1m0.01\033[0mmeV)]\n");
	   printf("            [-m well effective mass (\033[1m0.067\033[0mm0)][-n barrier mass (\033[1m0.067\033[0mm0)]\n");
	   printf("            [-V barrier height (\033[1m100\033[0mmeV)]\n");
	   exit(0);
 }
 argv++;
 argv++;
 argc--;
 argc--;
}

/* Calculate interfaces	*/

I1=0;	I2=L1;	I3=L1+L2;	I4=L1+L2+L3;

/* Initialise energy	*/

E=dE; 

/* Open output file	*/

FT=fopen("T.r","w");

do      /* loop increments energy */
{
 k=sqrt(2*m_w*E)/hbar;
 K=sqrt(2*m_b*(V-E))/hbar;

 /* Define matrices */

 M1.M[0][0].re=1;			M1.M[0][0].im=0;
 M1.M[0][1].re=1;			M1.M[0][1].im=0;
 M1.M[1][0].re=0;			M1.M[1][0].im=k/m_w;
 M1.M[1][1].re=0;			M1.M[1][1].im=-k/m_w;

 M2.M[0][0].re=1;			M2.M[0][0].im=0;
 M2.M[0][1].re=1;			M2.M[0][1].im=0;
 M2.M[1][0].re=K/m_b;			M2.M[1][0].im=0;
 M2.M[1][1].re=-K/m_b;			M2.M[1][1].im=0;

 M3.M[0][0].re=exp(K*I2);		M3.M[0][0].im=0;
 M3.M[0][1].re=exp(-K*I2);		M3.M[0][1].im=0;
 M3.M[1][0].re=K*exp(K*I2)/m_b;		M3.M[1][0].im=0;
 M3.M[1][1].re=-K*exp(-K*I2)/m_b;	M3.M[1][1].im=0;

 M4.M[0][0].re=cos(k*I2);		M4.M[0][0].im=sin(k*I2);
 M4.M[0][1].re=cos(k*I2);		M4.M[0][1].im=-sin(k*I2);
 M4.M[1][0].re=-k*sin(k*I2)/m_w;	M4.M[1][0].im=k*cos(k*I2)/m_w;
 M4.M[1][1].re=k*sin(-k*I2)/m_w;	M4.M[1][1].im=-k*cos(k*I2)/m_w;

 M5.M[0][0].re=cos(k*I3);		M5.M[0][0].im=sin(k*I3);
 M5.M[0][1].re=cos(k*I3);		M5.M[0][1].im=-sin(k*I3);
 M5.M[1][0].re=-k*sin(k*I3)/m_w;	M5.M[1][0].im=k*cos(k*I3)/m_w;
 M5.M[1][1].re=k*sin(-k*I3)/m_w;	M5.M[1][1].im=-k*cos(k*I3)/m_w;

 M6.M[0][0].re=exp(K*I3);		M6.M[0][0].im=0;
 M6.M[0][1].re=exp(-K*I3);		M6.M[0][1].im=0;
 M6.M[1][0].re=K*exp(K*I3)/m_b;		M6.M[1][0].im=0;
 M6.M[1][1].re=-K*exp(-K*I3)/m_b;	M6.M[1][1].im=0;

 M7.M[0][0].re=exp(K*I4);		M7.M[0][0].im=0;
 M7.M[0][1].re=exp(-K*I4);		M7.M[0][1].im=0;
 M7.M[1][0].re=K*exp(K*I4)/m_b;		M7.M[1][0].im=0;
 M7.M[1][1].re=-K*exp(-K*I4)/m_b;	M7.M[1][1].im=0;

 M8.M[0][0].re=cos(k*I4);		M8.M[0][0].im=sin(k*I4);
 M8.M[0][1].re=cos(k*I4);		M8.M[0][1].im=-sin(k*I4);
 M8.M[1][0].re=-k*sin(k*I4)/m_w;	M8.M[1][0].im=k*cos(k*I4)/m_w;
 M8.M[1][1].re=k*sin(-k*I4)/m_w;	M8.M[1][1].im=-k*cos(k*I4)/m_w;

 M=cmat2x2mult(invcmat2x2(M1),
    cmat2x2mult(M2,
     cmat2x2mult(invcmat2x2(M3),
      cmat2x2mult(M4,
       cmat2x2mult(invcmat2x2(M5),
        cmat2x2mult(M6,
         cmat2x2mult(invcmat2x2(M7),M8)
        )
       )
      )
     )
    )
   );

 T=1/(M.M[0][0].re*M.M[0][0].re+M.M[0][0].im*M.M[0][0].im);

 fprintf(FT,"%20.17le %20.17le\n",E/(1e-3*e_0),T);
 E+=dE;
}while(E<V);

fclose(FT);

}        /* end main */



cmat2x2
cmat2x2mult(cmat2x2 M1,cmat2x2 M2)
{
 /* Multiplies two complex 2x2 matrices together	*/

 cmat2x2	M;

 M.M[0][0]=cadd(cmult(M1.M[0][0],M2.M[0][0]),cmult(M1.M[0][1],M2.M[1][0]));
 M.M[0][1]=cadd(cmult(M1.M[0][0],M2.M[0][1]),cmult(M1.M[0][1],M2.M[1][1]));
 M.M[1][0]=cadd(cmult(M1.M[1][0],M2.M[0][0]),cmult(M1.M[1][1],M2.M[1][0]));
 M.M[1][1]=cadd(cmult(M1.M[1][0],M2.M[0][1]),cmult(M1.M[1][1],M2.M[1][1]));
 
 return(M);
}



cmat2x2
invcmat2x2(cmat2x2 M)
{
 /* Calculates the inverse of a complex 2x2	*/

 complex	detcmat2x2();

 cmat2x2	Minv;
 complex 	cminus1;

 cminus1.re=-1;
 cminus1.im=0;

 Minv.M[0][0]=cdiv(M.M[1][1],detcmat2x2(M));
 Minv.M[0][1]=cdiv(cmult(cminus1,M.M[0][1]),detcmat2x2(M));
 Minv.M[1][0]=cdiv(cmult(cminus1,M.M[1][0]),detcmat2x2(M));
 Minv.M[1][1]=cdiv(M.M[0][0],detcmat2x2(M));

 return(Minv);
 
}



complex
detcmat2x2(cmat2x2 M)
{
 /* Calculates the determinant of a complex 2x2	*/

 return(csub(cmult(M.M[0][0],M.M[1][1]),cmult(M.M[0][1],M.M[1][0])));
}
