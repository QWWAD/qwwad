/*=========================================================
     efkpsl  Envelope Function Kronig-Penney SuperLattice
  =========================================================*/

/* This program calculates the energy levels for a user 
   supplied wave vector along the growth (z-) axis of an 
   infinite Kronig-Penney superlattice.
   The code is based around:

   Different masses in well and barrier

   -hbar^2 d^2 Y + V Y = E Y,   0<z<a, V=0, a<z<b, V=V
   ------  ----                 
     2m*   dz^2                

   with the additional constraint of the boundary conditions

   Y and 1 dY , continuous
         - --
         m dz
 
   this represents the Hamiltonian P 1 P + V
                                     -
                                     m

   The system is solved by expressing the
   standard condition as a function f(x)=0 and
   implementing a Newton-Raphson iterative method
   where the independent variable x is the energy.
   The first estimate of the energy is found by 
   following the function along the x axis until it
   changes sign then using a midpoint rule.

   Paul Harrison, April 1998				*/

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
double	f();		/* function of energy			*/
double	a;		/* well width				*/
double	b;		/* barrier width			*/
double	dy;		/* derivative of y			*/
double	E;		/* uncorrelated energies		*/
double	dx;		/* small energy increment		*/
double	k;		/* wave vector				*/
double	m_w;		/* effective masses in the well		*/
double	m_b;		/* effective masses in the barrier	*/
double	V;		/* potential				*/
double	x;		/* energy estimate			*/
double	y;		/* current functional value		*/
double	y1;		/* temporary y value			*/
double	y2;		/*     "     "   " 			*/
int	i_state;	/* principal quantum number index	*/
int	state;		/* number of states			*/
char	p;		/* particle				*/
char	filename[9];	/* output filename			*/

FILE	*FE;		/* pointer to energy file 		*/


/* default values, appropriate to GaAs-GaAlAs */

a=100e-10;          
b=100e-10;
k=0;
m_w=0.067*m0;
m_b=0.067*m0;
p='e';
V=100*1e-3*e_0;   
state=1;

/* Computational values	*/

dx=1e-3*e_0;        /* arbitrarily small energy---0.1meV   */

x=dx;    /* first energy estimate */
 
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
  case 'k':
	   k=atof(argv[2]);
           break;
  case 'm':
	   m_w=atof(argv[2])*m0;
	   break;
  case 'n':
	   m_b=atof(argv[2])*m0;
	   break;
  case 'p':
           p=*argv[2];
           switch(p)
           {
            case 'e': break;
            case 'h': break;
            case 'l': break;
            default:  printf("Usage:  efkpsl [-p particle (\033[1me\033[0m, h, or l)]\n");
                      exit(0);
           }
           break;
  case 's':
	   state=atoi(argv[2]);
	   break;
  case 'V':
	   V=atof(argv[2])*1e-3*e_0;
	   break;
  default:
	   printf("Usage:  efkpsl [-a well width (\033[1m100\033[0mA)][-b barrier width (\033[1m100\033[0mA)]\n");
           printf("               [-k wave vector (\033[1m0\033[0m pi/L)\n");
	   printf("               [-m well mass (\033[1m0.067\033[0mm0)][-n barrier mass (\033[1m0.067\033[0mm0)\n");
	   printf("               [-p particle (\033[1me\033[0m, h, or l)][-s # states \033[1m1\033[0m]\n");
	   printf("               [-V barrier height (\033[1m100\033[0mmeV)]\n");
	   exit(0);
 }
 argv++;
 argv++;
 argc--;
 argc--;
}

/* Need to convert k into SI now sure of L (=a+b)	*/

k*=pi/(a+b);

sprintf(filename,"E%c.r",p);
FE=fopen(filename,"w");

for(i_state=1;i_state<=state;i_state++)
{
 do      /* loop increments energy until function changes sign */
 {
  y1=f(a,b,x,k,m_b,m_w,V);
  x+=dx;
  y2=f(a,b,x,k,m_b,m_w,V);
 }while((y1*y2>0)&&(x<10*V));

 /* Note the end stop to prevent infinite loop in absence of solution	*/ 

 x-=fabs(y2)/(fabs(y1)+fabs(y2))*dx;

 do      /* Newton-Raphson iterative method */
 {
  y=f(a,b,x,k,m_b,m_w,V);
  dy=(f(a,b,x+dx/100,k,m_b,m_w,V)-f(a,b,x-dx/100,k,m_b,m_w,V))/(2*dx/100);
  x-=y/dy;
 }while(fabs(y/dy)>1e-9*e_0);	
 

 E=x;

 fprintf(FE,"%i %20.17le\n",i_state,E/(1e-3*e_0));

 x+=dx;       /* clears x from solution */

}/* end i_state */



fclose(FE);

}        /* end main */



double
f(a,b,energy,k,m_b,m_w,V)

/* This function is the standard fw result =0 */

double a;
double b;
double energy;     /* local energy */
double k;
double m_b;
double m_w;
double V;
{
 double F;	    /* value of function	   */
 double k_w;        /* wave vector in well         */
 double k_b;        /* wave vector in barrier      */
 double K;	    /* wave function decay const.  */


 if(energy<V)
 {
  k_w=sqrt(2*m_w/hbar*energy/hbar);
  K=sqrt(2*m_b/hbar*(V-energy)/hbar);

  F=cos(k_w*a)*cosh(K*b)-sin(k_w*a)*sinh(K*b)*
        (sqr(m_b*k_w)-sqr(m_w*K))/(2*m_w*m_b*k_w*K)-cos(k*(a+b));
 }
 else
 {
  k_w=sqrt(2*m_w/hbar*energy/hbar);
  k_b=sqrt(2*m_b/hbar*(energy-V)/hbar);


  F=cos(k_w*a)*cos(k_b*b)-sin(k_w*a)*sin(k_b*b)*
        (sqr(m_b*k_w)-sqr(m_w*k_b))/(2*m_w*m_b*k_w*k_b)-cos(k*(a+b));
 }

 return(F);
}     
