/*=========================================================
                           efsqw
  =========================================================*/

/* One particle energies and wavefunctions in a single
   quantum well.  Three basic theoretical approaches are
   contained within this sourcecode

   (i) Constant mass

   -hbar^2 d^2 Y + V Y = E Y,   -lw < z < +lw,  V=0
   ------  ----                  --        --
     2m*   dz^2                   2         2

   (ii) Different masses in well and barrier

   -hbar^2 d^2 Y + V Y = E Y,   -lw < z < +lw,  V=0
   ------  ----                  --        --
     2m*   dz^2                   2         2

   with the additional constraint of the boundary conditions

   Y and 1 dY , continuous
         - --
         m dz
 
   this represents the Hamiltonian P 1 P + V
                                     -
                                     m

   (iii) Different masses in well and barrier

   -hbar^2 d^2 Y + V Y = E Y,   -lw < z < +lw,  V=0
   ------  ----                  --        --
     2m*   dz^2                   2         2

   with the additional constraint of the boundary conditions

   Y and dY , continuous
         --
         dz

   this represents the Hamiltonian 1    P P    1     + V
                                   -           -
                                 sqrt(m)     sqrt(m)

   The code is based around point (ii). Point (i) is 
   implemented by the user selecting m_b=m_w.  Point
   (iii) is implemented by allowing different m_b and m_w
   in the evaluation of k and K, but m_b is artificially
   forced to equal m_w for the boundary conditions.


   The system is solved by expressing the
   standard condition as a function f(x)=0 and
   implementing a Newton-Raphson iterative method
   where the independent variable x is the energy.
   The first estimate of the energy is found by 
   following the function along the x axis until it
   changes sign then using a midpoint rule.

   Paul Harrison 1992 

   Major modifications/simplifications

   Paul Harrison, March 1998				*/

#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <strings.h>
#include <math.h>
#include "struct.h"
#include "maths.h"
#include "const.h"

#define   N  1001      /* number of wavefunction sampling
                          points, must be odd to provide an
                          even number of strips for numerical 
                          integration.                            */


int main(int argc,char *argv[])
{
double	df_dx();	/* derivative of function wrt energy	*/
double	f();		/* function of energy			*/
void	wavef();	/* generates wave functions		*/
double	a;		/* well width				*/
double	b;		/* barrier width			*/
double	dy;		/* derivative of y			*/
double	E;		/* uncorrelated energies		*/
double	dx;		/* small energy increment		*/
double	m_w;		/* effective masses in the well		*/
double	m_b;		/* effective masses in the barrier	*/
double	m_B;		/* m_b for the boundary conditions	*/
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
bool	parity_flag;	/* true  => odd parity   
                           false => even parity			*/
bool	T_flag;		/* KE operator flag def.=1=>D(1/m)D  */


/* default values, appropriate to GaAs-GaAlAs */

a=100e-10;          
b=200e-10;
m_w=0.067*m0;
m_b=0.067*m0;
p='e';
V=100*1e-3*e_0;   
state=1;
T_flag=true;

/* Computational values	*/

dx=1e-4*e_0;        /* arbitrarily small energy---0.1meV   */

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
	   T_flag=false;
           argv--;
           argc++;
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
            default:  printf("Usage:  efsqw [-p particle (e, h, or l)]\n");
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
	   printf("Usage:  efsqw [-a well width (\033[1m100\033[0mA)][-b barrier width (\033[1m200\033[0mA)]\n");
           printf("              [-k use alternative KE operator (1/m)PP \033[1mP(1/m)P\033[0m]\n");
	   printf("              [-m well mass (\033[1m0.067\033[0mm0)][-n barrier mass (\033[1m0.067\033[0mm0)\n");
	   printf("              [-p particle (\033[1me\033[0m, h, or l)][-s # states \033[1m1\033[0m]\n");
	   printf("              [-V barrier height (\033[1m100\033[0mmeV)]\n");
	   exit(0);
 }
 argv++;
 argv++;
 argc--;
 argc--;
}

/* Setting the barrier mass equal to the well mass within the boundary 
   conditions produces T=(1/m)PP */

if(T_flag)m_B=m_b;	
else m_B=m_w;

sprintf(filename,"E%c.r",p);
FE=fopen(filename,"w");

for(i_state=1;i_state<=state;i_state++)
{

 /* deduce parity */

 if(i_state%2==1)
 {
  parity_flag=false;		/* even parity (state 1,3 etc.) */
 }
 else
 {
  parity_flag=true;		/* odd parity (state 2,4 etc.)  */
 }


  do      /* loop increments energy until function changes sign */
  {
   y1=f(a,x,m_B,m_b,m_w,V,parity_flag);
   x+=dx;
   y2=f(a,x,m_B,m_b,m_w,V,parity_flag);
  }while (y1*y2>0);

  x-=fabs(y2)/(fabs(y1)+fabs(y2))*dx;

  do      /* Newton-Raphson iterative method */
  {
   y=f(a,x,m_B,m_b,m_w,V,parity_flag);
   dy=df_dx(a,x,m_B,m_b,m_w,V,parity_flag);
   x-=y/dy;
  }while (fabs(y/dy)>1e-12*e_0);

  E=x;

  fprintf(FE,"%i %20.17le\n",i_state,E/(1e-3*e_0));

  wavef(a,b,E,m_b,m_w,V,i_state,p,parity_flag);

 x+=dx;       /* clears x from solution */

}/* end i_state */



fclose(FE);

return EXIT_SUCCESS;
}        /* end main */




/**
 * This function is the derivative of standard fw result =0
 */
double 
df_dx(a,energy,m_B,m_b,m_w,V,parity_flag)

double a;
double energy;     /* local energy */
double m_B;
double m_b;
double m_w;
double V;
bool   parity_flag;
{
 /* Find electron wavevector using Eq. 2.81, QWWAD3 */
 const double k=sqrt(2*m_w/hbar*energy/hbar);

 /* Energy-derivitive of wavevector (Eq. 2.88, QWWAD3) */
 const double dk_de=sqrt(2*m_w)/(2*hbar*sqrt(energy));

 /* Energy-derivitive of wavefunction decay constant (Eq. 2.88, QWWAD3) */
 const double dK_de=sqrt(2*m_b)/(-2*hbar*sqrt(V-energy));

 if(parity_flag) /* ODD parity */
 {
   /* df/dE for odd states (Eq. 2.87, QWWAD3) */
   return dk_de*cot(k*a/2)/m_w-k*a*sqr(cosec(k*a/2))*dk_de/(2*m_w)+dK_de/m_B;
 }
 else /* EVEN parity */
 {
   /* df/dE for even states (Eq. 2.86, QWWAD3) */
   return dk_de*tan(k*a/2)/m_w+k*a*sqr(sec(k*a/2))*dk_de/(2*m_w)-dK_de/m_B;
 }
}     



double
f(a,energy,m_B,m_b,m_w,V,parity_flag)

/* This function is the standard fw result =0 */

double a;
double energy;     /* local energy */
double m_B;
double m_b;
double m_w;
double V;
bool   parity_flag;
{
 double k;        /* electron wave vector        */
 double K;        /* wavefunction decay constant */

 k=sqrt(2*m_w/hbar*energy/hbar);
 K=sqrt(2*m_b/hbar*(V-energy)/hbar);

 if(parity_flag)
 {
  return(k*cot(k*a/2)/m_w+K/m_B);
 }
 else
 {
  return(k*tan(k*a/2)/m_w-K/m_B);
 }
}     
 


void
wavef(a,b,E,m_b,m_w,V,i_state,p,parity_flag)

/* This function calculates the uncorrelated one particle
   wavefunctions for the electron and hole and writes
   them to an external file.                              */

double a;              /* well width             */
double b;              /* barrier width          */
double E;
double m_b;
double m_w;
double V;
int    i_state;
char   p;
bool	parity_flag;
{
 double A;         /* In the well the wavefunction psi=Acoskz */
 double B;         /* and in the barrier  psi=Bexp(-Kz)       */
 double k;         /* wavevector in the well                  */
 double K;         /* decay constant in the barrier           */
 double norm_int;  /* integral over all space of psi*psi      */
 double psi[N];    /* wavefunction                            */
 double z[N];      /* displacement along growth direction     */
 int i_z;          /* z displacement index                    */
 char filename[9]; /* output filename                         */
 FILE *wf_out;

 /* Generate z co-ordinates	*/

 for (i_z=0;i_z<N;i_z++) z[i_z]=i_z*(a+2*b)/(N-1)-(b+a/2);

 /* Define k and K	*/

  k=sqrt(2*m_w/hbar*E/hbar);
  K=sqrt(2*m_b/hbar*(V-E)/hbar);

  if(parity_flag)	/* odd parity wavefunction */
  {
   B=1;              /* choose B arbitrarily                    */

   A=B*exp(-K*a/2)/sin(k*a/2);

   for (i_z=0;i_z<N;i_z++)		/* calculate wavefunctions */
   {
    if (z[i_z]<(-a/2))
    {
     psi[i_z]=-B*exp(-K*fabs(z[i_z]));		/* barrier */
    }
    if ((z[i_z]>=(-a/2))&&(z[i_z]<(a/2)))
    {
     psi[i_z]=A*sin(k*z[i_z]);			/* well    */
    }
    if (z[i_z]>=(a/2))
    {
     psi[i_z]=B*exp(-K*z[i_z]);			/* barrier */
    }
   }/* end for*/
  
    /* normalisation integral for odd parity type I */

    norm_int=sqr(A)*(a/2-sin(k*a)/(2*k))-
             sqr(B)*exp(-K*a)*(exp(-2*K*b)-1)/K;

  }
  else		/* even parity wavefunction */
  {
   B=1;              /* choose B arbitrarily                    */

   A=B*exp(-K*a/2)/cos(k*a/2);

   for (i_z=0;i_z<N;i_z++)           /* calculate wavefunctions */
   {
    if (z[i_z]<(-a/2))
    {
     psi[i_z]=B*exp(-K*fabs(z[i_z]));          /* barrier */
    }
    if ((z[i_z]>=(-a/2))&&(z[i_z]<(a/2)))
    {
     psi[i_z]=A*cos(k*z[i_z]);                 /* well    */
    }
    if (z[i_z]>=(a/2))
    {
     psi[i_z]=B*exp(-K*z[i_z]);                /* barrier */
    }
   }/* end for*/
  
    /* normalisation integral for even parity type I */

    norm_int=sqr(A)*(a/2+sin(k*a)/(2*k))+
             sqr(B)*exp(-K*a)*(1-exp(-2*K*b))/K;

   }

  /* normalise wavefunction */

  for(i_z=0;i_z<N;i_z++)psi[i_z]/=sqrt(norm_int);

 sprintf(filename,"wf_%c%i.r",p,i_state);
 wf_out=fopen(filename,"w");    /* open the external file  */

 for (i_z=0;i_z<N;i_z++)
 {
  fprintf(wf_out,"%e %e\n",z[i_z],psi[i_z]);
 }

 fclose(wf_out);              /* close the external file */
}

