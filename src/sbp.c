/*==================================================================
      sbp	SubBand Populations
  ==================================================================

   This program generates the Fermi-Dirac distribution function
   for an individual subband given its population and the lattice
   temperature.

	Input files:
			Ep.r	Subband energies file, p=e,h,l

	Output files:
			FDX.r	F-D distribution for subband X

   Paul Harrison, March 1997                                   */

#include <stdio.h>
#include <stdlib.h>
#include <strings.h>
#include <math.h>
#include <signal.h>
#include <malloc.h>
#include "struct.h"
#include "const.h"
#include "maths.h"
#include "bools.h"

typedef
struct	{
 double	E;		    /* total electron energy		  */
 double	SR;		    /* scattering rate            	  */
} data;

main(int argc,char *argv[])
{
double	calc_fermilevel();	/* calculates Fermi level		*/
double	*read_energies();	/* reads subband energies from Ep.r	*/
double	*read_populations();	/* reads subband populations from N.r	*/
void	calc_dist();		/* calculates electron distribution	*/

double	*E;			/* pointer to energy data		*/
double	Ef;			/* Fermi energy 			*/
double	m;			/* effective mass			*/
double	*N;			/* number of electrons/area 		*/
double	T;			/* temperature				*/
int	n;			/* number of lines in `filename'	*/
int	nE;			/* number of energies to output FD	*/
int	s;			/* index over subband states		*/
char	p;			/* particle, e, h or l			*/
boolean	FD_flag;		/* if true prints out all FD distrib	*/
FILE	*FEf;			/* file pointer to Fermi Energy file	*/

/* default values */

FD_flag=false;
m=0.067*m0;		/* GaAs electron effective mass		*/
nE=1000;
p='e';
T=300;

while((argc>1)&&(argv[1][0]=='-'))
{
 switch(argv[1][1])
 {
  case 'f':
           FD_flag=true;
           argv--;
           argc++;
           break;
  case 'm':
	   m=atof(argv[2])*m0;
	   break;
  case 'n':
	   nE=atoi(argv[2]);
	   break;
  case 'p':
           p=*argv[2];
           switch(p)
           {
            case 'e': break;
            case 'h': break;
            case 'l': break;
            default:  printf("Usage:  sbp [-p particle (\033[1me\033[0m, h, or l)]\n");
                      exit(0);
           }
           break;
  case 'T':
	   T=atof(argv[2]);
	   break;
  default :
	   printf("Usage:  sbp [-f (output distributions \033[1mfalse\033[0m)][-m mass (\033[1m0.067\033[0mm0)]\n");
	   printf("            [-n (number of output energies \033[1m1000\033[0m)][-p particle (\033[1me\033[0m, h, or l)]\n");
	   printf("            [-T temperature (\033[1m300\033[0mK)]\n");
	   exit(0);
 }
 argv++;
 argv++;
 argc--;
 argc--;
}

E=read_energies(p,&n);		/* reads subband energy file	*/

N=read_populations(n);		/* reads subband populations file	*/

if((FEf=fopen("Ef.r","w"))==0)
 {fprintf(stderr,"Error: Cannot open input file 'Ef.r'!\n");exit(0);}

for(s=0;s<n;s++)		/* s=0 => ground state		*/
{
 Ef=calc_fermilevel(*(E+s),m,*(N+s),T);

 fprintf(FEf,"%i %20.17le\n",s+1,Ef/(1e-3*e_0));

 if(FD_flag) calc_dist(*(E+s),Ef,m,T,nE,s);
}

fclose(FEf);

free(E);


} /* end main */



void
calc_dist(Emin,Ef,m,T,nE,s)

/* This function calculates the probability of occupation of the 
   subband energies.

                                                                        */
 

double	Emin;			/* subband minima			*/
double  Ef;			/* Fermi energy                         */
double  m;                      /* effective mass                       */
double  T;                      /* temperature                          */
int	nE;			/* number of energies to output FD	*/
int	s;			/* number of subband			*/
{
double	Vmax(); 

double	dE;			/* energy increment in integration	*/
double	E;			/* Energy				*/
double	Emax;
double	f;			/* probability of occupation		*/
double	Ne=0;			/* total electron density (for checking)*/
double	vmax;			/* maximum value of potential		*/
int     i;			/* general index			*/
char	filename[9];		/* output filename for FD distribs.	*/
FILE	*FFD;			/* file pointer to Fermi-Dirac probability
				   of occupation of subband file	*/

sprintf(filename,"FD%i.r",s+1);
FFD=fopen(filename,"w");

vmax=Vmax();
Emax=Ef+10*kb*T;if(Emax>vmax) Emax=vmax;

/* Implement trapezium rule integration, i.e. `ends+2middles' 	*/

f=1/(exp((Emin-Ef)/(kb*T))+1);
fprintf(FFD,"%20.17le %20.17le\n",Emin/(1e-3*e_0),f);
Ne+=f;

dE=(Emax-Emin)/(nE-1);
E=Emin;
for(i=1;i<nE-1;i++)
{
 E+=dE;
 f=1/(exp((E-Ef)/(kb*T))+1);
 fprintf(FFD,"%20.17le %20.17le\n",E/(1e-3*e_0),f);
 Ne+=2*f;
}

f=1/(exp((Emax-Ef)/(kb*T))+1);
fprintf(FFD,"%20.17le %20.17le\n",Emax/(1e-3*e_0),f);
Ne+=f;

Ne*=m/(pi*sqr(hbar))*0.5*dE;

printf("Ne=%20.17le\n",Ne/1e+14);

fclose(FFD);
}



double
calc_fermilevel(E,m,N,T)
 
/* Solutions are sought, using a stepping algorithm, a midpoint rule and 
   finally a Newton-Raphson method, to the equation f(E_F)=0.  This 
   function has been derived by integrating the total density of states
   multiplied by the Fermi-Dirac distribution function across the in-plane
   energy---thus giving the total electron density, i.e.

         oo
   Ne=  I  f(E)N(E) dE
         Eo
   
   where f(E) is the normal Fermi-Dirac distribution function and
   N(E)=m/(pi*sqr(hbar)) is the areal density of states in a QW, see
   Bastard p12.
									*/
   

double	E;			/* subband minima 			*/
double	m;			/* effective mass			*/
double	N;			/* number of electrons/unit area	*/
double	T;			/* temperature				*/
{
double  f();			/* function to be solved		*/
double	Vmax();			/* maximum value of potential		*/

double	delta_E=0.001*1e-3*e_0;	/* energy increment			*/
double 	Emin;			/* subband minimum			*/
double 	Emax;			/* subband maximum			*/
double	vmax;			/* potential maximum, i.e. top of well	*/
double	x;			/* independent variable			*/
double	y1;			/* dependent variable			*/
double	y2;			/* dependent variable			*/

Emin=E;				/* subband minimum 			*/

vmax=Vmax();			/* calulate potential maximum		*/

x=Emin-20*kb*T;			/* first value of x			*/

/* In this implementation, the upper limit of integration is set at the 
   Fermi level+10kT, limited at potential maximum			*/

Emax=x+10*kb*T;if(Emax>vmax) Emax=vmax;

y2=f(x,Emax,Emin,m,N,T);

do
{
 y1=y2;
 x+=delta_E;
 Emax=x+10*kb*T;if(Emax>vmax) Emax=vmax;
 y2=f(x,Emax,Emin,m,N,T);
}while(y1*y2>0);

/* improve estimate using midpoint rule */

x-=fabs(y2)/(fabs(y1)+fabs(y2))*delta_E;

return(x);
}



double
f(E_F,Emax,Emin,m,N,T)

/* Function to be solved						
									*/

double 	E_F;			/* Fermi energy				*/
double 	Emax;			/* subband maximum (top of QW)		*/
double 	Emin;			/* subband minimum			*/
double	m;			/* effective mass			*/
double	N;			/* number of electrons/unit area	*/
double	T;			/* temperature				*/
{
double	y;			/* dependent variable			*/

y=m/(pi*hbar)*(kb*T/hbar)*
  (
   ((Emax-E_F)/(kb*T)-log(1+exp((Emax-E_F)/(kb*T))))-
   ((Emin-E_F)/(kb*T)-log(1+exp((Emin-E_F)/(kb*T))))
  )
  -N;

return(y);
}



double
*read_energies(p,n)

/* This function reads the potential into memory and returns the start
   address of this block of memory and the number of lines	   */

char	p;
int	*n;

{
 double	*E;
 int	i=0;		/* index over the energies			*/
 char	filename[9];	/* filename string				*/
 FILE 	*FE;		/* file pointer to energy data 			*/
 data	*fdata;		/* temporary pointer to data			*/
 data	*data_start;	/* start address of data			*/


 sprintf(filename,"E%c.r",p);
 if((FE=fopen(filename,"r"))==0)
 {
   fprintf(stderr,"Error: Cannot open input file '%s'!\n",filename);
   exit(0);
 }

 *n=0;
 while(fscanf(FE,"%*i %*e")!=EOF)
  (*n)++;
 rewind(FE);

 E=(double *)calloc(*n,sizeof(double));
 if (E==0)  {
  fprintf(stderr,"Cannot allocate memory!\n");
  exit(0);
 }

 while(fscanf(FE,"%*i %le",E+i)!=EOF)
 {
  *(E+i)*=1e-3*e_0;		/*convert meV->J		*/
  i++;
 }

 fclose(FE);

 return(E);

}



double
*read_populations(n)

/* This function reads the populations into memory and returns the start
   address of this block of memory and the number of lines	   */

int	n;

{
 double	*N;
 int	i=0;		/* index over the energies			*/
 int	m;		/* counter over number of populations defined	*/
 FILE 	*FN;		/* file pointer to energy data 			*/

 if((FN=fopen("N.r","r"))==0)
 {
   fprintf(stderr,"Error: Cannot open input file 'N.r'!\n");
   exit(0);
 }

 m=0;
 while(fscanf(FN,"%*i %*f")!=EOF)
  m++;
 rewind(FN);

 if(m!=n)
  {printf("Subband populations not defined for all levels!\n");exit(0);}

 N=(double *)calloc(n,sizeof(double));
 if (N==0)  {
  fprintf(stderr,"Cannot allocate memory!\n");
  exit(0);
 }

 while(fscanf(FN,"%*i %lf",N+i)!=EOF)
 {
  *(N+i)*=1e+10*1e+4;	/*convert from units of 10^10cm^-2->m^-2 */
  i++;
 }

 fclose(FN);

 return(N);

}



double
Vmax()

/* This function scans the file v.r and returns the maximum value of the
   potential.
									*/
{
 double	max;			/* maximum value of potential energy	*/
 double	v;			/* potential				*/
 FILE	*Fv;			/* file pointer to v.r			*/

max=0;

if((Fv=fopen("v.r","r"))==0)
 {fprintf(stderr,"Error: Cannot open input file 'v.r'!\n");exit(0);}

while(fscanf(Fv,"%*e %le",&v)!=EOF)
 if(v>max) max=v;

fclose(Fv);

return(max);

}
