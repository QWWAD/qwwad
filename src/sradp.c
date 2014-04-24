/*==================================================================
           sradp Scattering Rate Acoustic Deformation Potential
  ==================================================================*/

/* This program calculates the carrier-acoustic deformation potential
   (acoustic phonon) scattering rate for parallel parabolic
   intra- and intersubband events.  The required rates are provided
   by the user in the file `rrp.r'.  The other necessary inputs are listed
   below.


	Input files:		rrp.r	contains required rates
				wf_xy.r	x=particle y=state
				Ex.r	x=particle, energies

	Output files:		AC[a,e]if.r	absorption and emission rates


    Paul Harrison, March 2000 
    									*/
#include <complex.h> 
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <signal.h>
#include "struct.h"
#include "qclsim-constants.h"
#include "maths.h"


int main(int argc,char *argv[])
{
double	Gsqr();		/* form factor calculation			*/
double	*read_E();	/* reads subband minima energies from file	*/
double	Vmax();		/* maximum value of the potential		*/
void    ff_output();    /* output formfactors                           */
data11	*ff_table();
data12	*read_wf();	/* reads wavefunctions into memory		*/

double	A0;		/* lattice constant				*/
double	alpha1;		/* solutions for phonon wavevector Kz		*/
double	alpha2;		/* solutions for phonon wavevector Kz		*/
double	arg;		/* argument of square root function		*/
double	Da;		/* acoustic deformation potential		*/
double	DeltaE;		/* Ef-Ei, see notes				*/
double	delta_z;	/* step length along growth (z-) axis		*/
double  dki;            /* step length for loop over ki                 */
double	dKz;		/* step length for integration over Kz		*/
double	dtheta;		/* step length for integration over theta	*/
double	*E;		/* pointer to subband minima			*/
double	Ephonon;	/* phonon energy				*/
double  kimax;          /* maximum value of ki                          */
double  ki;             /* carrier momentum (wave vector actually)	*/
double	Kz;		/* phonon wavevector				*/
double	m;		/* carrier effective mass			*/
double	N0;		/* number of phonons at LO phonon energy	*/
double	omega_0;	/* angular frequency of LO phonon		*/
double	rho;		/* density of material				*/
double	T;		/* temperature					*/
double	theta;		/* the angle theta				*/
double	Upsilon_a;	/* absorption scattering rate prefactor		*/
double	Upsilon_e;	/* emission scattering rate prefactor		*/
double	Vs;		/* velocity of sound				*/
double	Wif;		/* temporary sum (for integration)		*/
double	Waif;		/* the scattering rate				*/
double	Weif;		/* the scattering rate				*/
int	iKz;		/* index over Kz				*/
int     iki;            /* index over ki                                */
int	itheta;		/* index over theta 				*/
int	state[2];	/* electron state index				*/
int	n;		/* length of wavefunctions file			*/
int	nE;		/* number of subband minima in Ep.r		*/
int     nki;            /* number of ki calculations                    */
int	nKz;		/* number of Kz values for lookup table		*/
int	ntheta;		/* number of points in theta integration 	*/
char	filename[9];	/* character string for output filename		*/
char	p;		/* particle					*/
bool	ff_flag;	/* form factor flag, output to file if true	*/
data11	*Gifsqr;	/* particular values of form factor squared	*/
data12	*wf;		/* start address of wave function structure	*/
FILE	*FACa;		/* pointer to absorption output file		*/
FILE	*FACe;		/* pointer to emission   output file		*/
FILE	*Frrp;		/* scattering rate required c-c rates		*/

/* default values */

A0=5.65*1e-10;			/* lattice constant for GaAs		*/
Da=7.0*e;			/* for GaAs				*/
Ephonon=2*1e-3*e;	    	/* bulk AC phonon energy, assume small	*/
ff_flag=false;			/* don't output formfactors	*/
m=0.067*me;			/* GaAs electron value		*/
p='e';				/* electron			*/
rho=5317.5;			/* density, default GaAs	*/
T=300;				/* temperature			*/
Vs=5117.0;			/* velocity of sound, GaAs	*/

/* default values for numerical calculations	*/

nKz=300;
nki=300;
ntheta=100;

/* calculate step lengths	*/

dKz=2/(A0*(float)nKz);		/* range of phonon integration as 2/A0	*/
dtheta=pi/((float)ntheta);	/* theta integration from 0 to pi	*/

while((argc>1)&&(argv[1][0]=='-'))
{
 switch(argv[1][1])
 {
  case 'A':
           A0=atof(argv[2])*1e-10;
           break;
  case 'a':
	   ff_flag=true;
	   argv--;
           argc++;
           break;
  case 'E':
	   Ephonon=atof(argv[2])*1e-3*e;
	   break;
  case 'm':
	   m=atof(argv[2])*me;
	   break;
  case 'p':
	   p=*argv[2];
	   switch(p)
	   {
	    case 'e': break;
	    case 'h': break;
	    case 'l': break;
	    default:  printf("Usage:  srcc [-p particle (\033[1me\033[0m, h, or l)]\n");
                      exit(0);
	   }
	   break;
  case 'T':
	   T=atof(argv[2]);
	   break;
  case 'v':
	   Vs=atof(argv[2]);
	   break;
  default :
	   printf("Usage:  sradp [-A lattice constant (\033[1m5.65\033[0mA)][-a generate form factors \033[1mfalse\033[0m]\n");
	   printf("              [-d density (\033[1m5317.5\033[0mkgm^-3)][-D deformation potential (\033[1m7.0\033[0meV)]\n");
	   printf("              [-E phonon energy (\033[1m2\033[0mmeV)]\n");
	   printf("              [-m mass (\033[1m0.067\033[0mm0)][-p particle (\033[1me\033[0m, h, or l)]\n");
	   printf("              [-T temperature (\033[1m300\033[0mK)][-v velocity of sound (\033[1m5117.0\033[0mm/s)]\n");
	   exit(0);

 }
 argv++;
 argv++;
 argc--;
 argc--;
}

/* calculate often used constants	*/

omega_0=Ephonon/hBar;		/* phonon angular frequency	*/
N0=1/(exp(Ephonon/(kB*T))-1);	/* Bose-Einstein factor	*/

Upsilon_a=gsl_pow_2(Da)*m*N0/(rho*Vs*4*gsl_pow_2(pi*hBar));

Upsilon_e=gsl_pow_2(Da)*m*(N0+1)/(rho*Vs*4*gsl_pow_2(pi*hBar));

E=read_E(p,&nE);	/* read in subband minima	*/

if((Frrp=fopen("rrp.r","r"))==0)
 {fprintf(stderr,"Error: Cannot open input file 'rrp.r'!\n");exit(0);}

while(fscanf(Frrp,"%i %i",&state[0],&state[1])!=EOF)
{
 wf=read_wf(&n,state,p);		/* reads potential file	*/

 delta_z=((wf+1)->a)-(wf->a);		/* Assumes regular grid	*/

 Gifsqr=ff_table(dKz,delta_z,wf,n,nKz);		/* generates formfactor table	*/

 if(ff_flag)ff_output(Gifsqr,nKz,state);	/* Outputs formfactors	*/

 /* Generate filename for particular mechanism and open file	*/

 sprintf(filename,"ACa%i%i.r",state[0],state[1]);	/* absorption	*/
 FACa=fopen(filename,"w");			
 sprintf(filename,"ACe%i%i.r",state[0],state[1]);	/* emission	*/
 FACe=fopen(filename,"w");			

 /* calculate Delta E, as a zero energy phonon is assumed, no need to 
    consider emission and absorption processes as in e-LO scattering	*/

 DeltaE=*(E+state[1]-1)-*(E+state[0]-1);	/* Ef-Ei	*/

 /* calculate maximum value of ki and hence ki step length */

 kimax=sqrt(2*m*(Vmax()-*(E+state[0]-1)))/hBar; /* gsl_pow_2(hBar*kimax)/2m=Vmax-Ei */
 dki=kimax/((float)nki);

 for(iki=0;iki<nki;iki++)       /* calculate e-AC rate for all ki	*/
 {
  ki=dki*(float)iki+dki/100;	/* second term avoids ki=0 pole	*/
  Wif=0;			/* Initialize for integration   */

  /* Integral around angle theta	*/

  for(itheta=0;itheta<ntheta;itheta++)
  {
   theta=dtheta*(float)itheta;

   /* Integral over phonon wavevector Kz	*/

   for(iKz=0;iKz<nKz;iKz++)
   {
    Kz=(Gifsqr+iKz)->a;

    arg=gsl_pow_2(ki*cos(theta))-2*m*DeltaE/gsl_pow_2(hBar);	/* sqrt argument */
    if(arg>=0)
    {
     alpha1=-ki*cos(theta)+sqrt(arg);
     alpha2=-ki*cos(theta)-sqrt(arg);

     /* alpha1 and alpha2 represent solutions for the in-plane polar
        coordinate Kxy of the carrier momentum---they must be positive, hence
        use Heaviside unit step function to ignore other contributions	*/

     Wif+=((Gifsqr+iKz)->b)*
	   (alpha1*Theta(alpha1)*gsl_hypot(alpha1, Kz)+
	    alpha2*Theta(alpha2)*gsl_hypot(alpha2, Kz))/
	   (alpha1-alpha2);
    }

   } /* end integral over Kz	*/
  } /* end integral over theta	*/

  Waif=2*Upsilon_a*Wif*dKz*dtheta;
  Weif=2*Upsilon_e*Wif*dKz*dtheta;	

  /* Now check for energy conservation!, would be faster with a nasty `if'
     statement just after the beginning of the ki loop!                 */

  Weif*=Theta(gsl_pow_2(hBar*ki)/(2*m)-DeltaE-Ephonon)*
        Theta(Vmax()-*(E+state[0]-1)+Ephonon-gsl_pow_2(hBar*ki)/(2*m));
  Waif*=Theta(gsl_pow_2(hBar*ki)/(2*m)-DeltaE+Ephonon)*
        Theta(Vmax()-*(E+state[0]-1)-Ephonon-gsl_pow_2(hBar*ki)/(2*m));

  /* output scattering rate versus carrier energy=subband minima+in-plane
     kinetic energy						*/

  fprintf(FACa,"%20.17le %20.17le\n",(*(E+state[0]-1)+gsl_pow_2(hBar*ki)/(2*m))/
                                    (1e-3*e),Waif);

  fprintf(FACe,"%20.17le %20.17le\n",(*(E+state[0]-1)+gsl_pow_2(hBar*ki)/(2*m))/
                                    (1e-3*e),Weif);

 }
 fclose(FACa);	/* close output file for this mechanism	*/
 fclose(FACe);	/* close output file for this mechanism	*/

 free(wf);	/* free memory of wavefunctions */

} /* end while over states */

free(E);

fclose(Frrp);

return EXIT_SUCCESS;
} /* end main */



data11
*ff_table(dKz,delta_z,wf,n,nKz)

/* This function outputs the formfactors into files	*/

double	dKz;
double	delta_z;
data12	*wf;
int	n;
int	nKz;

{
 double	Gsqr();

 double	Kz;		/* phonon wave vector				*/
 int	iKz;
 data11	*Gifsqr;	/* pointer to formfactor versus Kz data	*/

 Gifsqr=(data11 *)calloc(nKz,sizeof(data11));
  if (Gifsqr==0)  {
   fprintf(stderr,"Cannot allocate memory!\n");
   exit(0);
  }

 for(iKz=0;iKz<nKz;iKz++)
 {
  Kz=(float)iKz*dKz;
 
  (Gifsqr+iKz)->a=Kz;
  (Gifsqr+iKz)->b=Gsqr(Kz,delta_z,n,wf);


 } /* end Kz	*/

 return(Gifsqr);
}



double
Gsqr(Kz,delta_z,n,wf)

/* This function calculates the overlap integral squared between the two
   states	*/

double	Kz;
double	delta_z;
int	n;
data12	*wf;

{
 complex double	G;	/* integral over z and hence form factor	*/
 int		iz;	/* index over z		*/

 G=0;
 for(iz=0;iz<n;iz++)		/* Integral of i(=0) and f(=2) over z	*/
 {
  G += cos(Kz*(wf+iz)->a)*((wf+iz)->b[0])*((wf+iz)->b[1]) + I * sin(Kz*(wf+iz)->a)*((wf+iz)->b[0])*((wf+iz)->b[1]);
 }
 G*=delta_z;

 /* TODO: Replace with std::norm(G) in C++ */
return gsl_pow_2(cabs(G));	/* square-modulus of complex number	*/

}



void
ff_output(Gifsqr,nKz,state)

/* This function outputs the formfactors into files	*/

data11	*Gifsqr;
int	nKz;
int	state[];


{
 int	iKz;
 char	filename[9];	/* output filename				*/
 FILE	*FGif;		/* output file for form factors versus Kz	*/ 

 /* First generate filename and then open file	*/

 sprintf(filename,"G%i%i.r",state[0],state[1]);	
 if((FGif=fopen(filename,"w"))==0)
  {fprintf(stderr,"Error: Cannot open input file '%s'!\n",filename);exit(0);}

for(iKz=0;iKz<nKz;iKz++)
{
 fprintf(FGif,"%le %le\n",(Gifsqr+iKz)->a,(Gifsqr+iKz)->b);

} /* end Kz	*/

 fclose(FGif);
}



double
*read_E(p,nE)

/* This function reads the potential into memory and returns the start
   address of this block of memory and the number of lines	   */

char	p;
int	*nE;

{
 double	*E;
 int	i=0;		/* index over the energies			*/
 char	filename[9];	/* filename string				*/
 FILE 	*FE;		/* file pointer to energy data 			*/

 sprintf(filename,"E%c.r",p);
 if((FE=fopen(filename,"r"))==0)
 {
   fprintf(stderr,"Error: Cannot open input file '%s'!\n",filename);
   exit(0);
 }

 *nE=0;
 while(fscanf(FE,"%*i %*le")!=EOF)
  (*nE)++;
 rewind(FE);

 E=(double *)calloc(*nE,sizeof(double));
 if (E==0)  {
  fprintf(stderr,"Cannot allocate memory!\n");
  exit(0);
 }

 while(fscanf(FE,"%*i %le",E+i)!=EOF)
 {
  *(E+i)*=1e-3*e;		/*convert meV->J		*/
  i++;
 }

 fclose(FE);

 return(E);

}



data12
*read_wf(n,state,p)

/* This function reads the potential into memory and returns the start
   address of this block of memory and the number of lines	   */

int	*n;
int	state[];
char	p;

{
 char	filename[9];	/* input filename			*/
 int	i;
 int	ii;
 FILE 	*Fwf;		/* file pointer to wavefunction file	*/
 FILE	*Fv;		/* file pointer to potential file	*/
 data12	*wf;		/* pointer to wave function structure	*/

 /* Use potential file for number of lines	*/
 
 if((Fv=fopen("v.r","r"))==0)
 {fprintf(stderr,"Error: Cannot open input file 'v.r'!\n");exit(0);}

 *n=0;	
 while(fscanf(Fv,"%*le %*le")!=EOF)
  (*n)++;
 fclose(Fv);

 /* Allocate memory for all four wavefunctions */

 wf=(data12 *)calloc(*n,sizeof(data12));
 if(wf==0){fprintf(stderr,"Cannot allocate memory!\n");exit(0);}

 /* Read in both wave functions	*/

 for(i=0;i<2;i++)
 {
  sprintf(filename,"wf_%c%i.r",p,state[i]);	/* Open each file	*/
  if((Fwf=fopen(filename,"r"))==0)
   {fprintf(stderr,"Error: Cannot open input file '%s'!\n",filename);exit(0);}

  ii=0;						/* Read in each file	*/
  while(fscanf(Fwf,"%le %le",&((wf+ii)->a),&((wf+ii)->b[i]))!=EOF)
   ii++;

 fclose(Fwf);					/* Close each file	*/
 }

 return(wf);

}



double
Vmax()

/* This function scans the file v.r and returns the maximum value of the
   potential.
                                                                        */
{
 double max;                    /* maximum value of potential energy    */
 double v;                      /* potential                            */
 FILE   *Fv;                    /* file pointer to v.r                  */

max=0;

if((Fv=fopen("v.r","r"))==0)
 {fprintf(stderr,"Error: Cannot open input file 'v.r'!\n");exit(0);}

while(fscanf(Fv,"%*le %le",&v)!=EOF)
 if(v>max) max=v;

fclose(Fv);

return(max);

}

