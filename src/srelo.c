/*==================================================================
              srelo Scattering Rate Electron-LO phonon
  ==================================================================*/

/* This program calculates the electron-LO phonon scattering rate for
   both intra- and intersubband events.  The required rates are provided
   by the user in the file `rrp.r'.  The other necessary inputs are listed
   below.


	Input files:		rrp.r	contains required rates
				wf_xy.r	x=particle y=state
				Ex.r	x=particle, energies

	Output files:		LO[a,e]if.r absorption and emission rates


    Paul Harrison, January 1998 
 
    Improvements January 1999						*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <signal.h>
#include <malloc.h>
#include "struct.h"
#include "const.h"
#include "maths.h"
#include "bools.h"
#include "io.h"


main(int argc,char *argv[])
{
double	Gsqr();		/* form factor calculation			*/
double	*read_E();	/* reads subband minima energies from file	*/
double	Vmax();		/* maximum value of the potential		*/
void    ff_output();    /* output formfactors                           */
data11	*ff_table();
data12	*read_wf();	/* reads wavefunctions into memory		*/

double	A0;		/* lattice constant				*/
double	Delta_a;	/* Ef-Ei-Ephonon, see notes			*/
double	Delta_e;	/* Ef-Ei+Ephonon, see notes			*/
double	delta_z;	/* step length along growth (z-) axis		*/
double  dki;            /* step length for loop over ki                 */
double	dKz;		/* step length for integration over Kz	*/
double	*E;		/* pointer to subband minima			*/
double	Ephonon;	/* phonon energy				*/
double	epsilon_s;	/* low frequency dielectric constant		*/
double	epsilon_inf;	/* high frequency dielectric constant		*/
double  kimax;          /* maximum value of ki                          */
double  ki;             /* carrier momentum (wave vector actually)	*/
double	Kz;		/* phonon wavevector				*/
double	m;		/* carrier effective mass			*/
double	N0;		/* number of phonons at LO phonon energy	*/
double	omega_0;	/* angular frequency of LO phonon		*/
double	T;		/* temperature					*/
double	Upsilon_star_a;	/* scattering rate prefactor			*/
double	Upsilon_star_e;	/* scattering rate prefactor			*/
double	Waif;		/* absorption scattering rate			*/
double	Weif;		/* emission scattering rate			*/
int	iKz;		/* index over Kz				*/
int     iki;            /* index over ki                                */
int	state[2];	/* electron state index				*/
int	n;		/* length of wavefunctions file			*/
int	nE;		/* number of subband minima in Ep.r		*/
int     nki;            /* number of ki calculations                    */
int	nKz;		/* number of Kz values for lookup table	*/
char	filename[9];	/* character string for output filename		*/
char	p;		/* particle					*/
boolean	ff_flag;	/* form factor flag, output to file if true	*/
data11	*Gifsqr;	/* particular values of form factor squared	*/
data12	*wf;		/* start address of wave function structure	*/
FILE	*FLOa;		/* pointer to absorption output file		*/
FILE	*FLOe;		/* pointer to emission   output file		*/
FILE	*Frrp;		/* scattering rate required c-c rates		*/

/* default values */

A0=5.65*1e-10;			/* lattice constant for GaAs		*/
Ephonon=36*1e-3*e_0;    	/* bulk LO phonon energy, 36 meV in GaAs*/
epsilon_s=13.18*epsilon_0;	/* low frequency dielectric constant for GaAs*/
epsilon_inf=10.89*epsilon_0;	/* high frequency dielectric constant for GaAs*/
ff_flag=false;			/* don't output formfactors	*/
m=0.067*m0;			/* GaAs electron value		*/
p='e';				/* electron			*/
T=300;				/* temperature			*/

/* default values for numerical calculations	*/

nKz=1000;
nki=1000;

/* calculate step lengths	*/

dKz=2/(A0*(float)nKz);	/* Taken range of phonon integration as 2/A0 */

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
	   Ephonon=atof(argv[2])*1e-3*e_0;
	   break;
  case 'e':
	   epsilon_s=atof(argv[2])*epsilon_0;
	   break;
  case 'f':
	   epsilon_inf=atof(argv[2])*epsilon_0;
	   break;
  case 'm':
	   m=atof(argv[2])*m0;
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
  default :
	   printf("Usage:  srelo [-A lattice constant (\033[1m5.65\033[0mA)][-a generate form factors \033[1mfalse\033[0m]\n");
	   printf("              [-e low frequency dielectric constant (\033[1m13.18\033[0mepsilon_0]\n");
	   printf("              [-f high frequency dielectric constant (\033[1m10.89\033[0mepsilon_0]\n");
	   printf("              [-m mass (\033[1m0.067\033[0mm0)][-p particle (\033[1me\033[0m, h, or l)]\n");
	   printf("              [-T temperature (\033[1m300\033[0mK)]\n");
	   exit(0);

 }
 argv++;
 argv++;
 argc--;
 argc--;
}

/* calculate often used constants	*/

omega_0=Ephonon/hbar;		/* phonon angular frequency	*/
N0=1/(exp(Ephonon/(kb*T))-1);	/* Bose-Einstein factor	*/

Upsilon_star_a=pi*e_0*e_0*omega_0/epsilon_s*(epsilon_s/epsilon_inf-1)*(N0)
               *2*m/sqr(hbar)*2/(8*pi*pi*pi);

Upsilon_star_e=pi*e_0*e_0*omega_0/epsilon_s*(epsilon_s/epsilon_inf-1)*(N0+1)
               *2*m/sqr(hbar)*2/(8*pi*pi*pi);

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

 sprintf(filename,"LOa%i%i.r",state[0],state[1]);	/* absorption	*/
 FLOa=fopen(filename,"w");			
 sprintf(filename,"LOe%i%i.r",state[0],state[1]);	/* emission	*/
 FLOe=fopen(filename,"w");			

 /* calculate Delta variables, constant for each mechanism	*/

 Delta_a=*(E+state[1]-1)-*(E+state[0]-1)-Ephonon;	/* Ef-Ei-Ephonon */

 Delta_e=*(E+state[1]-1)-*(E+state[0]-1)+Ephonon;	/* Ef-Ei+Ephonon */

 /* calculate maximum value of ki and hence ki step length */

 kimax=sqrt(2*m*(Vmax()-*(E+state[0]-1)))/hbar; /* sqr(hbar*kimax)/2m=Vmax-Ei */
 dki=kimax/((float)nki);

 for(iki=0;iki<nki;iki++)       /* calculate e-LO rate for all ki	*/
 {
  ki=dki*(float)iki;
  Waif=0;                       /* Initialize for integration   */
  Weif=0;                       /* Initialize for integration   */

  /* Integral over phonon wavevector Kz	*/

  for(iKz=0;iKz<nKz;iKz++)
  {
   Kz=(Gifsqr+iKz)->a;

   Waif+=((Gifsqr+iKz)->b)/
         sqrt(sqr(sqr(Kz))+
              2*sqr(Kz)*(2*sqr(ki)-2*m*Delta_a/sqr(hbar))+
              sqr(2*m*Delta_a/sqr(hbar))
             );

   Weif+=((Gifsqr+iKz)->b)/
         sqrt(sqr(sqr(Kz))+
              2*sqr(Kz)*(2*sqr(ki)-2*m*Delta_e/sqr(hbar))+
              sqr(2*m*Delta_e/sqr(hbar))
             );



  } /* end integral over Kz	*/

  Waif*=Upsilon_star_a*pi*dKz;	/* Note integral from 0->inf, hence *2	*/
  Weif*=Upsilon_star_e*pi*dKz;	/* Note integral from 0->inf, hence *2	*/

  /* Now check for energy conservation!, would be faster with a nasty `if'
     statement just after the beginning of the ki loop!			*/

  Weif*=Theta(sqr(hbar*ki)/(2*m)-Delta_e)*
        Theta(Vmax()-*(E+state[0]-1)+Ephonon-sqr(hbar*ki)/(2*m));

  Waif*=Theta(sqr(hbar*ki)/(2*m)-Delta_a)*
        Theta(Vmax()-*(E+state[0]-1)-Ephonon-sqr(hbar*ki)/(2*m));

  /* output scattering rate versus carrier energy=subband minima+in-plane
     kinetic energy						*/

  fprintf(FLOa,"%20.17le %20.17le\n",(*(E+state[0]-1)+sqr(hbar*ki)/(2*m))/
                                    (1e-3*e_0),Waif);

  fprintf(FLOe,"%20.17le %20.17le\n",(*(E+state[0]-1)+sqr(hbar*ki)/(2*m))/
                                    (1e-3*e_0),Weif);

 }
 fclose(FLOa);	/* close output file for this mechanism	*/
 fclose(FLOe);	/* close output file for this mechanism	*/

 free(wf);	/* free memory of wavefunctions */

} /* end while over states */

free(E);

fclose(Frrp);

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
 complex	G;	/* integral over z and hence form factor	*/
 int		iz;	/* index over z		*/

 G.re=0;G.im=0;
 for(iz=0;iz<n;iz++)		/* Integral of i(=0) and f(=2) over z	*/
 {
  G.re+=cos(Kz*(wf+iz)->a)*((wf+iz)->b[0])*((wf+iz)->b[1]);
  G.im+=sin(Kz*(wf+iz)->a)*((wf+iz)->b[0])*((wf+iz)->b[1]);
 }
 G.re*=delta_z;
 G.im*=delta_z;

return(cmod(G)*cmod(G));	/* cmod---modulus of complex number	*/

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
  *(E+i)*=1e-3*e_0;		/*convert meV->J		*/
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

