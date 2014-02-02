/*==================================================================
              srcc  Scattering Rate Carrier-Carrier
  ==================================================================*/

/* This program calculates the carrier-carrier scattering rate for
   both intra- and intersubband events.  The required rates are provided
   by the user in the file `rr.r'.  The other necessary inputs are listed
   below.

	Input files:	rr.r	contains required rates
			wf_xy.r	x=particle y=state
			N.r	subband populations
			Ex.r	x=particle, energies
			Ef.r	subband Fermi energies

	Output files:	ccABCD.r	cc rate versus Ei for each mechanism


    Paul Harrison, March 1997
    Modifications, February 1999					*/

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


typedef
struct	{
 double	a;		/* z value of files			*/
 double	b[4];		/* wave function			*/
} data14;

main(int argc,char *argv[])
{
double	A();		/* form factor calculation			*/
double	lookup_ff();	/* look up form factor in table			*/
double	lookup_PI();	/* look up screening function in table		*/
double	PI();		/* screening term				*/
double	*read_E();	/* reads subband minima energies from file	*/
double	*read_Ef();	/* reads Fermi energies from file		*/
double	*read_N();	/* reads subband populations from file		*/
double	Vmax();		/* maximum value of the potential		*/
void	output_ff();	/* output formfactors				*/
data11	*ff_table();	/* generates form factor v q_perp lookup table	*/
data11	*PI_table();	/* generates screening function lookup table	*/
data14	*read_wf();	/* reads wavefunctions into memory		*/

double	alpha;		/* angle between ki and kj			*/
double	Deltak0sqr;	/* twice the change in KE, see Smet (55)	*/
double	dalpha;		/* step length for alpha integration		*/
double	delta_z;	/* mesh length along growth (z-) axis		*/
double	dtheta;		/* step length for theta integration		*/
double	dki;		/* step length for loop over ki 		*/
double	dkj;		/* step length for kj integration		*/
double	*E;		/* pointer to subband minima			*/
double	*Ef;		/* pointer to Fermi energies			*/
double	epsilon;	/* low frequency dielectric constant		*/
double	ki;		/* carrier momentum				*/
double	kifermi;	/* Fermi wave vector of state i			*/
double	kij;		/* (vector)kj-(vector)(ki)			*/
double	kimax;		/* maximum value of ki				*/
double	kj;		/* carrier momentum				*/
double	kjmax;		/* maximum value of kj				*/
double	m;		/* carrier effective mass			*/
double	*N;		/* pointer to subband populations		*/
double	P;		/* probability factor, see Smet			*/
double	q_perp;		/* in-plane momentum, |ki-kf|			*/
double	q_perpsqr4;	/* 4*q_perp*q_perp, checks vailidity of q_perp	*/
double	T;		/* temperature					*/
double	theta;		/* angle between kij and kfg			*/
double	W;		/* arbitrary well width, soley for output	*/
double	Wbar;		/* FD weighted mean of Wijfg			*/
double	Wijfg;		/* the carrier-carrier scattering rate		*/
int	i;		/* general index				*/
int	ialpha;		/* index over alpha				*/
int	iki;		/* index over ki				*/
int	ikj;		/* index over kj				*/
int	itheta;		/* index over theta				*/
int	state[4];	/* electron state index				*/
int	n;		/* length of wavefunctions file			*/
int	nalpha;		/* number of strips in alpha integration	*/
int	nE;		/* number of subband minima in Ep.r		*/
int	nki;		/* number of ki calculations			*/
int	nkj;		/* number of strips in |kj| integration		*/
int	nq;		/* number of q_perp values for lookup table	*/
int	ntheta;		/* number of strips in theta integration	*/
char	filename[9];	/* character string for output filename		*/
char	p;		/* particle					*/
boolean	ff_flag;	/* form factor flag, output to file if true	*/
boolean	S_flag;		/* screening flag, include screening if true	*/
data11	*Aijfg;		/* form factor over all 4 wave functions	*/
data11	*PIii;		/* pointer to screening function versus q_perp	*/
data14	*wf;		/* start address of wave function structure	*/
FILE	*Fcc;		/* pointer to output file			*/
FILE	*FccABCD;	/* pointer to weighted mean output file		*/
FILE	*Frr;		/* scattering rate required c-c rates		*/

/* default values */

epsilon=13.18*epsilon_0;/* low frequency dielectric constant for GaAs	*/
ff_flag=false;		/* don't output formfactors	*/
m=0.067*m0;		/* GaAs electron value		*/
p='e';			/* electron			*/
T=300;			/* temperature			*/
W=250e-10;		/* a well width, same as Smet	*/
S_flag=true;		/* include screening by default	*/

/* default values for numerical calculations	*/

nalpha=100;
ntheta=100;
nki=100;
nkj=100;
nq=100;

/* calculate step lengths	*/

dalpha=2*pi/((float)nalpha);
dtheta=2*pi/((float)ntheta);

while((argc>1)&&(argv[1][0]=='-'))
{
 switch(argv[1][1])
 {
  case 'a':
	   ff_flag=true;
	   argv--;
           argc++;
           break;
  case 'e':
	   epsilon=atof(argv[2])*epsilon_0;
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
  case 'S':
	   S_flag=false;
	   argv--;
           argc++;
           break;
  case 'T':
	   T=atof(argv[2]);
	   break;
  case 'w':
	   W=atof(argv[2])*1e-10;
	   break;
  default :
	   printf("Usage:  srcc [-a generate form factors \033[1mfalse\033[0m][-e permittivity (\033[1m13.18\033[0mepsilon_0)]\n");
	   printf("             [-m mass (\033[1m0.067m0\033[0m)][-p particle (\033[1me\033[0m, h, or l)]\n");
	   printf("             [-S screening \033[1mtrue\033[0m]\n");
	   printf("             [-T temperature (\033[1m300\033[0mK)][-w a well width (\033[1m250\033[0mA)]\n");
	   exit(0);

 }
 argv++;
 argv++;
 argc--;
 argc--;
}

E=read_E(p,&nE);	/* read in subband minima	*/

Ef=read_Ef(p,nE);	/* read in Fermi energies	*/

N=read_N(nE);		/* read in subband populations	*/

if((Frr=fopen("rr.r","r"))==0)
 {fprintf(stderr,"Error: Cannot open input file 'rr.r'!\n");exit(0);}

FccABCD=fopen("ccABCD.r","w");	/* open file for output of weighted means */

while(fscanf(Frr,"%i %i %i %i",&state[0],&state[1],&state[2],&state[3])!=EOF)
{
 wf=read_wf(&n,state,p);		/* reads potential file	*/

 delta_z=((wf+1)->a)-(wf->a);		/* Assumes regular grid	*/

 if(ff_flag)output_ff(delta_z,W,wf,n,state);	/* Outputs formfactors	*/

 /* Generate filename for particular mechanism and open file	*/

 sprintf(filename,"cc%i%i%i%i.r",state[0],state[1],state[2],state[3]);
 Fcc=fopen(filename,"w");			

 /* Calculate Delta k0^2	*/

 if((state[0]+state[1])==(state[2]+state[3])) Deltak0sqr=0;
 else
  Deltak0sqr=4*m*(*(E+state[0]-1)+*(E+state[1]-1)
                 -*(E+state[2]-1)-*(E+state[3]-1))/(hbar*hbar);	

 Aijfg=ff_table(Deltak0sqr,delta_z,m,E,wf,n,nq,state);	/* generate Aijfg table	*/

 /* Calculate Fermi wave vector of state i, kifermi, note divisor `1'
    is a degeneracy factor representing number of valleys, GaAs=1	*/

 kifermi=sqrt(2*pi*(*(N+state[0]-1))/1);

 PIii=PI_table(Aijfg,m,E,Ef,kifermi,T,n,nq,state,S_flag);/* generate PIii table	*/

 /* calculate maximum value of ki & kj and hence kj step length	*/

 kimax=sqrt(2*m*(Vmax()-*(E+state[0]-1)))/hbar;	/* sqr(hbar*kimax)/2m=Vmax-Ei	*/
 dki=kimax/((float)nki);

 kjmax=sqrt(2*m*(Vmax()-*(E+state[1]-1)))/hbar;	/* sqr(hbar*kjmax)/2m=Vmax-Ej	*/
 dkj=kjmax/((float)nkj);

 Wbar=0;			/* initialise integral sum */

 for(iki=0;iki<nki;iki++)	/* calculate c-c rate for all ki	*/
 {
 ki=dki*(float)iki;
 Wijfg=0;			/* Initialize for integration	*/
 for(ikj=0;ikj<nkj;ikj++)	/* integrate over |kj|		*/
 {
  kj=dkj*(float)ikj;

  P=1/(exp((*(E+state[1]-1)+hbar*kj*hbar*kj/(2*m)-*(Ef+state[1]-1))/(kb*T))+1);

  for(ialpha=0;ialpha<nalpha;ialpha++)	/* Integral over alpha	*/
  {
   alpha=dalpha*(float)ialpha;
 
   kij=sqrt(ki*ki+kj*kj-2*ki*kj*cos(alpha));	/* calculate kij	*/

   for(itheta=0;itheta<ntheta;itheta++)		/* Integral over theta	*/
   {
    theta=dtheta*(float)itheta;	/* move origin to avoid 0*/
  
    /* calculate argument of sqrt function=4*q_perp*q_perp, see (8.179),
       to check for imaginary q_perp, if argument is positive, q_perp is
       real and hence calculate scattering rate, otherwise ignore and move
       onto next q_perp							*/

    q_perpsqr4=2*kij*kij+Deltak0sqr-2*kij*sqrt(kij*kij+Deltak0sqr)*cos(theta);

    /* recall areal element in plane polars is kj.dkj.dalpha */

    if(q_perpsqr4>=0) 
    {
     q_perp=sqrt(q_perpsqr4)/2;
     Wijfg+=sqr(lookup_ff(Aijfg,q_perp,nq)/
                (q_perp+2*pi*e_0*e_0*lookup_PI(PIii,q_perp,nq)
                 *lookup_ff(Aijfg,q_perp,nq)/(4*pi*epsilon)
                )
               )*P*kj;
    }

    /* Note screening term has been absorbed into the above equation
       to avoid pole at q_perp=0				*/

   } /* end theta */

  } /* end alpha */

 } /* end kj   */
 
 Wijfg*=dtheta*dalpha*dkj;	/* multiply by all step lengths	*/

 Wijfg*=sqr(e_0*e_0/(hbar*4*pi*epsilon))*m/(pi*hbar);

 /* output scattering rate versus carrier energy=subband minima+in-plane
    kinetic energy						*/

 fprintf(Fcc,"%20.17le %20.17le\n",(*(E+state[0]-1)+sqr(hbar*ki)/(2*m))/
                                   (1e-3*e_0),Wijfg);

 /* calculate Fermi-Dirac weighted mean of scattering rates over the 
    initial carrier states, note that the integral step length 
    dE=2*sqr(hbar)*ki*dki/(2m)					*/

 Wbar+=Wijfg*ki/(exp((*(E+state[0]-1)+sqr(hbar*ki)/(2*m)-*(Ef+state[0]-1))/(kb*T))+1);

 } /* end ki	*/

 Wbar*=dki/(pi*(*(N+state[0]-1)));

 fprintf(FccABCD,"%i %i %i %i %20.17le\n",state[0],state[1],state[2],state[3],
                                          Wbar);

 fclose(Fcc);	/* close output file for this mechanism	*/

 free(wf);	/* free memory of wavefunctions, formfactors and screening */
 free(Aijfg);
 free(PIii);

} /* end while over states */

free(E);
free(Ef);

fclose(FccABCD);	/* close weighted mean output file	*/
fclose(Frr);

} /* end main */



double
A(delta_z,q_perp,n,wf)

/* This function calculates the overlap integral over all four carrier
   states		*/

double	delta_z;
double	q_perp;
int	n;
data14	*wf;

{
 double	A;	/* integral over z and hence form factor	*/
 double	B;	/* integral over z'	*/
 int	iz;	/* index over z		*/
 int	izd;	/* index over z'	*/

 A=0;
 for(iz=0;iz<n;iz++)		/* Integral of i(=0) and f(=2) over z	*/
 {
  B=0;
  for(izd=0;izd<n;izd++)	/* Integral of j(=1) and g(=3) over z'	*/
  {
   B+=((wf+izd)->b[1])*((wf+izd)->b[3])
      *exp(-q_perp*fabs(((wf+iz)->a)-((wf+izd)->a)));
  }
  B*=delta_z;
  A+=(((wf+iz)->b[0])*((wf+iz)->b[2]))*B;
 }
 A*=delta_z;

return(A);

}



double
PI(E,Ef,kifermi,m,q_perp,T,state)

/* This function returns the screening factor, referred to by Smet as 
   e_sc
									*/

double	*E;		/* pointer to subband minima		*/
double	*Ef;		/* pointer to Fermi energies		*/
double	kifermi;	/* Fermi wave vector of state i		*/
double	m;		/* effective mass			*/
double	q_perp;		/* In-plane wave vector			*/
double	T;		/* temperature				*/
int	state[];	/* electron state index			*/
{
 double	dmu;		/* mu interval				*/
 double	dI;		/* intervals in I			*/
 double	mu;		/* integration variable			*/
 double	I;		/* temperature dependence integral	*/
 double	P;		/* polarizability			*/

 /* Equation 43 of Smet	*/

 if(q_perp<=2*kifermi)
  P=m/(pi*hbar*hbar);
 else
  P=m/(pi*hbar*hbar)*(1-sqrt(1-sqr(2*kifermi/q_perp)));

 /* Now perform the integration, equation 44 of Smet	*/

 mu=*(E+state[0]-1);dmu=1e-3*e_0;I=0;
 do
 {
  dI=1/(4*kb*T*sqr(cosh((*(Ef+state[0]-1)-mu)/(2*kb*T))));
  I+=dI*dmu;
  mu+=dmu;
 }while(dI>I/100);	/* continue until integral converges	*/

 P*=I;

return(P);

}



data11
*PI_table(Aijfg,m,E,Ef,kifermi,T,n,nq,state,S_flag)

/* This function creates the polarizability table */

data11	*Aijfg;
double	m;
double	*E;
double	*Ef;
double	kifermi;
double	T;
int	n;
int	nq;		/* Number of q_perp values for table	*/
int	state[];
boolean	S_flag;


{
 double	PI();

 double	q_perp;	/* In-plane wave vector			*/
 int	iq;
 data11	*PIii;	/* pointer to screening function data	*/

 PIii=(data11 *)calloc(nq,sizeof(data11));
 if (PIii==0)  {
  fprintf(stderr,"Cannot allocate memory!\n");
  exit(0);
 }

 for(iq=0;iq<nq;iq++)
 {
  q_perp=(Aijfg+iq)->a;
  (PIii+iq)->a=q_perp;

  /* Allow screening to be turned off	*/

  if(S_flag)(PIii+iq)->b=PI(E,Ef,kifermi,m,q_perp,T,state);
  else (PIii+iq)->b=0;
 }

 return(PIii);
}




data11
*ff_table(Deltak0sqr,delta_z,m,E,wf,n,nq,state)

/* This function creates a form factor look up table	*/

double	Deltak0sqr;
double	delta_z;
double	m;
double	*E;
data14	*wf;
int	n;
int	nq;		/* Number of q_perp values for table	*/
int	state[];


{
 double	A();
 double	Vmax();

 double	dq;	/* interval in q_perp			*/
 double	kimax;	/* maximum value of ki			*/
 double	kjmax;	/* maximum value of ki			*/
 double	q_perp;	/* In-plane wave vector			*/
 double	q_perp_max;	/* maximum in-plane wave vector	*/
 double	vmax;	/* maximum potential			*/
 int	iq;
 data11	*Aijfg;

 vmax=Vmax();
 kimax=sqrt(2*m*(vmax-*(E+state[0]-1)))/hbar;	/* sqr(hbar*kimax)/2m=Vmax-Ei	*/
 kjmax=sqrt(2*m*(vmax-*(E+state[1]-1)))/hbar;	/* sqr(hbar*kfmax)/2m=Vmax-Ef	*/

 Aijfg=(data11 *)calloc(nq,sizeof(data11));
  if (Aijfg==0)  {
   fprintf(stderr,"Cannot allocate memory!\n");
   exit(0);
  }

 q_perp_max=sqrt(2*sqr(kimax+kjmax)+Deltak0sqr+2*(kimax+kjmax)*
                 sqrt(sqr(kimax+kjmax)+Deltak0sqr))/2;

 dq=q_perp_max/((float)(nq-1));
 for(iq=0;iq<nq;iq++)
 {
  q_perp=dq*(float)iq;
  (Aijfg+iq)->a=q_perp;
  (Aijfg+iq)->b=A(delta_z,q_perp,n,wf);
 }

 return(Aijfg);
}



double
lookup_ff(Aijfg,q_perp,nq)

/* This function looks up the form factor Aijfg in the table generated
   by ff_table								*/

data11	*Aijfg;
double	q_perp;
int	nq;

{
 double	A;	/* The looked up form factor	*/
 int	iq=0;	/* index over q_perp in table	*/

 if(q_perp>((Aijfg+nq-1)->a))
 {printf("q_perp>maximum allowed q!\n");exit(0);}

 while(q_perp>=((Aijfg+iq)->a))
 {
  iq++;		/* retreive the ff above q_perp	*/
 }

 /* Linearly interpolate the form factor from between the values
    directly above and below q_perp				*/

 A=((Aijfg+iq-1)->b)+(((Aijfg+iq)->b)-((Aijfg+iq-1)->b))*
                   (q_perp-((Aijfg+iq-1)->a))/
                   (((Aijfg+iq)->a)-((Aijfg+iq-1)->a));

 return(A);

}



double
lookup_PI(PIii,q_perp,nq)

/* This function looks up the screening function PIii in the table 
   generated by PI_table					*/

data11	*PIii;
double	q_perp;
int	nq;

{
 double	P;	/* The looked up screening function	*/
 int	iq=0;	/* index over q_perp in table		*/

 if(q_perp>((PIii+nq-1)->a))
 {printf("q_perp>maximum allowed q!\n");exit(0);}

 while(q_perp>=((PIii+iq)->a))
 {
  iq++;		/* retreive the ff above q_perp	*/
 }

 /* Linearly interpolate the form factor from between the values
    directly above and below q_perp				*/

 P=((PIii+iq-1)->b)+(((PIii+iq)->b)-((PIii+iq-1)->b))*
                   (q_perp-((PIii+iq-1)->a))/
                   (((PIii+iq)->a)-((PIii+iq-1)->a));

 return(P);

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



double
*read_Ef(p,nE)

/* This function reads the potential into memory and returns the start
   address of this block of memory and the number of lines	   */

char	p;
int	nE;

{
 double	*Ef;
 int	i=0;		/* index over the energies			*/
 char	filename[9];	/* filename string				*/
 FILE 	*FEf;		/* file pointer to energy data 			*/

 sprintf(filename,"Ef.r",p);	/* Retain general structure for another day	*/
 if((FEf=fopen(filename,"r"))==0)
 {
   fprintf(stderr,"Error: Cannot open input file '%s'!\n",filename);
   exit(0);
 }

 Ef=(double *)calloc(nE,sizeof(double));
 if (Ef==0)  {
  fprintf(stderr,"Cannot allocate memory!\n");
  exit(0);
 }

 while(fscanf(FEf,"%*i %le",Ef+i)!=EOF)
 {
  *(Ef+i)*=1e-3*e_0;		/*convert meV->J		*/
  i++;
 }

 fclose(FEf);

 return(Ef);
}



double
*read_N(nE)

/* This function reads the subband populations into memory and returns 
   the start address of this block of memory	   */

int	nE;

{
 double	*N;
 int	i=0;		/* index over the energies			*/
 char	filename[9];	/* filename string				*/
 FILE 	*FN;		/* file pointer to subband populations data	*/

 if((FN=fopen("N.r","r"))==0)
 {
   fprintf(stderr,"Error: Cannot open input file 'N.r'!\n");
   exit(0);
 }

 N=(double *)calloc(nE,sizeof(double));
 if (N==0)  {
  fprintf(stderr,"Cannot allocate memory!\n");
  exit(0);
 }

 while(fscanf(FN,"%*i %le",N+i)!=EOF)
 {
  *(N+i)*=1e+10*1e+4;	/*convert from units of 10^10cm^-2->m^-2 */
  i++;
 }

 fclose(FN);

 return(N);
}



data14
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
 data14	*wf;		/* pointer to wave function structure	*/

 /* Use potential file for number of lines	*/
 
 if((Fv=fopen("v.r","r"))==0)
 {fprintf(stderr,"Error: Cannot open input file 'v.r'!\n");exit(0);}

 *n=0;	
 while(fscanf(Fv,"%*le %*le")!=EOF)
  (*n)++;
 fclose(Fv);

 /* Allocate memory for all four wavefunctions */

 wf=(data14 *)calloc(*n,sizeof(data14));
 if(wf==0){fprintf(stderr,"Cannot allocate memory!\n");exit(0);}

 /* Read in all four wave functions	*/

 for(i=0;i<4;i++)
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



void
output_ff(delta_z,W,wf,n,state)

/* This function outputs the formfactors into files	*/

double	delta_z;
double	W;		/* Arbitrary well width to generate q	*/
data14	*wf;
int	n;
int	state[];


{
 double	A();

 double	Aijfg;
 double	q_perp;		/* In-plane wave vector				*/
 int	iq;
 char	filename[9];	/* output filename				*/
 FILE	*FA;		/* output file for form factors versus q_perp	*/

 /* First generate filename and then open file	*/

 sprintf(filename,"A%i%i%i%i.r",state[0],state[1],state[2],state[3]);	
 if((FA=fopen(filename,"w"))==0)
  {fprintf(stderr,"Error: Cannot open input file '%s'!\n",filename);exit(0);}

 for(iq=0;iq<100;iq++)
 {
  q_perp=6*iq/(100*W);
  Aijfg=A(delta_z,q_perp,n,wf);
  fprintf(FA,"%le %le\n",q_perp*W,sqr(Aijfg));
 }

 fclose(FA);
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

