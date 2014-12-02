/*==================================================================
                                i0
  ==================================================================*/

/* This program implements a variational technique to calculate the
   uncorrelated one particle energies of an electron or hole attatched
   to a single donor or acceptor at any position, in any user supplied 
   potential.  The potential is read from the file v.r

   This version is a single variational parameter calculation---the 
   three-dimensional (3D) approximation trial wavefunction:

		Psi=psi(z) phi(r)

   where psi(z) is the one-particle wave function without the impurity
   and phi(r) is the hydorgenic-like term, which can take any form, i.e.
   1s, 2s, 2pz, 2px.

		Input files:
		r_i.r		donor (or acceptor positions)
		v.r		one-dimensional potential
		wf_pn.r		wave functions

		Output files:
		e.r		total energies for each r_i
		l.r		Bohr radii (lambda) for each r_i

   Paul Harrison, June 2001
   								*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <signal.h>
#include <malloc.h>
#include "const.h"
#include "struct.h"
#include "maths.h"
#include "bools.h"
#include "io.h"

main(int argc,char *argv[])
{
double	read_delta_z();
double	Energy();	/* expectation value of Hamiltonian	*/
data11 *read_v();	/* reads potential file into memory	*/
data11 *read_wf();	/* reads eigenvector into memory	*/

double	dz;		/* z separation of input potentials	*/
double	E;		/* electron (or hole) energies		*/
double	epsilon;	/* permitivity of material		*/
double	f,fdash;	/* function (dE/dlambda) and derivative 
			   to be solved				*/
double	lambda;		/* Bohr radius (variational)		*/
double	lambda_step;	/* Bohr radius increment		*/
double	lambda_0;	/* Bohr radius of bulk impurity		*/
double	m;		/* mass of particle			*/
double	r_i;		/* position of impurity			*/
double	y1,y2,y3;	/* E(lambda-), E(lambda), E(lambda+)	*/
int	i_i;		/* index of impurity			*/
int	n;		/* number of lines in potential file	*/
int	state;		/* principal quantum number		*/
int	S;		/* impurity level, i.e. `1s', `2px' etc	*/
char	p;		/* particle				*/
char	State[9];	/* string containing impurity level	*/
data11  *V;		/* start address of potential		*/
data11  *wf;		/* start address of wave function	*/
FILE   *fe;		/* file pointer for energies		*/
FILE   *fl;		/* file pointer for lambda_0		*/
FILE   *fr_i;		/* file pointer to donor positions	*/
 

/* default values */

epsilon=13.18*epsilon_0;
m=0.067*m0;
p='e';
state=1;
S=1;

/* computational default values */

i_i=0;          

while((argc>1)&&(argv[1][0]=='-'))
{
 switch(argv[1][1])
 {
  case 'e':
	   epsilon=atof(argv[2])*epsilon_0;
	   break;
  case 'm':
	   m=atof(argv[2])*m0;
	   break;
  case 's':
	   state=atoi(argv[2]);
	   break;
  case 'S':
	   strcpy(State,argv[2]);
           if(!strcmp(State,"1s"))S=9;
           if(!strcmp(State,"2s"))S=2;
           if(!strcmp(State,"2px"))S=3;
           if(!strcmp(State,"2pz"))S=4;
           switch(S)
           {
	    case 2 :break;
	    case 3 :break;
	    case 4 :break;
	    case 9 :S=1;break;
	    default:
		    printf("The `%s' impurity level wave function is not defined\n",State);
		    exit(0);
		    break;
	   }
	   break;
  case 'p':
	   p=*argv[2];
	   switch(p)
	   {
	    case 'e': break;
	    case 'h': break;
	    case 'l': break;
	    default:  printf("Usage:  i0 [-p particle (e, h, or l)]\n");
                      exit(0);
	   }
	   break;
  default :
           printf("Usage:  i0 [-e relative permittivity \033[1m13.18\033[0m]\n");
	   printf("           [-m mass (\033[1m0.067\033[0mm0)]\n");
	   printf("           [-s subband (\033[1m1\033[0m)][-S impurity level (\033[1m1s\033[0m)]\n");
	   exit(0);
 }
 argv++;
 argv++;
 argc--;
 argc--;
}


  V=read_v(&n);			/* reads potential file		*/
  wf=read_wf(n,state,p);	/* reads wavefunction into memory	*/

  dz=read_delta_z(V);		/* z- (growth) direction step length	*/

  lambda_0=4*pi*epsilon*(hbar/e_0)*(hbar/e_0)/m;/* Bohr	theory (1s)	*/

  /* Open files for output of data */

  fe=fopen("e.r","w");           /* E versus r_i	*/
  fl=fopen("l.r","w");           /* lambda versus r_i	*/

  /* Different donor (or acceptor) positions */

  if((fr_i=fopen("r_i.r","r"))==0)    /* open r_i data for reading	*/
  {
   fprintf(stderr,"Error: Cannot open input file 'r_i.r'!\n");
   exit(0);
  }

  /* read in each donor position from r_i.r and perform variational
     calculation for each one						*/

  while(fscanf(fr_i,"%lf\n",&r_i)!=EOF)
  {
   lambda=lambda_0;				/* initial lambda guess	*/
   if((S==2)||(S==3)||(S==4))lambda*=2;		/* for n=2 levels	*/

   lambda_step=lambda/100;

   /* Newton-Raphson iteration for solution of lambda, this occurs when
      dE/dlambda=0, hence the function f is dE/dlambda and f'=d2E/dlambda^2
    								*/

   do
   {
    y1=Energy(wf,V,dz,epsilon,m,lambda-lambda_step,lambda_0,r_i,n,S);
    y2=Energy(wf,V,dz,epsilon,m,lambda,lambda_0,r_i,n,S);
    y3=Energy(wf,V,dz,epsilon,m,lambda+lambda_step,lambda_0,r_i,n,S);

    f=(y3-y1)/(2*lambda_step);
    fdash=(y3-2*y2+y1)/(lambda_step*lambda_step);

    lambda-=f/fdash;

    printf("r_i %4.2f A lambda %4.2f A energy %4.3f meV\n",
            r_i/1e-10,lambda/1e-10,y2/(1e-3*e_0));

   }while(fabs(f/fdash)>1e-10);	/* obtain lambda_0 to 1 Angstrom	*/


   /* Output total energy (E) of impurity/heterostructure system 
      and Bohr radii (lambda), in meV and Angstrom respectively */

   fprintf(fe,"%le %le\n",r_i/1e-10,y2/(1e-3*e_0));
   fprintf(fl,"%le %le\n",r_i/1e-10,lambda/1e-10);

   i_i++;            /* index for impurity atoms */

  }/* end while r_i */

  fclose(fr_i);
  fclose(fe);
  fclose(fl);
  free(V);
  free(wf);


} /* end main */





double
Energy(wf,V,dz,epsilon,m,lambda,lambda_0,r_i,n,S)

data11	*wf;
data11	*V;
double	dz;
double	epsilon;
double	m;
double	lambda;
double	lambda_0;
double	r_i;
int	n;
int	S;

/* This function calculates the expectation value (the energy) of the
   Hamiltonian operator	*/

{
 double	dx;
 double	dy;
 double	Psi();	/* the wave function		*/
 double	Psixyz;	/* Psi(x,y,z)			*/
 double	d2Pdx2;	/* 2nd derivative of Psi wrt x	*/
 double	d2Pdy2;	/* 2nd derivative of Psi wrt y	*/
 double	d2Pdz2;	/* 2nd derivative of Psi wrt z	*/
 double	t=0;	/* intermediate value of `top'	*/
 double	b=0;	/* intermediate value of `bot'	*/
 double	top=0;	/* <Psi|H|Psi>			*/
 double	bot=0;	/* <Psi|Psi>			*/
 double	r;	/* distance from impurity	*/
 double	v;	/* the potential V(z)		*/
 double	x;	/* the spatial coordinate x	*/
 double	y;	/* the spatial coordinate y	*/
 int	iz;	/* index over z coordinates	*/
 
 dx=lambda_0/10;
 dy=lambda_0/10;

 for(iz=1;iz<(n-1);iz++)	/* integration along z-axis	*/
 {
  v=(V+iz)->b;			/* the potential V(z)		*/
  				
  /* integrate over the plane, include `+dx/2' to avoid `x=0'	*/

  for(x=-3*lambda_0+dx/8;x<(3*lambda_0);x+=dx)	
  {
   t=0;b=0;
   for(y=-3*lambda_0+dy/8;y<(3*lambda_0);y+=dy)	
   {
    Psixyz=Psi(wf,lambda,x,y,r_i,iz,S);	/* reused, so remember	*/

    /* Calculate the second derivatives along x, y and z	*/

    d2Pdx2=(Psi(wf,lambda,x+dx,y,r_i,iz,S)-2*Psixyz+Psi(wf,lambda,x-dx,y,r_i,iz,S))/(dx*dx);
    d2Pdy2=(Psi(wf,lambda,x,y+dy,r_i,iz,S)-2*Psixyz+Psi(wf,lambda,x,y-dy,r_i,iz,S))/(dy*dy);
    d2Pdz2=(Psi(wf,lambda,x,y,r_i,iz+1,S)-2*Psixyz+Psi(wf,lambda,x,y,r_i,iz-1,S))/(dz*dz);
	   
    /* Need distance from impurity for Coloumb term		*/

    r=sqrt(x*x+y*y+(((wf+iz)->a)-r_i)*(((wf+iz)->a)-r_i));

    t+=Psixyz*(-(hbar/(2*m))*hbar*(d2Pdx2+d2Pdy2+d2Pdz2)
	+(v-e_0*e_0/(4*pi*epsilon*r))*Psixyz)*dy;
    b+=Psixyz*Psixyz*dy;		/* variable step requires this	*/

    //if(y>(lambda_0/10))dy=lambda_0/100;	/* variable integration step	*/
    //if(y>lambda_0)dy=lambda_0/10;
   }
   //if(x>(lambda_0/10))dx=lambda_0/100;	/* variable integration step	*/
   //if(x>lambda_0)dx=lambda_0/10;		
   top+=t*dx;bot+=b*dx;			/* variable step requires this	*/
  }
 }

 return(top/bot);
}



double
Psi(wf,lambda,x,y,r_i,iz,S)

/* The wave function psi(z)phi(r)					*/

data11	*wf;
double	lambda;
double	x;
double	y;
double	r_i;
int	iz;
int	S;

{
 double	r;
 r=sqrt(x*x+y*y+(((wf+iz)->a)-r_i)*(((wf+iz)->a)-r_i));

 if(S==1)return(((wf+iz)->b)*exp(-r/lambda));			/* 1s	*/
 if(S==2)return(((wf+iz)->b)*(1-r/lambda)*exp(-r/lambda));	/* 2s	*/
 if(S==3)return(((wf+iz)->b)*x*exp(-r/lambda));			/* 2px	*/
 if(S==4)return(((wf+iz)->b)*(((wf+iz)->a)-r_i)*exp(-r/lambda));/* 2pz	*/
}



double 
read_delta_z(Vp)

/* This function calculates the separation along the z (growth) 
   direction of the user supplied potentials                                        */

data11 *Vp;
{
 double z[2];           /* displacement along growth direction     */

 z[0] = Vp->a;
 Vp++;
 z[1] = Vp->a;
 return(z[1]-z[0]);
}



data11
*read_v(n)

/* This function reads the potential into memory and returns the start
   address of this block of memory and the number of lines         */

int     *n;

{
 FILE   *fp;            /* file pointer to potential file          */
 data11  *Vp;           /* temporary pointer to potential          */
 data11  *V;            /* start address of potential              */

 if((fp=fopen("v.r","r"))==0)
 {
   fprintf(stderr,"Error: Cannot open input file 'v.r'!\n");
   exit(0);
 }
 *n=0;
 while(fscanf(fp,"%*le %*le")!=EOF)
  (*n)++;
 rewind(fp);


 V=(data11 *)calloc(*n,sizeof(data11));
 if (V==0)  {
  fprintf(stderr,"Cannot allocate memory!\n");
  exit(0);
 }
 Vp=V;

 while(fscanf(fp,"%le %le", &(Vp->a), &(Vp->b))!=EOF)
  Vp++;

 fclose(fp);
 return(V);

}



data11
*read_wf(n,state,p)

/* This function reads the potential into memory and returns the start
   address of this block of memory and the number of lines	   */

int	n;
int	state;
char	p;

{
 char	filename[9];	/* input filename			*/
 int	i;		/* index				*/
 FILE 	*Fwf;		/* file pointer to wavefunction file	*/
 data11	*wf;		/* pointer to wave function structure	*/

 /* Open wave function file for reading	*/
 
 sprintf(filename,"wf_%c%i.r",p,state);	
 if((Fwf=fopen(filename,"r"))==0)
 {fprintf(stderr,"Error: Cannot open input file '%s'!\n",filename);exit(0);}

 /* Allocate memory for wavefunction */

 wf=(data11 *)calloc(n,sizeof(data11));
 if(wf==0){fprintf(stderr,"Cannot allocate memory!\n");exit(0);}

 /* Read in wave function 	*/

 i=0;						/* Read in each file	*/
 while(fscanf(Fwf,"%le %le",&((wf+i)->a),&((wf+i)->b))!=EOF)
  i++;

 fclose(Fwf);					/* Close each file	*/

 return(wf);

}

