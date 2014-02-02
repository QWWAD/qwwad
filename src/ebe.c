/*==================================================================
                                 ebe
  ==================================================================*/

/* This program calculates the exciton binding energy from 
   user supplied uncorrelated (or correlated) wavefunctions,
   usually read in from the files wf_pX.r.

   Input files:
     wf_eX.r      electron wavefunction versus z  
     wf_hX.r      hole wavefunction versus z  

   Output files:
     ABC.r         terms A, B, Ct and Cv versus lambda
     beta.r        beta corresponding to minimum Eb versus lambda
     EX0.r         minimum binding energy, corresponding lambda and beta
     EX0-lambda.r  minimum binding energy versus lambda
     p.r           uncorrelated probaility of e-h separation

   Paul Harrison, March 1994                                  */

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

typedef
struct	{
 double	z;		    /* z value of files    		  */
 double	wf[2];		    /* electron and hole values   	  */
} files;
 
typedef
struct  {
 double a;                  /* electron and hole separation       */
 double p;                  /* probability of separation by a     */
} probs;

main(int argc,char *argv[])
{
double Eb_1S();             /* exciton binding energy of 1S      */
double read_delta_z();
boolean repeat_beta();      /* repeat beta variational loop      */
boolean repeat_lambda();    /* repeat lambda variational         */
files  *read_data();        /* reads data from external files    */
probs  *pP_calc();          /* calculates p(a), Pm(a) and
                                 Pmu(a)                          */
void   input();

double beta;                /* dimensionality parameter          */
double beta_start;          /* initial beta                      */
double beta_step;           /* beta increment                    */
double beta_stop;           /* final beta                        */
double beta_0;              /* beta for Eb_min                   */
double beta_0_lambda;       /* beta for Eb_min_beta              */
double delta_a;             /* separation of adjacent a values   */
double delta_z;             /* z separation of input potentials  */
double Eb;                  /* exciton binding energy (<0=bound) */
double Eb_min;              /* minimum Eb for lambda variation   */
double Eb_min_beta;         /* minimum Eb for beta variation     */
double epsilon;             /* relative permittivity of material */
double lambda;              /* Bohr radius                       */
double lambda_start;        /* initial Bohr radius               */
double lambda_step;         /* lambda increment                  */
double lambda_stop;         /* final lambda                      */
double lambda_0;            /* lambda for Eb_min                 */
double m[2];		    /* e and h z-axis masses		 */
double m_xy[2];		    /* e and h x-y plane masses		 */
double mu_xy;		    /* exciton reduced mass in x-y plane */
int    n;		    /* length of potential file		 */
int    N_x;                 /* number of points in x integration */
int    state[2];	    /* electron and hole states          */
boolean output_flag;        /* if set, write data to screen      */
boolean repeat_flag_beta;   /* repeat variational beta loop flag */
boolean repeat_flag_lambda; /* repeat variational lambda loop    */
char   c;                   /* general character                 */
FILE   *FABC;               /* file pointer to ABC.r             */
FILE   *Fbeta;              /* file pointer to beta.r            */
FILE   *FEX0l;              /* file pointer to EX0-lambda.r      */
FILE   *FEX0;               /* file pointer to EX0.r             */
files  *data_start;    	    /* start address of wavefunction	 */
probs  *pP_start;           /* start address of probabilities    */

/* default values */

state[0]=1;
state[1]=1;
beta_start=0.001;
beta_step=0.05;
beta_stop=-1.0;
epsilon=13.18*epsilon_0;
lambda_start=70e-10;
lambda_step=1e-10;
lambda_stop=-1e-10;
m[0]=0.067*m0;
m[1]=0.62*m0;
output_flag=false;
repeat_flag_beta=true;
repeat_flag_lambda=true;

/* computational default values */

N_x=100;

while((argc>1)&&(argv[1][0]=='-'))
{
 switch(argv[1][1])
 {
  case 'a':
	   state[0]=atoi(argv[2]);
	   break;
  case 'b':
	   state[1]=atoi(argv[2]);
	   break;
  case 'e':
           epsilon=atof(argv[2])*epsilon_0;
           break;
  case 'm':
           m[0]=atof(argv[2])*m0;
           break;
  case 'n':
           m[1]=atof(argv[2])*m0;
           break;
  case 'N':
	   N_x=atoi(argv[2]);
	   break;
  case 'o':
	   output_flag=true;
 	   argv--;
	   argc++;
           break;
  case 's':
           lambda_start=atof(argv[2])*1e-10;
           break;
  case 't':
           lambda_step=atof(argv[2])*1e-10;
           break;
  case 'u':
           lambda_stop=atof(argv[2])*1e-10;
           break;
  case 'w':
           beta_start=atof(argv[2]);
           break;
  case 'x':
           beta_step=atof(argv[2]);
           break;
  case 'y':
           beta_stop=atof(argv[2]);
           break;
  default :
           printf("Usage:  ebe [-a # electron eigenstate \033[1m1\033[0m][-b # hole eigenstate \033[1m1\033[0m]\n");
           printf("            [-e relative permittivity \033[1m13.18\033[0m]\n");
           printf("            [-m electron mass (\033[1m0.067\033[0mm0)][-n hole mass (\033[1m0.62\033[0mm0)]\n");
           printf("            [-N # points in w integration \033[1m100\033[0m][-o output data to screen \033[1mfalse\033[0m]\n");
           printf("            [-s starting lambda (\033[1m70\033[0mA)][-t lambda increment (\033[1m1\033[0mA)]\n");
           printf("            [-u final lambda (\033[1m-1\033[0mA)]\n");
           printf("            [-w starting beta \033[1m0.001\033[0m][-x beta increment \033[1m0.05\033[0m][-y final beta \033[1m-1\033[0m]\n");
           exit(0);
 }
 argv++;
 argv++;
 argc--;
 argc--;
}

data_start=read_data(state,&n);	/* reads wave functions */

delta_z=read_delta_z(data_start);
delta_a=delta_z;
Eb_min=1*e_0;			/* i.e. 1eV ! */
m_xy[0]=m[0];	m_xy[1]=m[1];	/* assumes isotropic mass for now	*/
mu_xy=1/(1/m_xy[0]+1/m_xy[1]);	/* calculate reduced mass in-plane	*/

FABC=fopen("ABC.r","w");
Fbeta=fopen("beta-lambda.r","w");
FEX0l=fopen("EX0-lambda.r","w");

pP_start=pP_calc(delta_z,n,data_start); /* calculates p and P's, returns start
                                           address of structure               */
 
if(output_flag)printf("  l/A   beta   Eb/meV  T/meV  V/meV   OS/arb.\n");

lambda=lambda_start;
do
{
 beta=beta_start;
 Eb_min_beta=1*e_0;           /* i.e. 1eV ! */
 do
 {

Eb=Eb_1S(data_start,pP_start,FABC,beta,delta_a,epsilon,lambda,m,mu_xy,N_x,n,output_flag);

  repeat_flag_beta=repeat_beta(&beta,&beta_0_lambda,&Eb,&Eb_min_beta);

  beta+=beta_step;

 }while((repeat_flag_beta&&(beta_stop<0))||(beta<beta_stop));

 fprintf(FEX0l,"%lf %lf\n",lambda/1e-10,Eb_min_beta/(1e-3*e_0));
 fprintf(Fbeta,"%lf %lf\n",lambda/1e-10,beta_0_lambda);

 repeat_flag_lambda=repeat_lambda(&beta_0,&beta_0_lambda,&Eb_min_beta,&Eb_min,
                                  &lambda,&lambda_0);

 lambda+=lambda_step;   /* increment Bohr radius */

}while((repeat_flag_lambda&&(lambda_stop<0))||(lambda<lambda_stop));
  

/* Write out final data to file	*/

FEX0=fopen("EX0.r","w");
fprintf(FEX0,"%6.3lf %6.2lf %6.3lf\n",Eb_min/(1e-3*e_0),lambda_0/1e-10,beta_0);
fclose(FEX0);

fclose(FABC);
fclose(Fbeta);
fclose(FEX0l);

free(data_start);
free(pP_start);


} /* end main */



double
Eb_1S(fdata,ppP,FABC,beta,delta_a,epsilon,lambda,m,mu_xy,N_x,n,output_flag)

files  *fdata;              /* pointer to data structure         */
probs  *ppP;                /* pointer to pP structure           */
FILE   *FABC;               /* file pointer to ABC.r             */
double beta;
double delta_a;
double epsilon;
double lambda;
double m[];
double mu_xy;
int    N_x;
int    n;
boolean	output_flag;
{
 double F();                 /* F(a)---see notes!                 */
 double G();                 /* G(a)---see notes!                 */
 double J();                 /* J(a)---see notes!                 */
 double K();                 /* K(a)---see notes!                 */
 double A;                   /* {\cal A}, see notes!              */
 double B;                   /* {\cal B}, see notes!              */
 double Ct;                  /* kinetic energy component of C     */
 double Cv;                  /* potential energy component of C   */
 double D;                   /* {\cal D}, see notes!              */
 double Eb;                  /* exciton binding energy (<0=bound) */
 double O;                   /* {\cal O}, overlap integral        */
 int    i_a;                 /* index through a values            */
 
 A=0;B=0;Ct=0;Cv=0;D=0;O=0;         /* initialize variables before integration */
 for(i_a=0;i_a<n;i_a++)
 {
  A+=(ppP->p)*G(ppP->a,beta,lambda,N_x)*delta_a;
  B+=(ppP->p)*G(ppP->a,beta,lambda,N_x)*delta_a;
  Ct+=(ppP->p)*J(ppP->a,beta,lambda,N_x)*delta_a;
  Cv+=(ppP->p)*K(ppP->a,beta,lambda,N_x)*delta_a;
  D+=(ppP->p)*F(ppP->a,beta,lambda)*delta_a;
  
  O+=(fdata->wf[0])*(fdata->wf[1])*delta_a;  /* strictly delta_z, but = */

  ppP++;
  fdata++;
 }

 A*=hbar*hbar/(2*m[0]);         /* multiply integrals by constant factors */
 B*=hbar*hbar/(2*m[1]);
 Ct*=-hbar*hbar/(2*mu_xy);
 Cv*=-e_0*e_0/(4*pi*epsilon);

 fprintf(FABC,"%6.2lf %6.3lf %6.3lf %6.3lf %6.3lf\n",lambda/1e-10,
         A/D/(1e-3*e_0),B/D/(1e-3*e_0),Ct/D/(1e-3*e_0),Cv/D/(1e-3*e_0));

 Eb=(A+B+Ct+Cv)/D;

 if(output_flag)
  printf("%6.2lf %6.3lf %6.3lf %6.3lf %6.3lf %6.3le\n",lambda/1e-10,beta,
         Eb/(1e-3*e_0),(A+B+Ct)/D/(1e-3*e_0),Cv/D/(1e-3*e_0),sqr(O)/D);

 return(Eb);
}



boolean
repeat_beta(beta,beta_0_lambda,Eb,Eb_min_beta)

double *beta;
double *beta_0_lambda;
double *Eb;
double *Eb_min_beta;
{
 boolean flag;

 if(*Eb<*Eb_min_beta)
 {
  *Eb_min_beta=*Eb;
  *beta_0_lambda=*beta;
  flag=true;
 }
 else
 {
  flag=false;
 }

 return(flag);
}



boolean
repeat_lambda(beta_0,beta_0_lambda,Eb_min_beta,Eb_min,lambda,lambda_0)

double *beta_0;
double *beta_0_lambda;
double *Eb_min_beta;
double *Eb_min;
double *lambda;
double *lambda_0;
{
 boolean flag;

 if(*Eb_min_beta<*Eb_min)
 {
  *Eb_min=*Eb_min_beta;
  *lambda_0=*lambda;
  *beta_0=*beta_0_lambda;
  flag=true;
 }
 else
 {
  flag=false;
 }

 return(flag);
}



double 
F(a,beta,lambda)

/* This function returns the value of F(a) */

double a;
double beta;
double lambda;
{
 double f;
 f=2*pi*lambda*(sqrt(1-sqr(beta))*a/2+lambda/4)*
   exp(-2*sqrt(1-sqr(beta))*a/lambda);

 return(f);
}



double 
G(a,beta,lambda,N_x)

/* This function returns the value of G(a),
   to overcome the problem of divergence when
   x=0, the integration is performed using a
   midpoint sum, the strip width being delta_x */

double a;
double beta;
double lambda;
int    N_x;
{
 double delta_x;
 double g;
 double x;  /* dummy variable---see notes! */

 delta_x=(1.0-0.0)/(float)N_x;

 g=0;       /* initialize variable */

 for(x=delta_x/2;x<1;x+=delta_x)
 {
  g+=1/(1/x+x)*exp(-sqrt(1-sqr(beta))*a*(1/x+x)/lambda)
     *(1/sqr(x)-1)*delta_x;
 }
 g*=2*pi*sqr(1-sqr(beta))*sqr(a)/sqr(lambda);

 return(g);
}



double 
J(a,beta,lambda,N_x)

/* This function returns the value of J(a),
   to overcome the problem of divergence when
   x=0, the integration is performed using a
   midpoint sum, the strip width being delta_x */

double a;
double beta;
double lambda;
int    N_x;
{
 double delta_x;
 double j13;     /* J1+J3---see notes! */
 double j24;     /* J2+J4---see notes! */
 double x;       /* dummy variable---see notes! */

 delta_x=(1.0-0.0)/(float)N_x;

 j13=2*pi*(sqrt(1-sqr(beta))*a/(2*lambda)-0.25)
     *exp(-2*sqrt(1-sqr(beta))*a/lambda);

 j24=0;        /* initialize variable */

 for(x=delta_x/2;x<1;x+=delta_x)
 {
  j24+=(
        -1/(lambda*sqr(1/x+x)/4)
        -sqrt(1-sqr(beta))*a/(sqr(lambda)*(1/x+x)/2)
       )
       *exp(-sqrt(1-sqr(beta))*a*(1/x+x)/lambda)
       *(1/sqr(x)-1)*delta_x;
 }
 j24*=2*pi*sqrt(1-sqr(beta))*a/2;

 return(j13+j24);
}



double 
K(a,beta,lambda,N_x)

/* This function returns the value of K(a),
   to overcome the problem of divergence when
   x=0, the integration is performed using a
   midpoint sum, the strip width being delta_x */

double a;
double beta;
double lambda;
int N_x;
{
 double delta_x;      /* step length of integration */
 double k;    
 double lower_limit;  /* lower limit of integration */
 double upper_limit;  /* upper limit of integration */
 double x;            /* dummy variable---see notes! */

 upper_limit=(1-sqrt(1-sqr(beta)))/beta;
 lower_limit=0;

 delta_x=(upper_limit-lower_limit)/(float)N_x;

 k=0;

 for(x=lower_limit+delta_x/2;x<upper_limit;x+=delta_x)
 {
  k+=exp(-beta*a*(1/x-x)/lambda)
     *(1/sqr(x)-1)*delta_x;
 }
 k*=2*pi*beta*a/2;

 return(k);
}


double 
read_delta_z(Vp)

/* This function calculates the separation along the z (growth) 
   direction of the user supplied wave functions, assumes regular
   one-dimensional mesh                                        */

files *Vp;
{
 double z[2];           /* displacement along growth direction     */

 z[0] = Vp->z;
 Vp++;
 z[1] = Vp->z;
 return(z[1]-z[0]);
}





files
*read_data(state,n)

/* This function reads the potential into memory and returns the start
   address of this block of memory and the number of lines	   */

int	state[];
int	*n;

{
 char	filenamee[9];	/* filename of electron wave function		*/
 char	filenameh[9];	/* filename of hole wave function		*/
 FILE 	*Fwfe;		/* file pointer to electron wavefunction file	*/
 FILE 	*Fwfh;		/* file pointer to hole wavefunction file	*/
 files  *ft;		/* temporary pointer to wavefunction		*/
 files  *data_start;	/* start address of wavefunction		*/

 /* Generate electron and hole filenames	*/

 sprintf(filenamee,"wf_e%i.r",state[0]);
 sprintf(filenameh,"wf_h%i.r",state[1]);

 /* open files	*/

 if((Fwfe=fopen(filenamee,"r"))==0)
  {fprintf(stderr,"Error: Cannot open input file '%s'!\n",filenamee);exit(0);}

 if((Fwfh=fopen(filenameh,"r"))==0)
  {fprintf(stderr,"Error: Cannot open input file '%s'!\n",filenameh);exit(0);}

 /* count number of lines in wave function files	*/

 *n=0;
 while(fscanf(Fwfe,"%*le %*le")!=EOF)
  (*n)++;
 rewind(Fwfe);

 /* allocate memory for wave functions	*/

 data_start=(files *)calloc(*n,sizeof(files));
 if (data_start==0)  
 {
  fprintf(stderr,"Cannot allocate memory!\n");
  exit(0);
 }
 ft=data_start;

 /* read wave functions into memory	*/

 while(fscanf(Fwfe,"%le %le",&(ft->z),&(ft->wf[0]))!=EOF)
 {
  fscanf(Fwfh,"%*le %le",&(ft->wf[1]));
  ft++;
 }

 fclose(Fwfe);
 fclose(Fwfh);

 return(data_start);

}



probs
*pP_calc(delta_a,n,data_start)

/* This function calculates the probabilities known as p(a), Pm(a)
   and Pmu(a) and returns the start address of the structure */

double  delta_a;        /* distance between adjacent points        */
int	n;              /* number of lines in wavefunction file    */
files   *data_start;    /* pointer to beginning of wavefunctions   */

{
 double temp_sum;
 int    i;
 int    j;
 FILE   *Fp;            /* file pointer to p(a) versus a, p.r      */
 files  *fdata;		/* pointer to wavefunctions                */
 probs  *pP_start;	/* start address of pP                     */
 probs  *ppP;           /* start address of pP                     */

 pP_start=(probs *)calloc(n,sizeof(probs));
 if (pP_start==0)  
 {
  fprintf(stderr,"Cannot allocate memory!\n");
  exit(0);
 }

 /* Note integrals are calculated with a simple sum of left hand
   ordinate times strip width.  Hence the total number of ordinates
   is simply n.                                                    */
  
 temp_sum=0;
 ppP=pP_start;
 for(i=0;i<n;i++)
 {
  ppP->a=delta_a*(float)i;
  ppP->p=0;                /* initialize variable to zero */

  fdata=data_start;      /* reset wavefunction pointer  */
  for(j=0;j<n-i;j++)
  {
   ppP->p+=(sqr((fdata+i)->wf[0])*sqr(fdata->wf[1])
            +sqr(fdata->wf[0])*sqr((fdata+i)->wf[1])
           )*delta_a;
   fdata++;
  }
  temp_sum+=ppP->p*delta_a;
  ppP++;
 }

 Fp=fopen("p.r","w");
 ppP=pP_start;
 for(i=0;i<n;i++)
 {
  fprintf(Fp,"%20.17le %20.17le\n",ppP->a,ppP->p);
  ppP++;
 }
 fclose(Fp);

 /*printf("temp_sum=%20.17le\n",temp_sum);*/
 return(pP_start);
}

