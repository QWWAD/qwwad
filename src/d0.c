/*==================================================================
                                d0  
  ==================================================================*/

/* This program implements a variational technique to calculate the
   uncorrelated one particle energies of an electron attatched to a 
   single donor at any position, in any user supplied potential.  
   The potential is read from the file v.r

   This version is a two variational parameter calculation---the variable
   symmetry trial wavefunction:

		Psi=chi(z) exp(-r"/lambda)

   where	r"=sqrt(x^2+y^2+zeta^2(z-r_d)^2)

		Input files:
		r_d.r		donor (or acceptor positions)
		v.r		one-dimensional potential

		Output files:
		e.r		total energies for each r_d
		l.r		Bohr radii (lambda) for each r_d
		wfn.r		wave functions, both Psi and chi, n=0,1,2..
		zeta.r		zeta value for each r_d

   Paul Harrison, February 1993 

   Substantially revised,

   Paul Harrison, February 1998					*/

#include <stdbool.h>
#include <stdio.h>
#include <math.h>
#include <signal.h>
#include <gsl/gsl_math.h>
#include "d0-helpers.h"
#include "maths.h"
#include "qclsim-constants.h"

static double I_1(const double lambda,
                  const double z_dash,
                  const double zeta);
static double I_2(const double lambda,
                  const double z_dash,
                  const double zeta);
static double I_3(const double lambda,
                  const double z_dash,
                  const double zeta,
                  const int    N_w);
static double I_4(const double lambda,
                  const double z_dash,
                  const double zeta,
                  const int    N_w);

int main(int argc,char *argv[])
{
double psi_at_inf();
double V_min();
void   wavefunctions();
bool    repeat_lambda();    
bool    repeat_zeta();    

double d_E;                 /* infinitesmal energy               */
double delta_E;             /* small but finite energy           */
double delta_z;             /* z separation of input potentials  */
double dy;                  /* derivative of function            */
double E;                   /* electron (or hole) energies       */
double epsilon;             /* permitivity of material           */
double lambda;              /* Bohr radius (variational)         */

/* TODO: lambda_0 is found iteratively. Check that this is a sensible initial value */
double lambda_0 = 0;        /* Bohr radius of electron (or hole) */
double lambda_start;        /* initial Bohr radius               */
double lambda_step;         /* Bohr radius increment             */
double lambda_stop;         /* final lambda                      */
double mstar;               /* electron mass  	                 */
double r_d;                 /* donor (or acceptor) position      */
double x;                   /* independent variable (energy)     */
double x_min;               /* variational calculation minimum x */
double x_min_zeta;          /* minimum x as zeta varies only     */
double y;                   /* function (psi at infinity)        */
double y1;                  /* temporary y value                 */
double y2;                  /* temporary y value                 */
double zeta;                /* value for zeta in wavefunction    */

/*
 * TODO: zeta_0 and zeta_0_lambda are found iteratively. 
 *       Check that this is a sensible initial value
 */
double zeta_0 = 0;          /* zeta of solution                  */
double zeta_0_lambda = 0;   /* zeta of minimum for each lambda   */
double zeta_start;          /* initial zeta                      */
double zeta_step;           /* zeta increment                    */
double zeta_stop;           /* final zeta                        */
int    i_d;                 /* donor (or acceptor) index         */
int    N_w;                 /* number of strips in w integration */
size_t n;		    /* number of lines of potential file */
bool   repeat_flag_zeta;   /* variational flag=>new zeta        */
bool   repeat_flag_lambda; /* variational flag=>new lambda      */
data11  *Vstart;             /* start address of potential        */
FILE   *fe;                 /* file pointer for energies         */
FILE   *fl;                 /* file pointer for lambda_0         */
FILE   *fzeta;              /* file pointer for zeta_0           */
FILE   *fr_d;               /* file pointer to donor positions   */
 

/* default values */

delta_E=1e-3*e;
epsilon=13.18*eps0;
lambda_start=50.0e-10;
lambda_step=1.0e-10;
lambda_stop=-1.0e-10;
zeta_start=0.001;
zeta_step=0.1;
zeta_stop=-1.0;
mstar=0.067*me;

/* computational default values */

i_d=0;          
d_E=delta_E/1e+6;
N_w=100;

while((argc>1)&&(argv[1][0]=='-'))
{
 switch(argv[1][1])
 {
  case 'd':
	   delta_E=atof(argv[2])*1e-3*e;
	   break;
  case 'e':
	   epsilon=atof(argv[2])*eps0;
	   break;
  case 'm':
	   mstar=atof(argv[2])*me;
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
	   zeta_start=atof(argv[2]);
	   break;
  case 'x':
	   zeta_step=atof(argv[2]);
	   break;
  case 'y':
	   zeta_stop=atof(argv[2]);
	   break;
  default :
           printf("Usage:  d0 [-d energy step (\033[1m1\033[0mmeV)][-e relative permittivity \033[1m13.18\033[0m]\n");
	   printf("           [-m mass (\033[1m0.067\033[0mm0)]\n");
	   printf("           [-s starting lambda (\033[1m50\033[0mA)][-t lambda increment (\033[1m1\033[0mA)]\n");
	   printf("           [-u final lambda (\033[1m-1\033[0mA)]\n");
	   printf("           [-w starting zeta \033[1m0.001\033[0m][-x zeta increment \033[1m0.1\033[0m][-y final zeta \033[1m-1\033[0m]\n");
	   exit(0);
 }
 argv++;
 argv++;
 argc--;
 argc--;
}

  Vstart = read_v(&n);                  /* reads potential file */

  delta_z=read_delta_z(Vstart);

  /* Open files for output of data */

  fe=fopen("e.r","w");           /* E versus r_d	*/
  fl=fopen("l.r","w");           /* lambda versus r_d	*/
  fzeta=fopen("zeta.r","w");     /* zeta versus r_d	*/

  /* Different donor (or acceptor) positions */

  if((fr_d=fopen("r_d.r","r"))==0)    /* open r_d data for reading	*/
  {
   fprintf(stderr,"Error: Cannot open input file 'r_d.r'!\n");
   exit(0);
  }

  /* read in each donor position from r_d.r and perform variational
     calculation for each one						*/

  while(fscanf(fr_d,"%lf\n",&r_d)!=EOF)
  {
   lambda=lambda_start;    /* initial lambda value               */
   x_min=e;              /* minimum energy of single donor 1eV */

   /* Variational calculation */

   do
   {
    zeta=zeta_start;
    x_min_zeta=1*e;
    do
    {

    /* initial energy estimate=minimum potential-binding energy
                               of particle to free ionised dopant */

    x=V_min(Vstart,n)-gsl_pow_2(e)/(4*pi*epsilon*lambda);   

    /* increment energy-search for f(x)=0 */

    y2=psi_at_inf(x,delta_z,epsilon,lambda,mstar,r_d,Vstart,n,zeta,N_w);

    do
    {
     y1=y2;
     x+=delta_E;
     y2=psi_at_inf(x,delta_z,epsilon,lambda,mstar,r_d,Vstart,n,zeta,N_w);
    }while(y1*y2>0);

   /* improve estimate using midpoint rule */

    x-=fabs(y2)/(fabs(y1)+fabs(y2))*delta_E;

   /* implement Newton-Raphson method */

    do
    {
     y=psi_at_inf(x,delta_z,epsilon,lambda,mstar,r_d,Vstart,n,zeta,N_w);
     dy=(psi_at_inf(x+d_E,delta_z,epsilon,lambda,mstar,r_d,Vstart,n,zeta,N_w)-
         psi_at_inf(x-d_E,delta_z,epsilon,lambda,mstar,r_d,Vstart,n,zeta,N_w))/
        (2.0*d_E);
     x-=y/dy;

    }while(fabs(y/dy)>1e-9*e);

    printf("r_d %le lambda %le zeta %le energy %le meV\n",
            r_d,lambda,zeta,x/(1e-3*e));       

    repeat_flag_zeta=repeat_zeta(zeta,&zeta_0_lambda,x,&x_min_zeta);
    zeta+=zeta_step;		/* increments zeta	*/
    }while((repeat_flag_zeta&&(zeta_stop<0))||(zeta<zeta_stop));

    repeat_flag_lambda=repeat_lambda(lambda,&lambda_0,x_min_zeta,&x_min,&zeta_0,zeta_0_lambda);

    lambda+=lambda_step;     /* increments Bohr radius */

   }while((repeat_flag_lambda&&(lambda_stop<0))||(lambda<lambda_stop));

   E=x_min;   /* assign the energy E to the minimum x_min */

   /* Output neutral dopant binding energies (E) and 
      Bohr radii (lambda) in meV and Angstrom respectively */

   fprintf(fe,"%le %le\n",r_d/1e-10,E/(1e-3*e));
   fprintf(fl,"%le %le\n",r_d/1e-10,lambda_0/1e-10);
   fprintf(fzeta,"%le %le\n",r_d/1e-10,zeta_0);

   wavefunctions(delta_z,E,epsilon,lambda_0,mstar,r_d,
                 zeta_0,i_d,Vstart,n,N_w);

   i_d++;            /* donor index */

  }/* end while r_d */

  fclose(fr_d);
  fclose(fe);
  fclose(fl);
  fclose(fzeta);
  free(Vstart);

  return EXIT_SUCCESS;
} /* end main */


/**
 * Compares minimum value of energy for this lambda,
 * if it is a new true minima, then repeat for new lambda
 *
 * \param[in]  lambda        Bohr radius (variational)
 * \param[out] lambda_0      Bohr radius of electron (or hole)
 * \param[in]  x_min_zeta    Minimum x as zeta varies only
 * \param[out] x_min         Variational calculation minimum x
 * \param[out] zeta_0        Zeta of solution
 * \param[in]  zeta_0_lambda Zeta of minimum for each lambda
 *
 * \returns True if x_min_zeta < x_min
 */
bool repeat_lambda(const double  lambda,
                   double       *lambda_0,
                   const double  x_min_zeta,
                   double       *x_min,
                   double       *zeta_0,
                   const double  zeta_0_lambda)
{
 bool flag = false;

 if(x_min_zeta<*x_min)    
 {      
  *x_min=x_min_zeta;                 /* set new minimum        */
  *lambda_0=lambda;
  *zeta_0=zeta_0_lambda;
  flag=true;                          /* repeat for new lambda  */
 }

 return flag;
}

/**
 * compares current energy value with the minima for this
 * value of lambda only, if a new minima for this lambda, then repeat
 */
bool repeat_zeta(const double  zeta,
                 double       *zeta_0_lambda,
                 const double  x,
                 double       *x_min_zeta)
{
 bool flag = false;
 
 if(x<*x_min_zeta)
 {
  *x_min_zeta=x;
  *zeta_0_lambda=zeta;
  flag=true;
 }

 return flag;
}



double
psi_at_inf(E,delta_z,epsilon,lambda,mstar,r_d,Vp,n,zeta,N_w)     

/* This function returns the value of the wavefunction (psi)
   at +infinity for a given value of the energy.  The solution
   to the energy occurs for psi(+infinity)=0.                      */

double E;
double delta_z;
double epsilon;
double lambda;
double mstar;
double r_d;
int n;
double zeta;
data11 *Vp;
int N_w;
{
 double alpha;		     /* coefficient of second derivative, see notes */
 double beta;                /* coefficient of first derivative            */
 double gamma;               /* coefficient of function                    */
 double delta_psi;           /* initial wavefunction value                 */
 double I1,I2,I3,I4;         /* particular values of the functions         */
 double kappa;     
 double psi[3];              /* wavefunction at z-delta_z, z, z+delta_z    */
 int i;			     /* index			                   */


 /* ignore potential V corresponding to psi[0] and psi[1] */

 Vp++;

 /* boundary conditions */

 kappa=sqrt(2*mstar/hBar*(Vp->b-E)/hBar);

 delta_psi=1.0e-10;

 psi[0]=delta_psi;                 /* arbitrary number close to zero */
 psi[1]=psi[0]*exp(kappa*delta_z); /* exponential growth produce by  */

 for(i=1;i<n;i++)
 {
  I1=I_1(lambda,(Vp->a)-r_d,zeta);
  I2=I_2(lambda,(Vp->a)-r_d,zeta);
  I3=I_3(lambda,(Vp->a)-r_d,zeta,N_w);
  I4=I_4(lambda,(Vp->a)-r_d,zeta,N_w);

  alpha=I1;
  beta=2*I2;
  gamma=I3+(2*mstar*gsl_pow_2(e/hBar)/(4*pi*epsilon))*I4
          -(2*mstar/hBar)*((Vp->b)-E)*I1/hBar;

  psi[2]=((-1+beta*delta_z/(2*alpha))*psi[0]
          +(2-gsl_pow_2(delta_z)*gamma/alpha)*psi[1]
         )/(1+beta*delta_z/(2*alpha));

  psi[0]=psi[1];
  psi[1]=psi[2];
  Vp++;
 } /* end for */

 return(psi[0]-delta_psi);
}

/* This function calculates and writes the wavefunctions
   both psi(z) and chi(z) to the external file wf(n).r.		*/
void wavefunctions(const double  delta_z,
                   const double  E,
                   const double  epsilon,
                   const double  lambda,
                   const double  mstar,
                   const double  r_d,
                   const double  zeta,
                   const int     i_d,
                   data11       *Vp,
                   const size_t  n,
                   const int     N_w)
{
 double        alpha;		 /* coefficient of second derivative, see notes */
 double        beta;            /* coefficient of first derivative             */
 double        gamma;           /* coefficient of function                     */
 double        delta_psi;       /* initial wavefunction value                  */
 double        I1,I2,I3,I4;     /* particular values of the functions	        */
 double        kappa;
 double        Npsi=0;          /* normalisation integral for psi              */
 double        Nchi=0;          /* normalisation integral for chi              */
 double        psi[3];          /* wavefunctions at z-d_z,z,z+d_z              */
 unsigned int  i;               /* index                                       */
 char          filename[9];     /* character string for wavefunction filename  */
 FILE         *fw;             /* file wf.r                                   */
 data12       *wf_start;      /* pointer to start of w.f                     */
 data12       *wf;		 /* wavefunction pointer, note b[0]=psi, 
                                 b[1]=chi, see notes                         */

 wf_start=(data12 *)calloc(n,sizeof(data12)); /* allocates memory for wavefunctions */
 if (wf_start==0)  {
  fprintf(stderr,"Cannot allocate memory!\n");
  exit(0);
 }
 wf=wf_start;

 /* boundary conditions */

 kappa=sqrt(2*mstar/hBar*((Vp->b)-E)/hBar);

 delta_psi=1.0e-10;

 psi[0]=delta_psi;                 /* arbitrary number close to zero */
 psi[1]=psi[0]*exp(kappa*delta_z); /* exponential growth produce by  */

 /* write first value of wavefunction */ 

 wf->a=Vp->a;
 wf->b[0]=psi[0]*exp(-zeta*fabs((Vp->a)-r_d)/lambda); 
 wf->b[1]=psi[0]; 

 wf++;
 Vp++;
 wf->a=Vp->a;
 wf->b[0]=psi[1]*exp(-zeta*fabs((Vp->a)-r_d)/lambda); 
 wf->b[1]=psi[1]; 
 
 /* calculate unnormalised wavefunction */

 for(i=2;i<n;i++) /* points 0 and 1 defined	*/
 {
  I1=I_1(lambda,(Vp->a)-r_d,zeta);
  I2=I_2(lambda,(Vp->a)-r_d,zeta);
  I3=I_3(lambda,(Vp->a)-r_d,zeta,N_w);
  I4=I_4(lambda,(Vp->a)-r_d,zeta,N_w);

  alpha=I1;
  beta=2*I2;
  gamma=I3+(2*mstar*gsl_pow_2(e/hBar)/(4*pi*epsilon))*I4
          -(2*mstar/hBar)*((Vp->b)-E)*I1/hBar;

  psi[2]=((-1+beta*delta_z/(2*alpha))*psi[0]
          +(2-gsl_pow_2(delta_z)*gamma/alpha)*psi[1]
         )/(1+beta*delta_z/(2*alpha));


  /* write wavefunction point corresponding to current z and V */

  wf++;		/* increment pointer ready for writing */
  Vp++;		/* increment pointer for next loop, note increment first
                   as wf derives z from Vp                       */

  wf->a=Vp->a;
  wf->b[0]=psi[1]*exp(-zeta*fabs(Vp->a-r_d)/lambda);
  wf->b[1]=psi[1];


  psi[0]=psi[1];       
  psi[1]=psi[2];

 }

 /* calculate normalisation integral	*/

 wf=wf_start;
 for(i=0;i<n;i++)
 {
  Npsi+=gsl_pow_2(wf->b[0])*delta_z;
  Nchi+=gsl_pow_2(wf->b[1])*delta_z;
  wf++;
 }

 /* divide unnormalised wavefunction by square root
    of normalisation integral                       */
 
 wf=wf_start;
 for(i=0;i<n;i++)
 {
  wf->b[0]=wf->b[0]/sqrt(Npsi);
  wf->b[1]=wf->b[1]/sqrt(Nchi);
  wf++;
 }

 /* generate output filename (and open file for writing) using 
    the basis wf%i.r where the integer %i is the donor index i_d  */

 sprintf(filename,"wf%i.r",i_d);
 fw=fopen(filename,"w");

 /* write wavefunction wf(n).r	*/

 wf=wf_start;

 for(i=0;i<n;i++)
 {
  fprintf(fw,"%20.17le %le %le\n",wf->a,wf->b[0],wf->b[1]);
  wf++;
 }
 fclose(fw);
}

static double I_1(const double lambda,
                  const double z_dash,
                  const double zeta)
{
 /* Eq. 5.100, QWWAD3 */
 return 2*pi*(zeta*fabs(z_dash)*lambda/2+gsl_pow_2(lambda)/4)*
        exp(-2*zeta*fabs(z_dash)/lambda);
}

static double I_2(const double lambda,
                  const double z_dash,
                  const double zeta)
{
  /* Eq. 5.106, QWWAD3 */
  return 2*pi*(-gsl_pow_2(zeta)*z_dash/2)*exp(-2*zeta*fabs(z_dash)/lambda);
}

static double I_3(const double lambda,
                  const double z_dash,
                  const double zeta,
                  const int    N_w)
{
double I_33=0.0;
double I_34=0.0;
double w;

/* Eq. 5.113, QWWAD3 */
const double I_31=((-1-gsl_pow_2(zeta))/2)*exp(-2*zeta*fabs(z_dash)/lambda);
const double I_32=(zeta*fabs(z_dash)/(2*lambda)+0.25)*exp(-2*zeta*fabs(z_dash)/lambda);

/* perform integrations over `w' for I_33 and I_34, area simply given by
   sum of (height of centre of strip * strip width) */
const double delta_w=(1.0-0.0)/(float)N_w;
for(w=delta_w/2;w<1;w+=delta_w)
{
 /* Eq. 5.116, QWWAD3 */
 I_33+=exp(-zeta*fabs(z_dash)*(1/w+w)/lambda)*(1-gsl_pow_2(w))/gsl_pow_2(1+gsl_pow_2(w))*delta_w;

 /* Eq. 5.117, QWWAD3 */
 I_34+=exp(-zeta*fabs(z_dash)*(1/w+w)/lambda)*(1-gsl_pow_2(w))/(w*(1+gsl_pow_2(w)))*delta_w;
}

I_33*=2*(gsl_pow_3(zeta)-zeta)*fabs(z_dash)/lambda;
I_34*=(gsl_pow_4(zeta)-gsl_pow_2(zeta))*gsl_pow_2(z_dash/lambda);

return(2*pi*(I_31+I_32+I_33+I_34));
}



static double I_4(const double lambda,
                  const double z_dash,
                  const double zeta,
                  const int    N_w)
{
double I_40=0.0;
double w;
const double delta_w=(1.0-0.0)/(float)N_w;

/* Eq. 5.118, QWWAD3 */
for (w=delta_w/2;w<1;w+=delta_w)
{
 I_40+=exp(-2*fabs(z_dash)*sqrt(gsl_pow_2((1-gsl_pow_2(w))/(2*w))+gsl_pow_2(zeta))/lambda)
       *fabs(z_dash)*(1-gsl_pow_2(w))/(2*gsl_pow_2(w))*delta_w;
}

return 2*pi*I_40;
}
