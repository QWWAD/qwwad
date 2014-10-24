/*==================================================================
                                d02D
  ==================================================================*/

/* This program implements a variational technique to calculate the
   uncorrelated one particle energies of an electron attatched to a 
   single donor at any position, in any user supplied potential.  
   The potential is read from the file v.r

   This version is a single variational parameter calculation---the 
   two-dimensional (2D) approximation trial wavefunction:

		Psi=chi(z) exp(-r"/lambda)

   where	r"=sqrt(x^2+y^2)

		Input files:
		r_d.r		donor (or acceptor positions)
		v.r		one-dimensional potential

		Output files:
		e.r		total energies for each r_d
		l.r		Bohr radii (lambda) for each r_d
		wfn.r		wave functions, both Psi and chi, n=0,1,2..

   Paul Harrison, February 1993 

   Substantially revised,

   Paul Harrison, February 1998					*/

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <gsl/gsl_math.h>
#include "qclsim-constants.h"
#include "qclsim-fileio.h"
#include "struct.h"

using namespace Leeds;
using namespace constants;

static double I_1        (const double  lambda);
static double I_2        (const double  lambda);
static double I_3        (const double  lambda);
static double I_4        (const double  lambda,
                          const double  z_dash,
                          const size_t  N_w);
static double psi_at_inf(const double                 E,
                         const std::valarray<double> &z,
                         const double                 epsilon,
                         const double                 lambda,
                         const double                 mstar,
                         const double                 r_d,
                         const std::valarray<double> &Vp,
                         const size_t                 N_w);
static bool repeat_lambda(double       *lambda,
                          double       *lambda_0,
                          const double  x,
                          double       *x_min);
static void wavefunctions(const std::valarray<double> &z,
                          const double                 E,
                          const double                 epsilon,
                          const double                 lambda,
                          const double                 mstar,
                          const double                 r_d,
                          const int                    i_d,
                          const std::valarray<double> &Vp,
                          const size_t                 N_w);

int main(int argc,char *argv[])
{
double d_E;                 /* infinitesmal energy               */
double delta_E;             /* small but finite energy           */
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
double y;                   /* function (psi at infinity)        */
double y1;                  /* temporary y value                 */
double y2;                  /* temporary y value                 */
int    i_d;                 /* donor (or acceptor) index         */
int    N_w;                 /* number of strips in w integration */
bool   repeat_flag_lambda; /* variational flag=>new lambda      */
FILE   *fe;                 /* file pointer for energies         */
FILE   *fl;                 /* file pointer for lambda_0         */
FILE   *fr_d;               /* file pointer to donor positions   */
 
/* default values */
delta_E=1e-3*e;
epsilon=13.18*eps0;
lambda_start=50.0e-10;
lambda_step=1.0e-10;
lambda_stop=-1.0e-10;
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
  default :
           printf("Usage:  d02D [-d energy step (\033[1m1\033[0mmeV)][-e relative permittivity \033[1m13.18\033[0m]\n");
	   printf("             [-m mass (\033[1m0.067\033[0mm0)]\n");
	   printf("             [-s starting lambda (\033[1m50\033[0mA)][-t lambda increment (\033[1m1\033[0mA)]\n");
	   printf("             [-u final lambda (\033[1m-1\033[0mA)]\n");
	   exit(0);
 }
 argv++;
 argv++;
 argc--;
 argc--;
}

  std::valarray<double> z; // Spatial location [m]
  std::valarray<double> V; // Confining potential [J]
  read_table_xy("v.r", z, V);

  /* Open files for output of data */
  fe=fopen("e.r","w");           /* E versus r_d	*/
  fl=fopen("l.r","w");           /* lambda versus r_d	*/

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

    /* initial energy estimate=minimum potential-binding energy
                               of particle to free ionised dopant */

    x = V.min() - e*e/(4*pi*epsilon*lambda);   

    /* increment energy-search for f(x)=0 */

    y2=psi_at_inf(x,z,epsilon,lambda,mstar,r_d,V,N_w);

    do
    {
     y1=y2;
     x+=delta_E;
     y2=psi_at_inf(x,z,epsilon,lambda,mstar,r_d,V,N_w);
    }while(y1*y2>0);

   /* improve estimate using midpoint rule */

    x-=fabs(y2)/(fabs(y1)+fabs(y2))*delta_E;

   /* implement Newton-Raphson method */

    do
    {
     y=psi_at_inf(x,z,epsilon,lambda,mstar,r_d,V,N_w);
     dy=(psi_at_inf(x+d_E,z,epsilon,lambda,mstar,r_d,V,N_w)-
         psi_at_inf(x-d_E,z,epsilon,lambda,mstar,r_d,V,N_w))/
        (2.0*d_E);
     x-=y/dy;

    }while(fabs(y/dy)>1e-9*e);

    printf("r_d %le lambda %le energy %le meV\n",r_d,lambda,x/(1e-3*e));       

    repeat_flag_lambda=repeat_lambda(&lambda,&lambda_0,x,&x_min);

    lambda+=lambda_step;     /* increments Bohr radius */

   }while((repeat_flag_lambda&&(lambda_stop<0))||(lambda<lambda_stop));

   E=x_min;   /* assign the energy E to the minimum x_min */

   /* Output neutral dopant binding energies (E) and 
      Bohr radii (lambda) in meV and Angstrom respectively */

   fprintf(fe,"%le %le\n",r_d/1e-10,E/(1e-3*e));
   fprintf(fl,"%le %le\n",r_d/1e-10,lambda_0/1e-10);

   wavefunctions(z,E,epsilon,lambda_0,mstar,r_d,i_d,V,N_w);

   i_d++;            /* donor index */

  }/* end while r_d */

  fclose(fr_d);
  fclose(fe);
  fclose(fl);

  return EXIT_SUCCESS;
} /* end main */

/**
 * \brief compares minimum value of energy for this lambda
 *
 * \details If it is a new true minimum, then repeat for new lambda
 */
static bool repeat_lambda(double       *lambda,
                          double       *lambda_0,
                          const double  x,
                          double       *x_min)
{
 bool flag;

 if(x<*x_min)    
 {      
  *x_min=x;                           /* set new minimum        */
  *lambda_0=*lambda;
  flag=true;                          /* repeat for new lambda  */
 }
 else 
 {
  flag=false;
 }    
 return(flag);
}

/**
 * \brief Finds the value of the wavefunction at +infinity for a given energy.
 *
 * \details The solution to the energy occurs for psi(+infinity)=0.
 *
 * \returns The wavefunction at \f$\psi(\infty)\f$
 */
static double psi_at_inf(const double                 E,
                         const std::valarray<double> &z,
                         const double                 epsilon,
                         const double                 lambda,
                         const double                 mstar,
                         const double                 r_d,
                         const std::valarray<double> &Vp,
                         const size_t                 N_w)
{
    const size_t nz = z.size();
    const double dz = z[1] - z[0];

    std::valarray<double> psi(3); // Wavefunction amplitude at 3 adjacent points

    /* boundary conditions */
    const double kappa=sqrt(2*mstar*(Vp[1]-E))/hBar;

    const double delta_psi = 1.e-10;
    psi[0] = delta_psi; // Initial wave function value (arbitrarily small)
    psi[1] = psi[0]*exp(kappa*dz); // exponential growth

    // Compute coefficients for Schroedinger equation [QWWAD3, 5.23]
    // Note that these are all constants so they can be computed outside
    // the loop
    const double I1=I_1(lambda);
    const double I3=I_3(lambda);
    const double alpha = I1;                // Coeff. of 2nd derivative
    const double beta  = 2.0 * I_2(lambda); // Coeff. of 1st deriviative

    // Loop over spatial points (ignore first)
    for(unsigned int iz = 1; iz < nz; ++iz)
    {
        const double I4 = I_4(lambda, z[iz] - r_d, N_w);

        const double gamma = I3 + (2*mstar*gsl_pow_2(e/hBar)/(4*pi*epsilon))*I4
                             - 2*mstar/hBar * (Vp[iz]-E) * I1/hBar;

        psi[2]=((-1+beta*dz/(2*alpha))*psi[0]
                +(2-dz*dz*gamma/alpha)*psi[1]
               )/(1+beta*dz/(2*alpha));

        psi[0]=psi[1];
        psi[1]=psi[2];
    }

    return psi[0]-delta_psi;
}

/**
 * \brief Calculates and writes the wavefunctions
 *
 * \details Both psi(z) and chi(z) are written to the external file wf(n).r
 */
static void wavefunctions(const std::valarray<double> &z,
                          const double                 E,
                          const double                 epsilon,
                          const double                 lambda,
                          const double                 mstar,
                          const double                 r_d,
                          const int                    i_d,
                          const std::valarray<double> &Vp,
                          const size_t                 N_w)
{
 double delta_psi;       /* initial wavefunction value                  */
 double kappa;
 double Npsi=0;          /* normalisation integral for psi              */
 double Nchi=0;          /* normalisation integral for chi              */
 double psi[3];          /* wavefunctions at z-d_z,z,z+d_z              */
 char   filename[9];     /* character string for wavefunction filename  */
 FILE   *fw;             /* file wf.r                                   */
 data12  *wf_start;      /* pointer to start of w.f                     */
 data12  *wf;		 /* wavefunction pointer, note b[0]=psi, 
                            b[1]=chi, see notes                         */

 const size_t nz = z.size();
 const size_t dz = z[1] - z[0];

 wf_start = new data12[nz]; // allocates memory for wavefunctions
 wf=wf_start;

 /* boundary conditions */
 kappa=sqrt(2*mstar/hBar*(Vp[0]-E)/hBar);

 delta_psi=1.0e-10;

 psi[0]=delta_psi;                 /* arbitrary number close to zero */
 psi[1]=psi[0]*exp(kappa*dz); /* exponential growth produce by  */

 /* write first value of wavefunction */ 
 wf->a=z[0];
 wf->b[0]=psi[0]; 
 wf->b[1]=psi[0]; 

 wf++;
 wf->a=z[1];
 wf->b[0]=psi[1]; 
 wf->b[1]=psi[1]; 
 
 const double I1=I_1(lambda);
 const double I2=I_2(lambda);
 const double I3=I_3(lambda);
 const double alpha = I1;   // Coefficient of second derivative, see notes
 const double beta  = 2*I2; // Coefficient of first derivative

 // calculate unnormalised wavefunction
 // Note that points 0 and 1 were already defined before the loop
 for(unsigned int iz = 2; iz<nz; ++iz)
 {
  const double I4=I_4(lambda, z[iz-1]-r_d, N_w);

  // Coefficient of function
  const double gamma=I3+(2*mstar*gsl_pow_2(e/hBar)/(4*pi*epsilon))*I4
          -(2*mstar/hBar)*(Vp[iz-1]-E)*I1/hBar;

  psi[2]=((-1+beta*dz/(2*alpha))*psi[0]
          +(2-dz*dz*gamma/alpha)*psi[1]
         )/(1+beta*dz/(2*alpha));

  // write wavefunction point corresponding to current z and V

  wf++;	// increment pointer ready for writing

  wf->a=z[iz];
  wf->b[0]=psi[1];
  wf->b[1]=psi[1];

  psi[0]=psi[1];       
  psi[1]=psi[2];
 }

 // calculate normalisation integral
 wf=wf_start;
 for(unsigned int iz=0; iz<nz; iz++)
 {
  Npsi+=gsl_pow_2(wf->b[0])*dz;
  Nchi+=gsl_pow_2(wf->b[1])*dz;
  wf++;
 }

 /* divide unnormalised wavefunction by square root
    of normalisation integral                       */
 
 wf=wf_start;
 for(unsigned int iz=0; iz<nz; iz++)
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

 for(unsigned int iz=0;iz<nz;iz++)
 {
  fprintf(fw,"%20.17le %le %le\n",wf->a,wf->b[0],wf->b[1]);
  wf++;
 }
 fclose(fw);
}

/**
 * Computes the binding energy integral \f$I_1\f$ for a 2D trial wavefunction
 *
 * \param[in] lambda variational parameter [m]
 *
 * \returns Binding energy integral [m^2]
 *
 * \details See Eq. 5.39, QWWAD3. The integral is solved analytically as
 *          \f[
 *            I_1 = 2\pi\frac{\lambda^2}{4}
 *          \f]
 */
static double I_1(const double lambda)
{
 return 2*pi*gsl_pow_2(lambda)/4;
}

/**
 * Computes the binding energy integral \f$I_2\f$ for a 2D trial wavefunction
 *
 * \param[in] lambda variational parameter [m]
 *
 * \returns Binding energy integral [m]
 *
 * \details See Eq. 5.42, QWWAD3. The integral evaluates to zero
 */
static double I_2(const double lambda)
{
 (void)lambda; /* Silence compiler warning about unused param */
 return 0;
}

/**
 * Computes the binding energy integral \f$I_3\f$ for a 2D trial wavefunction
 *
 * \param[in] lambda variational parameter [m]
 *
 * \returns Binding energy integral [dimensionless]
 *
 * \details See Eq. 5.52, QWWAD3. The integral is solved analytically as
 *          \f[
 *            I_3 = 2\pi \left(-\frac{1}{4}\right)
 *          \f]
 */
static double I_3(const double lambda)
{
 (void)lambda; /* Silence compiler warning about unused param */
 return 2*pi*(-0.25);
}

/**
 * Computes the binding energy integral \f$I_4\f$ for a 2D trial wavefunction
 *
 * \param[in] lambda variational parameter [m]
 * \param[in] z_dash displacement between electron and donor in z-direction [m]
 * \param[in] N_w    number of samples to take in integration
 *
 * \returns Binding energy integral [m]
 *
 * \details See Eq. 5.54, QWWAD3. The integral is given by
 *          \f[
 *            I_4 = 2\pi \int_{|z'|}^{\infty}\frac{\exp{-\frac{2\sqrt{r'^2-z'^2}}{\lambda}}}{r'}r'\mathrm{d}r'
 *          \f]
 */
static double I_4(const double lambda,
                  const double z_dash,
                  const size_t N_w)
{
    double I_40=0.0; /* Integral term */
    double w;
    double delta_w=(1.0-0.0)/(float)N_w;

    for (w=delta_w/2;w<1;w+=delta_w)
    {
        I_40+=exp(-fabs(z_dash)*(1/w-w)/lambda)*fabs(z_dash)
            *(1-gsl_pow_2(w))/(2*gsl_pow_2(w))*delta_w;
    }

    return 2*pi*I_40;
}
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :
