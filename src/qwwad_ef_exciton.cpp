/**
 * \file  qwwad_ef_exciton.cpp
 * \brief Find exciton binding energies
 * \author Paul Harrison  <p.harrison@shu.ac.uk>
 * \author Alex Valavanis <a.valavanis@leeds.ac.uk>

   Input files:
     wf_eX.r      electron wavefunction versus z  
     wf_hX.r      hole wavefunction versus z  

   Output files:
     ABC.r         terms A, B, Ct and Cv versus lambda
     beta.r        beta corresponding to minimum Eb versus lambda
     EX0.r         minimum binding energy, corresponding lambda and beta
     EX0-lambda.r  minimum binding energy versus lambda
     p.r           uncorrelated probaility of e-h separation
   */

#include <cstdio>
#include <cstdlib>
#include <gsl/gsl_math.h>

#include "qwwad/constants.h"
#include "qwwad/options.h"

using namespace QWWAD;
using namespace constants;

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

static bool repeat_beta(const double  beta,
                        double       *beta_0_lambda,
                        const double  Eb,
                        double       *Eb_min_beta);

static bool repeat_lambda(double       *beta_0,
                          const double  beta_0_lambda,
                          const double  Eb_min_beta,
                          double       *Eb_min,
                          const double  lambda,
                          double       *lambda_0);

static double Eb_1S(files        *fdata,
                    probs        *ppP,
                    FILE         *FABC,
                    const double  beta,
                    const double  delta_a,
                    const double  epsilon,
                    const double  lambda,
                    const double  m[],
                    const double  mu_xy,
                    const int     N_x,
                    const int     n,
                    const bool    output_flag);

files * read_data(unsigned int state[], int *n);

double read_delta_z(files *Vp);

probs * pP_calc(double  delta_a,
                int     n,
                files  *data_start);

double F(double a,
         double beta,
         double lambda);

double G(double a,
         double beta,
         double lambda,
         int    N_x);

double J(double a,
         double beta,
         double lambda,
         int    N_x);

double K(double a,
         double beta,
         double lambda,
         int N_x);

/**
 * \brief Configure command-line options
 */
Options configure_options(int argc, char* argv[])
{
    Options opt;
    std::string doc("Find exciton binding energy");

    opt.add_option<double>      ("dcpermittivity,e", 13.18,  "Bulk relative permittivity");
    opt.add_option<double>      ("betastart,w",       0.001, "Initial value for beta symmetry parameter search");
    opt.add_option<double>      ("betastep,x",        0.05,  "Increment for beta symmetry parameter search");
    opt.add_option<double>      ("betastop,y",       -1.0,   "Final value for beta symmetry parameter search");
    opt.add_option<double>      ("lambdastart,s",    70.0,   "Initial value for Bohr radius search [Angstrom]");
    opt.add_option<double>      ("lambdastep,t",      1.0,   "Increment for beta symmetry parameter search");
    opt.add_option<double>      ("lambdastop,u",     -1.0,   "Final value for beta symmetry parameter search");
    opt.add_option<double>      ("electronmass,m",    0.067, "Bulk electron effective mass (relative to free electron)");
    opt.add_option<double>      ("holemass,n",        0.62,  "Bulk hole effective mass (relative to free electron)");
    opt.add_option<unsigned int>("electronstate,a",   1,     "Index of electron state");
    opt.add_option<unsigned int>("holestate,b",       1,     "Index of hole state");
    opt.add_option<size_t>      ("nx,N",            100,     "Number of samples for x-integration");

    opt.add_prog_specific_options_and_parse(argc, argv, doc);

    return opt;
}

int main(int argc,char *argv[])
{
    const auto opt = configure_options(argc, argv);

    const auto epsilon      = opt.get_option<double>("dcpermittivity") * eps0; // Permittivity [F/m]
    const auto beta_start   = opt.get_option<double>("betastart");             // Initial beta
    const auto beta_step    = opt.get_option<double>("betastep");              // beta increment
    const auto beta_stop    = opt.get_option<double>("betastop");              // Final beta
    const auto lambda_start = opt.get_option<double>("lambdastart") * 1e-10;   // Initial Bohr radius [m]
    const auto lambda_step  = opt.get_option<double>("lambdastep")  * 1e-10;   // Bohr radius increment [m]
    const auto lambda_stop  = opt.get_option<double>("lambdastop")  * 1e-10;   // Final Bohr radius [m]
    const auto N_x          = opt.get_option<size_t>("nx");                    // Number of steps for x-integration

    // default values
    unsigned int state[2]; // electron and hole states
    state[0] = opt.get_option<unsigned int>("electronstate");
    state[1] = opt.get_option<unsigned int>("holestate");

    // e and h z-axis masses
    double m[2];
    m[0] = opt.get_option<double>("electronmass") * me;
    m[1] = opt.get_option<double>("holemass") * me;

    const auto output_flag = opt.get_verbose(); // if set, write data to screen

    int n; // length of potential file
    files *data_start = read_data(state,&n); // reads wave functions

    // z separation of input potentials 
    double delta_z = read_delta_z(data_start);
    double delta_a = delta_z; // separation of adjacent a values
    double Eb_min  = e;       // minimum Eb for lambda variation, i.e., 1 eV !

    double m_xy[2]; // e and h x-y plane masses
    m_xy[0]=m[0];	m_xy[1]=m[1];	// assumes isotropic mass for now

    const double mu_xy=1/(1/m_xy[0]+1/m_xy[1]);	/* calculate reduced mass in-plane	*/

    FILE *FABC  = fopen("ABC.r","w");
    FILE *Fbeta = fopen("beta-lambda.r","w");
    FILE *FEX0l=fopen("EX0-lambda.r","w");

    // calculates p and P's, returns start address of structure
    probs *pP_start=pP_calc(delta_z,n,data_start);

    if(output_flag)printf("  l/A   beta   Eb/meV  T/meV  V/meV   OS/arb.\n");

    double lambda=lambda_start; // Bohr radius
    double beta_0 = beta_start; // beta for Eb_min

    // TODO: lambda_0 is found iteratively. Check that this is a sensible initial value
    double lambda_0 = 0;            // lambda for Eb_min
    bool repeat_flag_lambda=true; // repeat variational lambda loop

    do
    {
        /* Find minimum binding energy for beta variation */
        double Eb_min_beta = e; /* Start with 1 eV as a huge initial estimate */
        bool   repeat_flag_beta = true;   /* repeat variational beta loop flag */
        double beta=beta_start;
        double beta_0_lambda = beta; /* Value of beta that gives the minimum binding energy */

        /* Loop through beta values and store the minimum binding energy as
         * Eb_min_beta.  The corresponding beta value is stored as beta_0_lambda */ 
        do
        {
            /* Find exciton binding energy (<0=bound) */
            const double Eb=Eb_1S(data_start,pP_start,FABC,beta,delta_a,epsilon,lambda,m,mu_xy,N_x,n,output_flag);

            repeat_flag_beta=repeat_beta(beta,&beta_0_lambda,Eb,&Eb_min_beta);

            beta+=beta_step;
        }while((repeat_flag_beta&&(beta_stop<0))||(beta<beta_stop));

        fprintf(FEX0l,"%lf %lf\n",lambda/1e-10,Eb_min_beta/(1e-3*e));
        fprintf(Fbeta,"%lf %lf\n",lambda/1e-10,beta_0_lambda);

        repeat_flag_lambda=repeat_lambda(&beta_0,beta_0_lambda,Eb_min_beta,&Eb_min,
                lambda,&lambda_0);

        lambda+=lambda_step;   /* increment Bohr radius */
    }while((repeat_flag_lambda&&(lambda_stop<0))||(lambda<lambda_stop));

    /* Write out final data to file	*/
    FILE *FEX0=fopen("EX0.r","w");
    fprintf(FEX0,"%6.3lf %6.2lf %6.3lf\n",Eb_min/(1e-3*e),lambda_0/1e-10,beta_0);
    fclose(FEX0);

    fclose(FABC);
    fclose(Fbeta);
    fclose(FEX0l);

    free(data_start);
    free(pP_start);

    return EXIT_SUCCESS;
}

/**
 * \brief Find binding energy of 1S exciton
 *
 * \param fdata pointer to data structure
 * \param probs pointer to pP structure
 * \param FABC  file pointer to ABC.r
 */
static double Eb_1S(files        *fdata,
                    probs        *ppP,
                    FILE         *FABC,
                    const double  beta,
                    const double  delta_a,
                    const double  epsilon,
                    const double  lambda,
                    const double  m[],
                    const double  mu_xy,
                    const int     N_x,
                    const int     n,
                    const bool    output_flag)
{
 double A  = 0;              /* {\cal A}, see notes!              */
 double B  = 0;              /* {\cal B}, see notes!              */
 double Ct = 0;              /* kinetic energy component of C     */
 double Cv = 0;              /* potential energy component of C   */
 double D  = 0;              /* {\cal D}, see notes!              */
 double Eb;                  /* exciton binding energy (<0=bound) */
 double O  = 0;              /* {\cal O}, overlap integral        */
 int    i_a;                 /* index through a values            */
 double C  = 0;              /* electron--hole interaction term */
 
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

 A*=hBar*hBar/(2*m[0]);         /* multiply integrals by constant factors */
 B*=hBar*hBar/(2*m[1]);
 Ct*=-hBar*hBar/(2*mu_xy);
 Cv*=-e*e/(4*pi*epsilon);

 fprintf(FABC,"%6.2lf %6.3lf %6.3lf %6.3lf %6.3lf\n",lambda/1e-10,
         A/D/(1e-3*e),B/D/(1e-3*e),Ct/D/(1e-3*e),Cv/D/(1e-3*e));

 C = Ct + Cv;  /* Find e--h term [QWWAD4, eq. 6.36] */
 Eb=(A+B+C)/D; /* Binding energy [QWWAD4, eq. 6.44] */

 if(output_flag)
  printf("%6.2lf %6.3lf %6.3lf %6.3lf %6.3lf %6.3le\n",lambda/1e-10,beta,
         Eb/(1e-3*e),(A+B+Ct)/D/(1e-3*e),Cv/D/(1e-3*e),gsl_pow_2(O)/D);

 return Eb;
}

/**
 * \brief Repeat beta variational
 */
static bool repeat_beta(const double  beta,
                        double       *beta_0_lambda,
                        const double  Eb,
                        double       *Eb_min_beta)
{
    bool flag = false;

    if(Eb < *Eb_min_beta)
    {
        *Eb_min_beta   = Eb;
        *beta_0_lambda = beta;
        flag           = true;
    }

    return flag;
}

/**
 * \brief Repeat lambda variational
 */
static bool repeat_lambda(double       *beta_0,
                          const double  beta_0_lambda,
                          const double  Eb_min_beta,
                          double       *Eb_min,
                          const double  lambda,
                          double       *lambda_0)
{
    bool flag = false;

    /* If this is the smallest binding energy
     * so far, copy the energy, lambda and beta
     * and flag the function for another iteration */
    if(Eb_min_beta < *Eb_min)
    {
        *Eb_min   = Eb_min_beta;
        *lambda_0 = lambda;
        *beta_0   = beta_0_lambda;
        flag      = true;
    }

    return flag;
}

/**
 * \brief returns the value of F(a)
 */
double F(double a,
         double beta,
         double lambda)
{
 double f;
 f=2*pi*lambda*(sqrt(1-gsl_pow_2(beta))*a/2+lambda/4)*
   exp(-2*sqrt(1-gsl_pow_2(beta))*a/lambda);

 return(f);
}

/* This function returns the value of G(a),
   to overcome the problem of divergence when
   x=0, the integration is performed using a
   midpoint sum, the strip width being delta_x */
double G(double a,
         double beta,
         double lambda,
         int    N_x)
{
 double delta_x;
 double g;
 double x;  /* dummy variable---see notes! */

 delta_x=(1.0-0.0)/(float)N_x;

 g=0;       /* initialize variable */

 for(x=delta_x/2;x<1;x+=delta_x)
 {
  g+=1/(1/x+x)*exp(-sqrt(1-gsl_pow_2(beta))*a*(1/x+x)/lambda)
     *(1/gsl_pow_2(x)-1)*delta_x;
 }
 g*=2*pi*gsl_pow_2(1-gsl_pow_2(beta))*gsl_pow_2(a)/gsl_pow_2(lambda);

 return(g);
}

/**
 * \brief returns the value of J(a)
 *
 * \details To overcome the problem of divergence when x=0, the integration is performed
 *          using a midpoint sum, the strip width being delta_x
 */
double J(double a,
         double beta,
         double lambda,
         int    N_x)
{
 double delta_x;
 double j13;     /* J1+J3---see notes! */
 double j24;     /* J2+J4---see notes! */
 double x;       /* dummy variable---see notes! */

 delta_x=(1.0-0.0)/(float)N_x;

 j13=2*pi*(sqrt(1-gsl_pow_2(beta))*a/(2*lambda)-0.25)
     *exp(-2*sqrt(1-gsl_pow_2(beta))*a/lambda);

 j24=0;        /* initialize variable */

 for(x=delta_x/2;x<1;x+=delta_x)
 {
  j24+=(
        -1/(lambda*gsl_pow_2(1/x+x)/4)
        -sqrt(1-gsl_pow_2(beta))*a/(gsl_pow_2(lambda)*(1/x+x)/2)
       )
       *exp(-sqrt(1-gsl_pow_2(beta))*a*(1/x+x)/lambda)
       *(1/gsl_pow_2(x)-1)*delta_x;
 }
 j24*=2*pi*sqrt(1-gsl_pow_2(beta))*a/2;

 return(j13+j24);
}

/**
 * \brief returns the value of K(a), to overcome the problem of divergence when
 *        x=0, the integration is performed using a midpoint sum, the strip width being delta_x
 */
double K(double a,
         double beta,
         double lambda,
         int N_x)
{
 double delta_x;      /* step length of integration */
 double k;    
 double lower_limit;  /* lower limit of integration */
 double upper_limit;  /* upper limit of integration */
 double x;            /* dummy variable---see notes! */

 upper_limit=(1-sqrt(1-gsl_pow_2(beta)))/beta;
 lower_limit=0;

 delta_x=(upper_limit-lower_limit)/(float)N_x;

 k=0;

 for(x=lower_limit+delta_x/2;x<upper_limit;x+=delta_x)
 {
  k+=exp(-beta*a*(1/x-x)/lambda)
     *(1/gsl_pow_2(x)-1)*delta_x;
 }
 k*=2*pi*beta*a/2;

 return(k);
}

/**
 * \brief Calculates the separation along the z (growth)
 *        direction of the user supplied wave functions, assumes regular
 *        one-dimensional mesh
 */
double read_delta_z(files *Vp)
{
 double z[2];           /* displacement along growth direction     */

 z[0] = Vp->z;
 Vp++;
 z[1] = Vp->z;
 return(z[1]-z[0]);
}

/**
 * \brief reads the potential into memory and returns the start
 *        address of this block of memory and the number of lines
 */
files * read_data(unsigned int state[], int *n)
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
 while(fscanf(Fwfe,"%*e %*e")!=EOF)
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
  int n_read = fscanf(Fwfh,"%*e %le",&(ft->wf[1]));

  if (n_read == 2)
    ft++;
 }

 fclose(Fwfe);
 fclose(Fwfh);

 return(data_start);

}

/**
 * \brief Calculate probabilities known as p(a), Pm(a) and Pmu(a)
 *
 * \param[in] delta_a    distance between adjacent points
 * \param[in] n          number of lines in wavefunction file
 * \param[in] data_start pointer to beginning of wavefunctions
 */
probs * pP_calc(double  delta_a,
                int     n,
                files  *data_start)
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
   ppP->p+=(gsl_pow_2((fdata+i)->wf[0])*gsl_pow_2(fdata->wf[1])
            +gsl_pow_2(fdata->wf[0])*gsl_pow_2((fdata+i)->wf[1])
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
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :
