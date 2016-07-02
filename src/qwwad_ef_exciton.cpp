/**
 * \file  qwwad_ef_exciton.cpp
 * \brief Find exciton binding energies
 * \author Paul Harrison  <p.harrison@shu.ac.uk>
 * \author Alex Valavanis <a.valavanis@leeds.ac.uk>
 */

//#include <cstdio>
#include <cstdlib>
#include <valarray>

#include <gsl/gsl_math.h>

#include "qwwad/constants.h"
#include "qwwad/file-io.h"
#include "qwwad/options.h"

using namespace QWWAD;
using namespace constants;

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

static double Eb_1S(const std::valarray<double> z,
                    const std::valarray<double> psi_e,
                    const std::valarray<double> psi_h,
                    const std::valarray<double> p,
                    const double  beta,
                    const double  epsilon,
                    const double  lambda,
                    const double  m[],
                    const double  mu_xy,
                    const int     N_x,
                    const bool    output_flag);

std::valarray<double> pP_calc(const std::valarray<double> &z,
                              const std::valarray<double> &psi_e,
                              const std::valarray<double> &psi_h);

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

    opt.add_option<double>      ("dcpermittivity,e",     13.18,  "Bulk relative permittivity");
    opt.add_option<double>      ("betastart,w",           0.001, "Initial value for beta symmetry parameter search");
    opt.add_option<double>      ("betastep,x",            0.05,  "Increment for beta symmetry parameter search");
    opt.add_option<double>      ("betastop,y",           -1.0,   "Final value for beta symmetry parameter search");
    opt.add_option<double>      ("lambdastart,s",        70.0,   "Initial value for Bohr radius search [Angstrom]");
    opt.add_option<double>      ("lambdastep,t",          1.0,   "Increment for beta symmetry parameter search");
    opt.add_option<double>      ("lambdastop,u",         -1.0,   "Final value for beta symmetry parameter search");
    opt.add_option<double>      ("electronmass,m",        0.067, "Bulk electron effective mass (relative to free electron)");
    opt.add_option<double>      ("holemass,n",            0.62,  "Bulk hole effective mass (relative to free electron)");
    opt.add_option<unsigned int>("electronstate,a",       1,     "Index of electron state");
    opt.add_option<unsigned int>("holestate,b",           1,     "Index of hole state");
    opt.add_option<size_t>      ("nx,N",                100,     "Number of samples for x-integration");
    opt.add_option<std::string> ("searchlogfile", "searchlog.r", "Filename for search log");

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

    const auto e_state = opt.get_option<unsigned int>("electronstate");
    const auto h_state = opt.get_option<unsigned int>("holestate");

    // e and h z-axis masses
    double m[2];
    m[0] = opt.get_option<double>("electronmass") * me;
    m[1] = opt.get_option<double>("holemass") * me;

    const auto output_flag = opt.get_verbose(); // if set, write data to screen

    // Construct filenames for wavefunction input and read data
    // TODO: We've assumed for now that the spatial positions are identical in each file
    std::ostringstream psi_e_file, psi_h_file;
    psi_e_file << "wf_e" << e_state << ".r";
    psi_h_file << "wf_h" << h_state << ".r";

    std::valarray<double> z;     // Spatial location [m]
    std::valarray<double> psi_e; // Electron wave function [m^{-0.5}]
    std::valarray<double> psi_h; // Hole wave function [m^{-0.5}]
    read_table(psi_e_file.str(), z, psi_e);
    read_table(psi_h_file.str(), z, psi_h);

    double Eb_min = e;       // minimum Eb for lambda variation, i.e., 1 eV !

    double m_xy[2]; // e and h x-y plane masses
    m_xy[0]=m[0];	m_xy[1]=m[1];	// assumes isotropic mass for now

    const double mu_xy=1/(1/m_xy[0]+1/m_xy[1]);	/* calculate reduced mass in-plane	*/

//    FILE *Fbeta = fopen("beta-lambda.r","w");
//    FILE *FEX0l=fopen("EX0-lambda.r","w");

    // calculate probabilities
    const auto p = pP_calc(z, psi_e, psi_h);

    // Output table of probabilities
    const std::valarray<double> p_Angstrom = p*1e-10;
    write_table("p.r", z, p_Angstrom);

    if(output_flag)printf("  l/A   beta   Eb/meV  T/meV  V/meV   OS/arb.\n");

    double lambda=lambda_start; // Bohr radius
    double beta_0 = beta_start; // beta for Eb_min

    // TODO: lambda_0 is found iteratively. Check that this is a sensible initial value
    double lambda_0 = 0;            // lambda for Eb_min
    bool repeat_flag_lambda=true; // repeat variational lambda loop

    // Tables for the search log
    std::vector<double> lambda_log; 
    std::vector<double> beta_log;
    std::vector<double> Eb_log;

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
            // Find exciton binding energy (<0=bound)
            const auto Eb = Eb_1S(z, psi_e, psi_h, p, beta, epsilon, lambda, m, mu_xy, N_x, output_flag);

            repeat_flag_beta=repeat_beta(beta,&beta_0_lambda,Eb,&Eb_min_beta);

            // Update search log
            lambda_log.push_back(lambda*1e10);
            beta_log.push_back(beta);
            Eb_log.push_back(Eb_min_beta*1000/e);

            beta+=beta_step;
        }while((repeat_flag_beta&&(beta_stop<0))||(beta<beta_stop));

//        fprintf(FEX0l,"%lf %lf\n",lambda/1e-10,Eb_min_beta/(1e-3*e));
//        fprintf(Fbeta,"%lf %lf\n",lambda/1e-10,beta_0_lambda);

        repeat_flag_lambda=repeat_lambda(&beta_0,beta_0_lambda,Eb_min_beta,&Eb_min,
                lambda,&lambda_0);

        lambda+=lambda_step;   /* increment Bohr radius */
    }while((repeat_flag_lambda&&(lambda_stop<0))||(lambda<lambda_stop));

    // Write out final data to file
    std::ofstream FEX0("EX0.r");
    FEX0 << Eb_min*1000/e << " " << lambda_0*1e10 << " " << beta_0 << std::endl;
    FEX0.close();
    
    const auto searchlogfile = opt.get_option<std::string>("searchlogfile");
    write_table(searchlogfile, lambda_log, beta_log, Eb_log);

//    fclose(Fbeta);
//    fclose(FEX0l);

    return EXIT_SUCCESS;
}

/**
 * \brief Find binding energy of 1S exciton
 *
 * \param probs pointer to pP structure
 */
static double Eb_1S(const std::valarray<double> z,
                    const std::valarray<double> psi_e,
                    const std::valarray<double> psi_h,
                    const std::valarray<double> p,
                    const double  beta,
                    const double  epsilon,
                    const double  lambda,
                    const double  m[],
                    const double  mu_xy,
                    const int     N_x,
                    const bool    output_flag)
{
    double A  = 0;              /* {\cal A}, see notes!              */
    double B  = 0;              /* {\cal B}, see notes!              */
    double Ct = 0;              /* kinetic energy component of C     */
    double Cv = 0;              /* potential energy component of C   */
    double D  = 0;              /* {\cal D}, see notes!              */
    double Eb;                  /* exciton binding energy (<0=bound) */
    double O  = 0;              /* {\cal O}, overlap integral        */
    double C  = 0;              /* electron--hole interaction term */

    const auto nz = z.size();
    const auto dz = z[1] - z[0];

    // Loop over spatial separations
    for(unsigned int ia=0; ia<nz; ++ia)
    {
        A  += p[ia] * G(z[ia], beta,lambda,N_x)*dz;
        B  += p[ia] * G(z[ia], beta,lambda,N_x)*dz;
        Ct += p[ia] * J(z[ia], beta,lambda,N_x)*dz;
        Cv += p[ia] * K(z[ia], beta,lambda,N_x)*dz;
        D  += p[ia] * F(z[ia], beta,lambda)*dz;

        O += psi_e[ia]*psi_h[ia]*dz;
    }

    A*=hBar*hBar/(2*m[0]);         /* multiply integrals by constant factors */
    B*=hBar*hBar/(2*m[1]);
    Ct*=-hBar*hBar/(2*mu_xy);
    Cv*=-e*e/(4*pi*epsilon);

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
 * \brief Calculate probabilities known as p(a), Pm(a) and Pmu(a)
 */
std::valarray<double> pP_calc(const std::valarray<double> &z,
                              const std::valarray<double> &psi_e,
                              const std::valarray<double> &psi_h)
{
    const auto nz = z.size();    // Number of spatial samples  
    const auto dz = z[1] - z[0]; // Separation between spatial samples [m]

    std::valarray<double> p(nz);

    // Loop over spatial separation
    for(unsigned int ia=0; ia<nz; ++ia)
    {
        p[ia] = 0.0;    // initialize probability to zero

        // Loop over hole location and find p(a) [QWWAD4, 6.23]
        for(unsigned int iz = 0; iz < nz-ia; ++iz)
        {
            p[ia] += (psi_e[iz+ia]*psi_e[iz+ia]*psi_h[iz]*psi_h[iz]
                     +psi_h[iz+ia]*psi_h[iz+ia]*psi_e[iz]*psi_e[iz]) * dz;
        }
    }

    return p;
}
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :
