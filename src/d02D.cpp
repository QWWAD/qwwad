/**
 * \file   d02D.cpp
 * \brief  Calculates state of electron attached to donor in a user-specified potential
 * \author Paul Harrison  <p.harrison@shu.ac.uk>
 * \author Alex Valavanis <a.valavanis@leeds.ac.uk>
 *
 * \details Implements a variational technique to calculate the
 *          uncorrelated one particle energies of an electron attatched to a 
 *          single donor at any position, in any user supplied potential.  
 *          The potential is read from the file v.r
 *
 *          This version is a single variational parameter calculation---the 
 *          two-dimensional (2D) approximation trial wavefunction:
 *
 *          \f[
 *            \Psi(z)=\chi(z) \exp(-r^{\prime\prime}/\lambda)
 *          \f]
 *
 *          where \f$r^{\prime\prime} = \sqrt(x^2+y^2)\f$
 *
 *   Input files:
 *  r_d.r		donor (or acceptor positions)
 *  v.r		one-dimensional potential
 *
 *  Output files:
 *  e.r		total energies for each r_d
 *  l.r		Bohr radii (lambda) for each r_d
 *  wfn.r		wave functions, both Psi and chi, n=0,1,2..
 */

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <gsl/gsl_math.h>
#include "qclsim-constants.h"
#include "qclsim-fileio.h"
#include "qwwad-options.h"

using namespace Leeds;
using namespace constants;

static double I_1        (const double                 lambda);
static double I_2        (const double                 lambda);
static double I_3        (const double                 lambda);
static double I_4        (const double                 lambda,
                          const double                 z_dash,
                          const size_t                 N_w);
static double psi_at_inf (const double                  E,
                          const std::valarray<double>  &z,
                          const double                  epsilon,
                          const double                  lambda,
                          const double                  mstar,
                          const double                  r_d,
                          const std::valarray<double>  &Vp,
                          const size_t                  N_w);
static bool repeat_lambda(const double                  lambda,
                          double                       &lambda_0,
                          const double                  E,
                          double                       &E0);
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
    Options opt;
    std::string doc("Find state of electron attached to a donor in a 2D system");

    opt.add_numeric_option("dE,d",              1, "Energy step for Shooting solver [meV]");
    opt.add_numeric_option("epsilon,e",     13.18, "Bulk relative permittivity");
    opt.add_numeric_option("mass,m",        0.067, "Bulk effective mass (relative to free electron)");
    opt.add_numeric_option("lambdastart,s",    50, "Initial value for Bohr radius search [Angstrom]");
    opt.add_numeric_option("lambdastep,t",      1, "Step size for Bohr radius search [Angstrom]");
    opt.add_numeric_option("lambdastop,u",     -1, "Final value for Bohr radius search [Angstrom]");

    opt.add_prog_specific_options_and_parse(argc, argv, doc);

    /* computational default values */
    int N_w=100;             // Number of strips in w integration

    const double delta_E = opt.get_numeric_option("dE") * 1e-3*e;    // Energy increment [J]
    const double epsilon = opt.get_numeric_option("epsilon") * eps0; // Permittivity [F/m]
    const double mstar   = opt.get_numeric_option("mass") * me;      // Effective mass [kg]

    const double lambda_start = opt.get_numeric_option("lambdastart") * 1e-10; // Initial Bohr radius [m]
    const double lambda_step  = opt.get_numeric_option("lambdastep")  * 1e-10; // Bohr radius increment [m]
    const double lambda_stop  = opt.get_numeric_option("lambdastop")  * 1e-10; // Final Bohr radius [m]

    std::valarray<double> z; // Spatial location [m]
    std::valarray<double> V; // Confining potential [J]
    read_table_xy("v.r", z, V);

    /* Open files for output of data */
    FILE *fe=fopen("e.r","w");           // E versus r_d
    FILE *fl=fopen("l.r","w");           // lambda versus r_d

    // Read list of donor (or acceptor) positions
    std::valarray<double> r_d; // [m]
    read_table_x("r_d.r", r_d);

    // Perform variational calculation for each donor/acceptor position
    for(unsigned int i_d = 0; i_d < r_d.size(); ++i_d)
    {
        double lambda=lambda_start; // Initial Bohr radius value [m]
        double x_min=e;             // Minimum energy of single donor 1eV

        // TODO: lambda_0 is found iteratively. Check that this is a sensible initial value
        double lambda_0 = 0;        /* Bohr radius of electron (or hole) */
        bool   repeat_flag_lambda;  /* variational flag=>new lambda      */

        /* Variational calculation */
        do
        {
            /* initial energy estimate=minimum potential-binding energy
               of particle to free ionised dopant */
            double x = V.min() - e*e/(4*pi*epsilon*lambda);   
            double y1;                  /* temporary y value                 */

            /* increment energy-search for f(x)=0 */
            double y2=psi_at_inf(x,z,epsilon,lambda,mstar,r_d[i_d],V,N_w);

            do
            {
                y1=y2;
                x+=delta_E;
                y2=psi_at_inf(x,z,epsilon,lambda,mstar,r_d[i_d],V,N_w);
            }while(y1*y2>0);

            /* improve estimate using midpoint rule */

            x-=fabs(y2)/(fabs(y1)+fabs(y2))*delta_E;

            /* implement Newton-Raphson method */
            double y;  // function (psi at infinity)
            double d_E=delta_E/1e+6; // Infinitesimal energy
            double dy; // Derivative
            do
            {
                y=psi_at_inf(x,z,epsilon,lambda,mstar,r_d[i_d],V,N_w);
                std::cout << y << std::endl;
                dy=(psi_at_inf(x+d_E,z,epsilon,lambda,mstar,r_d[i_d],V,N_w)-
                    psi_at_inf(x-d_E,z,epsilon,lambda,mstar,r_d[i_d],V,N_w))/
                    (2.0*d_E);
                x-=y/dy;
            }while(fabs(y/dy)>1e-9*e);

            printf("r_d %le lambda %le energy %le meV\n",r_d[i_d],lambda,x/(1e-3*e));

            repeat_flag_lambda=repeat_lambda(lambda,lambda_0,x,x_min);

            lambda+=lambda_step;     /* increments Bohr radius */
        }while((repeat_flag_lambda&&(lambda_stop<0))||(lambda<lambda_stop));

        const double E=x_min; // assign the energy E to the minimum x_min

        /* Output neutral dopant binding energies (E) and 
           Bohr radii (lambda) in meV and Angstrom respectively */

        fprintf(fe,"%le %le\n",r_d[i_d]/1e-10,E/(1e-3*e));
        fprintf(fl,"%le %le\n",r_d[i_d]/1e-10,lambda_0/1e-10);

        wavefunctions(z,E,epsilon,lambda_0,mstar,r_d[i_d],i_d,V,N_w);
    }/* end loop over r_d */

    fclose(fe);
    fclose(fl);

    return EXIT_SUCCESS;
}

/**
 * \brief compares minimum value of energy for this lambda
 *
 * \param[in]     lambda   The new Bohr radius to be tested [m]
 * \param[in,out] lambda_0 The previous estimate of the Bohr radius [m]
 * \param[in]     E        The new binding energy [J]
 * \param[in,out] E0       The previous estimate of the binding energy [J]
 *
 * \returns True is the binding energy is lower than the previous estimate
 *
 * \details The new binding energy is compared with the previous value. If it
 *          is lower, then the Bohr radius and binding energy are overwritten
 *          and the function returns "true" to indicate that we haven't yet
 *          found the absolute minimum.
 */
static bool repeat_lambda(const double  lambda,
                          double       &lambda_0,
                          const double  E,
                          double       &E0)
{
    bool flag = false; // True if we find a new minimum 

    // Set new minimum if required
    if(E < E0)
    {
        E0       = E;
        lambda_0 = lambda;
        flag=true;
    }

    return flag;
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

    // boundary conditions
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
    const size_t nz = z.size();
    const double dz = z[1] - z[0];

    std::valarray<double> psi(nz); // Complete expression for wavefunction [m^{-1/2}]
    std::valarray<double> chi(nz); // Wavefunction envelope [m^{-1/2}]

    // boundary conditions
    psi[0] = 1;
    chi[0] = 1;
    double chi_next = 1.0; 

    const double I1=I_1(lambda);
    const double I2=I_2(lambda);
    const double I3=I_3(lambda);
    const double alpha = I1;   // Coefficient of second derivative, see notes
    const double beta  = 2*I2; // Coefficient of first derivative

    // calculate unnormalised wavefunction
    // Note that points 0 and 1 were already defined before the loop
    for(unsigned int iz = 0; iz<nz; ++iz)
    {
        // Wave function amplitude at previous point
        double chi_prev = 0.0;

        if(iz != 0)
            chi_prev = chi[iz-1];

        const double I4=I_4(lambda, z[iz]-r_d, N_w);

        // Coefficient of function
        const double gamma = I3 + mstar*e*e*I4/(2.0*pi*epsilon*hBar*hBar)
                             - 2.0*mstar*(Vp[iz]-E)*I1/(hBar*hBar);

        chi_next = ((-1.0+beta*dz/(2.0*alpha))*chi_prev
                    +(2.0-dz*dz*gamma/alpha)*chi[iz]
                   )/(1.0+beta*dz/(2.0*alpha));

        if (iz != nz - 1)
        {
            chi[iz+1] = chi_next;

            // The complete wave function at (x,y) = 0
            // is just the same as the envelope when we're considering
            // the 2D symmetrical case
            psi[iz+1] = chi[iz+1];
        }
    }

    // calculate normalisation integral
    double Npsi=pow(psi,2.0).sum() * dz; // normalisation integral for psi
    double Nchi=pow(chi,2.0).sum() * dz; // normalisation integral for chi

    /* divide unnormalised wavefunction by square root
       of normalisation integral                       */
    psi /= sqrt(Npsi);
    chi /= sqrt(Nchi);

    /* generate output filename (and open file for writing) using 
       the basis wf%i.r where the integer %i is the donor index i_d  */
    char   filename[9];     /* character string for wavefunction filename  */
    sprintf(filename,"wf%i.r",i_d);

    write_table_xyz(filename, z, psi, chi);
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
 * \param[in] nw    number of samples to take in integration
 *
 * \returns Binding energy integral [m]
 *
 * \details See Eq. 5.59, QWWAD3. The integral is given by
 *          \f[
 *              I_4=2\pi\int_0^1 \exp{\left[-\frac{\vert z^\prime\vert\left(\frac{1}{w}-w\right) }{\lambda}\right]}\;\;\vert z^\prime\vert \frac{1-w^2}{2w^2}\;\;\text{d}w
 *          \f]
 */
static double I_4(const double lambda,
                  const double z_dash,
                  const size_t nw)
{
    const double dw = 1.0 /(float)(nw-1);   // Step size in w variable
    const double z_dash_abs = fabs(z_dash); // Magnitude of displacement [m]

    double I_40=0.0; // Integral term

    // Exclude w = 0 point to avoid singularity
    for (unsigned int iw = 1; iw < nw; ++iw)
    {
        const double w = iw * dw;

        I_40 += exp(-z_dash_abs * (1/w-w)/lambda) *
                (1-w*w)/(2*w*w)*dw;
    }

    return 2*pi*z_dash_abs * I_40;
}
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :
