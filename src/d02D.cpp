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
#include <gsl/gsl_errno.h>
#include <gsl/gsl_roots.h>
#include "qclsim-constants.h"
#include "qclsim-fileio.h"
#include "qclsim-maths.h"
#include "qwwad-options.h"

using namespace Leeds;
using namespace constants;

static double I_1        (const double                 lambda);
static double I_2        (const double                 lambda);
static double I_3        (const double                 lambda);
static double I_4        (const double                 lambda,
                          const double                 z_dash,
                          const size_t                 N_w);

/// Parameters used for Shooting method
struct shoot_params
{
    const std::valarray<double>  &z;
    const double                  epsilon;
    const double                  lambda;
    const double                  mstar;
    const double                  r_d;
    const std::valarray<double>  &Vp;
    const size_t                  N_w;
};

static double psi_at_inf (double  E,
                          void   *params);
static bool repeat_lambda(const double                  lambda,
                          double                       &lambda_0,
                          const double                  E,
                          double                       &E0);
static double shoot_wavefunction(std::valarray<double>       &psi,
                                 std::valarray<double>       &chi,
                                 const std::valarray<double> &z,
                                 const double                 E,
                                 const double                 epsilon,
                                 const double                 lambda,
                                 const double                 mstar,
                                 const double                 r_d,
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
    const size_t N_w=99;             // Number of strips in w integration

    const double delta_E = opt.get_numeric_option("dE") * 1e-3*e;    // Energy increment [J]
    const double epsilon = opt.get_numeric_option("epsilon") * eps0; // Permittivity [F/m]
    const double mstar   = opt.get_numeric_option("mass") * me;      // Effective mass [kg]

    const double lambda_start = opt.get_numeric_option("lambdastart") * 1e-10; // Initial Bohr radius [m]
    const double lambda_step  = opt.get_numeric_option("lambdastep")  * 1e-10; // Bohr radius increment [m]
    const double lambda_stop  = opt.get_numeric_option("lambdastop")  * 1e-10; // Final Bohr radius [m]

    std::valarray<double> z; // Spatial location [m]
    std::valarray<double> V; // Confining potential [J]
    read_table_xy("v.r", z, V);

    // Read list of donor (or acceptor) positions
    std::valarray<double> r_d; // [m]
    read_table_x("r_d.r", r_d);

    // Solutions for each donor position
    std::valarray<double> E0(r_d.size());       // Binding energy [J]
    std::valarray<double> lambda_0(r_d.size()); // Bohr radius [m]

    // Perform variational calculation for each donor/acceptor position
    for(unsigned int i_d = 0; i_d < r_d.size(); ++i_d)
    {
        double lambda=lambda_start; // Initial Bohr radius value [m]
        E0[i_d] = e;             // Minimum energy of single donor 1eV

        bool   repeat_flag_lambda;  /* variational flag=>new lambda      */

        // Variational calculation (search over lambda)
        do
        {
            shoot_params params = {z, epsilon, lambda, mstar, r_d[i_d], V, N_w};
            gsl_function f;
            f.function = &psi_at_inf;
            f.params   = &params;
            gsl_root_fsolver *solver = gsl_root_fsolver_alloc(gsl_root_fsolver_brent);

            /* initial energy estimate=minimum potential-binding energy
               of particle to free ionised dopant */
            double Elo = V.min() - e*e/(4*pi*epsilon*lambda);

            // Value for y=f(x) at bottom of search range
            const double y1 = GSL_FN_EVAL(&f,Elo);

            // Find the range in which the solution lies by incrementing the
            // upper limit of the search range until the function changes sign.
            // Since this coarse search can require many iterations, we keep the
            // lower limit fixed to minimise the amount of computation at each step.
            // The Brent algorithm is extremely fast, so it really doesn't matter that
            // the range we find here is large.
            //
            // Note the end stop to prevent infinite loop in absence of solution
            //
            // TODO: Make the cut-off configurable
            double y2 = y1;
            double Ehi = Elo;
            do
            {
                Ehi += delta_E;
                y2=GSL_FN_EVAL(&f, Ehi);
            }while(y1*y2>0);

            double E = (Elo + Ehi)/2;
            gsl_root_fsolver_set(solver, &f, Elo, Ehi);
            int status = 0;

            // Improve the estimate of the solution using the Brent algorithm
            // until we hit a desired level of precision
            do
            {
                status = gsl_root_fsolver_iterate(solver);
                E   = gsl_root_fsolver_root(solver);
                Elo = gsl_root_fsolver_x_lower(solver);
                Ehi = gsl_root_fsolver_x_upper(solver);
                status = gsl_root_test_interval(Elo, Ehi, 1e-12*e, 0);
            }while(status == GSL_CONTINUE);

            // Stop if we've exceeded the cut-off energy
            if(gsl_fcmp(E, V.max(), e*1e-12) == 1)
                std::cerr << "Exceeded Vmax" << std::endl;

            printf("r_d %le lambda %le energy %le meV\n",r_d[i_d],lambda,E/(1e-3*e));

            repeat_flag_lambda=repeat_lambda(lambda,lambda_0[i_d],E,E0[i_d]);

            lambda+=lambda_step; // increments Bohr radius
        }while((repeat_flag_lambda&&(lambda_stop<0))||(lambda<lambda_stop));

        /* Output neutral dopant binding energies (E) and 
           Bohr radii (lambda) in meV and Angstrom respectively */
        const std::valarray<double> E_out = E0*1000.0/e;
        const std::valarray<double> r_d_out = r_d/1e-10;
        const std::valarray<double> lambda_out = lambda_0*1.0e10;
        write_table_xy("e.r", r_d_out, E_out);
        write_table_xy("l.r", r_d_out, lambda_out);

        std::valarray<double> psi(z.size());
        std::valarray<double> chi(z.size());
        const double psi_inf = shoot_wavefunction(psi, chi, z, E0[i_d], epsilon,lambda_0[i_d],mstar,r_d[i_d],V,N_w);

        // Check that wavefunction is tightly bound
        // TODO: Implement a better check
        if(gsl_fcmp(fabs(psi_inf), 0, 1) == 1)
            throw "Warning: Wavefunction is not tightly bound";

        /* generate output filename (and open file for writing) using 
           the basis wf%i.r where the integer %i is the donor index i_d  */
        char   filename[9];     /* character string for wavefunction filename  */
        sprintf(filename,"wf%i.r",i_d);

        write_table_xyz(filename, z, psi, chi);
    }/* end loop over r_d */

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
static double psi_at_inf (double  E,
                          void   *params)
{
    const shoot_params *p = reinterpret_cast<shoot_params *>(params);
    std::valarray<double> psi(p->z.size()); // Wavefunction amplitude
    std::valarray<double> chi(p->z.size()); // Wavefunction envelope amplitude
    const double psi_inf = shoot_wavefunction(psi, chi, p->z, E, p->epsilon,p->lambda,p->mstar,p->r_d,p->Vp,p->N_w);

    return psi_inf;
}

/**
 * \brief Calculates a wavefunction iteratively from left to right of structure
 *
 * \details The value of the wavefunction is taken to be zero at the point
 *          immediately to the left of the potential profile. Subsequent
 *          values are computed using QWWAD3, Eq. 5.28.
 *
 * \param[out] psi     Array to which complete wavefunction will be written [m^{-1/2}]
 * \param[out] chi     Array to which wavefunction envelope will be written [m^{-1/2}]
 * \param[in]  E       Energy at which to compute wavefunction
 * \param[in]  z       Spatial positions [m]
 * \param[in]  V       Potential profile [J]
 * \param[in]  m0      Band-edge effective mass [kg]
 * \param[in]  alpha   Nonparabolicity parameter [1/J]
 *
 * \returns The wavefunction amplitude at the point immediately to the right of the structure
 */
static double shoot_wavefunction(std::valarray<double>       &psi,
                                 std::valarray<double>       &chi,
                                 const std::valarray<double> &z,
                                 const double                 E,
                                 const double                 epsilon,
                                 const double                 lambda,
                                 const double                 mstar,
                                 const double                 r_d,
                                 const std::valarray<double> &Vp,
                                 const size_t                 N_w)
{
    const size_t nz = z.size();
    const double dz = z[1] - z[0];

    psi.resize(nz);
    chi.resize(nz);

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
    double Npsi=integral(pow(psi,2.0),dz); // normalisation integral for psi
    double Nchi=integral(pow(chi,2.0),dz); // normalisation integral for chi

    /* divide unnormalised wavefunction by square root
       of normalisation integral                       */
    psi /= sqrt(Npsi);
    chi /= sqrt(Nchi);

    return chi_next/sqrt(Npsi);
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

    std::valarray<double> dI4_dw(nw); // Function to be integrated

    // Exclude w = 0 point to avoid singularity
    for (unsigned int iw = 1; iw < nw; ++iw)
    {
        const double w = iw * dw;

        dI4_dw[iw] = exp(-z_dash_abs * (1/w-w)/lambda) *
                (1-w*w)/(2*w*w);
    }

    // Compute integral [QWWAD3, Eq. 5.59]
    const double I4=2.0*pi*z_dash_abs * integral(dI4_dw, dw);

    return I4;
}
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :
