/**
 * \file    gde.cpp
 * \brief   General Diffusion Equation
 * \author  Paul Harrison <p.harrison@shu.ac.uk>
 * \author  Alex Valavanis <a.valavanis@leeds.ac.uk>
 *
 * \details Produces the general solution to the diffusion
 *          equation:
 *          \f[
 *            \frac{\partial n}{\partial t} = 
 *                \frac{\partial}{\partial x}
 *                D \frac{\partial n}{\partial x}
 *          \f]
 *
 *          for \f$n=n(x,t)\f$ and \f$D=D(x,t,n)\f$.
 *
 *  Input files:
 *    x.r           initial (t=0) concentration profile versus z  
 *
 *  Output files:
 *    X.r           final (diffused) concentration profile 
 */

#include <iostream>
#include <cstdlib>

#include <gsl/gsl_math.h>

#include "qwwad/file-io.h"
#include "qwwad/options.h"

using namespace QWWAD;

static void diffuse(const std::valarray<double> &z,
                    std::valarray<double>       &x,
                    const std::valarray<double> &D,
                    const double                delta_t);

static void check_stability(const double dt,
                            const double dz,
                            const double D)
{
    const double dt_max = dz*dz/(2*D);

    if (dt > dt_max)
    {
        std::cerr << "User-specified time step (dt = " << dt << " s) exceeds stability criterion (dt < " << dt_max << " s). "
                  << "You can fix this by choosing a lower value using the --dt option, or by increasing the spatial-step size in your "
                  << "input data files." << std::endl;
        exit(EXIT_FAILURE);
    }
}

int main(int argc,char *argv[])
{
    Options opt;
    std::string doc("Solve the generalised diffusion equation");

    opt.add_option<double>     ("dt,d",          0.01, "Time-step [s]");
    opt.add_option<double>     ("coeff,D",        1.0, "Diffusion coefficient [Angstrom^2/s]");
    opt.add_option<double>     ("time,t",         1.0, "End time for simulation [s]");
    opt.add_option<std::string>("mode,a",  "constant", "Form of diffusion coefficient");
    opt.add_option<std::string>("infile",       "x.r", "File from which input profile of diffusant will be read");
    opt.add_option<std::string>("outfile",      "X.r", "File to which output profile of diffusant will be written");

    opt.add_prog_specific_options_and_parse(argc, argv, doc);

    const auto t_final = opt.get_option<double>("time");          // [s]
    const auto dt      = opt.get_option<double>("dt");            // [s]
    const auto D0      = opt.get_option<double>("coeff") * 1e-20; // [m^2/s]
    const auto mode    = opt.get_option<std::string>("mode");

    std::valarray<double> z; // Spatial location [m]
    std::valarray<double> x; // Initial diffusant profile
    read_table(opt.get_option<std::string>("infile").c_str(), z, x);

    const size_t nz = z.size(); // Number of spatial points

    std::valarray<double> D(nz); // Diffusion coefficient

    for(double t=dt; t<=t_final; t+=dt)
    {
        if (mode == "constant")
        {
            D = D0;	// set constant diffusion coeff.
        }
        else if(mode == "concentration-dependent")
        {
            // TODO: Make this configurable
            const double k = 1e-20; // Concentration factor [m^2/s]

            // Find concentration-dependent diffusion coefficient
            // [4.14, QWWAD4]
            D = k*x*x;
        }
        else if(mode == "depth-dependent")
        {
            // TODO: Make this configurable
            const double D0    = 10*1e-20;   // Magnitude of distribution [m^2/s]
            const double z0    = 1800*1e-10; // Centre of diff. coeff. distribution [m]
            const double sigma = 600*1e-10;  // Width of distribution [m]

            // Find depth-dependent diffusion coefficient
            // [4.16, QWWAD4]
            D = D0*exp(-pow((z-z0)/sigma, 2)/2);
        }
        else if(mode == "time-dependent")
        {
            // TODO: Make this configurable
            const double D0    = 10*1e-20;   // Magnitude of distribution [m^2/s]
            const double z0    = 1800*1e-10; // Centre of diff. coeff. distribution [m]
            const double sigma = 600*1e-10;  // Width of distribution [m]
            const double tau   = 100;        // Decay time-constant for diffusion [s]

            // Find time and depth-dependent diffusion coefficient
            // [4.18, QWWAD4]
            D = D0*exp(-pow((z-z0)/sigma, 2)/2)*exp(-t/tau);
        }
        else
        {
            std::cerr << "Diffusion mode: " << mode << " not recognised" << std::endl;
            exit(EXIT_FAILURE);
        }

        diffuse(z, x, D, dt);
    }

    write_table(opt.get_option<std::string>("outfile").c_str(), z, x);

    return EXIT_SUCCESS;
}

/**
 * Projects the diffusant profile a short time interval delta_t into the future
 *
 * \param[in]     z        spatial profile [m]
 * \param[in,out] x        diffusant profile
 * \param[in]     D        Diffusion coefficient at each point [m^2/s]
 * \param[in]     delta_t  time step [s]
 */
static void diffuse(const std::valarray<double> &z,
                    std::valarray<double>       &x,
                    const std::valarray<double> &D,
                    const double                 delta_t)
{
    const double dz = z[1] - z[0];
    const size_t nz = z.size();

    check_stability(delta_t, dz, D.max());

    std::valarray<double> x_new(nz); // Modified diffusion profile

    for(unsigned int iz=1; iz<nz-1; ++iz)
    {
        x_new[iz]=delta_t*
            (
             (D[iz+1]-D[iz-1]) * (x[iz+1]-x[iz-1])/gsl_pow_2(2*dz)
             +D[iz] * (x[iz+1]-2*x[iz]+x[iz-1])/gsl_pow_2(dz)
            )
            + x[iz];
    }

    /* Impose `closed-system' boundary conditions. See section 4.3, QWWAD3 */
    x_new[0]    = x_new[1];
    x_new[nz-1] = x_new[nz-2];

    x = x_new; // Copy new profile
}
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :
