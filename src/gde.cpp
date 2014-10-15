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
 *
 *  Paul Harrison, July 1994
 *	
 *  Modifications June 1998
 */

#include <iostream>
#include <cstdlib>

#include "qclsim-fileio.h"
#include "qwwad-options.h"
#include "dox.h"

using namespace Leeds;

void          calculate_D(const std::valarray<double> &z,
                          const std::valarray<double> &x,
                          std::valarray<double>       &D,
                          double                       t);
static void   diffuse    (const std::valarray<double> &z,
                          std::valarray<double>       &x,
                          const std::valarray<double> &D,
                          const double                delta_t);

int main(int argc,char *argv[])
{
    Options opt;
    std::string doc("Solve the generalised diffusion equation");

    opt.add_numeric_option("dt,d",          0.01, "Time-step [s]");
    opt.add_numeric_option("coeff,D",        1.0, "Diffusion coefficient [Angstrom^2/s]");
    opt.add_numeric_option("time,t",         1.0, "End time for simulation [s]");
    opt.add_string_option ("mode,a",  "constant", "Form of diffusion coefficient");

    opt.add_prog_specific_options_and_parse(argc, argv, doc);

    const double t_final   = opt.get_numeric_option("time");          // [s]
    const double dt        = opt.get_numeric_option("dt");            // [s]
    const double D0        = opt.get_numeric_option("coeff") * 1e-20; // [m^2/s]
    std::string  mode      = opt.get_string_option("mode");

    std::valarray<double> z; // Spatial location [m]
    std::valarray<double> x; // Initial diffusant profile
    read_table_xy("x.r", z, x);

    const size_t n  = z.size(); // Number of spatial points

    std::valarray<double> D(n); // Diffusion coefficient

    if (mode == "constant")
    {
        for(unsigned int iz=0; iz<n; ++iz)
            D[iz] = D0;	// set constant diffusion coeff.

        for(double t=dt; t<=t_final; t+=dt)
            diffuse(z, x, D, dt);
    }
    else if(mode == "file")
    {
        read_table_xy("D.r", z, D); // read D from file

        for(double t=dt; t<=t_final; t+=dt)
            diffuse(z, x, D, dt);
    }
    else if(mode == "general")
    {
        for(double t=dt; t<=t_final; t+=dt)
        {
            calculate_D(z, x, D, t);	/* calculate D using `dox.c'	*/
            diffuse(z, x, D, dt);
        }
    }
    else if(mode == "time")
    {
        read_table_xy("D.r", z, D); // read D at t=0 from file

        for(double t=dt; t<=t_final; t+=dt)
        {
            calculate_D(z, x, D, t);	/* calculate subsequent D */
            diffuse(z, x, D, dt);
        }
    }
    else
    {
        std::cerr << "Diffusion mode: " << mode << " not recognised" << std::endl;
        exit(EXIT_FAILURE);
    }

    write_table_xy("X.r", z, x);

    return EXIT_SUCCESS;
}

/**
 * \brief recalculates the diffusion coefficient
 *
 * \param[in]     z spatial profile [m]
 * \param[in]     x diffusant profile
 * \param[in,out] D diffusion coefficient at each point [m^2/s]
 * \param[in]     t time
 *
 * \details recalculates D for all points along the z-axis when D is a function of the concentration
 */
void calculate_D(const std::valarray<double> &z,
                 const std::valarray<double> &x,
                 std::valarray<double>       &D,
                 double                       t)
{
    for(unsigned int iz=0; iz<z.size(); ++iz)
        D[iz] = D_of_x(D[iz], x[iz], z[iz], t);
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
