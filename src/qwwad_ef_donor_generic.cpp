/**
 * \file   qwwad_ef_donor_generic.cpp
 * \brief  Use variational technique to minimise Hamiltonian of donor state
 * \author Paul Harrison  <p.harrison@shu.ac.uk>
 * \author Alex Valavanis <a.valavanis@leeds.ac.uk>
 */

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <iostream>
#include <sstream>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_min.h>
#include "qwwad/constants.h"
#include "qwwad/file-io.h"
#include "qwwad/maths-helpers.h"
#include "qwwad/options.h"

using namespace QWWAD;
using namespace constants;

/**
 * \brief Identifiers for different states
 */
enum StateID
{
    STATE_1S,
    STATE_2S,
    STATE_2PX,
    STATE_2PZ
};

struct EnergyParams
{
    std::valarray<double> wf;
    std::valarray<double> V;
    std::valarray<double> z;
    double                epsilon;
    double                m;
    double                r_i;
    StateID               S;
};

double Energy(double lambda,
              void   *params);

double Psi(const double  psi,
           const double  lambda,
           const double  x,
           const double  y,
           const double  z,
           const StateID S);

/**
 * \brief Configure command-line options
 */
Options configure_options(int argc, char* argv[])
{
    Options opt;
    std::string doc("Find state of electron attached to a donor in a 2D system using a generic search");

    opt.add_option<double>     ("dcpermittivity,e",   13.18, "Bulk relative permittivity");
    opt.add_option<double>     ("mass,m",             0.067, "Bulk effective mass (relative to free electron)");
    opt.add_option<char>       ("particle,p",           'e', "ID of particle to be used: 'e', 'h' or 'l' for "
                                                             "electrons, heavy holes or light holes respectively.");
    opt.add_option<size_t>     ("subband",                1, "Principal quantum number of subband for which to find impurity state.");
    opt.add_option<std::string>("impuritystate",       "1s", "Symmetry of impurity state (1s, 2s, 2px, or 2pz)");
    opt.add_option<double>     ("donorposition,r",           "Location of donor ion [Angstrom]");

    opt.add_prog_specific_options_and_parse(argc, argv, doc);

    return opt;
}

int main(int argc,char *argv[])
{
    const auto opt = configure_options(argc, argv);

    const auto epsilon = opt.get_option<double>("dcpermittivity") * eps0; // Permittivity [F/m]
    const auto m       = opt.get_option<double>("mass") * me;             // Effective mass [kg]
    const auto p       = opt.get_option<char>  ("particle");              // Particle ID (e, h, or l)
    const auto subband = opt.get_option<size_t>("subband");               // Principal quantum number of state to find

    StateID S = STATE_1S;

    const auto impuritystate_string = opt.get_option<std::string>("impuritystate");

    if (impuritystate_string == "1s")
        S = STATE_1S;
    else if (impuritystate_string == "2s")
        S = STATE_2S;
    else if (impuritystate_string == "2px")
        S = STATE_2PX;
    else if (impuritystate_string == "2pz")
        S = STATE_2PZ;
    else
    {
        std::cerr << "Unknown impurity state ID: " << impuritystate_string << std::endl;
        exit(EXIT_FAILURE);
    }

    std::valarray<double> z; // Spatial location [m]
    std::valarray<double> V; // Confining potential [J]
    read_table("v.r", z, V);

    std::ostringstream filename; // input filename
    filename << "wf_" << p << subband << ".r";

    std::valarray<double> z_tmp; // Dummy file for unused spatial locations
    std::valarray<double> wf;    // Wave function samples at each point [m^{-1/2}]
    read_table(filename.str(), z_tmp, wf);

    const auto lambda_0=4*pi*epsilon*(hBar/e)*(hBar/e)/m;/* Bohr	theory (1s)	*/

    /* Open files for output of data */

    FILE *fe=fopen("e.r","w");           /* E versus r_i	*/
    FILE *fl=fopen("l.r","w");           /* lambda versus r_i	*/

    auto r_d = (z[z.size()-1] + z[0])/2.0;

    if (opt.get_argument_known("donorposition") > 0)
        r_d = opt.get_option<double>("donorposition") * 1e-10;

    // Perform variational calculation for each donor/acceptor position
    double lambda=lambda_0;	// initial lambda guess

    // Double the estimate of Bohr radius if we're in a second orbital
    // This isn't correct for 2pz, but it's still better than the 1s
    // estimate!
    if((S==STATE_2S)||(S==STATE_2PX)||(S==STATE_2PZ))lambda*=2;

    /* Newton-Raphson iteration for solution of lambda, this occurs when
       dE/dlambda=0, hence the function f is dE/dlambda and f'=d2E/dlambda^2
       */
    EnergyParams params = {wf, V, z, epsilon, m, r_d, S};

    // Set up the numerical solver using GSL
    gsl_function f;
    f.function = &Energy;
    f.params   = &params;

    gsl_min_fminimizer *s = gsl_min_fminimizer_alloc(gsl_min_fminimizer_brent);
    gsl_min_fminimizer_set(s, &f, lambda, lambda/5, lambda*10);

    size_t max_iter = 100; // Maximum number of iterations before giving up
    int status = 0;        // Error flag for GSL
    unsigned int iter=0;   // The number of iterations attempted so far

    double E = 1000*e; // Minimum energy of carrier [J]

    // Variational calculation (search over lambda)
    do
    {
        ++iter;
        status  = gsl_min_fminimizer_iterate(s);
        const double lambda_lo = gsl_min_fminimizer_x_lower(s);
        const double lambda_hi = gsl_min_fminimizer_x_upper(s);
        lambda = gsl_min_fminimizer_x_minimum(s);
        E      = gsl_min_fminimizer_f_minimum(s);
        status = gsl_min_test_interval(lambda_lo, lambda_hi, 0.1e-10, 0.0);
        printf("r_d %le lambda %le energy %le meV\n", r_d, lambda, E/(1e-3*e));
    }while((status == GSL_CONTINUE) && (iter < max_iter));

    gsl_min_fminimizer_free(s);

    /* Output total energy (E) of impurity/heterostructure system 
       and Bohr radii (lambda), in meV and Angstrom respectively */
    fprintf(fe,"%le %le\n", r_d/1e-10,E/(1e-3*e));
    fprintf(fl,"%le %le\n", r_d/1e-10,lambda/1e-10);

    fclose(fe);
    fclose(fl);

    return EXIT_SUCCESS;
} /* end main */

/* This function calculates the expectation value (the energy) of the
   Hamiltonian operator	*/
double Energy(double  lambda,
              void   *params)
{
    auto *p = reinterpret_cast<EnergyParams *>(params);
    const auto dz  = p->z[1] - p->z[0]; // z- (growth) direction step length [m]
    const auto nz  = p->z.size();       // Number of spatial samples in z direction
    const auto dxy = lambda/10;         // Step size for in-plane integration [m]
    const auto nxy = 31;                // Number of samples to use in integration over x and y

    // Integrands wrt z for calculating wavefunction overlap
    // and Hamiltonian
    std::valarray<double> PD_integrand_z(nz);
    std::valarray<double> H_integrand_z(nz);

    // Pre-calculate a couple of params to speed things up
    const auto hBar_sq_by_2m  = hBar*hBar/(2.0*p->m);
    const auto e_sq_by_4pieps = e*e/(4.0*pi*p->epsilon);

    // Compute integrand over the z-axis, skipping both end-points since we
    // need the 2nd derivatives
    for(unsigned int iz=1; iz < nz-1; ++iz)
    {
        // Integrands wrt (x,z) for calculating wavefunction overlap
        // and Hamiltonian
        std::valarray<double> PD_integrand_xz(nxy);
        std::valarray<double> H_integrand_xz(nxy);

        const auto z_dash = p->z[iz] - p->r_i; // Separation from donor in z-direction [m]

        for(unsigned int ix=0; ix<nxy; ++ix)	
        {
            const auto x = ix*dxy;

            // Integrands wrt (x,y,z) for calculating wavefunction overlap
            // and Hamiltonian
            std::valarray<double> PD_integrand_xyz(nxy);
            std::valarray<double> H_integrand_xyz(nxy);

            // Wavefunction at (x, y - dy, z_dash)
            auto Psixyz_last_y = Psi(p->wf[iz], lambda, x, -dxy, z_dash, p->S);

            // Wavefunction at (x,y,z_dash)
            auto Psixyz = Psi(p->wf[iz],lambda, x, 0, z_dash, p->S);

            const auto r_xz = hypot(x, z_dash);
            
            for(unsigned int iy=0; iy<nxy; ++iy)
            {
                const auto y = iy*dxy;

                // Wavefunction at (x, y+dy, z_dash)
                const auto Psixyz_next_y = Psi(p->wf[iz],lambda,x,y+dxy, z_dash, p->S);

                // Calculate the second derivatives along x, y and z
                const auto d2Pdx2=(Psi(p->wf[iz],lambda,x+dxy,y, z_dash, p->S)-
                        2*Psixyz+
                        Psi(p->wf[iz],lambda,x-dxy,y, z_dash, p->S))/(dxy*dxy);

                const auto d2Pdy2=(Psixyz_next_y - 2*Psixyz + Psixyz_last_y)/(dxy*dxy);

                const auto d2Pdz2=(Psi(p->wf[iz+1],lambda,x,y,p->z[iz+1]-p->r_i,p->S)-
                        2*Psixyz+
                        Psi(p->wf[iz-1],lambda,x,y,p->z[iz-1]-p->r_i,p->S))/(dz*dz);

                // Distance from impurity for Coloumb term [m]
                const auto r = hypot(y, r_xz);

                // The Laplacian of Psi
                const auto laplace_Psi = d2Pdx2 + d2Pdy2 + d2Pdz2;

                // The integrand for the Hamiltonian expectation value
                // QWWAD 3, Eq. 5.142
                H_integrand_xyz[iy] = Psixyz*(-hBar_sq_by_2m*laplace_Psi
                        +(p->V[iz]-e_sq_by_4pieps/r)*Psixyz);

                PD_integrand_xyz[iy] = Psixyz*Psixyz;

                // Reuse values for next iteration
                Psixyz_last_y = Psixyz;
                Psixyz        = Psixyz_next_y;
            }

            // Approximate the singularity at r=0 with value at (r=dy)
            // Note that this preserves the symmetry of the function around (x,y) = 0.
            if(ix==0 && abs(z_dash) < dz)
                H_integrand_xyz[0] = H_integrand_xyz[1];

            // Perform integration over y, noting that a factor of 2 is included
            // to account for even symmetry
            H_integrand_xz[ix]  = 2*simps(H_integrand_xyz, dxy);
            PD_integrand_xz[ix] = 2*simps(PD_integrand_xyz, dxy);
        }

        // Perform integration over x, noting that a factor of 2 is included
        // to account for even symmetry
        PD_integrand_z[iz] = 2*simps(PD_integrand_xz, dxy);
        H_integrand_z[iz]  = 2*simps(H_integrand_xz, dxy);
    }

    // Note that endpoints of the integral can keep their default value of zero, since
    // psi decays to zero at infinity

    // Compute the final value of the energy using Eq. 5.141, QWWAD3
    const auto H_exp = integral(H_integrand_z, dz);
    const auto norm  = integral(PD_integrand_z, dz);
    const auto E = H_exp/norm;

    return E;
}

/**
 * \brief The wave function psi(z)phi(r)
 */
double Psi(const double  psi,
           const double  lambda,
           const double  x,
           const double  y,
           const double  z,
           const StateID S)
{
    const auto r_xy = hypot(x,y);
    const auto r    = hypot(r_xy,z);

    double result = 0.0;

    switch(S)
    {
        case STATE_1S:
            result = psi*exp(-r/lambda);
            break;
        case STATE_2S:
            result = psi*(1-r/lambda)*exp(-r/lambda);
            break;
        case STATE_2PX:
            result = psi*fabs(x)*exp(-r/lambda);
            break;
        case STATE_2PZ:
            result = psi*fabs(z)*exp(-r/lambda);
            break;
        default:
            std::cerr << "Unrecognised orbital" << std::endl;
            exit(EXIT_FAILURE);
    }

    return result;
}
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :
