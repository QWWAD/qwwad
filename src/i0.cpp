/*==================================================================
                                i0
  ==================================================================*/

/* This program implements a variational technique to calculate the
   uncorrelated one particle energies of an electron or hole attatched
   to a single donor or acceptor at any position, in any user supplied 
   potential.  The potential is read from the file v.r

   This version is a single variational parameter calculation---the 
   three-dimensional (3D) approximation trial wavefunction:

		Psi=psi(z) phi(r)

   where psi(z) is the one-particle wave function without the impurity
   and phi(r) is the hydorgenic-like term, which can take any form, i.e.
   1s, 2s, 2pz, 2px.

		Input files:
		r_i.r		donor (or acceptor positions)
		v.r		one-dimensional potential
		wf_pn.r		wave functions

		Output files:
		e.r		total energies for each r_i
		l.r		Bohr radii (lambda) for each r_i

   Paul Harrison, June 2001
   								*/

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_min.h>
#include "struct.h"
#include "qclsim-constants.h"
#include "qclsim-fileio.h"
#include "qclsim-maths.h"

using namespace Leeds;
using namespace constants;

struct EnergyParams
{
    std::valarray<double> wf;
    std::valarray<double> V;
    std::valarray<double> z;
    double                epsilon;
    double                m;
    double                r_i;
    int	                  S;
};

double Energy(double lambda,
              void   *params);

double Psi(const double psi,
           const double lambda,
           const double x,
           const double y,
           const double z,
           const int    S);

int main(int argc,char *argv[])
{
char	State[9];	/* string containing impurity level	*/

/* default values */
double epsilon = 13.18*eps0; // Permittivity [F/m]
double m       = 0.067*me;   // Mass of particle [kg]
char   p       = 'e';        // Particle ID
int    state   = 1;          // Principal quantum number
int    S       = 1;          // Impurity level `1s', `2px', etc

while((argc>1)&&(argv[1][0]=='-'))
{
 switch(argv[1][1])
 {
  case 'e':
	   epsilon=atof(argv[2])*eps0;
	   break;
  case 'm':
	   m=atof(argv[2])*me;
	   break;
  case 's':
	   state=atoi(argv[2]);
	   break;
  case 'S':
	   strcpy(State,argv[2]);
           if(!strcmp(State,"1s"))S=9;
           if(!strcmp(State,"2s"))S=2;
           if(!strcmp(State,"2px"))S=3;
           if(!strcmp(State,"2pz"))S=4;
           switch(S)
           {
	    case 2 :break;
	    case 3 :break;
	    case 4 :break;
	    case 9 :S=1;break;
	    default:
		    printf("The `%s' impurity level wave function is not defined\n",State);
		    exit(0);
		    break;
	   }
	   break;
  case 'p':
	   p=*argv[2];
	   switch(p)
	   {
	    case 'e': break;
	    case 'h': break;
	    case 'l': break;
	    default:  printf("Usage:  i0 [-p particle (e, h, or l)]\n");
                      exit(0);
	   }
	   break;
  default :
           printf("Usage:  i0 [-e relative permittivity \033[1m13.18\033[0m]\n");
	   printf("           [-m mass (\033[1m0.067\033[0mm0)]\n");
	   printf("           [-s subband (\033[1m1\033[0m)][-S impurity level (\033[1m1s\033[0m)]\n");
	   exit(0);
 }
 argv++;
 argv++;
 argc--;
 argc--;
}

    std::valarray<double> z; // Spatial location [m]
    std::valarray<double> V; // Confining potential [J]
    read_table("v.r", z, V);

    char	filename[9];	/* input filename			*/
    sprintf(filename,"wf_%c%i.r",p,state);
    std::valarray<double> z_tmp; // Dummy file for unused spatial locations
    std::valarray<double> wf;    // Wave function samples at each point [m^{-1/2}]
    read_table(filename, z_tmp, wf);

  const double lambda_0=4*pi*epsilon*(hBar/e)*(hBar/e)/m;/* Bohr	theory (1s)	*/

  /* Open files for output of data */

  FILE *fe=fopen("e.r","w");           /* E versus r_i	*/
  FILE *fl=fopen("l.r","w");           /* lambda versus r_i	*/

  // Read list of donor (or acceptor) positions
  std::valarray<double> r_d; // [m]
  read_table("r_d.r", r_d);

  // Perform variational calculation for each donor/acceptor position
  for(unsigned int i_d = 0; i_d < r_d.size(); ++i_d)
  {
   double lambda=lambda_0;	// initial lambda guess

   // Double the estimate of Bohr radius if we're in a second orbital
   // This isn't correct for 2pz, but it's still better than the 1s
   // estimate!
   if((S==2)||(S==3)||(S==4))lambda*=2;

   /* Newton-Raphson iteration for solution of lambda, this occurs when
      dE/dlambda=0, hence the function f is dE/dlambda and f'=d2E/dlambda^2
    								*/
   EnergyParams params = {wf, V, z, epsilon, m, r_d[i_d], S};

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
       printf("r_d %le lambda %le energy %le meV\n", r_d[i_d], lambda, E/(1e-3*e));
   }while((status == GSL_CONTINUE) && (iter < max_iter));

   gsl_min_fminimizer_free(s);
   /* Output total energy (E) of impurity/heterostructure system 
      and Bohr radii (lambda), in meV and Angstrom respectively */
   fprintf(fe,"%le %le\n",r_d[i_d]/1e-10,E/(1e-3*e));
   fprintf(fl,"%le %le\n",r_d[i_d]/1e-10,lambda/1e-10);
  }/* end while r_i */

  fclose(fe);
  fclose(fl);

  return EXIT_SUCCESS;
} /* end main */

/* This function calculates the expectation value (the energy) of the
   Hamiltonian operator	*/
double Energy(double  lambda,
              void   *params)
{
    EnergyParams *p = reinterpret_cast<EnergyParams *>(params);
    const double dz  = p->z[1] - p->z[0]; // z- (growth) direction step length [m]
    const size_t nz  = p->z.size();       // Number of spatial samples in z direction
    const double dxy = lambda/10;         // Step size for in-plane integration [m]
    const size_t nxy = 31;                // Number of samples to use in integration over x and y

    // Integrands wrt z for calculating wavefunction overlap
    // and Hamiltonian
    std::valarray<double> PD_integrand_z(nz);
    std::valarray<double> H_integrand_z(nz);

    // Pre-calculate a couple of params to speed things up
    const double hBar_sq_by_2m  = hBar*hBar/(2.0*p->m);
    const double e_sq_by_4pieps = e*e/(4.0*pi*p->epsilon);

    // Compute integrand over the z-axis, skipping both end-points since we
    // need the 2nd derivatives
    for(unsigned int iz=1;iz < nz-1;iz++)
    {
        // Integrands wrt (x,z) for calculating wavefunction overlap
        // and Hamiltonian
        std::valarray<double> PD_integrand_xz(nxy);
        std::valarray<double> H_integrand_xz(nxy);

        const double z_dash = p->z[iz] - p->r_i; // Separation from donor in z-direction [m]

        for(unsigned int ix=0; ix<nxy; ++ix)	
        {
            const double x = ix*dxy;

            // Integrands wrt (x,y,z) for calculating wavefunction overlap
            // and Hamiltonian
            std::valarray<double> PD_integrand_xyz(nxy);
            std::valarray<double> H_integrand_xyz(nxy);

            // Wavefunction at (x, y - dy, z_dash)
            double Psixyz_last_y = Psi(p->wf[iz], lambda, x, -dxy, z_dash, p->S);

            // Wavefunction at (x,y,z_dash)
            double Psixyz = Psi(p->wf[iz],lambda, x, 0, z_dash, p->S);

            const double rsq_xz = x*x + z_dash*z_dash;

            for(unsigned int iy=0; iy<nxy; ++iy)
            {
                const double y = iy*dxy;

                // Wavefunction at (x, y+dy, z_dash)
                const double Psixyz_next_y = Psi(p->wf[iz],lambda,x,y+dxy, z_dash, p->S);

                // Calculate the second derivatives along x, y and z
                const double d2Pdx2=(Psi(p->wf[iz],lambda,x+dxy,y, z_dash, p->S)-
                        2*Psixyz+
                        Psi(p->wf[iz],lambda,x-dxy,y, z_dash, p->S))/(dxy*dxy);

                const double d2Pdy2=(Psixyz_next_y - 2*Psixyz + Psixyz_last_y)/(dxy*dxy);

                const double d2Pdz2=(Psi(p->wf[iz+1],lambda,x,y,p->z[iz+1]-p->r_i,p->S)-
                        2*Psixyz+
                        Psi(p->wf[iz-1],lambda,x,y,p->z[iz-1]-p->r_i,p->S))/(dz*dz);

                // Distance from impurity for Coloumb term [m]
                const double r = sqrt(y*y + rsq_xz);

                // The Laplacian of Psi
                const double laplace_Psi = d2Pdx2 + d2Pdy2 + d2Pdz2;

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
            if(ix==0 && z_dash == 0)
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
    const double H_exp = integral(H_integrand_z, dz);
    const double norm  = integral(PD_integrand_z, dz);
    const double E = H_exp/norm;

    return E;
}

/**
 * \brief The wave function psi(z)phi(r)
 */
double Psi(const double psi,
           const double lambda,
           const double x,
           const double y,
           const double z,
           const int    S)
{
    const double r = sqrt(x*x+y*y+z*z);

    double result = 0.0;

    switch(S)
    {
        case 1:
            result = psi*exp(-r/lambda); /* 1s	*/
            break;
        case 2:
            result = psi*(1-r/lambda)*exp(-r/lambda);	/* 2s	*/
            break;
        case 3:
            result = psi*fabs(x)*exp(-r/lambda); /* 2px	*/
            break;
        case 4:
            result = psi*fabs(z)*exp(-r/lambda); /* 2pz	*/
            break;
        default:
            fprintf(stderr, "Unrecognised orbital\n");
            exit(EXIT_FAILURE);
    }

    return result;
}
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :
