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
#include "struct.h"
#include "qclsim-constants.h"
#include "qclsim-fileio.h"
#include "qclsim-maths.h"

using namespace Leeds;
using namespace constants;

double Energy(const std::valarray<double> &wf,
              const std::valarray<double> &V,
              const std::valarray<double> &z,
              double	epsilon,
              double	m,
              double	lambda,
              double	r_i,
              int	S);
double Psi(const double psi,
           const double lambda,
           const double x,
           const double y,
           const double z,
           const int    S);

int main(int argc,char *argv[])
{
double	epsilon;	/* permitivity of material		*/
double	f,fdash;	/* function (dE/dlambda) and derivative 
			   to be solved				*/
double	lambda;		/* Bohr radius (variational)		*/
double	lambda_step;	/* Bohr radius increment		*/
double	lambda_0;	/* Bohr radius of bulk impurity		*/
double	m;		/* mass of particle			*/
double	y1,y2,y3;	/* E(lambda-), E(lambda), E(lambda+)	*/
int	state;		/* principal quantum number		*/
int	S;		/* impurity level, i.e. `1s', `2px' etc	*/
char	p;		/* particle				*/
char	State[9];	/* string containing impurity level	*/
FILE   *fe;		/* file pointer for energies		*/
FILE   *fl;		/* file pointer for lambda_0		*/

/* default values */
epsilon=13.18*eps0;
m=0.067*me;
p='e';
state=1;
S=1;

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
    read_table_xy("v.r", z, V);

    char	filename[9];	/* input filename			*/
    sprintf(filename,"wf_%c%i.r",p,state);
    std::valarray<double> z_tmp; // Dummy file for unused spatial locations
    std::valarray<double> wf;    // Wave function samples at each point [m^{-1/2}]
    read_table_xy(filename, z_tmp, wf);

  lambda_0=4*pi*epsilon*(hBar/e)*(hBar/e)/m;/* Bohr	theory (1s)	*/

  /* Open files for output of data */

  fe=fopen("e.r","w");           /* E versus r_i	*/
  fl=fopen("l.r","w");           /* lambda versus r_i	*/

  // Read list of donor (or acceptor) positions
  std::valarray<double> r_d; // [m]
  read_table_x("r_d.r", r_d);

  // Perform variational calculation for each donor/acceptor position
  for(unsigned int i_d = 0; i_d < r_d.size(); ++i_d)
  {
   lambda=lambda_0;				/* initial lambda guess	*/
   if((S==2)||(S==3))lambda*=2;			/* for 2s, 2px NOT 2pz	*/

   lambda_step=lambda/100;

   /* Newton-Raphson iteration for solution of lambda, this occurs when
      dE/dlambda=0, hence the function f is dE/dlambda and f'=d2E/dlambda^2
    								*/
   do
   {
    y1=Energy(wf,V,z,epsilon,m,lambda-lambda_step,r_d[i_d],S);
    y2=Energy(wf,V,z,epsilon,m,lambda,r_d[i_d],S);
    y3=Energy(wf,V,z,epsilon,m,lambda+lambda_step,r_d[i_d],S);

    f=(y3-y1)/(2*lambda_step);
    fdash=(y3-2*y2+y1)/(lambda_step*lambda_step);

    printf("r_i %4.2f A lambda %4.2f A energy %4.3f meV\n",
            r_d[i_d]/1e-10,lambda/1e-10,y2/(1e-3*e));

    lambda-=f/fdash;

   }while(fabs(f/fdash)>1e-10);	

   /* Output total energy (E) of impurity/heterostructure system 
      and Bohr radii (lambda), in meV and Angstrom respectively */
   fprintf(fe,"%le %le\n",r_d[i_d]/1e-10,y2/(1e-3*e));
   fprintf(fl,"%le %le\n",r_d[i_d]/1e-10,lambda/1e-10);
  }/* end while r_i */

  fclose(fe);
  fclose(fl);

  return EXIT_SUCCESS;
} /* end main */

/* This function calculates the expectation value (the energy) of the
   Hamiltonian operator	*/
double Energy(const std::valarray<double> &wf,
              const std::valarray<double> &V,
              const std::valarray<double> &z,
              double	epsilon,
              double	m,
              double	lambda,
              double	r_i,
              int	S)
{
 const double dz  = z[1] - z[0]; // z- (growth) direction step length [m]
 const size_t nz  = V.size();    // Number of spatial samples in z direction
 const double dxy = lambda/10;   // Step size for in-plane integration [m]
 const size_t nxy = 31;          // Number of samples to use in integration over x and y

 // Integrands wrt z for calculating wavefunction overlap
 // and Hamiltonian
 std::valarray<double> PD_integrand_z(nz);
 std::valarray<double> H_integrand_z(nz);

 // Compute integrand over the z-axis, skipping both end-points since we
 // need the 2nd derivatives
 for(unsigned int iz=1;iz < nz-1;iz++)
 {
     // Integrands wrt (x,z) for calculating wavefunction overlap
     // and Hamiltonian
     std::valarray<double> PD_integrand_xz(nxy);
     std::valarray<double> H_integrand_xz(nxy);

     // integrate over the plane, skipping singularity at x=0
     // Note that we correct for this later
     for(unsigned int ix=1; ix<nxy; ++ix)	
     {
         const double x = ix*dxy;

         // Integrands wrt (x,y,z) for calculating wavefunction overlap
         // and Hamiltonian
         std::valarray<double> PD_integrand_xyz(nxy);
         std::valarray<double> H_integrand_xyz(nxy);

         for(unsigned int iy=1; iy<nxy; ++iy)
         {
             const double y = iy*dxy;

             // Wavefunction at this point
             const double Psixyz = Psi(wf[iz],lambda,x,y,z[iz]-r_i,S);

             // Calculate the second derivatives along x, y and z
             const double d2Pdx2=(Psi(wf[iz],lambda,x+dxy,y,z[iz]-r_i,S)-
                     2*Psixyz+
                     Psi(wf[iz],lambda,x-dxy,y,z[iz]-r_i,S))/(dxy*dxy);

             const double d2Pdy2=(Psi(wf[iz],lambda,x,y+dxy,z[iz]-r_i,S)-
                     2*Psixyz+
                     Psi(wf[iz],lambda,x,y-dxy,z[iz]-r_i,S))/(dxy*dxy);

             const double d2Pdz2=(Psi(wf[iz+1],lambda,x,y,z[iz+1]-r_i,S)-
                     2*Psixyz+
                     Psi(wf[iz-1],lambda,x,y,z[iz-1]-r_i,S))/(dz*dz);

             // Distance from impurity for Coloumb term [m]
             const double r=gsl_hypot3(x, y, (z[iz]-r_i));

             // The Laplacian of Psi
             const double laplace_Psi = d2Pdx2 + d2Pdy2 + d2Pdz2;

             // The integrand for the Hamiltonian expectation value
             // QWWAD 3, Eq. 5.142
             H_integrand_xyz[iy] = Psixyz*(-hBar*hBar/(2*m)*laplace_Psi
                     +(V[iz]-e*e/(4*pi*epsilon*r))*Psixyz);

             PD_integrand_xyz[iy] = Psixyz*Psixyz;
         }
         // Approximate the singularities with values at neighbouring point
         // Note that this preserves the symmetry of the functions around 0.
         PD_integrand_xyz[0] = PD_integrand_xyz[1];
         H_integrand_xyz[0]  = H_integrand_xyz[1];

      // Perform integration over y, noting that a factor of 2 is included
      // to account for even symmetry
      H_integrand_xz[ix]  = 2*simps(H_integrand_xyz, dxy);
      PD_integrand_xz[ix] = 2*simps(PD_integrand_xyz, dxy);
  }

  // Approximate the singularities with values at neighbouring point
  // Note that this preserves the symmetry of the functions around 0.
  PD_integrand_xz[0] = PD_integrand_xz[1];
  H_integrand_xz[0]  = H_integrand_xz[1];

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
        default:
            fprintf(stderr, "Unrecognised orbital\n");
            exit(EXIT_FAILURE);
    }

    return result;
}
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :
