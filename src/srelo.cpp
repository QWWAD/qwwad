/*==================================================================
              srelo Scattering Rate Electron-LO phonon
  ==================================================================*/

/* This program calculates the electron-LO phonon scattering rate for
   both intra- and intersubband events.  The required rates are provided
   by the user in the file `rrp.r'.  The other necessary inputs are listed
   below.


	Input files:		rrp.r	contains required rates
				wf_xy.r	x=particle y=state
				Ex.r	x=particle, energies

	Output files:		LO[a,e]if.r absorption and emission rates


    Paul Harrison, January 1998 
 
    Improvements January 1999						*/

#include <complex>
#include <cstdio>
#include <cstdlib>
#include <math.h>
#include <gsl/gsl_math.h>
#include "struct.h"
#include "qclsim-constants.h"
#include "qclsim-fileio.h"
#include "qclsim-subband.h"
using namespace Leeds;
using namespace constants;

static void ff_table(const double   dKz,
                     const Subband &isb,
                     const Subband &fsb,
                     unsigned int   nKz,
                     std::valarray<double> &Kz,
                     std::valarray<double> &Gifsqr);

void ff_output(const std::valarray<double> &Kz,
               const std::valarray<double> &Gifsqr,
               unsigned int        i,
               unsigned int        f);

static double Gsqr(const double   Kz,
                   const Subband &isb,
                   const Subband &fsb);
static unsigned int Theta(const double x)
{
    if(x > 0) return 1;
    else      return 0;
}

int main(int argc,char *argv[])
{
double	A0;		/* lattice constant				*/
double	Delta_a;	/* Ef-Ei-Ephonon, see notes			*/
double	Delta_e;	/* Ef-Ei+Ephonon, see notes			*/
double  dki;            /* step length for loop over ki                 */
double	dKz;		/* step length for integration over Kz	*/
double	Ephonon;	/* phonon energy				*/
double	epsilon_s;	/* low frequency dielectric constant		*/
double	epsilon_inf;	/* high frequency dielectric constant		*/
double  kimax;          /* maximum value of ki                          */
double  ki;             /* carrier momentum (wave vector actually)	*/
double	m;		/* carrier effective mass			*/
double	N0;		/* number of phonons at LO phonon energy	*/
double	omega_0;	/* angular frequency of LO phonon		*/
double	T;		/* temperature					*/
double	Upsilon_star_a;	/* scattering rate prefactor			*/
double	Upsilon_star_e;	/* scattering rate prefactor			*/
double	Waif;		/* absorption scattering rate			*/
double	Weif;		/* emission scattering rate			*/
int	iKz;		/* index over Kz				*/
int     iki;            /* index over ki                                */
int	state[2];	/* electron state index				*/
int     nki;            /* number of ki calculations                    */
int	nKz;		/* number of Kz values for lookup table	*/
char	filename[9];	/* character string for output filename		*/
char	p;		/* particle					*/
bool	ff_flag;	/* form factor flag, output to file if true	*/
FILE	*FLOa;		/* pointer to absorption output file		*/
FILE	*FLOe;		/* pointer to emission   output file		*/

/* default values */

A0=5.65*1e-10;			/* lattice constant for GaAs		*/
Ephonon=36*1e-3*e;    	/* bulk LO phonon energy, 36 meV in GaAs*/
epsilon_s=13.18*eps0;	/* low frequency dielectric constant for GaAs*/
epsilon_inf=10.89*eps0;	/* high frequency dielectric constant for GaAs*/
ff_flag=false;			/* don't output formfactors	*/
m=0.067*me;			/* GaAs electron value		*/
p='e';				/* electron			*/
T=300;				/* temperature			*/

/* default values for numerical calculations	*/

nKz=1000;
nki=1000;

/* calculate step lengths	*/

dKz=2/(A0*(float)nKz);	/* Taken range of phonon integration as 2/A0 */

while((argc>1)&&(argv[1][0]=='-'))
{
 switch(argv[1][1])
 {
  case 'A':
           A0=atof(argv[2])*1e-10;
           break;
  case 'a':
	   ff_flag=true;
	   argv--;
           argc++;
           break;
  case 'E':
	   Ephonon=atof(argv[2])*1e-3*e;
	   break;
  case 'e':
	   epsilon_s=atof(argv[2])*eps0;
	   break;
  case 'f':
	   epsilon_inf=atof(argv[2])*eps0;
	   break;
  case 'm':
	   m=atof(argv[2])*me;
	   break;
  case 'p':
	   p=*argv[2];
	   switch(p)
	   {
	    case 'e': break;
	    case 'h': break;
	    case 'l': break;
	    default:  printf("Usage:  srcc [-p particle (\033[1me\033[0m, h, or l)]\n");
                      exit(0);
	   }
	   break;
  case 'T':
	   T=atof(argv[2]);
	   break;
  default :
	   printf("Usage:  srelo [-A lattice constant (\033[1m5.65\033[0mA)][-a generate form factors \033[1mfalse\033[0m]\n");
	   printf("              [-e low frequency dielectric constant (\033[1m13.18\033[0mepsilon_0]\n");
	   printf("              [-f high frequency dielectric constant (\033[1m10.89\033[0mepsilon_0]\n");
	   printf("              [-m mass (\033[1m0.067\033[0mm0)][-p particle (\033[1me\033[0m, h, or l)]\n");
	   printf("              [-T temperature (\033[1m300\033[0mK)]\n");
	   exit(0);

 }
 argv++;
 argv++;
 argc--;
 argc--;
}

/* calculate often used constants	*/

omega_0=Ephonon/hBar;		/* phonon angular frequency	*/
N0=1/(exp(Ephonon/(kB*T))-1);	/* Bose-Einstein factor	*/

Upsilon_star_a=pi*e*e*omega_0/epsilon_s*(epsilon_s/epsilon_inf-1)*(N0)
    *2*m/gsl_pow_2(hBar)*2/(8*pi*pi*pi);

Upsilon_star_e=pi*e*e*omega_0/epsilon_s*(epsilon_s/epsilon_inf-1)*(N0+1)
    *2*m/gsl_pow_2(hBar)*2/(8*pi*pi*pi);

    std::ostringstream E_filename; // Energy filename string
    E_filename << "E" << p << ".r";
    std::ostringstream wf_prefix;  // Wavefunction filename prefix
    wf_prefix << "wf_" << p;

    // Read data for all subbands from file
    std::vector<Subband> subbands = Subband::read_from_file(E_filename.str(),
            wf_prefix.str(),
            ".r",
            m);

// Read and set carrier distributions within each subband
std::valarray<double>       Ef;      // Fermi energies [J]
std::valarray<double>       N;       // Subband populations [m^{-2}]
std::valarray<unsigned int> indices; // Subband indices (garbage)
read_table("Ef.r", indices, Ef);
Ef *= e/1000.0; // Rescale to J
read_table("N.r", N);	// read populations

for(unsigned int isb = 0; isb < subbands.size(); ++isb)
subbands[isb].set_distribution(Ef[isb], N[isb]);

std::valarray<double> z;
std::valarray<double> V;
read_table("v.r", z, V);
const double Vmax = V.max();

// Read list of wanted transitions
std::valarray<unsigned int> i_indices;
std::valarray<unsigned int> f_indices;

read_table("rrp.r", i_indices, f_indices);

// Loop over all desired transitions
for(unsigned int itx = 0; itx < i_indices.size(); ++itx)
{
    // State indices for this transition (NB., these are indexed from 1)
    unsigned int i = i_indices[itx];
    unsigned int f = f_indices[itx];

    // Convenience labels for each subband (NB., these are indexed from 0)
    const Subband isb = subbands[i-1];
    const Subband fsb = subbands[f-1];

    // Subband minima
    const double Ei = isb.get_E();
    const double Ef = fsb.get_E();
    std::valarray<double> Kz(nKz);
    std::valarray<double> Gifsqr(nKz);

    ff_table(dKz,isb,fsb,nKz,Kz,Gifsqr); // generates formfactor table

    // Output form-factors if desired
    if(ff_flag)
        ff_output(Kz, Gifsqr, i,f);

 /* Generate filename for particular mechanism and open file	*/

 sprintf(filename,"LOa%i%i.r",state[0],state[1]);	/* absorption	*/
 FLOa=fopen(filename,"w");			
 sprintf(filename,"LOe%i%i.r",state[0],state[1]);	/* emission	*/
 FLOe=fopen(filename,"w");			

 /* calculate Delta variables, constant for each mechanism	*/
 Delta_a = Ef - Ei - Ephonon;
 Delta_e = Ef - Ei + Ephonon;

 /* calculate maximum value of ki and hence ki step length */
 kimax=sqrt(2*m* (Vmax - Ei))/hBar; /* sqr(hBar*kimax)/2m=Vmax-Ei */
 dki=kimax/((float)nki);

 for(iki=0;iki<nki;iki++)       /* calculate e-LO rate for all ki	*/
 {
  ki=dki*(float)iki;
  Waif=0;                       /* Initialize for integration   */
  Weif=0;                       /* Initialize for integration   */

  /* Integral over phonon wavevector Kz	*/

  for(iKz=0;iKz<nKz;iKz++)
  {
   Waif += Gifsqr[iKz] /
         sqrt(gsl_pow_4(Kz[iKz]) +
              2*gsl_pow_2(Kz[iKz]) * (2*gsl_pow_2(ki)-2*m*Delta_a/gsl_pow_2(hBar))+
              gsl_pow_2(2*m*Delta_a/gsl_pow_2(hBar))
             );

   Weif+= Gifsqr[iKz] /
         sqrt(gsl_pow_4(Kz[iKz])+
              2*gsl_pow_2(Kz[iKz])*(2*gsl_pow_2(ki)-2*m*Delta_e/gsl_pow_2(hBar))+
              gsl_pow_2(2*m*Delta_e/gsl_pow_2(hBar))
             );



  } /* end integral over Kz	*/

  Waif*=Upsilon_star_a*pi*dKz;	/* Note integral from 0->inf, hence *2	*/
  Weif*=Upsilon_star_e*pi*dKz;	/* Note integral from 0->inf, hence *2	*/

  /* Now check for energy conservation!, would be faster with a nasty `if'
     statement just after the beginning of the ki loop!			*/

  Weif*=Theta(gsl_pow_2(hBar*ki)/(2*m)-Delta_e)*
        Theta(Vmax-Ei+Ephonon-gsl_pow_2(hBar*ki)/(2*m));

  Waif*=Theta(gsl_pow_2(hBar*ki)/(2*m)-Delta_a)*
        Theta(Vmax-Ei-Ephonon-gsl_pow_2(hBar*ki)/(2*m));

  /* output scattering rate versus carrier energy=subband minima+in-plane
     kinetic energy						*/

  fprintf(FLOa,"%20.17le %20.17le\n",(Ei + gsl_pow_2(hBar*ki)/(2*m))/
                                    (1e-3*e),Waif);

  fprintf(FLOe,"%20.17le %20.17le\n",(Ei + gsl_pow_2(hBar*ki)/(2*m))/
                                    (1e-3*e),Weif);

 }
 fclose(FLOa);	/* close output file for this mechanism	*/
 fclose(FLOe);	/* close output file for this mechanism	*/
} /* end while over states */

return EXIT_SUCCESS;
} /* end main */

/* This function outputs the formfactors into files	*/
static void ff_table(const double   dKz,
                     const Subband &isb,
                     const Subband &fsb,
                     unsigned int   nKz,
                     std::valarray<double> &Kz,
                     std::valarray<double> &Gifsqr)
{
    for(unsigned int iKz=0;iKz<nKz;iKz++)
    {
        Kz[iKz]     = iKz*dKz;               // Magnitude of phonon wave vector
        Gifsqr[iKz] = Gsqr(Kz[iKz], isb, fsb); // Squared form-factor
    }
}

/* This function calculates the overlap integral squared between the two
   states	*/
static double Gsqr(const double   Kz,
                   const Subband &isb,
                   const Subband &fsb)
{
 const std::valarray<double> z = isb.z_array();
 const double dz = z[1] - z[0];
 const double nz = z.size();
 const std::valarray<double> psi_i = isb.psi_array();
 const std::valarray<double> psi_f = fsb.psi_array();

 std::complex<double> I(0,1); // Imaginary unit

 // Integral of i(=0) and f(=2) over z
 std::valarray< std::complex<double> > G_integrand_dz;
 for(unsigned int iz=0; iz<nz; ++iz)
 {
  G_integrand_dz[iz] = exp(Kz*z[iz]*I) * psi_i[iz] * psi_f[iz];
 }

 std::complex<double> G = integral(G_integrand_dz, dz);

 return norm(G);	/* cmod---modulus of complex number	*/
}

/* This function outputs the formfactors into files	*/
void ff_output(const std::valarray<double> &Kz,
               const std::valarray<double> &Gifsqr,
               unsigned int        i,
               unsigned int        f)
{
 char	filename[9];	/* output filename				*/
 sprintf(filename,"G%i%i.r",i,f);	
 write_table(filename, Kz, Gifsqr);
}
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :
