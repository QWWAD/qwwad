/**
 * \file  sbp.cpp
 * \brief Calculate subband populations
 *
 * \details This program generates the Fermi-Dirac distribution function
 *          for an individual subband given its population and the lattice
 *          temperature.
 *
 *          Input files:
 *            Ep.r  Subband energies file, p=e,h,l
 *
 *          Output files:
 *            FDX.r F-D distribution for subband X
 *
 * \author Paul Harrison  <p.harrison@shu.ac.uk>
 * \author Alex Valavanis <a.valavanis@leeds.ac.uk>
 */

#include <cstdio>
#include <cstdlib>
#include <gsl/gsl_math.h>
#include "qclsim-constants.h"
#include "qwwad-fileio.h"
#include "qwwad-options.h"

using namespace Leeds;
using namespace constants;

typedef
struct	{
 double	E;		    /* total electron energy		  */
 double	SR;		    /* scattering rate            	  */
} data;

double   f(double E_F, double Emax, double Emin, double m, int N, double T);
double * read_populations(int n);
void     calc_dist(double Emin, double Ef, double m, double T, int nE, int s);
double   calc_fermilevel(double E, double m, double N, double T);
double   Vmax();

/**
 * Handler for command-line options
 */
class SBPOptions : public Options
{
    public:
        SBPOptions(int argc, char* argv[])
        {
            try
            {
                program_specific_options->add_options()
                    ("fd,f", po::bool_switch()->default_value(false),
                     "Output Fermi-Dirac distribution.")

                    ("mass,m", po::value<double>()->default_value(0.067),
                     "Effective mass (relative to that of a free electron)")

                    ("nenergy,n", po::value<size_t>()->default_value(1000),
                     "Number of energy samples for output file")

                    ("particle,p", po::value<char>()->default_value('e'),
                     "Particle to be used: 'e', 'h' or 'l'")

                    ("temperature,T", po::value<double>()->default_value(300),
                     "Temperature of carrier distribution [K]")
                    ;

                std::string doc("Find the Fermi-Dirac distribution function for an individual "
                                "subband given its population and temperature."
                                "The subband minimum is read from the file \"E*.r\", "
                                "and the distribution is written to \"FDi.r\" where the '*' "
                                "is replaced by the particle ID in each case and the "
                                "'i' is replaced by the number of the state");

                add_prog_specific_options_and_parse(argc, argv, doc);	
            }
            catch(std::exception &e)
            {
                std::cerr << e.what() << std::endl;
                exit(EXIT_FAILURE);
            }
        }

        /// \returns the effective mass [kg]
        bool get_fd_flag() const {return vm["fd"].as<bool>();}

        /// \returns the effective mass [kg]
        double get_mass() const {return vm["mass"].as<double>()*me;}

        /// \returns the particle ID
        char get_particle() const {return vm["particle"].as<char>();}

        /// \returns the number of energy samples
        size_t get_n_energy() const {return vm["nenergy"].as<size_t>();}

        /// \returns the temperature of the carrier distributeion [K]
        double get_temperature() const {return vm["temperature"].as<double>();}
};

int main(int argc,char *argv[])
{
    SBPOptions opt(argc, argv);

FILE	*FEf;			/* file pointer to Fermi Energy file	*/

const bool   FD_flag = opt.get_fd_flag();
const double m       = opt.get_mass();
const size_t nE      = opt.get_n_energy();
const char   p       = opt.get_particle();
const double T       = opt.get_temperature();

std::valarray<double> E=read_E(p); // reads subband energy file
const size_t n = E.size();
double *N=read_populations(n);		/* reads subband populations file	*/

if((FEf=fopen("Ef.r","w"))==0)
 {fprintf(stderr,"Error: Cannot open input file 'Ef.r'!\n");exit(0);}

for(unsigned int s=0; s<n; s++) // s=0 => ground state
{
 const double Ef=calc_fermilevel(E[s],m,*(N+s),T);

 fprintf(FEf,"%i %20.17le\n",s+1,Ef/(1e-3*e));

 if(FD_flag) calc_dist(E[s],Ef,m,T,nE,s);
}

fclose(FEf);

return EXIT_SUCCESS;
} /* end main */

/**
 * \brief calculates the probability of occupation of the subband energies
 *
 * \param[in] Emin subband minima
 * \param[in] Ef   Fermi energy
 * \param[in] m    effective mass
 * \param[in] T    temperature
 * \param[in] nE   number of energies to output FD
 * \param[in] s    number of subband
 */
void calc_dist(double Emin, double Ef, double m, double T, int nE, int s)
{
double	dE;			/* energy increment in integration	*/
double	E;			/* Energy				*/
double	Emax;
double	f;			/* probability of occupation		*/
double	Ne=0;			/* total electron density (for checking)*/
double	vmax;			/* maximum value of potential		*/
int     i;			/* general index			*/
char	filename[9];		/* output filename for FD distribs.	*/
FILE	*FFD;			/* file pointer to Fermi-Dirac probability
				   of occupation of subband file	*/

sprintf(filename,"FD%i.r",s+1);
FFD=fopen(filename,"w");

vmax=Vmax();
Emax=Ef+10*kB*T;if(Emax>vmax) Emax=vmax;

/* Implement trapezium rule integration, i.e. `ends+2middles' 	*/

f=1/(exp((Emin-Ef)/(kB*T))+1);
fprintf(FFD,"%20.17le %20.17le\n",Emin/(1e-3*e),f);
Ne+=f;

dE=(Emax-Emin)/(nE-1);
E=Emin;
for(i=1;i<nE-1;i++)
{
 E+=dE;
 f=1/(exp((E-Ef)/(kB*T))+1);
 fprintf(FFD,"%20.17le %20.17le\n",E/(1e-3*e),f);
 Ne+=2*f;
}

f=1/(exp((Emax-Ef)/(kB*T))+1);
fprintf(FFD,"%20.17le %20.17le\n",Emax/(1e-3*e),f);
Ne+=f;

Ne*=m/(pi*gsl_pow_2(hBar))*0.5*dE;

printf("Ne=%20.17le\n",Ne/1e+14);

fclose(FFD);
}

/**
 * \brief calculates Fermi level
 *
 * \param[in] E subband minima
 * \param[in] m effective mass
 * \param[in] N number of electrons/unit area
 * \param[in] T temperature
 *
 * \details Solutions are sought, using a stepping algorithm, a midpoint rule and
 *          finally a Newton-Raphson method, to the equation f(E_F)=0.  This
 *          function has been derived by integrating the total density of states
 *          multiplied by the Fermi-Dirac distribution function across the in-plane
 *          energy---thus giving the total electron density, i.e.
 *
 *                 oo
 *           Ne=  I  f(E)N(E) dE
 *                 Eo
 *
 *         where f(E) is the normal Fermi-Dirac distribution function and
 *         N(E)=m/(pi*sqr(hbar)) is the areal density of states in a QW, see
 *         Bastard p12.
 */
double calc_fermilevel(double E, double m, double N, double T)
{
double	delta_E=0.001*1e-3*e;	/* energy increment			*/
double 	Emin;			/* subband minimum			*/
double 	Emax;			/* subband maximum			*/
double	vmax;			/* potential maximum, i.e. top of well	*/
double	x;			/* independent variable			*/
double	y1;			/* dependent variable			*/
double	y2;			/* dependent variable			*/

Emin=E;				/* subband minimum 			*/

vmax=Vmax();			/* calulate potential maximum		*/

x=Emin-20*kB*T;			/* first value of x			*/

/* In this implementation, the upper limit of integration is set at the 
   Fermi level+10kT, limited at potential maximum			*/

Emax=x+10*kB*T;if(Emax>vmax) Emax=vmax;

y2=f(x,Emax,Emin,m,N,T);

do
{
 y1=y2;
 x+=delta_E;
 Emax=x+10*kB*T;if(Emax>vmax) Emax=vmax;
 y2=f(x,Emax,Emin,m,N,T);
}while(y1*y2>0);

/* improve estimate using midpoint rule */

x-=fabs(y2)/(fabs(y1)+fabs(y2))*delta_E;

return(x);
}

/**
 * Function to be solved
 *
 * \param[in] E_F   Fermi energy
 * \param[in] Emax  subband maximum (top of QW)
 * \param[in] Emin  subband minimum
 * \param[in] m     effective mass
 * \param[in] N     number of electrons/unit area
 * \param[in] T     temperature
 */
double f(double E_F, double Emax, double Emin, double m, int N, double T)
{
double	y;			/* dependent variable			*/

y=m/(pi*hBar)*(kB*T/hBar)*
  (
   ((Emax-E_F)/(kB*T)-log(1+exp((Emax-E_F)/(kB*T))))-
   ((Emin-E_F)/(kB*T)-log(1+exp((Emin-E_F)/(kB*T))))
  )
  -N;

return(y);
}

/**
 * \brief reads subband populations from N.r
 *
 * \details This function reads the populations into memory and returns the start
 *          address of this block of memory and the number of lines
 */
double * read_populations(int n)
{
 double	*N;
 int	i=0;		/* index over the energies			*/
 int	m;		/* counter over number of populations defined	*/
 FILE 	*FN;		/* file pointer to energy data 			*/

 if((FN=fopen("N.r","r"))==0)
 {
   fprintf(stderr,"Error: Cannot open input file 'N.r'!\n");
   exit(0);
 }

 m=0;
 while(fscanf(FN,"%*i %*f")!=EOF)
  m++;
 rewind(FN);

 if(m!=n)
  {printf("Subband populations not defined for all levels!\n");exit(0);}

 N=(double *)calloc(n,sizeof(double));
 if (N==0)  {
  fprintf(stderr,"Cannot allocate memory!\n");
  exit(0);
 }

 while(fscanf(FN,"%*i %lf",N+i)!=EOF)
 {
  *(N+i)*=1e+10*1e+4;	/*convert from units of 10^10cm^-2->m^-2 */
  i++;
 }

 fclose(FN);

 return(N);

}

/**
 * Scans the file v.r and returns the maximum value of the potential.
 */
double Vmax()
{
 double	max;			/* maximum value of potential energy	*/
 double	v;			/* potential				*/
 FILE	*Fv;			/* file pointer to v.r			*/

max=0;

if((Fv=fopen("v.r","r"))==0)
 {fprintf(stderr,"Error: Cannot open input file 'v.r'!\n");exit(0);}

while(fscanf(Fv,"%*e %le",&v)!=EOF)
 if(v>max) max=v;

fclose(Fv);

return(max);

}
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :
