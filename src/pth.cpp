/*=========================================================
                           pth
  =========================================================*/

/* this programme calculates the confined energy levels of a 
   Poschl Teller potential hole and writes the potential to a 
   file (v.r) suitable for solution with the shooting method.

   Paul Harrison, May 1998				*/

#include <cstdio>
#include <cstdlib>
#include <gsl/gsl_math.h>
#include "qclsim-constants.h"

using namespace Leeds;
using namespace constants;

/**
 * Configure command-line options for the program
 */
Options configure_options(int argc, char* argv[])
{
    Options opt;

    opt.add_numeric_option("alpha,a",    0.1,   "Width parameter [1/angstrom].");
    opt.add_numeric_option("lambda,l",   2.0,   "Depth parameter.");
    opt.add_numeric_option("mass,m",     0.067, "Effective mass (relative to free electron).");
    opt.add_size_option   ("nz,N",       100,   "Number of spatial points for output file.");
    opt.add_size_option   ("nst,s",      1,     "Number of states to find.");
    opt.add_char_option   ("particle,p", 'e',   "ID of particle to be used: 'e', 'h' or 'l', for electrons, heavy holes or light holes respectively.");
    opt.add_numeric_option("vcb",        0.00,  "Band-edge potential [eV]");
    opt.add_numeric_option("alpha",      0.00,  "Non-parabolicity parameter [eV^{-1}]");
    opt.add_numeric_option("E-cutoff",          "Cut-off energy for solutions [meV]");

    std::string doc("Generate a Poeschl--Teller potential profile.");

    std::string details("The following output text files are created:\n"
                        "  'E*.r'   \tEnergy of each state:\n"
                        "           \tCOLUMN 1: state index.\n"
                        "           \tCOLUMN 2: energy [meV].\n"
                        "  'wf_*i.r'\tWave function amplitude at each position\n"
                        "           \tCOLUMN 1: position [m].\n"
                        "           \tCOLUMN 2: wave function amplitude [m^{-1/2}].\n"
                        "\n"
                        "\tIn each case, the '*' is replaced by the particle ID and the 'i' is replaced by the number of the state.\n"
                        "\n"
                        "Examples:\n"
                        "   Compute the ground state in a 150-angstrom well with effective mass = 0.1 m0:\n\n"
                        "   efiw --width 150 --mass 0.1\n"
                        "\n"
                        "   Compute the first three heavy-hole states in a 200-angstrom well, using effective mass = 0.62 m0:\n\n"
                        "   efiw --width 200 --mass 0.62 --particle h");

    opt.add_prog_specific_options_and_parse(argc, argv, doc, details);

    return opt;
};

int main(int argc,char *argv[])
{
    const Options opt = configure_options(argc, argv);

    const double alpha  = opt.get_numeric_option("alpha") * 1e10; // Width parameter [1/m]
    const double lambda = opt.get_numeric_option("lambda");       // Depth parameter

double	E;		/* energy				*/
double	L;		/* length of potential			*/
double	m;		/* effective mass			*/
double	V;		/* potential				*/
double	z;		/* displacement along growth axis	*/
int	n;		/* principal quantum number		*/
int	N;		/* number of z points per Angstrom	*/
char	p;		/* particle (e, h, or l)		*/
char    filename[9];    /* output filename                      */

FILE	*FE;		/* pointer to energy file 		*/
FILE	*FV;		/* pointer to potential file 		*/

/* default values, appropriate to GaAs-GaAlAs */
lambda=2.0;
L=300e-10;
m=0.067*me;
N=1;
p='e';

while((argc>1)&&(argv[1][0]=='-'))
{
 switch(argv[1][1])
 {
  case 'l':
	   lambda=atof(argv[2]);
	   break;
  case 'L':
	   L=atof(argv[2])*1e-10;
	   break;
  case 'm':
	   m=atof(argv[2])*me;
	   break;
  case 'N':
	   N=atoi(argv[2]);
	   break;
  case 'p':
           p=*argv[2];
           switch(p)
           {
            case 'e': break;
            case 'h': break;
            case 'l': break;
            default:  printf("Usage:  pth [-p particle (\033[1me\033[0m, h, or l)]\n");
                      exit(0);
           }
           break;
  default:
	   printf("Usage:  pth [-a width parameter alpha (\033[1m0.1\033[0m])/A\n");
           printf("            [-l depth parameter lambda \033[1m2.0\033[0m]\n");
           printf("            [-L extent of z- domain (\033[1m300\033[0mA)]\n");
	   printf("            [-m well mass (\033[1m0.067\033[0mm0)]\n");
	   printf("            [-N points/Angstrom (\033[1m1\033[0m/A)]\n");
           printf("            [-p particle (\033[1me\033[0m, h, or l)]\n");
	   exit(0);
 }
 argv++;
 argv++;
 argc--;
 argc--;
}

/* Write potential to file	*/

FV=fopen("v.r","w");

z=-L/2;		/* initialise z position	*/
do      
{
 V=-gsl_pow_2(hBar*alpha)*lambda*(lambda-1)/(2*m*gsl_pow_2(cosh(alpha*z)));
 fprintf(FV,"%20.17le %20.17le\n",z,V);
 z+=1/((float)N)*1e-10;
}while(z<=(L/2));

fclose(FV);

/* Write energy levels to file `Ee.r'	*/

sprintf(filename,"E%c.r",p);	/* define output filename	*/
FE=fopen(filename,"w");

n=0;
while((lambda-1-n)>=0)	
{
 E=-gsl_pow_2(hBar*alpha)*gsl_pow_2(lambda-1-(float)n)/(2*m);
 
 /* Write data to file, note n->n+1 to conform with other standards	*/

 fprintf(FE,"%i %20.17le\n",n+1,E/(1e-3*e));	

 n++;
}

fclose(FE);

return EXIT_SUCCESS;
}
