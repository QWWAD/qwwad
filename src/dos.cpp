/**
 * \file   dos.c Density of states calculator
 *
 * \author Paul Harrison  <p.harrison@shu.ac.uk>
 * \author Alex Valavanis <a.valavanis@leeds.ac.uk>
 *
 * \details This program calculates the density of states for bulk (3D),
 *          quantum wells (2D) and quantum wires (1D), for a series of subband
 *          minima which are read in from the external file `Ee.r', or `Eh.r'
 */

#include "qclsim-constants.h"
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <gsl/gsl_math.h>

using namespace Leeds;
using namespace constants;

double * read_E(char  p,
        int  *nE);

int main(int argc,char *argv[])
{
    double	dos_bulk;	/* bulk density of states			*/
    double	dos_2D;		/* quantum well (2D) density of states		*/
    double	dos_1D;		/* quantum wire (1D) density of states		*/
    double	energy;		/* the independent variable			*/
    double	*E;		/* pointer to subband minima			*/
    double	m;		/* effective mass				*/
    int	nE;		/* number of confined states			*/
    int	n;		/* number of output energies			*/
    int	i;		/* index over subband minima			*/
    int	ie;		/* index over energy				*/
    char	p;		/* electron, light- or heavy-hole		*/
    FILE	*Frho;		/* pointer to output file rho.r			*/

    /* default values */

    m=0.067*me;		/* GaAs electron value		*/
    p='e';			/* electron			*/

    /* default values for numerical calculations	*/

    n=1000;

    while((argc>1)&&(argv[1][0]=='-'))
    {
        switch(argv[1][1])
        {
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
                    default:  printf("Usage:  dos [-p particle (\033[1me\033[0m, h, or l)]\n");
                              exit(0);
                }
                break;
            default :
                printf("Usage:  dos [-m mass (\033[1m0.067\033[0mm0)][-p particle (\033[1me\033[0m, h, or l)]\n");
                exit(0);

        }
        argv++;
        argv++;
        argc--;
        argc--;
    }

    E=read_E(p,&nE);	/* read in subband minima	*/

    Frho=fopen("rho.r","w");

    for(ie=0;ie<=n;ie++)
    {
        energy=(float)ie*1e-3*e;		/* convert meV-> J	*/
        dos_bulk=gsl_pow_3(sqrt(2*m)/hBar)*sqrt(energy)/(2*gsl_pow_2(pi));

        dos_2D=0;		/* initialise before sum over subbands	*/
        dos_1D=0;		/* initialise before sum over subbands	*/
        for(i=0;i<nE;i++)
        {
            if(energy > E[i])
                dos_2D+=m/(pi*gsl_pow_2(hBar));

            if(energy>*(E+i))dos_1D+=sqrt(2*m)/hBar/(pi*sqrt(energy-*(E+i)));
        }
        fprintf(Frho,"%le %le %le %le\n",energy/(1e-3*e),dos_bulk,dos_2D,dos_1D);
    }

    fclose(Frho);
    return EXIT_SUCCESS;
} /* end main */

/**
 * Reads subband minima into memory and returns the start
 * address of this block of memory and the number of lines
 */
double * read_E(char  p,
        int  *nE)
{
    double	*E;
    int	i=0;		/* index over the energies			*/
    char	filename[9];	/* filename string				*/
    FILE 	*FE;		/* file pointer to energy data 			*/

    sprintf(filename,"E%c.r",p);
    if((FE=fopen(filename,"r"))==0)
    {
        fprintf(stderr,"Error: Cannot open input file '%s'!\n",filename);
        exit(0);
    }

    *nE=0;
    while(fscanf(FE,"%*i %*e")!=EOF)
        (*nE)++;
    rewind(FE);

    E=(double *)calloc(*nE,sizeof(double));
    if (E==0)  {
        fprintf(stderr,"Cannot allocate memory!\n");
        exit(0);
    }

    while(fscanf(FE,"%*i %le",E+i)!=EOF)
    {
        *(E+i)*=1e-3*e;		/*convert meV->J		*/
        i++;
    }

    fclose(FE);

    return(E);

}
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :
