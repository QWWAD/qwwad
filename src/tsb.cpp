/**
 * \file   tsb.cpp
 * \brief  Calculate transmission coefficient for single barrier
 * \author Paul Harrison  <p.harrison@shu.ac.uk>
 * \author Alex Valavanis <a.valavanis@leeds.ac.uk>
 */

#include <cstdio>
#include <cstdlib>
#include <strings.h>
#include <cmath>
#include <gsl/gsl_math.h>
#include "struct.h"
#include "maths.h"
#include "qclsim-constants.h"

using namespace Leeds;
using namespace constants;

int main(int argc,char *argv[])
{
    double	dE;		/* energy step				*/
    double	E;		/* energy				*/
    double	k;		/* wave vector in `well' material	*/
    double	kdash;		/* wave vector in barrier material	*/
    double	K;		/* decay constant in barrier material	*/
    double	L;		/* well width				*/
    double	m;		/* effective mass 			*/
    double	T;		/* transmission coefficient		*/
    double	V;		/* barrier height			*/
    FILE	*FT;		/* pointer to output file `T.r'		*/


    /* default values, appropriate to GaAs-GaAlAs */

    L=100e-10;          
    m=0.067*me;
    V=100*1e-3*e;   

    /* Computational values	*/

    dE=1e-4*e;        /* arbitrarily small energy---0.1meV   */


    while((argc>1)&&(argv[1][0]=='-'))
    {
        switch(argv[1][1])
        {
            case 'd':
                dE=atof(argv[2])*1e-3*e;
                break;
            case 'L':
                L=atof(argv[2])*1e-10;
                break;
            case 'm':
                m=atof(argv[2])*me;
                break;
            case 'V':
                V=atof(argv[2])*1e-3*e;
                break;
            default:
                printf("Usage:  tsb [-d energy step (\033[1m0.1\033[0mmeV)][-L barrier width (\033[1m100\033[0mA)]\n");
                printf("            [-m effective mass (\033[1m0.067\033[0mm0)][-V barrier height (\033[1m100\033[0mmeV)]\n");
                exit(0);
        }
        argv++;
        argv++;
        argc--;
        argc--;
    }

    /* Initialise energy	*/

    E=0; 

    /* Open output file	*/

    FT=fopen("T.r","w");

    do      /* loop increments energy */
    {
        k=sqrt(2*m*E)/hBar;
        K=sqrt(2*m*(V-E))/hBar;
        T=1/(1+gsl_pow_2((k*k+K*K)/(2*k*K))*gsl_pow_2(sinh(K*L)));
        fprintf(FT,"%20.17le %20.17le\n",E/(1e-3*e),T);
        E+=dE;
    }while(E<V);

    do      /* loop increments energy */
    {
        k=sqrt(2*m*E)/hBar;
        kdash=sqrt(2*m*(E-V))/hBar;
        T=1/(1+gsl_pow_2((k*k-kdash*kdash)/(2*k*kdash))*gsl_pow_2(sin(kdash*L)));
        fprintf(FT,"%20.17le %20.17le\n",E/(1e-3*e),T);
        E+=dE;
    }while(E<10*V);

    fclose(FT);

    return EXIT_SUCCESS;
}        /* end main */
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :
