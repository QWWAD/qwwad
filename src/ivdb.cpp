/*=========================================================
  ivdb I-V for Double Barrier
  =========================================================*/

/* This code reads in the transmission coefficient T(E) for a 
   double barrier and uses a very simple model to calculate 
   an I-V curve.

   Paul Harrison, May 1998		*/

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <valarray>
#include <gsl/gsl_math.h>
#include "qclsim-constants.h"
#include "qclsim-fileio.h"

using namespace Leeds;
using namespace constants;

int main(int argc,char *argv[])
{
    double	DeltaE;		/* change in energy of band minima	*/
    double	dE;		/* energy integration step length	*/
    double	Ef;		/* Fermi energy				*/
    double  F;              /* Electric field                       */
    double	f_FD;		/* Fermi Dirac distribution function	*/
    double	current;	/* current				*/
    double  L1;             /* left barrier width                   */
    double  L2;             /* central well width                   */
    double	L3;		/* right hand barrier width		*/
    double	m;		/* effective mass			*/
    double	rho;		/* bulk density of states		*/
    double	T;		/* temperature				*/
    double	V;		/* barrier height			*/

    FILE	*FIV;		/* pointer to output file `IV.r'	*/


    /* default values, appropriate to GaAs-GaAlAs */
    Ef=2*1e-3*e;		/* just set Ef=2meV to represent some fixed density */
    F=0.0;
    L1=100e-10;
    L2=100e-10;
    L3=100e-10;
    m=0.067*me;
    T=300;

    /* Computational values	*/


    while((argc>1)&&(argv[1][0]=='-'))
    {
        switch(argv[1][1])
        {
            case 'a':
                L1=atof(argv[2])*1e-10;
                break;
            case 'b':
                L2=atof(argv[2])*1e-10;
                break;
            case 'c':
                L3=atof(argv[2])*1e-10;
                break;
            case 'm':
                m=atof(argv[2])*me;
                break;
            case 'T':
                T=atof(argv[2]);
                break;
            default:
                printf("Usage:  ivdb [-a left hand barrier width (\033[1m100\033[0mA)][-b well width (\033[1m100\033[0mA)]\n");
                printf("             [-c right hand barrier width (\033[1m100\033[0mA)]\n");
                printf("             [-m effective mass (\033[1m0.067\033[0mm0)]\n");
                printf("             [-T temperature (\033[1m300\033[0mK)]\n");
                exit(0);
        }
        argv++;
        argv++;
        argc--;
        argc--;
    }

    std::valarray<double> Tx; // Transmission coefficient
    std::valarray<double> E;  // Energy [J]

    read_table_xy("T.r", E, Tx);
    E *= 1e-3*e; // Rescale energy to J
    const size_t n = E.size();

    FIV=fopen("IV.r","w");		/* Open output file for writing	*/

    dE=E[1] - E[0];	/* Calculate step length, assume constant */

    for(unsigned int iF=0; iF<100; ++iF)		/* Loop for different fields	*/
    {
        F=(double)iF*1e+5;		/* Convert kVcm^-1-->Vm^-1	*/
        DeltaE=e*F*(L1+0.5*L2);	/* Calculate DeltaE	*/
        V=F*(L1+L2+L3);		/* Calculate voltage	*/

        current=0;			/* Initialise current for integration	*/
        for(unsigned int iE=0;iE<n;iE++)
        {
            if(E[iE] > DeltaE)	/* only add contribution to current if E>DeltaE	*/
            {
                rho=gsl_pow_3(sqrt(2*m)/hBar)*sqrt(E[iE]-DeltaE)/(2*gsl_pow_2(pi));
                f_FD=1/(exp((E[iE]-(Ef+DeltaE))/(kB*T))+1);

                current+=Tx[iE]*f_FD*rho*dE;
                /* if(iF==0)printf("%20.17le %20.17le\n",E,rho*f_FD); output f(E)rho(E) */
            }
        }

        fprintf(FIV,"%20.17le %20.17le\n",V,current);	/* output data */
    } /* end loop over iE	*/

    fclose(FIV);

    return EXIT_SUCCESS;
}
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :
