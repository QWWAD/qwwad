/**
 * \file   ivdb.cpp
 * \brief  Calculate I-V for double barrier structure
 * \author Paul Harrison  <p.harrison@shu.ac.uk>
 * \author Alex Valavanis <a.valavanis@leeds.ac.uk>
 *
 * \details This code reads in the transmission coefficient T(E) for a 
 *          double barrier and uses a very simple model to calculate 
 *          an I-V curve.
 */

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <valarray>
#include <gsl/gsl_math.h>
#include "dos-functions.h"
#include "qclsim-constants.h"
#include "qclsim-fermi.h"
#include "qclsim-fileio.h"

using namespace Leeds;
using namespace constants;

int main(int argc,char *argv[])
{
    double L1;             /* left barrier width                   */
    double L2;             /* central well width                   */
    double L3;		/* right hand barrier width		*/
    double m;		/* effective mass			*/
    double T;		/* temperature				*/

    /* default values, appropriate to GaAs-GaAlAs */
    const double Ef=2*1e-3*e; // just set Ef=2meV to represent some fixed density
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

    const double dE = E[1] - E[0]; // Calculate energy step length, assume constant
    const size_t nF = 100; // Number of field points
    std::valarray<double> V(nF); // Voltage drop across structure [V]
    std::valarray<double> current(nF);

    // Loop over field
    for(unsigned int iF=0; iF<nF; ++iF)
    {
        const double F = iF*1e+5;            // Electric field [Vm^{-1}]
        const double DeltaE=e*F*(L1+0.5*L2); // Field induced shift in band energy
        V[iF] = F*(L1+L2+L3);

        current[iF] = 0; // Initialise current for integration
        for(unsigned int iE=0;iE<n;iE++)
        {
            if(E[iE] > DeltaE)	/* only add contribution to current if E>DeltaE	*/
            {
                const double rho   = calculate_dos_3D(m, E[iE] - DeltaE); // Bulk density of states
                const double _f_FD = f_FD(Ef + DeltaE, E[iE], T); // Fermi function

                current[iF] += Tx[iE]*_f_FD*rho*dE;
                /* if(iF==0)printf("%20.17le %20.17le\n",E,rho*f_FD); output f(E)rho(E) */
            }
        }
    }

    write_table_xy("IV.r", V, current);

    return EXIT_SUCCESS;
}
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :
