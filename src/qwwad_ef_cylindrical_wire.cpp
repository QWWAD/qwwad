/**
 * \file   qwwad_ef_cylindrical_wire.cpp
 * \brief  Envelope Function Circular WIRE
 * \author Paul Harrison  <p.harrison@shu.ac.uk>
 * \author Alex Valavanis <a.valavanis@leeds.ac.uk>
 *
 * \details This program uses a shooting technique to calculate the
 *          uncorrelated one particle energies of any user supplied
 *          radial potential.  The potential is read from the file v.r
 *
 *          This program has been butchered from `efshoot', there may be 
 *          a little redundant code lying around...
 */

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include "ef-helpers.h"
#include "struct.h"
#include "maths.h"
#include "qwwad/constants.h"
#include "qwwad/file-io-deprecated.h"

using namespace QWWAD;
using namespace constants;

static double psi_at_inf(const double  E,
                         const double  delta_z,
                         files        *fdata,
                         const data11 *data_m0Eg,
                         const int     n,
                         const bool    np_flag);

int main(int argc,char *argv[])
{
double d_E;		/* infinitesmal energy               */
double delta_E;		/* small but finite energy           */
double delta_z;		/* z separation of input potentials  */
double dy;		/* derivative of function            */
double E;		/* electron and hole energies        */
double E_start;		/* initial energies (if not default) */
double x;		/* independent variable (energy)     */
double y;		/* function (psi at infinity)        */
double y1;		/* temporary y value                 */
double y2;		/* temporary y value                 */
int    i_state;		/* electron/hole state index         */
int    n;		/* length of potential file		 */
int    state;		/* electron and hole output states   */
char   p;		/* particle				 */
char   filename[9];	/* output filename			 */
bool   np_flag;         /* Hamiltonian flag def.=1=>D(1/m)D  */
FILE   *FE;		/* outputfile for el. energy states  */
files  *data_start;	/* start address of potential	 */
data11 *data_m0Eg=NULL;	/* start address of m(0) and Eg		*/

/* default values */

delta_E=1e-3*e;
E_start=0e-3*e;
np_flag=false;
p='e';
state=1;

while((argc>1)&&(argv[1][0]=='-'))
{
 switch(argv[1][1])
 {
  case 'a':
	   np_flag=true;
	   argv--;
	   argc++;
	   break;
  case 'd':
	   delta_E=atof(argv[2])*1e-3*e;
	   break;
  case 'e':
	   E_start=atof(argv[2])*1e-3*e;
	   break;
  case 'p':
	   p=*argv[2];
	   switch(p)
	   {
	    case 'e': break;
	    case 'h': break;
	    case 'l': break;
	    default:  printf("Usage:  efcwire [-p particle (e, h, or l)]\n");
                      exit(0);
	   }
	   break;
  case 's':
	   state=atoi(argv[2]);
	   break;
  default :
           printf("Usage:  efcwire [-a include non-parabolicity \033[1mfalse\033[0m]\n");
	   printf("                [-d energy step (\033[1m1\033[0mmeV)][-e initial energy (\033[1m0\033[0mmeV)]\n");
	   printf("                [-p particle (\033[1me\033[0m, h, or l)]\n");
	   printf("                [-s number of states \033[1m1\033[0m]\n");
	   exit(0);

 }
 argv++;
 argv++;
 argc--;
 argc--;
}

d_E=1e-8*e;

sprintf(filename,"E%c.r",p);
FE=fopen(filename,"w");

data_start=read_data(&n);			/* reads potential file	*/

if(np_flag)
    data_m0Eg=read_Egdata(n,data_start);	/* reads bandgap data	*/	

delta_z=read_delta_z(data_start);

if(fabs(E_start)<1e-3*1e-3*e) x=V_min(data_start,n);
 else x=E_start;                       /* initial energy estimate */

for(i_state=1;i_state<=state;i_state++)  
{
 /* increment energy-search for f(x)=0 */
 y2=psi_at_inf(x,delta_z,data_start,data_m0Eg,n,np_flag);

 do
 {
  y1=y2;
  x+=delta_E;
  y2=psi_at_inf(x,delta_z,data_start,data_m0Eg,n,np_flag);
 }while(y1*y2>0);

/* improve estimate using midpoint rule */

 x-=fabs(y2)/(fabs(y1)+fabs(y2))*delta_E;

/* implement Newton-Raphson method */

 do
 {
  y=psi_at_inf(x,delta_z,data_start,data_m0Eg,n,np_flag);
  dy=(psi_at_inf(x+d_E,delta_z,data_start,data_m0Eg,n,np_flag)-
      psi_at_inf(x-d_E,delta_z,data_start,data_m0Eg,n,np_flag))/
     (2.0*d_E);
  x-=y/dy;
 }while(fabs(y/dy)>1e-12*e);

 E=x;

 fprintf(FE,"%i %24.17le\n",i_state,E/(1e-3*e));

 x+=delta_E;    /* clears x from solution */

 }/* end i_state */

fclose(FE);
free(data_start);
if(np_flag)free(data_m0Eg);

return EXIT_SUCCESS;
} /* end main */

/**
 * Finds the value of the wavefunction (\f$\psi\f$)
 * at \f$+\infty\f$ for a given energy. The solution
 * to the energy occurs for \f$\psi(+\infty)=0\f$.
 *
 * \param[in]      E         energy
 * \param[in]      delta_z
 * \param[in, out] fdata
 * \param[in]      data_m0Eg
 * \param[in]      n
 * \param[in]      np_flag   Whether to use non-parabolic effective mass
 *
 * \returns The value of the wavefunction at \f$\infty\f$
 */
static double psi_at_inf(const double  E,
                         const double  dr,
                         files        *fdata,
                         const data11 *data_m0Eg,
                         const int     n,
                         const bool    np_flag)
{
    double psi[3];              /* wavefunction at z-delta_z,
                                   z and z+delta_z              */
    int	i;		     /* index			     */

    /* account for non-parabolicity with mstar(E)	*/
    if(np_flag)
    {
        for(i=0;i<n;i++)
        {
            // Find nonparabolicity parameter using Eq. 3.77, QWWAD3
            const double alpha= pow(1-data_m0Eg[i].a/me, 2) / data_m0Eg[i].b;

            /* Find effective mass at the desired energy using Eq. 3.76, QWWAD3 */
            fdata[i].mstar=data_m0Eg[i].a*(1.0+alpha*(E-fdata[i].V));
        }
    }

    /* boundary conditions: Eq. 8.55, QWWAD3 */
    psi[0]=1.0;
    psi[1]=1.0;
    psi[2]=0;

    for(i=1;i<(n-1);i++)              /* last potential not used */
    {
        const double r = fdata[i].z;
        const double m = fdata[i].mstar;
        const double V = fdata[i].V;

        /* Find wavefunction at next point, using Eq. 8.54, QWWAD3 */
        psi[2]=(
                2.0*r*(2.0*m*dr*dr/(hBar*hBar)*(V-E)+2.0)*
                psi[1]+
                (-2.0*r+dr)*psi[0]
               )
            /(2.0*r+dr);

/*        printf("%e %e\n", V, psi[2]); */

        /* Shift previous values along ready for next iteration */
        psi[0]=psi[1];
        psi[1]=psi[2];
    } 

    return psi[2];
}
