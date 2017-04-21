/**
 * \file   qwwad_ef_cylindrical_wire_wf.cpp
 * \brief  Envelope Function Circular Wire Wave Functions
 * \author Paul Harrison <p.harrison@shu.ac.uk>
 * \author Alex Valavanis <a.valavanis@leeds.ac.uk>
 *
 * \details This program uses a shooting technique to calculate the
 *          uncorrelated one particle wavefunctions of any user supplied
 *          radial potential.  The potential is read from the file v.r
 *          The eigenenergies have been calculated previously with
 *          `qwwad_ef_cylindrical_wire'
 *          and are stored in the file E?.r
 */

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include "ef-helpers.h"
#include "qwwad/constants.h"
#include "maths.h"

using namespace QWWAD;
using namespace constants;

static double wf(const double  E,
                 const double  delta_z,
                 files        *fdata,
                 const data11 *data_m0Eg,
                 data11       *data_zwf,
                 const int     n,
                 const bool    np_flag);

int main(int argc,char *argv[])
{
double	delta_z;	/* z separation of input potentials	*/
double	E;		/* electron and hole energies		*/
double	N;		/* normalization integral		*/
int	i_state;	/* electron/hole state index		*/
int	i;		/* index				*/
int	n;		/* length of potential file		*/
char	p;		/* particle				*/
char	filename[9];	/* output filename			*/
bool	np_flag;	/* non-parabolicity flag		*/
FILE	*FE;		/* outputfile for el. energy states	*/
FILE	*Fwf;		/* outputfile for wavefunctions		*/
files	*data_start;	/* start address of potential		*/
data11	*data_m0Eg = NULL;	/* start address of m(0) and Eg		*/
data11	*data_zwf;	/* start address of z and wf		*/

/* default values */

np_flag=false;
p='e';

while((argc>1)&&(argv[1][0]=='-'))
{
 switch(argv[1][1])
 {
  case 'a':
	   np_flag=true;
	   argv--;
	   argc++;
	   break;
  case 'p':
	   p=*argv[2];
	   switch(p)
	   {
	    case 'e': break;
	    case 'h': break;
	    case 'l': break;
	    default:  printf("Usage:  efcwwf [-p particle (\033[1me\033[0m, h, or l)]\n");
                      exit(0);
	   }
	   break;
  default :
	   printf("Usage:  efcwwf [-a include non-parabolicity \033[1mfalse\033[0m]\n");
	   printf("               [-p particle (\033[1me\033[0m, h, or l)]\n");
	   exit(0);
 }
 argv++;
 argv++;
 argc--;
 argc--;
}


data_start=read_data(&n);			/* reads potential file	*/
if(np_flag)data_m0Eg=read_Egdata(n,data_start);	/* reads bandgap data	*/

delta_z=read_delta_z(data_start);

/* allocate memory for wavefunction */

data_zwf=(data11 *)calloc(n,sizeof(data11));
if(data_zwf==0){fprintf(stderr,"Cannot allocate memory!\n");exit(0);}

/* copy z data into wavefunction structure	*/

for(i=0;i<n;i++)
 (data_zwf+i)->a=(data_start+i)->z;

/* open appropriate energy file and read in energies sequentially	*/

sprintf(filename,"E%c.r",p);
FE=fopen(filename,"r");

while(fscanf(FE,"%i %le",&i_state,&E)!=EOF)
{
 E*=1e-3*e;	/* convert E from meV->J	*/
 N=wf(E,delta_z,data_start,data_m0Eg,data_zwf,n,np_flag);

 sprintf(filename,"wf_%c%i.r",p,i_state);
 Fwf=fopen(filename,"w");

 for(i=0;i<n;i++)
  fprintf(Fwf,"%20.17le %20.17le\n",(data_zwf+i)->a,((data_zwf+i)->b)/sqrt(N));

 fclose(Fwf);

}/* end while */

fclose(FE);
free(data_start);
if(np_flag)free(data_m0Eg);
free(data_zwf);

return EXIT_SUCCESS;
} /* end main */

/**
 * Finds the wavefunction (\f$\psi\f$) for a given energy
 *
 * \param[in]      E         energy
 * \param[in]      delta_z
 * \param[in, out] fdata
 * \param[in]      data_m0Eg
 * \param[out]     data_zwf
 * \param[in]      n
 * \param[in]      np_flag   Whether to use non-parabolic effective mass
 *
 * \returns The value of the wavefunction at \f$\infty\f$
 */
static double wf(const double  E,
                 const double  delta_z,
                 files        *fdata,
                 const data11 *data_m0Eg,
                 data11       *data_zwf,
                 const int     n,
                 const bool    np_flag)
{
 double N=0;		     /* normalization integral       */
 double psi[3];              /* wavefunction at z-delta_z,
                                z and z+delta_z              */
 int	i;		     /* index			     */

 /* account for non-parabolicity with mstar(E)	*/

 if(np_flag)
 {
  for(i=0;i<n;i++)
  {
   /* Find nonparabolicity parameter using Eq. 3.77, QWWAD3 */
   const double alpha=pow(1-data_m0Eg[i].a/me,2)/data_m0Eg[i].b;

   /* Find effective mass at the desired energy using Eq. 3.76, QWWAD3 */
   fdata[i].mstar=data_m0Eg[i].a*(1.0+alpha*(E-fdata[i].V));
  }
 }

 /* boundary conditions: Eq. 8.55. QWWAD3 */
 psi[0]=1.0;
 psi[1]=1.0;

 data_zwf[0].b=psi[0];
 data_zwf[1].b=psi[1];

 N+=pow(psi[0],2);
 N+=pow(psi[1],2);

 fdata++;                    /* ignore data corresponding to psi[0] */

  for(i=1;i<(n-1);i++)              /* last potential not used */
  {
   /* Find wavefunction at next point, using Eq. 8.54, QWWAD3 */
   psi[2]=(
           2*(fdata->z)*(2*(fdata->mstar)*pow(delta_z/hBar,2)*(fdata->V-E)+2)*
	   psi[1]+
           (-2*(fdata->z)+delta_z)*psi[0]
          )
           /(2*(fdata->z)+delta_z);
   data_zwf[i+1].b=psi[2];
   N+=pow(psi[2],2);
   psi[0]=psi[1];
   psi[1]=psi[2];
   fdata++;
  } 

 return N*delta_z;
}
