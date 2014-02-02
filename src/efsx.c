/*===============================================================
       efsx Envelope Function Structure to concentration (X)
  ===============================================================*/

/* This program produces the function of alloy concentration 
   x (and y for quaternaries) which defines the semiconductor
   heterostructure. The input data comes from a file 's.r',
   the output file 'x.r' is in a format suitable for conversion into
   electron and hole potentials and effective mass data, using the 
   programs `efxv' and `efxm'. 

   The output files 'm.r' and 'm_perp.r' are required by shoot.
   
   Paul Harrison, December 1996 (Major modifications, quaternaries,
                                 command line arguments)           

   Paul Harrison, May 1998, further minor modifications        
 
   Paul Harrison, June 1998, added dopant profiles		*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "struct.h" 
#include "maths.h" 
#include "const.h"
#include "bools.h"

main(int argc,char *argv[])
{
void    Create_files();	/* create files 's.r' and 'm.r'      */
int	No_of_el();	/* number of elements per line       */
int	N;		/* number of points Angstrom         */     
int	Nel;		/* Number of elements per line in s.r */
FILE    *Fp;		/* file pointer to input data        */

/* default values */

N=1;

while((argc>1)&&(argv[1][0]=='-'))
{
 switch(argv[1][1])
 {
  case 'N':
	   N=atoi(argv[2]);
	   break;
  default :
	   printf("Usage:  efsx [-N points/Angstrom \033[1m1\033[0m]\n\n");
	   printf("Format of input file `s.r' for material A(1-x-y)B(x)C(y)D or\n");
	   printf("A(1-x)B(x)C(1-y)D(y)---check material definitions with command\n");
	   printf("\n\033[1mefxm -M help\033[0m\n\n");
	   printf("=========================begin s.r===================================\n");
	   printf("layer thickness (\033[1mA\033[0m)\t x\t y\t doping level (\033[1m10^18cm^-3\033[0m)\n");
           exit(0);
 }
 argv++;
 argv++;
 argc--;
 argc--;
}

if((Fp=fopen("s.r","r"))==0)
 {printf("Error: Cannot open input file `s.r'!\n");exit(0);}  

/* In case file v0.r exists, remove it.  Thus ensuring each time a new
   structure is designed, existing files are handled correctly	*/

unlink("v0.r");	

Nel=No_of_el(Fp);

/* Check the number of columns returned is as expected.
   The default is to expect quaternaries, but if a 4th column exists too
   then produce a dopant profile					*/

switch(Nel)
{
 case 2:
	printf("Error: Expect quaternaries, add a null third column to `s.r'!\n");
        exit(0);
 case 3:break;
 case 4:break;
 default:
	printf("Error: Check input file 's.r'!\n");exit(0);
}

Create_files(Fp,N,Nel);

fclose(Fp);
}





int
No_of_el(fp)

/* This function returns the number of elements per line of the file
   associated with fp. -1 is returned if this number is not an 
   integer, the required format of 's.r'.  On exit the filepointer is 
   set to the beginning of the file.                                  */

FILE	*fp;
{
double	dummy;         /* dummy variable	            */
int     i1;            /* ASCII value of bytes read from s.r*/
int	i2=10;         /* ASCII value of previous element   */
int     lines=0;       /* number of lines                   */
int     elements=0;    /* number of elements                */

do {
  i1 = getc(fp);
		/* ASCII value 10 => EOL, doesn't count empty lines */
  lines += (i1==10 && i2!=10) ? 1 : 0;
  i2 = i1;
} while(i1 != EOF);
rewind(fp);

while(fscanf(fp,"%le", &dummy) != EOF) 
  elements++;

rewind(fp);

if (elements%lines==0) return(elements/lines);
else return(-1);
}



void
Create_files(fp,N,elements)

/* This function creates the output file 'x.r'.  
   The filepointer is changed on exit.                     */

FILE	*fp;
int	N;
int	elements;

{
double	Nda;		/* volume concentration of donors or acceptors	*/
double  w;		/* width of the region in A			*/
double  x;		/* alloy concentration				*/
double  y;		/* alloy concentration				*/
double  z_0=0;		/* z in multiples of 1./N			*/
double	z_1=0;		/* actual z coordinate				*/
int	n;		/* counter					*/
FILE    *Fd;		/* file pointer to d versus z data		*/
FILE    *Fx;		/* file pointer to x versus z data		*/

Fx=fopen("x.r","w");

if(elements==3)
{
 while(fscanf(fp,"%lf %lf %lf",&w,&x,&y)!=EOF)
 {
  n=0;
  z_1+=w*1e-10;

  while(Nint((z_0+n*1e-10/(float)N)*1e14)<Nint(z_1*1e14))
  {
   fprintf(Fx,"%20.17le %le %le\n",z_0+n*1e-10/(float)N,x,y);
   n++;
  }
 z_0+=n*1e-10/(float)N;
 }
}	/* end if	*/
else	/* If number of elements is 4 then output dopant data too 	*/
{
 Fd=fopen("d.r","w");

 while(fscanf(fp,"%lf %lf %lf %lf",&w,&x,&y,&Nda)!=EOF)
 {
  n=0;
  z_1+=w*1e-10;
  Nda*=1e+18*1e+6;		/* convert into cm^-3 and then into m^-3 */

  while(Nint((z_0+n*1e-10/(float)N)*1e14)<Nint(z_1*1e14))
  {
   fprintf(Fx,"%20.17le %le %le\n",z_0+n*1e-10/(float)N,x,y);
   fprintf(Fd,"%20.17le %20.17le\n",z_0+n*1e-10/(float)N,Nda);
   n++;
  }
 z_0+=n*1e-10/(float)N;
 }

 fclose(Fd);
}	/* end else	*/

fclose(Fx);
}
