/*===============================================================
        efxm Envelope Function concentration (X) to Mass
  ===============================================================*/

/* This program produces the function of the masses of the electrons
   and holes in the z direction and in the xy plane, the results are
   written to the files 'm.r' and 'm_perp.r' respectively. Input 
   is required from 'x.r'.
   
   Substantially revised, quadratic dependencies, multiple material 
   support, command line arguments

   Paul Harrison, September 1995                               */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "struct.h"
#include "maths.h"
#include "const.h"

int main(int argc,char *argv[])
{
double	m;		/* effective mass 			*/
double	mp;		/* effective mass perpendicular growth	*/
double	mstar;		/* constant effective mass		*/
double	x;		/* alloy concentration `x'		*/
double	y;		/* alloy concentration `y'		*/
double	z;		/* displacement along growth (z-) axis	*/
char	material[9];	/* material string			*/
char	Material;	/* material character			*/
char	p;		/* particle e, h, or l			*/
FILE	*Fx;		/* Mn concentraion file 'x.r'		*/
FILE    *Fm;		/* effective mass in z direction, 'm.r'	*/
FILE    *Fmp;		/* effective mass in the xy plane, 'm_perp.r'	*/

/* Define global defaults */

Material='a';
mstar=0;
p='e';

while((argc>1)&&(argv[1][0]=='-'))
{
 switch(argv[1][1])
 {
  case 'M':
           (void)strcpy(material,argv[2]);
           if(!strcmp(material,"gaalas"))Material='z';
           if(!strcmp(material,"cdmnte"))Material='b';
           if(!strcmp(material,"inalgaas"))Material='c';
	   switch(Material)
	   {
	    case 'b': break;
	    case 'c': break;
	    case 'z': Material='a';break;
	    default :
                      printf("The only materials defined in the database are\n");
                      printf("Ga(1-x)Al(x)As, Cd(1-x)Mn(x)Te, In(1-x-y)Al(x)Ga(y)As\n");
                       exit(0);
	   }

           break;
  case 'm':
	   mstar=atof(argv[2]);
	   break;
  case 'p':
           p=*argv[2];
           switch(p)
           {
            case 'e': break;
            case 'h': break;
            case 'l': break;
            default:  printf("Usage:  efxm [-p particle (\033[1me\033[0m, h, or l)]\n");
                      exit(0);
           }
           break;
  default :
	   printf("Usage:  efxm [-M material \033[1mGa(1-x)Al(x)As\033[0m][-p particle (\033[1me\033[0m, h, or l)]\n");
	   printf("             [-m mass (m0), if set overrides m(x,y,z) \033[1mfalse\033[0m]\n");
           exit(0);

 }
 argv++;
 argv++;
 argc--;
 argc--;
}

if((Fx=fopen("x.r","r"))==0)
 {fprintf(stderr,"Error: Cannot open input file 'x.r'!\n");exit(0);}  
   
Fm=fopen("m.r","w");
Fmp=fopen("m_perp.r","w");

/* if mstar is set to a value greater than zero, then this overrides the
   dependencies set below and creates a file containing constant mass */

if(mstar>0) 
{
 while(fscanf(Fx,"%lf %lf %lf",&z,&x,&y)!=EOF)
 {
 m=mstar*m0;
 mp=mstar*m0;
 fprintf(Fm,"%20.17le %20.17le\n",z,m);
 fprintf(Fmp,"%20.17le %20.17le\n",z,mp);
 }
}

/* if mstar left as a default of zero then invoke usual dependencies */

else 
switch(Material)
{
 case 'a':	/* Ga(1-x)Al(x)As, S. Adachi, `GaAs and related materials' */
	  while(fscanf(Fx,"%lf %lf %lf",&z,&x,&y)!=EOF)  
          {
	   switch(p)
	   {
	    case 'e':
		     m=(0.067+0.083*x)*m0;
		     mp=(0.067+0.083*x)*m0;
		     break;
	    case 'h':
		     m=(0.62+0.14*x)*m0;
		     mp=(0.62+0.14*x)*m0;
		     break;
	    case 'l':
		     printf("Data not defined for Ga(1-x)Al(x)As light-hole\n");
                     exit(0);
	   }
	   fprintf(Fm,"%20.17le %20.17le\n",z,m);
	   fprintf(Fmp,"%20.17le %20.17le\n",z,mp);
	  }
	  break;
 case 'b':	/* Cd(1-x)Mn(x)Te, Long, 23rd Phys. Semicond. p1819 */
	  while(fscanf(Fx,"%lf %lf %lf",&z,&x,&y)!=EOF)  
          {
           switch(p)
           {
            case 'e':
                     m=(0.11+0.067*x)*m0;
                     mp=(0.11+0.067*x)*m0;
                     break;
            case 'h':
                     m=(0.60+0.21*x+0.15*x*x)*m0;
                     mp=(0.60+0.21*x+0.15*x*x)*m0;
                     break;
            case 'l':
                     m=(0.18+0.14*x)*m0;
                     mp=(0.18+0.14*x)*m0;
                     break;
           }
           fprintf(Fm,"%20.17le %20.17le\n",z,m);
           fprintf(Fmp,"%20.17le %20.17le\n",z,mp);
          }
	  break;
 case 'c':	/* In(1-x-y)Al(x)Ga(y)As, Landolt&Bornstein III/22a p156 */
	  while(fscanf(Fx,"%lf %lf %lf",&z,&x,&y)!=EOF)  
          {
           switch(p)
           {
            case 'e':
                     m=(0.0427+0.0685*x)*m0;
                     mp=(0.0427+0.0685*x)*m0;
                     break;
            case 'h':
                     printf("Data not defined for In(1-x-y)Al(x)Ga(y)As heavy-hole\n");
                     exit(0);
            case 'l':
                     printf("Data not defined for In(1-x-y)Al(x)Ga(y)As light-hole\n");
                     exit(0);
           }
           fprintf(Fm,"%20.17le %20.17le\n",z,m);
           fprintf(Fmp,"%20.17le %20.17le\n",z,mp);
	  }
	  break;
}



fclose(Fx);
fclose(Fm);
fclose(Fmp);

return EXIT_SUCCESS;
} /* end main */
