/*=========================================================
             efxv Envelope Function x to v
  =========================================================*/

/* 
   This program converts the structure as defined in terms
   of alloy components into a potential profile for either 
   electron, light- or heavy-hole.  It has support for multiple 
   material systems, ternaries as well as quaternaries.

   In addition generation of the bandgap allows for band
   non-parabolicity in efshoot.

 
   Paul Harrison, December 1996				 */


#include <stdio.h>
#include <strings.h>
#include <math.h>
#include "struct.h"
#include "maths.h"
#include "bools.h"
#include "const.h"

main(int argc,char *argv[])
{
double	  dV;		   /* Delta V, total band discontinuity */
double	  Eg;		   /* Bandgap				*/
double	  V;		   /* CB and VB potentials		*/
double	  x;	           /* alloy concentration               */
double	  y;	           /* alloy concentration               */
double    z;               /* growth direction                  */
int       n;               /* return value of fopen             */
char	  material[9];	   /* material string			*/
char	  Material;	   /* material character		*/
char	  p;		   /* particle (e, h, or l)		*/
FILE      *Fx;
FILE      *Fv;
FILE      *FEg;		/* pointer to Eg.r file			*/
boolean	  Eg_flag;	/* print bandgap or not			*/

/* Define global defaults */

Material='a';
p='e';
Eg_flag=false;

while((argc>1)&&(argv[1][0]=='-'))
{
 switch(argv[1][1])
 {
  case 'g':
	   Eg_flag=true;
 	   argv--;
	   argc++;
	   break;
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
  case 'p':
           p=*argv[2];
           switch(p)
           {
            case 'e': break;
            case 'h': break;
            case 'l': break;
            default:  printf("Usage:  efxv [-p particle (\033[1me\033[0m, h, or l)]\n");
                      exit(0);
           }
           break;
  default:
	  printf("Usage:  efxv [-M material \033[1mGa(1-x)Al(x)As\033[0m][-p particle (\033[1me\033[0m, h, or l)]\n");
	  printf("             [-g output bandgap \033[1mfalse\033[0m]\n");
	  exit(0);

 }
 argv++;
 argv++;
 argc--;
 argc--;
}

/* If either of the reference potential files exist, i.e., v0.r---the zero
   electric field potential file, or v1.r---the zero dopant reference, then
   remove them.  Thus ensuring each time a new structure is designed, 
   existing files are handled correctly.	*/

unlink("v0.r");	unlink("v1.r");

if((Fx=fopen("x.r","r"))==NULL)
 {printf("Error: Cannot open input file 'x.r'!\n");exit(0);}

if((Fv=fopen("v.r","w"))==NULL)
 {printf("Error: Cannot open file 'v.r'!\n");exit(0);}

if(Eg_flag)
 if((FEg=fopen("Eg.r","w"))==NULL)
  {printf("Error: Cannot open file 'Eg.r'!\n");exit(0);}


switch(Material)
{
 case 'a':	/* Ga(1-x)Al(x)As	*/

          while((fscanf(Fx,"%lf %lf %lf\n",&z,&x,&y))!=EOF)
          {
           dV=(1.247*x)*e_0;
	    switch(p)
	    {
	     case 'e':V=0.67*dV;break;
	     case 'h':V=0.33*dV;break;
	     case 'l':printf("Data not defined for Ga(1-x)Al(x)As light-hole\n");
		      exit(0);
	    }
           fprintf(Fv,"%20.17e %20.17e\n",z,V);
	   if(Eg_flag){Eg=1.426*e_0+dV;fprintf(FEg,"%20.17e %20.17e\n",z,Eg);}
          }
          break;

 case 'b':	/* Cd(1-x)Mn(x)Te	*/

	  while((fscanf(Fx,"%lf %lf %lf\n",&z,&x,&y))!=EOF)
          {
           dV=(1.587*x)*e_0;
	    switch(p)
	    {
	     case 'e':V=0.70*dV;break;
	     case 'h':V=0.30*dV;break;
	     case 'l':printf("Data not defined for Cd(1-x)Mn(x)Te light-hole\n");
		      exit(0);
	    }
           fprintf(Fv,"%20.17e %20.17e\n",z,V);
	   if(Eg_flag){Eg=1.606*e_0+dV;fprintf(FEg,"%20.17e %20.17e\n",z,Eg);}
          }
          break;

 case 'c':	/* In(1-x-y)Al(x)Ga(y)As, Landolt & Bornstein, III/22a, p156 */

	  while((fscanf(Fx,"%lf %lf %lf\n",&z,&x,&y))!=EOF)
          {
           dV=(2.093*x+0.629*y+0.577*x*x+0.436*y*y+1.013*x*y
               -2.0*x*x*(1-x-y))*e_0;
	    switch(p)
	    {
	     /* 53% gives an offset with AlAs of 1.2 eV---close to that 
	        of Hirayama which takes account of strain */
	     case 'e':V=0.53*dV;break;  
	     case 'h':V=0.47*dV;break;
	     case 'l':printf("Data not defined for In(1-x-y)Al(x)Ga(y)As light-hole\n");
		      exit(0);
	    }
           fprintf(Fv,"%20.17e %20.17e\n",z,V);
	   if(Eg_flag){Eg=0.36*e_0+dV;fprintf(FEg,"%20.17e %20.17e\n",z,Eg);}
          }
          break;
}
   
fclose(Fx);
fclose(Fv);
if(Eg_flag)fclose(FEg);


}        /* end main */

