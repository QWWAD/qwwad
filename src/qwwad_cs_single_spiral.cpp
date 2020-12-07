/**
 * \file   qwwad_cs_single_spiral.cpp
 * \brief  Crystal Structure Single Spiral
 * \author Paul Harrison  <p.harrison@shu.ac.uk>
 * \author Alex Valavanis <a.valavanis@leeds.ac.uk>
 *
 * \details This program generates the atomic positions of a single spiral
 *          along the z-axis of a zinc blende crystal and writes them in 
 *          XYZ format to the file zb.xyz
 */

#include <array>
#include <cstdio>
#include <cmath>
#include <cstdlib>
#include "struct.h"

const static unsigned int n_basis = 4; // Number of basis vectors

static void write_ap(const double A0,
                     const int    n_z,
                     vector       a,
                     const std::array<vector, n_basis> &T,
                     const std::string                 &anion,
                     const std::string                 &cation);

auto main(int argc,char *argv[]) -> int
{
double	A0;		       /* the lattice constant				*/
int	n_z;		       /* number of lattice points along z-axis of cell*/
std::string cation("GA");      /* cation species				*/
std::string anion("AS");       /* anion species				*/
vector	a;		       /* lattice vectors 				*/
std::array<vector, n_basis> T; // basis vectors

/* default values	*/

n_z=1;
A0=5.65;		/* break all the rules and keep in Angstrom	*/

while((argc>1)&&(argv[1][0]=='-'))
{
 switch(argv[1][1])
 {
  case 'A':
           A0=atof(argv[2]);
           break;
  case 'a':
           anion = argv[2];
           break;
  case 'c':
           cation = argv[2];
           break;
  case 'z':
           n_z=atoi(argv[2]);
           break;
  default :
	   printf("Usage:  csss [-a anion \033[1mAS\033[0m][-c cation \033[1mGA\033[0m]\n");
	   printf("             [-z # cells \033[1m1\033[0m][-A lattice constant (\033[1m5.65\033[0mA)]\n");
	   exit(0);
 }
 argv++;
 argv++;
 argc--;
 argc--;
}

a.x=0;a.y=0;a.z=1.0;          /* in order to construct a */

T[0].x=-1.0/8;T[0].y=-1.0/8;T[0].z=-1.0/8;
T[1].x=+1.0/8;T[1].y=+1.0/8;T[1].z=+1.0/8;
T[2].x=-1.0/8;T[2].y=+3.0/8;T[2].z=+3.0/8;
T[3].x=-3.0/8;T[3].y=+1.0/8;T[3].z=+5.0/8;

write_ap(A0,n_z,a,T,anion,cation);
return EXIT_SUCCESS;
}/* end main */

/**
 * \param A0     lattice constant
 * \param n_z    number of lattice points along z-axis of cell
 * \param a      lattice vectors (a1, a2, a3 plus null vector)
 * \param T      basis vector
 * \param anion  anion species
 * \param cation cation species
 */
static void write_ap(const double A0,
                     const int    n_z,
                     vector       a,
                     const std::array<vector,n_basis> &T,
                     const std::string                &anion,
                     const std::string                &cation)
{
 int    i_n_z;  /* index to n_z */
 vector t;      /* general vector representing atom within cell  */
 FILE	*Fap;	/* pointer to output file	*/

 Fap=fopen("atoms.xyz","w");

 /* Write number of atoms to first line of file, then leave blank line	*/

 fprintf(Fap,"%i\n\n",4*n_z);

 for(i_n_z=0;i_n_z<n_z;i_n_z++)
 {
  t.x=T[0].x;
  t.y=T[0].y;
  t.z=i_n_z*a.z+T[0].z;
  fprintf(Fap,"%s %9.3f %9.3f %9.3f\n",cation.c_str(),t.x*A0,t.y*A0,t.z*A0);
  t.x=T[1].x;
  t.y=T[1].y;
  t.z=i_n_z*a.z+T[1].z;
  fprintf(Fap,"%s %9.3f %9.3f %9.3f\n",anion.c_str(),t.x*A0,t.y*A0,t.z*A0);
  t.x=T[2].x;
  t.y=T[2].y;
  t.z=i_n_z*a.z+T[2].z;
  fprintf(Fap,"%s %9.3f %9.3f %9.3f\n",cation.c_str(),t.x*A0,t.y*A0,t.z*A0);
  t.x=T[3].x;
  t.y=T[3].y;
  t.z=i_n_z*a.z+T[3].z;
  fprintf(Fap,"%s %9.3f %9.3f %9.3f\n",anion.c_str(),t.x*A0,t.y*A0,t.z*A0);
 }

 fclose(Fap);
}

