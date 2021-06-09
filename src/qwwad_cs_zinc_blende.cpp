/**
 * \file   qwwad_cs_zinc_blende.cpp
 * \brief  Crystal Structure Zinc Blende 
 * \author Paul Harrison  <p.harrison@shu.ac.uk>
 * \author Alex Valavanis <a.valavanis@leeds.ac.uk>
 *
 * \details This program generates the atomic positions of a zinc blende 
 *          crystal and writes them in XYZ format to the file atoms.xyz
 */

#include <array>
#include <cstdio>
#include <cmath>
#include <cstdlib>
#include <string>

struct vector
{
 double x;
 double y;
 double z;
};

static const size_t n_lattice = 4; // Number of lattice vectors
static const size_t n_basis   = 2; // Number of basis vectors

static void write_ap(double A0,
                     int    n_x,
                     int    n_y,
                     int    n_z,
                     const std::array<vector,n_lattice> &a,
                     const std::array<vector,n_basis>   &T,
                     const std::string                  &anion,
                     const std::string                  &cation);

auto main(int argc,char *argv[]) -> int
{
double	A0;		/* the lattice constant				*/
int	n_x;		/* number of lattice points along x-axis of cell*/
int	n_y;		/* number of lattice points along y-axis of cell*/
int	n_z;		/* number of lattice points along z-axis of cell*/
std::string cation("GA");	/* cation species				*/
std::string anion("AS");	/* anion species				*/
std::array<vector,n_lattice> a;	/* lattice vectors (a1, a2, a3 plus null vector)*/
std::array<vector,n_basis>   T;	/* basis vector					*/

/* default values	*/

n_x=1;n_y=1;n_z=1;
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
  case 'x':
           n_x=atoi(argv[2]);
           break;
  case 'y':
           n_y=atoi(argv[2]);
           break;
  case 'z':
           n_z=atoi(argv[2]);
           break;
  default :
	   printf("Usage:  cszb [-a anion \033[1mGA\033[0m][-c cation \033[1mAS\033[0m]\n");
	   printf("             [-x # cells along x-axis \033[1m1\033[0m][-y # cells \033[1m1\033[0m][-z # cells \033[1m1\033[0m]\n");
	   printf("             [-A lattice constant (\033[1m5.65\033[0mA)]\n");
	   exit(EXIT_SUCCESS);
 }
 argv++;
 argv++;
 argc--;
 argc--;
}

a[0].x=0;a[0].y=0;a[0].z=0;          /* in order to construct a */
a[1].x=1.0/2;a[1].y=1.0/2;a[1].z=0;    /* cubic crystal, this alg-*/
a[2].x=0;a[2].y=1.0/2;a[2].z=1.0/2;    /* orithm requires null ve-*/
a[3].x=1.0/2;a[3].y=0;a[3].z=1.0/2;    /* ctor-called a[0]        */

T[0].x=-1.0/8;T[0].y=-1.0/8;T[0].z=-1.0/8;
T[1].x=+1.0/8;T[1].y=+1.0/8;T[1].z=+1.0/8;

write_ap(A0,n_x,n_y,n_z,a,T,anion,cation);

return EXIT_SUCCESS;
}/* end main */

/**
 * \param A0     lattice constant
 * \param n_x    number of lattice points along x-axis of cell
 * \param n_y    number of lattice points along y-axis of cell
 * \param n_z    number of lattice points along z-axis of cell
 * \param a      lattice vectors (a1, a2, a3 plus null vector)
 * \param T      basis vector
 * \param anion  anion species
 * \param cation cation species
 */
static void write_ap(double A0,
                     int    n_x,
                     int    n_y,
                     int    n_z,
                     const std::array<vector,n_lattice> &a,
                     const std::array<vector,n_basis>   &T,
                     const std::string &anion,
                     const std::string &cation)
{
 int    i_a;    /* index across primitive vectors `a' */
 int    i_n_x;  /* index to n_x */
 int    i_n_y;  /* index to n_y */
 int    i_n_z;  /* index to n_z */
 vector t;      /* general vertor representing atom within cell  */
 FILE	*Fap;	/* pointer to output file	*/

 Fap=fopen("atoms.xyz","w");

 /* Write number of atoms to first line of file, then leave blank line	*/

 fprintf(Fap,"%i\n\n",8*n_x*n_y*n_z);

 for(i_n_x=0;i_n_x<n_x;i_n_x++) {
  for(i_n_y=0;i_n_y<n_y;i_n_y++) {
   for(i_n_z=0;i_n_z<n_z;i_n_z++) {
    for(i_a=0;i_a<=3;i_a++) {
      t.x=i_n_x+a[i_a].x+T[0].x;
      t.y=i_n_y+a[i_a].y+T[0].y;
      t.z=i_n_z+a[i_a].z+T[0].z;
      fprintf(Fap,"%s %9.3f %9.3f %9.3f\n",cation.c_str(),t.x*A0,t.y*A0,t.z*A0);
      t.x=i_n_x+a[i_a].x+T[1].x;
      t.y=i_n_y+a[i_a].y+T[1].y;
      t.z=i_n_z+a[i_a].z+T[1].z;
      fprintf(Fap,"%s %9.3f %9.3f %9.3f\n",anion.c_str(),t.x*A0,t.y*A0,t.z*A0);
     }
   }
  }
 }

 fclose(Fap);
}

// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :
