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

#include <cstdio>
#include <cstdlib>
#include "qclsim-constants.h"
#include "qclsim-fileio.h"
#include "qwwad-options.h"

using namespace Leeds;
using namespace constants;

/**
 * Handler for command-line options
 */
class EFXVOptions : public Options
{
    public:
        EFXVOptions(int argc, char* argv[])
        {
            try
            {
                program_specific_options->add_options()
                    ("material,M", po::value<std::string>()->default_value("gaalas"),
                     "Material ID: \"gaalas\" for Ga(1-x)Al(x)As, \"cdmnte\" "
                     "for Cd(1-x)Mn(x)Te, or \"inalgaas\" for In(1-x-y)Al(x)Ga(y)As")

                    ("particle,p", po::value<char>()->default_value('e'),
                     "Particle to be used: 'e', 'h' or 'l'")

                    ("print-bandgap,g", po::bool_switch()->default_value(false),
                     "Print the bandgap profile to Eg.r")
                    ;

                std::string doc("Find the band-edge profile for a heterostructure and print to. "
                                "v.r.  The alloy profile is read from the file alloy-profile.dat "
                                "which should contain two or three columns:\n"
                                " * position [m]\n"
                                " * alloy fraction 1\n"
                                " * alloy fraction 2 (OPTIONAL)\n" 
                                "For materials where the 2nd alloy component is not needed (E.g., AlGaAs), "
                                "just omit the last column");

                add_prog_specific_options_and_parse(argc, argv, doc);	
            }
            catch(std::exception &e)
            {
                std::cerr << e.what() << std::endl;
                exit(EXIT_FAILURE);
            }
        }

        /**
         * \returns True if we want to output the bandgap
         */
        bool print_bandgap() const {return vm["print-bandgap"].as<bool>();}

        /**
         * \returns The material identifier
         */
        char get_material() const {
            std::string mat_string(vm["material"].as<std::string>());
            char Material;
            if     (mat_string.compare("gaalas")   == 0) Material='a';
            else if(mat_string.compare("cdmnte")   == 0) Material='b';
            else if(mat_string.compare("inalgaas") == 0) Material='c';
            else
            {
                fprintf(stderr,"The only materials defined in the database are\n");
                fprintf(stderr,"Ga(1-x)Al(x)As, Cd(1-x)Mn(x)Te, In(1-x-y)Al(x)Ga(y)As\n");
                exit(EXIT_FAILURE);
            }
            return Material;
        }

        /**
         * \returns the particle ID
         */
        char get_particle() const {return vm["particle"].as<char>();}
};

int main(int argc,char *argv[])
{
    const EFXVOptions opt(argc, argv);

    char  Material = opt.get_material();  // material character
    char  p = opt.get_particle();   // particle (e, h, or l)

    /* If either of the reference potential files exist, i.e., v0.r---the zero
       electric field potential file, or v1.r---the zero dopant reference, then
       remove them.  Thus ensuring each time a new structure is designed, 
       existing files are handled correctly.	*/
    remove("v0.r");
    remove("v1.r");

    std::valarray<double> z;
    std::valarray<double> x;

    switch(Material)
    {
        case 'a':	/* Ga(1-x)Al(x)As	*/
            {
                read_table_xy("alloy-profile.dat", z, x);

                std::valarray<double> dV=(1.247*x)*e;
                std::valarray<double> V(z.size());
                switch(p)
                {
                    case 'e':V=0.67*dV;break;
                    case 'h':V=0.33*dV;break;
                    case 'l':printf("Data not defined for Ga(1-x)Al(x)As light-hole\n");
                             exit(EXIT_FAILURE);
                }
                write_table_xy("v.r", z, V);

                if(opt.print_bandgap())
                {
                    std::valarray<double> Eg = 1.426*e+dV;
                    write_table_xy("Eg.r", z, Eg);
                }
            }
            break;

        case 'b':	/* Cd(1-x)Mn(x)Te	*/
            {
                read_table_xy("alloy-profile.dat", z, x);

                std::valarray<double> dV=(1.587*x)*e;
                std::valarray<double> V(z.size());

                switch(p)
                {
                    case 'e':V=0.70*dV;break;
                    case 'h':V=0.30*dV;break;
                    case 'l':printf("Data not defined for Cd(1-x)Mn(x)Te light-hole\n");
                             exit(EXIT_FAILURE);
                }

                write_table_xy("v.r", z, V);

                if(opt.print_bandgap())
                {
                    std::valarray<double> Eg=1.606*e+dV;
                    write_table_xy("Eg.r", z, Eg);
                }
            }
            break;

        case 'c':	/* In(1-x-y)Al(x)Ga(y)As, Landolt & Bornstein, III/22a, p156 */
            {
                std::valarray<double> y;

                read_table_xyz("alloy-profile.dat", z, x, y);
                std::valarray<double> dV=(2.093*x+0.629*y+0.577*x*x+0.436*y*y+1.013*x*y
                        +2.0*x*x*(x+y-1))*e;
                std::valarray<double> V(z.size());

                switch(p)
                {
                    /* 53% gives an offset with AlAs of 1.2 eV---close to that 
                       of Hirayama which takes account of strain */
                    case 'e':V=0.53*dV;break;  
                    case 'h':V=0.47*dV;break;
                    case 'l':printf("Data not defined for In(1-x-y)Al(x)Ga(y)As light-hole\n");
                             exit(EXIT_FAILURE);
                }

                write_table_xy("v.r", z, V);

                if(opt.print_bandgap())
                {
                    std::valarray<double> Eg=0.36*e+dV;
                    write_table_xy("Eg.r",z,Eg);
                }
            }
            break;
    }

    return EXIT_SUCCESS;
}
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :
