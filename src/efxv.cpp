/**
 * \file   efxv.cpp
 * \author Paul Harrison  <p.harrison@shu.ac.uk>
 * \author Alex Valavanis <a.valavanis@leeds.ac.uk>
 *
 * \brief Envelope Function x to v
 *
 * \details Converts the structure as defined in terms
 *          of alloy components into a potential profile for either 
 *          electron, light- or heavy-hole.  It has support for multiple 
 *          material systems, ternaries as well as quaternaries.
 *
 *          In addition generation of the bandgap allows for band
 *          non-parabolicity in efshoot.
 */ 

#include <cstdio>
#include <cstdlib>
#include <iostream>
#include "qwwad/constants.h"
#include "qclsim-fileio.h"
#include "qwwad-options.h"

using namespace Leeds;
using namespace QWWAD;
using namespace constants;

/**
 * Handler for command-line options
 */
class EFXVOptions : public Options
{
    private:
        bool auto_mass; ///< Calculate effective mass automatically

        /**
         * \brief Constant value of effective mass.
         *
         * This is only used if the constant-mass approximation is chosen 
         */
        double mass;

    public:
        EFXVOptions(int argc, char* argv[])
        {
            add_string_option("mass,m",     "auto",   "Set a constant effective-mass across the structure "
                                                      "(relative to free electron). "
                                                      "If not specified, the mass is calculated automatically "
                                                      "for all positions in the material.");
            add_string_option("material,M", "gaalas", "Material ID: \"gaalas\" for Ga(1-x)Al(x)As, "
                                                      "\"cdmnte\" for Cd(1-x)Mn(x)Te, or "
                                                      "\"inalgaas\" for In(1-x-y)Al(x)Ga(y)As");
            add_char_option  ("particle,p",      'e', "Particle to be used: 'e', 'h' or 'l'");
            add_string_option("alloyfile",     "x.r", "File from which alloy is read");

            std::string doc("Find material parameters for a given heterostructure. "
                            "The alloy data is read for each point in the system and used to tabulate "
                            "the band-edge profile, effective mass and permittivity.");

            add_prog_specific_options_and_parse(argc, argv, doc);	

            // Parse the effective-mass calculation type            
            std::string mass_arg(vm["mass"].as<std::string>());

            if(!strcmp(mass_arg.c_str(), "auto"))
                auto_mass = true;
            else if(atof(mass_arg.c_str()) > 0.0)
            {
                auto_mass = false;
                mass = atof(mass_arg.c_str());
            }
            else
            {
                std::ostringstream oss;
                oss << "Cannot parse mass type: " << mass_arg;
                throw std::runtime_error(oss.str());
            }
        }

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

        /**
         * \returns true if we calculate effective mass automatically
         */
        bool compute_mass() const {return auto_mass;}

        /**
         * \returns The constant effective mass in the material
         */
        double get_mass() const {return mass;}
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
    std::valarray<double> V;      // Band-edge potential
    std::valarray<double> m;      // Effective mass
    std::valarray<double> mp;     // Effective mass perpendicular to growth
    std::valarray<double> Eg;     // Bandgap
    std::valarray<double> eps_dc; // Low-frequency permittivity [F/m]

    const char *alloyfile = opt.get_string_option("alloyfile").c_str();

    switch(Material)
    {
        case 'a':	/* Ga(1-x)Al(x)As	*/
            {
                read_table(alloyfile, z, x);
                V.resize(z.size());
                Eg.resize(z.size());
                m.resize(z.size());
                mp.resize(z.size());
                eps_dc.resize(z.size());

                std::valarray<double> dV=(1.247*x)*e; // Total band discontinuity

                switch(p)
                {
                    case 'e':
                        {
                            V=0.67*dV;

                            // Mass data: S. Adachi, `GaAs and related materials' */
                            if(opt.compute_mass())
                            {
                                m=(0.067+0.083*x)*me;
                                mp=(0.067+0.083*x)*me;
                            }
                        }
                        break;
                    case 'h':
                        {
                            V=0.33*dV;

                            if(opt.compute_mass())
                            {
                                m=(0.62+0.14*x)*me;
                                mp=(0.62+0.14*x)*me;
                            }
                        }
                        break;
                    case 'l':printf("Data not defined for Ga(1-x)Al(x)As light-hole\n");
                             exit(EXIT_FAILURE);
                }

                Eg     = 1.426*e+dV;
                eps_dc = ((x-1.0)*(-12.9) + x*10.06)*eps0;
            }
            break;

        case 'b':	/* Cd(1-x)Mn(x)Te	*/
            {
                read_table(alloyfile, z, x);
                V.resize(z.size());
                Eg.resize(z.size());
                m.resize(z.size());
                mp.resize(z.size());
                eps_dc.resize(z.size());

                std::valarray<double> dV=(1.587*x)*e;

                switch(p)
                {
                    case 'e':
                        {
                            V=0.70*dV;

                            // Mass data: Long, 23rd Phys. Semicond. p1819
                            if(opt.compute_mass())
                            {
                                m=(0.11+0.067*x)*me;
                                mp=(0.11+0.067*x)*me;
                            }
                        }
                        break;
                    case 'h':
                        {
                            V=0.30*dV;

                            if(opt.compute_mass())
                            {
                                m=(0.60+0.21*x+0.15*x*x)*me;
                                mp=(0.60+0.21*x+0.15*x*x)*me;
                            }
                        }
                        break;
                    case 'l':
                        {
                            if(opt.compute_mass())
                            {
                                m=(0.18+0.14*x)*me;
                                mp=(0.18+0.14*x)*me;
                            }
                            fprintf(stderr, "Warning: Potential data not defined for Cd(1-x)Mn(x)Te light-hole\n");
                        }
                }

                Eg=1.606*e+dV;
                eps_dc = 10.2*eps0; // Just use CdTe value - can't immediately find MnTe in literature (AV)
            }
            break;

        case 'c':	/* In(1-x-y)Al(x)Ga(y)As, Landolt & Bornstein, III/22a, p156 */
            {
                std::valarray<double> y;

                read_table(alloyfile, z, x, y);
                V.resize(z.size());
                Eg.resize(z.size());
                m.resize(z.size());
                mp.resize(z.size());
                eps_dc.resize(z.size());

                std::valarray<double> dV=(2.093*x+0.629*y+0.577*x*x+0.436*y*y+1.013*x*y
                        +2.0*x*x*(x+y-1))*e;

                switch(p)
                {
                    /* 53% gives an offset with AlAs of 1.2 eV---close to that 
                       of Hirayama which takes account of strain */
                    case 'e':
                        {
                            V=0.53*dV;
                            if(opt.compute_mass())
                            {
                                m=(0.0427+0.0685*x)*me;
                                mp=(0.0427+0.0685*x)*me;
                            }
                        }
                        break;  
                    case 'h':
                        {
                            V=0.47*dV;
                            if(opt.compute_mass())
                                fprintf(stderr, "Warning: Mass data not defined for In(1-x-y)Al(x)Ga(y)As light-hole\n");
                        }
                        break;
                    case 'l':
                        {
                            printf("Data not defined for In(1-x-y)Al(x)Ga(y)As light-hole\n");
                            exit(EXIT_FAILURE);
                        }
                }

                Eg=0.36*e+dV;
                std::cerr << "Warning: Unknown dc permittivity for InAlGaAs" << std::endl;
            }
            break;
    }

    write_table("v.r", z, V);
    write_table("Eg.r", z, Eg);
    write_table("eps-dc.r", z, eps_dc);

    if(!opt.compute_mass())
    {
        m = z*0.0 + opt.get_mass()*me;
        mp = z*0.0 + opt.get_mass()*me;
    }

    write_table("m.r", z, m);
    write_table("m_perp.r", z, mp);

    // Find nonparabolicity parameter
    std::valarray<double> alpha = 1.0/Eg;
    write_table("alpha.r", z, alpha);

    return EXIT_SUCCESS;
}
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :
