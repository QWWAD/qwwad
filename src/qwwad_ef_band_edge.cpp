/**
 * \file   qwwad_ef_band_edge.cpp
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

#include <cstdlib>
#include <iostream>
#include <valarray>

#include "qwwad/constants.h"
#include "qwwad/file-io.h"
#include "qwwad/options.h"

using namespace QWWAD;
using namespace constants;

/**
 * Handler for command-line options
 */
class BandEdgeOptions : public Options
{
    public:
        BandEdgeOptions(int argc, char** argv)
        {
            add_option<double>     ("mass,m",               "Set a constant effective-mass across the structure "
                                                            "(relative to free electron). "
                                                            "If not specified, the mass is calculated automatically "
                                                            "for all positions in the material.");
            add_option<std::string>("material,M", "gaalas", "Material ID: \"gaalas\" for Ga(1-x)Al(x)As, "
                                                            "\"cdmnte\" for Cd(1-x)Mn(x)Te, or "
                                                            "\"inalgaas\" for In(1-x-y)Al(x)Ga(y)As");
            add_option<char>       ("particle,p",            'e',        "Particle to be used: 'e', 'h' or 'l'");
            add_option<std::string>("dcpermittivityfile",    "eps_dc.r", "File containing the dc permittivity");
            add_option<std::string>("alloyfile",             "x.r",      "File from which alloy is read");
            add_option<std::string>("bandedgepotentialfile", "v_b.r",    "File to which band-edge potential is written");

            std::string doc("Find band-edge parameters for a heterostructure");

            add_prog_specific_options_and_parse(argc, argv, doc);	
        }

        /**
         * \returns The material identifier
         */
        [[nodiscard]] char get_material() const {
            const auto mat_string = get_option<std::string>("material");
            char Material;
            if     (mat_string.compare("gaalas")   == 0) Material='a';
            else if(mat_string.compare("cdmnte")   == 0) Material='b';
            else if(mat_string.compare("inalgaas") == 0) Material='c';
            else
            {
                std::cerr << "The only materials defined in the database are "
                             "Ga(1-x)Al(x)As, Cd(1-x)Mn(x)Te and In(1-x-y)Al(x)Ga(y)As" << std::endl;
                exit(EXIT_FAILURE);
            }
            return Material;
        }
};

/**
 * \brief dc relative permittivity values
 *
 * \todo Read from material library
 */
static double const eps_dc_GaAs = 12.9;
static double const eps_dc_AlAs = 10.06;
static double const eps_dc_InAs = 15.15;

/**
 * \brief Heavy-hole effective masses
 *
 * \todo Read from material library
 */
static double const m_hh_GaAs = 0.51;
static double const m_hh_AlAs = 0.76;
static double const m_hh_InAs = 0.41;

/**
 * \brief Light-hole effective masses
 *
 * \todo Read from material library
 */
static double const m_lh_GaAs = 0.082;
static double const m_lh_AlAs = 0.15;
static double const m_lh_InAs = 0.026;

int main(int argc,char *argv[])
{
    const BandEdgeOptions opt(argc, argv);

    const auto Material = opt.get_material(); // material character
    const auto p        = opt.get_option<char>("particle"); // particle (e, h, or l)


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

    // Flag whether effective mass needs computing or setting manually
    bool compute_mass = true;

    if (opt.get_argument_known("mass"))
        compute_mass = false;

    const auto alloyfile = opt.get_option<std::string>("alloyfile");

    switch(Material)
    {
        case 'a':	/* Ga(1-x)Al(x)As	*/
            {
                read_table(alloyfile.c_str(), z, x);
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
                            if(compute_mass)
                            {
                                m=(0.067+0.083*x)*me;
                                mp=(0.067+0.083*x)*me;
                            }
                        }
                        break;
                    case 'h':
                        {
                            V=0.33*dV;

                            if(compute_mass)
                            {
                                m=(0.62+0.14*x)*me;
                                mp=(0.62+0.14*x)*me;
                            }
                        }
                        break;
                    case 'l':
                        {
                            std::cerr << "Data not defined for Ga(1-x)Al(x)As light-hole" << std::endl;
                            exit(EXIT_FAILURE);
                        }
                }

                Eg     = 1.426*e+dV;
                eps_dc = ((1.0-x)*eps_dc_GaAs + x*eps_dc_AlAs)*eps0;
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
                            if(compute_mass)
                            {
                                m=(0.11+0.067*x)*me;
                                mp=(0.11+0.067*x)*me;
                            }
                        }
                        break;
                    case 'h':
                        {
                            V=0.30*dV;

                            if(compute_mass)
                            {
                                m=(0.60+0.21*x+0.15*x*x)*me;
                                mp=(0.60+0.21*x+0.15*x*x)*me;
                            }
                        }
                        break;
                    case 'l':
                        {
                            if(compute_mass)
                            {
                                m=(0.18+0.14*x)*me;
                                mp=(0.18+0.14*x)*me;
                            }
                            std::cerr << "Warning: Potential data not defined for Cd(1-x)Mn(x)Te light-hole" << std::endl;
                        }
                }

                Eg=1.606*e+dV;
                eps_dc = 10.2*eps0; // Just use CdTe value - can't immediately find MnTe in literature (AV)
            }
            break;

        case 'c':	/* In(1-x-y)Al(x)Ga(y)As, Landolt & Bornstein, III/22a, p156 */
            {
                std::valarray<double> y;

                read_table(alloyfile.c_str(), z, x, y);
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
                            if(compute_mass)
                            {
                                m=(0.0427+0.0685*x)*me;
                                mp=(0.0427+0.0685*x)*me;
                            }
                        }
                        break;  
                    case 'h':
                        {
                            V=0.47*dV;
                            if(compute_mass) {
                                // Linearly interpolate between InAs, AlAs and GaAs for now
                                // TODO: Find a more accurate interpolation
                                m = (m_hh_InAs * (1.0-x-y) + m_hh_GaAs * y + m_hh_AlAs * x)*me;
                            }
                        }
                        break;
                    case 'l':
                        {
                            std::cerr << "Warning: Potential data not defined for In(1-x-y)Al(x)Ga(y)As light-hole" << std::endl;

                            if(compute_mass) {
                                // Linearly interpolate between InAs, AlAs and GaAs for now
                                // TODO: Find a more accurate interpolation
                                m = (m_lh_InAs * (1.0-x-y) + m_lh_GaAs * y + m_lh_AlAs * x)*me;
                            }

                        }
                }

                Eg=0.36*e+dV;

                // Linearly interpolate between InAs, AlAs and GaAs for now
                // TODO: Find a more accurate interpolation
                eps_dc = (eps_dc_InAs * (1.0-x-y) + eps_dc_GaAs * y + eps_dc_AlAs * x)*eps0;
            }
            break;
    }

    const auto bandedgepotentialfile = opt.get_option<std::string>("bandedgepotentialfile");
    const auto dcpermittivityfile    = opt.get_option<std::string>("dcpermittivityfile");
    write_table(bandedgepotentialfile.c_str(), z, V);
    write_table("Eg.r", z, Eg);
    write_table(dcpermittivityfile.c_str(), z, eps_dc);

    if(!compute_mass)
    {
        const auto mass = opt.get_option<double>("mass")*me;
        m = z*0.0 + mass;
        mp = z*0.0 + mass;
    }

    write_table("m.r", z, m);
    write_table("m_perp.r", z, mp);

    // Find nonparabolicity parameter
    std::valarray<double> alpha = 1.0/Eg;
    write_table("alpha.r", z, alpha);

    return EXIT_SUCCESS;
}
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :
