/**
 * \file     find_heterostructure.cpp
 * \author   Alex Valavanis <a.valavanis@leeds.ac.uk>
 * \date     2012-08-03
 * 
 * \brief    Front-end for Heterostructure class
 *
 * \details  Reads data from a three-column file containing a description of 
 *           each layer in a heterostructure:
 *           - Layer width (in nm, or angstroms)
 *           - Alloy fraction (between 0 and 1)
 *           - Volume doping density (in cm^{-3})
 *
 *           This program creates a Heterostructure class using the data in the
 *           input file, and then outputs tables of alloy fraction and
 *           doping density at each point in the structure, using a 
 *           user-specified spatial profile.
 */

#include <iostream>
#include "qclsim-fileio.h"
#include "qwwad-heterostructure.h"
#include "qwwad-options.h"

using namespace Leeds;

/** 
 * \brief Argument values read from the command-line.
 *
 * \details The values are set by the option parser
 */
class HeterostructureOptions : public Options
{
    private:
        Unit unit;       ///< Length unit for calculation
        double Ldiff;    ///< Diffusion length [m]
        void scale_diffusion_length(); // Scale the diffusion length to metres 
    public:
        HeterostructureOptions(int argc, char* argv[]);

        /**
         * \brief Returns the filename for the input data
         *
         * \returns The filename for the input data
         */
        std::string get_input_filename() const {return vm["infile"].as<std::string>();}

        /**
         * \brief Returns the filename for alloy profile
         *
         * \returns The filename for the alloy profile
         */
        std::string get_alloy_filename() const {return vm["alloy-file"].as<std::string>();}

        /**
         * \brief Returns the filename for doping profile
         *
         * \returns The filename for the doping profile
         */
        std::string get_doping_filename() const {return vm["doping-file"].as<std::string>();}

        /**
         * \brief Returns the filename for interface locations
         *
         * \returns The filename for the interface locations
         */
        std::string get_interfaces_filename() const {return vm["interfaces-file"].as<std::string>();}

        /**
         * \brief Returns the unit of measurement for lengths
         *
         * \returns The unit of measurement for lengths
         */
        Unit get_unit() const {return unit;}

        /**
         * \brief Returns the diffusion length [m]
         *
         * \returns The diffusion length in metres
         */
        double get_Ldiff() const {return Ldiff;}

        /**
         * \returns The number of spatial points (per period) to output
         */
        size_t get_nz_1per() const {return vm["nz-1per"].as<size_t>();}

        /**
         * \brief Returns the number of periods
         */
        size_t get_nper() const {return vm["nper"].as<size_t>();}

        void print() const;
};

/**
 * \brief Scales the diffusion length according to user-selected measurement unit
 */
void HeterostructureOptions::scale_diffusion_length() 
{
    // Now scale the diffusion length
    switch(unit)
    {
        case UNIT_NM: 
            Ldiff*=1.0e-9; 
            break;
        case UNIT_ANGSTROM:  
            Ldiff*=1.0e-10;
    }
}

/**
 * \brief Constructor: Define and parse all user options
 *
 * \param[in] argc Number of command-line arguments
 * \param[in] argv Array of command-line arguments
 */
HeterostructureOptions::HeterostructureOptions(int argc, char* argv[]) :
    unit(UNIT_ANGSTROM),
    Ldiff(0.0)
{
    try
    {
        // Specific configuration options for this program
        program_specific_options->add_options()
            ("unit,u", po::value<std::string>()->default_value("angstrom"), 
             "Set length unit.  Acceptable values are 'A': "
             "Ångstroms or 'n': nanometres.")

            ("ldiff,l", po::value(&Ldiff)->default_value(0),
             "Set diffusion length.")

            ("nz-1per", po::value<size_t>()->default_value(1000),
             "Number of points (per period) within the structure")

            ("nper,p", po::value<size_t>()->default_value(1),
             "Number of periods to output")

            ("infile,i", po::value<std::string>()->default_value("s.r"),
             "Set filename from which to read input data.")

            ("interfaces-file,f", 
             po::value<std::string>()->default_value("interfaces.dat"),
             "Set filename for interface locations.")

            ("alloy-file,x",
             po::value<std::string>()->default_value("alloy-profile.dat"),
             "Set filename for alloy profile.")

            ("doping-file,d",
             po::value<std::string>()->default_value("doping-profile.dat"),
             "Set filename for doping profile.")
            ;

        std::string doc = "Create alloy and doping profiles for "
            "a specified semiconductor heterostructure.  Interdiffusion "
            "of alloys across interfaces may also be specified.";

        add_prog_specific_options_and_parse(argc, argv, doc);	

        // Perform a bit of post-processing on the options
        {
            if (vm.count("unit"))
            {
                const char* arg = vm["unit"].as<std::string>().c_str();
                if(!strcmp (arg, "A") ||
                        !strcmp(arg, "angstrom") ||
                        !strcmp(arg, "Angstrom") ||
                        !strcmp(arg, "angstroms") ||
                        !strcmp(arg, "Angstroms"))
                    unit = UNIT_ANGSTROM; // Sets unit to [Å]
                else if(!strcmp(arg, "n") ||
                        !strcmp(arg, "nm") ||
                        !strcmp(arg, "nanometre") ||
                        !strcmp(arg, "nanometer") ||
                        !strcmp(arg, "nanometres") ||
                        !strcmp(arg, "nanometers"))
                    unit = UNIT_NM; // Sets unit to [nm]
            }

            if (Ldiff < 0)
                throw std::domain_error("Diffusion length must be positive.");

            if (vm["nper"].as<size_t>() < 1)
                throw std::domain_error("Number of periods must be positive.");
        }
    }
    catch(std::exception& e)
    {
        std::cerr << e.what() << std::endl;
        exit(EXIT_FAILURE);
    }

    if(get_verbose())
        print();

    scale_diffusion_length();
}

/**
 * \brief Prints a list of user options to screen
 */
void HeterostructureOptions::print() const
{
    const char* unit_string = NULL;

    printf("Program options have been set as follows...\n");

    switch(unit)
    {
        case UNIT_NM:
            unit_string="nm";
            break;
        case UNIT_ANGSTROM:
            unit_string="angstroms";
            break;
        default:
            unit_string="unknown";
    }

    printf(" * Unit of length for input: %s\n", unit_string);
    printf(" * Diffusion length: %f %s\n",Ldiff,unit_string);
    printf(" * Number of points per period: %i\n", (int)get_nz_1per());
    printf(" * Number of periods to output: %i\n", (int)get_nper());
    printf(" * Filename of input structure: %s\n", get_input_filename().c_str());
    printf(" * Filename for interface locations: %s\n",get_interfaces_filename().c_str());
    printf(" * Filename for alloy profile: %s\n", get_alloy_filename().c_str());
    printf(" * Filename for doping profile: %s\n", get_doping_filename().c_str());
    printf("\n");
}

int main(int argc, char* argv[])
{
    // Read command-line options
    const HeterostructureOptions opt(argc,argv);

    // Create a new heterostructure using input data
    const Heterostructure *het = Heterostructure::read_from_file(opt.get_input_filename(),
                                                                 opt.get_unit(),
                                                                 opt.get_nz_1per(),
                                                                 opt.get_nper(),
                                                                 opt.get_Ldiff());

    if(opt.get_verbose())
    {
        std::cout << "Period length:                       " << het->get_period_length() << std::endl
                  << "Number of spatial points per period: " << het->get_nz_1per()       << std::endl
                  << "Actual spatial resolution:           " << het->get_dz() << " m"    << std::endl;

        for(unsigned int iL = 0; iL < het->get_n_layers_total(); iL++)
            printf("Top of layer %u is %e\n", iL, het->get_height_at_top_of_layer(iL));
    }
    
    // Output the index of each interface to file
    write_table_x(opt.get_interfaces_filename().c_str(), het->get_layer_top_indices());

    std::ofstream stream(opt.get_alloy_filename().c_str());
    for(unsigned int iz = 0; iz < het->get_z().size(); ++iz)
    {
        stream << std::setprecision(20) << std::scientific << het->get_z()[iz] << "\t";

        for(unsigned int ialloy = 0; ialloy < het->get_x_diffuse_array().at(0).size(); ++ialloy)
            stream << std::setprecision(20) << std::scientific << het->get_x_diffuse_array().at(iz)[ialloy] << "\t";

        stream << std::endl;
    }

    write_table_xy(opt.get_doping_filename().c_str(), het->get_z(), het->get_n3D_array());

    delete het;

    return EXIT_SUCCESS;
}
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :
