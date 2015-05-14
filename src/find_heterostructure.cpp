/**
 * \file     find_heterostructure.cpp
 * \author   Alex Valavanis <a.valavanis@leeds.ac.uk>
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
#include "qwwad/file-io.h"
#include "qwwad-heterostructure.h"
#include "qwwad-options.h"

using namespace QWWAD;

/** 
 * \brief Argument values read from the command-line.
 *
 * \details The values are set by the option parser
 */
class HeterostructureOptions : public Options
{
    private:
        Unit unit;       ///< Length unit for calculation
    public:
        HeterostructureOptions(int argc, char* argv[]);

        /**
         * \brief Returns the unit of measurement for lengths
         *
         * \returns The unit of measurement for lengths
         */
        Unit get_unit() const {return unit;}

        double get_dz_max() const
        {
            if(vm.count("dz-max") == 0)
                throw std::runtime_error("Spatial separation not specified");

            double result = get_numeric_option("dz-max");

            // Override result use spatial resolution setting if specified
            if (vm.count("res-min") == 1)
                result = 1.0 / get_numeric_option("res-min");

            if (result < 0)
                throw std::domain_error("Spatial separation must be positive.");

            switch(unit)
            {
                case UNIT_NM:
                    result*=1.0e-9;
                    break;
                case UNIT_ANGSTROM:
                    result*=1.0e-10;
            }

            return result;
        }

        void print() const;
};

/**
 * \brief Constructor: Define and parse all user options
 *
 * \param[in] argc Number of command-line arguments
 * \param[in] argv Array of command-line arguments
 */
HeterostructureOptions::HeterostructureOptions(int argc, char* argv[]) :
    unit(UNIT_ANGSTROM)
{
    std::string description("This program generates a mesh (i.e., a set of spatial points)\n"
                            "and tabulates the composition of a heterostructure at each point.\n"
                            "Normally, this will be the very first step in performing a numerical\n"
                            "simulation of a system. The user should create an input file first,\n"
                            "containing a table describing the system, and then run this program.\n"
                            "Multi-period structures can be generated if desired.");

    add_numeric_option("dz-max",              0.1,          "Maximum separation between spatial points.");
    add_numeric_option("res-min",                           "Minimum spatial resolution. Overrides the --dz-max option");
    add_size_option   ("nz-1per",               0,          "Number of points (per period) within the structure. "
                                                            "If specified, this overrides the --dz-max and "
                                                            "--res-min options");
    add_size_option   ("nper,p",                1,          "Number of periods to output");
    add_string_option ("infile,i",          "s.r",          "Filename from which to read input data.");
    add_string_option ("interfaces-file,f", "interfaces.r", "Filename to which interface locations are written.");
    add_string_option ("alloy-file,x",      "x.r",          "Filename to which alloy profile is written.");
    add_string_option ("doping-file,d",     "d.r",          "Filename to which doping profile is written.");
    add_string_option ("unit,u",            "angstrom",     "Set length unit.  Acceptable values are 'A': "
                                                            "Angstroms or 'n': nanometres.");

    add_prog_specific_options_and_parse(argc, argv, description);

    // Parse which unit we're using
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

    if (get_size_option("nper") < 1)
        throw std::domain_error("Number of periods must be positive.");

    if(get_verbose())
        print();
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
    }

    std::cout << " * Unit of length for input: " << unit_string << std::endl;
    std::cout << " * Number of points per period: " << get_size_option("nz-1per") << std::endl;
    std::cout << " * Number of periods to output: " << get_size_option("nper") << std::endl;
    std::cout << " * Filename of input structure: " << get_string_option("infile") << std::endl;
    std::cout << " * Filename for interface locations: " << get_string_option("interfaces-file") << std::endl;
    std::cout << " * Filename for alloy profile: " << get_string_option("alloy-file") << std::endl;
    std::cout << " * Filename for doping profile: " << get_string_option("doping-file") << std::endl;
    std::cout << std::endl;
}

int main(int argc, char* argv[])
{
    // Read command-line options
    const HeterostructureOptions opt(argc,argv);

    // Create a new heterostructure using input data
    const Heterostructure *het = (opt.get_size_option("nz-1per") != 0) ?  // Force the number of points per period if specified
                                 Heterostructure::create_from_file(opt.get_string_option("infile"),
                                                                   opt.get_unit(),
                                                                   opt.get_size_option("nz-1per"),
                                                                   opt.get_size_option("nper"))
                                 :
                                 Heterostructure::create_from_file_auto_nz(opt.get_string_option("infile"),
                                                                           opt.get_unit(),
                                                                           opt.get_size_option("nper"),
                                                                           opt.get_dz_max());

    if(opt.get_verbose())
    {
        std::cout << "Period length:                       " << het->get_period_length() << std::endl
                  << "Number of spatial points per period: " << het->get_nz_1per()       << std::endl
                  << "Actual spatial resolution:           " << het->get_dz() << " m"    << std::endl;

        for(unsigned int iL = 0; iL < het->get_n_layers_total(); iL++)
            printf("Top of layer %u is %e\n", iL, het->get_height_at_top_of_layer(iL));
    }
    
    // Output the index of each interface to file
    write_table(opt.get_string_option("interfaces-file").c_str(), het->get_layer_top_indices());

    std::ofstream stream(opt.get_string_option("alloy-file").c_str());
    for(unsigned int iz = 0; iz < het->get_z().size(); ++iz)
    {
        stream << std::setprecision(20) << std::scientific << het->get_z()[iz] << "\t";

        for(unsigned int ialloy = 0; ialloy < het->get_x_array().at(0).size(); ++ialloy)
            stream << std::setprecision(20) << std::scientific << het->get_x_array().at(iz)[ialloy] << "\t";

        stream << std::endl;
    }

    write_table(opt.get_string_option("doping-file").c_str(), het->get_z(), het->get_n3D_array());

    delete het;

    return EXIT_SUCCESS;
}
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :
