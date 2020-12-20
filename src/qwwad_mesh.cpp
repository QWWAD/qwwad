/**
 * \file   qwwad_mesh.cpp
 * \author Alex Valavanis <a.valavanis@leeds.ac.uk>
 * 
 * \brief  Reads data from a three-column file containing a description of 
 *         each layer in a heterostructure:
 *         - Layer width (in nm, or angstroms)
 *         - Alloy fraction (between 0 and 1)
 *         - Volume doping density (in cm^{-3})
 *
 *         This program creates a Mesh class using the data in the
 *         input file, and then outputs tables of alloy fraction and
 *         doping density at each point in the structure, using a 
 *         user-specified spatial profile.
 */

#include <iostream>
#include "qwwad/file-io.h"
#include "qwwad/mesh.h"
#include "qwwad/options.h"

using namespace QWWAD;

/** 
 * \brief Argument values read from the command-line.
 *
 * \details The values are set by the option parser
 */
class MeshOptions : public Options
{
    public:
        MeshOptions(int argc, char** argv);

        [[nodiscard]] auto get_dz_max() const -> double;

        void print() const;
};

auto
MeshOptions::get_dz_max() const -> double
{
    if(!get_argument_known("dzmax")) {
        throw std::runtime_error("Spatial separation not specified");
    }

    auto result = get_option<double>("dzmax");

    // Override result use spatial resolution setting if specified
    if (get_argument_known("zresmin"))
        result = 1.0 / get_option<double>("zresmin");

    if (result < 0)
        throw std::domain_error("Spatial separation must be positive.");

    result*=1.0e-10; // Convert from angstrom -> m

    return result;
}

/**
 * \brief Constructor: Define and parse all user options
 *
 * \param[in] argc Number of command-line arguments
 * \param[in] argv Array of command-line arguments
 */
MeshOptions::MeshOptions(int argc, char** argv)
{
    std::string description("Generate a mesh of samples of structural data.");

    add_option<double>     ("dzmax",               0.1,         "Maximum separation between spatial points.");
    add_option<double>     ("zresmin",                          "Minimum spatial resolution. Overrides the --dz-max option");
    add_option<size_t>     ("nz1per",                0,         "Number of points (per period) within the structure. "
                                                                "If specified, this overrides the --dzmax and "
                                                                "--zresmin options");
    add_option<size_t>     ("nper,p",                1,         "Number of periods to output");
    add_option<std::string>("layerfile,i",          "s.r",      "Filename from which to read input data.");
    add_option<std::string>("interfacesfile,f", "interfaces.r", "Filename to which interface locations are written.");
    add_option<std::string>("alloyfile,x",      "x.r",          "Filename to which alloy profile is written.");
    add_option<std::string>("dopingfile,d",     "d.r",          "Filename to which doping profile is written.");

    add_prog_specific_options_and_parse(argc, argv, description);

    if (get_option<size_t>("nper") == 0)
        throw std::domain_error("Number of periods must be positive.");

    if(get_verbose())
        print();
}

/**
 * \brief Prints a list of user options to screen
 */
void MeshOptions::print() const
{
    printf("Program options have been set as follows...\n");

    std::cout << " * Number of points per period:      " << get_option<size_t>("nz1per")              << std::endl;
    std::cout << " * Number of periods to output:      " << get_option<size_t>("nper")                <<  std::endl;
    std::cout << " * Filename of input structure:      " << get_option<std::string>("layerfile")      << std::endl;
    std::cout << " * Filename for interface locations: " << get_option<std::string>("interfacesfile") << std::endl;
    std::cout << " * Filename for alloy profile:       " << get_option<std::string>("alloyfile")      << std::endl;
    std::cout << " * Filename for doping profile:      " << get_option<std::string>("dopingfile")     << std::endl;
    std::cout << std::endl;
}

auto main(int argc, char* argv[]) -> int
{
    // Read command-line options
    const MeshOptions opt(argc,argv);

    // Create a new Mesh using input data
    const auto nz_1per = opt.get_option<size_t>("nz1per");
    const auto het = (nz_1per != 0) ?  // Force the number of points per period if specified
                     Mesh::create_from_file(opt.get_option<std::string>("layerfile"),
                                                       opt.get_option<size_t>("nz1per"),
                                                       opt.get_option<size_t>("nper"))
                     :
                     Mesh::create_from_file_auto_nz(opt.get_option<std::string>("layerfile"),
                                                               opt.get_option<size_t>("nper"),
                                                               opt.get_dz_max());

    if(opt.get_verbose())
    {
        std::cout << "Period length:                   " << het->get_period_length() << std::endl
                  << "Number of mesh-cells per period: " << het->get_ncell_1per()    << std::endl
                  << "Actual spatial resolution:       " << het->get_dz() << " m"    << std::endl;

        for(unsigned int iL = 0; iL < het->get_n_layers_total(); iL++)
            printf("Top of layer %u is %e\n", iL, het->get_height_at_top_of_layer(iL));
    }
    
    // Output the index of each interface to file
    write_table(opt.get_option<std::string>("interfacesfile").c_str(), het->get_layer_top_indices());

    // Now output each alloy fraction to file
    std::ofstream stream(opt.get_option<std::string>("alloyfile").c_str());

    const auto ncell  = het->get_ncell();
    const auto z      = het->get_z();
    const auto nalloy = het->get_n_alloy();

    const auto alloy = het->get_x_array();

    for(unsigned int iz = 0; iz < ncell; ++iz)
    {
        stream << std::setprecision(20) << std::scientific << z[iz] << "\t";

        for(unsigned int ialloy = 0; ialloy < nalloy; ++ialloy)
            stream << std::setprecision(20) << std::scientific << alloy.at(iz)[ialloy] << "\t";

        stream << std::endl;
    }

    write_table(opt.get_option<std::string>("dopingfile").c_str(), het->get_z(), het->get_n3D_array());
    delete het;

    return EXIT_SUCCESS;
}
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :
