/**
 * \file     qwwad_material_property.cpp
 * \brief    Queries parameters from the material library
 * \author   Alex Valavanis <a.valavanis@leeds.ac.uk>
 */

#include <iostream>
#include <glibmm/ustring.h>

#include "qwwad/material.h"
#include "qwwad/material-library.h"
#include "qwwad/material-property-numeric.h"
#include "qwwad/material-property-string.h"
#include "qwwad/options.h"

using namespace QWWAD;

/** 
 * \brief Argument values read from the command-line.
 *
 * \details The values are set by the option parser
 */
class MatLibOptions : public Options
{
    public:
        MatLibOptions(int argc, char** argv)
        {
            try
            {
                add_option<std::string>("filename",    "", "Material library file to read. If this is not specified, "
                                                           "the default material library for the system will be used.");
                add_option<std::string>("property,p",      "Name of property to look up.");
                add_option<std::string>("material",        "Name of material to look up.");
                add_option<bool>       ("show-unit,u",     "Show the unit for the property rather than just its value");
                add_option<double>     ("variable,x",   0, "Optional input parameter for properties of the form y=f(x)");
                add_option<bool>       ("list-properties", "List all known property names for material, and then exit");

                make_option_positional("material");
                make_option_positional("property");

                std::string doc = "Queries the value of a property from the material "
                                  "database.";

                add_prog_specific_options_and_parse(argc, argv, doc);	
            }
            catch(std::exception& e)
            {
                std::cerr << e.what() << std::endl;
                exit(EXIT_FAILURE);
            }

            if(get_verbose()) {
                print();
            }
        }

        /**
         * \returns Dump user-options to screen
         */
        void print() const
        {
            std::cout << "Querying property: " << get_option<std::string>("property")
                      << " for " << get_option<std::string>("material")
                      << "." << std::endl;
        }
};

auto main(int argc, char* argv[]) -> int
{
    MatLibOptions opt(argc, argv);
    const auto filename = opt.get_option<std::string>("filename");

    MaterialLibrary lib(filename);

    if(!opt.get_argument_known("material")) {
        std::cerr << "Material was not specified" << std::endl;
        std::cerr << "Usage: " << argv[0] << " material_name property_name" << std::endl;
        exit(EXIT_FAILURE);
    }

    // Read the desired material from the library
    const auto material_name = opt.get_option<std::string>("material");
    const Material *mat;

    try {
        mat = lib.get_material(material_name);
    } catch(std::exception &e) {
        std::cerr << e.what() << std::endl;
        exit(EXIT_FAILURE);
    }

    if(opt.get_option<bool>("list-properties")) {
        auto properties = mat->get_all_properties();

        for(auto property : properties) {
            std::cout << property.first << ":" << std::endl
                                               << "\t|-->\t" << property.second->get_description() << std::endl
                                               << "\t|-->\t" << property.second->get_reference() << std::endl;
        }

        exit(EXIT_SUCCESS);
    }

    const auto property_name = opt.get_option<std::string>("property");

    MaterialProperty const * prop;

    try {
        prop = mat->get_property(property_name);
    } catch(std::exception &e) {
        std::cerr << e.what() << std::endl;
        exit(EXIT_FAILURE);
    }

    const auto * const text_property = dynamic_cast<MaterialPropertyString const *>(prop);

    if(text_property != nullptr) {
        std::cout << text_property->get_text() << std::endl;
    } else {
        const auto * const numeric_property = dynamic_cast<MaterialPropertyNumeric const *>(prop);

        const auto x = opt.get_option<double>("variable");
        std::cout << numeric_property->get_val(x);

        if(opt.get_option<bool>("show-unit")) {
            std::cout << " " << numeric_property->get_unit();
        }

        std::cout << std::endl;
    }

    return EXIT_SUCCESS;
}
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :
