/**
 * \file     matlib-query.cpp
 * \brief    Queries parameters from the material library
 * \author   Alex Valavanis <a.valavanis@leeds.ac.uk>
 * \date     2014-02-07
 */

#include <iostream>
#include "qclsim-material.h"
#include "qclsim-material-property.h"
#include "material_library.h"
#include "qwwad-options.h"
#include <glibmm/ustring.h>

/** 
 * \brief Argument values read from the command-line.
 *
 * \details The values are set by the option parser
 */
class MatLibOptions : public Options
{
    public:
        MatLibOptions(int argc, char* argv[])
        {
            try
            {
                // Specific configuration options for this program
                program_specific_options->add_options()
                    ("filename", po::value<std::string>()->default_value(""), 
                     "Material library file to read. If this is not specified, "
                     "the default material library for the system will be used.")

                    ("property,p",
                     po::value<std::string>()->default_value(""),
                     "Name of property to look up.")

                    ("material",
                     po::value<std::string>()->default_value(""),
                     "Name of material to look up.")

                    ("text,t",
                     po::bool_switch()->default_value(false),
                     "Use this flag if the property is a string of text")

                    ("show-unit,u",
                     po::bool_switch()->default_value(false),
                     "Show the unit for the property rather than just its value")

                    ("variable,x",
                     po::value<double>()->default_value(0),
                     "Optional input parameter for properties of the form y=f(x)")
                    ;

                std::string doc = "Queries the value of a property from the material "
                                  "database.";

                add_prog_specific_options_and_parse(argc, argv, doc);	
            }
            catch(std::exception& e)
            {
                std::cerr << e.what() << std::endl;
                exit(EXIT_FAILURE);
            }

            if(get_verbose())
                print();
        }

        /**
         * \brief Returns the filename for the material database
         *
         * \returns The filename for the input data
         */
        Glib::ustring get_filename() const {return vm["filename"].as<std::string>();}

        /**
         * \returns True if the property is text; false if it's a number
         */
        bool is_text() const {return vm["text"].as<bool>();}

        /**
         * \returns True if the user wants to see the unit for the property
         */
        bool show_unit() const {return vm["show-unit"].as<bool>();}

        /**
         * \returns the name of the property we're inspecting
         */
        std::string get_property() const {return vm["property"].as<std::string>();}

        /**
         * \returns the name of the property we're inspecting
         */
        double get_var() const {return vm["variable"].as<double>();}

        /**
         * \returns Dump user-options to screen
         */
        void print() const
        {
            std::cout << "Querying property: " << get_property()
                      << " for " << get_string_option("material")
                      << "." << std::endl;
        }
};

int main(int argc, char* argv[])
{
    MatLibOptions opt(argc, argv);
    MaterialLibrary lib(opt.get_filename());

    Material          mat  = lib.get_material(opt.get_string_option("material"));
    MaterialProperty *prop = mat.get_property(opt.get_property().c_str());

    if(opt.is_text())
        std::cout << prop->get_text() << std::endl;
    else
    {
        std::cout << prop->get_val(opt.get_var());

        if(opt.show_unit())
            std::cout << prop->get_unit();

        std::cout << std::endl;
    }

    return EXIT_SUCCESS;
}
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :
