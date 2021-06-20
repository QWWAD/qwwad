/**
 * \file     material-library.cpp
 * \brief    Fundamental material data library
 *
 * \details  This is a set of functions for determining properties of a number
 *           of semiconductor materials.  The properties are supposed to be
 *           only things that can be "looked up".  More complicated derived 
 *           properties really belong somewhere else.
 *
 * \author   Alex Valavanis <a.valavanis@leeds.ac.uk>
 * 	     Andrew Grier <el09a2g@leeds.ac.uk>
 */

#include "material-library.h"
#include "material.h"
#include "material-property-numeric.h"
#include <cstdlib>
#include <stdexcept>
#include <sstream>
#include <cassert>
#include <cstring>
#include <gsl/gsl_math.h>
#include "maths-helpers.h"
#include "constants.h"

namespace QWWAD {
using namespace constants;

auto MaterialLibrary::get_property_unit(Glib::ustring &mat_name,
                                        Glib::ustring &property_name) const -> const Glib::ustring &
{
    const auto *property = dynamic_cast<MaterialPropertyNumeric const *>(get_property(mat_name, property_name));
    return property->get_unit();
}

/**
 * Constructor loads material data from XML file
 *
 * param[in] filename Name of input file
 */
MaterialLibrary::MaterialLibrary(const Glib::ustring &filename)
{
    std::string fname(filename);
    // If no filename was specified, read from default data file
    if(fname.empty())
    {
        std::stringstream fname_str;
        fname_str << QWWAD_PKGDATADIR << "/material-library.xml";
        fname_str >> fname;
    }

    xmlpp::DomParser parser(fname, true);

    auto *doc          = parser.get_document();
    auto *root_element = doc->get_root_node();

    if(root_element != nullptr)
    {
        // Get a list of all known materials from the XML file
        auto material_nodes = root_element->get_children("material");

        // Iterate through all material nodes and add materials to the list
        for(auto *node : material_nodes)
        {
            // Check that the node is really an element
            auto *elem = dynamic_cast<xmlpp::Element *>(node);

            if(elem != nullptr)
            {
                // Add the material to the list
                auto name = elem->get_attribute_value("name");
                materials.insert(name, new Material(elem));
            }
        }
    }
    else
    {
        throw std::runtime_error("Couldn't find root element");
    }
}

/**
 * Find a material with a given name in the library
 *
 * \param[in] mat_name The name of the material
 *
 * \return The material from the library
 *
 * \throws boost::bad_ptr_container_operation if the material could not be found
 */
auto MaterialLibrary::get_material(const Glib::ustring &mat_name) const -> Material const *
{
    Material const * mat;

    try
    {
        mat = &materials.at(mat_name);
    }
    catch(std::exception &e)
    {
        std::ostringstream oss;
        oss << "Could not find material: " << mat_name << " in the material library" << std::endl;
        throw std::runtime_error(oss.str());
    }

    return mat;
}

/**
 * Get a property for a given element
 *
 * \param[in] mat_name      Name of the material
 * \param[in] property_name Name of the property
 *
 * \return The material property object
 */
auto MaterialLibrary::get_property(Glib::ustring &mat_name,
                                                       Glib::ustring &property_name) const -> MaterialProperty const *
{
    return materials.at(mat_name).get_property(property_name);
}

auto MaterialLibrary::get_val(Glib::ustring &mat_name,
                                Glib::ustring &property_name) -> double
{
    const auto * const property = materials.at(mat_name).get_property(property_name);
    const auto * const numeric_property = dynamic_cast<MaterialPropertyNumeric const *>(property);

    return numeric_property->get_val();
}

/**
 * \param[in] mat_name The name of a material to look up
 *
 * \returns the material with the given name
 */
auto MaterialLibrary::get_material(const char  *mat_name) const -> Material const *
{
    std::string str(mat_name);
    return &materials.at(str);
}
}
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :
