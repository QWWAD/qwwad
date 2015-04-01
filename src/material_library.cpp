/**
 * \file     material_library.cpp
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

#include "material_library.h"
#include "qclsim-material.h"
#include "qclsim-material-property.h"
#include <cstdlib>
#include <stdexcept>
#include <sstream>
#include <assert.h>
#include <cstring>
#include <gsl/gsl_math.h>
#include "qclsim-maths.h"

#include "qclsim-constants.h"
using namespace Leeds;
using namespace Leeds::constants;

typedef xmlpp::Node::NodeList::iterator NodeListIter;

const Glib::ustring & MaterialLibrary::get_property_unit(Glib::ustring &mat_name,
                                                         Glib::ustring &property_name)
{
    MaterialProperty *property = get_property(mat_name, property_name);
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
    if(fname == "")
    {
        std::stringstream fname_str;
        fname_str << QWWAD_PKGDATADIR << "/material-library.xml";
        fname_str >> fname;
    }

    parser = new xmlpp::DomParser(fname, true);

    doc = parser->get_document();
    root_element = doc->get_root_node();

    if(root_element)
    {
        // Get a list of all known materials from the XML file
        material_nodes = root_element->get_children("material");

        // Iterate through all material nodes and add materials to the list
        for(NodeListIter ielem = material_nodes.begin(); ielem != material_nodes.end(); ++ielem)
        {
            // Check that the node is really an element
            xmlpp::Element *elem = dynamic_cast<xmlpp::Element *>(*ielem);

            if(elem)
            {
                // Add the material to the list
                Glib::ustring name = elem->get_attribute_value("name");
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
 */
Material * MaterialLibrary::get_material(const Glib::ustring &mat_name)
{
    // Look through all materials that we have already parsed
    Material mat;
    return &materials.at(mat_name);
}

MaterialLibrary::~MaterialLibrary()
{
    // TODO: Commented out to hack-fix a double free. Should investigate leaks!
    // delete parser;
}

/**
 * Get a property for a given element
 *
 * \param[in] mat_name      Name of the material
 * \param[in] property_name Name of the property
 *
 * \return The material property object
 */
MaterialProperty * MaterialLibrary::get_property(Glib::ustring &mat_name,
                                                 Glib::ustring &property_name)
{
    return materials.at(mat_name).get_property(property_name);
}

double MaterialLibrary::get_val(Glib::ustring &mat_name,
                                Glib::ustring &property_name)
{
    MaterialProperty *property = materials.at(mat_name).get_property(property_name);

    return property->get_val();
}

/**
 * \param[in] mat_name The name of a material to look up
 *
 * \returns the material with the given name
 */
Material * MaterialLibrary::get_material(const char  *mat_name)
{
    std::string str(mat_name);
    return &materials.at(str);
}
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :
