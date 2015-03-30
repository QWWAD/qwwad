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
#include <boost/ptr_container/ptr_vector.hpp>
#include <libxml++/libxml++.h>

#include "qclsim-constants.h"
using namespace Leeds;
using namespace Leeds::constants;

typedef xmlpp::Node::NodeList::iterator NodeListIter;

class MaterialLibraryImpl {
public:
    MaterialLibraryImpl(const Glib::ustring &filename);
    ~MaterialLibraryImpl();

    boost::ptr_vector<Material>  materials;
    xmlpp::Node::NodeList        material_nodes;
    xmlpp::DomParser            *parser;
    xmlpp::Document             *doc;
    xmlpp::Element              *root_element;
};

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
    priv = new MaterialLibraryImpl(filename);
}

MaterialLibraryImpl::MaterialLibraryImpl(const Glib::ustring &filename)
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
    for(boost::ptr_vector<Material>::iterator imat = priv->materials.begin(); imat != priv->materials.end(); ++imat)
    {
        if(imat->get_name() == mat_name)
        {
            return &(*imat);
        }
    }

    // If we haven't already seen the material then look in the XML data
    xmlpp::Element *elem = 0;

    // Look through the XML material elements until we find one whose
    // name matches the one we want
    NodeListIter ielem = priv->material_nodes.begin();

    while(ielem != priv->material_nodes.end() and !elem) 
    {
        elem = dynamic_cast<xmlpp::Element *>(*ielem);

        if(elem)
        {
            const Glib::ustring name = elem->get_attribute_value("name");

            if(name != mat_name)
                elem = 0;
        }

        ++ielem;
    }

    if(elem)
    {
        // Add the material to the cache and return it
        priv->materials.push_back(new Material(elem));
        return &(priv->materials.back());
    }
    else
    {
        std::ostringstream oss;
        oss << "Couldn't find material " << mat_name;
        throw std::runtime_error(oss.str());
        return 0;
    }
}

MaterialLibraryImpl::~MaterialLibraryImpl()
{
    delete parser;
}

MaterialLibrary::~MaterialLibrary()
{
    delete priv;
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
    Material *mat = get_material(mat_name);
    return mat->get_property(property_name);
}

double MaterialLibrary::get_val(Glib::ustring &mat_name,
                                Glib::ustring &property_name)
{
    Material         *material = get_material(mat_name);
    MaterialProperty *property = material->get_property(property_name);

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
    return get_material(str);
}
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :
