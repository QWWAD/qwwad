/**
 * \file   qclsim-material.cpp
 * \brief  Class to describe a material
 * \author Alex Valavanis <a.valavanis@leeds.ac.uk>
 */
#include <stdexcept>
#include <iostream>
#include "material_library.h"
#include "qclsim-material.h"
#include "qclsim-material-property.h"
#include "qclsim-material-property-interp.h"
#include "qclsim-material-property-poly.h"
#include "qwwad-material-property-constant.h"
#include "qwwad-material-property-string.h"

typedef xmlpp::Node::NodeList::iterator NodeListIter;

/** Return the name of the material */
const Glib::ustring & Material::get_name() const {
    return name;
}
    
/** Return the underlying XML representation */
xmlpp::Element * Material::get_elem() const {
    return elem;
}

/*
 * Return the description of the material
 *
 * \details If a description attribute is not provided, then the name
 *          attribute is used instead.
 *
 * \returns The description as a string
 */
const Glib::ustring & Material::get_description() const
{
    return description;
}

Material::~Material()
{
//    delete priv;
}

Material::Material(const Material *mat)
    : elem(mat->get_elem())
{
    read_properties_from_xml();
}

Material::Material(xmlpp::Element *elem)
    : elem(elem)
{
    read_properties_from_xml();
}

void Material::read_properties_from_xml()
{
    if(elem)
    {
        // Set name and description of this material
        name           = elem->get_attribute_value("name");
        description    = elem->get_attribute_value("description");

        if(description == "")
            description = name;

        // Get a list of all known properties for this material
        property_nodes = elem->get_children("property");

        // Loop through all property nodes and add them to the cache
        for(NodeListIter iprop = property_nodes.begin(); iprop != property_nodes.end(); ++iprop)
        {
            xmlpp::Element *prop = dynamic_cast<xmlpp::Element *>(*iprop);

            if(prop)
            {
                Glib::ustring prop_name = prop->get_attribute_value("name");

                // Figure out what type of property it is and add an appropriate object to the tree
                if(!prop->get_children("interp").empty())
                    properties.insert(prop_name, new MaterialPropertyInterp(prop));
                else if(!prop->get_children("poly").empty())
                    properties.insert(prop_name, new MaterialPropertyPoly(prop));
                else if(!prop->get_children("string").empty())
                    properties.insert(prop_name, new MaterialPropertyString(prop));
                else if(prop->has_child_text())
                    properties.insert(prop_name, new MaterialPropertyConstant(prop));
                else
                {
                    std::ostringstream oss;
                    oss << "Property type could not be recognised for " << prop_name << " in material " << name << std::endl;
                    throw std::runtime_error(oss.str());
                }
            }
        }
    }
    else
        throw std::runtime_error("Invalid XML element");
}

MaterialProperty * Material::get_property(const char* name)
{
    Glib::ustring prop_name(name);
    return get_property(prop_name);
}

/**
 * Get a property node
 *
 * \param[in] property_name Name of the property
 *
 * \return The XML node for this material property
 */
MaterialProperty * Material::get_property(Glib::ustring &property_name)
{
    return &properties.at(property_name);
}
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :
