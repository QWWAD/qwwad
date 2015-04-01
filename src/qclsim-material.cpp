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

//                std::cout << name << ":" << prop_name << std::endl;

                // Figure out what type of property it is
                if(!prop->get_children("interp").empty())
                {
  //                  std::cout << "...adding as interp node" << std::endl;
                    properties.insert(prop_name, new MaterialPropertyInterp(prop));
                }
                else
                    properties.insert(prop_name, new MaterialProperty(prop));
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
