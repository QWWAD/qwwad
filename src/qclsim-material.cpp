/**
 * \file   qclsim-material.cpp
 * \brief  Class to describe a material
 * \author Alex Valavanis <a.valavanis@leeds.ac.uk>
 */
#include <stdexcept>
#include <boost/ptr_container/ptr_vector.hpp>
#include <libxml++/libxml++.h>

#include "material_library.h"
#include "qclsim-material.h"
#include "qclsim-material-property.h"

typedef xmlpp::Node::NodeList::iterator NodeListIter;

/**
 * Private members of \c Material class
 */
struct MaterialImpl {
    MaterialImpl(xmlpp::Element *elem);

    /// Cached set of material properties
    boost::ptr_vector<MaterialProperty>   properties;

    xmlpp::Element        *elem;           ///< Underlying XML data
    xmlpp::Node::NodeList  property_nodes; ///< Set of material properties
    Glib::ustring          name;           ///< The name of the material
    Glib::ustring          description;    ///< The description of the material
};

/** Return the name of the material */
const Glib::ustring & Material::get_name() const {
    return priv->name;
}
    
/** Return the underlying XML representation */
xmlpp::Element * Material::get_elem() const {
    return priv->elem;
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
    return priv->description;
}

Material::Material(const Material *mat)
    : priv(new MaterialImpl(mat->get_elem()))
{}

Material::Material(xmlpp::Element *elem)
    : priv(new MaterialImpl(elem))
{}

Material::~Material()
{
//    delete priv;
}

// \todo Add all the property nodes into a map here, using appropriate sub-classes of
// MaterialProperty
MaterialImpl::MaterialImpl(xmlpp::Element *elem)
    : elem(elem)
{
    if(elem)
    {
        // Get a list of all known properties for this material
        property_nodes = elem->get_children("property");
        name           = elem->get_attribute_value("name");
        description    = elem->get_attribute_value("description");

        if(description == "")
            description = name;
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
    // Look through all properties that we have already parsed
    for(boost::ptr_vector<MaterialProperty>::iterator iprop = priv->properties.begin(); iprop != priv->properties.end(); ++iprop)
    {
        if(iprop->get_name() == property_name)
            return &(*iprop);
    }

    // If we haven't already seen this property, then look in the XML data
    xmlpp::Element *elem = 0;

    NodeListIter ielem = priv->property_nodes.begin();

    while(ielem != priv->property_nodes.end() and !elem)
    {
        elem = dynamic_cast<xmlpp::Element *>(*ielem);

        if(elem)
        {
            const Glib::ustring name = elem->get_attribute_value("name");

            if(name != property_name)
                elem = 0;
        }

        ++ielem;
    }

    if(elem)
    {
        // Add the property to the cache and return it
        priv->properties.push_back(new MaterialProperty(elem));
        return &(priv->properties.back());
    }
    else
    {
        std::ostringstream oss;
        oss << "Couldn't find property " << property_name << " for material " << get_name();
        throw std::runtime_error(oss.str());
        return 0;
    }
}
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :
