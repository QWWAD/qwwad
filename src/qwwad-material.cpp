/**
 * \file   qwwad-material.cpp
 * \brief  Class to describe a material
 * \author Alex Valavanis <a.valavanis@leeds.ac.uk>
 */
#include <stdexcept>
#include <iostream>
#include "material_library.h"
#include "qwwad-material.h"
#include "qwwad-material-property-interp.h"
#include "qwwad-material-property-poly.h"
#include "qwwad-material-property-constant.h"
#include "qwwad-material-property-string.h"

namespace QWWAD {
/** Return the name of the material */
const Glib::ustring & Material::get_name() const {
    return name;
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

Material::Material(const Material *mat)
{
    properties = mat->properties.clone();
}

Material::Material(xmlpp::Element *elem)
{
    if(elem)
    {
        // Set name and description of this material
        name           = elem->get_attribute_value("name");
        description    = elem->get_attribute_value("description");

        if(description == "")
            description = name;

        // Get a list of all known properties for this material
        auto property_nodes = elem->get_children("property");

        // Loop through all property nodes and add them to the cache
        for(auto node : property_nodes)
        {
            auto prop = dynamic_cast<xmlpp::Element *>(node);

            if(prop)
            {
                auto prop_name = prop->get_attribute_value("name");

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

MaterialProperty const * Material::get_property(const char *name) const
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
MaterialProperty const * Material::get_property(const Glib::ustring &property_name) const
{
    MaterialProperty const * prop;

    try
    {
        prop = &properties.at(property_name);
    }
    catch(std::exception &e)
    {
        std::ostringstream oss;
        oss << "Could not find property: " << property_name << " in the material library" << std::endl;
        throw std::runtime_error(oss.str());
    }

    return &properties.at(property_name);
}

MaterialPropertyNumeric const * Material::get_numeric_property(const char *name) const
{
    Glib::ustring prop_name(name);
    return get_numeric_property(prop_name);
}

MaterialPropertyNumeric const * Material::get_numeric_property(Glib::ustring &name) const
{
    auto prop = dynamic_cast<MaterialPropertyNumeric const *>(get_property(name));
    return prop;
}

double Material::get_property_value(const char   *name,
                                    const double  x) const
{
    Glib::ustring prop_name(name);
    return get_property_value(prop_name, x);
}

/**
 * Get the value of a numerical property
 */
double Material::get_property_value(Glib::ustring &name,
                                    const double   x) const
{
    auto prop = get_numeric_property(name);
    return prop->get_val(x);
}
} // end namespace
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :
