/**
 * \file   qwwad-material-property-numeric.cpp
 * \brief  A class to handle an individual property of a material
 * \author Alex Valavanis <a.valavanis@leeds.ac.uk>
 */

#include <glibmm/ustring.h>
#include <libxml++/libxml++.h>
#include "qwwad-material-property-numeric.h"

/**
 * Create a material property object from a given XML element
 *
 * \param[in] elem The XML element to translate
 *
 * \details The type of the property is determined by looking
 *          at the contents of the XML element
 */
MaterialPropertyNumeric::MaterialPropertyNumeric(xmlpp::Element *elem) :
    MaterialProperty(elem)
{
    // Read the unit of the parameter
    _unit = elem->get_attribute_value("unit");
}

/// Return the unit for the property
const Glib::ustring & MaterialPropertyNumeric::get_unit() const {
    return _unit;
}
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :
