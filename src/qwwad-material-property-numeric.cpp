/**
 * \file   qwwad-material-property-numeric.cpp
 * \brief  A class to handle an individual property of a material
 * \author Alex Valavanis <a.valavanis@leeds.ac.uk>
 */

#include <glibmm/ustring.h>
#include <libxml++/libxml++.h>
#include "qwwad-material-property-numeric.h"

namespace QWWAD {
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

/**
 * Create a material property object using specified values
 *
 * \param[in] name        The name of the property
 * \param[in] description A description of the property
 * \param[in] reference   A literature reference for the property
 * \param[in] unit        The unit associated with the property
 */
MaterialPropertyNumeric::MaterialPropertyNumeric(decltype(_name)        name,
                                                 decltype(_description) description,
                                                 decltype(_reference)   reference,
                                                 decltype(_unit)        unit) :
    MaterialProperty(name, description, reference),
    _unit(unit)
{}

/// Return the unit for the property
const decltype(MaterialPropertyNumeric::_unit) &
MaterialPropertyNumeric::get_unit() const
{
    return _unit;
}
} // end namespace
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :
