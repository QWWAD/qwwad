/**
 * \file   material-property.cpp
 * \brief  A class to handle an individual property of a material
 * \author Alex Valavanis <a.valavanis@leeds.ac.uk>
 */

#include <stdexcept>
#include <utility>

#include <glibmm/ustring.h>
#include <libxml++/libxml++.h>
#include "constants.h"
#include "material-property.h"
#include "maths-helpers.h"

namespace QWWAD {
using namespace constants;

/**
 * Create a material property object from a given XML element
 *
 * \param[in] elem The XML element to translate
 *
 * \details The name and description attributes are read from a <property> XML element
 */
MaterialProperty::MaterialProperty(xmlpp::Element *elem)
{
    // Read the name of the parameter
    auto name_str = elem->get_attribute_value("name");

    if(name_str.empty()) {
        throw std::runtime_error("Material property found with no name");
    }

    _name = name_str;

    _description = elem->get_attribute_value("description");
}

/**
 * Create a material property object using specified values
 *
 * \param[in] name        The name of the property
 * \param[in] description A description of the property
 * \param[in] reference   A literature reference for the property
 */
MaterialProperty::MaterialProperty(const decltype(_name)  &name,
                                   decltype(_description)  description,
                                   decltype(_reference)    reference)
    :
        _name(""),
        _description(std::move(description)),
        _reference(std::move(reference))
{
    if(name.empty()) {
        throw std::runtime_error("Material property must have a name");
    }

    _name = name;
}

auto MaterialProperty::clone() const -> MaterialProperty*
{
    return new MaterialProperty(_name, _description, _reference);
}

/**
 * Return the name of the property
 *
 * \return The name as a string
 */
auto
MaterialProperty::get_name() const -> const decltype(MaterialProperty::_name) &
{
    return _name;
}

/**
 * Return the description of the property
 *
 * \return The description as a string
 */
auto
MaterialProperty::get_description() const -> const decltype(MaterialProperty::_description) &
{
    return _description;
}

/**
 * Return the literature reference for the property
 *
 * \return The reference as a string
 */
auto
MaterialProperty::get_reference() const -> const decltype(MaterialProperty::_reference) &
{
    return _reference;
}
} // end namespace
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :
