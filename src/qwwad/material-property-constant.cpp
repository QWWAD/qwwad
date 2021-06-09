/**
 * \file   material-property-constant.cpp
 * \brief  An individual property of a material, represented by a constant value
 * \author Alex Valavanis <a.valavanis@leeds.ac.uk>
 */

#include <stdexcept>
#include <libxml++/libxml++.h>
#include "material-property-constant.h"

namespace QWWAD {
/**
 * Initialise a MaterialPropertyConstant from an XML element
 *
 * \param[in] elem An XML element containing the constant data value.
 *
 * \details The XML element must contain a single double-precision number as its child text
 */
MaterialPropertyConstant::MaterialPropertyConstant(xmlpp::Element *elem) :
    MaterialPropertyNumeric(elem),
    _constant(0.0)
{
    if(elem->has_child_text()) // Parse a constant value
    {
        auto *val_node = elem->get_child_text();

        // Read raw text value
        auto _text = val_node->get_content().raw();

        // Parse text as a number
        std::stringstream s(_text);
        s >> _constant;
    }
    else
    {
        std::ostringstream oss;
        oss << "Couldn't parse constant value for " << _name << "." << std::endl;
        throw std::runtime_error(oss.str());
    }
}

/**
 * Create a material property object using specified values
 *
 * \param[in] name        The name of the property
 * \param[in] description A description of the property
 * \param[in] reference   A literature reference for the property
 * \param[in] unit        The unit associated with the property
 * \param[in] value       The constant value of the property
 */
MaterialPropertyConstant::MaterialPropertyConstant(decltype(_name)        name,
                                                   decltype(_description) description,
                                                   decltype(_reference)   reference,
                                                   decltype(_unit)        unit,
                                                   decltype(_constant)    value) :
    MaterialPropertyNumeric(name, description, reference, unit),
    _constant(value)
{}

/**
 * \returns a copy of the current object
 */
auto MaterialPropertyConstant::clone() const -> MaterialPropertyConstant *
{
    return new MaterialPropertyConstant(_name, _description, _reference, _unit, _constant);
}

/**
 * \brief Return the numerical value of the property
 *
 * \param[in] x (unused) This is just a placeholder, since the value is constant
 *
 * \returns The parameter value
 */
auto
MaterialPropertyConstant::get_val(const double /* x */) const -> decltype(_constant)
{
    return _constant;
}
} // end namespace
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :
