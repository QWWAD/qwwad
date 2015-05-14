/**
 * \file   qclsim-material-property.cpp
 * \brief  A class to handle an individual property of a material
 * \author Alex Valavanis <a.valavanis@leeds.ac.uk>
 */

#include <stdexcept>
#include <glibmm/ustring.h>
#include <libxml++/libxml++.h>
#include "qwwad/constants.h"
#include "qclsim-material-property.h"
#include "qclsim-maths.h"

using namespace Leeds;
using namespace QWWAD;
using namespace constants;

/**
 * Return the property as a string
 *
 * \return the property as a raw string of text
 *
 * \details This only works if the property has a single
 *          child text element
 */
const Glib::ustring & MaterialProperty::get_text() const
{
    return _text;
}

/**
 * Create a material property object from a given XML element
 *
 * \param[in] elem The XML element to translate
 *
 * \details The type of the property is determined by looking
 *          at the contents of the XML element
 */
MaterialProperty::MaterialProperty(xmlpp::Element *elem) :
    elem(elem),
    _constant(0)
{
    // Read the name and unit of the parameter
    _name = elem->get_attribute_value("name");
    _unit = elem->get_attribute_value("unit");

    if(elem->has_child_text()) // Parse a constant value
    {
        xmlpp::TextNode *val_node = elem->get_child_text();

        // Store raw text value
        _text = val_node->get_content().raw();

        // Parse text as a number
        std::stringstream s(_text);
        s >> _constant;
    }
    else
    {
        std::ostringstream oss;
        oss << "Couldn't parse the property " << _name;
        throw std::runtime_error(oss.str());
    }
}

/**
 * Return the numerical value of the property
 *
 * \param[in] x If the property is a function of x, then you can
 *              specify x here.  If you leave it blank, x=0 is
 *              assumed.  For constant-valued properties, x is
 *              ignored.
 *
 * \returns The value in whichever units are specified in the XML file
 *
 * \todo Make the library unit-aware (i.e., enable unit conversions etc)
 *
 * \todo Make the library type-aware (i.e., "know" that a property
 *       represents an energy/length/time etc) and check that the
 *       unit makes sense.
 */
double MaterialProperty::get_val(const double /* x */) const
{
    return _constant;
}

/// Return the name of the property
const Glib::ustring & MaterialProperty::get_name() const {
    return _name;
}

/// Return the unit for the property
const Glib::ustring & MaterialProperty::get_unit() const {
    return _unit;
}
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :
