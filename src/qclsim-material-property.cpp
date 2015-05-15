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
#include "qwwad/maths-helpers.h"

using namespace QWWAD;
using namespace constants;

/**
 * Create a material property object from a given XML element
 *
 * \param[in] elem The XML element to translate
 */
MaterialProperty::MaterialProperty(xmlpp::Element *elem) :
    elem(elem)
{
    // Read the name of the parameter
    _name = elem->get_attribute_value("name");
}

/// Return the name of the property
const Glib::ustring & MaterialProperty::get_name() const {
    return _name;
}
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :
