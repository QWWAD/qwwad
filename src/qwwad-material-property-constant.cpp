/**
 * \file   qwwad-material-property-constant.cpp
 * \brief  An individual property of a material, represented by a constant value
 * \author Alex Valavanis <a.valavanis@leeds.ac.uk>
 */

#include <stdexcept>
#include <libxml++/libxml++.h>
#include "qwwad-material-property-constant.h"

typedef xmlpp::Node::NodeList::iterator NodeListIter;

MaterialPropertyConstant::MaterialPropertyConstant(xmlpp::Element *elem) :
    MaterialPropertyNumeric(elem),
    _constant(0.0)
{
    if(elem->has_child_text()) // Parse a constant value
    {
        auto val_node = elem->get_child_text();

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
 * \brief Return the numerical value of the property
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
double MaterialPropertyConstant::get_val(const double /* x */) const
{
    return _constant;
}
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :
