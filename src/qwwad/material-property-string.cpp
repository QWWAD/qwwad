/**
 * \file   material-property-string.cpp
 * \brief  An individual property of a material, represented by a text string
 * \author Alex Valavanis <a.valavanis@leeds.ac.uk>
 */

#include <stdexcept>
#include <libxml++/libxml++.h>
#include "material-property-string.h"

namespace QWWAD {
MaterialPropertyString::MaterialPropertyString(xmlpp::Element *elem) :
    MaterialProperty(elem)
{
    auto string_nodes = elem->get_children("string");

    if(string_nodes.size() == 1) // Parse a constant value
    {
        auto node = dynamic_cast<xmlpp::Element *>(*string_nodes.begin());

        if(node->has_child_text())
            _text = node->get_child_text()->get_content();
        else
        {
            std::ostringstream oss;
            oss << "Couldn't find text for " << _name << "." << std::endl;
            throw std::runtime_error(oss.str());
        }
    }
    else
    {
        std::ostringstream oss;
        oss << "Couldn't parse text string for " << _name << "." << std::endl;
        throw std::runtime_error(oss.str());
    }
}

/**
 * Return the property as a string
 *
 * \return the property value as a raw string
 */
auto MaterialPropertyString::get_text() const -> const Glib::ustring &
{
    return _text;
}
} // end namespace
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :
