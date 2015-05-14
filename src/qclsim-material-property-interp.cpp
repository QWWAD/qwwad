/**
 * \file   qclsim-material-property-interp.cpp
 * \brief  An individual property of a material that is interpolated between two values
 * \author Alex Valavanis <a.valavanis@leeds.ac.uk>
 */

#include <iostream>
#include <stdexcept>
#include <libxml++/libxml++.h>
#include "qclsim-material-property-interp.h"
#include "qwwad/maths-helpers.h"

namespace QWWAD {
/**
 * \brief Initialise a new interpolated material property
 *
 * \param[in] elem An XML element containing the property data
 */
MaterialPropertyInterp::MaterialPropertyInterp(xmlpp::Element *elem) :
    MaterialProperty(elem),
    _xmin(0),
    _xmax(1)
{
    xmlpp::Node::NodeList all_interp_nodes = elem->get_children("interp");

    // Check for interpolated values
    if(all_interp_nodes.size() == 1)
    {
        xmlpp::Element *interp_node = dynamic_cast<xmlpp::Element *>(all_interp_nodes.front());

        if(interp_node)
        {
            xmlpp::Node::NodeList y0_nodes = interp_node->get_children("y0");
            xmlpp::Node::NodeList y1_nodes = interp_node->get_children("y1");
            xmlpp::Element *y0_node = dynamic_cast<xmlpp::Element *>(y0_nodes.front());
            xmlpp::Element *y1_node = dynamic_cast<xmlpp::Element *>(y1_nodes.front());

            if(y0_node and y1_node)
            {
                xmlpp::TextNode *y0_txtnode = y0_node->get_child_text();
                xmlpp::TextNode *y1_txtnode = y1_node->get_child_text();
                std::stringstream y0_str(y0_txtnode->get_content().raw());
                std::stringstream y1_str(y1_txtnode->get_content().raw());

                y0_str >> _interp_y0;
                y1_str >> _interp_y1;

                // Read optional bowing factor
                xmlpp::Node::NodeList b_nodes = interp_node->get_children("bow");

                if(!b_nodes.empty())
                {
                    xmlpp::Element *b_node = dynamic_cast<xmlpp::Element *>(b_nodes.front());
                    _interp_is_bowed = true;
                    xmlpp::TextNode *b_txtnode = b_node->get_child_text();
                    std::stringstream b_str(b_txtnode->get_content().raw());
                    b_str >> _interp_b;
                }

                // Read validity limits
                std::stringstream xmin_str(interp_node->get_attribute_value("xmin").raw());
                std::stringstream xmax_str(interp_node->get_attribute_value("xmax").raw());
                xmin_str >> _xmin;
                xmax_str >> _xmax;
            }
        }
        else
        {
            std::ostringstream oss;
            oss << "Couldn't parse interpolation list for " << _name << "." << std::endl;
            oss << "Expected 1 'interp' element; found " << all_interp_nodes.size() << std::endl;
            throw std::runtime_error(oss.str());
        }
    } // Finish parsing <interp> tag
}

double MaterialPropertyInterp::get_val(const double x) const
{
    if (x < _xmin or x > _xmax)
    {
        std::ostringstream oss;
        oss << "x-value " << x << " is outside the permitted range (" << _xmin << "," << _xmax << ") for property " << _name << std::endl;
        throw std::domain_error(oss.str());
    }

    return lin_interp(_interp_y0, _interp_y1, x, _interp_b);
}

double MaterialPropertyInterp::get_interp_y0() const
{
    return _interp_y0;
}

double MaterialPropertyInterp::get_interp_y1() const
{
    return _interp_y1;
}

double MaterialPropertyInterp::get_interp_b() const
{
    return _interp_b;
}

bool MaterialPropertyInterp::interp_is_bowed() const {
    return _interp_is_bowed;
}
} // namespace
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :
