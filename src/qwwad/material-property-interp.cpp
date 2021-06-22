/**
 * \file   material-property-interp.cpp
 * \brief  An individual property of a material that is interpolated between two values
 * \author Alex Valavanis <a.valavanis@leeds.ac.uk>
 */

#include <stdexcept>
#include <libxml++/libxml++.h>
#include "material-property-interp.h"
#include "maths-helpers.h"

namespace QWWAD {
/**
 * \brief Initialise a new interpolated material property
 *
 * \param[in] elem An XML element containing the property data
 */
MaterialPropertyInterp::MaterialPropertyInterp(xmlpp::Element *elem) :
    MaterialPropertyNumeric(elem),
    _xmin(0),
    _xmax(1)
{
    auto all_interp_nodes = elem->get_children("interp");

    // Check for interpolated values
    if(all_interp_nodes.size() == 1)
    {
        auto *interp_node = dynamic_cast<xmlpp::Element *>(all_interp_nodes.front());

        if(interp_node != nullptr)
        {
            auto y0_nodes = interp_node->get_children("y0");
            auto y1_nodes = interp_node->get_children("y1");
            auto *y0_node  = dynamic_cast<xmlpp::Element *>(y0_nodes.front());
            auto *y1_node  = dynamic_cast<xmlpp::Element *>(y1_nodes.front());

            if(y0_node != nullptr && y1_node != nullptr)
            {
                auto *y0_txtnode = y0_node->get_child_text();
                auto *y1_txtnode = y1_node->get_child_text();
                std::stringstream y0_str(y0_txtnode->get_content().raw());
                std::stringstream y1_str(y1_txtnode->get_content().raw());

                y0_str >> _y0;
                y1_str >> _y1;

                // Read optional bowing factor if specified
                auto b_nodes = interp_node->get_children("bow");

                if(!b_nodes.empty())
                {
                    auto *b_node = dynamic_cast<xmlpp::Element *>(b_nodes.front());
                    auto *b_txtnode = b_node->get_child_text();
                    std::stringstream b_str(b_txtnode->get_content().raw());
                    b_str >> _b;
                }
                else {
                    _b = 0.0;
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

/**
 * Create a material property object using specified values
 *
 * \param[in] name        The name of the property
 * \param[in] description A description of the property
 * \param[in] reference   A literature reference for the property
 * \param[in] unit        The unit associated with the property
 * \param[in] y0          The value of the property when \f$x=0\f$ 
 * \param[in] y1          The value of the property when \f$x=1\f$ 
 * \param[in] b           The bowing parameter
 *
 * \details The limits for interpolation are set as \f$x = (0.0, 1.0)\f$
 */
MaterialPropertyInterp::MaterialPropertyInterp(const decltype(_name)        &name,
                                               const decltype(_description) &description,
                                               const decltype(_reference)   &reference,
                                               decltype(_unit)               unit,
                                               decltype(_y0)                 y0,
                                               decltype(_y1)                 y1,
                                               decltype(_b)                  b)
    :
    MaterialPropertyNumeric(name, description, reference, std::move(unit)),
    _y0(y0),
    _y1(y1),
    _b(b),
    _xmin(0.0),
    _xmax(1.0)
{}

auto MaterialPropertyInterp::clone() const -> MaterialPropertyInterp*
{
    auto *mat = new MaterialPropertyInterp(_name, _description, _reference, _unit, _y0, _y1, _b);
    mat->set_limits(_xmin, _xmax);
    return mat;
}

auto
MaterialPropertyInterp::get_val(const double x) const -> decltype(_y0)
{
    if (x < _xmin or x > _xmax)
    {
        std::ostringstream oss;
        oss << "x-value " << x << " is outside the permitted range (" << _xmin << "," << _xmax << ") for property " << _name << std::endl;
        throw std::domain_error(oss.str());
    }

    return lin_interp(_y0, _y1, x, _b);
}

/**
 * Set the validity limits for the interpolation
 *
 * \param[in] xmin The lower limit for interpolation
 * \param[in] xmax The upper limit for interpolation
 */
void MaterialPropertyInterp::set_limits(const decltype(_xmin) xmin,
                                        const decltype(_xmax) xmax)
{
    if (xmin >= xmax)
    {
        std::ostringstream oss;
        oss << "Lower limit: " << xmin << " must be lower than upper limit: " << xmax << std::endl;
        throw std::runtime_error(oss.str());
    }

    _xmin = xmin;
    _xmax = xmax;
}

/**
 * Set the validity limits for the interpolation
 *
 * \param[out] xmin The lower limit for interpolation
 * \param[out] xmax The upper limit for interpolation
 */
void MaterialPropertyInterp::get_limits(decltype(_xmin) &xmin,
                                        decltype(_xmax) &xmax) const
{
    xmin = _xmin;
    xmax = _xmax;
}
} // namespace
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :
