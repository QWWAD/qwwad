/**
 * \file   material-property-poly.cpp
 * \brief  An individual property of a material, represented by a polynomial fitting
 * \author Alex Valavanis <a.valavanis@leeds.ac.uk>
 */

#include <stdexcept>
#include <libxml++/libxml++.h>
#include "material-property-poly.h"
#include <cmath>

namespace QWWAD {
MaterialPropertyPoly::MaterialPropertyPoly(xmlpp::Element *elem) :
    MaterialPropertyNumeric(elem)
{
    auto poly_nodes = elem->get_children("poly");

    if(poly_nodes.size() == 1) // Parse a polynomial value
    {
        auto polynomial_node  = dynamic_cast<xmlpp::Element *>(poly_nodes.front());

        if(polynomial_node)
        {
            auto poly_terms = polynomial_node->get_children("ai");

            // Loop through all terms of the polynomial
            for(auto term : poly_terms)
            {
                // TODO: Probably need some error checking in here

                auto ai_elem = dynamic_cast<xmlpp::Element *>(term);
                std::stringstream i_str(ai_elem->get_attribute_value("i").raw());
                std::stringstream ai_str(ai_elem->get_child_text()->get_content().raw());

                int    i  = 0;   // Index of polynomial term
                double ai = 0.0; // Polynomial coefficient

                i_str  >> i;
                ai_str >> ai;
                _poly_coeffs[i] = ai;
            }
        }
    }
    else
    {
            std::ostringstream oss;
            oss << "Couldn't parse polynomial fitting for " << _name << "." << std::endl;
            oss << "Expected 1 'poly' element; found " << poly_nodes.size() << std::endl;
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
 * \param[in] poly_coeffs The set of polynomial coefficients
 *
 * \details The coefficients are stored in a map, with the index being used to specify the
 *          order of the polynomial term, and the value being used to store the coefficient.
 */
MaterialPropertyPoly::MaterialPropertyPoly(decltype(_name)        name,
                                           decltype(_description) description,
                                           decltype(_reference)   reference,
                                           decltype(_unit)        unit,
                                           decltype(_poly_coeffs) poly_coeffs) :
    MaterialPropertyNumeric(name, description, reference, unit),
    _poly_coeffs(poly_coeffs)
{}

/**
 * \return a clone of the current object
 */
MaterialPropertyPoly * MaterialPropertyPoly::clone() const
{
    return new MaterialPropertyPoly(_name, _description, _reference, _unit, _poly_coeffs);
}

double MaterialPropertyPoly::get_val(const double x) const
{
    double val = 0; // Output value

    for(auto term : _poly_coeffs)
    {
        const auto  i = term.first;
        const auto ai = term.second;

        // TODO: This fallback mechanism is provided for systems that
        //       don't have new versions of GSL (1.16?). Ultimately,
        //       we should drop support.
#if defined(HAVE_GSL_POW_UINT)
        val += ai * gsl_pow_uint(x, i);
#else
        val += ai * pow(x, i);
#endif
    }

    return val;
}
} // end namespace
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :
