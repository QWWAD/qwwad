/**
 * \file   qclsim-material-property-poly.cpp
 * \brief  An individual property of a material, represented by a polynomial fitting
 * \author Alex Valavanis <a.valavanis@leeds.ac.uk>
 */

#include <stdexcept>
#include <libxml++/libxml++.h>
#include "qclsim-material-property-poly.h"
#include <cmath>

typedef xmlpp::Node::NodeList::iterator NodeListIter;

MaterialPropertyPoly::MaterialPropertyPoly(xmlpp::Element *elem) :
    MaterialPropertyNumeric(elem),
    _poly_index(std::vector<int>(0)),
    _poly_coeff(std::vector<double>(0))
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
                _poly_index.push_back(i);
                _poly_coeff.push_back(ai);
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

double MaterialPropertyPoly::get_val(const double x) const
{
    double val = 0; // Output value

    for(unsigned int iterm = 0; iterm < _poly_index.size(); ++iterm)
    {
        const auto ai = _poly_coeff[iterm];
        const auto  i = _poly_index[iterm];

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
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :
