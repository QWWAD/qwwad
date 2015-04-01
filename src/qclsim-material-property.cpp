/**
 * \file   qclsim-material-property.cpp
 * \brief  A class to handle an individual property of a material
 * \author Alex Valavanis <a.valavanis@leeds.ac.uk>
 */

#include <stdexcept>
#include <glibmm/ustring.h>
#include <libxml++/libxml++.h>
#include "qclsim-constants.h"
#include "qclsim-material-property.h"
#include "qclsim-maths.h"

using namespace Leeds;
using namespace constants;

typedef xmlpp::Node::NodeList::iterator NodeListIter;

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
    if(type != MATERIAL_PROPERTY_CONSTANT)
    {
        std::ostringstream oss;
        oss << "Couldn't read a text value for the property " << get_name();
        throw std::runtime_error(oss.str());
    }

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
    _constant(0),
    _poly_index(std::vector<int>(0)),
    _poly_coeff(std::vector<double>(0))
{
    // Read the name and unit of the parameter
    _name = elem->get_attribute_value("name");
    _unit = elem->get_attribute_value("unit");
    bool parsing_complete = false;

    if(!elem->get_children("poly").empty()) // Parse a polynomial value
    {
        type = MATERIAL_PROPERTY_POLY;
        xmlpp::Node::NodeList poly_nodes = elem->get_children("poly");
        xmlpp::Element *polynomial_node  = dynamic_cast<xmlpp::Element *>(poly_nodes.front());

        if(polynomial_node)
        {
            xmlpp::Node::NodeList poly_terms = polynomial_node->get_children("ai");

            // Loop through all terms of the polynomial
            for(NodeListIter term = poly_terms.begin(); term != poly_terms.end(); ++term)
            {
                // TODO: Probably need some error checking in here

                xmlpp::Element *ai_elem = dynamic_cast<xmlpp::Element *>(*term);
                std::stringstream i_str(ai_elem->get_attribute_value("i").raw());
                std::stringstream ai_str(ai_elem->get_child_text()->get_content().raw());

                int    i  = 0;   // Index of polynomial term
                double ai = 0.0; // Polynomial coefficient

                i_str  >> i;
                ai_str >> ai;
                _poly_index.push_back(i);
                _poly_coeff.push_back(ai);

                parsing_complete = true;
            }
        }
    }

    if(!parsing_complete and elem->has_child_text()) // Parse a constant value
    {
        type = MATERIAL_PROPERTY_CONSTANT;
        xmlpp::TextNode *val_node = elem->get_child_text();
        _text = val_node->get_content().raw();
        std::stringstream s(_text);
        s >> _constant;
        parsing_complete = true;
    }

    if(!parsing_complete)
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
 * \returns The value in SI units
 *
 * \details At the moment, the parser can interpret data in the following
 *          forms:
 *
 *          - Constant values
 *
 * \todo Allow non-SI output?
 * \todo Make the library type-aware (i.e., "know" that a property
 *       represents an energy/length/time etc) and check that the
 *       unit makes sense.
 */
double MaterialProperty::get_val(const double x) const
{
    double val = 0;

    switch(type) {
        case MATERIAL_PROPERTY_CONSTANT:
            val = _constant;
            break;
        case MATERIAL_PROPERTY_POLY:
            for(unsigned int iterm = 0; iterm < _poly_index.size(); ++iterm)
            {
                const double ai = _poly_coeff[iterm];
                unsigned int  i = _poly_index[iterm];

                // TODO: This fallback mechanism is provided for systems that
                //       don't have new versions of GSL (1.16?). Ultimately,
                //       we should drop support.
#if defined(HAVE_GSL_POW_UINT)
                val += ai * gsl_pow_uint(x, i);
#else
                val += ai * pow(x, i);
#endif
            }
            break;
    }

    return val;
}

/// Return the name of the property
const Glib::ustring & MaterialProperty::get_name() const {
    return _name;
}

/// Return the unit for the property
const Glib::ustring & MaterialProperty::get_unit() const {
    return _unit;
}

MaterialPropertyType MaterialProperty::get_type() const {
    return type;
}
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :
