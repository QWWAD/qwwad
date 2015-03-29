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

static double convert_to_SI(const double y, const std::string &unit);

struct MaterialPropertyImpl {
    MaterialPropertyImpl(xmlpp::Element *elem);

    double get_val(const double x) const;

    xmlpp::Element       *elem; ///< Underlying XML data
    MaterialPropertyType  type; ///< The way in which the property is specified

    double        _constant; ///< The value of a constant property
    Glib::ustring _text;     ///< Text representation of a constant property

    // Values for interpolated properties
    double _interp_y0; ///< The value when input variable is 0
    double _interp_y1; ///< The value when input variable is 1

    bool   _interp_is_bowed; ///< Whether a bowing factor is used
    double _interp_b;  ///< The optional bowing factor

    // Values for polynomial properties
    std::vector<int>    _poly_index; ///< The index of each term
    std::vector<double> _poly_coeff; ///< The coefficient of each term

    double _xmin; ///< The minimum value of x that can be used in value lookup
    double _xmax; ///< The maximum value of x that can be used in value lookup

    Glib::ustring _name; ///< The name of the property
    Glib::ustring _unit; ///< The unit of the property
};

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
    if(priv->type != MATERIAL_PROPERTY_CONSTANT)
    {
        std::ostringstream oss;
        oss << "Couldn't read a text value for the property " << get_name();
        throw std::runtime_error(oss.str());
    }

    return priv->_text;
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
    priv(new MaterialPropertyImpl(elem))
{}

MaterialProperty::~MaterialProperty()
{
    delete priv;
}

MaterialPropertyImpl::MaterialPropertyImpl(xmlpp::Element *elem) :
    elem(elem),
    _constant(0),
    _interp_y0(0),
    _interp_y1(0),
    _interp_is_bowed(false),
    _interp_b(0),
    _poly_index(std::vector<int>(0)),
    _poly_coeff(std::vector<double>(0)),
    _xmin(0.0),
    _xmax(DBL_MAX)
{
    _name = elem->get_attribute_value("name");
    _unit = elem->get_attribute_value("unit");
    bool parsing_complete = false;

    xmlpp::Node::NodeList all_interp_nodes = elem->get_children("interp");

    // Check for interpolated values
    if(!all_interp_nodes.empty())
    {
        type = MATERIAL_PROPERTY_INTERP;

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

                parsing_complete = true;
            }
        }
        else
        {
            std::ostringstream oss;
            oss << "Couldn't parse interpolation list for " << _name;
            throw std::runtime_error(oss.str());
        }
    } // Finish parsing <interp> tag

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
 *          - Interpolation between two data points (ust <interp> element)
 *
 * \todo Allow non-SI output?
 * \todo Make the library type-aware (i.e., "know" that a property
 *       represents an energy/length/time etc) and check that the
 *       unit makes sense.
 */
double MaterialProperty::get_val(const double x) const
{
    return priv->get_val(x);
}

double MaterialPropertyImpl::get_val(const double x) const
{
    if (x < _xmin or x > _xmax)
    {
        std::ostringstream oss;
        oss << "x-value " << x << " is outside the permitted range (" << _xmin << "," << _xmax << ") for property " << _name << std::endl;
        throw std::runtime_error(oss.str());
    }

    double val = 0;

    switch(type) {
        case MATERIAL_PROPERTY_CONSTANT:
            val = _constant;
            break;
        case MATERIAL_PROPERTY_INTERP:
            val = lin_interp(_interp_y0, _interp_y1, x, _interp_b);
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

    const std::string unit = _unit;
    return convert_to_SI(val, unit);
}

double MaterialProperty::get_interp_y0() const
{
    const std::string unit = get_unit();
    return convert_to_SI(priv->_interp_y0, unit);
}

double MaterialProperty::get_interp_y1() const
{
    const std::string unit = get_unit();
    return convert_to_SI(priv->_interp_y1, unit);
}

double MaterialProperty::get_interp_b() const
{
    const std::string unit = get_unit();
    return convert_to_SI(priv->_interp_b, unit);
}

/**
 * Convert a given parameter to SI units
 *
 * \todo Make this smarter...
 *       * We don't always want SI
 *       * We shouldn't have to mess around with this list whenever we
 *         introduce a new unit!
 */
static double convert_to_SI(const double y, const std::string &unit)
{
    double val = y;

    if(unit == "") // Dimensionless quantity
        val *= 1.0;    // Bit of a hack... just a null operator

    // Other SI quantities
    else if(unit == "K" or unit == "kg m^{-3}" or unit == "J kg^{-1} K^{-1}")
        val *= 1.0;

    // Convert energies to joules
    else if(unit == "eV")
        val *= e;
    else if(unit == "meV")
        val *= 0.001 * e;

    // Convert lengths to metres
    else if(unit == "angstrom")
        val *= 1e-10;

    // Convert other composite units
    else if(unit == "eV/K")
        val *= e;

//    else
  //      std::cerr << "Unrecognised unit: " << unit << std::endl;

    return val;
}
    
/// Return the name of the property
const Glib::ustring & MaterialProperty::get_name() const {
    return priv->_name;
}

/// Return the unit for the property
const Glib::ustring & MaterialProperty::get_unit() const {
    return priv->_unit;
}

MaterialPropertyType MaterialProperty::get_type() const {
    return priv->type;
}

bool MaterialProperty::interp_is_bowed() const {
    return priv->_interp_is_bowed;
}
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :
