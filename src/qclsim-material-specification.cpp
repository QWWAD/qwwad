#include "material_library.h"
#include "qclsim-material.h"
#include "qclsim-material-property.h"
#include "qwwad-material-property-numeric.h"
#include "qwwad-material-property-string.h"
#include "qclsim-material-specification.h"
#include <glibmm/ustring.h>
#include <stdexcept>

/**
 * Default constructor
 */
MaterialSpecification::MaterialSpecification() :
    lib(new MaterialLibrary("")),
    alloy(0.0),
    n3d(0.0)
{
    xml = new Material(lib->get_material("AlGaAs"));
}

/**
 * Specify a material by providing a list of properties
 */
MaterialSpecification::MaterialSpecification(const char          *mat_name,
                                             const double         alloy,
                                             const double         doping) :
    lib(new MaterialLibrary("")),
    alloy(alloy),
    n3d(doping)
{
    xml = new Material(lib->get_material(mat_name));
}
MaterialSpecification::MaterialSpecification(const Glib::ustring &mat_name,
                                             const double         alloy,
                                             const double         doping) :
    lib(new MaterialLibrary("")),
    alloy(alloy),
    n3d(doping)
{
    xml = new Material(lib->get_material(mat_name));
}

/**
 * Copy constructor
 *
 * \param[in] spec The material specification to be copied
 */
MaterialSpecification::MaterialSpecification(const MaterialSpecification &spec) :
    lib(new MaterialLibrary("")),
    alloy(spec.alloy),
    n3d(spec.n3d)
{
    std::string mat_name(spec.xml->get_name());
    xml = new Material(lib->get_material(mat_name));
}

MaterialSpecification::~MaterialSpecification()
{
    delete xml;
    delete lib;
}

MaterialSpecification & MaterialSpecification::operator=(const MaterialSpecification &rhs)
{
    lib    = new MaterialLibrary("");
    alloy  = rhs.alloy;
    n3d    = rhs.n3d;
    std::string mat_name(rhs.xml->get_name());
    xml    = new Material(lib->get_material(mat_name));

    return *this;
}

/**
 * Get the value of a property, using zero as the argument
 *
 * \param[in] prop_name The name of the property to look up
 *
 * \returns The value of the property
 */
double MaterialSpecification::get_prop_val_0(const char *prop_name) const
{
	Glib::ustring str(prop_name);
	return get_prop_val_0(str);
}

double MaterialSpecification::get_prop_val_0(Glib::ustring &prop_name) const
{
    return get_prop_val(prop_name, 0.0);
}

/**
 * Gets the text from a text property
 *
 * \param[in] prop_name The name of the property to look up
 *
 * \returns The text stored in the property
 */
const Glib::ustring & MaterialSpecification::get_prop_text(Glib::ustring &prop_name) const
{
    const auto property = dynamic_cast<MaterialPropertyString *>(xml->get_property(prop_name));

    if(!property)
    {
        std::ostringstream oss;
        oss << "Could not read text data for property: " << prop_name << std::endl;
        throw std::runtime_error(oss.str());
    }

    return property->get_text();
}

const Glib::ustring & MaterialSpecification::get_prop_text(const char *prop_name) const
{
    Glib::ustring str(prop_name);
    return get_prop_text(str);
}

double MaterialSpecification::get_prop_val(const char   *prop_name,
                                           const double  x) const
{
    Glib::ustring str(prop_name);
    return get_prop_val(str, x);
}

double MaterialSpecification::get_prop_val(Glib::ustring &prop_name,
                                           const double   x) const
{
    const auto property = dynamic_cast<MaterialPropertyNumeric *>(xml->get_property(prop_name));

    if(!property)
    {
        std::ostringstream oss;
        oss << "Could not read numerical data for property: " << prop_name << std::endl;
        throw std::runtime_error(oss.str());
    }

    return property->get_val(x);
}

/**
 * Get the value of a property, using the alloy concentration
 * as the input parameter
 *
 * \param[in] prop_name The name of the property to look up
 *
 * \returns The value of the property
 */
double MaterialSpecification::get_prop_val_x(const char *prop_name) const
{
    Glib::ustring str(prop_name);
    return get_prop_val_x(str);
}

double MaterialSpecification::get_prop_val_x(Glib::ustring &prop_name) const
{
    return get_prop_val(prop_name, alloy);
}
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :
