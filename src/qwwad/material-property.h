#ifndef QCLSIM_MATERIAL_PROPERTY_H
#define QCLSIM_MATERIAL_PROPERTY_H

#include <glibmm/ustring.h>

namespace xmlpp {
class Element;
}

namespace QWWAD {
/**
 * A base class describing a single physical property of a material
 */
class MaterialProperty {
protected:
    Glib::ustring _name;        ///< The name of the property
    Glib::ustring _description; ///< A description of the property
    Glib::ustring _reference;   ///< A literature reference

public:
    MaterialProperty() = default;
    MaterialProperty(xmlpp::Element *elem);
    MaterialProperty(decltype(_name)        name,
                     decltype(_description) description,
                     decltype(_reference)   reference);

    [[nodiscard]] virtual MaterialProperty * clone() const;

    [[nodiscard]] virtual const decltype(_name)        & get_name()        const;
    [[nodiscard]] virtual const decltype(_description) & get_description() const;
    [[nodiscard]] virtual const decltype(_reference)   & get_reference()   const;
};

inline MaterialProperty * new_clone(const MaterialProperty & mat)
{
    return mat.clone();
}
} // end namespace
#endif
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :
