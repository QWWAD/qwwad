#ifndef QWWAD_MATERIAL_PROPERTY_NUMERIC_H
#define QWWAD_MATERIAL_PROPERTY_NUMERIC_H

#include "material-property.h"

namespace QWWAD {
/**
 * A physical property of a material that can be described by a numerical value
 *
 * \details The numerical value may be obtained using the get_val() function.
 *          The relevant unit for the property may be obtained using get_unit().
 */
class MaterialPropertyNumeric : public MaterialProperty {
protected:
    Glib::ustring _unit; ///< The unit of the property

public:
    MaterialPropertyNumeric() = default;
    MaterialPropertyNumeric(xmlpp::Element *elem);
    MaterialPropertyNumeric(const decltype(_name)        &name,
                            const decltype(_description) &description,
                            const decltype(_reference)   &reference,
                            decltype(_unit)               unit);

    [[nodiscard]] auto get_unit() const -> const decltype(_unit) &;

    [[nodiscard]] virtual auto get_val(double x = 0) const -> double = 0;
};
} // end namespace
#endif
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :
