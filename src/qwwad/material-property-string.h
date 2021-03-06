#ifndef QWWAD_MATERIAL_PROPERTY_STRING_H
#define QWWAD_MATERIAL_PROPERTY_STRING_H

#include "material-property.h"

namespace QWWAD {
class MaterialPropertyString : public MaterialProperty {
public:
    MaterialPropertyString() = default;
    MaterialPropertyString(xmlpp::Element *elem);

    [[nodiscard]] auto get_text() const -> const Glib::ustring &;

private:
    Glib::ustring _text; ///< The text value of the property
};
} // end namespace
#endif
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :
