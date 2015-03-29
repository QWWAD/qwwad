#ifndef QCLSIM_MATERIAL
#define QCLSIM_MATERIAL

class MaterialProperty;
struct MaterialImpl;

namespace xmlpp {
class Element;
}

namespace Glib {
class ustring;
}

/** Wrapper for XML data for a material */
class Material {
public:
    Material(const Material *mat);
    Material(xmlpp::Element *elem);
    ~Material();

    const Glib::ustring & get_name() const;
    const Glib::ustring & get_description() const;

    MaterialProperty * get_property(const char    *property_name);
    MaterialProperty * get_property(Glib::ustring &property_name);

    xmlpp::Element * get_elem() const;

private:
    MaterialImpl *priv;
};

#endif // QCLSIM_MATERIAL
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :
