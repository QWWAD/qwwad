#ifndef CONFIG_H
#define CONFIG_H

#include <string_view>

constexpr std::string_view PACKAGE_NAME      = "${PROJECT_NAME}";
constexpr std::string_view PACKAGE_VERSION   = "${qwwad_VERSION}";
constexpr std::string_view PACKAGE_URL       = "${QWWAD_URL}";
constexpr std::string_view PACKAGE_BUGREPORT = "${QWWAD_BUGREPORT}";

#endif // CONFIG_H

// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :
