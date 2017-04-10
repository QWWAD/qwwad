message(STATUS "Creating build files in ${CMAKE_CURRENT_BINARY_DIR}")

# Set path for shared resource files (e.g., example scripts)
if(NOT SHARE_INSTALL)
  set(SHARE_INSTALL
      "share"
      CACHE STRING
      "Data file install path. Must be a relative path (from CMAKE_INSTALL_PREFIX), with no trailing slash.")
endif()

set(QWWAD_SHARE_INSTALL "${SHARE_INSTALL}/qwwad")
mark_as_advanced(QWWAD_SHARE_INSTALL)
