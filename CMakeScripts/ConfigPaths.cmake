# Define useful directory variables

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

set(INSTALL_LIB_DIR  "${CMAKE_INSTALL_PREFIX}/lib"     CACHE PATH "Installation directory for libraries")
set(INSTALL_INC_DIR  "${CMAKE_INSTALL_PREFIX}/include" CACHE PATH "Installation directory for headers")
add_definitions(-DQWWAD_PKGDATADIR="${CMAKE_INSTALL_PREFIX}/${QWWAD_SHARE_INSTALL}")
