prefix=@CMAKE_INSTALL_PREFIX@
exec_prefix=@CMAKE_INSTALL_PREFIX@
libdir=@INSTALL_LIB_DIR@
includedir=@INSTALL_INC_DIR@

Name: QWWAD
Description: Quantum Wells, Wires and Dots
Version: @QWWAD_VERSION@

Requires: gsl libxml++-2.6 >= @LIBXMLPP_REQUIRED_VERSION@
Cflags: -I${includedir}
Libs: -L${libdir} -lqwwad
