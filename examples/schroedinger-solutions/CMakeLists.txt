# Define installation path
set(QWWAD_EX_SS_INSTALL "${QWWAD_SHARE_INSTALL}/examples/schroedinger-solutions")

message("Installing scripts to ${QWWAD_EX_SS_INSTALL}")

# Install index file
install(FILES index DESTINATION ${QWWAD_EX_SS_INSTALL})

# Find and install all scripts
file(GLOB _SS_SCRIPTS "*.sh")
install(PROGRAMS ${_SS_SCRIPTS} DESTINATION ${QWWAD_EX_SS_INSTALL})
