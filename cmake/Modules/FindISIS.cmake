# (Keep your ISISROOT logic at the top...)

# 1. Find ISIS Core
find_path(ISIS_INCLUDE_DIR
  NAMES    Isis.h
  HINTS    "${ISISROOT}/include/isis" "${ISISROOT}/include"
)

# Search for the cspice dir
find_path(CSPICE_INCLUDE_DIR
  NAMES  SpiceUsr.h SpiceCK.h
  HINTS  "${ISISROOT}/include/naif" "${ISISROOT}/include/cspice"
)
if(ISIS_INCLUDE_DIR AND CSPICE_INCLUDE_DIR)
    set(ISIS_INCLUDE_DIRS ${ISIS_INCLUDE_DIR} ${CSPICE_INCLUDE_DIR})
    # Export to parent scope
    set(ISIS_INCLUDE_DIRS ${ISIS_INCLUDE_DIRS} PARENT_SCOPE)
endif()

set(ISIS_LIBRARIES ${ISIS_LIBRARIES} CACHE INTERNAL "")
set(ISIS_FOUND TRUE)