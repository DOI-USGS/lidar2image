
# - FindGFortran.cmake
#
# The following variables are set when GFORTRAN is found:
# GFORTRAN_FOUND
# GFORTRAN_LIBS

set(TEMP_CMAKE_FIND_LIBRARY_SUFFIXES ${CMAKE_FIND_LIBRARY_SUFFIXES})
set(CMAKE_FIND_LIBRARY_SUFFIXES .a .so.3)
set(CMAKE_FIND_LIBRARY_PREFIXES lib)

find_library (GFORTRAN_LIBS gfortran
    PATHS /sw /usr/lib /usr/lib64 /usr/local /usr/local/gfortran /opt/local ${CMAKE_INSTALL_PREFIX}
    PATH_SUFFIXES lib gcc41 gcc42 gcc43 gcc44 gcc45 gcc46 gcc48
)
  
if (GFORTRAN_LIBS)
    set (GFORTRAN_FOUND TRUE)
    message("-- GFortran Found: ${GFORTRAN_LIBS}")
else (GFORTRAN_LIBS)
    set (GFORTRAN_FOUND FALSE)
    message("-- GFortran Not Found")
endif (GFORTRAN_LIBS)
  
set(CMAKE_FIND_LIBRARY_SUFFIXES ${TEMP_CMAKE_FIND_LIBRARY_SUFFIXES})
