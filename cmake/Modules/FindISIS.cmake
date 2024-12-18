# Find ISIS and its dependencies, such as CSPICE, Qt, and patchelf

message(STATUS "ISISROOT is $ENV{ISISROOT}")

find_path(ISIS_INCLUDE_DIR
  NAME   Isis.h
  HINTS  "${ISISROOT}/include/isis"
         "$ENV{ISISROOT}/include/isis"
         "$ENV{HOME}/projects/isis/inc"
         "$ENV{HOME}/packages/isis/inc"
)

if(ISIS_INCLUDE_DIR)
  message(STATUS "Found the ISIS include directory: ${ISIS_INCLUDE_DIR}")
else()
  message(SEND_ERROR "Could not find the ISIS include directory: ${ISIS_INCLUDE_DIR}")
endif()

if ("$ENV{ISISROOT}" STREQUAL "")
  get_filename_component(ISISROOT ${ISIS_INCLUDE_DIR} PATH)
endif()
message(STATUS "ISISROOT is ${ISISROOT}")

# Searching for each library that makes up a component
foreach(COMPONENT ${ISIS_FIND_COMPONENTS})
  string(TOUPPER ${COMPONENT} UPPERCOMPONENT)
  set( ISIS_${UPPERCOMPONENT}_LIBRARY "ISIS_${UPPERCOMPONENT}_LIBRARY-NOTFOUND")

  find_library(ISIS_${UPPERCOMPONENT}_LIBRARY
    NAMES ${COMPONENT}
    HINTS $ENV{ISISROOT}
    PATH_SUFFIXES lib lib64
    )
  
  mark_as_advanced( ISIS_${UPPERCOMPONENT}_LIBRARY )

  if(ISIS_${UPPERCOMPONENT}_LIBRARY)
    set(ISIS_${UPPERCOMPONENT}_FOUND TRUE CACHE INTERNAL "If the ISIS ${UPPERCOMPONENT} library was found")
  endif()

endforeach(COMPONENT)

set(ISIS_INCLUDE_DIRS ${ISIS_INCLUDE_DIR})
set(ISIS_3RD "$ENV{ISISROOT}/lib")

if(ISIS_INCLUDE_DIR)
  set(ISIS_FOUND TRUE)
else(ISIS_INCLUDE_DIR)
  set(ISIS_FOUND FALSE)
endif(ISIS_INCLUDE_DIR)

# Summarize
if(ISIS_FOUND)
  message(STATUS "Found the following ISIS libraries:")
  foreach( COMPONENT ${ISIS_FIND_COMPONENTS})
    string(TOUPPER ${COMPONENT} UPPERCOMPONENT )
    if ( ISIS_${UPPERCOMPONENT}_FOUND )
      message(STATUS "  ${ISIS_${UPPERCOMPONENT}_LIBRARY}")
      set(ISIS_LIBRARIES ${ISIS_LIBRARIES} ${ISIS_${UPPERCOMPONENT}_LIBRARY})
    else ( ISIS_${UPPERCOMPONENT}_FOUND )
      message(SEND_ERROR "Unable to find ${COMPONENT}")
    endif()
  endforeach()
else(ISIS_FOUND)
  message(SEND_ERROR "Unable to find requested ISIS libraries")
endif(ISIS_FOUND)

# Search for the cspice dir
find_path(CSPICE_INCLUDE_DIR
  NAMES  SpiceCK.h
  HINTS  "${CSPICE_ROOT}/naif"
         "$ENV{ISISROOT}/include/cspice"
         "$ENV{HOME}/packages/cspice/naif"
         "$ENV{HOME}/projects/cspice/naif"
         "/usr/local/packages/cspice/naif"
         "/usgs/pkgs/local/v004/include/naif"
)
if(CSPICE_INCLUDE_DIR)
  message(STATUS "Found the CSPICE include directory: ${CSPICE_INCLUDE_DIR}")
else()
  message(SEND_ERROR "Not found the CSPICE include directory: ${CSPICE_INCLUDE_DIR}")
endif()

if ("${CSPICE_ROOT}" STREQUAL "")
  get_filename_component(CSPICE_ROOT ${CSPICE_INCLUDE_DIR} PATH)
endif()
message(STATUS "CSPICE_ROOT is ${CSPICE_ROOT}")

  find_library(lib_tiff
               NAMES tiff
               PATHS $ENV{ISISROOT}/lib
	      )
  find_library(lib_cspice
               NAMES "cspice"
               PATHS $ENV{ISISROOT}/lib
	      )
  find_library(lib_geos
               NAMES "geos"
               PATHS $ENV{ISISROOT}/lib
	       )
  find_library(lib_gsl
               NAMES "gsl"
               PATHS $ENV{ISISROOT}/lib
	      )	      
  find_library(lib_gslcblas
               NAMES "gslcblas"
               PATHS $ENV{ISISROOT}/lib
	       )
  find_library(lib_protobuf
               NAMES "protobuf"
               PATHS $ENV{ISISROOT}/lib)

  find_library(lib_cholmod
               NAMES "cholmod"
               PATHS $ENV{ISISROOT}/lib
	       )
  find_library(lib_amd
               NAMES "amd"
               PATHS $ENV{ISISROOT}/lib
	       )
  find_library(lib_colamd
               NAMES "colamd"
               PATHS $ENV{ISISROOT}/lib
	       )
  find_library(lib_superlu
               NAMES "superlu"
               PATHS $ENV{ISISROOT}/lib
	       )
  find_library(lib_kdu_a63R
               NAMES "kdu_a63R"
               PATHS $ENV{ISISROOT}/lib
	       )
  find_library(lib_blas
               NAMES "blas"
               PATHS $ENV{ISISROOT}/lib
	       )
  find_library(lib_z
               NAMES "z"
               PATHS $ENV{ISISROOT}/lib
               )
  find_library(lib_xerces-c
               NAMES "xerces-c"
               PATHS $ENV{ISISROOT}/lib
	       )
	       
  if(NOT lib_tiff)
      message (STATUS "lib_tiff NOT FOUND in $ENV{ISISROOT}/lib")
  else(NOT lib_tiff)
      message(STATUS "Found libtiff in: $ENV{ISISROOT}/lib")
  endif(NOT lib_tiff)
  
  if(NOT lib_cspice)
      message (STATUS "lib_cspice not found in $ENV{ISISROOT}/lib")
  else(NOT lib_cspice)
      message (STATUS "lib_cspice found in $ENV{ISISROOT}/lib")
  endif(NOT lib_cspice)
  
  if(NOT lib_protobuf)
      message (STATUS "lib_protobuf not found in " $ENV{ISISROOT}/lib)
  else(NOT lib_protobuf)
      message (STATUS "lib_protobuf found in " $ENV{ISISROOT}/lib)
  endif(NOT lib_protobuf)

 if (NOT APPLE)

  find_library(lib_isis
               NAMES "isis"
               PATHS $ENV{ISISROOT}/lib)
  find_library(lib_gfortran
               NAMES "gfortran"
               PATHS $ENV{ISISROOT}/lib)
  find_library(lib_lapack
               NAMES "lapack"
               PATHS $ENV{ISISROOT}/lib)
  find_library(lib_QtXmlPatterns
               NAMES "QtXmlPatterns"
               PATHS $ENV{ISISROOT}/lib)
  find_library(lib_QtXml
               NAMES "QtXml"
               PATHS $ENV{ISISROOT}/lib)
  find_library(lib_QtNetwork
               NAMES "QtNetwork"
               PATHS $ENV{ISISROOT}/lib)
  find_library(lib_QtSql
               NAMES "QtSql"
               PATHS $ENV{ISISROOT}/lib)
  find_library(lib_QtCore
               NAMES "QtCore"
               PATHS $ENV{ISISROOT}/lib)
  find_library(lib_QtSvg
               NAMES "QtSvg"
               PATHS $ENV{ISISROOT}/lib)
  find_library(lib_QtTest
               NAMES "QtTest"
               PATHS $ENV{ISISROOT}/lib)
  find_library(lib_QtOpenGL
               NAMES "QtOpenGL"
               PATHS $ENV{ISISROOT}/lib)
  find_library(lib_QtDBus
               NAMES "QtDBus"
               PATHS $ENV{ISISROOT}/lib)
  find_library(lib_QtWebKit
               NAMES "QtWebKit"
               PATHS $ENV{ISISROOT}/lib)
  find_library(lib_qwt
               NAMES "qwt"
               PATHS $ENV{ISISROOT}/lib)	       
 endif()


if (NOT APPLE)

  # Search for patchelf
  get_filename_component(PATCHELF_GUESS_DIR "${PATCHELF}" PATH)
  find_path(PATCHELF_DIR
    NAMES  patchelf
    HINTS  "${PATCHELF_GUESS_DIR}"
    "$ENV{ISISROOT}/bin"
    )
  set(PATCHELF "${PATCHELF_DIR}/patchelf")
  if (EXISTS ${PATCHELF})
    message(STATUS "Found patchelf: ${PATCHELF}")
  else()
    message(SEND_ERROR "Invalid path to patchelf: ${PATCHELF}")
  endif()
  
  # Search for Qt
  find_path(QT_DIR
    NAMES  qt
    HINTS "$ENV{ISISROOT}/include"
    NO_DEFAULT_PATH
    )
  set(QT_ROOT "${QT_DIR}")
  if (EXISTS "${QT_ROOT}/")
      message(STATUS "Found Qt in: ${QT_ROOT}")
  else()
    message(SEND_ERROR "Invalid path to Qt: ${QT_ROOT}")
  endif()
  
endif()

if (APPLE)
if (EXISTS "${ISIS_3RD}/libcamd.dylib")
message(STATUS "Found ${ISIS_3RD}/libcamd.dylib")
set( CHOLMODLIB_DEPS "-lcamd")
else()
set( CHOLMODLIB_DEPS "")
endif()
endif()


