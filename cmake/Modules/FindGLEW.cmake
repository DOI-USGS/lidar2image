#IF(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
#
# Try to find GLEW library and include path.
# Once done this will define
#
# GLEW_FOUND
# GLEW_INCLUDE_PATH
# GLEW_LIBRARY
#

FIND_PATH( GLEW_INCLUDE_PATH GL/glew.h
	/usr/include
	/usr/local/include
	DOC "The directory where GL/glew.h resides")
FIND_LIBRARY( GLEW_LIBRARY
	NAMES GLEW glew
	PATHS
	/usr/lib64
	/usr/lib
	/usr/local/lib64
	/usr/local/lib
        /opt/local/lib
	DOC "The GLEW library")

IF (GLEW_INCLUDE_PATH)
	SET( GLEW_FOUND 1 CACHE STRING "Set to 1 if GLEW is found, 0 otherwise")
ELSE (GLEW_INCLUDE_PATH)
	SET( GLEW_FOUND 0 CACHE STRING "Set to 1 if GLEW is found, 0 otherwise")
ENDIF (GLEW_INCLUDE_PATH)

MARK_AS_ADVANCED( GLEW_FOUND )

if( GLEW_FOUND )
  message( STATUS "  Found GLEW version in ${GLEW_LIBRARY}" )
endif( GLEW_FOUND )

#ENDIF(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
