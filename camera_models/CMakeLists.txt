project(CAMERA_MODELS)
cmake_minimum_required(VERSION 2.6 FATAL_ERROR)

set( CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/tests )
set( CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/../cmake/Modules/ ${CMAKE_MODULE_PATH} )

include(SetupInstall)
include(SetupRPATH)

find_package(Eigen3)

include_directories(${EIGEN3_INCLUDE_DIR} ../common)

setup_rpath()

#set(CMAKE_VERBOSE_MAKEFILE ON)
add_executable (cahv_to_pinhole cahv_to_pinhole.cc 
		CAHV.cc  pinhole.cc
		../common/GeometricTransf.cc
		../common/StringUtils.cc)
target_link_libraries (cahv_to_pinhole)

add_executable (pinhole_to_cahv pinhole_to_cahv.cc 
		CAHV.cc  pinhole.cc
		../common/GeometricTransf.cc
		../common/StringUtils.cc)
target_link_libraries (pinhole_to_cahv)





