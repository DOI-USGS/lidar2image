cmake_minimum_required(VERSION 3.10 FATAL_ERROR)
project(CAMERA_MODELS LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(NOT APPLE)
    add_definitions(-D_GLIBCXX_USE_CXX11_ABI=0)
endif()

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/tests)
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/../cmake/Modules/ ${CMAKE_MODULE_PATH})

include(SetupInstall)
include(SetupRPATH)

find_package(Eigen3 REQUIRED)

include_directories(
    ${EIGEN3_INCLUDE_DIR} 
    ${PROJECT_SOURCE_DIR}/../common
)

setup_rpath()

add_executable(cahv_to_pinhole 
    cahv_to_pinhole.cc 
    CAHV.cc 
    pinhole.cc
    ../common/GeometricTransf.cc
    ../common/StringUtils.cc
)
target_link_libraries(cahv_to_pinhole Eigen3::Eigen)

add_executable(pinhole_to_cahv 
    pinhole_to_cahv.cc 
    CAHV.cc 
    pinhole.cc
    ../common/GeometricTransf.cc
    ../common/StringUtils.cc
)
target_link_libraries(pinhole_to_cahv Eigen3::Eigen)

install(TARGETS cahv_to_pinhole pinhole_to_cahv DESTINATION bin)