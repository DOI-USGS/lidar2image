cmake_minimum_required(VERSION 3.10)
project(lidar_tracks LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(NOT APPLE)
    add_definitions(-D_GLIBCXX_USE_CXX11_ABI=0)
endif()

set( CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/../cmake/Modules/ ${CMAKE_MODULE_PATH} )
set( CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/)

include(SetupInstall)
include(SetupCTest)
include(SetupRPATH)

find_package(Eigen3 REQUIRED)
find_package(GDAL REQUIRED)
find_package(ISIS REQUIRED COMPONENTS isis)

if( CMAKE_COMPILER_IS_GNUCXX )
	set( WARN_FLAG_ADDITIONS "-fno-strict-aliasing" )
endif( CMAKE_COMPILER_IS_GNUCXX )


## We must add the in-source and out-of-source
## include paths 
##--------------------------------------------
include_directories(
  ${ISIS_INCLUDE_DIRS}
  ${EIGEN3_INCLUDE_DIR}
  ${GDAL_INCLUDE_DIRS}
  ${PROJECT_SOURCE_DIR}/..
)

setup_rpath()
setup_ctest()

add_executable(track_splitter track_splitter.cpp
               TracksLOLA.cpp
               ../common/ReflectanceModels.cpp
    ../common/FileListUtils.cpp ../common/StringUtils.cpp
    ../geotif_processing/CoordTransform.cpp)
target_link_libraries(track_splitter 
    ${GDAL_LIBRARY}
    ${ISIS_LIBRARIES}
)

install(TARGETS track_splitter DESTINATION bin)

if(CTEST_ENABLED)
    add_test(NAME track_splitter_test COMMAND bash ./track_splitter_test.sh)
endif()