cmake_minimum_required(VERSION 3.10)
project(GEOTIF_PROCESSING_MODULE LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(NOT APPLE)
    add_definitions(-D_GLIBCXX_USE_CXX11_ABI=0)
endif()

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -Wall")

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/../cmake/Modules ${CMAKE_MODULE_PATH})
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/tests)

if(CMAKE_COMPILER_IS_GNUCXX)
    set(WARN_FLAG_ADDITIONS "-fno-strict-aliasing")
endif()

include(SetupInstall)
include(SetupRPATH)

find_package(GDAL REQUIRED)
find_package(Eigen3 REQUIRED)
find_package(PCL REQUIRED)

## include paths 
##--------------------------------------------
include_directories(
    ${PROJECT_SOURCE_DIR}/..
    ${GDAL_INCLUDE_DIR}
    ${EIGEN3_INCLUDE_DIR}
    ${PCL_INCLUDE_DIRS}
)

## Set up RPATH after all find_package calls
##-----------------------------------------------
setup_rpath()

add_executable(geo_assembler 
    geo_assembler.cc
    CoordTransform.cc 
    GeoAssembler.cc
    GeoTiling.cc 
    GeoUtils.cc
    ../common/StringUtils.cc 
    ../common/FileListUtils.cc
    ../common/tiling.cc 
    ../common/ImageProcessing.cc
    ../common/FileSystemUtils.cc
)
target_link_libraries(geo_assembler ${GDAL_LIBRARY})

add_executable(geo_pyramid 
    geo_pyramid.cc 
    GeoTiling.cc 
    GeoUtils.cc 
    CoordTransform.cc 
    ../common/StringUtils.cc 
    ../common/tiling.cc
)
target_link_libraries(geo_pyramid ${GDAL_LIBRARY})

add_executable(geo_pyramid_assembler 
    geo_pyramid_assembler.cc
    GeoAssembler.cc 
    GeoPyramidAssembler.cc 
    CoordTransform.cc
    GeoTiling.cc 
    GeoUtils.cc
    ../common/StringUtils.cc 
    ../common/FileListUtils.cc
    ../common/FileSystemUtils.cc 
    ../common/tiling.cc
    ../common/ImageProcessing.cc
)
target_link_libraries(geo_pyramid_assembler ${PCL_LIBRARIES} ${GDAL_LIBRARY})

add_executable(geo_diff geo_diff.cc ../common/StringUtils.cc)
target_link_libraries(geo_diff ${GDAL_LIBRARY})

add_executable(geo_translate 
    geo_translate.cc
    GeoUtils.cc 
    CoordTransform.cc
    ../common/tiling.cc
)
target_link_libraries(geo_translate ${GDAL_LIBRARY})

add_executable(geo_resize 
    geo_resize.cc 
    GeoUtils.cc  
    CoordTransform.cc
    ../common/StringUtils.cc 
    ../common/tiling.cc
)
target_link_libraries(geo_resize ${GDAL_LIBRARY})

install(TARGETS 
    geo_assembler 
    geo_pyramid 
    geo_pyramid_assembler 
    geo_diff 
    geo_translate 
    geo_resize 
    DESTINATION bin
)