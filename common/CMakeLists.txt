project(COMMON_MODULE)
cmake_minimum_required(VERSION 2.6 FATAL_ERROR)

set( CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/tests )
set( CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/../cmake/Modules/ ${CMAKE_MODULE_PATH} )
set(CMAKE_INSTALL_PREFIX ../../Install_dir CACHE INTERNAL "" )

include(SetupInstall)
include(SetupRPATH)

find_package( OpenCV REQUIRED )
find_package(PCL 1.6 REQUIRED)

include_directories(${OPENCV_INCLUDE_DIRS})
include_directories(${PCL_INCLUDE_DIRS})
link_directories(${PCL_LIBRARY_DIRS})
add_definitions(${PCL_DEFINITIONS})

set (LABLIB3_INCLUDE_PATH ${PROJECT_SOURCE_DIR}/../thirdparty/PDS-tools-v49/include)
set (LABLIB3_LIB_PATH ${PROJECT_SOURCE_DIR}/../thirdparty/PDS-tools-v49/lib)
set (LABLIB3_LIB lablib3)
include_directories(${OPENCV_INCLUDE_DIRS} ${LABLIB3_INCLUDE_PATH})
link_directories( ${LABLIB3_LIB_PATH} )

set( DIR_NAME "common" )
set( HEADERS
  GeometricTransf.h
  PdsToOpencv.h
  pc.h
  PdsRead.h
  StringUtils.h
  Tiling.h
  FileListUtils.h
  NotAssignable.h
  NotCopyable.h
  RotMat3d.h
  VectorUtil.h
  ImageProcessing.h
)

## Set up RPATH after all find_package calls
##-----------------------------------------------
setup_rpath()

project(tile_maker)
add_library(tile_maker SHARED tiling.cc FileListUtils.cc StringUtils.cc)

add_executable(image_tiling image_tiling.cc)
target_link_libraries(image_tiling tile_maker ${OpenCV_LIBS})

add_executable(pds_read pds_read.cc PdsToOpencv.cc PdsRead.cc)
target_link_libraries(pds_read ${OpenCV_LIBS} ${LABLIB3_LIB})

add_executable(gaussian_blur gaussian_blur.cc ImageProcessing.cc)
target_link_libraries(gaussian_blur ${OpenCV_LIBS})

#install libraries
install_library(tile_maker)

#install executables
install_targets(/bin image_tiling)

#install header files
install_files(/include/${DIR_NAME} FILES ${HEADERS})
