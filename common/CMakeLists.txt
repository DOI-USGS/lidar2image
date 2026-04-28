cmake_minimum_required(VERSION 3.10 FATAL_ERROR)
project(COMMON_MODULE LANGUAGES CXX)

# Set C++ Standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(NOT APPLE)
    add_definitions(-D_GLIBCXX_USE_CXX11_ABI=0)
endif()

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/tests)
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/../cmake/Modules/ ${CMAKE_MODULE_PATH})


include(SetupInstall)
include(SetupRPATH)

find_package(OpenCV REQUIRED)
find_package(PCL 1.6 REQUIRED)

include_directories(
    ${OPENCV_INCLUDE_DIRS} 
    ${PCL_INCLUDE_DIRS}
    ${PROJECT_SOURCE_DIR}/../thirdparty/PDS-tools-v49/include
)
link_directories(${PCL_LIBRARY_DIRS})
link_directories(${PROJECT_SOURCE_DIR}/../thirdparty/PDS-tools-v49/lib)
add_definitions(${PCL_DEFINITIONS})

set(LABLIB3_LIB lablib3)
set(DIR_NAME "common")
set(HEADERS
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

add_library(tile_maker SHARED tiling.cc FileListUtils.cc StringUtils.cc)
target_link_libraries(tile_maker ${OpenCV_LIBS})

add_executable(image_tiling image_tiling.cc)
target_link_libraries(image_tiling tile_maker ${OpenCV_LIBS})

add_executable(pds_read pds_read.cc PdsToOpencv.cc PdsRead.cc)
target_link_libraries(pds_read ${OpenCV_LIBS} ${LABLIB3_LIB})

add_executable(gaussian_blur gaussian_blur.cc ImageProcessing.cc)
target_link_libraries(gaussian_blur ${OpenCV_LIBS})

install(TARGETS tile_maker 
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib)

install(TARGETS image_tiling pds_read gaussian_blur 
        DESTINATION bin)

install(FILES ${HEADERS} 
        DESTINATION include/${DIR_NAME})