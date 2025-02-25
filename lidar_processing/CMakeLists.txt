project(lidar_tracks)
cmake_minimum_required(VERSION 2.6)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set( CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/../cmake/Modules/ ${CMAKE_MODULE_PATH} )
set( CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/)

include(SetupInstall)
include(SetupCTest)
include(SetupRPATH)

#find_package( OpenCV REQUIRED )
find_package( Eigen3 )
find_package( GDAL )

if( CMAKE_COMPILER_IS_GNUCXX )
	set( WARN_FLAG_ADDITIONS "-fno-strict-aliasing" )
endif( CMAKE_COMPILER_IS_GNUCXX )


## We must add the in-source and out-of-source
## include paths 
##--------------------------------------------
include_directories(
  ${PROJECT_SOURCE_DIR}/..
  ${EIGEN3_INCLUDE_DIR}
 #${OpenCV_INCLUDE_DIR}
  ${GDAL_INCLUDE_DIR}
)

setup_rpath()
setup_ctest()

add_executable(track_splitter track_splitter.cpp
               TracksLOLA.cpp
               ../common/ReflectanceModels.cpp
               ../common/FileListUtils.cpp ../common/StringUtils.cpp
               ../geotif_processing/CoordTransform.cpp)
#target_link_libraries(track_splitter ${OpenCV_LIBS} ${GDAL_LIBRARY})
target_link_libraries(track_splitter ${GDAL_LIBRARY})

install_targets(/bin track_splitter)

if (CTEST_ENABLED)
  add_test(track_splitter_test bash ./track_splitter_test.sh)
endif(CTEST_ENABLED)
