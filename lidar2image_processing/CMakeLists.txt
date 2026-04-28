cmake_minimum_required(VERSION 3.10)
project(lidar2image_processing VERSION 1.0.1 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set( CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/../cmake/Modules/ ${CMAKE_MODULE_PATH} )
set( CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/)
if(NOT DEFINED ENV{PREFIX})
    set(CONDA_PREFIX $ENV{ISISROOT})
else()
    set(CONDA_PREFIX $ENV{PREFIX})
endif()

if( CMAKE_COMPILER_IS_GNUCXX )
	set( WARN_FLAG_ADDITIONS "-fno-strict-aliasing" )
endif( CMAKE_COMPILER_IS_GNUCXX )

include(SetupInstall)
include(SetupCTest)
include(SetupRPATH)

find_package(Boost 1.5 REQUIRED COMPONENTS system filesystem)
find_package(Cholmod REQUIRED)
find_package(Eigen3)
find_package(GDAL)
find_package(ISIS REQUIRED COMPONENTS isis)
find_package(Qt5 COMPONENTS
                Core
                Concurrent
                Gui
                Multimedia
                MultimediaWidgets
                Network
                OpenGL 
                PrintSupport
                Qml
                Quick
                Script
                ScriptTools
                Sql
                Svg
                Test
                WebChannel
                Widgets
                Xml
                XmlPatterns REQUIRED)
find_package(Qwt)

find_path(CSPICE_INCLUDE_DIR SpiceUsr.h
    PATHS 
    "${CONDA_PREFIX}/include/naif"
    "${CONDA_PREFIX}/include/cspice"
)

list( APPEND CMAKE_INSTALL_RPATH ${CONDA_PREFIX}/lib )

if (APPLE)   
    set(COMPILE_FLAGS  "-std=c17 -pedantic -F${CONDA_PREFIX} -Wall -ansi -mmacosx-version-min=10.6 -DISIS_LITTLE_ENDIAN=1 -D_LARGEFILE64_SOURCE -D_LARGEFILE_SOURCE -fPIC -DGMM_USES_SUPERLU -O1 -DQT_NO_DEBUG -DQT_GUI_LIB -DQT_CORE_LIB -DQT_SHARED -DENABLEJP2K=0")
    set(LINK_FLAGS "-L${CONDA_PREFIX}/lib -F${CONDA_PREFIX}/ -headerpad_max_install_names -arch x86_64 -Xarch_x86_64 -mmacosx-version-min=10.6 -bind_at_load -Wl,-rpath,${CONDA_PREFIX} -Wl,-w -Wl,-rpath,@loader_path/..  -lm -framework ApplicationServices  ${CHOLMODLIB_LIBLIST} -framework Accelerate ")
else()
    set(COMPILE_FLAGS " -Wall -ansi -DISIS_LITTLE_ENDIAN=1 -fPIC -DGMM_USES_SUPERLU -O1 -DENABLEJP2K=0")
    set(LINK_FLAGS "-L${CONDA_PREFIX}/lib -Wl,-E -Xlinker -z -Xlinker origin -Xlinker -rpath -Xlinker ${CONDA_PREFIX}/lib -pthread -lX11")
endif()
    
## We must add the in-source and out-of-source
## include paths 
##--------------------------------------------
 include_directories(
  ${CMAKE_SOURCE_DIR}
  ${PROJECT_SOURCE_DIR}/..
  ${EIGEN3_INCLUDE_DIR}
  ${ISIS_INCLUDE_DIRS}
  ${CSPICE_INCLUDE_DIR}           
  ${CONDA_PREFIX}/include
  ${CONDA_PREFIX}/include/qt
)
 setup_rpath()
 setup_ctest()


 #lidar_image_align
 add_executable(lidar_image_align
                lidar_image_align.cpp
                LidarImageAlign.cpp
                ../common/StringUtils.cpp
                ../common/FileListUtils.cpp
                ../common/ImageProcessing.cpp
                ../common/ReflectanceModels.cpp
                ../common/tiling.cpp
                ../camera_models/PVLRead.cpp
                ../camera_models/IsisInterfaceATK.cpp
                ../geotif_processing/CoordTransform.cpp
                ../geotif_processing/GeoUtils.cpp
                ../lidar_processing/TracksGCP.cpp
                ../lidar_processing/TracksLOLA.cpp
                ../lidar_processing/TracksFeatures.cpp)
  
target_compile_features(lidar_image_align PRIVATE cxx_std_17)

set_target_properties(lidar_image_align PROPERTIES
                       LINK_FLAGS ${LINK_FLAGS}
                     )

target_link_directories(lidar_image_align 
                        PRIVATE 
                        ${CONDA_PREFIX}/lib
                        )

target_link_libraries(lidar_image_align
                       ${GDAL_LIBRARY}
                       ${Boost_LIBRARIES}
                       ${ISIS_LIBRARIES}
                       isis
                       cspice
                       ${CHOLMOD_LIBRARY}
                       ${CAMD_LIBRARY}
                       ${AMD_LIBRARY}
                       ${COLAMD_LIBRARY}
                       ${BLAS_LIBRARY}
                       Qt5::XmlPatterns
                       Qt5::Xml
                       Qt5::Network
                       Qt5::Sql
                       Qt5::Gui
                       Qt5::Core
                       Qt5::Svg
                       Qt5::Test
                       Qt5::OpenGL
                       ${QWT_LIBRARY}
                       )

install_targets(/bin lidar_image_align)

if (CTEST_ENABLED)
  add_test(lidar_image_align_test bash ./lidar_image_align_test.sh)
endif(CTEST_ENABLED)
