cmake_minimum_required(VERSION 3.10)
project(lidar2image_processing VERSION 1.0.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set( CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/../cmake/Modules/ ${CMAKE_MODULE_PATH} )
set( CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/tests )

if( CMAKE_COMPILER_IS_GNUCXX )
	set( WARN_FLAG_ADDITIONS "-fno-strict-aliasing" )
endif( CMAKE_COMPILER_IS_GNUCXX )

include(SetupInstall)
include(SetupCTest)
include(SetupRPATH)

set(BOOST_ROOT /usr/local)
find_package(Boost 1.5 REQUIRED COMPONENTS system filesystem)
find_package(Cholmod REQUIRED)
find_package(Eigen3)
#find_package(OpenCV)
find_package(GDAL)
find_package(ISIS REQUIRED COMPONENTS isis)
find_package(Qt5 COMPONENTS
                Core
                Concurrent
                Gui
                Multimedia
                MultimediaWidgets
                Network
                OpenGL # Needed to install mesa-common-dev for this!
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
                XmlPatterns
                # Search this path explicitly for MacOS OpenGL Framework
                PATHS /System/Library/Frameworks/ REQUIRED)

find_package(Qwt)

list( APPEND CMAKE_INSTALL_RPATH $ENV{ISISROOT}/lib )

if (APPLE)   
    #set(COMPILE_FLAGS  "-std=c17 -pedantic -F${ISIS_3RD} -Wall -ansi -arch x86_64 -Xarch_x86_64 -mmacosx-version-min=10.6 -DISIS_LITTLE_ENDIAN=1 -D_LARGEFILE64_SOURCE -D_LARGEFILE_SOURCE -fPIC -DGMM_USES_SUPERLU -O1 -DQT_NO_DEBUG -DQT_GUI_LIB -DQT_CORE_LIB -DQT_SHARED -DQT_GUI_LIB -DQT_CORE_LIB -DQT_SHARED -DENABLEJP2K=0")
    set(COMPILE_FLAGS  "-std=c17 -pedantic -F${ISIS_ROOT} -Wall -ansi -mmacosx-version-min=10.6 -DISIS_LITTLE_ENDIAN=1 -D_LARGEFILE64_SOURCE -D_LARGEFILE_SOURCE -fPIC -DGMM_USES_SUPERLU -O1 -DQT_NO_DEBUG -DQT_GUI_LIB -DQT_CORE_LIB -DQT_SHARED -DQT_GUI_LIB -DQT_CORE_LIB -DQT_SHARED -DENABLEJP2K=0")
    #set(LINK_FLAGS "-L${ISIS_3RD} -F${ISIS_3RD} -headerpad_max_install_names -arch x86_64 -Xarch_x86_64 -mmacosx-version-min=10.6 -bind_at_load -Wl,-rpath,${ISISROOT} -Wl,-w -Wl,-rpath,@loader_path/..  -lm -framework ApplicationServices -framework QtXmlPatterns -framework QtXml -framework QtNetwork -framework QtSql -framework QtGui -framework QtCore -framework QtSvg -framework QtTest -framework QtWebKit -framework QtOpenGL -framework qwt  ${CHOLMODLIB_DEPS} -framework Accelerate ")
    set(LINK_FLAGS "-L$ENV{ISISROOT}/lib -F$ENV{ISISROOT}/ -headerpad_max_install_names -arch x86_64 -Xarch_x86_64 -mmacosx-version-min=10.6 -bind_at_load -Wl,-rpath,$ENV{ISISROOT} -Wl,-w -Wl,-rpath,@loader_path/..  -lm -framework ApplicationServices  ${CHOLMODLIB_LIBLIST} -framework Accelerate ")
else()
    set(COMPILE_FLAGS " -Wall -ansi -DISIS_LITTLE_ENDIAN=1 -fPIC -DGMM_USES_SUPERLU -O1 -DENABLEJP2K=0")
    set(LINK_FLAGS "-L$ENV{ISISROOT}/3rdParty/lib -Wl,-E -Xlinker -z -Xlinker origin -Xlinker -rpath -Xlinker $ENV{ISISROOT}/lib -pthread -lX11")
endif()
    
## We must add the in-source and out-of-source
## include paths 
##--------------------------------------------
 include_directories(
  ${CMAKE_SOURCE_DIR}
  ${PROJECT_SOURCE_DIR}/..
  ${EIGEN3_INCLUDE_DIR}
  #${OPENCV_INCLUDE_DIRS}
  ${PCL_INCLUDE_DIRS}
  ${ISIS_INCLUDE_DIR}
  ${CSPICE_ROOT}/naif
  ${CSPICE_ROOT}/cspice
  ${CSPICE_ROOT}/Qt/QtCore
  ${CSPICE_ROOT}/Qt/QtXml
  ${CSPICE_ROOT}/Qt/QtXmlPatterns
  ${CSPICE_ROOT}/Qt
  ${CSPICE_ROOT}/gdal
  ${CSPICE_ROOT}/
  ${CSPICE_ROOT}/Qt/QtCore.framework/Headers
  ${CSPICE_ROOT}/Qt/QtAssistant.framework/Headers
  ${CSPICE_ROOT}/Qt/QtGui.framework/Headers
  ${CSPICE_ROOT}/Qt/QtNetwork.framework/Headers
  ${CSPICE_ROOT}/Qt/QtOpenGL.framework/Headers
  ${CSPICE_ROOT}/Qt/QtScript.framework/Headers
  ${CSPICE_ROOT}/Qt/QtScriptTools.framework/Headers
  ${CSPICE_ROOT}/Qt/QtSql.framework/Headers
  ${CSPICE_ROOT}/Qt/QtSvg.framework/Headers
  ${CSPICE_ROOT}/Qt/QtTest.framework/Headers
  ${CSPICE_ROOT}/Qt/QtWebKit.framework/Headers
  ${CSPICE_ROOT}/Qt/QtXml.framework/Headers
  ${CSPICE_ROOT}/Qt/QtXmlPatterns.framework/Headers
  ${CSPICE_ROOT}/qwt.framework/Headers
  ${CSPICE_ROOT}/xercesc
  ${CSPICE_ROOT}/tiff
  ${CSPICE_ROOT}/tnt
  ${CSPICE_ROOT}/jama
  ${CSPICE_ROOT}/geos
  ${CSPICE_ROOT}/gsl
  ${CSPICE_ROOT}/google
  ${CSPICE_ROOT}/superlu
)

 ## Set up RPATH after all find_package calls
 ##-----------------------------------------------
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
                ../camera_models/PVLRead.cpp	../camera_models/IsisInterfaceATK.cpp
                ../geotif_processing/CoordTransform.cpp
                ../geotif_processing/GeoUtils.cpp
                ../lidar_processing/TracksGCP.cpp
                ../lidar_processing/TracksLOLA.cpp
                ../lidar_processing/TracksFeatures.cpp)
  
target_compile_features(lidar_image_align PRIVATE cxx_std_17)

set_target_properties(lidar_image_align PROPERTIES
                       #COMPILE_FLAGS ${COMPILE_FLAGS}
                       LINK_FLAGS ${LINK_FLAGS}
                     )

target_link_directories(lidar_image_align
                        PRIVATE
                       $ENV{ISISROOT}
                       )

target_link_libraries(lidar_image_align
                       ${GDAL_LIBRARY}
                       ${Boost_LIBRARIES}
                       ${ISIS_LIBRARIES}
                       ${lib_cspice}
                       ${lib_geos}
                       ${lib_gsl}
                       ${lib_gslcblas}
                       ${lib_protobuf}
                       ${lib_tiff}
                       ${CHOLMOD_LIBRARY}
                       ${CAMD_LIBRARY}
                       ${AMD_LIBRARY}
                       ${COLAMD_LIBRARY}
                       ${lib_superlu}
                       ${BLAS_LIBRARY}
                       ${lib_z}
                       ${lib_xerces-c}
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

# Here we use ISIS's own tool to fix the paths to dynamic libraries
if(NOT APPLE)
    MESSAGE("NOT APPLE")
    add_custom_command(
        TARGET lidar_image_align 
        POST_BUILD
        COMMAND ${PATCHELF} --set-rpath $ENV{ISISROOT}/lib lidar_image_align
        COMMENT "Editing the run-time paths"
        VERBATIM
    )
endif()

install_targets(/bin lidar_image_align)

if (CTEST_ENABLED)
  add_test(lidar_image_align_test bash ./lidar_image_align_test.sh)
endif(CTEST_ENABLED)
