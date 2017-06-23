# Copyright (c) 2012-2017, Association Scientifique pour la Geologie et ses Applications (ASGA)
# All rights reserved.
# 
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#     * Redistributions of source code must retain the above copyright
#       notice, this list of conditions and the following disclaimer.
#     * Redistributions in binary form must reproduce the above copyright
#       notice, this list of conditions and the following disclaimer in the
#       documentation and/or other materials provided with the distribution.
#     * Neither the name of ASGA nor the
#       names of its contributors may be used to endorse or promote products
#       derived from this software without specific prior written permission.
# 
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
# ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL ASGA BE LIABLE FOR ANY
# DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
# (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
# ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
# 
#     http://www.ring-team.org
# 
#     RING Project
#     Ecole Nationale Superieure de Geologie - GeoRessources
#     2 Rue du Doyen Marcel Roubault - TSA 70605
#     54518 VANDOEUVRE-LES-NANCY
#     FRANCE

# Root CMakeList of RINGMesh project
cmake_minimum_required(VERSION 2.8)

#------------------------------------------------------------------------------------------------
# Options of the RINGMesh project 
# They can be set using CMake GUI, or in command line: cmake -D<var>:<type>=<value> MyProjectFolder

# Configuration modes - Debug and Release by default + RelWithDebInfo for Windows.
# RelWithDebInfo is mandatory to debug some plugins in Windows in debug mode (e.g., a Gocad plugin) to
# avoid issues in particular in the STL (e.g., bad transfer of std::string from Gocad to RINGMesh).
set(CMAKE_CONFIGURATION_TYPES Debug Release RelWithDebInfo CACHE CMAKE_CONFIGURATION_TYPES "Configuration types")

# If a path is specified RINGMesh will be compiled to work with - empty by default
set(MG_TETRA "" CACHE PATH "Path to MG-Tetra")
if(MG_TETRA)
    set(USE_MG_TETRA ON)
    message(STATUS "Using MG-Tetra directory =  ${MG_TETRA}")
endif(MG_TETRA)
  

# By default RINGMesh will be compiled with OpenMP flags
option(USE_OPENMP "Use Open MP" ON)

# By default RINGMesh is compiled with Tetgen  
# Check Tetgen licences to ensure you have the right to use it freely
option(RINGMESH_WITH_TETGEN "Use Tetgen tetrahedral mesher" ON)

# Optional components of RINGMesh - creation of specific targets
option(RINGMESH_WITH_GRAPHICS "Compile viewer" OFF)
include(CMakeDependentOption)
# Option to build ringmesh-view. Disable when RINGMESH_WITH_GRAPHICS
# is OFF. The fist time RINGMESH_WITH_GRAPHICS is ON, this
# option is ON.
CMAKE_DEPENDENT_OPTION(BUILD_RINGMESH_VIEW "Build ringmesh-view" ON
    "RINGMESH_WITH_GRAPHICS" OFF)

option(RINGMESH_WITH_UTILITIES "Compile utility executables" OFF)
option(RINGMESH_WITH_TUTORIALS "Compile API trainings and tuturials" OFF)
option(RINGMESH_WITH_TESTS "Compile test projects" OFF)
option(RINGMESH_WITH_GEOLOGYJS "Enable export to html" OFF)

# Option to test ringmesh-view. Disable when BUILD_RINGMESH_VIEW
# is OFF. The fist time BUILD_RINGMESH_VIEW is ON, this
# option is ON.
CMAKE_DEPENDENT_OPTION(RINGMESH_TEST_GRAPHICS "Test ringmesh-view" ON
    "BUILD_RINGMESH_VIEW;RINGMESH_WITH_TESTS" OFF)

# Optional custom steps 
option(BUILD_DOCUMENTATION "Create and install the HTML documentation (requires Doxygen)")

# Optional custom configuration of data directories
if(RINGMESH_WITH_TESTS)
    set(TEST_DIRECTORY ${CMAKE_SOURCE_DIR}/tests)
    set(TEST_SOURCE_DIRECTORY ${TEST_DIRECTORY}/src)
    set(TEST_DATA_INPUT_DIRECTORY "${TEST_DIRECTORY}/data")
    # PROJECT_BINARY_DIR not defined at this point
    set(TEST_DATA_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/build/tests/output" CACHE PATH "Output data location")
    if(UNIX)    
      list(APPEND CMAKE_CONFIGURATION_TYPES Coverage)
    endif()    
endif()

# Optional custom configuration of tutorials directories
if(RINGMESH_WITH_TUTORIALS)
    set(TUTORIALS_DIRECTORY ${CMAKE_SOURCE_DIR}/doc/tutorials)
    set(TUTORIALS_SOURCE_DIRECTORY ${TUTORIALS_DIRECTORY}/src)
    set(TUTORIALS_DATA_INPUT_DIRECTORY "${TUTORIALS_DIRECTORY}/data")
    # PROJECT_BINARY_DIR not defined at this point
    set(TUTORIALS_DATA_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/build/tutorials/output" CACHE PATH "Output data location")
    if(UNIX)    
      list(APPEND CMAKE_CONFIGURATION_TYPES Coverage)
    endif() 
endif()

if(RINGMESH_WITH_GEOLOGYJS)
    set(GEOLOGY_JS_DIRECTORY "${CMAKE_SOURCE_DIR}/../geology.js/cpp" CACHE PATH "geology.js cpp directory location")
    message(STATUS ${GEOLOGY_JS_DIRECTORY})
endif(RINGMESH_WITH_GEOLOGYJS)

# Turn on the ability to create folders to organize projects and files
# It creates "CMakePredefinedTargets" folder by default and adds CMake
# defined projects like INSTALL.vcproj and ZERO_CHECK.vcproj
set_property(GLOBAL PROPERTY USE_FOLDERS ON)

#------------------------------------------------------------------------------------------------
# Additional cmake modules
include(ExternalProject)
# build export macros
include (GenerateExportHeader) 

# Enable testing with CTest
enable_testing()

#------------------------------------------------------------------------------------------------
# Define the project RINGMesh      
project(RINGMesh)

# Define version number 
# It is then exported to the configuration file
set (RINGMesh_VERSION_MAJOR 4)
set (RINGMesh_VERSION_MINOR 0)
set (RINGMesh_VERSION_PATCH 0)
set (RINGMesh_VERSION ${RINGMesh_VERSION_MAJOR}.${RINGMesh_VERSION_MINOR}.${RINGMesh_VERSION_PATCH})


#------------------------------------------------------------------------------------------------
# Generate configuration directories for single-configuration generators (Make)
# and run cmake configuration command in each one of them
if(UNIX)    
    # If the CMAKE_BUILD_TYPE is no defined
    # i.e. at first run of that file
    if(CMAKE_BUILD_TYPE STREQUAL "")
    
      # For all configuration types  
      foreach(config IN LISTS CMAKE_CONFIGURATION_TYPES)
          # The binary directory for this configuration
          set(project_binary_dir_config  ${PROJECT_BINARY_DIR}/ringmesh/${config})
    
          # Create the directory
          file(MAKE_DIRECTORY ${project_binary_dir_config})
          
          # Launch cmake for this configuration
          # by specifying the CMAKE_BUILD_TYPE and the PROJECT _BINARY_DIR, 
          # and pass on all the options
          execute_process(
             COMMAND ${CMAKE_COMMAND} ${CMAKE_SOURCE_DIR}
                -DCMAKE_BUILD_TYPE:STRING=${config}
                -DPROJECT_BINARY_DIR:STRING=${project_binary_dir_config}
                -DGLOBAL_BINARY_DIR:STRING=${PROJECT_BINARY_DIR}
                -DMG_TETRA:STRING=${MG_TETRA}
                -DUSE_OPENMP:BOOL=${USE_OPENMP}
                -DRINGMESH_WITH_TETGEN:BOOL=${RINGMESH_WITH_TETGEN} 
                -DRINGMESH_WITH_GRAPHICS:BOOL=${RINGMESH_WITH_GRAPHICS}
                -DRINGMESH_TEST_GRAPHICS:BOOL=${RINGMESH_TEST_GRAPHICS}
                -DBUILD_RINGMESH_VIEW:BOOL=${BUILD_RINGMESH_VIEW}
                -DRINGMESH_WITH_UTILITIES:BOOL=${RINGMESH_WITH_UTILITIES}
                -DRINGMESH_WITH_TESTS:BOOL=${RINGMESH_WITH_TESTS}
                -DTEST_DATA_INPUT_DIRECTORY:STRING=${TEST_DATA_INPUT_DIRECTORY}
                -DTEST_DATA_OUTPUT_DIRECTORY:STRING=${TEST_DATA_OUTPUT_DIRECTORY}        
                -DRINGMESH_WITH_TUTORIALS:BOOL=${RINGMESH_WITH_TUTORIALS}
                -DTUTORIALS_DATA_INPUT_DIRECTORY:STRING=${TUTORIALS_DATA_INPUT_DIRECTORY}
                -DTUTORIALS_DATA_OUTPUT_DIRECTORY:STRING=${TUTORIALS_DATA_OUTPUT_DIRECTORY}  
                -DRINGMESH_WITH_GEOLOGYJS:BOOL=${RINGMESH_WITH_GEOLOGYJS}
                -DGEOLOGY_JS_DIRECTORY:STRING=${GEOLOGY_JS_DIRECTORY} 
                -DBUILD_DOCUMENTATION:BOOL=${BUILD_DOCUMENTATION}                     
                WORKING_DIRECTORY ${project_binary_dir_config})                
      endforeach()
      
      # Get out when all configurations have been run
      # We do not want to execute the rest of the file in that case
      return()
    endif()
else(UNIX)
    set(GLOBAL_BINARY_DIR  ${PROJECT_BINARY_DIR}/..)
endif(UNIX)

message(STATUS "RINGMesh binary directory is: ${PROJECT_BINARY_DIR}") 
message(STATUS "RINGMesh source directory is: ${PROJECT_SOURCE_DIR}") 

include(cmake/utils.cmake)

#------------------------------------------------------------------------------------------------
# Platform dependent settings
if(WIN32)
    # Select the CMAKE_CONFIGURATION_TYPES from Visual Studio
    set(COMPILATION_OPTION --config ${CMAKE_CFG_INTDIR})
elseif(UNIX)
    set(CMAKE_POSITION_INDEPENDENT_CODE TRUE)
    add_compile_options(-Wall -Wextra -Wno-long-long -Wconversion)
    
    # Determine gcc version and activate additional warnings available in latest versions
    execute_process(COMMAND ${CMAKE_CXX_COMPILER} -dumpversion OUTPUT_VARIABLE G++_VERSION)

    if(G++_VERSION VERSION_LESS 4.8 )
        message(FATAL_ERROR "RINGMesh require G++ version >= 4.8")
    endif()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
    add_compile_options(-Wsign-conversion)
    add_compile_options(-Wdouble-promotion)
    add_compile_options(-Wno-attributes)

     if(CMAKE_BUILD_TYPE STREQUAL "Coverage")
       set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${PROJECT_SOURCE_DIR}/cmake)
       set(COVERALLS ON)
       set(COVERALLS_UPLOAD ON)
       include(Coveralls)
       coveralls_turn_on_coverage()
    endif()
    
    # Enable parallel computation
    set(cpuinfo_file "/proc/cpuinfo")
    if(EXISTS "${cpuinfo_file}")
        file(STRINGS "${cpuinfo_file}" procs REGEX "^processor.: [0-9]+$")
        list(LENGTH procs PROCESSOR_COUNT)
    endif()
    set(COMPILATION_OPTION -- -j${PROCESSOR_COUNT})
else(WIN32)
    # Nothing else really supported
    message(FATAL_ERROR "Your platform is not supported. Please modify the cmakelist file.")    
endif(WIN32)

#------------------------------------------------------------------------------------------------
# Configure a header file to pass on some of the CMake settings
# to the source code
set(generated_configure_file ${PROJECT_BINARY_DIR}/ringmesh/ringmesh_config.h)
set(input_configure_file ${PROJECT_SOURCE_DIR}/include/ringmesh/basic/ringmesh_config.h.in)
configure_file(${input_configure_file} ${generated_configure_file})

# We want to be able to include this file
include_directories(${PROJECT_BINARY_DIR})

#------------------------------------------------------------------------------------------------
# GEOGRAM 
# Set the path to Geogram code
set(GEOGRAM_PATH ${PROJECT_SOURCE_DIR}/third_party/geogram)

# Geogram platform dependent settings
if(WIN32)
    set(GEOGRAM_PATH_BIN ${GLOBAL_BINARY_DIR}/geogram)
    set(geoplatform Win-vs-dynamic-generic)  
    # Extra lib 
    set (EXTRA_LIBS ${EXTRA_LIBS} psapi)
    # TODO check that it is really necessary [JP]
    add_compile_options(-DGEO_DYNAMIC_LIBS) 
else(WIN32)
    set(GEOGRAM_PATH_BIN ${GLOBAL_BINARY_DIR}/geogram/${CMAKE_BUILD_TYPE})
    set(geoplatform Linux64-gcc-dynamic)
endif(WIN32)

# Define Geogram as an external project that we know how to
# configure and compile
ExternalProject_Add(geogram_ext
  PREFIX ${GEOGRAM_PATH_BIN}
  
  #--Download step--------------
  DOWNLOAD_COMMAND ""
  
  #--Update/Patch step----------
  UPDATE_COMMAND ""

  #--Configure step-------------
  SOURCE_DIR ${GEOGRAM_PATH}
  CONFIGURE_COMMAND ${CMAKE_COMMAND} ${GEOGRAM_PATH}
        -G ${CMAKE_GENERATOR} 
        -DVORPALINE_PLATFORM:STRING=${geoplatform}
        -DGEOGRAM_WITH_TETGEN:BOOL=${RINGMESH_WITH_TETGEN} 
        -DGEOGRAM_WITH_GRAPHICS:BOOL=${RINGMESH_WITH_GRAPHICS}
        -DGEOGRAM_WITH_EXPLORAGRAM:BOOL=OFF
        -DGEOGRAM_WITH_LUA:BOOL=OFF
        -DCMAKE_BUILD_TYPE:STRING=${CMAKE_BUILD_TYPE}
  
  #--Build step-----------------
  BINARY_DIR ${GEOGRAM_PATH_BIN}
  #-- Command to build geogram
  BUILD_COMMAND ${CMAKE_COMMAND} --build ${GEOGRAM_PATH_BIN} ${COMPILATION_OPTION}

  #--Install step---------------
  INSTALL_COMMAND "" 
)

ExternalProject_Add_Step(geogram_ext forcebuild
    DEPENDERS build
    ALWAYS 1
  )

# Add geogram include directories to the current ones
include_directories(SYSTEM ${GEOGRAM_PATH}/src/lib)

# Add geogram project libs to the libs with which RINGMesh will link
set(EXTRA_LIBS ${EXTRA_LIBS} geogram)
if(RINGMESH_WITH_GRAPHICS)
    include(${GEOGRAM_PATH}/cmake/opengl.cmake)
    set(EXTRA_LIBS ${EXTRA_LIBS} geogram_gfx ${OPENGL_LIBRARIES})
endif(RINGMESH_WITH_GRAPHICS)
    
# Add geogram bin directories to the current ones 
# It would be preferable to set the imported library location [JP]
link_directories(${GEOGRAM_PATH_BIN}/lib)

#------------------------------------------------------------------------------------------------    
# Get MG-Tetra if a path is given
if(MG_TETRA)   
  include_directories(SYSTEM ${MG_TETRA}/include)
 
  if(CMAKE_SIZEOF_VOID_P EQUAL 8)
    set(PLATFORM _64)
  else(CMAKE_SIZEOF_VOID_P EQUAL 8)
    SET(PLATFORM)
  endif(CMAKE_SIZEOF_VOID_P EQUAL 8)
  
  # Is there not a nice way to import this library? [JP]
  if(WIN32)
    link_directories(${MG_TETRA}/lib/WIN7${PLATFORM}_VC10)
  else(WIN32)
    link_directories(${MG_TETRA}/lib/Linux${PLATFORM})
  endif(WIN32)
   
  set (EXTRA_LIBS ${EXTRA_LIBS} meshgems mg-tetra)
endif()

#------------------------------------------------------------------------------------------------
# Add compiler options for OpenMP
if(USE_OPENMP)
    find_package(OpenMP)
    if (OPENMP_FOUND)
        set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
        set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    endif()
endif(USE_OPENMP)

#--------------------------------------------------------
# geology.js library
if(RINGMESH_WITH_GEOLOGYJS)
    message(STATUS "Using geology.js directory = ${GEOLOGY_JS_DIRECTORY}")
    # geology.js platform dependent settings. CMAKE_BUILD_TYPE is empty on Windows.
    set(GEOLOGY_JS_PATH_BIN ${GEOLOGY_JS_DIRECTORY}/build/${CMAKE_BUILD_TYPE})
    # Add geology.js include directories to the current ones
    include_directories(SYSTEM ${GEOLOGY_JS_DIRECTORY}/include)

    # Add geology.js project libs to the libs with which RINGMesh will link
    set(EXTRA_LIBS ${EXTRA_LIBS} GEOLOGYJS)
    # Add geology.js bin directories to the current ones 
    # It would be preferable to set the imported library location [JP]
    link_directories(${GEOLOGY_JS_PATH_BIN}/lib)
endif(RINGMESH_WITH_GEOLOGYJS)

#------------------------------------------------------------------------------------------------
# Collect the library files
FILE(GLOB directories RELATIVE ${PROJECT_SOURCE_DIR}/src/ringmesh ${PROJECT_SOURCE_DIR}/src/ringmesh/*)
FOREACH(directory ${directories})
    IF(IS_DIRECTORY ${PROJECT_SOURCE_DIR}/src/ringmesh/${directory})
        source_file_directory(ringmesh_src ${directory})
        include_file_directory(ringmesh_include ${directory})
    ENDIF()
ENDFOREACH()

file(GLOB zlib_src "${PROJECT_SOURCE_DIR}/third_party/zlib/*.c")
file(GLOB zlib_include "${PROJECT_SOURCE_DIR}/third_party/zlib/*.h")
# Add configure file to ringmesh includes
set(ringmesh_include ${ringmesh_include} ${generated_configure_file})

# Group files of a project in directories
# All the files should be added to the sources of the project
# and the global property USE_FOLDERS must be ON
source_group("Zlib" FILES ${zlib_src} ${zlib_include})

set(tinyxml2_src "${PROJECT_SOURCE_DIR}/third_party/tinyxml2/tinyxml2.cpp")
set(tinyxml2_include "${PROJECT_SOURCE_DIR}/third_party/tinyxml2/tinyxml2.h")

# Group files of a project in directories
# All the files should be added to the sources of the project
# and the global property USE_FOLDERS must be ON
source_group("tinyxml2" FILES ${tinyxml2_src} ${tinyxml2_include})

#------------------------------------------------------------------------------------------------
# Build configuration
# Organize outputs in a Linux way
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)

# Add include directories of RINGMesh
include_directories(${PROJECT_SOURCE_DIR}/include)
include_directories(SYSTEM ${PROJECT_SOURCE_DIR}/third_party)

#add generated export file to include directories
set(generated_export_file ${PROJECT_BINARY_DIR}/ringmesh/ringmesh_export.h)
set(ringmesh_include ${ringmesh_include} ${generated_export_file})

# Add the RINGMesh target as a shared library
add_library(RINGMesh SHARED
    ${ringmesh_include} 
    ${ringmesh_src} 
    ${zlib_src}
    ${zlib_include}
    ${tinyxml2_src}
    ${tinyxml2_include})
    
# Function for generation of export macros for libraries (dll)
generate_export_header(RINGMesh EXPORT_MACRO_NAME RINGMESH_API EXPORT_FILE_NAME ${generated_export_file} )

# RINGMesh depends on Geogram
add_dependencies(RINGMesh geogram_ext)

# Libraries to which RINGMesh should link
target_link_libraries(RINGMesh ${EXTRA_LIBS})

#------------------------------------------------------------------------------------------------
# Optional modules configuration
if(BUILD_RINGMESH_VIEW)
    message(STATUS "Configure ringmesh-view")
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)
    message(STATUS "ringmesh-view directory ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}")
    include(${CMAKE_SOURCE_DIR}/third_party/geogram/cmake/opengl.cmake)
    add_ringmesh_binary(ringmesh-view geogram_gfx ${OPENGL_LIBRARIES})
    copy_for_windows()
endif(BUILD_RINGMESH_VIEW)

if(RINGMESH_WITH_UTILITIES)
    message(STATUS "Configuring RINGMesh with utilities")
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin/utilities)
    message(STATUS "RINGMesh utilities directory ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}")
    add_ringmesh_binary(ringmesh-convert)
    add_ringmesh_binary(ringmesh-surface-convert)
    add_ringmesh_binary(ringmesh-stats)
    add_ringmesh_binary(ringmesh-repair)
    add_ringmesh_binary(ringmesh-mesh-quality)
    add_ringmesh_binary(ringmesh-rotate)
    add_ringmesh_binary(ringmesh-translate)
    copy_for_windows()
endif(RINGMESH_WITH_UTILITIES)

if(RINGMESH_WITH_TUTORIALS)
    message(STATUS "Configuring RINGMesh with tutorials")
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin/tutorials)
    message(STATUS "RINGMesh tutorials directory ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}")
    add_subdirectory(doc/tutorials)
    copy_for_windows()
endif(RINGMESH_WITH_TUTORIALS)

if(RINGMESH_WITH_TESTS)
    message(STATUS "Configuring RINGMesh with tests")
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin/tests)
    message(STATUS "RINGMesh tests directory ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}")
    add_subdirectory(tests)
    copy_for_windows()
endif(RINGMESH_WITH_TESTS)


# Documentation
if(BUILD_DOCUMENTATION)
    message(STATUS "Configuring RINGMesh with doxygen")
    add_subdirectory(doc)
endif()


if (COVERALLS)
set(ringmesh_files ${ringmesh_include} ${ringmesh_src})
    coveralls_setup(
    "${ringmesh_files}" 
        ${COVERALLS_UPLOAD}                 
        "${PROJECT_SOURCE_DIR}/cmake")
endif()

#------------------------------------------------------------------------------------------------
# Configure CPack

set(CPACK_PACKAGE_NAME ${CMAKE_PROJECT_NAME})
set(CPACK_PACKAGE_VERSION_MAJOR ${RINGMesh_VERSION_MAJOR})
set(CPACK_PACKAGE_VERSION_MINOR ${RINGMesh_VERSION_MINOR})
set(CPACK_PACKAGE_VERSION_PATCH ${RINGMesh_VERSION_PATCH})
set(CPACK_PACKAGE_VENDOR "RING-TEAM (www.ring-team.org)")
set(CPACK_SOURCE_GENERATOR "ZIP")

set(CPACK_SOURCE_IGNORE_FILES "/build/;/.hg/;/_CPack_Packages/")

# This must always be last!
include(CPack)
