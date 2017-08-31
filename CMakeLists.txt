# Copyright (c) 2012-2017, Association Scientifique pour la Geologie et ses
# Applications (ASGA). All rights reserved.
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
#------------------------------------------------------------------------------------------------
# Root CMakeList of RINGMesh project
#------------------------------------------------------------------------------------------------

cmake_minimum_required(VERSION 2.8)
set(USER_CONFIG "${CMAKE_SOURCE_DIR}/cmake/UserConfig.cmake" CACHE PATH
  "Path to optional user configuration file.")

#------------------------------------------------------------------------------------------------
# Options of the RINGMesh project
# They can be set
# 1 - using CMake GUI,
# 2 - in command line: cmake -D<var>:<type>=<value> MyProjectFolder
# 3 - by creating a UserConfig.cmake file in the cmake directory (see DefaultConfig.cmake)

# Load user settings before the defaults - this way the defaults will not
# overwrite the user set options. If the user has not set all options, we still
# have the defaults.
message(STATUS "(optional) USER_CONFIG = ${USER_CONFIG}")
include("${USER_CONFIG}" OPTIONAL)
include("${CMAKE_SOURCE_DIR}/cmake/DefaultConfig.cmake")

#------------------------------------------------------------------------------------------------
#ensure consistent choice in configuration options
include(CMakeDependentOption)
# Option to build ringmesh-view. Disable when RINGMESH_WITH_GRAPHICS
# is OFF. The fist time RINGMESH_WITH_GRAPHICS is ON, this
# option is ON.
CMAKE_DEPENDENT_OPTION(BUILD_RINGMESH_VIEW "Build ringmesh-view" ON
    "RINGMESH_WITH_GRAPHICS" OFF)
# Option to test ringmesh-view. Disable when BUILD_RINGMESH_VIEW
# is OFF. The fist time BUILD_RINGMESH_VIEW is ON, this
# option is ON.
CMAKE_DEPENDENT_OPTION(RINGMESH_TEST_GRAPHICS "Test ringmesh-view" ON
    "BUILD_RINGMESH_VIEW;RINGMESH_WITH_TESTS" OFF)

#------------------------------------------------------------------------------------------------
# Turn on the ability to create folders to organize projects and files
# It creates "CMakePredefinedTargets" folder by default and adds CMake
# defined projects like INSTALL.vcproj and ZERO_CHECK.vcproj
set_property(GLOBAL PROPERTY USE_FOLDERS ON)

#------------------------------------------------------------------------------------------------
# Additional cmake modules
include(ExternalProject)


# Enable testing with CTest
enable_testing()

#------------------------------------------------------------------------------------------------
# Define the project RINGMesh
project(RINGMesh)

# Define version number
# It is then exported to the configuration file
set (RINGMesh_VERSION_MAJOR 4)
set (RINGMesh_VERSION_MINOR 1)
set (RINGMesh_VERSION_PATCH 0)
set (RINGMesh_VERSION ${RINGMesh_VERSION_MAJOR}.${RINGMesh_VERSION_MINOR}.${RINGMesh_VERSION_PATCH})


#------------------------------------------------------------------------------------------------
# Generate configuration directories for single-configuration generators (Make)
# and run cmake configuration command in each one of them
if(CMAKE_GENERATOR STREQUAL "Unix Makefiles")
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
                -DPROPAGATE_COMPILER_TO_THIRD_PARTIES:BOOL=${PROPAGATE_COMPILER_TO_THIRD_PARTIES}
                WORKING_DIRECTORY ${project_binary_dir_config})
      endforeach()

      # Get out when all configurations have been run
      # We do not want to execute the rest of the file in that case
      return()
    endif()
else(CMAKE_GENERATOR STREQUAL "Unix Makefiles")
    set(GLOBAL_BINARY_DIR  ${PROJECT_BINARY_DIR}/..)
endif(CMAKE_GENERATOR STREQUAL "Unix Makefiles")

message(STATUS "RINGMesh binary directory is: ${PROJECT_BINARY_DIR}")
message(STATUS "RINGMesh source directory is: ${PROJECT_SOURCE_DIR}")

include(cmake/utils.cmake)

#------------------------------------------------------------------------------------------------
# Platform dependent settings
if(UNIX)
    set(CMAKE_POSITION_INDEPENDENT_CODE TRUE)
    add_compile_options(-Wall -Wextra -Wno-long-long -Wconversion 
        -Wsign-conversion -Wdouble-promotion -Wno-attributes)

    # Determine gcc version and activate additional warnings available in latest versions
    execute_process(COMMAND ${CMAKE_CXX_COMPILER} -dumpversion OUTPUT_VARIABLE COMPILER_VERSION)

    # C++11 fully works on Linux from gcc/g++ 4.8 (before not fully operational).
    # On MacOS, it seems that C++11 is working since gcc/g++ 4.2.
    if(APPLE)
        if (NOT "${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
            message(FATAL_ERROR "RINGMesh on Apple requires Clang compiler")
        endif()
        set(CMAKE_MACOSX_RPATH 1)
        set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
    endif(APPLE)
    if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
        set(supported_clang_version 3.8)
        if(COMPILER_VERSION VERSION_LESS ${supported_clang_version})
            message(FATAL_ERROR "RINGMesh require Clang version >= ${supported_clang_version}")
        endif()
    elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
        set(supported_gcc_version 4.8)
        if(COMPILER_VERSION VERSION_LESS ${supported_gcc_version})
            message(FATAL_ERROR "RINGMesh require G++ version >= ${supported_gcc_version}")
        endif()
    endif()
    # In Cmake > 3, there is the option: set (CMAKE_CXX_STANDARD 11). To see in the future.
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -pthread")

    if(CMAKE_BUILD_TYPE STREQUAL "Coverage")
       set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${PROJECT_SOURCE_DIR}/cmake)
       set(COVERALLS ON)
       set(COVERALLS_UPLOAD ON)
       include(Coveralls)
       coveralls_turn_on_coverage()
    endif()
endif(UNIX)

if(CMAKE_GENERATOR STREQUAL "Unix Makefiles")
    # Enable parallel computation
    set(cpuinfo_file "/proc/cpuinfo")
    if(EXISTS "${cpuinfo_file}")
        file(STRINGS "${cpuinfo_file}" procs REGEX "^processor.: [0-9]+$")
        list(LENGTH procs PROCESSOR_COUNT)
    endif()
    set(COMPILATION_OPTION -- -j${PROCESSOR_COUNT})
else(CMAKE_GENERATOR STREQUAL "Unix Makefiles")
    # Select the CMAKE_CONFIGURATION_TYPES from Visual Studio or Xcode.
    # This enables in particular that dependencies of RINGMesh (e.g., geogram, zlib, etc.)
    # are automatically built with the same build type as RINGMesh.
    set(COMPILATION_OPTION --config ${CMAKE_CFG_INTDIR})
endif(CMAKE_GENERATOR STREQUAL "Unix Makefiles")

#------------------------------------------------------------------------------------------------
# Compil Third Parties
include("${CMAKE_SOURCE_DIR}/cmake/CompilThirdParties.cmake")


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
file(GLOB directories RELATIVE ${PROJECT_SOURCE_DIR}/src/ringmesh ${PROJECT_SOURCE_DIR}/src/ringmesh/*)
foreach(directory ${directories})
    if(IS_DIRECTORY ${PROJECT_SOURCE_DIR}/src/ringmesh/${directory})
        source_file_directory(ringmesh_src ${directory})
        include_file_directory(ringmesh_include ${directory})
    endif()
endforeach()

#------------------------------------------------------------------------------------------------
# Build configuration
# Organize outputs in a Linux way
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)

# Add include directories of RINGMesh
include_directories(${PROJECT_SOURCE_DIR}/include)

# Add the RINGMesh target as a shared library
add_library(RINGMesh SHARED
    ${ringmesh_include}
    ${ringmesh_src})

# RINGMesh depends on Geogram Tinyxml2 and zlib
add_dependencies(RINGMesh geogram_ext tinyxml2_ext zlib_ext minizip_ext)

# Libraries to which RINGMesh should link
target_link_libraries(RINGMesh ${EXTRA_LIBS})

#------------------------------------------------------------------------------------------------
# file automatically generated by cmake

# generate configure file to pass on some of the CMake settings
# to the source code
# Add configure file to ringmesh includes

set(generated_configure_file ${PROJECT_BINARY_DIR}/ringmesh/ringmesh_config.h)
set(input_configure_file ${PROJECT_SOURCE_DIR}/include/ringmesh/basic/ringmesh_config.h.in)
configure_file(${input_configure_file} ${generated_configure_file})
set(ringmesh_include ${ringmesh_include} ${generated_configure_file})

#generate macros for API export and add itbto include directories
include (GenerateExportHeader)
set(generated_export_file ${PROJECT_BINARY_DIR}/ringmesh/ringmesh_export.h)
set(ringmesh_include ${ringmesh_include} ${generated_export_file})
generate_export_header(RINGMesh EXPORT_MACRO_NAME RINGMESH_API EXPORT_FILE_NAME ${generated_export_file} )

# We want to be able to include these file
include_directories(${PROJECT_BINARY_DIR})

#------------------------------------------------------------------------------------------------
# Optional modules configuration
set(binary_source_dir ${PROJECT_SOURCE_DIR}/src/bin)
if(BUILD_RINGMESH_VIEW)
    message(STATUS "Configure ringmesh-view")
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)
    message(STATUS "ringmesh-view directory ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}")
    include(${CMAKE_SOURCE_DIR}/third_party/geogram/cmake/opengl.cmake)
    add_ringmesh_binary(${binary_source_dir}/ringmesh-view geogram_gfx ${OPENGL_LIBRARIES})
    copy_for_windows()
endif(BUILD_RINGMESH_VIEW)

if(RINGMESH_WITH_UTILITIES)
    message(STATUS "Configuring RINGMesh with utilities")
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin/utilities)
    message(STATUS "RINGMesh utilities directory ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}")

    # Get the paths of the utility files
    file(GLOB utility_sources "${binary_source_dir}/utilities/*.cpp")
    foreach(utility_src ${utility_sources})
        add_ringmesh_binary(${utility_src})
    endforeach(utility_src ${utility_sources})
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


set(ringmesh_files ${ringmesh_include} ${ringmesh_src})
if(COVERALLS)
    coveralls_setup(
    "${ringmesh_files}"
        ${COVERALLS_UPLOAD}
        "${PROJECT_SOURCE_DIR}/cmake")
endif()

# additional target to perform clang-tidy run, requires clang-tidy
find_program(CLANG_TIDY "clang-tidy")
if(CLANG_TIDY)
    message(STATUS "Configuring RINGMesh with clang-tidy")
    get_directory_property(dirs INCLUDE_DIRECTORIES)
    foreach(dir ${dirs})
        set(include_dirs -I${dir} ${include_dirs})
    endforeach()
    add_custom_target(
        tidy
        COMMAND ${CLANG_TIDY}
        ${ringmesh_files}
        -config=''
        --
        -std=c++11
        ${include_dirs}
    )
endif()

# additional target to perform clang-format run, requires clang-format
find_program(CLANG_FORMAT "clang-format")
if(CLANG_FORMAT)
    message(STATUS "Configuring RINGMesh with clang-format")
    add_custom_target(
        format
        COMMAND ${CLANG_FORMAT}
        -style=file
        -i
        ${ringmesh_files}
    )
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
