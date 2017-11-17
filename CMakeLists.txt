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
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
# THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
# PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL ASGA BE LIABLE FOR ANY DIRECT,
# INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
# (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
# ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF
# THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
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

cmake_minimum_required(VERSION 3.1)
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

set(RINGMESH_EXTRA_ARGS
    -DPROJECT_BINARY_DIR:STRING=${project_binary_dir_config}
    -DMG_TETRA:STRING=${MG_TETRA}
    -DRINGMESH_WITH_TETGEN:BOOL=${RINGMESH_WITH_TETGEN}
    -DRINGMESH_WITH_GRAPHICS:BOOL=${RINGMESH_WITH_GRAPHICS}
    -DRINGMESH_TEST_GRAPHICS:BOOL=${RINGMESH_TEST_GRAPHICS}
    -DBUILD_RINGMESH_VIEW:BOOL=${BUILD_RINGMESH_VIEW}
    -DRINGMESH_WITH_UTILITIES:BOOL=${RINGMESH_WITH_UTILITIES}
    -DRINGMESH_WITH_TESTS:BOOL=${RINGMESH_WITH_TESTS}
    -DRINGMESH_WITH_TUTORIALS:BOOL=${RINGMESH_WITH_TUTORIALS}
    -DBUILD_DOCUMENTATION:BOOL=${BUILD_DOCUMENTATION}
    -DBUILD_GEOGRAM_WITHOUT_EXE:BOOL=${BUILD_GEOGRAM_WITHOUT_EXE})


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
# Generate configuration directories for single-configuration generators (Make)
# and run cmake configuration command in each one of them
if(CMAKE_GENERATOR STREQUAL "Unix Makefiles")
    # If the CMAKE_BUILD_TYPE is no defined
    # i.e. at first run of that file
    if(NOT CMAKE_BUILD_TYPE)
      # For all configuration types
      foreach(config IN LISTS CMAKE_CONFIGURATION_TYPES)
          project(MULTI_CONFIG NONE)
          # The binary directory for this configuration
          set(project_binary_dir_config  ${PROJECT_BINARY_DIR}/${config})

          # Create the directory
          file(MAKE_DIRECTORY ${project_binary_dir_config})

          # Launch cmake for this configuration
          # by specifying the CMAKE_BUILD_TYPE and the PROJECT _BINARY_DIR,
          # and pass on all the options
          execute_process(
             COMMAND ${CMAKE_COMMAND} ${CMAKE_SOURCE_DIR}
                -G ${CMAKE_GENERATOR}
                ${RINGMESH_EXTRA_ARGS} 
                -DCMAKE_BUILD_TYPE=${config}
                -DGLOBAL_BINARY_DIR=${PROJECT_BINARY_DIR}
                WORKING_DIRECTORY ${project_binary_dir_config})
      endforeach()

      # Get out when all configurations have been run
      # We do not want to execute the rest of the file in that case
      return()
    endif()
endif(CMAKE_GENERATOR STREQUAL "Unix Makefiles")

if(USE_SUPERBUILD)
    # execute the superbuild
    project(SUPERBUILD)
    if(NOT GLOBAL_BINARY_DIR)
        set(GLOBAL_BINARY_DIR ${PROJECT_BINARY_DIR})
    endif()
    include("${CMAKE_SOURCE_DIR}/cmake/SuperBuild.cmake")
    return()
else()
    # Define the project RINGMesh
    project(RINGMesh)
endif()

#------------------------------------------------------------------------------------------------
# Turn on the ability to create folders to organize projects and files
# It creates "CMakePredefinedTargets" folder by default and adds CMake
# defined projects like INSTALL.vcproj and ZERO_CHECK.vcproj
set_property(GLOBAL PROPERTY USE_FOLDERS ON)


# Define version number
# It is then exported to the configuration file
set(RINGMesh_VERSION_MAJOR 5)
set(RINGMesh_VERSION_MINOR 0)
set(RINGMesh_VERSION_PATCH 0)
set(RINGMesh_VERSION ${RINGMesh_VERSION_MAJOR}.${RINGMesh_VERSION_MINOR}.${RINGMesh_VERSION_PATCH})

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
    if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
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
    set(CMAKE_CXX_STANDARD 11)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread")

    if(CMAKE_BUILD_TYPE STREQUAL "Coverage")
       include(${PROJECT_SOURCE_DIR}/cmake/Coverage.cmake)
    endif()
else(UNIX)
    if(NOT "${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
        message(FATAL_ERROR "RINGMesh on Windows requires Microsoft Visual C++")
    endif()
    # 1800 corresponds to VS 12 2013.
    # See https://en.wikipedia.org/wiki/Microsoft_Visual_C%2B%2B for VS version numbers.
    if(MSVC_VERSION LESS 1800)
        message(FATAL_ERROR "RINGMesh requires Visual Studio greater or equal to VS 12 2013 WIN64")
    endif(MSVC_VERSION LESS 1800)
    
    add_definitions(-DGEO_DYNAMIC_LIBS)
    set(VORPALINE_BUILD_DYNAMIC ON)
endif(UNIX)

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

  set(EXTRA_LIBS ${EXTRA_LIBS} meshgems mg-tetra meshgems_stubs)
endif()

#------------------------------------------------------------------------------------------------
# Collect the library files
add_folder(ringmesh_src ringmesh_include basic)
add_folder(ringmesh_src ringmesh_include geogram_extension)
add_folder(ringmesh_src ringmesh_include geomodel/builder)
add_folder(ringmesh_src ringmesh_include geomodel/core)
add_folder(ringmesh_src ringmesh_include geomodel/tools)
add_folder(ringmesh_src ringmesh_include io)
add_folder(ringmesh_src ringmesh_include mesh)
add_folder(ringmesh_src ringmesh_include tetrahedralize)
add_folder(ringmesh_src ringmesh_include visualize)

# Add the RINGMesh target as a shared library
add_library(RINGMesh SHARED
    ${ringmesh_include}
    ${ringmesh_src})

#------------------------------------------------------------------------------------------------
# Build configuration
# Organize outputs in a Linux way
set_target_properties(RINGMesh PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
set_target_properties(RINGMesh PROPERTIES ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)

# Add include directories of RINGMesh
target_include_directories(RINGMesh PUBLIC ${PROJECT_SOURCE_DIR}/include)

message(STATUS "ZLIB_ROOT = ${ZLIB_ROOT}")
find_package(ZLIB REQUIRED)
list(APPEND CMAKE_MODULE_PATH ${GEOGRAM_INSTALL_PREFIX}/lib/cmake/modules)
message ("GEOGRAM_INSTALL_PREFIX " ${GEOGRAM_INSTALL_PREFIX} )
find_package(Geogram REQUIRED)
find_package(tinyxml2 REQUIRED PATHS ${TINYXML2_INSTALL_PREFIX})
include(${MINIZIP_INSTALL_PREFIX}/cmake/minizip-exports.cmake)
target_link_libraries(RINGMesh PUBLIC Geogram::geogram PRIVATE ZLIB::ZLIB tinyxml2 MINIZIP::minizip)

get_target_property(TEST1 RINGMesh LINK_LIBRARIES)
get_target_property(TEST2 RINGMesh IMPORTED_LINK_DEPENDENT_LIBRARIES)
message(STATUS "TEST1 = ${TEST1}")
message(STATUS "TEST2 = ${TEST2}")

if(RINGMESH_WITH_GRAPHICS)
    target_link_libraries(RINGMesh PUBLIC Geogram::geogram_gfx)
endif(RINGMESH_WITH_GRAPHICS)

# Libraries to which RINGMesh should link
target_link_libraries(RINGMesh LINK_PRIVATE ${EXTRA_LIBS})

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
target_include_directories(RINGMesh PUBLIC ${PROJECT_BINARY_DIR})

#------------------------------------------------------------------------------------------------
# Optional modules configuration

# ringmesh_files definition must be before the binaries (utilities, ringmesh-viewm, etc.).
# Thus this variable can be updated by the called cmake files.
# This line is for clang utilities.
set(ringmesh_files ${ringmesh_include} ${ringmesh_src})

# Collect the IO files
file(GLOB directories RELATIVE ${PROJECT_SOURCE_DIR}/src/ringmesh/io ${PROJECT_SOURCE_DIR}/src/ringmesh/io/*)
foreach(directory ${directories})
    if(IS_DIRECTORY ${PROJECT_SOURCE_DIR}/src/ringmesh/io/${directory})
        source_file_directory(ringmesh_files io/${directory})
    endif()
endforeach()

set(binary_source_dir ${PROJECT_SOURCE_DIR}/src/bin)
if(BUILD_RINGMESH_VIEW)
    message(STATUS "Configure ringmesh-view")
    add_ringmesh_binary(${binary_source_dir}/ringmesh-view.cpp)
    copy_for_windows(${PROJECT_BINARY_DIR}/bin)
endif()

if(RINGMESH_WITH_UTILITIES)
    message(STATUS "Configuring RINGMesh with utilities")
    # Get the paths of the utility files
    file(GLOB utility_sources "${binary_source_dir}/utilities/*.cpp")
    foreach(utility_src ${utility_sources})
        add_ringmesh_utility(${utility_src})
    endforeach()
    copy_for_windows(${PROJECT_BINARY_DIR}/bin/utilities)
endif()

if(RINGMESH_WITH_TUTORIALS)
    message(STATUS "Configuring RINGMesh with tutorials")
    add_subdirectory(doc/tutorials)
    copy_for_windows(${PROJECT_BINARY_DIR}/bin/tutorials)
endif()

if(RINGMESH_WITH_TESTS)
    # Enable testing with CTest
    enable_testing()
    message(STATUS "Configuring RINGMesh with tests")
    add_subdirectory(tests)
    copy_for_windows(${PROJECT_BINARY_DIR}/bin/tests)
endif()


# Documentation
if(BUILD_DOCUMENTATION)
    message(STATUS "Configuring RINGMesh with doxygen")
    add_subdirectory(doc)
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
