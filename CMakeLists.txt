# Copyright (c) 2012-2016, Association Scientifique pour la Geologie et ses Applications (ASGA)
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
# Author Jeanne Pellerin
cmake_minimum_required(VERSION 2.8)

#------------------------------------------------------------------------------------------------
# Options of the RINGMesh project 
# They can be set using CMake GUI, or in command line: cmake -D<var>:<type>=<value> MyProjectFolder

# Configuration modes - Debug and Release by default 
set(CMAKE_CONFIGURATION_TYPES 
      Debug Release CACHE CMAKE_CONFIGURATION_TYPES "Configuration types")

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
option(RINGMESH_WITH_UTILITIES "Compile utility executables" OFF)
option(RINGMESH_WITH_TESTS "Compile test projects" OFF)

# Optional custom steps 
option(BUILD_DOCUMENTATION "Create and install the HTML documentation (requires Doxygen)")

# Optional custom configuration of data directories
if(RINGMESH_WITH_TESTS)
    set(TEST_DIRECTORY ${CMAKE_SOURCE_DIR}/tests)
    set(TEST_SOURCE_DIRECTORY ${TEST_DIRECTORY}/src)
    set(TEST_DATA_INPUT_DIRECTORY "${TEST_DIRECTORY}/data" CACHE PATH "Input data location")
set(TEST_DATA_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/tests/output")
endif()

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
set (RINGMesh_VERSION_MAJOR 2)
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
                -DRINGMESH_WITH_UTILITIES:BOOL=${RINGMESH_WITH_UTILITIES}
                -DRINGMESH_WITH_TESTS:BOOL=${RINGMESH_WITH_TESTS}
			    -DBUILD_DOCUMENTATION:BOOL=${BUILD_DOCUMENTATION}
				-DTEST_DATA_INPUT_DIRECTORY:STRING=${TEST_DATA_INPUT_DIRECTORY}
			    -DTEST_DATA_OUTPUT_DIRECTORY:STRING=${TEST_DATA_OUTPUT_DIRECTORY}                
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

#------------------------------------------------------------------------------------------------
# Platform dependent settings
if(WIN32)
    add_compile_options(-DRINGMESH_EXPORTS) 
	
	# Select the CMAKE_CONFIGURATION_TYPES from Visual Studio
    set(COMPILATION_OPTION --config ${CMAKE_CFG_INTDIR})
elseif(UNIX)
    set(CMAKE_POSITION_INDEPENDENT_CODE TRUE)
    add_compile_options(-Wall -Wextra -Wno-long-long -Wconversion)
    
    # Determine gcc version and activate additional warnings available in latest versions
    execute_process(COMMAND ${CMAKE_C_COMPILER} -dumpversion OUTPUT_VARIABLE GCC_VERSION)

    if (GCC_VERSION VERSION_GREATER 4.3 OR GCC_VERSION VERSION_EQUAL 4.3)
        message(STATUS "GCC version >= 4.3, activating sign conversion warnings")
        add_compile_options(-Wsign-conversion)
    endif()

    if (GCC_VERSION VERSION_GREATER 4.6 OR GCC_VERSION VERSION_EQUAL 4.6)
        message(STATUS "GCC version >= 4.6, activating double promotion warnings")
        add_compile_options(-Wdouble-promotion)
    endif()
    
    # Enable parallel computation
    set(COMPILATION_OPTION -- -j4)    
else(WIN32)
    # Nothing else really supported
    message(FATAL_ERROR "Your platform is not supported. Please modify the cmakelist file.")    
endif(WIN32)

#------------------------------------------------------------------------------------------------
# Configure a header file to pass on some of the CMake settings
# to the source code
set(generated_configure_file ${PROJECT_BINARY_DIR}/ringmesh/ringmesh_config.h)
set(input_configure_file ${PROJECT_SOURCE_DIR}/include/ringmesh/ringmesh_config.h.in)
configure_file(${input_configure_file} ${generated_configure_file})

# We want to be able to include this file
include_directories(${PROJECT_BINARY_DIR})
 
#------------------------------------------------------------------------------------------------
# GEOGRAM 
# Set the path to Geogram code
set(GEOGRAM_PATH ${PROJECT_SOURCE_DIR}/src/third_party/geogram)

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
        -DGEOGRAM_WITH_MEDIT:BOOL=${RINGMESH_WITH_GRAPHICS}
        -DGEOGRAM_WITH_GRAPHICS:BOOL=${RINGMESH_WITH_GRAPHICS}
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
include_directories(${GEOGRAM_PATH}/src/lib)

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
  include_directories(${MG_TETRA}/include)
 
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

#------------------------------------------------------------------------------------------------
# Collect the library files
file(GLOB ringmesh_src "${PROJECT_SOURCE_DIR}/src/ringmesh/*.cpp")
file(GLOB ringmesh_include "${PROJECT_SOURCE_DIR}/include/ringmesh/*.h")
file(GLOB zlib_src "${PROJECT_SOURCE_DIR}/src/third_party/zlib/*.c")
file(GLOB zlib_include "${PROJECT_SOURCE_DIR}/include/third_party/zlib/*.h")

# Add configure file to ringmesh includes
set(ringmesh_include ${ringmesh_include} ${generated_configure_file})

# Group files of a project in directories
# For this to work, all the files should be added to the sources of the project
# and the global property USE_FOLDERS must be ON
source_group("Zlib" FILES ${zlib_src} ${zlib_include})
source_group("Source Files" FILES ${ringmesh_src})
source_group("Header Files" FILES ${ringmesh_include})

#------------------------------------------------------------------------------------------------
# Build configuration
# Organize outputs in a Linux way
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)

# Add include directories of RINGMesh
include_directories(${PROJECT_SOURCE_DIR}/include)

# Add the RINGMesh target as a shared library
add_library(RINGMesh SHARED
    ${ringmesh_include} 
    ${ringmesh_src} 
    ${zlib_src}
    ${zlib_include})

# RINGMesh depends on Geogram
add_dependencies(RINGMesh geogram_ext)

# Libraries to which RINGMesh should link
target_link_libraries(RINGMesh ${EXTRA_LIBS})

if(WIN32)
    # On windows, without proper installation steps, we need to 
    # copy of Geogram dll and pdb information to RINGMesh
    # to be able to launch RINGMesh utilities and tests from the debugger
    add_custom_command(TARGET RINGMesh POST_BUILD 
       COMMAND  "${CMAKE_COMMAND}" -E copy_directory 
           "${GEOGRAM_PATH_BIN}/bin/$<CONFIGURATION>"
           "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/$<CONFIGURATION>" 
       COMMENT "Copy geogram binaries")
endif(WIN32)


#------------------------------------------------------------------------------------------------
# Optional modules configuration
if(RINGMESH_WITH_GRAPHICS)
    message(STATUS "Configuring RINGMesh with graphics")
    add_subdirectory(src/ringmeshview)
endif(RINGMESH_WITH_GRAPHICS)

if(RINGMESH_WITH_UTILITIES)
    message(STATUS "Configuring RINGMesh with utilities")
    add_subdirectory(src/ringmeshconvert)
    add_subdirectory(src/ringmesh_surface_convert)
    add_subdirectory(src/ringmeshstats)
endif(RINGMESH_WITH_UTILITIES)

if(RINGMESH_WITH_TESTS)
    message(STATUS "Configuring RINGMesh with tests")
    add_subdirectory(tests)
endif(RINGMESH_WITH_TESTS)


# Documentation
if(BUILD_DOCUMENTATION)
    message(STATUS "Configuring RINGMesh with doxygen")
    add_subdirectory(doc)
endif()

#------------------------------------------------------------------------------------------------
# Configure CPack

set(CPACK_PACKAGE_NAME ${CMAKE_PROJECT_NAME})
set(CPACK_PACKAGE_VERSION_MAJOR ${RINGMesh_VERSION_MAJOR})
set(CPACK_PACKAGE_VERSION_MINOR ${RINGMesh_VERSION_MINOR})
set(CPACK_PACKAGE_VERSION_PATCH ${RINGMesh_VERSION_PATCH})
set(CPACK_PACKAGE_VENDOR "RING-TEAM (www.ring-team.org)")

if(WIN32)
    set(CPACK_SOURCE_GENERATOR "ZIP")
else()
    set(CPACK_SOURCE_GENERATOR "TGZ")
endif()

set(CPACK_SOURCE_IGNORE_FILES "/build/;/.hg/;/_CPack_Packages/")

# This must always be last!
include(CPack)
