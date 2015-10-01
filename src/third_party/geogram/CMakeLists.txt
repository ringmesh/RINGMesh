##############################################################################
# Geogram/Vorpaline root CMakeList
##############################################################################

# Make sure that VORPALINE_PLATFORM is defined

if(NOT DEFINED VORPALINE_PLATFORM)
     if(WIN32) 
        message( 
           STATUS
           " Using Win-vs-generic (default),\n"
           " (if need be, use CMake variable VORPALINE_PLATFORM to override)."
        )
        set(VORPALINE_PLATFORM Win-vs-generic)
     else()
        message(FATAL_ERROR
           " CMake variable VORPALINE_PLATFORM is not defined.\n"
           " Please run configure.{sh,bat} to setup the build tree."
        )
     endif()
endif()

# CMake 2.8.11 is required for 2 reasons:
# - it is the first version that fully supports the specification of Visual
# Studio toolsets (v110_xp).
# - it is the version that supports the command string(TIMESTAMP ...)

cmake_minimum_required(VERSION 2.8.11)

##############################################################################

# Determine whether Geogram is built with Vorpaline
if(IS_DIRECTORY ${CMAKE_SOURCE_DIR}/src/lib/vorpalib)
   message(STATUS "Configuring build for Geogram + Vorpaline")
   set(GEOGRAM_WITH_VORPALINE TRUE)
else()
   message(STATUS "Configuring build for standalone Geogram (without Vorpaline)")
   set(GEOGRAM_WITH_VORPALINE FALSE)   
endif()

if(GEOGRAM_WITH_VORPALINE)
   project(Vorpaline)
else()
   project(Geogram)
endif()

set(VORPALINE_VERSION_MAJOR 1)
set(VORPALINE_VERSION_MINOR 0)
set(VORPALINE_VERSION_PATCH 0)
set(VORPALINE_VERSION ${VORPALINE_VERSION_MAJOR}.${VORPALINE_VERSION_MINOR}.${VORPALINE_VERSION_PATCH})

# Optional modules
# So why are they set to TRUE ??? [JP]
# set(GEOGRAM_WITH_MEDIT TRUE)
# set(GEOGRAM_WITH_GRAPHICS TRUE)

# Determine the current Build-OS (Build-platform without the compiler info)
string(REGEX REPLACE "-[^-]+$" "" VORPALINE_OS ${VORPALINE_PLATFORM})

# Determine the current build date
string(TIMESTAMP VORPALINE_BUILD_DATE "%Y-%m-%d %H:%M:%S")
string(TIMESTAMP YEAR "%Y")

# Determine the current build number
# This is set by Jenkins in environment variable BUILD_NUMBER
set(VORPALINE_BUILD_NUMBER $ENV{BUILD_NUMBER})


##############################################################################
# Get SVN revision info

if(GEOGRAM_WITH_VORPALINE)
   find_package(Subversion QUIET)
   if(NOT SUBVERSION_FOUND)
       message(WARNING "Subversion executable not found - cannot determine current revision")
   else()
       Subversion_WC_INFO(${PROJECT_SOURCE_DIR} Vorpaline)
       message(STATUS "Vorpaline revision is ${Vorpaline_WC_REVISION}")
       set(VORPALINE_SVN_REVISION ${Vorpaline_WC_REVISION})
   endif()
endif()   

##############################################################################
# Build configuration

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)

if(EXISTS ${CMAKE_SOURCE_DIR}/CMakeOptions.txt)
   message(STATUS "Using local options file: ${CMAKE_SOURCE_DIR}/CMakeOptions.txt")
   include(${CMAKE_SOURCE_DIR}/CMakeOptions.txt)
endif()


include(${CMAKE_SOURCE_DIR}/cmake/utilities.cmake)
include(${CMAKE_SOURCE_DIR}/cmake/platforms/${VORPALINE_PLATFORM}/config.cmake)

# This test is there to keep CMake happy about unused variable CMAKE_BUILD_TYPE
if(CMAKE_BUILD_TYPE STREQUAL "")
endif()

# Static versus dynamic builds
if(VORPALINE_BUILD_DYNAMIC)
    set(BUILD_SHARED_LIBS TRUE)
    # Object files in OBJECT libraries are compiled in static mode, event if
    # BUILD_SHARED_LIBS is true! We must set CMAKE_POSITION_INDEPENDENT_CODE
    # to force compilation in dynamic mode.
    set(CMAKE_POSITION_INDEPENDENT_CODE TRUE)
    add_definitions(-DGEO_DYNAMIC_LIBS)
else()
    set(BUILD_SHARED_LIBS FALSE)
endif()

##############################################################################
# Vorpaline sources

include_directories(${CMAKE_SOURCE_DIR}/src/lib)

add_subdirectory(src/lib/geogram)
add_subdirectory(src/bin/vorpastat)
add_subdirectory(src/bin/vorpacomp)

if(GEOGRAM_WITH_FPG)
   add_subdirectory(src/bin/fpg)
endif()

if(GEOGRAM_WITH_VORPALINE)
   add_subdirectory(src/lib/vorpalib)
   add_subdirectory(src/bin/vorpaline)
else()
   add_subdirectory(src/bin/vorpalite)
endif()

if(GEOGRAM_WITH_GRAPHICS)
   add_subdirectory(src/lib/geogram_gfx)
   add_subdirectory(src/bin/vorpaview)

   if(GEOGRAM_WITH_MEDIT)
       add_subdirectory(src/bin/medit)
   endif()

endif()

add_subdirectory(src/test)
add_subdirectory(doc)
add_subdirectory(tests)

##############################################################################
# Cleanup from previous builds

file(REMOVE ${CMAKE_BINARY_DIR}/doc/LICENSE.txt)


##############################################################################
# Vorpaline installation

# Configure CPack

SET(CPACK_PACKAGE_NAME ${CMAKE_PROJECT_NAME})
set(CPACK_SYSTEM_NAME ${VORPALINE_OS})
SET(CPACK_PACKAGE_VENDOR "INRIA - ALICE")
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "A flexible mesh generator")

set(CPACK_PACKAGE_VERSION_MAJOR ${VORPALINE_VERSION_MAJOR})
set(CPACK_PACKAGE_VERSION_MINOR ${VORPALINE_VERSION_MINOR})
set(CPACK_PACKAGE_VERSION_PATCH ${VORPALINE_VERSION_PATCH})
set(CPACK_PACKAGE_VERSION ${VORPALINE_VERSION})
set(CPACK_PACKAGE_FILE_NAME ${CPACK_PACKAGE_NAME}-${CPACK_PACKAGE_VERSION})
set(CPACK_COMPONENT_INCLUDE_TOPLEVEL_DIRECTORY true)

if(WIN32)
    set(CPACK_GENERATOR ZIP)
else()
    set(CPACK_GENERATOR TGZ)
endif()

# Enable component-based packaging for archive generators (TGZ, ZIP)
set(CPACK_ARCHIVE_COMPONENT_INSTALL ON)

set(CPACK_COMPONENTS_ALL runtime devkit devkit-full doc-devkit doc-devkit-full)
set(CPACK_COMPONENTS_GROUPING "IGNORE")

set(CPACK_COMPONENT_RUNTIME_DISPLAY_NAME "Vorpaline Application")
set(CPACK_COMPONENT_RUNTIME_GROUP "Runtime")

set(CPACK_COMPONENT_DEVKIT_DISPLAY_NAME "Vorpaline Developer Kit")
set(CPACK_COMPONENT_DEVKIT_GROUP "Development")

set(CPACK_COMPONENT_DEVKIT-FULL_DISPLAY_NAME "Vorpaline Full Developer Kit")
set(CPACK_COMPONENT_DEVKIT-FULL_GROUP "Development")

set(CPACK_COMPONENT_DOC-DEVKIT_DISPLAY_NAME "Vorpaline API Developer Kit Documentation")
set(CPACK_COMPONENT_DOC-DEVKIT_GROUP "Documentation")

set(CPACK_COMPONENT_DOC-DEVKIT-FULL_DISPLAY_NAME "Vorpaline Full Developer Kit Documentation")
set(CPACK_COMPONENT_DOC-DEVKIT-FULL_GROUP "Documentation")

set(CPACK_COMPONENT_DOC-DEVKIT-INTERNAL_DISPLAY_NAME "Vorpaline Internal Developer Kit Documentation")
set(CPACK_COMPONENT_DOC-DEVKIT-INTERNAL_GROUP "Documentation")

# Copy the helper script to build individual packages to the binary directory
configure_file(
    tools/make_package.pl.in
    make_package.pl
    @ONLY
)

# This must always be last!
include(CPack)

