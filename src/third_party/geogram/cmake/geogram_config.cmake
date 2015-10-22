
if(NOT DEFINED GEOGRAM_SOURCE_DIR)
   set(GEOGRAM_SOURCE_DIR ${CMAKE_SOURCE_DIR})
endif()

if(EXISTS ${GEOGRAM_SOURCE_DIR}/CMakeOptions.txt)
   message(STATUS "Using local options file: ${GEOGRAM_SOURCE_DIR}/CMakeOptions.txt")
   include(${GEOGRAM_SOURCE_DIR}/CMakeOptions.txt)
endif()

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

# Determine whether Geogram is built with Vorpaline
if(IS_DIRECTORY ${GEOGRAM_SOURCE_DIR}/src/lib/vorpalib)
   message(STATUS "Configuring build for Geogram + Vorpaline")
   set(GEOGRAM_WITH_VORPALINE TRUE)
else()
   message(STATUS "Configuring build for standalone Geogram (without Vorpaline)")
   set(GEOGRAM_WITH_VORPALINE FALSE)   
endif()

# This test is there to keep CMake happy about unused variable CMAKE_BUILD_TYPE
if(CMAKE_BUILD_TYPE STREQUAL "")
endif()
