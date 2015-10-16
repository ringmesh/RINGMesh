
include(${GEOGRAM_SOURCE_DIR}/cmake/utilities.cmake)
include(${GEOGRAM_SOURCE_DIR}/cmake/platforms/${VORPALINE_PLATFORM}/config.cmake)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)

# Static versus dynamic builds
if(VORPALINE_BUILD_DYNAMIC)
    set(BUILD_SHARED_LIBS TRUE)
    # Object files in OBJECT libraries are compiled in static mode, even if
    # BUILD_SHARED_LIBS is true! We must set CMAKE_POSITION_INDEPENDENT_CODE
    # to force compilation in dynamic mode.
    set(CMAKE_POSITION_INDEPENDENT_CODE TRUE)
    add_definitions(-DGEO_DYNAMIC_LIBS)
else()
    set(BUILD_SHARED_LIBS FALSE)
endif()

##############################################################################

# RELATIVE_OUTPUT_DIR: Binary build directory relative to source directory
# Since the build tree is a subdirectory of the source tree, it is
#  found by replacing the source dir with an empty string in the bin dir.

string(REPLACE ${CMAKE_SOURCE_DIR} "" RELATIVE_OUTPUT_DIR ${PROJECT_BINARY_DIR})

# RELATIVE_BIN_DIR
# RELATIVE_LIB_DIR
#  Directories where libraries and executables are generated, relative to
# source directory.

if(WIN32)
    set(MSVC_CONFIG \$\(Configuration\))
 	set(RELATIVE_BIN_DIR ${RELATIVE_OUTPUT_DIR}/bin/${MSVC_CONFIG}/)
	set(RELATIVE_LIB_DIR ${RELATIVE_OUTPUT_DIR}/lib/${MSVC_CONFIG}/)
else()
	set(RELATIVE_BIN_DIR ${RELATIVE_OUTPUT_DIR}/bin/)
	set(RELATIVE_LIB_DIR ${RELATIVE_OUTPUT_DIR}/lib/)
endif()

##############################################################################

include_directories(${GEOGRAM_SOURCE_DIR}/src/lib)
link_directories(${GEOGRAM_SOURCE_DIR}/${RELATIVE_LIB_DIR})

##############################################################################
