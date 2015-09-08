cmake_minimum_required(VERSION 2.8)

project(RINGMesh)


##############################################################################
# Options
IF(EXISTS ${CMAKE_SOURCE_DIR}/CMakeOptions.txt)
   MESSAGE(INFO: "Using local options file: " .${CMAKE_SOURCE_DIR}/CMakeOptions.txt)
   INCLUDE(${CMAKE_SOURCE_DIR}/CMakeOptions.txt)
ENDIF(EXISTS ${CMAKE_SOURCE_DIR}/CMakeOptions.txt)


##############################################################################
# Files
AUX_SOURCE_DIRECTORY(src/ringmesh SOURCES)
AUX_SOURCE_DIRECTORY(src/third_party/zlib SOURCES)
INCLUDE_DIRECTORIES(include)


function(aux_source_directories var folder)
    set(sources)
    foreach(dir ${ARGN})
        file(GLOB _sources RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} "${dir}/*.[ch]" "${dir}/*.[ch]pp")
        list(APPEND sources ${_sources})
    endforeach()
    source_group(${folder} FILES ${sources})
    set(${var} ${${var}} ${sources} PARENT_SCOPE)
endfunction()

IF(WIN32)
    aux_source_directories(SOURCES "Header Files" include/ringmesh)
ENDIF(WIN32)
##############################################################################
# Geogram
INCLUDE_DIRECTORIES(src/third_party/geogram/src/lib)
IF(WIN32)
    ADD_DEFINITIONS(-DRINGMESH_EXPORTS)
	ADD_DEFINITIONS(-DGEO_DYNAMIC_LIBS)

    link_directories(${CMAKE_SOURCE_DIR}/build/geogram/${PLATFORM}/lib)
ELSE(WIN32)
    link_directories(${CMAKE_SOURCE_DIR}/build/geogram/${PLATFORM}-${CMAKE_BUILD_TYPE}/lib)
ENDIF(WIN32)

##############################################################################
# Geogram_gfx
IF(NOT DEFINED RINGMESH_WITH_GRAPHICS)
  SET(RINGMESH_WITH_GRAPHICS FALSE)
ENDIF(NOT DEFINED RINGMESH_WITH_GRAPHICS)
if(${RINGMESH_WITH_GRAPHICS})
    add_definitions(-DRINGMESH_WITH_GRAPHICS)
endif()


##############################################################################
# MG-Tetra
IF(NOT DEFINED USE_MG_TETRA)
  SET(USE_MG_TETRA FALSE)
ENDIF(NOT DEFINED USE_MG_TETRA)
IF(${USE_MG_TETRA})	
  add_definitions(-DUSE_MG_TETRA)
  IF(${MG_TETRA} MATCHES .+)
    MESSAGE(STATUS "Using MG-Tetra directory =  ${MG_TETRA}")
    SET(MG_TETRA_DIR  ${MG_TETRA})
  ELSEIF( $ENV{MG_TETRA} MATCHES .+)
    MESSAGE(STATUS "Using MG-Tetra directory =  $ENV{MG_TETRA}")
    SET(MG_TETRA_DIR  $ENV{MG_TETRA})
  ELSE(${MG_TETRA} MATCHES .+)
    MESSAGE(STATUS "No MG-Tetra directory found")
  ENDIF(${MG_TETRA} MATCHES .+)
  INCLUDE_DIRECTORIES(${MG_TETRA_DIR}/include)
 
  if( CMAKE_SIZEOF_VOID_P EQUAL 8 )
    SET(PLATFORM _64)
  ELSE(CMAKE_SIZEOF_VOID_P EQUAL 8 )
    SET(PLATFORM)
  ENDIF(CMAKE_SIZEOF_VOID_P EQUAL 8 )
  
  IF(WIN32)
    LINK_DIRECTORIES(${MG_TETRA_DIR}/lib/WIN7${PLATFORM}_VC10)
  ELSE(WIN32)
    LINK_DIRECTORIES(${MG_TETRA_DIR}/lib/Linux${PLATFORM})
  ENDIF(WIN32)
ENDIF(${USE_MG_TETRA})

##############################################################################
#To use OpenMP 
if(NOT DEFINED USE_OPENMP)
    SET(USE_OPENMP TRUE)
endif()
IF(${USE_OPENMP})
    IF(WIN32) 
        SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /wd4275 /wd4996 /MP")
        SET(CMAKE_C_FLAGS "${CMAKE_CXX_FLAGS} /wd4275 /wd4996 /MP")
    ELSE(WIN32)
        SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
        SET(CMAKE_C_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
    ENDIF(WIN32)
    add_definitions(-DRINGMESH_USE_OPENMP)
ENDIF()

IF(UNIX)
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")
    SET(CMAKE_C_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")
ENDIF(UNIX)


##############################################################################
# Build configuration

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)

add_library(RINGMesh SHARED ${SOURCES})
target_link_libraries(RINGMesh geogram)

##############################################################################
# Additional optional link libraries
IF(${USE_MG_TETRA})
  target_link_libraries(RINGMesh meshgems mg-tetra)
ENDIF(${USE_MG_TETRA})
if(${RINGMESH_WITH_GRAPHICS})
  TARGET_LINK_LIBRARIES(RINGMesh geogram_gfx)
endif()
IF(WIN32)
    target_link_libraries(RINGMesh psapi)
ENDIF(WIN32)


# Output binary
if(${RINGMESH_WITH_GRAPHICS})
    add_subdirectory(src/bin/ringmeshview)
endif()
add_subdirectory(src/bin/ringmeshconvert)
include(tests/tests.cmake)
