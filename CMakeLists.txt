cmake_minimum_required(VERSION 2.8)


message("CMAKE_BUILD_TYPE is equal to " ${CMAKE_BUILD_TYPE})
# debug suffix for libraries

if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    set(SUFF "_d")
    MESSAGE("DEBUG! SUFFIX = ${SUFF}")
endif(CMAKE_BUILD_TYPE STREQUAL "Debug")

set(GRGMESH_PROJECT GRGMesh${SUFF})
project(${GRGMESH_PROJECT})


function(aux_source_directories var folder)
    set(sources)
    foreach(dir ${ARGN})
        file(GLOB _sources RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} "${dir}/*.[ch]" "${dir}/*.[ch]pp")
        list(APPEND sources ${_sources})
    endforeach()
    source_group(${folder} FILES ${sources})
    set(${var} ${${var}} ${sources} PARENT_SCOPE)
endfunction()

# Options
IF(EXISTS ${CMAKE_SOURCE_DIR}/CMakeOptions.txt)
   MESSAGE(INFO: "Using local options file: " .${CMAKE_SOURCE_DIR}/CMakeOptions.txt)
   INCLUDE(${CMAKE_SOURCE_DIR}/CMakeOptions.txt)
ENDIF(EXISTS ${CMAKE_SOURCE_DIR}/CMakeOptions.txt)



# Files
AUX_SOURCE_DIRECTORY(src/grgmesh SOURCES)
AUX_SOURCE_DIRECTORY(src/third_party/zip SOURCES)
INCLUDE_DIRECTORIES(include)

# Geogram
INCLUDE_DIRECTORIES(src/third_party/geogram/src/lib)
IF(WIN32)
    ADD_DEFINITIONS(-DGRGMESH_EXPORTS)
	ADD_DEFINITIONS(-DGEO_DYNAMIC_LIBS)
    aux_source_directories(SOURCES "Header Files" include/grgmesh)
    link_directories(${CMAKE_SOURCE_DIR}/build/geogram/${GEOGRAM_PLATFORM}/lib/${CMAKE_BUILD_TYPE})
ELSE(WIN32)
    link_directories(${CMAKE_SOURCE_DIR}/build/geogram/${GEOGRAM_PLATFORM}-${CMAKE_BUILD_TYPE}/lib)
ENDIF(WIN32)

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


#To use OpenMP 
IF(WIN32) 
       SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /wd4275 /wd4996 /MP")
       SET(CMAKE_C_FLAGS "${CMAKE_CXX_FLAGS} /wd4275 /wd4996 /MP")
ELSE(WIN32)
       SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
       SET(CMAKE_C_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
ENDIF(WIN32)

IF(UNIX)
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")
    SET(CMAKE_C_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")
ENDIF(UNIX)

# Output library
SET(LIBRARY_OUTPUT_PATH lib)
add_library(${GRGMESH_PROJECT} SHARED ${SOURCES})
target_link_libraries(${GRGMESH_PROJECT} geogram)

# Additional link libraries
IF(${USE_MG_TETRA})
  TARGET_LINK_LIBRARIES(${GRGMESH_PROJECT} meshgems mg-tetra)
ENDIF(${USE_MG_TETRA})

IF(WIN32)
    target_link_libraries(${GRGMESH_PROJECT} psapi)
ENDIF(WIN32)


if(UNIX) 
	install(DIRECTORY include/grgmesh DESTINATION include)
	install(DIRECTORY include/geogram DESTINATION include)
		
	install(TARGETS ${GRGMESH_PROJECT}
		LIBRARY DESTINATION lib
		RUNTIME DESTINATION lib 
		CONFIGURATIONS All)	
endif(UNIX)
