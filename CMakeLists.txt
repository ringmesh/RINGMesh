cmake_minimum_required(VERSION 2.8)
project(GRGMesh)

IF (EXISTS CMakeOptions.txt)
   MESSAGE(INFO: "Using local options file: " CMakeOptions.txt)
   INCLUDE(CMakeOptions.txt)
ENDIF(EXISTS CMakeOptions.txt)

IF(WIN32)
  ADD_DEFINITIONS(-DGRGMESH_EXPORTS)
ENDIF(WIN32)

AUX_SOURCE_DIRECTORY(src SOURCES)
INCLUDE_DIRECTORIES(include)

SET(EXECUTABLE_OUTPUT_PATH lib)
add_library(GRGMesh ${SOURCES})

