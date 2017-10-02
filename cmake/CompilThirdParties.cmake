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
# GEOGRAM
# Set the path to Geogram code
set(GEOGRAM_PATH ${PROJECT_SOURCE_DIR}/third_party/geogram)

# Geogram platform dependent settings
if(WIN32)
    set(GEOGRAM_PATH_BIN ${GLOBAL_BINARY_DIR}/third_party/geogram)
    set(geoplatform Win-vs-dynamic-generic)
    # Extra lib
    set (EXTRA_LIBS ${EXTRA_LIBS} psapi)
    # TODO check that it is really necessary [JP]
    add_compile_options(-DGEO_DYNAMIC_LIBS)
    add_definitions(-D_CRT_SECURE_NO_WARNINGS)
else(WIN32)
    set(GEOGRAM_PATH_BIN ${GLOBAL_BINARY_DIR}/third_party/geogram/${CMAKE_BUILD_TYPE})

    if(APPLE)
        set(geoplatform Darwin-clang-dynamic)
    else(APPLE)
    # Linux
        if(${PROPAGATE_COMPILER_TO_THIRD_PARTIES})
            if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
                message(STATUS "Using Clang compiler to compile third parties")
                set(geoplatform Linux64-clang-dynamic)
            elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
                message(STATUS "Using gcc compiler to compile third parties")
                set(geoplatform Linux64-gcc-dynamic)
            endif()
        else(${PROPAGATE_COMPILER_TO_THIRD_PARTIES})
            message(STATUS "Using gcc default compiler to compile third parties")
            set(geoplatform Linux64-gcc-dynamic)
            find_program(CMAKE_C_COMPILER NAMES $ENV{CC} gcc PATHS ENV PATH NO_DEFAULT_PATH)
            find_program(CMAKE_CXX_COMPILER NAMES $ENV{CXX} g++ PATHS ENV PATH NO_DEFAULT_PATH)
        endif()
    endif(APPLE)
endif(WIN32)

# Special case for emscripten
if(EMSCRIPTEN_DIRECTORY)
    set(geoplatform Emscripten-clang)
endif()

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
        -DGEOGRAM_WITH_LUA:BOOL=OFF
        -DGEOGRAM_WITH_TETGEN:BOOL=${RINGMESH_WITH_TETGEN} 
        -DGEOGRAM_WITH_GRAPHICS:BOOL=${RINGMESH_WITH_GRAPHICS}
        -DGEOGRAM_WITH_EXPLORAGRAM:BOOL=OFF
        -DCMAKE_BUILD_TYPE:STRING=${CMAKE_BUILD_TYPE}
        -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
        -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
        -DEMSCRIPTEN_DIR=${EMSCRIPTEN_DIRECTORY}
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

# Add geogram bin directories to the current ones.
# It would be preferable to set the imported library location [JP].
# CMAKE_CFG_INTDIR is needed for Xcode (in MacOS) because the executables
# need the complete path to geogram libraries (with the configuration:
# Release or Debug).
link_directories(${GEOGRAM_PATH_BIN}/lib/${CMAKE_CFG_INTDIR})

#------------------------------------------------------------------------------------------------
# tinyxml2
# Set the path to tinyxml2 code
set(TINYXML2_PATH ${PROJECT_SOURCE_DIR}/third_party/tinyxml2)

# tinyxml2 platform dependent settings
if(CMAKE_GENERATOR STREQUAL "Unix Makefiles")
    set(TINYXML2_PATH_BIN ${GLOBAL_BINARY_DIR}/third_party/tinyxml2/${CMAKE_BUILD_TYPE})
else(CMAKE_GENERATOR STREQUAL "Unix Makefiles")
    set(TINYXML2_PATH_BIN ${GLOBAL_BINARY_DIR}/third_party/tinyxml2)
endif(CMAKE_GENERATOR STREQUAL "Unix Makefiles")

# Define tinyxml2 as an external project that we know how to
# configure and compile
# Define Tinyxml2 as an external project that we know how to
# configure and compile
ExternalProject_Add(tinyxml2_ext
  PREFIX ${TINYXML2_PATH_BIN}

  #--Download step--------------
  DOWNLOAD_COMMAND ""

  #--Update/Patch step----------
  UPDATE_COMMAND ""

  #--Configure step-------------
  SOURCE_DIR ${TINYXML2_PATH}
      CONFIGURE_COMMAND ${CMAKE_COMMAND} ${TINYXML2_PATH}
          -G ${CMAKE_GENERATOR} 
          -DCMAKE_BUILD_TYPE:STRING=${CMAKE_BUILD_TYPE}
          -DBUILD_TESTING:BOOL=OFF
          -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
          -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}

  #--Build step-----------------
  BINARY_DIR ${TINYXML2_PATH_BIN}
  #-- Command to build tinyxml2
  BUILD_COMMAND ${CMAKE_COMMAND} --build ${TINYXML2_PATH_BIN} ${COMPILATION_OPTION}

  #--Install step---------------
  INSTALL_COMMAND ""
)

ExternalProject_Add_Step(tinyxml2_ext forcebuild
    DEPENDERS build
    ALWAYS 1
  )

# Add tinyxml2 include directories to the current ones
include_directories(SYSTEM ${PROJECT_SOURCE_DIR}/third_party)

# Add tinyxml2 project libs to the libs with which RINGMesh will link
set(EXTRA_LIBS ${EXTRA_LIBS} debug tinyxml2d optimized tinyxml2)

# Add tinyxml2 bin directories to the current ones 
# It would be preferable to set the imported library location [JP]
link_directories(${TINYXML2_PATH_BIN})

#------------------------------------------------------------------------------------------------
# zlib
# Set the path to zlib code
set(ZLIB_PATH ${PROJECT_SOURCE_DIR}/third_party/zlib)

# zib platform dependent settings
if(CMAKE_GENERATOR STREQUAL "Unix Makefiles")
    set(ZLIB_PATH_BIN ${GLOBAL_BINARY_DIR}/third_party/zlib/${CMAKE_BUILD_TYPE})
else(CMAKE_GENERATOR STREQUAL "Unix Makefiles")
    set(ZLIB_PATH_BIN ${GLOBAL_BINARY_DIR}/third_party/zlib)
endif(CMAKE_GENERATOR STREQUAL "Unix Makefiles")

# Define zlib as an external project that we know how to
# configure and compile
ExternalProject_Add(zlib_ext
  PREFIX ${ZLIB_PATH_BIN}

  #--Download step--------------
  DOWNLOAD_COMMAND ""

  #--Update/Patch step----------
  UPDATE_COMMAND ""
  #--Configure step-------------
  SOURCE_DIR ${ZLIB_PATH}
      CONFIGURE_COMMAND ${CMAKE_COMMAND} ${ZLIB_PATH}
          -G ${CMAKE_GENERATOR}
          -DCMAKE_INSTALL_PREFIX:PATH=${ZLIB_PATH_BIN}/install
          -DCMAKE_BUILD_TYPE:STRING=${CMAKE_BUILD_TYPE}
          -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
          -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}

  #--Build step-----------------
  BINARY_DIR ${ZLIB_PATH_BIN}
  #-- Command to build zlib
  BUILD_COMMAND ${CMAKE_COMMAND} --build ${ZLIB_PATH_BIN} ${COMPILATION_OPTION}

  #--Install step---------------
  INSTALL_DIR ${ZLIB_PATH_BIN}/install
)

ExternalProject_Add_Step(zlib_ext forcebuild
    DEPENDERS build
    ALWAYS 1
  )

# Add zlib include directories to the current ones
include_directories(SYSTEM ${ZLIB_PATH_BIN}/install/include)

# Add zlib project libs to the libs with which RINGMesh will link
if(UNIX)
    set(EXTRA_LIBS ${EXTRA_LIBS} z)
else()
    set(EXTRA_LIBS ${EXTRA_LIBS} debug zlibd optimized zlib)
endif()


# Add zlib bin directories to the current ones
# It would be preferable to set the imported library location [JP]
link_directories(${ZLIB_PATH_BIN}/install/lib)


#------------------------------------------------------------------------------------------------
# minizip
# Set the path to minizip code
set(MINIZIP_PATH ${PROJECT_SOURCE_DIR}/third_party/minizip)

# zib platform dependent settings
if(CMAKE_GENERATOR STREQUAL "Unix Makefiles")
    set(MINIZIP_PATH_BIN ${GLOBAL_BINARY_DIR}/third_party/minizip/${CMAKE_BUILD_TYPE})
else(CMAKE_GENERATOR STREQUAL "Unix Makefiles")
    set(MINIZIP_PATH_BIN ${GLOBAL_BINARY_DIR}/third_party/minizip)
endif(CMAKE_GENERATOR STREQUAL "Unix Makefiles")

# Define minizip as an external project that we know how to
# configure and compile
ExternalProject_Add(minizip_ext
  PREFIX ${MINIZIP_PATH_BIN}

  #--Download step--------------
  DOWNLOAD_COMMAND ""

  #--Update/Patch step----------
  UPDATE_COMMAND ""

  #--Configure step-------------
  SOURCE_DIR ${MINIZIP_PATH}
      CONFIGURE_COMMAND ${CMAKE_COMMAND} ${MINIZIP_PATH}
          -G ${CMAKE_GENERATOR}
          -DCMAKE_BUILD_TYPE:STRING=${CMAKE_BUILD_TYPE}
          -DUSE_AES:BOOL=OFF
          -DZLIB_ROOT:PATH=${ZLIB_PATH_BIN}/install
          -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
          -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}

  #--Build step-----------------
  BINARY_DIR ${MINIZIP_PATH_BIN}
  #-- Command to build minizip
  BUILD_COMMAND ${CMAKE_COMMAND} --build ${MINIZIP_PATH_BIN} ${COMPILATION_OPTION}

  #--Install step---------------
  INSTALL_COMMAND ""
)

ExternalProject_Add_Step(minizip_ext forcebuild
    DEPENDERS build
    ALWAYS 1
  )

add_dependencies(minizip_ext zlib_ext)

# Add minizip include directories to the current ones
# same as minizip

# Add minizip project libs to the libs with which RINGMesh will link
set(EXTRA_LIBS ${EXTRA_LIBS} debug minizipd optimized minizip)

# Add minizip bin directories to the current ones
# It would be preferable to set the imported library location [JP]
link_directories(${MINIZIP_PATH_BIN})
