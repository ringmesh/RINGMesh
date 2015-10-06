cmake_minimum_required(VERSION 2.8)

project(RINGMesh)
# include some cmake modules
include(ExternalProject)


set(CMAKE_CONFIGURATION_TYPES Debug Release)

# define version number 
set (RINGMesh_VERSION_MAJOR 1)
set (RINGMesh_VERSION_MINOR 0)

# To check that a given function exists we can use 
#include (${CMAKE_ROOT}/Modules/CheckFunctionExists.cmake)
#check_function_exists (log HAVE_LOG)
# and so HAVE_LOG stores true of false depending on the system


# Turn on the ability to create folders to organize projects (.vcproj)
# It creates "CMakePredefinedTargets" folder by default and adds CMake
# defined projects like INSTALL.vcproj and ZERO_CHECK.vcproj
set_property(GLOBAL PROPERTY USE_FOLDERS ON)

 
# Options can be set with the CMakeOptions file, but also in Cmake GUI
# or in command line cmake -D <var>:<type>=<value> MyProjectFolder
option(USE_MG_TETRA "Compile with mesher MG-Tetra" OFF)
option(MG_TETRA "Path to MG-Tetra" "")

option(RINGMESH_WITH_GRAPHICS "Compile viewer" OFF)
option(RINGMESH_WITH_UTILITIES "Compile utility executables" OFF)
option(RINGMESH_WITH_TESTS "Compile test projects" OFF)

# To disable OpenMP
option(USE_OPENMP "Use Open MP" OFF)
 
# if(WIN32)
    # add_definitions(-DRINGMESH_EXPORTS -DGEO_DYNAMIC_LIBS)
# endif(WIN32)

# Configure a header file to pass some of the CMake settings
# to the source code
configure_file (
  "${PROJECT_SOURCE_DIR}/include/ringmesh/ringmesh_config.h.in"
  "${PROJECT_BINARY_DIR}/ringmesh_config.h"
)
  
  
##############################################################################
# Read a file for the global setting of the options ??
# Do we really need it ?
# if(EXISTS ${CMAKE_SOURCE_DIR}/CMakeOptions.txt)
   # MESSAGE(INFO: "Using local options file: " .${CMAKE_SOURCE_DIR}/CMakeOptions.txt)
   # INCLUDE(${CMAKE_SOURCE_DIR}/CMakeOptions.txt)
# endif(EXISTS ${CMAKE_SOURCE_DIR}/CMakeOptions.txt)


##############################################################################

# Let's define properly that Geogram is another project we rely on
# A FINIR  PAS FINI
ExternalProject_Add(geogram_ext
  PREFIX ${PROJECT_SOURCE_DIR}/src/third_party/geogram
  DOWNLOAD_COMMAND ""
  UPDATE_COMMAND ""
  
  SOURCE_DIR ${PROJECT_SOURCE_DIR}/src/third_party/geogram/src/lib  
  CONFIGURE_COMMAND ${CMAKE_COMMAND} . -DVORPALINE_PLATFORM:STRING=Win-vs-dynamic-generic
  
  BINARY_DIR ${PROJECT_SOURCE_DIR}/build/geogram
  
  BUILD_COMMAND ${CMAKE_MAKE_PROGRAM} 
  
  INSTALL ""
)

# Je voudrai copier la lib de geogram là ou je veux
# ou alors step INSTALL ? ça doit marcher ....
# ExternalProject_Add_Step(geogram_ext Copying
  # COMMAND  "${CMAKE_COMMAND}" -E copy_if_different  $<TARGET_FILE:geogram_ext>
  # ${PROJECT_BINARY_DIR}/$<CONFIGURATION>/$<TARGET_FILE_NAME:geogram_ext>
# )

# besoin de ça ?
# add_library(geogram SHARED IMPORTED)

ExternalProject_Get_Property(geogram_ext SOURCE_DIR)
include_directories(${SOURCE_DIR})


# Add Geogram to the libs with which RINGMesh 
# will be linked ? ça marche ou il faut être plus précis sur le chemin ?
set (EXTRA_LIBS ${EXTRA_LIBS} geogram)

# on devrait les avoir automatiquement si on fait ça correctement
# if(WIN32)
    # link_directories(${CMAKE_SOURCE_DIR}/build/geogram/${PLATFORM}/lib)
# ELSE(WIN32)
    # link_directories(${CMAKE_SOURCE_DIR}/build/geogram/${PLATFORM}-${CMAKE_BUILD_TYPE}/lib)
# ENDIF(WIN32)


if(USE_MG_TETRA)
  if(${MG_TETRA} MATCHES .+)
    message(STATUS "Using MG-Tetra directory =  ${MG_TETRA}")
    set(MG_TETRA_DIR  ${MG_TETRA})
  elseif( $ENV{MG_TETRA} MATCHES .+)
    message(STATUS "Using MG-Tetra directory =  $ENV{MG_TETRA}")
    set(MG_TETRA_DIR  $ENV{MG_TETRA})
  else(${MG_TETRA} MATCHES .+)
    message(STATUS "No MG-Tetra directory found")
  endif(${MG_TETRA} MATCHES .+)
  include_directories(${MG_TETRA_DIR}/include)
 
  if( CMAKE_SIZEOF_VOID_P EQUAL 8 )
    set(PLATFORM _64)
  else(CMAKE_SIZEOF_VOID_P EQUAL 8 )
    SET(PLATFORM)
  endif(CMAKE_SIZEOF_VOID_P EQUAL 8 )
  
  # Mouais NON trouver mieux - importer cette librairie ? comment ? 
  if(WIN32)
    link_directories(${MG_TETRA_DIR}/lib/WIN7${PLATFORM}_VC10)
  else(WIN32)
    link_directories(${MG_TETRA_DIR}/lib/Linux${PLATFORM})
  endif(WIN32)
  
  set (EXTRA_LIBS ${EXTRA_LIBS} meshgems mg-tetra)
endif(USE_MG_TETRA)


if(RINGMESH_WITH_GRAPHICS)
    # faire pareil qu'avec geogram - definir le projet ? recuperer la target ?
    add_subdirectory(src/bin/ringmeshview)
    set (EXTRA_LIBS ${EXTRA_LIBS} geogram_gfx)
endif(RINGMESH_WITH_GRAPHICS)

if(RINGMESH_WITH_UTILITIES)
    add_subdirectory(src/bin/ringmeshconvert)
endif(RINGMESH_WITH_UTILITIES)

if(RINGMESH_WITH_TESTS)
    add_subdirectory(tests)
endif(RINGMESH_WITH_TESTS)


if(USE_OPENMP)
# Vérifier ça ?? vraiment utile ?
    if(WIN32) 
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /wd4275 /wd4996 /MP")
        set(CMAKE_C_FLAGS "${CMAKE_CXX_FLAGS} /wd4275 /wd4996 /MP")
    else(WIN32)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
        set(CMAKE_C_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
    endif(WIN32)
endif(USE_OPENMP)



# CMake documentation states "We do not recommend using GLOB to collect a list of
# source files from your source tree. If no CMakeLists.txt file changes when a 
# source is added or removed then the generated build system cannot know when to 
# ask CMake to regenerate".
file(GLOB ringmesh_src "${PROJECT_SOURCE_DIR}/src/ringmesh/*.cpp")
file(GLOB ringmesh_include "${PROJECT_SOURCE_DIR}/include/ringmesh/*.h")
file(GLOB zlib_src "${PROJECT_SOURCE_DIR}/src/third_party/zlib/*.c")
file(GLOB zlib_include "${PROJECT_SOURCE_DIR}/include/third_party/zlib/*.h")


# Group files of a project in directories
# For this to work all the files should be added to the sources of the project
source_group("Zlib" FILES ${zlib_src} ${zlib_include})
source_group("Source Files" FILES ${ringmesh_src})
source_group("Header Files" FILES ${ringmesh_include})

# tiens on peut faire des dossiers de vcproj
# il nous en faut un pour les tests, ce sera plus joli
# Creates folder "libraries" and adds target project (math.vcproj)
# set_property(TARGET math PROPERTY FOLDER "libraries")


if(UNIX)
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")
    SET(CMAKE_C_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")
endif(UNIX)
if(WIN32)
    set (EXTRA_LIBS ${EXTRA_LIBS} psapi)
endif(WIN32)

##############################################################################
# Build configuration

# set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)
# set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
# set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)


include_directories(${PROJECT_SOURCE_DIR}/include)

add_library(RINGMesh SHARED
    ${ringmesh_include} 
    ${ringmesh_src} 
    ${zlib_src}
    ${zlib_include}
  )


#RINGMesh depends on Geogram. Modification of Geogram -> modification of RINGMesh
add_dependencies(RINGMesh geogram_ext)


# Add all libraries to which RINGMesh should link
target_link_libraries(RINGMesh ${EXTRA_LIBS})



# Check how it is possible to launch Doxygen documentation generation
# It seems to be possible

############################################################################
# On Windows we need to copy Geogram.dll in the RINGMesh bin directories
# to launch RINGMesh executables

# Optional copy of RINGMesh dll in a dependency 


