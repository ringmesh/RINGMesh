# Root CMakeList of RINGMesh project
cmake_minimum_required(VERSION 2.8)

# To configure only Debug and Release modes only
# Otherwise 2 additional configs on Windows
set(CMAKE_CONFIGURATION_TYPES 
      Debug Release CACHE CMAKE_CONFIGURATION_TYPES INTERNAL)


# Turn on the ability to create folders to organize projects and files
# It creates "CMakePredefinedTargets" folder by default and adds CMake
# defined projects like INSTALL.vcproj and ZERO_CHECK.vcproj
set_property(GLOBAL PROPERTY USE_FOLDERS ON)

# Needed cmake modules
include(ExternalProject)


# Define the project RINGMesh      
project(RINGMesh)

# Define version number 
# It is then exported to the configuration file
set (RINGMesh_VERSION_MAJOR 1)
set (RINGMesh_VERSION_MINOR 0)

# To check that a given function exists we can use 
#include(CheckFunctionExists)
#check_function_exists (log HAVE_LOG)
# and so HAVE_LOG stores true of false depending on the system

 
# Options can be set with the CMakeOptions file, but also in Cmake GUI
# or in command line cmake -D <var>:<type>=<value> MyProjectFolder

# Path to MG_TETRA is asked for to the user in the interface
set(MG_TETRA "" CACHE FILEPATH "Path to MG-Tetra")

option(USE_OPENMP "Use Open MP" ON)

# Optional components of RINGMesh - creation of specific targets
option(RINGMESH_WITH_GRAPHICS "Compile viewer" OFF)
option(RINGMESH_WITH_UTILITIES "Compile utility executables" OFF)
option(RINGMESH_WITH_TESTS "Compile test projects" OFF)

# Optional custom steps 
option(BUILD_DOCUMENTATION "Create and install the HTML documentation (requires Doxygen)")


if(UNIX)    
    if(CMAKE_BUILD_TYPE STREQUAL "")   
      foreach(i IN LISTS CMAKE_CONFIGURATION_TYPES)
          set(project_binary_dir_config  ${CMAKE_SOURCE_DIR}/build/ringmesh/${i})
          file(MAKE_DIRECTORY ${project_binary_dir_config})
          execute_process(
             COMMAND ${CMAKE_COMMAND} ${CMAKE_SOURCE_DIR}
                -DSKIP:BOOL=ON
                -DCMAKE_BUILD_TYPE=${i}
                -DPROJECT_BINARY_DIR:STRING=${project_binary_dir_config}
                -DMG_TETRA:STRING=${MG_TETRA}
                -DUSE_OPENMP:BOOL=${USE_OPENMP}
                -DRINGMESH_WITH_GRAPHICS:BOOL=${RINGMESH_WITH_GRAPHICS}
                -DRINGMESH_WITH_UTILITIES:BOOL=${RINGMESH_WITH_UTILITIES}
                -DRINGMESH_WITH_TESTS:BOOL=${RINGMESH_WITH_TESTS}
                -DBUILD_DOCUMENTATION:BOOL=${BUILD_DOCUMENTATION}
                WORKING_DIRECTORY ${project_binary_dir_config})
      endforeach()
      return()
    endif()
endif(UNIX)

# Platform dependent settings
if(WIN32)
    set(RINGMESH_EXPORTS ON)
    set(GEO_DYNAMIC_LIBS ON)
    # Extra library needed with geogram I think [JP]
    set (EXTRA_LIBS ${EXTRA_LIBS} psapi)
    set(geoplatform Win-vs-dynamic-generic)    
else(WIN32)
    set(geoplatform Linux64-gcc-dynamic)
endif(WIN32)
if(UNIX)
    # If you really need it, please comment [JP]
    # Maybe rather test the compiler than the system
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")
    set(CMAKE_C_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")
endif(UNIX)


# Configure a header file to pass some of the CMake settings
# to the source code
set(generated_configure_file ${PROJECT_BINARY_DIR}/ringmesh/ringmesh_config.h)
configure_file(
  "${PROJECT_SOURCE_DIR}/include/ringmesh/ringmesh_config.h.in"
   ${generated_configure_file}
)
# We want to be able to include this file
include_directories(${PROJECT_BINARY_DIR})
 
#------------------------------------------------------------------------------------------------
# Define Geogram as an external project that we know how to
# configure and compile

set(GEOGRAM_PATH ${PROJECT_SOURCE_DIR}/src/third_party/geogram)
# Maybe change it for Linux [JP]
if(WIN32)
    set(GEOGRAM_PATH_BIN ${PROJECT_SOURCE_DIR}/build/geogram)
else(WIN32)
    set(GEOGRAM_PATH_BIN ${PROJECT_SOURCE_DIR}/build/geogram/${CMAKE_BUILD_TYPE})
endif(WIN32)

ExternalProject_Add(geogram_ext
  PREFIX ${GEOGRAM_PATH}
  
  #--Download step--------------
  DOWNLOAD_COMMAND ""
  
  #--Update/Patch step----------
  UPDATE_COMMAND ""

  #--Configure step-------------
  SOURCE_DIR ${GEOGRAM_PATH}
  CONFIGURE_COMMAND ${CMAKE_COMMAND} ${GEOGRAM_PATH}
        -G ${CMAKE_GENERATOR} -DVORPALINE_PLATFORM:STRING=${geoplatform}
        -DGEOGRAM_WITH_TETGEN:BOOL=ON
        -DGEOGRAM_WITH_MEDIT:BOOL=ON 
        -DGEOGRAM_WITH_GRAPHICS=${RINGMESH_WITH_GRAPHICS}
        -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
  # maybe missing options for linux [JP]
  
  #--Build step-----------------
  BINARY_DIR ${GEOGRAM_PATH_BIN}  
  # Build all configurations if I understood well [JP]
  BUILD_COMMAND ${CMAKE_COMMAND} --build ${GEOGRAM_PATH_BIN}

  #--Install step---------------
  INSTALL_COMMAND "" 
)

# Add geogram include directories to the current ones
include_directories(${GEOGRAM_PATH}/src/lib)

# Add geogram project libs to the libs with which RINGMesh will link
set(EXTRA_LIBS ${EXTRA_LIBS} geogram)
if(RINGMESH_WITH_GRAPHICS)
    set(EXTRA_LIBS ${EXTRA_LIBS} geogram_gfx)
endif(RINGMESH_WITH_GRAPHICS)
    
# Add geogram bin directories to the current ones 
# It would be preferable to set the imported library location
# but I did not managed to do it [JP]
# Check this on Linux [JP]
link_directories(${GEOGRAM_PATH_BIN}/lib)


# On Windows I desperately need to copy the dll of geogram 
# with those of RINGMesh but still did not manage [JP]
# Maybe it would be cleaner with a clean installation step for geogram

#------------------------------------------------------------------------------------------------    

 # I do not use it, so please clean this section yourself [JP]
if(not ${MG_TETRA} STREQUAL "")   
  message(STATUS "Using MG-Tetra directory =  ${MG_TETRA}")
  include_directories(${MG_TETRA}/include)
 
  if(CMAKE_SIZEOF_VOID_P EQUAL 8)
    set(PLATFORM _64)
  else(CMAKE_SIZEOF_VOID_P EQUAL 8)
    SET(PLATFORM)
  endif(CMAKE_SIZEOF_VOID_P EQUAL 8)
  
  # Is there not a nice way to import this library [JP]
  if(WIN32)
    link_directories(${MG_TETRA}/lib/WIN7${PLATFORM}_VC10)
  else(WIN32)
    link_directories(${MG_TETRA}/lib/Linux${PLATFORM})
  endif(WIN32)
  
  # Additional libraries with which RINGMesh should link
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

# Collect files of the project
file(GLOB ringmesh_src "${PROJECT_SOURCE_DIR}/src/ringmesh/*.cpp")
file(GLOB ringmesh_include "${PROJECT_SOURCE_DIR}/include/ringmesh/*.h")
file(GLOB zlib_src "${PROJECT_SOURCE_DIR}/src/third_party/zlib/*.c")
file(GLOB zlib_include "${PROJECT_SOURCE_DIR}/include/third_party/zlib/*.h")

# Add configure file to ringmesh includes
set(ringmesh_include ${ringmesh_include} ${generated_configure_file})

# Group files of a project in directories
# For this to work all the files should be added to the sources of the project
# and the global property USE_FOLDERS must be ON
source_group("Zlib" FILES ${zlib_src} ${zlib_include})
source_group("Source Files" FILES ${ringmesh_src})
source_group("Header Files" FILES ${ringmesh_include})

#------------------------------------------------------------------------------------------------
# Build configuration

# Organize cmake outputs in a Linux way
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)

# Include directories of RINGMesh
include_directories(${PROJECT_SOURCE_DIR}/include)

# Add the RINGMesh target as a shared library
add_library(RINGMesh SHARED
    ${ringmesh_include} 
    ${ringmesh_src} 
    ${zlib_src}
    ${zlib_include})

#RINGMesh depends on Geogram. Modification of Geogram -> modification of RINGMesh
#This does not happen ? snif snif...
add_dependencies(RINGMesh geogram_ext)

if(WIN32)
    # Optional copy of RINGMesh dll in a dependency 
    add_custom_command(TARGET RINGMesh PRE_BUILD 
     COMMAND  "${CMAKE_COMMAND}" -E copy_if_different 
         "${GEOGRAM_PATH_BIN}/bin/$<CONFIGURATION>/geogram.dll"
         "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/$<CONFIGURATION>/geogram.dll" 
     COMMENT "Copying to output directory")
endif()

# Add all libraries to which RINGMesh should link
target_link_libraries(RINGMesh ${EXTRA_LIBS})


# Optional modules configuration
if(RINGMESH_WITH_GRAPHICS)
    add_subdirectory(src/bin/ringmeshview)
endif(RINGMESH_WITH_GRAPHICS)

if(RINGMESH_WITH_UTILITIES)
    add_subdirectory(src/bin/ringmeshconvert)
endif(RINGMESH_WITH_UTILITIES)

if(RINGMESH_WITH_TESTS)
    add_subdirectory(tests)
endif(RINGMESH_WITH_TESTS)

#------------------------------------------------------------------------------------------------


if(BUILD_DOCUMENTATION)
    find_package(Doxygen)
    if(NOT DOXYGEN_FOUND)
        message(FATAL_ERROR "Doxygen is needed to build the documentation.")
    endif()
    
    set(doxyfile ${PROJECT_SOURCE_DIR}/doc/Doxyfile)

    # To set properly project version and paths
    # use an input file and use configure_file to populate the doxyfile
    # set(doxyfile_in ${PROJECT_SOURCE_DIR}/Doxyfile.in)
    # configure_file(${doxyfile_in} ${doxyfile} @ONLY)  
    
    # Create a target that can be build from the IDE 
    add_custom_target(doc
        COMMAND ${DOXYGEN_EXECUTABLE} ${doxyfile}
        WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/doc
        COMMENT "Generating API documentation with Doxygen"
        VERBATIM)

    # We should think about installation steps fo the whole project [JP]
    # install(DIRECTORY ${PROJECT_SOURCE_DIR}/doc/html DESTINATION share/doc)
endif()
