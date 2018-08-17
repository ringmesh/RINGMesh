# Define the project RINGMesh
project(RINGMesh)

#------------------------------------------------------------------------------------------------
# Turn on the ability to create folders to organize projects and files
# It creates "CMakePredefinedTargets" folder by default and adds CMake
# defined projects like INSTALL.vcproj and ZERO_CHECK.vcproj
set_property(GLOBAL PROPERTY USE_FOLDERS ON)

# Define version number
# It is then exported to the configuration file
set (RINGMesh_VERSION_MAJOR 5)
set (RINGMesh_VERSION_MINOR 1)
set (RINGMesh_VERSION_PATCH 0)
set(RINGMesh_VERSION ${RINGMesh_VERSION_MAJOR}.${RINGMesh_VERSION_MINOR}.${RINGMesh_VERSION_PATCH})

message(STATUS "RINGMesh binary directory is: ${PROJECT_BINARY_DIR}")
message(STATUS "RINGMesh source directory is: ${PROJECT_SOURCE_DIR}")

include(cmake/utils.cmake)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
    
#------------------------------------------------------------------------------------------------
# Platform dependent settings
if(UNIX)
    add_compile_options(-Wall -Wextra -Wno-long-long -Wconversion
        -Wsign-conversion -Wdouble-promotion -Wno-attributes)

    if(APPLE)
        if (NOT "${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
            message(WARNING "RINGMesh on Apple is only tested with Clang compiler")
        endif()
    else(APPLE)
        # pthread is ignored on MacOS
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread")
    endif(APPLE)

    if(CMAKE_BUILD_TYPE STREQUAL "Coverage")
       include(tools/Coverage.cmake)
    endif()
else(UNIX)
    if(NOT "${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
        message(WARNING "RINGMesh on Windows is only tested with Microsoft Visual C++")
    endif()
    add_custom_target(copy_dll ALL)
endif(UNIX)

# RINGMesh depends on Geogram Tinyxml2 Minizip and zlib

#list(APPEND CMAKE_MODULE_PATH ${GEOGRAM_INSTALL_PREFIX}/lib/cmake/modules)
#set(VORPALINE_BUILD_DYNAMIC TRUE CACHE BOOL "")
#find_package(Geogram REQUIRED)


# already done in RINGMeshConfig.cmake.in
# Path to geogram
#################

set(VORPALINE_BUILD_DYNAMIC TRUE CACHE BOOL "")

if(IS_DIRECTORY ${CMAKE_SOURCE_DIR}/geogram/)
   set(
      GEOGRAM_SOURCE_DIR "${CMAKE_SOURCE_DIR}/geogram/"
      CACHE PATH "full path to the Geogram (or Vorpaline) installation"
   )
   set(USE_BUILTIN_GEOGRAM TRUE)
else()
   set(
      GEOGRAM_SOURCE_DIR "${CMAKE_SOURCE_DIR}/../geogram/"
      CACHE PATH "full path to the Geogram (or Vorpaline) installation"
   )
   set(USE_BUILTIN_GEOGRAM FALSE)
endif()

message(STATUS "GEOGRAM Path ${GEOGRAM_SOURCE_DIR}")
include(${GEOGRAM_SOURCE_DIR}/cmake/geogram.cmake)


find_package(ZLIB REQUIRED)
find_package(tinyxml2 REQUIRED PATHS ${TINYXML2_INSTALL_PREFIX})
include(${MINIZIP_INSTALL_PREFIX}/cmake/minizip-exports.cmake)

install(
    DIRECTORY
        ${TINYXML2_INSTALL_PREFIX}/
        ${ZLIB_ROOT}/
#        ${GEOGRAM_INSTALL_PREFIX}/
    DESTINATION
        .
)

include(GenerateExportHeader)

if(RINGMESH_WITH_GUI)
    message(STATUS "Configure RINGMesh with GUI")
    find_program(NPM NAMES npm.cmd npm)
    if(NOT NPM)
        message(FATAL_ERROR "npm is needed to create the GUI")
    endif()
    execute_process(COMMAND ${NPM} install)
    include(node_modules/node-cmake/NodeJS.cmake)
    nodejs_init()
    
    function(add_nodejs_module NAME)
        _add_nodejs_module(${NAME} ${ARGN})
        set_target_properties(${NAME} 
            PROPERTIES 
                C_VISIBILITY_PRESET default
                CXX_VISIBILITY_PRESET default
        )
    endfunction()
    
    set(NBIND_DIR ${PROJECT_SOURCE_DIR}/third_party/nbind)
    set(NBIND_SOURCE_FILES
        ${NBIND_DIR}/src/common.cc
        ${NBIND_DIR}/src/reflect.cc
        ${NBIND_DIR}/src/v8/Binding.cc
        ${NBIND_DIR}/src/v8/Buffer.cc
    )
    add_nodejs_module(nbind ${NBIND_SOURCE_FILES})
    target_include_directories(nbind 
        SYSTEM PUBLIC ${NBIND_DIR}/include)
    target_compile_definitions(nbind
        PUBLIC
            -DBUILDING_NODE_EXTENSION
            -DUSING_V8_SHARED
            -DUSING_UV_SHARED
    )
    generate_export_header(nbind
        EXPORT_MACRO_NAME nbind_api 
        EXPORT_FILE_NAME ${NBIND_DIR}/include/nbind/export.h
    )
endif()

#------------------------------------------------------------------------------------------------
# Build configuration
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
#------------------------------------------------------------------------------------------------
# Install configuration
set(CMAKE_MACOSX_RPATH ON)
set(CMAKE_INSTALL_RPATH ".")
# add the automatically determined parts of the RPATH
# which point to directories outside the build tree to the install RPATH
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

#------------------------------------------------------------------------------------------------
# file automatically generated by cmake

# generate configure file to pass on some of the CMake settings
# to the source code
configure_file(
    ${PROJECT_SOURCE_DIR}/cmake/ringmesh_config.h.in 
    ${PROJECT_BINARY_DIR}/ringmesh/ringmesh_config.h
)
install(
    FILES ${PROJECT_BINARY_DIR}/ringmesh/ringmesh_config.h 
    DESTINATION include/ringmesh
)

# Exports RINGMesh target
include(CMakePackageConfigHelpers)
include(InstallRequiredSystemLibraries)
configure_package_config_file(
    cmake/RINGMeshConfig.cmake.in 
    ${CMAKE_BINARY_DIR}/lib/cmake/RINGMesh/RINGMeshConfig.cmake
    INSTALL_DESTINATION lib/cmake/RINGMesh
    PATH_VARS GEOGRAM_INSTALL_PREFIX
)
install(
    FILES ${CMAKE_BINARY_DIR}/lib/cmake/RINGMesh/RINGMeshConfig.cmake 
    DESTINATION lib/cmake/RINGMesh
)
install(DIRECTORY include/ringmesh DESTINATION include)

#------------------------------------------------------------------------------------------------
# Configure the ringmesh libraries

add_ringmesh_library(basic)
add_ringmesh_library(geogram_extension)
add_ringmesh_library(geomodel/builder)
add_ringmesh_library(geomodel/core)
add_ringmesh_library(geomodel/tools)
add_ringmesh_library(io)
add_ringmesh_library(mesh)
add_ringmesh_library(tetrahedralize)
if(RINGMESH_WITH_GRAPHICS)
    add_ringmesh_library(visualize)
endif(RINGMESH_WITH_GRAPHICS)

#Copy dll from RINGMesh third parties to make it accessible for plugins
copy_deps_dll_window()

#------------------------------------------------------------------------------------------------
# Optional modules configuration

set(binary_source_dir ${PROJECT_SOURCE_DIR}/src/bin)
if(BUILD_RINGMESH_VIEW)
    message(STATUS "Configure ringmesh-view")
    add_ringmesh_binary(${binary_source_dir}/ringmesh-view.cpp visualize)
    copy_for_all_ringmesh_dlls(${PROJECT_BINARY_DIR}/bin)
endif()

if(RINGMESH_WITH_UTILITIES)
    message(STATUS "Configuring RINGMesh with utilities")
    # Get the paths of the utility files
    file(GLOB utility_sources "${binary_source_dir}/utilities/*.cpp")
    foreach(utility_src ${utility_sources})
        add_ringmesh_binary(${utility_src} geomodel_tools io)
    endforeach()
    copy_for_all_ringmesh_dlls(${PROJECT_BINARY_DIR}/bin)
endif()

if(RINGMESH_WITH_TUTORIALS)
    message(STATUS "Configuring RINGMesh with tutorials")
    add_subdirectory(doc/tutorials)
    copy_for_all_ringmesh_dlls(${PROJECT_BINARY_DIR}/bin/tutorials)
endif()

if(RINGMESH_WITH_TESTS)
    # Enable testing with CTest
    enable_testing()
    message(STATUS "Configuring RINGMesh with tests")
    add_subdirectory(tests)
    copy_for_all_ringmesh_dlls(${PROJECT_BINARY_DIR}/bin/tests)
endif()

#------------------------------------------------------------------------------------------------
# Configure CPack

set(CPACK_PACKAGE_NAME ${CMAKE_PROJECT_NAME})
set(CPACK_PACKAGE_VERSION_MAJOR ${RINGMesh_VERSION_MAJOR})
set(CPACK_PACKAGE_VERSION_MINOR ${RINGMesh_VERSION_MINOR})
set(CPACK_PACKAGE_VERSION_PATCH ${RINGMesh_VERSION_PATCH})
set(CPACK_PACKAGE_VENDOR "RING-TEAM (www.ring-team.org)")

set(CPACK_SOURCE_GENERATOR "ZIP")
set(CPACK_SOURCE_IGNORE_FILES "/build/;/.git/;/_CPack_Packages/")

if(WIN32)
    set(CPACK_GENERATOR "ZIP")
else()
    set(CPACK_GENERATOR "TGZ")
endif()


# This must always be last!
include(CPack)
