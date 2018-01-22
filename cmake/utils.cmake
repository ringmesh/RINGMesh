# Copyright (c) 2012-2018, Association Scientifique pour la Geologie et ses
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

function(add_ringmesh_library directory)
    string(REPLACE "/" "_" target_name ${directory})
    add_library(${target_name} SHARED "")
    set_target_properties(${target_name} 
        PROPERTIES 
            OUTPUT_NAME RINGMesh_${target_name}
            FOLDER "Libraries"
    )
    target_include_directories(${target_name} 
        PUBLIC   
            $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
            $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}>
            $<INSTALL_INTERFACE:include>
    )
    target_link_libraries(${target_name} PUBLIC Geogram::geogram)
    if(WIN32)
        target_compile_definitions(${target_name} 
            PUBLIC 
                -DGEO_DYNAMIC_LIBS 
                # Following are meant for geogram
                -D_CRT_SECURE_NO_WARNINGS
        )
        add_dependencies(copy_dll ${target_name})
    endif()
    set(lib_include_dir ${PROJECT_SOURCE_DIR}/include/ringmesh/${directory})
    set(lib_source_dir ${PROJECT_SOURCE_DIR}/src/ringmesh/${directory})
    include(${PROJECT_SOURCE_DIR}/src/ringmesh/${directory}/CMakeLists.txt)
    
    export(TARGETS ${target_name} 
        NAMESPACE RINGMesh:: 
        FILE lib/cmake/RINGMesh/RINGMesh_${target_name}_target.cmake
    )
    generate_export_header(${target_name} 
        EXPORT_MACRO_NAME ${target_name}_api 
        EXPORT_FILE_NAME ${PROJECT_BINARY_DIR}/ringmesh/${directory}/export.h
    )
    install(FILES ${PROJECT_BINARY_DIR}/ringmesh/${directory}/export.h
        DESTINATION include/ringmesh/${directory}
    )
    install(TARGETS ${target_name} 
        EXPORT ${target_name}
        RUNTIME DESTINATION bin
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib
    )
    install(EXPORT ${target_name}
        FILE RINGMesh_${target_name}_target.cmake
        NAMESPACE RINGMesh::
        DESTINATION lib/cmake/RINGMesh
    )
endfunction()

function(add_js_target target src)
    if(NOT RINGMESH_WITH_GUI)
        return()
    endif()
    set(js_target ${target}_js)
    set(NBIND_FILE ${NBIND_DIR}/src/v8/Binding.cc)
    add_nodejs_module(${js_target} ${src} ${NBIND_FILE})
    set(target_node_name RINGMesh_${target_name})
    configure_file(${PROJECT_BINARY_DIR}/nbind.js.in
        ${PROJECT_BINARY_DIR}/lib/${target_node_name}.js)
    set_target_properties(${js_target}
        PROPERTIES 
            OUTPUT_NAME ${target_node_name}
    )
    target_link_libraries(${js_target} LINK_PUBLIC ${target} nbind)
endfunction()

macro(copy_for_windows directory)

    # On windows, without proper installation steps, we need to
    # copy of Geogram dll and pdb information to RINGMesh
    # to be able to launch RINGMesh utilities and tests from the debugger

    # The dll and debug info of RINGMesh are in
    # build/ringmesh/Debug or build/ringmesh/Release.
if(WIN32)
    add_custom_command(TARGET copy_dll POST_BUILD
        COMMAND  "${CMAKE_COMMAND}" -E copy_directory
            "${PROJECT_BINARY_DIR}/$<CONFIGURATION>"
            "${directory}/$<CONFIGURATION>"
            COMMENT "Copy RINGMesh dll")
    add_custom_command(TARGET copy_dll POST_BUILD
        COMMAND  "${CMAKE_COMMAND}" -E copy_directory
            "${GEOGRAM_INSTALL_PREFIX}/bin"
            "${directory}/$<CONFIGURATION>"
            COMMENT "Copy geogram binaries")
    add_custom_command(TARGET copy_dll POST_BUILD
        COMMAND  "${CMAKE_COMMAND}" -E copy_directory
            "${GEOGRAM_INSTALL_PREFIX}/lib"
            "${directory}/$<CONFIGURATION>"
            COMMENT "Copy geogram visualization libraries")
    add_custom_command(TARGET copy_dll POST_BUILD
        COMMAND  "${CMAKE_COMMAND}" -E copy_directory
            "${ZLIB_ROOT}/bin"
            "${directory}/$<CONFIGURATION>"
            COMMENT "Copy zlib binaries")
    add_custom_command(TARGET copy_dll POST_BUILD
        COMMAND  "${CMAKE_COMMAND}" -E copy_directory
            "${TINYXML2_INSTALL_PREFIX}/bin"
            "${directory}/$<CONFIGURATION>"
            COMMENT "Copy tinyxml2 binaries")
endif(WIN32)
endmacro()

macro(add_ringmesh_executable exe_path folder_name)
    get_filename_component(exe_name ${exe_path} NAME_WE)

    # Set the target as an executable
    add_executable(${exe_name} ${exe_path})    
    foreach(dependency ${ARGN})
        target_link_libraries(${exe_name} PRIVATE ${dependency})
    endforeach()
    
    # Add the project to a folder of projects for the tests
    if(APPLE)
        set(OS_RPATH "@executable_path")
    else()
        set(OS_RPATH "$ORIGIN")
    endif()
    set_target_properties(${exe_name} 
        PROPERTIES 
            FOLDER ${folder_name}
            INSTALL_RPATH "${OS_RPATH}/../lib"
    )
endmacro()

function(add_ringmesh_binary bin_path)
    add_ringmesh_executable(${bin_path} "Utilities" ${ARGN})
    install(TARGETS ${exe_name} RUNTIME DESTINATION bin)
    set_target_properties(${exe_name} PROPERTIES 
        RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)
endfunction()

function(add_ringmesh_test cpp_file_path)
    add_ringmesh_executable(${cpp_file_path} "Tests" ${ARGN})
    # Add the test to CTest
    add_test(NAME ${exe_name} COMMAND ${exe_name})
    set_target_properties(${exe_name} PROPERTIES 
        RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin/tests)
endfunction()


function(add_ringmesh_tutorial cpp_file_path)
    add_ringmesh_executable(${cpp_file_path} "Tutorials" ${ARGN})
    # Add the test to CTest
    add_test(NAME ${exe_name} COMMAND ${exe_name})
    set_target_properties(${exe_name} PROPERTIES 
        RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin/tutorials)
endfunction()
