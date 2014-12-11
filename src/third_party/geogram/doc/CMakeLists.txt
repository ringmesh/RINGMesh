#
# Build rules for the documentation
#

#! The minimum doxygen version number allowed here
set(DOXYGEN_MINIMUM_VERSION 1.7.0)

# Check for doxygen and dot

find_package(Doxygen ${DOXYGEN_MINIMUM_VERSION} QUIET)
if(NOT DOXYGEN_FOUND)
    message(WARNING "Doxygen >= ${DOXYGEN_MINIMUM_VERSION} not found, cannot generate documentation")
    unset(DOXYGEN_EXECUTABLE CACHE)
    return()
endif()

if(DOXYGEN_DOT_FOUND)
    set(DOXYGEN_HAVE_DOT YES)
else()
    set(DOXYGEN_HAVE_DOT NO)
endif()


# Generate a target for generating a specific documentation type

function(add_doc_target doc_type)

    set(doc_output_dir ${CMAKE_CURRENT_BINARY_DIR})

    configure_file(${doc_type}.dox.in ${doc_type}.dox @ONLY)

    add_custom_target(
        doc-${doc_type}
        COMMAND ${CMAKE_COMMAND} -E remove_directory ${doc_output_dir}/${doc_type}
        COMMAND ${DOXYGEN_EXECUTABLE} ${doc_type}.dox # > ${doc_type}.log 2>&1
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
        COMMENT "Generate ${doc_type} documentation for ${PROJECT_NAME}"
    )

    set_property(
        DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
        APPEND PROPERTY ADDITIONAL_MAKE_CLEAN_FILES ${doc_output_dir}/${doc_type}
    )

endfunction()


# Add a custom documentation targets

add_doc_target(devkit)
add_doc_target(devkit-full)
add_doc_target(devkit-internal)
add_doc_target(devkit-internal-light)

# TODO: add end-user documentation

# Copy the Version info file to the build directory
configure_file(VERSION.txt.in VERSION.txt @ONLY)

# Install documentation

install(FILES README.txt DESTINATION doc COMPONENT runtime)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/LICENSE.txt DESTINATION doc COMPONENT runtime OPTIONAL)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/VERSION.txt DESTINATION doc COMPONENT runtime OPTIONAL)
install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/devkit/html DESTINATION doc/devkit COMPONENT doc-devkit OPTIONAL)
install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/devkit-full/html DESTINATION doc/devkit COMPONENT doc-devkit-full OPTIONAL)
install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/devkit-internal/html DESTINATION doc/devkit COMPONENT doc-devkit-internal OPTIONAL)

