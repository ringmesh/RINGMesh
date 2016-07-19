# CMakeList for the documentation generation of the RINGMesh project
# Author 

# Requested minimal doxygen version 
set(DOXYGEN_MINIMUM_VERSION 1.6.0)

# Check for doxygen and dot

find_package(Doxygen ${DOXYGEN_MINIMUM_VERSION} QUIET)
if(NOT DOXYGEN_FOUND)
    message(WARNING "Doxygen >= ${DOXYGEN_MINIMUM_VERSION} not found, cannot generate documentation")
    unset(DOXYGEN_EXECUTABLE CACHE)
    return()
else()
  # Generate a target for generating a specific documentation type
  function(add_doc_target doc_type)
    set(doc_output_dir ${CMAKE_CURRENT_BINARY_DIR})

    configure_file(${doc_type}.dox.in ${doc_type}.dox)
    
    set(doc_target doc-${doc_type})
    
    add_custom_target(
        ${doc_target}
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
  add_doc_target(devkit-lite)
endif()
