include(${CMAKE_SOURCE_DIR}/cmake/opengl.cmake)

add_subdirectory(third_party)

aux_source_directories(SOURCES "Source Files"              .)
aux_source_directories(SOURCES "Source Files\\glut_viewer" glut_viewer)

add_library(geogram_gfx ${SOURCES} $<TARGET_OBJECTS:geogram_gfx_third_party>)

if(VORPALINE_BUILD_DYNAMIC)
    # Must provide all symbols at link time
    target_link_libraries(geogram_gfx ${OPENGL_LIBRARIES})
endif()

# Install the library
install_devkit_targets(geogram_gfx)

# There are no public header files to install so far
