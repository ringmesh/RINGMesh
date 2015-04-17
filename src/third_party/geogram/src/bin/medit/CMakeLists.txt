include(${CMAKE_SOURCE_DIR}/cmake/opengl.cmake)

# Reset the warning level (too may warnings that won't be fixed)
vor_reset_warning_level()

if(WIN32)
    # Disable warning C4244: '=' : conversion from 'double' to 'float', possible loss of data
    # There are far too many conversion warnings in medit, fixing them would be a loss of time,
    # better disable the warning

    add_definitions(/wd4244)

    # Disable warning C4702: unreachable code
    # Some parts of the code have been purposely made unreachable
    # by inserting an early return()

    add_definitions(/wd4702)
endif()

aux_source_directories(SOURCES "" .)
add_executable(medit ${SOURCES})
target_link_libraries(medit geogram_gfx ${OPENGL_LIBRARIES})
if(UNIX)
    target_link_libraries(medit m)
endif()

install_runtime_targets(medit)
