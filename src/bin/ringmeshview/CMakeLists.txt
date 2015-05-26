include(${CMAKE_SOURCE_DIR}/src/third_party/geogram/cmake/opengl.cmake)

aux_source_directories(SRC "" .)
add_executable(ringmeshview ${SRC})
target_link_libraries(ringmeshview RINGMesh geogram geogram_gfx ${OPENGL_LIBRARIES})


