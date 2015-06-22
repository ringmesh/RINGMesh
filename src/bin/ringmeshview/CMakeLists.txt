include(${CMAKE_SOURCE_DIR}/src/third_party/geogram/cmake/opengl.cmake)
IF(WIN32)
    link_directories(${CMAKE_SOURCE_DIR}/build/ringmesh/${GEOGRAM_PLATFORM}/lib/${CMAKE_BUILD_TYPE})
ELSE(WIN32)
    link_directories(${CMAKE_SOURCE_DIR}/build/ringmesh/${GEOGRAM_PLATFORM}-${CMAKE_BUILD_TYPE}/lib)
ENDIF(WIN32)
aux_source_directories(SRC "" .)
add_executable(ringmeshview ${SRC})
target_link_libraries(ringmeshview RINGMesh geogram geogram_gfx ${OPENGL_LIBRARIES})


