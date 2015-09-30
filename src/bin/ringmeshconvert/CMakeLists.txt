IF(WIN32)
    link_directories(${CMAKE_SOURCE_DIR}/build/ringmesh/${GEOGRAM_PLATFORM}/lib/${CMAKE_BUILD_TYPE})
ELSE(WIN32)
    link_directories(${CMAKE_SOURCE_DIR}/build/ringmesh/${GEOGRAM_PLATFORM}-${CMAKE_BUILD_TYPE}/lib)
ENDIF(WIN32)
aux_source_directories(SRC "" .)
add_executable(ringmeshconvert ${SRC})
target_link_libraries(ringmeshconvert RINGMesh geogram)

