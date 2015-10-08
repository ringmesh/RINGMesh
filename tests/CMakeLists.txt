# CMakeLists for tests 

IF(WIN32)
    link_directories(${CMAKE_SOURCE_DIR}/build/ringmesh/${GEOGRAM_PLATFORM}/lib/${CMAKE_BUILD_TYPE})
ELSE()
    link_directories(${CMAKE_SOURCE_DIR}/build/ringmesh/${GEOGRAM_PLATFORM}-${CMAKE_BUILD_TYPE}/lib)
ENDIF()

# One subdirectory per test. 
add_subdirectory(test_io_ml)
add_subdirectory(test_io_bm)
add_subdirectory(test_make_unique)
add_subdirectory(test_build_bm_surface)
add_subdirectory(test_rotation)
add_subdirectory(test_get_dependent_elements)
add_subdirectory(test_tetrahedralize)
