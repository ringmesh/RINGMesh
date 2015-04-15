if(CMAKE_BUILD_TYPE STREQUAL "Debug")
function(add_test folder)
    add_subdirectory(tests/${folder})
endfunction()

IF(WIN32)
    link_directories(build/ringmesh/Win64-vs2012/lib/${CMAKE_BUILD_TYPE})
ELSE()
    link_directories(build/ringmesh/Linux64-gcc-${CMAKE_BUILD_TYPE}/lib)
ENDIF()
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/tests/bin)

add_test(test_io_ml)
add_test(test_io_bm)
add_test(test_intersection)

ENDIF()
