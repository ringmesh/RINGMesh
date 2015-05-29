aux_source_directories(SRC "" .)
add_executable(test_io_bm ${SRC})
target_link_libraries(test_io_bm RINGMesh geogram)

