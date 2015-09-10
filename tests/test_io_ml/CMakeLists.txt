aux_source_directories(SRC "" .)
add_executable(test_io_ml ${SRC})
target_link_libraries(test_io_ml RINGMesh geogram)

