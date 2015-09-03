aux_source_directories(SRC "" .)
add_executable(test_make_unique ${SRC})
target_link_libraries(test_make_unique RINGMesh geogram)

