aux_source_directories(SRC "" .)
add_executable(test_rotation ${SRC})
target_link_libraries(test_rotation RINGMesh geogram)

