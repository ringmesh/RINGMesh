aux_source_directories(SRC "" .)
add_executable(test_get_dependent_elements ${SRC})
target_link_libraries(test_get_dependent_elements RINGMesh geogram)

