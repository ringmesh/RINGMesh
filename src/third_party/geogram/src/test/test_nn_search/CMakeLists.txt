aux_source_directories(SOURCES "" .)
vor_add_executable(test_nn_search ${SOURCES})
target_link_libraries(test_nn_search geogram)


