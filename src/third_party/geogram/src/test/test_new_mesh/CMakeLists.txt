aux_source_directories(SOURCES "" .)
vor_add_executable(test_new_mesh ${SOURCES})
target_link_libraries(test_new_mesh geogram)


