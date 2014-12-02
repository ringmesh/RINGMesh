aux_source_directories(SOURCES "" .)
vor_add_executable(compute_OTM ${SOURCES})
target_link_libraries(compute_OTM geogram)
