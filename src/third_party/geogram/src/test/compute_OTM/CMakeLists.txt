aux_source_directories(SOURCES "" .)
vor_add_executable(compute_OTM ${SOURCES})

if(GEOGRAM_WITH_VORPALINE)
  target_link_libraries(compute_OTM geogram vorpalib)
else()
  target_link_libraries(compute_OTM geogram)
endif()