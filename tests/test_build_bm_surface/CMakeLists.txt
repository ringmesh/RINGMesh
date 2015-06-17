aux_source_directories(SRC "" .)
add_executable(test_build_bm_surface ${SRC})
target_link_libraries(test_build_bm_surface RINGMesh geogram)

