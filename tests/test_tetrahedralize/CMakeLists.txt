aux_source_directories(SRC "" .)
add_executable(test_tetrahedralize ${SRC})
target_link_libraries(test_tetrahedralize RINGMesh geogram)

