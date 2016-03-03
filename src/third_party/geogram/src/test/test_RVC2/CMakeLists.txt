aux_source_directories(SOURCES "" .)
vor_add_executable(test_RVC2 ${SOURCES})

## Probably bad practice ##
add_library(LSN least_square_normal.cpp)
target_link_libraries(LSN geogram)
target_link_libraries(test_RVC2 geogram LSN)
###########################



