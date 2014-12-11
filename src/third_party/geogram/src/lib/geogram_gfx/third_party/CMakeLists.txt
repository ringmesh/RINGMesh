vor_reset_warning_level()

aux_source_directories(SOURCES "Source Files"              .)
aux_source_directories(SOURCES "Source Files\\glew"        glew)
aux_source_directories(SOURCES "Source Files\\freeglut"    freeglut)
aux_source_directories(SOURCES "Source Files\\AntTweakBar" AntTweakBar)

add_library(geogram_gfx_third_party OBJECT ${SOURCES})

