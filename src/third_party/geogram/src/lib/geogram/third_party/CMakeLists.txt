vor_reset_warning_level()

aux_source_directories(SOURCES "Source Files"           .)
aux_source_directories(SOURCES "Source Files\\HLBFGS"   HLBFGS)
aux_source_directories(SOURCES "Source Files\\LM6"      LM6)
aux_source_directories(SOURCES "Source Files\\rply"     rply)
aux_source_directories(SOURCES "Source Files\\ANN"      ANN)
aux_source_directories(SOURCES "Source Files\\shewchuk" shewchuk)
aux_source_directories(SOURCES "Source Files\\NL"       NL)

if(GEOGRAM_WITH_TETGEN)
    aux_source_directories(SOURCES "Source Files\\tetgen" tetgen)
endif()

add_library(geogram_third_party OBJECT ${SOURCES})

