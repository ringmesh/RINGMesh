vor_reset_warning_level()

# Since we are in a CMake "object library" (i.e.,
# a set of object files compiled with same options),
# the automatically generated geogram_gfx_EXPORTS
# pre-processor definition of the parent library
# is not automatically inherited (this is why
# I set it manually)
if(BUILD_SHARED_LIBS)
   add_definitions(-Dgeogram_gfx_EXPORTS)
endif()

aux_source_directories(SOURCES "Source Files"              .)
aux_source_directories(SOURCES "Source Files\\glew"        glew)
aux_source_directories(SOURCES "Source Files\\freeglut"    freeglut)
aux_source_directories(SOURCES "Source Files\\AntTweakBar" AntTweakBar)

add_library(geogram_gfx_third_party OBJECT ${SOURCES})

