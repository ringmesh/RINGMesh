MESSAGE(STATUS "RINGMesh home ${RINGMesh_HOME}")
MESSAGE(STATUS "Build type ${CMAKE_BUILD_TYPE}")

find_path(RINGMesh_INCLUDE_DIR NAMES ringmesh PATHS ${RINGMesh_HOME}/include)
if(WIN32)
    find_path(RINGMesh_CONFIG_INCLUDE_DIR NAMES ringmesh PATHS ${RINGMesh_HOME}/build/ringmesh)
else(WIN32)
    find_path(RINGMesh_CONFIG_INCLUDE_DIR NAMES ringmesh PATHS ${RINGMesh_HOME}/build/ringmesh/${CMAKE_BUILD_TYPE})
endif(WIN32)
find_path(GEOGRAM_INCLUDE_DIR NAMES geogram PATHS ${RINGMesh_HOME}/third_party/geogram/src/lib)
find_path(THIRD_PARTY_INCLUDE_DIR NAMES zlib PATHS ${RINGMesh_HOME}/third_party)
if(WIN32)
    find_library(RINGMesh_LIBRARY NAMES RINGMesh PATHS ${RINGMesh_HOME}/build/ringmesh/lib)
else(WIN32)
    find_library(RINGMesh_LIBRARY NAMES RINGMesh PATHS ${RINGMesh_HOME}/build/ringmesh/${CMAKE_BUILD_TYPE}/lib)
endif(WIN32)
if(WIN32)
    find_library(GEOGRAM_LIBRARY NAMES geogram PATHS ${RINGMesh_HOME}/build/geogram/lib)
else(WIN32)
    find_library(GEOGRAM_LIBRARY NAMES geogram PATHS ${RINGMesh_HOME}/build/geogram/${CMAKE_BUILD_TYPE}/lib)
endif(WIN32)

include(FindPackageHandleStandardArgs)

# handle the QUIETLY and REQUIRED arguments and set LOGGING_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args(RINGMesh DEFAULT_MSG RINGMesh_INCLUDE_DIR RINGMesh_CONFIG_INCLUDE_DIR GEOGRAM_INCLUDE_DIR RINGMesh_LIBRARY GEOGRAM_LIBRARY THIRD_PARTY_INCLUDE_DIR)


set(RINGMesh_LIBRARIES ${RINGMesh_LIBRARY}  ${GEOGRAM_LIBRARY} )
set(RINGMesh_INCLUDE_DIRS ${RINGMesh_INCLUDE_DIR} ${RINGMesh_CONFIG_INCLUDE_DIR} ${GEOGRAM_INCLUDE_DIR} ${THIRD_PARTY_INCLUDE_DIR})
#set(RINGMesh_DEFINITIONS )

# Tell cmake GUIs to ignore the "local" variables.
mark_as_advanced(RINGMesh_INCLUDE_DIR RINGMesh_CONFIG_INCLUDE_DIR GEOGRAM_INCLUDE_DIR THIRD_PARTY_INCLUDE_DIR RINGMesh_LIBRARY GEOGRAM_LIBRARY)
