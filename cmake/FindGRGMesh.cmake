# - Find GRGMesh
# Find the native GRGMesh includes and library
#
#   GRGMesh_FOUND       - True if GRGMesh found.
#   GRGMesh_INCLUDE_DIR - where to find includes
#   GRGMesh_LIBRARIES   - List of libraries when using GRGMesh.
#

if(GRGMesh_INCLUDE_DIR)
    # Already in cache, be silent
    set(GRGMesh_FIND_QUIETLY TRUE)
endif(GRGMesh_INCLUDE_DIR)

find_path(GRGMesh_INCLUDE_DIR "grgmesh_assert.h"
          PATH_SUFFIXES "GRGMesh")

find_library(GRGMesh_LIB "GRGMesh")

set(GRGMesh_INCLUDE_DIRS "${GRGMesh_INCLUDE_DIR}" "${GRGMesh_INCLUDE_DIR}/Modules" "${GRGMesh_INCLUDE_DIR}/Noise")
set(GRGMesh_LIBRARIES GRGMesh) 

# handle the QUIETLY and REQUIRED arguments and set GRGMesh_FOUND to TRUE if
# all listed variables are TRUE
include("FindPackageHandleStandardArgs")
find_package_handle_standard_args("GRGMesh" DEFAULT_MSG GRGMesh_INCLUDE_DIR GRGMesh_LIBRARIES)

mark_as_advanced(GRGMesh_INCLUDE_DIR GRGMesh_LIB)

