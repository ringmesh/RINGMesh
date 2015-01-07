# - Find GRGMesh
# Find the native GRGMesh includes and library
#
#  GRGMesh_FOUND        - True if GRGMesh found.
#  GRGMesh_INCLUDE_DIR  - where to find includes
#  GRGMesh_LIBRARY_DBG  - link these to use GRGMesh with DEBUG config
#  GRGMesh_LIBRARY_REL  - link these to use GRGMesh with RELEASE config
#  GRGMesh_LIBRARY      - set of libraries, normally GRGMesh_LIBRARY_DBG:GRGMesh_LIBRARY_REL


include(FindPkgMacros)
findpkg_begin(GRGMesh)

# Get path, convert backslashes as ${ENV_${var}}
getenv_path(GRGMesh_HOME)

# construct search paths
set(GRGMesh_PREFIX_PATH 
    ${GRGMesh_HOME} 
    ${ENV_GRGMesh_HOME} 
)


find_path(GRGMesh_INCLUDE_DIR 
	NAMES "grgmesh/"
	PATH_SUFFIXES "include" 
	PATHS ${GRGMesh_PREFIX_PATH}
)
#"include/grgmesh" grgmesh_assert.h

find_library(GRGMesh_LIBRARY_REL 
	NAMES "GRGMesh"
	PATH_SUFFIXES "lib"
	PATHS ${GRGMesh_PREFIX_PATH}
)


find_library(GRGMesh_LIBRARY_DBG  
	NAMES "GRGMesh_d"
	PATH_SUFFIXES "lib"
	PATHS ${GRGMesh_PREFIX_PATH}
)

# setting options as advanced to hide them in GUI
mark_as_advanced(GRGMesh_INCLUDE_DIR GRGMesh_LIBRARY_DBG GRGMesh_LIBRARY_REL)
make_library_set(GRGMesh_LIBRARY)

# Add verbosity if needed
if (NOT DEFINED SILENT_FIND_PACKAGE)
	SET(SILENT_FIND_PACKAGE TRUE)
endif (NOT DEFINED SILENT_FIND_PACKAGE)
if (NOT SILENT_FIND_PACKAGE)
	message("GRGMesh_PREFIX_PATH is " ${GRGMesh_PREFIX_PATH})
	message("GRGMesh_INCLUDE_DIR is " ${GRGMesh_INCLUDE_DIR})
	message("GRGMesh_LIBRARY_REL is " ${GRGMesh_LIBRARY_REL})
	message("GRGMesh_LIBRARY_DBG is " ${GRGMesh_LIBRARY_DBG})
endif (NOT SILENT_FIND_PACKAGE)

findpkg_finish(GRGMesh)
