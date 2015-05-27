# - Find RINGMesh
# Find the native RINGMesh includes and library
#
#  RINGMesh_FOUND        - True if RINGMesh found.
#  RINGMesh_INCLUDE_DIR  - where to find includes
#  RINGMesh_LIBRARY_DBG  - link these to use RINGMesh with DEBUG config
#  RINGMesh_LIBRARY_REL  - link these to use RINGMesh with RELEASE config
#  RINGMesh_LIBRARY      - set of libraries, normally RINGMesh_LIBRARY_DBG:RINGMesh_LIBRARY_REL


include(FindPkgMacros)
findpkg_begin(RINGMesh)

# Get path, convert backslashes as ${ENV_${var}}
getenv_path(RINGMesh_HOME)

# construct search paths
set(RINGMesh_PREFIX_PATH 
    ${RINGMesh_HOME} 
    ${ENV_RINGMesh_HOME} 
)


find_path(RINGMesh_INCLUDE_DIR 
	NAMES "ringmesh/"
	PATH_SUFFIXES "include" 
	PATHS ${RINGMesh_PREFIX_PATH}
)
#"include/ringmesh" ringmesh_assert.h

find_library(RINGMesh_LIBRARY_REL 
	NAMES "RINGMesh"
	PATH_SUFFIXES "lib" "build/ringmesh/Linux64-gcc-Release/lib" "build/ringmesh/Win64-vs2010/lib/Release"
	PATHS ${RINGMesh_PREFIX_PATH}
)


find_library(RINGMesh_LIBRARY_DBG  
	NAMES "RINGMesh"
	PATH_SUFFIXES "lib" "build/ringmesh/Linux64-gcc-Debug/lib" "build/ringmesh/Win64-vs2010/lib/Debug"
	PATHS ${RINGMesh_PREFIX_PATH}
)

# setting options as advanced to hide them in GUI
mark_as_advanced(RINGMesh_INCLUDE_DIR RINGMesh_LIBRARY_DBG RINGMesh_LIBRARY_REL)
make_library_set(RINGMesh_LIBRARY)

# Add verbosity if needed
if (NOT DEFINED SILENT_FIND_PACKAGE)
	SET(SILENT_FIND_PACKAGE TRUE)
endif (NOT DEFINED SILENT_FIND_PACKAGE)
if (NOT SILENT_FIND_PACKAGE)
	message("RINGMesh_PREFIX_PATH is " ${RINGMesh_PREFIX_PATH})
	message("RINGMesh_INCLUDE_DIR is " ${RINGMesh_INCLUDE_DIR})
	message("RINGMesh_LIBRARY_REL is " ${RINGMesh_LIBRARY_REL})
	message("RINGMesh_LIBRARY_DBG is " ${RINGMesh_LIBRARY_DBG})
endif (NOT SILENT_FIND_PACKAGE)

findpkg_finish(RINGMesh)
