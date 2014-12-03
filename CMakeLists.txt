cmake_minimum_required(VERSION 2.8)
project(GRGMesh)


# Macro GatherProjectFiles definition
IF(WIN32)
	# visual studio specific section
	# controls inclusion of header files in solution
	# not really pretty but working COPIED FROM GRAPHITE
	SET(VISUAL_STUDIO_INCLUDE_HEADERS_IN_SOLUTION TRUE)
	# list of extensions added to the header folder
	SET(VISUAL_STUDIO_INCLUDE_HEADERS_EXTENSIONS_LIST *.h *.hpp *.hxx *.inl)	
	IF(VISUAL_STUDIO_INCLUDE_HEADERS_IN_SOLUTION)
		MACRO ( GatherProjectFiles ProjectName ProjectDir ProjectSources )		 
			# get the file lists
			SET ( Dir ${CMAKE_HOME_DIRECTORY}/${ProjectDir} )			
			FOREACH(var ${VISUAL_STUDIO_INCLUDE_HEADERS_EXTENSIONS_LIST})
				LIST(APPEND VS_EXTLIST_INCLUSION ${Dir}/${var})
			ENDFOREACH(var)
			FILE ( GLOB_RECURSE DirSources ${VS_EXTLIST_INCLUSION})			
		   SET ( "${ProjectSources}" ${DirSources} )
		ENDMACRO ( GatherProjectFiles )		
	ELSE(VISUAL_STUDIO_INCLUDE_HEADERS_IN_SOLUTION)	
		MACRO ( GatherProjectFiles ProjectName ProjectDir ProjectSources )
		ENDMACRO ( GatherProjectFiles )		
	ENDIF(VISUAL_STUDIO_INCLUDE_HEADERS_IN_SOLUTION)
ENDIF(WIN32)

IF (EXISTS ${CMAKE_SOURCE_DIR}/CMakeOptions.txt)
   MESSAGE(INFO: "Using local options file: " .${CMAKE_SOURCE_DIR}/CMakeOptions.txt)
   INCLUDE(${CMAKE_SOURCE_DIR}/CMakeOptions.txt)
ENDIF(EXISTS ${CMAKE_SOURCE_DIR}/CMakeOptions.txt)

IF(WIN32)
  ADD_DEFINITIONS(-DGRGMESH_EXPORTS)
ELSE(WIN32)
	IF(CMAKE_BUILD_TYPE MATCHES "Debug")
        # commented-out for now, since it causes problems when linking with
        # C++ libraries that are not compiled with this flag (changes the
        # size/definition of STL objects)
       # SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D_GLIBCXX_DEBUG")
	ENDIF(CMAKE_BUILD_TYPE MATCHES "Debug")
ENDIF(WIN32)

AUX_SOURCE_DIRECTORY(src/grgmesh SOURCES)
AUX_SOURCE_DIRECTORY(src/third_party SOURCES)
AUX_SOURCE_DIRECTORY(src/third_party/ANN SOURCES)

INCLUDE_DIRECTORIES(include)
INCLUDE_DIRECTORIES(include/third_party)

# Geogram
INCLUDE_DIRECTORIES(src/third_party/geogram/src/lib)
IF(WIN32)
    link_directories(${CMAKE_SOURCE_DIR}/build/geogram/${GEOGRAM_PLATFORM}/lib/${CMAKE_BUILD_TYPE})
ELSE(WIN32)
    link_directories(${CMAKE_SOURCE_DIR}/build/geogram/${GEOGRAM_PLATFORM}-${CMAKE_BUILD_TYPE}/lib)
ENDIF(WIN32)

# MG-Tetra
IF(NOT DEFINED USE_MG_TETRA)
  SET(USE_MG_TETRA FALSE)
ENDIF(NOT DEFINED USE_MG_TETRA)
IF(USE_MG_TETRA)	
  add_definitions(-DUSE_MG_TETRA)
  IF(${MG_TETRA} MATCHES .+)
    MESSAGE(STATUS "Using MG-Tetra directory =  ${MG_TETRA}")
    SET(MG_TETRA_DIR  ${MG_TETRA})
  ELSEIF( $ENV{MG_TETRA} MATCHES .+)
    MESSAGE(STATUS "Using MG-Tetra directory =  $ENV{MG_TETRA}")
    SET(MG_TETRA_DIR  $ENV{MG_TETRA})
  ELSE(${MG_TETRA} MATCHES .+)
    MESSAGE(STATUS "No MG-Tetra directory found")
  ENDIF(${MG_TETRA} MATCHES .+)


  INCLUDE_DIRECTORIES(${MG_TETRA_DIR}/include)
  
  if( CMAKE_SIZEOF_VOID_P EQUAL 8 )
    SET(PLATFORM _64)
  ELSE(CMAKE_SIZEOF_VOID_P EQUAL 8 )
    SET(PLATFORM)
  ENDIF(CMAKE_SIZEOF_VOID_P EQUAL 8 )
  
  IF(WIN32)
    LINK_DIRECTORIES(${MG_TETRA_DIR}/lib/WIN7${PLATFORM}_VC10)
  ELSE(WIN32)
    LINK_DIRECTORIES(${MG_TETRA_DIR}/lib/Linux${PLATFORM})
  ENDIF(WIN32)
ENDIF(USE_MG_TETRA)


#To use OpenMP 
IF(WIN32) 
       SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /wd4275 /wd4996 /MP")
       SET(CMAKE_C_FLAGS "${CMAKE_CXX_FLAGS} /wd4275 /wd4996 /MP")
ELSE(WIN32)
       SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
       SET(CMAKE_C_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
ENDIF(WIN32)

IF(UNIX)
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")
    SET(CMAKE_C_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")
ENDIF(UNIX)

SET(LIBRARY_OUTPUT_PATH lib)
IF(WIN32)
    GatherProjectFiles(VorteX . HeaderFiles)
    IF(GRGMESH_STATIC)
        add_library(GRGMesh STATIC ${SOURCES} ${HeaderFiles})
    ELSE(GRGMESH_STATIC)
        add_library(GRGMesh SHARED ${SOURCES} ${HeaderFiles})
    ENDIF(GRGMESH_STATIC)
ELSE(WIN32)
    IF(GRGMESH_STATIC)
        add_library(GRGMesh STATIC ${SOURCES})
    ELSE(GRGMESH_STATIC)
        add_library(GRGMesh SHARED ${SOURCES})
    ENDIF(GRGMESH_STATIC)
ENDIF(WIN32)

IF(USE_MG_TETRA)
  TARGET_LINK_LIBRARIES(GRGMesh meshgems mg-tetra)
ENDIF(USE_MG_TETRA)

SET(LIB_GeoGram geogram)
if(UNIX)
    SET(LIB_GeoGram -Wl,-whole-archive -lgeogram -Wl,-no-whole-archive)
ENDIF(UNIX)
target_link_libraries(GRGMesh ${LIB_GeoGram})

#add_executable(test src/main.cpp)
#target_link_libraries(test ${CMAKE_SOURCE_DIR}/build/geogram/Linux64-gcc-Debug/lib/libgeogram.a)
