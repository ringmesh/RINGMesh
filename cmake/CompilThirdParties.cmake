#------------------------------------------------------------------------------------------------
# GEOGRAM 
# Set the path to Geogram code
set(GEOGRAM_PATH ${PROJECT_SOURCE_DIR}/third_party/geogram)

# Geogram platform dependent settings
if(WIN32)
    set(GEOGRAM_PATH_BIN ${GLOBAL_BINARY_DIR}/geogram)
    set(geoplatform Win-vs-dynamic-generic)  
    # Extra lib 
    set (EXTRA_LIBS ${EXTRA_LIBS} psapi)
    # TODO check that it is really necessary [JP]
    add_compile_options(-DGEO_DYNAMIC_LIBS) 
else(WIN32)
    set(GEOGRAM_PATH_BIN ${GLOBAL_BINARY_DIR}/geogram/${CMAKE_BUILD_TYPE})
    set(geoplatform Linux64-gcc-dynamic)
endif(WIN32)

# Define Geogram as an external project that we know how to
# configure and compile
ExternalProject_Add(geogram_ext
  PREFIX ${GEOGRAM_PATH_BIN}
  
  #--Download step--------------
  DOWNLOAD_COMMAND ""
  
  #--Update/Patch step----------
  UPDATE_COMMAND ""

  #--Configure step-------------
  SOURCE_DIR ${GEOGRAM_PATH}
  CONFIGURE_COMMAND ${CMAKE_COMMAND} ${GEOGRAM_PATH}
        -G ${CMAKE_GENERATOR} 
        -DVORPALINE_PLATFORM:STRING=${geoplatform}
        -DGEOGRAM_WITH_TETGEN:BOOL=${RINGMESH_WITH_TETGEN} 
        -DGEOGRAM_WITH_GRAPHICS:BOOL=${RINGMESH_WITH_GRAPHICS}
        -DGEOGRAM_WITH_EXPLORAGRAM:BOOL=OFF
        -DCMAKE_BUILD_TYPE:STRING=${CMAKE_BUILD_TYPE}
  
  #--Build step-----------------
  BINARY_DIR ${GEOGRAM_PATH_BIN}
  #-- Command to build geogram
  BUILD_COMMAND ${CMAKE_COMMAND} --build ${GEOGRAM_PATH_BIN} ${COMPILATION_OPTION}

  #--Install step---------------
  INSTALL_COMMAND "" 
)

ExternalProject_Add_Step(geogram_ext forcebuild
    DEPENDERS build
    ALWAYS 1
  )

# Add geogram include directories to the current ones
include_directories(SYSTEM ${GEOGRAM_PATH}/src/lib)

# Add geogram project libs to the libs with which RINGMesh will link
set(EXTRA_LIBS ${EXTRA_LIBS} geogram)
if(RINGMESH_WITH_GRAPHICS)
    include(${GEOGRAM_PATH}/cmake/opengl.cmake)
    set(EXTRA_LIBS ${EXTRA_LIBS} geogram_gfx ${OPENGL_LIBRARIES})
endif(RINGMESH_WITH_GRAPHICS)
    
# Add geogram bin directories to the current ones 
# It would be preferable to set the imported library location [JP]
link_directories(${GEOGRAM_PATH_BIN}/lib)

#------------------------------------------------------------------------------------------------
# tinyxml2 
# Set the path to tinyxml2 code
set(TINYXML2_PATH ${PROJECT_SOURCE_DIR}/third_party/tinyxml2)

# tinyxml2 platform dependent settings
if(WIN32)
    # Select the CMAKE_CONFIGURATION_TYPES from Visual Studio
    set(COMPILATION_OPTION --config ${CMAKE_CFG_INTDIR})
	set(TINYXML2_PATH_BIN  ${TINYXML2_PATH_BIN}/tinyxml2)
elseif(UNIX)
	set(TINYXML2_PATH_BIN  ${TINYXML2_PATH_BIN}/tinyxml2/${config})
    set(CMAKE_POSITION_INDEPENDENT_CODE TRUE)
    add_compile_options(-Wall -Wextra -Wno-long-long -Wconversion)
    
    # Determine gcc version and activate additional warnings available in latest versions
    execute_process(COMMAND ${CMAKE_CXX_COMPILER} -dumpversion OUTPUT_VARIABLE G++_VERSION)

    if(G++_VERSION VERSION_LESS 4.8 )
        message(FATAL_ERROR "RINGMesh require G++ version >= 4.8")
    endif()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
    add_compile_options(-Wsign-conversion)
    add_compile_options(-Wdouble-promotion)
else(WIN32)
    # Nothing else really supported
    message(FATAL_ERROR "Your platform is not supported. Please modify the cmakelist file.")    
endif(WIN32)

# Define tinyxml2 as an external project that we know how to
# configure and compile
ExternalProject_Add(tinyxml2_ext
  PREFIX ${TINYXML2_PATH_BIN}

  #--Configure step-------------
  SOURCE_DIR ${TINYXML2_PATH}
  CONFIGURE_COMMAND ${CMAKE_COMMAND} ${TINYXML2_PATH}
        -G ${CMAKE_GENERATOR} 
  
  #--Build step-----------------
  BINARY_DIR ${TINYXML2_PATH_BIN}
  #-- Command to build geogram
  BUILD_COMMAND ${CMAKE_COMMAND} --build ${TINYXML2_PATH_BIN}

)

ExternalProject_Add_Step(tinyxml2_ext forcebuild
    DEPENDERS build
    ALWAYS 1
  )

# Add geogram include directories to the current ones
include_directories(SYSTEM ${TINYXML2_PATH})

# Add geogram project libs to the libs with which RINGMesh will link
set(EXTRA_LIBS ${EXTRA_LIBS} tinyxml2)
    
# Add geogram bin directories to the current ones 
# It would be preferable to set the imported library location [JP]
link_directories(${TINYXML2_PATH_BIN}/lib)

#------------------------------------------------------------------------------------------------
# zlib 
# Set the path to zlib code
#set(ZLIB_PATH ${PROJECT_SOURCE_DIR}/third_party/zlib)