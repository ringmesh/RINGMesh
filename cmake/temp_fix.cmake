
# Add geogram include directories to the current ones
include_directories(SYSTEM ${PROJECT_SOURCE_DIR}/third_party/geogram/src/lib)

# Add geogram project libs to the libs with which RINGMesh will link
set(EXTRA_LIBS ${EXTRA_LIBS} geogram)
if(RINGMESH_WITH_GRAPHICS)
    include(${PROJECT_SOURCE_DIR}/third_party/geogram/cmake/opengl.cmake)
    set(EXTRA_LIBS ${EXTRA_LIBS} geogram_gfx ${OPENGL_LIBRARIES})
endif(RINGMESH_WITH_GRAPHICS)

# Add geogram bin directories to the current ones.
# It would be preferable to set the imported library location [JP].
# CMAKE_CFG_INTDIR is needed for Xcode (in MacOS) because the executables
# need the complete path to geogram libraries (with the configuration:
# Release or Debug).
link_directories(${GLOBAL_BINARY_DIR}/third_party/geogram/${CMAKE_BUILD_TYPE}/lib/${CMAKE_CFG_INTDIR})


# Add zlib project libs to the libs with which RINGMesh will link
if(UNIX)
    set(EXTRA_LIBS ${EXTRA_LIBS} z)
else()
    set(EXTRA_LIBS ${EXTRA_LIBS} debug zlibd optimized zlib)
endif()


# Add minizip project libs to the libs with which RINGMesh will link
set(EXTRA_LIBS ${EXTRA_LIBS} debug minizipd optimized minizip)


# Add tinyxml2 include directories to the current ones
include_directories(SYSTEM ${PROJECT_SOURCE_DIR}/third_party)


# Add tinyxml2 bin directories to the current ones 
# It would be preferable to set the imported library location [JP]
link_directories(${GLOBAL_BINARY_DIR}/third_party/tinyxml2/${CMAKE_BUILD_TYPE})

# Add tinyxml2 project libs to the libs with which RINGMesh will link
set(EXTRA_LIBS ${EXTRA_LIBS} debug tinyxml2d optimized tinyxml2)



# Add zlib include directories to the current ones
include_directories(SYSTEM ${GLOBAL_BINARY_DIR}/third_party/zlib/install/include)


# Add zlib bin directories to the current ones
# It would be preferable to set the imported library location [JP]
link_directories(${GLOBAL_BINARY_DIR}/third_party/zlib/${CMAKE_BUILD_TYPE}/install/lib)


# Add minizip bin directories to the current ones
# It would be preferable to set the imported library location [JP]
link_directories(${GLOBAL_BINARY_DIR}/third_party/minizip/${CMAKE_BUILD_TYPE})