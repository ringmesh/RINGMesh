# CMakeList for RINGMesh viewer subproject

include(${CMAKE_SOURCE_DIR}/src/third_party/geogram/cmake/opengl.cmake)

# Set the target as an executable 
add_executable(ringmeshview main.cpp)

target_link_libraries(ringmeshview RINGMesh geogram geogram_gfx ${OPENGL_LIBRARIES})

# Add the project to a folder of projects for the tests 
set_property(TARGET ringmeshview PROPERTY FOLDER "Utilities")

