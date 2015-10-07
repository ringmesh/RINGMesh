# CMakeLists for RINGMesh tests
# This file is included by the root cmake list of the project 

set(TEST_DIRECTORY ${PROJECT_SOURCE_DIR}/tests/src)

# Create the directory for the outputs of the tests
set(test_output_directory ${PROJECT_BINARY_DIR}/tests/output)
file(MAKE_DIRECTORY ${test_output_directory}) 

# Configure a header file to pass some of the CMake settings
# to the source code
set(test_config_file ${PROJECT_BINARY_DIR}/ringmesh/ringmesh_tests_config.h)
configure_file(
  "${TEST_DIRECTORY}/ringmesh_tests_config.h.in"
   ${test_config_file}
)

# The libraries with which the test are supposed to link
set(RINGMESH_LIBS geogram ringmesh)

# Get the name of the test files 
file(GLOB TESTS RELATIVE ${TEST_DIRECTORY} "${TEST_DIRECTORY}/*.cpp")

# For all the test files found
foreach(testsourcefile ${TESTS})
    # Get the name of the file without extension
    get_filename_component(testname ${testsourcefile} NAME_WE)
    
    # Set the include directories 
    include_directories(${PROJECT_SOURCE_DIR}/include)
    include_directories(${GEOGRAM_PATH}/src/lib)
    include_directories(${PROJECT_BINARY_DIR})
        
    # Add the target for that file
    add_executable(${testname} ${TEST_DIRECTORY}/${testsourcefile})
    
    # Link RINGMesh Make sure YourLib is linked to each app
    target_link_libraries(${testname} geogram RINGMesh)
    
    add_dependencies(${testname} RINGMesh)
    
    # Set the same output directories than for RINGMesh
    # So that the tests can be run smoothly on Windows
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)
    set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
    set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
    
    # Add the project to a folder of projects for the tests 
    set_property(TARGET ${testname} PROPERTY FOLDER "Tests")
endforeach(testsourcefile ${TESTS})