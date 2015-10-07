# CMakeLists for RINGMesh tests

link_directories(${PROJECT_BINARY_DIR}/lib)

# Look at the test with cmake
# DO A LOOP here on a list of test name at least !!

set(TEST_DIRECTORY ${PROJECT_SOURCE_DIR}/tests/src)

# If necessary, use the RELATIVE flag, otherwise each source file may be listed 
# with full pathname. RELATIVE may makes it easier to extract an executable name
# automatically.
# file( GLOB APP_SOURCES RELATIVE app/*.cxx )

# file(GLOB TESTS RELATIVE ${TEST_DIRECTORY} *.cpp)
# tiens on peut faire des dossiers de vcproj
# il nous en faut un pour les tests, ce sera plus joli
# Creates folder "libraries" and adds target project (math.vcproj)

# Create folder of projects
# set_property(TARGET math PROPERTY FOLDER "libraries")

# foreach(testsourcefile ${TESTS} )
    # I used a simple string replace, to cut off .cpp.
    # string( REPLACE ".cpp" "" testname ${testsourcefile} )
    # add_executable( ${testname} ${testsourcefile} )
    # Make sure YourLib is linked to each app
    # target_link_libraries( ${testname} YourLib )
# endforeach( testsourcefile ${APP_SOURCES} )

include_directories(${PROJECT_SOURCE_DIR}/include)
include_directories(${GEOGRAM_PATH}/src/lib)
add_executable(test_io_bm ${TEST_DIRECTORY}/test_io_bm/main.cpp)
target_link_libraries(test_io_bm RINGMesh geogram)

# add_executable(test_io_ml ${PROJECT_SOURCE_DIR}/main.cpp)
# target_link_libraries(test_io_ml RINGMesh geogram)

# add_executable(test_make_unique ${PROJECT_SOURCE_DIR}/main.cpp)
# target_link_libraries(test_make_unique RINGMesh geogram)

# add_executable(test_rotation ${PROJECT_SOURCE_DIR}/main.cpp)
# target_link_libraries(test_rotation RINGMesh geogram)

# add_executable(test_tetrahedralize ${PROJECT_SOURCE_DIR}/main.cpp)
# target_link_libraries(test_tetrahedralize RINGMesh geogram)

# add_executable(test_build_bm_surface ${PROJECT_SOURCE_DIR}/main.cpp)
# target_link_libraries(test_build_bm_surface RINGMesh geogram)

# add_executable(test_get_dependent_elements ${TEST_DIRECTORY}//main.cpp)
# target_link_libraries(test_get_dependent_elements RINGMesh geogram)

