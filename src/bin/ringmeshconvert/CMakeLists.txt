# CMakeList for RINGMesh file format converter subproject

# Set the target as an executable 
add_executable(ringmeshconvert main.cpp)
   
target_link_libraries(ringmeshconvert RINGMesh geogram)

# Add the project to a folder of projects for the tests 
set_property(TARGET ringmeshconvert PROPERTY FOLDER "Utilities")


