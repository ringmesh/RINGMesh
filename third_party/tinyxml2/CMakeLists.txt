cmake_minimum_required(VERSION 2.8)
project (tinyxml2)

# Organize outputs in a Linux way
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)

if(UNIX)
	set(CMAKE_POSITION_INDEPENDENT_CODE TRUE)
endif(UNIX)


include_directories(${PROJECT_SOURCE_DIR})
set(tinyxml2_src "${PROJECT_SOURCE_DIR}/tinyxml2.cpp")
set(tinyxml2_include "${PROJECT_SOURCE_DIR}/tinyxml2.h")

# Add the tinyxml2 target as a shared library
add_library(tinyxml2 STATIC
    ${tinyxml2_src}
    ${tinyxml2_include})
	
target_link_libraries(tinyxml2)