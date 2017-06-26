cmake_minimum_required(VERSION 2.8)
project (zlib)

# Organize outputs in a Linux way
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)

include_directories(${PROJECT_SOURCE_DIR})
file(GLOB zlib_src "${PROJECT_SOURCE_DIR}/zlib/*.c")
file(GLOB zlib_include "${PROJECT_SOURCE_DIR}/zlib/*.h")

# Add the tinyxml2 target as a shared library
add_library(zlib STATIC
    ${zlib_src}
    ${zlib_include})
	
target_link_libraries( zlib )