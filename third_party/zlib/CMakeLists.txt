cmake_minimum_required(VERSION 2.8)
project (zlib)

include_directories(${PROJECT_SOURCE_DIR})
file(GLOB zlib_src "${PROJECT_SOURCE_DIR}/zlib/*.c")
file(GLOB zlib_include "${PROJECT_SOURCE_DIR}/zlib/*.h")

# Add the tinyxml2 target as a shared library
add_library(zlib SHARED
    ${zlib_src}
    ${zlib_include})
	
target_link_libraries(zlib)