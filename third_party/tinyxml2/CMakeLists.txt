cmake_minimum_required(VERSION 2.8)
project (tinyxml2)
#------------------------------------------------------------------------------------------------
# Platform dependent settings
if(WIN32)
    # Select the CMAKE_CONFIGURATION_TYPES from Visual Studio
    set(COMPILATION_OPTION --config ${CMAKE_CFG_INTDIR})
elseif(UNIX)
    set(CMAKE_POSITION_INDEPENDENT_CODE TRUE)
    add_compile_options(-Wall -Wextra -Wno-long-long -Wconversion)
    
    # Determine gcc version and activate additional warnings available in latest versions
    execute_process(COMMAND ${CMAKE_CXX_COMPILER} -dumpversion OUTPUT_VARIABLE G++_VERSION)

    if(G++_VERSION VERSION_LESS 4.8 )
        message(FATAL_ERROR "tinyxml2 require G++ version >= 4.8")
    endif()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
    add_compile_options(-Wsign-conversion)
    add_compile_options(-Wdouble-promotion)

else(WIN32)
    # Nothing else really supported
    message(FATAL_ERROR "Your platform is not supported. Please modify the cmakelist file.")    
endif(WIN32)

# Organize outputs in a Linux way
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)

include_directories(${CMAKE_SOURCE_DIR})

set(tinyxml2_src "tinyxml2.cpp")
set(tinyxml2_include "tinyxml2.h")

# Add the tinyxml2 target as a shared library
add_library(tinyxml2 SHARED
    ${tinyxml2_src}
    ${tinyxml2_include})
	
target_link_libraries(tinyxml2)