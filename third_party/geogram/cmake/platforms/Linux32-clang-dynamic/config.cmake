set(VORPALINE_ARCH_64 false)
set(VORPALINE_BUILD_DYNAMIC true)
include(${GEOGRAM_SOURCE_DIR}/cmake/platforms/Linux-clang.cmake)
add_flags(CMAKE_CXX_FLAGS -m32)
add_flags(CMAKE_C_FLAGS -m32)

