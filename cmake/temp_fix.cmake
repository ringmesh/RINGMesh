
# Add geogram include directories to the current ones
include_directories(SYSTEM ${PROJECT_SOURCE_DIR}/third_party/geogram/src/lib)


# Add geogram bin directories to the current ones.
# It would be preferable to set the imported library location [JP].
# CMAKE_CFG_INTDIR is needed for Xcode (in MacOS) because the executables
# need the complete path to geogram libraries (with the configuration:
# Release or Debug).
link_directories(${GLOBAL_BINARY_DIR}/third_party/geogram/lib/${CMAKE_CFG_INTDIR})


# Add tinyxml2 include directories to the current ones
include_directories(SYSTEM ${PROJECT_SOURCE_DIR}/third_party)


# Add tinyxml2 bin directories to the current ones 
# It would be preferable to set the imported library location [JP]
link_directories(${GLOBAL_BINARY_DIR}/third_party/tinyxml2)


# Add zlib include directories to the current ones
include_directories(SYSTEM ${GLOBAL_BINARY_DIR}/third_party/zlib/install/include)


# Add zlib bin directories to the current ones
# It would be preferable to set the imported library location [JP]
link_directories(${GLOBAL_BINARY_DIR}/third_party/zlib//install/lib)


# Add minizip bin directories to the current ones
# It would be preferable to set the imported library location [JP]
link_directories(${GLOBAL_BINARY_DIR}/third_party/minizip)