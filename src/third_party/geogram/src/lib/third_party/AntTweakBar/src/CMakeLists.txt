
include(${CMAKE_SOURCE_DIR}/cmake/opengl.cmake)

set(AntTweakBar_SOURCE_DIR ${CMAKE_SOURCE_DIR}/src/lib/third_party/AntTweakBar)

include_directories("${AntTweakBar_SOURCE_DIR}/src"
                    "${AntTweakBar_SOURCE_DIR}/include")

set(common_SOURCES TwEventGLFW.c TwEventGLUT.c 
                   LoadOGL.cpp TwBar.cpp TwFonts.cpp 
                   TwOpenGL.cpp TwPrecomp.cpp LoadOGLCore.cpp 
                   TwColors.cpp  TwMgr.cpp TwOpenGLCore.cpp)

set(common_HEADERS MiniGLFW.h MiniGLUT.h 
                   AntPerfTimer.h LoadOGL.h TwColors.h TwGraph.h 
                   TwOpenGL.h TwPrecomp.h LoadOGLCore.h 
                   TwBar.h TwFonts.h TwMgr.h TwOpenGLCore.h resource.h)

add_definitions(-DTW_EXPORTS -DTW_NO_LIB_PRAGMA)
                 
if(WIN32)
  set(AntTweakBar_HEADERS ${common_HEADERS})                   
  set(AntTweakBar_SOURCES ${common_SOURCES} TwEventWin.c)
elseif(APPLE)
  set(AntTweakBar_HEADERS ${common_HEADERS})                       
  set(AntTweakBar_SOURCES ${common_SOURCES})
else()
  set(AntTweakBar_HEADERS ${common_HEADERS})                     
  set(AntTweakBar_SOURCES ${common_SOURCES} TwEventX11.c)
  add_definitions(-D_UNIX -D__PLACEMENT_NEW_INLINE)
endif()

add_library(AntTweakBar ${AntTweakBar_SOURCES} ${AntTweakBar_HEADERS})

if(VORPALINE_BUILD_DYNAMIC)
    target_link_libraries(AntTweakBar ${OPENGL_LIBRARIES})
endif()  
                  
