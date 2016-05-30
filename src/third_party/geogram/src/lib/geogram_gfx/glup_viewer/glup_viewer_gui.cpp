/*
 *    _____   _       _   _     ____   
 *   /  ___| | |     | | | |   /  _ \  
 *   | |     | |     | | | |   | |_\ \ 
 *   | |  _  | |     | | | |   |  __ / 
 *   | |_| | | |___  | |_| |   | |     
 *   \_____/ |_____| \_____/   |_|     
 *
 *    _     _   _   _____   _          __  _____   _____
 *   | |   / / | | | ____| | |        / / | ____| |  _  \
 *   | |  / /  | | | |__   | |  __   / /  | |__   | |_| |
 *   | | / /   | | |  __|  | | /  | / /   |  __|  |  _  /
 *   | |/ /    | | | |___  | |/   |/ /    | |___  | | \ \
 *   |___/     |_| |_____| |___/|___/     |_____| |_|  \_\
 *
 *  Version 1.0
 *  Bruno Levy, April 2016
 *  INRIA, Project ALICE
 *
 *  Used internally by GLUP Viewer for interfacing with ImGUI
 *
 */

#include <geogram_gfx/basic/GL.h>
#include <geogram_gfx/glup_viewer/glup_viewer_gui.h>
#include <geogram_gfx/glup_viewer/glup_viewer.h>
#include <geogram_gfx/third_party/ImGui/imgui.h>
#include <geogram_gfx/third_party/ImGui/imgui_impl_glfw_gl3.h>
#include <geogram_gfx/third_party/ImGui/imgui_impl_glfw.h>
#include <geogram_gfx/third_party/quicktext/glQuickText.h>
#include <geogram/basic/logger.h>

#ifdef GEO_OS_EMSCRIPTEN
#include <GLFW/glfw3.h>
#include <emscripten.h>
#else
#include <third_party/glfw/include/GLFW/glfw3.h>
#endif

#include <geogram_gfx/GLUP/GLUP.h>

#include <string.h>
#include <iostream>

#ifdef GEO_GL_LEGACY
static bool vanillaGL = false;
#endif

void glup_viewer_gui_init(GLFWwindow* w) {
#ifdef GEO_GL_LEGACY        
    vanillaGL = (strcmp(glupCurrentProfileName(), "VanillaGL") == 0);
    if(vanillaGL) {
        GEO::Logger::out("ImGUI") << "Viewer GUI init (Vanilla)"
                                  << std::endl;        
        ImGui_ImplGlfw_Init(w, false);        
    } else
#endif
    {
        GEO::Logger::out("ImGUI") << "Viewer GUI init (GL3)"
                                  << std::endl;        
        ImGui_ImplGlfwGL3_Init(w, false);        
    }
}

void glup_viewer_gui_cleanup() {
#ifdef GEO_GL_LEGACY        
    if(vanillaGL) {    
        ImGui_ImplGlfw_Shutdown();
    } else
#endif
    {
        ImGui_ImplGlfwGL3_Shutdown();        
    }
}

void glup_viewer_gui_begin_frame() {
    glfwPollEvents();
#ifdef GEO_GL_LEGACY            
    if(vanillaGL) {
        ImGui_ImplGlfw_NewFrame();        
    } else
#endif        
    {
        ImGui_ImplGlfwGL3_NewFrame();
    }
}

void glup_viewer_gui_end_frame() {
    GlupViewerDisplayFunc overlay_func = glup_viewer_get_overlay_func();
    if(overlay_func != nil) {
        overlay_func();
        ImGui::Render();
    }
}

int glup_viewer_gui_takes_input() {
    if(!glup_viewer_is_enabled(GLUP_VIEWER_TWEAKBARS)) {
        return 0;
    }
    return (
        ImGui::GetIO().WantCaptureMouse ||
        ImGui::GetIO().WantCaptureKeyboard
    ) ? 1 : 0;
}

void glup_viewer_gui_mouse_button_callback(GLFWwindow* window, int button, int action, int mods) {
#ifdef GEO_GL_LEGACY            
    if(vanillaGL) {
        ImGui_ImplGlfw_MouseButtonCallback(window, button, action, mods);
    } else
#endif        
    {
        ImGui_ImplGlfwGL3_MouseButtonCallback(window, button, action, mods);
    }
}

void glup_viewer_gui_scroll_callback(GLFWwindow* window, double xoffset, double yoffset) {
#ifdef GEO_GL_LEGACY            
    if(vanillaGL) {
        ImGui_ImplGlfw_ScrollCallback(window, xoffset, yoffset);
    } else
#endif        
    {
        ImGui_ImplGlfwGL3_ScrollCallback(window, xoffset, yoffset);
    }
}

void glup_viewer_gui_key_callback(GLFWwindow* window, int key, int scancode, int action, int mods) {
#ifdef GEO_GL_LEGACY            
    if(vanillaGL) {
        ImGui_ImplGlFw_KeyCallback(window, key, scancode, action, mods);
    } else
#endif        
    {
        ImGui_ImplGlfwGL3_KeyCallback(window, key, scancode, action, mods);
    }
}

void glup_viewer_gui_char_callback(GLFWwindow* window, unsigned int c) {
#ifdef GEO_GL_LEGACY            
    if(vanillaGL) {
        ImGui_ImplGlfw_CharCallback(window, c);
    } else
#endif        
    {
        ImGui_ImplGlfwGL3_CharCallback(window, c);
    }
}

void glup_viewer_gui_resize(int width, int height) {
    ImGui::GetIO().DisplaySize.x = float(width);
    ImGui::GetIO().DisplaySize.y = float(height);    
}
