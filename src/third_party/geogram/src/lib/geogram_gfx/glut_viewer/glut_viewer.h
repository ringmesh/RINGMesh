/*
 *    _____   _       _   _   _____        _     _   _   _____   _          __  _____   _____
 *   /  ___| | |     | | | | |_   _|      | |   / / | | | ____| | |        / / | ____| |  _  \
 *   | |     | |     | | | |   | |        | |  / /  | | | |__   | |  __   / /  | |__   | |_| |
 *   | |  _  | |     | | | |   | |        | | / /   | | |  __|  | | /  | / /   |  __|  |  _  /
 *   | |_| | | |___  | |_| |   | |        | |/ /    | | | |___  | |/   |/ /    | |___  | | \ \
 *   \_____/ |_____| \_____/   |_|        |___/     |_| |_____| |___/|___/     |_____| |_|  \_\
 *
 *  Version 1.0
 *  Bruno Levy, August 2009
 *  INRIA, Project ALICE
 *
 */

#ifndef __GLUT_VIEWER__
#define __GLUT_VIEWER__

/* #define WITH_HDR */
#define WITH_ANTTWEAKBAR
/* #define WITH_PNG  */
/* #define WITH_GEEX */

#include <geogram_gfx/third_party/glew/glew.h>


#if defined(_MSC_VER) && defined(GEO_DYNAMIC_LIBS)
#ifdef geogram_gfx_EXPORTS
#define GLUT_VIEWER_API __declspec(dllexport) 
#else
#define GLUT_VIEWER_API __declspec(dllimport) 
#endif
#else
#define GLUT_VIEWER_API
#endif

#ifdef WITH_HDR
#include "glut_viewer_hdr.h"
#endif

#include <malloc.h>

#ifdef __cplusplus
#include <string>
#endif

#ifdef __cplusplus
extern "C" {
#endif

typedef void (* GlutViewerDisplayFunc)();
typedef int (* GlutViewerKeyboardFunc)(char key);
typedef void (* GlutViewerKeyFunc)();
typedef void (* GlutViewerInitFunc)();
typedef void (* GlutViewerDragDropFunc)(char*);

enum GlutViewerEvent {
    GLUT_VIEWER_DOWN,
    GLUT_VIEWER_MOVE,
    GLUT_VIEWER_UP
};

typedef GLboolean (* GlutViewerMouseFunc)(
    float x, float y, int button, enum GlutViewerEvent event
);

#define GLUT_VIEWER_IDLE_REDRAW 0
#define GLUT_VIEWER_DRAW_SCENE 1
#define GLUT_VIEWER_SHOW_HELP 2
#define GLUT_VIEWER_BACKGROUND 3
#define GLUT_VIEWER_HDR 4
#define GLUT_VIEWER_ROTATE_LIGHT 5
#define GLUT_VIEWER_3D 6
#define GLUT_VIEWER_TWEAKBARS 7
#define GLUT_VIEWER_HDR_VIGNETTE 8
#define GLUT_VIEWER_HDR_UNSHARP_MASKING 9
#define GLUT_VIEWER_HDR_POSITIVE_UNSHARP_MASKING 10
#define GLUT_VIEWER_STEREOSCOPIC_DISPLAY 11
#define GLUT_VIEWER_CLIP 12
#define GLUT_VIEWER_SHOW_CLIP 13
#define GLUT_VIEWER_EDIT_CLIP 14
#define GLUT_VIEWER_FIXED_CLIP 15
#define GLUT_VIEWER_FULL_SCREEN 16
   
extern GLUT_VIEWER_API void glut_viewer_enable(int cap);
extern GLUT_VIEWER_API void glut_viewer_disable(int cap);
extern GLUT_VIEWER_API GLboolean glut_viewer_is_enabled(int cap);
extern GLUT_VIEWER_API GLboolean* glut_viewer_is_enabled_ptr(int cap);
extern GLUT_VIEWER_API void glut_viewer_toggle(int cap);

#define GLUT_VIEWER_HDR_EXPOSURE 0
#define GLUT_VIEWER_HDR_BLUR_AMOUNT 1
#define GLUT_VIEWER_HDR_BLUR_WIDTH 2
#define GLUT_VIEWER_HDR_GAMMA 3
#define GLUT_VIEWER_HDR_UNSHARP_MASKING_GAMMA 4
#define GLUT_VIEWER_STEREOSCOPIC_EYE_DISTANCE 5
#define GLUT_VIEWER_ZOOM 6
extern GLUT_VIEWER_API void glut_viewer_set_float(int param, GLfloat value);
extern GLUT_VIEWER_API GLfloat glut_viewer_get_float(int param);
extern GLUT_VIEWER_API GLfloat* glut_viewer_float_ptr(int param);

extern GLUT_VIEWER_API void glut_viewer_main_loop(int argc, char** argv);
extern GLUT_VIEWER_API void glut_viewer_exit_main_loop();
extern GLUT_VIEWER_API void glut_viewer_set_window_title(char* title);
extern GLUT_VIEWER_API void glut_viewer_set_display_func(GlutViewerDisplayFunc f);
extern GLUT_VIEWER_API void glut_viewer_set_overlay_func(GlutViewerDisplayFunc f);
extern GLUT_VIEWER_API void glut_viewer_set_keyboard_func(GlutViewerKeyboardFunc f);
extern GLUT_VIEWER_API void glut_viewer_set_mouse_func(GlutViewerMouseFunc f);
extern GLUT_VIEWER_API void glut_viewer_set_init_func(GlutViewerInitFunc f);
extern GLUT_VIEWER_API void glut_viewer_set_drag_drop_func(GlutViewerDragDropFunc f);
extern GLUT_VIEWER_API void glut_viewer_add_toggle(
    char key, GLboolean* pointer, const char* description
);
extern GLUT_VIEWER_API void glut_viewer_add_key_func(
    char key, GlutViewerKeyFunc f, const char* description
);
extern GLUT_VIEWER_API void glut_viewer_unbind_key(char key);
extern GLUT_VIEWER_API void glut_viewer_set_region_of_interest(
    float xmin, float ymin, float zmin, float xmax, float ymax, float zmax
);
extern GLUT_VIEWER_API void glut_viewer_set_screen_size(int w, int h);
extern GLUT_VIEWER_API void glut_viewer_get_screen_size(int* w, int* h);

extern GLUT_VIEWER_API void glut_viewer_clear_text();
extern GLUT_VIEWER_API void glut_viewer_printf(char* format, ...);
extern GLUT_VIEWER_API void glut_viewer_set_skybox(int cube_texture);
extern GLUT_VIEWER_API int glut_viewer_get_skybox();
extern GLUT_VIEWER_API void glut_viewer_set_background(int texture);
extern GLUT_VIEWER_API int glut_viewer_get_background();
extern GLUT_VIEWER_API void glut_viewer_set_background_color(GLfloat r, GLfloat g, GLfloat b);
extern GLUT_VIEWER_API void glut_viewer_set_background_color2(GLfloat r, GLfloat g, GLfloat b);
extern GLUT_VIEWER_API GLfloat* glut_viewer_get_background_color();
extern GLUT_VIEWER_API GLfloat* glut_viewer_get_background_color2();

extern GLUT_VIEWER_API float* glut_viewer_get_scene_quaternion();
extern GLUT_VIEWER_API float* glut_viewer_get_light_matrix();
extern GLUT_VIEWER_API float* glut_viewer_get_light_quaternion();
extern GLUT_VIEWER_API float* glut_viewer_get_clip_quaternion();

extern GLUT_VIEWER_API void glTexImage2DXPM(char const** xpm_data);
extern GLUT_VIEWER_API void glTexImage2Dfile(const char* filename);

extern GLUT_VIEWER_API GLboolean glut_viewer_load_image(
    const char* filename, GLuint* width, GLuint* height, GLuint* bpp, GLvoid** pixels
);

extern GLUT_VIEWER_API int glut_viewer_fps();

extern GLUT_VIEWER_API void glut_viewer_redraw();

extern GLUT_VIEWER_API void glut_viewer_save_transform_for_picking();
extern GLUT_VIEWER_API void glut_viewer_get_picked_ray(GLdouble* p, GLdouble* v);
extern GLUT_VIEWER_API void glut_viewer_get_picked_point(
    GLdouble* p, GLboolean* hit_background
);

extern GLUT_VIEWER_API void glut_viewer_keyboard_callback(unsigned char c, int x, int y);

#ifdef __cplusplus
}
#endif

#ifdef __cplusplus
inline void glut_viewer_printf(const std::string& s) {
    glut_viewer_printf((char*) s.c_str());
}

#endif

#endif

