/*
 *    _____   _       _   _   _____
 *   /  ___| | |     | | | | |_   _|
 *   | |     | |     | | | |   | |
 *   | |  _  | |     | | | |   | |
 *   | |_| | | |___  | |_| |   | |
 *   \_____/ |_____| \_____/   |_|
 *    _     _   _   _____   _          __  _____   _____
 *   | |   / / | | | ____| | |        / / | ____| |  _  \
 *   | |  / /  | | | |__   | |  __   / /  | |__   | |_| |
 *   | | / /   | | |  __|  | | /  | / /   |  __|  |  _  /
 *   | |/ /    | | | |___  | |/   |/ /    | |___  | | \ \
 *   |___/     |_| |_____| |___/|___/     |_____| |_|  \_\
 *
 *  Version 1.0
 *  Bruno Levy, August 2009
 *  INRIA, Project ALICE
 *
 */

#include "glut_viewer.h"

#ifdef WITH_HDR
#include "glut_viewer_hdr.h"
#endif

#include <stdlib.h>
#include <stdarg.h>

#include <geogram_gfx/third_party/freeglut/glut.h>
#include <geogram_gfx/third_party/freeglut/freeglut_ext.h>

#include <stdio.h>

#include <math.h>
#include <string.h>
#include <ctype.h>

#ifdef WIN32
#include <windows.h>
#else
#include <time.h>
#endif

#ifdef WITH_PNG
#include <png/png.h>
#endif

#include <setjmp.h>

#ifdef WITH_ANTTWEAKBAR
#include <geogram_gfx/third_party/AntTweakBar/AntTweakBar.h>
#endif

#define glut_viewer_assert(x)                                                   \
    if(!(x)) {                                                                   \
        fprintf(stderr, "%s: assertion fail %s -- %d\n", #x, __FILE__, __LINE__); \
        abort();                                                                \
    }

#define glut_viewer_assert_not_reached()                                        \
    {                                                                            \
        fprintf(stderr, "should not go there!! %s -- %d\n", __FILE__, __LINE__); \
        abort();                                                                \
    }

#define glut_viewer_argused(x) (void) x

static GlutViewerDisplayFunc display_func = NULL;
static GlutViewerDisplayFunc overlay_func = NULL;
static GlutViewerKeyboardFunc keyboard_func = NULL;
static GlutViewerMouseFunc mouse_func = NULL;
static GlutViewerInitFunc init_func = NULL;
static GlutViewerDragDropFunc drag_drop_func = NULL;

static GLboolean* toggle[256];
static char* key_description[256];
static GlutViewerKeyFunc key_func[256];

static GLboolean caps[17] = {
    GL_FALSE, /* IDLE_REDRAW  */
    GL_FALSE, /* DRAW_SCENE   */
    GL_FALSE, /* SHOW_HELP    */
    GL_FALSE, /* BACKGROUND   */
    GL_FALSE, /* HDR          */
    GL_FALSE, /* ROTATE_LIGHT */
    GL_TRUE,  /* 3D           */
    GL_TRUE,  /* TWEAKBARS    */
    GL_TRUE,  /* HDR_VIGNETTE */
    GL_FALSE, /* HDR_UNSHARP_MASKING          */
    GL_FALSE, /* HDR_POSITIVE_UNSHARP_MASKING */
    GL_FALSE, /* STEREOSCOPIC_DISPLAY         */
    GL_FALSE, /* CLIP         */
    GL_TRUE,  /* SHOW_CLIP    */
    GL_FALSE, /* EDIT_CLIP    */
    GL_FALSE, /* FIXED_CLIP   */
    GL_FALSE  /* FULL_SCREEN  */   
};
static int nb_caps = 17;

static GLfloat params[7] = {
    0.5f,  /* HDR_EXPOSURE              */
    -1.0f, /* HDR_BLUR_AMOUNT           */
    0.5f,  /* HDR_BLUR_WIDTH            */
    0.55f, /* HDR_GAMMA                 */
    0.8f,  /* HDR_UNSHARP_MASKING_GAMMA */
    0.09f, /* STEREOSCOPIC_EYE_DISTANCE */
    1.0f   /* ZOOM                      */
};
static int nb_params = 7;

static char* title = "g33>|< Viewer";

static GLuint skybox_tex = 0;
static GLuint background_tex = 0;

/* ========================== Trackball prototypes ========================= */
/* (from SGI, see copyright below) */

/*
 * trackball.h
 * A virtual trackball implementation
 * Written by Gavin Bell for Silicon Graphics, November 1988.
 */

/*
 * Pass the x and y coordinates of the last and current positions of
 * the mouse, scaled so they are from (-1.0 ... 1.0).
 *
 * The resulting rotation is returned as a quaternion rotation in the
 * first paramater.
 */
void trackball(float q[4], float p1x, float p1y, float p2x, float p2y);

/*
 * Given two quaternions, add them together to get a third quaternion.
 * Adding quaternions to get a compound rotation is analagous to adding
 * translations to get a compound translation.  When incrementally
 * adding rotations, the first argument here should be the new
 * rotation, the second and third the total rotation (which will be
 * over-written with the resulting new total rotation).
 */
void add_quats(float* q1, float* q2, float* dest);

/*
 * A useful function, builds a rotation matrix in Matrix based on
 * given quaternion.
 */
void build_rotmatrix(float m[4][4], float q[4]);

/*
 * This function computes a quaternion based on an axis (defined by
 * the given vector) and an angle about which to rotate.  The angle is
 * expressed in radians.  The result is put into the third argument.
 */
void axis_to_quat(float a[3], float phi, float q[4]);

/* ==================== Timing ========================================= */

double now() {
#ifdef WIN32
    return (double) (GetTickCount()) / 1000.0;
#else
    return (double) clock() / (double) (CLOCKS_PER_SEC);
#endif
}

int glut_viewer_fps() {
    static int init = 1;
    static double ref = 0;
    static int frame = 0;
    static double result = 0;
    if(init) {
        ref = now();
        init = 0;
    }
    frame++;
    if((frame % 20) == 0) {
        double new_t = now();
        result = frame / (new_t - ref);
        ref = new_t;
        frame = 0;
    }
    return (int) result;
}

/* ==================== GlutViewer implementation ========================= */

int glut_viewer_W = 800, glut_viewer_H = 800;
static int last_x, last_y;
static float cur_rot[4] = {0.0, 0.0, 0.0, 0.0};
static float cur_xlat[3] = {0.0, 0.0, 0.0};
static float cur_rot_light[4] = {0.0, 0.0, 0.0, 0.0};
static float cur_rot_clip[4] = {0.0, 0.0, 0.0, 0.0};
static float cur_xlat_clip[3] = {0.0, 0.0, 0.0};

static enum {
    NONE, ROTATE, PAN, ZOOM
} mode = NONE;
static float xmin = -1, ymin = -1, zmin = -1, xmax = 1, ymax = 1, zmax = 1;
static float roi_radius = 1.4f;
static int window_w = 0, window_h = 0;

static void idle(void);

float* glut_viewer_get_scene_quaternion() {
    return cur_rot;
}

float* glut_viewer_get_light_quaternion() {
    return cur_rot_light;
}

float* glut_viewer_get_clip_quaternion() {
    return cur_rot_clip;
}

void glut_viewer_enable(int cap) {
#ifndef WITH_HDR
    if(cap == GLUT_VIEWER_HDR) {
        fprintf(stderr, "glut_viewer was not compiled with support for HDR\n");
        return;
    }
#endif
    glut_viewer_assert(cap < nb_caps);
    caps[cap] = GL_TRUE;
    if(cap == GLUT_VIEWER_IDLE_REDRAW) {
        glutIdleFunc(idle);
    }
}

void glut_viewer_disable(int cap) {
    glut_viewer_assert(cap < nb_caps);
    caps[cap] = GL_FALSE;
    if(cap == GLUT_VIEWER_IDLE_REDRAW) {
        glutIdleFunc(NULL);
    }
}

GLboolean glut_viewer_is_enabled(int cap) {
    glut_viewer_assert(cap < nb_caps);
    return caps[cap];
}

GLboolean* glut_viewer_is_enabled_ptr(int cap) {
    glut_viewer_assert(cap < nb_caps);
    return &(caps[cap]);
}

void glut_viewer_toggle(int cap) {
    glut_viewer_assert(cap < nb_caps);
    caps[cap] = !caps[cap];
}

void glut_viewer_set_float(int param, GLfloat value) {
    glut_viewer_assert(param < nb_params);
    params[param] = value;
}

GLfloat glut_viewer_get_float(int param) {
    glut_viewer_assert(param < nb_params);
    return params[param];
}

GLfloat* glut_viewer_float_ptr(int param) {
    glut_viewer_assert(param < nb_params);
    return &(params[param]);
}

void glut_viewer_get_screen_size(int* w, int* h) {
    *w = glut_viewer_W;
    *h = glut_viewer_H;
}

void glut_viewer_set_screen_size(int w, int h) {
    glut_viewer_W = w;
    glut_viewer_H = h;
}

static void reshape(int w, int h) {
#ifdef WITH_HDR
    if(glut_viewer_is_enabled(GLUT_VIEWER_HDR)) {
        glut_viewer_hdr_resize(w, h);
        window_w = w;
        window_h = h;
    }
    else
#endif
    {
        glViewport(0, 0, w, h);
        glut_viewer_W = w;
        glut_viewer_H = h;
        window_w = w;
        window_h = h;
    }
#ifdef WITH_ANTTWEAKBAR
    TwWindowSize(glut_viewer_W, glut_viewer_H);
#endif
}

static GLboolean transform_saved = GL_FALSE;
static GLdouble modelview_save[16];
static GLdouble project_save[16];
static GLint viewport_save[4];

void glut_viewer_save_transform_for_picking() {
    transform_saved = GL_TRUE;
    glGetDoublev(GL_MODELVIEW_MATRIX, modelview_save);
    glGetDoublev(GL_PROJECTION_MATRIX, project_save);
    glGetIntegerv(GL_VIEWPORT, viewport_save);
}

static GLdouble ray_p_save[3];
static GLdouble ray_v_save[3];
static GLdouble ray_p3d_save[3];
static GLboolean ray_hit_background;

void glut_viewer_get_picked_ray(GLdouble* p, GLdouble* v) {
    unsigned int i;
    for(i = 0; i < 3; i++) {
        p[i] = ray_p_save[i];
        v[i] = ray_v_save[i];
    }
}

void glut_viewer_get_picked_point(GLdouble* p, GLboolean* hit_background) {
    unsigned int i;
    for(i = 0; i < 3; i++) {
        p[i] = ray_p3d_save[i];
    }
    *hit_background = ray_hit_background;
}

static void save_picked_ray(int x_in, int y_in) {
    unsigned int i;
    GLdouble x, y;
    GLfloat z;
    GLdouble temp[3];
    if(transform_saved) {
        x = (double) x_in;
        y = (double) (viewport_save[3] - y_in);

        /*
         * Get the correct value of Z
         * (since even 2D mode uses perspective xform, if we do not do
         * that, we cannot get correct x,y).
         */
        gluProject(
            x, y, 0.0, modelview_save, project_save, viewport_save,
            &(temp[0]), &(temp[1]), &(temp[2])
        );

        gluUnProject(
            x, y, temp[2], modelview_save, project_save, viewport_save,
            &(ray_p_save[0]), &(ray_p_save[1]), &(ray_p_save[2])
        );
        gluUnProject(
            x, y, 1.0, modelview_save, project_save, viewport_save,
            &(ray_v_save[0]), &(ray_v_save[1]), &(ray_v_save[2])
        );
        for(i = 0; i < 3; i++) {
            ray_v_save[i] -= ray_p_save[i];
        }
        glReadPixels((GLint) x, (GLint) y, 1, 1, GL_DEPTH_COMPONENT, GL_FLOAT, &z);
        gluUnProject(
            x, y, z, modelview_save, project_save, viewport_save,
            &(ray_p3d_save[0]), &(ray_p3d_save[1]), &(ray_p3d_save[2])
        );
        ray_hit_background = (z == 1.0f);
    }
}

static GLboolean call_mouse_func(int x, int y, int button, enum GlutViewerEvent event) {
    double l = (double) (window_w > window_h ? window_w : window_h) / 2.0;
    double fx = x - (double) window_w / 2.0;
    double fy = y - (double) window_h / 2.0;
    int mx = (int) (3000.0 * fx / l);
    int my = -(int) (3000.0 * fy / l);
    GLboolean result = GL_FALSE;
    save_picked_ray(x, y);
    if(mouse_func == NULL) {
        return result;
    }
    result = mouse_func((float) mx, (float) my, button, event);
    glutPostRedisplay();
    return result;
}

static void mouse(int button, int state, int x, int y) {
#ifdef WITH_ANTTWEAKBAR
    if(
        glut_viewer_is_enabled(GLUT_VIEWER_TWEAKBARS) &&
        TwEventMouseButtonGLUT(button, state, x, y)
    ) {
        glutPostRedisplay();
        return;
    }
#endif
    if(call_mouse_func(x, y, button, state ? GLUT_VIEWER_UP : GLUT_VIEWER_DOWN)) {
        return;
    }
    mode = NONE;
    if(state != GLUT_DOWN) {
        return;
    }
    switch(button) {
        case GLUT_LEFT_BUTTON:
            mode = PAN;
            break;
        case GLUT_MIDDLE_BUTTON:
            mode = ZOOM;
            break;
        case GLUT_RIGHT_BUTTON:
            mode = ROTATE;
            break;
        case 3:
            *glut_viewer_float_ptr(GLUT_VIEWER_ZOOM) /= 1.1f;
            glutPostRedisplay();
            break;
        case 4:
            *glut_viewer_float_ptr(GLUT_VIEWER_ZOOM) *= 1.1f;
            glutPostRedisplay();
            break;
    }
    last_x = x;
    last_y = y;
}

static void passive_mouse(int x, int y) {
#ifdef WITH_ANTTWEAKBAR
    if(
        glut_viewer_is_enabled(GLUT_VIEWER_TWEAKBARS) &&
        TwEventMouseMotionGLUT(x, y)
    ) {
        glutPostRedisplay();
        return;
    }
#endif
}

static void motion(int x, int y) {

    float delta_rot[4];
    int W = window_w;
    int H = window_h;

#ifdef WITH_ANTTWEAKBAR
    if(
        glut_viewer_is_enabled(GLUT_VIEWER_TWEAKBARS) &&
        TwEventMouseMotionGLUT(x, y)
    ) {
        glutPostRedisplay();
        return;
    }
#endif

    if(call_mouse_func(x, y, -1, GLUT_VIEWER_MOVE)) {
        return;
    }

    switch(mode) {
        case ROTATE:
        {
            trackball(delta_rot,
                (float) (2 * last_x - W) / (float) W,
                (float) (H - 2 * last_y) / (float) H,
                (float) (2 * x - W) / (float) W,
                (float) (H - 2 * y) / (float) H
            );
            if(glut_viewer_is_enabled(GLUT_VIEWER_ROTATE_LIGHT)) {
                add_quats(delta_rot, cur_rot_light, cur_rot_light);
            } else {
                if(
                    glut_viewer_is_enabled(GLUT_VIEWER_EDIT_CLIP) ||
                    !glut_viewer_is_enabled(GLUT_VIEWER_FIXED_CLIP)
                ) {
                    add_quats(delta_rot, cur_rot_clip, cur_rot_clip);
                }
                if(!glut_viewer_is_enabled(GLUT_VIEWER_EDIT_CLIP)) {
                    add_quats(delta_rot, cur_rot, cur_rot);
                }
            }
        } break;
        case PAN:
        {
            if(!glut_viewer_is_enabled(GLUT_VIEWER_ROTATE_LIGHT)) {
                float delta_x = (float) (last_x - x) / (float) W;
                float delta_y = (float) (y - last_y) / (float) H;
                if(!glut_viewer_is_enabled(GLUT_VIEWER_EDIT_CLIP)) {
                    cur_xlat[0] -= 2.0f * delta_x / params[GLUT_VIEWER_ZOOM];
                    cur_xlat[1] -= 2.0f * delta_y / params[GLUT_VIEWER_ZOOM];
                }
                if(
                    glut_viewer_is_enabled(GLUT_VIEWER_EDIT_CLIP) ||
                    !glut_viewer_is_enabled(GLUT_VIEWER_FIXED_CLIP)
                ) {
                    cur_xlat_clip[0] -= 2.0f * delta_x / params[GLUT_VIEWER_ZOOM];
                    cur_xlat_clip[1] -= 2.0f * delta_y / params[GLUT_VIEWER_ZOOM];
                }
            }
        } break;
        case ZOOM:
        {
            if(!glut_viewer_is_enabled(GLUT_VIEWER_ROTATE_LIGHT)) {
                params[GLUT_VIEWER_ZOOM] *= (1.0f + (float) (y - last_y) / (float) H);
            }
        } break;
        default:
            break;
    }

    last_x = x;
    last_y = y;
    glutPostRedisplay();
}

static void drag_drop(const char* p) {
    if(drag_drop_func != NULL) {
        drag_drop_func((char*) (p));
        glutPostRedisplay();
    }
}

static float cur_text_x = 0;
static float cur_text_y = 0;

void glut_viewer_clear_text() {
    if(glut_viewer_H < 600 || glut_viewer_W < 600) {
        glLineWidth(1);
    } else {
        glLineWidth(2);
    }
    cur_text_x = -2800;
    cur_text_y = 2800;
    if(glut_viewer_W > glut_viewer_H) {
        cur_text_y *= ((float) glut_viewer_H / (float) glut_viewer_W);
    } else {
        cur_text_x *= ((float) glut_viewer_W / (float) glut_viewer_H);
    }
}

void glut_viewer_printf(char* format, ...) {
    va_list args;
    char buffer[1024], * p;
    va_start(args, format);
    vsprintf(buffer, format, args);
    va_end(args);
    glPushMatrix();
    glTranslatef(cur_text_x, cur_text_y, 0);
    for(p = buffer; *p; p++) {
        glutStrokeCharacter(GLUT_STROKE_MONO_ROMAN, *p);
    }
    glPopMatrix();
    cur_text_y -= 150;
}

static void draw_foreground() {
    int i;
    float aspect;

    glPointSize(1.0);
    glDisable(GL_POINT_SMOOTH);
    glut_viewer_clear_text();
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    aspect = (float) (glut_viewer_W) / (float) (glut_viewer_H);
    if(aspect < 1) {
        gluOrtho2D(-3000 * aspect, 3000 * aspect, -3000, 3000);
    } else {
        gluOrtho2D(-3000, 3000, -3000 / aspect, 3000 / aspect);
    }
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    glDisable(GL_DEPTH_TEST);
    glDisable(GL_LIGHTING);

    if(glut_viewer_is_enabled(GLUT_VIEWER_SHOW_HELP)) {
        glDisable(GL_CLIP_PLANE0);        
        glEnable(GL_BLEND);
        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
        glColor4f(0, 0, 0.5, 0.5);
        glBegin(GL_POLYGON);
        glVertex2f(-3000, -3000);
        glVertex2f(-3000, 3000);
        glVertex2f(3000, 3000);
        glVertex2f(3000, -3000);
        glEnd();
        glDisable(GL_BLEND);
        glColor4f(5, 5, 5, 1);

        glLineWidth(1);

        glColor3f(5, 5, 5);

        if(glut_viewer_is_enabled(GLUT_VIEWER_IDLE_REDRAW)) {
            glut_viewer_printf(" --- %s help [%4d FPS ] --- ", title, glut_viewer_fps());
        } else {
            glut_viewer_printf(" --- %s help --- ", title);
        }
        glut_viewer_printf("");

        for(i = 0; i < 256; i++) {
            if(!isprint(i)) {
                continue;
            }
            if(key_description[i] != NULL) {
                if(toggle[i] != NULL) {
                    glut_viewer_printf(
                        "%c:toggle %s (%s)\n",
                        (char) i, key_description[i],
                        *toggle[i] ? "on" : "off"
                    );
                }
                if(key_func[i] != NULL) {
                    glut_viewer_printf("%c:%s\n", (char) i, key_description[i]);
                }
            }
        }
        if(glut_viewer_is_enabled(GLUT_VIEWER_HDR)) {
            glut_viewer_printf(
                "Exposure: [<]%f[>]",
                (double) glut_viewer_get_float(GLUT_VIEWER_HDR_EXPOSURE)
            );
        }
    }

    if(
        glut_viewer_is_enabled(GLUT_VIEWER_TWEAKBARS) &&
        overlay_func != NULL
    ) {
        overlay_func();
    }
#ifdef WITH_ANTTWEAKBAR
    if(glut_viewer_is_enabled(GLUT_VIEWER_TWEAKBARS)) {
        TwDraw();
    }
#endif
}

static void display(void);

static void hdr_display(void) {
#ifdef WITH_HDR
    if(glut_viewer_is_enabled(GLUT_VIEWER_HDR)) {
        glut_viewer_hdr_begin_frame();
    }
#endif
    display();
#ifdef WITH_HDR
    if(glut_viewer_is_enabled(GLUT_VIEWER_HDR)) {
        glut_viewer_hdr_end_frame();
    }
#endif
    glutSwapBuffers();
}

static GLboolean in_main_loop_ = GL_FALSE;
void glut_viewer_redraw() {
    if(in_main_loop_) {
        hdr_display();
    }
}

#ifndef M_PI
#define M_PI 3.14159
#endif

static float to_radians(float x) {
    return x * 2.0f * (float) M_PI / 360.0f;
}

GLfloat* glut_viewer_get_light_matrix() {
    static GLfloat m[4][4];
    build_rotmatrix(m, cur_rot_light);
    return &m[0][0];
}

static float bkg1[3] = {1.0, 1.0, 1.0};
static float bkg2[3] = {0.0, 0.0, 0.5};

void glut_viewer_set_background_color(GLfloat r, GLfloat g, GLfloat b) {
    bkg1[0] = r;
    bkg1[1] = g;
    bkg1[2] = b;
}

void glut_viewer_set_background_color2(GLfloat r, GLfloat g, GLfloat b) {
    bkg2[0] = r;
    bkg2[1] = g;
    bkg2[2] = b;
}

GLfloat* glut_viewer_get_background_color() {
    return bkg1;
}

GLfloat* glut_viewer_get_background_color2() {
    return bkg2;
}

static void draw_background() {
    float z = 1.0f;
    float w = 1.0f;
    float h = 1.0f;

    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glDisable(GL_DEPTH_TEST);
    glDisable(GL_POLYGON_OFFSET_FILL);
    glClear((GLbitfield) (GL_DEPTH_BUFFER_BIT));
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    glShadeModel(GL_SMOOTH);
    glDisable(GL_LIGHTING);
    if(background_tex == 0) {
        glBegin(GL_POLYGON);
        glColor3fv(bkg1);
        glVertex3f(-1, -1, z);
        glVertex3f(1, -1, z);
        glColor3fv(bkg2);
        glVertex3f(1, 1, z);
        glVertex3f(-1, 1, z);
        glEnd();
    } else {
        glEnable(GL_TEXTURE_2D);
        glBindTexture(GL_TEXTURE_2D, background_tex);
        glColor4f(1, 1, 1, 1);
        glBegin(GL_POLYGON);
        glTexCoord2f(0, 0);
        glVertex3f(-w, -h, z);
        glTexCoord2f(1, 0);
        glVertex3f(w, -h, z);
        glTexCoord2f(1, 1);
        glVertex3f(w, h, z);
        glTexCoord2f(0, 1);
        glVertex3f(-w, h, z);
        glEnd();
        glDisable(GL_TEXTURE_2D);
    }
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_POLYGON_OFFSET_FILL);
    glEnable(GL_LIGHTING);
}

GLboolean gluInvertMatrix(float m[16], float invOut[16]) {
    float inv[16], det;
    int i;

    inv[0] = m[5] * m[10] * m[15] -
        m[5] * m[11] * m[14] -
        m[9] * m[6] * m[15] +
        m[9] * m[7] * m[14] +
        m[13] * m[6] * m[11] -
        m[13] * m[7] * m[10];

    inv[4] = -m[4] * m[10] * m[15] +
        m[4] * m[11] * m[14] +
        m[8] * m[6] * m[15] -
        m[8] * m[7] * m[14] -
        m[12] * m[6] * m[11] +
        m[12] * m[7] * m[10];

    inv[8] = m[4] * m[9] * m[15] -
        m[4] * m[11] * m[13] -
        m[8] * m[5] * m[15] +
        m[8] * m[7] * m[13] +
        m[12] * m[5] * m[11] -
        m[12] * m[7] * m[9];

    inv[12] = -m[4] * m[9] * m[14] +
        m[4] * m[10] * m[13] +
        m[8] * m[5] * m[14] -
        m[8] * m[6] * m[13] -
        m[12] * m[5] * m[10] +
        m[12] * m[6] * m[9];

    inv[1] = -m[1] * m[10] * m[15] +
        m[1] * m[11] * m[14] +
        m[9] * m[2] * m[15] -
        m[9] * m[3] * m[14] -
        m[13] * m[2] * m[11] +
        m[13] * m[3] * m[10];

    inv[5] = m[0] * m[10] * m[15] -
        m[0] * m[11] * m[14] -
        m[8] * m[2] * m[15] +
        m[8] * m[3] * m[14] +
        m[12] * m[2] * m[11] -
        m[12] * m[3] * m[10];

    inv[9] = -m[0] * m[9] * m[15] +
        m[0] * m[11] * m[13] +
        m[8] * m[1] * m[15] -
        m[8] * m[3] * m[13] -
        m[12] * m[1] * m[11] +
        m[12] * m[3] * m[9];

    inv[13] = m[0] * m[9] * m[14] -
        m[0] * m[10] * m[13] -
        m[8] * m[1] * m[14] +
        m[8] * m[2] * m[13] +
        m[12] * m[1] * m[10] -
        m[12] * m[2] * m[9];

    inv[2] = m[1] * m[6] * m[15] -
        m[1] * m[7] * m[14] -
        m[5] * m[2] * m[15] +
        m[5] * m[3] * m[14] +
        m[13] * m[2] * m[7] -
        m[13] * m[3] * m[6];

    inv[6] = -m[0] * m[6] * m[15] +
        m[0] * m[7] * m[14] +
        m[4] * m[2] * m[15] -
        m[4] * m[3] * m[14] -
        m[12] * m[2] * m[7] +
        m[12] * m[3] * m[6];

    inv[10] = m[0] * m[5] * m[15] -
        m[0] * m[7] * m[13] -
        m[4] * m[1] * m[15] +
        m[4] * m[3] * m[13] +
        m[12] * m[1] * m[7] -
        m[12] * m[3] * m[5];

    inv[14] = -m[0] * m[5] * m[14] +
        m[0] * m[6] * m[13] +
        m[4] * m[1] * m[14] -
        m[4] * m[2] * m[13] -
        m[12] * m[1] * m[6] +
        m[12] * m[2] * m[5];

    inv[3] = -m[1] * m[6] * m[11] +
        m[1] * m[7] * m[10] +
        m[5] * m[2] * m[11] -
        m[5] * m[3] * m[10] -
        m[9] * m[2] * m[7] +
        m[9] * m[3] * m[6];

    inv[7] = m[0] * m[6] * m[11] -
        m[0] * m[7] * m[10] -
        m[4] * m[2] * m[11] +
        m[4] * m[3] * m[10] +
        m[8] * m[2] * m[7] -
        m[8] * m[3] * m[6];

    inv[11] = -m[0] * m[5] * m[11] +
        m[0] * m[7] * m[9] +
        m[4] * m[1] * m[11] -
        m[4] * m[3] * m[9] -
        m[8] * m[1] * m[7] +
        m[8] * m[3] * m[5];

    inv[15] = m[0] * m[5] * m[10] -
        m[0] * m[6] * m[9] -
        m[4] * m[1] * m[10] +
        m[4] * m[2] * m[9] +
        m[8] * m[1] * m[6] -
        m[8] * m[2] * m[5];

    det = m[0] * inv[0] + m[1] * inv[4] + m[2] * inv[8] + m[3] * inv[12];

    if(det == 0) {
        return GL_FALSE;
    }

    det = 1.0f / det;

    for(i = 0; i < 16; i++) {
        invOut[i] = inv[i] * det;
    }

    return GL_TRUE;
}

static void actually_render_display(double offset) {
    GLfloat m[4][4];

    /* White diffuse light. */
    static GLfloat light_diffuse[] = {1.0, 1.0, 1.0, 1.0};

    /* Infinite light location. */
    static GLfloat light_position[] = {1.0, 1.0, 1.0, 0.0};

    /* field of view of the larger dimension in degrees */
    float camera_aperture = 9.0;

    float zoom = params[GLUT_VIEWER_ZOOM];
    float zNear = 1.0;                /* near clipping plane     */
    float zFar = 10.0;                /* far clipping plane      */
    float zScreen = 5.0;              /* screen projection plane */

    /* half of the distance between the eyes, if in stereo mode */
    float eye_offset = (float) offset;

    /* aspect ratio of the window */
    float aspect = (float) glut_viewer_W / (float) glut_viewer_H;

    /* half the width of the screen */
    float vue_max_size = (float) (
        zScreen * (float) tan(to_radians(camera_aperture) / 2.0f)
    );

    /* from the central point of view
       shift of the vue from the current point of view */
    float vue_shift = eye_offset * zNear / zScreen;

    float right;
    float top;
    float sq_w;

    double clip_eqn[4];
    float* background_color;

    if(glut_viewer_is_enabled(GLUT_VIEWER_BACKGROUND)) {
        draw_background();
    } else {
        glClearColor(bkg1[0], bkg1[1], bkg1[2], 1.0f);
        glClear((GLbitfield) (GL_DEPTH_BUFFER_BIT | GL_COLOR_BUFFER_BIT));
    }

    glEnable(GL_DEPTH_TEST);

    /* Setup the view. */
    glViewport(0, 0, glut_viewer_W, glut_viewer_H);

    /* Setup the projection matrix */
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();

    if(glut_viewer_is_enabled(GLUT_VIEWER_3D)) {
        if(aspect < 1) {
            top = vue_max_size;
            right = top * aspect;
        } else {
            right = vue_max_size;
            top = right / aspect;
        }
        right /= zoom;
        top /= zoom;
        glFrustum(-right - vue_shift, right - vue_shift, -top, top, zNear, zFar);
        glTranslatef(-eye_offset, 0.0, 0.0);
    } else {
        float x = 1.0f / zoom;
        float y = 1.0f / zoom;
        if(aspect > 1.0f) {
            x *= aspect;
        } else {
            y /= aspect;
        }
        glOrtho(-x, x, -y, y, zNear, zFar);
    }

    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();

    build_rotmatrix(m, cur_rot_light);
    glMultMatrixf(&m[0][0]);

    /* Enable a single OpenGL light. */
    glLightfv(GL_LIGHT0, GL_DIFFUSE, light_diffuse);
    glLightfv(GL_LIGHT0, GL_POSITION, light_position);

    glLoadIdentity();

    glDisable(GL_CLIP_PLANE0);
    if(glut_viewer_is_enabled(GLUT_VIEWER_CLIP)) {
        glPushMatrix();

        /* translate the world of the distance between eye and center */
        glTranslatef(0.0, 0.0, -zScreen);
        glTranslatef(cur_xlat_clip[0], cur_xlat_clip[1], cur_xlat_clip[2]);
        build_rotmatrix(m, cur_rot_clip);
        glMultMatrixf(&m[0][0]);

        background_color = glut_viewer_get_background_color();
        glColor3f(
            1.0f - background_color[0],
            1.0f - background_color[1],
            1.0f - background_color[2]
        );

        if(glut_viewer_is_enabled(GLUT_VIEWER_SHOW_CLIP)) {
            sq_w = 1.25f / params[GLUT_VIEWER_ZOOM];
            glLineWidth(4);
            glBegin(GL_LINE_LOOP);
            glVertex3f(sq_w, -sq_w, 0.0f);
            glVertex3f(sq_w, sq_w, 0.0f);
            glVertex3f(-sq_w, sq_w, 0.0f);
            glVertex3f(-sq_w, -sq_w, 0.0f);
            glEnd();
            glLineWidth(1);
            glBegin(GL_LINES);
            glVertex3f(-sq_w, 0.0f, 0.0f);
            glVertex3f(sq_w, 0.0f, 0.0f);
            glVertex3f(0.0f, -sq_w, 0.0f);
            glVertex3f(0.0f, sq_w, 0.0f);
            glEnd();
        }

        clip_eqn[0] = 0.0;
        clip_eqn[1] = 0.0;
        clip_eqn[2] = -1.0;
        clip_eqn[3] = 0.0;
        glEnable(GL_CLIP_PLANE0);
        glClipPlane(GL_CLIP_PLANE0, clip_eqn);
        glPopMatrix();
    } else {
        clip_eqn[0] = 0.0;
        clip_eqn[1] = 0.0;
        clip_eqn[2] = 0.0;
        clip_eqn[3] = 0.0;
        glClipPlane(GL_CLIP_PLANE0, clip_eqn);
    }

    /* translate the world of the distance between eye and center */
    glTranslatef(0.0, 0.0, -zScreen);

    /* apply pan */
    glTranslatef(cur_xlat[0], cur_xlat[1], cur_xlat[2]);

    /* if enabled, apply rotation */
    if(glut_viewer_is_enabled(GLUT_VIEWER_3D)) {
        build_rotmatrix(m, cur_rot);
        glMultMatrixf(&m[0][0]);
    }

    /* restrict view to the Region Of Interest */
    glScalef(1.0f / roi_radius, 1.0f / roi_radius, 1.0f / roi_radius);
    glTranslatef(
        -0.5f * (xmin + xmax),
        -0.5f * (ymin + ymax),
        -0.5f * (zmin + zmax)
    );
    if(skybox_tex != 0) {
#ifdef WITH_HDR
        GLfloat m[4][4];
        glMatrixMode(GL_MODELVIEW);
        glPushMatrix();
        glLoadIdentity();
        build_rotmatrix(m, cur_rot_light);
        glMultMatrixf(&m[0][0]);
        glut_viewer_draw_skybox(skybox_tex);
        glPopMatrix();
#else
        fprintf(stderr, "glut_viewer was compiled without HDR support, cannot draw skybox\n");
#endif
    }

    if(display_func != NULL) {
        glut_viewer_save_transform_for_picking();
        display_func();
    }

    draw_foreground();
}

static void display(void) {
    double eye_dist = glut_viewer_get_float(GLUT_VIEWER_STEREOSCOPIC_EYE_DISTANCE);

    if(init_func != NULL) {
        init_func();
        init_func = NULL;
    }

    if(glut_viewer_is_enabled(GLUT_VIEWER_STEREOSCOPIC_DISPLAY)) {
        /* left eye  */
        glDrawBuffer(GL_BACK_LEFT);
        glClear(GL_COLOR_BUFFER_BIT);
        actually_render_display(-eye_dist);

        /* right eye */
        glDrawBuffer(GL_BACK_RIGHT);
        glClear(GL_COLOR_BUFFER_BIT);
        actually_render_display(eye_dist);
    } else {
        glDrawBuffer(GL_BACK);
        glClear(GL_COLOR_BUFFER_BIT);
        actually_render_display(0.0);
    }
}

static void idle(void) {
    if(glut_viewer_is_enabled(GLUT_VIEWER_IDLE_REDRAW)) {
        hdr_display();
    }
}

static void copy_image_to_clipboard();

void glut_viewer_keyboard_callback(unsigned char c, int x, int y) {

#ifdef WITH_ANTTWEAKBAR
    if(glut_viewer_is_enabled(GLUT_VIEWER_TWEAKBARS) && TwEventKeyboardGLUT(c, x, y)) {
        glutPostRedisplay();
        return;
    }
#endif
    if(c == 3) {   /* 3 = <Ctrl> C */
        printf("copying image to clipboard\n");
        copy_image_to_clipboard();
        return;
    }
    if(keyboard_func != NULL) {
        keyboard_func((char) c);
    }
    if(toggle[c] != NULL) {
        *toggle[c] = !*toggle[c];
        glutPostRedisplay();
    }
    if(key_func[c] != NULL) {
        key_func[c]();
        glutPostRedisplay();
    }
    glutPostRedisplay();
}

static void toggle_clip() {
    glut_viewer_toggle(GLUT_VIEWER_CLIP);
    if(!glut_viewer_is_enabled(GLUT_VIEWER_CLIP)) {
        glut_viewer_disable(GLUT_VIEWER_EDIT_CLIP);
    }
}

static void toggle_edit_clip() {
    if(glut_viewer_is_enabled(GLUT_VIEWER_CLIP)) {
        glut_viewer_toggle(GLUT_VIEWER_EDIT_CLIP);
    } else {
        glut_viewer_disable(GLUT_VIEWER_EDIT_CLIP);
    }
}

static void toggle_fixed_clip() {
    if(glut_viewer_is_enabled(GLUT_VIEWER_CLIP)) {
        glut_viewer_toggle(GLUT_VIEWER_FIXED_CLIP);
        glut_viewer_disable(GLUT_VIEWER_EDIT_CLIP);
    }
}

void glut_viewer_special_callback(int c, int x, int y) {
    glut_viewer_argused(x);
    glut_viewer_argused(y);
    switch(c) {
        case GLUT_KEY_F1:
            toggle_clip();
            break;
        case GLUT_KEY_F2:
            toggle_edit_clip();
            break;
        case GLUT_KEY_F3:
            toggle_fixed_clip();
            break;
        case GLUT_KEY_F4:
            glut_viewer_toggle(GLUT_VIEWER_SHOW_CLIP);
            break;
    }
    glutPostRedisplay();
}

static void quit() {
    glut_viewer_exit_main_loop();
}

static void home() {
    cur_xlat[0] = cur_xlat[1] = cur_xlat[2] = 0.0;
    params[GLUT_VIEWER_ZOOM] = 1.0;
    trackball(cur_rot, 0.0, 0.0, 0.0, 0.0);
    trackball(cur_rot_light, 0.0, 0.0, 0.0, 0.0);
    trackball(cur_rot_clip, 0.0, 0.0, 0.0, 0.0);
}

#ifdef WITH_HDR
static void dec_exposure() {
    *glut_viewer_float_ptr(GLUT_VIEWER_HDR_EXPOSURE) /= 2.0;
}

static void inc_exposure() {
    *glut_viewer_float_ptr(GLUT_VIEWER_HDR_EXPOSURE) *= 2.0;
}

#endif

static void init() {

#ifdef WITH_HDR
    if(glut_viewer_is_enabled(GLUT_VIEWER_HDR)) {
        glut_viewer_hdr_init();
    }
#endif

    glut_viewer_disable(GLUT_VIEWER_IDLE_REDRAW);
    glut_viewer_enable(GLUT_VIEWER_DRAW_SCENE);
    if(!glut_viewer_is_enabled(GLUT_VIEWER_HDR)) {
        glut_viewer_enable(GLUT_VIEWER_BACKGROUND);
    }

    glut_viewer_add_key_func('q', quit, "quit");
    glut_viewer_add_key_func(27, quit, "quit");
    glut_viewer_add_key_func('H', home, "home");

#ifdef WITH_HDR
    if(glut_viewer_is_enabled(GLUT_VIEWER_HDR)) {
        glut_viewer_add_key_func('<', dec_exposure, NULL);
        glut_viewer_add_key_func('>', inc_exposure, NULL);
    }
#endif

    glut_viewer_add_toggle('h', &caps[GLUT_VIEWER_SHOW_HELP], "help");
    glut_viewer_add_toggle('l', &caps[GLUT_VIEWER_ROTATE_LIGHT], "light rotation");

    glEnable(GL_LIGHT0);
    glEnable(GL_LIGHTING);
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_NORMALIZE);

    glPolygonOffset(1.0, 2.0);
    glEnable(GL_POLYGON_OFFSET_FILL);

    home();
}

void glut_viewer_set_window_title(char* s) {
    title = s;
}

void glut_viewer_set_display_func(GlutViewerDisplayFunc f) {
    display_func = f;
}

void glut_viewer_set_overlay_func(GlutViewerDisplayFunc f) {
    overlay_func = f;
}

void glut_viewer_set_keyboard_func(GlutViewerKeyboardFunc f) {
    keyboard_func = f;
}

void glut_viewer_set_mouse_func(GlutViewerMouseFunc f) {
    mouse_func = f;
}

void glut_viewer_set_init_func(GlutViewerInitFunc f) {
    init_func = f;
}

void glut_viewer_set_drag_drop_func(GlutViewerDragDropFunc f) {
    drag_drop_func = f;
}

static void init_keys_if_needed() {
    static int first = 1;
    int i;
    if(first) {
        first = 0;
        for(i = 0; i < 256; i++) {
            toggle[i] = NULL;
            key_description[i] = NULL;
            key_func[i] = NULL;
        }
    }
}

static void set_key_description(char key, const char* description) {
    if(key_description[(int) key] != NULL) {
        printf(
            "Warning: overriding key %c with %s (was %s)\n",
            key, description ? description : "NULL", key_description[(int) key]
        );
    }
    free(key_description[(int) key]);
    if(description == NULL) {
        key_description[(int) key] = NULL;
    } else {
        key_description[(int) key] = malloc(strlen(description) + 1);
        strcpy(key_description[(int) key], description);
    }
}

void glut_viewer_add_toggle(char key, GLboolean* pointer, const char* description) {
    init_keys_if_needed();
    toggle[(int) key] = pointer;
    key_func[(int) key] = NULL;
    set_key_description(key, description);
}

void glut_viewer_add_key_func(char key, GlutViewerKeyFunc f, const char* description) {
    init_keys_if_needed();
    key_func[(int) key] = f;
    toggle[(int) key] = NULL;
    set_key_description(key, description);
}

void glut_viewer_unbind_key(char key) {
    key_func[(int) key] = NULL;
    toggle[(int) key] = NULL;
    set_key_description(key, NULL);
}

static void cleanup(void) {
    int i;
#ifdef WITH_ANTTWEAKBAR
    TwTerminate();
#endif
    for(i = 0; i < 256; ++i) {
        if(key_description[i] != NULL) {
            free(key_description[i]);
        }
    }
}

static jmp_buf env;

void glut_viewer_exit_main_loop() {
    /* longjmp(env, 1) ; */
    /* TODO: find a clean way of exiting the loop
       without exiting the program. longjmp is
       not compatible with C++...
     */
    exit(-1);
}

#ifdef WITH_ANTTWEAKBAR

static void TW_CALL SetQuatCallback(const void* value, void* clientData) {
    const float* dvalue = (const float*) (value);
    float* dclientdata = (float*) (clientData);
    dclientdata[0] = -dvalue[0];
    dclientdata[1] = -dvalue[1];
    dclientdata[2] = dvalue[2];
    dclientdata[3] = dvalue[3];
}

static void TW_CALL GetQuatCallback(void* value, void* clientData) {
    float* dvalue = (float*) (value);
    float* dclientdata = (float*) (clientData);
    dvalue[0] = -dclientdata[0];
    dvalue[1] = -dclientdata[1];
    dvalue[2] = dclientdata[2];
    dvalue[3] = dclientdata[3];
}

static void create_tw_gui() {
    TwBar* viewer_bar = TwNewBar("Viewer");
    TwDefine(" Viewer position='232 16' size='200 200' iconify");

    TwAddVarCB(viewer_bar, "Rotate", TW_TYPE_QUAT4F,
        SetQuatCallback, GetQuatCallback,
        glut_viewer_get_scene_quaternion(), "open"
    );

    TwAddVarCB(viewer_bar, "Light", TW_TYPE_QUAT4F,
        SetQuatCallback, GetQuatCallback,
        glut_viewer_get_light_quaternion(), "open"
    );

    if(glut_viewer_is_enabled(GLUT_VIEWER_HDR)) {
        TwAddVarRW(
            viewer_bar, "Exposure", TW_TYPE_FLOAT,
            glut_viewer_float_ptr(GLUT_VIEWER_HDR_EXPOSURE),
            "min=0.001 max=3.0 step=0.1"
        );
        TwAddVarRW(
            viewer_bar, "Gamma", TW_TYPE_FLOAT,
            glut_viewer_float_ptr(GLUT_VIEWER_HDR_GAMMA),
            "min=0.2 max=1.5 step=0.1"
        );
        TwAddVarRW(
            viewer_bar, "Vignette", TW_TYPE_BOOL8,
            glut_viewer_is_enabled_ptr(GLUT_VIEWER_HDR_VIGNETTE),
            ""
        );
        TwAddVarRW(
            viewer_bar, "Blur amount", TW_TYPE_FLOAT,
            glut_viewer_float_ptr(GLUT_VIEWER_HDR_BLUR_AMOUNT),
            "min=0.0 max=1.0 step=0.1"
        );
        TwAddVarRW(
            viewer_bar, "Blur width", TW_TYPE_FLOAT,
            glut_viewer_float_ptr(GLUT_VIEWER_HDR_BLUR_WIDTH),
            "min=0.0 max=20.0 step=1.0"
        );
        TwAddVarRW(
            viewer_bar, "Unshrp. Msk.", TW_TYPE_BOOL8,
            glut_viewer_is_enabled_ptr(GLUT_VIEWER_HDR_UNSHARP_MASKING),
            ""
        );
        TwAddVarRW(
            viewer_bar, "Unshrp. Msk. +", TW_TYPE_BOOL8,
            glut_viewer_is_enabled_ptr(GLUT_VIEWER_HDR_POSITIVE_UNSHARP_MASKING),
            ""
        );
        TwAddVarRW(
            viewer_bar, "Unshrp. Msk. Gamma", TW_TYPE_FLOAT,
            glut_viewer_float_ptr(GLUT_VIEWER_HDR_UNSHARP_MASKING_GAMMA),
            "min=0.2 max=1.5 step=0.1"
        );
    }

    if(glut_viewer_is_enabled(GLUT_VIEWER_STEREOSCOPIC_DISPLAY)) {
        TwAddVarRW(
            viewer_bar, "Stereo Vision", TW_TYPE_BOOL8,
            glut_viewer_is_enabled_ptr(GLUT_VIEWER_STEREOSCOPIC_DISPLAY),
            ""
        );
        TwAddVarRW(
            viewer_bar, "Stereo Eye dist.", TW_TYPE_FLOAT,
            glut_viewer_float_ptr(GLUT_VIEWER_STEREOSCOPIC_EYE_DISTANCE),
            "min=0.0 max=0.25 step=0.005"
        );
    }
}

#endif

void glut_viewer_main_loop(int argc, char** argv) {
    int i;
    glutInit(&argc, argv);
    in_main_loop_ = GL_TRUE;
    /*
       Note: we do not need GLUT_DEPTH in HDR mode (depth
       buffer is created with the offscreen buffer)
     */
    glutInitDisplayMode(
        GLUT_DOUBLE
        | GLUT_RGB
        | (unsigned int) (glut_viewer_is_enabled(GLUT_VIEWER_HDR) ? 0 : GLUT_DEPTH)
        | (unsigned int) (glut_viewer_is_enabled(GLUT_VIEWER_STEREOSCOPIC_DISPLAY) ? GLUT_STEREO : 0)
    );
    glutInitWindowSize(glut_viewer_W, glut_viewer_H);
    glutCreateWindow(title);
    for(i = 1; i < argc; i++) {
        if(!strcmp(argv[i], "-fs")) {
            glutFullScreen();
        }
    } 
    if(glut_viewer_is_enabled(GLUT_VIEWER_FULL_SCREEN)) {
       glutFullScreen();
    } 

#ifdef WITH_ANTTWEAKBAR
    TwInit(TW_OPENGL, NULL);
    TwWindowSize(glut_viewer_W, glut_viewer_H);
    glutSpecialFunc((GLUTspecialfun) TwEventSpecialGLUT);
    TwGLUTModifiersFunc(glutGetModifiers);
    create_tw_gui();
#endif

    glutReshapeFunc(reshape);
    glutDisplayFunc(hdr_display);
    glutMouseFunc(mouse);
    glutPassiveMotionFunc(passive_mouse);
    glutMotionFunc(motion);
    glutKeyboardFunc(glut_viewer_keyboard_callback);
    glutSpecialFunc(glut_viewer_special_callback);
#if FREEGLUT_HAS_DRAGDROP
    glutDragDropFunc(drag_drop);
#endif
    init();

    if(glut_viewer_is_enabled(GLUT_VIEWER_IDLE_REDRAW)) {
        glutIdleFunc(idle);
    }

    atexit(cleanup);
    if(!setjmp(env)) {
        glutMainLoop();
    }

    in_main_loop_ = GL_FALSE;
}

void glut_viewer_set_region_of_interest(
    float xm, float ym, float zm, float xM, float yM, float zM
) {
    xmin = xm;
    ymin = ym;
    zmin = zm;
    xmax = xM;
    ymax = yM;
    zmax = zM;
    roi_radius = (float) sqrt(
        0.25f * (xmax - xmin) * (xmax - xmin) +
        0.25f * (ymax - ymin) * (ymax - ymin) +
        0.25f * (zmax - zmin) * (zmax - zmin)
    );
}

void glut_viewer_set_skybox(int cube_texture) {
    skybox_tex = (GLuint) cube_texture;
}

int glut_viewer_get_skybox() {
    return (int) skybox_tex;
}

void glut_viewer_set_background(int tex) {
    background_tex = (GLuint) tex;
}

int glut_viewer_get_background() {
    return (int) background_tex;
}

/* =============================== Trackball code ========================= */

/*
 * (c) Copyright 1993, 1994, Silicon Graphics, Inc.
 * ALL RIGHTS RESERVED
 * Permission to use, copy, modify, and distribute this software for
 * any purpose and without fee is hereby granted, provided that the above
 * copyright notice appear in all copies and that both the copyright notice
 * and this permission notice appear in supporting documentation, and that
 * the name of Silicon Graphics, Inc. not be used in advertising
 * or publicity pertaining to distribution of the software without specific,
 * written prior permission.
 *
 * THE MATERIAL EMBODIED ON THIS SOFTWARE IS PROVIDED TO YOU "AS-IS"
 * AND WITHOUT WARRANTY OF ANY KIND, EXPRESS, IMPLIED OR OTHERWISE,
 * INCLUDING WITHOUT LIMITATION, ANY WARRANTY OF MERCHANTABILITY OR
 * FITNESS FOR A PARTICULAR PURPOSE.  IN NO EVENT SHALL SILICON
 * GRAPHICS, INC.  BE LIABLE TO YOU OR ANYONE ELSE FOR ANY DIRECT,
 * SPECIAL, INCIDENTAL, INDIRECT OR CONSEQUENTIAL DAMAGES OF ANY
 * KIND, OR ANY DAMAGES WHATSOEVER, INCLUDING WITHOUT LIMITATION,
 * LOSS OF PROFIT, LOSS OF USE, SAVINGS OR REVENUE, OR THE CLAIMS OF
 * THIRD PARTIES, WHETHER OR NOT SILICON GRAPHICS, INC.  HAS BEEN
 * ADVISED OF THE POSSIBILITY OF SUCH LOSS, HOWEVER CAUSED AND ON
 * ANY THEORY OF LIABILITY, ARISING OUT OF OR IN CONNECTION WITH THE
 * POSSESSION, USE OR PERFORMANCE OF THIS SOFTWARE.
 *
 * US Government Users Restricted Rights
 * Use, duplication, or disclosure by the Government is subject to
 * restrictions set forth in FAR 52.227.19(c)(2) or subparagraph
 * (c)(1)(ii) of the Rights in Technical Data and Computer Software
 * clause at DFARS 252.227-7013 and/or in similar or successor
 * clauses in the FAR or the DOD or NASA FAR Supplement.
 * Unpublished-- rights reserved under the copyright laws of the
 * United States.  Contractor/manufacturer is Silicon Graphics,
 * Inc., 2011 N.  Shoreline Blvd., Mountain View, CA 94039-7311.
 *
 * OpenGL(TM) is a trademark of Silicon Graphics, Inc.
 */
/*
 * Trackball code:
 *
 * Implementation of a virtual trackball.
 * Implemented by Gavin Bell, lots of ideas from Thant Tessman and
 *   the August '88 issue of Siggraph's "Computer Graphics," pp. 121-129.
 *
 * Vector manip code:
 *
 * Original code from:
 * David M. Ciemiewicz, Mark Grossman, Henry Moreton, and Paul Haeberli
 *
 * Much mucking with by:
 * Gavin Bell
 */

/*
 * This size should really be based on the distance from the center of
 * rotation to the point on the object underneath the mouse.  That
 * point would then track the mouse as closely as possible.  This is a
 * simple example, though, so that is left as an Exercise for the
 * Programmer.
 */
#define TRACKBALLSIZE (0.8f)

/*
 * Local function prototypes (not defined in trackball.h)
 */
static float tb_project_to_sphere(float, float, float);
static void normalize_quat(float[4]);

void
vzero(float* v)
{
    v[0] = 0.0;
    v[1] = 0.0;
    v[2] = 0.0;
}

void
vset(float* v, float x, float y, float z)
{
    v[0] = x;
    v[1] = y;
    v[2] = z;
}

void
vsub(const float* src1, const float* src2, float* dst)
{
    dst[0] = src1[0] - src2[0];
    dst[1] = src1[1] - src2[1];
    dst[2] = src1[2] - src2[2];
}

void
vcopy(const float* v1, float* v2)
{
    register int i;
    for(i = 0; i < 3; i++) {
        v2[i] = v1[i];
    }
}

void
vcross(const float* v1, const float* v2, float* cross)
{
    float temp[3];

    temp[0] = (v1[1] * v2[2]) - (v1[2] * v2[1]);
    temp[1] = (v1[2] * v2[0]) - (v1[0] * v2[2]);
    temp[2] = (v1[0] * v2[1]) - (v1[1] * v2[0]);
    vcopy(temp, cross);
}

float
vlength(const float* v)
{
    return (float) sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
}

void
vscale(float* v, float div)
{
    v[0] *= div;
    v[1] *= div;
    v[2] *= div;
}

void
vnormal(float* v)
{
    vscale(v, 1.0f / vlength(v));
}

float
vdot(const float* v1, const float* v2)
{
    return v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2];
}

void
vadd(const float* src1, const float* src2, float* dst)
{
    dst[0] = src1[0] + src2[0];
    dst[1] = src1[1] + src2[1];
    dst[2] = src1[2] + src2[2];
}

/*
 * Ok, simulate a track-ball.  Project the points onto the virtual
 * trackball, then figure out the axis of rotation, which is the cross
 * product of P1 P2 and O P1 (O is the center of the ball, 0,0,0)
 * Note:  This is a deformed trackball-- is a trackball in the center,
 * but is deformed into a hyperbolic sheet of rotation away from the
 * center.  This particular function was chosen after trying out
 * several variations.
 *
 * It is assumed that the arguments to this routine are in the range
 * (-1.0 ... 1.0)
 */
void
trackball(float q[4], float p1x, float p1y, float p2x, float p2y)
{
    float a[3]; /* Axis of rotation */
    float phi;  /* how much to rotate about axis */
    float p1[3], p2[3], d[3];
    float t;

    if(p1x == p2x && p1y == p2y) {
        /* Zero rotation */
        vzero(q);
        q[3] = 1.0;
        return;
    }

    /*
     * First, figure out z-coordinates for projection of P1 and P2 to
     * deformed sphere
     */
    vset(p1, p1x, p1y, tb_project_to_sphere(TRACKBALLSIZE, p1x, p1y));
    vset(p2, p2x, p2y, tb_project_to_sphere(TRACKBALLSIZE, p2x, p2y));

    /*
     *  Now, we want the cross product of P1 and P2
     */
    vcross(p2, p1, a);

    /*
     *  Figure out how much to rotate around that axis.
     */
    vsub(p1, p2, d);
    t = vlength(d) / (2.0f * TRACKBALLSIZE);

    /*
     * Avoid problems with out-of-control values...
     */
    if(t > 1.0f) {
        t = 1.0;
    }
    if(t < -1.0f) {
        t = -1.0;
    }
    phi = 2.0f * (float) asin(t);

    axis_to_quat(a, phi, q);
}

/*
 *  Given an axis and angle, compute quaternion.
 */
void
axis_to_quat(float a[3], float phi, float q[4])
{
    vnormal(a);
    vcopy(a, q);
    vscale(q, (float) sin(phi / 2.0f));
    q[3] = (float) cos(phi / 2.0f);
}

/*
 * Project an x,y pair onto a sphere of radius r OR a hyperbolic sheet
 * if we are away from the center of the sphere.
 */
static float
tb_project_to_sphere(float r, float x, float y)
{
    float d, t, z;

    d = (float) sqrt(x * x + y * y);
    if(d < r * 0.70710678118654752440f) {     /* Inside sphere */
        z = (float) sqrt(r * r - d * d);
    } else {           /* On hyperbola */
        t = r / 1.41421356237309504880f;
        z = t * t / d;
    }
    return z;
}

/*
 * Given two rotations, e1 and e2, expressed as quaternion rotations,
 * figure out the equivalent single rotation and stuff it into dest.
 *
 * NOTE: This routine is written so that q1 or q2 may be the same
 * as dest (or each other).
 */

void
add_quats(float q1[4], float q2[4], float dest[4])
{
    float t1[4], t2[4], t3[4];
    float tf[4];

    vcopy(q1, t1);
    vscale(t1, q2[3]);

    vcopy(q2, t2);
    vscale(t2, q1[3]);

    vcross(q2, q1, t3);
    vadd(t1, t2, tf);
    vadd(t3, tf, tf);
    tf[3] = q1[3] * q2[3] - vdot(q1, q2);

    dest[0] = tf[0];
    dest[1] = tf[1];
    dest[2] = tf[2];
    dest[3] = tf[3];

    normalize_quat(dest);
}

/*
 * Quaternions always obey:  a^2 + b^2 + c^2 + d^2 = 1.0
 * If they don't add up to 1.0, dividing by their magnitued will
 * renormalize them.
 *
 * Note: See the following for more information on quaternions:
 *
 * - Shoemake, K., Animating rotation with quaternion curves, Computer
 *   Graphics 19, No 3 (Proc. SIGGRAPH'85), 245-254, 1985.
 * - Pletinckx, D., Quaternion calculus as a basic tool in computer
 *   graphics, The Visual Computer 5, 2-13, 1989.
 */
static void
normalize_quat(float q[4])
{
    int i;
    float mag;

    mag = (q[0] * q[0] + q[1] * q[1] + q[2] * q[2] + q[3] * q[3]);
    for(i = 0; i < 4; i++) {
        q[i] /= mag;
    }
}

/*
 * Build a rotation matrix, given a quaternion rotation.
 *
 */
void
build_rotmatrix(float m[4][4], float q[4])
{
    m[0][0] = 1.0f - 2.0f * (q[1] * q[1] + q[2] * q[2]);
    m[0][1] = 2.0f * (q[0] * q[1] - q[2] * q[3]);
    m[0][2] = 2.0f * (q[2] * q[0] + q[1] * q[3]);
    m[0][3] = 0.0f;

    m[1][0] = 2.0f * (q[0] * q[1] + q[2] * q[3]);
    m[1][1] = 1.0f - 2.0f * (q[2] * q[2] + q[0] * q[0]);
    m[1][2] = 2.0f * (q[1] * q[2] - q[0] * q[3]);
    m[1][3] = 0.0f;

    m[2][0] = 2.0f * (q[2] * q[0] - q[1] * q[3]);
    m[2][1] = 2.0f * (q[1] * q[2] + q[0] * q[3]);
    m[2][2] = 1.0f - 2.0f * (q[1] * q[1] + q[0] * q[0]);
    m[2][3] = 0.0f;

    m[3][0] = 0.0f;
    m[3][1] = 0.0f;
    m[3][2] = 0.0f;
    m[3][3] = 1.0f;
}

/*---------------------------------------------------------------------------*/

#ifdef WIN32

static void copy_image_to_clipboard() {
    int w = glut_viewer_W;
    int h = glut_viewer_H;
    /*    int row_len = w * 3 ;  */
    int image_size, size;
    HANDLE handle;
    char* base_mem;
    char* pData;
    BITMAPINFOHEADER* header;

    /* Step 1: Try to open Window's clipboard
     *   NULL -> bind to current process
     */
    if(!OpenClipboard(NULL)) {
        return;
    }

    /* Step 2: create a shared memory segment, with
     * a DIB (Device Independent Bitmap) in it.
     */

    image_size = 3 * w * h;
    size = sizeof(BITMAPINFOHEADER) + image_size;

    handle = (HANDLE) GlobalAlloc(GHND, size);
    if(handle == NULL) {
        return;
    }

    pData = GlobalLock((HGLOBAL) handle);
    header = (BITMAPINFOHEADER*) pData;

    header->biSize = sizeof(BITMAPINFOHEADER);
    header->biWidth = w;
    header->biHeight = h;
    header->biPlanes = 1;
    header->biBitCount = 24;
    header->biCompression = BI_RGB;
    header->biSizeImage = 0;
    header->biXPelsPerMeter = 1000000;
    header->biYPelsPerMeter = 1000000;
    header->biClrUsed = 0;
    header->biClrImportant = 0;

    base_mem = pData + sizeof(BITMAPINFOHEADER);
    glReadPixels(0, 0, w, h, GL_BGR_EXT, GL_UNSIGNED_BYTE, base_mem);

    /* Step 3: Detach the clipboard from current process. */
    GlobalUnlock((HGLOBAL) handle);
    EmptyClipboard();
    SetClipboardData(CF_DIB, handle);
    CloseClipboard();
}

#else
static void copy_image_to_clipboard() {
}

#endif

/* ---------------------------------------------------------------------------------------------- */

int htoi(char digit) {
    if(digit >= '0' && digit <= '9') {
        return digit - '0';
    }
    if(digit >= 'a' && digit <= 'f') {
        return digit - 'a' + 10;
    }
    if(digit >= 'A' && digit <= 'F') {
        return digit - 'A' + 10;
    }
    fprintf(stderr, "xpm: unknown digit\n");
    return 0;
}

/* The colormap. */
static unsigned char i2r[1024];
static unsigned char i2g[1024];
static unsigned char i2b[1024];
static unsigned char i2a[1024];

/* Converts a two-digit XPM color code into
 *  a color index.
 */
static int char_to_index[256][256];

void glTexImage2DXPM(const char** xpm_data) {
    int width, height, nb_colors, chars_per_pixel;
    int line = 0;
    int color = 0;
    int key1 = 0, key2 = 0;
    char* colorcode;
    int r, g, b;
    int none;
    int x, y;
    unsigned char* rgba;
    unsigned char* pixel;
    sscanf(xpm_data[line], "%d%d%d%d", &width, &height, &nb_colors, &chars_per_pixel);
    line++;
    if(nb_colors > 1024) {
        fprintf(stderr, "xpm with more than 1024 colors\n");
        return;
    }
    if(chars_per_pixel != 1 && chars_per_pixel != 2) {
        fprintf(stderr, "xpm with more than 2 chars per pixel\n");
        return;
    }
    for(color = 0; color < nb_colors; color++) {
        key1 = xpm_data[line][0];
        key2 = (chars_per_pixel == 2) ? xpm_data[line][1] : 0;
        colorcode = strstr(xpm_data[line], "c #");
        none = 0;
        if(colorcode == NULL) {
            colorcode = "c #000000";
            if(strstr(xpm_data[line], "None") != NULL) {
                none = 1;
            } else {
                fprintf(stderr, "unknown xpm color entry (replaced with black)\n");
            }
        }
        colorcode += 3;

        r = 16 * htoi(colorcode[0]) + htoi(colorcode[1]);
        g = 16 * htoi(colorcode[2]) + htoi(colorcode[3]);
        b = 16 * htoi(colorcode[4]) + htoi(colorcode[5]);

        i2r[color] = (unsigned char) r;
        i2g[color] = (unsigned char) g;
        i2b[color] = (unsigned char) b;
        i2a[color] = none ? 0 : 255;
        char_to_index[key1][key2] = color;
        line++;
    }
    rgba = (unsigned char*) malloc((size_t) (width * height * 4));
    pixel = rgba;
    for(y = 0; y < height; y++) {
        for(x = 0; x < width; x++) {
            if(chars_per_pixel == 2) {
                key1 = xpm_data[line][2 * x];
                key2 = xpm_data[line][2 * x + 1];
            } else {
                key1 = xpm_data[line][x];
                key2 = 0;
            }
            color = char_to_index[key1][key2];
            pixel[0] = i2r[color];
            pixel[1] = i2g[color];
            pixel[2] = i2b[color];
            pixel[3] = i2a[color];
            pixel += 4;
        }
        line++;
    }
    gluBuild2DMipmaps(
        GL_TEXTURE_2D, GL_RGBA, width, height, GL_RGBA, GL_UNSIGNED_BYTE, rgba
    );
    free(rgba);
}

/*--------------------------------------------------------------------------------------------*/

#ifdef WITH_PNG
static GLboolean glut_viewer_load_image_PNG(
    const char* filename, GLuint* width_in, GLuint* height_in, GLuint* bpp, GLvoid** pixels_in
) {

    FILE* in;
    png_structp png_ptr;
    png_infop info_ptr;
    unsigned long width, height;
    int bit_depth, color_type, interlace_type;
    png_bytep row_pointer;
    unsigned int row;
    char** pixels = (char**) pixels_in;
    png_uint_32 w, h;

    in = fopen(filename, "rb");
    if(in == NULL) {
        fprintf(stderr, "PNG loader: failed opening file %s\n", filename);
        return GL_FALSE;
    }

    png_ptr = png_create_read_struct(
        PNG_LIBPNG_VER_STRING,
        (png_voidp) NULL, (png_error_ptr) NULL, (png_error_ptr) NULL
    );
    if(png_ptr == NULL) {
        fclose(in);
        fprintf(stderr, "PNG loader: error while loading %s\n", filename);
        return GL_FALSE;
    }

    info_ptr = png_create_info_struct(png_ptr);
    if(info_ptr == NULL) {
        png_destroy_read_struct(
            &png_ptr, (png_infopp) NULL, (png_infopp) NULL
        );
        fprintf(stderr, "PNG loader: error while loading %s\n", filename);
        fclose(in);
        return GL_FALSE;
    }

    png_init_io(png_ptr, in);

    /* read header  */
    png_read_info(png_ptr, info_ptr);
    png_get_IHDR(
        png_ptr, info_ptr, &w, &h, &bit_depth, &color_type,
        &interlace_type, NULL, NULL
    );
    width = w;
    height = h;

    if(color_type == PNG_COLOR_TYPE_GRAY) {
        *bpp = 1;
    } else if(color_type == PNG_COLOR_TYPE_RGB_ALPHA) {
        *bpp = 4;
    } else {
        *bpp = 3;
    }

    if(color_type == PNG_COLOR_TYPE_PALETTE) {
        png_set_expand(png_ptr);
    }

    *pixels = malloc(width * height * *bpp);
    if(*pixels == NULL) {
        png_destroy_read_struct(
            &png_ptr, (png_infopp) NULL, (png_infopp) NULL
        );
        fprintf(stderr, "PNG loader: error while loading %s, unable to allocate sufficient memory\n", filename);
        fclose(in);
        return GL_FALSE;
    }

    /* Read the image one line at a time */
    row_pointer = (png_bytep) malloc(width * *bpp);

    for(row = 0; row < height; row++) {
        png_read_rows(png_ptr, &row_pointer, NULL, 1);
        memcpy(
            *pixels + (height - 1 - row) * width * *bpp, row_pointer, width * *bpp
        );
    }
    free(row_pointer);

    png_read_end(png_ptr, info_ptr);
    png_destroy_read_struct(&png_ptr, &info_ptr, (png_infopp) NULL);

    *width_in = (GLuint) width;
    *height_in = (GLuint) height;

    fclose(in);
    return GL_TRUE;
}

#endif

/*--------------------------------------------------------------------------------------------*/

#define GLUT_VIEWER_BMP_RLE_8 1
#define GLUT_VIEWER_BMP_RLE_4 2

typedef unsigned int bmp_int32;
typedef unsigned short bmp_int16;

#define readH(x) fread(&(header.x), sizeof(header.x), 1, f);

typedef struct {
    bmp_int16 sType;
    bmp_int32 iSizeOfFile;
    bmp_int16 sReserved1;
    bmp_int16 sReserved2;
    bmp_int32 iOffBits;
    bmp_int32 iSize;
    bmp_int32 iWidth;
    bmp_int32 iHeight;
    bmp_int16 sPlanes;
    bmp_int16 sBitCount;
    bmp_int32 iCompression;
    bmp_int32 iSizeImage;
    bmp_int32 iXpelsPerMeter;
    bmp_int32 iYpelsPerMeter;
    bmp_int32 iClrUsed;
    bmp_int32 iClrImportant;
} GlutViewerBMPHeader;

static void rgb_to_bgr(GLuint width, GLuint height, GLuint bpp, GLvoid* pixels) {
    char* p = (char*) pixels;
    int i;
    char tmp;
    if(bpp != 3 && bpp != 4) {
        return;
    }
    for(i = 0; i < (int) (width * height); i++) {
        tmp = p[0];
        p[0] = p[2];
        p[2] = tmp;
        p += bpp;
    }
}

static GLboolean glut_viewer_load_image_BMP(
    const char* filename, GLuint* width, GLuint* height, GLuint* bpp, GLvoid** pixels_in
) {
    GlutViewerBMPHeader header;
    char padding[3];
    size_t rowlen;
    int i;
    char* pixels;
    FILE* f = fopen(filename, "rb");
    int ok = 1;

    if(f == NULL) {
        fprintf(stderr, "cannot open %s\n", filename);
        return GL_FALSE;
    }


   
    /*
     * Argh, I need to read the header field by field, seems
     * that the compiler aligns the fields in a way that makes
     * the memory layout mismatch the file layout.
     */
   
    memset(&header, 0, sizeof(GlutViewerBMPHeader));
    ok = ok && readH(sType);
    ok = ok && readH(iSizeOfFile);
    ok = ok && readH(sReserved1);
    ok = ok && readH(sReserved2);
    ok = ok && readH(iOffBits);
    ok = ok && readH(iSize);
    ok = ok && readH(iWidth);
    ok = ok && readH(iHeight);
    ok = ok && readH(sPlanes);
    ok = ok && readH(sBitCount);
    ok = ok && readH(iCompression);
    ok = ok && readH(iSizeImage);
    ok = ok && readH(iXpelsPerMeter);
    ok = ok && readH(iYpelsPerMeter);
    ok = ok && readH(iClrUsed);
    ok = ok && readH(iClrImportant);

    if(!ok) {
        fprintf(stderr, "Error while reading BMP header of %s\n", filename);
        fclose(f);
        return GL_FALSE;
    }

    if((header.iCompression & GLUT_VIEWER_BMP_RLE_8) || (header.iCompression & GLUT_VIEWER_BMP_RLE_4)) {
        fprintf(stderr, "BMP loader: cannot load compressed BMP files, sorry\n");
        fclose(f);
        return GL_FALSE;
    }
    *width = header.iWidth;
    *height = header.iHeight;
    *bpp = header.sBitCount / 8;

    pixels = malloc(*width * *height * *bpp);

    rowlen = (size_t) (*width * *bpp);
    for(i = 0; i < (int) *height; i++) {
        ok = ok && fread(pixels + i * (int) rowlen, rowlen, (size_t) 1, f);
        /* BMP rows are padded to 4-bytes multiples */
        ok = ok && fread(padding, rowlen % 4, (size_t) 1, f);
    }

    rgb_to_bgr(*width, *height, *bpp, pixels);
    *pixels_in = pixels;

    fclose(f);

    if(!ok) {
        fprintf(stderr, "Error while reading BMP data from %s\n", filename);
        free(pixels);
        return GL_FALSE;
    }

    return GL_TRUE;
}

/*--------------------------------------------------------------------------------------------*/

static const char* extension(const char* filename) {
    const char* point_loc = strrchr(filename, '.');
    return (point_loc == NULL) ? NULL : (point_loc + 1);
}

GLboolean glut_viewer_load_image(
    const char* filename, GLuint* width, GLuint* height, GLuint* bpp, GLvoid** pixels
) {
    const char* ext = extension(filename);
#ifdef WITH_PNG
    if(!strcmp(ext, "png") || !strcmp(ext, "PNG")) {
        return glut_viewer_load_image_PNG(filename, width, height, bpp, pixels);
    }
#endif
    if(!strcmp(ext, "bmp") || !strcmp(ext, "BMP")) {
        return glut_viewer_load_image_BMP(filename, width, height, bpp, pixels);
    }
    return GL_FALSE;
}

void glTexImage2Dfile(const char* filename) {
    GLuint width, height, bpp;
    GLvoid* pixels;
    if(!glut_viewer_load_image(filename, &width, &height, &bpp, &pixels)) {
        return;
    }
    switch(bpp) {
        case 1:
            gluBuild2DMipmaps(
                GL_TEXTURE_2D, GL_LUMINANCE,
                (GLsizei) width, (GLsizei) height, GL_LUMINANCE, GL_UNSIGNED_BYTE, pixels
            );
            break;
        case 3:
            gluBuild2DMipmaps(
                GL_TEXTURE_2D, GL_RGB, (GLsizei) width, (GLsizei) height, GL_RGB, GL_UNSIGNED_BYTE, pixels
            );
            break;
        case 4:
            gluBuild2DMipmaps(
                GL_TEXTURE_2D, GL_RGBA, (GLsizei) width, (GLsizei) height, GL_RGBA, GL_UNSIGNED_BYTE, pixels
            );
            break;
        default:
            glut_viewer_assert_not_reached();
            break;
    }
    free(pixels);
}

