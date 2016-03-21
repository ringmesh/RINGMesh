/*
 *  Copyright (c) 2012-2014, Bruno Levy
 *  All rights reserved.
 *
 *  Redistribution and use in source and binary forms, with or without
 *  modification, are permitted provided that the following conditions are met:
 *
 *  * Redistributions of source code must retain the above copyright notice,
 *  this list of conditions and the following disclaimer.
 *  * Redistributions in binary form must reproduce the above copyright notice,
 *  this list of conditions and the following disclaimer in the documentation
 *  and/or other materials provided with the distribution.
 *  * Neither the name of the ALICE Project-Team nor the names of its
 *  contributors may be used to endorse or promote products derived from this
 *  software without specific prior written permission.
 * 
 *  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 *  AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 *  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 *  ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 *  LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 *  CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 *  SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 *  INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 *  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 *  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 *  POSSIBILITY OF SUCH DAMAGE.
 *
 *  If you modify this software, you should include a notice giving the
 *  name of the person performing the modification, the date of modification,
 *  and the reason for such modification.
 *
 *  Contact: Bruno Levy
 *
 *     Bruno.Levy@inria.fr
 *     http://www.loria.fr/~levy
 *
 *     ALICE Project
 *     LORIA, INRIA Lorraine, 
 *     Campus Scientifique, BP 239
 *     54506 VANDOEUVRE LES NANCY CEDEX 
 *     FRANCE
 *
 */

#include <geogram_gfx/GLUP/GLUP.h>
#include <geogram_gfx/GLUP/GLUP_context_GLSL.h>
#include <geogram_gfx/GLUP/GLUP_context_VanillaGL.h>
#include <geogram_gfx/basic/GLSL.h>
#include <geogram/basic/logger.h>
#include <geogram/basic/command_line.h>

/*****************************************************************************/

namespace GLUP {
    using namespace GEO;
    static Context* current_context_ = nil;
    static std::set<Context*> all_contexts_;
    static bool initialized_ = false;
    void cleanup() {
        for(
            std::set<Context*>::iterator it=all_contexts_.begin();
            it != all_contexts_.end(); ++it
        ) {
            delete *it;
        }
        all_contexts_.clear();
    }
    
}

/*****************************************************************************/

const char* glupUniformStateDeclaration() {
    return GLUP::current_context_->uniform_state_declaration();
}

void GLUP_API glupBindUniformState(GLUPuint program) {
    GLUP::current_context_->bind_uniform_state(program);
}


/**
 * \brief Tests whether tessellation shaders are supported by OpenGL.
 * \details Some drivers may declare to be OpenGL 4.5 compliant whereas
 *  they do not have tesselation shader (for instance, I have an old
 *  NVidia quadro that does that...)
 * \retval true if tessellation shaders are supported
 * \retval false otherwise
 */
bool supports_tessellation_shader() {
    bool result = true;
    GLuint s_handle = glCreateShader(GL_TESS_CONTROL_SHADER);
    result = result && (s_handle != 0);
    if (s_handle != 0) {
        glDeleteShader(s_handle);
    }

    // Clear OpenGL error flag.
    while(glGetError() != GL_NO_ERROR) {
    }

    return result;
}

GLUPcontext glupCreateContext() {

    if(!GLUP::initialized_) {
        GLUP::initialized_ = true;
        atexit(GLUP::cleanup);
    }
    
    std::string GLUP_profile = GEO::CmdLine::get_arg("gfx:GLUP_profile");
    GLUP::Context* result = nil;

    if(GLUP_profile == "auto") {
        double GLSL_version = GEO::GLSL::supported_language_version();
        const GLubyte* vendor = glGetString(GL_VENDOR);
        if(!GEO::String::string_starts_with(
               std::string((const char*)vendor), "NVIDIA")
        ) {
            GEO::Logger::out("GLUP") << "Non-NVIDIA GPU" << std::endl;

            if(GEO::CmdLine::get_arg("gfx:GL_profile") == "compatibility") {
                GEO::Logger::out("GLUP")
                    << "Switching to VanillaGL" << std::endl;
                GEO::Logger::out("GLUP")
                    << "Use gfx:GLUP_profile to override"
                    << std::endl;
                GLSL_version = 0.0;
            } else {
                GEO::Logger::warn("GLUP")
                    << "Cannot switch to VanillaGL" << std::endl;
                GEO::Logger::warn("GLUP")
                    << "Needs gfx:GL_profile=compatibility" << std::endl;
                GEO::Logger::warn("GLUP")
                    << "(trying anyway with GLUP150/GLUP440)" << std::endl;
            }
        }
        if (GLSL_version >= 4.4) {
            if (!supports_tessellation_shader()) {
                GEO::Logger::out("GLUP")
                    << "GLSL version >= 4.4 but tessellation unsupported"
                    << std::endl;
                GEO::Logger::out("GLUP") << "Downgrading to GLUP 150..."
                                         << std::endl;
                GLSL_version = 1.5;
            }
        }
        if(GLSL_version < 1.5) {
            GLUP_profile = "VanillaGL";
        } else if(GLSL_version < 4.4) {
            GLUP_profile = "GLUP150";
        } else {
            GLUP_profile = "GLUP440";
        }
    }

    GEO::Logger::out("GLUP") << "Using " << GLUP_profile << " profile"
                        << std::endl;
    
    if(GLUP_profile == "GLUP150") {
        result = new GLUP::Context_GLSL150;        
    } else if(GLUP_profile == "GLUP440") {
        result = new GLUP::Context_GLSL440;
    } else if(GLUP_profile == "VanillaGL") {
        result = new GLUP::Context_VanillaGL;
    } else {
        GEO::Logger::warn("GLUP")
            << GLUP_profile << "unknown profile, falling back to VanillaGL"
            << std::endl;
        result = new GLUP::Context_VanillaGL;        
    }
    
    
    try {
        result->setup();
    } catch(...) {
        GEO::Logger::warn("GLUP")
            << "Caught an exception, downgrading to VanillaGL"
            << std::endl;
        delete result;
        result = new GLUP::Context_VanillaGL;
        result->setup();
    }

    GLUP::all_contexts_.insert(result);
    
    return result;
}

void glupDeleteContext(GLUPcontext context_in) {
    GLUP::Context* context =
        reinterpret_cast<GLUP::Context*>(context_in);

    std::set<GLUP::Context*>::iterator it = GLUP::all_contexts_.find(context);
    geo_assert(it != GLUP::all_contexts_.end());
    GLUP::all_contexts_.erase(it);
    
    if(GLUP::current_context_ == context) {
        GLUP::current_context_ = nil;
    }
    delete context;
}


GLUPcontext glupCurrentContext() {
    return GLUP::current_context_;
}

const char* glupCurrentProfileName() {
    return GLUP::current_context_->profile_name();
}

void glupMakeCurrent(GLUPcontext context) {
    GLUP::current_context_ = reinterpret_cast<GLUP::Context*>(context);
}

void glupCopyFromGLState(GLUPbitfield which_attributes) {
    GLUP::current_context_->copy_from_GL_state(which_attributes);
}

void glupCopyToGLState(GLUPbitfield which_attributes) {
    GLUP::current_context_->copy_to_GL_state(which_attributes);
}

GLUPboolean glupPrimitiveSupportsArrayMode(GLUPprimitive prim) {
    return GLUP::current_context_->primitive_supports_array_mode(prim) ?
        GL_TRUE : GL_FALSE ;
}

/****************** Enable / Disable ***************************/

void glupEnable(GLUPtoggle toggle) {
    GLUP::current_context_->uniform_state().toggle[toggle].set(GL_TRUE);
}

void glupDisable(GLUPtoggle toggle) {
    GLUP::current_context_->uniform_state().toggle[toggle].set(GL_FALSE);
}

GLUPboolean glupIsEnabled(GLUPtoggle toggle) {
    return GLUP::current_context_->uniform_state().toggle[toggle].get();
}

/********************** Texturing ******************************/

void GLUP_API glupTextureType(GLUPtextureType type) {
    GLUP::current_context_->uniform_state().texture_type.set(type);
}

GLUPtextureType glupGetTextureType() {
    return GLUPtextureType(
        GLUP::current_context_->uniform_state().texture_type.get()
    );
}

void GLUP_API glupTextureMode(GLUPtextureMode mode) {
    GLUP::current_context_->uniform_state().texture_mode.set(mode);    
}

GLUPtextureMode glupGetTextureMode() {
    return GLUPtextureMode(
        GLUP::current_context_->uniform_state().texture_mode.get()
    );
}

/****************** Drawing state ******************************/

void glupSetColor4fv(GLUPcolor color, const GLUPfloat* rgba) {
    if(color == GLUP_FRONT_AND_BACK_COLOR) {
        glupSetColor4fv(GLUP_FRONT_COLOR, rgba);
        glupSetColor4fv(GLUP_BACK_COLOR, rgba);        
    } else {
        GLUP::current_context_->uniform_state().color[color].set(rgba);
    }
}

void glupGetColor4fv(GLUPcolor color, float* rgba) {
    geo_assert(color != GLUP_FRONT_AND_BACK_COLOR);
    GLUP::current_context_->uniform_state().color[color].get(rgba);
}

void glupSetColor3fv(GLUPcolor color, const GLUPfloat* rgba) {
    glupSetColor4f(color, rgba[0], rgba[1], rgba[2], 1.0);
}

void glupSetColor4f(
    GLUPcolor color, GLUPfloat r, GLUPfloat g, GLUPfloat b, GLUPfloat a
) {
    if(color == GLUP_FRONT_AND_BACK_COLOR) {
        glupSetColor4f(GLUP_FRONT_COLOR, r, g, b, a);
        glupSetColor4f(GLUP_BACK_COLOR, r, g, b, a);
    } else {
        GLUPfloat* ptr =
            GLUP::current_context_->uniform_state().color[color].get_pointer();
        ptr[0] = r;
        ptr[1] = g;
        ptr[2] = b;
        ptr[3] = a;
    }
}

void glupSetColor3f(GLUPcolor color, GLUPfloat r, GLUPfloat g, GLUPfloat b) {
    glupSetColor4f(color, r, g, b, 1.0f);
}

void glupSetColor4dv(GLUPcolor color, const GLUPdouble* rgba) {
    glupSetColor4f(
        color,
        GLUPfloat(rgba[0]),
        GLUPfloat(rgba[1]),
        GLUPfloat(rgba[2]),
        GLUPfloat(rgba[3])
    );
}

void glupSetColor3dv(GLUPcolor color, const GLUPdouble* rgba) {
    glupSetColor4f(
        color,
        GLUPfloat(rgba[0]),
        GLUPfloat(rgba[1]),
        GLUPfloat(rgba[2]),
        1.0f
    );
}

void glupSetColor4d(
    GLUPcolor color, GLUPdouble r, GLUPdouble g, GLUPdouble b, GLUPdouble a
) {
    glupSetColor4f(
        color,
        GLUPfloat(r),
        GLUPfloat(g),
        GLUPfloat(b),
        GLUPfloat(a)
    );
}

void glupSetColor3d(
    GLUPcolor color, GLUPdouble r, GLUPdouble g, GLUPdouble b
) {
    glupSetColor4f(
        color,
        GLUPfloat(r),
        GLUPfloat(g),
        GLUPfloat(b),
        1.0f
    );
}

void glupLightVector3f(GLUPfloat x, GLUPfloat y, GLUPfloat z) {
    GLUPfloat* ptr =
        GLUP::current_context_->uniform_state().light_vector.get_pointer();
    ptr[0] = x;
    ptr[1] = y;
    ptr[2] = z;
    GLUP::current_context_->flag_lighting_as_dirty();
}

void glupLightVector3fv(GLUPfloat* xyz) {
    GLUP::current_context_->uniform_state().light_vector.set(xyz);
    GLUP::current_context_->flag_lighting_as_dirty();
}

void glupSetMeshWidth(GLUPint width) {
    GLUP::current_context_->uniform_state().mesh_width.set(GLfloat(width));
}

GLUPint glupGetMeshWidth() {
    return GLUPint(GLUP::current_context_->uniform_state().mesh_width.get());
}

void glupSetCellsShrink(GLUPfloat x) {
    x = GEO::geo_min(x, 1.0f);
    x = GEO::geo_max(x, 0.0f);
    GLUP::current_context_->uniform_state().cells_shrink.set(x);
}

GLUPfloat glupGetCellsShrink() {
    return GLUP::current_context_->uniform_state().cells_shrink.get();    
}

/****************** Picking ******************************/

void glupPickingMode(GLUPpickingMode mode) {
    GLUP::current_context_->uniform_state().picking_mode.set(mode);
}

GLUPpickingMode glupGetPickingMode() {
    return GLUPpickingMode(
        GLUP::current_context_->uniform_state().picking_mode.get()
    );    
}

void glupPickingId(GLUPuint64 id) {
    // TODO: uint64
    GLUP::current_context_->uniform_state().picking_id.set(GLint(id));
}

GLUPuint64 glupGetPickingId() {
    // TODO: uint64    
    return GLUPuint64(
        GLUP::current_context_->uniform_state().picking_id.get()
    );
}

void glupBasePickingId(GLUPuint64 id) {
    // TODO: uint64
    GLUP::current_context_->uniform_state().base_picking_id.set(GLint(id));    
}

GLUPuint64 glupGetBasePickingId() {
    // TODO: uint64
    return GLUPuint64(
        GLUP::current_context_->uniform_state().base_picking_id.get()
    );
}

/****************** Clipping ******************************/

void glupClipMode(GLUPclipMode mode) {
    GLUP::current_context_->uniform_state().clipping_mode.set(mode);
}

GLUPclipMode glupGetClipMode() {
    return GLUPclipMode(
        GLUP::current_context_->uniform_state().clipping_mode.get()
    );
}

void glupClipPlane(const GLUPdouble* eqn_in) {
    const GLfloat* modelview =
        GLUP::current_context_->get_matrix(GLUP_MODELVIEW_MATRIX);
    GLfloat modelview_invert[16];
    if(!GLUP::invert_matrix(modelview_invert,modelview)) {
        GEO::Logger::warn("GLUP") << "Singular ModelView matrix"
                             << std::endl;
        GLUP::show_matrix(modelview);
    }
    GLfloat* state_world_clip_plane =
        GLUP::current_context_->uniform_state().world_clip_plane.get_pointer();
    GLfloat* state_clip_plane =
        GLUP::current_context_->uniform_state().clip_plane.get_pointer();
    for(GEO::index_t i=0; i<4; ++i) {
        state_world_clip_plane[i] = float(eqn_in[i]);
    }
    GLUP::mult_matrix_vector(
        state_clip_plane,modelview_invert,state_world_clip_plane
    );
}

void glupGetClipPlane(GLUPdouble* eqn) {
    const GLfloat* ptr =
        GLUP::current_context_->uniform_state().clip_plane.get_pointer();
    eqn[0] = GLdouble(ptr[0]);
    eqn[1] = GLdouble(ptr[1]);
    eqn[2] = GLdouble(ptr[2]);
    eqn[3] = GLdouble(ptr[3]);    
}

/******************* Matrices ***************************/


void glupMatrixMode(GLUPmatrix matrix) {
    GLUP::current_context_->set_matrix_mode(matrix);
}

GLUPmatrix glupGetMatrixMode() {
    return GLUP::current_context_->get_matrix_mode();
}

void glupPushMatrix() {
    GLUP::current_context_->push_matrix();
}

void glupPopMatrix() {
    GLUP::current_context_->pop_matrix();    
}

void glupGetMatrixdv(GLUPmatrix matrix, GLUPdouble* ptr) {
    for(GEO::index_t i=0; i<16; ++i) {
        ptr[i] = GLUPdouble(
            GLUP::current_context_->get_matrix(matrix)[i]
        );
    }
}

void glupGetMatrixfv(GLUPmatrix matrix, GLUPfloat* ptr) {
    for(GEO::index_t i=0; i<16; ++i) {
        GLUP::copy_vector(
            ptr, GLUP::current_context_->get_matrix(matrix), 16
        );
    }
}

void glupLoadIdentity() {
    GLUP::current_context_->load_identity();
}

void glupLoadMatrixf(const GLUPfloat* M) {
    GLUP::current_context_->load_matrix(M);
}

void glupLoadMatrixd(const GLUPdouble* M) {
    GLfloat Mf[16];
    for(GEO::index_t i=0; i<16; ++i) {
        Mf[i] = GLfloat(M[i]);
    }
    glupLoadMatrixf(Mf);
}    

void glupMultMatrixf(const GLUPfloat* M) {
    GLUP::current_context_->mult_matrix(M);
}

void glupMultMatrixd(const GLUPdouble* M) {
    GLfloat Mf[16];
    for(GEO::index_t i=0; i<16; ++i) {
        Mf[i] = GLfloat(M[i]);
    }
    glupMultMatrixf(Mf);
}    

void glupTranslatef(GLUPfloat x, GLUPfloat y, GLUPfloat z) {
    GLfloat M[16];

    M[4*0+0] = 1.0f;
    M[4*0+1] = 0.0f;
    M[4*0+2] = 0.0f;
    M[4*0+3] = x;

    M[4*1+0] = 0.0f;
    M[4*1+1] = 1.0f;
    M[4*1+2] = 0.0f;
    M[4*1+3] = y;

    M[4*2+0] = 0.0f;
    M[4*2+1] = 0.0f;
    M[4*2+2] = 1.0f;
    M[4*2+3] = z;

    M[4*3+0] = 0.0f;
    M[4*3+1] = 0.0f;
    M[4*3+2] = 0.0f;
    M[4*3+3] = 1.0f;

    GLUP::transpose_matrix(M);
    
    glupMultMatrixf(M);
}

void glupTranslated(GLUPdouble x, GLUPdouble y, GLUPdouble z) {
    glupTranslatef(GLfloat(x), GLfloat(y), GLfloat(z));
}

void glupScalef(GLUPfloat sx, GLUPfloat sy, GLUPfloat sz) {
    GLfloat M[16];

    M[4*0+0] = sx;
    M[4*0+1] = 0.0f;
    M[4*0+2] = 0.0f;
    M[4*0+3] = 0.0f;

    M[4*1+0] = 0.0f;
    M[4*1+1] = sy;
    M[4*1+2] = 0.0f;
    M[4*1+3] = 0.0f;

    M[4*2+0] = 0.0f;
    M[4*2+1] = 0.0f;
    M[4*2+2] = sz;
    M[4*2+3] = 0.0f;

    M[4*3+0] = 0.0f;
    M[4*3+1] = 0.0f;
    M[4*3+2] = 0.0f;
    M[4*3+3] = 1.0f;

    glupMultMatrixf(M);
}

void glupScaled(GLUPdouble sx, GLUPdouble sy, GLUPdouble sz) {
    glupScalef(GLfloat(sx), GLfloat(sy), GLfloat(sz));    
}

void glupRotatef(
    GLUPfloat angle, GLUPfloat x, GLUPfloat y, GLUPfloat z
) {
    GLfloat l = 1.0f / ::sqrtf(x*x+y*y+z*z);
    x *= l;
    y *= l;
    z *= l;
    GLfloat s = ::sinf(angle);
    GLfloat c = ::cosf(angle);
    GLfloat M[16];

    M[4*0+0] = x*x*(1.0f-c)+c;
    M[4*0+1] = x*y*(1.0f-c)-z*s;
    M[4*0+2] = x*z*(1.0f-c)+y*s;
    M[4*0+3] = 0.0f;

    M[4*1+0] = y*x*(1.0f-c)+z*s;
    M[4*1+1] = y*y*(1.0f-c)+c;
    M[4*1+2] = y*z*(1.0f-c)-x*s;
    M[4*1+3] = 0.0f;

    M[4*2+0] = z*x*(1.0f-c)-y*s;
    M[4*2+1] = z*y*(1.0f-c)+x*s;
    M[4*2+2] = z*z*(1.0f-c)+c;
    M[4*2+3] = 0.0f;

    M[4*3+0] = 0.0f;
    M[4*3+1] = 0.0f;
    M[4*3+2] = 0.0f;
    M[4*3+3] = 1.0f;

    GLUP::transpose_matrix(M);
    
    glupMultMatrixf(M);
}

void glupRotated(
    GLUPdouble angle, GLUPdouble x, GLUPdouble y, GLUPdouble z
) {
    glupRotatef(GLfloat(angle), GLfloat(x), GLfloat(y), GLfloat(z));
}


void glupOrtho(
    GLUPdouble left, GLUPdouble right,
    GLUPdouble bottom, GLUPdouble top,
    GLUPdouble nearVal, GLUPdouble farVal
) {
    GLfloat M[16];

    GLdouble tx = -(right+left)/(right-left);
    GLdouble ty = -(top+bottom)/(top-bottom);
    GLdouble tz = -(farVal+nearVal)/(farVal-nearVal);
    
    M[4*0+0] = GLfloat(2.0 / (right-left));
    M[4*0+1] = 0.0f;
    M[4*0+2] = 0.0f;
    M[4*0+3] = GLfloat(tx);

    M[4*1+0] = 0.0f;
    M[4*1+1] = GLfloat(2.0 / (top-bottom));
    M[4*1+2] = 0.0f;
    M[4*1+3] = GLfloat(ty);

    M[4*2+0] = 0.0f;
    M[4*2+1] = 0.0f;
    M[4*2+2] = GLfloat(-2.0 / (farVal - nearVal));
    M[4*2+3] = GLfloat(tz);

    M[4*3+0] = 0.0f;
    M[4*3+1] = 0.0f;
    M[4*3+2] = 0.0f;
    M[4*3+3] = 1.0f;
    
    GLUP::transpose_matrix(M);
    glupMultMatrixf(M);
}

void glupOrtho2D(
    GLUPdouble left, GLUPdouble right, GLUPdouble bottom, GLUPdouble top
) {
    glupOrtho(left, right, bottom, top, -1.0, 1.0);
}

void glupFrustum(
    GLUPdouble left, GLUPdouble right,
    GLUPdouble bottom, GLUPdouble top,
    GLUPdouble nearVal, GLUPdouble farVal
) {
    GLfloat M[16];

    GLdouble A = (right + left) / (right - left);
    GLdouble B = (top + bottom) / (top - bottom);
    GLdouble C = -(farVal + nearVal) / (farVal - nearVal);
    GLdouble D = -2.0*farVal*nearVal / (farVal - nearVal);
    
    M[4*0+0] = GLfloat(2.0 * nearVal / (right - left));
    M[4*0+1] = 0.0f;
    M[4*0+2] = GLfloat(A);
    M[4*0+3] = 0.0f;

    M[4*1+0] = 0.0f;
    M[4*1+1] = GLfloat(2.0 * nearVal / (top - bottom));
    M[4*1+2] = GLfloat(B);
    M[4*1+3] = 0.0f;

    M[4*2+0] = 0.0f;
    M[4*2+1] = 0.0f;
    M[4*2+2] = GLfloat(C);
    M[4*2+3] = GLfloat(D);

    M[4*3+0] =  0.0f;
    M[4*3+1] =  0.0f;
    M[4*3+2] = -1.0f;
    M[4*3+3] =  0.0f;

    GLUP::transpose_matrix(M);
    glupMultMatrixf(M);
} 

void glupPerspective(
    GLUPdouble fovy, GLUPdouble aspect,
    GLUPdouble zNear, GLUPdouble zFar
) {
    GLfloat M[16];
    
    double f = 1.0 / tan(fovy * M_PI / 180.0);

    M[4*0+0] = GLfloat(f / aspect);
    M[4*0+1] = 0.0f;
    M[4*0+2] = 0.0f;
    M[4*0+3] = 0.0f;

    M[4*1+0] = 0.0f;
    M[4*1+1] = GLfloat(f);
    M[4*1+2] = 0.0f;
    M[4*1+3] = 0.0f;

    M[4*2+0] = 0.0f;
    M[4*2+1] = 0.0f;
    M[4*2+2] = GLfloat((zFar+zNear)/(zNear-zFar));
    M[4*2+3] = GLfloat(2.0*zFar*zNear/(zNear-zFar));

    M[4*3+0] =  0.0f;
    M[4*3+1] =  0.0f;
    M[4*3+2] = -1.0f;
    M[4*3+3] =  0.0f;

    GLUP::transpose_matrix(M);    
    glupMultMatrixf(M);    
}

GLUPboolean glupUnProject(
    GLUPdouble winx, GLUPdouble winy, GLUPdouble winz,
    const GLUPdouble modelMatrix[16],
    const GLUPdouble projMatrix[16],
    const GLUPint viewport[4],
    GLUPdouble *objx, GLUPdouble *objy, GLUPdouble *objz
) {
    double modelviewproject[16];
    double modelviewproject_inv[16];
    GLUP::mult_matrices(modelviewproject, modelMatrix, projMatrix);
    if(!GLUP::invert_matrix(modelviewproject_inv, modelviewproject)) {
        return GL_FALSE;
    }

    double in[4];
    in[0] = winx;
    in[1] = winy;
    in[2] = winz;
    in[3] = 1.0;

    // Invert viewport transform
    in[0] = (in[0] - double(viewport[0])) / double(viewport[2]);
    in[1] = (in[1] - double(viewport[1])) / double(viewport[3]);

    // Map to [-1, 1]
    in[0] = in[0] * 2.0 - 1.0;
    in[1] = in[1] * 2.0 - 1.0;
    in[2] = in[2] * 2.0 - 1.0;

    GLUP::transpose_matrix(modelviewproject_inv);
    
    double out[4];
    GLUP::mult_matrix_vector(out, modelviewproject_inv, in);

    if(out[3] == 0.0) {
        return GL_FALSE;
    }

    *objx = out[0] / out[3];
    *objy = out[1] / out[3];
    *objz = out[2] / out[3];
    
    return GL_TRUE;
}

GLUPboolean glupInvertMatrixfv(
    GLUPfloat Minvert[16],    
    const GLUPfloat M[16]
) {
    return GLUP::invert_matrix(Minvert, M);
}

GLUPboolean glupInvertMatrixdv(
    GLUPdouble Minvert[16],    
    const GLUPdouble M[16]
) {
    return GLUP::invert_matrix(Minvert, M);
}



/******************* Drawing ***************************/

void glupDrawArrays(
    GLUPprimitive primitive, GLUPint first, GLUPsizei count
) {
    GLUP::current_context_->draw_arrays(
        primitive, first, count
    );
}
    
void glupDrawElements(
    GLUPprimitive primitive, GLUPsizei count,
    GLUPenum type, const GLUPvoid* indices
) {
    GLUP::current_context_->draw_elements(
        primitive, count, type, indices
    );
}

void glupBegin(GLUPprimitive primitive) {
    GLUP::current_context_->begin(primitive);    
}

void glupEnd() {
    GLUP::current_context_->end();
}

void glupVertex2fv(const GLUPfloat* xy) {
    GLUP::current_context_->immediate_vertex(xy[0], xy[1]);
}

void glupVertex3fv(const GLUPfloat* xyz) {
    GLUP::current_context_->immediate_vertex(xyz[0], xyz[1], xyz[2]);    
}

void glupVertex4fv(const GLUPfloat* xyzw) {
    GLUP::current_context_->immediate_vertex(
        xyzw[0], xyzw[1], xyzw[2], xyzw[3]
    );        
}

void glupVertex2dv(const GLUPdouble* xy) {
    GLUP::current_context_->immediate_vertex(
        GLfloat(xy[0]),
        GLfloat(xy[1])
    );
}

void glupVertex3dv(const GLUPdouble* xyz) {
    GLUP::current_context_->immediate_vertex(
        GLfloat(xyz[0]),
        GLfloat(xyz[1]),
        GLfloat(xyz[2])        
    );
}

void glupVertex4dv(const GLUPdouble* xyzw) {
    GLUP::current_context_->immediate_vertex(
        GLfloat(xyzw[0]),
        GLfloat(xyzw[1]),
        GLfloat(xyzw[2]),
        GLfloat(xyzw[3])                
    );
}

void glupVertex2f(GLUPfloat x, GLUPfloat y) {
    GLUP::current_context_->immediate_vertex(x,y);    
}        

void glupVertex3f(GLUPfloat x, GLUPfloat y, GLUPfloat z) {
    GLUP::current_context_->immediate_vertex(x,y,z);        
}    

void glupVertex4f(GLUPfloat x, GLUPfloat y, GLUPfloat z, GLUPfloat w) {
    GLUP::current_context_->immediate_vertex(x,y,z,w);            
}

void glupVertex2d(GLUPdouble x, GLUPdouble y) {
    GLUP::current_context_->immediate_vertex(
        GLfloat(x),
        GLfloat(y)
    );
}        

void glupVertex3d(GLUPdouble x, GLUPdouble y, GLUPdouble z) {
    GLUP::current_context_->immediate_vertex(
        GLfloat(x),
        GLfloat(y),
        GLfloat(z)
    );
}    

void glupVertex4d(GLUPdouble x, GLUPdouble y, GLUPdouble z, GLUPdouble w) {
    GLUP::current_context_->immediate_vertex(
        GLfloat(x),
        GLfloat(y),
        GLfloat(z),
        GLfloat(w)                
    );
}

void glupColor3fv(const GLUPfloat* rgb) {
    GLUP::current_context_->immediate_color(rgb[0], rgb[1], rgb[2]);
}

void glupColor4fv(const GLUPfloat* rgba) {
    GLUP::current_context_->immediate_color(rgba[0], rgba[1], rgba[2], rgba[3]);
}

void glupColor3dv(const GLUPdouble* rgb) {
    GLUP::current_context_->immediate_color(
        GLfloat(rgb[0]),
        GLfloat(rgb[1]),
        GLfloat(rgb[2])
    );    
}

void glupColor4dv(const GLUPdouble* rgba) {
    GLUP::current_context_->immediate_color(
        GLfloat(rgba[0]),
        GLfloat(rgba[1]),
        GLfloat(rgba[2]),
        GLfloat(rgba[3])        
    );    
}

void glupColor3f(GLUPfloat r, GLUPfloat g, GLUPfloat b) {
    GLUP::current_context_->immediate_color(r, g, b);    
}    

void glupColor4f(GLUPfloat r, GLUPfloat g, GLUPfloat b, GLUPfloat a) {
    GLUP::current_context_->immediate_color(r, g, b, a);        
}

void glupColor3d(GLUPdouble r, GLUPdouble g, GLUPdouble b) {
    GLUP::current_context_->immediate_color(
        GLfloat(r),
        GLfloat(g),
        GLfloat(b)
    );    
}    

void glupColor4d(GLUPdouble r, GLUPdouble g, GLUPdouble b, GLUPdouble a) {
    GLUP::current_context_->immediate_color(
        GLfloat(r),
        GLfloat(g),
        GLfloat(b),
        GLfloat(a)
    );    
}

void glupTexCoord2fv(const GLUPfloat* st) {
    GLUP::current_context_->immediate_tex_coord(st[0], st[1]);    
}

void glupTexCoord3fv(const GLUPfloat* stu) {
    GLUP::current_context_->immediate_tex_coord(stu[0], stu[1], stu[2]);        
}

void glupTexCoord4fv(const GLUPfloat* stuv) {
    GLUP::current_context_->immediate_tex_coord(
        stuv[0], stuv[1], stuv[2], stuv[3]
    );            
}

void glupTexCoord2dv(const GLUPdouble* st) {
    GLUP::current_context_->immediate_tex_coord(
        GLfloat(st[0]),
        GLfloat(st[1])
    );            
}

void glupTexCoord3dv(const GLUPdouble* stu) {
    GLUP::current_context_->immediate_tex_coord(
        GLfloat(stu[0]),
        GLfloat(stu[1]),
        GLfloat(stu[2])
    );            
}

void glupTexCoord4dv(const GLUPdouble* stuv) {
    GLUP::current_context_->immediate_tex_coord(
        GLfloat(stuv[0]),
        GLfloat(stuv[1]),
        GLfloat(stuv[2]),
        GLfloat(stuv[3])
    );            
}

void glupTexCoord1f(GLUPfloat s) {
    GLUP::current_context_->immediate_tex_coord(s);        
}

void glupTexCoord2f(GLUPfloat s, GLUPfloat t) {
    GLUP::current_context_->immediate_tex_coord(s,t);        
}        

void glupTexCoord3f(GLUPfloat s, GLUPfloat t, GLUPfloat u) {
    GLUP::current_context_->immediate_tex_coord(s,t,u);        
}    

void glupTexCoord4f(GLUPfloat s, GLUPfloat t, GLUPfloat u, GLUPfloat v) {
    GLUP::current_context_->immediate_tex_coord(s,t,u,v);    
}

void glupTexCoord1d(GLUPdouble s) {
    GLUP::current_context_->immediate_tex_coord(
        GLfloat(s)
    );        
}            

void glupTexCoord2d(GLUPdouble s, GLUPdouble t) {
    GLUP::current_context_->immediate_tex_coord(
        GLfloat(s),
        GLfloat(t)
    );        
}        

void glupTexCoord3d(GLUPdouble s, GLUPdouble t, GLUPdouble u) {
    GLUP::current_context_->immediate_tex_coord(
        GLfloat(s),
        GLfloat(t),
        GLfloat(u)
    );        
}    

void glupTexCoord4d(GLUPdouble s, GLUPdouble t, GLUPdouble u, GLUPdouble v) {
    GLUP::current_context_->immediate_tex_coord(
        GLfloat(s),
        GLfloat(t),
        GLfloat(u),
        GLfloat(v)
    );        
}

void glupUseProgram(GLUPuint program) {
    GLUP::current_context_->set_user_program(program);
}

/****************************************************************************/
