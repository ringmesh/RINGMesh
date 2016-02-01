/*
 *  Copyright (c) 2012-2016, Bruno Levy
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

#include <geogram_gfx/basic/GLUP_private.h>
#include <geogram_gfx/basic/GLSL.h>
#include <geogram/basic/command_line.h>
#include <geogram/basic/progress.h>

// TODO: uber shader framework: bug sur combinaison suivante:
// (meme avec GLSL440):
//      couleurs sommets off / texture on:
//           me fait du garbage (et ce sur toutes les primitives).
//  -> j'ai un "workaround" (dans ce cas l\`a, on laisse le shader
//  determiner dynamiquement que les couleurs ne sont pas utilis\'ees),
// voir si on peut faire mieux que \c{c}a...

// TODO: documenter un peu tout ca (en particulier le "vertex_gather_mode")

// http://github.prideout.net/modern-opengl-prezo/

// TODO: for points: early fragment tests,   depth greater
// (+ vertex shader qui "ramene le point devant" comme \c{c}a le
//  fragment shader "pousse le point" et on a l'early fragment test).
//
//  layout(early_fragment_tests) in;
//  layout(depth_greater) out float gl_FragDepth;

// TODO: implanter GLUP_POLYGON? (ou pas..., peut-etre plutot GLUP_TRIANGLE_FAN)
// TODO: implanter glup clip mode standard et slice
// TODO: finir le profile VanillaGL (du boulot...)
// TODO: Je me demande toujours si glBegin() / glEnd() n'existe pas encore,
//  avec les glVertexAttrib(). Peut-etre que j'ai fait mes ImmediateBuffers
//  pour rien... (mais je pense que de toute fa\c{c}on glBegin()/glEnd()
//  n'existe pas sous OpenGL ES, donc \c{c}a servira quand m\^eme...)


namespace GLUP {
    using namespace GEO;

    
    void show_matrix(const GLfloat M[16]) {
        for(index_t i=0; i<4; ++i) {
            for(index_t j=0; j<4; ++j) {
                std::cerr << M[4*i+j] << ' ';
            }
            std::cerr << std::endl;
        }
    }
    
    // Used to determine maximum vertex index, that needs
    // to be an integer multiple of the number of vertices
    // per primitive.
    index_t nb_vertices_per_primitive[] = {
        1, // GLUP_POINTS     =0,
        2, // GLUP_LINES      =1,
        3, // GLUP_TRIANGLES  =2,
        4, // GLUP_QUADS      =3,
        4, // GLUP_TETRAHEDRA =4,
        8, // GLUP_HEXAHEDRA  =5,
        6, // GLUP_PRISMS     =6,
        5  // GLUP_PYRAMIDS   =7
    };

    static index_t nb_vertices_per_GL_primitive(GLenum primitive) {
        index_t result = 0;
        switch(primitive) {
        case GL_POINTS:
            result = 1;
            break;
        case GL_LINES:
            result = 2;
            break;
        case GL_TRIANGLES:
            result = 3;
            break;
        case GL_LINES_ADJACENCY:
            result = 4;
            break;
        case GL_TRIANGLES_ADJACENCY:
            result = 6;
            break;
        default:
            geo_assert_not_reached;
            break;
        }
        return result;
    }
    
    const char* primitive_name[] = {
        "GLUP_POINTS",
        "GLUP_LINES",
        "GLUP_TRIANGLES",
        "GLUP_QUADS",
        "GLUP_TETRAHEDRA",
        "GLUP_HEXAHEDRA",
        "GLUP_PRISMS",
        "GLUP_PYRAMIDS"
    };
    
    
    static const char* GLUP_uniform_state_source = 
        "  layout(shared)                              \n"
        "  uniform GLUPStateBlock {                    \n"
        
        "     bool vertex_colors_enabled;              \n"        
        
        "     vec4  front_color;                       \n"
        "     vec4  back_color;                        \n"
        
        "     bool draw_mesh_enabled;                  \n"
        "     vec4  mesh_color;                        \n"
        "     float mesh_width;                        \n"
        
        "     bool lighting_enabled;                   \n"
        "     vec3 light_vector;                       \n"
        "     vec3 light_half_vector;                  \n"
        
        "     bool texturing_enabled;                  \n"
        "     int  texture_mode;                       \n"
        "     int  texture_type;                       \n"        
        
        "     float cells_shrink;                      \n"
        
        "     bool picking_enabled;                    \n"         
        "     int   picking_mode;                      \n"
        "     int   picking_id;                        \n" 
        "     int   base_picking_id;                   \n" 

        "     bool clipping_enabled;                   \n"
        "     int   clipping_mode;                     \n"
        "     vec4  clip_plane;                        \n"
        
        "     mat4 modelviewprojection_matrix;         \n"
        "     mat4 modelview_matrix;                   \n"
        "     mat3 normal_matrix;                      \n"
        "     mat4 texture_matrix;                     \n"
        "  } GLUP;                                     \n"
        "                                              \n"
        "  const int GLUP_CLIP_STANDARD         = 1;   \n"
        "  const int GLUP_CLIP_WHOLE_CELLS      = 2;   \n"
        "  const int GLUP_CLIP_STRADDLING_CELLS = 3;   \n"
        "  const int GLUP_CLIP_SLICE_CELLS      = 4;   \n"
        "                                              \n"
        "  const int GLUP_TEXTURE_1D = 1;              \n"
        "  const int GLUP_TEXTURE_2D = 2;              \n"
        "  const int GLUP_TEXTURE_3D = 3;              \n"
        "                                              \n"
        "  const int GLUP_TEXTURE_REPLACE  = 0;        \n"
        "  const int GLUP_TEXTURE_MODULATE = 1;        \n"
        "  const int GLUP_TEXTURE_ADD      = 2;        \n"
        ;
    
    
    /*
    ** Invert 4x4 matrix.
    ** Contributed by David Moore (See Mesa bug #6748)
    */
    GLboolean invert_matrix(const GLfloat m[16], GLfloat inv[16]) {
        
        inv[0]  = m[5]*m[10]*m[15] - m[5]*m[11]*m[14] - m[9]*m[6]*m[15]
                + m[9]*m[7]*m[14] + m[13]*m[6]*m[11] - m[13]*m[7]*m[10];
        inv[4]  = -m[4]*m[10]*m[15] + m[4]*m[11]*m[14] + m[8]*m[6]*m[15]
                - m[8]*m[7]*m[14] - m[12]*m[6]*m[11] + m[12]*m[7]*m[10];
        inv[8]  = m[4]*m[9]*m[15] - m[4]*m[11]*m[13] - m[8]*m[5]*m[15]
                + m[8]*m[7]*m[13] + m[12]*m[5]*m[11] - m[12]*m[7]*m[9];
        inv[12] = -m[4]*m[9]*m[14] + m[4]*m[10]*m[13] + m[8]*m[5]*m[14]
                - m[8]*m[6]*m[13] - m[12]*m[5]*m[10] + m[12]*m[6]*m[9];
        inv[1]  = -m[1]*m[10]*m[15] + m[1]*m[11]*m[14] + m[9]*m[2]*m[15]
                - m[9]*m[3]*m[14] - m[13]*m[2]*m[11] + m[13]*m[3]*m[10];
        inv[5]  = m[0]*m[10]*m[15] - m[0]*m[11]*m[14] - m[8]*m[2]*m[15]
                + m[8]*m[3]*m[14] + m[12]*m[2]*m[11] - m[12]*m[3]*m[10];
        inv[9]  = -m[0]*m[9]*m[15] + m[0]*m[11]*m[13] + m[8]*m[1]*m[15]
                - m[8]*m[3]*m[13] - m[12]*m[1]*m[11] + m[12]*m[3]*m[9];
        inv[13] = m[0]*m[9]*m[14] - m[0]*m[10]*m[13] - m[8]*m[1]*m[14]
                + m[8]*m[2]*m[13] + m[12]*m[1]*m[10] - m[12]*m[2]*m[9];
        inv[2]  = m[1]*m[6]*m[15] - m[1]*m[7]*m[14] - m[5]*m[2]*m[15]
                + m[5]*m[3]*m[14] + m[13]*m[2]*m[7] - m[13]*m[3]*m[6];
        inv[6]  = -m[0]*m[6]*m[15] + m[0]*m[7]*m[14] + m[4]*m[2]*m[15]
                - m[4]*m[3]*m[14] - m[12]*m[2]*m[7] + m[12]*m[3]*m[6];
        inv[10] = m[0]*m[5]*m[15] - m[0]*m[7]*m[13] - m[4]*m[1]*m[15]
                + m[4]*m[3]*m[13] + m[12]*m[1]*m[7] - m[12]*m[3]*m[5];
        inv[14] = -m[0]*m[5]*m[14] + m[0]*m[6]*m[13] + m[4]*m[1]*m[14]
                - m[4]*m[2]*m[13] - m[12]*m[1]*m[6] + m[12]*m[2]*m[5];
        inv[3]  = -m[1]*m[6]*m[11] + m[1]*m[7]*m[10] + m[5]*m[2]*m[11]
                - m[5]*m[3]*m[10] - m[9]*m[2]*m[7] + m[9]*m[3]*m[6];
        inv[7]  = m[0]*m[6]*m[11] - m[0]*m[7]*m[10] - m[4]*m[2]*m[11]
                + m[4]*m[3]*m[10] + m[8]*m[2]*m[7] - m[8]*m[3]*m[6];
        inv[11] = -m[0]*m[5]*m[11] + m[0]*m[7]*m[9] + m[4]*m[1]*m[11]
                - m[4]*m[3]*m[9] - m[8]*m[1]*m[7] + m[8]*m[3]*m[5];
        inv[15] = m[0]*m[5]*m[10] - m[0]*m[6]*m[9] - m[4]*m[1]*m[10]
                + m[4]*m[2]*m[9] + m[8]*m[1]*m[6] - m[8]*m[2]*m[5];
        
        GLfloat det =
            m[0]*inv[0] + m[1]*inv[4] + m[2]*inv[8] + m[3]*inv[12];
        
        if (det == 0.0f) {
            return GL_FALSE;
        }
        
        det = 1.0f / det;
        
        for (index_t i = 0; i < 16; ++i) {
            inv[i] *= det;
        }
        
        return GL_TRUE;
    }
    
    void mult_matrices(
        const GLfloat m1[16], const GLfloat m2[16], GLfloat out[16]
    ) {
        Memory::clear(out, sizeof(GLfloat)*16);
        for(index_t i=0; i<4; ++i) {
            for(index_t j=0; j<4; ++j) {
                for(index_t k=0; k<4; ++k) {
                    out[i*4+j] += m1[i*4+k]*m2[k*4+j];
                }
            }
        }
    }

    void mult_matrix_vector(
        const GLfloat m[16], const GLfloat v[4], GLfloat out[4]
    ) {
        Memory::clear(out, sizeof(GLfloat)*4);
        for(index_t i=0; i<4; ++i) {
            for(index_t j=0; j<4; ++j) {
                out[i] += v[j] * m[4*i+j];
            }
        }
    }

    
    void load_identity_matrix(GLfloat out[16]) {
        for(index_t i=0; i<4; ++i) {
            for(index_t j=0; j<4; ++j) {
                out[i*4+j] = (i == j) ? 1.0f : 0.0f;
            }
        }
    }

    /***********************************************************/
    
    void StateVariableBase::initialize(
        Context* context, const char* name
    ) {
        context_ = context;
        name_ = name;
        address_ = context_->get_state_variable_address(name);
    }

    void StateVariableBase::flag_uniform_buffer_as_dirty() {
        context_->flag_uniform_buffer_as_dirty();
    }

    /***********************************************************/

    const char* Context::uniform_state_declaration() {
        return GLUP_uniform_state_source;
    }
    
    Context::Context() {
        matrices_dirty_=true;        
        default_program_ = 0;
        uniform_buffer_=0;
        uniform_binding_point_=0;
        uniform_buffer_size_=0;
        uniform_buffer_data_=nil;
        uniform_buffer_dirty_=true;
        
        matrix_mode_ = GLUP_MODELVIEW_MATRIX;
        matrices_dirty_ = true;

        precompile_shaders_ =
            CmdLine::get_arg_bool("gfx:GLUP_precompile_shaders");

        use_core_profile_ =
            CmdLine::get_arg_bool("gfx:GLUP_use_core_profile");
    }
    
    Context::~Context() {
        if(default_program_ != 0) {
            glDeleteProgram(default_program_);
        }
        glDeleteBuffers(1, &uniform_buffer_);
        delete[] uniform_buffer_data_;                
    }

    bool Context::primitive_supports_array_mode(GLUPprimitive prim) const {
        return primitive_info_[prim].implemented &&
            !primitive_info_[prim].vertex_gather_mode;
    }
    
    void Context::setup() {
        
        uniform_buffer_dirty_=true;
        matrices_dirty_=true;        
        uniform_binding_point_=1;
        
        // A minimalistic GLSL program that uses the GLUP context.
        // It is there just to use GLSL introspection API to lookup
        // the offsets of GLUP context state variables.

        static const char* shader_source_header_ =
            "#version 150 core \n" ;
        
        static const char* vertex_shader_source_ =
            "in vec3 position;                                        \n"
            "                                                         \n"
            "void main() {                                            \n"
            "   gl_Position =                                         \n"
            "     GLUP.modelviewprojection_matrix*vec4(position,1.0); \n"
            "}                                                        \n";

        static const char* fragment_shader_source_ =
            "out vec4 colorOut;                      \n"
            "void main() {                           \n"
            "   colorOut = vec4(1.0, 1.0, 1.0, 1.0); \n" 
            "}                                       \n";

        // Create a GLSL program that uses the UBO

        GLuint GLUP_vertex_shader = GLSL::compile_shader(
            GL_VERTEX_SHADER,
            shader_source_header_,
            GLUP_uniform_state_source,
            vertex_shader_source_,
            0
        );

        GLuint GLUP_fragment_shader = GLSL::compile_shader(
            GL_FRAGMENT_SHADER,
            shader_source_header_,
            fragment_shader_source_,
            0
        );

        default_program_ = GLSL::create_program_from_shaders(
            GLUP_vertex_shader,
            GLUP_fragment_shader,
            0
        );

        // Get UBO size

        GLuint UBO_index =
            glGetUniformBlockIndex(default_program_, "GLUPStateBlock");
        
        glUniformBlockBinding(
            default_program_, UBO_index, uniform_binding_point_
        );

        glGetActiveUniformBlockiv(
            default_program_, UBO_index,
            GL_UNIFORM_BLOCK_DATA_SIZE,
            &uniform_buffer_size_
        );

        // Create UBO

        uniform_buffer_data_ = new Memory::byte[uniform_buffer_size_];
        Memory::clear(uniform_buffer_data_, size_t(uniform_buffer_size_));
        glGenBuffers(1, &uniform_buffer_);
        glBindBuffer(GL_UNIFORM_BUFFER, uniform_buffer_);
        glBufferData(
            GL_UNIFORM_BUFFER,
            uniform_buffer_size_,
            uniform_buffer_data_,
            GL_DYNAMIC_DRAW
        );
                
        glBindBufferBase(
            GL_UNIFORM_BUFFER,
            uniform_binding_point_,
            uniform_buffer_
        );
        glBindBuffer(GL_UNIFORM_BUFFER, 0);

        setup_state_variables();
        setup_immediate_buffers();
        
        //  Indicate that all toggle states should be
        // read from the state (default mode, superseded later).
        setup_shaders_source_for_toggles(0,~0);
        setup_primitives();
    }

    void Context::setup_state_variables() {
        uniform_state_.toggle.push_back(
            StateVariable<GLboolean>(this,"lighting_enabled",GL_TRUE)
        );
        uniform_state_.toggle.push_back(
            StateVariable<GLboolean>(this,"vertex_colors_enabled",GL_FALSE)
        );
        uniform_state_.toggle.push_back(
            StateVariable<GLboolean>(this,"texturing_enabled",GL_FALSE)
        );
        uniform_state_.toggle.push_back(
            StateVariable<GLboolean>(this,"draw_mesh_enabled",GL_FALSE)
        );
        uniform_state_.toggle.push_back(
            StateVariable<GLboolean>(this,"clipping_enabled",GL_FALSE)
        );
        uniform_state_.toggle.push_back(
            StateVariable<GLboolean>(this,"picking_enabled",GL_FALSE)
        );

        uniform_state_.color.push_back(
            VectorStateVariable(this, "front_color", 4)
        );
        uniform_state_.color.push_back(
            VectorStateVariable(this, "back_color", 4)
        );
        uniform_state_.color.push_back(
            VectorStateVariable(this, "mesh_color", 4)
        );
        
        uniform_state_.light_vector.initialize(this, "light_vector", 3);
        uniform_state_.light_half_vector.initialize(
            this, "light_half_vector", 3
        );
        
        uniform_state_.mesh_width.initialize(this, "mesh_width", 1.0f);
        uniform_state_.cells_shrink.initialize(this, "cells_shrink", 0.0f);

        uniform_state_.picking_mode.initialize(
            this, "picking_mode", GLUP_PICK_PRIMITIVE
        );
        uniform_state_.picking_id.initialize(this, "picking_id", 0);
        uniform_state_.base_picking_id.initialize(this, "base_picking_id", 0);

        uniform_state_.clipping_mode.initialize(
            this, "clipping_mode", GLUP_CLIP_STANDARD
        );
        uniform_state_.clip_plane.initialize(this, "clip_plane", 4);

        uniform_state_.texture_mode.initialize(
            this, "texture_mode", GLUP_TEXTURE_MODULATE
        );

        uniform_state_.texture_type.initialize(
            this, "texture_type", GLUP_TEXTURE_2D
        );
        
        uniform_state_.modelview_matrix.initialize(this, "modelview_matrix");
        uniform_state_.modelviewprojection_matrix.initialize(
            this, "modelviewprojection_matrix"
        );
        uniform_state_.normal_matrix.initialize(this, "normal_matrix");
        uniform_state_.texture_matrix.initialize(this, "texture_matrix");
        
        matrix_mode_ = GLUP_MODELVIEW_MATRIX;
        update_matrices();
    }

    void Context::setup_immediate_buffers() {

        glGenVertexArrays(1,&immediate_state_.VAO());
        glBindVertexArray(immediate_state_.VAO());

        for(index_t i=0; i<immediate_state_.buffer.size(); ++i) {
            update_buffer_object(
                immediate_state_.buffer[i].VBO(),
                GL_ARRAY_BUFFER,
                immediate_state_.buffer[i].size_in_bytes(),
                nil // no need to copy the buffer, it will be overwritten after.
            );
            // Note: the VAO is still bound.
            glVertexAttribPointer(
                i,
                GLint(immediate_state_.buffer[i].dimension()),
                GL_FLOAT,
                GL_FALSE,
                0,  // stride
                0   // pointer (relative to bound VBO beginning)
            );
        }
        
        glBindVertexArray(0);
        glBindBuffer(GL_ARRAY_BUFFER,0);
    }

    void Context::stream_immediate_buffers() {
        for(index_t i=0; i<immediate_state_.buffer.size(); ++i) {
            if(
                immediate_state_.buffer[i].is_enabled() &&
                immediate_state_.buffer[i].VBO() != 0
            ) {
                update_buffer_object(
                    immediate_state_.buffer[i].VBO(),
                    GL_ARRAY_BUFFER,
                    immediate_state_.buffer[i].size_in_bytes(),
                    immediate_state_.buffer[i].data(),
                    true // streaming enabled.
                );
            }
        }
        glBindBuffer(GL_ARRAY_BUFFER,0);        
    }
    
    Memory::pointer Context::get_state_variable_address(const char* name_in) {
        std::string name = std::string("GLUPStateBlock.") + name_in;
        GLuint index = GL_INVALID_INDEX;
        const char* pname = name.c_str();
        glGetUniformIndices(default_program_, 1, &pname, &index);
        geo_assert(index != GL_INVALID_INDEX);
        GLint offset;
        glGetActiveUniformsiv(
            default_program_, 1, &index, GL_UNIFORM_OFFSET, &offset
        );
        geo_assert(offset != -1);
        return uniform_buffer_data_ + offset;
    }

    void Context::copy_from_GL_state(GLUPbitfield which_attributes) {
        // No "GL state" in core profile.
        if(use_core_profile_) {
            return;
        }
        
        uniform_buffer_dirty_ = true;
        
        //  I observed that without these instructions,
        // glGet() does not always
        // returns the latest set values under Windows.
        glFlush();
        glFinish();

        if(which_attributes & GLUP_MATRICES_ATTRIBUTES_BIT) {
            matrices_dirty_ = true;            
            glGetFloatv(
                GL_PROJECTION_MATRIX,
                matrix_stack_[GLUP_PROJECTION_MATRIX].top()
            );
            glGetFloatv(
                GL_MODELVIEW_MATRIX,
                matrix_stack_[GLUP_MODELVIEW_MATRIX].top()
            );
            glGetFloatv(
                GL_TEXTURE_MATRIX,
                matrix_stack_[GLUP_TEXTURE_MATRIX].top()
            );
        }

        if(which_attributes & GLUP_CLIPPING_ATTRIBUTES_BIT) {
            uniform_state_.toggle[GLUP_CLIPPING].set(
                glIsEnabled(GL_CLIP_PLANE0)
            );
            GLdouble clip_plane_d[4];
            glGetClipPlane(GL_CLIP_PLANE0, clip_plane_d);
            copy_vector(
                uniform_state_.clip_plane.get_pointer(), clip_plane_d, 4
            );
        }
        
        if(which_attributes & GLUP_LIGHTING_ATTRIBUTES_BIT) {
            lighting_dirty_ = true;
            uniform_state_.toggle[GLUP_LIGHTING].set(glIsEnabled(GL_LIGHTING));
            GLfloat light[4];
            glGetLightfv(GL_LIGHT0, GL_POSITION, light);
            copy_vector(uniform_state_.light_vector.get_pointer(), light, 3);
        }            

        if(which_attributes & GLUP_COLORS_ATTRIBUTES_BIT) {
            glGetMaterialfv(
                GL_FRONT, GL_DIFFUSE,
                uniform_state_.color[GLUP_FRONT_COLOR].get_pointer()
            );
            glGetMaterialfv(
                GL_BACK, GL_DIFFUSE,
                uniform_state_.color[GLUP_BACK_COLOR].get_pointer()
            );
        }
    }

    void Context::copy_to_GL_state(GLUPbitfield which_attributes) {

        // No "GL state" in core profile.
        if(use_core_profile_) {
            return;
        }
        
        if(which_attributes & GLUP_MATRICES_ATTRIBUTES_BIT) {
            GLint mode_save;
            glGetIntegerv(GL_MATRIX_MODE, &mode_save);
            glMatrixMode(GL_PROJECTION);
            glLoadMatrixf(matrix_stack_[GLUP_PROJECTION_MATRIX].top());
            glMatrixMode(GL_MODELVIEW);
            glLoadMatrixf(matrix_stack_[GLUP_MODELVIEW_MATRIX].top());
            glMatrixMode(GL_TEXTURE);
            glLoadMatrixf(matrix_stack_[GLUP_TEXTURE_MATRIX].top());
            glMatrixMode(GLenum(mode_save));
        }

        if(which_attributes & GLUP_CLIPPING_ATTRIBUTES_BIT) {
            if(uniform_state_.toggle[GLUP_CLIPPING].get()) {
                glEnable(GL_CLIP_PLANE0);
            } else {
                glDisable(GL_CLIP_PLANE0);                
            }
            GLdouble clip_plane_d[4];
            glGetClipPlane(GL_CLIP_PLANE0, clip_plane_d);
            copy_vector(
                clip_plane_d, uniform_state_.clip_plane.get_pointer(), 4
            );            
            glClipPlane(GL_CLIP_PLANE0, clip_plane_d);
        }

        if(which_attributes & GLUP_LIGHTING_ATTRIBUTES_BIT) {
            if(uniform_state_.toggle[GLUP_LIGHTING].get()) {
                glEnable(GL_LIGHTING);
            } else {
                glDisable(GL_LIGHTING);
            }
            GLfloat light[4];
            copy_vector(light, uniform_state_.light_vector.get_pointer(), 3);
            light[3] = 0.0f;
            glLightfv(GL_LIGHT0, GL_POSITION, light);
        }

        if(which_attributes & GLUP_COLORS_ATTRIBUTES_BIT) {
            glMaterialfv(
                GL_FRONT, GL_DIFFUSE,
                uniform_state_.color[GLUP_FRONT_COLOR].get_pointer()
            );
            glMaterialfv(
                GL_BACK, GL_DIFFUSE,
                uniform_state_.color[GLUP_BACK_COLOR].get_pointer()
            );
        }
    }
    
    void Context::bind_uniform_state(GLuint program) {
        GLuint UBO_index = glGetUniformBlockIndex(
            program, "GLUPStateBlock"
        );
        if(UBO_index != GL_INVALID_INDEX) {
            glUniformBlockBinding(
                program, UBO_index, uniform_binding_point_
            );
        }
    }
    
    void Context::begin(GLUPprimitive primitive) {
        update_toggles_config();
        create_program_if_needed(primitive);
        if(!primitive_info_[primitive].implemented) {
            Logger::warn("GLUP")
                << "glupBegin(): "
                << primitive_name[primitive]
                << " not implemented in this profile" << std::endl;
        }
        
        update_uniform_buffer();

        //   If the primitive has a special VAO to be used for immediate
        // mode, then bind it.
        if(primitive_info_[primitive].vertex_gather_mode) {
            glBindVertexArray(
                primitive_info_[primitive].vertex_gather_mode_VAO
            );            
        } else {
            // Else use the regular VAO used by all immediate-mode primitives
            // (if there is one).
            if(immediate_state_.VAO() != 0) {
                glBindVertexArray(immediate_state_.VAO());
            }
        }

        if(uniform_state_.toggle[GLUP_VERTEX_COLORS].get()) {
            immediate_state_.buffer[GLUP_COLOR_ATTRIBUTE].enable();
        } else {
            immediate_state_.buffer[GLUP_COLOR_ATTRIBUTE].disable();
        }
        
        if(uniform_state_.toggle[GLUP_TEXTURING].get()) {
            immediate_state_.buffer[GLUP_TEX_COORD_ATTRIBUTE].enable();
        } else {
            immediate_state_.buffer[GLUP_TEX_COORD_ATTRIBUTE].disable();
        }
        
        immediate_state_.begin(primitive);
        
        if(primitive_info_[primitive].vertex_gather_mode) {
            index_t n = nb_vertices_per_primitive[primitive];
            GLenum GL_primitive = primitive_info_[primitive].GL_primitive;
            n /= nb_vertices_per_GL_primitive(GL_primitive);
            for(index_t i=0; i<immediate_state_.buffer.size(); ++i) {
                if(immediate_state_.buffer[i].is_enabled()) {                
                    for(index_t j=0; j<n; ++j) {
                        glEnableVertexAttribArray(i*n+j);
                    }
                }
            }
        } else {
            for(index_t i=0; i<immediate_state_.buffer.size(); ++i) {
                if(immediate_state_.buffer[i].is_enabled()) {
                    glEnableVertexAttribArray(i);
                    
                    //   If not using VBO, specify vertex attrib pointer
                    // using old style API.
                    if(immediate_state_.buffer[i].VBO() == 0) {                
                        glVertexAttribPointer(
                            i,
                            GLint(immediate_state_.buffer[i].dimension()),
                            GL_FLOAT, GL_FALSE, 0,
                            immediate_state_.buffer[i].data()
                        );
                    }
                }
            }
        }

        prepare_to_draw(primitive);
        
        glUseProgram(primitive_info_[primitive].program[toggles_config_]);
    }

    void Context::end() {
        flush_immediate_buffers();
        glUseProgram(0);

        if(primitive_info_[immediate_state_.primitive()].vertex_gather_mode) {
            index_t n = nb_vertices_per_primitive[immediate_state_.primitive()];
            GLenum GL_primitive =
                primitive_info_[immediate_state_.primitive()].GL_primitive;
            n /= nb_vertices_per_GL_primitive(GL_primitive);
            for(index_t i=0; i<immediate_state_.buffer.size(); ++i) {
                for(index_t j=0; j<n; ++j) {
                    glDisableVertexAttribArray(i*n+j);
                }
            }
        } else {
            for(index_t i=0; i<immediate_state_.buffer.size(); ++i) {
                glDisableVertexAttribArray(i);
            }
        }

        if(
            uniform_state_.toggle[GLUP_PICKING].get() &&
            uniform_state_.picking_mode.get() == GLUP_PICK_PRIMITIVE
        ) {
            update_base_picking_id(0);
        }

        if(
            primitive_info_[immediate_state_.primitive()].vertex_gather_mode_VAO
            != 0 || immediate_state_.VAO() != 0
        ) {
            glBindVertexArray(0);
        }
    }

    void Context::draw_arrays(
        GLUPprimitive primitive, GLUPint first, GLUPsizei count
    ) {
        update_toggles_config();
        create_program_if_needed(primitive);
        if(!primitive_info_[primitive].implemented) {
            Logger::warn("GLUP")
                << "glupDrawArrays(): "
                << primitive_name[primitive]
                << " not implemented in this profile" << std::endl;
            return;
        }
        prepare_to_draw(primitive);
        update_uniform_buffer();
        glUseProgram(primitive_info_[primitive].program[toggles_config_]);
        glDrawArrays(primitive_info_[primitive].GL_primitive, first, count);
        glUseProgram(0);
    }

    void Context::draw_elements(
        GLUPprimitive primitive, GLUPsizei count,
        GLUPenum type, const GLUPvoid* indices
    ) {
        update_toggles_config();
        create_program_if_needed(primitive);
        if(!primitive_info_[primitive].implemented) {
            Logger::warn("GLUP")
                << "glupDrawElements(): "
                << primitive_name[primitive]
                << " not implemented in this profile" << std::endl;
            return;
        }
        prepare_to_draw(primitive);
        update_uniform_buffer();
        glUseProgram(primitive_info_[primitive].program[toggles_config_]);
        glDrawElements(
            primitive_info_[primitive].GL_primitive, count, type, indices
        );
        glUseProgram(0);
    }

    void Context::prepare_to_draw(GLUPprimitive primitive) {
        if(primitive_info_[primitive].GL_primitive == GL_PATCHES) {
            // We generate an isoline for each patch, with the
            // minimum tesselation level. This generates two
            // vertices (we discard one of them in the geometry
            // shader).
            static float levels[4] = {1.0, 1.0, 0.0, 0.0};
            glPatchParameterfv(GL_PATCH_DEFAULT_OUTER_LEVEL, levels);
                
            // Specify number of vertices for GL_PATCH.
            glPatchParameteri(
                GL_PATCH_VERTICES, GLint(nb_vertices_per_primitive[primitive])
            );
        }
    }
    
    void Context::update_matrices() {
        
        if(!matrices_dirty_) {
            return;
        }

        GLfloat* modelview = matrix_stack_[GLUP_MODELVIEW_MATRIX].top();
        GLfloat* projection = matrix_stack_[GLUP_PROJECTION_MATRIX].top();

        copy_vector(
            uniform_state_.modelview_matrix.get_pointer(), modelview, 16
        );
        mult_matrices(
            modelview, projection,
            uniform_state_.modelviewprojection_matrix.get_pointer()
        );
        GLfloat modelview_invert[16];
        GLboolean OK = invert_matrix(modelview, modelview_invert);
        if(!OK) {
            Logger::warn("GLUP") << "Singular ModelView matrix"
                                 << std::endl;
            show_matrix(modelview);
        }
        GLfloat* normal_matrix = uniform_state_.normal_matrix.get_pointer();
        //   Copy the upper leftmost 3x3 part of the transpose of
        // modelview_invert_matrix to normal_matrix
        for(index_t i=0; i<3; ++i) {
            for(index_t j=0; j<3; ++j) {
                // Yes, it is i*4
                //   (with a '4' though it is a mat3 with a '3'),
                // mat3 rows are padded-aligned in UBOs !
                // TODO: query padding in UBO using introspection
                // (is it possible ? does not seem to work, so
                //  is padding with 4 universal ??)
                normal_matrix[i*4+j] = modelview_invert[j*4+i];
            }
        }

        copy_vector(
            uniform_state_.texture_matrix.get_pointer(),
            matrix_stack_[GLUP_TEXTURE_MATRIX].top(),
            16
        );
        
        matrices_dirty_ = false;
    }

    void Context::update_lighting() {
        if(!lighting_dirty_) {
            return;
        }

        GLfloat* light_vector =
            uniform_state_.light_vector.get_pointer();
        
        GLfloat* light_half_vector =
            uniform_state_.light_half_vector.get_pointer();

        // Normalize light vector
        
        normalize_vector(light_vector);

        // Compute half vector
        
        copy_vector(light_half_vector, light_vector, 3);
        light_half_vector[2] += 1.0f;
        normalize_vector(light_half_vector);        
        
        lighting_dirty_ = false;
    }

    void Context::update_base_picking_id(GLint new_value) {
        Memory::pointer address = uniform_state_.base_picking_id.address();
        index_t offset = index_t(address - uniform_buffer_data_);
        uniform_state_.base_picking_id.set(new_value);
        glBindBuffer(GL_UNIFORM_BUFFER, uniform_buffer_);
        glBufferSubData(
            GL_UNIFORM_BUFFER,
            offset,
            sizeof(int),
            address
        );
        glBindBuffer(GL_UNIFORM_BUFFER, 0);
    }
    
    void Context::update_uniform_buffer() {
        if(!uniform_buffer_dirty_) {
            return;
        }
        update_matrices();
        update_lighting();
        glBindBuffer(GL_UNIFORM_BUFFER, uniform_buffer_);
        glBufferSubData(
            GL_UNIFORM_BUFFER,
            0,
            uniform_buffer_size_,
            uniform_buffer_data_
        );
        glBindBuffer(GL_UNIFORM_BUFFER, 0);
        uniform_buffer_dirty_ = false;
    }

    void Context::flush_immediate_buffers() {
        if(immediate_state_.nb_vertices() == 0) {
            return;
        }
        
        // Sends the data from the buffers to OpenGL if VBO are used.
        stream_immediate_buffers();

        GLsizei nb_vertices = GLsizei(immediate_state_.nb_vertices());

        if(primitive_info_[immediate_state_.primitive()].vertex_gather_mode) {
            nb_vertices /= GLsizei(
                nb_vertices_per_primitive[immediate_state_.primitive()]
            );
            nb_vertices *= GLsizei(
                nb_vertices_per_GL_primitive(
                    primitive_info_[immediate_state_.primitive()].GL_primitive
                )
            );
        }
        
        glDrawArrays(
            primitive_info_[immediate_state_.primitive()].GL_primitive,
            0,
            nb_vertices
        );

        // Picking mode uses GLSL primitive_id variable, and add
        // GLUP state base_picking_id to it. This code updates
        // GLUP state base_picking_id and adds the number of drawn
        // primitives to it, so that the next batch will start with
        // the correct base_picking_id
        if(
            uniform_state_.toggle[GLUP_PICKING].get() &&
            uniform_state_.picking_mode.get() == GLUP_PICK_PRIMITIVE
        ) {
            update_base_picking_id(
                uniform_state_.base_picking_id.get() +
                GLint(immediate_state_.nb_primitives())
            );
        }

        immediate_state_.reset();
    }

    /***********************************************************************/

    void Context::set_primitive_info(
        GLUPprimitive glup_primitive, GLenum gl_primitive, GLuint program
    ) {
        primitive_info_[glup_primitive].implemented = true;
        primitive_info_[glup_primitive].GL_primitive = gl_primitive;
        primitive_info_[glup_primitive].program[toggles_config_] = program;
        primitive_info_[glup_primitive].program_initialized[toggles_config_] =
            true;
        
        bind_uniform_state(program);

        glBindAttribLocation(program, 0, "vertex_in");
        glBindAttribLocation(program, 1, "color_in");
        glBindAttribLocation(program, 2, "tex_coord_in");


        // Bind all textures to texture unit zero.
        GLSL::set_program_uniform_by_name(program, "texture1D", 0);
        GLSL::set_program_uniform_by_name(program, "texture2D", 0);
        GLSL::set_program_uniform_by_name(program, "texture3D", 0);        
    }

    void Context::set_primitive_info_vertex_gather_mode(
        GLUPprimitive glup_primitive, GLenum GL_primitive, GLuint program
    ) {
        set_primitive_info(glup_primitive, GL_primitive, program);
        index_t n = nb_vertices_per_primitive[glup_primitive];
        n /= nb_vertices_per_GL_primitive(GL_primitive);
        
        // Attribute location are bound here, programatically,
        // since their index depends on the number of vertices,
        // and GLSL does not like to have that in the declaration
        // (when saying layout(binding = nb_vertices), the GLSL
        // compiler does not "see" that nb_vertices is a constant).
        glBindAttribLocation(program, 0, "vertex_in");
        glBindAttribLocation(program, n, "color_in");
        glBindAttribLocation(program, 2*n, "tex_coord_in");

        primitive_info_[glup_primitive].vertex_gather_mode = true;
    
        //   We need a special VAO: memory layout is different since
        // we use a single vertex with an array of n attributes...
        //   Note: since the function can be called several times (one
        // per toggles configuration), make sute the VAO does not already
        // exists.
        if(primitive_info_[glup_primitive].vertex_gather_mode_VAO == 0) {
            glGenVertexArrays(
                1, &(primitive_info_[glup_primitive].vertex_gather_mode_VAO)
            );
            glBindVertexArray(
                primitive_info_[glup_primitive].vertex_gather_mode_VAO
            );

            for(index_t i=0; i<immediate_state_.buffer.size(); ++i) {
                glBindBuffer(GL_ARRAY_BUFFER,immediate_state_.buffer[i].VBO());
                for(index_t j=0; j<n; ++j) {
                    glVertexAttribPointer(
                        i*n+j,
                        4,
                        GL_FLOAT,
                        GL_FALSE,
                        GLsizei(sizeof(GL_FLOAT)*4*n),        // stride
                        (const GLvoid*)(sizeof(GL_FLOAT)*4*j) // pointer   
                    );
                }
            }

            glBindBuffer(GL_ARRAY_BUFFER,0);        
            glBindVertexArray(0);
        }
    }

    void Context::setup_primitives() {
        primitive_info_.resize(GLUP_NB_PRIMITIVES);

        if(!precompile_shaders_) {
            return;
        }
        
        GEO::Logger::out("GLUP compile") << "Optimizing shaders..."
                                         << std::endl;

        {
            GEO::ProgressTask progress(
                "GLUP compile", PrimitiveInfo::nb_toggles_configs
            );
            for(index_t i=0; i<PrimitiveInfo::nb_toggles_configs; ++i) {
                setup_shaders_source_for_toggles_config(i);
                toggles_config_ = i;
                setup_GLUP_POINTS();
                setup_GLUP_LINES();
                setup_GLUP_TRIANGLES();
                setup_GLUP_QUADS();
                setup_GLUP_TETRAHEDRA();
                setup_GLUP_HEXAHEDRA();
                setup_GLUP_PRISMS();
                setup_GLUP_PYRAMIDS();
                progress.next();
            }
        }
        
        GEO::Logger::out("GLUP compile") << "Shaders ready."
                                         << std::endl;
    }

    void Context::setup_shaders_source_for_toggles(
        GLUPbitfield toggles_state,
        GLUPbitfield toggles_undetermined
    ) {
        //   I think that the GLSL compiler has a bug whenever
        // the first vertex attribute (color) is not used,
        // whereas the second one (tex coord) is used (by 'used',
        // as mean 'determined statically to be potentially used
        // by the shader').... (or it's me who has a bug).
        //  Therefore, whenever texture coordinates are used,
        // and colors are not used, we don't tell the GLSL compiler
        // and let the shader figure out dynamically from the state.
        if(
            ((toggles_state & (1 << GLUP_TEXTURING)) != 0) &&
            ((toggles_state & (1 << GLUP_VERTEX_COLORS)) == 0)            
        ) {
            toggles_undetermined |= (1 << GLUP_VERTEX_COLORS);
        }
        
        toggles_shader_source_ = "";
        for(index_t i=0; i<uniform_state_.toggle.size(); ++i) {
            toggles_shader_source_ +=
                ("bool " + uniform_state_.toggle[i].name() + "() {\n");
            if(toggles_undetermined & (1 << i)) {
                toggles_shader_source_ +=
                    "   return GLUP." +
                    uniform_state_.toggle[i].name() +
                    ";\n";
            } else {
                if(toggles_state & (1 << i)) {
                    toggles_shader_source_ +=
                        "   return true; \n";
                } else {
                    toggles_shader_source_ +=
                        "   return false; \n";
                }
            }
            toggles_shader_source_ += "}\n";
        }
    }

    void Context::update_toggles_config() {
        if(uniform_state_.toggle[GLUP_PICKING].get()) {
            toggles_config_ = (1u << GLUP_PICKING);
        } else {
            toggles_config_ = 0;
            for(index_t i=0; i<uniform_state_.toggle.size(); ++i) {
                if(uniform_state_.toggle[i].get()) {
                    toggles_config_ |= (1u << i);
                }
            }
        }
    }

    void Context::create_program_if_needed(GLUPprimitive primitive) {
        if(!primitive_info_[primitive].program_initialized[toggles_config_]) {
            setup_shaders_source_for_toggles_config(toggles_config_);
            switch(primitive) {
            case GLUP_POINTS:
                setup_GLUP_POINTS();
                break;
            case GLUP_LINES:
                setup_GLUP_LINES();
                break;
            case GLUP_TRIANGLES:
                setup_GLUP_TRIANGLES();
                break;
            case GLUP_QUADS:
                setup_GLUP_QUADS();
                break;
            case GLUP_TETRAHEDRA:
                setup_GLUP_TETRAHEDRA();
                break;
            case GLUP_HEXAHEDRA:
                setup_GLUP_HEXAHEDRA();
                break;
            case GLUP_PRISMS:
                setup_GLUP_PRISMS();
                break;
            case GLUP_PYRAMIDS:
                setup_GLUP_PYRAMIDS();
                break;
            default:
                break;
            }
        }
    }
    
    void Context::setup_GLUP_POINTS() {
    }

    void Context::setup_GLUP_LINES() {
    }

    void Context::setup_GLUP_TRIANGLES() {
    }

    void Context::setup_GLUP_QUADS() {
    }

    void Context::setup_GLUP_TETRAHEDRA() {
    }

    void Context::setup_GLUP_HEXAHEDRA() {
    }

    void Context::setup_GLUP_PRISMS() {
    }

    void Context::setup_GLUP_PYRAMIDS() {
    }
    
    /***********************************************************************/
    /***********************************************************************/
    /**** GLUP implementation using GLSL 1.5                             ***/
    /***********************************************************************/
    /***********************************************************************/
    
    const char* GLUP150_shader_source_header =
        "#version 150 core                          \n"
        ;

    const char* GLUP150_vshader_in_out_declaration =
        "in vec4 vertex_in;                         \n"
        "in vec4 color_in;                          \n"
        "in vec4 tex_coord_in;                      \n"
        "out VertexData {                           \n"
        "   vec4 color;                             \n"
        "   vec4 tex_coord;                         \n"
        "} VertexOut;                               \n";


    // The geometry shader gets the input geometry and
    // attributes (color, tex_coords) through the functions
    // vertex_in(), color_in() (and tex_coord_in() in the
    // future) and the predicate prim_is_discarded().
    // There can be other ways of implementing these
    // functions, see vertex_gather_mode and gather
    // tesselation evaluation shader.
    // The function prim_is_discarded() returns true if
    // the primitive should be ignored by the geometry
    // shader (this happens with the gather tesselation
    // evaluation shader).
    
    const char* GLUP150_gshader_in_out_declaration =
        "in VertexData {                            \n"
        "    vec4 color;                            \n"
        "    vec4 tex_coord;                        \n"
        "} VertexIn[];                              \n"
        "                                           \n"
        "vec4 vertex_in(int i) {                    \n"
        "    return gl_in[i].gl_Position;           \n"
        "}                                          \n"
        "                                           \n"
        "vec4 color_in(int i) {                     \n"
        "    return VertexIn[i].color;              \n"
        "}                                          \n"
        "                                           \n"
        "vec4 tex_coord_in(int i) {                 \n"
        "    return VertexIn[i].tex_coord;          \n"
        "}                                          \n"
        "                                           \n"        
        "out FragmentData {                         \n"
        "   vec4 color;                             \n"
        "   vec4 tex_coord;                         \n"
        "   flat float diffuse;                     \n"
        "   flat float specular;                    \n"
        "   vec2 bary;                              \n"
        "   flat vec3 edge_mask;                    \n"
        "} FragmentOut;                             \n"
        "                                           \n"
        "bool prim_is_discarded() {                 \n"
        "   return false;                           \n"
        "}                                          \n";

   
    const char* GLUP150_fshader_in_out_declaration =
        "out vec4 frag_color ;                      \n"
        "                                           \n"
        "in FragmentData {                          \n"
        "   vec4 color;                             \n"
        "   vec4 tex_coord;                         \n"
        "   flat float diffuse;                     \n"
        "   flat float specular;                    \n"
        "   vec2 bary;                              \n"
        "   flat vec3 edge_mask;                    \n"
        "} FragmentIn;                              \n"
        "                                           \n"
        "uniform sampler1D texture1D;               \n"
        "uniform sampler2D texture2D;               \n"
        "uniform sampler3D texture3D;               \n";

    /**
     * \brief Declaration of input/output for simple fragment shaders
     *   (for points and lines).
     * \details This one is used when the fragment shader is directly
     *  plugged to the vertex shader (i.e. without geometry shader).
     */
    const char* GLUP150_simple_fshader_in_out_declaration =
        "out vec4 frag_color ;                      \n"
        "                                           \n"
        "in VertexData {                            \n"
        "   vec4 color;                             \n"
        "   vec4 tex_coord;                         \n"
        "} FragmentIn;                              \n"
        "                                           \n"
        "uniform sampler1D texture1D;               \n"
        "uniform sampler2D texture2D;               \n"
        "uniform sampler3D texture3D;               \n";

    // Note There is packUnorm4x8() and unpackUnorm4x8() that does what
    // we want, but it is only supported in GLSL 4.1...
    const char* GLUP150_fshader_utils =
        "vec4 int_to_vec4(int x) {                                           \n"
        "  return vec4(                                                      \n"
        "     float(x         & 255)/255.0,                                  \n"
        "     float((x >>  8) & 255)/255.0,                                  \n"
        "     float((x >> 16) & 255)/255.0,                                  \n"
        "     float((x >> 24) & 255)/255.0                                   \n"
        "  );                                                                \n"
        "}                                                                   \n"
        "                                                                    \n"
        "void output_picking_id() {                                          \n"
        "   if(GLUP.picking_mode == 1) {                                     \n"
        "      frag_color = int_to_vec4(gl_PrimitiveID+GLUP.base_picking_id);\n"
        "   } else {                                                         \n"
        "      frag_color = int_to_vec4(GLUP.picking_id);                    \n"
        "   }                                                                \n"
        "}                                                                   \n"
        "                                                                    \n"
        "void get_color() {                                                  \n"
        "   if(vertex_colors_enabled()) {                                    \n"
        "        frag_color = FragmentIn.color;                              \n"
        "   } else {                                                         \n"
        "        frag_color = gl_FrontFacing ?                               \n"
        "                       GLUP.front_color : GLUP.back_color;          \n"
        "   }                                                                \n"
        "   if(texturing_enabled()) {                                        \n"
        "       vec4 tex_color;                                              \n"
        "       switch(GLUP.texture_type) {                                  \n"
        "          case GLUP_TEXTURE_1D:                                     \n"
        "            tex_color = texture(texture1D, FragmentIn.tex_coord.x); \n"
        "            break;                                                  \n"
        "          case GLUP_TEXTURE_2D:                                     \n"
        "           tex_color = texture(texture2D, FragmentIn.tex_coord.xy); \n"
        "           break;                                                   \n"
        "          case GLUP_TEXTURE_3D:                                     \n"
        "           tex_color = texture(texture3D, FragmentIn.tex_coord.xyz);\n"
        "           break;                                                   \n"
        "       }                                                            \n"
        "       switch(GLUP.texture_mode) {                                  \n"
        "          case GLUP_TEXTURE_REPLACE:                                \n"
        "             frag_color = tex_color;                                \n"
        "             break;                                                 \n"
        "          case GLUP_TEXTURE_MODULATE:                               \n"
        "             frag_color *= tex_color;                               \n"
        "             break;                                                 \n"
        "          case GLUP_TEXTURE_ADD:                                    \n"
        "             frag_color += tex_color;                               \n"
        "             break;                                                 \n"
        "       }                                                            \n"
        "   }                                                                \n"
        "}                                                                   \n"
        "                                                                    \n"
        "void output_lighting(float diff, float spec) {                      \n"
        "   float s = gl_FrontFacing ? 1.0 : -1.0 ;                          \n"
        "   float sdiffuse = s * diff ;                                      \n"
        "   if(sdiffuse > 0.0) {                                             \n"
        "       vec3 vspec = spec*vec3(1.0,1.0,1.0);                         \n"
        "       frag_color = sdiffuse*frag_color + vec4(vspec,1.0);          \n"
        "       frag_color.rgb += vec3(0.2, 0.2, 0.2);                       \n"
        "   } else {                                                         \n"
        "       frag_color = vec4(0.2, 0.2, 0.2, 1.0);                       \n"
        "   }                                                                \n"
        "}                                                                   \n"
        ;
    
#define GLUP150_std                       \
        GLUP150_shader_source_header,     \
        GLUP_uniform_state_source,        \
        toggles_declaration()

    // GLUP_POINTS ********************************************************

    const char* GLUP150_points_and_lines_vshader_source =
        "void main() {                                              \n"
        "    if(clipping_enabled()) {                               \n"
        "       gl_ClipDistance[0] = dot(                           \n"
        "          GLUP.modelview_matrix*vertex_in,                 \n"
        "          GLUP.clip_plane                                  \n"
        "       );                                                  \n"
        "   } else {                                                \n"
        "      gl_ClipDistance[0] = 0.0;                            \n"
        "   }                                                       \n"
        "   if(vertex_colors_enabled()) {                           \n"
        "      VertexOut.color = color_in;                          \n"
        "   }                                                       \n"
        "   if(texturing_enabled()) {                               \n"
        "      VertexOut.tex_coord =                                \n"
        "                       GLUP.texture_matrix * tex_coord_in; \n"
        "   }                                                       \n"
        "   gl_Position =                                           \n"
        "                GLUP.modelviewprojection_matrix*vertex_in; \n"
        "}                                                          \n";
    
    // Note: depth update is not correct, it should be something like:
    // (to be checked...)
    // gl_FragDepth = gl_FragCoord.z +
    //   (pt_size*0.0001)/3.0 * gl_ProjectionMatrix[2].z * sqrt(1.0 - r2);

    const char* GLUP150_points_fshader_source =
        "#extension GL_ARB_conservative_depth : enable                      \n"
        "layout (depth_less) out float gl_FragDepth;                        \n"
        "                                                                   \n"
        "void main() {                                                      \n"
        "   vec2 V = 2.0*(gl_PointCoord - vec2(0.5, 0.5));                  \n"
        "   float one_minus_r2 = 1.0 - dot(V,V);                            \n"
        "   if(one_minus_r2 < 0.0) {                                        \n"
        "      discard;                                                     \n"
        "   }                                                               \n"
        "   vec3 W = vec3(V.x, -V.y, sqrt(one_minus_r2));                   \n"
        "   gl_FragDepth = gl_FragCoord.z - 0.001 * W.z;                    \n"
        "   if(picking_enabled()) {                                         \n"
        "        output_picking_id();                                       \n"
        "   } else {                                                        \n"
        "        get_color();                                               \n"
        "        if(lighting_enabled()) {                                   \n"
        "            float diff = dot(W,GLUP.light_vector);                 \n"
        "            float spec = dot(W,GLUP.light_half_vector);            \n"
        "            spec = pow(spec,30.0);                                 \n"
        "            output_lighting(diff,spec);                            \n"
        "        }                                                          \n"
        "    }                                                              \n"
        "}                                                                  \n";

    void Context_GLSL150::setup_GLUP_POINTS() {
        
        glEnable(GL_POINT_SPRITE);
        glTexEnvi(GL_POINT_SPRITE, GL_COORD_REPLACE, GL_TRUE);
        glEnable(GL_POINT_SMOOTH);
        
        GLuint vshader = GLSL::compile_shader(
            GL_VERTEX_SHADER,
            GLUP150_std,
            GLUP150_vshader_in_out_declaration,
            GLUP150_points_and_lines_vshader_source,
            0
        );

        GLuint fshader = GLSL::compile_shader(
            GL_FRAGMENT_SHADER,
            GLUP150_std,
            GLUP150_simple_fshader_in_out_declaration,
            GLUP150_fshader_utils,
            GLUP150_points_fshader_source,
            0
        );

        GLuint program = GLSL::create_program_from_shaders(
            vshader, fshader, 0 
        );
        
        set_primitive_info(GLUP_POINTS, GL_POINTS, program);

        glDeleteShader(vshader);
        glDeleteShader(fshader);
    }

    // GLUP_LINES *******************************************************

    const char* GLUP150_lines_fshader_source =
        "void main() {                                              \n"
        "   if(picking_enabled()) {                                 \n"
        "      output_picking_id();                                 \n"
        "   } else {                                                \n"
        "      get_color();                                         \n"
        "   }                                                       \n"
        "}                                                          \n";
    
    void Context_GLSL150::setup_GLUP_LINES() {
        
        GLuint vshader = GLSL::compile_shader(
            GL_VERTEX_SHADER,
            GLUP150_std,
            GLUP150_vshader_in_out_declaration,
            GLUP150_points_and_lines_vshader_source,
            0
        );
        
        GLuint fshader = GLSL::compile_shader(
            GL_FRAGMENT_SHADER,
            GLUP150_std,
            GLUP150_simple_fshader_in_out_declaration,
            GLUP150_fshader_utils,                        
            GLUP150_lines_fshader_source,
            0
        );

        GLuint program = GLSL::create_program_from_shaders(
            vshader, fshader, 0 
        );
        
        set_primitive_info(GLUP_LINES, GL_LINES, program);

        glDeleteShader(vshader);
        glDeleteShader(fshader);
    }

    // GLUP_TRIANGLES *******************************************************

    /**
     * \brief The fragment shader for polygons if GLSL version is 1.50.
     */
    const char* GLUP150_triangle_fshader_source = 
        "float edge_factor() {                                              \n"
        "    vec3 bary3 = vec3(                                             \n"
        "       FragmentIn.bary.x,                                          \n"
        "       FragmentIn.bary.y,                                          \n"
        "       1.0-FragmentIn.bary.x-FragmentIn.bary.y                     \n"
        "    ) ;                                                            \n"
        "    vec3 d = fwidth(bary3);                                        \n"
        "    vec3 a3 = smoothstep(                                          \n"
        "                  vec3(0.0,0.0,0.0), d*GLUP.mesh_width, bary3      \n"
        "    );                                                             \n"
        "    a3 = vec3(1.0, 1.0, 1.0)                                       \n"
        "           - FragmentIn.edge_mask + FragmentIn.edge_mask*a3;       \n"
        "    return min(min(a3.x, a3.y), a3.z);                             \n"
        "}                                                                  \n"
        "                                                                   \n"
        "void main() {                                                      \n"
        "    if(picking_enabled()) {                                        \n"
        "        output_picking_id();                                       \n"
        "    } else {                                                       \n"
        "        get_color();                                               \n"
        "        if(lighting_enabled()) {                                   \n"
        "          output_lighting(FragmentIn.diffuse, FragmentIn.specular);\n"
        "        }                                                          \n"
        "        if(draw_mesh_enabled()) {                                  \n"
        "            frag_color = mix(                                      \n"
        "                  GLUP.mesh_color,frag_color,edge_factor()         \n"
        "            );                                                     \n"
        "        }                                                          \n"
        "    }                                                              \n"
        "}                                                                  \n";
    
    /** 
     * \brief Some utility functions for the geometry shaders.
     * \details Provides functions for clipping, projection, and
     *  for generating shaded polygons.
     *  - flat_shaded_triangle(p1,p2,p3,pp1,pp2,pp3,do_clip) where
     *   (p1,p2,p3) are the coordinates in world space, (pp1,pp2,pp3) the
     *   transformed coordinates in clip space and do_clip specifies whether
     *   the triangle should be clipped.
     *  - flat_shaded_quad(p1,p2,p3,p4,pp1,pp2,pp3,pp4,do_clip,edges)
     */
    const char* GLUP150_gshader_utils_source =
        "out float gl_ClipDistance[];                                       \n"
        "vec4 projected[nb_vertices];                                       \n"
        "                                                                   \n"
        "void project_vertices() {                                          \n"
        "   for(int i=0; i<nb_vertices; ++i) {                              \n"
        "     projected[i] = GLUP.modelviewprojection_matrix*vertex_in(i);  \n"
        "   }                                                               \n"
        "   if(GLUP.cells_shrink != 0.0) {                                  \n"
        "       vec4 g = vec4(0.0, 0.0, 0.0, 0.0);                          \n"
        "       for(int i=0; i<nb_vertices; ++i) {                          \n"
        "            g += projected[i];                                     \n"
        "       }                                                           \n"
        "       g /= float(nb_vertices);                                    \n"
        "       float s = GLUP.cells_shrink;                                \n"
        "       for(int i=0; i<nb_vertices; ++i) {                          \n"
        "            projected[i] = mix(projected[i], g, s);                \n"
        "       }                                                           \n"
        "   }                                                               \n"
        "}                                                                  \n"
        "                                                                   \n"
        "float clip(vec4 V, bool do_clip) {                                 \n"
        "  return do_clip?dot(GLUP.modelview_matrix*V,GLUP.clip_plane):1.0; \n"
        "}                                                                  \n"
        "                                                                   \n"
        "bool cell_is_clipped() {                                           \n"
        "  if(prim_is_discarded()) {                                        \n"
        "     return true;                                                  \n"
        "  }                                                                \n"
        "  if(!clipping_enabled()) {                                        \n"
        "     return false;                                                 \n"
        "  }                                                                \n"
        "  int count = 0;                                                   \n"
        "  for(int i=0; i<nb_vertices; ++i) {                               \n"
        "       count += int(clip(vertex_in(i),true) >= 0.0);               \n"
        "  }                                                                \n"
        "  if(GLUP.clipping_mode==GLUP_CLIP_WHOLE_CELLS && count == 0) {    \n"
        "    return true;                                                   \n"
        "  }                                                                \n"
        "  if(                                                              \n"
        "      GLUP.clipping_mode==GLUP_CLIP_STRADDLING_CELLS &&            \n"
        "      (count==0 || count==nb_vertices)                             \n"
        "  ) {                                                              \n"
        "    return true;                                                   \n"
        "  }                                                                \n"
        "  return false;                                                    \n"
        "}                                                                  \n"
        "                                                                   \n"
        " /* L is supposed to be normalized */                              \n"
        "float cosangle(vec3 N, vec3 L) {                                   \n"
        "   float s = inversesqrt(dot(N,N)) ;                               \n"
        "   return s*dot(N,L) ;                                             \n"
        "}                                                                  \n"
        "                                                                   \n"
        "void compute_lighting(vec3 N) {                                    \n"
        "     FragmentOut.diffuse = cosangle(N,GLUP.light_vector) ;         \n"
        "     FragmentOut.specular = abs(cosangle(N,GLUP.light_half_vector));\n"
        "     FragmentOut.specular = pow(FragmentOut.specular,30.0);        \n"
        "}                                                                  \n"
        "                                                                   \n"
        "void emit_vertex(int i, bool do_clip) {                            \n"
        "   gl_ClipDistance[0] = clip(vertex_in(i),do_clip);                \n"
        "   gl_Position = projected[i];                                     \n"
        "   if(vertex_colors_enabled()) {                                   \n"
        "      FragmentOut.color = color_in(i);                             \n"
        "   }                                                               \n"
        "   if(texturing_enabled()) {                                       \n"
        "      FragmentOut.tex_coord = tex_coord_in(i);                     \n"
        "   }                                                               \n"
        "   EmitVertex();                                                   \n"
        "}                                                                  \n"
        "                                                                   \n"
        "void flat_shaded_triangle(int i1, int i2, int i3, bool do_clip) {  \n"
        "   if(lighting_enabled() && !picking_enabled()) {                  \n"
        "      vec4 p1 = vertex_in(i1);                                     \n"
        "      vec4 p2 = vertex_in(i2);                                     \n"
        "      vec4 p3 = vertex_in(i3);                                     \n"
        "      vec3 N = GLUP.normal_matrix *                                \n"
        "            cross((p2-p1).xyz,(p3-p1).xyz) ;                       \n"
        "      compute_lighting(N);                                         \n"
        "   }                                                               \n"
        "   FragmentOut.edge_mask = vec3(1.0,1.0,1.0);                      \n"
        "   FragmentOut.bary = vec2(0.0,0.0);                               \n"
        "   emit_vertex(i1,do_clip);                                        \n"
        "   FragmentOut.bary = vec2(1.0,0.0);                               \n"
        "   emit_vertex(i2,do_clip);                                        \n"
        "   FragmentOut.bary = vec2(0.0,1.0);                               \n"
        "   emit_vertex(i3,do_clip);                                        \n"
        "   EndPrimitive();                                                 \n"
        "}                                                                  \n"
        "                                                                   \n"
        "void flat_shaded_quad(                                             \n"
        "     int i1, int i2, int i3, int i4,                               \n"
        "     bool do_clip                                                  \n"
        "  ) {                                                              \n"
        "   if(lighting_enabled() && !picking_enabled()) {                  \n"
        "      vec4 p1 = vertex_in(i1);                                     \n"
        "      vec4 p2 = vertex_in(i2);                                     \n"
        "      vec4 p3 = vertex_in(i3);                                     \n"
        "      vec4 p4 = vertex_in(i4);                                     \n"
        "      vec3 N = GLUP.normal_matrix * (                              \n" 
        "           cross((p2-p1).xyz,(p4-p1).xyz) -                        \n"
        "           cross((p4-p3).xyz,(p2-p3).xyz)                          \n"
        "      );                                                           \n"
        "      compute_lighting(N);                                         \n"
        "   }                                                               \n"
        "   FragmentOut.edge_mask = vec3(0.0, 1.0, 1.0);                    \n"
        "   FragmentOut.bary=vec2(1.0,0.0);                                 \n"
        "   emit_vertex(i1,do_clip);                                        \n"
        "   FragmentOut.bary=vec2(0.0,1.0);                                 \n"
        "   emit_vertex(i2,do_clip);                                        \n"
        "   FragmentOut.bary=vec2(0.0,0.0);                                 \n"
        "   emit_vertex(i3,do_clip);                                        \n"
        "   FragmentOut.edge_mask = vec3(0.0, 1.0, 1.0);                    \n"
        "   FragmentOut.bary=vec2(1.0,0.0);                                 \n"
        "   emit_vertex(i4,do_clip);                                        \n"
        "   EndPrimitive();                                                 \n"
        "}                                                                  \n"
        ;

    /**
     * \brief The pass-through vertex shader.
     * \details Used by points, quads, tets, prisms
     */
    const char* GLUP150_vshader_pass_through_source =
        " void main(void) {                                                 \n"
        "     if(vertex_colors_enabled()) {                                 \n"
        "        VertexOut.color = color_in;                                \n"
        "     }                                                             \n"
        "     if(texturing_enabled()) {                                     \n"
        "        VertexOut.tex_coord = GLUP.texture_matrix * tex_coord_in;  \n"
        "     }                                                             \n"
        "     gl_Position = vertex_in;                                      \n"
        " }                                                                 \n";


    /**
     * \brief The geometry shader for triangles.
     * \details Uses vshader_pass_through and gshader_utils.
     */
    const char* GLUP150_gshader_tri_source =
        "layout(triangles) in;                                              \n"
        "layout(triangle_strip, max_vertices = 3) out;                      \n"
        "                                                                   \n"
        "void main() {                                                      \n"
        "    gl_PrimitiveID = gl_PrimitiveIDIn;                             \n"
        "    project_vertices();                                            \n"
        "    flat_shaded_triangle(0,1,2,true);                              \n"
        "}                                                                  \n";
    
    void Context_GLSL150::setup_GLUP_TRIANGLES() {

        GLuint vshader = GLSL::compile_shader(
            GL_VERTEX_SHADER,
            GLUP150_std,
            GLUP150_vshader_in_out_declaration,
            GLUP150_vshader_pass_through_source,
            0
        );
        
        GLuint gshader = GLSL::compile_shader(
            GL_GEOMETRY_SHADER,
            GLUP150_std,            
            GLUP150_gshader_in_out_declaration,
            "const int nb_vertices = 3;",
            GLUP150_gshader_utils_source,
            GLUP150_gshader_tri_source,
            0
        );

        GLuint fshader = GLSL::compile_shader(
            GL_FRAGMENT_SHADER,
            GLUP150_std,
            GLUP150_fshader_in_out_declaration,
            GLUP150_fshader_utils,
            GLUP150_triangle_fshader_source,
            0
        );
                
        GLuint program = 
            GLSL::create_program_from_shaders(
                vshader, gshader, fshader, 0 
            );

        set_primitive_info(GLUP_TRIANGLES, GL_TRIANGLES, program);
        
        glDeleteShader(vshader);
        glDeleteShader(gshader);
        glDeleteShader(fshader);        
    }

    // GLUP_QUADS ***********************************************************

    /**
     * \brief The geometry shader for quads.
     * \details Uses vshader_pass_through and gshader_utils.
     */
    const char* GLUP150_gshader_quad_source =
        "layout(lines_adjacency) in;                                        \n"
        "layout(triangle_strip, max_vertices = 4) out;                      \n"
        "                                                                   \n"
        "void main() {                                                      \n"
        "    gl_PrimitiveID = gl_PrimitiveIDIn;                             \n"
        "    project_vertices();                                            \n"
        "    flat_shaded_quad(0,1,3,2,true);                                \n"
        "}                                                                  \n";
    
    void Context_GLSL150::setup_GLUP_QUADS() {

        GLuint vshader = GLSL::compile_shader(
            GL_VERTEX_SHADER,
            GLUP150_std,
            GLUP150_vshader_in_out_declaration,
            GLUP150_vshader_pass_through_source,
            0
        );

        
        GLuint gshader = GLSL::compile_shader(
            GL_GEOMETRY_SHADER,
            GLUP150_std,
            GLUP150_gshader_in_out_declaration,
            "const int nb_vertices = 4;",
            GLUP150_gshader_utils_source,
            GLUP150_gshader_quad_source,
            0
        );

        GLuint fshader = GLSL::compile_shader(
            GL_FRAGMENT_SHADER,
            GLUP150_std,
            GLUP150_fshader_in_out_declaration,
            GLUP150_fshader_utils,                        
            GLUP150_triangle_fshader_source,
            0
        );

        GLuint program = 
            GLSL::create_program_from_shaders(
                vshader, gshader, fshader, 0
            );
        
        set_primitive_info(GLUP_QUADS, GL_LINES_ADJACENCY, program);
        
        glDeleteShader(vshader);
        glDeleteShader(gshader);
        glDeleteShader(fshader);        
    }

    // GLUP_TETRAHEDRA ******************************************************

    /**
     * \brief The geometry shader for tetrahedra.
     * \details Uses v_shader_pass_through and gshader_utils.
     */
    const char* GLUP150_gshader_tet_source =
        "layout(lines_adjacency) in;                                        \n"
        "layout(triangle_strip, max_vertices = 12) out;                     \n"
        "                                                                   \n"
        "void main() {                                                      \n"
        "    if(cell_is_clipped()) {                                        \n"
        "        return;                                                    \n"
        "    }                                                              \n"
        "    gl_PrimitiveID = gl_PrimitiveIDIn;                             \n"
        "    project_vertices();                                            \n"
        "    flat_shaded_triangle(0,1,2,false);                             \n"
        "    flat_shaded_triangle(1,0,3,false);                             \n"
        "    flat_shaded_triangle(0,2,3,false);                             \n"
        "    flat_shaded_triangle(2,1,3,false);                             \n"
        "}                                                                  \n";
    
    
    void Context_GLSL150::setup_GLUP_TETRAHEDRA() {

        GLuint vshader = GLSL::compile_shader(
            GL_VERTEX_SHADER,
            GLUP150_std,            
            GLUP150_vshader_in_out_declaration,
            GLUP150_vshader_pass_through_source,
            0
        );

        GLuint gshader = GLSL::compile_shader(
            GL_GEOMETRY_SHADER,
            GLUP150_std,                        
            GLUP150_gshader_in_out_declaration,
            "const int nb_vertices = 4; ",
            GLUP150_gshader_utils_source,
            GLUP150_gshader_tet_source,
            0
        );

        GLuint fshader = GLSL::compile_shader(
            GL_FRAGMENT_SHADER,
            GLUP150_std,
            GLUP150_fshader_in_out_declaration,
            GLUP150_fshader_utils,                        
            GLUP150_triangle_fshader_source,
            0
        );

        GLuint program = 
            GLSL::create_program_from_shaders(
                vshader, gshader, fshader, 0 
            );
        
        set_primitive_info(GLUP_TETRAHEDRA, GL_LINES_ADJACENCY, program);
        
        glDeleteShader(vshader);
        glDeleteShader(gshader);
        glDeleteShader(fshader);        
    }


    // GLUP_PRISMS **********************************************************

    /**
     * \brief The geometry shader for prisms
     * \details Uses v_shader_pass_through and gshader_utils.
     */
    const char* GLUP150_gshader_prism_source =
        "layout(triangles_adjacency) in;                                    \n"
        "layout(triangle_strip, max_vertices = 18) out;                     \n"
        "                                                                   \n"
        "void main() {                                                      \n"
        "    if(cell_is_clipped()) {                                        \n"
        "        return;                                                    \n"
        "    }                                                              \n"
        "    gl_PrimitiveID = gl_PrimitiveIDIn;                             \n"
        "    project_vertices();                                            \n"
        "    flat_shaded_triangle(0,1,2,false);                             \n"
        "    flat_shaded_triangle(5,4,3,false);                             \n"
        "    flat_shaded_quad(0,3,1,4,false);                               \n"
        "    flat_shaded_quad(0,2,3,5,false);                               \n"
        "    flat_shaded_quad(1,4,2,5,false);                               \n" 
        "}                                                                  \n";

    
    void Context_GLSL150::setup_GLUP_PRISMS() {
        
        GLuint vshader = GLSL::compile_shader(
            GL_VERTEX_SHADER,
            GLUP150_std,            
            GLUP150_vshader_in_out_declaration,
            GLUP150_vshader_pass_through_source,
            0
        );

        GLuint gshader = GLSL::compile_shader(
            GL_GEOMETRY_SHADER,
            GLUP150_std,                        
            GLUP150_gshader_in_out_declaration,
            "const int nb_vertices = 6; ",
            GLUP150_gshader_utils_source,
            GLUP150_gshader_prism_source,
            0
        );

        GLuint fshader = GLSL::compile_shader(
            GL_FRAGMENT_SHADER,
            GLUP150_std,
            GLUP150_fshader_in_out_declaration,
            GLUP150_fshader_utils,                        
            GLUP150_triangle_fshader_source,
            0
        );

        GLuint program = 
            GLSL::create_program_from_shaders(
                vshader, gshader, fshader, 0 
            );

        set_primitive_info(GLUP_PRISMS, GL_TRIANGLES_ADJACENCY, program);
        
        glDeleteShader(vshader);
        glDeleteShader(gshader);
        glDeleteShader(fshader);        
    }

    // GLUP_HEXAHEDRA *******************************************************

    const char* GLUP150_vshader_gather_source =
        " const int nb_vertices_per_GL =                      \n"
        "                      nb_vertices / nb_vertices_GL;  \n"
        " in vec4 vertex_in[nb_vertices_per_GL];              \n"
        " in vec4 color_in[nb_vertices_per_GL];               \n"
        " in vec4 tex_coord_in[nb_vertices_per_GL];           \n"
        "                                                     \n"
        "out GVertexData {                                    \n"
        "    vec4 other_vertex[nb_vertices_per_GL-1];         \n"     
        "    vec4 color[nb_vertices_per_GL];                  \n"
        "    vec4 tex_coord[nb_vertices_per_GL];              \n"
        "} VertexOut;                                         \n"
        "                                                     \n"
        "void main() {                                        \n"
        "   for(int i=1; i<nb_vertices_per_GL; ++i) {         \n"
        "       VertexOut.other_vertex[i-1] = vertex_in[i];   \n"
        "   }                                                 \n"
        "   if(texturing_enabled()) {                         \n"
        "       for(int i=0; i<nb_vertices_per_GL; ++i) {     \n"
        "           VertexOut.tex_coord[i] =                  \n"
        "              GLUP.texture_matrix * tex_coord_in[i]; \n"
        "       }                                             \n"
        "   }                                                 \n"
        "   if(vertex_colors_enabled()) {                     \n"
        "       for(int i=0; i<nb_vertices_per_GL; ++i) {     \n"
        "           VertexOut.color[i] = color_in[i];         \n"
        "       }                                             \n"
        "   }                                                 \n"
        "   gl_Position = vertex_in[0];                       \n"
        "}                                                    \n"
        ;

    // To be used for primitives that have a number of vertices
    // that does not match existing OpenGL primitives. In that
    // case, the other vertices are gathered by the vertex shader,
    // and all vertices are merged into the attributes of a single
    // vertex. The geometry shader then expands this single vertex
    // into the primitive.
    const char* GLUP150_gshader_gather_in_out_declaration =
        " const int nb_vertices_per_GL =            \n"
        "            nb_vertices / nb_vertices_GL;  \n"
        "in GVertexData {                           \n"
        "   vec4 other_vertex[nb_vertices_per_GL-1];\n"     
        "   vec4 color[nb_vertices_per_GL];         \n"
        "   vec4 tex_coord[nb_vertices_per_GL];     \n"
        "} VertexIn[];                              \n"
        "                                           \n"
        "vec4 vertex_in(int i) {                    \n"
        "   int i0 = i / nb_vertices_per_GL;        \n"
        "   int i1 = i % nb_vertices_per_GL;        \n"
        "   return (i1==0) ? gl_in[i0].gl_Position :\n"
        "          VertexIn[i0].other_vertex[i1-1]; \n"
        "}                                          \n"
        "                                           \n"
        "vec4 color_in(int i) {                     \n"
        "    int i0 = i / nb_vertices_per_GL;       \n"
        "    int i1 = i % nb_vertices_per_GL;       \n"
        "    return VertexIn[i0].color[i1];         \n"
        "}                                          \n"
        "                                           \n"
        "vec4 tex_coord_in(int i) {                 \n"  
        "    int i0 = i / nb_vertices_per_GL;       \n"
        "    int i1 = i % nb_vertices_per_GL;       \n"
        "    return VertexIn[i0].tex_coord[i1];     \n"
        "}                                          \n"
        "                                           \n"        
        "out FragmentData {                         \n"
        "   vec4 color;                             \n"
        "   vec4 tex_coord;                         \n"
        "   flat float diffuse;                     \n"
        "   flat float specular;                    \n"
        "   vec2 bary;                              \n"
        "   flat vec3 edge_mask;                    \n"
        "} FragmentOut;                             \n"
        "                                           \n"
        "bool prim_is_discarded() {                 \n"
        "   return false;                           \n"
        "}                                          \n";

    /**
     * \brief The geometry shader for hexahedra
     * \details Uses v_shader_gather and gshader_utils.
     */
    const char* GLUP150_gshader_hex_source =
        "layout(triangle_strip, max_vertices = 24) out;                     \n"
        "                                                                   \n"
        "void main() {                                                      \n"
        "    if(cell_is_clipped()) {                                        \n"
        "        return;                                                    \n"
        "    }                                                              \n"
        "    gl_PrimitiveID = gl_PrimitiveIDIn;                             \n"
        "    project_vertices();                                            \n"
        "    flat_shaded_quad(0,2,4,6,false);                               \n"
        "    flat_shaded_quad(3,1,7,5,false);                               \n"
        "    flat_shaded_quad(1,0,5,4,false);                               \n"
        "    flat_shaded_quad(2,3,6,7,false);                               \n"
        "    flat_shaded_quad(1,3,0,2,false);                               \n"
        "    flat_shaded_quad(4,6,5,7,false);                               \n" 
        "}                                                                  \n";
    
    void Context_GLSL150::setup_GLUP_HEXAHEDRA() {

        GLuint vshader = GLSL::compile_shader(
            GL_VERTEX_SHADER,
            GLUP150_std,
            "const int nb_vertices = 8;",
            "const int nb_vertices_GL = 4;",
            GLUP150_vshader_gather_source,
            0
        );

        GLuint gshader = GLSL::compile_shader(
            GL_GEOMETRY_SHADER,
            GLUP150_std,
            "const int nb_vertices = 8;",
            "const int nb_vertices_GL = 4;",
            "layout(lines_adjacency) in;",
            GLUP150_gshader_gather_in_out_declaration,
            GLUP150_gshader_utils_source,
            GLUP150_gshader_hex_source,
            0
        );

        GLuint fshader = GLSL::compile_shader(
            GL_FRAGMENT_SHADER,
            GLUP150_std,
            GLUP150_fshader_in_out_declaration,
            GLUP150_fshader_utils,                        
            GLUP150_triangle_fshader_source,
            0
        );

        GLuint program = 
            GLSL::create_program_from_shaders(
                vshader, gshader, fshader, 0 
            );

        set_primitive_info_vertex_gather_mode(
            GLUP_HEXAHEDRA, GL_LINES_ADJACENCY, program
        );
        
        glDeleteShader(vshader);
        glDeleteShader(gshader);
        glDeleteShader(fshader);                
    }
    
    // GLUP_PYRAMIDS ********************************************************

    /**
     * \brief The geometry shader for pyramids
     * \details Uses v_shader_gather and gshader_utils.
     */
    const char* GLUP150_gshader_pyramid_source =
        "layout(points) in;                                                 \n"
        "layout(triangle_strip, max_vertices = 28) out;                     \n"
        "                                                                   \n"
        "void main() {                                                      \n"
        "    if(cell_is_clipped()) {                                        \n"
        "        return;                                                    \n"
        "    }                                                              \n"
        "    gl_PrimitiveID = gl_PrimitiveIDIn;                             \n"
        "    project_vertices();                                            \n"
        "    flat_shaded_quad(0,1,3,2,false);                               \n"
        "    flat_shaded_triangle(0,4,1,false);                             \n"
        "    flat_shaded_triangle(0,3,4,false);                             \n"
        "    flat_shaded_triangle(2,4,3,false);                             \n"
        "    flat_shaded_triangle(2,1,4,false);                             \n"
        "}                                                                  \n";
    
    void Context_GLSL150::setup_GLUP_PYRAMIDS() {

        GLuint vshader = GLSL::compile_shader(
            GL_VERTEX_SHADER,
            GLUP150_std,
            "const int nb_vertices = 5;",
            "const int nb_vertices_GL = 1;",            
            GLUP150_vshader_gather_source,
            0
        );

        GLuint gshader = GLSL::compile_shader(
            GL_GEOMETRY_SHADER,
            GLUP150_std,
            "const int nb_vertices = 5;",
            "const int nb_vertices_GL = 1;",            
            GLUP150_gshader_gather_in_out_declaration,
            GLUP150_gshader_utils_source,
            GLUP150_gshader_pyramid_source,
            0
        );

        GLuint fshader = GLSL::compile_shader(
            GL_FRAGMENT_SHADER,
            GLUP150_std,
            GLUP150_fshader_in_out_declaration,
            GLUP150_fshader_utils,                        
            GLUP150_triangle_fshader_source,
            0
        );

        GLuint program = 
            GLSL::create_program_from_shaders(
                vshader, gshader, fshader, 0 
            );

        set_primitive_info_vertex_gather_mode(
            GLUP_PYRAMIDS, GL_POINTS, program
        );
        
        glDeleteShader(vshader);
        glDeleteShader(gshader);
        glDeleteShader(fshader);                
    }

    /***********************************************************************/
    /***********************************************************************/
    /**** GLUP implementation using GLSL 4.40                            ***/
    /***********************************************************************/
    /***********************************************************************/
    
    const char* GLUP440_shader_source_header =
        "#version 440 core                                     \n"
        ;

#define GLUP440_std                       \
        GLUP440_shader_source_header,     \
        GLUP_uniform_state_source,        \
        toggles_declaration()
    
     const char* GLUP440_teshader_gather_source =
        "layout(isolines, point_mode) in;                                   \n"
        "                                                                   \n"
        "in VertexData {                                                    \n"
        "   vec4 color;                                                     \n"
        "   vec4 tex_coord;                                                 \n"
        "} VertexIn[];                                                      \n"
        "                                                                   \n"
        "out GTVertexData {                                                 \n"
        "    vec4 vertex[nb_vertices];                                      \n"
        "    vec4 color[nb_vertices];                                       \n"
        "    vec4 tex_coord[nb_vertices];                                   \n"
        "    bool discard_me;                                               \n"
        "} VertexOut;                                                       \n"
        "                                                                   \n"
        "void main() {                                                      \n"
        "   VertexOut.discard_me = (gl_TessCoord.x > 0.5);                  \n"
        "   if(VertexOut.discard_me) {                                      \n"
        "        return;                                                    \n"
        "   }                                                               \n"
        "   for(int i=0; i<nb_vertices; ++i) {                              \n"
        "        VertexOut.vertex[i] = gl_in[i].gl_Position;                \n"
        "   }                                                               \n"
        "   if(vertex_colors_enabled()) {                                   \n"
        "     for(int i=0; i<nb_vertices; ++i) {                            \n"
        "        VertexOut.color[i] = VertexIn[i].color;                    \n"
        "     }                                                             \n"
        "   }                                                               \n"
        "   if(texturing_enabled()) {                                       \n"
        "     for(int i=0; i<nb_vertices; ++i) {                            \n"
        "        VertexOut.tex_coord[i] =                                   \n"
        "                      GLUP.texture_matrix * VertexIn[i].tex_coord; \n"
        "     }                                                             \n"
        "   }                                                               \n"
        "}                                                                  \n"
        ;

    // To be used for primitives that have a number of vertices
    // that does not match existing OpenGL primitives. In that
    // case, a tessellation shader fetches the vertices, using
    // GL_PATCHES primitive and generates two vertices (with a
    // lot of attributes). The second one needs to be discarded
    // (discard_me = true).
    
    const char* GLUP440_gshader_tegather_in_out_declaration =
        "in GTVertexData {                          \n"
        "    vec4 vertex[nb_vertices];              \n"     
        "    vec4 color[nb_vertices];               \n"
        "    vec4 tex_coord[nb_vertices];           \n"
        "    bool discard_me;                       \n"
        "} VertexIn[];                              \n"
        "                                           \n"
        "vec4 vertex_in(int i) {                    \n"
        "   return VertexIn[0].vertex[i];           \n"
        "}                                          \n"
        "                                           \n"
        "vec4 color_in(int i) {                     \n"
        "    return VertexIn[0].color[i];           \n"
        "}                                          \n"
        "                                           \n"
        "vec4 tex_coord_in(int i) {                 \n"
        "    return VertexIn[0].tex_coord[i];       \n"
        "}                                          \n"
        "                                           \n"        
        "out FragmentData {                         \n"
        "   vec4 color;                             \n"
        "   vec4 tex_coord;                         \n"        
        "   flat float diffuse;                     \n"
        "   flat float specular;                    \n"
        "   vec2 bary;                              \n"
        "   flat vec3 edge_mask;                    \n"
        "} FragmentOut;                             \n"
        "                                           \n"        
        "bool prim_is_discarded() {                 \n"
        "     return VertexIn[0].discard_me;        \n"
        "}                                          \n";
    
    void Context_GLSL440::setup_GLUP_HEXAHEDRA() {

        if(!GEO::CmdLine::get_arg_bool("gfx:GLSL_tesselation")) {
            Context_GLSL150::setup_GLUP_HEXAHEDRA();
            return;
        }

        GLuint vshader = GLSL::compile_shader(
            GL_VERTEX_SHADER,
            GLUP150_std,
            GLUP150_vshader_in_out_declaration,
            GLUP150_vshader_pass_through_source,
            0
        );

        GLuint teshader = GLSL::compile_shader(
            GL_TESS_EVALUATION_SHADER,
            GLUP440_std,
            "const int nb_vertices = 8;",
            GLUP440_teshader_gather_source,
            0
        );
        
        GLuint gshader = GLSL::compile_shader(
            GL_GEOMETRY_SHADER,
            GLUP440_std,
            "const int nb_vertices = 8;",
            "layout(points) in;",
            GLUP440_gshader_tegather_in_out_declaration,            
            GLUP150_gshader_utils_source,
            GLUP150_gshader_hex_source,
            0
        );

        GLuint fshader = GLSL::compile_shader(
            GL_FRAGMENT_SHADER,
            GLUP150_std,
            GLUP150_fshader_in_out_declaration,
            GLUP150_fshader_utils,                        
            GLUP150_triangle_fshader_source,
            0
        );

        GLuint program = 
            GLSL::create_program_from_shaders(
                vshader, teshader, gshader, fshader, 0 
            );

        set_primitive_info(
            GLUP_HEXAHEDRA, GL_PATCHES, program
        );
        
        glDeleteShader(teshader);
        glDeleteShader(vshader);        
        glDeleteShader(gshader);
        glDeleteShader(fshader);                
    }
    
    void Context_GLSL440::setup_GLUP_PYRAMIDS() {

        if(!GEO::CmdLine::get_arg_bool("gfx:GLSL_tesselation")) {
            Context_GLSL150::setup_GLUP_PYRAMIDS();
            return;
        }

        GLuint teshader = GLSL::compile_shader(
            GL_TESS_EVALUATION_SHADER,
            GLUP440_std,
            "const int nb_vertices = 5;",
            GLUP440_teshader_gather_source,
            0
        );

        GLuint vshader = GLSL::compile_shader(
            GL_VERTEX_SHADER,
            GLUP150_std,
            GLUP150_vshader_in_out_declaration,
            GLUP150_vshader_pass_through_source,
            0
        );
        
        GLuint gshader = GLSL::compile_shader(
            GL_GEOMETRY_SHADER,
            GLUP440_std,
            "const int nb_vertices = 5;",
            GLUP440_gshader_tegather_in_out_declaration,            
            GLUP150_gshader_utils_source,
            GLUP150_gshader_pyramid_source,
            0
        );

        GLuint fshader = GLSL::compile_shader(
            GL_FRAGMENT_SHADER,
            GLUP150_std,
            GLUP150_fshader_in_out_declaration,
            GLUP150_fshader_utils,                        
            GLUP150_triangle_fshader_source,
            0
        );

        GLuint program = 
            GLSL::create_program_from_shaders(
                vshader, teshader, gshader, fshader, 0 
            );

        set_primitive_info(
            GLUP_PYRAMIDS, GL_PATCHES, program
        );
        
        glDeleteShader(teshader);
        glDeleteShader(vshader);        
        glDeleteShader(gshader);
        glDeleteShader(fshader);                
    }

    /***********************************************************************/

    void Context_VanillaGL::setup() {
        uniform_buffer_dirty_=true;
        matrices_dirty_=true;        
        uniform_binding_point_=1;

        std::map<std::string, GLsizei> type_to_size;
        type_to_size["bool"] = sizeof(int);
        type_to_size["vec4"] = 4*sizeof(GLfloat);
        type_to_size["vec3"] = 3*sizeof(GLfloat);
        type_to_size["mat4"] = 4*4*sizeof(GLfloat);
        type_to_size["mat3"] = 4*3*sizeof(GLfloat); // yes, 4, there is padding
        type_to_size["float"] = sizeof(GLfloat);
        type_to_size["int"] = sizeof(GLint);

        // Parse uniform state declaration in order to "emulate" it...
        uniform_buffer_size_ = 0;
        std::istringstream input(GLUP_uniform_state_source);
        std::string line;
        while(std::getline(input,line)) {
            std::vector<std::string> words;
            GEO::String::split_string(line, ' ', words);
            if(
                (words.size() == 2 && words[1][words[1].length()-1] == ';') ||
                (words.size() == 3 && words[2] == ";")
            ) {
                std::string vartype = words[0];
                std::string varname = words[1];
                if(varname[varname.length()-1] == ';') {
                    varname = varname.substr(0, varname.length()-1);
                }
                if(type_to_size.find(vartype) != type_to_size.end()) {
                    variable_to_offset_[varname] = uniform_buffer_size_;
                    uniform_buffer_size_ += type_to_size[vartype];
                }
            }
        }
        
        uniform_buffer_data_ = new Memory::byte[uniform_buffer_size_];
        Memory::clear(uniform_buffer_data_, size_t(uniform_buffer_size_));

        setup_state_variables();
        setup_immediate_buffers();
        setup_primitives();
    }

    void Context_VanillaGL::flush_immediate_buffers() {
        switch(immediate_state_.primitive()) {
        case GLUP_POINTS:
            draw_immediate_buffer_GLUP_POINTS();
            break;
        case GLUP_LINES:
            draw_immediate_buffer_GLUP_LINES();
            break;
        default:
            Logger::warn("GLUP VanillaGL")
                << primitive_name[immediate_state_.primitive()]                
                << ":not implemented yet"
                << std::endl;
            break;
        }
        immediate_state_.reset();
    }

    void Context_VanillaGL::setup_immediate_buffers() {
    }

    void Context_VanillaGL::setup_primitives() {
        primitive_info_.resize(GLUP_NB_PRIMITIVES);
        primitive_info_[GLUP_POINTS].implemented = true;
        primitive_info_[GLUP_LINES].implemented = true;        
    }
    
    Memory::pointer Context_VanillaGL::get_state_variable_address(
        const char* name
    ) {
        geo_assert(variable_to_offset_.find(name) != variable_to_offset_.end());
        return uniform_buffer_data_ + variable_to_offset_[name];
    }
    
    bool Context_VanillaGL::primitive_supports_array_mode(
        GLUPprimitive prim
    ) const {
        geo_argused(prim);
        return false;
    }

    void Context_VanillaGL::draw_immediate_buffer_GLUP_POINTS() {
        glDisable(GL_LIGHTING);
        glEnable(GL_POINT_SMOOTH);
        glBegin(GL_POINTS);
        for(index_t i=0; i<immediate_state_.nb_vertices(); ++i) {
            if(immediate_state_.buffer[1].is_enabled()) {
                glColor4fv(immediate_state_.buffer[1].element_ptr(i)); 
            }
            if(immediate_state_.buffer[2].is_enabled()) {
                glTexCoord4fv(immediate_state_.buffer[2].element_ptr(i)); 
            }
            glVertex4fv(immediate_state_.buffer[0].element_ptr(i));
        }
        glEnd();
    }

    void Context_VanillaGL::draw_immediate_buffer_GLUP_LINES() {
        glDisable(GL_LIGHTING);
        glBegin(GL_LINES);
        for(index_t i=0; i<immediate_state_.nb_vertices(); ++i) {
            if(immediate_state_.buffer[1].is_enabled()) {
                glColor4fv(immediate_state_.buffer[1].element_ptr(i)); 
            }
            if(immediate_state_.buffer[2].is_enabled()) {
                glTexCoord4fv(immediate_state_.buffer[2].element_ptr(i)); 
            }
            glVertex4fv(immediate_state_.buffer[0].element_ptr(i));
        }
        glEnd();
    }
}

