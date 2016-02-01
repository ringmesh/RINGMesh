/*
 *  Copyright (c) 2012-2015, Bruno Levy
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

#include <geogram_gfx/basic/GL.h>
#include <geogram_gfx/basic/GLUP.h>
#include <geogram/basic/logger.h>

namespace {
    using namespace GEO;
 
    /** 
     * \brief Converts a mat4 into a matrix for OpenGL.
     * \param[in] m a const reference to the matrix
     * \return a const pointer to the converted matrix,
     *  stored as a static array of 16 doubles.
     * \note The matrix is transposed before being sent to OpenGL
     *  because Geogram uses the convention with column
     *  vectors and OpenGL the convention with row vectors
     *  to represent the transformed points.
     */
    const GLdouble* convert_matrix(
        const mat4& m
    ) {
        static double result[16] ;
        index_t k = 0 ;
        for(index_t i=0; i<4; i++) {
            for(index_t j=0; j<4; j++) {
                result[k] = m(i,j) ;
                k++ ;
            }
        }
        return result ;
    }
}

namespace GEO {
    
    void glLoadMatrix(const mat4& m) {
        glLoadMatrixd(convert_matrix(m));
    }

    void glMultMatrix(const mat4& m) {
        glMultMatrixd(convert_matrix(m));
    }

    void glupLoadMatrix(const mat4& m) {
        glupLoadMatrixd(convert_matrix(m));
    }

    void glupMultMatrix(const mat4& m) {
        glupMultMatrixd(convert_matrix(m));
    }
    

    GLint64 get_size_of_bound_buffer_object(GLenum target) {
        static bool init = false;
        static bool use_glGetBufferParameteri64v = true;

        // Note: there is a version of glGetBufferParameteriv that uses
        // 64 bit parameters. Since array data larger than 4Gb will be
        // common place, it is this version that should be used. However,
        // it is not supported by the Intel driver (therefore we fallback
        // to the standard 32 bits version if such a driver is detected).
        if(!init) {
            init = true;
            const char* vendor = (const char*)glGetString(
                GL_VENDOR
            );
            use_glGetBufferParameteri64v = (
                strlen(vendor) < 5 || strncmp(vendor, "Intel", 5) != 0
            );
            if(!use_glGetBufferParameteri64v) {
                Logger::warn("GLSL")
                    << "Buggy Intel driver detected (working around...)"
                    << std::endl;
            }
        }
        GLint64 result=0;
        if(use_glGetBufferParameteri64v) {
            glGetBufferParameteri64v(target,GL_BUFFER_SIZE,&result);
        } else {
            GLint result32=0;
            glGetBufferParameteriv(target,GL_BUFFER_SIZE,&result32);
            result = GLint64(result32);
        }
        return result;
    }


    void update_buffer_object(
        GLuint& buffer_id, GLenum target, size_t new_size, const void* data,
        bool streaming
    ) {
        if(new_size == 0) {
            if(buffer_id != 0) {
                glDeleteBuffers(1, &buffer_id);
                buffer_id = 0;
            }
            return;
        }

        GLint64 size = 0;        
        if(buffer_id == 0) {
            glGenBuffers(1, &buffer_id);
            glBindBuffer(target, buffer_id);            
        } else {
            glBindBuffer(target, buffer_id);
            size = get_size_of_bound_buffer_object(target);
        }
        
        if(new_size == size_t(size)) {
            if(streaming) {
                //   Binding nil makes the GPU-side allocated buffer "orphan",
                // if there was a rendering operation currently using it, then
                // it can safely continue.
                glBufferData(target, GLsizeiptr(size), nil, GL_STREAM_DRAW);
                //   And here we bind a fresh new block of GPU-side memory.
                // See https://www.opengl.org/wiki/Buffer_Object_Streaming
                glBufferData(target, GLsizeiptr(size), data, GL_STREAM_DRAW);
            } else {
                glBufferSubData(target, 0, GLsizeiptr(size), data);
            }
        } else {
            glBufferData(
                target, GLsizeiptr(new_size), data, GL_STATIC_DRAW
            );
        }
    }

    void update_or_check_buffer_object(
        GLuint& buffer_id, GLenum target, size_t new_size, const void* data,
        bool update
    ) {
        if(update) {
            update_buffer_object(buffer_id, target, new_size, data);
        } else {
            glBindBuffer(target, buffer_id);
            if(new_size != size_t(get_size_of_bound_buffer_object(target))) {
                Logger::warn("OpenGL")
                    << "Buffer Object does not have the expected size."
                    << std::endl;
                Logger::warn("OpenGL")
                    << "An object was probably changed "
                    << "without notifying/updating the graphics."
                    << std::endl;
                Logger::warn("OpenGL")
                    << "Forcing Buffer Object update."
                    << std::endl;
                update_buffer_object(buffer_id, target, new_size, data);
            }
        }
    }

    void check_gl(const char* file, int line) {
        GLenum error_code = glGetError() ;
        bool has_opengl_errors = false ;
        while(error_code != GL_NO_ERROR) {
            has_opengl_errors = true ;
            Logger::err("OpenGL")
                << file << ":" << line << " " 
                << (char*)(gluErrorString(error_code)) << std::endl ;
            error_code = glGetError() ;
        }
        geo_argused(has_opengl_errors);
        // geo_debug_assert(!has_opengl_errors);
    }
}
