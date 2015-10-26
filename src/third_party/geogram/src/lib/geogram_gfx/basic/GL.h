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

#ifndef __GEOGRAM_GFX_BASIC_GL__
#define __GEOGRAM_GFX_BASIC_GL__

#include <geogram_gfx/basic/common.h>
#include <geogram_gfx/api/defs.h>
#include <geogram_gfx/third_party/glew/glew.h>
#include <geogram/basic/geometry.h>

/**
 * \file geogram_gfx/basic/common.h
 * \brief Some utility functions for OpenGL graphics.
 */

namespace GEO {
    /**
     * \brief Sends a vertex to OpenGL.
     * \param[in] v a const reference to the vertex to be sent.
     * \note This uses the old pipeline (glBegin() / glEnd() calls).
     */
    inline void glVertex(const vec3& v) {
        glVertex3dv(v.data());
    }

    /**
     * \brief Sends a vertex to OpenGL.
     * \param[in] v a const reference to the vertex to be sent, in
     *  homogeneous coordinates (4d).
     * \note This uses the old pipeline (glBegin() / glEnd() calls).
     */
    inline void glVertex(const vec4& v) {
        glVertex4dv(v.data());
    }

    /**
     * \brief Sends a RGB color to OpenGL.
     * \param[in] v a const reference to the color to be sent.
     * \note This uses the old pipeline (glBegin() / glEnd() calls).
     */
    inline void glColor(const vec3& v) {
        glColor3dv(v.data());
    }

    /**
     * \brief Sends a RGBA color to OpenGL.
     * \param[in] v a const reference to the color to be sent.
     * \note This uses the old pipeline (glBegin() / glEnd() calls).
     */
    inline void glColor(const vec4& v) {
        glColor4dv(v.data());
    }

    /**
     * \brief Sends a normal to OpenGL.
     * \param[in] v a const reference to the normal to be sent.
     * \note This uses the old pipeline (glBegin() / glEnd() calls).
     */
    inline void glNormal(const vec3& v) {
        glNormal3dv(v.data());
    }

    /**
     * \brief Multiplies the current OpenGL matrix
     *   with another one.
     * \param[in] m a const reference to the matrix.
     * \note m is transposed before being sent to OpenGL
     *  because Geogram uses the convention with column
     *  vectors and OpenGL the convention with row vectors
     *  to represent the transformed points.
     */
    void GEOGRAM_API glMultMatrix(const mat4& m);

    /**
     * \brief Replaces the current OpenGL matrix
     *   with a user defined one.
     * \param[in] m a const reference to the matrix.
     * \note m is transposed before being sent to OpenGL
     *  because Geogram uses the convention with column
     *  vectors and OpenGL the convention with row vectors
     *  to represent the transformed points.
     */
    void GEOGRAM_API glLoadMatrix(const mat4& m);


    /**
     * \brief Gets the size (in bytes) of the OpenGL buffer 
     *  bound to a specified target.
     * \param[in] target buffer object target 
     *   (GL_ARRAY_BUFFER, GL_INDEX_BUFFER ...)
     * \return the size in bytes of the buffer object bound 
     *  to \p target.
     */
    GLint64 GEOGRAM_API get_size_of_bound_buffer_object(GLenum target);

    /**
     * \brief Updates the content of an OpenGL buffer object, 
     *   and resizes it if need be.
     * \param[in,out] buffer_id OpenGL opaque id of the buffer object. 
     *   0 means uninitialized.
     *   may be changed on exit if the buffer needed to be created or
     *   destroyed.
     * \param[in] target buffer object target 
     *   (GL_ARRAY_BUFFER, GL_INDEX_BUFFER ...)
     * \param[in] new_size of the buffer data, in bytes
     * \param[in] data pointer to the data to be copied into the buffer, 
     *  of length new_size
     */
    void GEOGRAM_API update_buffer_object(
        GLuint& buffer_id, GLenum target, size_t new_size, const void* data
    );


    /**
     * \brief Updates the content of an OpenGL buffer object, 
     *   and resizes it if need be, or tests whether it has the
     *   size it should have.
     * \param[in,out] buffer_id OpenGL opaque id of the buffer object. 
     *   0 means uninitialized.
     *   may be changed on exit if the buffer needed to be created or
     *   destroyed.
     * \param[in] target buffer object target 
     *   (GL_ARRAY_BUFFER, GL_INDEX_BUFFER ...)
     * \param[in] new_size of the buffer data, in bytes
     * \param[in] data pointer to the data to be copied into the buffer, 
     *  of length new_size
     * \param[in] update 
     *  - if false, the buffer will be updated, and resized
     *    if need be. 
     *  - if true, the size of the buffer will be tested, and an error 
     *    message will be displayed in the logger if it does not match
     *    the specified size (and update will be forced).
     */
    void GEOGRAM_API update_or_check_buffer_object(
        GLuint& buffer_id, GLenum target, size_t new_size, const void* data,
        bool update
    );
}

#endif

