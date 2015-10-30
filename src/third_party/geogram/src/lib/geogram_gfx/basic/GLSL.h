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

#include <geogram_gfx/basic/common.h>
#include <geogram/basic/numeric.h>

/**
 * \file geogram_gfx/basic/GLSL.h
 * \brief Utilities for manipulating GLSL shaders.
 */

namespace GEO {

    namespace GLSL {

        /**
         * \brief Exception thrown when a GLSL shader fails to
         *  compiled.
         * \details Can occur when OpenGL driver or hardware
         *  does not support some features.
         */
        struct GEOGRAM_GFX_API GLSLCompileError : std::exception {
            /**
             * \brief Gets the string identifying the exception
             */
            virtual const char* what() const throw();
        };

        
        /**
         * \brief Compiles a shader for a specific target.
         * \details One can split the source of the shader into
         *  different strings, one of them being used for library
         *  functions common to different shaders.
         *  It may seem more natural to generate a shader object with library 
         *  functions, but OpenGL documentation does not recommend
         *  to do so (and it did not seem to work). Errors are detected and 
         *  displayed to std::err.
         * \param[in] target the OpenGL shader target ()
         * \param[in] sources an array of pointer to ASCII strings that contain 
         *    the source of the shader 
         * \param[in] nb_sources number of strings in \p sources
         * \return the OpenGL opaque Id of the created shader object
         * \throw GLSLCompileError
         */
        GLuint GEOGRAM_GFX_API compile_shader(
            GLenum target, const char** sources, index_t nb_sources
        );

        /**
         * \brief Compiles a shader for a specific target.
         * \param[in] target the OpenGL shader target ()
         * \param[in] source the source of the shader as an ASCII string
         * \return the OpenGL opaque Id of the created shader object
         * \throw GLSLCompileError
         */
        inline GLuint GEOGRAM_GFX_API compile_shader(
            GLenum target, const char* source
        ) {
            return compile_shader(target, &source, 1);
        }


        /**
         * \brief Compiles a shader for a specific target.
         * \details One can split the source of the shader into
         *  different strings, one of them being used for library
         *  functions common to different shaders.
         *  It may seem more natural to generate a shader object with library 
         *  functions, but OpenGL documentation does not recommend
         *  to do so (and it did not seem to work). Errors are detected and 
         *  displayed to std::err.
         * \param[in] target the OpenGL shader target ()
         * \param[in] source1, source2 ASCII strings that will be 
         *  concatened to form the source of the shader.
         * \return the OpenGL opaque Id of the created shader object
         * \throw GLSLCompileError
         */
        GLuint GEOGRAM_GFX_API compile_shader(
            GLenum target, const char* source1, const char* source2
        );

        /**
         * \brief Compiles a shader for a specific target.
         * \details One can split the source of the shader into
         *  different strings, one of them being used for library
         *  functions common to different shaders.
         *  It may seem more natural to generate a shader object with library 
         *  functions, but OpenGL documentation does not recommend
         *  to do so (and it did not seem to work). Errors are detected and 
         *  displayed to std::err.
         * \param[in] target the OpenGL shader target ()
         * \param[in] source1, source2, source3 ASCII strings that will be 
         *  concatened to form the source of the shader.
         * \return the OpenGL opaque Id of the created shader object
         * \throw GLSLCompileError
         */
        GLuint GEOGRAM_GFX_API compile_shader(
            GLenum target,
            const char* source1, const char* source2, const char* source3
        );
        
        /**
         * \brief Creates a GLSL program from a zero-terminated list of shaders
         * \details Errors are detected and displayed to the Logger.
         * \param[in] shader the first shader of the list
         * \return the OpenGL opaque Id of the created program
         */
        GLuint GEOGRAM_GFX_API create_program_from_shaders(GLuint shader, ...);


        /**
         * \brief Creates a GLSL program from a string.
         * \details The string may contain several shaders. Each shader
         *   is delimited by begin-end statements: #begin(SHADER_TYPE) / #end(SHADER_TYPE)
         *   where SHADER_TYPE is one of GL_VERTEX_SHADER, GL_FRAGMENT_SHADER, GL_GEOMETRY_SHADER,
         *   GL_TESS_CONTROL_SHADER, GL_TESS_EVALUATION_SHADER.
         * \param[in,out] string the combined shaders that constitute the program. 
         * \param[in] copy_string if true, the input string is copied internally. The function
         *   temporarily modifies the input string (and then restores it on exit). This may
         *   be forbidden when input string is a constant char array (string litteral in source
         *   code). In this case, the input string is copied to a temporary buffer.
         * \return the OpenGL opaque Id of the created shader object
         * \throw GLSLCompileError
         */
        GLuint GEOGRAM_GFX_API create_program_from_string(
            const char* string, bool copy_string = true
        );

        /**
         * \brief Creates a GLSL program from a file.
         * \details The file contains a list of shaders, delimited by
         *   begin-end statements (see setup_program_from_string()).
         * \param[in] filename the name of the file
         * \throw GLSLCompileError
         */
        GLuint GEOGRAM_GFX_API create_program_from_file(const std::string& filename);


        /**
         * \brief Sets a uniform variable in a shader by name.
         * \param[in] shader_id the handle to the GLSL shader
         * \param[in] name the name of the uniform variable,
         *   as specified in the GLSL source of the shader.
         * \param[in] val the value of the parameter
         * \tparam T the type of the parameter. Needs to match
         *  the type of the uniform parameter in the GLSL source.
         */
        template <class T> inline bool set_program_uniform_by_name(
            GLuint shader_id, const char* name, T val
        ) {
            geo_argused(shader_id);
            geo_argused(name);
            geo_argused(val);
            geo_assert_not_reached;
            return false;
        }

        template<> inline bool set_program_uniform_by_name(
            GLuint shader_id, const char* name, bool val
        ) {
            GLint location = glGetUniformLocation(shader_id, name) ;
            if(location < 0) {
                return false ;
            }
            glUseProgram(shader_id);
            glUniform1i(location, val ? 1 : 0) ;
            glUseProgram(0);
            return true;
        }

        template<> inline bool set_program_uniform_by_name(
            GLuint shader_id, const char* name, float val
        ) {
            GLint location = glGetUniformLocation(shader_id, name) ;
            if(location < 0) {
                return false ;
            }
            glUseProgram(shader_id);
            glUniform1f(location, val) ;
            glUseProgram(0);
            return true;
        }

        template<> inline bool set_program_uniform_by_name(
            GLuint shader_id, const char* name, double val
        ) {
            GLint location = glGetUniformLocation(shader_id, name) ;
            if(location < 0) {
                return false ;
            }
            glUseProgram(shader_id);
            glUniform1f(location, float(val)) ;
            glUseProgram(0);
            return true;
        }
        
        template<> inline bool set_program_uniform_by_name(
            GLuint shader_id, const char* name, int val
        ) {
            GLint location = glGetUniformLocation(shader_id, name) ;
            if(location < 0) {
                return false ;
            }
            glUseProgram(shader_id);
            glUniform1i(location, val) ;
            glUseProgram(0);
            return true;
        }

        /**
         * \brief Sets an array of uniform variables in a shader by name.
         * \param[in] shader_id the handle to the GLSL shader
         * \param[in] name the name of the uniform variable,
         *   as specified in the GLSL source of the shader.
         * \param(in] count number of values
         * \param[in] values a pointer to an array of values of size count
         */
        inline bool set_program_uniform_by_name(
            GLuint shader_id, const char* name, index_t count, float* values
        ) {
            GLint location = glGetUniformLocation(shader_id, name) ;
            if(location < 0) {
                return false ;
            }
            glUseProgram(shader_id);
            glUniform1fv(location, GLsizei(count), values) ;
            glUseProgram(0);
            return true;
        }

        
    }
}
