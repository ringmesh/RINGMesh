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

#include <geogram_gfx/basic/GLSL.h>
#include <geogram/basic/logger.h>
#include <cstdarg>

namespace GEO {

        namespace GLSL {

            /*************************************************************/

            const char* GLSLCompileError::what() const throw() {
                return "GLSL Compile Error";
            }
            
            /*************************************************************/

            GLuint compile_shader(
                GLenum target, const char* source1, const char* source2
            ) {
                const char* sources[2];
                sources[0] = source1;
                sources[1] = source2;
                return compile_shader(target, &sources[0], 2);
            }

            GLuint compile_shader(
                GLenum target,
                const char* source1, const char* source2, const char* source3
            ) {
                const char* sources[3];
                sources[0] = source1;
                sources[1] = source2;
                sources[2] = source3;
                return compile_shader(target, &sources[0], 2);
            }
            
            GLuint compile_shader(
                GLenum target, const char** sources, index_t nb_sources
            ) {
                GLuint s_handle = glCreateShader(target);
                if(s_handle == 0) {
                    Logger::err("GLSL") << "Could not create shader"
                                        << std::endl;
                    exit(1);
                }
                glShaderSource(s_handle, (GLsizei)nb_sources, sources, 0);
                glCompileShader(s_handle);
                GLint compile_status;
                glGetShaderiv(s_handle, GL_COMPILE_STATUS, &compile_status);
                if(!compile_status) {
                    GLchar compiler_message[4096];
                    glGetShaderInfoLog(
                        s_handle, sizeof(compiler_message), 0, compiler_message
                        );
                    Logger::err("GLSL")
                        << "compiler status :"
                        << compile_status << std::endl;
                    Logger::err("GLSL")
                        << "compiler message:"
                        << compiler_message << std::endl;
                    glDeleteShader(s_handle);
                    s_handle = 0;
                    throw GLSLCompileError();
                }
                return s_handle;
            }


            GLuint setup_program(GLuint shader, ...) {
                GLuint program = glCreateProgram();
                va_list args;
                va_start(args,shader);
                while(shader != 0) {
                    glAttachShader(program, shader);
                    shader = va_arg(args, GLuint);
                }
                va_end(args);
                glLinkProgram(program);
                GLint link_status;
                glGetProgramiv(program, GL_LINK_STATUS, &link_status);
                if(!link_status) {
                    GLchar linker_message[4096];
                    glGetProgramInfoLog(
                        program, sizeof(linker_message), 0, linker_message
                    );
                    Logger::err("GLSL") << "linker status :"
                                        << link_status << std::endl;
                    Logger::err("GLSL") << "linker message:"
                                        << linker_message << std::endl;
                    glDeleteProgram(program);
                    program = 0;
                }
                return program;
            }

            /*****************************************************************/
            
            
        }
}

