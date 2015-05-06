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

/**
 * \file geogram_gfx/basic/GLSL.h
 * \brief Utilities for manipulating GLSL shaders.
 */

namespace GEO {

    namespace GLSL {

        /**
         * \brief Exception thrown when a task is canceled
         * \see Progress::cancel()
         */
        struct GEOGRAM_GFX_API GLSLCompileError : std::exception {
            /**
             * \brief Gets the string identifying the exception
             */
            virtual const char* what() const throw();
        };

        
        /**
         * \brief Compiles a shader for a specific target.
         * \details One may use \p source1 for library functions common 
         *  to different shaders (then the code of the shader is in 
         *   \p source2). 
         *  It may seem more natural to generate a shader object with library 
         *  functions, but OpenGL documentation does not recommend
         *  to do so (and it did not seem to work). Errors are detected and 
         *  displayed to std::err.
         * \param[in] target the OpenGL shader target ()
         * \param[in] source1 the source of the shader (ASCII string)
         * \param[in] source2 an optional additional source string or 0 
         *  if unused
         * \return the OpenGL opaque Id of the created shader object
         * \throw GLSLCompileError
         */
        GLuint GEOGRAM_GFX_API compile_shader(
            GLenum target, const char* source1, const char* source2=0
        );


        /**
         * \brief Creates a program from a zero-terminated list of shaders
         * \details Errors are detected and displayed to the Logger.
         * \param[in] shader the first shader of the list
         * \return the OpenGL opaque Id of the created program
         */
        GLuint GEOGRAM_GFX_API setup_program(GLuint shader, ...);
    }
}
