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

#ifndef __GEOGRAM_BASIC_STACKTRACE__
#define __GEOGRAM_BASIC_STACKTRACE__

#include <geogram/basic/common.h>
#include <iosfwd>
#ifdef GEO_OS_WINDOWS
#include <windows.h>
#endif

/**
 * \file geogram/basic/stacktrace.h
 * \brief Functions to dump a stacktrace
 */

namespace GEO {

    /**
     * \brief Stack trace dumper
     * \details This class provides several functions to dump the current
     * stacktrace. It is especially useful when the program terminates abnormally.
     *
     * \internal
     * The windows implementation is borrowed from a very good article about the
     * Visual C++ Exception Model. It has extended to support both I386 and X64
     * architectures. The symbol extraction function write_function_name() did not
     * work and has been replaced by a working one.
     *
     * For more details, see: http://members.gamedev.net/sicrane/articles/exception.html
     */
    struct GEOGRAM_API StackTrace {

        /**
         * \brief Initializes the symbol table
         * \details This function must be called once before printing a stack
         * trace in order to perform necessary initializations (preload
         * application symbols, ...)
         * \note This function is only available on Windows
         */
        static void initialize();

        /**
         * \brief Prints the current stack trace to the standard output
         * \param[in] skip numbers of frames to skip at the beginning of the stack
         * trace (default is 0)
         */
        static void print_stack_trace(int skip = 0);

        /**
         * \brief Prints the current stack trace to a stream
         * \param[in] os stream to send the stack trace to
         * \param[in] skip numbers of frames to skip at the beginning of the stack
         * trace (default is 0)
         */
        static void print_stack_trace(std::ostream& os, int skip = 0);

#ifdef GEO_OS_WINDOWS

        /**
         * \brief Prints the current stack trace to a stream with context
         * \param[in] os stream to send the stack trace to
         * \param[in] pctx current execution context. In the scope of exception
         * handlers, the context is given by the exception record. When NULL is
         * passed, the execution context is determined automatically.
         * \param[in] skip numbers of frames to skip at the beginning of the stack
         * trace (default is 0)
         * \note This function is only available on Windows
         */
        static void print_stack_trace(std::ostream& os, CONTEXT* pctx, int skip = 0);

#endif
    };
}

#endif

