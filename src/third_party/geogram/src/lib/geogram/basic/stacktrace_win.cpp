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

#include <geogram/basic/common.h>

#ifdef GEO_OS_WINDOWS

#include <geogram/basic/stacktrace.h>

// Copyright (c) 2007 Howard Jeng
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following condition:
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
// THE SOFTWARE.

#include <iostream>
#include <string>
#include <windows.h>
#include <dbghelp.h>

#pragma comment(lib, "dbghelp.lib")

namespace {

    /**
     * \brief Gets a module path
     * \details Retrieves the fully qualified path for the file that contains
     * the specified module \p module.
     * \param[in] module handle to a module
     * \return the full path to the module file
     */
    std::string get_module_path(HMODULE module = 0) {
        char path_name[MAX_PATH] = {};
        DWORD size = GetModuleFileNameA(module, path_name, MAX_PATH);
        return std::string(path_name, size);
    }

    /**
     * \brief Prints the module at a given address
     * \details Prints the fully qualified path of the module file that
     * contains address \p program_counter in process \p process. The path
     * is printed to output stream \p os.
     * \param[in] os output stream to print the path to
     * \param[in] process handle to the current process
     * \param[in] program_counter current program counter in the process
     */
    void write_module_name(std::ostream& os, HANDLE process, DWORD64 program_counter) {
        DWORD64 module_base = SymGetModuleBase64(process, program_counter);
        if(module_base) {
            std::string module_name = get_module_path(reinterpret_cast<HMODULE>(module_base));
            if(!module_name.empty()) {
                os << module_name << "|";
            }
            else {
                os << "Unknown module|";
            }
        } else {
            os << "Unknown module|";
        }
    }

    /** Maximum size of a file path */
    const size_t MAX_NAMELEN = 1024;

    /**
     * \brief Data structure used to retrieve symbol names
     */
    struct Symbol {
        IMAGEHLP_SYMBOL64 is;
        char buffer[MAX_NAMELEN + 1];
    };

    /**
     * \brief Prints the symbol at a given address
     * \details Prints the name of the symbol (function name) that contains
     * address \p program_counter in process \p process. The symbol name
     * is printed to output stream \p os.
     * \param[in] os output stream to print the path to
     * \param[in] process handle to the current process
     * \param[in] program_counter current program counter in the process
     */
    void write_function_name(std::ostream& os, HANDLE process, DWORD64 program_counter) {
        Symbol sym;
        memset(&sym, 0, sizeof(Symbol));
        sym.is.SizeOfStruct = sizeof(IMAGEHLP_SYMBOL64);
        sym.is.MaxNameLength = MAX_NAMELEN;

        DWORD64 dummy = 0;
        if(SymGetSymFromAddr64(process, program_counter, &dummy, &sym.is) != FALSE) {
            os << sym.is.Name << "()";
        } else {
            os << "Unknown function";
        }
    }

    /**
     * \brief Prints the source file and line at a given address
     * \details Prints the path to the source file and the line number that
     * corresponds to address \p program_counter in process \p process. The
     * path and line number are printed to output stream \p os.
     * \param[in] os output stream to print the path to
     * \param[in] process handle to the current process
     * \param[in] program_counter current program counter in the process
     */
    void write_file_and_line(std::ostream& os, HANDLE process, DWORD64 program_counter) {
        IMAGEHLP_LINE64 ih_line = {sizeof(IMAGEHLP_LINE64)};
        DWORD dummy = 0;
        if(SymGetLineFromAddr64(process, program_counter, &dummy, &ih_line)) {
            os << "|" << ih_line.FileName << ":" << ih_line.LineNumber;
        }
    }
}

namespace GEO {

    void StackTrace::initialize() {
        SymInitialize(GetCurrentProcess(), 0, TRUE);
        DWORD options = SymGetOptions();
        options |= SYMOPT_LOAD_LINES;
        SymSetOptions(options);
    }

    void StackTrace::print_stack_trace(int skip) {
        print_stack_trace(std::cout, skip);
    }

    void StackTrace::print_stack_trace(std::ostream& os, int skip) {
        print_stack_trace(os, 0, skip);
    }

    void StackTrace::print_stack_trace(std::ostream& os, CONTEXT* pctx, int skip) {

        // Normally, the context is passed by the exception handler
        // If not, take the current execution context

        CONTEXT ctx;
        if(pctx != NULL) {
            ctx = *pctx;
        } else {
            memset(&ctx, 0, sizeof(CONTEXT));
            ctx.ContextFlags = CONTEXT_FULL;
            RtlCaptureContext(&ctx);
        }

        // Prepare the initial stack frame from the context

        STACKFRAME64 sf;
        memset(&sf, 0, sizeof(sf));

#if _M_X64
        os << "Stack trace (x64 architecture)" << std::endl;
        DWORD imageType = IMAGE_FILE_MACHINE_AMD64;
        sf.AddrPC.Offset = ctx.Rip;
        sf.AddrPC.Mode = AddrModeFlat;
        sf.AddrStack.Offset = ctx.Rsp;
        sf.AddrStack.Mode = AddrModeFlat;
        sf.AddrFrame.Offset = ctx.Rbp;
        // NOTE: Some stacktrace implementations take the Rsp register
        // It does not seem to impact the result -- TV
        // sf.AddrFrame.Offset = ctx.Rsp;
        sf.AddrFrame.Mode = AddrModeFlat;
#elif _M_IX86
        os << "Stack trace (i386 architecture)" << std::endl;
        DWORD imageType = IMAGE_FILE_MACHINE_I386;
        sf.AddrPC.Offset = ctx.Eip;
        sf.AddrPC.Mode = AddrModeFlat;
        sf.AddrStack.Offset = ctx.Esp;
        sf.AddrStack.Mode = AddrModeFlat;
        sf.AddrFrame.Offset = ctx.Ebp;
        sf.AddrFrame.Mode = AddrModeFlat;
#else
#error "Platform not supported!"
#endif

        HANDLE process = GetCurrentProcess();
        HANDLE thread = GetCurrentThread();

        os << std::uppercase;
        for(;;) {
            SetLastError(0);
            BOOL stack_walk_ok = StackWalk64(
                imageType, process, thread,
                &sf, &ctx,
                NULL,
                &SymFunctionTableAccess64,
                &SymGetModuleBase64,
                NULL
            );
            if(!stack_walk_ok || !sf.AddrFrame.Offset) {
                return;
            }

            if(skip) {
                --skip;
            } else {
                // write the address
                os << std::hex << reinterpret_cast<void*>(sf.AddrPC.Offset) << "|" << std::dec;

                write_module_name(os, process, sf.AddrPC.Offset);
                write_function_name(os, process, sf.AddrPC.Offset);
                write_file_and_line(os, process, sf.AddrPC.Offset);

                os << "\n";
            }
        }
    }
}

#endif

