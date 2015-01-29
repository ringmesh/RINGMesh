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

#if defined(GEO_OS_UNIX) && !defined(GEO_OS_ANDROID)

#include <geogram/basic/stacktrace.h>

/** \internal
 * By default, print_stack_trace() executes program addr2line to extract
 * symbol information. Setting STACKTRACE_WITH_BFD to 1 uses the BFD API to
 * extract symbol information. Both solutions have their own advantages and
 * drawbacks, see below for more details.
 */
#define STACKTRACE_WITH_BFD 0

#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif

#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <unistd.h>
#include <execinfo.h>
#include <sys/wait.h>

#if STACKTRACE_WITH_BFD
#include <bfd.h>
#include <cxxabi.h>
#endif

namespace {

    /** Path top the current executable */
    char executable_path[1024] = {-1};

    /** Maximum number of frames to display */
    const int MAX_FRAMES = 200;

    /** Table of stackframe addresses */
    void* addresses[MAX_FRAMES + 1];

    /**
     * Returns the path to the current executable
     */
    const char* current_executable() {

        if(executable_path[0] == -1) {
            ssize_t l = readlink(
                "/proc/self/exe",
                executable_path,
                sizeof(executable_path)
            );
            if(l == -1) {
                std::cerr << "failed to find executable: " << strerror(errno) << std::endl;
                return 0;
            }
            executable_path[l] = 0;
        }

        return executable_path;
    }

#if STACKTRACE_WITH_BFD

    // Globals retained across calls to resolve.
    bfd* abfd = 0;
    asymbol** syms = 0;
    asection* text = 0;

    /**
     * Symbol information returned by resolve_symbol()
     */
    struct symbol_info {
        /** Demangled function name */
        std::string file_name;
        /** Path to the source file that defines the symbol (debug mode only) */
        std::string func_name;
        /** Line number in the source file where the symbol is defined (debug mode only) */
        size_t line;
    };

    /**
     * Retrieves function name, file and line for a given address.
     *
     * Searches the given \p address in the executable BFD data structure  and
     * returns the demangled function name, filename and line number (in debug
     * mode only).
     *
     * \note Implementation notes:
     * To make this code compile the following Linux packages must be
     * available on the system:
     * - binutils-devel
     * - zlib-static (for static builds)
     *
     * \note In addition, the executable must be linked with the following
     * libraries: -lbfd -liberty -ldl -lz
     */
    bool resolve_symbol(const void* address, symbol_info& info) {
        if(abfd == 0) {
            const char* ename = current_executable();

            bfd_init();

            abfd = bfd_openr(ename, 0);
            if(!abfd) {
                perror("bfd_openr failed");
                return false;
            }

            /* oddly, this is required for it to work... */
            bfd_check_format(abfd, bfd_object);

            long storage_needed = bfd_get_symtab_upper_bound(abfd);
            syms = (asymbol**) malloc(storage_needed);
            /*long cSymbols =*/ (void) bfd_canonicalize_symtab(abfd, syms);

            text = bfd_get_section_by_name(abfd, ".text");
        }

        ssize_t offset = ((bfd_vma) address) - text->vma;
        if(offset > 0) {
            const char* file;
            const char* func;
            unsigned int line;
            if(bfd_find_nearest_line(abfd, text, syms, offset, &file, &func, &line)) {

                info.file_name = file ? file : "??";

                if(func == 0) {
                    info.line = 0;
                    info.func_name = "??";
                } else {
                    info.line = line;

                    // Demangle the function name
                    int status;
                    char* demangled = abi::__cxa_demangle(func, 0, 0, &status);
                    if(demangled == 0) {
                        info.func_name = func;
                    } else {
                        info.func_name = demangled;
                        free(demangled);
                    }
                }

                return true;
            }
        }

        return false;
    }

#else

    /**
     * Reads standard output from a command
     * \param[in] command command to execute
     * \param[in] output array of lines that will contain the command standard output
     * \retval true if the command successfully executed
     * \retval false if the command failed or exited with a non zero status code
     */
    bool read_from_pipe(const std::string& command, std::vector<std::string>& output) {
        FILE* fp = popen(command.c_str(), "r");
        if(fp == NULL) {
            std::cerr << "Failed to open pipe: " << strerror(errno) << std::endl;
            return false;
        }

        std::string line;
        char buf[1024];
        while(fgets(buf, sizeof(buf), fp)) {
            line += buf;
            size_t n = strlen(buf);
            if(n > 1 && buf[n - 1] == '\n') {
                output.push_back(line);
                line.clear();
            }
        }

        if(!line.empty()) {
            output.push_back(line);
        }

        int status = pclose(fp);
        if(status == -1) {
            std::cerr << "Child close error: " << strerror(errno) << std::endl;
            return false;
        }

        if(WIFEXITED(status)) {
            if(WEXITSTATUS(status) == 0) {
                return true;
            } else {
                std::cerr << "Child exited with status " << WEXITSTATUS(status) << std::endl;
                return false;
            }
        }
        if(WIFSIGNALED(status)) {
            std::cerr << "Child killed - signal " << WTERMSIG(status) << std::endl;
            return false;
        }
        if(WIFSTOPPED(status)) {
            std::cerr << "Child stopped - signal " << WSTOPSIG(status) << std::endl;
            return false;
        }
#ifdef WIFCONTINUED
        else if(WIFCONTINUED(status)) {
            /* Not all implementations support this */
            std::cerr << "Child continued" << std::endl;
            return false;
        }
#endif

        /* Non-standard case -- may never happen */
        std::cerr << "Unexpected status 0x" << std::hex << status << std::dec << std::endl;
        return false;
    }

#endif
}

namespace GEO {

    void StackTrace::initialize() {
    }

    void StackTrace::print_stack_trace(int skip) {
        print_stack_trace(std::cout, skip);
    }

    void StackTrace::print_stack_trace(std::ostream& os, int skip) {

        int nb_frames = backtrace(addresses, MAX_FRAMES);
        if(nb_frames == 0) {
            std::cerr << "backtrace: stack frame unavailable: " << strerror(errno) << std::endl;
            return;
        }

        // Adjust the value of frames to skip

        if(skip <= 0) {
            os << "Stacktrace (" << nb_frames << " stack frames)\n";
            skip = 0;
        } else {
            os << "Stacktrace (" << nb_frames << " stack frames, " << skip << " skipped)\n";
        }

#if STACKTRACE_WITH_BFD

        // Use the BFD API to extract symbol information from the executable.
        // The advantage of this solution compared to calling addr2line (see
        // below) is that we do not not need to call an external program
        // through a pipe. The drawbacks is that it requires to link the
        // executable with libraries: -lbfd -liberty -ldl -lz

        for(int i = skip; i < nb_frames; ++i) {
            symbol_info info;
            if(resolve_symbol(addresses[i], info)) {
                os << addresses[i] << ": " << info.func_name << " at " << info.file_name << ":" << info.line << std::endl;
            } else {
                os << addresses[i] << ": symbol information not available" << std::endl;
            }
        }

#else

        // Use program addr2line to extract symbol information from the
        // executable. The advantage of this solution is that it gives detailed
        // information about inlined functions that cannot be obtained with
        // resolve_symbol(). The drawback is that we must execute an external
        // program from the current executable which can be in an unstable state.

        std::ostringstream cmd;
        cmd << "/usr/bin/addr2line -a -C -f -i -p -e '" << current_executable() << "'";

        for(int i = skip; i < nb_frames; ++i) {
            cmd << " " << addresses[i];
        }

        // os << "Executing: " << cmd.str() << std::endl;

        std::vector<std::string> output;
        if(read_from_pipe(cmd.str(), output)) {
            for(std::vector<std::string>::const_iterator i = output.begin(); i != output.end(); ++i) {
                os << "= " << *i;
            }
            return;
        }

        // If the call to addr2line failed, we call backtrace_symbols() to get
        // information about the stack trace symbols. NOTE: backtrace_symbols()
        // does not work with static executables. Moreover it REQUIRES to link the
        // executable with -rdynamic to export all symbols. If this is not the
        // case, symbols will not be visible.

        char** symbols = backtrace_symbols(addresses, nb_frames);
        if(symbols == 0) {
            std::cerr << "backtrace_symbols: stack frame symbols unavailable: " << strerror(errno) << std::endl;
            return;
        }
        for(int i = skip; i < nb_frames; ++i) {
            os << symbols[i] << std::endl;
        }
        free(symbols);

#endif
    }
}

#endif

