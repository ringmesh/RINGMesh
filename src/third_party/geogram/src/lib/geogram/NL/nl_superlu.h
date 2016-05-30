/*
 *  Copyright (c) 2004-2010, Bruno Levy
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
 *     levy@loria.fr
 *
 *     ALICE Project
 *     LORIA, INRIA Lorraine, 
 *     Campus Scientifique, BP 239
 *     54506 VANDOEUVRE LES NANCY CEDEX 
 *     FRANCE
 *
 */

#ifndef OPENNL_SUPERLU_H
#define OPENNL_SUPERLU_H

#include "nl_private.h"

/**
 * \file geogram/NL/nl_superlu.h
 * \brief Internal OpenNL functions that interfaces SuperLU.
 */

/**
 * \brief Solves the system in the current OpenNL 
 *   context using SUPERLU.
 * \details This function should not be called directly by client code.
 *  To use SUPERLU, first call nlInitExtension("SUPERLU")
 *  then specify:
 *   - nlSolverParameteri(NL_SOLVER, NL_SUPERLU_EXT) 
 *     if no pre-ordering should be used
 *   - nlSolverParameteri(NL_SOLVER, NL_PERM_SUPERLU_EXT) 
 *     to use pre-ordering for general matrices
 *   - nlSolverParameteri(NL_SOLVER, NL_SYMMETRIC_SUPERLU_EXT) 
 *     to use pre-ordering for symmetric matrices
 * \retval NL_TRUE if solve was successful
 * \retval NL_FALSE otherwise
 */
NLboolean nlSolve_SUPERLU(void);

/**
 * \brief Initializes the SUPERLU extension
 * \details This dynamically loads the SuperLU 
 *  library available in the system (if available) and
 *  retreives the symbols in there. It supports SuperLU 3.x
 *  and SuperLU 4.x. 
 * \retval NL_TRUE if SUPERLU could be successfully
 *   dynamically loaded and all functions could be
 *   found in it.
 * \retval NL_FALSE otherwise.
 * \note For now, only implemented under Linux in 
 *  dynamic libraries mode
 */
NLboolean nlInitExtension_SUPERLU(void);


#endif
