/*
 * Copyright (c) 2012-2016, Association Scientifique pour la Geologie et ses Applications (ASGA)
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *     * Redistributions in binary form must reproduce the above copyright
 *       notice, this list of conditions and the following disclaimer in the
 *       documentation and/or other materials provided with the distribution.
 *     * Neither the name of the <organization> nor the
 *       names of its contributors may be used to endorse or promote products
 *       derived from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL <COPYRIGHT HOLDER> BE LIABLE FOR ANY
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 *  Contacts:
 *
 *
 *
 *
 *     http://www.ring-team.org
 *
 *     RING Project
 *     Ecole Nationale Superieure de Geologie - GeoRessources
 *     2 Rue du Doyen Marcel Roubault - TSA 70605
 *     54518 VANDOEUVRE-LES-NANCY
 *     FRANCE
 */

#ifndef __RINGMESH_ASSERT__
#define __RINGMESH_ASSERT__

#include <ringmesh/common.h>

#include <string>

#include <geogram/basic/assert.h>

namespace RINGMesh {
    static void ringmesh_assertion_failed(
        const std::string& condition_string,
        const std::string& file,
        int line )
    {
#if WIN32
        DebugBreak() ;
#endif
        GEO::geo_assertion_failed( condition_string, file, line ) ;
    }

    static void ringmesh_should_not_have_reached(
        const std::string& file,
        int line )
    {
#if WIN32
        DebugBreak() ;
#endif
        GEO::geo_should_not_have_reached( file, line ) ;
    }
}

#define ringmesh_assert( x ) \
    {                                        \
        if( !( x ) ) {                                                 \
            ::RINGMesh::ringmesh_assertion_failed( # x, __FILE__, __LINE__ ) ;   \
        }                                                          \
    }

#define ringmesh_assert_not_reached \
    {                               \
        ::RINGMesh::ringmesh_should_not_have_reached( __FILE__, __LINE__ ) ;   \
    }

#ifdef RINGMESH_DEBUG
  #define ringmesh_debug_assert( x ) ringmesh_assert( x )
  #define ringmesh_debug_assert_not_reached ringmesh_assert_not_reached
#else
  #define ringmesh_debug_assert( x )
  #define ringmesh_debug_assert_not_reached
#endif

#endif
