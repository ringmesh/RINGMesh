/*
 * Copyright (c) 2012-2015, Association Scientifique pour la Geologie et ses Applications (ASGA)
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
 *     Arnaud.Botella@univ-lorraine.fr 
 *     Antoine.Mazuyer@univ-lorraine.fr 
 *     Jeanne.Pellerin@wias-berlin.de
 *
 *     http://www.ring-team.org
 *
 *     RING Project
 *     Ecole Nationale Superieure de Geologie - GeoRessources
 *     2 Rue du Doyen Marcel Roubault - TSA 70605
 *     54518 VANDOEUVRE-LES-NANCY 
 *     FRANCE
*/


#ifndef __RINGMESH_TYPES__
#define __RINGMESH_TYPES__

#include <geogram/basic/geometry.h>
#include <geogram/basic/numeric.h>

namespace RINGMesh {
    /* If you need interger of 8bits of any other one
     * it is sufficient to write using GEO::Numeric::uint8 in your file. 
     * 
     * Dummy variables were removed, the pollute the namespace and 
     * it is quite easy to do without them.
     */

    // Basic types used inn RINGMesh
    // Using definitions of Geogram/basic/numeric.h   
    using GEO::Numeric::float32 ;
    using GEO::Numeric::float64 ;

    using GEO::Numeric::max_float32 ;
    using GEO::Numeric::min_float32 ;
    using GEO::Numeric::max_float64 ;
    using GEO::Numeric::min_float64 ;

    const float64 epsilon = 1E-8 ;
    const float64 epsilon_sq = epsilon*epsilon ;

    // This is an array of 3 doubles
    using GEO::vec3 ;
    // This is an unsigned int
    using GEO::index_t ;
    // This is an int
    using GEO::signed_index_t ;

    // This is the value used in RINGMesh for a invalid index
    const static index_t NO_ID(-1) ;
}

#endif
