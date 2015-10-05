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
 *     Ecole Nationale Superieure de Geologie - Georessources
 *     2 Rue du Doyen Marcel Roubault - TSA 70605
 *     54518 VANDOEUVRE-LES-NANCY 
 *     FRANCE
*/


#ifndef __RINGMESH_TYPES__
#define __RINGMESH_TYPES__

#include <geogram/basic/geometry.h>

namespace RINGMesh {

    enum Sign {
        NEGATIVE = -1, ZERO = 0, POSITIVE = 1
    } ;  

    template< class T >
    inline Sign sign( T x )
    {
        return ( x > 0 ) ? POSITIVE : ( ( x < 0 ) ? NEGATIVE : ZERO ) ;
    }


    /**********************************************************/
    /* @todo Why do we need to redefine those ? [JP]
     */


    typedef unsigned char byte ;
    /*! Generic pointer type */
    typedef byte* pointer ;

    /*! Null pointer type */
    //#define nil NULL // already in Geogram

    /*! Integer type with a width of 8 bits */
    typedef char int8 ;
    /*
     * @todo Review : I do not understand the point of these uncommented 
     * dummy variables. I would love to see them removed. To avoid warnings
     * it is sufficient and clearer for the reader to use the default vec3() [JP]
     */
    static int8 dummy_int8 ;

    /*! Integer type with a width of 16 bits */
    typedef short int16 ;
    static int16 dummy_int16 ;

    /*! Integer type with a width of 32 bits */
    typedef int int32 ;
    static int32 dummy_int32 ;

    /*! Integer type with a width of 64 bits */
    typedef long int int64 ;
    static int64 dummy_int64 ;

    /*! Unsigned integer type with a width of 8 bits */
    typedef unsigned char uint8 ;
    static uint8 dummy_uint8 ;

    /*! Unsigned integer type with a width of 16 bits */
    typedef unsigned short uint16 ;
    static uint16 dummy_uint16 ;

    /*! Unsigned integer type with a width of 32 bits */
    typedef unsigned int uint32 ;
    static uint32 dummy_uint32 ;

    /*! Unsigned integer type with a width of 64 bits */
    typedef unsigned long int uint64 ;
    static uint64 dummy_uint64 ;

    /*! Floating point type with a width of 32 bits */
    typedef float float32 ;
    static float32 dummy_float32 ;

    /*! Floating point type with a width of 64 bits */
    typedef double float64 ;
    static float64 dummy_float64 ;

    static bool dummy_bool ;

    const float32 big_float32 = 1e10f ;
    const float32 small_float32 = 1e-10f ;
    const float64 big_float64 = 1e20 ;
    const float64 small_float64 = 1e-20 ;
    const float64 epsilon = 1E-8 ;
    const float64 epsilon_sq = epsilon*epsilon ;

    typedef GEO::vec3 vec3 ;
    static vec3 dummy_vec3 ;

    typedef GEO::index_t index_t ;
    static index_t dummy_index_t ;
    typedef GEO::signed_index_t signed_index_t ;
    static signed_index_t dummy_signed_index_t ;

    const static index_t NO_ID = index_t( -1 ) ;

    const std::string surface_att_name = "region" ;
    const std::string region_att_name = "region" ;
    const std::string order_att_name = "order" ;

}

#endif
