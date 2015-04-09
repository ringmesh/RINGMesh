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
 *     http://www.gocad.org
 *
 *     GOCAD Project
 *     Ecole Nationale Superieure de Geologie - Georessources
 *     2 Rue du Doyen Marcel Roubault - TSA 70605
 *     54518 VANDOEUVRE-LES-NANCY 
 *     FRANCE
*/


#ifndef __RINGMESH_TYPES__
#define __RINGMESH_TYPES__

#include <geogram/basic/geometry.h>

namespace RINGMesh {

    //    _____
    //   |_   _|  _ _ __  ___ ___
    //     | || || | '_ \/ -_|_-<
    //     |_| \_, | .__/\___/__/
    //         |__/|_|
    typedef unsigned char byte ;
    /** Generic pointer type */
    typedef byte* pointer ;

    /** Null pointer type */
    //#define nil NULL // already in Geogram

    /** Integer type with a width of 8 bits */
    typedef char int8 ;
    static int8 dummy_int8 ;

    /** Integer type with a width of 16 bits */
    typedef short int16 ;
    static int16 dummy_int16 ;

    /** Integer type with a width of 32 bits */
    typedef int int32 ;
    static int32 dummy_int32 ;

    /** Integer type with a width of 64 bits */
    typedef long int int64 ;
    static int64 dummy_int64 ;

    /** Unsigned integer type with a width of 8 bits */
    typedef unsigned char uint8 ;
    static uint8 dummy_uint8 ;

    /** Unsigned integer type with a width of 16 bits */
    typedef unsigned short uint16 ;
    static uint16 dummy_uint16 ;

    /** Unsigned integer type with a width of 32 bits */
    typedef unsigned int uint32 ;
    static uint32 dummy_uint32 ;

    /** Unsigned integer type with a width of 64 bits */
    typedef unsigned long int uint64 ;
    static uint64 dummy_uint64 ;

    /** Floating point type with a width of 32 bits */
    typedef float float32 ;
    static float32 dummy_float32 ;

    /** Floating point type with a width of 64 bits */
    typedef double float64 ;
    static float64 dummy_float64 ;

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
    //    ___
    //   | __|_ _ _  _ _ __  ___
    //   | _|| ' \ || | '  \(_-<
    //   |___|_||_\_,_|_|_|_/__/
    //

    enum ElementType {
        VERTEX  = -1,
        LINE    =  0,
        TRGL    =  1,
        QUAD    =  2,
        TETRA   =  3,
        PYRAMID =  4,
        PRISM   =  5,
        HEXA    =  6
    } ;


    enum TetraMethod {
        TetGen, MG_Tetra
    } ;

    //     ___     _ _      _               _      _
    //    / __|___| | |  __| |___ _____ _ _(_)_ __| |_ ___ _ _ ___
    //   | (__/ -_) | | / _` / -_|_-< _| '_| | '_ \  _/ _ \ '_(_-<
    //    \___\___|_|_| \__,_\___/__|__|_| |_| .__/\__\___/_| /__/
    //                                       |_|
    struct CellDescriptor {
        uint8 nb_vertices ;
        uint8 nb_facets ;
        uint8 nb_vertices_in_facet[6] ;
        uint8 facet[6][4] ;
        uint8 nb_edges ;
        uint8 edge[12][2] ;
        uint8 edges_in_facet[6][4] ;
    } ;

    static CellDescriptor line_descriptor = {
        2,              // nb_vertices
        1,              // nb_facets
        { 2 }, // nb_vertices in facet
        { { 0, 1 } },  // facets
        1,              // nb_edges
        { { 0, 1 } } , // edges
        { { 0 } } //edges_in_facet
    } ;
    static CellDescriptor trgl_descriptor = {
        3,              // nb_vertices
        1,              // nb_facets
        { 3 }, // nb_vertices in facet
        { { 0, 1, 2 } },  // facets
        3,              // nb_edges
        { { 0, 1 }, { 1, 2 }, { 2, 0 } } , // edges
        { { 0, 1, 2 } } //edges_in_facet
    } ;
    static CellDescriptor quad_descriptor = {
        4,              // nb_vertices
        1,              // nb_facets
        { 4 }, // nb_vertices in facet
        { { 0, 1, 2, 3 } },  // facets
        4,              // nb_edges
        { { 0, 1 }, { 1, 2 }, { 2, 3 }, { 3, 0} } , // edges
        { { 0, 1, 2, 3 } } //edges_in_facet
    } ;
    static CellDescriptor tetra_descriptor = {
        4,              // nb_vertices
        4,              // nb_facets
        { 3, 3, 3, 3 }, // nb_vertices in facet
        { { 3, 1, 2 }, { 3, 2, 0 }, { 0, 1, 3 }, { 1, 0, 2 } },  // facets
        6,              // nb_edges
        { { 1, 3 }, { 3, 2 }, { 2, 1 }, { 0, 2 }, { 3, 0 }, { 1, 0 } } , // edges
        { { 2, 0, 1 }, { 1, 4, 3 }, { 4, 0, 5 }, { 3, 2, 5 } } //edges_in_facet
    } ;
    static CellDescriptor hexa_descriptor = {
        8,                      // nb_vertices
        6,                      // nb_facets
        { 4, 4, 4, 4, 4, 4 },   // nb_vertices in facet
        { { 0, 2, 6, 4 }, { 3, 1, 5, 7 }, { 1, 0, 4, 5 }, { 2, 3, 7, 6 },
          { 1, 3, 2, 0 }, { 4, 6, 7, 5 } },  // facets
        12,                     // nb_edges
        { { 0, 2 }, { 2, 6 }, { 6, 4 }, { 4, 0 }, { 3, 1 }, { 1, 5 }, { 5, 7 },
          { 7, 3 }, { 1, 0 }, { 4, 5 }, { 2, 3 }, { 7, 6 } },   // edges
        { { 0, 1, 2, 3 }, { 4, 5, 6, 7 }, { 3, 8, 5, 9 },
          { 1, 11, 7, 10 }, { 0, 10, 4, 8 }, { 2, 9, 6, 11 } }  //edges_in_facet
    } ;
    static CellDescriptor prism_descriptor = {
        6,                  // nb_vertices
        5,                  // nb_facets
        { 3, 3, 4, 4, 4 },  // nb_vertices in facet
        { { 0, 1, 2 }, { 3, 5, 4 }, { 0, 3, 4, 1 }, { 0, 2, 5, 3 }, { 1, 4, 5, 2 } },  // facets
        9,                  // nb_edges
        { { 0, 1 }, { 1, 2 }, { 2, 0 }, { 3, 5 }, { 5, 4 },
          { 4, 3 }, { 0, 3 }, { 4, 1 }, { 2, 5 } },  // edges
        { { 0, 1, 2 }, { 5, 4, 3 }, { 6, 5, 7, 0 },
          { 2, 8, 3, 6 }, { 7, 4, 8, 1 } }  //edges_in_facet
    } ;
    static CellDescriptor pyramid_descriptor = {
        5,                  // nb_vertices
        5,                  // nb_facets
        { 4, 3, 3, 3, 3 },  // nb_vertices in facet
        { { 0, 1, 2, 3 }, { 0, 4, 1 }, { 0, 3, 4 }, { 2, 4, 3 }, { 2, 1, 4 } },    // facets
        8,                  // nb_edges
        { { 0, 1 }, { 1, 2 }, { 2, 3 }, { 3, 0 },
          { 0, 4 }, { 4, 1 }, { 3, 4 }, { 2, 4 } },  // edges
        { { 0, 1, 2, 3 }, { 4, 5, 0 }, { 3, 6, 4 },
          { 6, 7, 2 }, { 1, 5, 7 } }  //edges_in_facet
    } ;

    static ElementType cell_types[9] = {
        VERTEX,
        VERTEX,
        VERTEX,
        VERTEX,
        TETRA,
        PYRAMID,
        PRISM,
        VERTEX,
        HEXA
    } ;

    static CellDescriptor* cell_descriptors[7] = {
       &line_descriptor,
       &trgl_descriptor,
       &quad_descriptor,
       &tetra_descriptor,
       &pyramid_descriptor,
       &prism_descriptor,
       &hexa_descriptor
    } ;

}

#endif
