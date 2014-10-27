/*[
* Association Scientifique pour la Geologie et ses Applications (ASGA)
* Copyright (c) 1993-2013 ASGA. All Rights Reserved.
*
* This program is a Trade Secret of the ASGA and it is not to be:
* - reproduced, published, or disclosed to other,
* - distributed or displayed,
* - used for purposes or on Sites other than described
*   in the GOCAD Advancement Agreement,
* without the prior written authorization of the ASGA. Licencee
* agrees to attach or embed this Notice on all copies of the program,
* including partial copies or modified versions thereof.
]*/

#ifndef __GRGMESH_TYPES__
#define __GRGMESH_TYPES__

#include <grgmeshlib/common.h>
#include <stdint.h>

namespace GRGMesh {

    /** Generic pointer type */
    typedef void* pointer ;
    /** Small pointer type */
    typedef uint8_t* small_pointer ;

    /** Null pointer type */
    #define nil 0

    /** Integer type with a width of 8 bits */
    typedef int8_t int8 ;
    static int8 dummy_int8 ;

    /** Integer type with a width of 16 bits */
    typedef int16_t int16 ;
    static int16 dummy_int16 ;

    /** Integer type with a width of 32 bits */
    typedef int32_t int32 ;
    static int32 dummy_int32 ;

    /** Integer type with a width of 64 bits */
    typedef int64_t int64 ;
    static int64 dummy_int64 ;

    /** Unsigned integer type with a width of 8 bits */
    typedef uint8_t uint8 ;
    static uint8 dummy_uint8 ;

    /** Unsigned integer type with a width of 16 bits */
    typedef uint16_t uint16 ;
    static uint16 dummy_uint16 ;

    /** Unsigned integer type with a width of 32 bits */
    typedef uint32_t uint32 ;
    static uint32 dummy_uint32 ;

    /** Unsigned integer type with a width of 64 bits */
    typedef uint64_t uint64 ;
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

    enum Tetra_method {
        TetGen, MG_Tetra
    } ;

    enum CellType {
        LINE    = 0,
        TRGL    = 1,
        QUAD    = 2,
        TETRA   = 3,
        PYRAMID = 4,
        PRISM   = 5,
        HEXA    = 6
    } ;

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

}

#endif
