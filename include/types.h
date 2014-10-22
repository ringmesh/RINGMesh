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

#include <stdint.h>

namespace GRGMesh {

    /** Generic pointer type */
    typedef void* pointer ;

    /** Integer type with a width of 8 bits */
    typedef int8_t int8 ;

    /** Integer type with a width of 16 bits */
    typedef int16_t int16 ;

    /** Integer type with a width of 32 bits */
    typedef int32_t int32 ;

    /** Integer type with a width of 64 bits */
    typedef int64_t int64 ;

    /** Unsigned integer type with a width of 8 bits */
    typedef uint8_t uint8 ;

    /** Unsigned integer type with a width of 16 bits */
    typedef uint16_t uint16 ;

    /** Unsigned integer type with a width of 32 bits */
    typedef uint32_t uint32 ;

    /** Unsigned integer type with a width of 64 bits */
    typedef uint64_t uint64 ;

    /** Floating point type with a width of 32 bits */
    typedef float float32 ;

    /** Floating point type with a width of 64 bits */
    typedef double float64 ;

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
}

#endif
