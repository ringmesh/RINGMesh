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

#ifndef __GRGMESH_MIXED_MESH__
#define __GRGMESH_MIXED_MESH__

#include <common.h>
#include <mesh.h>

#include <vector> 

namespace GRGMesh {

    #define MixedMesh_TEMPLATE template< int DIM >
    #define MixedMesh_DECL MixedMesh< int DIM >

    static std::vector< int > vector_init = std::vector< int >( 12, -1 ) ;
    /*!
     * Mesh which can handle different type of elements.
     */
    template< int DIM=3 > class GRGMESH_API MixedMesh: public Mesh< DIM > {

    public:
        MixedMesh()
        {
        }
        ~MixedMesh()
        {
        }

        unsigned int nb_vertices_in_cell( unsigned int c ) const
        {
            grgmesh_debug_assert( c < nb_c() );
            return cell_end( c ) - cell_begin( c ) ;
        }


        CellType cell_type( unsigned int c ) const
        {
            int result = nb_v_to_cell_type[nb_vertices_in_cell( c )] ;
            grgmesh_debug_assert( result >= 0 );
            return CellType( result ) ;
        }

        CellDescriptor* cell_descriptor( unsigned int c ) const
        {
            CellDescriptor* result =
                nb_v_to_cell_descriptor[nb_vertices_in_cell( c )] ;
            grgmesh_debug_assert( result != 0 );
            return result ;
        }

        unsigned int nb_facets_in_cell( unsigned int c ) const
        { return cell_descriptor( c )->nb_facets ; }

        unsigned int nb_vertices_in_cell_facet(
            unsigned int c,
            unsigned int f ) const
        {
            CellDescriptor* desc = cell_descriptor( c ) ;
            grgmesh_debug_assert( f < desc->nb_facets );
            return desc->nb_vertices_in_facet[f] ;
        }

        unsigned int cell_facet_vertex(
            unsigned int c,
            unsigned int f,
            unsigned int v ) const
        {
            CellDescriptor* desc = cell_descriptor( c ) ;
            grgmesh_debug_assert( f < desc->nb_facets );
            grgmesh_debug_assert( v < desc->nb_vertices_in_facet[f] );
            return cell_vertex( c, desc->facet[f][v] ) ;
        }

    private:
        void copy( const MixedMesh& rhs )
        {
        }

    private:
        MixedMesh( const MixedMesh& rhs ) ;
        MixedMesh& operator=( const MixedMesh& rhs ) ;

        std::vector< uint64 > line_ ;
        std::vector< uint64 > trgl_ ;
        std::vector< uint64 > quad_ ;
        std::vector< uint64 > tetra_ ;
        std::vector< uint64 > pyramid_ ;
        std::vector< uint64 > prism_ ;
        std::vector< uint64 > hexa_ ;


    public:
        static int nb_v_to_cell_type[9] ;
        static CellDescriptor* nb_v_to_cell_descriptor[9] ;
    } ;

    MixedMesh_TEMPLATE
    static MixedMesh_DECL::CellDescriptor tet_descriptor = {
        4,              // nb_vertices
        4,              // nb_facets
        { 3, 3, 3, 3 }, // nb_vertices in facet
        { { 3, 1, 2 }, { 3, 2, 0 }, { 0, 1, 3 }, { 1, 0, 2 } },  // facets
        6,              // nb_edges
        { { 1, 3 }, { 3, 2 }, { 2, 1 }, { 0, 2 }, { 3, 0 }, { 1, 0 } } , // edges
        { { 2, 0, 1 }, { 1, 4, 3 }, { 4, 0, 5 }, { 3, 2, 5 } } //edges_in_facet
    } ;

    MixedMesh_TEMPLATE
    static MixedMesh_DECL::CellDescriptor hex_descriptor = {
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

    MixedMesh_TEMPLATE
    static MixedMesh_DECL::CellDescriptor prism_descriptor = {
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

    MixedMesh_TEMPLATE
    static MixedMesh_DECL::CellDescriptor pyramid_descriptor = {
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

    MixedMesh_TEMPLATE
    MixedMesh_DECL::CellDescriptor* MixedMesh_DECL::nb_v_to_cell_descriptor[9] = {
        0, 0, 0, 0, &tet_descriptor, &pyramid_descriptor, &prism_descriptor, 0,
        &hex_descriptor } ;

}

#endif
