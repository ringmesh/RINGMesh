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

#include <grgmeshlib/common.h>
#include <grgmeshlib/mesh.h>

#include <vector> 

namespace GRGMesh {

    /*!
     * Mesh which can handle different type of elements.
     */
    class GRGMESH_API MixedMesh: public Mesh {

    public:
        MixedMesh()
        {
        }
        virtual ~MixedMesh()
        {
        }

        uint64 nb_vertices_in_cell( uint64 c ) const
        {
            grgmesh_debug_assert( c < nb_c() ) ;
            return cell_end( c ) - cell_begin( c ) ;
        }

        CellType cell_type( uint64 c ) const ;

        CellDescriptor* cell_descriptor( uint64 c ) const
        {
            CellDescriptor* result = cell_descriptor_[cell_type( c )] ;
            grgmesh_debug_assert( result != 0 ) ;
            return result ;
        }

        uint64 nb_facets_in_cell( uint64 c ) const
        {
            return cell_descriptor( c )->nb_facets ;
        }

        uint64 nb_vertices_in_cell_facet( uint64 c, uint64 f ) const
        {
            CellDescriptor* desc = cell_descriptor( c ) ;
            grgmesh_debug_assert( f < desc->nb_facets ) ;
            return desc->nb_vertices_in_facet[f] ;
        }

        uint64 cell_facet_vertex( uint64 c, uint64 f, uint64 v ) const
        {
            CellDescriptor* desc = cell_descriptor( c ) ;
            grgmesh_debug_assert( f < desc->nb_facets ) ; grgmesh_debug_assert( v < desc->nb_vertices_in_facet[f] ) ;
            return cell_vertex_index( c, desc->facet[f][v] ) ;
        }

    private:
        void copy( const MixedMesh& rhs )
        {
        }

    private:
        MixedMesh( const MixedMesh& rhs ) ;
        MixedMesh& operator=( const MixedMesh& rhs ) ;

        std::vector< uint64 > cells_[7] ;
        std::vector< CellDescriptor* > cell_descriptor_ ;
    } ;

}

#endif
