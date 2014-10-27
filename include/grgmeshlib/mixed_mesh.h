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
        friend class MixedMeshMutator ;

    public:
        MixedMesh()
        {
        }
        virtual ~MixedMesh()
        {
        }
        MixedMesh( const MixedMesh& rhs ) { copy( rhs ) ; }
        MixedMesh& operator=( const MixedMesh& rhs ) { copy( rhs ) ; return *this ; }

        virtual uint64 nb_cells() const
        {
            return nb_lines() + nb_triangles() + nb_quad() + nb_tetra()
                + nb_pyramids() + nb_prisms() + nb_hexa() ;
        }
        uint64 nb_lines() const { return std::max( uint64(0), cells_[LINE].size()-1 ) ; }
        uint64 nb_triangles() const { return std::max( uint64(0), cells_[TRGL].size()-1 ) ; }
        uint64 nb_quad() const { return std::max( uint64(0), cells_[QUAD].size()-1 ) ; }
        uint64 nb_tetra() const { return std::max( uint64(0), cells_[TETRA].size()-1 ) ; }
        uint64 nb_pyramids() const { return std::max( uint64(0), cells_[PYRAMID].size()-1 ) ; }
        uint64 nb_prisms() const { return std::max( uint64(0), cells_[PRISM].size()-1 ) ; }
        uint64 nb_hexa() const { return std::max( uint64(0), cells_[HEXA].size()-1 ) ; }
        virtual uint8 nb_vertices_in_cell( uint64 c ) const
        {
            grgmesh_debug_assert( c < nb_cells() ) ;
            return cell_descriptor( c )->nb_vertices ;
        }
        virtual uint8 nb_facets_in_cell( uint64 c ) const
        {
            grgmesh_debug_assert( c < nb_cells() ) ;
            return cell_descriptor( c )->nb_facets ;
        }
        virtual uint8 nb_vertices_in_cell_facet( uint64 c, uint8 f ) const
        {
            grgmesh_debug_assert( c < nb_cells() ) ;
            return cell_descriptor( c )->nb_vertices_in_facet[f] ;
        }
        virtual uint64 cell_begin( uint64 c ) const {
            uint64 real_index ;
            return cells_[cell_type( c, real_index)][real_index] ;
        }
        virtual uint64 cell_end( uint64 c ) const {
            uint64 real_index ;
            return cells_[cell_type( c+1, real_index)][real_index] ;
        }

        CellType cell_type( uint64 c, uint64& c_index = dummy_uint64 ) const ;

        const CellDescriptor* cell_descriptor( uint64 c ) const
        {
            const CellDescriptor* result = cell_descriptor_[cell_type( c )] ;
            grgmesh_debug_assert( result != 0 ) ;
            return result ;
        }

        uint64 cell_facet_vertex( uint64 c, uint8 f, uint8 v ) const
        {
            const CellDescriptor* desc = cell_descriptor( c ) ;
            grgmesh_debug_assert( f < desc->nb_facets ) ; grgmesh_debug_assert( v < desc->nb_vertices_in_facet[f] ) ;
            return cell_vertex_index( c, desc->facet[f][v] ) ;
        }

        uint64 tetra_vertex_index( uint64 c, uint8 f, uint8 v ) const {
            return vertex_index( cells_[TETRA][4*c + cell_descriptor_[TETRA]->facet[f][v]] ) ;
        }

    private:
        void copy( const MixedMesh& rhs )
        {
            for( uint8 i = 0; i < 7; i++) cells_[i] = rhs.cells_[i] ;
        }

    private:
        std::vector< uint64 > cells_[7] ;
        static const CellDescriptor* cell_descriptor_[7] ;
    } ;


    class GRGMESH_API MixedMeshMutator: public MeshMutator {
    public:
        MixedMeshMutator( MixedMesh& mixed_mesh )
            : MeshMutator( mixed_mesh ), mixed_mesh_( mixed_mesh )
        {
        }
        MixedMeshMutator( const MixedMesh& mixed_mesh )
            :
                MeshMutator( const_cast< MixedMesh& >( mixed_mesh ) ),
                mixed_mesh_( const_cast< MixedMesh& >( mixed_mesh ) )
        {
        }
        virtual ~MixedMeshMutator() {}

        std::vector< uint64 >* cells() { return mixed_mesh_.cells_ ; }
        std::vector< uint64 >& lines() { return mixed_mesh_.cells_[LINE] ; }
        std::vector< uint64 >& triangles() { return mixed_mesh_.cells_[TRGL] ; }
        std::vector< uint64 >& quad() { return mixed_mesh_.cells_[QUAD] ; }
        std::vector< uint64 >& tetra() { return mixed_mesh_.cells_[TETRA] ; }
        std::vector< uint64 >& pyramids() { return mixed_mesh_.cells_[PYRAMID] ; }
        std::vector< uint64 >& prisms() { return mixed_mesh_.cells_[PRISM] ; }
        std::vector< uint64 >& hexa() { return mixed_mesh_.cells_[HEXA] ; }

    private:
        MixedMesh& mixed_mesh_ ;
    } ;

}

#endif
