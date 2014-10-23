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

#ifndef __GRGMESH_MESH__
#define __GRGMESH_MESH__

#include <grgmeshlib/common.h>
#include <grgmeshlib/vecn.h>

#include <vector>

namespace GRGMesh {

    class GRGMESH_API Mesh {

    public:
        void clear()
        {
            vertices_.clear() ;
            vertex_indices_.clear() ;
        }

        uint64 nb_v() const
        {
            return vertices_.size() ;
        }
        virtual uint64 nb_c() const = 0 ;

        /*!
         * Get the vertex first index of a cell
         * @param c Cell index
         */
        virtual uint64 cell_begin( uint64 c ) const = 0 ;
        /*!
         * Get the vertex last index of a cell (the first of the next cell)
         * @param c Cell index
         */
        virtual uint64 cell_end( uint64 c ) const = 0 ;
        virtual uint64 nb_vertices_in_cell( uint64 c ) const = 0 ;
        virtual uint64 nb_facets_in_cell( uint64 c ) const = 0 ;
        virtual uint64 nb_vertices_in_cell_facet( uint64 c, uint64 f ) const = 0 ;

        uint64 cell_vertex_index( uint64 c, uint64 v ) const
        {
            grgmesh_debug_assert( v < nb_vertices_in_cell( c ) ) ;
            return vertex_indices_[cell_begin( c ) + v] ;
        }
        const vec3& cell_vertex( uint64 c, uint64 v ) const
        {
            grgmesh_debug_assert( v < nb_vertices_in_cell( c ) ) ;
            return vertices_[vertex_indices_[cell_begin( c ) + v]] ;
        }


        virtual CellType cell_type( uint64 c ) const = 0 ;
        virtual CellDescriptor* cell_descriptor( uint64 c ) const = 0 ;
        uint64 cell_facet_vertex( uint64 c, uint64 f, uint64 v ) const
        {
            CellDescriptor* desc = cell_descriptor( c ) ;
            grgmesh_debug_assert( f < desc->nb_facets ) ;
            grgmesh_debug_assert( v < desc->nb_vertices_in_facet[f] ) ;
            return cell_vertex_index( c, desc->facet[f][v] ) ;
        }

        vec3 cell_centroid( uint64 c ) const
        {
            vec3 result ;
            for( uint64 i = cell_begin( c ); i < cell_end( c ); i++ ) {
                result += cell_vertex( c, i ) ;
            }
            return result / nb_vertices_in_cell( c );
        }

    protected:
        Mesh()
        {
        }
        ~Mesh()
        {
        }
        void copy( const Mesh& rhs )
        {
        }

    private:
        ///List of the vertices
        std::vector< vec3 > vertices_ ;
        ///Mapping between the list of the vertices in a cell and the actual vertices
        std::vector< uint64 > vertex_indices_ ;

    } ;

}

#endif
