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

#include <grgmesh/aabb_tree.h>
#include <grgmesh/permutation.h>
#include <grgmesh/reorder.h>

#include <fstream>

namespace GRGMesh {

    static index_t max_node_index(
        index_t node_index,
        index_t b,
        index_t e )
    {
        grgmesh_debug_assert( e > b ) ;
        if( b + 1 == e ) {
            return node_index ;
        }
        index_t m = b + ( e - b ) / 2 ;
        index_t childl = 2 * node_index ;
        index_t childr = 2 * node_index + 1 ;
        return std::max( max_node_index( childl, b, m ),
            max_node_index( childr, m, e ) ) ;
    }
    void FacetAABBTree::init_bboxes_recursive(
        index_t node_index,
        index_t b,
        index_t e )
    {
        grgmesh_debug_assert( node_index < bboxes_.size() ) ;
        grgmesh_debug_assert( b != e ) ;
        if( b + 1 == e ) {
            Box3d bbox ;
            bbox.add_point( mesh_.vertex( b, 0 ) ) ;
            bbox.add_point( mesh_.vertex( b, 1 ) ) ;
            bbox.add_point( mesh_.vertex( b, 2 ) ) ;
            bboxes_[node_index] = bbox ;
            return ;
        }
        index_t m = b + ( e - b ) / 2 ;
        index_t childl = 2 * node_index ;
        index_t childr = 2 * node_index + 1 ;
        grgmesh_debug_assert( childl < bboxes_.size() ) ;
        grgmesh_debug_assert( childr < bboxes_.size() ) ;
        init_bboxes_recursive( childl, b, m ) ;
        init_bboxes_recursive( childr, m, e ) ;
        grgmesh_debug_assert( childl < bboxes_.size() ) ;
        grgmesh_debug_assert( childr < bboxes_.size() ) ;
        bboxes_[node_index] = bboxes_[childl].bbox_union( bboxes_[childr] ) ;
    }
    void FacetAABBTree::reorder_morton()
    {
        std::vector< int32 > sorted_indices ;
        SurfaceMutator mutator( mesh_ ) ;

        // Step 1: reorder vertices
        morton_vertex_sort( mesh_, sorted_indices ) ;

        Permutation::apply( mutator.vertices(), sorted_indices ) ;

        Surface::VertexAttributeManager* vertex_manager =
            mesh_.vertex_attribute_manager() ;
        std::vector< std::string > vertex_attribute_names ;
        vertex_manager->list_named_attributes( vertex_attribute_names ) ;
        for( index_t i = 0; i < vertex_attribute_names.size(); i++ ) {
            AttributeStore_var attribute =
                vertex_manager->resolve_named_attribute_store(
                    vertex_attribute_names[i] ) ;
            Permutation::apply( attribute->data( 0 ), sorted_indices,
                attribute->item_size() ) ;
        }

        Permutation::invert(sorted_indices) ;
        std::vector< index_t >& facets = mutator.facets() ;
        for( index_t t = 0; t < mesh_.nb_cells(); t++ ) {
            for( index_t p = 0; p < 3; p++ ) {
                facets[3*t+p] = sorted_indices[facets[3*t+p]] ;
            }
        }

        sorted_indices.clear() ;

        // Step 2: reorder facets
        morton_cell_sort( mesh_, sorted_indices ) ;

        Permutation::apply( &mutator.facets()[0], sorted_indices, sizeof(index_t) * 3 ) ;
        Permutation::apply( &mutator.adjacents()[0], sorted_indices, sizeof(signed_index_t) * 3 ) ;
        Surface::FacetAttributeManager* facet_manager =
            mesh_.facet_attribute_manager() ;
        std::vector< std::string > facet_attribute_names ;
        facet_manager->list_named_attributes( facet_attribute_names ) ;
        for( index_t i = 0; i < facet_attribute_names.size(); i++ ) {
            AttributeStore_var attribute =
                facet_manager->resolve_named_attribute_store(
                    facet_attribute_names[i] ) ;
            Permutation::apply( attribute->data( 0 ), sorted_indices,
                attribute->item_size() ) ;
        }

        Permutation::invert(sorted_indices) ;
        std::vector< index_t >& adjacents = mutator.adjacents() ;
        for( index_t t = 0; t < mesh_.nb_cells(); t++ ) {
            for( index_t p = 0; p < 3; p++ ) {
                if( !mesh_.is_on_border(t,p) ) {
                    adjacents[3*t+p] = sorted_indices[adjacents[3*t+p]] ;
                }
            }
        }
    }

    float64 FacetAABBTree::get_nearest_point(
        const Surface& M,
        const vec3& p,
        int t,
        vec3& nearest_p )
    {
        const vec3& p1 = M.vertex( t, 0 ) ;
        const vec3& p2 = M.vertex( t, 1 ) ;
        const vec3& p3 = M.vertex( t, 2 ) ;
        float64 distance = Utils::point_triangle_distance( p, p1, p2, p3,
            nearest_p ) ;
        if( Utils::point_inside_triangle( p, p1, p2, p3 ) ) {
            nearest_p = p ;
            distance = -distance ;
        }

        return distance ;
    }

    FacetAABBTree::FacetAABBTree( Surface& M )
        : mesh_( M )
    {
        reorder_morton() ;
        bboxes_.resize( max_node_index( 1, 0, mesh_.nb_cells() ) + 1 ) ; // <-- this is because size == max_index + 1 !!!
        init_bboxes_recursive( 1, 0, mesh_.nb_cells() ) ;
    }

    void FacetAABBTree::get_nearest_facet_hint(
        const vec3& p,
        index_t& nearest_t,
        vec3& nearest_point,
        float64& sq_dist ) const
    {

        // Find a good initial value for nearest_f by traversing
        // the boxes and selecting the child such that the center
        // of its bounding box is nearer to the query point.
        // For a large mesh (20M facets) this gains up to 10%
        // performance as compared to picking nearest_f randomly.
        uint32 b = 0 ;
        uint32 e = mesh_.nb_cells() - 1 ;
        if( e > 0 ) {
            index_t n = 1 ;
            while( e != b + 1 ) {
                index_t m = b + ( e - b ) / 2 ;
                index_t childl = 2 * n ;
                index_t childr = 2 * n + 1 ;
                if( bboxes_[childl].distance_to_center( p )
                    < bboxes_[childr].distance_to_center( p ) ) {
                    e = m ;
                    n = childl ;
                } else {
                    b = m ;
                    n = childr ;
                }
            }
        }
        nearest_t = b ;

        nearest_point = mesh_.vertex( nearest_t, 0 ) ;
        sq_dist = length2( p - nearest_point ) ;
    }

    void FacetAABBTree::nearest_facet_recursive(
        const vec3& p,
        index_t& nearest_t,
        vec3& nearest_point,
        float64& sq_dist,
        index_t n,
        index_t b,
        index_t e ) const
    {
        grgmesh_debug_assert( e > b ) ;

        // If node is a leaf: compute point-facet distance
        // and replace current if nearer
        if( b + 1 == e ) {
            vec3 cur_nearest_point ;
            float64 cur_sq_dist = get_nearest_point( mesh_, p, b, cur_nearest_point ) ;
            if( cur_sq_dist < sq_dist ) {
                nearest_t = b ;
                nearest_point = cur_nearest_point ;
                sq_dist = cur_sq_dist ;
            }
            return ;
        }
        index_t m = b + ( e - b ) / 2 ;
        index_t childl = 2 * n ;
        index_t childr = 2 * n + 1 ;

        float64 dl = bboxes_[childl].signed_distance( p ) ;
        float64 dr = bboxes_[childr].signed_distance( p ) ;

        // Traverse the "nearest" child first, so that it has more chances
        // to prune the traversal of the other child.
        if( dl < dr ) {
            if( dl < sq_dist ) {
                nearest_facet_recursive( p, nearest_t, nearest_point, sq_dist,
                    childl, b, m ) ;
            }
            if( dr < sq_dist ) {
                nearest_facet_recursive( p, nearest_t, nearest_point, sq_dist,
                    childr, m, e ) ;
            }
        } else {
            if( dr < sq_dist ) {
                nearest_facet_recursive( p, nearest_t, nearest_point, sq_dist,
                    childr, m, e ) ;
            }
            if( dl < sq_dist ) {
                nearest_facet_recursive( p, nearest_t, nearest_point, sq_dist,
                    childl, b, m ) ;
            }
        }
    }

}
