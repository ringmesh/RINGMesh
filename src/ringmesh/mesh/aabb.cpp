/*
 * Copyright (c) 2012-2016, Association Scientifique pour la Geologie et ses Applications (ASGA)
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *     * Redistributions in binary form must reproduce the above copyright
 *       notice, this list of conditions and the following disclaimer in the
 *       documentation and/or other materials provided with the distribution.
 *     * Neither the name of ASGA nor the
 *       names of its contributors may be used to endorse or promote products
 *       derived from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL ASGA BE LIABLE FOR ANY
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 *     http://www.ring-team.org
 *
 *     RING Project
 *     Ecole Nationale Superieure de Geologie - GeoRessources
 *     2 Rue du Doyen Marcel Roubault - TSA 70605
 *     54518 VANDOEUVRE-LES-NANCY
 *     FRANCE
 */

#include <ringmesh/mesh/aabb.h>

#include <ringmesh/mesh/mesh.h>

/// Copied and adapted from Geogram

namespace {

    using namespace RINGMesh ;

    /**
     * \brief Computes the axis-aligned bounding box of a mesh cell
     * \param[in] mesh the mesh
     * \param[in] element_id the index of the element in mesh \p mesh
     * \param[out] bbox the bounding box of the facet
     */
    void get_element_bbox( const MeshBase& mesh, index_t element_id, Box3d& bbox )
    {
        for( index_t v = 0; v < mesh.nb_mesh_element_vertices( element_id ); v++ ) {
            const vec3& point = mesh.mesh_element_vertex( element_id, v ) ;
            bbox.add_point( point ) ;
        }
    }
    
    /**
     * \brief Computes the maximum node index in a subtree
     * \param[in] node_index node index of the root of the subtree
     * \param[in] element_begin first element index in the subtree
     * \param[in] element_end one position past the last element index in the subtree
     * \return the maximum node index in the subtree rooted at \p node_index
     */
    index_t max_node_index(
        index_t node_index,
        index_t element_begin,
        index_t element_end )
    {
        ringmesh_assert( element_end > element_begin ) ;
        if( element_begin + 1 == element_end ) {
            return node_index ;
        }
        index_t element_middle = element_begin
            + ( element_end - element_begin ) / 2 ;
        index_t child_left = 2 * node_index ;
        index_t child_right = 2 * node_index + 1 ;
        return std::max( max_node_index( child_left, element_begin, element_middle ),
            max_node_index( child_right, element_middle, element_end ) ) ;
    }

    /**
     * \brief Computes the hierarchy of bounding boxes recursively.
     * \param[in] mesh the mesh
     * \param[in] bboxes the array of bounding boxes
     * \param[in] node_index the index of the root of the subtree
     * \param[in] element_begin first element index in the subtree
     * \param[in] element_end one position past the last element index in the subtree
     */
    void init_bboxes_recursive(
        const MeshBase& mesh,
        std::vector< Box3d >& bboxes,
        index_t node_index,
        index_t element_begin,
        index_t element_end )
    {
        ringmesh_assert( node_index < bboxes.size() ) ;
        ringmesh_assert( element_begin != element_end ) ;
        if( element_begin + 1 == element_end ) {
            get_element_bbox( mesh, element_begin, bboxes[node_index] ) ;
            return ;
        }
        index_t element_middle = element_begin
            + ( element_end - element_begin ) / 2 ;
        index_t child_left = 2 * node_index ;
        index_t child_right = 2 * node_index + 1 ;
        ringmesh_assert( child_left < bboxes.size() ) ;
        ringmesh_assert( child_right < bboxes.size() ) ;
        init_bboxes_recursive( mesh, bboxes, child_left, element_begin,
            element_middle ) ;
        init_bboxes_recursive( mesh, bboxes, child_right, element_middle,
            element_end ) ;
        ringmesh_assert( child_left < bboxes.size() ) ;
        ringmesh_assert( child_right < bboxes.size() ) ;
        bboxes[node_index] = bboxes[child_left].bbox_union( bboxes[child_right] ) ;
    }

}

/****************************************************************************/

namespace RINGMesh {

    AABBTree::AABBTree( const MeshBase& mesh )
    {
        bboxes_.resize( max_node_index( 1, 0, mesh.nb_mesh_elements() ) + 1 ) ;
        init_bboxes_recursive( mesh, bboxes_, 1, 0, mesh.nb_mesh_elements() ) ;
    }
        
}

