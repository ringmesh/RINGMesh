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
 *
 *
 *
 *
 *
 *     http://www.ring-team.org
 *
 *     RING Project
 *     Ecole Nationale Superieure de Geologie - GeoRessources
 *     2 Rue du Doyen Marcel Roubault - TSA 70605
 *     54518 VANDOEUVRE-LES-NANCY
 *     FRANCE
 */

#ifndef __RINGMESH_GEO_MESH__
#define __RINGMESH_GEO_MESH__

#include <ringmesh/common.h>
#include <geogram/mesh/mesh.h>

namespace RINGMesh {
    /* 
     * @brief class to encapsulate mesh structure in order to provide an API on which we base the RINGMesh algorithm 
     * @note For now, we encapsulate the GEO::Mesh class. We can develop the concept using a factory to build several encapsulating classes. 
     */
    class RINGMESH_API GeoMesh {
    ringmesh_disable_copy( GeoMesh ) ;

    public:
    GeoMesh(index_t dimension=3, bool single_precision=false) ;
    ~GeoMesh() ;

    /* 
     * @brief Gets a point.
     * @param[in] v_id the vertex, in 0.. @function nb_vetices()-1.
     * @return a modifiable reference to the point that corresponds to the vertex.
     */
    const vec3& vertex( index_t v_id ) const ;

    /* 
     * @brief Gets the number of point in the GeoMesh.
     */
    index_t nb_vetices( ) const ;

    /* 
     * @brief Gets the index of an edge vertex.
     * @param[in] edge_id index of the edge.
     * @param[in] vertex_id local index of the vertex, in {0,1} 
     * @return the global index of vertex \param vertex_id in edge \param edge_id.
     */
    index_t edge_vertex( index_t edge_id, index_t vertex_id ) const ;

    /* 
     * @return the index of the adjacent edge of the one starting at \param vertex_id
     * @TODO: implement... for now we are doing the implicit hypothese that edges are ordered with vertex.
     */
    index_t edge_adjacent ( index_t vertex_id ) const { return NO_ID; }

    /* 
     * @brief Gets the number of edges in the GeoMesh.
     */
    index_t nb_edges() const ;

    /* 
     * @brief Gets a vertex by facet and local vertex index. 
     * @param[in] facet_id the facet index.
     * @param[in] vertex_id the local vertex index in \param facet_id.
     * @return the global vertex index.
     * @precondition vertex_id < nomber of vertices of the facet.
     */
    index_t facet_vertex (index_t facet_id, index_t vertex_id) const ;

    /* 
     * @return the index of the adjacent facet of \param facet_id
     * along the edge starting at \param vertex_id
     */
    index_t facet_adjacent ( index_t facet_id, index_t vertex_id ) const ;

    /* 
     * @brief Gets the number of facets in the GeoMesh.
     */
    index_t nb_facets () const ;

    /* 
     * @brief Gets a vertex by cell and local vertex index.
     * @param[in] cell_id the cell index.
     * @param[in] vertex_id the local vertex index in \param cell_id.
     * @return the global vertex index.
     * @precondition vertex_id<nomber of vertices of the cell.
     */
    index_t cell_vertex ( index_t cell_id, index_t vertex_id ) const ;

    /* 
     * @return the index of the adjacent cell of \param cell_id along the facet \param facet_id
     */
    index_t cell_adjacent (index_t cell_id, index_t facet_id) const ;

    /* 
     * @brief Gets the number of cells in the GeoMesh.
     */
    index_t nb_cells () const ;

    private:
    GEO::Mesh* mesh_ ;

    } ;

}

#endif