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

#ifndef __RINGMESH_MESH__
#define __RINGMESH_MESH__

#include <ringmesh/common.h>
#include <ringmesh/geo_model_element.h>

#include <geogram/mesh/mesh.h>
#include <geogram/mesh/mesh_io.h>

namespace GEO {
    class MeshFacetsAABB ;
    class MeshCellsAABB ;
}

namespace RINGMesh {

    /* 
     * @brief class to encapsulate mesh structure in order to provide an API on which we base the RINGMesh algorithm 
     * @note For now, we encapsulate the GEO::Mesh class. We can develop the concept using a factory to build several encapsulating classes. 
     */
    class RINGMESH_API Mesh {
    ringmesh_disable_copy( Mesh ) ;
        friend class MeshBuilder ;
    public:

        /*
         * @brief Mesh constructor.
         * @param[in] dimension dimension of the vertices.
         * @parm[in] single_precision if true, vertices are stored in single precision (float), else they are stored as double precision (double)..
         */
        Mesh(
            const GeoModelElement& geo_model_elment,
            index_t dimension,
            bool single_precision )
            :
                geo_model_elment_( geo_model_elment ),
                facets_aabb_( nil ),
                cells_aabb_( nil )
        {
            mesh_ = new GEO::Mesh( dimension, single_precision ) ;
        }
        ~Mesh()
        {
            delete mesh_ ;
        }

        /*
         * @brief Copy a mesh into this one.
         * @param[in] rhs a const reference to the mesh to be copied.
         * @parm[in] copy_attributes if true, all attributes are copied.
         * @param[in] what a combination of MESH_VERTICES, MESH_EDGES, MESH_FACETS, MESH_CELLS flags. Set to MESH_ALL_ELEMENTS to copy everything (default). If MESH_VERTICES is not set, then the mesh is cleared.
         * @return a modifiable reference to the point that corresponds to the vertex.
         */
        void copy(
            const Mesh& rhs,
            bool copy_attributes,
            GEO::MeshElementsFlags what )
        {
            mesh_->copy( *rhs.mesh_, copy_attributes, what ) ;
        }

        void save_mesh( const std::string filename )
        {
            GEO::mesh_save( *mesh_, filename ) ;
        }

        GEO::AttributesManager& vertex_attribute_manager() const
        {
            return mesh_->vertices.attributes() ;
        }
        GEO::AttributesManager& facet_attribute_manager() const
        {
            return mesh_->facets.attributes() ;
        }
        GEO::AttributesManager& cell_attribute_manager() const
        {
            return mesh_->cells.attributes() ;
        }

        /*!
         * @brief Create an AABB tree for a Mesh facets
         * @pre The GeoModelElement must be simplicial
         * @warning SIDE EFFECTS: The mesh vertices are reordered.
         * That is why the global Mesh vertices are deleted (This is BAD)
         */

        const GEO::MeshFacetsAABB& facets_aabb() const ;

        /*!
         * @brief Create an AABB tree for a Mesh cells
         * @pre The GeoModelElement must be simplicial
         * @warning SIDE EFFECTS: The mesh vertices are reordered.
         * That is why the global Mesh vertices are deleted (This is BAD)
         */
        const GEO::MeshCellsAABB& cells_aabb() const ;

        /*
         * @brief Gets a dimension of points in the mesh.
         */
        index_t vertices_dimension() const
        {
            return mesh_->vertices.dimension() ;
        }
        /*
         * @brief Gets a point.
         * @param[in] v_id the vertex, in 0.. @function nb_vetices()-1.
         * @return reference to the point that corresponds to the vertex.
         */
        const vec3& vertex( index_t v_id ) const
        {
            return mesh_->vertices.point( v_id ) ;
        }

        /*
         * @brief Gets the number of point in the Mesh.
         */
        index_t nb_vertices() const
        {
            return mesh_->vertices.nb() ;
        }

        /*
         * @brief Gets the index of an edge vertex.
         * @param[in] edge_id index of the edge.
         * @param[in] vertex_id local index of the vertex, in {0,1}
         * @return the global index of vertex \param vertex_id in edge \param edge_id.
         */
        index_t edge_vertex( index_t edge_id, index_t vertex_id ) const
        {
            return mesh_->edges.vertex( edge_id, vertex_id ) ;
        }

        /*
         * @return the index of the adjacent edge of the one starting at \param vertex_id
         * @TODO: implement... for now we are doing the implicit hypothese that edges are ordered with vertex.
         */
        index_t edge_adjacent( index_t vertex_id ) const
        {
            return NO_ID ;
        }

        /*
         * @brief Gets the number of all the edges in the whole Mesh.
         */
        index_t nb_edges() const
        {
            return mesh_->edges.nb() ;
        }

        /*
         * @brief Gets the vertex index by facet index and local vertex index.
         * @param[in] facet_id the facet index.
         * @param[in] vertex_id the local edge index in \param facet_id.
         * @return the global facet index adjacent to the \param edge_id of the facet \param facet_id.
         * @precondition  \param edge_id < number of edge of the facet \param facet_id .
         */
        index_t facet_vertex( index_t facet_id, index_t vertex_id ) const
        {
            return mesh_->facets.vertex( facet_id, vertex_id ) ;
        }
        /*
         * @brief Gets an adjacent facet index by facet index and local edge index.
         * @param[in] facet_id the facet index.
         * @param[in] edge_id the local edge index in \param facet_id.
         * @return the global facet index adjacent to the \param edge_id of the facet \param facet_id.
         * @precondition  \param edge_id < number of edge of the facet \param facet_id .
         */
        index_t facet_adjacent( index_t facet_id, index_t edge_id ) const
        {
            return mesh_->facets.adjacent( facet_id, edge_id ) ;
        }
        /*
         * @brief Gets the number of vertices in the facet \param facet_id.
         * @param[in] facet_id facet index
         */
        index_t facet_nb_vertices( index_t facet_id ) const
        {
            return mesh_->facets.nb_vertices( facet_id ) ;
        }

        /*
         * @brief Gets the number of all facets in the whole Mesh.
         */
        index_t nb_facets() const
        {
            return mesh_->facets.nb() ;
        }

        /*
         * @brief Gets a vertex index by cell and local vertex index.
         * @param[in] cell_id the cell index.
         * @param[in] vertex_id the local vertex index in \param cell_id.
         * @return the global vertex index.
         * @precondition vertex_id<number of vertices of the cell.
         */
        index_t cell_vertex( index_t cell_id, index_t vertex_id ) const
        {
            return mesh_->cells.vertex( cell_id, vertex_id ) ;
        }

        /*
         * @brief Gets a vertex by cell facet and local vertex index.
         * @param[in] cell_id index of the cell
         * @param[in] facet_id index of the facet in the cell \param cell_id
         * @param[in] vertex_id index of the vertex in the facet \param facet_id
         * @return the global vertex index.
         * @precondition vertex_id < number of vertices in the facet \param facet_id and facet_id number of facet in th cell \param cell_id
         */
        index_t cell_facet_vertex(
            index_t cell_id,
            index_t facet_id,
            index_t vertex_id ) const
        {
            return mesh_->cells.facet_vertex( cell_id, facet_id, vertex_id ) ;
        }

        /*
         * @return the index of the adjacent cell of \param cell_id along the facet \param facet_id
         */
        index_t cell_adjacent( index_t cell_id, index_t facet_id ) const
        {
            return mesh_->cells.adjacent( cell_id, facet_id ) ;
        }

        /*
         * @brief Gets the number of facet in a cell
         * @param[in] cell_id index of the cell
         * @return the number of facet of the cell \param cell_id
         */
        index_t nb_cell_facets( index_t cell_id ) const
        {
            return mesh_->cells.nb_facets( cell_id ) ;
        }

        /*
         * @brief Gets the number of vertices of a facet in a cell
         * @param[in] cell_id index of the cell
         * @param[in] facet_id index of the facet in the cell \param cell_id
         * @return the number of vertices in the facet \param facet_id in the cell \param cell_id
         */
        index_t nb_cell_facet_vertices( index_t cell_id, index_t facet_id ) const
        {
            return mesh_->cells.facet_nb_vertices( cell_id, facet_id ) ;
        }

        /*
         * @brief Gets the number of vertices of a cell
         * @param[in] cell_id index of the cell
         * @return the number of vertices in the cell \param cell_id
         */
        index_t nb_cell_vertices( index_t cell_id ) const
        {
            return mesh_->cells.nb_vertices( cell_id ) ;
        }
        /*
         * @brief Gets the number of cells in the Mesh.
         */
        index_t nb_cells() const
        {
            return mesh_->cells.nb() ;
        }

        /*
         * @brief Gets the type of a cell.
         * @param[in] cell_id the cell index, in 0..nb()-1
         */
        GEO::MeshCellType cell_type( index_t cell_id ) const
        {
            return mesh_->cells.type( cell_id ) ;
        }

    private:
        GEO::Mesh* mesh_ ;
        const GeoModelElement& geo_model_elment_ ;

        GEO::MeshFacetsAABB* facets_aabb_ ;
        GEO::MeshCellsAABB* cells_aabb_ ;

    } ;

    class RINGMESH_API MeshBuilder {
    ringmesh_disable_copy( MeshBuilder ) ;

    public:
        MeshBuilder( Mesh& mesh )
            : mesh_( mesh )
        {
        }
        ~MeshBuilder() ;

        index_t create_vertex( const double* coords )
        {
            return mesh_.mesh_->vertices.create_vertex( coords ) ;
        }

    private:
        Mesh& mesh_ ;

    } ;
}

#endif
