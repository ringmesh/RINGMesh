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

#include <ringmesh/basic/common.h>

#include <geogram/basic/command_line.h>
#include <geogram/mesh/mesh_io.h>

#include <ringmesh/basic/geometry.h>
#include <ringmesh/geogram_extension/geogram_extension.h>
#include <ringmesh/mesh/aabb.h>

namespace RINGMesh {
    class GeoModel ;
    class MeshBaseBuilder ;
    class Mesh0DBuilder ;
    class Mesh1DBuilder ;
    class Mesh2DBuilder ;
    class Mesh3DBuilder ;
    class MeshAllDBuilder ;
}

namespace RINGMesh {

    typedef std::string MeshType ;

    /*!
     * class base class for encapsulating Mesh structure
     * @brief encapsulate adimensional mesh functionalities in order to provide an API
     * on which we base the RINGMesh algorithm
     * @note For now, we encapsulate the GEO::Mesh class.
     */
    class RINGMESH_API MeshBase: public GEO::Counted {
    ringmesh_disable_copy( MeshBase ) ;
        friend class MeshBaseBuilder ;

    public:

        virtual ~MeshBase() ;

        virtual void save_mesh( const std::string& filename ) const = 0 ;

        /*!
         * get access to GEO::MESH... only for GFX..
         * @todo Remove this function as soon as the GEO::MeshGFX is encapsulated
         */
        virtual const GEO::Mesh& gfx_mesh() const = 0 ;

        //TODO maybe reimplement the function with a RINGMesh::Mesh??
        virtual void print_mesh_bounded_attributes() const = 0 ;
        /*!
         * \name Vertex methods
         * @{
         */
        /*!
         * @brief Gets a point.
         * @param[in] v_id the vertex, in 0.. @function nb_vetices()-1.
         * @return const reference to the point that corresponds to the vertex.
         */
        virtual const vec3& vertex( index_t v_id ) const = 0 ;
        /*
         * @brief Gets the number of vertices in the Mesh.
         */
        virtual index_t nb_vertices() const = 0 ;

        virtual GEO::AttributesManager& vertex_attribute_manager() const = 0 ;

        /*!
         * @brief return the ColocaterANN at vertices
         * @warning the ColocaterANN is destroy when calling the Mesh::facets_aabb() and Mesh::cells_aabb()
         */
        const ColocaterANN& vertices_colocater_ann() const
        {
            if( vertices_ann_ == nil ) {
                std::vector< vec3 > vec_vertices( nb_vertices() ) ;
                for( index_t v = 0; v < nb_vertices(); ++v ) {
                    vec_vertices[v] = vertex( v ) ;
                }
                vertices_ann_ = new ColocaterANN( vec_vertices, true ) ;
            }
            return *vertices_ann_ ;
        }

        virtual const MeshType type_name() const = 0 ;

        virtual const std::string default_extension() const = 0 ;

        /*!
         * @}
         */
    protected:
        /*!
         * @brief MeshBase constructor.
         * @param[in] dimension dimension of the vertices.
         * @param[in] single_precision if true, vertices are stored in single precision (float),
         * else they are stored as double precision (double)..
         */
        MeshBase()
            : vertices_ann_( nil )
        {
        }

    protected:
        mutable ColocaterANN* vertices_ann_ ;
    } ;

    /*!
     * class for encapsulating isolated vertices structure
     */
    class RINGMESH_API Mesh0D: public virtual MeshBase {
    ringmesh_disable_copy( Mesh0D ) ;
        friend class Mesh0DBuilder ;

    public:
        virtual ~Mesh0D()
        {
        }

        static Mesh0D* create_mesh( const MeshType type ) ;
    protected:
        /*!
         * @brief Mesh0D constructor.
         */
        Mesh0D()
            : MeshBase()
        {
        }
    } ;
    typedef GEO::SmartPointer< Mesh0D > Mesh0D_var ;
    typedef GEO::Factory0< Mesh0D > Mesh0DFactory ;
#define ringmesh_register_mesh_0d(type) \
    geo_register_creator(RINGMesh::Mesh0DFactory, type, type::type_name_static())


    /*!
     * class for encapsulating 1D mesh component
     */
    class RINGMESH_API Mesh1D: public virtual MeshBase {
    ringmesh_disable_copy( Mesh1D ) ;
        friend class Mesh1DBuilder ;
        friend class GeogramMeshBuilder ;

    public:
        virtual ~Mesh1D()
        {
            if( edges_ann_ != nil ) delete edges_ann_ ;
            if( edges_aabb_ != nil ) delete edges_aabb_ ;
        }

        static Mesh1D* create_mesh( const MeshType type ) ;

        /*
         * @brief Gets the index of an edge vertex.
         * @param[in] edge_id index of the edge.
         * @param[in] vertex_id local index of the vertex, in {0,1}
         * @return the global index of vertex \param vertex_id in edge \param edge_id.
         */
        virtual index_t edge_vertex( index_t edge_id, index_t vertex_id ) const = 0 ;

        /*!
         * @brief Gets the number of all the edges in the whole Mesh.
         */
        virtual index_t nb_edges() const = 0 ;

        /*!
         * @brief Gets the length of the edge \param edge_id
         */
        double edge_length( index_t edge_id ) const
        {
            const vec3& e0 = vertex( edge_vertex( edge_id, 0 ) ) ;
            const vec3& e1 = vertex( edge_vertex( edge_id, 1 ) ) ;
            return ( e1 - e0 ).length() ;
        }

        vec3 edge_barycenter( index_t edge_id ) const
        {
            const vec3& e0 = vertex( edge_vertex( edge_id, 0 ) ) ;
            const vec3& e1 = vertex( edge_vertex( edge_id, 1 ) ) ;
            return ( e1 + e0 ) / 2. ;
        }
        
        /*!
         * @brief return the ColocaterANN at edges
         * @warning the ColocaterANN is destroy when calling the Mesh::facets_aabb() and Mesh::cells_aabb()
         */
        const ColocaterANN& edges_colocater_ann() const
        {
            if( edges_ann_ == nil ) {
                std::vector< vec3 > edge_centers( nb_edges() ) ;
                for( index_t e = 0; e < nb_edges(); ++e ) {
                    edge_centers[e] = edge_barycenter( e ) ;
                }
                edges_ann_ = new ColocaterANN( edge_centers, true ) ;
            }
            return *edges_ann_ ;
        }
        /*!
         * @brief Creates an AABB tree for a Mesh edges
         */
        const AABBTree1D& edges_aabb() const
        {
            if( edges_aabb_ == nil ) {
                edges_aabb_ = new AABBTree1D( *this ) ;
            }
            return *edges_aabb_ ;
        }

        virtual GEO::AttributesManager& edge_attribute_manager() const = 0 ;
    protected:
        Mesh1D()
            : MeshBase(), edges_ann_( nil ), edges_aabb_( nil )
        {
        }

    protected:
        mutable ColocaterANN* edges_ann_ ;
        mutable AABBTree1D* edges_aabb_ ;
    } ;
    typedef GEO::SmartPointer< Mesh1D > Mesh1D_var ;
    typedef GEO::Factory0< Mesh1D > Mesh1DFactory ;
#define ringmesh_register_mesh_1d(type) \
    geo_register_creator(RINGMesh::Mesh1DFactory, type, type::type_name_static())

    /*!
     * class for encapsulating 2D mesh component
     */
    class RINGMESH_API Mesh2D: public virtual MeshBase {
    ringmesh_disable_copy( Mesh2D ) ;
        friend class Mesh2DBuilder ;

    public:
        virtual ~Mesh2D()
        {
            if( facets_ann_ != nil ) delete facets_ann_ ;
            if( facets_aabb_ != nil ) delete facets_aabb_ ;
        }

        static Mesh2D* create_mesh( const MeshType type ) ;

        /*!
         * @brief Gets the vertex index by facet index and local vertex index.
         * @param[in] facet_id the facet index.
         * @param[in] vertex_id the local edge index in \param facet_id.
         * @return the global facet index adjacent to the \param edge_id of the facet \param facet_id.
         * @precondition  \param edge_id < number of edge of the facet \param facet_id .
         */
        virtual index_t facet_vertex( index_t facet_id, index_t vertex_id ) const = 0 ;
        /*!
         * @brief Gets the number of all facets in the whole Mesh.
         */
        virtual index_t nb_facets() const = 0 ;
        /*!
         * @brief Gets the number of vertices in the facet \param facet_id.
         * @param[in] facet_id facet index
         */
        virtual index_t nb_facet_vertices( index_t facet_id ) const = 0 ;
        /*!
         * @brief Gets the next vertex index in the facet \param facet_id.
         * @param[in] facet_id facet index
         * @param[in] vertex_id current index
         */
        index_t next_facet_vertex( index_t facet_id, index_t vertex_id ) const
        {
            ringmesh_assert( vertex_id < nb_facet_vertices( facet_id ) ) ;
            if( vertex_id != nb_facet_vertices( facet_id ) - 1 ) {
                return vertex_id + 1 ;
            } else {
                return 0 ;
            }
        }
        /*!
         * @brief Gets the previous vertex index in the facet \param facet_id.
         * @param[in] facet_id facet index
         * @param[in] vertex_id current index
         */
        index_t prev_facet_vertex( index_t facet_id, index_t vertex_id ) const
        {
            ringmesh_assert( vertex_id < nb_facet_vertices( facet_id ) ) ;
            if( vertex_id > 0 ) {
                return vertex_id - 1 ;
            } else {
                return nb_facet_vertices( facet_id ) - 1 ;
            }
        }
        /*!
         * @brief Gets an adjacent facet index by facet index and local edge index.
         * @param[in] facet_id the facet index.
         * @param[in] edge_id the local edge index in \param facet_id.
         * @return the global facet index adjacent to the \param edge_id of the facet \param facet_id.
         * @precondition  \param edge_id < number of edge of the facet \param facet_id .
         */
        virtual index_t facet_adjacent(
            index_t facet_id,
            index_t edge_id ) const = 0 ;
        virtual GEO::AttributesManager& facet_attribute_manager() const = 0 ;
        /*!
         * @brief Tests whether all the facets are triangles. when all the facets are triangles, storage and access is optimized.
         * @return True if all facets are triangles and False otherwise.
         */
        virtual bool facets_are_simplicies() const = 0 ;
        /*!
         * return true if the facet \param facet_id is a triangle
         */
        bool is_triangle( index_t facet_id ) const
        {
            return nb_facet_vertices( facet_id ) == 3 ;
        }

        /*!
         * Computes the Mesh facet normal
         * @param[in] facet_id the facet index
         * @return the facet normal
         */
        vec3 facet_normal( index_t facet_id ) const
        {
            const vec3& p1 = vertex( facet_vertex( facet_id, 0 ) ) ;
            const vec3& p2 = vertex( facet_vertex( facet_id, 1 ) ) ;
            const vec3& p3 = vertex( facet_vertex( facet_id, 2 ) ) ;
            vec3 norm = cross( p2 - p1, p3 - p1 ) ;
            return normalize( norm ) ;
        }
        /*!
         * Computes the Mesh facet barycenter
         * @param[in] facet_id the facet index
         * @return the facet center
         */
        vec3 facet_barycenter( index_t facet_id ) const
        {
            vec3 result( 0.0, 0.0, 0.0 ) ;
            double count = 0.0 ;
            for( index_t v = 0; v < nb_facet_vertices( facet_id ); ++v ) {
                result += vertex( facet_vertex( facet_id, v ) ) ;
                count += 1.0 ;
            }
            return ( 1.0 / count ) * result ;
        }
        /*!
         * Computes the Mesh facet area
         * @param[in] facet_id the facet index
         * @return the facet area
         */
        double facet_area( index_t facet_id ) const
        {
            double result = 0.0 ;
            if( nb_facet_vertices( facet_id ) == 0 ) {
                return result ;
            }
            const vec3& p1 = vertex( facet_vertex( facet_id, 0 ) ) ;
            for( index_t i = 1; i + 1 < nb_facet_vertices( facet_id ); i++ ) {
                const vec3& p2 = vertex( facet_vertex( facet_id, i ) ) ;
                const vec3& p3 = vertex( facet_vertex( facet_id, i + 1 ) ) ;
                result += 0.5 * length( cross( p2 - p1, p3 - p1 ) ) ;
            }
            return result ;
        }

        
        /*!
         * @brief return the ColocaterANN at facets
         */
        const ColocaterANN& facets_colocater_ann() const
        {
            if( facets_ann_ == nil ) {
                std::vector< vec3 > facet_centers( nb_facets() ) ;
                for( index_t f = 0; f < nb_facets(); ++f ) {
                    facet_centers[f] = facet_barycenter( f ) ;
                }
                facets_ann_ = new ColocaterANN( facet_centers, true ) ;
            }
            return *facets_ann_ ;
        }
        /*!
         * @brief Creates an AABB tree for a Mesh facets
         */
        const AABBTree2D& facets_aabb() const
        {
            if( facets_aabb_ == nil ) {
                facets_aabb_ = new AABBTree2D( *this ) ;
            }
            return *facets_aabb_ ;
        }
    protected:
        Mesh2D()
            : MeshBase(), facets_ann_( nil ), facets_aabb_( nil )
        {
        }

    protected:
        mutable ColocaterANN* facets_ann_ ;
        mutable AABBTree2D* facets_aabb_ ;
    } ;
    typedef GEO::SmartPointer< Mesh2D > Mesh2D_var ;
    typedef GEO::Factory0< Mesh2D > Mesh2DFactory ;
#define ringmesh_register_mesh_2d(type) \
    geo_register_creator(RINGMesh::Mesh2DFactory, type, type::type_name_static())

    /*!
     * class for encapsulating 3D mesh component
     */
    class RINGMESH_API Mesh3D: public virtual MeshBase {
    ringmesh_disable_copy( Mesh3D ) ;
        friend class Mesh3DBuilder ;

    public:
        virtual ~Mesh3D()
        {
            if( cell_facets_ann_ != nil ) delete cell_facets_ann_ ;
            if( cell_ann_ != nil ) delete cell_ann_ ;
            if( cell_aabb_ != nil ) delete cell_aabb_ ;
        }

        static Mesh3D* create_mesh( const MeshType type ) ;

        /*!
         * @brief Gets a vertex index by cell and local vertex index.
         * @param[in] cell_id the cell index.
         * @param[in] vertex_id the local vertex index in \param cell_id.
         * @return the global vertex index.
         * @precondition vertex_id<number of vertices of the cell.
         */
        virtual index_t cell_vertex( index_t cell_id, index_t vertex_id ) const = 0 ;

        /*!
         * @brief Gets a vertex index by cell and local edge and local vertex index.
         * @param[in] cell_id the cell index.
         * @param[in] edge_id the local edge index in \param cell_id.
         * @param[in] vertex_id the local vertex index in \param cell_id.
         * @return the global vertex index.
         * @precondition vertex_id<number of vertices of the cell.
         */
        virtual index_t cell_edge_vertex(
            index_t cell_id,
            index_t edge_id,
            index_t vertex_id ) const = 0 ;

        /*!
         * @brief Gets a vertex by cell facet and local vertex index.
         * @param[in] cell_id index of the cell
         * @param[in] facet_id index of the facet in the cell \param cell_id
         * @param[in] vertex_id index of the vertex in the facet \param facet_id
         * @return the global vertex index.
         * @precondition vertex_id < number of vertices in the facet \param facet_id
         * and facet_id number of facet in th cell \param cell_id
         */
        virtual index_t cell_facet_vertex(
            index_t cell_id,
            index_t facet_id,
            index_t vertex_id ) const = 0 ;

        /*!
         * @brief Gets a facet index by cell and local facet index.
         * @param[in] cell_id index of the cell
         * @param[in] facet_id index of the facet in the cell \param cell_id
         * @return the global facet index.
         */
        virtual index_t cell_facet( index_t cell_id, index_t facet_id ) const = 0 ;

        /*!
         * @brief Gets the number of facet in a cell
         * @param[in] cell_id index of the cell
         * @return the number of facet of the cell \param cell_id
         */
        virtual index_t nb_cell_facets( index_t cell_id ) const = 0 ;
        /*!
         * @brief Gets the total number of facet in a all cells
         */
        virtual index_t nb_cell_facets() const = 0 ;

        /*!
         * @brief Gets the number of edges in a cell
         * @param[in] cell_id index of the cell
         * @return the number of facet of the cell \param cell_id
         */
        virtual index_t nb_cell_edges( index_t cell_id ) const = 0 ;

        /*!
         * @brief Gets the number of vertices of a facet in a cell
         * @param[in] cell_id index of the cell
         * @param[in] facet_id index of the facet in the cell \param cell_id
         * @return the number of vertices in the facet \param facet_id in the cell \param cell_id
         */
        virtual index_t nb_cell_facet_vertices(
            index_t cell_id,
            index_t facet_id ) const = 0 ;

        /*!
         * @brief Gets the number of vertices of a cell
         * @param[in] cell_id index of the cell
         * @return the number of vertices in the cell \param cell_id
         */
        virtual index_t nb_cell_vertices( index_t cell_id ) const = 0 ;

        /*!
         * @brief Gets the number of cells in the Mesh.
         */
        virtual index_t nb_cells() const = 0 ;

        virtual index_t cell_begin( index_t cell_id ) const = 0 ;

        virtual index_t cell_end( index_t cell_id ) const = 0 ;

        /*!
         * @return the index of the adjacent cell of \param cell_id along the facet \param facet_id
         */
        virtual index_t cell_adjacent(
            index_t cell_id,
            index_t facet_id ) const = 0 ;

        virtual GEO::AttributesManager& cell_attribute_manager() const = 0 ;

        virtual GEO::AttributesManager& cell_facet_attribute_manager() const = 0 ;

        /*!
         * @brief Gets the type of a cell.
         * @param[in] cell_id the cell index, in 0..nb()-1
         */
        virtual GEO::MeshCellType cell_type( index_t cell_id ) const = 0 ;

        /*!
         * @brief Tests whether all the cells are tetrahedra. when all the cells are tetrahedra, storage and access is optimized.
         * @return True if all cells are tetrahedra and False otherwise.
         */
        virtual bool cells_are_simplicies() const = 0 ;

        /*!
         * Computes the Mesh cell facet barycenter
         * @param[in] cell_id the cell index
         * @param[in] facet_id the facet index in the cell
         * @return the cell facet center
         */
        vec3 cell_facet_barycenter( index_t cell_id, index_t facet_id ) const
        {
            vec3 result( 0., 0., 0. ) ;
            index_t nb_vertices = nb_cell_facet_vertices( cell_id, facet_id ) ;
            for( index_t v = 0; v < nb_vertices; ++v ) {
                result += vertex( cell_facet_vertex( cell_id, facet_id, v ) ) ;
            }
            ringmesh_assert( nb_vertices > 0 ) ;

            return result / nb_vertices ;
        }
        /*!
         * Compute the non weighted barycenter of the \param cell_id
         */
        vec3 cell_barycenter( index_t cell_id ) const
        {
            vec3 result( 0.0, 0.0, 0.0 ) ;
            double count = 0.0 ;
            for( index_t v = 0; v < nb_cell_vertices( cell_id ); ++v ) {
                result += vertex( cell_vertex( cell_id, v ) ) ;
                count += 1.0 ;
            }
            return ( 1.0 / count ) * result ;
        }
        /*!
         * Computes the Mesh cell facet normal
         * @param[in] cell_id the cell index
         * @param[in] facet_id the facet index in the cell
         * @return the cell facet normal
         */
        vec3 cell_facet_normal( index_t cell_id, index_t facet_id ) const
        {
            ringmesh_assert( cell_id < nb_cells() ) ;
            ringmesh_assert( facet_id < nb_cell_facets( cell_id ) ) ;

            const vec3& p1 = vertex( cell_facet_vertex( cell_id, facet_id, 0 ) ) ;
            const vec3& p2 = vertex( cell_facet_vertex( cell_id, facet_id, 1 ) ) ;
            const vec3& p3 = vertex( cell_facet_vertex( cell_id, facet_id, 2 ) ) ;

            return cross( p2 - p1, p3 - p1 ) ;
        }

        /*!
         * @brief compute the volume of the cell \param cell_id.
         */
        virtual double cell_volume( index_t cell_id ) const = 0 ;

        index_t find_cell_corner( index_t cell_id, index_t vertex_id ) const
        {
            for( index_t v = 0; v < nb_cell_vertices( cell_id ); ++v ) {
                if( cell_vertex( cell_id, v ) == vertex_id ) {
                    return cell_vertex( cell_id, v ) ;
                }
            }
            return NO_ID ;
        }

        /*!
         * @brief return the ColocaterANN at cell facets
         * @warning the ColocaterANN is destroy when calling the Mesh::facets_aabb() and Mesh::cells_aabb()
         */
        const ColocaterANN& cell_facets_colocater_ann() const
        {
            if( cell_facets_ann_ == nil ) {
                std::vector< vec3 > cell_facet_centers( nb_cell_facets() ) ;
                index_t cf = 0 ;
                for( index_t c = 0; c < nb_cells(); ++c ) {
                    for( index_t f = 0; f < nb_cell_facets( c ); ++f ) {
                        cell_facet_centers[cf] = cell_facet_barycenter( c, f ) ;
                        ++cf ;
                    }
                }
                cell_facets_ann_ = new ColocaterANN( cell_facet_centers, true ) ;
            }
            return *cell_facets_ann_ ;
        }
        /*!
         * @brief return the ColocaterANN at cells
         */
        const ColocaterANN& cells_colocater_ann() const
        {
            if( cell_ann_ == nil ) {
                std::vector< vec3 > cell_centers( nb_cells() ) ;
                for( index_t c = 0; c < nb_cells(); ++c ) {
                    cell_centers[c] = cell_barycenter( c ) ;
                }
                cell_ann_ = new ColocaterANN( cell_centers, true ) ;
            }
            return *cell_ann_ ;
        }
        /*!
         * @brief Creates an AABB tree for a Mesh cells
         */
        const AABBTree3D& cells_aabb() const
        {
            if( cell_aabb_ == nil ) {
                cell_aabb_ = new AABBTree3D( *this ) ;
            }
            return *cell_aabb_ ;
        }
    protected:
        Mesh3D()
            :
                MeshBase(),
                cell_facets_ann_( nil ),
                cell_ann_( nil ),
                cell_aabb_( nil )
        {
        }

    protected:
        mutable ColocaterANN* cell_facets_ann_ ;
        mutable ColocaterANN* cell_ann_ ;
        mutable AABBTree3D* cell_aabb_ ;
    } ;
    typedef GEO::SmartPointer< Mesh3D > Mesh3D_var ;
    typedef GEO::Factory0< Mesh3D > Mesh3DFactory ;
#define ringmesh_register_mesh_3d(type) \
    geo_register_creator(RINGMesh::Mesh3DFactory, type, type::type_name_static())

    class RINGMESH_API MeshAllD: public virtual Mesh0D,
        public virtual Mesh1D,
        public virtual Mesh2D,
        public virtual Mesh3D {
    ringmesh_disable_copy( MeshAllD ) ;
        friend class MeshAllDBuilder ;

    public:
        virtual ~MeshAllD()
        {
        }

        static MeshAllD* create_mesh( const MeshType type ) ;
    protected:
        MeshAllD()
            : Mesh0D(), Mesh1D(), Mesh2D(), Mesh3D()
        {
        }
    } ;
    typedef GEO::SmartPointer< MeshAllD > MeshAllD_var ;
    typedef GEO::Factory0< MeshAllD > MeshAllDFactory ;
#define ringmesh_register_mesh_alld(type) \
    geo_register_creator(RINGMesh::MeshAllDFactory, type, type::type_name_static())
}

#endif
