/*
 * Copyright (c) 2012-2015, Association Scientifique pour la Geologie et ses Applications (ASGA)
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
 *  Contacts:
 *     Arnaud.Botella@univ-lorraine.fr
 *     Antoine.Mazuyer@univ-lorraine.fr
 *     Jeanne.Pellerin@wias-berlin.de
 *
 *     http://www.ring-team.org
 *
 *     RING Project
 *     Ecole Nationale Supérieure de Géologie - Georessources
 *     2 Rue du Doyen Marcel Roubault - TSA 70605
 *     54518 VANDOEUVRE-LES-NANCY
 *     FRANCE
 */


#include <ringmesh/geo_model_mesh.h>
#include <ringmesh/geo_model.h>
#include <ringmesh/utils.h>

#include <geogram/basic/algorithm.h>
#include <geogram/mesh/mesh_repair.h>
#include <geogram/points/colocate.h>

namespace {
    using namespace RINGMesh ;

    inline GeoModelMeshElement& cast_gmm_element(
        const GeoModel& M,
        GME::TYPE T,
        index_t i )
    {
        return dynamic_cast< GeoModelMeshElement& >( const_cast< GME& >( M.element(
            GME::gme_t( T, i ) ) ) ) ;
    }
}

namespace RINGMesh {

    GeoModelMeshVertices::GeoModelMeshVertices( GeoModelMesh& gmm, GEO::Mesh& mesh )
        : gmm_( gmm ), gm_( gmm.model() ), mesh_( mesh ), kdtree_( nil )
    {
    }

    GeoModelMeshVertices::~GeoModelMeshVertices()
    {
        clear_kdtree() ;
    }

    bool GeoModelMeshVertices::is_initialized() const {
        return mesh_.vertices.nb() > 0 ;
    }

    void GeoModelMeshVertices::test_and_initialize() const
    {
        if( !is_initialized() ) {
            const_cast< GeoModelMeshVertices* >( this )->initialize() ;
        }
    }

    void GeoModelMeshVertices::test_kdtree_and_initialize() const
    {
        test_and_initialize() ;
        if( !kdtree_ ) {
            const_cast< GeoModelMeshVertices* >( this )->initialize_kdtree() ;
        }
    }

    void GeoModelMeshVertices::initialize()
    {
        mesh_.clear() ;

        // Total number of vertices in the
        // Corners, Lines, and Surfaces of the GeoModel
        index_t nb = 0 ;
        for( index_t t = GME::CORNER; t <= GME::REGION; ++t ) {
            GME::TYPE T = static_cast< GME::TYPE >( t ) ;
            for( index_t e = 0; e < gm_.nb_elements( T ); ++e ) {
                nb += gm_.mesh_element( GME::gme_t( T, e ) ).nb_vertices() ;
            }
        }
        // Get out if no vertices
        if( nb == 0 ) {
            return ;
        }

        // Fill the vertices
        mesh_.vertices.create_vertices( nb ) ;
        gme_vertices_.resize( nb ) ;

        index_t index = 0 ;
        for( index_t t = GME::CORNER; t <= GME::REGION; ++t ) {
            GME::TYPE T = static_cast< GME::TYPE >( t ) ;
            for( index_t e = 0; e < gm_.nb_elements( T ); ++e ) {
                GeoModelMeshElement& E = cast_gmm_element( gm_, T, e ) ;
                if( E.nb_vertices() == 0 ) continue ;
                GEO::Memory::copy( mesh_.vertices.point_ptr( index ),
                    E.vertex( 0 ).data(), 3 * E.nb_vertices() * sizeof(double) ) ;
                for( index_t v = 0; v < E.nb_vertices(); v++ ) {
                    // Global index stored at BME level
                    E.set_model_vertex_id( v, index ) ;
                    // Index in the BME stored at global level
                    gme_vertices_[index].push_back( VertexInGME( E.gme_id(), v ) ) ;
                    // Global vertex index increment
                    index++ ;
                }
            }
        }
        // Remove colocated vertices
        remove_colocated() ;
    }

    void GeoModelMeshVertices::clear()
    {
//        GEO::Process::acquire_spinlock( lock_ ) ;

        mesh_.clear( true, true ) ;
        gme_vertices_.clear() ;
        clear_kdtree() ;

        // Clear the model vertex id information
        // for the Corner - Line - Surface - REGION
        for( index_t t = GME::CORNER; t <= GME::REGION; ++t ) {
            GME::TYPE T = static_cast< GME::TYPE >( t ) ;
            RINGMESH_PARALLEL_LOOP_DYNAMIC
            for( index_t e = 0; e < gm_.nb_elements( T ); ++e ) {
                GeoModelMeshElement& E = cast_gmm_element( gm_, T, e ) ;
                for( index_t v = 0; v < E.nb_vertices(); v++ ) {
                    E.set_model_vertex_id( v, NO_ID ) ;
                }
            }
        }
//        GEO::Process::release_spinlock( lock_ ) ;
    }

    void GeoModelMeshVertices::clear_kdtree()
    {
        if( kdtree_ ) {
            delete kdtree_ ;
            kdtree_ = nil ;
        }
    }

    void GeoModelMeshVertices::initialize_kdtree()
    {
        kdtree_ = new ColocaterANN( mesh_ );
#ifdef RINGMESH_DEBUG
        // Paranoia
        GEO::vector< index_t > old2new ;
        ringmesh_debug_assert(
            GEO::Geom::colocate( mesh_.vertices.point_ptr( 0 ), 3,
                mesh_.vertices.nb(), old2new, epsilon, mesh_.vertices.dimension() )
                == mesh_.vertices.nb() ) ;
#endif
    }

    index_t GeoModelMeshVertices::nb() const
    {
        test_and_initialize() ;
        ringmesh_debug_assert( gme_vertices_.size() == mesh_.vertices.nb() ) ;
        return mesh_.vertices.nb() ;
    }

    const vec3& GeoModelMeshVertices::vertex( index_t v ) const
    {
        test_and_initialize() ;
        ringmesh_debug_assert( v < nb() ) ;
        return mesh_.vertices.point( v ) ;
    }

    index_t GeoModelMeshVertices::index( const vec3& p ) const
    {
        test_kdtree_and_initialize() ;
        std::vector< index_t > vertices ;
        kdtree_->get_colocated( p, vertices ) ;
        if( vertices.empty() ) {
            return NO_ID ;
        } else {
            ringmesh_debug_assert( vertices.size() == 1 ) ;
            return vertices[0] ;
        }
    }

    const std::vector< GeoModelMeshVertices::VertexInGME >&
    GeoModelMeshVertices::gme_vertices( index_t v ) const
    {
        test_and_initialize() ;
        return gme_vertices_[v] ;
    }

    index_t GeoModelMeshVertices::add_vertex( const vec3& point )
    {
        clear_kdtree() ;
        gme_vertices_.push_back( std::vector< VertexInGME >() ) ;
        return mesh_.vertices.create_vertex( point.data() ) ;
    }

    void GeoModelMeshVertices::add_to_bme(
        index_t v,
        const VertexInGME& v_gme )
    {
        test_and_initialize() ;
        ringmesh_debug_assert( v < nb() ) ;
        ringmesh_debug_assert( gme_vertices_.size() == nb() ) ;
        // Assert if adding twice the same thing - not a normal behavior
        ringmesh_debug_assert(
            std::find( gme_vertices_[v].begin(), gme_vertices_[v].end(), v_gme )
                == gme_vertices_[v].end() ) ;

        gme_vertices_[v].push_back( v_gme ) ;
    }

    void GeoModelMeshVertices::set_gme(
        index_t unique_id,
        index_t k,
        const VertexInGME& v )
    {
        test_and_initialize() ;
        ringmesh_debug_assert( unique_id < nb() ) ;
        ringmesh_debug_assert( k < gme_vertices( unique_id ).size() ) ;
        gme_vertices_[unique_id][k] = v ;
    }

    void GeoModelMeshVertices::update_point( index_t v, const vec3& point )
    {
        test_and_initialize() ;
        ringmesh_debug_assert( v < nb() ) ;
        // Change the position of the unique_vertex
        mesh_.vertices.point( v ) = point ;
        clear_kdtree() ;

        const std::vector< VertexInGME >& gme_v = gme_vertices( v ) ;
        for( index_t i = 0; i < gme_v.size(); i++ ) {
            const VertexInGME& info = gme_v[i] ;
            const_cast< GMME& >( gm_.mesh_element( GME::gme_t( info.gme_id ) ) ).set_vertex(
                info.v_id, point, false ) ;
        }
    }

    void GeoModelMeshVertices::remove_colocated()
    {
        // Get out if nothing to do
        // and compute the points if they are not initialized yet
        if( nb() == 0 ) {
            return ;
        }
        // Identify and invalidate colocated vertices
        GEO::vector< index_t > old2new ;
        if( GEO::Geom::colocate( mesh_.vertices.point_ptr( 0 ), 3,
            mesh_.vertices.nb(), old2new, epsilon, mesh_.vertices.dimension() )
            != mesh_.vertices.nb() ) {
            std::vector< index_t > stupid_copy( old2new.begin(), old2new.end() ) ;
            erase_vertices( stupid_copy ) ;
        }
    }

    void GeoModelMeshVertices::erase_invalid_vertices()
    {

        index_t nb_todelete = 0 ;
        std::vector< index_t > to_delete( nb() ) ; // Here nb() represents the number of vertices before removal of the elements

        for( index_t v = 0; v < nb(); ++v ) {
            std::vector< VertexInGME >& related = gme_vertices_[v] ;
            index_t nb_invalid = 0 ;

            // Get the invalid BMEVertices for the current global vertex
            for( index_t i = 0; i < related.size(); ++i ) {

                if( !related[i].is_defined() ) {
                    // To ease removal of invalid BMEVertices
                    related[i] = VertexInGME() ;
                    nb_invalid++ ;
                }
            }

            if( nb_invalid < related.size() ) {
                to_delete[v] = v ;
                related.erase(
                    std::remove( related.begin(), related.end(), VertexInGME() ),
                    related.end() ) ;
            } else {
                // This vertex must be deleted
                to_delete[v] = NO_ID ;
                nb_todelete++ ;
                // std::erase of all elements has an undefined behavior
                related.clear() ;
            }
        }

        if( nb_todelete > 0 ) {
            erase_vertices( to_delete ) ;
        }
    }

    void GeoModelMeshVertices::erase_vertices( std::vector< index_t >& to_delete )
    {
        ringmesh_debug_assert( to_delete.size() == nb() ) ;

        // For mesh vertices deletion
        GEO::vector< index_t > to_delete_geo( nb(), 0 ) ;

        // Fill the delete information for geogram
        // Recycle the to_delete vertex to get the mapping between
        // new and old points. This is implemented to be the same
        // as what is done in the delete_elements function in geogram
        index_t nb_todelete = 0 ;
        index_t cur = 0 ;
        for( index_t v = 0; v < nb(); ++v ) {
            if( to_delete[v] != v ) {
                to_delete_geo[v] = 1 ;
                nb_todelete++ ;
                if( to_delete[v] != NO_ID ) {
                    ringmesh_debug_assert( to_delete[v] < v ) ;
                    to_delete[v] = to_delete[to_delete[v]] ;
                }
            } else {
                to_delete[v] = cur ;
                ++cur ;
            }
        }
        if( nb_todelete == 0 ) {
            return ;
        }
        if( nb_todelete == nb() ) {
            // Clear everything
            clear() ;
            return ;
        }

        // Empty the gme_vertices_ of the deleted vertices and erase them
        for( index_t v = 0; v < nb(); ++v ) {
            if( to_delete_geo[v] == 1 ) {
                gme_vertices_[v].clear() ;
            }
        }

        gme_vertices_.erase(
            std::remove( gme_vertices_.begin(), gme_vertices_.end(),
                std::vector< VertexInGME >() ), gme_vertices_.end() ) ;

        // Delete the vertices - false is to not remove
        // isolated vertices (here all the vertices)
        mesh_.vertices.delete_elements( to_delete_geo, false ) ;

#ifdef RINGMESH_DEBUG
        // Paranoia - check that we have the same mapping than the
        // delete_elements function in Geogram
        for( index_t v = 0; v < nb(); ++v ) {
            ringmesh_debug_assert(
                to_delete_geo[v] == NO_ID || to_delete_geo[v] == to_delete[v] ) ;
        }
#endif

        // Update model_vertex_ids in BMME
        for( index_t t = GME::CORNER; t <= GME::REGION; ++t ) {
            GME::TYPE T = static_cast< GME::TYPE >( t ) ;

            for( index_t e = 0; e < gm_.nb_elements( T ); ++e ) {
                GeoModelMeshElement& E = cast_gmm_element( gm_, T, e ) ;

                for( index_t v = 0; v < E.nb_vertices(); v++ ) {
                    index_t old_id = E.model_vertex_id( v ) ;
                    index_t new_id = to_delete[old_id] ;
                    // If new_id is NO_ID the vertex should be removed afterwards
                    // from the BMME
                    ringmesh_debug_assert( new_id != NO_ID ) ;
                    E.set_model_vertex_id( v, new_id ) ;

                    /*!
                     * @todo Review: I don't understand this for and what it does...
                     * When we remove a region, this for add stupid vertices inside the
                     * vector... [AB]
                     */
                    // Merge gme_vertices_ information
                    if( std::find( gme_vertices_[new_id].begin(),
                        gme_vertices_[new_id].end(), VertexInGME( E.gme_id(), v ) )
                        == gme_vertices_[new_id].end() ) {
                        gme_vertices_[new_id].push_back(
                            VertexInGME( E.gme_id(), v ) ) ;
                    }
                }
            }
        }

        // The Kd-tree should be updated next time we need it
        clear_kdtree() ;
    }

    /*******************************************************************************/

    GeoModelMeshFacets::GeoModelMeshFacets( GeoModelMesh& gmm, GEO::Mesh& mesh )
        :
            gmm_( gmm ),
            gm_( gmm.model() ),
            mesh_( mesh ),
            nb_triangle_( 0 ),
            nb_quad_( 0 )
    {
    }

    GeoModelMeshFacets::~GeoModelMeshFacets()
    {
        unbind_attribute() ;
    }

    void GeoModelMeshFacets::bind_attribute()
    {
        if( !surface_id_.is_bound() ) {
            surface_id_.bind( gmm_.facet_attribute_manager(), surface_att_name ) ;
        }
    }

    void GeoModelMeshFacets::unbind_attribute()
    {
        if( surface_id_.is_bound() ) {
            surface_id_.unbind() ;
        }
    }

    bool GeoModelMeshFacets::is_initialized() const
    {
        return mesh_.facets.nb() > 0 ;
    }

    index_t GeoModelMeshFacets::nb() const
    {
        test_and_initialize() ;
        return mesh_.facets.nb() ;
    }

    index_t GeoModelMeshFacets::vertex( index_t f, index_t v ) const
    {
        test_and_initialize() ;
        ringmesh_debug_assert( f < mesh_.facets.nb() ) ;
        ringmesh_debug_assert( v < mesh_.facets.nb_vertices( f ) ) ;
        return mesh_.facets.vertex( f, v ) ;
    }

    index_t GeoModelMeshFacets::adjacent( index_t f, index_t e ) const
    {
        test_and_initialize() ;
        ringmesh_debug_assert( f < mesh_.facets.nb() ) ;
        ringmesh_debug_assert( e < mesh_.facets.nb_vertices( f ) ) ;
        return mesh_.facets.adjacent( f, e ) ;
    }

    index_t GeoModelMeshFacets::surface( index_t f ) const
    {
        test_and_initialize() ;
        ringmesh_debug_assert( f < mesh_.facets.nb() ) ;
        return surface_id_[f] ;
    }

    index_t GeoModelMeshFacets::facet_in_surface( index_t f ) const
    {
        test_and_initialize() ;
        ringmesh_debug_assert( f < mesh_.facets.nb() ) ;
        return f - surface_facet_ptr_[ALL*surface( f )] ;
    }

    GeoModelMeshFacets::FacetType GeoModelMeshFacets::facet_type(
        index_t f,
        index_t& index ) const
    {
        test_and_initialize() ;
        ringmesh_debug_assert( f < mesh_.facets.nb() ) ;
        index_t facet = facet_in_surface( f ) ;
        index_t s = surface( f ) ;
        for( index_t t = TRIANGLE; t < ALL; t++ ) {
            FacetType T = static_cast< FacetType >( t ) ;
            if( facet < nb_facets( s, T ) ) {
                index = facet ;
                return T ;
            }
            facet -= nb_facets( s, T ) ;
        }
        index = NO_ID ;
        ringmesh_assert_not_reached ;
        return NO_FACET ;
    }

    index_t GeoModelMeshFacets::nb_facets( FacetType type ) const
    {
        test_and_initialize() ;
        switch( type ) {
            case TRIANGLE:
                return nb_triangle() ;
            case QUAD:
                return nb_quad() ;
            case POLYGON:
                return nb_polygon() ;
            case ALL:
                return nb() ;
            default:
                ringmesh_assert_not_reached ;
                return 0 ;
        }
    }

    index_t GeoModelMeshFacets::nb_facets( index_t s, FacetType type ) const
    {
        test_and_initialize() ;
        ringmesh_debug_assert( s < gm_.nb_surfaces() ) ;
        switch( type ) {
            case TRIANGLE:
                return nb_triangle( s ) ;
            case QUAD:
                return nb_quad( s ) ;
            case POLYGON:
                return nb_polygon( s ) ;
            case ALL:
                return surface_facet_ptr_[ALL * ( s + 1 )]
                    - surface_facet_ptr_[ALL * s] ;
            default:
                ringmesh_assert_not_reached ;
                return 0 ;
        }
    }

    index_t GeoModelMeshFacets::facet( index_t s, index_t f, FacetType type ) const
    {
        test_and_initialize() ;
        ringmesh_debug_assert( s < gm_.nb_surfaces() ) ;
        switch( type ) {
            case TRIANGLE:
                return triangle( s, f ) ;
            case QUAD:
                return quad( s, f ) ;
            case POLYGON:
                return polygon( s, f ) ;
            case ALL:
                return surface_facet_ptr_[ALL * s] + f ;
            default:
                ringmesh_assert_not_reached ;
                return 0 ;
        }
    }

    index_t GeoModelMeshFacets::nb_triangle() const
    {
        test_and_initialize() ;
        return nb_triangle_ ;
    }

    index_t GeoModelMeshFacets::nb_triangle( index_t s ) const
    {
        test_and_initialize() ;
        ringmesh_debug_assert( s < gm_.nb_surfaces() ) ;
        return surface_facet_ptr_[ALL * s + (TRIANGLE+1)]
            - surface_facet_ptr_[ALL * s + TRIANGLE] ;
    }

    index_t GeoModelMeshFacets::triangle( index_t s, index_t t ) const
    {
        test_and_initialize() ;
        ringmesh_debug_assert( s < gm_.nb_surfaces() ) ;
        return surface_facet_ptr_[ALL * s + TRIANGLE] + t ;
    }

    index_t GeoModelMeshFacets::nb_quad() const
    {
        test_and_initialize() ;
        return nb_quad_ ;
    }

    index_t GeoModelMeshFacets::nb_quad( index_t s ) const
    {
        test_and_initialize() ;
        ringmesh_debug_assert( s < gm_.nb_surfaces() ) ;
        return surface_facet_ptr_[ALL * s + (QUAD+1)]
            - surface_facet_ptr_[ALL * s + QUAD] ;
    }

    index_t GeoModelMeshFacets::quad( index_t s, index_t q ) const
    {
        test_and_initialize() ;
        ringmesh_debug_assert( s < gm_.nb_surfaces() ) ;
        return surface_facet_ptr_[ALL * s + QUAD] + q ;
    }

    index_t GeoModelMeshFacets::nb_polygon() const
    {
        test_and_initialize() ;
        return nb_polygon_ ;
    }

    index_t GeoModelMeshFacets::nb_polygon( index_t s ) const
    {
        test_and_initialize() ;
        ringmesh_debug_assert( s < gm_.nb_surfaces() ) ;
        return surface_facet_ptr_[ALL * s + (POLYGON+1)]
            - surface_facet_ptr_[ALL * s + POLYGON] ;
    }

    index_t GeoModelMeshFacets::polygon( index_t s, index_t p ) const
    {
        test_and_initialize() ;
        ringmesh_debug_assert( s < gm_.nb_surfaces() ) ;
        return surface_facet_ptr_[ALL * s + POLYGON] + p ;
    }

    void GeoModelMeshFacets::clear()
    {
        surface_facet_ptr_.clear() ;
        nb_triangle_ = 0 ;
        nb_quad_ = 0 ;
        mesh_.facets.clear() ;
    }

    void GeoModelMeshFacets::test_and_initialize() const
    {
        if( !is_initialized() ) {
            const_cast< GeoModelMeshFacets* >( this )->initialize() ;
        }
    }

    class GeoModelMeshFacetsSort {
    public:
        GeoModelMeshFacetsSort(
            const GEO::Mesh& mesh,
            const GEO::Attribute< index_t >& surface_id )
            : mesh_( mesh ), surface_id_( surface_id )
        {
        }

        bool operator()( index_t i, index_t j ) const
        {
            if( surface_id_[i] != surface_id_[j] ) {
                return surface_id_[i] < surface_id_[j] ;
            } else {
                return mesh_.facets.nb_vertices( i ) < mesh_.facets.nb_vertices( j ) ;
            }
        }
    private:
        const GEO::Mesh& mesh_ ;
        const GEO::Attribute< index_t >& surface_id_ ;
    } ;

    void GeoModelMeshFacets::initialize()
    {
        gmm_.vertices.test_and_initialize() ;

        // Compute the total number of facet per type
        std::vector< index_t > nb_facet_per_type( ALL+1, 0 ) ;
        for( index_t s = 0; s < gm_.nb_surfaces(); s++ ) {
            const Surface& surface = gm_.surface( s ) ;
            if( surface.is_triangulated() ) {
                nb_facet_per_type[TRIANGLE+1] += surface.nb_cells() ;
            } else {
                for( index_t f = 0; f < surface.nb_cells(); f++ ) {
                    switch( surface.nb_vertices_in_facet( f ) ) {
                        case 3:
                            nb_facet_per_type[TRIANGLE+1]++ ;
                            break ;
                        case 4:
                            nb_facet_per_type[QUAD+1]++ ;
                            break ;
                        default:
                            nb_facet_per_type[POLYGON+1]++ ;
                            break ;
                    }
                }
            }
        }

        // Create triangles and quads, the polygons will be handle later
        if( nb_facet_per_type[TRIANGLE+1] ) {
            mesh_.facets.create_triangles( nb_facet_per_type[TRIANGLE+1] ) ;
        }
        if( nb_facet_per_type[QUAD+1] ) {
            mesh_.facets.create_quads( nb_facet_per_type[QUAD+1] ) ;
        }
        index_t nb_polygons = nb_facet_per_type[POLYGON+1] ;

        // Compute the facet offset
        for( index_t t = TRIANGLE+1; t < ALL; t++ ) {
            nb_facet_per_type[t+1] += nb_facet_per_type[t] ;
        }

        // Fill the triangles and quads created above
        std::vector< index_t > cur_facet_per_type( ALL, 0 ) ;
        for( index_t s = 0; s < gm_.nb_surfaces(); s++ ) {
            const Surface& surface = gm_.surface( s ) ;
            for( index_t f = 0; f < surface.nb_cells(); f++ ) {
                index_t nb_vertices = surface.nb_vertices_in_facet( f ) ;
                if( nb_vertices > 4 ) continue ;
                FacetType T = static_cast< FacetType >( nb_vertices - 3 ) ;
                index_t cur_facet = nb_facet_per_type[T] + cur_facet_per_type[T]++ ;
                for( index_t v = 0; v < nb_vertices; v++ ) {
                    mesh_.facets.set_vertex( cur_facet, v,
                        surface.model_vertex_id( f, v ) ) ;
                }
                surface_id_[cur_facet] = s ;
            }
        }

        // Create and fill polygons
        if( nb_polygons ) {
            for( index_t s = 0; s < gm_.nb_surfaces(); s++ ) {
                const Surface& surface = gm_.surface( s ) ;
                for( index_t f = 0; f < surface.nb_cells(); f++ ) {
                    index_t nb_vertices = surface.nb_vertices_in_facet( f ) ;
                    if( nb_vertices < 5 ) continue ;
                    GEO::vector< index_t > vertices( nb_vertices ) ;
                    for( index_t v = 0; v < nb_vertices; v++ ) {
                        vertices[v] = surface.model_vertex_id( f, v ) ;
                    }
                    index_t cur_facet = mesh_.facets.create_polygon( vertices ) ;
                    surface_id_[cur_facet] = s ;
                }
            }
        }

        // Compute facet adjacencies
        mesh_.facets.connect() ;

        // Permute facets to sort them per surface and per type
        // Example for a mesh with two surfaces and only triangles and quads
        // [TRGL,TRGL, .. , QUAD, QUAD .. , TRGL, TRGL, ... , QUAD, QUAD ..]
        // |          surface 0           |             surface 1           |
        GEO::vector< index_t > sorted_indices( mesh_.vertices.nb() ) ;
        for( index_t i = 0; i < mesh_.vertices.nb(); i++ ) {
            sorted_indices[i] = i ;
        }
        GeoModelMeshFacetsSort action( mesh_, surface_id_ ) ;
        GEO::sort( sorted_indices.begin(), sorted_indices.end(), action ) ;
        mesh_.facets.permute_elements( sorted_indices ) ;
    }


    /*******************************************************************************/


    GeoModelMesh::GeoModelMesh( const GeoModel& gm )
        :
            gm_( gm ),
            mesh_( new GEO::Mesh ),
            vertices( *this, *mesh_ ),
            facets( *this, *mesh_ )
    {
    }

    GeoModelMesh::~GeoModelMesh()
    {
        facets.unbind_attribute() ;
        delete mesh_ ;
    }

    void GeoModelMesh::remove_colocated_vertices()
    {
        vertices.remove_colocated() ;
    }

    void GeoModelMesh::erase_vertices( std::vector< index_t >& to_delete )
    {
        vertices.erase_vertices( to_delete ) ;
    }


    void GeoModelMesh::erase_invalid_vertices()
    {
        vertices.erase_invalid_vertices() ;
    }

}
