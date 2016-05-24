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

/*!
 * @file Implementation of all GeoModelElements classes
 * @author Jeanne Pellerin and Arnaud Botella 
 */

#include <ringmesh/geo_model_element.h>

#include <algorithm>
#include <fstream>
#include <set>
#include <stack>

#include <geogram/basic/logger.h>
#include <geogram/basic/geometry_nd.h>

#include <geogram/mesh/mesh_AABB.h>
#include <geogram/mesh/mesh_geometry.h>
#include <geogram/mesh/mesh_intersection.h>
#include <geogram/mesh/mesh_repair.h>
#include <geogram/mesh/mesh_topology.h>
#include <geogram/mesh/triangle_intersection.h>

#include <ringmesh/algorithm.h>
#include <ringmesh/geo_model.h>
#include <ringmesh/geo_model_validity.h>
#include <ringmesh/geogram_extension.h>
#include <ringmesh/geometry.h>

namespace {
    /* Definition of functions that we do not want exported in the interface */
    using namespace RINGMesh ;

    typedef GeoModelElement::gme_t gme_t ;
    typedef GeoModelMeshElement GMME ;

    /*!
     * @brief Checks that the model vertex indices of @param E
     *       are in a valid range
     */
    bool check_range_model_vertex_ids( const GMME& E )
    {
        /// Check that the stored model vertex indices are in a valid range
        for( index_t i = 0; i < E.nb_vertices(); ++i ) {
            if( E.model_vertex_id( i ) == NO_ID
                && E.model_vertex_id( i ) >= E.model().mesh.vertices.nb() ) {
                GEO::Logger::warn( "GeoModelElement" )
                    << "Invalid model vertex index in " << E.gme_id() << std::endl ;
                return false ;
            }
        }
        return true ;
    }

    index_t compute_nb_volume_connected_components( const Region& M )
    {
        static const index_t NO_COMPONENT = index_t( -1 ) ;
        std::vector< index_t > component( M.nb_polytope(), NO_COMPONENT ) ;
        index_t nb_components = 0 ;
        for( index_t cell = 0; cell < M.nb_polytope(); cell++ ) {
            if( component[cell] == NO_COMPONENT ) {
                std::stack< index_t > S ;
                S.push( cell ) ;
                component[cell] = nb_components ;
                do {
                    index_t cur_cell = S.top() ;
                    S.pop() ;
                    for( index_t facet = 0; facet < M.nb_cell_facets( cur_cell );
                        facet++ ) {
                        index_t adj_cell = M.cell_adjacent_index( cur_cell, facet ) ;
                        if( adj_cell != GEO::NO_CELL
                            && component[adj_cell] == NO_COMPONENT ) {
                            S.push( adj_cell ) ;
                            component[adj_cell] = nb_components ;
                        }
                    }
                } while( !S.empty() ) ;
                nb_components++ ;
            }
        }
        return nb_components ;
    }

    /*!
     * @brief Count the number of times each vertex is in an edge or facet
     *
     * @param[in] gmme The GeoModelMeshElement
     * @param[out] nb Resized to the number of vertices of the mesh.
     *      Number of times one vertex appear in an polytope collection of the GeoModelMeshElement edge or facet of the mesh.
     */
    void count_vertex_occurences(
        const GeoModelMeshElement& E,
        std::vector< index_t >& nb )
    {
        nb.resize( E.nb_vertices(), 0 ) ;
        for( index_t polytope_index = 0; polytope_index < E.nb_polytope();
            ++polytope_index ) {
            for( index_t vertex = 0;
                vertex < E.nb_polytope_vertices( polytope_index ); ++vertex ) {
                ++nb[E.polytope_vertex_index( polytope_index, vertex )] ;
            }
        }
    }

    index_t count_nb_isolated_vertices( const GeoModelMeshElement& mesh )
    {
        std::vector< index_t > nb ;
        count_vertex_occurences( mesh, nb ) ;
        return static_cast< index_t >( std::count( nb.begin(), nb.end(), 0 ) ) ;
    }


    bool check_mesh_element_vertices_are_different(
        std::vector< index_t >& vertices,
        std::vector< index_t >& vertices_global )
    {
        ringmesh_assert(
            std::count( vertices.begin(), vertices.end(), NO_ID ) == 0 ) ;
        ringmesh_assert(
            std::count( vertices_global.begin(), vertices_global.end(), NO_ID ) == 0 ) ;
        // 0 is the default value of the model_vertex_id
        // If we have only 0 either this is a degenerate facets, but most certainly
        // model vertex ids are not good
        ringmesh_assert(
            std::count( vertices_global.begin(), vertices_global.end(), 0 )
            != vertices_global.size() ) ;

        std::sort( vertices.begin(), vertices.end() ) ;
        std::sort( vertices_global.begin(), vertices_global.end() ) ;
        return std::unique( vertices.begin(), vertices.end() ) != vertices.end()
            || std::unique( vertices_global.begin(), vertices_global.end() )
                != vertices_global.end() ;
    }

    /*!
     * @brief Returns true if the surface facet is incident twice to the same vertex
     */
    bool facet_is_degenerate( const Surface& S, index_t f )
    {
        index_t nb_facet_vertices = S.nb_polytope_vertices( f ) ;
        std::vector< index_t > corners( nb_facet_vertices, NO_ID ) ;
        std::vector< index_t > corners_global( nb_facet_vertices, NO_ID ) ;
        index_t v = 0 ;
        for( index_t c = S.facet_begin( f ); c < S.facet_end( f ); ++c ) {
            corners[v] = c ;
            corners_global[v] = S.model_vertex_id( f, v ) ;
            v++ ;
        }
        return check_mesh_element_vertices_are_different( corners, corners_global ) ;
    }    
    /*!
     * @brief Returns true if the region cell is incident twice to the same vertex
     * or if the cell volume is negative or inferior to epsilon
     */
    bool cell_is_degenerate( const Region& region, index_t cell_index )
    {
        index_t nb_vertices_in_cell = region.nb_polytope_vertices( cell_index ) ;
        std::vector< index_t > vertices( nb_vertices_in_cell, NO_ID ) ;
        std::vector< index_t > vertices_global( nb_vertices_in_cell, NO_ID ) ;
        for( index_t v = 0; v < nb_vertices_in_cell; v++ ) {
            vertices[v] = region.polytope_vertex_index( cell_index, v ) ;
            vertices_global[v] = region.model_vertex_id( cell_index, v ) ;
        }
        double volume = region.polytope_size( cell_index ) ;
        return check_mesh_element_vertices_are_different( vertices, vertices_global )
            || volume < epsilon ;
	}

    /*!
     * @brief Debug: Save a Surface of the model in the file OBJ format is used
     * @todo Move this function to an API providing utility functions on a
     * GeoModel and its Elements ? [JP]
     */
    void save_surface_as_obj_file( const Surface& S, const std::string& file_name )
    {
        std::ofstream out( file_name.c_str() ) ;
        if( out.bad() ) {
            GEO::Logger::err( "I/O" ) << "Error when opening the file: "
                << file_name.c_str() << std::endl ;
            return ;
        }
        out.precision( 16 ) ;
        for( index_t p = 0; p < S.nb_vertices(); p++ ) {
            const vec3& V = S.vertex( p ) ;
            out << "v" << " " << V.x << " " << V.y << " " << V.z << std::endl ;
        }
        for( index_t f = 0; f < S.nb_polytope(); f++ ) {
            out << "f" << " " ;
            for( index_t v = 0; v < S.nb_polytope_vertices( f ); v++ ) {
                out << S.polytope_vertex_index( f, v ) + 1 << " " ;
            }
            out << std::endl ;
        }
    }
}

namespace RINGMesh {
    /*!
     * @brief Map the name of a geological type with a value of GEOL_FEATURE
     *
     * @param[in] in Name of the feature. Can be
     * \li "reverse_fault"
     * \li "normal_fault"
     * \li "fault"
     * \li "top"
     * \li "none"
     * \li "topographic"
     * \li "unconformity"
     * \li "boundary"
     * Other strings will end up in \p NO_GEOL
     * @return The geological feature index
     * @todo Add other types of unconformity, see RINGMesh::GeoModelElement::TYPE. --GC
     */
    GeoModelElement::GEOL_FEATURE GeoModelElement::determine_geological_type(
        const std::string& in )
    {
        if( in == "reverse_fault" ) {
            return REVERSE_FAULT ;
        } else if( in == "normal_fault" ) {
            return NORMAL_FAULT ;
        } else if( in == "fault" ) {
            return FAULT ;
        } else if( in == "top" ) {
            return STRATI ;
        } else if( in == "none" ) {
            // This might seem strange - but it seems that what's
            // Gocad is doing
            return STRATI ;
        } else if( in == "topographic" ) {
            return STRATI ;
        } else if( in == "unconformity" ) {
            return UNCONFORMITY ;
        } else if( in == "boundary" ) {
            return VOI ;
        } else {
            // Default case - no information
            return NO_GEOL ;
        }
    }

    /*!
     * \return the (uppercase) string associated to a GeoModelELement::TYPE
     */
    std::string GeoModelElement::type_name( GME::TYPE t )
    {
        switch( t ) {
            case CORNER:
                return "CORNER" ;
            case LINE:
                return "LINE" ;
            case SURFACE:
                return "SURFACE" ;
            case REGION:
                return "REGION" ;
            case CONTACT:
                return "CONTACT" ;
            case INTERFACE:
                return "INTERFACE" ;
            case LAYER:
                return "LAYER" ;
            default:
                return "NO_TYPE_NAME" ;
        }
    }

    /*!
     * \return the (lowercase) string associated to a
     * GeoModelELement::GEOL_FEATURE
     */
    std::string GeoModelElement::geol_name( GME::GEOL_FEATURE t )
    {
        switch( t ) {
            case STRATI:
                return "top" ;
            case FAULT:
                return "fault" ;
            case REVERSE_FAULT:
                return "reverse_fault" ;
            case NORMAL_FAULT:
                return "normal_fault" ;
            case UNCONFORMITY:
                return "unconformity" ;
            case VOI:
                return "boundary" ;
            case NO_GEOL:
                return "no_geological_feature" ;
            default:
                return "no_geological_feature" ;
                break ;
        }
    }

    /*!
     * @brief Defines the type of the parent of an element of type @param t
     *        If no parent is allowed returns NO_TYPE
     * @details The elements that can have a parent are LINE, SURFACE, and REGION
     */
    GeoModelElement::TYPE GeoModelElement::parent_type( GME::TYPE t )
    {
        switch( t ) {
            case LINE:
                return CONTACT ;
            case SURFACE:
                return INTERFACE ;
            case REGION:
                return LAYER ;
            default:
                // The others have no parent
                return NO_TYPE ;
        }
    }

    /*!
     * @brief Defines the type of a child of an element of type @param t
     *        If no child is allowed returns NO_TYPE
     * @details The elements that can have a parent are CONTACT, INTERFACE, and LAYER
     */
    GeoModelElement::TYPE GeoModelElement::child_type( GME::TYPE t )
    {
        switch( t ) {
            case CONTACT:
                return LINE ;
            case INTERFACE:
                return SURFACE ;
            case LAYER:
                return REGION ;
            default:
                return NO_TYPE ;
        }
    }

    /*!
     * @brief Defines the type of an element on the boundary of an element of type @param t
     *        If no boundary is allowed returns NO_TYPE
     * @details The elements that can have a boundary are LINE, SURFACE, and REGION
     */
    GeoModelElement::TYPE GeoModelElement::boundary_type( GeoModelElement::TYPE t )
    {
        switch( t ) {
            case LINE:
                return CORNER ;
            case SURFACE:
                return LINE ;
            case REGION:
                return SURFACE ;
            default:
                return NO_TYPE ;
        }
    }

    /*!
     * @brief Defines the type of an element into which boundary an element of type @param t can be
     *        If no in_boundary is allowed returns NO_TYPE
     * @details The elements that can be in the boundary of another are CORNER, LINE, and SURFACE
     */
    GeoModelElement::TYPE GeoModelElement::in_boundary_type(
        GeoModelElement::TYPE t )
    {
        switch( t ) {
            case CORNER:
                return LINE ;
            case LINE:
                return SURFACE ;
            case SURFACE:
                return REGION ;
            default:
                return NO_TYPE ;
        }
    }

    /*!
     * @brief Dimension 0, 1, 2, or 3 of an element of type @param t
     */
    index_t GeoModelElement::dimension( GME::TYPE t )
    {
        switch( t ) {
            case CORNER:
                return 0 ;
            case LINE:
                return 1 ;
            case SURFACE:
                return 2 ;
            case REGION:
                return 3 ;
            case CONTACT:
                return 1 ;
            case INTERFACE:
                return 2 ;
            case LAYER:
                return 3 ;
            default:
                return NO_ID ;
        }
    }

    /*!
     * @brief Return true if \param type is a CORNER, LINE or SURFACE
     */
    bool GeoModelElement::has_mesh( GME::TYPE type )
    {
        return type <= REGION ;
    }

    bool GeoModelElement::is_connectivity_valid() const
    {
        bool valid = true ;

        /// 1. Check the validity of identification information
        ///    in the model - Universe has no index, but a TYPE
        if( gme_id() == gme_t() ) {
            GEO::Logger::err( "GeoModelElement" ) << " Element associated to model "
                << model().name() << "has no type and no index " << std::endl ;
            valid = false ;
        }

        if( !valid ) {
            // If previous information are not valid
            // No further checks are possible 
            // This really should not happen 
            ringmesh_assert( valid ) ;
            return valid ;
        }

        if( index() >= model().nb_elements( type() ) ) {
            GEO::Logger::warn( "GeoModelElement" ) << " Element index " << gme_id()
                << " is not valid " << " There are " << model().nb_elements( type() )
                << " element of that type in model " << model().name() << std::endl ;
            // This really should not happen
            valid = false ;
            ringmesh_assert( valid ) ;
            return valid ;
        }
        if( &( model().element( gme_id() ) ) != this ) {
            GEO::Logger::err( "GeoModelElement" ) << " Element " << gme_id()
                << " in model " << model().name() << " does not match this element"
                << std::endl ;
            // This really should not happen
            ringmesh_assert( valid ) ;
            valid = false ;
            return valid ;
        }

        /// 2. Check that required information for the TYPE is defined
        ///    and that reverse information is stored by the corresponding
        ///    elements
        TYPE T = type() ;

        // Boundaries
        if( boundary_allowed( T ) ) {
            if( T == REGION ) {
                if( nb_boundaries() == 0 ) {
                    GEO::Logger::warn( "GeoModelElement" ) << gme_id()
                        << " has no boundaries " << std::endl ;
                    valid = false ;
                }
            }
            // A Line must have 2 corners - they are identical if the Line is closed
            if( T == LINE && nb_boundaries() != 2 ) {
                GEO::Logger::warn( "GeoModelElement" ) << gme_id()
                    << " does not have 2 corners" << std::endl ;
                valid = false ;
            }
            // No requirement on Surface - it may have no boundary - bubble

            // All elements in the boundary must have this in their
            // in_boundary vector
            for( index_t i = 0; i < nb_boundaries(); ++i ) {
                const GME& E = boundary( i ) ;
                bool found = false ;
                index_t j = 0 ;
                while( !found && j < E.nb_in_boundary() ) {
                    if( E.in_boundary_gme( j ) == gme_id() ) {
                        found = true ;
                    }
                    j++ ;
                }
                if( !found ) {
                    GEO::Logger::warn( "GeoModelElement" )
                        << "Inconsistency boundary-in_boundary between " << gme_id()
                        << " and " << E.gme_id() << std::endl ;
                    valid = false ;
                }
            }
        }

        // In_boundary
        if( in_boundary_allowed( T ) ) {
            // Fix for a .ml for which VOI Surface are only on the boundary of Universe
            // Can we keep this ? Or should we compute the Region
            if( nb_in_boundary() == 0 ) {
                GEO::Logger::warn( "GeoModelElement" ) << gme_id()
                    << " is in the boundary of no element " << std::endl ;
                valid = false ;
            }

            // All elements in the in_boundary must have this in their
            // boundary vector
            for( index_t i = 0; i < nb_in_boundary(); ++i ) {
                const GME& E = in_boundary( i ) ;
                bool found = false ;
                index_t j = 0 ;
                while( !found && j < E.nb_boundaries() ) {
                    if( E.boundary_gme( j ) == gme_id() ) {
                        found = true ;
                    }
                    j++ ;
                }
                if( !found ) {
                    GEO::Logger::warn( "GeoModelElement" )
                        << "Inconsistency in_boundary-boundary between " << gme_id()
                        << " and " << E.gme_id() << std::endl ;
                    valid = false ;
                }
            }
        }

        // Parent - High level elements are not mandatory
        // But if the model has elements of the parent type, the element must have a parent
        if( parent_allowed( T ) ) {
            bool model_has_parent_elements(
                model().nb_elements( parent_type( T ) ) > 0 ) ;
            if( model_has_parent_elements ) {
                if( has_parent() ) {
                    const GME& E = parent() ;
                    // The parent must have this element in its children
                    bool found = false ;
                    index_t j = 0 ;
                    while( !found && j < E.nb_children() ) {
                        if( E.child_id( j ) == gme_id() ) {
                            found = true ;
                        }
                        j++ ;
                    }
                    if( !found ) {
                        GEO::Logger::warn( "GeoModelElement" )
                            << "Inconsistency parent-child between " << gme_id()
                            << " and " << E.gme_id() << std::endl ;
                        valid = false ;
                    }
                } else {
                    GEO::Logger::warn( "GeoModelElement" ) << gme_id()
                        << " has no geological parent element " << std::endl ;
                    valid = false ;
                }
            }
        }

        // Children
        if( child_allowed( T ) ) {
            if( nb_children() == 0 ) {
                GEO::Logger::warn( "GeoModelElement" ) << gme_id()
                    << " has no children mesh element, so no geometry "
                    << std::endl ;
                valid = false ;
            }

            // All children must have this element as a parent
            for( index_t i = 0; i < nb_children(); ++i ) {
                if( child( i ).parent_id() != gme_id() ) {
                    GEO::Logger::warn( "GeoModelElement" )
                        << "Inconsistency child-parent between " << gme_id()
                        << " and " << child( i ).gme_id() << std::endl ;
                    valid = false ;
                }
            }
        }
        return valid ;
    }

    /*!
     * @return Assert that the parent exists and returns it.
     */
    const GeoModelElement& GeoModelElement::parent() const
    {
        ringmesh_assert( parent_id().is_defined() ) ;
        return model().element( parent_id() ) ;
    }

    /*!
     *
     * @param[in] x Index of the boundary element
     * @return Asserts that is exists and returns the element on the boundary
     */
    const GeoModelElement& GeoModelElement::boundary( index_t x ) const
    {
        ringmesh_assert( x < nb_boundaries() ) ;
        return model().element( boundary_gme( x ) ) ;
    }

    /*!
     *
     * @param[in] x Index of the in_boundary element
     * @return Asserts that it exists and returns the element in in_boundary.
     */
    const GeoModelElement& GeoModelElement::in_boundary( index_t x ) const
    {
        ringmesh_assert( x < nb_in_boundary() ) ;
        return model().element( in_boundary_gme( x ) ) ;
    }

    /*!
     *
     * @param[in] x Index of the child
     * @return Asserts that the child exists and returns it.
     */
    const GeoModelElement& GeoModelElement::child( index_t x ) const
    {
        ringmesh_assert( x < nb_children() ) ;
        return model().element( child_id( x ) ) ;
    }

    /*!
     * @brief Checks if this element define the model external boundary
     * @details Test if the element is in the Surfaces defining the universe 
     */
    bool GeoModelElement::is_on_voi() const
    {
        TYPE T = type() ;
        if( T == SURFACE ) {
            for( index_t i = 0; i < model().universe().nb_boundaries(); ++i ) {
                if( model().universe().boundary_gme( i ) == gme_id() ) {
                    return true ;
                }
            }
            return false ;
        } else if( T == LINE || T == CORNER ) {
            // True if one of the incident surface define the universe
            for( index_t i = 0; i < nb_in_boundary(); ++i ) {
                if( in_boundary( i ).is_on_voi() ) {
                    return true ;
                }
            }
            return false ;
        } else if( T == REGION || T == LAYER ) {
            return false ;
        } else if( T == INTERFACE || T == CONTACT ) {
            // Check that all children are on the voi
            if( nb_children() > 0 ) {
                for( index_t i = 0; i < nb_children(); ++i ) {
                    if( !child( i ).is_on_voi() ) {
                        return false ;
                    }
                }
                return true ;
            } else {
                return false ;
            }
        } else {
            ringmesh_assert( false ) ;
            return false ;
        }
    }

    /*!
     * @brief Check if this element an inside border of rhs
     * @details That can be Surface stopping in a Region, or Line stopping in a Surface.
     * @param[in] rhs The element to test
     */
    bool GeoModelElement::is_inside_border( const GeoModelElement& rhs ) const
    {
        // Find out if this surface is twice in the in_boundary vector
        return std::count( in_boundary_.begin(), in_boundary_.end(), rhs.gme_id() )
            > 1 ;
    }

    /*!
     * @brief Check if one element is twice in the boundary
     */
    bool GeoModelElement::has_inside_border() const
    {
        for( index_t i = 0; i < nb_boundaries(); ++i ) {
            if( boundary( i ).is_inside_border( *this ) ) {
                return true ;
            }
        }
        return false ;
    }

    /*********************************************************************/

    const std::string GeoModelMeshElement::model_vertex_id_att_name()
    {
        return "model_vertex_id" ;
    }

    GeoModelMeshElement::~GeoModelMeshElement()
    {
        unbind_attributes() ;
#ifdef RINGMESH_DEBUG
        mesh_.print_mesh_bounded_attributes() ;
#endif
    }

    /*!
     * @brief Binds attributes stored by the GME on the Mesh
     */
    void GeoModelMeshElement::bind_attributes()
    {
        model_vertex_id_.bind( mesh_.vertex_attribute_manager(),
            model_vertex_id_att_name() ) ;
    }
    /*!
     * @brief Unbinds attributes stored by the GME on the Mesh
     */
    void GeoModelMeshElement::unbind_attributes()
    {
        model_vertex_id_.unbind() ;
    }

    bool GeoModelMeshElement::are_model_vertex_indices_valid() const
    {
        bool valid = true ;
        // For all vertices
        // Check that the global vertex has an index backward to 
        // the vertex of this element
        for( index_t v = 0; v < nb_vertices(); ++v ) {
            index_t model_v = model_vertex_id( v ) ;

            const std::vector< GMEVertex >& backward =
                model().mesh.vertices.gme_vertices( model_v ) ;

            GMEVertex cur_v( gme_id(), v ) ;
            index_t count_v = static_cast< index_t >( std::count( backward.begin(),
                backward.end(), cur_v ) ) ;

            if( count_v != 1 ) {
                GEO::Logger::warn( "GeoModelElement" ) << gme_id() << " vertex " << v
                    << " appears " << count_v
                    << " in the related global model vertex " << model_v
                    << std::endl ;
                valid = false ;
            }
        }
        return valid ;
    }

    index_t GeoModelMeshElement::gmme_vertex_index_from_model(
        index_t model_vertex_id ) const
    {
        const std::vector< GMEVertex >& gme_vertices =
            model().mesh.vertices.gme_vertices( model_vertex_id ) ;

        for( index_t i = 0; i < gme_vertices.size(); i++ ) {
            const GMEVertex& info = gme_vertices[i] ;
            if( info.gme_id == gme_id() ) {
                return info.v_id ;
            }
        }
        return NO_ID ;
    }

    std::vector< index_t > GeoModelMeshElement::gme_vertex_indices(
        index_t model_vertex_id ) const
    {
        const std::vector< GMEVertex >& all_vertices =
            model().mesh.vertices.gme_vertices( model_vertex_id ) ;

        std::vector< index_t > this_gme_vertices ;
        for( index_t i = 0; i < all_vertices.size(); i++ ) {
            const GMEVertex& gme_vertex = all_vertices[i] ;
            if( gme_vertex.gme_id == gme_id() ) {
                this_gme_vertices.push_back( gme_vertex.v_id ) ;
            }
        }
        return this_gme_vertices ;
    }

    /**************************************************************/

    /*!
     * @brief Check that the Corner mesh is a unique point
     */
    bool Corner::is_mesh_valid() const
    {
        bool valid = true ;
        if( nb_vertices() != 1 ) {
            GEO::Logger::err( "GeoModelElement" ) << "Corner " << index()
                << " mesh has " << mesh_.nb_vertices() << " vertices " << std::endl ;
            valid = false ;
        }
        if( mesh_.nb_edges() != 0 ) {
            GEO::Logger::err( "GeoModelElement" ) << "Corner " << index()
                << " mesh has " << mesh_.nb_edges() << " edges " << std::endl ;
            valid = false ;
        }
        if( mesh_.nb_facets() != 0 ) {
            GEO::Logger::err( "GeoModelElement" ) << "Corner " << index()
                << " mesh has " << mesh_.nb_facets() << " facets " << std::endl ;
            valid = false ;
        }
        if( mesh_.nb_cells() != 0 ) {
            GEO::Logger::warn( "GeoModelElement" ) << "Corner " << index()
                << " mesh has " << mesh_.nb_cells() << " cells " << std::endl ;
            valid = false ;
        }
        // The default point is (0., 0., 0.) and there might be a valid
        // Corner at this position.
        /*if( mesh_.vertices.point( 0 ) == vec3() ) {
         GEO::Logger::warn( "GeoModelElement" )
         << "Corner " << index()
         << " point is default " << std::endl ;
         valid = false ;
         }*/
        return valid ;
    }

    /***************************************************************/

    /*!
     * @brief Construct a Line
     *
     * @param[in] model The parent model
     * @param[in] id The index of the line in the lines_ vector of the parent model
     */
    Line::Line( const GeoModel& model, index_t id )
        : GeoModelMeshElement( model, LINE, id )
    {
    }

    /*!
     * @brief Check that the mesh of the Line is valid
     * @details Check that 
     *  - the GEO::Mesh has more than 1 vertex - more than 1 edge - no facets - no cells.
     *  - global indices of vertices in the model are in a valid range 
     *  - each vertex is in 2 edges except extremities that are in 1 edge
     * 
     * Does not check:
     *  - Self-intersection - I suppose there are no segment - segment intersection (JP)
     *  - Duplicated edge - most probably ruled out with the duplicated vertex test (JP)
     *  - Duplicated vertex (verified at GeoModel level)
     */
    bool Line::is_mesh_valid() const
    {
        bool valid = true ;

        // Check that the GEO::Mesh has the expected elements
        if( nb_vertices() < 2 ) {
            GEO::Logger::err( "GeoModelElement" ) << "Line " << index() << " has "
                << mesh_.nb_vertices() << " vertices " << std::endl ;
            valid = false ;
        }
        if( mesh_.nb_edges() == 0 ) {
            GEO::Logger::err( "GeoModelElement" ) << "Line " << index()
                << " mesh has " << mesh_.nb_edges() << " edges " << std::endl ;
            valid = false ;
        }
        if( mesh_.nb_facets() != 0 ) {
            GEO::Logger::err( "GeoModelElement" ) << "Line " << index()
                << " mesh has " << mesh_.nb_facets() << " facets " << std::endl ;
            valid = false ;
        }
        if( mesh_.nb_cells() != 0 ) {
            GEO::Logger::err( "GeoModelElement" ) << "Line " << index()
                << " mesh has " << mesh_.nb_cells() << " cells " << std::endl ;
            valid = false ;
        }

        // Model indices must be valid
        valid = check_range_model_vertex_ids( *this ) && valid ;

        if( nb_vertices() > 1 ) {
            // Count the number of edges in which each vertex is
            std::vector< index_t > nb ;
            count_vertex_occurences( *this, nb ) ;
            index_t nb0 = 0 ;
            index_t nb1 = 0 ;
            index_t nb2 = 0 ;
            for( index_t i = 0; i < nb.size(); ++i ) {
                if( nb[i] == 0 )
                    ++nb0 ;
                else if( nb[i] == 1 )
                    ++nb1 ;
                else if( nb[i] == 2 ) ++nb2 ;
            }

            // Vertices at extremitites must be in only one edge
            if( nb.front() != 1 || nb.back() != 1 ) {
                GEO::Logger::err( "GeoModelElement" )
                    << "Invalid extremity points in " << gme_id() << std::endl ;
                valid = false ;
            }
            // No isolated vertices are allowed
            if( nb0 > 0 ) {
                GEO::Logger::warn( "GeoModelElement" ) << nb0
                    << " isolated vertices in " << gme_id() << std::endl ;
                valid = false ;
            }
            // Only the two extremities are in only 1 edge 
            // One connected component condition
            if( nb1 != 2 ) {
                GEO::Logger::warn( "GeoModelElement" )
                    << "More than one connected component for " << gme_id()
                    << std::endl ;
                valid = false ;
            }
            // All the others must be in 2 edges and 2 edges only
            // Manifold condition
            if( nb2 != nb.size() - 2 ) {
                GEO::Logger::warn( "GeoModelElement" ) << "Non-manifold element"
                    << gme_id() << std::endl ;
                valid = false ;
            }
        }

        // No zero edge length
        index_t nb_degenerated = 0 ;
        for( index_t e = 0; e < nb_polytope(); ++e ) {
            double l = length( polytope_vertex( e, 1 ) - polytope_vertex( e, 0 ) ) ;
            if( l < epsilon ) {
                nb_degenerated++ ;
            }
        }
        if( nb_degenerated > 0 ) {
            GEO::Logger::warn( "GeoModelElement" ) << nb_degenerated
                << " degenerated edges in " << gme_id() << std::endl ;
            valid = false ;
        }
        return valid ;
    }

    /********************************************************************/

    /*!
     * @brief Check that the mesh of the Surface is valid
     * @details Check that
     *  - the GEO::Mesh has more than 2 vertices, at least 1 facet, no cells.
     *  - global indices of vertices in the model are in a valid range
     *  - no degenerate facet 
     *  - one connected component 
     *
     *  Some tests are not performed here but globally on the GeoModel
     *  - intersection of facets 
     *  - non-manifold edges 
     *  - duplicated vertices are on a boundary Line ending in the Surface 
     * 
     *
     *  Some tests are not performed     
     *  - non-manifold points
     *  - surface orientability
     *  - planarity of polygonal facets 
     *
     * @todo Check that there is no duplicated facet 
     */
    bool Surface::is_mesh_valid() const
    {
        bool valid = true ;
        // Check that the GEO::Mesh has the expected elements
        // at least 3 vertices and one facet.
        if( nb_vertices() < 3 ) {
            GEO::Logger::warn( "GeoModelElement" ) << gme_id()
                << " has less than 3 vertices " << std::endl ;
            valid = false ;
        }
        // Is it important to have edges or not ?
        // I would say we do not care (JP) - so no check on that 
        if( mesh_.nb_facets() == 0 ) {
            GEO::Logger::warn( "GeoModelElement" ) << gme_id() << " has no facets "
                << std::endl ;
            valid = false ;
        }
        if( mesh_.nb_cells() != 0 ) {
            GEO::Logger::warn( "GeoModelElement" ) << gme_id() << " has "
                << mesh_.nb_cells() << " cells " << std::endl ;
            valid = false ;
        }

        // No isolated vertices
        index_t nb_isolated_vertices = count_nb_isolated_vertices( *this ) ;
        if( nb_isolated_vertices > 0 ) {
            GEO::Logger::warn( "GeoModelElement" ) << gme_id() << " mesh has "
                << nb_isolated_vertices << " isolated vertices " << std::endl ;
            valid = false ;
        }

        // No zero area facet
        // No facet incident to the same vertex check local and global indices
        index_t nb_degenerate = 0 ;
        for( index_t f = 0; f < mesh_.nb_facets(); f++ ) {
            if( facet_is_degenerate( *this, f ) ) {
                nb_degenerate++ ;
            }
        }
        if( nb_degenerate != 0 ) {
            GEO::Logger::warn( "GeoModelElement" ) << gme_id() << " mesh has "
                << nb_degenerate << " degenerate facets " << std::endl ;
            valid = false ;
        }

        // No duplicated facet
        GEO::vector< index_t > colocated ;
        // GEO::mesh_detect_duplicated_facets( mesh_, colocated ) ; // not implemented yet 
        index_t nb_duplicated_f = 0 ;
        for( index_t f = 0; f < colocated.size(); ++f ) {
            if( colocated[f] != f ) {
                nb_duplicated_f++ ;
            }
        }
        if( nb_duplicated_f > 0 ) {
            GEO::Logger::warn( "GeoModelElement" ) << gme_id() << " mesh has "
                << nb_duplicated_f << " duplicated facets " << std::endl ;
            valid = false ;
        }

        // One connected component  
        index_t cc = mesh_.nb_connected_components() ;
        if( cc != 1 ) {
            GEO::Logger::warn( "GeoModelElement" ) << gme_id() << " mesh has " << cc
                << " connected components " << std::endl ;
            valid = false ;
#ifdef RINGMESH_DEBUG
            std::ostringstream file ;
            file << validity_errors_directory << "/" << "invalid_surf_" << index()
                << ".obj" ;
            save_surface_as_obj_file( *this, file.str() ) ;

#endif  
        }
        return valid ;
    }

    /*!
     * @brief Traversal of a surface border
     * @details From the input facet f, get the facet that share vertex v and
     * get the indices of vertex v and of the following vertex in this new facet.
     * The next facet next_f may be the same, and from is required to avoid going back.
     *
     * @param[in] f Index of the facet
     * @param[in] from Index in the facet of the previous point on the border - gives the direction
     * @param[in] v Index in the facet of the point for which we want the next point on border
     * @param[out] next_f Index of the facet containing the next point on border
     * @param[out] v_in_next Index of vertex v in facet next_f
     * @param[out] next_in_next Index of the next vertex on border in facet v_in_next
     */
    void Surface::next_on_border(
        index_t f,
        index_t from,
        index_t v,
        index_t& next_f,
        index_t& v_in_next,
        index_t& next_in_next ) const
    {
        ringmesh_assert( v < nb_polytope_vertices( f ) ) ;
        ringmesh_assert( is_on_border( f, v ) || is_on_border( f, from ) ) ;

        index_t V = polytope_vertex_index( f, v ) ;

        // We want the next triangle that is on the boundary and share V
        // If there is no such triangle, the next vertex on the boundary
        // is the vertex of F neighbor of V that is not from

        // Get the facets around the shared vertex that are on the boundary
        // There must be one (the current one) or two (the next one on boundary)
        std::vector< index_t > facets ;
        index_t nb_around = facets_around_vertex( V, facets, true, f ) ;
        ringmesh_assert( nb_around < 3 && nb_around > 0 ) ;

        next_f = facets[0] ;

        if( nb_around == 2 ) {
            if( next_f == f ) {
                next_f = facets[1] ;
            }
            ringmesh_assert( next_f != NO_ID ) ;

            // Now get the other vertex that is on the boundary opposite to p1
            v_in_next = vertex_index_in_facet( next_f, V ) ;
            ringmesh_assert( v_in_next != NO_ID ) ;

            // The edges containing V in next_f are
            // the edge starting at v_in_next and the one ending there
            index_t prev_v_in_next = prev_facet_vertex_index( next_f, v_in_next ) ;

            bool e0_on_boundary = is_on_border( next_f, v_in_next ) ;

            // Only one must be on the boundary otherwise there is a corner missing
            ringmesh_assert( e0_on_boundary != is_on_border( next_f, prev_v_in_next ) ) ;

            // From the edge that is on boundary get the next vertex on this boundary
            // If the edge starting at p_in_next is on boundary, new_vertex is its next
            // If the edge ending at p_in_next is on boundary, new vertex is its prev
            next_in_next =
                e0_on_boundary ?
                    next_facet_vertex_index( next_f, v_in_next ) : prev_v_in_next ;
        } else if( nb_around == 1 ) {
            // V must be in two border edges of facet f
            // Get the id in the facet of the vertex neighbor of v1 that is not v0
            v_in_next = v ;
            if( prev_facet_vertex_index( f, v ) == from ) {
                ringmesh_assert( is_on_border( f, v ) ) ;
                next_in_next = next_facet_vertex_index( f, v ) ;
            } else {
                ringmesh_assert( is_on_border( f, prev_facet_vertex_index( f, v ) ) ) ;
                next_in_next = prev_facet_vertex_index( f, v ) ;
            }
        }
    }

    /*!
     * @brief Get the next edge on the border
     * @param[in] f Input facet index
     * @param[in] e Edge index in the facet
     * @param[out] next_f Next facet index
     * @param[out] next_e Next edge index in the facet
     */
    void Surface::next_on_border(
        index_t f,
        index_t e,
        index_t& next_f,
        index_t& next_e ) const
    {
        index_t v = next_facet_vertex_index( f, e ) ;
        index_t next_in_next( NO_ID ) ;
        return next_on_border( f, e, v, next_f, next_e, next_in_next ) ;
    }

    /*!
     * @brief Get the first facet of the surface that has an edge linking the two vertices (ids in the surface)
     *
     * @param[in] in0 Index of the first vertex in the surface
     * @param[in] in1 Index of the second vertex in the surface
     * @return NO_ID or the index of the facet
     */
    index_t Surface::facet_from_surface_vertex_ids( index_t in0, index_t in1 ) const
    {
        ringmesh_assert(
            in0 < nb_vertices() && in1 <nb_vertices() ) ;

        // Another possible, probably faster, algorithm is to check if the 2 indices
        // are neighbors in facets_ and check that they are in the same facet

        // Check if the edge is in one of the facet
        for( index_t f = 0; f < nb_polytope(); ++f ) {
            bool found = false ;
            index_t prev = polytope_vertex_index( f,
                nb_polytope_vertices( f ) - 1 ) ;
            for( index_t v = 0; v < nb_polytope_vertices( f ); ++v ) {
                index_t p = polytope_vertex_index( f, v ) ;
                if( ( prev == in0 && p == in1 ) || ( prev == in1 && p == in0 ) ) {
                    found = true ;
                    break ;
                }
                prev = p ;
            }
            if( found ) {
                return f ;
            }
        }
        return NO_ID ;
    }

    /*!
     * @brief Get the first facet of the surface that has an edge linking the
     * two vertices (ids in the model)
     *
     * @param[in] i0 Index of the first vertex in the model
     * @param[in] i1 Index of the second vertex in the model
     * @return NO_ID or the index of the facet
     */
    index_t Surface::facet_from_model_vertex_ids( index_t i0, index_t i1 ) const
    {
        index_t facet = NO_ID ;
        index_t edge = NO_ID ;
        edge_from_model_vertex_ids( i0, i1, facet, edge ) ;
        return facet ;
    }

    /*!
     * @brief Determine the facet and the edge in this facet linking the 2 vertices
     * @details There might be two pairs facet-edge. This only gets the first.
     *
     * @param[in] i0 First vertex index in the model
     * @param[in] i1 Second vertex index in the model
     * @param[out] facet NO_ID or facet index in the surface
     * @param[out] edge NO_ID or edge index in the facet
     */
    void Surface::edge_from_model_vertex_ids(
        index_t i0,
        index_t i1,
        index_t& facet,
        index_t& edge ) const
    {
        edge = NO_ID ;

        // If a facet is given, look for the edge in this facet only
        if( facet != NO_ID ) {
            for( index_t v = 0; v < nb_polytope_vertices( facet ); ++v ) {
                index_t prev = model_vertex_id( facet,
                    prev_facet_vertex_index( facet, v ) ) ;
                index_t p = model_vertex_id( facet, v ) ;
                if( ( prev == i0 && p == i1 ) || ( prev == i1 && p == i0 ) ) {
                    edge = prev_facet_vertex_index( facet, v ) ;
                    return ;
                }
            }
        } else {
            for( index_t f = 0; f < nb_polytope(); ++f ) {
                facet = f ;
                edge_from_model_vertex_ids( i0, i1, facet, edge ) ;
                if( edge != NO_ID ) {
                    return ;
                }
            }
        }

        // If we get here, no facet was found get out
        facet = NO_ID ;
        edge = NO_ID ;
    }

    /*!
     * @brief Determine the facet and the edge linking the 2 vertices with the same orientation
     * @details There might be two pairs facet-edge. This only gets the first.
     *
     * @param[in] i0 First vertex index in the model
     * @param[in] i1 Second vertex index in the model
     * @param[out] facet NO_ID or facet index in the surface
     * @param[out] edge NO_ID or edge index in the facet
     */
    void Surface::oriented_edge_from_model_vertex_ids(
        index_t i0,
        index_t i1,
        index_t& facet,
        index_t& edge ) const
    {
        // Copy from above .. tant pis
        edge = NO_ID ;

        // If a facet is given, look for the oriented edge in this facet only
        if( facet != NO_ID ) {
            for( index_t v = 0; v < nb_polytope_vertices( facet ); ++v ) {
                index_t p = model_vertex_id( facet, v ) ;
                index_t next = model_vertex_id( facet,
                    next_facet_vertex_index( facet, v ) ) ;

                if( p == i0 && next == i1 ) {
                    edge = v ;
                    return ;
                }
            }
        } else {
            for( index_t f = 0; f < nb_polytope(); ++f ) {
                facet = f ;
                oriented_edge_from_model_vertex_ids( i0, i1, facet, edge ) ;
                if( edge != NO_ID ) {
                    return ;
                }
            }
        }
        facet = NO_ID ;
    }

    /*!
     * @brief Comparator of two vec3
     */
    struct comp_vec3bis {
        bool operator()( const vec3& l, const vec3& r ) const
        {
            if( l.x != r.x ) {
                return l.x < r.x ;
            }
            if( l.y != r.y ) {
                return l.y < r.y ;
            }
            return l.z < r.z ;
        }
    } ;

    /*!
     * @brief Determines the facets around a vertex
     *
     * @param[in] v Index ot the vertex in the surface
     * @param[in] result Indices of the facets containing @param v
     * @param[in] border_only If true only facets on the border are considered
     * @return The number of facet found
     */
    index_t Surface::facets_around_vertex(
        index_t v,
        std::vector< index_t >& result,
        bool border_only ) const
    {
        index_t f = NO_ID ;

        // I tried using an AABB tree to accelerate the function
        // but apparently this does not do exactly the same than brute force
        // I do not understand why (JP)
        // I have problem with closed line in model A6. No idea why !! (JP)

        /*   // What should be an adequate limit on the number of
         // facets under which we do not use the AABB tree ?
         // When building an AABB tree the Mesh is triangulated
         // We do not want that
         if( mesh().facets.are_simplices() && mesh().facets.nb() > 10 ) {
         double dist = DBL_MAX ;
         vec3 nearest ;
         f = tools.aabb().nearest_facet( vertex( v ), nearest, dist ) ;
         // Check that the point is indeed a vertex of the facet
         if( facet_vertex_id( f, v ) == NO_ID ) {
         f = NO_ID ;
         }
         }
         */
        // So, we are back to the brute force stupid approach             
        for( index_t i = 0; i < nb_polytope(); ++i ) {
            for( index_t lv = 0; lv < nb_polytope_vertices( i ); lv++ ) {
                if( polytope_vertex_index( i, lv ) == v ) {
                    f = i ;
                    break ;
                }
            }
        }
        return facets_around_vertex( v, result, border_only, f ) ;
    }

    /*!
     * @brief Determines the facets around a vertex
     *
     * @param[in] P Index ot the vertex in the surface
     * @param[in] result Indices of the facets containing @param P
     * @param[in] border_only If true only facets on the border are considered
     * @param[in] f0 Index of one facet containing the vertex @param P
     * @return The number of facet found
     *
     * @todo Evaluate if this is fast enough !!
     */
    index_t Surface::facets_around_vertex(
        index_t P,
        std::vector< index_t >& result,
        bool border_only,
        index_t f0 ) const
    {
        result.resize( 0 ) ;

        if( f0 == NO_ID ) {
            return 0 ;
        }

        // Flag the visited facets
        std::vector< index_t > visited ;
        visited.reserve( 10 ) ;

        // Stack of the adjacent facets
        std::stack< index_t > S ;
        S.push( f0 ) ;
        visited.push_back( f0 ) ;

        do {
            index_t f = S.top() ;
            S.pop() ;

            for( index_t v = 0; v < nb_polytope_vertices( f ); ++v ) {
                if( polytope_vertex_index( f, v ) == P ) {
                    index_t adj_P = facet_adjacent_index( f, v ) ;
                    index_t prev = prev_facet_vertex_index( f, v ) ;
                    index_t adj_prev = facet_adjacent_index( f, prev ) ;

                    if( adj_P != NO_ADJACENT ) {
                        // The edge starting at P is not on the boundary
                        if( !contains( visited, adj_P ) ) {
                            S.push( adj_P ) ;
                            visited.push_back( adj_P ) ;
                        }
                    }
                    if( adj_prev != NO_ADJACENT ) {
                        // The edge ending at P is not on the boundary
                        if( !contains( visited, adj_prev ) ) {
                            S.push( adj_prev ) ;
                            visited.push_back( adj_prev ) ;
                        }
                    }

                    if( border_only ) {
                        if( adj_P == NO_ADJACENT || adj_prev == NO_ADJACENT ) {
                            result.push_back( f ) ;
                        }
                    } else {
                        result.push_back( f ) ;
                    }

                    // We are done with this facet
                    break ;
                }
            }
        } while( !S.empty() ) ;

        return static_cast< index_t >( result.size() ) ;
    }

    /*!
     * @brief Compute closest vertex in a facet to a point
     * @param[in] f Facet index
     * @param[in] v Coordinates of the point to which distance is measured
     * @return Index of the vertex of @param f closest to @param v
     */
    index_t Surface::closest_vertex_in_facet( index_t f, const vec3& v ) const
    {
        index_t result = 0 ;
        double dist = DBL_MAX ;
        for( index_t p = 0; p < nb_polytope_vertices( f ); p++ ) {
            double distance = length2( v - polytope_vertex( f, p ) ) ;
            if( dist > distance ) {
                dist = distance ;
                result = p ;
            }
        }
        return result ;
    }

    /********************************************************************/

    bool Region::is_mesh_valid() const
    {
        if( !is_meshed() ) {
            return true ;
        } else {
            bool valid = true ;
            // Check that the GEO::Mesh has the expected elements
            // at least 4 vertices and one cell.
            if( mesh_.nb_vertices() < 4 ) {
                GEO::Logger::warn( "GeoModelElement" ) << gme_id()
                    << " has less than 4 vertices " << std::endl ;
                valid = false ;
            }

            // No isolated vertices
            index_t nb_isolated_vertices = count_nb_isolated_vertices( *this ) ;
            if( nb_isolated_vertices > 0 ) {
                GEO::Logger::warn( "GeoModelElement" ) << gme_id() << " mesh has "
                    << nb_isolated_vertices << " isolated vertices " << std::endl ;
                valid = false ;
            }

            // No cell with negative volume
            // No cell incident to the same vertex check local and global indices
            index_t nb_degenerate = 0 ;
            for( index_t c = 0; c < mesh_.nb_cells(); c++ ) {
                if( cell_is_degenerate( *this, c ) ) {
                    nb_degenerate++ ;
                }
            }
            if( nb_degenerate != 0 ) {
                GEO::Logger::warn( "GeoModelElement" ) << gme_id() << " mesh has "
                    << nb_degenerate << " degenerate cells " << std::endl ;
                valid = false ;
            }

            // One connected component
            index_t cc = compute_nb_volume_connected_components( *this ) ;
            if( cc != 1 ) {
                GEO::Logger::warn( "GeoModelElement" ) << gme_id() << " mesh has "
                    << cc << " connected components " << std::endl ;
                valid = false ;
            }
            return valid ;
        }
    }
    void Region::compute_region_volumes_per_cell_type(
        double& tet_volume,
        double& pyramid_volume,
        double& prism_volume,
        double& hex_volume,
        double& poly_volume ) const
    {
        for( index_t c = 0; c < nb_polytope(); c++ ) {
            index_t nb_vertices = nb_polytope_vertices( c ) ;
            double volume = mesh_.cell_volume( c ) ;
            switch( nb_vertices ) {
                case 4:
                    tet_volume += volume ;
                    break ;
                case 5:
                    pyramid_volume += volume ;
                    break ;
                case 6:
                    prism_volume += volume ;
                    break ;
                case 8:
                    hex_volume += volume ;
                    break ;
                default:
                    poly_volume += volume ;
                    break ;
            }
        }
    }

}
