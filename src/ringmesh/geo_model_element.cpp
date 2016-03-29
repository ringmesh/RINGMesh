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

#include <geogram/mesh/mesh.h>
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

    /*!
     * @brief Count the number of times each vertex is in an edge or facet
     *
     * @param[in] M The mesh
     * @param[out] nb Resized to the number of vertices of the mesh.
     *      Number of times one vertex appear in an edge or facet of the mesh.
     */
    void count_vertex_occurences( const GEO::Mesh& M, std::vector< index_t >& nb )
    {
        nb.resize( M.vertices.nb(), 0 ) ;
        for( index_t e = 0; e < M.edges.nb(); ++e ) {
            ++nb[M.edges.vertex( e, 0 )] ;
            ++nb[M.edges.vertex( e, 1 )] ;
        }
        for( index_t f = 0; f < M.facets.nb(); ++f ) {
            for( index_t co = M.facets.corners_begin( f );
                co < M.facets.corners_end( f ); ++co ) {
                ++nb[M.facet_corners.vertex( co )] ;
            }
        }
    }

    /*!
     * @brief Returns true if the surface facet is incident twice to the same vertex
     */
    bool facet_is_degenerate( const Surface& S, index_t f )
    {
        std::vector< index_t > corners( S.nb_vertices_in_facet( f ), NO_ID ) ;
        std::vector< index_t > corners_global( S.nb_vertices_in_facet( f ), NO_ID ) ;
        index_t v = 0 ;
        for( index_t c = S.facet_begin( f ); c < S.facet_end( f ); ++c ) {
            corners[v] = c ;
            corners_global[v] = S.model_vertex_id( f, v ) ;
            v++ ;
        }
        ringmesh_assert(
            std::count( corners.begin(), corners.end(), NO_ID ) == 0 ) ;
        ringmesh_assert(
            std::count( corners_global.begin(), corners_global.end(), NO_ID ) == 0 ) ;
        // 0 is the default value of the model_vertex_id
        // If we have only 0 either this is a degenerate facets, but most certainly
        // model vertex ids are not good 
        ringmesh_assert(
            std::count( corners_global.begin(), corners_global.end(), 0 )
                != corners_global.size() ) ;

        std::sort( corners.begin(), corners.end() ) ;
        std::sort( corners_global.begin(), corners_global.end() ) ;
        return std::unique( corners.begin(), corners.end() ) != corners.end()
            || std::unique( corners_global.begin(), corners_global.end() )
                != corners_global.end() ;
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
        for( index_t f = 0; f < S.nb_cells(); f++ ) {
            out << "f" << " " ;
            for( index_t v = 0; v < S.nb_vertices_in_facet( f ); v++ ) {
                out << S.surf_vertex_id( f, v ) + 1 << " " ;
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
            bool model_has_parent_elements( model().nb_elements( parent_type( T ) ) > 0 ) ;
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
        print_bounded_attributes( mesh_ ) ;
#endif
    }

    /*!
     * @brief Binds attributes stored by the GME on the Mesh
     */
    void GeoModelMeshElement::bind_attributes()
    {
        model_vertex_id_.bind( mesh_.vertices.attributes(),
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
        if( mesh_.vertices.nb() != 1 ) {
            GEO::Logger::err( "GeoModelElement" ) << "Corner " << index()
                << " mesh has " << mesh_.vertices.nb() << " vertices " << std::endl ;
            valid = false ;
        }
        if( mesh_.edges.nb() != 0 ) {
            GEO::Logger::err( "GeoModelElement" ) << "Corner " << index()
                << " mesh has " << mesh_.edges.nb() << " edges " << std::endl ;
            valid = false ;
        }
        if( mesh_.facets.nb() != 0 ) {
            GEO::Logger::err( "GeoModelElement" ) << "Corner " << index()
                << " mesh has " << mesh_.facets.nb() << " facets " << std::endl ;
            valid = false ;
        }
        if( mesh_.cells.nb() != 0 ) {
            GEO::Logger::warn( "GeoModelElement" ) << "Corner " << index()
                << " mesh has " << mesh_.cells.nb() << " cells " << std::endl ;
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
        if( mesh_.vertices.nb() < 2 ) {
            GEO::Logger::err( "GeoModelElement" ) << "Line " << index() << " has "
                << mesh_.vertices.nb() << " vertices " << std::endl ;
            valid = false ;
        }
        if( mesh_.edges.nb() == 0 ) {
            GEO::Logger::err( "GeoModelElement" ) << "Line " << index()
                << " mesh has " << mesh_.edges.nb() << " edges " << std::endl ;
            valid = false ;
        }
        if( mesh_.facets.nb() != 0 ) {
            GEO::Logger::err( "GeoModelElement" ) << "Line " << index()
                << " mesh has " << mesh_.facets.nb() << " facets " << std::endl ;
            valid = false ;
        }
        if( mesh_.cells.nb() != 0 ) {
            GEO::Logger::err( "GeoModelElement" ) << "Line " << index()
                << " mesh has " << mesh_.cells.nb() << " cells " << std::endl ;
            valid = false ;
        }

        // Model indices must be valid
        valid = check_range_model_vertex_ids( *this ) && valid ;

        if( mesh_.vertices.nb() > 1 ) {
            // Count the number of edges in which each vertex is
            std::vector< index_t > nb ;
            count_vertex_occurences( mesh(), nb ) ;
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
        for( index_t e = 0; e < nb_cells(); ++e ) {
            double l = length( vertex( e, 1 ) - vertex( e, 0 ) ) ;
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
        if( mesh_.vertices.nb() < 3 ) {
            GEO::Logger::warn( "GeoModelElement" ) << gme_id()
                << " has less than 3 vertices " << std::endl ;
            valid = false ;
        }
        // Is it important to have edges or not ?
        // I would say we do not care (JP) - so no check on that 
        if( mesh_.facets.nb() == 0 ) {
            GEO::Logger::warn( "GeoModelElement" ) << gme_id() << " has no facets "
                << std::endl ;
            valid = false ;
        }
        if( mesh_.cells.nb() != 0 ) {
            GEO::Logger::warn( "GeoModelElement" ) << gme_id() << " has "
                << mesh_.cells.nb() << " cells " << std::endl ;
            valid = false ;
        }

        // No isolated vertices
        std::vector< index_t > nb ;
        count_vertex_occurences( mesh(), nb ) ;
        index_t nb0 = static_cast< index_t >( std::count( nb.begin(), nb.end(), 0 ) ) ;
        if( nb0 > 0 ) {
            GEO::Logger::warn( "GeoModelElement" ) << gme_id() << " mesh has " << nb0
                << " isolated vertices " << std::endl ;
            valid = false ;
        }

        // No zero area facet
        // No facet incident to the same vertex check local and global indices
        index_t nb_degenerate = 0 ;
        for( index_t f = 0; f < mesh_.facets.nb(); f++ ) {
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
        index_t cc = GEO::mesh_nb_connected_components( mesh_ ) ;
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
        ringmesh_assert( v < nb_vertices_in_facet( f ) ) ;
        ringmesh_assert( is_on_border( f, v ) || is_on_border( f, from ) ) ;

        index_t V = surf_vertex_id( f, v ) ;

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
            v_in_next = facet_vertex_id( next_f, V ) ;
            ringmesh_assert( v_in_next != NO_ID ) ;

            // The edges containing V in next_f are
            // the edge starting at v_in_next and the one ending there
            index_t prev_v_in_next = prev_in_facet( next_f, v_in_next ) ;

            bool e0_on_boundary = is_on_border( next_f, v_in_next ) ;

            // Only one must be on the boundary otherwise there is a corner missing
            ringmesh_assert( e0_on_boundary != is_on_border( next_f, prev_v_in_next ) ) ;

            // From the edge that is on boundary get the next vertex on this boundary
            // If the edge starting at p_in_next is on boundary, new_vertex is its next
            // If the edge ending at p_in_next is on boundary, new vertex is its prev
            next_in_next = e0_on_boundary ?
                    next_in_facet( next_f, v_in_next ) : prev_v_in_next ;
        } else if( nb_around == 1 ) {
            // V must be in two border edges of facet f
            // Get the id in the facet of the vertex neighbor of v1 that is not v0
            v_in_next = v ;
            if( prev_in_facet( f, v ) == from ) {
                ringmesh_assert( is_on_border( f, v ) ) ;
                next_in_next = next_in_facet( f, v ) ;
            } else {
                ringmesh_assert( is_on_border( f, prev_in_facet( f, v ) ) ) ;
                next_in_next = prev_in_facet( f, v ) ;
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
        index_t v = next_in_facet( f, e ) ;
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
            in0 < mesh_.vertices.nb() && in1 < mesh_.vertices.nb() ) ;

        // Another possible, probably faster, algorithm is to check if the 2 indices
        // are neighbors in facets_ and check that they are in the same facet

        // Check if the edge is in one of the facet
        for( index_t f = 0; f < nb_cells(); ++f ) {
            bool found = false ;
            index_t prev = surf_vertex_id( f, nb_vertices_in_facet( f ) - 1 ) ;
            for( index_t v = 0; v < nb_vertices_in_facet( f ); ++v ) {
                index_t p = surf_vertex_id( f, v ) ;
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
            for( index_t v = 0; v < nb_vertices_in_facet( facet ); ++v ) {
                index_t prev = model_vertex_id( facet, prev_in_facet( facet, v ) ) ;
                index_t p = model_vertex_id( facet, v ) ;
                if( ( prev == i0 && p == i1 ) || ( prev == i1 && p == i0 ) ) {
                    edge = prev_in_facet( facet, v ) ;
                    return ;
                }
            }
        } else {
            for( index_t f = 0; f < nb_cells(); ++f ) {
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
            for( index_t v = 0; v < nb_vertices_in_facet( facet ); ++v ) {
                index_t p = model_vertex_id( facet, v ) ;
                index_t next = model_vertex_id( facet, next_in_facet( facet, v ) ) ;

                if( p == i0 && next == i1 ) {
                    edge = v ;
                    return ;
                }
            }
        } else {
            for( index_t f = 0; f < nb_cells(); ++f ) {
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
     * @brief Convert vertex surface index to an index in a facet
     * @param[in] f Index of the facet
     * @param[in] surf_vertex_id_in Index of the vertex in the surface
     * @return NO_ID or index of the vertex in the facet
     */
    index_t Surface::facet_vertex_id( index_t f, index_t surf_vertex_id_in ) const
    {
        for( index_t v = 0; v < nb_vertices_in_facet( f ); v++ ) {
            if( surf_vertex_id( f, v ) == surf_vertex_id_in ) {
                return v ;
            }
        }
        return NO_ID ;
    }

    /*!
     * @brief Convert model vertex index to an index in a facet
     * @param[in] f Index of the facet
     * @param[in] model_v_id Index of the vertex in the GeoModel
     * @return NO_ID or index of the vertex in the facet
     */
    index_t Surface::facet_id_from_model( index_t f, index_t model_v_id ) const
    {
        for( index_t v = 0; v < nb_vertices_in_facet( f ); v++ ) {
            if( model_vertex_id( f, v ) == model_v_id ) {
                return v ;
            }
        }
        return NO_ID ;
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
        for( index_t i = 0; i < nb_cells(); ++i ) {
            for( index_t lv = 0; lv < nb_vertices_in_facet( i ); lv++ ) {
                if( surf_vertex_id( i, lv ) == v ) {
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

            for( index_t v = 0; v < nb_vertices_in_facet( f ); ++v ) {
                if( surf_vertex_id( f, v ) == P ) {
                    index_t adj_P = adjacent( f, v ) ;
                    index_t prev = prev_in_facet( f, v ) ;
                    index_t adj_prev = adjacent( f, prev ) ;

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

        return static_cast<index_t> ( result.size() ) ;
    }

    /*!
     * Get the facet normal
     * @param[in] f Facet index
     * @return Normal to the facet
     */
    vec3 Surface::facet_normal( index_t facet_index ) const
    {
        return normalize( GEO::Geom::mesh_facet_normal( mesh_, facet_index ) ) ;
    }

    vec3 Surface::facet_barycenter( index_t facet_index ) const
    {
        return GEO::Geom::mesh_facet_center( mesh_, facet_index ) ;
    }

    double Surface::facet_area( index_t facet_index ) const
    {
        return GEO::Geom::mesh_facet_area( mesh_, facet_index ) ;
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
        for( index_t p = 0; p < nb_vertices_in_facet( f ); p++ ) {
            double distance = length2( v - vertex( f, p ) ) ;
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
            GEO::Logger::warn( "GeoModel" )
                << "TO DO : Mesh validity function on Regions is to implement "
                << std::endl ;
            return true ;
        }
    }

    /********************************************************************/

    SurfaceTools::SurfaceTools( const Surface& surface )
        : surface_( surface ), aabb_( nil ), ann_( nil )
    {
    }

    SurfaceTools::~SurfaceTools()
    {
        if( aabb_ ) delete aabb_ ;
        if( ann_ ) delete ann_ ;
    }

    /*! 
     * @brief Create an AABB tree for a Surface
     * @pre The Surface mesh must be simplicial
     * @warning SIDE EFFECTS: The Surface mesh vertices are reordered.
     * That is why the global Mesh vertices are deleted (This is BAD)
     */
    const GEO::MeshFacetsAABB& SurfaceTools::aabb() const
    {
        GeoModel& M = const_cast< GeoModel& >( surface_.model() ) ;
        if( M.mesh.vertices.is_initialized() ) {
            GEO::Logger::warn( "AABB" )
                << "Creation of AABB results in deletion of the GeoModelMeshVertices"
                << std::endl ;
            M.mesh.vertices.clear() ;
        }
        if( aabb_ == nil ) {
            // Geogram triangulates the Mesh when creating the AABB tree
            ringmesh_assert( surface_.mesh().facets.are_simplices() ) ;

            // Very bad side effect
            // The root cause of the problem is the duplication of many things 
            // in our GeoModel structure [JP]
            M.mesh.vertices.clear() ;

            aabb_ = new GEO::MeshFacetsAABB( surface_.mesh() ) ;
        }
        return *aabb_ ;
    }

    const ColocaterANN& SurfaceTools::ann() const
    {
        if( ann_ == nil ) {
            ann_ = new ColocaterANN( surface_.mesh(), ColocaterANN::VERTICES ) ;
        }
        return *ann_ ;
    }

    /********************************************************************/

    RegionTools::RegionTools( const Region& region )
        : region_( region ), aabb_( nil ), ann_( nil )
    {
    }

    RegionTools::~RegionTools()
    {
        if( aabb_ ) delete aabb_ ;
        if( ann_ ) delete ann_ ;
    }

    const GEO::MeshCellsAABB& RegionTools::aabb() const
    {
        GeoModel& M = const_cast< GeoModel& >( region_.model() ) ;
        if( M.mesh.vertices.is_initialized() ) {
            GEO::Logger::warn( "AABB" )
                << "Creation of AABB results in deletion of the GeoModelMeshVertices"
                << std::endl ;
            M.mesh.vertices.clear() ;
        }
        if( aabb_ == nil ) {
            aabb_ = new GEO::MeshCellsAABB( region_.mesh() ) ;
            /// @todo Et pourquoi creer AABB me fait vider les sommets ?
            /// @todo Il faut un mecanisme update de ces RegionTools correct.
            // if( ann_ ) {
            //     delete ann_ ;
            //     this_not_const->ann_ = nil ;
            // }

            // Building an AABB reorders the mesh vertices and facets
            // Very annoying if model_vertex_ids are set because we need
            // to update the model vertices.
            /* GeoModel& M = const_cast< GeoModel& >( region_.model() ) ;
             if( M.mesh.vertices.is_initialized() ) {
             for( index_t sv = 0; sv < region_.nb_vertices(); ++sv ) {
             index_t v = region_.model_vertex_id( sv ) ;
             const std::vector< GMEVertex >& to_update =
             M.mesh.vertices.gme_vertices( v ) ;

             index_t count_skipped = 0 ;
             for( index_t i = 0; i < to_update.size(); ++i ) {
             if( to_update[i].gme_id == region_.gme_id() ) {
             M.mesh.vertices.set_gme( v, i,
             GMEVertex( region_.gme_id(), sv ) ) ;
             break ;
             }
             }
             }
             } */
        }
        return *aabb_ ;
    }

    const ColocaterANN& RegionTools::ann() const
    {
        if( ann_ == nil ) {
            ann_ = new ColocaterANN( region_.mesh(), ColocaterANN::VERTICES ) ;
        }
        return *ann_ ;
    }
}
