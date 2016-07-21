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
 * @file Implementation of all GeoModelEntities classes
 * @author Jeanne Pellerin and Arnaud Botella 
 */

#include <ringmesh/geo_model_entity.h>

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

    typedef GeoModelEntity::gme_t gme_t ;
    typedef GeoModelMeshEntity GMME ;
    typedef GeoModelGeologicalEntity GMGE ;

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
                Logger::warn( "GeoModelEntity" )
                    << "Invalid model vertex index in " << E.gme_id() << std::endl ;
                return false ;
            }
        }
        return true ;
    }

    index_t compute_nb_volume_connected_components( const Region& M )
    {
        static const index_t NO_COMPONENT = index_t( -1 ) ;
        std::vector< index_t > component( M.nb_mesh_elements(), NO_COMPONENT ) ;
        index_t nb_components = 0 ;
        for( index_t cell = 0; cell < M.nb_mesh_elements(); cell++ ) {
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
     * @param[in] gmme The GeoModelMeshEntity
     * @param[out] nb Resized to the number of vertices of the mesh.
     *      Number of times one vertex appear in an mesh_element collection of 
     *      the GeoModelMeshEntity edge or facet of the mesh.
     */
    void count_vertex_occurences(
        const GeoModelMeshEntity& E,
        std::vector< index_t >& nb )
    {
        nb.resize( E.nb_vertices(), 0 ) ;
        for( index_t mesh_element_index = 0; mesh_element_index < E.nb_mesh_elements();
            ++mesh_element_index ) {
            for( index_t vertex = 0;
                vertex < E.nb_mesh_element_vertices( mesh_element_index ); ++vertex ) {
                ++nb[E.mesh_element_vertex_index( mesh_element_index, vertex )] ;
            }
        }
    }

    index_t count_nb_isolated_vertices( const GeoModelMeshEntity& mesh )
    {
        std::vector< index_t > nb ;
        count_vertex_occurences( mesh, nb ) ;
        return static_cast< index_t >( std::count( nb.begin(), nb.end(), 0 ) ) ;
    }

    bool check_mesh_entity_vertices_are_different(
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
        index_t nb_facet_vertices = S.nb_mesh_element_vertices( f ) ;
        std::vector< index_t > corners( nb_facet_vertices, NO_ID ) ;
        std::vector< index_t > corners_global( nb_facet_vertices, NO_ID ) ;
        index_t v = 0 ;
        for( index_t c = S.facet_begin( f ); c < S.facet_end( f ); ++c ) {
            corners[v] = c ;
            corners_global[v] = S.model_vertex_id( f, v ) ;
            v++ ;
        }
        return check_mesh_entity_vertices_are_different( corners, corners_global ) ;
    }    
    /*!
     * @brief Returns true if the region cell is incident twice to the same vertex
     * or if the cell volume is negative or inferior to epsilon
     */
    bool cell_is_degenerate( const Region& region, index_t cell_index )
    {
        index_t nb_vertices_in_cell = region.nb_mesh_element_vertices( cell_index ) ;
        std::vector< index_t > vertices( nb_vertices_in_cell, NO_ID ) ;
        std::vector< index_t > vertices_global( nb_vertices_in_cell, NO_ID ) ;
        for( index_t v = 0; v < nb_vertices_in_cell; v++ ) {
            vertices[v] = region.mesh_element_vertex_index( cell_index, v ) ;
            vertices_global[v] = region.model_vertex_id( cell_index, v ) ;
        }
        double volume = region.mesh_element_size( cell_index ) ;
        return check_mesh_entity_vertices_are_different( vertices, vertices_global )
            || volume < epsilon ;
	}

    /*!
     * @brief Debug: Save a Surface of the model in the file OBJ format is used
     * @todo Move this function to an API providing utility functions on a
     * GeoModel and its Entities ? [JP]
     */
    void save_surface_as_obj_file( const Surface& S, const std::string& file_name )
    {
        std::ofstream out( file_name.c_str() ) ;
        if( out.bad() ) {
            Logger::err( "I/O" ) << "Error when opening the file: "
                << file_name.c_str() << std::endl ;
            return ;
        }
        out.precision( 16 ) ;
        for( index_t p = 0; p < S.nb_vertices(); p++ ) {
            const vec3& V = S.vertex( p ) ;
            out << "v" << " " << V.x << " " << V.y << " " << V.z << std::endl ;
        }
        for( index_t f = 0; f < S.nb_mesh_elements(); f++ ) {
            out << "f" << " " ;
            for( index_t v = 0; v < S.nb_mesh_element_vertices( f ); v++ ) {
                out << S.mesh_element_vertex_index( f, v ) + 1 << " " ;
            }
            out << std::endl ;
        }
    }
}

namespace RINGMesh {
    Universe::Universe( const GeoModel& model )
        : GeoModelEntity( model, NO_ID, type_name_ )
    {
    }

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
     * @todo Add other types of unconformity, see RINGMesh::GeoModelEntity::TYPE. --GC
     */
    GeoModelEntity::GEOL_FEATURE GeoModelEntity::determine_geological_type(
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
     * \return the (lowercase) string associated to a
     * GeoModelELement::GEOL_FEATURE
     */
    std::string GeoModelEntity::geol_name( GME::GEOL_FEATURE t )
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
     * @brief Return true if \param type is a CORNER, LINE or SURFACE
     *
    bool GeoModelEntity::has_mesh( GME::TYPE type )
    {
        return type <= REGION ;
    } */

    bool GeoModelEntity::is_connectivity_valid() const
    {
        bool valid = true ;

/*        /// 1. Check the validity of identification information
        ///    in the model - Universe has no index, but a TYPE
        if( gme_id() == gme_t() ) {
            Logger::err( "GeoModelEntity" ) << " Entity associated to model "
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

        if( index() >= model().nb_entities( type() ) ) {
            Logger::warn( "GeoModelEntity" ) << " Entity index " << gme_id()
                << " is not valid " << " There are " << model().nb_entities( type() )
                << " entity of that type in model " << model().name() << std::endl ;
            // This really should not happen
            valid = false ;
            ringmesh_assert( valid ) ;
            return valid ;
        }
        if( &( model().entity( gme_id() ) ) != this ) {
            Logger::err( "GeoModelEntity" ) << " Entity " << gme_id()
                << " in model " << model().name() << " does not match this entity"
                << std::endl ;
            // This really should not happen
            ringmesh_assert( valid ) ;
            valid = false ;
            return valid ;
        }

        /// 2. Check that required information for the TYPE is defined
        ///    and that reverse information is stored by the corresponding
        ///    entities
        TYPE T = type() ;

        // Boundaries
        if( boundary_allowed( T ) ) {
            if( T == REGION ) {
                if( nb_boundaries() == 0 ) {
                    Logger::warn( "GeoModelEntity" ) << gme_id()
                        << " has no boundaries " << std::endl ;
                    valid = false ;
                }
            }
            // A Line must have 2 corners - they are identical if the Line is closed
            if( T == LINE && nb_boundaries() != 2 ) {
                Logger::warn( "GeoModelEntity" ) << gme_id()
                    << " does not have 2 corners" << std::endl ;
                valid = false ;
            }
            // No requirement on Surface - it may have no boundary - bubble

            // All entities in the boundary must have this in their
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
                    Logger::warn( "GeoModelEntity" )
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
                Logger::warn( "GeoModelEntity" ) << gme_id()
                    << " is in the boundary of no entity " << std::endl ;
                valid = false ;
            }

            // All entities in the in_boundary must have this in their
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
                    Logger::warn( "GeoModelEntity" )
                        << "Inconsistency in_boundary-boundary between " << gme_id()
                        << " and " << E.gme_id() << std::endl ;
                    valid = false ;
                }
            }
        }

        // Parent - High level entities are not mandatory
        // But if the model has entities of the parent type, the entity must have a parent
        if( parent_allowed( T ) ) {
            bool model_has_parent_entities(
                model().nb_entities( parent_type( T ) ) > 0 ) ;
            if( model_has_parent_entities ) {
                if( has_parent() ) {
                    const GME& E = parent() ;
                    // The parent must have this entity in its children
                    bool found = false ;
                    index_t j = 0 ;
                    while( !found && j < E.nb_children() ) {
                        if( E.child_id( j ) == gme_id() ) {
                            found = true ;
                        }
                        j++ ;
                    }
                    if( !found ) {
                        Logger::warn( "GeoModelEntity" )
                            << "Inconsistency parent-child between " << gme_id()
                            << " and " << E.gme_id() << std::endl ;
                        valid = false ;
                    }
                } else {
                    Logger::warn( "GeoModelEntity" ) << gme_id()
                        << " has no geological parent entity " << std::endl ;
                    valid = false ;
                }
            }
        }

        // Children
        if( child_allowed( T ) ) {
            if( nb_children() == 0 ) {
                Logger::warn( "GeoModelEntity" ) << gme_id()
                    << " has no children mesh entity, so no geometry "
                    << std::endl ;
                valid = false ;
            }

            // All children must have this entity as a parent
            for( index_t i = 0; i < nb_children(); ++i ) {
                if( child( i ).parent_id() != gme_id() ) {
                    Logger::warn( "GeoModelEntity" )
                        << "Inconsistency child-parent between " << gme_id()
                        << " and " << child( i ).gme_id() << std::endl ;
                    valid = false ;
                }
            }
        }
        return valid ;
   */

        return false ; // TO SPLIT between classes
    }

}
