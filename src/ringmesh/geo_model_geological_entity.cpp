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

#include <ringmesh/geo_model_geological_entity.h>

#include <geogram/basic/algorithm.h> // for validity - that may be moved 

#include <ringmesh/io.h>
#include <ringmesh/geo_model.h>
#include <ringmesh/geo_model_mesh_entity.h>


namespace RINGMesh {
    
    typedef std::string EntityType ;

    const GeoModelMeshEntity& GeoModelGeologicalEntity::child( index_t x ) const
    {
        return model().mesh_entity( child_gme( x ) ) ;
    }

    bool GeoModelGeologicalEntity::is_on_voi() const
    {
        for( index_t i = 0; i < nb_children(); i++ ) {
            if( !child( i ).is_on_voi() ) return false ;
        }
        return true ;
    }

    bool GeoModelGeologicalEntity::is_connectivity_valid() const
    {
        bool valid = true ;
        if( nb_children() == 0 ) {
            Logger::warn( "GeologicalEntity" ) << gme_id()
                << " is undefined. No child. "
                << std::endl ;
            valid = false ;
        } else {
            // All children must have this entity as a parent
            const EntityType entity_type = type_name() ;
            for( index_t i = 0; i < nb_children(); ++i ) {
                const GeoModelMeshEntity& one_child = child( i ) ;
                if( one_child.parent_gme( entity_type ) != gme_id() ) {
                    Logger::warn( "GeoModelEntity" )
                        << "Inconsistency child-parent between " << gme_id()
                        << " and " << one_child.gme_id() << std::endl ;
                    valid = false ;
                }
            }
        }
        return valid ;
    }

    void GeoModelGeologicalEntity::initialize()
    {
        ringmesh_register_GeoModelGeologicalEntity_creator( Contact ) ;
        ringmesh_register_GeoModelGeologicalEntity_creator( Interface ) ;
        ringmesh_register_GeoModelGeologicalEntity_creator( Layer ) ;
    }

    const std::string Contact::child_type_name() const
    {
        return Line::type_name_static() ;
    }

    const std::string Interface::child_type_name() const
    {
        return Surface::type_name_static() ;
    }

    const std::string Layer::child_type_name() const
    {
        return Region::type_name_static() ;
    }

      
    /*!
     * @brief Get the entities in the boundary of which @param E is
     * @details For GMME, get the contents of the in_boundary vector
     *          For high level entities, determine in_boundary high level entities
     */
    void in_boundary_gme(
        const GeoModelGeologicalEntity& E,
        std::vector< gme_t >& in_boundary )
    {
        in_boundary.clear() ;

        // We are dealing with high level entities
        // Need to go through the children to get information
        for( index_t i = 0; i < E.nb_children(); ++i ) {
            for( index_t j = 0; j < E.child( i ).nb_in_boundary(); ++j ) {
                in_boundary.push_back(
                    E.child( i ).in_boundary( j ).parent_gme(
                        Layer::type_name_static() ) ) ;
            }
        }
        // Remove duplicates
        GEO::sort_unique( in_boundary ) ;
    }
    bool is_geomodel_geology_valid( const GeoModel& GM )
    {
        bool valid = true ;
        for( index_t l = 0; l < GM.nb_lines(); ++l ) {
            if( GM.line( l ).nb_in_boundary() == 1 ) {
                const GME& S = GM.line( l ).in_boundary( 0 ) ;
                if( !GME::is_fault( S.geological_feature() ) ) {
                    Logger::warn( "GeoModel" ) << " Invalid free border: "
                        << GM.line( l ).gme_id() << " is in the boundary of Surface "
                        << S.gme_id() << " that is not a FAULT " << std::endl
                        << std::endl ;
                    valid = false ;
                }
            }
        }

        for( index_t i = 0; i < GM.nb_geological_entities( Interface::type_name_static() ); ++i ) {
            std::vector< gme_t > layers ;
            const GeoModelGeologicalEntity& entity = GM.geological_entity(
                Interface::type_name_static(), i ) ;
            in_boundary_gme( entity, layers ) ;
            if( layers.empty() ) {
                Logger::warn( "GeoModel" ) << " Invalid interface: "
                    << entity.gme_id()
                    << " is in the boundary of no Layer " << std::endl ;
                valid = false ;
            }
            if( entity.geological_feature() == GME::STRATI
                && layers.size() > 2 ) {
                Logger::warn( "GeoModel" ) << " Invalid horizon: "
                    << entity.gme_id() << " is in the boundary of "
                    << layers.size() << " Layers: " ;
                for( index_t j = 0; j < layers.size(); ++j ) {
                    Logger::warn( "GeoModel" ) << layers[ j ] << " ; " ;
                }
                Logger::warn( "GeoModel" ) << std::endl ;
                valid = false ;
            }
        }
        return valid ;
    }

}
