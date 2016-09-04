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
 * @file Implementation of all GeoModelGeologicalEntities classes
 * @author Jeanne Pellerin and Arnaud Botella 
 */

#include <ringmesh/geomodel/geo_model_geological_entity.h>

#include <geogram/basic/algorithm.h> // for validity - that may be moved 

#include <ringmesh/geomodel/geo_model.h>
#include <ringmesh/geomodel/geo_model_mesh_entity.h>
#include <ringmesh/io/io.h>


namespace {
    using namespace RINGMesh ;

    bool check_has_children( const GeoModelGeologicalEntity& E )
    {
        if( E.nb_children() == 0 ) {
            Logger::warn( "GeoModel" ) << E.type_name() << " " << E.index()
                << " has no children" << std::endl ;
            return false ;
        }
        return true ;
    }
}

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
    bool GeoModelGeologicalEntity::is_index_valid() const
    {
    	return index() < model().nb_geological_entities( type_name() ) ;
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

    bool Contact::is_valid() const
    {
        return check_has_children( *this ) ;
    }

    const std::string Interface::child_type_name() const
    {
        return Surface::type_name_static() ;
    }

    bool Interface::is_valid() const
    {
        return check_has_children( *this ) ;
    }

    const std::string Layer::child_type_name() const
    {
        return Region::type_name_static() ;
    }

    bool Layer::is_valid() const
    {
        return check_has_children( *this ) ;
    }
      
}
