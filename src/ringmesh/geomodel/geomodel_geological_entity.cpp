/*
 * Copyright (c) 2012-2017, Association Scientifique pour la Geologie et ses Applications (ASGA)
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

#include <ringmesh/geomodel/geomodel_geological_entity.h>

#include <ringmesh/geomodel/geomodel.h>
#include <ringmesh/geomodel/geomodel_mesh_entity.h>

namespace {
    using namespace RINGMesh;

    bool check_has_children( const GeoModelGeologicalEntity& E )
    {
        if( E.nb_children() == 0 ) {
            Logger::warn( "GeoModel", E.type_name(), " ", E.index(),
                " has no children" );
            return false;
        }
        return true;
    }
}

namespace RINGMesh {

    GeoModelGeologicalEntity::GEOL_FEATURE GeoModelGeologicalEntity::determine_geological_type(
        const std::string& in )
    {
        if( in == "reverse_fault" ) {
            return GEOL_FEATURE::REVERSE_FAULT;
        } else if( in == "normal_fault" ) {
            return GEOL_FEATURE::NORMAL_FAULT;
        } else if( in == "fault" ) {
            return GEOL_FEATURE::FAULT;
        } else if( in == "top" ) {
            return GEOL_FEATURE::STRATI;
        } else if( in == "none" ) {
            // This might seem strange - but it seems that what's
            // Gocad is doing
            return GEOL_FEATURE::STRATI;
        } else if( in == "topographic" ) {
            return GEOL_FEATURE::STRATI;
        } else if( in == "unconformity" ) {
            return GEOL_FEATURE::UNCONFORMITY;
        } else if( in == "boundary" ) {
            return GEOL_FEATURE::VOI;
        } else {
            // Default case - no information
            return GEOL_FEATURE::NO_GEOL;
        }
    }

    std::string GeoModelGeologicalEntity::geol_name(
        GeoModelGeologicalEntity::GEOL_FEATURE t )
    {
        switch( t ) {
            case GEOL_FEATURE::STRATI:
                return "top";
            case GEOL_FEATURE::FAULT:
                return "fault";
            case GEOL_FEATURE::REVERSE_FAULT:
                return "reverse_fault";
            case GEOL_FEATURE::NORMAL_FAULT:
                return "normal_fault";
            case GEOL_FEATURE::UNCONFORMITY:
                return "unconformity";
            case GEOL_FEATURE::VOI:
                return "boundary";
            case GEOL_FEATURE::NO_GEOL:
                return "no_geological_feature";
            default:
                return "no_geological_feature";
                break;
        }
    }

    const gmme_id& GeoModelGeologicalEntity::child_gmme( index_t x ) const
    {
        ringmesh_assert( x < nb_children() );
        return geomodel().entity_type_manager().relationship_manager.child_of_gmge(
            children_[x] );
    }
    const GeoModelMeshEntity& GeoModelGeologicalEntity::child( index_t x ) const
    {
        return geomodel().mesh_entity( child_gmme( x ) );
    }

    bool GeoModelGeologicalEntity::is_on_voi() const
    {
        for( index_t i = 0; i < nb_children(); i++ ) {
            if( !child( i ).is_on_voi() ) return false;
        }
        return true;
    }

    bool GeoModelGeologicalEntity::is_index_valid() const
    {
        return index() < geomodel().nb_geological_entities( type_name() );
    }

    bool GeoModelGeologicalEntity::is_connectivity_valid() const
    {
        bool valid = true;
        if( nb_children() == 0 ) {
            Logger::warn( "GeologicalEntity", gmge(), " is undefined. No child. " );
            valid = false;
        } else {
            // All children must have this entity as a parent
            const GeologicalEntityType entity_type = type_name();
            for( index_t i = 0; i < nb_children(); ++i ) {
                const GeoModelMeshEntity& one_child = child( i );
                if( one_child.parent_gmge( entity_type ) != gmge() ) {
                    Logger::warn( "GeoModelEntity",
                        "Inconsistency child-parent between ", gmge(), " and ",
                        one_child.gmme() );
                    valid = false;
                }
            }
        }
        return valid;
    }

    bool GeoModelGeologicalEntity::is_identification_valid() const
    {
        bool is_valid = true;
        if( !gmge().is_defined() ) {
            Logger::err( "GeoModelGeologicalEntity",
                " Entity associated to geomodel ", geomodel().name(),
                "has no type and/or no index " );
            is_valid = false;
            // No further checks are possible - This really should not happen
            ringmesh_assert_not_reached;
        }
        if( !is_index_valid() ) {
            Logger::warn( "GeoModelGeologicalEntity", " Entity index ", gmge(),
                " is not valid. " );
            // This really should not happen
            is_valid = false;
            ringmesh_assert_not_reached;
        }
        return is_valid;
    }

    bool GeoModelGeologicalEntity::is_valid() const
    {
        return check_has_children( *this );
    }

    void GeoModelGeologicalEntity::initialize()
    {
        ringmesh_register_GeoModelGeologicalEntity_creator( Contact );
        ringmesh_register_GeoModelGeologicalEntity_creator( Interface );
        ringmesh_register_GeoModelGeologicalEntity_creator( Layer );
    }

    MeshEntityType Contact::child_type_name() const
    {
        return Line::type_name_static();
    }

    MeshEntityType Interface::child_type_name() const
    {
        return Surface::type_name_static();
    }

    MeshEntityType Layer::child_type_name() const
    {
        return Region::type_name_static();
    }

    std::unique_ptr< GeoModelGeologicalEntity > GeoModelGeologicalEntityAccess::create_geological_entity(
        const GeologicalEntityType& type,
        const GeoModel& geomodel,
        index_t index_in_geomodel )
    {
        GeoModelGeologicalEntity* GMGE =
            GeoModelGeologicalEntityFactory::create_object( type, geomodel );
        GMGE->id_ = index_in_geomodel;
        return std::unique_ptr< GeoModelGeologicalEntity >( GMGE );
    }

}
