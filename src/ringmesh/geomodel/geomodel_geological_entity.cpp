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

    template< index_t DIMENSION >
    bool check_has_children( const GeoModelGeologicalEntity< DIMENSION >& E )
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

    template< index_t DIMENSION >
    typename GeoModelGeologicalEntity< DIMENSION >::GEOL_FEATURE GeoModelGeologicalEntity<
        DIMENSION >::determine_geological_type( const std::string& in )
    {
        if( in == "reverse_fault" ) {
            return REVERSE_FAULT;
        } else if( in == "normal_fault" ) {
            return NORMAL_FAULT;
        } else if( in == "fault" ) {
            return FAULT;
        } else if( in == "top" ) {
            return STRATI;
        } else if( in == "none" ) {
            // This might seem strange - but it seems that what's
            // Gocad is doing
            return STRATI;
        } else if( in == "topographic" ) {
            return STRATI;
        } else if( in == "unconformity" ) {
            return UNCONFORMITY;
        } else if( in == "boundary" ) {
            return VOI;
        } else {
            // Default case - no information
            return NO_GEOL;
        }
    }

    template< index_t DIMENSION >
    std::string GeoModelGeologicalEntity< DIMENSION >::geol_name(
        GeoModelGeologicalEntity< DIMENSION >::GEOL_FEATURE t )
    {
        switch( t ) {
            case STRATI:
                return "top";
            case FAULT:
                return "fault";
            case REVERSE_FAULT:
                return "reverse_fault";
            case NORMAL_FAULT:
                return "normal_fault";
            case UNCONFORMITY:
                return "unconformity";
            case VOI:
                return "boundary";
            case NO_GEOL:
                return "no_geological_feature";
            default:
                return "no_geological_feature";
                break;
        }
    }

    template< index_t DIMENSION >
    const gmme_id& GeoModelGeologicalEntity< DIMENSION >::child_gmme(
        index_t x ) const
    {
        ringmesh_assert( x < nb_children() );
        return this->geomodel().entity_type_manager().relationship_manager.child_of_gmge(
            children_[x] );
    }

    template< index_t DIMENSION >
    const GeoModelMeshEntity< DIMENSION >& GeoModelGeologicalEntity< DIMENSION >::child(
        index_t x ) const
    {
        return this->geomodel().mesh_entity( child_gmme( x ) );
    }

    template< index_t DIMENSION >
    bool GeoModelGeologicalEntity< DIMENSION >::is_on_voi() const
    {
        for( index_t i = 0; i < nb_children(); i++ ) {
            if( !child( i ).is_on_voi() ) return false;
        }
        return true;
    }

    template< index_t DIMENSION >
    bool GeoModelGeologicalEntity< DIMENSION >::is_index_valid() const
    {
        return this->index() < this->geomodel().nb_geological_entities( type_name() );
    }

    template< index_t DIMENSION >
    bool GeoModelGeologicalEntity< DIMENSION >::is_connectivity_valid() const
    {
        bool valid = true;
        if( nb_children() == 0 ) {
            Logger::warn( "GeologicalEntity", gmge(), " is undefined. No child. " );
            valid = false;
        } else {
            // All children must have this entity as a parent
            const GeologicalEntityType entity_type = type_name();
            for( index_t i = 0; i < nb_children(); ++i ) {
                const GeoModelMeshEntity< DIMENSION >& one_child = child( i );
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

    template< index_t DIMENSION >
    bool GeoModelGeologicalEntity< DIMENSION >::is_identification_valid() const
    {
        bool is_valid = true;
        if( !gmge().is_defined() ) {
            Logger::err( "GeoModelGeologicalEntity",
                " Entity associated to geomodel ", this->geomodel().name(),
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

    template< index_t DIMENSION >
    bool GeoModelGeologicalEntity< DIMENSION >::is_valid() const
    {
        return check_has_children( *this );
    }

    template< index_t DIMENSION >
    void GeoModelGeologicalEntity< DIMENSION >::initialize()
    {
        ringmesh_register_GeoModelGeologicalEntity_creator( Contact< DIMENSION > );
        ringmesh_register_GeoModelGeologicalEntity_creator( Interface< DIMENSION > );
        ringmesh_register_GeoModelGeologicalEntity_creator( Layer< DIMENSION > );
    }

    template< index_t DIMENSION >
    MeshEntityType Contact< DIMENSION >::child_type_name() const
    {
        return Line< DIMENSION >::type_name_static();
    }

    template< index_t DIMENSION >
    MeshEntityType Interface< DIMENSION >::child_type_name() const
    {
        return Surface< DIMENSION >::type_name_static();
    }

    template< index_t DIMENSION >
    MeshEntityType Layer< DIMENSION >::child_type_name() const
    {
        return Region< DIMENSION >::type_name_static();
    }

    template< index_t DIMENSION >
    std::unique_ptr< GeoModelGeologicalEntity< DIMENSION > > GeoModelGeologicalEntityAccess<
        DIMENSION >::create_geological_entity(
        const GeologicalEntityType& type,
        const GeoModel& geomodel,
        index_t index_in_geomodel )
    {
        GeoModelGeologicalEntity< DIMENSION >* GMGE =
            GeoModelGeologicalEntityFactory::create_object( type, geomodel );
        GMGE->id_ = index_in_geomodel;
        return std::unique_ptr< GeoModelGeologicalEntity< DIMENSION > >( GMGE );
    }

//    template class GeoModelGeologicalEntity< 2 > ;
//    template class GeoModelGeologicalEntityAccess< 2 > ;

    template class GeoModelGeologicalEntity< 3 > ;
    template class GeoModelGeologicalEntityAccess< 3 > ;

}
