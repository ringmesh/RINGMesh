/*
 * Copyright (c) 2012-2018, Association Scientifique pour la Geologie et ses
 * Applications (ASGA). All rights reserved.
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
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
 * THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 * PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL ASGA BE LIABLE FOR ANY DIRECT,
 * INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF
 * THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
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

#include <ringmesh/geomodel/core/entity_type.h>
#include <ringmesh/geomodel/core/geomodel.h>
#include <ringmesh/geomodel/core/geomodel_geological_entity.h>
#include <ringmesh/geomodel/core/geomodel_mesh_entity.h>

namespace
{
    using namespace RINGMesh;

    template < index_t DIMENSION >
    bool check_has_children( const GeoModelGeologicalEntity< DIMENSION >& E )
    {
        if( E.nb_children() == 0 )
        {
            Logger::warn(
                "GeoModel", E.type_name(), " ", E.index(), " has no children" );
            return false;
        }
        return true;
    }
} // namespace

namespace RINGMesh
{
    template < index_t DIMENSION >
    typename GeoModelGeologicalEntity< DIMENSION >::GEOL_FEATURE
        GeoModelGeologicalEntity< DIMENSION >::determine_geological_type(
            const std::string& in )
    {
        if( in == "reverse_fault" )
        {
            return GEOL_FEATURE::REVERSE_FAULT;
        }
        if( in == "normal_fault" )
        {
            return GEOL_FEATURE::NORMAL_FAULT;
        }
        if( in == "fault" )
        {
            return GEOL_FEATURE::FAULT;
        }
        if( in == "top" )
        {
            return GEOL_FEATURE::STRATI;
        }
        if( in == "none" )
        {
            // This might seem strange - but it seems that what's
            // Gocad is doing
            return GEOL_FEATURE::STRATI;
        }
        if( in == "topographic" )
        {
            return GEOL_FEATURE::STRATI;
        }
        if( in == "unconformity" )
        {
            return GEOL_FEATURE::UNCONFORMITY;
        }
        if( in == "boundary" )
        {
            return GEOL_FEATURE::VOI;
        }
        if( in == "lease" )
        {
            return GEOL_FEATURE::VOI;
        }
        // Default case - no information
        return GEOL_FEATURE::NO_GEOL;
    }

    template < index_t DIMENSION >
    std::string GeoModelGeologicalEntity< DIMENSION >::geol_name(
        typename GeoModelGeologicalEntity< DIMENSION >::GEOL_FEATURE feature )
    {
        switch( feature )
        {
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

    template < index_t DIMENSION >
    bool GeoModelGeologicalEntity< DIMENSION >::is_fault(
        typename GeoModelGeologicalEntity< DIMENSION >::GEOL_FEATURE feature )
    {
        return feature == GEOL_FEATURE::FAULT
               || feature == GEOL_FEATURE::REVERSE_FAULT
               || feature == GEOL_FEATURE::NORMAL_FAULT;
    }

    template < index_t DIMENSION >
    const gmme_id& GeoModelGeologicalEntity< DIMENSION >::child_gmme(
        index_t x ) const
    {
        ringmesh_assert( x < nb_children() );
        return this->geomodel()
            .entity_type_manager()
            .relationship_manager.child_of_gmge( children_[x] );
    }

    template < index_t DIMENSION >
    gmge_id GeoModelGeologicalEntity< DIMENSION >::gmge() const
    {
        return gmge_id( type_name(), this->index() );
    }

    template < index_t DIMENSION >
    GeologicalEntityType
        GeoModelGeologicalEntity< DIMENSION >::entity_type() const
    {
        return gmge().type();
    }

    template < index_t DIMENSION >
    GeologicalEntityType
        GeoModelGeologicalEntity< DIMENSION >::type_name_static()
    {
        return ForbiddenGeologicalEntityType::type_name_static();
    }

    template < index_t DIMENSION >
    GeologicalEntityType
        GeoModelGeologicalEntity< DIMENSION >::type_name() const
    {
        return type_name_static();
    }

    template < index_t DIMENSION >
    const GeoModelMeshEntity< DIMENSION >&
        GeoModelGeologicalEntity< DIMENSION >::child( index_t x ) const
    {
        return this->geomodel().mesh_entity( child_gmme( x ) );
    }

    template < index_t DIMENSION >
    bool GeoModelGeologicalEntity< DIMENSION >::is_on_voi() const
    {
        for( auto i : range( nb_children() ) )
        {
            if( !child( i ).is_on_voi() )
            {
                return false;
            }
        }
        return true;
    }

    template < index_t DIMENSION >
    bool GeoModelGeologicalEntity< DIMENSION >::is_index_valid() const
    {
        return this->index()
               < this->geomodel().nb_geological_entities( type_name() );
    }

    template < index_t DIMENSION >
    void GeoModelGeologicalEntity< DIMENSION >::copy_geological_entity(
        const GeoModelGeologicalEntity< DIMENSION >& from )
    {
        this->copy_name( from );
        geol_feature_ = from.geol_feature_;
        children_ = from.children_;
    }

    template < index_t DIMENSION >
    bool GeoModelGeologicalEntity< DIMENSION >::is_connectivity_valid() const
    {
        bool valid = true;
        if( nb_children() == 0 )
        {
            Logger::warn(
                "GeologicalEntity", gmge(), " is undefined. No child. " );
            valid = false;
        }
        else
        {
            // All children must have this entity as a parent
            const GeologicalEntityType entity_type = type_name();
            for( auto i : range( nb_children() ) )
            {
                const GeoModelMeshEntity< DIMENSION >& one_child = child( i );
                if( one_child.parent_gmge( entity_type ) != gmge() )
                {
                    Logger::warn( "GeoModelEntity",
                        "Inconsistency child-parent between ", gmge(), " and ",
                        one_child.gmme() );
                    valid = false;
                }
            }
        }
        return valid;
    }

    template < index_t DIMENSION >
    bool GeoModelGeologicalEntity< DIMENSION >::is_identification_valid() const
    {
        bool is_valid = true;
        if( !gmge().is_defined() )
        {
            Logger::err( "GeoModelGeologicalEntity",
                " Entity associated to geomodel ", this->geomodel().name(),
                "has no type and/or no index " );
            is_valid = false;
            // No further checks are possible - This really should not happen
            ringmesh_assert_not_reached;
        }
        if( !is_index_valid() )
        {
            Logger::warn( "GeoModelGeologicalEntity", " Entity index ", gmge(),
                " is not valid. " );
            // This really should not happen
            is_valid = false;
            ringmesh_assert_not_reached;
        }
        return is_valid;
    }

    template < index_t DIMENSION >
    bool GeoModelGeologicalEntity< DIMENSION >::is_valid() const
    {
        return check_has_children( *this );
    }

    template <>
    void GeoModelGeologicalEntity< 2 >::initialize()
    {
        GeoModelGeologicalEntityFactory2D::register_creator< Contact2D >(
            Contact2D::type_name_static() );
        GeoModelGeologicalEntityFactory2D::register_creator< Interface2D >(
            Interface2D::type_name_static() );
        GeoModelGeologicalEntityFactory2D::register_creator< Layer2D >(
            Layer2D::type_name_static() );
    }

    template <>
    void GeoModelGeologicalEntity< 3 >::initialize()
    {
        GeoModelGeologicalEntityFactory3D::register_creator< Contact3D >(
            Contact3D::type_name_static() );
        GeoModelGeologicalEntityFactory3D::register_creator< Interface3D >(
            Interface3D::type_name_static() );
        GeoModelGeologicalEntityFactory3D::register_creator< Layer3D >(
            Layer3D::type_name_static() );
    }

    template <>
    MeshEntityType Contact< 2 >::child_type_name() const
    {
        return Corner2D::type_name_static();
    }

    template <>
    MeshEntityType Interface< 2 >::child_type_name() const
    {
        return Line2D::type_name_static();
    }

    template <>
    MeshEntityType Layer< 2 >::child_type_name() const
    {
        return Surface2D::type_name_static();
    }

    template < index_t DIMENSION >
    GeologicalEntityType Contact< DIMENSION >::type_name_static()
    {
        return GeologicalEntityType( "Contact" );
    }

    template < index_t DIMENSION >
    GeologicalEntityType Contact< DIMENSION >::type_name() const
    {
        return type_name_static();
    }

    template <>
    MeshEntityType Contact< 3 >::child_type_name() const
    {
        return Line3D::type_name_static();
    }

    template < index_t DIMENSION >
    GeologicalEntityType Interface< DIMENSION >::type_name_static()
    {
        return GeologicalEntityType( "Interface" );
    }

    template < index_t DIMENSION >
    GeologicalEntityType Interface< DIMENSION >::type_name() const
    {
        return type_name_static();
    }

    template <>
    MeshEntityType Interface< 3 >::child_type_name() const
    {
        return Surface3D::type_name_static();
    }

    template < index_t DIMENSION >
    GeologicalEntityType Layer< DIMENSION >::type_name_static()
    {
        return GeologicalEntityType( "Layer" );
    }

    template < index_t DIMENSION >
    GeologicalEntityType Layer< DIMENSION >::type_name() const
    {
        return type_name_static();
    }

    template <>
    MeshEntityType Layer< 3 >::child_type_name() const
    {
        return Region3D::type_name_static();
    }

    template class geomodel_core_api GeoModelGeologicalEntity< 2 >;
    template class geomodel_core_api Contact< 2 >;
    template class geomodel_core_api Interface< 2 >;
    template class geomodel_core_api Layer< 2 >;

    template class geomodel_core_api GeoModelGeologicalEntity< 3 >;
    template class geomodel_core_api Contact< 3 >;
    template class geomodel_core_api Interface< 3 >;
    template class geomodel_core_api Layer< 3 >;

} // namespace RINGMesh
