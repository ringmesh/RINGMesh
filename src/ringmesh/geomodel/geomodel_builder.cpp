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

#include <ringmesh/geomodel/geomodel_builder.h>

#include <stack>

#include <ringmesh/basic/geometry.h>

#include <ringmesh/geomodel/geomodel_api.h>

/*!
 * @file ringmesh/geomodel/geomodel_builder.cpp
 * @brief Implementation of the classes to build GeoModel from various inputs
 * @author Jeanne Pellerin
 */

namespace RINGMesh {

    GeoModelBuilderFile::GeoModelBuilderFile(
        GeoModel< 3 >& geomodel,
        const std::string& filename )
        : GeoModelBuilder( geomodel ), filename_( filename )
    {

    }

    template< index_t DIMENSION >
    GeoModelBuilderCopy< DIMENSION >::GeoModelBuilderCopy(
        GeoModelBuilder< DIMENSION >& builder,
        GeoModel< DIMENSION >& geomodel )
        : builder_( builder ), geomodel_( geomodel ), geomodel_access_( geomodel )
    {
    }

    template< index_t DIMENSION >
    void GeoModelBuilderCopy< DIMENSION >::copy_geomodel(
        const GeoModel< DIMENSION >& from )
    {
        builder_.topology.copy_topology( from );
        builder_.geometry.copy_meshes( from );
        builder_.geology.copy_geology( from );
    }

    template< index_t DIMENSION >
    GeoModelBuilderInfo< DIMENSION >::GeoModelBuilderInfo(
        GeoModelBuilder< DIMENSION >& builder,
        GeoModel< DIMENSION >& geomodel )
        : builder_( builder ), geomodel_( geomodel ), geomodel_access_( geomodel )
    {
    }

    template< index_t DIMENSION >
    GeoModelBuilder< DIMENSION >::GeoModelBuilder( GeoModel< DIMENSION >& geomodel )
        :
            topology( *this, geomodel ),
            geometry( *this, geomodel ),
            geology( *this, geomodel ),
            removal( *this, geomodel ),
            repair( *this, geomodel ),
            copy( *this, geomodel ),
            info( *this, geomodel ),
            from_surfaces( *this, geomodel ),
            geomodel_( geomodel ),
            geomodel_access_( geomodel )
    {
    }

    template< index_t DIMENSION >
    void GeoModelBuilder< DIMENSION >::end_geomodel()
    {
        if( geomodel_.name().empty() ) {
            info.set_geomodel_name( "model_default_name" );
        }

        geometry.cut_surfaces_by_internal_lines();
        geometry.cut_regions_by_internal_surfaces();
        topology.compute_universe();

        // Deliberate clear of the geomodel vertices used for geomodel building
        geomodel_.mesh.vertices.clear();
    }

    template< index_t DIMENSION >
    GeoModelBuilderGeology< DIMENSION >::GeoModelBuilderGeology(
        GeoModelBuilder< DIMENSION >& builder,
        GeoModel< DIMENSION >& geomodel )
        : builder_( builder ), geomodel_( geomodel ), geomodel_access_( geomodel )
    {
    }

    template< index_t DIMENSION >
    void GeoModelBuilderGeology< DIMENSION >::copy_geology(
        const GeoModel< DIMENSION >& from )
    {
        for( index_t t = 0; t < from.nb_geological_entity_types(); t++ ) {
            builder_.geology.copy_geological_entity_topology( from,
                from.geological_entity_type( t ) );
        }
    }

    template< index_t DIMENSION >
    bool GeoModelBuilderGeology< DIMENSION >::create_geological_entities(
        const GeologicalEntityType& type,
        index_t nb_additional_entities )
    {
        find_or_create_geological_entity_type( type );
        std::vector< std::unique_ptr< GeoModelGeologicalEntity< DIMENSION > > >& store =
            geomodel_access_.modifiable_geological_entities( type );
        index_t old_size = static_cast< index_t >( store.size() );
        index_t new_size = old_size + nb_additional_entities;
        store.reserve( new_size );
        for( index_t i = old_size; i < new_size; i++ ) {
            store.emplace_back(
                GeoModelGeologicalEntityAccess< DIMENSION >::create_geological_entity(
                    type, geomodel_, i ) );
        }
        return true;
    }

    template< index_t DIMENSION >
    bool GeoModelBuilderGeology< DIMENSION >::check_if_boundary_incident_entity_relation_already_exists(
        const gmge_id& parent,
        const gmme_id& children )
    {
        const GeoModelMeshEntity< DIMENSION >& children_mesh_entity =
            geomodel_.mesh_entity( children );
        if( children_mesh_entity.has_parent( parent.type() ) ) {
            ringmesh_assert(
                parent == children_mesh_entity.parent_gmge( parent.type() ) );
            return true;
        } else {
            return false;
        }
    }

    template< index_t DIMENSION >
    void GeoModelBuilderGeology< DIMENSION >::add_parent_children_relation(
        const gmge_id& parent,
        const gmme_id& children )
    {
        GeoModelGeologicalEntity< DIMENSION >& parent_entity =
            geomodel_access_.modifiable_geological_entity( parent );
        const std::vector< GeologicalEntityType >& parent_entity_types =
            geomodel_.entity_type_manager().relationship_manager.parent_types(
                children.type() );
        if( !contains( parent_entity_types, parent.type() ) ) {
            std::ostringstream message;
            message << "Wrong Parent type in the parent children relation between "
                << parent << " and " << children;
            throw RINGMeshException( "Entity", message.str() );
        }

        if( check_if_boundary_incident_entity_relation_already_exists( parent,
            children ) ) {
            return;
        }

        GeoModelMeshEntity< DIMENSION >& children_entity =
            geomodel_access_.modifiable_mesh_entity( children );
        const MeshEntityType& children_type =
            geomodel_.entity_type_manager().relationship_manager.child_type(
                parent.type() );

        if( children_type != children.type() ) {
            std::ostringstream message;
            message << "Wrong children type in the parent children relation between "
                << parent << " and " << children;
            throw RINGMeshException( "Entity", message.str() );
        }

        RelationshipManager& manager =
            geomodel_access_.modifiable_entity_type_manager().relationship_manager;
        index_t relation_id = manager.add_parent_child_relationship( parent,
            children );
        GeoModelGeologicalEntityAccess< DIMENSION > parent_access( parent_entity );
        parent_access.modifiable_children().push_back( relation_id );
        GeoModelMeshEntityAccess< DIMENSION > children_access( children_entity );
        children_access.modifiable_parents().push_back( relation_id );
    }

    template< index_t DIMENSION >
    void GeoModelBuilderGeology< DIMENSION >::remove_parent_children_relation(
        const gmge_id& parent,
        const gmme_id& children )
    {
        RelationshipManager& manager =
            geomodel_access_.modifiable_entity_type_manager().relationship_manager;
        index_t relation_id = manager.find_parent_child_relationship( parent,
            children );
        if( relation_id == NO_ID ) {
            std::ostringstream message;
            message << "No parent children relation found between " << parent
                << " and " << children;
            throw RINGMeshException( "Entity", message.str() );
        }
        GeoModelGeologicalEntityAccess< DIMENSION > parent_access(
            geomodel_access_.modifiable_geological_entity( parent ) );
        std::vector< index_t >& childs = parent_access.modifiable_children();
        std::remove_if( childs.begin(), childs.end(),
            [relation_id](index_t relation) {return relation == relation_id;} );
        GeoModelMeshEntityAccess< DIMENSION > children_access(
            geomodel_access_.modifiable_mesh_entity( children ) );
        std::vector< index_t >& parents = children_access.modifiable_parents();
        std::remove_if( parents.begin(), parents.end(),
            [relation_id](index_t relation) {return relation == relation_id;} );
    }

    template< index_t DIMENSION >
    void GeoModelBuilderGeology< DIMENSION >::delete_geological_entity(
        const GeologicalEntityType& type,
        index_t index )
    {
        geomodel_access_.modifiable_geological_entities( type )[index].reset();
    }

    template< index_t DIMENSION >
    gmge_id GeoModelBuilderGeology< DIMENSION >::create_geological_entity(
        const GeologicalEntityType& type )
    {
        index_t index = find_or_create_geological_entity_type( type );
        index_t id =
            static_cast< index_t >( geomodel_.nb_geological_entities( type ) );
        geomodel_access_.modifiable_geological_entities()[index].emplace_back(
            GeoModelGeologicalEntityAccess< DIMENSION >::create_geological_entity(
                type, geomodel_, id ) );
        return geomodel_access_.modifiable_geological_entities()[index].back()->gmge();
    }

    template< index_t DIMENSION >
    index_t GeoModelBuilderGeology< DIMENSION >::find_or_create_geological_entity_type(
        const GeologicalEntityType& type )
    {
        index_t type_index =
            geomodel_.entity_type_manager().geological_entity_manager.geological_entity_type_index(
                type );
        if( type_index == NO_ID ) {
            type_index = create_geological_entity_type( type );
        }
        return type_index;
    }

    template< index_t DIMENSION >
    index_t GeoModelBuilderGeology< DIMENSION >::create_geological_entity_type(
        const GeologicalEntityType& type )
    {
        ringmesh_assert( GeoModelGeologicalEntityFactory3D::has_creator( type ) );

        geomodel_access_.modifiable_entity_type_manager().geological_entity_manager.geological_entity_types_.push_back(
            type );
        geomodel_access_.modifiable_geological_entities().push_back(
            std::vector< std::unique_ptr< GeoModelGeologicalEntity< DIMENSION > > >() );
        std::unique_ptr< GeoModelGeologicalEntity< DIMENSION > > E(
            GeoModelGeologicalEntityFactory3D::create_object( type, geomodel_ ) );

        const MeshEntityType child_type = E->child_type_name();
        RelationshipManager& parentage =
            geomodel_access_.modifiable_entity_type_manager().relationship_manager;
        parentage.register_geology_relationship( type, child_type );

        return geomodel_.entity_type_manager().geological_entity_manager.nb_geological_entity_types()
            - 1;
    }

    template< index_t DIMENSION >
    void GeoModelBuilderGeology< DIMENSION >::copy_geological_entity_topology(
        const GeoModel< DIMENSION >& from,
        const GeologicalEntityType& type )
    {
        create_geological_entities( type, from.nb_geological_entities( type ) );

        RINGMESH_PARALLEL_LOOP
        for( index_t e = 0; e < geomodel_.nb_geological_entities( type ); ++e ) {
            gmge_id id( type, e );
            GeoModelGeologicalEntityAccess< DIMENSION > gmge_access(
                geomodel_access_.modifiable_geological_entity( id ) );
            gmge_access.copy( from.geological_entity( id ) );
        }
    }

    template< index_t DIMENSION >
    void GeoModelBuilderGeology< DIMENSION >::build_contacts()
    {
        if( geomodel_.entity_type_manager().geological_entity_manager.is_valid_type(
            Contact< DIMENSION >::type_name_static() )
            && geomodel_.nb_geological_entities(
                Contact< DIMENSION >::type_name_static() ) > 0 ) {
            return;
        }

        std::vector< std::set< gmge_id > > interfaces;
        for( index_t i = 0; i < geomodel_.nb_lines(); ++i ) {
            const Line< DIMENSION >& L = geomodel_.line( i );
            std::set< gmge_id > cur_interfaces;
            for( index_t j = 0; j < L.nb_incident_entities(); ++j ) {
                const GeoModelMeshEntity< DIMENSION >& S = L.incident_entity( j );
                gmge_id parent_interface = S.parent_gmge(
                    Interface< DIMENSION >::type_name_static() );
                cur_interfaces.insert( parent_interface );
            }
            gmge_id contact_id;
            for( index_t j = 0; j < interfaces.size(); ++j ) {
                if( cur_interfaces.size() == interfaces[j].size()
                    && std::equal( cur_interfaces.begin(), cur_interfaces.end(),
                        interfaces[j].begin() ) ) {
                    contact_id = gmge_id( Contact< DIMENSION >::type_name_static(),
                        j );
                    break;
                }
            }
            if( !contact_id.is_defined() ) {
                contact_id = create_geological_entity(
                    Contact< DIMENSION >::type_name_static() );
                ringmesh_assert( contact_id.index() == interfaces.size() );
                interfaces.push_back( cur_interfaces );
                // Create a name for this contact
                std::string name = "contact";
                for( const gmge_id& it : cur_interfaces ) {
                    name += "_";
                    name += geomodel_.geological_entity( it ).name();
                }
                builder_.info.set_geological_entity_name( contact_id, name );
            }
            add_parent_children_relation( contact_id, L.gmme() );
        }
    }

    template< index_t DIMENSION >
    void GeoModelBuilderGeology< DIMENSION >::set_geological_entity_geol_feature(
        const gmge_id& gmge_id,
        typename GeoModelGeologicalEntity< DIMENSION >::GEOL_FEATURE geol_feature )
    {
        GeoModelGeologicalEntityAccess< DIMENSION > gmge_access(
            geomodel_access_.modifiable_geological_entity( gmge_id ) );
        gmge_access.modifiable_geol_feature() = geol_feature;
    }
    //    template class RINGMESH_API GeoModelBuilderInfo< 2 > ;
    //    template class RINGMESH_API GeoModelBuilderGeology< 2 > ;
    //    template class RINGMESH_API GeoModelBuilderCopy< 2 > ;
    //    template class RINGMESH_API GeoModelBuilderFromSurfaces< 2 > ;

    template class RINGMESH_API GeoModelBuilder< 3 > ;
    template class RINGMESH_API GeoModelBuilderInfo< 3 > ;
    template class RINGMESH_API GeoModelBuilderGeology< 3 > ;
    template class RINGMESH_API GeoModelBuilderCopy< 3 > ;
//    template class RINGMESH_API GeoModelBuilderFromSurfaces< 3 > ;

} // namespace
