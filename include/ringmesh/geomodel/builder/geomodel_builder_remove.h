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

#pragma once

#include <ringmesh/geomodel/builder/common.h>
#include <ringmesh/geomodel/builder/geomodel_builder_access.h>
#include <set>

/*!
 * @brief Builder tools to remove entities from a GeoModel.
 * @author Antoine Mazuyer
 * @author Pierre Anquez
 */

namespace RINGMesh
{
    FORWARD_DECLARATION_DIMENSION_CLASS( GeoModelBuilderBase );
    FORWARD_DECLARATION_DIMENSION_CLASS( GeoModelBuilder );
    FORWARD_DECLARATION_DIMENSION_CLASS( GeoModelMeshEntity );
    FORWARD_DECLARATION_DIMENSION_CLASS( Region );

    ALIAS_3D( GeoModel );
    ALIAS_3D( GeoModelBuilder );
    ALIAS_3D( GeoModelMeshEntity );
    ALIAS_3D( Region );
} // namespace RINGMesh

namespace RINGMesh
{
    /*!
     * @brief Builder tools to remove entities from a GeoModel
     */
    template < index_t DIMENSION >
    class geomodel_builder_api GeoModelBuilderRemoveBase
    {
        ringmesh_disable_copy_and_move( GeoModelBuilderRemoveBase );
        ringmesh_template_assert_2d_or_3d( DIMENSION );

    public:
        virtual ~GeoModelBuilderRemoveBase();

        /*!
         * @brief Remove a list of mesh entities of the geomodel
         * @details No check is done on the consistency of this removal
         *          The entities and all references to them are removed.
         *          All dependent entities should be in the set of entities to
         * remove,
         *          with a prior call to get_dependent_entities function.
         *
         */
        void remove_mesh_entities( const std::set< gmme_id >& entities );

        /*!
         * @brief Remove a list of geological entities of the geomodel
         * @details No check is done on the consistency of this removal
         *          The entities and all references to them are removed.
         *          All dependent entities should be in the set of entities to
         * remove,
         *          with a prior call to get_dependent_entities function.
         *
         */
        void remove_geological_entities( const std::set< gmge_id >& entities );

        /*!
         * Should be rewritten. Put as it was before someone removed it...
         */
        void remove_entities_and_dependencies(
            const std::set< gmme_id >& entities_to_remove );

    protected:
        GeoModelBuilderRemoveBase( GeoModelBuilder< DIMENSION >& builder,
            GeoModel< DIMENSION >& geomodel );

        virtual void update_mesh_entity( GeoModelMeshEntity< DIMENSION >& ME );

        // ---  High level functions ----------
        void initialize_for_removal(
            const std::set< gmme_id >& mesh_entities_to_remove );
        void do_delete_flagged_mesh_entities();
        void do_delete_flagged_geological_entities();

        void check_if_entities_are_meshed(
            const std::set< gmme_id >& mesh_entities_to_remove );

        void initialize_costly_storage();

        void delete_mesh_entity( index_t type, index_t index );

        void clear_null_mesh_entities( index_t type );

        void clear_null_geological_entities( index_t type );

        void update_mesh_entity_connectivity();

        void update_geological_entity_connectivity();

        //        void remove_dependencies()
        //        {
        //            std::set< gme_id > new_gmme_to_remove ;
        //            for( index_t me = 0;
        //                me < geomodel().nb_mesh_entities( starting_dependency_
        //                ); me++ ) {
        //                const GeoModelMeshEntity& cur_gmme =
        //                geomodel().mesh_entity(
        //                    starting_dependency_, me ) ;
        //                if( cur_gmme.incident_entity( 0 ).index() == NO_ID
        //                    && cur_gmme.nb_incident_entities() == 1 ) {
        //                    new_gmme_to_remove.insert( cur_gmme.gme_id() ) ;
        //                }
        //            }
        //            if( starting_dependency_ != Corner::type_name_static() ) {
        //                starting_dependency_ =
        //                EntityTypeManager::boundary_type(
        //                    starting_dependency_ ) ;
        //                remove_mesh_entities_with_dependencies(
        //                new_gmme_to_remove ) ;
        //            }
        //        }

        //------  Initialization -------
        void fill_removed_entities_and_mapping();

        void fill_to_erase_vectors(
            const std::set< gmme_id >& mesh_entities_to_remove );

        void fill_nb_children_vector();

        void fill_nb_initial_entities();

        void fill_entity_type_to_index_map();

        // ---- Easier access to relationships between EntityTypes
        index_t mesh_entity_type_index(
            const GeoModelMeshEntity< DIMENSION >& E ) const;

        index_t geological_entity_type_index(
            const GeoModelGeologicalEntity< DIMENSION >& E ) const;

        index_t children_type_index( const GeologicalEntityType& type ) const;

        const MeshEntityType children_type(
            const GeologicalEntityType& type ) const;

        index_t boundary_type_index( const MeshEntityType& type ) const;

        const MeshEntityType& boundary_entity_type(
            const MeshEntityType& type ) const;

        /// TODO unused function. To handle during removal refactoring BC.
        index_t incident_entity_type_to_index(
            const MeshEntityType& type ) const;

        const MeshEntityType& incident_entity_type(
            const MeshEntityType& type ) const;

        bool is_mesh_entity( index_t i ) const
        {
            return i < nb_mesh_entity_types_;
        }

        bool is_geological_entity( index_t i ) const
        {
            return !is_mesh_entity( i );
        }

        // ----  Update connectivity functions  ------

        void flag_geological_entities_without_children();

        void set_mesh_entity_index(
            GeoModelMeshEntity< DIMENSION >& mesh_entity,
            index_t new_index_in_geomodel );
        void set_geological_entity_index(
            GeoModelGeologicalEntity< DIMENSION >& geological_entity,
            index_t new_index_in_geomodel );

        void update_mesh_entity_index(
            GeoModelMeshEntity< DIMENSION >& mesh_entity );
        void update_geological_entity_index(
            GeoModelGeologicalEntity< DIMENSION >& geological_entity );
        void update_mesh_entity_boundaries(
            GeoModelMeshEntity< DIMENSION >& mesh_entity );

        void set_boundary_side(
            Region3D& region, index_t boundary_index, bool new_side );

        void update_mesh_entity_incident_entity(
            GeoModelMeshEntity< DIMENSION >& mesh_entity );
        void update_mesh_entity_parents(
            GeoModelMeshEntity< DIMENSION >& mesh_entity );
        void update_geological_entity_children(
            GeoModelGeologicalEntity< DIMENSION >& geological_entity );

        // --- Deletion of some values the GeoModel storage
        template < typename TEST, typename THINGS_TO_DELETE >
        void remove_invalid_values(
            std::vector< THINGS_TO_DELETE >& vector, const TEST& test );

        void delete_invalid_children(
            GeoModelGeologicalEntity< DIMENSION >& E );

        void delete_invalid_boundaries( GeoModelMeshEntity< DIMENSION >& E );

        void delete_invalid_incident_entity(
            GeoModelMeshEntity< DIMENSION >& E );

        void delete_invalid_parents( GeoModelMeshEntity< DIMENSION >& E );

        index_t mesh_entity_type_to_index( const MeshEntityType& type ) const;

        index_t geological_entity_type_to_index(
            const GeologicalEntityType& type ) const;
        const MeshEntityType& index_to_mesh_entity_type( index_t index ) const;

        const GeologicalEntityType& index_to_geological_entity_type(
            index_t index ) const;

        void reinitialize_utility_vectors();

    protected:
        GeoModelBuilder< DIMENSION >& builder_;
        GeoModel< DIMENSION >& geomodel_;
        GeoModelAccess< DIMENSION > geomodel_access_;

        index_t nb_entity_types_;
        index_t nb_geological_entity_types_;
        index_t nb_mesh_entity_types_;

        std::vector< index_t > nb_initial_mesh_entities_;
        std::vector< index_t > nb_initial_geological_entities_;

        std::vector< index_t > nb_removed_mesh_entities_;
        std::vector< index_t > nb_removed_geological_entities_;

        /*! For each type of entity, store a vector of where the
         * entities to remove from the geomodel are flagged with NO_ID. */
        std::vector< std::vector< bool > > mesh_entity_to_erase_;
        /*! Stores the mapping table between indices for each type of
         *  element before and after the removal of entities */
        std::vector< std::vector< index_t > > old_2_new_mesh_entity_;
        std::vector< std::vector< index_t > > nb_childs_;

        std::vector< std::vector< index_t > > old_2_new_geological_entity_;

        // std::map< EntityType, index_t > entity_type_to_index_ ;
        // std::map< index_t, EntityType > index_to_entity_type_ ;
        //        std::vector< EntityType > all_entity_types_ ;

        std::vector< MeshEntityType > mesh_entity_types_;
        std::vector< GeologicalEntityType > geological_entity_types_;
    };

    ALIAS_2D_AND_3D( GeoModelBuilderRemoveBase );

    template < index_t DIMENSION >
    class geomodel_builder_api GeoModelBuilderRemove
        : public GeoModelBuilderRemoveBase< DIMENSION >
    {
        friend class GeoModelBuilderBase< DIMENSION >;
        friend class GeoModelBuilder< DIMENSION >;

    public:
        GeoModelBuilderRemove( GeoModelBuilder< DIMENSION >& builder,
            GeoModel< DIMENSION >& geomodel );
    };

    template <>
    class geomodel_builder_api GeoModelBuilderRemove< 3 >
        : public GeoModelBuilderRemoveBase< 3 >
    {
        friend class GeoModelBuilderBase< 3 >;
        friend class GeoModelBuilder< 3 >;

    public:
        GeoModelBuilderRemove(
            GeoModelBuilder3D& builder, GeoModel3D& geomodel );

    private:
        void update_mesh_entity( GeoModelMeshEntity3D& ME ) override;

        void set_boundary_side(
            Region3D& R, index_t boundary_index, bool new_side );

        void update_region_boundary_signs( Region3D& R );

        void delete_invalid_signs( Region3D& R );
    };
} // namespace RINGMesh
