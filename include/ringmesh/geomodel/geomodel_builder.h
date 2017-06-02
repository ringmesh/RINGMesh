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

#pragma once

#include <ringmesh/basic/common.h>

#include <ringmesh/geomodel/geomodel.h>
#include <ringmesh/geomodel/geomodel_entity.h>
#include <ringmesh/geomodel/geomodel_mesh_entity.h>
#include <ringmesh/geomodel/geomodel_geological_entity.h>

#include <ringmesh/geomodel/geomodel_builder_geometry.h>
#include <ringmesh/geomodel/geomodel_builder_remove.h>
#include <ringmesh/geomodel/geomodel_builder_repair.h>
#include <ringmesh/geomodel/geomodel_builder_topology.h>

/*!
 * @brief Classes to build GeoModel from various inputs
 * @author Jeanne Pellerin
 */

namespace RINGMesh {
    /*!
     * @brief First draft of flags to build a GeoModel
     * @todo Implements functions to set, access the values, depending on what ?
     * To check the consistency of the options. What do we do about the other entities ? [JP] 
     * 
     * @todo We need to keep track of the status of the GeoModel when building it:
     * same flags or some others ?    
     *
     * @todo To separate in two classes ? One providing the low level functions set, assign etc,
     * and the other one some high level functions. [JP]
     */
    class GeoModelBuildingFlags {
    public:
        GeoModelBuildingFlags()
        {
            compute_corners = false;
            compute_lines = false;
            compute_surfaces = false;
            compute_regions_brep = false;
            compute_regions_mesh = false;
        }
        bool compute_corners;
        bool compute_lines;
        bool compute_surfaces;
        bool compute_regions_brep;
        bool compute_regions_mesh;
    };

    // Implementation details
    class GeoModelRegionFromSurfaces;
}

namespace RINGMesh {

    class RINGMESH_API GeoModelBuilderInfo {
    ringmesh_disable_copy( GeoModelBuilderInfo );
        friend class GeoModelBuilder;

    public:
        /*!
         *@brief Set the name of the geomodel
         */
        void set_geomodel_name( const std::string& name )
        {
            geomodel_access_.modifiable_name() = name;
        }

        /*!
         *@brief Set the name of a geomodel mesh entity
         */
        void set_mesh_entity_name( const gmme_id& gmme_id, const std::string& name )
        {
            GeoModelMeshEntityAccess< 3 > gmme_access(
                geomodel_access_.modifiable_mesh_entity( gmme_id ) );
            gmme_access.modifiable_name() = name;

        }

        /*!
         *@brief Set the name of a geomodel geological entity
         */
        void set_geological_entity_name(
            const gmge_id& gmge_id,
            const std::string& name )
        {
            GeoModelGeologicalEntityAccess< 3 > gmge_access(
                geomodel_access_.modifiable_geological_entity( gmge_id ) );
            gmge_access.modifiable_name() = name;
        }

    protected:
        GeoModelBuilderInfo( GeoModelBuilder& builder, GeoModel< 3 >& geomodel );

    private:
        GeoModelBuilder& builder_;
        GeoModel< 3 >& geomodel_;
        GeoModelAccess< 3 > geomodel_access_;

    };

    template< index_t DIMENSION >
    class GeoModelBuilderGeology {
    ringmesh_disable_copy( GeoModelBuilderGeology );
        ringmesh_template_assert_2d_or_3d( DIMENSION );
        friend class GeoModelBuilder;

    public:
        void copy_geology( const GeoModel< DIMENSION >& from );

        /*!
         * @brief Create and store a geological entity of the given type
         * @return The index of the created geological entity
         */
        gmge_id create_geological_entity( const GeologicalEntityType& type );

        bool create_geological_entities(
            const GeologicalEntityType& type,
            index_t nb );

        void set_geological_entity_geol_feature(
            const gmge_id& gmge_id,
            typename GeoModelGeologicalEntity< DIMENSION >::GEOL_FEATURE geol_feature );

        void set_mesh_entity_parent(
            const gmme_id& child_gmme,
            index_t id,
            const gmge_id& parent_gmge )
        {
            /// No check on the validity of the index of the entity parents_
            /// NO_ID is used to flag entities to delete
            GeoModelMeshEntity< DIMENSION >& mesh_entity =
                geomodel_access_.modifiable_mesh_entity( child_gmme );
            ringmesh_assert( id < mesh_entity.nb_parents() );
            GeoModelMeshEntityAccess< DIMENSION > gmme_access( mesh_entity );
            index_t relationship_id = gmme_access.modifiable_parents()[id];
            RelationshipManager& manager =
                geomodel_access_.modifiable_entity_type_manager().relationship_manager;
            manager.set_parent_to_parent_child_relationship( relationship_id,
                parent_gmge );
        }

        void add_parent_children_relation(
            const gmge_id& parent,
            const gmme_id& children );

        void remove_parent_children_relation(
            const gmge_id& parent,
            const gmme_id& children );
        void set_geological_entity_child(
            const gmge_id& parent_gmge,
            index_t id,
            index_t child_id )
        {
            /// No check on the validity of the index of the entity child_index
            /// NO_ID is used to flag entities to delete
            GeoModelGeologicalEntity< DIMENSION >& geol_entity =
                geomodel_access_.modifiable_geological_entity( parent_gmge );
            const MeshEntityType& child_type =
                geomodel_.entity_type_manager().relationship_manager.child_type(
                    parent_gmge.type() );
            gmme_id child( child_type, child_id );
            GeoModelGeologicalEntityAccess< DIMENSION > gmge_access( geol_entity );
            index_t relationship_id = gmge_access.modifiable_children()[id];
            RelationshipManager& manager =
                geomodel_access_.modifiable_entity_type_manager().relationship_manager;
            manager.set_child_to_parent_child_relationship( relationship_id, child );
        }

        void delete_geological_entity(
            const GeologicalEntityType& type,
            index_t index );

        /*!
         * @brief Build the Contacts
         * @details One contact is a group of lines shared by the same Interfaces
         */
        void build_contacts();

    protected:
        GeoModelBuilderGeology(
            GeoModelBuilder& builder,
            GeoModel< DIMENSION >& geomodel );

    private:
        index_t create_geological_entity_type( const GeologicalEntityType& type );

        index_t find_or_create_geological_entity_type(
            const GeologicalEntityType& type );

        void copy_geological_entity_topology(
            const GeoModel< DIMENSION >& from,
            const GeologicalEntityType& type );

        bool check_if_boundary_incident_entity_relation_already_exists(
            const gmge_id& parent,
            const gmme_id& children );

    private:
        GeoModelBuilder& builder_;
        GeoModel< DIMENSION >& geomodel_;
        GeoModelAccess< DIMENSION > geomodel_access_;
    };

    class RINGMESH_API GeoModelBuilderCopy {
    ringmesh_disable_copy( GeoModelBuilderCopy );
        friend class GeoModelBuilder;
    public:
        void copy_geomodel( const GeoModel< 3 >& from );

    private:
        GeoModelBuilderCopy( GeoModelBuilder& builder, GeoModel< 3 >& geomodel );

    private:
        GeoModelBuilder& builder_;
        GeoModel< 3 >& geomodel_;
        GeoModelAccess< 3 > geomodel_access_;
    };

    class RINGMESH_API GeoModelBuilderFromSurfaces {
    ringmesh_disable_copy( GeoModelBuilderFromSurfaces );
        friend class GeoModelBuilder;

    public:
        /*!
         * Delete all GeoModelRegionFromSurfaces owned by the builder
         */
        virtual ~GeoModelBuilderFromSurfaces();
        /*
         * @brief From a GeoModel in which only Surfaces are defined,
         * create Corners, Lines and Regions depending on the building flags
         * @note Validity is not checked
         */
        void build();

        /*!
         * @brief From the Surfaces of the GeoModel, build its Lines and Corners
         */
        bool build_lines_and_corners_from_surfaces();

        /*!
         * @brief Build the regions of the GeoModel from the Surfaces
         * @pre Function build_lines_and_corners_from_surfaces
         * must have been called before
         */
        bool build_brep_regions_from_surfaces();

    private:
        GeoModelBuilderFromSurfaces(
            GeoModelBuilder& builder,
            GeoModel< 3 >& geomodel );

    public:
        /*! Options to toggle the building of entities from the available entities */
        GeoModelBuildingFlags options_;

    private:
        GeoModelBuilder& builder_;
        GeoModel< 3 >& geomodel_;
        GeoModelAccess< 3 > geomodel_access_;

        /*! Internal information */
        std::vector< GeoModelRegionFromSurfaces* > regions_info_;
    };

    /*!
     * @brief Base class to build or edit a GeoModel
     * @details All needed functions are organized in several specific builder
     * in accordance with the kind of edition operation (copy, repair, ...) or
     * with the GeoModel part which is edited (topology, geometry, geology, info)
     */
    class RINGMESH_API GeoModelBuilder {
    ringmesh_disable_copy( GeoModelBuilder );

    public:
        GeoModelBuilder( GeoModel< 3 >& geomodel );
        virtual ~GeoModelBuilder() = default;

        /*!
         * @brief Finish up geomodel building and complete missing information.
         */
        void end_geomodel();

    public:
        GeoModelBuilderTopology< 3 > topology;
        GeoModelBuilderGeometry< 3 > geometry;
        GeoModelBuilderGeology< 3 > geology;
        GeoModelBuilderRemoval< 3 > removal;
        GeoModelBuilderRepair repair;
        GeoModelBuilderCopy copy;
        GeoModelBuilderInfo info;
        GeoModelBuilderFromSurfaces from_surfaces;

    protected:
        GeoModel< 3 >& geomodel_;
        GeoModelAccess< 3 > geomodel_access_;
    };

    /*!
     * @brief Abstract interface class to load and build GeoModels from files
     */
    class RINGMESH_API GeoModelBuilderFile: public GeoModelBuilder {
    public:
        GeoModelBuilderFile( GeoModel< 3 >& geomodel, const std::string& filename );

        virtual ~GeoModelBuilderFile() = default;

        void build_geomodel()
        {
            load_file();
            end_geomodel();
        }

    private:
        virtual void load_file() = 0;

    protected:
        std::string filename_;
    };
}
