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
            compute_corners = false ;
            compute_lines = false ;
            compute_surfaces = false ;
            compute_regions_brep = false ;
            compute_regions_mesh = false ;
        }
        bool compute_corners ;
        bool compute_lines ;
        bool compute_surfaces ;
        bool compute_regions_brep ;
        bool compute_regions_mesh ;
    } ;

    // Implementation details
    class GeoModelRegionFromSurfaces ;
}

namespace RINGMesh {

    class RINGMESH_API GeoModelBuilderInfo {
    ringmesh_disable_copy( GeoModelBuilderInfo ) ;
        friend class GeoModelBuilder ;

    public:
        /*!
         *@brief Set the name of the geomodel
         */
        void set_geomodel_name( const std::string& name )
        {
            geomodel_access_.modifiable_name() = name ;
        }

        /*!
         *@brief Set the name of a geomodel mesh entity
         */
        void set_mesh_entity_name( const gmme_t& gmme_id, const std::string& name )
        {
            GeoModelMeshEntityAccess gmme_access(
                geomodel_access_.modifiable_mesh_entity( gmme_id ) ) ;
            gmme_access.modifiable_name() = name ;

        }

        /*!
         *@brief Set the name of a geomodel geological entity
         */
        void set_geological_entity_name( const gmge_t& gmge_id, const std::string& name )
        {
            GeoModelGeologicalEntityAccess gmge_access(
                geomodel_access_.modifiable_geological_entity( gmge_id ) ) ;
            gmge_access.modifiable_name() = name ;
        }

    protected:
        GeoModelBuilderInfo( GeoModelBuilder& builder, GeoModel& geomodel ) ;

    private:
        GeoModelBuilder& builder_ ;
        GeoModel& geomodel_ ;
        GeoModelAccess geomodel_access_ ;

    } ;

    class RINGMESH_API GeoModelBuilderGeology {
    ringmesh_disable_copy( GeoModelBuilderGeology ) ;
        friend class GeoModelBuilder ;

    public:
        void copy_geology( const GeoModel& from ) ;

        /*!
         * @brief Create and store a geological entity of the given type
         * @return The index of the created geological entity
         */
        gmge_t create_geological_entity( const GeologicalEntityType& type ) ;

        bool create_geological_entities( const GeologicalEntityType& type, index_t nb ) ;

        /*!
         * @brief Fill the parent of all entities of the given type
         * @details If the parents do not have any child nothing is done.
         */
        void fill_mesh_entities_parent( const MeshEntityType& type ) ;

        /*!
         * @brief Fill the children of all entities of the given type
         * @details If the children entities do not have any parent information
         * nothing is done.
         */
        void fill_geological_entities_children( const GeologicalEntityType& type ) ;

        void complete_mesh_entities_geol_feature_from_first_parent(
            const MeshEntityType& type ) ;
        void complete_geological_entities_geol_feature_from_first_child(
            const GeologicalEntityType& type ) ;

        void set_mesh_entity_geol_feature(
            const gmme_t& gmme_id,
            GME::GEOL_FEATURE geol_feature )
        {
            GeoModelMeshEntityAccess gmme_access(
                geomodel_access_.modifiable_mesh_entity( gmme_id ) ) ;
            gmme_access.modifiable_geol_feature() = geol_feature ;

        }

        void set_geological_entity_geol_feature(
            const gmge_t& gmge_id,
            GME::GEOL_FEATURE geol_feature )
        {
            GeoModelGeologicalEntityAccess gmge_access(
                geomodel_access_.modifiable_geological_entity( gmge_id ) ) ;
            gmge_access.modifiable_geol_feature() = geol_feature ;
        }

        void add_mesh_entity_parent( const gmme_t& child_gmme, const gmge_t& parent_gmge )
        {
            GeoModelMeshEntity& mesh_entity =
                geomodel_access_.modifiable_mesh_entity( child_gmme ) ;
            GeoModelMeshEntityAccess gmme_access( mesh_entity ) ;
            DEBUG(child_gmme) ;
            DEBUG(parent_gmge) ;
            DEBUG("la") ;

            gmme_access.modifiable_parents().push_back( parent_gmge ) ;
        }

        void set_mesh_entity_parent(
            const gmme_t& child_gmme,
            index_t id,
            const gmge_t& parent_gmge )
        {
            /// No check on the validity of the index of the entity parents_
            /// NO_ID is used to flag entities to delete
            GeoModelMeshEntity& mesh_entity =
                geomodel_access_.modifiable_mesh_entity( child_gmme ) ;
            ringmesh_assert( id < mesh_entity.nb_parents() ) ;
            GeoModelMeshEntityAccess gmme_access( mesh_entity ) ;
            DEBUG(parent_gmge) ;
            DEBUG("ici") ;
            gmme_access.modifiable_parents()[id] = parent_gmge ;
        }

        void add_geological_entity_child( const gmge_t& parent_gmge, index_t child_id )
        {
            GeoModelGeologicalEntity& geol_entity =
                geomodel_access_.modifiable_geological_entity( parent_gmge ) ;
            const MeshEntityType& child_type =
                geomodel_.entity_type_manager().relationship_manager.child_type( parent_gmge.type() ) ;
            gmme_t child( child_type, child_id ) ;
            GeoModelGeologicalEntityAccess gmge_access( geol_entity ) ;
            gmge_access.modifiable_children().push_back( child ) ;
        }

        void set_geological_entity_child(
            const gmge_t& parent_gmge,
            index_t id,
            index_t child_id )
        {
            /// No check on the validity of the index of the entity child_index
            /// NO_ID is used to flag entities to delete
            GeoModelGeologicalEntity& geol_entity =
                geomodel_access_.modifiable_geological_entity( parent_gmge ) ;
            const MeshEntityType& child_type =
                geomodel_.entity_type_manager().relationship_manager.child_type( parent_gmge.type() ) ;
            gmme_t child( child_type, child_id ) ;
            GeoModelGeologicalEntityAccess gmge_access( geol_entity ) ;
            gmge_access.modifiable_children()[id] = child ;
        }

        void delete_geological_entity( const GeologicalEntityType& type, index_t index ) ;

    protected:
        GeoModelBuilderGeology( GeoModelBuilder& builder, GeoModel& geomodel ) ;

    private:
        index_t create_geological_entity_type( const GeologicalEntityType& type ) ;

        index_t find_or_create_geological_entity_type( const GeologicalEntityType& type ) ;

        void copy_geological_entity_topology(
            const GeoModel& from,
            const GeologicalEntityType& type ) ;

    private:
        GeoModelBuilder& builder_ ;
        GeoModel& geomodel_ ;
        GeoModelAccess geomodel_access_ ;
    } ;

    class RINGMESH_API GeoModelBuilderCopy {
    ringmesh_disable_copy( GeoModelBuilderCopy ) ;
        friend class GeoModelBuilder ;
    public:
        void copy_geomodel( const GeoModel& from ) ;

    private:
        GeoModelBuilderCopy( GeoModelBuilder& builder, GeoModel& geomodel ) ;

    private:
        GeoModelBuilder& builder_ ;
        GeoModel& geomodel_ ;
        GeoModelAccess geomodel_access_ ;
    } ;

    class RINGMESH_API GeoModelBuilderFromSurfaces {
    ringmesh_disable_copy( GeoModelBuilderFromSurfaces ) ;
        friend class GeoModelBuilder ;

    public:
        virtual ~GeoModelBuilderFromSurfaces() ;
        /*
         * @brief From a GeoModel in which only Surfaces are defined,
         * create Corners, Lines and Regions depending on the building flags
         * @note Validity is not checked
         */
        void build() ;

        /*!
         * @brief From the Surfaces of the GeoModel, build its Lines and Corners
         */
        bool build_lines_and_corners_from_surfaces() ;

        /*!
         * @brief Build the regions of the GeoModel from the Surfaces
         * @pre Function build_lines_and_corners_from_surfaces
         * must have been called before
         */
        bool build_brep_regions_from_surfaces() ;

    private:
        GeoModelBuilderFromSurfaces( GeoModelBuilder& builder, GeoModel& geomodel ) ;

    public:
        /*! Options to toggle the building of entities from the available entities */
        GeoModelBuildingFlags options_ ;

    private:
        GeoModelBuilder& builder_ ;
        GeoModel& geomodel_ ;
        GeoModelAccess geomodel_access_ ;

        /*! Internal information */
        std::vector< GeoModelRegionFromSurfaces* > regions_info_ ;
    } ;

    /*!
     * @brief Base class to build or edit a GeoModel
     * @details All needed functions are organized in several specific builder
     * in accordance with the kind of edition operation (copy, repair, ...) or
     * with the GeoModel part which is edited (topology, geometry, geology, info)
     */
    class RINGMESH_API GeoModelBuilder {
    ringmesh_disable_copy( GeoModelBuilder ) ;

    public:
        GeoModelBuilder( GeoModel& geomodel ) ;
        virtual ~GeoModelBuilder()
        {
        }

        /*!
         * @brief Finish up geomodel building and complete missing information.
         */
        void end_geomodel() ;

    public:
        GeoModelBuilderTopology topology ;
        GeoModelBuilderGeometry geometry ;
        GeoModelBuilderGeology geology ;
        GeoModelBuilderRemoval removal ;
        GeoModelBuilderRepair repair ;
        GeoModelBuilderCopy copy ;
        GeoModelBuilderInfo info ;
        GeoModelBuilderFromSurfaces from_surfaces ;

    protected:
        GeoModel& geomodel_ ;
        GeoModelAccess geomodel_access_ ;
    } ;

    /*!
     * @brief Abstract interface class to load and build GeoModels from files
     */
    class RINGMESH_API GeoModelBuilderFile: public GeoModelBuilder {
    public:
        GeoModelBuilderFile( GeoModel& geomodel, const std::string& filename ) ;

        virtual ~GeoModelBuilderFile()
        {
        }
        void build_geomodel()
        {
            load_file() ;
            end_geomodel() ;
        }

    private:
        virtual void load_file() = 0 ;

    protected:
        std::string filename_ ;
    } ;
}
