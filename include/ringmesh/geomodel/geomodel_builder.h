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

#ifndef __RINGMESH_GEOMODEL_BUILDER__
#define __RINGMESH_GEOMODEL_BUILDER__

#include <ringmesh/basic/common.h>

#include <ringmesh/geomodel/geomodel_editor.h>

/*!
 * @file ringmesh/geomodel_builder.h
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

    /*!
     * @brief Base class for all classes building a GeoModel.
     * @details Derive from this class to build or modify a GeoModel. 
     * @note NON Geometry related modifications are in GeoModelEditor class.
     * @todo To refactor and rename. We need a GeoModelTopologyEditor 
     * and a GeoModelGeometryEditor
     */
    class RINGMESH_API GeoModelBuilder: public GeoModelEditor {
    public:
        GeoModelBuilder( GeoModel& geomodel )
            : GeoModelEditor( geomodel ), options_()
        {
        }
        virtual ~GeoModelBuilder() ;

        /*!
         * @todo Implement it so that it returns true if the input options are consistent
         */
        void set_options( const GeoModelBuildingFlags& options )
        {
            options_ = options ;
        }

        void copy( const GeoModel& from )
        {
            copy_macro_topology( from ) ;
            copy_meshes( from ) ;
        }
        /*!
         * @brief Copy all entity meshes from the input geomodel
         * @pre The geomodel under construction has exaclty the same number of entities
         * than the input geomodel.
         */
        void copy_meshes( const GeoModel& from ) ;
        void copy_meshes( const GeoModel& from, const std::string& entity_type ) ;
        void copy_mesh( const GeoModel& from, const gme_t& mesh_entity ) ;

        void assign_mesh_to_entity( const MeshBase& mesh, const gme_t& to ) ;

        /*!
         * \name Set entity geometry from geometrical positions
         * @{
         */
        /*!
         * @brief Sets a vertex coordinates of a GeoModelMeshEntity
         * @param[in] entity_id the entity to edit
         * @param[in] v the index of the vertex in the entity
         * @param[in] point the coordinates to set
         * @param[in] update if true, updates all the colocated vertices
         * to the new coordinates (ie if edit a Corner coordinates, it will updates
         * its Lines, Surfaces...)
         */
        void set_mesh_entity_vertex(
            const gme_t& entity_id,
            index_t v,
            const vec3& point,
            bool update ) ;

        void set_mesh_entity_vertices(
            const gme_t& entity_id,
            const std::vector< vec3 >& points,
            bool clear ) ;

        /*!
         * @brief Sets the coordinates of a given existing Corner
         * @param[in] corner_id the index of the corner in the GeoModel
         * @param[in] point the coordinates to set
         */
        void set_corner( index_t corner_id, const vec3& point ) ;
        /*!
         * @brief Sets the mesh of a given existing Line
         * @param[in] line_id the index of the line in the GeoModel
         * @param[in] vertices the coordinates to set
         * @warning the vertices should be ordered from the first boundary
         * corner to the second one
         */
        void set_line( index_t line_id, const std::vector< vec3 >& vertices ) ;
        /*!
         * @brief Sets the mesh of a given existing Surface
         * @param[in] surface_id the index of the surface in the GeoModel
         * @param[in] surface_vertices the coordinates to set
         * @param[in] surface_facets the vertex indices of the facets
         * corresponding to \p surface_vertices
         * @param[in] surface_facet_ptr the index of each new facet start in \p surface_facets
         */
        void set_surface_geometry(
            index_t surface_id,
            const std::vector< vec3 >& surface_vertices,
            const std::vector< index_t >& surface_facets,
            const std::vector< index_t >& surface_facet_ptr ) ;
        /*!
         * @brief Sets the tetrahedral mesh of a given existing Region
         * @param[in] region_id the index of the region in the GeoModel
         * @param[in] points the coordinates to set
         * @param[in] tetras the vertex indices of the cells (to read 4 by 4)
         * corresponding to \p points
         */
        void set_region_geometry(
            index_t region_id,
            const std::vector< vec3 >& points,
            const std::vector< index_t >& tetras ) ;

        /*! @}
         * \name Set entity geometry using global GeoModel vertices
         * @{
         */
        void set_mesh_entity_vertex(
            const gme_t& id,
            index_t v,
            index_t geomodel_vertex ) ;

        void set_mesh_entity_vertices(
            const gme_t& entity_id,
            const std::vector< index_t >& geomodel_vertices,
            bool clear ) ;

        void set_corner( index_t corner_id, index_t unique_vertex ) ;

        void set_line( index_t id, const std::vector< index_t >& unique_vertices ) ;

        void set_surface_geometry(
            index_t surface_id,
            const std::vector< index_t >& surface_vertices,
            const std::vector< index_t >& surface_facets,
            const std::vector< index_t >& surface_facet_ptr ) ;

        void set_surface_geometry(
            index_t surface_id,
            const std::vector< index_t >& corners,
            const std::vector< index_t >& facet_ptr ) ;

        void set_surface_geometry(
            index_t surface_id,
            const std::vector< index_t >& triangle_corners ) ;

        void set_surface_geometry_with_adjacencies(
            index_t surface_id,
            const std::vector< index_t >& triangle_corners,
            const std::vector< index_t >& adjacent_triangles ) ;

        void set_surface_element_geometry(
            index_t surface_id,
            index_t facet_id,
            const std::vector< index_t >& corners ) ;

        void set_surface_element_adjacency(
            index_t surface_id,
            index_t facet_id,
            const std::vector< index_t >& adjacents ) ;

        void set_region_geometry(
            index_t region_id,
            const std::vector< index_t >& tet_corners ) ;

        void set_region_element_geometry(
            index_t region_id,
            index_t cell_id,
            const std::vector< index_t >& corners ) ;

        /*! @}
         * \name Create entity element
         * @{
         */

        index_t create_mesh_entity_vertices(
            const gme_t& entity_id,
            index_t nb_vertices ) ;

        index_t create_surface_facet(
            index_t surface_id,
            const std::vector< index_t >& vertex_indices ) ;

        index_t create_region_cell(
            index_t region_id,
            GEO::MeshCellType type,
            const std::vector< index_t >& vertex_indices ) ;

        index_t create_region_cells(
            index_t region_id,
            GEO::MeshCellType type,
            index_t nb_cells ) ;

        /*! @}
         * \name Delete mesh element entities
         * @{
         */

        void delete_mesh_entity_mesh( const gme_t& E_id ) ;
        void delete_mesh_entity_isolated_vertices( const gme_t& E_id ) ;
        void delete_mesh_entity_vertices(
            const gme_t& E_id,
            const std::vector< bool >& to_delete ) ;
        void delete_corner_vertex( index_t corner_id ) ;
        void delete_line_edges(
            index_t line_id,
            const std::vector< bool >& to_delete,
            bool remove_isolated_vertices ) ;
        void delete_surface_facets(
            index_t surface_id,
            const std::vector< bool >& to_delete,
            bool remove_isolated_vertices ) ;
        void delete_region_cells(
            index_t region_id,
            const std::vector< bool >& to_delete,
            bool remove_isolated_vertices ) ;

        /*! @}
         * \name Misc
         * @{
         */

        void compute_surface_adjacencies(
            index_t surface_id,
            bool recompute_adjacency = true ) ;
        void compute_region_adjacencies(
            index_t region_id,
            bool recompute_adjacency = true ) ;
        void triangulate_surface(
            const RINGMesh::Surface& surface_in,
            index_t surface_out ) ;

        gme_t find_or_create_corner( const vec3& point ) ;
        gme_t find_or_create_corner( index_t geomodel_point_id ) ;
        gme_t find_or_create_line( const std::vector< vec3 >& vertices ) ;
        gme_t find_or_create_line(
            const std::vector< index_t >& incident_surfaces,
            const gme_t& first_corner,
            const gme_t& second_corner ) ;

        void recompute_geomodel_mesh() ;

        /*!
         * @}
         * \name Model building functions
         */

        /*!
         * @brief From the Surfaces of the GeoModel, build its Lines and Corners
         */
        bool build_lines_and_corners_from_surfaces() ;

        /*!
         * @brief Build the regions of the GeoModel from the Surfaces
         * @pre Function build_lines_and_corners_from_surfaces must have been called before
         */
        bool build_brep_regions_from_surfaces() ;

        /*
         * @brief From a GeoModel in which only Surface are defined, create corners, contacts
         * and regions depending on the building flags
         * @note Valdity is not checked
         */
        void build_geomodel_from_surfaces() ;

        /*!
         * @brief Finish up geomodel building and complete missing information.
         */
        void end_geomodel() ;

    protected:
        void set_surface_facet_adjacencies(
            index_t surface_id,
            const std::vector< index_t >& facets_ids,
            const std::vector< index_t >& edges_ids,
            const std::vector< index_t >& adjacent_triangles ) ;

    protected:
        /*! Options to toggle the building of entities from the available entities */
        GeoModelBuildingFlags options_ ;

        /*! Internal information */
        std::vector< GeoModelRegionFromSurfaces* > regions_info_ ;

    private:
        void assign_surface_mesh_facets(
            index_t surface_id,
            const std::vector< index_t >& facets,
            const std::vector< index_t >& facet_ptr ) ;

        void assign_surface_triangle_mesh(
            index_t surface_id,
            const std::vector< index_t >& triangle_vertices ) ;
        void update_facet_vertices_around_facet_vertex(
            Surface& surface,
            index_t facet,
            index_t old_vertex,
            index_t new_vertex ) ;
        void update_facet_vertex(
            Surface& surface,
            const std::vector< index_t >& facets,
            index_t old_vertex,
            index_t new_vertex ) ;
        void update_cell_vertex(
            Region& region,
            const std::vector< index_t >& cells,
            index_t old_vertex,
            index_t new_vertex ) ;
        void assign_surface_triangle_mesh(
            index_t surface_id,
            const std::vector< index_t >& triangle_vertices,
            const std::vector< index_t >& adjacent_triangles ) ;

        void assign_region_tet_mesh(
            index_t region_id,
            const std::vector< index_t >& tet_vertices ) ;

        void compute_universe() ;

        void cut_surfaces_by_internal_lines() ;
        void cut_regions_by_internal_surfaces() ;

        void cut_surface_by_line( index_t surface_id, index_t line_id ) ;
        void cut_region_by_surface( index_t region_id, index_t surface_id ) ;
        void duplicate_surface_vertices_along_line(
            index_t surface_id,
            index_t line_id ) ;
        void duplicate_region_vertices_along_surface(
            index_t region_id,
            index_t surface_id ) ;
        index_t disconnect_surface_facets_along_line_edges(
            index_t surface_id,
            index_t line_id ) ;
        index_t disconnect_region_cells_along_surface_facets(
            index_t region_id,
            index_t surface_id ) ;
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

namespace RINGMesh {
    class GeoModelBuilder2 ;
}

namespace RINGMesh {

    class UniverseAccess {
    ringmesh_disable_copy( UniverseAccess ) ;
        friend class GeoModelBuilderTopology ;

    private:
        UniverseAccess( Universe& universe )
            : universe_( universe )
        {
        }

        std::vector< gme_t >& modifiable_boundaries()
        {
            return universe_.boundary_surfaces_ ;
        }

        std::vector< bool >& modifiable_sides()
        {
            return universe_.boundary_surface_sides_ ;
        }

        void copy( const Universe& from )
        {
            universe_.copy( from ) ;
        }

    private:
        Universe& universe_ ;
    } ;

    class GeoModelMeshEntityConstAccess {
    ringmesh_disable_copy( GeoModelMeshEntityConstAccess ) ;
        friend class GeoModelBuilderGeometry ;

    private:
        GeoModelMeshEntityConstAccess( const GeoModelMeshEntity& gme )
            : gmme_( gme )
        {
        }

        const MeshBase* mesh() const
        {
            return gmme_.mesh_ ;
        }

    private:
        const GeoModelMeshEntity& gmme_ ;
    } ;

    class GeoModelMeshEntityAccess {
    ringmesh_disable_copy( GeoModelMeshEntityAccess ) ;
        friend class GeoModelBuilderTopology ;
        friend class GeoModelBuilderGeometry ;
        friend class GeoModelBuilderGeology ;
        friend class GeoModelBuilderInfo ;

    private:
        GeoModelMeshEntityAccess( GeoModelMeshEntity& gme )
            : gmme_( gme )
        {
        }

        std::string& modifiable_name()
        {
            return gmme_.name_ ;
        }

        GME::GEOL_FEATURE& modifiable_geol_feature()
        {
            return gmme_.geol_feature_ ;
        }

        std::vector< gme_t >& modifiable_boundaries()
        {
            return gmme_.boundaries_ ;
        }

        std::vector< gme_t >& modifiable_in_boundaries()
        {
            return gmme_.in_boundary_ ;
        }

        std::vector< bool >& modifiable_sides()
        {
            ringmesh_assert( gmme_.type_name() == Region::type_name_static() ) ;
            return dynamic_cast< Region& >( gmme_ ).sides_ ;
        }

        std::vector< gme_t >& modifiable_parents()
        {
            return gmme_.parents_ ;
        }

        MeshBase* modifiable_mesh()
        {
            return gmme_.mesh_ ;
        }

        template< typename ENTITY >
        static ENTITY* create_entity(
            const GeoModel& geomodel,
            index_t id,
            const MeshType type )
        {
            return new ENTITY( geomodel, id, type ) ;
        }

        void copy( const GeoModelMeshEntity& from )
        {
            gmme_.copy( from ) ;
        }

    private:
        GeoModelMeshEntity& gmme_ ;
    } ;

    class GeoModelGeologicalEntityAccess {
    ringmesh_disable_copy( GeoModelGeologicalEntityAccess ) ;
        friend class GeoModelBuilderTopology ;
        friend class GeoModelBuilderGeology ;
        friend class GeoModelBuilderInfo ;

    private:
        GeoModelGeologicalEntityAccess( GeoModelGeologicalEntity& gmge )
            : gmge_( gmge )
        {
        }

        std::string& modifiable_name()
        {
            return gmge_.name_ ;
        }

        GME::GEOL_FEATURE& modifiable_geol_feature()
        {
            return gmge_.geol_feature_ ;
        }

        std::vector< gme_t >& modifiable_children()
        {
            return gmge_.children_ ;
        }

        static GeoModelGeologicalEntity* create_geological_entity(
            const EntityType& type,
            const GeoModel& geomodel,
            index_t index_in_geomodel )
        {
            GeoModelGeologicalEntity* E =
                GeoModelGeologicalEntityFactory::create_object( type, geomodel ) ;
            E->id_.index = index_in_geomodel ;
            return E ;
        }

        void copy( const GeoModelGeologicalEntity& from )
        {
            gmge_.copy( from ) ;
        }

    private:
        GeoModelGeologicalEntity& gmge_ ;
    } ;

    class GeoModelAccess {
    ringmesh_disable_copy( GeoModelAccess ) ;
        friend class GeoModelBuilder2 ;
        friend class GeoModelBuilderTopology ;
        friend class GeoModelBuilderGeometry ;
        friend class GeoModelBuilderGeology ;
        friend class GeoModelBuilderRemoval ;
        friend class GeoModelBuilderCopy ;
        friend class GeoModelBuilderInfo ;

    private:
        GeoModelAccess( GeoModel& geomodel )
            : geomodel_( geomodel )
        {
        }

//        const EntityTypeManager& entity_type_manager() const
//        {
//            return geomodel_.entity_type_manager_
//        }

        std::string& modifiable_name()
        {
            return geomodel_.geomodel_name_ ;
        }

        EntityTypeManager& modifiable_entity_type_manager()
        {
            return geomodel_.entity_type_manager_ ;
        }

        std::vector< GeoModelMeshEntity* >& modifiable_mesh_entities(
            const EntityType& type )
        {
            return const_cast< std::vector< GeoModelMeshEntity* >& >( geomodel_.mesh_entities(
                type ) ) ;
        }

        GeoModelMeshEntity& modifiable_mesh_entity( const gme_t& id )
        {
            return *modifiable_mesh_entities( id.type )[id.index] ;
        }

        std::vector< std::vector< GeoModelGeologicalEntity* > > modifiable_geological_entities()
        {
            return geomodel_.geological_entities_ ;
        }

        std::vector< GeoModelGeologicalEntity* >& modifiable_geological_entities(
            const EntityType& type )
        {
            return const_cast< std::vector< GeoModelGeologicalEntity* >& >( geomodel_.geological_entities(
                type ) ) ;
        }

        GeoModelGeologicalEntity& modifiable_geological_entity( const gme_t& id )
        {
            return *modifiable_geological_entities( id.type )[id.index] ;
        }

        Universe& modifiable_universe()
        {
            return geomodel_.universe_ ;
        }

        double& modifiable_epsilon()
        {
            return geomodel_.epsilon_ ;
        }

    private:
        GeoModel& geomodel_ ;
    } ;

    class RINGMESH_API GeoModelBuilderTopology {
    ringmesh_disable_copy( GeoModelBuilderTopology ) ;
        friend class GeoModelBuilder2 ;

    public:
        /*!
         * @brief Copy topological information from a geomodel
         * @details Copy all the geomodel entities and their relationship
         * ignoring their geometry
         * @param[in] from Model to copy the information from
         */
        void copy_topology( const GeoModel& from ) ;

        template< typename T >
        gme_t create_mesh_entity( const MeshType type = "" )
        {
            const EntityType entity_type = T::type_name_static() ;
            index_t nb_entities( geomodel_.nb_mesh_entities( entity_type ) ) ;
            index_t new_id( nb_entities ) ;
            T* new_entity = new T( geomodel_, new_id, type ) ;
            geomodel_access_.modifiable_mesh_entities( entity_type ).push_back(
                new_entity ) ;
            return new_entity->gme_id() ;
        }

        /*!
         * @brief Complete missing information in GeoModelEntities
         * boundaries - in_boundary - parent - children
         * @details For all 7 types of entities, check what information is available
         * for the first one and fill the entities of the same type accordingly
         * THIS MEANS that the all the entities of the same type have been initialized with
         * the same information
         */
        void complete_entity_connectivity() ;

        /*!
         * @brief Fill the boundaries of all entities of the given type
         * @details If the boundary entities do not have any in_boundary
         * information, nothing is done.
         */
        void fill_mesh_entities_boundaries( const EntityType& type ) ;

        /*!
         * @brief Fill the in_boundary vector of all entities of the given type
         * @details If the in_boundary entities do not have any boundary
         * information, nothing is done, and geomodel construction will eventually fail.
         */
        void fill_mesh_entities_in_boundaries( const EntityType& type ) ;

        void add_mesh_entity_boundary(
            const gme_t& gme_id,
            index_t boundary_id,
            bool side = false )
        {
            GeoModelMeshEntity& mesh_entity =
                geomodel_access_.modifiable_mesh_entity( gme_id ) ;
            const EntityType& b_type = geomodel_.entity_type_manager().boundary_type(
                gme_id.type ) ;
            gme_t boundary( b_type, boundary_id ) ;
            GeoModelMeshEntityAccess gme_access( mesh_entity ) ;
            gme_access.modifiable_boundaries().push_back( boundary ) ;

            if( gme_id.type == Region::type_name_static() ) {
                gme_access.modifiable_sides().push_back( side ) ;
            }
        }

        void set_mesh_entity_boundary(
            const gme_t& gme_id,
            index_t id,
            index_t boundary_id,
            bool side = false )
        {
            ringmesh_assert( id < geomodel_.mesh_entity( gme_id ).nb_boundaries() ) ;
            GeoModelMeshEntity& mesh_entity =
                geomodel_access_.modifiable_mesh_entity( gme_id ) ;
            const EntityType& b_type = geomodel_.entity_type_manager().boundary_type(
                gme_id.type ) ;
            gme_t boundary( b_type, boundary_id ) ;
            GeoModelMeshEntityAccess gme_access( mesh_entity ) ;
            gme_access.modifiable_boundaries()[id] = boundary ;

            if( gme_id.type == Region::type_name_static() ) {
                gme_access.modifiable_sides()[id] = side ;
            }
        }

        void add_universe_boundary( index_t boundary_id, bool side )
        {
            gme_t boundary( Surface::type_name_static(), boundary_id ) ;
            UniverseAccess universe_access(
                geomodel_access_.modifiable_universe() ) ;
            universe_access.modifiable_boundaries().push_back( boundary ) ;
            universe_access.modifiable_sides().push_back( side ) ;
        }

        void set_universe_boundary( index_t id, index_t boundary_id, bool side )
        {
            ringmesh_assert( id < geomodel_.universe().nb_boundaries() ) ;
            gme_t boundary( Surface::type_name_static(), boundary_id ) ;
            UniverseAccess universe_access(
                geomodel_access_.modifiable_universe() ) ;
            universe_access.modifiable_boundaries()[id] = boundary ;
            universe_access.modifiable_sides()[id] = side ;
        }

        void add_mesh_entity_in_boundary( const gme_t& t, index_t in_boundary_id )
        {
            GeoModelMeshEntity& mesh_entity =
                geomodel_access_.modifiable_mesh_entity( t ) ;
            const EntityType& in_b_type =
                geomodel_.entity_type_manager().in_boundary_type( t.type ) ;
            gme_t in_boundary( in_b_type, in_boundary_id ) ;
            GeoModelMeshEntityAccess gme_access( mesh_entity ) ;
            gme_access.modifiable_in_boundaries().push_back( in_boundary ) ;
        }

        void set_mesh_entity_in_boundary(
            const gme_t& gme_id,
            index_t id,
            index_t in_boundary_id )
        {
            /// No check on the validity of the index of the entity in_boundary
            /// NO_ID is used to flag entities to delete
            GeoModelMeshEntity& mesh_entity =
                geomodel_access_.modifiable_mesh_entity( gme_id ) ;
            ringmesh_assert( id < mesh_entity.nb_in_boundary() ) ;
            const EntityType& in_b_type =
                geomodel_.entity_type_manager().in_boundary_type( gme_id.type ) ;
            gme_t in_boundary( in_b_type, in_boundary_id ) ;
            GeoModelMeshEntityAccess gme_access( mesh_entity ) ;
            gme_access.modifiable_in_boundaries()[id] = in_boundary ;
        }

    private:
        GeoModelBuilderTopology( GeoModel& builder ) ;

        bool create_mesh_entities(
            const EntityType& type,
            index_t nb_additional_entities ) ;

        template< typename ENTITY >
        bool create_mesh_entities(
            index_t nb_additionnal_entities,
            const MeshType type = "" )
        {
            const EntityType entity_type = ENTITY::type_name_static() ;
            std::vector< GeoModelMeshEntity* >& store =
                geomodel_access_.modifiable_mesh_entities( entity_type ) ;
            index_t old_size = static_cast< index_t >( store.size() ) ;
            index_t new_size = old_size + nb_additionnal_entities ;
            store.resize( new_size, nil ) ;
            for( index_t i = old_size; i < new_size; i++ ) {
                store[i] = GeoModelMeshEntityAccess::create_entity< ENTITY >(
                    geomodel_, i, type ) ;
            }
            return true ;
        }
        bool create_geological_entities( const EntityType& type, index_t nb ) ;
        index_t find_or_create_geological_entity_type( const EntityType& type )
        {
            index_t type_index =
                geomodel_.entity_type_manager().geological_entity_type_index(
                    type ) ;
            if( type_index == NO_ID ) {
                type_index = create_geological_entity_type( type ) ;
            }
            return type_index ;
        }
        index_t create_geological_entity_type( const EntityType& type ) ;

        template< typename T >
        void complete_mesh_entity_connectivity() ;

        template< typename T >
        void copy_mesh_entity_topology( const GeoModel& from ) ;
        void copy_geological_entity_topology(
            const GeoModel& from,
            const EntityType& type ) ;

    private:
        GeoModel& geomodel_ ;
        GeoModelAccess geomodel_access_ ;
    } ;

    class RINGMESH_API GeoModelBuilderRemoval {
    ringmesh_disable_copy( GeoModelBuilderRemoval ) ;
        friend class GeoModelBuilder2 ;

    public:
        void remove_mesh_entities( const std::set< gme_t >& entities ) ;

        void remove_geological_entities( const std::set< gme_t >& entities ) ;

        /*!
         * @todo Could be moved in the API [JP]
         */
        bool get_dependent_entities( std::set< gme_t >& entities ) const ;

        /*!
         * Should be rewritten. Put as it was before someone removed it...
         */
        void remove_entities_and_dependencies(
            const std::set< gme_t >& entities_to_remove ) ;

    protected:
        GeoModelBuilderRemoval( GeoModel& builder ) ;

    private:
        GeoModel& geomodel_ ;
        GeoModelAccess geomodel_access_ ;
    } ;

    class RINGMESH_API GeoModelBuilderGeometry {
    ringmesh_disable_copy( GeoModelBuilderGeometry ) ;
        friend class GeoModelBuilder2 ;

    public:
        /*!
         * @brief Transfer general mesh information from one mesh
         * data structure to another one
         * @param[in] id the GeoModelMeshEntity id to operate on
         * @param[in] type the new mesh data structure type
         */
        void change_mesh_data_structure( const gme_t& id, const MeshType type ) ;

        /*!
         * @brief Copy all entity meshes from the input geomodel
         * @pre The geomodel under construction has exaclty the same number of entities
         * than the input geomodel.
         */
        void copy_meshes( const GeoModel& from ) ;

        /*!
         * \name Set entity geometry from geometrical positions
         * @{
         */
        /*!
         * @brief Sets a vertex coordinates of a GeoModelMeshEntity
         * @param[in] entity_id the entity to edit
         * @param[in] v the index of the vertex in the entity
         * @param[in] point the coordinates to set
         * @param[in] update if true, updates all the colocated vertices
         * to the new coordinates (ie if edit a Corner coordinates, it will updates
         * its Lines, Surfaces...)
         */
        void set_mesh_entity_vertex(
            const gme_t& entity_id,
            index_t v,
            const vec3& point,
            bool update ) ;

        void set_mesh_entity_vertices(
            const gme_t& entity_id,
            const std::vector< vec3 >& points,
            bool clear ) ;

        /*!
         * @brief Sets the coordinates of a given existing Corner
         * @param[in] corner_id the index of the corner in the GeoModel
         * @param[in] point the coordinates to set
         */
        void set_corner( index_t corner_id, const vec3& point ) ;
        /*!
         * @brief Sets the mesh of a given existing Line
         * @param[in] line_id the index of the line in the GeoModel
         * @param[in] vertices the coordinates to set
         * @warning the vertices should be ordered from the first boundary
         * corner to the second one
         */
        void set_line( index_t line_id, const std::vector< vec3 >& vertices ) ;
        /*!
         * @brief Sets the mesh of a given existing Surface
         * @param[in] surface_id the index of the surface in the GeoModel
         * @param[in] surface_vertices the coordinates to set
         * @param[in] surface_facets the vertex indices of the facets
         * corresponding to \p surface_vertices
         * @param[in] surface_facet_ptr the index of each new facet start in \p surface_facets
         */
        void set_surface_geometry(
            index_t surface_id,
            const std::vector< vec3 >& surface_vertices,
            const std::vector< index_t >& surface_facets,
            const std::vector< index_t >& surface_facet_ptr ) ;
        /*!
         * @brief Sets the tetrahedral mesh of a given existing Region
         * @param[in] region_id the index of the region in the GeoModel
         * @param[in] points the coordinates to set
         * @param[in] tetras the vertex indices of the cells (to read 4 by 4)
         * corresponding to \p points
         */
        void set_region_geometry(
            index_t region_id,
            const std::vector< vec3 >& points,
            const std::vector< index_t >& tetras ) ;

        /*! @}
         * \name Set entity geometry using global GeoModel vertices
         * @{
         */
        void set_mesh_entity_vertex(
            const gme_t& id,
            index_t v,
            index_t geomodel_vertex ) ;

        void set_mesh_entity_vertices(
            const gme_t& entity_id,
            const std::vector< index_t >& geomodel_vertices,
            bool clear ) ;

        void set_corner( index_t corner_id, index_t unique_vertex ) ;

        void set_line( index_t id, const std::vector< index_t >& unique_vertices ) ;

        void set_surface_geometry(
            index_t surface_id,
            const std::vector< index_t >& surface_vertices,
            const std::vector< index_t >& surface_facets,
            const std::vector< index_t >& surface_facet_ptr ) ;

        void set_surface_geometry(
            index_t surface_id,
            const std::vector< index_t >& corners,
            const std::vector< index_t >& facet_ptr ) ;

        void set_surface_geometry(
            index_t surface_id,
            const std::vector< index_t >& triangle_corners ) ;

        void set_surface_geometry_with_adjacencies(
            index_t surface_id,
            const std::vector< index_t >& triangle_corners,
            const std::vector< index_t >& adjacent_triangles ) ;

        void set_surface_element_geometry(
            index_t surface_id,
            index_t facet_id,
            const std::vector< index_t >& corners ) ;

        void set_surface_element_adjacency(
            index_t surface_id,
            index_t facet_id,
            const std::vector< index_t >& adjacents ) ;

        void set_region_geometry(
            index_t region_id,
            const std::vector< index_t >& tet_corners ) ;

        void set_region_element_geometry(
            index_t region_id,
            index_t cell_id,
            const std::vector< index_t >& corners ) ;

        /*! @}
         * \name Create entity element
         * @{
         */

        index_t create_mesh_entity_vertices(
            const gme_t& entity_id,
            index_t nb_vertices ) ;

        index_t create_surface_facet(
            index_t surface_id,
            const std::vector< index_t >& vertex_indices ) ;

        index_t create_region_cell(
            index_t region_id,
            GEO::MeshCellType type,
            const std::vector< index_t >& vertex_indices ) ;

        index_t create_region_cells(
            index_t region_id,
            GEO::MeshCellType type,
            index_t nb_cells ) ;

        /*! @}
         * \name Delete mesh element entities
         * @{
         */

        void delete_mesh_entity_mesh( const gme_t& E_id ) ;
        void delete_mesh_entity_isolated_vertices( const gme_t& E_id ) ;
        void delete_mesh_entity_vertices(
            const gme_t& E_id,
            const std::vector< bool >& to_delete ) ;
        void delete_corner_vertex( index_t corner_id ) ;
        void delete_line_edges(
            index_t line_id,
            const std::vector< bool >& to_delete,
            bool remove_isolated_vertices ) ;
        void delete_surface_facets(
            index_t surface_id,
            const std::vector< bool >& to_delete,
            bool remove_isolated_vertices ) ;
        void delete_region_cells(
            index_t region_id,
            const std::vector< bool >& to_delete,
            bool remove_isolated_vertices ) ;

        /*! @}
         * \name Misc
         * @{
         */

        void compute_surface_adjacencies(
            index_t surface_id,
            bool recompute_adjacency = true ) ;
        void compute_region_adjacencies(
            index_t region_id,
            bool recompute_adjacency = true ) ;
        void triangulate_surface(
            const RINGMesh::Surface& surface_in,
            index_t surface_out ) ;

    protected:
        GeoModelBuilderGeometry( GeoModel& builder ) ;

    private:
        void copy_meshes( const GeoModel& from, const std::string& entity_type ) ;
        void copy_mesh( const GeoModel& from, const gme_t& mesh_entity ) ;

        void assign_mesh_to_entity( const MeshBase& mesh, const gme_t& to ) ;

    private:
        GeoModel& geomodel_ ;
        GeoModelAccess geomodel_access_ ;
    } ;

    class RINGMESH_API GeoModelBuilderInfo {
    ringmesh_disable_copy( GeoModelBuilderInfo ) ;
        friend class GeoModelBuilder2 ;

    public:
        /*!
         *@brief Set the name of the geomodel
         */
        void set_geomodel_name( const std::string& name )
        {
            geomodel_access_.modifiable_name() = name ;
        }

        void set_entity_name( const gme_t& gme_id, const std::string& name )
        {
            if( geomodel_.is_mesh_entity_type( gme_id.type ) ) {
                GeoModelMeshEntityAccess gmme_access(
                    geomodel_access_.modifiable_mesh_entity( gme_id ) ) ;
                gmme_access.modifiable_name() = name ;
            } else {
                GeoModelGeologicalEntityAccess gmge_access(
                    geomodel_access_.modifiable_geological_entity( gme_id ) ) ;
                gmge_access.modifiable_name() = name ;
            }
        }

    protected:
        GeoModelBuilderInfo( GeoModel& builder ) ;

    private:
        GeoModel& geomodel_ ;
        GeoModelAccess geomodel_access_ ;

    } ;

    class RINGMESH_API GeoModelBuilderGeology {
    ringmesh_disable_copy( GeoModelBuilderGeology ) ;
        friend class GeoModelBuilder2 ;

    public:
        gme_t create_geological_entity( const EntityType& type ) ;

        /*!
         * @brief Fill the parent of all entities of the given type
         * @details If the parents do not have any child nothing is done.
         */
        void fill_mesh_entities_parent( const EntityType& type ) ;

        /*!
         * @brief Fill the children of all entities of the given type
         * @details If the children entities do not have any parent information
         * nothing is done.
         */
        void fill_geological_entities_children( const EntityType& type ) ;

        void complete_mesh_entities_geol_feature_from_first_parent(
            const EntityType& type ) ;
        void complete_geological_entities_geol_feature_from_first_child(
            const EntityType& type ) ;

        void set_entity_geol_feature(
            const gme_t& gme_id,
            GME::GEOL_FEATURE geol_feature )
        {
            if( geomodel_.is_mesh_entity_type( gme_id.type ) ) {
                GeoModelMeshEntityAccess gmme_access(
                    geomodel_access_.modifiable_mesh_entity( gme_id ) ) ;
                gmme_access.modifiable_geol_feature() = geol_feature ;
            } else {
                GeoModelGeologicalEntityAccess gmge_access(
                    geomodel_access_.modifiable_geological_entity( gme_id ) ) ;
                gmge_access.modifiable_geol_feature() = geol_feature ;
            }
        }

        void add_mesh_entity_parent( const gme_t& gme_id, const gme_t& parent_index )
        {
            GeoModelMeshEntity& mesh_entity =
                geomodel_access_.modifiable_mesh_entity( gme_id ) ;
            GeoModelMeshEntityAccess gmme_access( mesh_entity ) ;
            gmme_access.modifiable_parents().push_back( parent_index ) ;
        }

        void set_mesh_entity_parent(
            const gme_t& gme_id,
            index_t id,
            const gme_t& parent_index )
        {
            /// No check on the validity of the index of the entity parents_
            /// NO_ID is used to flag entities to delete
            GeoModelMeshEntity& mesh_entity =
                geomodel_access_.modifiable_mesh_entity( gme_id ) ;
            ringmesh_assert( id < mesh_entity.nb_parents() ) ;
            GeoModelMeshEntityAccess gmme_access( mesh_entity ) ;
            gmme_access.modifiable_parents()[id] = parent_index ;
        }

        void add_geological_entity_child( const gme_t& gme_id, index_t child_id )
        {
            GeoModelGeologicalEntity& geol_entity =
                geomodel_access_.modifiable_geological_entity( gme_id ) ;
            const EntityType& child_type =
                geomodel_.entity_type_manager().child_type( gme_id.type ) ;
            gme_t child( child_type, child_id ) ;
            GeoModelGeologicalEntityAccess gmge_access( geol_entity ) ;
            gmge_access.modifiable_children().push_back( child ) ;
        }

        void set_geological_entity_child(
            const gme_t& gme_id,
            index_t id,
            index_t child_id )
        {
            /// No check on the validity of the index of the entity child_index
            /// NO_ID is used to flag entities to delete
            GeoModelGeologicalEntity& geol_entity =
                geomodel_access_.modifiable_geological_entity( gme_id ) ;
            const EntityType& child_type =
                geomodel_.entity_type_manager().child_type( gme_id.type ) ;
            gme_t child( child_type, child_id ) ;
            GeoModelGeologicalEntityAccess gmge_access( geol_entity ) ;
            gmge_access.modifiable_children()[id] = child ;
        }

    protected:
        GeoModelBuilderGeology( GeoModel& builder ) ;

    private:
        GeoModel& geomodel_ ;
        GeoModelAccess geomodel_access_ ;

    } ;

    class RINGMESH_API GeoModelBuilderCopy {
    ringmesh_disable_copy( GeoModelBuilderCopy ) ;
        friend class GeoModelBuilder2 ;

    private:
        GeoModelBuilderCopy( GeoModelBuilder2& builder, GeoModel& geomodel ) ;

        void copy_geomodel( const GeoModel& from ) ;

    private:
        GeoModelBuilder2& builder_ ;
        GeoModel& geomodel_ ;
        GeoModelAccess geomodel_access_ ;

    } ;

    class RINGMESH_API GeoModelBuilder2 {
    ringmesh_disable_copy( GeoModelBuilder2 ) ;

    public:
        GeoModelBuilder2( GeoModel& geomodel ) ;
        virtual ~GeoModelBuilder2()
        {
        }

//        const GeoModel& geomodel() const
//        {
//            return geomodel_ ;
//        }

    private:
//        /*! The geomodel edited
//         */
//        GeoModel& geomodel_ ;

//        /*! Parameter to forbid element creation. Crucial to control
//         *  building of the geomodel and detect errors in find_or_create functions
//         */
//        bool create_entity_allowed_ ;

    public:
        GeoModelBuilderTopology topology ;
        GeoModelBuilderGeometry geometry ;
        GeoModelBuilderGeology geology ;
        GeoModelBuilderRemoval removal ;
        GeoModelBuilderCopy copy ;
        GeoModelBuilderInfo info ;
    } ;
}

#endif
