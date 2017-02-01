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
#include <ringmesh/geomodel/geomodel_builder_repair.h>

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
//    class RINGMESH_API GeoModelBuilder: public GeoModelEditor {
//    public:
//        GeoModelBuilder( GeoModel& geomodel )
//            : GeoModelEditor( geomodel ), options_()
//        {
//        }
//        virtual ~GeoModelBuilder() ;
//
//        /*!
//         * @todo Implement it so that it returns true if the input options are consistent
//         */
//        void set_options( const GeoModelBuildingFlags& options )
//        {
//            options_ = options ;
//        }
//
//        void copy( const GeoModel& from )
//        {
//            copy_macro_topology( from ) ;
//            copy_meshes( from ) ;
//        }
//        /*!
//         * @brief Copy all entity meshes from the input geomodel
//         * @pre The geomodel under construction has exaclty the same number of entities
//         * than the input geomodel.
//         */
//        void copy_meshes( const GeoModel& from ) ;
//        void copy_meshes( const GeoModel& from, const std::string& entity_type ) ;
//        void copy_mesh( const GeoModel& from, const gme_t& mesh_entity ) ;
//
//        void assign_mesh_to_entity( const MeshBase& mesh, const gme_t& to ) ;
//
//        /*!
//         * \name Set entity geometry from geometrical positions
//         * @{
//         */
//        /*!
//         * @brief Sets a vertex coordinates of a GeoModelMeshEntity
//         * @param[in] entity_id the entity to edit
//         * @param[in] v the index of the vertex in the entity
//         * @param[in] point the coordinates to set
//         * @param[in] update if true, updates all the colocated vertices
//         * to the new coordinates (ie if edit a Corner coordinates, it will updates
//         * its Lines, Surfaces...)
//         */
//        void set_mesh_entity_vertex(
//            const gme_t& entity_id,
//            index_t v,
//            const vec3& point,
//            bool update ) ;
//
//        void set_mesh_entity_vertices(
//            const gme_t& entity_id,
//            const std::vector< vec3 >& points,
//            bool clear ) ;
//
//        /*!
//         * @brief Sets the coordinates of a given existing Corner
//         * @param[in] corner_id the index of the corner in the GeoModel
//         * @param[in] point the coordinates to set
//         */
//        void set_corner( index_t corner_id, const vec3& point ) ;
//        /*!
//         * @brief Sets the mesh of a given existing Line
//         * @param[in] line_id the index of the line in the GeoModel
//         * @param[in] vertices the coordinates to set
//         * @warning the vertices should be ordered from the first boundary
//         * corner to the second one
//         */
//        void set_line( index_t line_id, const std::vector< vec3 >& vertices ) ;
//        /*!
//         * @brief Sets the mesh of a given existing Surface
//         * @param[in] surface_id the index of the surface in the GeoModel
//         * @param[in] surface_vertices the coordinates to set
//         * @param[in] surface_facets the vertex indices of the facets
//         * corresponding to \p surface_vertices
//         * @param[in] surface_facet_ptr the index of each new facet start in \p surface_facets
//         */
//        void set_surface_geometry(
//            index_t surface_id,
//            const std::vector< vec3 >& surface_vertices,
//            const std::vector< index_t >& surface_facets,
//            const std::vector< index_t >& surface_facet_ptr ) ;
//        /*!
//         * @brief Sets the tetrahedral mesh of a given existing Region
//         * @param[in] region_id the index of the region in the GeoModel
//         * @param[in] points the coordinates to set
//         * @param[in] tetras the vertex indices of the cells (to read 4 by 4)
//         * corresponding to \p points
//         */
//        void set_region_geometry(
//            index_t region_id,
//            const std::vector< vec3 >& points,
//            const std::vector< index_t >& tetras ) ;
//
//        /*! @}
//         * \name Set entity geometry using global GeoModel vertices
//         * @{
//         */
//        void set_mesh_entity_vertex(
//            const gme_t& id,
//            index_t v,
//            index_t geomodel_vertex ) ;
//
//        void set_mesh_entity_vertices(
//            const gme_t& entity_id,
//            const std::vector< index_t >& geomodel_vertices,
//            bool clear ) ;
//
//        void set_corner( index_t corner_id, index_t unique_vertex ) ;
//
//        void set_line( index_t id, const std::vector< index_t >& unique_vertices ) ;
//
//        void set_surface_geometry(
//            index_t surface_id,
//            const std::vector< index_t >& surface_vertices,
//            const std::vector< index_t >& surface_facets,
//            const std::vector< index_t >& surface_facet_ptr ) ;
//
//        void set_surface_geometry(
//            index_t surface_id,
//            const std::vector< index_t >& corners,
//            const std::vector< index_t >& facet_ptr ) ;
//
//        void set_surface_geometry(
//            index_t surface_id,
//            const std::vector< index_t >& triangle_corners ) ;
//
//        void set_surface_geometry_with_adjacencies(
//            index_t surface_id,
//            const std::vector< index_t >& triangle_corners,
//            const std::vector< index_t >& adjacent_triangles ) ;
//
//        void set_surface_element_geometry(
//            index_t surface_id,
//            index_t facet_id,
//            const std::vector< index_t >& corners ) ;
//
//        void set_surface_element_adjacency(
//            index_t surface_id,
//            index_t facet_id,
//            const std::vector< index_t >& adjacents ) ;
//
//        void set_region_geometry(
//            index_t region_id,
//            const std::vector< index_t >& tet_corners ) ;
//
//        void set_region_element_geometry(
//            index_t region_id,
//            index_t cell_id,
//            const std::vector< index_t >& corners ) ;
//
//        /*! @}
//         * \name Create entity element
//         * @{
//         */
//
//        index_t create_mesh_entity_vertices(
//            const gme_t& entity_id,
//            index_t nb_vertices ) ;
//
//        index_t create_surface_facet(
//            index_t surface_id,
//            const std::vector< index_t >& vertex_indices ) ;
//
//        index_t create_region_cell(
//            index_t region_id,
//            GEO::MeshCellType type,
//            const std::vector< index_t >& vertex_indices ) ;
//
//        index_t create_region_cells(
//            index_t region_id,
//            GEO::MeshCellType type,
//            index_t nb_cells ) ;
//
//        /*! @}
//         * \name Delete mesh element entities
//         * @{
//         */
//
//        void delete_mesh_entity_mesh( const gme_t& E_id ) ;
//        void delete_mesh_entity_isolated_vertices( const gme_t& E_id ) ;
//        void delete_mesh_entity_vertices(
//            const gme_t& E_id,
//            const std::vector< bool >& to_delete ) ;
//        void delete_corner_vertex( index_t corner_id ) ;
//        void delete_line_edges(
//            index_t line_id,
//            const std::vector< bool >& to_delete,
//            bool remove_isolated_vertices ) ;
//        void delete_surface_facets(
//            index_t surface_id,
//            const std::vector< bool >& to_delete,
//            bool remove_isolated_vertices ) ;
//        void delete_region_cells(
//            index_t region_id,
//            const std::vector< bool >& to_delete,
//            bool remove_isolated_vertices ) ;
//
//        /*! @}
//         * \name Misc
//         * @{
//         */
//
//        void compute_surface_adjacencies(
//            index_t surface_id,
//            bool recompute_adjacency = true ) ;
//        void compute_region_adjacencies(
//            index_t region_id,
//            bool recompute_adjacency = true ) ;
//        void triangulate_surface(
//            const RINGMesh::Surface& surface_in,
//            index_t surface_out ) ;
//
//        gme_t find_or_create_corner( const vec3& point ) ;
//        gme_t find_or_create_corner( index_t geomodel_point_id ) ;
//        gme_t find_or_create_line( const std::vector< vec3 >& vertices ) ;
//        gme_t find_or_create_line(
//            const std::vector< index_t >& incident_surfaces,
//            const gme_t& first_corner,
//            const gme_t& second_corner ) ;
//
//        void recompute_geomodel_mesh() ;
//
//        /*!
//         * @}
//         * \name Model building functions
//         */
//
//        /*!
//         * @brief From the Surfaces of the GeoModel, build its Lines and Corners
//         */
//        bool build_lines_and_corners_from_surfaces() ;
//
//        /*!
//         * @brief Build the regions of the GeoModel from the Surfaces
//         * @pre Function build_lines_and_corners_from_surfaces must have been called before
//         */
//        bool build_brep_regions_from_surfaces() ;
//
//        /*
//         * @brief From a GeoModel in which only Surface are defined, create corners, contacts
//         * and regions depending on the building flags
//         * @note Valdity is not checked
//         */
//        void build_geomodel_from_surfaces() ;
//
//        /*!
//         * @brief Finish up geomodel building and complete missing information.
//         */
//        void end_geomodel() ;
//
//    protected:
//        void set_surface_facet_adjacencies(
//            index_t surface_id,
//            const std::vector< index_t >& facets_ids,
//            const std::vector< index_t >& edges_ids,
//            const std::vector< index_t >& adjacent_triangles ) ;
//
//    protected:
//        /*! Options to toggle the building of entities from the available entities */
//        GeoModelBuildingFlags options_ ;
//
//        /*! Internal information */
//        std::vector< GeoModelRegionFromSurfaces* > regions_info_ ;
//
//    private:
//        void assign_surface_mesh_facets(
//            index_t surface_id,
//            const std::vector< index_t >& facets,
//            const std::vector< index_t >& facet_ptr ) ;
//
//        void assign_surface_triangle_mesh(
//            index_t surface_id,
//            const std::vector< index_t >& triangle_vertices ) ;
//        void update_facet_vertices_around_facet_vertex(
//            Surface& surface,
//            index_t facet,
//            index_t old_vertex,
//            index_t new_vertex ) ;
//        void update_facet_vertex(
//            Surface& surface,
//            const std::vector< index_t >& facets,
//            index_t old_vertex,
//            index_t new_vertex ) ;
//        void update_cell_vertex(
//            Region& region,
//            const std::vector< index_t >& cells,
//            index_t old_vertex,
//            index_t new_vertex ) ;
//        void assign_surface_triangle_mesh(
//            index_t surface_id,
//            const std::vector< index_t >& triangle_vertices,
//            const std::vector< index_t >& adjacent_triangles ) ;
//
//        void assign_region_tet_mesh(
//            index_t region_id,
//            const std::vector< index_t >& tet_vertices ) ;
//
//        void compute_universe() ;
//
//        void cut_surfaces_by_internal_lines() ;
//        void cut_regions_by_internal_surfaces() ;
//
//        void cut_surface_by_line( index_t surface_id, index_t line_id ) ;
//        void cut_region_by_surface( index_t region_id, index_t surface_id ) ;
//        void duplicate_surface_vertices_along_line(
//            index_t surface_id,
//            index_t line_id ) ;
//        void duplicate_region_vertices_along_surface(
//            index_t region_id,
//            index_t surface_id ) ;
//        index_t disconnect_surface_facets_along_line_edges(
//            index_t surface_id,
//            index_t line_id ) ;
//        index_t disconnect_region_cells_along_surface_facets(
//            index_t region_id,
//            index_t surface_id ) ;
//    } ;

}

namespace RINGMesh {
    class GeoModelBuilder2 ;
}

namespace RINGMesh {

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
        friend class GeoModelBuilderRemoval ;

    private:
        GeoModelMeshEntityAccess( GeoModelMeshEntity& gme )
            : gmme_( gme )
        {
        }

        std::string& modifiable_name()
        {
            return gmme_.name_ ;
        }

        index_t& modifiable_index()
        {
            return gmme_.id_.index ;
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

        void change_mesh_data_structure( const MeshType type )
        {
            if( EntityTypeManager::is_corner( gmme_.type_name() ) ) {
                Corner& corner = dynamic_cast< Corner& >( gmme_ ) ;
                Mesh0D* old_mesh = corner.mesh0d_ ;
                if( old_mesh->type_name() == type ) {
                    return ;
                }
                corner.update_mesh_storage_type( Mesh0D::create_mesh( type ) ) ;
                Mesh0DBuilder_var builder = Mesh0DBuilder::create_builder(
                    *corner.mesh0d_ ) ;
                builder->copy( *old_mesh, true ) ;
                delete old_mesh ;
            } else if( EntityTypeManager::is_line( gmme_.type_name() ) ) {
                Line& line = dynamic_cast< Line& >( gmme_ ) ;
                Mesh1D* old_mesh = line.mesh1d_ ;
                if( old_mesh->type_name() == type ) {
                    return ;
                }
                line.update_mesh_storage_type( Mesh1D::create_mesh( type ) ) ;
                Mesh1DBuilder_var builder = Mesh1DBuilder::create_builder(
                    *line.mesh1d_ ) ;
                builder->copy( *old_mesh, true ) ;
                delete old_mesh ;
            } else if( EntityTypeManager::is_surface( gmme_.type_name() ) ) {
                Surface& surface = dynamic_cast< Surface& >( gmme_ ) ;
                Mesh2D* old_mesh = surface.mesh2d_ ;
                if( old_mesh->type_name() == type ) {
                    return ;
                }
                surface.update_mesh_storage_type( Mesh2D::create_mesh( type ) ) ;
                Mesh2DBuilder_var builder = Mesh2DBuilder::create_builder(
                    *surface.mesh2d_ ) ;
                builder->copy( *old_mesh, true ) ;
                delete old_mesh ;
            } else if( EntityTypeManager::is_region( gmme_.type_name() ) ) {
                Region& region = dynamic_cast< Region& >( gmme_ ) ;
                Mesh3D* old_mesh = region.mesh3d_ ;
                if( old_mesh->type_name() == type ) {
                    return ;
                }
                region.update_mesh_storage_type( Mesh3D::create_mesh( type ) ) ;
                Mesh3DBuilder_var builder = Mesh3DBuilder::create_builder(
                    *region.mesh3d_ ) ;
                builder->copy( *old_mesh, true ) ;
                delete old_mesh ;
            } else {
                ringmesh_assert_not_reached ;
            }
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
        friend class GeoModelBuilderRemoval ;

    private:
        GeoModelGeologicalEntityAccess( GeoModelGeologicalEntity& gmge )
            : gmge_( gmge )
        {
        }

        std::string& modifiable_name()
        {
            return gmge_.name_ ;
        }

        index_t& modifiable_index()
        {
            return gmge_.id_.index ;
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

        /*!
         * @brief Add to the vector the entities which cannot exist if
         *        an entity in the set does not exist.
         * @return True if at least one entity was added.
         */
        bool get_dependent_entities( std::set< gme_t >& entities ) const ;

        template< typename T >
        gme_t create_mesh_entity( const MeshType type = "" )
        {
            const EntityType entity_type = T::type_name_static() ;
            index_t nb_entities( geomodel_.nb_mesh_entities( entity_type ) ) ;
            index_t new_id( nb_entities ) ;
            T* new_entity = GeoModelMeshEntityAccess::create_entity< T >( geomodel_,
                new_id, type ) ;
            geomodel_access_.modifiable_mesh_entities( entity_type ).push_back(
                new_entity ) ;
            return new_entity->gme_id() ;
        }

        gme_t create_geological_entity( const EntityType& type )
        {
            index_t index = find_or_create_geological_entity_type( type ) ;
            index_t id = static_cast< index_t >( geomodel_.nb_geological_entities(
                type ) ) ;
            GeoModelGeologicalEntity* E =
                GeoModelGeologicalEntityAccess::create_geological_entity( type,
                    geomodel_, id ) ;
            geomodel_access_.modifiable_geological_entities()[index].push_back( E ) ;
            return E->gme_id() ;
        }

        bool create_mesh_entities(
            const EntityType& type,
            index_t nb_additional_entities )
        {
            if( EntityTypeManager::is_corner( type ) ) {
                return create_mesh_entities< Corner >( nb_additional_entities ) ;
            } else if( EntityTypeManager::is_line( type ) ) {
                return create_mesh_entities< Line >( nb_additional_entities ) ;
            } else if( EntityTypeManager::is_surface( type ) ) {
                return create_mesh_entities< Surface >( nb_additional_entities ) ;
            } else if( EntityTypeManager::is_region( type ) ) {
                return create_mesh_entities< Region >( nb_additional_entities ) ;
            } else {
                ringmesh_assert_not_reached ;
                return false ;
            }
        }

        bool create_geological_entities( const EntityType& type, index_t nb ) ;

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

        void delete_mesh_entity( const EntityType& type, index_t index )
        {
            std::vector< GeoModelMeshEntity* >& store =
                geomodel_access_.modifiable_mesh_entities( type ) ;
            delete store[index] ;
            store[index] = nil ;
        }
        void delete_geological_entity( const EntityType& type, index_t index )
        {
            std::vector< GeoModelGeologicalEntity* >& store =
                geomodel_access_.modifiable_geological_entities( type ) ;
            delete store[index] ;
            store[index] = nil ;
        }

        /*!
         * @brief Finds or creates a corner at given coordinates.
         * @param[in] point Geometric location of the Corner
         * @return Index of the Corner
         */
        gme_t find_or_create_corner( const vec3& point ) ;
        gme_t find_or_create_corner( index_t geomodel_point_id ) ;
        gme_t find_or_create_line( const std::vector< vec3 >& vertices ) ;
        gme_t find_or_create_line(
            const std::vector< index_t >& incident_surfaces,
            const gme_t& first_corner,
            const gme_t& second_corner ) ;

        void compute_universe() ;

    private:
        GeoModelBuilderTopology( GeoModelBuilder2& builder, GeoModel& geomodel ) ;

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
        index_t create_geological_entity_type( const EntityType& type ) ;

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

        template< typename T >
        void complete_mesh_entity_connectivity() ;

        template< typename T >
        void copy_mesh_entity_topology( const GeoModel& from ) ;

        void copy_geological_entity_topology(
            const GeoModel& from,
            const EntityType& type ) ;

    private:
        GeoModelBuilder2& builder_ ;
        GeoModel& geomodel_ ;
        GeoModelAccess geomodel_access_ ;
    } ;

    class RINGMESH_API GeoModelBuilderRemoval {
    ringmesh_disable_copy( GeoModelBuilderRemoval ) ;
        friend class GeoModelBuilder2 ;

    public:

        /*!
         * @brief Remove a list of mesh entities of the geomodel
         * @details No check is done on the consistency of this removal
         *          The entities and all references to them are removed.
         *          All dependent entities should be in the set of entities to remove,
         *          with a prior call to get_dependent_entities function.
         *
         */
        void remove_mesh_entities( const std::set< gme_t >& entities ) ;

        void remove_geological_entities( const std::set< gme_t >& entities )
        {
            check_if_entities_are_not_meshed_entities( entities ) ;
            std::set< gme_t > mesh_entities ;
            for( std::set< gme_t >::const_iterator it = entities.begin();
                it != entities.end(); ++it ) {
                const GeoModelGeologicalEntity& cur_gmge =
                    geomodel_.geological_entity( *it ) ;
                for( index_t i = 0; i < cur_gmge.nb_children(); i++ ) {
                    mesh_entities.insert( cur_gmge.child( i ).gme_id() ) ;
                }
            }
            remove_mesh_entities( mesh_entities ) ;
        }

        /*!
         * Should be rewritten. Put as it was before someone removed it...
         */
        void remove_entities_and_dependencies(
            const std::set< gme_t >& entities_to_remove ) ;

    protected:
        GeoModelBuilderRemoval( GeoModelBuilder2& builder, GeoModel& geomodel ) ;

    private:
    private:
        // ---  High level functions ----------
        void initialize_for_removal(
            const std::set< gme_t >& mesh_entities_to_remove )
        {
            nb_mesh_entity_types_ = EntityTypeManager::nb_mesh_entity_types() ;
            nb_geological_entity_types_ = geomodel_.nb_geological_entity_types() ;
            nb_entity_types_ = nb_geological_entity_types_ + nb_mesh_entity_types_ ;
            nb_removed_mesh_entities_.resize( nb_mesh_entity_types_, 0 ) ;
            nb_removed_geological_entities_.resize( nb_geological_entity_types_,
                0 ) ;
            fill_entity_type_to_index_map() ;
            fill_nb_initial_entities() ;
            initialize_costly_storage() ;
            fill_nb_children_vector() ;

            check_if_entities_are_meshed( mesh_entities_to_remove ) ;
            fill_to_erase_vectors( mesh_entities_to_remove ) ;
            fill_removed_entities_and_mapping() ;
        }
        void do_delete_flagged_mesh_entities()
        {
            for( index_t i = 0; i < nb_mesh_entity_types_; ++i ) {
                for( index_t j = 0; j < nb_initial_mesh_entities_[i]; ++j ) {
                    if( mesh_entity_to_erase_[i][j] ) {
                        const EntityType& type_name = index_to_mesh_entity_type(
                            i ) ;
                        for( index_t p = 0;
                            p < geomodel_.mesh_entity( type_name, j ).nb_parents();
                            p++ ) {
                            gme_t parent =
                                geomodel_.mesh_entity( type_name, j ).parent_gme(
                                    p ) ;
                            nb_childs_[geological_entity_type_to_index( parent.type )][parent.index]-- ;
                        }

                        delete_mesh_entity( i, j ) ;
                    }
                }
                clear_null_mesh_entities( i ) ;
            }
        }
        void do_delete_flagged_geological_entities() ;

        void check_if_entities_are_meshed(
            const std::set< gme_t >& mesh_entities_to_remove )
        {
            for( std::set< gme_t >::const_iterator it =
                mesh_entities_to_remove.begin(); it != mesh_entities_to_remove.end();
                ++it ) {
                if( !RINGMesh::EntityTypeManager::is_mesh_entity_type( it->type ) ) {
                    throw RINGMeshException( "REMOVE",
                        "You try to remove a Geological Entity using mesh removal." ) ;
                }
            }
        }

        void check_if_entities_are_not_meshed_entities(
            const std::set< gme_t >& mesh_entities_to_remove )
        {
            for( std::set< gme_t >::const_iterator it =
                mesh_entities_to_remove.begin(); it != mesh_entities_to_remove.end();
                ++it ) {
                if( RINGMesh::EntityTypeManager::is_mesh_entity_type( it->type ) ) {
                    throw RINGMeshException( "REMOVE",
                        "You try to remove a Mesh Entity using geological removal" ) ;
                }
            }
        }
        void initialize_costly_storage()
        {
            mesh_entity_to_erase_.resize( nb_mesh_entity_types_ ) ;

            old_2_new_mesh_entity_.resize( nb_mesh_entity_types_ ) ;
            old_2_new_geological_entity_.resize( nb_geological_entity_types_ ) ;
            nb_childs_.resize( nb_geological_entity_types_ ) ;
            for( index_t i = 0; i < nb_mesh_entity_types_; ++i ) {
                index_t size = geomodel_.nb_mesh_entities(
                    index_to_mesh_entity_type( i ) ) ;
                mesh_entity_to_erase_[i].resize( size, false ) ;
                old_2_new_mesh_entity_[i].resize( size, 0 ) ;
            }

            for( index_t i = 0; i < nb_geological_entity_types_; ++i ) {
                index_t size = geomodel_.nb_geological_entities(
                    index_to_geological_entity_type( i ) ) ;
                old_2_new_geological_entity_[i].resize( size, 0 ) ;

                nb_childs_[i].resize( size, 0 ) ;

            }

        }
        void delete_mesh_entity( index_t type, index_t index ) ;

        void clear_null_mesh_entities( index_t type )
        {
            const EntityType& type_name = index_to_mesh_entity_type( type ) ;
            std::vector< GeoModelMeshEntity* >& store =
                geomodel_access_.modifiable_mesh_entities( type_name ) ;
            store.erase(
                std::remove( store.begin(), store.end(),
                    static_cast< GeoModelMeshEntity* >( nil ) ), store.end() ) ;

            // QC
            ringmesh_assert( geomodel_.nb_mesh_entities( type_name )
                == nb_initial_mesh_entities_[type] - nb_removed_mesh_entities_[type] ) ;
        }

        void clear_null_geological_entities( index_t type )
        {
            const EntityType& type_name = index_to_geological_entity_type( type ) ;
            std::vector< GeoModelGeologicalEntity* >& store =
                geomodel_access_.modifiable_geological_entities( type_name ) ;
            store.erase(
                std::remove( store.begin(), store.end(),
                    static_cast< GeoModelGeologicalEntity* >( nil ) ),
                store.end() ) ;

            // QC
            ringmesh_assert( geomodel_.nb_geological_entities( type_name )
                == nb_initial_geological_entities_[type] - nb_removed_geological_entities_[type] ) ;
        }
        void update_mesh_entity_connectivity()
        {
            for( index_t i = 0; i < nb_mesh_entity_types_; ++i ) {
                const EntityType& entity_type = index_to_mesh_entity_type( i ) ;

                for( index_t j = 0; j < geomodel_.nb_mesh_entities( entity_type );
                    ++j ) {
                    gme_t new_id( entity_type, j ) ;
                    GeoModelMeshEntity& ME = geomodel_access_.modifiable_mesh_entity(
                        new_id ) ;
                    update_mesh_entity_index( ME ) ;
                    ringmesh_assert( new_id == ME.gme_id() ) ;
                    update_mesh_entity_boundaries( ME ) ;
                    delete_invalid_boundaries( ME ) ;

                    update_mesh_entity_in_boundary( ME ) ;
                    delete_invalid_in_boundary( ME ) ;

                    if( ME.entity_type() == Region::type_name_static() ) {
                        Region& R = dynamic_cast< Region& >( ME ) ;
                        update_region_boundary_signs( R ) ;
                        delete_invalid_signs( R ) ;
                    }
                }
            }
        }

        void update_geological_entity_connectivity()
        {

            for( index_t i = 0; i < nb_geological_entity_types_; ++i ) {
                const EntityType& entity_type = index_to_geological_entity_type(
                    i ) ;
                for( index_t j = 0;
                    j < geomodel_.nb_geological_entities( entity_type ); ++j ) {
                    gme_t new_id( entity_type, j ) ;
                    GeoModelGeologicalEntity& GE =
                        geomodel_access_.modifiable_geological_entity( new_id ) ;
                    update_geological_entity_index( GE ) ;
                    update_geological_entity_children( GE ) ;
                    delete_invalid_children( GE ) ;
                }
            }

            for( index_t i = 0; i < nb_mesh_entity_types_; ++i ) {
                const EntityType& entity_type = index_to_mesh_entity_type( i ) ;

                for( index_t j = 0; j < geomodel_.nb_mesh_entities( entity_type );
                    ++j ) {
                    gme_t new_id( entity_type, j ) ;
                    GeoModelMeshEntity& ME = geomodel_access_.modifiable_mesh_entity(
                        new_id ) ;
                    update_mesh_entity_parents( ME ) ;
                    delete_invalid_parents( ME ) ;
                }
            }
        }

        void update_universe()
        {
            Universe& U = geomodel_access_.modifiable_universe() ;
            update_universe_sided_boundaries( U ) ;
            delete_invalid_universe_sided_boundaries( U ) ;
        }

        //        void remove_dependencies()
        //        {
        //            std::set< gme_t > new_gmme_to_remove ;
        //            for( index_t me = 0;
        //                me < geomodel().nb_mesh_entities( starting_dependency_ ); me++ ) {
        //                const GeoModelMeshEntity& cur_gmme = geomodel().mesh_entity(
        //                    starting_dependency_, me ) ;
        //                if( cur_gmme.in_boundary( 0 ).index() == NO_ID
        //                    && cur_gmme.nb_in_boundary() == 1 ) {
        //                    new_gmme_to_remove.insert( cur_gmme.gme_id() ) ;
        //                }
        //            }
        //            if( starting_dependency_ != Corner::type_name_static() ) {
        //                starting_dependency_ = EntityTypeManager::boundary_type(
        //                    starting_dependency_ ) ;
        //                remove_mesh_entities_with_dependencies( new_gmme_to_remove ) ;
        //            }
        //        }

        //------  Initialization -------
        void fill_removed_entities_and_mapping()
        {
            for( index_t i = 0; i < nb_mesh_entity_types_; ++i ) {
                for( index_t j = 0; j < nb_initial_mesh_entities_[i]; ++j ) {
                    if( mesh_entity_to_erase_[i][j] ) {
                        nb_removed_mesh_entities_[i]++ ;
                        old_2_new_mesh_entity_[i][j] = NO_ID ;
                    } else {
                        old_2_new_mesh_entity_[i][j] = j
                            - nb_removed_mesh_entities_[i] ;
                    }
                }
            }
        }
        void fill_to_erase_vectors(
            const std::set< gme_t >& mesh_entities_to_remove )
        {
            for( std::set< gme_t >::const_iterator it =
                mesh_entities_to_remove.begin(); it != mesh_entities_to_remove.end();
                ++it ) {
                gme_t cur = *it ;

                index_t type_index = mesh_entity_type_to_index( cur.type ) ;
                mesh_entity_to_erase_[type_index][cur.index] = true ;
            }
        }

        void fill_nb_children_vector()
        {
            for( index_t i = 0; i < nb_childs_.size(); i++ ) {
                for( index_t j = 0; j < nb_childs_[i].size(); j++ ) {
                    nb_childs_[i][j] = geomodel_.geological_entity(
                        index_to_geological_entity_type( i ), j ).nb_children() ;
                }
            }
        }

        void fill_nb_initial_entities()
        {
            nb_initial_mesh_entities_.resize( nb_mesh_entity_types_, 0 ) ;
            for( index_t i = 0; i < nb_mesh_entity_types_; ++i ) {
                const EntityType& type = index_to_mesh_entity_type( i ) ;
                nb_initial_mesh_entities_[i] = geomodel_.nb_mesh_entities( type ) ;
            }

            nb_initial_geological_entities_.resize( nb_geological_entity_types_,
                0 ) ;
            for( index_t i = 0; i < nb_geological_entity_types_; ++i ) {
                const EntityType& type = index_to_geological_entity_type( i ) ;
                nb_initial_geological_entities_[i] =
                    geomodel_.nb_geological_entities( type ) ;
            }
        }
        void fill_entity_type_to_index_map()
        {
            const EntityTypeManager& manager = geomodel_.entity_type_manager() ;
            mesh_entity_types_.insert( mesh_entity_types_.end(),
                manager.mesh_entity_types().begin(),
                manager.mesh_entity_types().end() ) ;

            if( nb_geological_entity_types_ != 0 ) {
                geological_entity_types_.insert( geological_entity_types_.end(),
                    manager.geological_entity_types().begin(),
                    manager.geological_entity_types().end() ) ;
            }

        }

        // ---- Easier access to relationships between EntityTypes
        index_t mesh_entity_type_index( const GeoModelMeshEntity& E ) const
        {
            const EntityType& type = E.type_name() ;
            return mesh_entity_type_to_index( type ) ;
        }
        index_t geological_entity_type_index(
            const GeoModelGeologicalEntity& E ) const
        {
            const EntityType& type = E.type_name() ;
            return geological_entity_type_to_index( type ) ;
        }
        index_t children_type_index( const EntityType& type ) const
        {
            const EntityType& child_type = children_type( type ) ;
            return mesh_entity_type_to_index( child_type ) ;
        }
        const EntityType children_type( const EntityType& type ) const
        {
            const EntityTypeManager& family = geomodel_.entity_type_manager() ;
            return family.child_type( type ) ;
        }
        index_t boundary_type_index( const EntityType& type ) const
        {
            const EntityType& b_type = boundary_type( type ) ;
            if( !EntityTypeManager::is_defined_type( b_type ) ) {
                return NO_ID ;
            } else {
                return mesh_entity_type_to_index( b_type ) ;
            }
        }
        const EntityType& boundary_type( const EntityType& type ) const
        {
            const EntityTypeManager& family = geomodel_.entity_type_manager() ;
            return family.boundary_type( type ) ;
        }
        index_t in_boundary_type_index( const EntityType& type ) const
        {
            const EntityType& in_b_type = in_boundary_type( type ) ;
            if( !EntityTypeManager::is_defined_type( in_b_type ) ) {
                return NO_ID ;
            } else {
                return mesh_entity_type_to_index( in_b_type ) ;
            }
        }
        const EntityType& in_boundary_type( const EntityType& type ) const
        {
            const EntityTypeManager& family = geomodel_.entity_type_manager() ;
            return family.in_boundary_type( type ) ;
        }
        bool is_mesh_entity( index_t i ) const
        {
            return i < nb_mesh_entity_types_ ;
        }
        bool is_geological_entity( index_t i ) const
        {
            return !is_mesh_entity( i ) ;
        }

        // ----  Update connectivity functions  ------

        void flag_geological_entities_without_children()
        {
            for( index_t i = 0; i < nb_childs_.size(); i++ ) {
                for( index_t j = 0; j < nb_childs_[i].size(); j++ ) {
                    if( nb_childs_[i][j] == 0 ) {
                        nb_removed_geological_entities_[i]++ ;
                        old_2_new_geological_entity_[i][j] = NO_ID ;
                    } else {
                        old_2_new_geological_entity_[i][j] = j
                            - nb_removed_geological_entities_[i] ;
                    }

                }
            }
        }

        void set_entity_index( GeoModelEntity& E, index_t new_index_in_geomodel ) ;

        void update_mesh_entity_index( GeoModelMeshEntity& ME ) ;
        void update_geological_entity_index( GeoModelGeologicalEntity& GE ) ;
        void update_mesh_entity_boundaries( GeoModelMeshEntity& ME ) ;

        void set_boundary_side( Region& R, index_t boundary_index, bool new_side )
        {
            ringmesh_assert( boundary_index < R.nb_boundaries() ) ;
            GeoModelMeshEntityAccess region_access(
                geomodel_access_.modifiable_mesh_entity( R.gme_id() ) ) ;
            region_access.modifiable_sides()[boundary_index] = new_side ;
        }

        void update_region_boundary_signs( Region& R )
        {
            const EntityType& surface_type = boundary_type( R.entity_type() ) ;
            gme_t invalid_value( surface_type, NO_ID ) ;

            index_t offset = 0 ;
            for( index_t i = 0; i + offset < R.nb_boundaries(); ++i ) {
                if( R.boundary_gme( i ) == invalid_value ) {
                    offset++ ;
                } else {
                    bool new_side = R.side( i + offset ) ;
                    set_boundary_side( R, i, new_side ) ;
                }
            }
        }
        void update_mesh_entity_in_boundary( GeoModelMeshEntity& E ) ;
        void update_mesh_entity_parents( GeoModelMeshEntity& E ) ;
        void update_geological_entity_children( GeoModelGeologicalEntity& E ) ;
        void update_universe_sided_boundaries( Universe& U ) ;

        // --- Deletion of some values the GeoModel storage
        void remove_invalid_values(
            std::vector< gme_t >& vector,
            const gme_t& invalid_value )
        {
            std::vector< gme_t >::iterator new_end = std::remove( vector.begin(),
                vector.end(), invalid_value ) ;
            if( new_end == vector.begin() ) {
                // Clear instead of erase, because the behavior would be undefined.
                vector.clear() ;
            } else if( new_end < vector.end() ) {
                vector.erase( new_end, vector.end() ) ;
            }
        }
        void delete_invalid_children( GeoModelGeologicalEntity& E )
        {
            if( E.nb_children() == 0 ) {
                return ;
            } else {
                const EntityType& child_type = children_type( E.entity_type() ) ;
                gme_t invalid_child( child_type, NO_ID ) ;
                GeoModelGeologicalEntityAccess gmge_access( E ) ;
                remove_invalid_values( gmge_access.modifiable_children(),
                    invalid_child ) ;
            }
        }
        void delete_invalid_boundaries( GeoModelMeshEntity& E )
        {
            const EntityType& b_type = boundary_type( E.entity_type() ) ;
            gme_t invalid( b_type, NO_ID ) ;
            if( !EntityTypeManager::is_defined_type( b_type ) ) {
                return ;
            } else {
                GeoModelMeshEntityAccess gmme_access( E ) ;
                remove_invalid_values( gmme_access.modifiable_boundaries(),
                    invalid ) ;
            }
        }
        void delete_invalid_in_boundary( GeoModelMeshEntity& E )
        {
            const EntityType& in_b_type = in_boundary_type( E.entity_type() ) ;
            gme_t invalid( in_b_type, NO_ID ) ;
            if( !EntityTypeManager::is_defined_type( in_b_type ) ) {
                return ;
            } else {
                GeoModelMeshEntityAccess gmme_access( E ) ;
                remove_invalid_values( gmme_access.modifiable_in_boundaries(),
                    invalid ) ;
            }
        }
        void delete_invalid_parents( GeoModelMeshEntity& E ) ;
        void delete_invalid_signs( Region& R )
        {
            GeoModelMeshEntityAccess region_access(
                geomodel_access_.modifiable_mesh_entity( R.gme_id() ) ) ;
            region_access.modifiable_sides().resize( R.nb_boundaries() ) ;
        }
        void delete_invalid_universe_sided_boundaries( Universe& U )
        {
            const EntityType& b_type = Surface::type_name_static() ;
            gme_t invalid( b_type, NO_ID ) ;
            UniverseAccess universe_access( U ) ;
            remove_invalid_values( universe_access.modifiable_boundaries(),
                invalid ) ;
            universe_access.modifiable_sides().resize( U.nb_boundaries() ) ;
        }

        index_t mesh_entity_type_to_index( const EntityType& type ) const
        {
            return find( mesh_entity_types_, type ) ;
        }

        index_t geological_entity_type_to_index( const EntityType& type ) const
        {
            return find( geological_entity_types_, type ) ;
        }
        const EntityType& index_to_mesh_entity_type( index_t index ) const
        {
            return mesh_entity_types_.at( index ) ;
        }

        const EntityType& index_to_geological_entity_type( index_t index ) const
        {
            return geological_entity_types_.at( index ) ;
        }

    private:
        GeoModelBuilder2& builder_ ;
        GeoModel& geomodel_ ;
        GeoModelAccess geomodel_access_ ;

        index_t nb_entity_types_ ;
        index_t nb_geological_entity_types_ ;
        index_t nb_mesh_entity_types_ ;

        std::vector< index_t > nb_initial_mesh_entities_ ;
        std::vector< index_t > nb_initial_geological_entities_ ;

        std::vector< index_t > nb_removed_mesh_entities_ ;
        std::vector< index_t > nb_removed_geological_entities_ ;

        /*! For each type of entity, store a vector of where the
         * entities to remove from the geomodel are flagged with NO_ID. */
        std::vector< std::vector< bool > > mesh_entity_to_erase_ ;
        /*! Stores the mapping table between indices for each type of
         *  element before and after the removal of entities */
        std::vector< std::vector< index_t > > old_2_new_mesh_entity_ ;
        std::vector< std::vector< index_t > > nb_childs_ ;

        std::vector< std::vector< index_t > > old_2_new_geological_entity_ ;

        //std::map< EntityType, index_t > entity_type_to_index_ ;
        //std::map< index_t, EntityType > index_to_entity_type_ ;
//        std::vector< EntityType > all_entity_types_ ;

        std::vector< EntityType > mesh_entity_types_ ;
        std::vector< EntityType > geological_entity_types_ ;

        EntityType starting_dependency_ ;
    } ;

    class RINGMESH_API GeoModelBuilderGeometry {
    ringmesh_disable_copy( GeoModelBuilderGeometry ) ;
        friend class GeoModelBuilder2 ;

    public:
        void recompute_geomodel_mesh() ;
        /*!
         * @brief Transfer general mesh information from one mesh
         * data structure to another one
         * @param[in] id the GeoModelMeshEntity id to operate on
         * @param[in] type the new mesh data structure type
         */
        void change_mesh_data_structure( const gme_t& id, const MeshType type )
        {
            GeoModelMeshEntityAccess gmme_access(
                geomodel_access_.modifiable_mesh_entity( id ) ) ;
            gmme_access.change_mesh_data_structure( type ) ;
        }

        /*!
         * @brief Copy all entity meshes from the input geomodel
         * @pre The geomodel under construction has exaclty the same number of entities
         * than the input geomodel.
         */
        void copy_meshes( const GeoModel& from ) ;

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

        void set_surface_facet_adjacencies(
            index_t surface_id,
            const std::vector< index_t >& facets_ids,
            const std::vector< index_t >& edges_ids,
            const std::vector< index_t >& adjacent_triangles ) ;

        /*!
         * @brief Computes and sets the adjacencies between the facets
         * @details The adjacent facet is given for each vertex of each facet for the edge
         * starting at this vertex.
         * If there is no neighbor inside the same Surface adjacent is set to NO_ID
         *
         * @param[in] surface_id Index of the surface
         * @param[in] recompute_adjacency If true, recompute the existing adjacencies
         */
        void compute_surface_adjacencies(
            index_t surface_id,
            bool recompute_adjacency = true ) ;
        void compute_region_adjacencies(
            index_t region_id,
            bool recompute_adjacency = true ) ;

        void cut_surfaces_by_internal_lines() ;

        void cut_regions_by_internal_surfaces() ;

        /*!
         * @brief Cuts a Surface along a Line assuming that the edges of the Line are edges of the Surface
         * @pre Surface is not already cut. Line L does not cut the Surface S into 2 connected components.
         * @todo Add a test for this function.
         */
        void cut_surface_by_line( index_t surface_id, index_t line_id ) ;

        void cut_region_by_surface( index_t region_id, index_t surface_id ) ;

    protected:
        GeoModelBuilderGeometry( GeoModelBuilder2& builder, GeoModel& geomodel ) ;

    private:
        void copy_meshes( const GeoModel& from, const std::string& entity_type ) ;
        void copy_mesh( const GeoModel& from, const gme_t& mesh_entity ) ;

        void assign_surface_mesh_facets(
            index_t surface_id,
            const std::vector< index_t >& facets,
            const std::vector< index_t >& facet_ptr ) ;

        void assign_surface_triangle_mesh(
            index_t surface_id,
            const std::vector< index_t >& triangle_vertices ) ;

        void assign_surface_triangle_mesh(
            index_t surface_id,
            const std::vector< index_t >& triangle_vertices,
            const std::vector< index_t >& adjacent_triangles ) ;

        void assign_region_tet_mesh(
            index_t region_id,
            const std::vector< index_t >& tet_vertices ) ;

        /*!
         * @brief Duplicates the surface vertices along the fake boundary
         * (NO_ID adjacencies but shared vertices) and duplicate the vertices
         */
        void duplicate_surface_vertices_along_line(
            index_t surface_id,
            index_t line_id ) ;
        void duplicate_region_vertices_along_surface(
            index_t region_id,
            index_t surface_id ) ;
        /*
         * @brief Resets the adjacencies for all Surface facets adjacent to the Line
         * @return The number of disconnection done
         * @pre All the edges of the Line are edges of at least one facet of the Surface
         */
        index_t disconnect_surface_facets_along_line_edges(
            index_t surface_id,
            index_t line_id ) ;
        index_t disconnect_region_cells_along_surface_facets(
            index_t region_id,
            index_t surface_id ) ;

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

    private:
        GeoModelBuilder2& builder_ ;
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
        GeoModelBuilderInfo( GeoModelBuilder2& builder, GeoModel& geomodel ) ;

    private:
        GeoModelBuilder2& builder_ ;
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
        GeoModelBuilderGeology( GeoModelBuilder2& builder, GeoModel& geomodel ) ;

    private:
        GeoModelBuilder2& builder_ ;
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

    class RINGMESH_API GeoModelBuilderFromSurfaces {
    ringmesh_disable_copy( GeoModelBuilderFromSurfaces ) ;
        friend class GeoModelBuilder2 ;

    public:
        /*
         * @brief From a GeoModel in which only Surfaces are defined,
         * create Corners, Lines and Regions depending on the building flags
         * @note Valdity is not checked
         */
        void build() ;

        /*!
         * @brief From the Surfaces of the GeoModel, build its Lines and Corners
         */
        bool build_lines_and_corners_from_surfaces() ;

    private:
        GeoModelBuilderFromSurfaces(
            GeoModelBuilder2& builder,
            GeoModel& geomodel ) ;

        /*!
         * @brief Build the regions of the GeoModel from the Surfaces
         * @pre Function build_lines_and_corners_from_surfaces
         * must have been called before
         */
        bool build_brep_regions_from_surfaces() ;

    public:
        /*! Options to toggle the building of entities from the available entities */
        GeoModelBuildingFlags options_ ;

    private:
        GeoModelBuilder2& builder_ ;
        GeoModel& geomodel_ ;
        GeoModelAccess geomodel_access_ ;

        /*! Internal information */
        std::vector< GeoModelRegionFromSurfaces* > regions_info_ ;
    } ;

    class RINGMESH_API GeoModelBuilder2 {
    ringmesh_disable_copy( GeoModelBuilder2 ) ;

    public:
        GeoModelBuilder2( GeoModel& geomodel ) ;
        virtual ~GeoModelBuilder2()
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
    class RINGMESH_API GeoModelBuilderFile: public GeoModelBuilder2 {
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

#endif
