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

#include <geogram/basic/attributes.h>

#include <ringmesh/basic/geometry.h>
#include <ringmesh/geomodel/core/geomodel.h>
#include <ringmesh/geomodel/core/geomodel_api.h>

#include <ringmesh/geomodel/core/geomodel_mesh_entity.h>
#include <ringmesh/io/geomodel_builder_gocad.h>

#include <ringmesh/mesh/mesh_index.h>
//#include <ringmesh/mesh/point_set_mesh.h>
//#include <ringmesh/mesh/line_mesh.h>
#include <ringmesh/mesh/surface_mesh.h>
#include <ringmesh/mesh/volume_mesh.h>
/*!
 * @brief Implementation of the class to build GeoModel from input
 * Gocad TSolid .so file
 * @author Pierre Anquez
 */

namespace
{
    using namespace RINGMesh;

    gmme_id find_corner( const GeoModel3D& geomodel, const vec3& point )
    {
        for( const auto& corner : geomodel.corners() )
        {
            if( corner.vertex( 0 ) == point )
            {
                return corner.gmme();
            }
        }
        return gmme_id();
    }

    std::string read_name_with_spaces(
        index_t field_id, const GEO::LineInput& line )
    {
        std::ostringstream oss;
        do
        {
            oss << line.field( field_id++ );
        } while( field_id < line.nb_fields() );
        return oss.str();
    }

    vec3 read_vertex_coordinates(
        GEO::LineInput& in, index_t start_field, int z_sign )
    {
        vec3 vertex;
        vertex.x = in.field_as_double( start_field++ );
        vertex.y = in.field_as_double( start_field++ );
        vertex.z = z_sign * in.field_as_double( start_field );
        return vertex;
    }

    std::vector< double > read_vertex_attributes(
        GEO::LineInput& in, index_t start_field, index_t nb_attribute_fields )
    {
        std::vector< double > vertex( nb_attribute_fields );
        for( auto& cur_attribute : vertex )
        {
            ringmesh_assert( !in.field_matches( start_field, "CNXYZ" ) );
            ringmesh_assert( !in.field_matches( start_field, "XYZ" ) );
            cur_attribute = in.field_as_double( start_field++ );
        }
        return vertex;
    }

    std::vector< double > read_cell_attributes(
        GEO::LineInput& in, index_t start_field, index_t nb_attribute_fields )
    {
        std::vector< double > cell( nb_attribute_fields );
        for( auto& cur_attribute : cell )
        {
            cur_attribute = in.field_as_double( start_field++ );
        }
        return cell;
    }

    /*!
     * \name Building surface
     * @{
     */

    /*!
     * @brief Gets the index of an Interface from its name
     * @param[in] geomodel GeoModel to consider
     * @param[in] interface_name Name of the interface to find
     * @return Index of the interface in the geomodel, NO_ID if not found.
     */
    template < index_t DIMENSION >
    gmge_id find_geological_entity( const GeoModel< DIMENSION >& geomodel,
        const GeologicalEntityType& geol_entity_type,
        const std::string& geol_entity_name )
    {
        for( auto& geol_entity : geomodel.geol_entities( geol_entity_type ) )
        {
            if( geol_entity.name() == geol_entity_name )
            {
                return geol_entity.gmge();
            }
        }
        return gmge_id();
    }

    /*!
     * @brief Structure used to build Line by GeoModelBuilderGocad
     */
    struct Border
    {
        Border( index_t part, index_t corner, index_t p0, index_t p1 )
            : part_id_( part ), corner_id_( corner ), p0_( p0 ), p1_( p1 )
        {
        }

        // Id of the Surface owning this Border
        index_t part_id_;

        // Id of p0 in the GeoModel corner vector
        index_t corner_id_;

        // Ids of the starting corner and second vertex on the border in the
        // Surface
        // to which this Border belong
        index_t p0_;
        index_t p1_;
    };

    /*!
     * @brief Gets the coordinates of the point from gocad index
     * @param[in] geomodel GeoModel to consider
     * @param[in] vertex_map Map between Gocad and GeoModel vertex indices
     * @param[in] point_gocad_id Gocad index of the point to get
     * @return Coordinates of the point
     */
    vec3 get_point_from_gocad_id( const GeoModel3D& geomodel,
        const VertexMap& vertex_map,
        index_t point_gocad_id )
    {
        index_t point_local_id = vertex_map.local_id( point_gocad_id );
        index_t point_region = vertex_map.region( point_gocad_id );

        return geomodel.region( point_region ).vertex( point_local_id );
    }

    /*!
     * @brief Gets the point and the index in the points vector to
     * build the polygons for one read gocad vertex
     * @param[in] vertex_gocad_id Gocad index of the vertex
     * @param[in] geomodel GeoModel to consider
     * @param[in] load_storage Set of tools useful for loading a GeoModel
     * @param[in] gocad_vertices2cur_surf_points Map between vertices with
     * gocad indices and position of the corresponding point
     * in the points vector
     * @param[out] cur_surf_points Vector of unique point coordinates
     * belonging to the surface
     * @param[out] cur_surf_polygons Vector of each polygon corner indices in
     * the cur_surf_points vector to build polygons
     */
    void get_surface_point_and_polygon_from_gocad_index(
        index_t vertex_gocad_id,
        const GeoModel3D& geomodel,
        const TSolidLoadingStorage& load_storage,
        std::vector< index_t >& gocad_vertices2cur_surf_points,
        std::vector< vec3 >& cur_surf_points,
        std::vector< index_t >& cur_surf_polygons )
    {
        if( vertex_gocad_id >= gocad_vertices2cur_surf_points.size() )
        {
            gocad_vertices2cur_surf_points.resize( vertex_gocad_id + 1, NO_ID );
        }

        if( gocad_vertices2cur_surf_points[vertex_gocad_id] == NO_ID )
        {
            // First time this polygon corner is met in polygon_corners
            vec3 point = get_point_from_gocad_id(
                geomodel, load_storage.vertex_map_, vertex_gocad_id );
            auto index = static_cast< index_t >( cur_surf_points.size() );
            cur_surf_polygons.push_back( index );
            gocad_vertices2cur_surf_points[vertex_gocad_id] = index;
            cur_surf_points.push_back( point );
        }
        else
        {
            // If this polygon corner has already been met in polygon_corners
            cur_surf_polygons.push_back(
                gocad_vertices2cur_surf_points[vertex_gocad_id] );
        }
    }

    /*!
     * @brief Gets the points and the indices in the points vector to
     * build the polygons
     * @param[in] geomodel GeoModel to consider
     * @param[in] load_storage Set of tools useful for loading a GeoModel
     * @param[out] cur_surf_points Vector of unique point coordinates
     * belonging to the surface
     * @param[out] cur_surf_polygons Vector of each polygon corner indices in
     * the cur_surf_points vector to build polygons
     */
    void get_surface_points_and_polygons_from_gocad_indices(
        const GeoModel3D& geomodel,
        const TSolidLoadingStorage& load_storage,
        std::vector< vec3 >& cur_surf_points,
        std::vector< index_t >& cur_surf_polygons )
    {
        std::vector< index_t > gocad_vertices2cur_surf_points;
        for( auto corner_gocad_id :
            load_storage.cur_surf_polygon_corners_gocad_id_ )
        {
            get_surface_point_and_polygon_from_gocad_index( corner_gocad_id,
                geomodel, load_storage, gocad_vertices2cur_surf_points,
                cur_surf_points, cur_surf_polygons );
        }
    }

    /*!
     * Builds surface by setting the points and polygons of the surface
     * @param[in] geomodel GeoModel to consider
     * @param[in] load_storage Set of tools useful for loading a GeoModel
     */
    void build_surface( GeoModelBuilderGocad& builder,
        GeoModel3D& geomodel,
        TSolidLoadingStorage& load_storage )
    {
        std::vector< vec3 > cur_surf_points;
        std::vector< index_t > cur_surf_polygons;
        get_surface_points_and_polygons_from_gocad_indices(
            geomodel, load_storage, cur_surf_points, cur_surf_polygons );
        builder.geometry.set_surface_geometry( load_storage.cur_surface_,
            cur_surf_points, cur_surf_polygons,
            load_storage.cur_surf_polygon_ptr_ );
        load_storage.cur_surf_polygon_corners_gocad_id_.clear();
        load_storage.cur_surf_polygon_ptr_.clear();
        load_storage.cur_surf_polygon_ptr_.push_back( 0 );
    }

    /*! @}
     * \name Linking surfaces and region boundaries
     * @{
     */

    /*!
     * @brief Builds a vector with the center of the cell
     * facets of a given region
     * @param[in] geomodel GeoModel to consider
     * @param[in] region_id Index of the region
     * @param[out] cell_facet_centers Vector of cell facet centers
     */
    void compute_region_cell_facet_centers( const GeoModel3D& geomodel,
        index_t region_id,
        std::vector< vec3 >& cell_facet_centers )
    {
        const Region3D& region = geomodel.region( region_id );
        const VolumeMesh3D& mesh = region.mesh();
        const index_t nb_cells = region.nb_mesh_elements();
        cell_facet_centers.reserve( 4 * nb_cells );
        for( auto c : range( nb_cells ) )
        {
            for( auto f : range( region.nb_cell_facets( c ) ) )
            {
                cell_facet_centers.push_back(
                    mesh.cell_facet_barycenter( CellLocalFacet( c, f ) ) );
            }
        }
    }

    /*!
     * @brief Computes the NNSearchs of the centers of cell facets for
     * each region
     * @param[in] geomodel GeoModel to consider
     * @return Pointers to the NNSearchs of regions
     */
    std::vector< std::unique_ptr< NNSearch3D > >
        compute_cell_facet_centers_region_nn_searchs(
            const GeoModel3D& geomodel )
    {
        std::vector< std::unique_ptr< NNSearch3D > > region_nn_searchs(
            geomodel.nb_regions() );
        for( auto r : range( geomodel.nb_regions() ) )
        {
            std::vector< vec3 > cell_facet_centers;
            compute_region_cell_facet_centers(
                geomodel, r, cell_facet_centers );
            region_nn_searchs[r].reset(
                new NNSearch3D( cell_facet_centers, true ) );
        }
        return region_nn_searchs;
    }

    /*!
     * @brief Tests if a surface is a boundary of a region.
     * @details If it is the case, add the surface to the boundaries of
     * the region and the region to the incident_entities of the surface
     * @param[in] surface Surface to test
     * @param[in] region_nn_search NNSearch of the region to test
     * @return a tuple containing:
     * - The number of surface sides bounding the region.
     * - Vector of colocated cell facet centers.
     */
    std::tuple< index_t, std::vector< index_t > >
        are_surface_sides_region_boundaries(
            const Surface3D& surface, const NNSearch3D& region_nn_search )
    {
        vec3 first_polygon_center = surface.mesh_element_barycenter( 0 );
        std::vector< index_t > colocated_cell_facet_centers;
        colocated_cell_facet_centers = region_nn_search.get_neighbors(
            first_polygon_center, surface.geomodel().epsilon() );
        return std::make_tuple(
            static_cast< index_t >( colocated_cell_facet_centers.size() ),
            colocated_cell_facet_centers );
    }

    /*!
     * @brief Determines which side of the surface is to be added in the
     * region boundaries
     * @param[in] geomodel GeoModel to consider
     * @param[in] region_id Index of the region
     * @param[in] surface_id Index of the surface
     * @param[in] cell_facet_center_id Index of the cell facet center
     * (i.e., cell_id * 4 + local_facet_id)
     * @return The side of the surface to add in the region boundaries,
     * i.e. true for the '+' side (normal size) and false for the
     * '-' side (other side)
     */
    bool determine_surface_side_to_add( const GeoModel3D& geomodel,
        index_t region_id,
        index_t surface_id,
        index_t cell_facet_center_id )
    {
        index_t local_facet_id = cell_facet_center_id % 4;
        index_t cell_id = ( cell_facet_center_id - local_facet_id ) / 4;
        vec3 cell_facet_normal =
            geomodel.region( region_id )
                .mesh()
                .cell_facet_normal( CellLocalFacet( cell_id, local_facet_id ) );
        vec3 first_polygon_normal =
            geomodel.surface( surface_id ).mesh().polygon_normal( 0 );
        return dot( first_polygon_normal, cell_facet_normal ) > 0;
    }

    /*!
     * @brief Both adds the surface in the boundaries of a region and
     * adds the region to the incident_entities of the surface
     * @param[in] region_id Index of the region
     * @param[in] surface_id Index of the surface
     * @param[in] surf_side Side of the surface bounding the region
     * @param[in] geomodel_builder Builder of the GeoModel to consider
     */
    void fill_region_and_surface_boundaries_links( index_t region_id,
        index_t surface_id,
        bool surf_side,
        GeoModelBuilderTSolid& geomodel_builder )
    {
        geomodel_builder.topology.add_region_surface_boundary_relation(
            region_id, surface_id, surf_side );
    }

    /*!
     * @brief Adds the both surface sides in the boundaries of a region
     * (internal boundary) and add twice the region to the incident_entities
     * of the surface
     * @param[in] region_id Index of the region
     * @param[in] surface_id Index of the surface
     * @param[in,out] geomodel_builder Builder of the GeoModel to consider
     */
    void add_both_surface_sides_to_region_boundaries( index_t region_id,
        index_t surface_id,
        GeoModelBuilderTSolid& geomodel_builder )

    {
        fill_region_and_surface_boundaries_links(
            region_id, surface_id, true, geomodel_builder );
        fill_region_and_surface_boundaries_links(
            region_id, surface_id, false, geomodel_builder );
    }

    /*!
     * @brief Adds one surface side in the boundaries of a region
     * and add the region to the incident_entities of the surface
     * @details The index of the cell facet center is used for the
     * determination of the side to add.
     * @param[in] region_id Index of the region
     * @param[in] surface_id Index of the surface
     * @param[in] cell_facet_center_id Index of the cell facet center
     * (i.e., cell_id * 4 + local_facet_id)
     * @param[in,out] geomodel_builder Builder of the GeoModel to consider
     */
    void add_one_surface_side_to_region_boundaries( index_t region_id,
        index_t surface_id,
        index_t cell_facet_center_id,
        GeoModelBuilderTSolid& geomodel_builder,
        const GeoModel3D& geomodel )
    {
        bool side = determine_surface_side_to_add(
            geomodel, region_id, surface_id, cell_facet_center_id );
        fill_region_and_surface_boundaries_links(
            region_id, surface_id, side, geomodel_builder );
    }

    /*!
     * @brief Adds the surface sides which bound the region to the
     * boundaries of the region (and add the region to incident entities
     * of the surface)
     * @param[in] surface_id Index of the surface
     * @param[in] region_id Index of the region
     * @param[in] colocated_cell_facet_centers Vector of colocated cell
     * facet centers
     * @param[in,out] geomodel_builder Builder of the GeoModel to consider
     */
    void add_surface_sides_to_region_boundaries( index_t surface_id,
        index_t region_id,
        const std::vector< index_t >& colocated_cell_facet_centers,
        const GeoModel3D& geomodel,
        GeoModelBuilderTSolid& geomodel_builder )
    {
        switch( colocated_cell_facet_centers.size() )
        {
        case 1:
            add_one_surface_side_to_region_boundaries( region_id, surface_id,
                colocated_cell_facet_centers[0], geomodel_builder, geomodel );
            break;
        case 2:
            add_both_surface_sides_to_region_boundaries(
                region_id, surface_id, geomodel_builder );
            break;
        default:
            ringmesh_assert_not_reached;
        }
    }

    /*!
     * @brief Sets the given surface as regions boundaries
     * @details Based on NNSearch, retrieves the regions bounded by the
     * given surface. One side or the both sides of the surface
     * could bound geomodel regions.
     * @param[in] surface_id Index of the surface
     * @param[in] region_nn_searchs Vector of NNSearchs of the geomodel regions
     * @param[in,out] geomodel_builder Builder of the GeoModel to consider
     */
    void add_surface_to_region_boundaries( index_t surface_id,
        const std::vector< std::unique_ptr< NNSearch3D > >& region_nn_searchs,
        const GeoModel3D& geomodel,
        GeoModelBuilderTSolid& geomodel_builder )
    {
        index_t cur_region{ 0 };
        index_t nb_added_surf_sides{ 0 };
        // Maximum 2 regions could be bounded by a single surface
        while( cur_region < geomodel.nb_regions() && nb_added_surf_sides < 2 )
        {
            index_t nb_surf_sides_are_boundary = NO_ID;
            std::vector< index_t > colocated_cell_facet_centers;
            std::tie(
                nb_surf_sides_are_boundary, colocated_cell_facet_centers ) =
                are_surface_sides_region_boundaries(
                    geomodel.surface( surface_id ),
                    *region_nn_searchs[cur_region] );
            if( nb_surf_sides_are_boundary > 0 )
            {
                add_surface_sides_to_region_boundaries( surface_id, cur_region,
                    colocated_cell_facet_centers, geomodel, geomodel_builder );
                nb_added_surf_sides += nb_surf_sides_are_boundary;
            }
            ++cur_region;
        }
        if( nb_added_surf_sides == 0 )
        {
            ringmesh_assert_not_reached;
        }
    }

    /*!
     * @brief Sets the boundaries of the GeoModel regions
     * @param[in,out] geomodel_builder Builder of the GeoModel to consider
     */
    void compute_boundaries_of_geomodel_regions(
        GeoModelBuilderTSolid& geomodel_builder, const GeoModel3D& geomodel )
    {
        std::vector< std::unique_ptr< NNSearch3D > > reg_nn_searchs =
            compute_cell_facet_centers_region_nn_searchs( geomodel );
        for( const auto& surface : geomodel.surfaces() )
        {
            add_surface_to_region_boundaries(
                surface.index(), reg_nn_searchs, geomodel, geomodel_builder );
        }
    }

    /*! @}
     * \name Surface internal borders determination
     * @{
     */

    /*!
     * @brief Finds if a surface polygon edge is an internal border
     * (i.e. shared by at least two surfaces)
     * @param[in] geomodel GeoModel to consider
     * @param[in] surface_id Index of the surface
     * @param[in] polygon Index of the polygon in the surface
     * @param[in] edge Index of the edge in the polygon
     * @param[in] surface_nns Unique pointers to the NNSearchs of surfaces
     * @param[in] surface_boxes Bounding Box of surfaces
     * @return True is the edge is found in at least another surface
     */
    bool is_edge_in_several_surfaces( const GeoModel3D& geomodel,
        index_t surface_id,
        index_t polygon,
        index_t edge,
        const std::vector< std::unique_ptr< NNSearch3D > >& surface_nns,
        const std::vector< Box3D >& surface_boxes )
    {
        const Surface3D& surface = geomodel.surface( surface_id );
        const SurfaceMesh3D& mesh = surface.mesh();
        const vec3 barycenter =
            mesh.polygon_edge_barycenter( PolygonLocalEdge( polygon, edge ) );
        std::vector< index_t > result;
        index_t tested_surf{ 0 };
        while( result.empty() && tested_surf < surface_nns.size() )
        {
            if( surface_boxes[tested_surf].contains( barycenter ) )
            {
                result = surface_nns[tested_surf]->get_neighbors(
                    barycenter, geomodel.epsilon() );
            }
            ++tested_surf;
        }
        return !result.empty();
    }

    /*!
     * @brief Gets the polygon edge barycenters of a given surface
     * @param[in] geomodel GeoModel to consider
     * @param[in] surface_id Index of the surface
     * @return Vector of all the border
     * edge barycenters of the surface
     */
    std::vector< vec3 > get_surface_border_edge_barycenters(
        const GeoModel3D& geomodel, index_t surface_id )
    {
        std::vector< vec3 > border_edge_barycenters;
        const Surface3D& surface = geomodel.surface( surface_id );
        const SurfaceMesh3D& mesh = surface.mesh();
        for( auto p : range( surface.nb_mesh_elements() ) )
        {
            for( auto e : range( surface.nb_mesh_element_vertices( p ) ) )
            {
                if( mesh.is_edge_on_border( PolygonLocalEdge( p, e ) ) )
                {
                    border_edge_barycenters.push_back(
                        mesh.polygon_edge_barycenter(
                            PolygonLocalEdge( p, e ) ) );
                }
            }
        }
        return border_edge_barycenters;
    }

    void assign_mesh_surface(
        GeoModelBuilderGocad& builder, MLLoadingStorage& load_storage )
    {
        std::vector< vec3 > vertices(
            load_storage.vertices_.begin() + load_storage.tface_vertex_ptr_,
            load_storage.vertices_.end() );
        for( index_t& id : load_storage.cur_surf_polygon_corners_gocad_id_ )
        {
            id -= load_storage.tface_vertex_ptr_;
        }
        builder.geometry.set_surface_geometry( load_storage.cur_surface_,
            vertices, load_storage.cur_surf_polygon_corners_gocad_id_,
            load_storage.cur_surf_polygon_ptr_ );
        load_storage.cur_surf_polygon_corners_gocad_id_.clear();
        load_storage.cur_surf_polygon_ptr_.clear();
        load_storage.cur_surf_polygon_ptr_.push_back( 0 );
        load_storage.cur_surface_++;
    }

    void assign_attributes_to_mesh( const Region3D& region,
        TSolidLoadingStorage& load_storage,
        const std::vector< std::vector< double > >& region_attributes )
    {
        index_t read_fields{ 0 };
        for( auto attrib_itr :
            range( load_storage.vertex_attribute_names_.size() ) )
        {
            std::string name = load_storage.vertex_attribute_names_[attrib_itr];
            if( region.vertex_attribute_manager().is_defined( name ) )
            {
                Logger::warn( "Transfer attribute", "The attribute ", name,
                    " already exists on the ", region.gmme() );
                continue;
            }
            GEO::Attribute< double > attr;
            index_t nb_dimensions =
                load_storage.vertex_attribute_dims_[attrib_itr];

            attr.create_vector_attribute( region.vertex_attribute_manager(),
                load_storage.vertex_attribute_names_[attrib_itr],
                nb_dimensions );
            // Does it resize all the past attributes to the size of the current
            // attribute?
            // Problematic, isn't it?
            region.vertex_attribute_manager().resize(
                static_cast< index_t >( region_attributes.size() )
                    * nb_dimensions
                + nb_dimensions );

            for( auto v_itr : range( region_attributes.size() ) )
            {
                for( auto attrib_dim_itr : range( nb_dimensions ) )
                {
                    attr[v_itr * nb_dimensions + attrib_dim_itr] =
                        region_attributes[v_itr][read_fields + attrib_dim_itr];
                }
            }
            read_fields += nb_dimensions;
        }
    }

    void assign_cell_attributes_to_mesh( const Region3D& region,
        TSolidLoadingStorage& load_storage,
        const std::vector< std::vector< double > >& region_attributes )
    {
        index_t read_fields{ 0 };
        for( auto attrib_itr :
            range( load_storage.cell_attribute_names_.size() ) )
        {
            std::string name = load_storage.cell_attribute_names_[attrib_itr];

            if( region.cell_attribute_manager().is_defined( name ) )
            {
                Logger::warn( "Transfer attribute", "The attribute ", name,
                    " already exists on the ", region.gmme() );
                continue;
            }

            GEO::Attribute< double > attr;
            index_t nb_dimensions =
                load_storage.cell_attribute_dims_[attrib_itr];
            attr.create_vector_attribute( region.cell_attribute_manager(),
                load_storage.cell_attribute_names_[attrib_itr], nb_dimensions );
            // Does it resize all the past attributes to the size of the current
            // attribute?
            // Problematic, isn't it?

            region.cell_attribute_manager().resize(
                static_cast< index_t >( region_attributes.size() )
                    * nb_dimensions
                + nb_dimensions );

            for( auto c_itr : range( region_attributes.size() ) )
            {
                for( auto attrib_dim_itr : range( nb_dimensions ) )
                {
                    attr[c_itr * nb_dimensions + attrib_dim_itr] =
                        region_attributes[c_itr][read_fields + attrib_dim_itr];
                }
            }
            read_fields += nb_dimensions;
        }
    }

    // Indices begin to 1 in Gocad
    index_t GOCAD_OFFSET{ 1 };

    class LoadZSign final : public GocadLineParser
    {
    public:
        LoadZSign( GeoModelBuilderGocad& gm_builder, GeoModel3D& geomodel )
            : GocadLineParser( gm_builder, geomodel )
        {
        }

    private:
        void execute(
            GEO::LineInput& line, GocadLoadingStorage& load_storage ) final
        {
            if( line.field_matches( 1, "Elevation" ) )
            {
                load_storage.z_sign_ = 1;
            }
            else if( line.field_matches( 1, "Depth" ) )
            {
                load_storage.z_sign_ = -1;
            }
            else
            {
                ringmesh_assert_not_reached;
            }
        }
    };

    class LoadTSurf final : public MLLineParser
    {
    public:
        LoadTSurf( GeoModelBuilderML& gm_builder, GeoModel3D& geomodel )
            : MLLineParser( gm_builder, geomodel )
        {
        }

    private:
        void execute(
            GEO::LineInput& line, MLLoadingStorage& load_storage ) final
        {
            ringmesh_unused( load_storage );
            std::string interface_name = read_name_with_spaces( 1, line );
            // Create an interface and set its name
            gmge_id interface_id = builder().geology.create_geological_entity(
                Interface3D::type_name_static() );
            builder().info.set_geological_entity_name(
                interface_id, interface_name );
        }
    };

    class LoadMLSurface final : public MLLineParser
    {
    public:
        LoadMLSurface( GeoModelBuilderML& gm_builder, GeoModel3D& geomodel )
            : MLLineParser( gm_builder, geomodel )
        {
        }

    private:
        void execute(
            GEO::LineInput& line, MLLoadingStorage& load_storage ) final
        {
            if( !load_storage.is_header_read_ )
            {
                /// Create Surface from the name of its parent Interface
                /// and its geological feature
                std::string geol = line.field( 2 );
                std::string interface_name = read_name_with_spaces( 3, line );
                create_surface( interface_name, geol );
            }
            else
            {
                if( !load_storage.vertices_.empty() )
                {
                    assign_mesh_surface( builder(), load_storage );
                    auto nb_vertices =
                        static_cast< index_t >( load_storage.vertices_.size() );
                    load_storage.tface_vertex_ptr_ = nb_vertices;
                }
            }
        }

        void create_surface(
            const std::string& interface_name, const std::string& type )
        {
            gmge_id parent = find_geological_entity(
                geomodel(), Interface3D::type_name_static(), interface_name );
            if( !interface_name.empty() )
            {
                ringmesh_assert( parent.is_defined() );
            }

            gmme_id children = builder().topology.create_mesh_entity(
                Surface3D::type_name_static() );
            builder().geology.add_parent_children_relation( parent, children );
            builder().geology.set_geological_entity_geol_feature( parent,
                GeoModelGeologicalEntity3D::determine_geological_type( type ) );
        }
    };

    class LoadLayer final : public MLLineParser
    {
    public:
        LoadLayer( GeoModelBuilderML& gm_builder, GeoModel3D& geomodel )
            : MLLineParser( gm_builder, geomodel )
        {
        }

    private:
        void execute(
            GEO::LineInput& line, MLLoadingStorage& load_storage ) final
        {
            ringmesh_unused( load_storage );
            /// Build the volumetric layers from their name and
            /// the ids of the regions they contain
            gmge_id layer_id = builder().geology.create_geological_entity(
                Layer3D::type_name_static() );
            builder().info.set_geological_entity_name(
                layer_id, line.field( 1 ) );
            bool end_layer = false;
            while( !end_layer )
            {
                line.get_line();
                line.get_fields();
                for( auto i : range( 5 ) )
                {
                    index_t region_id = line.field_as_uint( i );
                    if( region_id == 0 )
                    {
                        end_layer = true;
                        break;
                    }
                    // Remove Universe region
                    region_id -= geomodel().nb_surfaces() + 1;
                    // Correction because ids begin at 1 in the file
                    builder().geology.add_parent_children_relation(
                        layer_id, gmme_id( Region3D::type_name_static(),
                                      region_id - GOCAD_OFFSET ) );
                }
            }
        }
    };

    class MLEndSection final : public MLLineParser
    {
    public:
        MLEndSection( GeoModelBuilderML& gm_builder, GeoModel3D& geomodel )
            : MLLineParser( gm_builder, geomodel )
        {
        }

    private:
        void execute(
            GEO::LineInput& line, MLLoadingStorage& load_storage ) final
        {
            ringmesh_unused( line );
            if( !load_storage.is_header_read_ )
            {
                load_storage.is_header_read_ = true;
            }
            else
            {
                assign_mesh_surface( builder(), load_storage );
                load_storage.vertices_.clear();
                load_storage.tface_vertex_ptr_ = 0;
            }
        }
    };

    class LoadCorner final : public MLLineParser
    {
    public:
        LoadCorner( GeoModelBuilderML& gm_builder, GeoModel3D& geomodel )
            : MLLineParser( gm_builder, geomodel )
        {
        }

    private:
        void execute(
            GEO::LineInput& line, MLLoadingStorage& load_storage ) final
        {
            index_t v_id = line.field_as_uint( 1 ) - GOCAD_OFFSET;
            if( !find_corner( geomodel(), load_storage.vertices_[v_id] )
                     .is_defined() )
            {
                // Create the corner
                gmme_id corner_gme = builder().topology.create_mesh_entity(
                    Corner3D::type_name_static() );
                builder().geometry.set_corner(
                    corner_gme.index(), load_storage.vertices_[v_id] );
            }
        }
    };

    class LoadMLRegion final : public MLLineParser
    {
    public:
        LoadMLRegion( GeoModelBuilderML& gm_builder, GeoModel3D& geomodel )
            : MLLineParser( gm_builder, geomodel )
        {
        }

    private:
        void execute(
            GEO::LineInput& line, MLLoadingStorage& load_storage ) final
        {
            ringmesh_unused( load_storage );
            /// Read Region information and create them from their name,
            /// and the surfaces on their boundary
            std::string name = read_name_with_spaces( 2, line );

            std::vector< std::pair< index_t, bool > > region_boundaries =
                get_region_boundaries( line );

            // Create the entity if it is not the universe
            // Set the region name and boundaries
            if( name != "Universe" )
            {
                const gmme_id region_id = builder().topology.create_mesh_entity(
                    Region3D::type_name_static() );
                builder().info.set_mesh_entity_name( region_id, name );
                for( const std::pair< index_t, bool >& info :
                    region_boundaries )
                {
                    const index_t surface_id{ info.first };
                    builder().topology.add_region_surface_boundary_relation(
                        region_id.index(), surface_id, info.second );
                }
            }
        }

        std::vector< std::pair< index_t, bool > > get_region_boundaries(
            GEO::LineInput& line )
        {
            std::vector< std::pair< index_t, bool > > region_boundaries;
            bool end_region = false;
            while( !end_region )
            {
                line.get_line();
                line.get_fields();
                for( auto i : range( 5 ) )
                {
                    signed_index_t signed_id = line.field_as_int( i );
                    if( signed_id == 0 )
                    {
                        end_region = true;
                        break;
                    }
                    bool side = signed_id > 0;
                    auto id = static_cast< index_t >( std::abs( signed_id ) )
                              - GOCAD_OFFSET;
                    region_boundaries.emplace_back( id, side );
                }
            }
            return region_boundaries;
        }
    };

    class LoadRegion : public TSolidLineParser
    {
    public:
        LoadRegion( GeoModelBuilderTSolid& gm_builder, GeoModel3D& geomodel )
            : TSolidLineParser( gm_builder, geomodel )
        {
        }

    protected:
        /*!
         * @brief Creates an empty entity of type GeoModelEntity::REGION and
         * sets
         * its name from .so file
         * @param[in] region_name Name of the new region
         * @param[in] geomodel_builder Builder of the geomodel
         * @return The index of the initialized region
         */
        index_t initialize_region( const std::string& region_name,
            GeoModelBuilderGocad& geomodel_builder ) const
        {
            gmme_id cur_region = geomodel_builder.topology.create_mesh_entity(
                Region3D::type_name_static() );
            geomodel_builder.info.set_mesh_entity_name(
                cur_region, region_name );
            return cur_region.index();
        }
    };

    class LoadTSolidRegion final : public LoadRegion
    {
    public:
        LoadTSolidRegion(
            GeoModelBuilderTSolid& gm_builder, GeoModel3D& geomodel )
            : LoadRegion( gm_builder, geomodel )
        {
        }

    private:
        void execute(
            GEO::LineInput& line, TSolidLoadingStorage& load_storage ) final
        {
            if( !load_storage.vertices_.empty() )
            {
                builder().geometry.set_region_geometry(
                    load_storage.cur_region_, load_storage.vertices_,
                    load_storage.tetra_corners_ );
                assign_attributes_to_mesh(
                    geomodel().region( load_storage.cur_region_ ), load_storage,
                    load_storage.attributes_ );
                assign_cell_attributes_to_mesh(
                    geomodel().region( load_storage.cur_region_ ), load_storage,
                    load_storage.cell_attributes_ );
            }

            std::string region_name;
            if( line.nb_fields() == 1 )
            {
                region_name = "Unnamed";
            }
            else
            {
                region_name = line.field( 1 );
            }

            load_storage.attributes_.clear();
            load_storage.cell_attributes_.clear();
            load_storage.cur_region_ =
                initialize_region( region_name, builder() );
            load_storage.vertices_.clear();
            load_storage.tetra_corners_.clear();
        }
        void execute_light(
            GEO::LineInput& line, TSolidLoadingStorage& load_storage ) final
        {
            ringmesh_unused( line );
            ringmesh_unused( load_storage );
            // Nothing
        }
    };

    class LoadLightTSolidRegion final : public LoadRegion
    {
    public:
        LoadLightTSolidRegion(
            GeoModelBuilderTSolid& gm_builder, GeoModel3D& geomodel )
            : LoadRegion( gm_builder, geomodel )
        {
        }

    private:
        void execute(
            GEO::LineInput& line, TSolidLoadingStorage& load_storage ) final
        {
            ringmesh_unused( line );
            ringmesh_unused( load_storage );
            // Nothing
        }
        void execute_light(
            GEO::LineInput& line, TSolidLoadingStorage& load_storage ) final
        {
            std::string region_name = line.field( 2 );

            // Record new regions
            bool found;
            index_t region_id;
            std::tie( found, region_id ) =
                load_storage.vertex_map_.find_region_id_from_name(
                    region_name );
            if( !found )
            {
                region_id = initialize_region( region_name, builder() );
                load_storage.vertex_map_.add_new_region(
                    region_id, region_name );
            }

            load_storage.vertex_map_.record_vertex_with_its_region(
                load_storage.cur_gocad_vrtx_id1_, region_id );
            load_storage.vertex_map_.record_vertex_with_its_region(
                load_storage.cur_gocad_vrtx_id2_, region_id );
            load_storage.vertex_map_.record_vertex_with_its_region(
                load_storage.cur_gocad_vrtx_id3_, region_id );
            load_storage.vertex_map_.record_vertex_with_its_region(
                load_storage.cur_gocad_vrtx_id4_, region_id );

            }
    };

    class LoadVertex final : public GocadLineParser
    {
    public:
        LoadVertex( GeoModelBuilderGocad& gm_builder, GeoModel3D& geomodel )
            : GocadLineParser( gm_builder, geomodel )
        {
        }

    private:
        void execute(
            GEO::LineInput& line, GocadLoadingStorage& load_storage ) final
        {
            vec3 vertex =
                read_vertex_coordinates( line, 2, load_storage.z_sign_ );
            load_storage.vertices_.push_back( vertex );
            if( load_storage.nb_attribute_fields_ > 0 )
            {
                std::vector< double > attribute = read_vertex_attributes(
                    line, 5, load_storage.nb_attribute_fields_ );
                load_storage.attributes_.push_back( attribute );
            }
        }
    };

    class LoadMLAtom final : public MLLineParser
    {
    public:
        LoadMLAtom( GeoModelBuilderML& gm_builder, GeoModel3D& geomodel )
            : MLLineParser( gm_builder, geomodel )
        {
        }

    private:
        void execute(
            GEO::LineInput& line, MLLoadingStorage& load_storage ) final
        {
            index_t vertex_id = line.field_as_uint( 2 ) - GOCAD_OFFSET;
            const vec3& vertex = load_storage.vertices_[vertex_id];
            load_storage.vertices_.push_back( vertex );
        }
    };

    class LoadAttributeTSolid final : public TSolidLineParser
    {
    public:
        LoadAttributeTSolid(
            GeoModelBuilderTSolid& gm_builder, GeoModel3D& geomodel )
            : TSolidLineParser( gm_builder, geomodel )
        {
        }

    private:
        void execute(
            GEO::LineInput& line, TSolidLoadingStorage& load_storage ) final
        {
            load_storage.vertex_attribute_names_.reserve(
                line.nb_fields() - 1 );
            for( auto attrib_name_itr : range( 1, line.nb_fields() ) )
            {
                load_storage.vertex_attribute_names_.emplace_back(
                    line.field( attrib_name_itr ) );
            }
        }
        void execute_light(
            GEO::LineInput& line, TSolidLoadingStorage& load_storage ) final
        {
            execute( line, load_storage );
        }
    };

    class LoadCellAttributeTSolid final : public TSolidLineParser
    {
    public:
        LoadCellAttributeTSolid(
            GeoModelBuilderTSolid& gm_builder, GeoModel3D& geomodel )
            : TSolidLineParser( gm_builder, geomodel )
        {
        }

    private:
        void execute(
            GEO::LineInput& line, TSolidLoadingStorage& load_storage ) final
        {
            load_storage.cell_attribute_names_.reserve( line.nb_fields() - 1 );
            for( auto attrib_name_itr : range( 1, line.nb_fields() ) )
            {
                load_storage.cell_attribute_names_.emplace_back(
                    line.field( attrib_name_itr ) );
            }
        }
        void execute_light(
            GEO::LineInput& line, TSolidLoadingStorage& load_storage ) final
        {
            execute( line, load_storage );
        }
    };

    class LoadAttributeDimensionTSolid final : public TSolidLineParser
    {
    public:
        LoadAttributeDimensionTSolid(
            GeoModelBuilderTSolid& gm_builder, GeoModel3D& geomodel )
            : TSolidLineParser( gm_builder, geomodel )
        {
        }

    private:
        void execute(
            GEO::LineInput& line, TSolidLoadingStorage& load_storage ) final
        {
            load_storage.vertex_attribute_dims_.reserve( line.nb_fields() - 1 );
            for( auto attrib_size_itr : range( 1, line.nb_fields() ) )
            {
                load_storage.vertex_attribute_dims_.push_back(
                    line.field_as_uint( attrib_size_itr ) );
                load_storage.nb_attribute_fields_ +=
                    line.field_as_uint( attrib_size_itr );
            }
        }
        void execute_light(
            GEO::LineInput& line, TSolidLoadingStorage& load_storage ) final
        {
            execute( line, load_storage );
        }
    };

    class LoadCellAttributeDimensionTSolid final : public TSolidLineParser
    {
    public:
        LoadCellAttributeDimensionTSolid(
            GeoModelBuilderTSolid& gm_builder, GeoModel3D& geomodel )
            : TSolidLineParser( gm_builder, geomodel )
        {
        }

    private:
        void execute(
            GEO::LineInput& line, TSolidLoadingStorage& load_storage ) final
        {
            load_storage.cell_attribute_dims_.reserve( line.nb_fields() - 1 );
            for( auto attrib_size_itr : range( 1, line.nb_fields() ) )
            {
                load_storage.cell_attribute_dims_.push_back(
                    line.field_as_uint( attrib_size_itr ) );
                load_storage.nb_cell_attribute_fields_ +=
                    line.field_as_uint( attrib_size_itr );
            }
        }
        void execute_light(
            GEO::LineInput& line, TSolidLoadingStorage& load_storage ) final
        {
            execute( line, load_storage );
        }
    };

    class LoadTSolidVertex final : public TSolidLineParser
    {
    public:
        LoadTSolidVertex(
            GeoModelBuilderTSolid& gm_builder, GeoModel3D& geomodel )
            : TSolidLineParser( gm_builder, geomodel )
        {
        }

    private:
        void execute(
            GEO::LineInput& line, TSolidLoadingStorage& load_storage ) final
        {
            auto vertex_id =
                static_cast< index_t >( load_storage.vertices_.size() );
            load_storage.vertex_map_.add_vertex(
                vertex_id, load_storage.cur_region_ );
            GocadLineFactory::create( "VRTX", builder(), geomodel() )
                ->execute( line, load_storage );
        }
        void execute_light(
            GEO::LineInput& line, TSolidLoadingStorage& load_storage ) final
        {
            load_storage.vertex_map_.add_vertex(
                line.field_as_uint( 1 ) - GOCAD_OFFSET,
                load_storage.cur_region_ );
            GocadLineFactory::create( "VRTX", builder(), geomodel() )
                ->execute( line, load_storage );
        }
    };

    class LoadTSAtomic final : public TSolidLineParser
    {
    public:
        LoadTSAtomic( GeoModelBuilderTSolid& gm_builder, GeoModel3D& geomodel )
            : TSolidLineParser( gm_builder, geomodel )
        {
        }

    private:
        void execute(
            GEO::LineInput& line, TSolidLoadingStorage& load_storage ) final
        {
            read_and_add_atom_to_region_vertices( geomodel(), line,
                load_storage, load_storage.vertices_, load_storage.attributes_,
                load_storage.vertex_map_ );
        }
        void execute_light(
            GEO::LineInput& line, TSolidLoadingStorage& load_storage ) final
        {
            load_storage.lighttsolid_atom_map_.emplace(
                line.field_as_uint( 1 ) - GOCAD_OFFSET,
                line.field_as_uint( 2 ) - GOCAD_OFFSET );
            load_storage.vertex_map_.add_vertex(
                line.field_as_uint( 1 ) - GOCAD_OFFSET,
                load_storage.cur_region_ );
            load_storage.vertices_.emplace_back( vec3{} );
            if( load_storage.nb_attribute_fields_ > 0 )
            {
                std::vector< double > null_attrib(
                    load_storage.nb_attribute_fields_, 0 );
                load_storage.attributes_.push_back( null_attrib );
            }
        }

        /*!
         * @brief Reads atom information and adds it in the list
         * of region vertices only if it refers to a vertex of another region
         * @param[in] geomodel GeoModel
         * @param[in] line ACSII file reader
         * @param[in] load_storage Load informations
         * @param[in,out] region_vertices Vector of the coordinates of the
         * vertices of the region
         * @param[in,out] vertex_map Map between Gocad and GeoModel vertex
         * indices
         * @param[in,out] region_attributes Vector of the attributes of the
         * vertices of the region
         */
        void read_and_add_atom_to_region_vertices( const GeoModel3D& geomodel,
            const GEO::LineInput& line,
            const TSolidLoadingStorage& load_storage,
            std::vector< vec3 >& region_vertices,
            std::vector< std::vector< double > >& region_attributes,
            VertexMap& vertex_map )
        {
            const index_t referring_vertex =
                line.field_as_uint( 2 ) - GOCAD_OFFSET;
            const index_t referred_vertex_local_id =
                vertex_map.local_id( referring_vertex );
            const index_t referred_vertex_region_id =
                vertex_map.region( referring_vertex );
            if( referred_vertex_region_id < load_storage.cur_region_ )
            {
                // If the atom referred to a vertex of another region,
                // acting like for a vertex
                auto index = static_cast< index_t >( region_vertices.size() );
                vertex_map.add_vertex( index, load_storage.cur_region_ );

                region_vertices.push_back(
                    geomodel.region( referred_vertex_region_id )
                        .vertex( referred_vertex_local_id ) );

                std::vector< double > attribute_v;
                for( auto attrib_itr :
                    range( load_storage.vertex_attribute_names_.size() ) )
                {
                    std::string name =
                        load_storage.vertex_attribute_names_[attrib_itr];
                    index_t dim =
                        load_storage.vertex_attribute_dims_[attrib_itr];
                    GEO::Attribute< double > attr(
                        geomodel.region( referred_vertex_region_id )
                            .vertex_attribute_manager(),
                        name );
                    for( auto dim_itr : range( dim ) )
                    {
                        attribute_v.push_back(
                            attr[referred_vertex_local_id * dim + dim_itr] );
                    }
                }
                region_attributes.push_back( attribute_v );
            }
            else
            {
                // If the atom referred to an atom of the same region
                vertex_map.add_vertex(
                    referred_vertex_local_id, referred_vertex_region_id );
            }
        }
    };

    class LoadTetra final : public TSolidLineParser
    {
    public:
        LoadTetra( GeoModelBuilderTSolid& gm_builder, GeoModel3D& geomodel )
            : TSolidLineParser( gm_builder, geomodel )
        {
        }

    private:
        void execute(
            GEO::LineInput& line, TSolidLoadingStorage& load_storage ) final
        {
            std::vector< index_t > corners =
                read_tetraedra( line, load_storage.vertex_map_ );
            load_storage.tetra_corners_.insert(
                load_storage.tetra_corners_.end(), corners.begin(),
                corners.end() );
            if( load_storage.nb_cell_attribute_fields_ > 0 )
            {
                std::vector< double > attribute = read_cell_attributes(
                    line, 5, load_storage.nb_cell_attribute_fields_ );
                load_storage.cell_attributes_.push_back( attribute );
            }
        }
        void execute_light(
            GEO::LineInput& line, TSolidLoadingStorage& load_storage ) final
        {
            load_storage.cur_gocad_vrtx_id1_ =
                line.field_as_uint( 1 ) - GOCAD_OFFSET;
            load_storage.cur_gocad_vrtx_id2_ =
                line.field_as_uint( 2 ) - GOCAD_OFFSET;
            load_storage.cur_gocad_vrtx_id3_ =
                line.field_as_uint( 3 ) - GOCAD_OFFSET;
            load_storage.cur_gocad_vrtx_id4_ =
                line.field_as_uint( 4 ) - GOCAD_OFFSET;
            if( load_storage.nb_cell_attribute_fields_ > 0 ) {
                std::vector< double > attribute = read_cell_attributes(
                    line, 5, load_storage.nb_cell_attribute_fields_ );
                load_storage.cell_attributes_.push_back( attribute );
            }
        }

        /*!
         * @brief Reads the four vertices index of a tetrahedron
         * @details Reads gocad indices (from .so file) and transforms
         * them to vertex local (region) indices
         * @param[in] in ACSII file reader
         * @param[out] gocad_vertices2region_vertices Vector which maps the
         * indices
         * of vertices from Gocad .so file to the local (in region) indices of
         * vertices
         * @return Indices of the four vertices
         */
        std::vector< index_t > read_tetraedra(
            const GEO::LineInput& in, const VertexMap& vertex_map )
        {
            std::vector< index_t > corners_id( 4 );
            ringmesh_assert( corners_id.size() == 4 );
            corners_id[0] =
                vertex_map.local_id( in.field_as_uint( 1 ) - GOCAD_OFFSET );
            corners_id[1] =
                vertex_map.local_id( in.field_as_uint( 2 ) - GOCAD_OFFSET );
            corners_id[2] =
                vertex_map.local_id( in.field_as_uint( 3 ) - GOCAD_OFFSET );
            corners_id[3] =
                vertex_map.local_id( in.field_as_uint( 4 ) - GOCAD_OFFSET );
            return corners_id;
        }
    };

    class LoadName final : public GocadLineParser
    {
    public:
        LoadName( GeoModelBuilderGocad& gm_builder, GeoModel3D& geomodel )
            : GocadLineParser( gm_builder, geomodel )
        {
        }

    private:
        void execute(
            GEO::LineInput& line, GocadLoadingStorage& load_storage ) final
        {
            ringmesh_unused( load_storage );
            // Set to the GeoModel name if empty
            if( geomodel().name().empty() )
            {
                builder().info.set_geomodel_name(
                    read_name_with_spaces( 1, line ) );
            }
        }
    };

    class LoadLastRegion final : public TSolidLineParser
    {
    public:
        LoadLastRegion(
            GeoModelBuilderTSolid& gm_builder, GeoModel3D& geomodel )
            : TSolidLineParser( gm_builder, geomodel )
        {
        }

    private:
        void execute(
            GEO::LineInput& line, TSolidLoadingStorage& load_storage ) final
        {
            ringmesh_unused( line );
            if( !load_storage.vertices_.empty() )
            {
                builder().geometry.set_region_geometry(
                    load_storage.cur_region_, load_storage.vertices_,
                    load_storage.tetra_corners_ );
                assign_attributes_to_mesh(
                    geomodel().region( load_storage.cur_region_ ), load_storage,
                    load_storage.attributes_ );
                assign_cell_attributes_to_mesh(
                    geomodel().region( load_storage.cur_region_ ), load_storage,
                    load_storage.cell_attributes_ );

                load_storage.attributes_.clear();
                load_storage.cell_attributes_.clear();
                load_storage.vertices_.clear();
                load_storage.tetra_corners_.clear();
            }
        }
        void execute_light(
            GEO::LineInput& line, TSolidLoadingStorage& load_storage ) final
        {
            ringmesh_unused( line );
            get_light_tsolid_workflow_to_catch_up_with_tsolid_workflow(
                load_storage );
            // End of LightTSolid peculiar processing
        }

        std::vector< index_t > get_region_local_indices(
            const std::vector< VertexMap::RegionLocalVertex >&
                region_local_vertices )
        {
            std::vector< index_t > result;
            result.reserve( region_local_vertices.size() );
            for( auto& vertex : region_local_vertices )
            {
                result.push_back( vertex.local_id );
            }
            return result;
        }
        std::vector< vec3 > get_region_vertices(
            const std::vector< VertexMap::RegionLocalVertex >&
                region_local_vertices )
        {
            std::vector< vec3 > result;
            result.reserve( region_local_vertices.size() );
            for( auto& vertex : region_local_vertices )
            {
                result.push_back( vertex.tetra_vertex );
            }
            return result;
        }

        void get_light_tsolid_workflow_to_catch_up_with_tsolid_workflow(
            TSolidLoadingStorage& load_storage )
        {
            load_storage.vertex_map_.fill_with_lighttsolid_region_ids();

            std::vector< std::vector< index_t > > region_tetra_corners_local;
            std::vector< std::vector< vec3 > > region_vertices;
            std::vector< std::vector< std::vector< double > > >
                region_attributes;
            std::vector< std::vector< std::vector< double > > >
                region_cell_attributes;

            region_tetra_corners_local.resize(
                load_storage.vertex_map_.nb_regions() );
            region_vertices.resize( load_storage.vertex_map_.nb_regions() );
            region_attributes.resize( load_storage.vertex_map_.nb_regions() );
            region_cell_attributes.resize( load_storage.vertex_map_.nb_regions() );
            load_storage.vertex_map_.local_ids_.resize(
                load_storage.vertex_map_.nb_regions() );

            for( auto region_id : load_storage.vertex_map_.get_regions() )
            {
                ringmesh_assert( !load_storage.vertices_.empty() );
                /// Fill the region_vertices and local_ids
                std::vector< VertexMap::RegionLocalVertex >
                    region_local_indices =
                        load_storage.vertex_map_
                            .get_vertices_list_and_local_ids_from_gocad_ids(
                                load_storage.vertices_, region_id,
                                load_storage.lighttsolid_atom_map_ );
                region_vertices[region_id] =
                    get_region_vertices( region_local_indices );
                load_storage.vertex_map_.local_ids_[region_id] =
                    get_region_local_indices( region_local_indices );
                if( load_storage.nb_attribute_fields_ > 0 )
                {
                    load_storage.vertex_map_
                        .get_vertices_attributes_list_from_gocad_ids(
                            load_storage.attributes_, region_id,
                            load_storage.lighttsolid_atom_map_,
                            region_attributes[region_id] );
                    load_storage.vertex_map_
                        .get_cells_attributes_list_from_gocad_ids(
                            load_storage.cell_attributes_, region_id,
                            load_storage.lighttsolid_atom_map_,
                            region_cell_attributes[region_id] );
                }
            }

            load_storage.vertex_map_.fill_with_lighttsolid_local_ids();
            load_storage.vertex_map_.deal_with_same_region_atoms(
                load_storage.lighttsolid_atom_map_ );

            for( auto region_id : load_storage.vertex_map_.get_regions() )
            {
                /// Fill the region_tetra_corners
                load_storage.vertex_map_.get_tetra_corners_with_this_region_id(
                    region_id, region_tetra_corners_local[region_id] );

                builder().geometry.set_region_geometry( region_id,
                    region_vertices[region_id],
                    region_tetra_corners_local[region_id] );
                
                assign_attributes_to_mesh( geomodel().region( region_id ),
                    load_storage, region_attributes[region_id] );

                assign_cell_attributes_to_mesh(
                    geomodel().region( region_id ), load_storage,
                    region_cell_attributes[region_id] );
            }
            load_storage.tetra_corners_.clear();
            load_storage.attributes_.clear();
            load_storage.cell_attributes_.clear();
            load_storage.vertices_.clear();
            load_storage.lighttsolid_atom_map_.clear();
        }
    };

    class LoadInterface final : public TSolidLineParser
    {
    public:
        LoadInterface( GeoModelBuilderTSolid& gm_builder, GeoModel3D& geomodel )
            : TSolidLineParser( gm_builder, geomodel )
        {
        }

    private:
        void execute(
            GEO::LineInput& line, TSolidLoadingStorage& load_storage ) final
        {
            gmge_id created_interface =
                builder().geology.create_geological_entity(
                    Interface3D::type_name_static() );
            load_storage.cur_interface_ = created_interface.index();
            builder().info.set_geological_entity_name(
                created_interface, line.field( 1 ) );
        }
        void execute_light(
            GEO::LineInput& line, TSolidLoadingStorage& load_storage ) final
        {
            // LightTSolid Interface processing : same as TSolid processing
            execute( line, load_storage );
        }
    };

    class LoadSurface final : public TSolidLineParser
    {
    public:
        LoadSurface( GeoModelBuilderTSolid& gm_builder, GeoModel3D& geomodel )
            : TSolidLineParser( gm_builder, geomodel )
        {
        }

    private:
        void execute(
            GEO::LineInput& line, TSolidLoadingStorage& load_storage ) final
        {
            ringmesh_unused( line );
            // Compute the surface
            if( !load_storage.cur_surf_polygon_corners_gocad_id_.empty() )
            {
                build_surface( builder(), geomodel(), load_storage );
            }
            // Create a new surface
            gmme_id new_surface = builder().topology.create_mesh_entity(
                Surface3D::type_name_static() );
            load_storage.cur_surface_ = new_surface.index();
            builder().geology.add_parent_children_relation(
                gmge_id( Interface3D::type_name_static(),
                    load_storage.cur_interface_ ),
                new_surface );
        }
        void execute_light(
            GEO::LineInput& line, TSolidLoadingStorage& load_storage ) final
        {
            // LightTSolid Surface processing : same as TSolid processing
            execute( line, load_storage );
        }
    };

    class LoadLastSurface final : public TSolidLineParser
    {
    public:
        LoadLastSurface(
            GeoModelBuilderTSolid& gm_builder, GeoModel3D& geomodel )
            : TSolidLineParser( gm_builder, geomodel )
        {
        }

    private:
        void execute(
            GEO::LineInput& line, TSolidLoadingStorage& load_storage ) final
        {
            ringmesh_unused( line );
            // Compute the last surface
            if( !load_storage.cur_surf_polygon_corners_gocad_id_.empty() )
            {
                build_surface( builder(), geomodel(), load_storage );
            }
        }
        void execute_light(
            GEO::LineInput& line, TSolidLoadingStorage& load_storage ) final
        {
            // LightTSolid LastSurface processing : same as TSolid processing
            execute( line, load_storage );
        }
    };

    class LoadTriangle final : public GocadLineParser
    {
    public:
        LoadTriangle( GeoModelBuilderGocad& gm_builder, GeoModel3D& geomodel )
            : GocadLineParser( gm_builder, geomodel )
        {
        }

    private:
        void execute(
            GEO::LineInput& line, GocadLoadingStorage& load_storage ) final
        {
            read_triangle(
                line, load_storage.cur_surf_polygon_corners_gocad_id_ );
            load_storage.end_polygon();
        }

        /*!
         * @brief Reads the three vertices index of a triangle and adds
         * them to the polygon corners of the currently built surface
         * @details Reads gocad indices
         * @param[in] in ACSII file reader
         * @param[out] cur_surf_polygons Vector of each polygon corner indices
         * to build polygons
         */
        void read_triangle( const GEO::LineInput& in,
            std::vector< index_t >& cur_surf_polygons )
        {
            cur_surf_polygons.push_back( in.field_as_uint( 1 ) - GOCAD_OFFSET );
            cur_surf_polygons.push_back( in.field_as_uint( 2 ) - GOCAD_OFFSET );
            cur_surf_polygons.push_back( in.field_as_uint( 3 ) - GOCAD_OFFSET );
        }
    };

    void tsolid_import_factory_initialize()
    {
        TSolidLineFactory::register_creator< LoadTSolidRegion >( "TVOLUME" );
        TSolidLineFactory::register_creator< LoadLightTSolidRegion >( "#" );
        TSolidLineFactory::register_creator< LoadAttributeTSolid >(
            "PROPERTIES" );
        TSolidLineFactory::register_creator< LoadCellAttributeTSolid >(
            "TETRA_PROPERTIES" );
        TSolidLineFactory::register_creator< LoadAttributeDimensionTSolid >(
            "ESIZES" );
        TSolidLineFactory::register_creator< LoadCellAttributeDimensionTSolid >(
            "TETRA_ESIZES" );
        TSolidLineFactory::register_creator< LoadTSolidVertex >( "VRTX" );
        TSolidLineFactory::register_creator< LoadTSolidVertex >( "PVRTX" );
        TSolidLineFactory::register_creator< LoadTSAtomic >( "SHAREDVRTX" );
        TSolidLineFactory::register_creator< LoadTSAtomic >( "SHAREDPVRTX" );
        TSolidLineFactory::register_creator< LoadTSAtomic >( "ATOM" );
        TSolidLineFactory::register_creator< LoadTSAtomic >( "PATOM" );
        TSolidLineFactory::register_creator< LoadTetra >( "TETRA" );
        TSolidLineFactory::register_creator< LoadLastRegion >( "MODEL" );
        TSolidLineFactory::register_creator< LoadInterface >( "SURFACE" );
        TSolidLineFactory::register_creator< LoadSurface >( "TFACE" );
        TSolidLineFactory::register_creator< LoadLastSurface >( "END" );
    }

    void ml_import_factory_initialize()
    {
        MLLineFactory::register_creator< LoadTSurf >( "TSURF" );
        MLLineFactory::register_creator< LoadMLSurface >( "TFACE" );
        MLLineFactory::register_creator< LoadMLRegion >( "REGION" );
        MLLineFactory::register_creator< LoadLayer >( "LAYER" );
        MLLineFactory::register_creator< MLEndSection >( "END" );
        MLLineFactory::register_creator< LoadMLAtom >( "ATOM" );
        MLLineFactory::register_creator< LoadMLAtom >( "PATOM" );
    }

} // namespace
namespace RINGMesh
{
    Factory< std::string, MLLineParser, GeoModelBuilderML&, GeoModel3D& >
        ml_factory;
    Factory< std::string,
        TSolidLineParser,
        GeoModelBuilderTSolid&,
        GeoModel3D& >
        tsolid_factory;

    TSolidLineParser::TSolidLineParser(
        GeoModelBuilderTSolid& gm_builder, GeoModel3D& geomodel )
        : GocadBaseParser( gm_builder, geomodel )
    {
    }
    MLLineParser::MLLineParser(
        GeoModelBuilderML& gm_builder, GeoModel3D& geomodel )
        : GocadBaseParser( gm_builder, geomodel )
    {
    }

    void GeoModelBuilderGocad::read_file()
    {
        while( !file_line().eof() && file_line().get_line() )
        {
            file_line().get_fields();
            if( file_line().nb_fields() > 0 )
            {
                read_line();
            }
        }
    }

    GocadLoadingStorage::GocadLoadingStorage()
    {
        cur_surf_polygon_ptr_.push_back( 0 );
    }

    GeoModelBuilderTSolid::GeoModelBuilderTSolid(
        GeoModel3D& geomodel, std::string filename )
        : GeoModelBuilderGocad( geomodel, std::move( filename ) )
    {
        type_impl_[0].reset( new GeoModelBuilderTSolidImpl_TSolid(
            *this, geomodel, file_line(), tsolid_load_storage_ ) );
        type_impl_[1].reset( new GeoModelBuilderTSolidImpl_LightTSolid(
            *this, geomodel, file_line(), tsolid_load_storage_ ) );
    }

    void GeoModelBuilderTSolid::read_number_of_vertices()
    {
        GEO::LineInput line( filename() );
        tsolid_load_storage_.nb_vertices_ = 0;
        while( !line.eof() && line.get_line() )
        {
            line.get_fields();
            if( line.nb_fields() > 0 )
            {
                if( line.field_matches( 0, "VRTX" )
                    || line.field_matches( 0, "PVRTX" )
                    || line.field_matches( 0, "PATOM" )
                    || line.field_matches( 0, "ATOM" )
                    || line.field_matches( 0, "SHAREDPVRTX" )
                    || line.field_matches( 0, "SHAREDVRTX" ) )
                {
                    tsolid_load_storage_.nb_vertices_++;
                }
            }
        }
        tsolid_load_storage_.attributes_.reserve(
            tsolid_load_storage_.nb_vertices_ );
    }

    void GeoModelBuilderTSolid::load_file()
    {
        read_file();
        // Compute internal borders (by removing adjacencies on
        // triangle edges common to at least two surfaces)
        compute_surfaces_internal_borders();
        build_lines_and_corners_from_surfaces();
        compute_boundaries_of_geomodel_regions( *this, this->geomodel_ );
        geology.build_contacts();
    }

    void GeoModelBuilderTSolid::read_type()
    {
        GEO::LineInput line( filename() );
        while( !line.eof() && line.get_line() )
        {
            line.get_fields();
            if( line.nb_fields() > 0 )
            {
                if( line.field_matches( 0, "GOCAD" ) )
                {
                    if( line.field_matches( 1, "TSolid" ) )
                    {
                        file_type_ = TSolidType::TSOLID;
                    }
                    else
                    {
                        ringmesh_assert(
                            line.field_matches( 1, "LightTSolid" ) );
                        file_type_ = TSolidType::LIGHT_TSOLID;
                    }
                    break;
                }
            }
        }
    }

    void GeoModelBuilderTSolid::read_file()
    {
        read_type();
        read_number_of_vertices();
        GeoModelBuilderGocad::read_file();
    }

    void GeoModelBuilderTSolid::read_line()
    {
        tsolid_load_storage_.vertex_map_.reserve(
            tsolid_load_storage_.nb_vertices_ ); // Exactly
        tsolid_load_storage_.vertex_map_.reserve_nb_vertices(
            tsolid_load_storage_.nb_vertices_ ); // At least
        type_impl_[static_cast< index_t >( file_type_ )]->read_line();
    }

    void GeoModelBuilderTSolid::compute_surface_internal_borders(
        index_t surface_id,
        const std::vector< std::unique_ptr< NNSearch3D > >& surface_nns,
        const std::vector< Box3D >& surface_boxes )
    {
        const Surface3D& surface = geomodel_.surface( surface_id );
        const SurfaceMesh3D& mesh = surface.mesh();

        for( index_t p : range( surface.nb_mesh_elements() ) )
        {
            std::vector< index_t > adjacent_polygons_id( 3 );
            for( index_t e : range( 3 ) )
            {
                adjacent_polygons_id[e] =
                    surface.polygon_adjacent_index( PolygonLocalEdge( p, e ) );
                if( !mesh.is_edge_on_border( PolygonLocalEdge( p, e ) ) )
                {
                    bool internal_border =
                        is_edge_in_several_surfaces( geomodel_, surface_id, p,
                            e, surface_nns, surface_boxes );
                    if( internal_border )
                    {
                        adjacent_polygons_id[e] = NO_ID;
                    }
                }
            }
            geometry.set_surface_element_adjacency(
                surface_id, p, adjacent_polygons_id );
        }
    }

    void GeoModelBuilderTSolid::
        compute_polygon_edge_centers_nn_and_surface_boxes(
            std::vector< std::unique_ptr< NNSearch3D > >& surface_nns,
            std::vector< Box3D >& surface_boxes ) const
    {
        surface_nns.resize( geomodel_.nb_surfaces() );
        surface_boxes.resize( geomodel_.nb_surfaces() );

        for( const auto& surface : geomodel_.surfaces() )
        {
            for( index_t v : range( surface.nb_vertices() ) )
            {
                surface_boxes[surface.index()].add_point( surface.vertex( v ) );
            }
            std::vector< vec3 > border_edge_barycenters =
                get_surface_border_edge_barycenters(
                    geomodel_, surface.index() );
            surface_nns[surface.index()].reset(
                new NNSearch3D( border_edge_barycenters, true ) );
        }
    }

    void GeoModelBuilderTSolid::compute_surfaces_internal_borders()
    {
        std::vector< std::unique_ptr< NNSearch3D > > nn_searchs;
        std::vector< Box3D > boxes;
        compute_polygon_edge_centers_nn_and_surface_boxes( nn_searchs, boxes );
        for( const auto& surface : geomodel_.surfaces() )
        {
            compute_surface_internal_borders(
                surface.index(), nn_searchs, boxes );
        }
    }

    /*************************************************************************/

    MLLoadingStorage::MLLoadingStorage()
    {
        cur_surface_ = 0;
    }

    void GeoModelBuilderML::load_file()
    {
        read_file();
        geomodel_.mesh.vertices.test_and_initialize();
        build_lines_and_corners_from_surfaces();
        geology.build_contacts();
    }

    void GeoModelBuilderML::read_line()
    {
        std::string keyword = file_line().field( 0 );
        std::unique_ptr< MLLineParser > tsolid_parser(
            MLLineFactory::create( keyword, *this, geomodel_ ) );
        if( tsolid_parser )
        {
            tsolid_parser->execute( file_line(), ml_load_storage_ );
        }
        else
        {
            std::unique_ptr< GocadLineParser > gocad_parser =
                GocadLineFactory::create( keyword, *this, geomodel_ );
            if( gocad_parser )
            {
                gocad_parser->execute( file_line(), ml_load_storage_ );
            }
        }
    }

    void initialize_gocad_import_factories()
    {
        GocadLineFactory::register_creator< LoadZSign >( "ZPOSITIVE" );
        GocadLineFactory::register_creator< LoadVertex >( "VRTX" );
        GocadLineFactory::register_creator< LoadVertex >( "PVRTX" );
        GocadLineFactory::register_creator< LoadName >( "name:" );
        GocadLineFactory::register_creator< LoadTriangle >( "TRGL" );
        tsolid_import_factory_initialize();
        ml_import_factory_initialize();
    }
} // namespace RINGMesh
