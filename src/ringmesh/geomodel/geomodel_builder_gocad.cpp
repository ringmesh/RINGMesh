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

#include <ringmesh/geomodel/geomodel_builder_gocad.h>

#include <ringmesh/basic/geometry.h>
#include <ringmesh/geomodel/geomodel_api.h>

/*!
 * @brief Implementation of the class to build GeoModel from input
 * Gocad TSolid .so file
 * @author Pierre Anquez
 */

namespace {
    using namespace RINGMesh;

    gmme_id find_corner( const GeoModel< 3 >& geomodel, const vec3& point )
    {
        for( index_t i : range( geomodel.nb_corners() ) ) {
            if( geomodel.corner( i ).vertex( 0 ) == point ) {
                return gmme_id( Corner< 3 >::type_name_static(), i );
            }
        }
        return gmme_id();
    }

    std::string read_name_with_spaces( index_t field_id, const GEO::LineInput& line )
    {
        std::ostringstream oss;
        do {
            oss << line.field( field_id++ );
        } while( field_id < line.nb_fields() );
        return oss.str();
    }

    vec3 read_vertex_coordinates(
        GEO::LineInput& in,
        index_t start_field,
        int z_sign )
    {
        vec3 vertex;
        vertex.x = in.field_as_double( start_field++ );
        vertex.y = in.field_as_double( start_field++ );
        vertex.z = z_sign * in.field_as_double( start_field );
        return vertex;
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
    gmge_id find_interface(
        const GeoModel< 3 >& geomodel,
        const std::string& interface_name )
    {
        GeologicalEntityType type = Interface< 3 >::type_name_static();
        for( index_t i : range( geomodel.nb_geological_entities( type ) ) ) {
            if( geomodel.geological_entity( type, i ).name() == interface_name ) {
                return geomodel.geological_entity( type, i ).gmge();
            }
        }
        return gmge_id();
    }

    /*!
     * @brief Structure used to build Line by GeoModelBuilderGocad
     */
    struct Border {
        Border( index_t part, index_t corner, index_t p0, index_t p1 )
            : part_id_( part ), corner_id_( corner ), p0_( p0 ), p1_( p1 )
        {
        }

        // Id of the Surface owning this Border
        index_t part_id_;

        // Id of p0 in the GeoModel corner vector
        index_t corner_id_;

        // Ids of the starting corner and second vertex on the border in the Surface
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
    vec3 get_point_from_gocad_id(
        const GeoModel< 3 >& geomodel,
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
        const GeoModel< 3 >& geomodel,
        const TSolidLoadingStorage& load_storage,
        std::vector< index_t >& gocad_vertices2cur_surf_points,
        std::vector< vec3 >& cur_surf_points,
        std::vector< index_t >& cur_surf_polygons )
    {
        if( vertex_gocad_id >= gocad_vertices2cur_surf_points.size() ) {
            gocad_vertices2cur_surf_points.resize( vertex_gocad_id + 1, NO_ID );
        }

        if( gocad_vertices2cur_surf_points[vertex_gocad_id] == NO_ID ) {
            // First time this polygon corner is met in polygon_corners
            vec3 point = get_point_from_gocad_id( geomodel, load_storage.vertex_map_,
                vertex_gocad_id );
            index_t index = static_cast< index_t >( cur_surf_points.size() );
            cur_surf_polygons.push_back( index );
            gocad_vertices2cur_surf_points[vertex_gocad_id] = index;
            cur_surf_points.push_back( point );
        } else {
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
        const GeoModel< 3 >& geomodel,
        const TSolidLoadingStorage& load_storage,
        std::vector< vec3 >& cur_surf_points,
        std::vector< index_t >& cur_surf_polygons )
    {
        std::vector< index_t > gocad_vertices2cur_surf_points;
        for( index_t corner_gocad_id : load_storage.cur_surf_polygon_corners_gocad_id_ ) {
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
    void build_surface(
        GeoModelBuilderGocad& builder,
        GeoModel< 3 >& geomodel,
        TSolidLoadingStorage& load_storage )
    {
        std::vector< vec3 > cur_surf_points;
        std::vector< index_t > cur_surf_polygons;
        get_surface_points_and_polygons_from_gocad_indices( geomodel, load_storage,
            cur_surf_points, cur_surf_polygons );
        builder.geometry.set_surface_geometry( load_storage.cur_surface_,
            cur_surf_points, cur_surf_polygons, load_storage.cur_surf_polygon_ptr_ );
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
     * @param[in,out] cell_facet_centers Vector of cell facet centers
     */
    void compute_region_cell_facet_centers(
        const GeoModel< 3 >& geomodel,
        index_t region_id,
        std::vector< vec3 >& cell_facet_centers )
    {
        const Region< 3 >& region = geomodel.region( region_id );
        const VolumeMesh< 3 >& mesh = region.low_level_mesh_storage();
        const index_t nb_cells = region.nb_mesh_elements();
        cell_facet_centers.reserve( 4 * nb_cells );
        for( index_t c : range( nb_cells ) ) {
            for( index_t f : range( region.nb_cell_facets( c ) ) ) {
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
    std::vector< std::unique_ptr< NNSearch<3> > > compute_cell_facet_centers_region_nn_searchs(
        const GeoModel< 3 >& geomodel )
    {
        std::vector< std::unique_ptr< NNSearch<3> > > region_nn_searchs(
            geomodel.nb_regions() );
        for( index_t r : range( geomodel.nb_regions() ) ) {
            std::vector< vec3 > cell_facet_centers;
            compute_region_cell_facet_centers( geomodel, r, cell_facet_centers );
            region_nn_searchs[r].reset(
                new NNSearch< 3 >( cell_facet_centers, true ) );
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
    std::tuple< index_t, std::vector< index_t > > are_surface_sides_region_boundaries(
        const Surface< 3 >& surface,
        const NNSearch< 3 >& region_nn_search )
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
    bool determine_surface_side_to_add(
        const GeoModel< 3 >& geomodel,
        index_t region_id,
        index_t surface_id,
        index_t cell_facet_center_id )
    {
        index_t local_facet_id = cell_facet_center_id % 4;
        index_t cell_id = ( cell_facet_center_id - local_facet_id ) / 4;
        vec3 cell_facet_normal =
            geomodel.region( region_id ).low_level_mesh_storage().cell_facet_normal(
                CellLocalFacet( cell_id, local_facet_id ) );
        vec3 first_polygon_normal =
            geomodel.surface( surface_id ).low_level_mesh_storage().polygon_normal(
                0 );
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
    void fill_region_and_surface_boundaries_links(
        index_t region_id,
        index_t surface_id,
        bool surf_side,
        GeoModelBuilderTSolid& geomodel_builder )
    {
        geomodel_builder.topology.add_mesh_entity_boundary_relation(
            gmme_id( Region< 3 >::type_name_static(), region_id ),
            gmme_id( Surface< 3 >::type_name_static(), surface_id ), surf_side );
    }

    /*!
     * @brief Adds the both surface sides in the boundaries of a region
     * (internal boundary) and add twice the region to the incident_entities
     * of the surface
     * @param[in] region_id Index of the region
     * @param[in] surface_id Index of the surface
     * @param[in,out] geomodel_builder Builder of the GeoModel to consider
     */
    void add_both_surface_sides_to_region_boundaries(
        index_t region_id,
        index_t surface_id,
        GeoModelBuilderTSolid& geomodel_builder )

    {
        fill_region_and_surface_boundaries_links( region_id, surface_id, true,
            geomodel_builder );
        fill_region_and_surface_boundaries_links( region_id, surface_id, false,
            geomodel_builder );
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
    void add_one_surface_side_to_region_boundaries(
        index_t region_id,
        index_t surface_id,
        index_t cell_facet_center_id,
        GeoModelBuilderTSolid& geomodel_builder,
        const GeoModel< 3 >& geomodel )
    {
        bool side = determine_surface_side_to_add( geomodel, region_id, surface_id,
            cell_facet_center_id );
        fill_region_and_surface_boundaries_links( region_id, surface_id, side,
            geomodel_builder );
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
    void add_surface_sides_to_region_boundaries(
        index_t surface_id,
        index_t region_id,
        const std::vector< index_t >& colocated_cell_facet_centers,
        const GeoModel< 3 >& geomodel,
        GeoModelBuilderTSolid& geomodel_builder )
    {
        switch( colocated_cell_facet_centers.size() ) {
            case 1:
                add_one_surface_side_to_region_boundaries( region_id, surface_id,
                    colocated_cell_facet_centers[0], geomodel_builder, geomodel );
                break;
            case 2:
                add_both_surface_sides_to_region_boundaries( region_id, surface_id,
                    geomodel_builder );
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
    void add_surface_to_region_boundaries(
        index_t surface_id,
        const std::vector< std::unique_ptr< NNSearch< 3 > > >& region_nn_searchs,
        const GeoModel< 3 >& geomodel,
        GeoModelBuilderTSolid& geomodel_builder )
    {
        index_t cur_region = 0;
        index_t nb_added_surf_sides = 0;
        // Maximum 2 regions could be bounded by a single surface
        while( cur_region < geomodel.nb_regions() && nb_added_surf_sides < 2 ) {
            index_t nb_surf_sides_are_boundary = NO_ID;
            std::vector< index_t > colocated_cell_facet_centers;
            std::tie( nb_surf_sides_are_boundary, colocated_cell_facet_centers ) =
                are_surface_sides_region_boundaries( geomodel.surface( surface_id ),
                    *region_nn_searchs[cur_region] );
            if( nb_surf_sides_are_boundary > 0 ) {
                add_surface_sides_to_region_boundaries( surface_id, cur_region,
                    colocated_cell_facet_centers, geomodel, geomodel_builder );
                nb_added_surf_sides += nb_surf_sides_are_boundary;
            }
            ++cur_region;
        }
        if( nb_added_surf_sides == 0 ) {
            ringmesh_assert_not_reached;
        }
    }

    /*!
     * @brief Sets the boundaries of the GeoModel regions
     * @param[in,out] geomodel_builder Builder of the GeoModel to consider
     */
    void compute_boundaries_of_geomodel_regions(
        GeoModelBuilderTSolid& geomodel_builder,
        const GeoModel< 3 >& geomodel )
    {
        std::vector< std::unique_ptr< NNSearch< 3 > > > reg_nn_searchs =
            compute_cell_facet_centers_region_nn_searchs( geomodel );
        for( index_t s : range( geomodel.nb_surfaces() ) ) {
            add_surface_to_region_boundaries( s, reg_nn_searchs, geomodel,
                geomodel_builder );
        }
    }

    /*
     * @brief Adds the right surface sides in universe boundaries
     * @param[in] surf_side_minus Vector indicating if the '-' side of
     * surfaces are in the boundaries of geomodel regions
     * @param[in] surf_side_plus Vector indicating if the '+' side of
     * surfaces are in the boundaries of geomodel regions
     * @param[in,out] geomodel_builder Builder of the GeoModel to consider
     */
    void add_surfaces_to_universe_boundaries(
        const std::vector< bool >& surface_sides,
        index_t nb_surfaces,
        GeoModelBuilderTSolid& geomodel_builder )
    {
        for( index_t s : range( nb_surfaces ) ) {
            if( surface_sides[2 * s] && !surface_sides[2 * s + 1] ) {
                geomodel_builder.topology.add_universe_boundary( s, false );
            } else if( !surface_sides[2 * s] && surface_sides[2 * s + 1] ) {
                geomodel_builder.topology.add_universe_boundary( s, true );
            }
        }
    }

    /*
     * @brief Determines if each side of the surfaces are
     * in the boundaries of geomodel regions
     * @param[in] geomodel GeoModel to consider
     * @param[out] surf_side_minus Vector indicating if the '-' side of
     * surfaces are in the boundaries of geomodel regions
     * @return Vector indicating if the '+' side of
     * surfaces are in the boundaries of geomodel regions
     */
    std::vector< bool > determine_if_surface_sides_bound_regions(
        const GeoModel< 3 >& geomodel )
    {
        std::vector< bool > surface_sides( 2 * geomodel.nb_surfaces(), false );
        for( index_t r : range( geomodel.nb_regions() ) ) {
            for( index_t s : range( geomodel.region( r ).nb_boundaries() ) ) {
                if( geomodel.region( r ).side( s ) ) {
                    surface_sides[2 * geomodel.region( r ).boundary( s ).index() + 1] =
                        true;
                } else if( !geomodel.region( r ).side( s ) ) {
                    surface_sides[2 * geomodel.region( r ).boundary( s ).index()] =
                        true;
                } else {
                    ringmesh_assert_not_reached;
                }
            }
        }
        return surface_sides;
    }

    /*!
     * @brief Sets the boundaries of region Universe
     * @details A surface is set in the boundaries of region Universe if
     * only one of its sides belongs to the boundaries of other regions.
     * @param[in,out] geomodel_builder Builder of the GeoModel to consider
     */
    void compute_universe_boundaries(
        const GeoModel< 3 >& geomodel,
        GeoModelBuilderTSolid& geomodel_builder )
    {
        // The universe boundaries are the surfaces with only one side in all
        // the boundaries of the other regions
        std::vector< bool > surface_sides = determine_if_surface_sides_bound_regions(
            geomodel );
        add_surfaces_to_universe_boundaries( surface_sides, geomodel.nb_surfaces(),
            geomodel_builder );
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
    bool is_edge_in_several_surfaces(
        const GeoModel< 3 >& geomodel,
        index_t surface_id,
        index_t polygon,
        index_t edge,
        const std::vector< std::unique_ptr< NNSearch< 3 > > >& surface_nns,
        const std::vector< Box< 3 > >& surface_boxes )
    {
        const Surface< 3 >& surface = geomodel.surface( surface_id );
        const SurfaceMesh< 3 >& mesh = surface.low_level_mesh_storage();
        const vec3 barycenter = mesh.polygon_edge_barycenter(
            PolygonLocalEdge( polygon, edge ) );
        std::vector< index_t > result;
        index_t tested_surf = 0;
        while( result.empty() && tested_surf < surface_nns.size() ) {
            if( surface_boxes[tested_surf].contains( barycenter ) ) {
                result = surface_nns[tested_surf]->get_neighbors( barycenter,
                    geomodel.epsilon() );
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
        const GeoModel< 3 >& geomodel,
        index_t surface_id )
    {
        std::vector< vec3 > border_edge_barycenters;
        const Surface< 3 >& surface = geomodel.surface( surface_id );
        const SurfaceMesh< 3 >& mesh = surface.low_level_mesh_storage();
        for( index_t p : range( surface.nb_mesh_elements() ) ) {
            for( index_t e : range( surface.nb_mesh_element_vertices( p ) ) ) {
                if( mesh.is_edge_on_border( PolygonLocalEdge( p, e ) ) ) {
                    border_edge_barycenters.push_back(
                        mesh.polygon_edge_barycenter( PolygonLocalEdge( p, e ) ) );
                }
            }
        }
        return border_edge_barycenters;
    }

    void assign_mesh_surface(
        GeoModelBuilderGocad& builder,
        MLLoadingStorage& load_storage )
    {
        std::vector< vec3 > vertices(
            load_storage.vertices_.begin() + load_storage.tface_vertex_ptr_,
            load_storage.vertices_.end() );
        for( index_t& id : load_storage.cur_surf_polygon_corners_gocad_id_ ) {
            id -= load_storage.tface_vertex_ptr_;
        }
        builder.geometry.set_surface_geometry( load_storage.cur_surface_, vertices,
            load_storage.cur_surf_polygon_corners_gocad_id_,
            load_storage.cur_surf_polygon_ptr_ );
        load_storage.cur_surf_polygon_corners_gocad_id_.clear();
        load_storage.cur_surf_polygon_ptr_.clear();
        load_storage.cur_surf_polygon_ptr_.push_back( 0 );
        load_storage.cur_surface_++;
    }

    // Indices begin to 1 in Gocad
    index_t GOCAD_OFFSET = 1;

    class LoadZSign final : public GocadLineParser {
    private:
        void execute( GEO::LineInput& line, GocadLoadingStorage& load_storage ) final
        {
            if( line.field_matches( 1, "Elevation" ) ) {
                load_storage.z_sign_ = 1;
            } else if( line.field_matches( 1, "Depth" ) ) {
                load_storage.z_sign_ = -1;
            } else {
                ringmesh_assert_not_reached;
            }
        }
    };

    class LoadTSurf final : public MLLineParser {
    private:
        void execute( GEO::LineInput& line, MLLoadingStorage& load_storage ) final
        {
            ringmesh_unused( load_storage );
            std::string interface_name = read_name_with_spaces( 1, line );
            // Create an interface and set its name
            gmge_id interface_id = builder().geology.create_geological_entity(
                Interface< 3 >::type_name_static() );
            builder().info.set_geological_entity_name( interface_id,
                interface_name );
        }
    };

    class LoadMLSurface final : public MLLineParser {
    private:
        void execute( GEO::LineInput& line, MLLoadingStorage& load_storage ) final
        {
            if( !load_storage.is_header_read_ ) {
                /// Create Surface from the name of its parent Interface
                /// and its geological feature
                std::string geol = line.field( 2 );
                std::string interface_name = read_name_with_spaces( 3, line );
                create_surface( interface_name, geol );
            } else {
                if( !load_storage.vertices_.empty() ) {
                    assign_mesh_surface( builder(), load_storage );
                    index_t nb_vertices =
                        static_cast< index_t >( load_storage.vertices_.size() );
                    load_storage.tface_vertex_ptr_ = nb_vertices;
                }
            }
        }

        void create_surface(
            const std::string& interface_name,
            const std::string& type )
        {
            gmge_id parent = find_interface( geomodel(), interface_name );
            if( interface_name != "" ) {
                ringmesh_assert( parent.is_defined() );
            }

            gmme_id children = builder().topology.create_mesh_entity< Surface >();
            builder().geology.add_parent_children_relation( parent, children );
            builder().geology.set_geological_entity_geol_feature( parent,
                GeoModelGeologicalEntity< 3 >::determine_geological_type( type ) );
        }
    };

    class LoadLayer final : public MLLineParser {
    private:
        void execute( GEO::LineInput& line, MLLoadingStorage& load_storage ) final
        {
            ringmesh_unused( load_storage );
            /// Build the volumetric layers from their name and
            /// the ids of the regions they contain
            gmge_id layer_id = builder().geology.create_geological_entity(
                Layer< 3 >::type_name_static() );
            builder().info.set_geological_entity_name( layer_id, line.field( 1 ) );
            bool end_layer = false;
            while( !end_layer ) {
                line.get_line();
                line.get_fields();
                for( index_t i : range( 5 ) ) {
                    index_t region_id = line.field_as_uint( i );
                    if( region_id == 0 ) {
                        end_layer = true;
                        break;
                    } else {
                        // Remove Universe region
                        region_id -= geomodel().nb_surfaces() + 1;
                        // Correction because ids begin at 1 in the file
                        builder().geology.add_parent_children_relation( layer_id,
                            gmme_id( Region< 3 >::type_name_static(),
                                region_id - GOCAD_OFFSET ) );
                    }
                }
            }
        }
    };

    class MLEndSection final : public MLLineParser {
    private:
        void execute( GEO::LineInput& line, MLLoadingStorage& load_storage ) final
        {
            ringmesh_unused( line );
            if( !load_storage.is_header_read_ ) {
                load_storage.is_header_read_ = true;
            } else {
                assign_mesh_surface( builder(), load_storage );
                load_storage.vertices_.clear();
                load_storage.tface_vertex_ptr_ = 0;
            }
        }
    };

    class LoadCorner final : public MLLineParser {
    private:
        void execute( GEO::LineInput& line, MLLoadingStorage& load_storage ) final
        {
            index_t v_id = line.field_as_uint( 1 ) - GOCAD_OFFSET;
            if( !find_corner( geomodel(), load_storage.vertices_[v_id] ).is_defined() ) {
                // Create the corner
                gmme_id corner_gme =
                    builder().topology.create_mesh_entity< Corner >();
                builder().geometry.set_corner( corner_gme.index(),
                    load_storage.vertices_[v_id] );
            }
        }
    };

    class LoadMLRegion final : public MLLineParser {
    private:
        void execute( GEO::LineInput& line, MLLoadingStorage& load_storage ) final
        {
            ringmesh_unused( load_storage );
            /// Read Region information and create them from their name,
            /// and the surfaces on their boundary
            std::string name = read_name_with_spaces( 2, line );

            std::vector< std::pair< index_t, bool > > region_boundaries =
                get_region_boundaries( line );

            // Create the entity if it is not the universe
            // Set the region name and boundaries
            if( name
                != static_cast< std::string >( Universe< 3 >::universe_type_name() ) ) {
                gmme_id region_id =
                    builder().topology.create_mesh_entity< Region >();
                builder().info.set_mesh_entity_name( region_id, name );
                for( const std::pair< index_t, bool >& info : region_boundaries ) {
                    gmme_id surface_id( Surface< 3 >::type_name_static(),
                        info.first );
                    builder().topology.add_mesh_entity_boundary_relation( region_id,
                        surface_id, info.second );
                }
            } else {
                for( const std::pair< index_t, bool >& info : region_boundaries ) {
                    builder().topology.add_universe_boundary( info.first,
                        info.second );
                }
            }
        }

        std::vector< std::pair< index_t, bool > > get_region_boundaries(
            GEO::LineInput& line )
        {
            std::vector< std::pair< index_t, bool > > region_boundaries;
            bool end_region = false;
            while( !end_region ) {
                line.get_line();
                line.get_fields();
                for( index_t i : range( 5 ) ) {
                    signed_index_t signed_id = line.field_as_int( i );
                    if( signed_id == 0 ) {
                        end_region = true;
                        break;
                    }
                    bool side = signed_id > 0;
                    index_t id = static_cast< index_t >( std::abs( signed_id ) )
                        - GOCAD_OFFSET;
                    region_boundaries.push_back(
                        std::pair< index_t, bool >( id, side ) );
                }
            }
            return region_boundaries;
        }
    };

    class LoadRegion final: public TSolidLineParser {
    private:
        void execute( GEO::LineInput& line, TSolidLoadingStorage& load_storage ) final
        {
            if( !load_storage.vertices_.empty() ) {
                builder().geometry.set_region_geometry( load_storage.cur_region_,
                    load_storage.vertices_, load_storage.tetra_corners_ );
            }
            load_storage.cur_region_ = initialize_region( line.field( 1 ),
                builder() );
            load_storage.vertices_.clear();
            load_storage.tetra_corners_.clear();
        }

        /*!
         * @brief Creates an empty entity of type GeoModelEntity::REGION and sets
         * its name from .so file
         * @param[in] region_name Name of the new region
         * @param[in] geomodel_builder Builder of the geomodel
         * @return The index of the initialized region
         */
        index_t initialize_region(
            const std::string& region_name,
            GeoModelBuilderGocad& geomodel_builder )
        {
            gmme_id cur_region =
                geomodel_builder.topology.create_mesh_entity< Region >();
            geomodel_builder.info.set_mesh_entity_name( cur_region, region_name );
            return cur_region.index();
        }
    };

    class LoadVertex final : public GocadLineParser {
    private:
        void execute( GEO::LineInput& line, GocadLoadingStorage& load_storage ) final
        {
            vec3 vertex = read_vertex_coordinates( line, 2, load_storage.z_sign_ );
            load_storage.vertices_.push_back( vertex );
        }
    };

    class LoadMLAtom final : public MLLineParser {
    private:
        void execute( GEO::LineInput& line, MLLoadingStorage& load_storage ) final
        {
            index_t vertex_id = line.field_as_uint( 2 ) - GOCAD_OFFSET;
            const vec3& vertex = load_storage.vertices_[vertex_id];
            load_storage.vertices_.push_back( vertex );
        }
    };

    class LoadTSolidVertex final: public TSolidLineParser {
    private:
        void execute( GEO::LineInput& line, TSolidLoadingStorage& load_storage ) final
        {
            index_t vertex_id =
                static_cast< index_t >( load_storage.vertices_.size() );
            load_storage.vertex_map_.add_vertex( vertex_id,
                load_storage.cur_region_ );
            GocadLineParser::create( "VRTX", builder(), geomodel() )->execute( line,
                load_storage );
        }
    };

    class LoadTSAtomic final : public TSolidLineParser {
    private:
        void execute( GEO::LineInput& line, TSolidLoadingStorage& load_storage ) final
        {
            read_and_add_atom_to_region_vertices( geomodel(), line,
                load_storage.cur_region_, load_storage.vertices_,
                load_storage.vertex_map_ );
        }

        /*!
         * @brief Reads atom information and adds it in the list
         * of region vertices only if it refers to a vertex of another region
         * @param[in] line ACSII file reader
         * @param[in] region_id Index of the region
         * @param[in,out] region_vertices Vector of the coordinates of the
         * vertices of the region
         * @param[in] vertex_map Map between Gocad and GeoModel vertex indices
         */
        void read_and_add_atom_to_region_vertices(
            const GeoModel< 3 >& geomodel,
            GEO::LineInput& line,
            index_t region_id,
            std::vector< vec3 >& region_vertices,
            VertexMap& vertex_map )
        {
            const index_t referring_vertex = line.field_as_uint( 2 ) - GOCAD_OFFSET;
            const index_t referred_vertex_local_id = vertex_map.local_id(
                referring_vertex );
            const index_t referred_vertex_region_id = vertex_map.region(
                referring_vertex );
            if( referred_vertex_region_id < region_id ) {
                // If the atom referred to a vertex of another region,
                // acting like for a vertex
                index_t index = static_cast< index_t >( region_vertices.size() );
                vertex_map.add_vertex( index, region_id );
                region_vertices.push_back(
                    geomodel.region( referred_vertex_region_id ).vertex(
                        referred_vertex_local_id ) );
            } else {
                // If the atom referred to an atom of the same region
                vertex_map.add_vertex( referred_vertex_local_id,
                    referred_vertex_region_id );
            }
        }
    };

    class LoadTetra final : public TSolidLineParser {
    private:
        void execute( GEO::LineInput& line, TSolidLoadingStorage& load_storage ) final
        {
            std::vector< index_t > corners = read_tetraedra( line,
                load_storage.vertex_map_ );
            load_storage.tetra_corners_.insert( load_storage.tetra_corners_.end(),
                corners.begin(), corners.end() );
        }

        /*!
         * @brief Reads the four vertices index of a tetrahedron
         * @details Reads gocad indices (from .so file) and transforms
         * them to vertex local (region) indices
         * @param[in] in ACSII file reader
         * @param[out] gocad_vertices2region_vertices Vector which maps the indices
         * of vertices from Gocad .so file to the local (in region) indices of vertices
         * @return Indices of the four vertices
         */
        std::vector< index_t > read_tetraedra(
            GEO::LineInput& in,
            const VertexMap& vertex_map )
        {
            std::vector< index_t > corners_id( 4 );
            ringmesh_assert( corners_id.size() == 4 );
            corners_id[0] = vertex_map.local_id(
                in.field_as_uint( 1 ) - GOCAD_OFFSET );
            corners_id[1] = vertex_map.local_id(
                in.field_as_uint( 2 ) - GOCAD_OFFSET );
            corners_id[2] = vertex_map.local_id(
                in.field_as_uint( 3 ) - GOCAD_OFFSET );
            corners_id[3] = vertex_map.local_id(
                in.field_as_uint( 4 ) - GOCAD_OFFSET );
            return corners_id;
        }
    };

    class LoadName final : public GocadLineParser {
    private:
        void execute( GEO::LineInput& line, GocadLoadingStorage& load_storage ) final
        {
            ringmesh_unused( load_storage );
            // Set to the GeoModel name if empty
            if( geomodel().name().empty() ) {
                builder().info.set_geomodel_name( read_name_with_spaces( 1, line ) );
            }
        }
    };

    class LoadLastRegion final : public TSolidLineParser {
    private:
        void execute( GEO::LineInput& line, TSolidLoadingStorage& load_storage ) final
        {
            ringmesh_unused( line );
            if( !load_storage.vertices_.empty() ) {
                builder().geometry.set_region_geometry( load_storage.cur_region_,
                    load_storage.vertices_, load_storage.tetra_corners_ );
                load_storage.vertices_.clear();
                load_storage.tetra_corners_.clear();
            }
        }
    };

    class LoadInterface final : public TSolidLineParser {
    private:
        void execute( GEO::LineInput& line, TSolidLoadingStorage& load_storage ) final
        {
            gmge_id created_interface = builder().geology.create_geological_entity(
                Interface< 3 >::type_name_static() );
            load_storage.cur_interface_ = created_interface.index();
            builder().info.set_geological_entity_name( created_interface,
                line.field( 1 ) );
        }
    };

    class LoadSurface final : public TSolidLineParser {
    private:
        void execute( GEO::LineInput& line, TSolidLoadingStorage& load_storage ) final
        {
            ringmesh_unused( line );
            // Compute the surface
            if( !load_storage.cur_surf_polygon_corners_gocad_id_.empty() ) {
                build_surface( builder(), geomodel(), load_storage );
            }
            // Create a new surface
            gmme_id new_surface = builder().topology.create_mesh_entity< Surface >();
            load_storage.cur_surface_ = new_surface.index();
            builder().geology.add_parent_children_relation(
                gmge_id( Interface< 3 >::type_name_static(),
                    load_storage.cur_interface_ ), new_surface );
        }
    };

    class LoadLastSurface final : public TSolidLineParser {
    private:
        void execute( GEO::LineInput& line, TSolidLoadingStorage& load_storage ) final
        {
            ringmesh_unused( line );
            // Compute the last surface
            if( !load_storage.cur_surf_polygon_corners_gocad_id_.empty() ) {
                build_surface( builder(), geomodel(), load_storage );
            }
        }
    };

    class LoadTriangle final : public GocadLineParser {
    private:
        void execute( GEO::LineInput& line, GocadLoadingStorage& load_storage ) final
        {
            read_triangle( line, load_storage.cur_surf_polygon_corners_gocad_id_ );
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
        void read_triangle(
            GEO::LineInput& in,
            std::vector< index_t >& cur_surf_polygons )
        {
            cur_surf_polygons.push_back( in.field_as_uint( 1 ) - GOCAD_OFFSET );
            cur_surf_polygons.push_back( in.field_as_uint( 2 ) - GOCAD_OFFSET );
            cur_surf_polygons.push_back( in.field_as_uint( 3 ) - GOCAD_OFFSET );
        }
    };

    void tsolid_import_factory_initialize()
    {
        ringmesh_register_TSolidLineParser_creator( LoadRegion, "TVOLUME" );
        ringmesh_register_TSolidLineParser_creator( LoadTSolidVertex, "VRTX" );
        ringmesh_register_TSolidLineParser_creator( LoadTSolidVertex, "PVRTX" );
        ringmesh_register_TSolidLineParser_creator( LoadTSAtomic, "ATOM" );
        ringmesh_register_TSolidLineParser_creator( LoadTSAtomic, "PATOM" );
        ringmesh_register_TSolidLineParser_creator( LoadTetra, "TETRA" );
        ringmesh_register_TSolidLineParser_creator( LoadLastRegion, "MODEL" );
        ringmesh_register_TSolidLineParser_creator( LoadInterface, "SURFACE" );
        ringmesh_register_TSolidLineParser_creator( LoadSurface, "TFACE" );
        ringmesh_register_TSolidLineParser_creator( LoadLastSurface, "END" );
    }

    void ml_import_factory_initialize()
    {
        ringmesh_register_MLLineParser_creator( LoadTSurf, "TSURF" );
        ringmesh_register_MLLineParser_creator( LoadMLSurface, "TFACE" );
        ringmesh_register_MLLineParser_creator( LoadMLRegion, "REGION" );
        ringmesh_register_MLLineParser_creator( LoadLayer, "LAYER" );
        ringmesh_register_MLLineParser_creator( MLEndSection, "END" );
        ringmesh_register_MLLineParser_creator( LoadMLAtom, "ATOM" );
        ringmesh_register_MLLineParser_creator( LoadMLAtom, "PATOM" );
    }

}
// anonymous namespace

namespace RINGMesh {

    void GeoModelBuilderGocad::read_file()
    {
        while( !file_line_.eof() && file_line_.get_line() ) {
            file_line_.get_fields();
            if( file_line_.nb_fields() > 0 ) {
                read_line();
            }
        }
    }

    GocadLoadingStorage::GocadLoadingStorage()
    {
        cur_surf_polygon_ptr_.push_back( 0 );
    }

    std::unique_ptr< GocadLineParser > GocadLineParser::create(
        const std::string& keyword,
        GeoModelBuilderGocad& gm_builder,
        GeoModel< 3 >& geomodel )
    {
        std::unique_ptr< GocadLineParser > parser(
            GocadLineParserFactory::create_object( keyword ) );
        if( parser ) {
            parser->set_builder( gm_builder );
            parser->set_geomodel( geomodel );
        }
        return parser;
    }

    void GeoModelBuilderTSolid::load_file()
    {
        read_file();

        // Compute internal borders (by removing adjacencies on
        // triangle edges common to at least two surfaces)
        compute_surfaces_internal_borders();

        geomodel_.mesh.vertices.test_and_initialize();
        build_lines_and_corners_from_surfaces();

        compute_boundaries_of_geomodel_regions( *this, ( *this ).geomodel_ );

        compute_universe_boundaries( ( *this ).geomodel_, *this );

        geology.build_contacts();
    }

    void GeoModelBuilderTSolid::read_line()
    {
        std::string keyword = file_line_.field( 0 );
        std::unique_ptr< TSolidLineParser > tsolid_parser = TSolidLineParser::create(
            keyword, *this, geomodel_ );
        if( tsolid_parser ) {
            tsolid_parser->execute( file_line_, tsolid_load_storage_ );
        } else {
            std::unique_ptr< GocadLineParser > gocad_parser =
                GocadLineParser::create( keyword, *this, geomodel_ );
            if( gocad_parser ) {
                gocad_parser->execute( file_line_, tsolid_load_storage_ );
            }
        }
    }

    void GeoModelBuilderTSolid::compute_surface_internal_borders(
        index_t surface_id,
        const std::vector< std::unique_ptr< NNSearch< 3 > > >& surface_nns,
        const std::vector< Box< 3 > >& surface_boxes )
    {
        const Surface< 3 >& surface = geomodel_.surface( surface_id );
        const SurfaceMesh< 3 >& mesh = surface.low_level_mesh_storage();

        for( index_t p : range( surface.nb_mesh_elements() ) ) {
            std::vector< index_t > adjacent_polygons_id( 3 );
            for( index_t e : range( 3 ) ) {
                adjacent_polygons_id[e] = surface.polygon_adjacent_index(
                    PolygonLocalEdge( p, e ) );
                if( !mesh.is_edge_on_border( PolygonLocalEdge( p, e ) ) ) {
                    bool internal_border = is_edge_in_several_surfaces( geomodel_,
                        surface_id, p, e, surface_nns, surface_boxes );
                    if( internal_border ) {
                        adjacent_polygons_id[e] = NO_ID;
                    }
                }
            }
            geometry.set_surface_element_adjacency( surface_id, p,
                adjacent_polygons_id );
        }
    }

    void GeoModelBuilderTSolid::compute_polygon_edge_centers_nn_and_surface_boxes(
        std::vector< std::unique_ptr< NNSearch< 3 > > >& surface_nns,
        std::vector< Box< 3 > >& surface_boxes ) const
    {
        surface_nns.resize( geomodel_.nb_surfaces() );
        surface_boxes.resize( geomodel_.nb_surfaces() );

        for( index_t s : range( geomodel_.nb_surfaces() ) ) {
            const Surface< 3 >& surface = geomodel_.surface( s );
            for( index_t v : range( surface.nb_vertices() ) ) {
                surface_boxes[s].add_point( surface.vertex( v ) );
            }
            std::vector< vec3 > border_edge_barycenters =
                get_surface_border_edge_barycenters( geomodel_, s );
            surface_nns[s].reset(
                new NNSearch< 3 >( border_edge_barycenters, true ) );
        }
    }

    void GeoModelBuilderTSolid::compute_surfaces_internal_borders()
    {
        std::vector< std::unique_ptr< NNSearch< 3 > > > nn_searchs;
        std::vector< Box< 3 > > boxes;
        compute_polygon_edge_centers_nn_and_surface_boxes( nn_searchs, boxes );
        for( index_t s : range( geomodel_.nb_surfaces() ) ) {
            compute_surface_internal_borders( s, nn_searchs, boxes );
        }
    }

    std::unique_ptr< TSolidLineParser > TSolidLineParser::create(
        const std::string& keyword,
        GeoModelBuilderTSolid& gm_builder,
        GeoModel< 3 >& geomodel )
    {
        std::unique_ptr< TSolidLineParser > parser(
            TSolidLineParserFactory::create_object( keyword ) );
        if( parser ) {
            parser->set_builder( gm_builder );
            parser->set_geomodel( geomodel );
        }
        return parser;
    }

    /*************************************************************************/

    std::unique_ptr< MLLineParser > MLLineParser::create(
        const std::string& keyword,
        GeoModelBuilderML& gm_builder,
        GeoModel< 3 >& geomodel )
    {
        std::unique_ptr< MLLineParser > parser(
            MLLineParserFactory::create_object( keyword ) );
        if( parser ) {
            parser->set_builder( gm_builder );
            parser->set_geomodel( geomodel );
        }
        return parser;
    }

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
        std::string keyword = file_line_.field( 0 );
        std::unique_ptr< MLLineParser > tsolid_parser(
            MLLineParser::create( keyword, *this, geomodel_ ) );
        if( tsolid_parser ) {
            tsolid_parser->execute( file_line_, ml_load_storage_ );
        } else {
            std::unique_ptr< GocadLineParser > gocad_parser =
                GocadLineParser::create( keyword, *this, geomodel_ );
            if( gocad_parser ) {
                gocad_parser->execute( file_line_, ml_load_storage_ );
            }
        }
    }

    void initialize_gocad_import_factories()
    {
        ringmesh_register_GocadLineParser_creator( LoadZSign, "ZPOSITIVE" );
        ringmesh_register_GocadLineParser_creator( LoadVertex, "VRTX" );
        ringmesh_register_GocadLineParser_creator( LoadVertex, "PVRTX" );
        ringmesh_register_GocadLineParser_creator( LoadName, "name:" );
        ringmesh_register_GocadLineParser_creator( LoadTriangle, "TRGL" );
        tsolid_import_factory_initialize();
        ml_import_factory_initialize();
    }
}
// RINGMesh namespace
