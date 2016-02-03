/*
 * Copyright (c) 2012-2016, Association Scientifique pour la Geologie et ses Applications (ASGA)
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *     * Redistributions in binary form must reproduce the above copyright
 *       notice, this list of conditions and the following disclaimer in the
 *       documentation and/or other materials provided with the distribution.
 *     * Neither the name of the <organization> nor the
 *       names of its contributors may be used to endorse or promote products
 *       derived from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL <COPYRIGHT HOLDER> BE LIABLE FOR ANY
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 *
 *
 *
 *
 *
 *     http://www.ring-team.org
 *
 *     RING Project
 *     Ecole Nationale Superieure de Geologie - GeoRessources
 *     2 Rue du Doyen Marcel Roubault - TSA 70605
 *     54518 VANDOEUVRE-LES-NANCY
 *     FRANCE
 */

#include <ringmesh/geo_model_builder_so.h>

#include <iostream>
#include <iomanip>

#include <geogram/basic/logger.h>

#include <ringmesh/geometry.h>
#include <ringmesh/geogram_extension.h>
#include <ringmesh/utils.h>

/*!
* @file ringmesh/geo_model_builder_so.cpp
* @brief Implementation of the class to build GeoModel from input
* Gocad TSolid .so file
* @author Pierre Anquez
*/

namespace RINGMesh {

    class VertexMap {
    public:
        VertexMap(
            std::vector< index_t > gocad_vertices2region_vertices,
            std::vector< index_t > gocad_vertices2region_id ) :
                gocad_vertices2region_vertices_( gocad_vertices2region_vertices ),
                gocad_vertices2region_id_( gocad_vertices2region_id )
        {}

        const index_t local_id( index_t gocad_vertex_id ) const
        {
            return gocad_vertices2region_vertices_[ gocad_vertex_id ] ;
        }

        const index_t region( index_t gocad_vertex_id ) const
        {
            return gocad_vertices2region_id_[ gocad_vertex_id ] ;
        }

        void add_vertex( index_t local_vertex_id, index_t region_id )
        {
            gocad_vertices2region_vertices_.push_back( local_vertex_id ) ;
            gocad_vertices2region_id_.push_back( region_id ) ;
        }

        const index_t nb_vertex() const
        {
            ringmesh_debug_assert( gocad_vertices2region_vertices_.size() ==
                gocad_vertices2region_id_.size() )
            return gocad_vertices2region_vertices_.size() ;
        }

    private:

        /*!
         * Mapping the indices of vertices from Gocad .so file
         * to the local (in region) indices of vertices
         */
        std::vector< index_t > gocad_vertices2region_vertices_ ;
        /*!
         * Mapping the indices of vertices from Gocad .so file
         * to the region containing them
         */
        std::vector< index_t > gocad_vertices2region_id_ ;
    } ;
} // RINGMesh namespace

namespace {
    using namespace RINGMesh ;

    /*! @}
     * \name Information on the number of mesh elements from .so file
     * @{
     */

    /*!
     * @brief Counts number of vertices and tetras in each region
     * @details Reserve also space from the attributes which maps vertex
     * indices between gocad and region local indices.
     * @param[in] filename Path to the input .so file
     * @param[out] Vector built from number of vertices and tetras region
     * after region (i.e. [nb_v1, nb_t1, nb_v2, nb_t2, nb_v3, nb_t3, ...]
     * @param[out] gocad_vertices2region_vertices Vector which maps the indices
     * of vertices from Gocad .so file to the local (in region) indices of vertices
     * @param[out] gocad_vertices2region_id Vector which maps the indices of
     * vertices from Gocad .so file to the index of the region they belong to
     */
    void count_nb_vertices_and_tetras_per_region(
        const std::string& filename,
        std::vector< index_t >& nb_elements_par_region,
        index_t& nb_vertices_in_model )
    {
        nb_elements_par_region.clear() ;

        // Define a new LineInput counting number of elements
        GEO::LineInput line_input( filename ) ;

        // Initialize counters
        index_t cur_region = NO_ID ;
        index_t nb_vertices_in_region = 0 ;
        index_t nb_tetras_in_region = 0 ;
        index_t nb_surfaces_in_bmodel = 0 ;
        index_t nb_triangles_in_bmodel = 0 ;

        // Reading file
        while( !line_input.eof() && line_input.get_line() ) {
            line_input.get_fields() ;
            if( line_input.nb_fields() > 0 ) {
                if( line_input.field_matches( 0, "TVOLUME" ) ||
                        line_input.field_matches( 0, "MODEL" ) ) {
                    if( cur_region != NO_ID ) {
                        nb_elements_par_region.push_back(
                            nb_vertices_in_region ) ;
                        nb_elements_par_region.push_back( nb_tetras_in_region ) ;
                        nb_vertices_in_model += nb_vertices_in_region ;
                        nb_vertices_in_region = 0 ;
                        nb_tetras_in_region = 0 ;

                    }
                    ++cur_region ;
                } else if( line_input.field_matches( 0, "VRTX" ) ||
                           line_input.field_matches( 0, "PVRTX" ) ||
                           line_input.field_matches( 0, "ATOM" ) ||
                           line_input.field_matches( 0, "PATOM" ) ) {
                    ++nb_vertices_in_region ;
                } else if( line_input.field_matches( 0, "TETRA" ) ) {
                    ++nb_tetras_in_region ;
                }
            }
        }
    }

    /*!
     * @brief Shows number of vertices and tetras in each region
     * @param[in] nb_elements_per_region Vector built from number of vertices
     * and tetras region after region
     * (i.e. [nb_v1, nb_t1, nb_v2, nb_t2, nb_v3, nb_t3, ...]
     */
    void print_nb_vertices_and_tetras_per_region(
        const std::vector< index_t >& nb_elements_per_region )
    {
        const index_t nb_regions = 0.5 * nb_elements_per_region.size() ;
        GEO::Logger::out( "Mesh" )
            << "Mesh has " << nb_regions << " regions "
            << std::endl ;
        for( index_t i = 0 ; i < nb_regions ; ++i ) {
            GEO::Logger::out( "Mesh" )
                << "Region " << i << " has"
                << std::endl
                << std::setw( 10 ) << std::left
                << nb_elements_per_region.at( 2*i ) << " vertices "
                << std::endl
                << std::setw( 10 ) << std::left
                << nb_elements_per_region.at( 2*i + 1 ) << " tetras "
                << std::endl ;
        }
    }

    /*! @}
     * \name Properties import
     * @{
     */

    /*!
     * @brief Adds a new property to the mesh
     * @param[in] prop_name Name of the property to add
     * @param[in] attribute_manager Attribut manager in which add the new property
     */
    void add_new_property(
        const std::string& prop_name,
        GEO::AttributesManager& attribute_manager )
    {
        /// @todo All the property types are double.
        /// Change in order to have the good type for each property.
        GEO::Attribute< double > property( attribute_manager, prop_name ) ;
    }

    /*! @}
     * \name Volume mesh import
     * @{
     */

    /*!
     * @brief Creates an empty element of type GME::REGION and sets
     * its name from .so file
     * @param[in] region_name Name of the new region
     * @param[in] geomodel_builder Builder of the geomodel
     * @return The index of the initialized region
     */
    index_t initialize_region(
        const std::string& region_name,
        GeoModelBuilderTSolid& geomodel_builder )
    {
        GME::gme_t cur_region = geomodel_builder.create_element( GME::REGION ) ;
        geomodel_builder.set_element_name( cur_region, region_name ) ;
        return cur_region.index ;
    }

    /*!
     * @brief Clears the vectors region_vertices and tetra_corners and reserves
     * enough space for the next region elements
     * @param[in] nb_vertices_in_next_region Number of vertices in the
     * next region (to reverse space)
     * @param[in] nb_tetras_in_next_region Number of tetrahedra in the
     * @param[out] region_vertices Vector of the coordinates of the
     * vertices of the region to re-initialized.
     * @param[out] tetra_corners Vector of the region tetrahedra corner
     * indices to re-initialized.
     */
    void reinitialize_region_vertices_and_teras(
        const index_t nb_vertices_in_next_region,
        const index_t nb_tetras_in_next_region,
        std::vector < vec3 >& region_vertices,
        std::vector < index_t >& tetra_corners )
    {
        region_vertices.clear() ;
        region_vertices.reserve( nb_vertices_in_next_region ) ;
        tetra_corners.clear() ;
        tetra_corners.reserve( 4 * nb_tetras_in_next_region ) ;
    }

    /*! @}
     * \name Read mesh elements (points, triangles, tetrehedra)
     * @{
     */

    /*!
     * Reads the coordinates of a vertex from file
     * @param[in] in ACSII file reader
     * @param[in] z_sign FActor for z value in order to have z increasing upwards
     * @param[out] vertex Vertex
     */
    void read_vertex_coordinates(
        const GEO::LineInput& in,
        const int z_sign,
        vec3& vertex )
    {
        vertex.x = in.field_as_double( 2 ) ;
        vertex.y = in.field_as_double( 3 ) ;
        vertex.z = z_sign * in.field_as_double( 4 ) ;
    }

    /*!
     * @brief Reads the three vertices index of a triangle and adds
     * them to the facet corners of the currently built surface
     * @details Reads gocad indices
     * @param[in] in ACSII file reader
     * @param[out] cur_surf_facets Vector of each facet corner indices
     * to build facets
     */
    void read_triangle(
        const GEO::LineInput& in,
        std::vector< index_t >& cur_surf_facets )
    {
        cur_surf_facets.push_back( in.field_as_uint( 1 ) - 1 ) ;
        cur_surf_facets.push_back( in.field_as_uint( 2 ) - 1 ) ;
        cur_surf_facets.push_back( in.field_as_uint( 3 ) - 1 ) ;
    }

    /*!
     * @brief Reads the four vertices index of a tetrahedron
     * @details Reads gocad indices (from .so file) and transforms
     * them to vertex local (region) indices
     * @param[in] in ACSII file reader
     * @param[out] gocad_vertices2region_vertices Vector which maps the indices
     * of vertices from Gocad .so file to the local (in region) indices of vertices
     * @param[out] corners_id Indices of the four vertices
     */
    void read_tetraedra(
        const GEO::LineInput& in,
        const VertexMap& vertex_map,
        std::vector< index_t >& corners_id )
    {
        ringmesh_debug_assert( corners_id.size() == 4 ) ;
        corners_id[0] = vertex_map.local_id( in.field_as_uint( 1 ) - 1 ) ;
        corners_id[1] = vertex_map.local_id( in.field_as_uint( 2 ) - 1 ) ;
        corners_id[2] = vertex_map.local_id( in.field_as_uint( 3 ) - 1 ) ;
        corners_id[3] = vertex_map.local_id( in.field_as_uint( 4 ) - 1 ) ;
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
    void compute_region_cell_facet_centers(
        const GeoModel& geomodel,
        const index_t region_id,
        std::vector< vec3 >& cell_facet_centers )
    {
        const Region& region = geomodel.region( region_id ) ;
        const index_t nb_cells = region.nb_cells() ;
        cell_facet_centers.reserve( 4 * nb_cells ) ;
        for( index_t c = 0 ; c < nb_cells ; ++c ) {
            for( index_t f = 0 ; f <= 3 ; ++f ) {
                cell_facet_centers.push_back(
                    mesh_cell_facet_center( region.mesh(), c, f ) ) ;
            }
        }
    }

    /*!
     * @brief Computes the colocaters of the centers of cell facets for
     * each region
     * * @param[in] geomodel GeoModel to consider
     * @param[out] region_anns Pointers to the ColocaterANNs of regions
     */
    void compute_cell_facet_centers_region_anns(
        const GeoModel& geomodel,
        std::vector< ColocaterANN* >& region_anns )
    {
        for( index_t r = 0 ; r < geomodel.nb_regions() ; ++r ) {
            std::vector< vec3 > cell_facet_centers ;
            compute_region_cell_facet_centers( geomodel, r, cell_facet_centers ) ;
            region_anns[r] = new ColocaterANN( cell_facet_centers, true ) ;
        }
    }

    /*!
     * @brief Tests if a surface is a boundary of a region.
     * @details If it is the case, add the surface to the boundaries of
     * the region and the region to the in_boundaries of the surface
     * @param[in] surface Surface to test
     * @param[in] region_ann Vector of ColocaterANN of the region to test
     * @param[out] colocated_cell_facet_centers Vector of colocated cell
     * facet centers
     * @return The number of surface sides bounding the region
     */
    index_t are_surface_sides_region_boundaries(
        const Surface& surface,
        const ColocaterANN& region_ann,
        std::vector< index_t >& colocated_cell_facet_centers )
    {
        vec3 first_facet_center = surface.facet_barycenter( 0 ) ;
        region_ann.get_colocated( first_facet_center,
            colocated_cell_facet_centers ) ;
        return colocated_cell_facet_centers.size() ;
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
        const GeoModel& geomodel,
        const index_t region_id,
        const index_t surface_id,
        const index_t cell_facet_center_id )
    {
        index_t local_facet_id = cell_facet_center_id % 4 ;
        index_t cell_id =
            0.25 * ( cell_facet_center_id - local_facet_id ) ;
        vec3 cell_facet_normal =
            mesh_cell_facet_normal(
                geomodel.region( region_id ).mesh(),
                cell_id,
                local_facet_id ) ;
        vec3 first_facet_normal =
            geomodel.surface( surface_id ).facet_normal( 0 ) ;
        return dot( first_facet_normal, cell_facet_normal ) > 0 ;
    }

    /*!
     * @brief Both adds the surface in the boundaries of a region and
     * adds the region to the in_boundaries of the surface
     * @param[in] region_id Index of the region
     * @param[in] surface_id Index of the surface
     * @param[in] surf_side Side of the surface bounding the region
     * @param[in] geomodel_builder Builder of the GeoModel to consider
     */
    void fill_region_and_surface_boundaries_links(
        const index_t region_id,
        const index_t surface_id,
        const bool surf_side,
        GeoModelBuilderTSolid& geomodel_builder )
    {
        geomodel_builder.add_element_boundary(
            GME::gme_t( GME::REGION, region_id ),
            GME::gme_t( GME::SURFACE, surface_id ),
            surf_side ) ;
        geomodel_builder.add_element_in_boundary(
            GME::gme_t( GME::SURFACE, surface_id ),
            GME::gme_t( GME::REGION, region_id ) ) ;
    }

    /*!
     * @brief Adds the both surface sides in the boundaries of a region
     * (internal boundary) and add twice the region to the in_boundaries
     * of the surface
     * @param[in] region_id Index of the region
     * @param[in] surface_id Index of the surface
     * @param[in,out] geomodel_builder Builder of the GeoModel to consider
     */
    void add_both_surface_sides_to_region_boundaries(
        const index_t region_id,
        const index_t surface_id,
        GeoModelBuilderTSolid& geomodel_builder )

    {
        fill_region_and_surface_boundaries_links(
            region_id,
            surface_id,
            true,
            geomodel_builder ) ;
        fill_region_and_surface_boundaries_links(
            region_id,
            surface_id,
            false,
            geomodel_builder ) ;
    }

    /*!
     * @brief Adds one surface side in the boundaries of a region
     * and add the region to the in_boundaries of the surface
     * @details The index of the cell facet center is used for the
     * determination of the side to add.
     * @param[in] region_id Index of the region
     * @param[in] surface_id Index of the surface
     * @param[in] cell_facet_center_id Index of the cell facet center
     * (i.e., cell_id * 4 + local_facet_id)
     * @param[in,out] geomodel_builder Builder of the GeoModel to consider
     */
    void add_one_surface_side_to_region_boundaries(
        const index_t region_id,
        const index_t surface_id,
        const index_t cell_facet_center_id,
        GeoModelBuilderTSolid& geomodel_builder )
    {
        bool side = determine_surface_side_to_add(
            geomodel_builder.model(),
            region_id,
            surface_id,
            cell_facet_center_id ) ;
        fill_region_and_surface_boundaries_links(
            region_id,
            surface_id,
            side,
            geomodel_builder ) ;
    }

    /*!
     * @brief Adds the surface sides which bound the region to the
     * boundaries of the region (and add the region to in boundaries
     * of the surface)
     * @param[in] surface_id Index of the surface
     * @param[in] region_id Index of the region
     * @param[in] colocated_cell_facet_centers Vector of colocated cell
     * facet centers
     * @param[in,out] geomodel_builder Builder of the GeoModel to consider
     */
    void add_surface_sides_to_region_boundaries(
        const index_t surface_id,
        const index_t region_id,
        const std::vector< index_t >& colocated_cell_facet_centers,
        GeoModelBuilderTSolid& geomodel_builder )
    {
        switch( colocated_cell_facet_centers.size() ) {
            case 1 :
                add_one_surface_side_to_region_boundaries(
                    region_id,
                    surface_id,
                    colocated_cell_facet_centers[0],
                    geomodel_builder ) ;
                break ;
            case 2 :
                add_both_surface_sides_to_region_boundaries(
                    region_id,
                    surface_id,
                    geomodel_builder ) ;
                break ;
            default :
                ringmesh_assert_not_reached ;
        }
    }

    /*!
     * @brief Sets the given surface as regions boundaries
     * @details Based on ColocaterANN, retrieves the regions bounded by the
     * given surface. One side or the both sides of the surface
     * could bound model regions.
     * @param[in] surface_id Index of the surface
     * @param[in] region_anns Vector of ColocaterANN of the model regions
     * @param[in,out] geomodel_builder Builder of the GeoModel to consider
     */
    void add_surface_to_region_boundaries(
        const index_t surface_id,
        const std::vector< ColocaterANN* >& region_anns,
        GeoModelBuilderTSolid& geomodel_builder )
    {
        const GeoModel& model = geomodel_builder.model() ;
        index_t cur_region = 0 ;
        index_t nb_added_surf_sides = 0 ;
        // Maximum 2 regions could be bounded by a single surface
        while ( cur_region < model.nb_regions() && nb_added_surf_sides < 2 ) {
            std::vector< index_t > colocated_cell_facet_centers ;
            index_t nb_surf_sides_are_boundary =
                are_surface_sides_region_boundaries(
                    model.surface( surface_id ),
                    *region_anns[cur_region],
                    colocated_cell_facet_centers ) ;
            if ( nb_surf_sides_are_boundary > 0 ) {
                add_surface_sides_to_region_boundaries(
                    surface_id,
                    cur_region,
                    colocated_cell_facet_centers,
                    geomodel_builder ) ;
                nb_added_surf_sides += nb_surf_sides_are_boundary ;
            }
            ++cur_region ;
        }
        if ( nb_added_surf_sides == 0 ) {
            ringmesh_assert_not_reached ;
        }
    }

    /*!
     * @brief Sets the boundaries of the GeoModel regions
     * @param[in,out] geomodel_builder Builder of the GeoModel to consider
     */
    void compute_boundaries_of_geomodel_regions(
        GeoModelBuilderTSolid& geomodel_builder )
    {
        const GeoModel& model = geomodel_builder.model() ;
        std::vector< ColocaterANN* > reg_anns( model.nb_regions(), nil ) ;
        compute_cell_facet_centers_region_anns( model, reg_anns ) ;
        for( index_t s = 0 ; s < model.nb_surfaces() ; ++s ) {
            add_surface_to_region_boundaries( s, reg_anns, geomodel_builder ) ;
        }
        for( index_t r = 0 ; r < model.nb_regions() ; ++r ) {
            delete reg_anns[r] ;
        }
    }

    /*
     * @brief Adds the right surface sides in universe boundaries
     * @param[in] surf_side_minus Vector indicating if the '-' side of
     * surfaces are in the boundaries of model regions
     * @param[in] surf_side_plus Vector indicating if the '+' side of
     * surfaces are in the boundaries of model regions
     * @param[in,out] geomodel_builder Builder of the GeoModel to consider
     */
    void add_surfaces_to_universe_boundaries(
        const std::vector< bool >& surf_minus_side,
        const std::vector< bool >& surf_plus_side,
        GeoModelBuilderTSolid& geomodel_builder )
    {
        const GeoModel& model = geomodel_builder.model() ;
        for( index_t s = 0 ; s < model.nb_surfaces() ; ++s ) {
            if( surf_minus_side[s] && !surf_plus_side[s] ) {
                geomodel_builder.add_element_boundary(
                    GME::gme_t( GME::REGION, NO_ID ),
                    GME::gme_t( GME::SURFACE, s ),
                    false ) ;
            } else if( !surf_minus_side[s] && surf_plus_side[s] ) {
                geomodel_builder.add_element_boundary(
                    GME::gme_t( GME::REGION, NO_ID ),
                    GME::gme_t( GME::SURFACE, s ),
                    true ) ;
            }
        }
    }

    /*
     * @brief Determines if each side of the surfaces are
     * in the boundaries of model regions
     * @param[in] geomodel GeoModel to consider
     * @param[out] surf_side_minus Vector indicating if the '-' side of
     * surfaces are in the boundaries of model regions
     * @param[out] surf_side_plus Vector indicating if the '+' side of
     * surfaces are in the boundaries of model regions
     */
    void determine_if_surface_sides_bound_regions(
        const GeoModel& geomodel,
        std::vector< bool >& surf_minus_side,
        std::vector< bool >& surf_plus_side )
    {
        for( index_t r = 0 ; r < geomodel.nb_regions() ; ++r ) {
            for( index_t s = 0 ; s < geomodel.region(r).nb_boundaries() ; ++s ) {
                if( geomodel.region(r).side(s) ) {
                    surf_plus_side[ geomodel.region(r).boundary(s).index() ] = true ;
                } else if( !geomodel.region(r).side(s) ) {
                    surf_minus_side[ geomodel.region(r).boundary(s).index() ] = true ;
                } else {
                    ringmesh_assert_not_reached
                }
            }
        }
    }

    /*!
     * @brief Sets the boundaries of region Universe
     * @details A surface is set in the boundaries of region Universe if
     * only one of its sides belongs to the boundaries of other regions.
     * @param[in,out] geomodel_builder Builder of the GeoModel to consider
     */
    void compute_universe_boundaries( GeoModelBuilderTSolid& geomodel_builder )
    {
        const GeoModel& model = geomodel_builder.model() ;
        // The universe boundaries are the surfaces with only one side in all
        // the boundaries of the other regions
        std::vector< bool > surf_minus_side( model.nb_surfaces(), false ) ;
        std::vector< bool > surf_plus_side( model.nb_surfaces(), false ) ;
        determine_if_surface_sides_bound_regions(
            model,
            surf_minus_side,
            surf_plus_side ) ;
        add_surfaces_to_universe_boundaries(
            surf_minus_side,
            surf_plus_side,
            geomodel_builder ) ;
    }

    /*! @}
     * \name Surface internal borders determination
     * @{
     */

    /*!
     * @brief Finds if a surface facet edge is an internal border
     * (i.e. shared by at least two surfaces)
     * @param[in] geomodel GeoModel to consider
     * @param[in] surface_id Index of the surface
     * @param[in] facet Index of the facet in the surface
     * @param[in] edge Index of the edge in the facet
     * @param[in] surface_anns Pointers to the ColocaterANNs of surfaces
     * @param[in] surface_boxes Bounding Box of surfaces
     * @return True is the edge is found in at least another surface
     */
    bool is_edge_in_several_surfaces(
        const GeoModel& geomodel,
        const index_t surface_id,
        const index_t facet,
        const index_t edge,
        const std::vector< ColocaterANN* >& surface_anns,
        const std::vector< Box3d >& surface_boxes )
    {
        const Surface& S = geomodel.surface( surface_id ) ;
        const vec3 barycenter = GEO::Geom::barycenter(
            S.vertex( facet, edge ),
            S.vertex( facet, ( edge+1 ) % 3 ) ) ;
       std::vector< index_t > result ;
       index_t tested_surf = 0 ;
       while( result.empty() &&
           tested_surf < surface_anns.size() ) {
           if( surface_boxes[tested_surf].contains( barycenter ) ) {
               surface_anns[tested_surf]->
                   get_colocated( barycenter, result ) ;
           }
           ++tested_surf ;
       }
       return !result.empty() ;
    }

    /*!
     * @brief Computes internal borders of a given surface
     * @details A surface facet edge is an internal border if it is shared
     * by at least two surfaces. Adjacency of such a facet edge is set to
     * GEO::NO_FACET.
     * @param[in] geomodel GeoModel to consider
     * @param[in] surface_id Index of the surface
     * @param[in] surface_anns Pointers to the ColocaterANNs of surfaces
     */
    void compute_surface_internal_borders(
        const GeoModel& geomodel,
        const index_t surface_id,
        const std::vector< ColocaterANN* >& surface_anns,
        const std::vector< Box3d >& surface_boxes )
    {
        const Surface& S = geomodel.surface( surface_id ) ;
        for( index_t f = 0 ; f < S.nb_cells() ; ++f ) {
            for( index_t e = 0 ; e < 3 ; ++e ) {
               if( !S.is_on_border(f,e) ) {
                   bool internal_border = is_edge_in_several_surfaces(
                       geomodel, surface_id, f, e, surface_anns, surface_boxes ) ;
                   if( internal_border ) {
                       S.mesh().facets.set_adjacent( f, e, GEO::NO_FACET ) ;
                   }
               }
            }
        }
    }

    /*!
     * @brief Gets the facet edge barycenters of a given surface
     * @param[in] geomodel GeoModel to consider
     * @param[in] surface_id Index of the surface
     * @param[out] border_edge_barycenters Vector of all the border
     * edge barycenters of the surface
     */
    void get_surface_border_edge_barycenters(
        const GeoModel& geomodel,
        const index_t surface_id,
        std::vector< vec3 >& border_edge_barycenters )
    {
        const Surface& S = geomodel.surface( surface_id ) ;
        for( index_t f = 0 ; f < S.nb_cells() ; ++f ) {
            for( index_t e = 0 ; e < 3 ; ++e ) {
                if(S.is_on_border(f,e)) {
                    const vec3 barycenter = GEO::Geom::barycenter(
                        S.vertex( f, e ),
                        S.vertex( f, ( e+1 ) % 3 ) ) ;
                    border_edge_barycenters.push_back( barycenter ) ;
                }
            }
        }
    }

    /*!
     * @brief Computes the colocaters of the centers of facet edges for
     * each surface and their Box3d
     * @param[in] geomodel GeoModel to consider
     * @param[out] surface_anns Pointers to the ColocaterANNs of surfaces
     * @param[out] surface_boxes Bounding Box of surfaces
     */
    void compute_facet_edge_centers_anns_and_surface_boxes(
        const GeoModel& geomodel,
        std::vector< ColocaterANN* >& surface_anns,
        std::vector< Box3d >& surface_boxes )
    {
        for( index_t s = 0 ; s < geomodel.nb_surfaces() ; ++s ) {
            const Surface& S = geomodel.surface(s) ;
            for( index_t p = 0; p < S.nb_vertices(); p++ ) {
                surface_boxes[s].add_point( S.vertex( p ) ) ;
            }
            std::vector< vec3 > border_edge_barycenters ;
            get_surface_border_edge_barycenters(
                geomodel, s, border_edge_barycenters ) ;
            surface_anns[s] =
                new ColocaterANN( border_edge_barycenters, true ) ;
        }
    }

    /*!
     * @brief Computes internal borders of the model surfaces
     * @details An surface facet edge is an internal border if it is shared
     * by at least two surfaces. Adjacency of such a facet edge is set to
     * GEO::NO_FACET.
     * @param[in] geomodel GeoModel to consider
     */
    void compute_surfaces_internal_borders( const GeoModel& geomodel )
    {
        std::vector< ColocaterANN* > anns( geomodel.nb_surfaces(), nil ) ;
        std::vector< Box3d > boxes( geomodel.nb_surfaces() ) ;
        compute_facet_edge_centers_anns_and_surface_boxes( geomodel, anns, boxes ) ;
        for( index_t s = 0 ; s < geomodel.nb_surfaces() ; ++s ) {
            compute_surface_internal_borders( geomodel, s, anns, boxes ) ;
        }
        for( index_t s = 0 ; s < geomodel.nb_surfaces() ; ++s ) {
            delete anns[s] ;
        }
    }

    /*!
     * \name Reads and sets Gocad Coordinate System information from .so file
     * @{
     */

    /*!
     * @brief Sets the Gocad coordinate system axis name
     * @param[in] in ACSII file reader
     */
    void set_gocad_coordinate_system_axis_name(
        GEO::LineInput& in,
        std::vector< std::string >& gocad_coordinates_system_axis_name )
    {
        gocad_coordinates_system_axis_name.push_back( in.field(1) ) ;
        gocad_coordinates_system_axis_name.push_back( in.field(2) ) ;
        gocad_coordinates_system_axis_name.push_back( in.field(3) ) ;
    }

    /*!
     * @brief Sets the Gocad coordinate system axis units (meters, feet, ...)
     * @param[in] in ACSII file reader
     */
    void set_gocad_coordinate_system_axis_unit(
        GEO::LineInput& in,
        std::vector< std::string >& gocad_coordinates_system_axis_unit )
    {
        gocad_coordinates_system_axis_unit.push_back( in.field(1) ) ;
        gocad_coordinates_system_axis_unit.push_back( in.field(2) ) ;
        gocad_coordinates_system_axis_unit.push_back( in.field(3) ) ;
    }

    /*!
     * @brief Sets the Gocad coordinate system Z sign
     * @details If the Z values increase upwards, z_sign_ is positive,
     * else (Z values increasing downwards), z_sign_ is negative.
     * @param[in] in ACSII file reader
     */
    int set_gocad_coordinate_system_z_sign(
        GEO::LineInput& in )
    {
        if( in.field_matches( 1, "Elevation" ) ) {
            return 1 ;
        } else if( in.field_matches( 1, "Depth" ) ) {
            return -1 ;
        } else {
            ringmesh_assert_not_reached ;
        }
    }

    /*! @}
     */

} // anonymous namespace

namespace RINGMesh {

    bool GeoModelBuilderTSolid::load_file()
    {
        if( !in_.OK() ) {
            return false ;
        }

        index_t cur_region ;

        index_t nb_vertices_in_model ;

        // First : count the number of vertex and tetras
        // in each region
        std::vector< index_t > nb_elements_per_region ;
        count_nb_vertices_and_tetras_per_region(
            filename_,
            nb_elements_per_region,
            nb_vertices_in_model ) ;

        std::vector< index_t > gocad_vertices2region_vertices ;
        gocad_vertices2region_vertices.reserve( nb_vertices_in_model ) ;
        std::vector< index_t > gocad_vertices2region_id ;
        gocad_vertices2region_id.reserve( nb_vertices_in_model ) ;

        VertexMap vertex_map ( gocad_vertices2region_vertices,
            gocad_vertices2region_id ) ;

        // Region vertices
        std::vector< vec3 > region_vertices ;
        region_vertices.reserve( nb_elements_per_region[0] ) ;
        // Region tetrahedron corners
        std::vector< index_t > tetra_corners ;
        tetra_corners.reserve( 4 * nb_elements_per_region[0] ) ;

        GME::gme_t current_interface ;
        GME::gme_t current_surface ;
        std::vector< index_t > cur_surf_facets ;
        std::vector< index_t > cur_surf_facet_ptr (1, 0) ;

        index_t nb_vertex_properties = 0 ;
        index_t nb_cell_properties = 0 ;

        // Then : Reading .so file
        while( !in_.eof() && in_.get_line() ) {
            in_.get_fields() ;
            if( in_.nb_fields() > 0 ) {
                if( in_.field_matches( 0, "GOCAD_ORIGINAL_COORDINATE_SYSTEM" ) ) {
                    read_and_set_gocad_coordinate_system() ;
                } else if( in_.field_matches( 0, "PROPERTIES" ) ) {
                    nb_vertex_properties = in_.nb_fields() - 1 ;
                } else if( in_.field_matches( 0, "PROPERTY_CLASS_HEADER" ) ) {
                    add_new_property( in_.field( 1 ),
                        model_.mesh.vertex_attribute_manager() ) ;
                } else if( in_.field_matches( 0, "TETRA_PROPERTIES" ) ) {
                    nb_cell_properties = in_.nb_fields() - 1 ;
                } else if( in_.field_matches( 0, "TETRA_PROPERTY_CLASS_HEADER" ) ) {
                    add_new_property( in_.field( 1 ),
                        model_.mesh.cell_attribute_manager() ) ;
                } else if( in_.field_matches( 0, "TVOLUME" ) ) {
                    if( region_vertices.size() > 0 ) {
                        set_region_geometry(
                            cur_region,
                            region_vertices,
                            tetra_corners ) ;
                        reinitialize_region_vertices_and_teras(
                            nb_elements_per_region[ 2*cur_region ],
                            nb_elements_per_region[ 2*cur_region + 1 ],
                            region_vertices,
                            tetra_corners ) ;
                    }
                    cur_region = initialize_region( in_.field( 1 ), *this ) ;
                } else if( in_.field_matches( 0, "VRTX" ) ||
                    in_.field_matches( 0, "PVRTX" ) ) {
                    read_and_add_vertex_to_region_vertices(
                        cur_region,
                        region_vertices,
                        vertex_map ) ;
                } else if( in_.field_matches( 0, "ATOM" ) ||
                    in_.field_matches( 0, "PATOM" ) ) {
                    read_and_add_atom_to_region_vertices(
                        cur_region,
                        region_vertices,
                        vertex_map ) ;
                } else if( in_.field_matches( 0, "TETRA" ) ) {
                    // Reading and create a tetra
                    std::vector< index_t > corners(4) ;
                    read_tetraedra(
                        in_, vertex_map, corners ) ;
                    tetra_corners.insert( tetra_corners.end(),
                        corners.begin(),
                        corners.end() ) ;
                } else if( in_.field_matches( 0, "name:" ) ) {
                    // GeoModel name is set to the TSolid name.
                    set_model_name( in_.field( 1 ) ) ;
                } else if( in_.field_matches( 0, "MODEL" ) ) {
                    if( region_vertices.size() > 0 ) {
                        set_region_geometry(
                            cur_region,
                            region_vertices,
                            tetra_corners ) ;
                        reinitialize_region_vertices_and_teras(
                            0,
                            0,
                            region_vertices,
                            tetra_corners ) ;
                    }
                } else if( in_.field_matches( 0, "SURFACE" ) ) {
                    current_interface = create_element( GME::INTERFACE ) ;
                    set_element_name( current_interface, in_.field( 1 ) ) ;
                } else if( in_.field_matches( 0, "TFACE" ) ) {
                    // Compute the surface
                    if( cur_surf_facets.size() > 0 ) {
                        build_surface( vertex_map,
                            current_surface.index,
                            cur_surf_facets,
                            cur_surf_facet_ptr ) ;
                    }
                    // Create a new surface
                    current_surface = create_element( GME::SURFACE ) ;
                    set_element_parent( current_surface, current_interface ) ;
                    add_element_child( current_interface, current_surface ) ;
                } else if( in_.field_matches( 0, "TRGL" ) ) {
                    read_triangle( in_, cur_surf_facets ) ;
                    cur_surf_facet_ptr.push_back( cur_surf_facets.size() ) ;
                } else if( in_.field_matches( 0, "MODEL_REGION" ) ) {
                    // Compute the last surface
                    if( cur_surf_facets.size() > 0 ) {
                        build_surface( vertex_map,
                            current_surface.index,
                            cur_surf_facets,
                            cur_surf_facet_ptr ) ;
                    }
                }
            }
        }

        // Compute internal borders (by removing adjacencies on
        // triangle edges common to at least two surfaces)
        compute_surfaces_internal_borders( (*this).model() ) ;

        // Build GeoModel Lines and Corners from the surfaces
        model_.mesh.vertices.test_and_initialize() ;
        build_lines_and_corners_from_surfaces() ;

        // Regions boundaries
        compute_boundaries_of_geomodel_regions( *this ) ;

        // Universe boundaries
        compute_universe_boundaries( *this ) ;

        // Contacts building
        build_contacts() ;

        return true ;

    }

    void GeoModelBuilderTSolid::read_and_set_gocad_coordinate_system()
    {
        while( !in_.eof() && in_.get_line() ) {
            in_.get_fields() ;
            if( in_.field_matches( 0, "END_ORIGINAL_COORDINATE_SYSTEM" ) ) {
                return ;
            } else if( in_.field_matches( 0, "NAME" ) ) {
                // Useless for the moment
            } else if( in_.field_matches( 0, "PROJECTION" ) ) {
                // Useless for the moment
            } else if( in_.field_matches( 0, "DATUM" ) ) {
                // Useless for the moment
            } else if( in_.field_matches( 0, "AXIS_NAME" ) ) {
                set_gocad_coordinate_system_axis_name(
                    in_, gocad_coordinates_system_axis_name_ ) ;
            } else if( in_.field_matches( 0, "AXIS_UNIT" ) ) {
                set_gocad_coordinate_system_axis_unit(
                    in_, gocad_coordinates_system_axis_unit_ ) ;
            } else if( in_.field_matches( 0, "ZPOSITIVE" ) ) {
                z_sign_ = set_gocad_coordinate_system_z_sign( in_ ) ;
            }
        }
    }

    void GeoModelBuilderTSolid::read_and_add_vertex_to_region_vertices(
        const index_t region_id,
        std::vector < vec3 >& region_vertices,
        VertexMap& vertex_map )
    {
        vertex_map.add_vertex( region_vertices.size(), region_id ) ;
        vec3 vertex ;
        read_vertex_coordinates( in_, z_sign_, vertex ) ;
        region_vertices.push_back( vertex ) ;
    }

    void GeoModelBuilderTSolid::read_and_add_atom_to_region_vertices(
        const index_t region_id,
        std::vector < vec3 >& region_vertices,
        VertexMap& vertex_map )
    {
        const index_t referring_vertex = in_.field_as_double( 2 ) - 1 ;
        const index_t referred_vertex_local_id =
            vertex_map.local_id( referring_vertex ) ;
        const index_t referred_vertex_region_id =
            vertex_map.region( referring_vertex ) ;
        if( referred_vertex_region_id < region_id ) {
            // If the atom referred to a vertex of another region,
            // acting like for a vertex
            vertex_map.add_vertex( region_vertices.size(), region_id ) ;
            region_vertices.push_back(
                    model_.region( referred_vertex_region_id ).vertex(
                            referred_vertex_local_id ) ) ;
        } else {
            // If the atom referred to an atom of the same region
            vertex_map.add_vertex(
                referred_vertex_local_id, referred_vertex_region_id ) ;
        }
    }

    void GeoModelBuilderTSolid::build_surface(
        const VertexMap& vertex_map,
        const index_t surface_id,
        std::vector< index_t >& facet_corners,
        std::vector< index_t >& facet_ptr )
    {
        std::vector< vec3 > cur_surf_points ;
        std::vector< index_t > cur_surf_facets ;
        get_surface_points_and_facets_from_gocad_indices(
            vertex_map,
            facet_corners,
            cur_surf_points,
            cur_surf_facets ) ;
        set_surface_geometry( surface_id,
            cur_surf_points,
            cur_surf_facets,
            facet_ptr ) ;
        cur_surf_points.clear() ;
        facet_corners.clear() ;
        facet_ptr.clear() ;
        facet_ptr.push_back( 0 ) ;
    }

    void GeoModelBuilderTSolid::get_surface_points_and_facets_from_gocad_indices(
        const VertexMap& vertex_map,
        const std::vector< index_t >& facet_corners,
        std::vector< vec3 >& cur_surf_points,
        std::vector< index_t >& cur_surf_facets ) const
    {
        std::vector< index_t > gocad_vertices2cur_surf_points(
            vertex_map.nb_vertex(), NO_ID ) ;
        for( index_t co = 0 ; co < facet_corners.size() ; ++co ) {
            const index_t corner_gocad_id = facet_corners[ co ] ;
            get_surface_point_and_facet_from_gocad_index(
                corner_gocad_id,
                vertex_map,
                gocad_vertices2cur_surf_points,
                cur_surf_points,
                cur_surf_facets ) ;

        }
    }
    void GeoModelBuilderTSolid::get_surface_point_and_facet_from_gocad_index(
        const index_t vertex_gocad_id,
        const VertexMap& vertex_map,
        std::vector< index_t >& gocad_vertices2cur_surf_points,
        std::vector< vec3 >& cur_surf_points,
        std::vector< index_t >& cur_surf_facets ) const
    {
        if( gocad_vertices2cur_surf_points[ vertex_gocad_id ] == NO_ID ) {
            // First time this facet corner is met in facet_corners
            vec3 point ;
            get_point_from_gocad_id( vertex_gocad_id, vertex_map, point ) ;
            cur_surf_facets.push_back( cur_surf_points.size() ) ;
            gocad_vertices2cur_surf_points[ vertex_gocad_id ] =
                cur_surf_points.size() ;
            cur_surf_points.push_back( point ) ;
        } else {
            // If this facet corner has already been met in facet_corners
            cur_surf_facets.push_back(
                gocad_vertices2cur_surf_points[ vertex_gocad_id ] ) ;
        }
    }

    void GeoModelBuilderTSolid::get_point_from_gocad_id(
        const index_t point_gocad_id,
        const VertexMap& vertex_map,
        vec3& point ) const
    {
        const index_t point_local_id = vertex_map.local_id( point_gocad_id ) ;
        const index_t point_region = vertex_map.region( point_gocad_id ) ;
        point =
            model_.region( point_region ).vertex( point_local_id ) ;
    }
} // RINGMesh namespace
