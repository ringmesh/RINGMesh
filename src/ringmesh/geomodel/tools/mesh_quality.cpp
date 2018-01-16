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
#include <ringmesh/geomodel/core/geomodel.h>
#include <ringmesh/geomodel/core/geomodel_mesh_entity.h>
#include <ringmesh/geomodel/tools/mesh_quality.h>
#include <ringmesh/mesh/mesh.h>
#include <ringmesh/mesh/mesh_builder.h>
#include <ringmesh/mesh/mesh_index.h>

/*!
 * @author Benjamin Chauvin
 * This code is inspired from
 * http://people.sc.fsu.edu/~jburkardt/cpp_src/tet_mesh_quality/tet_mesh_quality.html
 */

namespace
{
    using namespace RINGMesh;

    /*!
     * @brief Computes the radius of the tetrahedron insphere.
     *
     * The tetrahedron insphere is the sphere inside the tetrahedron which
     * is tangent to each tetrahedron facet.
     *
     * @param[in] v0 first vertex of the tetrahedron.
     * @param[in] v1 second vertex of the tetrahedron.
     * @param[in] v2 third vertex of the tetrahedron.
     * @param[in] v3 fourth vertex of the tetrahedron.
     * @return the radius of the tetrahedron insphere.
     */
    double tetra_insphere_radius(
        const vec3& v0, const vec3& v1, const vec3& v2, const vec3& v3 )
    {
        double tet_volume = GEO::Geom::tetra_volume( v0, v1, v2, v3 );
        double A1 = GEO::Geom::triangle_area( v0, v1, v2 );
        double A2 = GEO::Geom::triangle_area( v1, v2, v3 );
        double A3 = GEO::Geom::triangle_area( v2, v3, v0 );
        double A4 = GEO::Geom::triangle_area( v3, v0, v1 );
        ringmesh_assert( A1 + A2 + A3 + A4 > global_epsilon );
        return ( 3 * tet_volume ) / ( A1 + A2 + A3 + A4 );
    }

    /*!
     * @brief Tetrahedron quality based on the insphere and circumsphere radii.
     *
     * The quality is between 0 and 1. 0 corresponds to a bad tetrahedron, and
     * 1 to a good tetrahedron (equilaterality).
     * For more information, see
     * <a
     * href="http://people.sc.fsu.edu/~jburkardt/cpp_src/tet_mesh_quality/tet_mesh_quality.html">
     * TET_MESH_QUALITY Interactive Program for Tet Mesh Quality</a>
     *
     * @param[in] v0 first vertex of the tetrahedron.
     * @param[in] v1 second vertex of the tetrahedron.
     * @param[in] v2 third vertex of the tetrahedron.
     * @param[in] v3 fourth vertex of the tetrahedron.
     * @return 3 * the insphere radius divided by the circumsphere radius.
     */
    double tet_quality_insphere_radius_by_circumsphere_radius(
        const vec3& v0, const vec3& v1, const vec3& v2, const vec3& v3 )
    {
        const vec3 tetra_circum_center =
            GEO::Geom::tetra_circum_center( v0, v1, v2, v3 );
        const double tetra_circum_radius =
            vec3( tetra_circum_center - v0 ).length();
        ringmesh_assert( std::abs( tetra_circum_radius
                                   - ( tetra_circum_center - v1 ).length() )
                         < global_epsilon );
        ringmesh_assert( std::abs( tetra_circum_radius
                                   - ( tetra_circum_center - v2 ).length() )
                         < global_epsilon );
        ringmesh_assert( std::abs( tetra_circum_radius
                                   - ( tetra_circum_center - v3 ).length() )
                         < global_epsilon );

        // insphere computation
        double in_radius = tetra_insphere_radius( v0, v1, v2, v3 );
        return 3. * in_radius / tetra_circum_radius;
    }

    /*!
     * @param[in] v0 first vertex of the tetrahedron.
     * @param[in] v1 second vertex of the tetrahedron.
     * @param[in] v2 third vertex of the tetrahedron.
     * @param[in] v3 fourth vertex of the tetrahedron.
     * @return the maximum of the tetrahedron edge length.
     */
    double max_tet_edge_length(
        const vec3& v0, const vec3& v1, const vec3& v2, const vec3& v3 )
    {
        double l1 = ( v1 - v0 ).length();
        double l2 = ( v2 - v0 ).length();
        double l3 = ( v3 - v0 ).length();
        double l4 = ( v2 - v1 ).length();
        double l5 = ( v3 - v1 ).length();
        double l6 = ( v3 - v2 ).length();
        return GEO::geo_max(
            GEO::geo_max(
                GEO::geo_max( GEO::geo_max( GEO::geo_max( l1, l2 ), l3 ), l4 ),
                l5 ),
            l6 );
    }

    /*!
     * @brief Tetrahedron quality based on the insphere radius and the maximum
     * edge length.
     *
     * The quality is between 0 and 1. 0 corresponds to a bad tetrahedron, and
     * 1 to a good tetrahedron (equilaterality).
     * For more information see
     * Du, Q., and D. Wang, 2005,
     * The optimal centroidal Voronoi tessellations and the gersho's conjecture
     * in the three-dimensional space,
     * Computers & Mathematics with Applications, v. 49, no. 9, p. 1355-1373,
     * <a href="http://doi.org/10.1016/j.camwa.2004.12.008">doi</a>,
     * <a
     * href="http://people.sc.fsu.edu/~jburkardt/cpp_src/tet_mesh_quality/tet_mesh_quality.html">
     * TET_MESH_QUALITY Interactive Program for Tet Mesh Quality</a>
     *
     * @param[in] v0 first vertex of the tetrahedron.
     * @param[in] v1 second vertex of the tetrahedron.
     * @param[in] v2 third vertex of the tetrahedron.
     * @param[in] v3 fourth vertex of the tetrahedron.
     * @return 2 * sqrt( 6 ) * the insphere radius divided by the maximum of the
     * tetrhedron edge length.
     */
    double tet_quality_insphere_radius_by_max_edge_length(
        const vec3& v0, const vec3& v1, const vec3& v2, const vec3& v3 )
    {
        double in_radius = tetra_insphere_radius( v0, v1, v2, v3 );
        double edge_length = max_tet_edge_length( v0, v1, v2, v3 );

        return 2 * std::sqrt( 6. ) * in_radius / edge_length;
    }

    /*!
     * @param[in] v0 first vertex of the tetrahedron.
     * @param[in] v1 second vertex of the tetrahedron.
     * @param[in] v2 third vertex of the tetrahedron.
     * @param[in] v3 fourth vertex of the tetrahedron.
     * @return the sum of the square tetrahedron edge length.
     */
    double sum_square_edge_length(
        const vec3& v0, const vec3& v1, const vec3& v2, const vec3& v3 )
    {
        double l1 = vec3( v1 - v0 ).length2();
        double l2 = vec3( v2 - v0 ).length2();
        double l3 = vec3( v3 - v0 ).length2();
        double l4 = vec3( v2 - v1 ).length2();
        double l5 = vec3( v3 - v1 ).length2();
        double l6 = vec3( v3 - v2 ).length2();
        return l1 + l2 + l3 + l4 + l5 + l6;
    }

    /*!
     * @brief Tetrahedron quality based on the tetrahedron volume and
     * the sum of the square edges.
     *
     * The quality is between 0 and 1. 0 corresponds to a bad tetrahedron, and
     * 1 to a good tetrahedron (equilaterality).
     * For more information see
     * <a
     * href="http://people.sc.fsu.edu/~jburkardt/cpp_src/tet_mesh_quality/tet_mesh_quality.html">
     * TET_MESH_QUALITY Interactive Program for Tet Mesh Quality</a>
     *
     * @param[in] v0 first vertex of the tetrahedron.
     * @param[in] v1 second vertex of the tetrahedron.
     * @param[in] v2 third vertex of the tetrahedron.
     * @param[in] v3 fourth vertex of the tetrahedron.
     *
     * @return 12. * (3 * volume)^(2/3) / sum of the square edge.
     */
    double tet_quality_volume_by_sum_square_edges(
        const vec3& v0, const vec3& v1, const vec3& v2, const vec3& v3 )
    {
        const double tet_volume = GEO::Geom::tetra_volume( v0, v1, v2, v3 );
        double sum_square_edge = sum_square_edge_length( v0, v1, v2, v3 );
        ringmesh_assert( sum_square_edge > global_epsilon );
        return 12. * std::pow( 3. * tet_volume, 2. / 3. ) / sum_square_edge;
    }

    /*!
     * @bried Computes the sinus of the half solid angle relatively to a
     * tetrahedron vertex.
     * @param[in] v0 first vertex of the tetrahedron. The solid angle is
     * computed
     * relatively to this vertex.
     * @param[in] v1 second vertex of the tetrahedron.
     * @param[in] v2 third vertex of the tetrahedron.
     * @param[in] v3 fourth vertex of the tetrahedron.
     * @return the sinus of the half solid angle on the vertex \p v0.
     */
    double sin_half_solid_angle(
        const vec3& v0, const vec3& v1, const vec3& v2, const vec3& v3 )
    {
        double tet_volume = GEO::Geom::tetra_volume( v0, v1, v2, v3 );
        double l01 = ( v1 - v0 ).length();
        double l02 = ( v2 - v0 ).length();
        double l03 = ( v3 - v0 ).length();
        double l12 = ( v2 - v1 ).length();
        double l13 = ( v3 - v1 ).length();
        double l23 = ( v3 - v2 ).length();
        double denominator = ( l01 + l02 + l12 ) * ( l01 + l02 - l12 )
                             * ( l02 + l03 + l23 ) * ( l02 + l03 - l23 )
                             * ( l03 + l01 + l13 ) * ( l03 + l01 - l13 );

        ringmesh_assert( denominator > global_epsilon_sq );

        denominator = std::sqrt( denominator );
        return 12 * tet_volume / denominator;
    }

    /*!
     * @brief Tetrahedron quality based on the solid angles.
     *
     * This metrics is based on the sinus of the half solid angle on each
     * vertex. It was shown on the literature that the minimum of the four
     * values provides an estimate of the tetrahedron quality.
     * For more information, see:
     * <a
     * href="http://people.eecs.berkeley.edu/~jrs/meshpapers/robnotes.pdf">robnotes.pdf</a>
     * p15,
     * Liu, A., and B. Joe, 1994, Relationship between tetrahedron shape
     * measures,
     * BIT, v. 34, no. 2, p. 268-287, <a
     * href="http://doi.org/10.1007/BF01955874">doi</a> and
     * <a
     * href="http://people.sc.fsu.edu/~jburkardt/cpp_src/tet_mesh_quality/tet_mesh_quality.html">
     * TET_MESH_QUALITY Interactive Program for Tet Mesh Quality</a>
     *
     * @param[in] v0 first vertex of the tetrahedron.
     * @param[in] v1 second vertex of the tetrahedron.
     * @param[in] v2 third vertex of the tetrahedron.
     * @param[in] v3 fourth vertex of the tetrahedron.
     *
     * @return 1.5 * sqrt( 6 ) * the minimun of the sinus of the half
     * solid angles. 1.5 * sqrt( 6 ) is a factor to scale the metrics between
     * 0 and 1.
     * 0 corresponds to a bad tetrahedron, and
     * 1 to a good tetrahedron (equilaterality).
     */
    double tet_quality_min_solid_angle(
        const vec3& v0, const vec3& v1, const vec3& v2, const vec3& v3 )
    {
        double min_sin_half_solid_angle = std::min(
            std::min( std::min( sin_half_solid_angle( v0, v1, v2, v3 ),
                          sin_half_solid_angle( v1, v0, v2, v3 ) ),
                sin_half_solid_angle( v2, v0, v1, v3 ) ),
            sin_half_solid_angle( v3, v0, v1, v2 ) );
        return 1.5 * std::sqrt( 6. ) * min_sin_half_solid_angle;
    }

    /*!
     * @param[in] mesh_qual_mode mesh quality number.
     * @return the property name associated to the mesh quality number
     * \p mesh_qual_mode.
     */
    std::string mesh_qual_mode_to_prop_name( MeshQualityMode mesh_qual_mode )
    {
        std::string quality_name;
        switch( mesh_qual_mode )
        {
        case INSPHERE_RADIUS_BY_CIRCUMSPHERE_RADIUS:
            quality_name = "INSPHERE_RADIUS_BY_CIRCUMSPHERE_RADIUS";
            break;
        case INSPHERE_RADIUS_BY_MAX_EDGE_LENGTH:
            quality_name = "INSPHERE_RADIUS_BY_MAX_EDGE_LENGTH";
            break;
        case VOLUME_BY_SUM_SQUARE_EDGE:
            quality_name = "VOLUME_BY_SUM_SQUARE_EDGE";
            break;
        case MIN_SOLID_ANGLE:
            quality_name = "MIN_SOLID_ANGLE";
            break;
        default:
            ringmesh_assert_not_reached;
        }
        ringmesh_assert( !quality_name.empty() );
        return quality_name;
    }

    /*!
     * @brief Gets the quality for one tetrahedron.
     *
     * The quality is between 0 and 1. 0 corresponds to a bad tetrahedron, and
     * 1 to a good tetrahedron (equilaterality).
     *
     * @param[in] v0 first vertex of the tetrahedron.
     * @param[in] v1 second vertex of the tetrahedron.
     * @param[in] v2 third vertex of the tetrahedron.
     * @param[in] v3 fourth vertex of the tetrahedron.
     * @param[in] mesh_qual_mode tetrahedron quality to get.
     * @return the tetrahedron quality.
     */
    double get_tet_quality( const vec3& v0,
        const vec3& v1,
        const vec3& v2,
        const vec3& v3,
        MeshQualityMode mesh_qual_mode )
    {
        double quality = -1;
        switch( mesh_qual_mode )
        {
        case INSPHERE_RADIUS_BY_CIRCUMSPHERE_RADIUS:
            quality = tet_quality_insphere_radius_by_circumsphere_radius(
                v0, v1, v2, v3 );
            break;
        case INSPHERE_RADIUS_BY_MAX_EDGE_LENGTH:
            quality = tet_quality_insphere_radius_by_max_edge_length(
                v0, v1, v2, v3 );
            break;
        case VOLUME_BY_SUM_SQUARE_EDGE:
            quality = tet_quality_volume_by_sum_square_edges( v0, v1, v2, v3 );
            break;
        case MIN_SOLID_ANGLE:
            quality = tet_quality_min_solid_angle( v0, v1, v2, v3 );
            break;
        default:
            ringmesh_assert_not_reached;
        }
        ringmesh_assert(
            quality > -1 * global_epsilon && quality < 1 + global_epsilon );
        return quality;
    }
} // namespace

namespace RINGMesh
{
    void compute_prop_tet_mesh_quality(
        MeshQualityMode mesh_qual_mode, const GeoModel3D& geomodel )
    {
        ringmesh_assert( geomodel.nb_regions() != 0 );
        for( const auto& region : geomodel.regions() )
        {
            ringmesh_assert( region.is_meshed() );
            ringmesh_assert( region.is_simplicial() );
            GEO::AttributesManager& reg_attr_mgr =
                region.cell_attribute_manager();
            GEO::Attribute< double > attr(
                reg_attr_mgr, mesh_qual_mode_to_prop_name( mesh_qual_mode ) );
            for( auto cell_itr : range( region.nb_mesh_elements() ) )
            {
                attr[cell_itr] = get_tet_quality(
                    region.mesh_element_vertex( { cell_itr, 0 } ),
                    region.mesh_element_vertex( { cell_itr, 1 } ),
                    region.mesh_element_vertex( { cell_itr, 2 } ),
                    region.mesh_element_vertex( { cell_itr, 3 } ),
                    mesh_qual_mode );
            }
        }
    }

    double fill_mesh_with_low_quality_cells( MeshQualityMode mesh_qual_mode,
        double min_quality,
        const GeoModel3D& geomodel,
        VolumeMesh3D& output_mesh )
    {
        ringmesh_assert( geomodel.nb_regions() != 0 );
        auto mesh_builder = VolumeMeshBuilder3D::create_builder( output_mesh );
        double min_qual_value{ max_float64() };
        for( const auto& region : geomodel.regions() )
        {
            ringmesh_assert( region.is_meshed() );
            ringmesh_assert( region.is_simplicial() );
            GEO::Attribute< double > quality_attribute(
                region.cell_attribute_manager(),
                mesh_qual_mode_to_prop_name( mesh_qual_mode ) );
            for( auto cell_id : range( region.nb_mesh_elements() ) )
            {
                if( quality_attribute[cell_id] < min_quality )
                {
                    auto first_new_vertex_id = output_mesh.nb_vertices();
                    for( auto v : range( 4 ) )
                    {
                        mesh_builder->create_vertex(
                            region.mesh_element_vertex( { cell_id, v } ) );
                    }
                    auto new_cell_id =
                        mesh_builder->create_cells( 1, CellType::TETRAHEDRON );
                    for( auto v_id : range( 4 ) )
                    {
                        mesh_builder->set_cell_vertex(
                            { new_cell_id, v_id }, first_new_vertex_id + v_id );
                    }
                    if( quality_attribute[cell_id] < min_qual_value )
                    {
                        min_qual_value = quality_attribute[cell_id];
                    }
                }
            }
        }
        return min_qual_value;
    }
} // namespace RINGMesh
