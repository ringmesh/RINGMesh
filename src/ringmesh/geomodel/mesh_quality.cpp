#include <ringmesh/geomodel/mesh_quality.h>

#include <ringmesh/geomodel/geo_model.h>
#include <ringmesh/geomodel/geo_model_mesh_entity.h>

/*!
 * @author Benjamin Chauvin
 */

namespace {

    using namespace RINGMesh ;

    double tetra_in_sphere_radius(
        const vec3& v0,
        const vec3& v1,
        const vec3& v2,
        const vec3& v3 )
    {
        double tet_volume = GEO::Geom::tetra_volume( v0, v1, v2, v3 ) ;
        double A1 = GEO::Geom::triangle_area( v0, v1, v2 ) ;
        double A2 = GEO::Geom::triangle_area( v1, v2, v3 ) ;
        double A3 = GEO::Geom::triangle_area( v2, v3, v0 ) ;
        double A4 = GEO::Geom::triangle_area( v3, v0, v1 ) ;
        ringmesh_assert( A1 + A2 + A3 + A4 > global_epsilon ) ;
        return ( 3 * tet_volume ) / ( A1 + A2 + A3 + A4 ) ;
    }

    double tet_quality_insphere_radius_by_circumsphere_radius(
        const vec3& v0,
        const vec3& v1,
        const vec3& v2,
        const vec3& v3 )
    {
        const vec3 tetra_circum_center = GEO::Geom::tetra_circum_center( v0, v1, v2,
            v3 ) ;
        const double tetra_circum_radius =
            vec3( tetra_circum_center - v0 ).length() ;
        ringmesh_assert(
            std::abs( tetra_circum_radius - ( tetra_circum_center - v1 ).length() )
            < global_epsilon ) ;
        ringmesh_assert(
            std::abs( tetra_circum_radius - ( tetra_circum_center - v2 ).length() )
            < global_epsilon ) ;
        ringmesh_assert(
            std::abs( tetra_circum_radius - ( tetra_circum_center - v3 ).length() )
            < global_epsilon ) ;

        // insphere computation
        double in_radius = tetra_in_sphere_radius( v0, v1, v2, v3 ) ;
        return 3. * in_radius / tetra_circum_radius ;
    }

    double max_tet_edge_length(
        const vec3& v0,
        const vec3& v1,
        const vec3& v2,
        const vec3& v3 )
    {
        double l1 = ( v1 - v0 ).length() ;
        double l2 = ( v2 - v0 ).length() ;
        double l3 = ( v3 - v0 ).length() ;
        double l4 = ( v2 - v1 ).length() ;
        double l5 = ( v3 - v1 ).length() ;
        double l6 = ( v3 - v2 ).length() ;
        return GEO::geo_max(
            GEO::geo_max(
                GEO::geo_max( GEO::geo_max( GEO::geo_max( l1, l2 ), l3 ), l4 ), l5 ),
            l6 ) ;

    }

    double tet_quality_insphere_radius_by_max_edge_length(
        const vec3& v0,
        const vec3& v1,
        const vec3& v2,
        const vec3& v3 )
    {
        double in_radius = tetra_in_sphere_radius( v0, v1, v2, v3 ) ;
        double edge_length = max_tet_edge_length( v0, v1, v2, v3 ) ;

        return 2 * sqrt( 6 ) * in_radius / edge_length ;
    }

    double sum_square_edge_length(
        const vec3& v0,
        const vec3& v1,
        const vec3& v2,
        const vec3& v3 )
    {
        double l1 = vec3( v1 - v0 ).length2() ;
        double l2 = vec3( v2 - v0 ).length2() ;
        double l3 = vec3( v3 - v0 ).length2() ;
        double l4 = vec3( v2 - v1 ).length2() ;
        double l5 = vec3( v3 - v1 ).length2() ;
        double l6 = vec3( v3 - v2 ).length2() ;
        return l1 + l2 + l3 + l4 + l5 + l6 ;

    }

    double tet_quality_volume_sum_square_edges(
        const vec3& v0,
        const vec3& v1,
        const vec3& v2,
        const vec3& v3 )
    {
        const double tet_volume = GEO::Geom::tetra_volume( v0, v1, v2, v3 ) ;
        double sum_square_edge = sum_square_edge_length( v0, v1, v2, v3 ) ;
        ringmesh_assert( sum_square_edge > global_epsilon ) ;
        return 12. * std::pow( 3. * tet_volume, 2. / 3. ) / sum_square_edge ;
    }

    double sin_half_solid_angle(
        const vec3& v0,
        const vec3& v1,
        const vec3& v2,
        const vec3& v3 )
    {
        double tet_volume = GEO::Geom::tetra_volume( v0, v1, v2, v3 ) ;
        double l01 = ( v1 - v0 ).length() ;
        double l02 = ( v2 - v0 ).length() ;
        double l03 = ( v3 - v0 ).length() ;
        double l12 = ( v2 - v1 ).length() ;
        double l13 = ( v3 - v1 ).length() ;
        double l23 = ( v3 - v2 ).length() ;
        double denominator = ( l01 + l02 + l12 ) * ( l01 + l02 - l12 )
            * ( l02 + l03 + l23 ) * ( l02 + l03 - l23 ) * ( l03 + l01 + l13 )
            * ( l03 + l01 - l13 ) ;

        ringmesh_assert( denominator > global_epsilon_sq ) ;

        denominator = sqrt( denominator ) ;
        return 12 * tet_volume / denominator ;
    }

    double tet_quality_min_solid_angle(
        const vec3& v0,
        const vec3& v1,
        const vec3& v2,
        const vec3& v3 )
    {
        double min_sin_half_solid_angle = std::min(
            std::min(
                std::min( sin_half_solid_angle( v0, v1, v2, v3 ),
                    sin_half_solid_angle( v1, v0, v2, v3 ) ),
                sin_half_solid_angle( v2, v0, v1, v3 ) ),
            sin_half_solid_angle( v3, v0, v1, v2 ) ) ;
        return 1.5 * sqrt( 6 ) * min_sin_half_solid_angle ;
    }

    std::string mesh_qual_mode_to_prop_name( MeshQualityMode mesh_qual_mode )
    {
        std::string quality_name = "" ;
        switch( mesh_qual_mode ) {
            case INSPHERE_RADIUS_BY_CIRCUMSPHERE_RADIUS:
                quality_name = "INSPHERE_RADIUS_BY_CIRCUMSPHERE_RADIUS" ;
                break ;
            case INSPHERE_RADIUS_BY_MAX_EDGE_LENGTH:
                quality_name = "INSPHERE_RADIUS_BY_MAX_EDGE_LENGTH" ;
                break ;
            case VOLUME_BY_SUM_SQUARE_EDGE:
                quality_name = "VOLUME_BY_SUM_SQUARE_EDGE" ;
                break ;
            case MIN_SOLID_ANGLE:
                quality_name = "MIN_SOLID_ANGLE" ;
                break ;
            default:
                ringmesh_assert_not_reached ;
        }
        ringmesh_assert( !quality_name.empty() ) ;
        return quality_name ;
    }

    double get_tet_quality(
        const vec3& v0,
        const vec3& v1,
        const vec3& v2,
        const vec3& v3,
        MeshQualityMode mesh_qual_mode )
    {
        double quality = -1 ;
        switch( mesh_qual_mode ) {
            case INSPHERE_RADIUS_BY_CIRCUMSPHERE_RADIUS:
                quality = tet_quality_insphere_radius_by_circumsphere_radius( v0, v1,
                    v2, v3 ) ;
                break ;
            case INSPHERE_RADIUS_BY_MAX_EDGE_LENGTH:
                quality = tet_quality_insphere_radius_by_max_edge_length( v0, v1, v2,
                    v3 ) ;
                break ;
            case VOLUME_BY_SUM_SQUARE_EDGE:
                quality = tet_quality_volume_sum_square_edges( v0, v1, v2, v3 ) ;
                break ;
            case MIN_SOLID_ANGLE:
                quality = tet_quality_min_solid_angle( v0, v1, v2, v3 ) ;
                break ;
            default:
                ringmesh_assert_not_reached ;
        }
        ringmesh_assert( quality > -1 * global_epsilon
            && quality < 1 + global_epsilon ) ;
        return quality ;
    }
}
namespace RINGMesh {

    void compute_prop_tet_mesh_quality(
        MeshQualityMode mesh_qual_mode,
        const GeoModel& geo_model )
    {
        ringmesh_assert( geo_model.nb_regions() != 0 ) ;
        for( index_t reg_itr = 0; reg_itr < geo_model.nb_regions(); ++reg_itr ) {
            const Region& cur_reg = geo_model.region( reg_itr ) ;
            ringmesh_assert( cur_reg.is_meshed() ) ;
            ringmesh_assert( cur_reg.is_simplicial() ) ;
            GEO::AttributesManager& reg_attr_mgr = cur_reg.cell_attribute_manager() ;
            GEO::Attribute< double > attr( reg_attr_mgr,
                mesh_qual_mode_to_prop_name( mesh_qual_mode ) ) ;
            for( index_t cell_itr = 0; cell_itr < cur_reg.nb_mesh_elements();
                ++cell_itr ) {
                attr[cell_itr] = get_tet_quality(
                    cur_reg.mesh_element_vertex( cell_itr, 0 ),
                    cur_reg.mesh_element_vertex( cell_itr, 1 ),
                    cur_reg.mesh_element_vertex( cell_itr, 2 ),
                    cur_reg.mesh_element_vertex( cell_itr, 3 ), mesh_qual_mode ) ;
            }
        }
    }

}
