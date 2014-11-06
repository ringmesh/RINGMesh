/*[
 * Association Scientifique pour la Geologie et ses Applications (ASGA)
 * Copyright (c) 1993-2013 ASGA. All Rights Reserved.
 *
 * This program is a Trade Secret of the ASGA and it is not to be:
 * - reproduced, published, or disclosed to other,
 * - distributed or displayed,
 * - used for purposes or on Sites other than described
 *   in the GOCAD Advancement Agreement,
 * without the prior written authorization of the ASGA. Licencee
 * agrees to attach or embed this Notice on all copies of the program,
 * including partial copies or modified versions thereof.
 ]*/

#include <grgmeshlib/mixed_mesh.h>
#include <grgmeshlib/utils.h>

namespace GRGMesh {

    const CellDescriptor* MixedMesh::cell_descriptor_[7] = {
        &line_descriptor, &trgl_descriptor, &quad_descriptor, &tetra_descriptor,
        &pyramid_descriptor, &prism_descriptor, &hexa_descriptor } ;

    ElementType MixedMesh::cell_type( uint32 c, uint32& c_index ) const
    {
        uint8 result = 0 ;
        uint32 size = 0 ;
        for( ; result < 7; result++ ) {
            size += cells_[result].size() ;
            if( c > size ) break ;
        }
        c_index = c - size ;
        return ElementType( result ) ;
    }

    uint32 MixedMesh::global_index(
        const ElementType& type,
        const uint32 index ) const
    {
        uint32 result = 0 ;
        for( uint8 i = 0; i < type; i++ ) {
            result += cells_[i].size() ;
        }
        result += index ;
        return result ;
    }

    float64 MixedMesh::get_nearest_point_in_cell(
        const vec3& p,
        uint32 c,
        vec3& nearest_p ) const
    {
        uint32 c_index ;
        switch( cell_type( c, c_index ) ) {
            case LINE:
                return get_nearest_point_in_line( p, c_index, nearest_p ) ;
            case TRGL:
                return get_nearest_point_in_triangle( p, c_index, nearest_p ) ;
            case QUAD:
                return get_nearest_point_in_quad( p, c_index, nearest_p ) ;
            case TETRA:
                return get_nearest_point_in_tetra( p, c_index, nearest_p ) ;
            case PYRAMID:
                return get_nearest_point_in_pyramid( p, c_index, nearest_p ) ;
            case PRISM:
                return get_nearest_point_in_prism( p, c_index, nearest_p ) ;
            case HEXA:
                return get_nearest_point_in_hexa( p, c_index, nearest_p ) ;
            default:
                grgmesh_assert_not_reached;
                return 0 ;
            }
        }

    float64 MixedMesh::get_nearest_point_in_line(
        const vec3& p,
        uint32 c,
        vec3& nearest_p ) const
    {
        return Utils::nearest_point_segment( p, line_vertex( c, 0 ),
            line_vertex( c, 1 ), nearest_p ) ;
    }

    float64 MixedMesh::get_nearest_point_in_triangle(
        const vec3& p,
        uint32 c,
        vec3& nearest_p ) const
    {
        const vec3& p1 = triangle_vertex( c, 0 ) ;
        const vec3& p2 = triangle_vertex( c, 1 ) ;
        const vec3& p3 = triangle_vertex( c, 2 ) ;
        float64 distance = Utils::point_triangle_distance( p, p1, p2, p3,
            nearest_p ) ;
        if( Utils::point_inside_triangle( p, p1, p2, p3 ) ) {
            nearest_p = p ;
            distance = -distance ;
        }

        return distance ;
    }

    float64 MixedMesh::get_nearest_point_in_quad(
        const vec3& p,
        uint32 c,
        vec3& nearest_p ) const
    {
        const vec3& p1 = quad_vertex( c, 0 ) ;
        const vec3& p2 = quad_vertex( c, 1 ) ;
        const vec3& p3 = quad_vertex( c, 2 ) ;
        const vec3& p4 = quad_vertex( c, 3 ) ;
        float64 distance = Utils::point_quad_distance( p, p1, p2, p3, p4,
            nearest_p ) ;
        if( Utils::point_inside_quad( p, p1, p2, p3, p4 ) ) {
            nearest_p = p ;
            distance = -distance ;
        }

        return distance ;
    }

    float64 MixedMesh::get_nearest_point_in_tetra(
        const vec3& p,
        uint32 c,
        vec3& nearest_p ) const
    {
        const vec3& p1 = tetra_vertex( c, 0 ) ;
        const vec3& p2 = tetra_vertex( c, 1 ) ;
        const vec3& p3 = tetra_vertex( c, 2 ) ;
        const vec3& p4 = tetra_vertex( c, 3 ) ;
        float64 distance = Utils::point_tetra_distance( p, p1, p2, p3, p4,
            nearest_p ) ;
        if( Utils::point_inside_tetra( p, p1, p2, p3, p4 ) ) {
            nearest_p = p ;
            distance = -distance ;
        }

        return distance ;
    }

    float64 MixedMesh::get_nearest_point_in_pyramid(
        const vec3& p,
        uint32 c,
        vec3& nearest_p ) const
    {
        const vec3& p1 = pyramid_vertex( c, 0 ) ;
        const vec3& p2 = pyramid_vertex( c, 1 ) ;
        const vec3& p3 = pyramid_vertex( c, 2 ) ;
        const vec3& p4 = pyramid_vertex( c, 3 ) ;
        const vec3& p5 = pyramid_vertex( c, 4 ) ;
        float64 distance = Utils::point_pyramid_distance( p, p1, p2, p3, p4, p5,
            nearest_p ) ;
        if( Utils::point_inside_pyramid( p, p1, p2, p3, p4, p5 ) ) {
            nearest_p = p ;
            distance = -distance ;
        }

        return distance ;
    }

    float64 MixedMesh::get_nearest_point_in_prism(
        const vec3& p,
        uint32 c,
        vec3& nearest_p ) const
    {
        const vec3& p1 = prism_vertex( c, 0 ) ;
        const vec3& p2 = prism_vertex( c, 1 ) ;
        const vec3& p3 = prism_vertex( c, 2 ) ;
        const vec3& p4 = prism_vertex( c, 3 ) ;
        const vec3& p5 = prism_vertex( c, 4 ) ;
        const vec3& p6 = prism_vertex( c, 5 ) ;
        float64 distance = Utils::point_prism_distance( p, p1, p2, p3, p4, p5, p6,
            nearest_p ) ;
        if( Utils::point_inside_prism( p, p1, p2, p3, p4, p5, p6 ) ) {
            nearest_p = p ;
            distance = -distance ;
        }

        return distance ;
    }

    float64 MixedMesh::get_nearest_point_in_hexa(
        const vec3& p,
        uint32 c,
        vec3& nearest_p ) const
    {
        const vec3& p1 = hexa_vertex( c, 0 ) ;
        const vec3& p2 = hexa_vertex( c, 1 ) ;
        const vec3& p3 = hexa_vertex( c, 2 ) ;
        const vec3& p4 = hexa_vertex( c, 3 ) ;
        const vec3& p5 = hexa_vertex( c, 4 ) ;
        const vec3& p6 = hexa_vertex( c, 5 ) ;
        const vec3& p7 = hexa_vertex( c, 6 ) ;
        const vec3& p8 = hexa_vertex( c, 7 ) ;
        float64 distance = Utils::point_hexa_distance( p, p1, p2, p3, p4, p5, p6, p7,
            p8, nearest_p ) ;
        if( Utils::point_inside_hexa( p, p1, p2, p3, p4, p5, p6, p7, p8 ) ) {
            nearest_p = p ;
            distance = -distance ;
        }

        return distance ;
    }

}

