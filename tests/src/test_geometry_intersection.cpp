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

#include <ringmesh/ringmesh_tests_config.h>

#include <vector>

#include <geogram/basic/vecg.h>

#include <ringmesh/basic/geometry.h>

/*!
 * @author Pierre Anquez
 */

using namespace RINGMesh ;

bool are_almost_equal( const vec3& vec0, const vec3& vec1 )
{
    return ( vec0 - vec1 ).length2() < global_epsilon_sq ;
}

void test_line_plane_intersection()
{
    Logger::out( "TEST" ) << "Test Line-Plane intersections" << std::endl ;
    vec3 O_plane( 1., 4., -2 ) ;
    vec3 N_plane( 1., 0., 0. ) ;

    // Intersection is the origin point of the line
    vec3 O_line1( 1., 5., 0. ) ;
    vec3 D_line1( -2., 4., 1. ) ;
    vec3 result1 ;
    vec3 answer1 = O_line1 ;
    bool intersect1 = line_plane_intersection( O_line1, D_line1, O_plane, N_plane,
        result1 ) ;
    if( !intersect1 || result1 != answer1 ) {
        throw RINGMeshException( "TEST", "True intersection1: KO" ) ;
    } else {
        Logger::out( "TEST" ) << "True intersection1: OK" << std::endl ;
    }

    // Intersection is a point
    vec3 O_line2( -41., 7., -28. ) ;
    vec3 D_line2( -2., 4., 1. ) ;
    vec3 result2 ;
    vec3 answer2( 1., -77., -49. ) ;
    bool intersect2 = line_plane_intersection( O_line2, D_line2, O_plane, N_plane,
        result2 ) ;
    if( !intersect2 || result2 != answer2 ) {
        throw RINGMeshException( "TEST", "True intersection2: KO" ) ;
    } else {
        Logger::out( "TEST" ) << "True intersection2: OK" << std::endl ;
    }

    // The line is parallel to the plane
    vec3 O_line3( 0., 1., 8. ) ;
    vec3 D_line3( 0., 2., 1. ) ;
    vec3 result3 ;
    bool intersect3 = line_plane_intersection( O_line3, D_line3, O_plane, N_plane,
        result3 ) ;
    if( intersect3 ) {
        throw RINGMeshException( "TEST", "Line parallel to the plane: KO" ) ;
    } else {
        Logger::out( "TEST" ) << "Line parallel to the plane: OK" << std::endl ;
    }

    // The line is included into the plane
    vec3 O_line4( 1., 1., 8. ) ;
    vec3 D_line4( 0., 2., 1. ) ;
    vec3 result4 ;
    bool intersect4 = line_plane_intersection( O_line4, D_line4, O_plane, N_plane,
        result4 ) ;
    if( intersect4 ) {
        throw RINGMeshException( "TEST", "Line included into the plane: KO" ) ;
    } else {
        Logger::out( "TEST" ) << "Line included into the plane: OK" << std::endl ;
    }

    Logger::out( "TEST" ) << " " << std::endl ;
}

void test_segment_plane_intersection()
{
    Logger::out( "TEST" ) << "Test Segment-Plane intersections" << std::endl ;
    vec3 O_plane( 1., -2., 3. ) ;
    vec3 N_plane( 1., 1.25, -2 ) ;

    // Intersection is a point
    vec3 seg10( 3., -1., 6.75 ) ;
    vec3 seg11( -3., -10., -8.25 ) ;
    vec3 result1 ;
    vec3 answer1( 1., -4., 1.75 ) ;
    bool intersect1 = segment_plane_intersection( seg10, seg11, O_plane, N_plane,
        result1 ) ;
    if( !intersect1 || result1 != answer1 ) {
        throw RINGMeshException( "TEST", "True intersection: KO" ) ;
    } else {
        Logger::out( "TEST" ) << "True intersection: OK" << std::endl ;
    }

    // Intersection is an extremal point of the segment
    vec3 seg20( 3., -1., 6.75 ) ;
    vec3 seg21( 1., -4., 1.75 ) ;
    vec3 result2 ;
    vec3 answer2 = seg21 ;
    bool intersect2 = segment_plane_intersection( seg20, seg21, O_plane, N_plane,
        result2 ) ;
    if( !intersect2 || result2 != answer2 ) {
        throw RINGMeshException( "TEST", "Intersection at segment extremity: KO" ) ;
    } else {
        Logger::out( "TEST" ) << "Intersection at segment extremity: OK"
            << std::endl ;
    }

    // Line intersects but not the segment
    vec3 seg30( 3., -1., 6.75 ) ;
    vec3 seg31( 2., -2.5, 4.25 ) ;
    vec3 result3 ;
    bool intersect3 = segment_plane_intersection( seg30, seg31, O_plane, N_plane,
        result3 ) ;
    if( intersect3 ) {
        throw RINGMeshException( "TEST", "Intersection outside segment: KO" ) ;
    } else {
        Logger::out( "TEST" ) << "Intersection outside segment: OK" << std::endl ;
    }

    // The segment is parallel to the plane
    vec3 seg40( 0., 0., 6.75 ) ;
    vec3 seg41( 1., 0., 7.25 ) ;
    vec3 result4 ;
    bool intersect4 = segment_plane_intersection( seg40, seg41, O_plane, N_plane,
        result4 ) ;
    if( intersect4 ) {
        throw RINGMeshException( "TEST", "Segment parallel to plane: KO" ) ;
    } else {
        Logger::out( "TEST" ) << "Segment parallel to plane: OK" << std::endl ;
    }

    // The segment is included into the plane
    vec3 seg50( 5., -2., 5 ) ;
    vec3 seg51( 1., -2., 3 ) ;
    vec3 result5 ;
    bool intersect5 = segment_plane_intersection( seg50, seg51, O_plane, N_plane,
        result5 ) ;
    if( intersect5 ) {
        throw RINGMeshException( "TEST", "Segment included into the plane: KO" ) ;
    } else {
        Logger::out( "TEST" ) << "Segment included into the plane: OK" << std::endl ;
    }

    Logger::out( "TEST" ) << " " << std::endl ;

}

void test_segment_triangle_intersection()
{
    Logger::out( "TEST" ) << "Test Segment-Triangle intersections" << std::endl ;

    vec3 trgl0( 1., 1., 0. ) ;
    vec3 trgl1( 2., 3., 0. ) ;
    vec3 trgl2( 4., -1., 0. ) ;

    // True intersection inside the triangle
    vec3 seg10( 2., 2., 3. ) ;
    vec3 seg11( 2., 2., -3. ) ;
    vec3 result1 ;
    vec3 answer1( 2., 2., 0. ) ;
    bool intersect1 = segment_triangle_intersection( seg10, seg11, trgl0, trgl1,
        trgl2, result1 ) ;
    if( !intersect1 || result1 != answer1 ) {
        throw RINGMeshException( "TEST", "Test1: KO" ) ;
    } else {
        Logger::out( "TEST" ) << "Test1: OK" << std::endl ;
    }

    // No intersection
    vec3 seg20( 2., 2., 3. ) ;
    vec3 seg21( 20., 2., -1. ) ;
    vec3 result2 ;
    vec3 answer2( 2., 2., 0. ) ;
    bool intersect2 = segment_triangle_intersection( seg20, seg21, trgl0, trgl1,
        trgl2, result2 ) ;
    if( intersect2 ) {
        throw RINGMeshException( "TEST", "Test2: KO" ) ;
    } else {
        Logger::out( "TEST" ) << "Test2: OK" << std::endl ;
    }

    // Intersection at a triangle vertex
    vec3 seg30( 1., 4., 3. ) ;
    vec3 seg31( 3., 2., -3. ) ;
    vec3 result3 ;
    vec3 answer3 = trgl1 ;
    bool intersect3 = segment_triangle_intersection( seg30, seg31, trgl0, trgl1,
        trgl2, result3 ) ;
    if( !intersect3 || result3 != answer3 ) {
        throw RINGMeshException( "TEST", "Test3: KO" ) ;
    } else {
        Logger::out( "TEST" ) << "Test3: OK" << std::endl ;
    }

    // Intersection on a triangle edge
    vec3 seg40( 1., -4., 3. ) ;
    vec3 seg41( 4., 4., -3. ) ;
    vec3 result4 ;
    vec3 answer4( 2.5, 0., 0. ) ;
    bool intersect4 = segment_triangle_intersection( seg40, seg41, trgl0, trgl1,
        trgl2, result4 ) ;
    if( !intersect4 || result4 != answer4 ) {
        throw RINGMeshException( "TEST", "Test4: KO" ) ;
    } else {
        Logger::out( "TEST" ) << "Test4: OK" << std::endl ;
    }

    // Segment included inside the triangle
    vec3 seg50( 2., 2., 0. ) ;
    vec3 seg51( 3., 1., 0. ) ;
    vec3 result5 ;
    bool intersect5 = segment_triangle_intersection( seg50, seg51, trgl0, trgl1,
        trgl2, result5 ) ;
    if( intersect5 ) {
        throw RINGMeshException( "TEST", "Test5: KO" ) ;
    } else {
        Logger::out( "TEST" ) << "Test5: OK" << std::endl ;
    }

    // Segment is a triangle edge
    vec3 seg60 = trgl1 ;
    vec3 seg61 = trgl0 ;
    vec3 result6 ;
    bool intersect6 = segment_triangle_intersection( seg60, seg61, trgl0, trgl1,
        trgl2, result6 ) ;
    if( intersect6 ) {
        throw RINGMeshException( "TEST", "Test6: KO" ) ;
    } else {
        Logger::out( "TEST" ) << "Test6: OK" << std::endl ;
    }

    // Segment in the same plane than triangle, one point inside the other outside
    vec3 seg70( 2., 2., 0. ) ;
    vec3 seg71( 4., 1., -0. ) ;
    vec3 result7 ;
    bool intersect7 = segment_triangle_intersection( seg70, seg71, trgl0, trgl1,
        trgl2, result7 ) ;
    if( intersect7 ) {
        throw RINGMeshException( "TEST", "Test7: KO" ) ;
    } else {
        Logger::out( "TEST" ) << "Test7: OK" << std::endl ;
    }

    Logger::out( "TEST" ) << " " << std::endl ;

}

void test_circle_plane_intersection()
{
    Logger::out( "TEST" ) << "Test Circle-Plane intersections" << std::endl ;
    vec3 O_plane( 3., 1., -1. ) ;
    vec3 N_plane( 0., 0., -2. ) ;

    // Circle parallel to the plane
    vec3 O_circle1( 2., 3., 4. ) ;
    vec3 N_circle1( 0., 0., 1. ) ;
    double r1 = 4. ;
    std::vector< vec3 > results1 ;
    bool intersect1 = circle_plane_intersection( O_plane, N_plane, O_circle1,
        N_circle1, r1, results1 ) ;
    if( intersect1 ) {
        throw RINGMeshException( "TEST", "Test circle parallel to plane: KO" ) ;
    } else {
        Logger::out( "TEST" ) << "Test circle parallel to plane: OK" << std::endl ;
    }

    // Circle adjacent to the plane
    vec3 O_circle2( 2., 3., 4. ) ;
    vec3 N_circle2( -1., 2., 0. ) ;
    double r2 = 5. ;
    std::vector< vec3 > results2 ;
    vec3 answer2( 2., 3., -1. ) ;
    bool intersect2 = circle_plane_intersection( O_plane, N_plane, O_circle2,
        N_circle2, r2, results2 ) ;
    if( !intersect2 || results2.size() != 1 ) {
        throw RINGMeshException( "TEST",
            "Test circle adjacent to the plane: KO (wrong number of points)" ) ;
    } else {
        if( !are_almost_equal( results2[0], answer2 ) ) {
            throw RINGMeshException( "TEST",
                "Test circle adjacent to the plane: KO (wrong point coordinates)" ) ;
        }
        Logger::out( "TEST" ) << "Test circle adjacent to the plane: OK"
            << std::endl ;
    }

    // Circle crossing the plane
    vec3 O_circle3( 2., 3., 0. ) ;
    vec3 N_circle3( 1., -1., 0. ) ;
    double r3 = 2. ;
    std::vector< vec3 > results3 ;
    vec3 answer31 = O_circle3
        + vec3( sqrt( 2 ) * cos( M_PI / 6 ), sqrt( 2 ) * cos( M_PI / 6 ), -1. ) ;
    vec3 answer32 = O_circle3
        + vec3( -sqrt( 2 ) * cos( M_PI / 6 ), -sqrt( 2 ) * cos( M_PI / 6 ), -1. ) ;
    bool intersect3 = circle_plane_intersection( O_plane, N_plane, O_circle3,
        N_circle3, r3, results3 ) ;
    if( !intersect3 || results3.size() != 2 ) {
        throw RINGMeshException( "TEST",
            "Test circle crossing the plane: KO (wrong number of points)" ) ;
    } else {
        if( !( are_almost_equal( results3[0], answer31 )
            || ( are_almost_equal( results3[0], answer32 ) ) )
            || !( are_almost_equal( results3[1], answer31 )
                || ( are_almost_equal( results3[1], answer32 ) ) ) ) {

            throw RINGMeshException( "TEST",
                "Test circle crossing the plane: KO (wrong point coordinates)" ) ;
        }
        Logger::out( "TEST" ) << "Test circle crossing the plane: OK" << std::endl ;
    }

    Logger::out( "TEST" ) << " " << std::endl ;
}

void test_disk_segment_intersection()
{
    Logger::out( "TEST" ) << "Test Disk-Segment intersections" << std::endl ;
    vec3 O_disk( 2., 2., 2. ) ;
    vec3 N_disk( 0., 4., 0. ) ;
    double disk_radius = 4. ;

    // Segment in the disk plane
    vec3 seg_10( 1., 2., 3. ) ;
    vec3 seg_11( 3., 2., 1. ) ;
    vec3 result1 ;
    bool intersect1 = disk_segment_intersection( seg_10, seg_11, O_disk, N_disk,
        disk_radius, result1 ) ;
    if( intersect1 ) {
        throw RINGMeshException( "TEST", "Test segment inside disk: KO" ) ;
    } else {
        Logger::out( "TEST" ) << "Test segment inside disk: OK" << std::endl ;
    }

    // Segment adjacent to the disk
    vec3 seg_20( -2., 0., -2. ) ;
    vec3 seg_21( -2., 4., 4. ) ;
    vec3 result2 ;
    bool intersect2 = disk_segment_intersection( seg_20, seg_21, O_disk, N_disk,
        disk_radius, result2 ) ;
    if( intersect2 ) {
        throw RINGMeshException( "TEST", "Test segment tangent to the disk: KO" ) ;
    } else {
        Logger::out( "TEST" ) << "Test segment tangent to the disk: OK"
            << std::endl ;
    }

    // Circle crossing the disk
    vec3 seg_30( 1., 1., 3. ) ;
    vec3 seg_31( 3., 3., 1. ) ;
    vec3 answer3( 2., 2., 2. ) ;
    vec3 result3 ;
    bool intersect3 = disk_segment_intersection( seg_30, seg_31, O_disk, N_disk,
        disk_radius, result3 ) ;
    if( !intersect3 ) {
        throw RINGMeshException( "TEST",
            "Test circle adjacent to the plane: KO (wrong number of points)" ) ;
    } else {
        if( !are_almost_equal( result3, answer3 ) ) {
            throw RINGMeshException( "TEST",
                "Test circle adjacent to the plane: KO (wrong point coordinates)" ) ;
        }
        Logger::out( "TEST" ) << "Test circle adjacent to the plane: OK"
            << std::endl ;
    }

    Logger::out( "TEST" ) << " " << std::endl ;
}

void test_circle_triangle_intersection()
{
    Logger::out( "TEST" ) << "Test Circle-Triangle intersections" << std::endl ;
    vec3 trgl0( 1., 1., 0. ) ;
    vec3 trgl1( 2., 3., 0. ) ;
    vec3 trgl2( 4., -1., 0. ) ;

    // The circle is adjacent to the triangle plane exactly at a triangle vertex
    vec3 O_circle1( 2., 3., 4. ) ;
    vec3 N_circle1( 1., 1., 0. ) ;
    double r1 = 4. ;
    std::vector< vec3 > results1 ;
    vec3 answer1 = trgl1 ;
    bool intersect1 = circle_triangle_intersection( trgl0, trgl1, trgl2, O_circle1,
        N_circle1, r1, results1 ) ;
    if( !intersect1 || results1.size() != 1 ) {
        throw RINGMeshException( "TEST", "Test1: KO (wrong number of points)" ) ;
    } else {
        if( results1[0] != answer1 ) {
            throw RINGMeshException( "TEST",
                "Test1: KO (wrong point coordinates)" ) ;
        }
        Logger::out( "TEST" ) << "Test1: OK" << std::endl ;
    }

    // One point inside triangle the other one outside
    vec3 O_circle2( 2., 2., 0. ) ;
    vec3 N_circle2( 1.5, 1., 0. ) ;
    double r2 = 1. ;
    std::vector< vec3 > results2 ;
    vec3 answer2( 2. + std::sqrt( 2. ), 2. - std::sqrt( 2. ), 0. ) ;
    bool intersect2 = circle_triangle_intersection( trgl0, trgl1, trgl2, O_circle2,
        N_circle2, r2, results2 ) ;
    if( !intersect2 || results2.size() != 1 ) {
        throw RINGMeshException( "TEST", "Test2: KO (wrong number of points)" ) ;
    } else {
        if( results1[0] != answer1 ) {
            throw RINGMeshException( "TEST",
                "Test2: KO (wrong point coordinates)" ) ;
        }
        Logger::out( "TEST" ) << "Test2: OK" << std::endl ;
    }

    Logger::out( "TEST" ) << " " << std::endl ;
}

void test_plane_plane_intersection()
{
    Logger::out( "TEST" ) << "Test Plane-Plane intersections" << std::endl ;

    vec3 O_P0( 4., -2., 0. ) ;
    vec3 N_P0( 1., -2., 4. ) ;

    // Two parallel planes
    vec3 O_P11( 6., 0., 1.52 ) ;
    vec3 N_P11( -2., 4., -8. ) ;
    vec3 O_inter_result1 ;
    vec3 D_inter_result1 ;
    bool intersect1 = plane_plane_intersection( O_P0, N_P0, O_P11, N_P11,
        O_inter_result1, D_inter_result1 ) ;
    if( intersect1 ) {
        throw RINGMeshException( "TEST", "Test parallel planes: KO" ) ;
    } else {
        Logger::out( "TEST" ) << "Test parallel planes: OK" << std::endl ;
    }

    // Two times the same plane
    vec3 O_P12( 4., -2., 0. ) ;
    vec3 N_P12( -1., 2., -4. ) ;
    vec3 O_inter_result2 ;
    vec3 D_inter_result2 ;
    bool intersect2 = plane_plane_intersection( O_P0, N_P0, O_P12, N_P12,
        O_inter_result2, D_inter_result2 ) ;
    if( intersect2 ) {
        throw RINGMeshException( "TEST", "Test same plane: KO" ) ;
    } else {
        Logger::out( "TEST" ) << "Test same plane: OK" << std::endl ;
    }

    // Two intersecting plane
    vec3 O_P13( -3., -1., 1. ) ;
    vec3 N_P13( 0., 2., -4. ) ;
    vec3 O_inter_result3 ;
    vec3 D_inter_result3 ;
    vec3 O_inter_answer3( 2., -3., 0. ) ;
    vec3 D_inter_answer3( 0., 2., 1. ) ;
    bool intersect3 = plane_plane_intersection( O_P0, N_P0, O_P13, N_P13,
        O_inter_result3, D_inter_result3 ) ;
    if( !intersect3
        || !are_almost_equal( normalize( D_inter_answer3 ),
            normalize( D_inter_result3 ) )
        || !are_almost_equal( normalize( D_inter_answer3 ),
            normalize( O_inter_result3 - O_inter_answer3 ) ) ) {
        throw RINGMeshException( "TEST", "Test intersecting planes: KO" ) ;
    } else {
        Logger::out( "TEST" ) << "Test intersecting planes: OK" << std::endl ;
    }

    Logger::out( "TEST" ) << " " << std::endl ;
}

int main()
{
    try {
        default_configure() ;

        Logger::out( "TEST" ) << "Test intersection algorithms" << std::endl ;

        test_line_plane_intersection() ;
        test_segment_plane_intersection() ;
        test_segment_triangle_intersection() ;
        test_circle_plane_intersection() ;
        test_disk_segment_intersection() ;
        test_circle_triangle_intersection() ;
        test_plane_plane_intersection() ;

    } catch( const RINGMeshException& e ) {
        Logger::err( e.category() ) << e.what() << std::endl ;
        return 1 ;
    } catch( const std::exception& e ) {
        Logger::err( "Exception" ) << e.what() << std::endl ;
        return 1 ;
    }
    Logger::out( "TEST" ) << "SUCCESS" << std::endl ;
    return 0 ;
}
