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

using namespace RINGMesh;

void verdict( bool condition, std::string test_name )
{
    if( !condition ) {
        throw RINGMeshException( "TEST", test_name, ": KO" );
    } else {
        Logger::out( "TEST", test_name + ": OK" );
    }
}

template< index_t DIMENSION >
bool are_almost_equal( const vecn< DIMENSION >& vec0, const vecn< DIMENSION >& vec1 )
{
    return ( vec0 - vec1 ).length2() < global_epsilon_sq;
}

void test_line_plane_intersection()
{
    Logger::out( "TEST", "Test Line-Plane intersections" );
    vec3 O_plane( 1., 4., -2 );
    vec3 N_plane( 1., 0., 0. );

    // Intersection is the origin point of the line
    vec3 O_line1( 1., 5., 0. );
    vec3 D_line1( -2., 4., 1. );
    vec3 result1;
    vec3 answer1 = O_line1;
    verdict(
        Intersection::line_plane( O_line1, D_line1, O_plane, N_plane, result1 )
            && result1 == answer1, "True intersection1" );

    // Intersection is a point
    vec3 O_line2( -41., 7., -28. );
    vec3 D_line2( -2., 4., 1. );
    vec3 result2;
    vec3 answer2( 1., -77., -49. );
    verdict(
        Intersection::line_plane( O_line2, D_line2, O_plane, N_plane, result2 )
            && result2 == answer2, "True intersection2" );

    // The line is parallel to the plane
    vec3 O_line3( 0., 1., 8. );
    vec3 D_line3( 0., 2., 1. );
    vec3 result3;
    verdict(
        !Intersection::line_plane( O_line3, D_line3, O_plane, N_plane, result3 ),
        "Line parallel to the plane" );

    // The line is included into the plane
    vec3 O_line4( 1., 1., 8. );
    vec3 D_line4( 0., 2., 1. );
    vec3 result4;
    verdict(
        !Intersection::line_plane( O_line4, D_line4, O_plane, N_plane, result4 ),
        "Line included into the plane" );

    Logger::out( "TEST", " " );
}

void test_segment_plane_intersection()
{
    Logger::out( "TEST", "Test Segment-Plane intersections" );
    vec3 O_plane( 1., -2., 3. );
    vec3 N_plane( 1., 1.25, -2 );

    // Intersection is a point
    vec3 seg10( 3., -1., 6.75 );
    vec3 seg11( -3., -10., -8.25 );
    vec3 result1;
    vec3 answer1( 1., -4., 1.75 );
    verdict(
        Intersection::segment_plane( seg10, seg11, O_plane, N_plane, result1 )
            && result1 == answer1, "True intersection" );

    // Intersection is an extremal point of the segment
    vec3 seg20( 3., -1., 6.75 );
    vec3 seg21( 1., -4., 1.75 );
    vec3 result2;
    vec3 answer2 = seg21;
    verdict(
        Intersection::segment_plane( seg20, seg21, O_plane, N_plane, result2 )
            && result2 == answer2, "Intersection at segment extremity" );

    // Line intersects but not the segment
    vec3 seg30( 3., -1., 6.75 );
    vec3 seg31( 2., -2.5, 4.25 );
    vec3 result3;
    verdict( !Intersection::segment_plane( seg30, seg31, O_plane, N_plane, result3 ),
        "Intersection on line but but into segment" );

    // The segment is parallel to the plane
    vec3 seg40( 0., 0., 6.75 );
    vec3 seg41( 1., 0., 7.25 );
    vec3 result4;
    verdict( !Intersection::segment_plane( seg40, seg41, O_plane, N_plane, result4 ),
        "Segment parallel to plane" );

    // The segment is included into the plane
    vec3 seg50( 5., -2., 5 );
    vec3 seg51( 1., -2., 3 );
    vec3 result5;
    verdict( !Intersection::segment_plane( seg50, seg51, O_plane, N_plane, result5 ),
        "Segment included into the plane" );

    Logger::out( "TEST", " " );

}

void test_segment_triangle_intersection()
{
    Logger::out( "TEST", "Test Segment-Triangle intersections" );

    vec3 trgl0( 1., 1., 0. );
    vec3 trgl1( 2., 3., 0. );
    vec3 trgl2( 4., -1., 0. );

    // True intersection inside the triangle
    vec3 seg10( 2., 2., 3. );
    vec3 seg11( 2., 2., -3. );
    vec3 result1;
    vec3 answer1( 2., 2., 0. );
    verdict(
        Intersection::segment_triangle( seg10, seg11, trgl0, trgl1, trgl2, result1 )
            && result1 == answer1, "Test1" );

    // No intersection
    vec3 seg20( 2., 2., 3. );
    vec3 seg21( 20., 2., -1. );
    vec3 result2;
    verdict(
        !Intersection::segment_triangle( seg20, seg21, trgl0, trgl1, trgl2,
            result2 ), "Test2" );

    // Intersection at a triangle vertex
    vec3 seg30( 1., 4., 3. );
    vec3 seg31( 3., 2., -3. );
    vec3 result3;
    vec3 answer3 = trgl1;
    verdict(
        Intersection::segment_triangle( seg30, seg31, trgl0, trgl1, trgl2, result3 )
            && result3 == answer3, "Test3" );

    // Intersection on a triangle edge
    vec3 seg40( 1., -4., 3. );
    vec3 seg41( 4., 4., -3. );
    vec3 result4;
    vec3 answer4( 2.5, 0., 0. );
    verdict(
        Intersection::segment_triangle( seg40, seg41, trgl0, trgl1, trgl2, result4 )
            && result4 == answer4, "Test4" );

    // Segment included inside the triangle
    vec3 seg50( 2., 2., 0. );
    vec3 seg51( 3., 1., 0. );
    vec3 result5;
    verdict(
        !Intersection::segment_triangle( seg50, seg51, trgl0, trgl1, trgl2,
            result5 ), "Test5" );

    // Segment is a triangle edge
    vec3 seg60 = trgl1;
    vec3 seg61 = trgl0;
    vec3 result6;
    verdict(
        !Intersection::segment_triangle( seg60, seg61, trgl0, trgl1, trgl2,
            result6 ), "Test6" );

    // Segment in the same plane than triangle, one point inside the other outside
    vec3 seg70( 2., 2., 0. );
    vec3 seg71( 4., 1., -0. );
    vec3 result7;
    verdict(
        !Intersection::segment_triangle( seg70, seg71, trgl0, trgl1, trgl2,
            result7 ), "Test7" );

    Logger::out( "TEST", " " );

}

void test_circle_plane_intersection()
{
    Logger::out( "TEST", "Test Circle-Plane intersections" );
    vec3 O_plane( 3., 1., -1. );
    vec3 N_plane( 0., 0., -2. );

    // Circle parallel to the plane
    vec3 O_circle1( 2., 3., 4. );
    vec3 N_circle1( 0., 0., 1. );
    double r1 = 4.;
    std::vector< vec3 > results1;
    verdict(
        !Intersection::circle_plane( O_plane, N_plane, O_circle1, N_circle1, r1,
            results1 ), "Test circle parallel to plane" );

    // Circle adjacent to the plane
    vec3 O_circle2( 2., 3., 4. );
    vec3 N_circle2( -1., 2., 0. );
    double r2 = 5.;
    std::vector< vec3 > results2;
    vec3 answer2( 2., 3., -1. );
    verdict(
        Intersection::circle_plane( O_plane, N_plane, O_circle2, N_circle2, r2,
            results2 ) && results2.size() == 1,
        "Test circle adjacent to the plane (number of points)" );
    verdict( are_almost_equal( results2[0], answer2 ),
        "Test circle adjacent to the plane (intersection coordinates)" );

    // Circle crossing the plane
    vec3 O_circle3( 2., 3., 0. );
    vec3 N_circle3( 1., -1., 0. );
    double r3 = 2.;
    std::vector< vec3 > results3;
    vec3 answer31 = O_circle3
        + vec3( std::sqrt( 2 ) * std::cos( M_PI / 6 ),
            std::sqrt( 2 ) * std::cos( M_PI / 6 ), -1. );
    vec3 answer32 = O_circle3
        + vec3( -std::sqrt( 2 ) * std::cos( M_PI / 6 ),
            -std::sqrt( 2 ) * std::cos( M_PI / 6 ), -1. );
    verdict(
        Intersection::circle_plane( O_plane, N_plane, O_circle3, N_circle3, r3,
            results3 ) && results3.size() == 2,
        "Test circle crossing the plane (number of points)" );
    verdict(
        ( are_almost_equal( results3[0], answer31 )
            || ( are_almost_equal( results3[0], answer32 ) ) )
            && ( are_almost_equal( results3[1], answer31 )
                || ( are_almost_equal( results3[1], answer32 ) ) ),
        "Test circle crossing the plane (intersection coordinates)" );

    Logger::out( "TEST", " " );
}

void test_disk_segment_intersection()
{
    Logger::out( "TEST", "Test Disk-Segment intersections" );
    vec3 O_disk( 2., 2., 2. );
    vec3 N_disk( 0., 4., 0. );
    double disk_radius = 4.;

    // Segment in the disk plane
    vec3 seg_10( 1., 2., 3. );
    vec3 seg_11( 3., 2., 1. );
    vec3 result1;
    verdict(
        !Intersection::disk_segment( seg_10, seg_11, O_disk, N_disk, disk_radius,
            result1 ), "Test segment inside disk" );

    // Segment adjacent to the disk
    vec3 seg_20( -2., 0., -2. );
    vec3 seg_21( -2., 4., 4. );
    vec3 result2;
    verdict(
        !Intersection::disk_segment( seg_20, seg_21, O_disk, N_disk, disk_radius,
            result2 ), "Test segment tangent to the disk" );

    // Circle crossing the disk
    vec3 seg_30( 1., 1., 3. );
    vec3 seg_31( 3., 3., 1. );
    vec3 answer3( 2., 2., 2. );
    vec3 result3;
    verdict(
        Intersection::disk_segment( seg_30, seg_31, O_disk, N_disk, disk_radius,
            result3 ), "Test circle adjacent to the plane (intersection exists)" );
    verdict( are_almost_equal( result3, answer3 ),
        "Test circle adjacent to the plane (intersection coordinates)" );

    Logger::out( "TEST", " " );
}

void test_circle_triangle_intersection()
{
    Logger::out( "TEST", "Test Circle-Triangle intersections" );
    vec3 trgl0( 1., 1., 0. );
    vec3 trgl1( 2., 3., 0. );
    vec3 trgl2( 4., -1., 0. );

    // The circle is adjacent to the triangle plane on a triangle border
    vec3 O_circle1( 1.5, 2., 4. );
    vec3 N_circle1( 1., 1., 0. );
    double r1 = 4.;
    std::vector< vec3 > results1;
    vec3 answer1( 1.5, 2., 0. );
    verdict(
        Intersection::circle_triangle( trgl0, trgl1, trgl2, O_circle1, N_circle1, r1,
            results1 ) && results1.size() == 1, "Test1 (number of points)" );
    verdict( are_almost_equal( results1[0], answer1 ),
        "Test1 (intersection coordinates)" );

    // One point inside triangle the other one outside
    vec3 O_circle2( 2., 2., 0. );
    vec3 N_circle2( 1., 1., 0. );
    double r2 = 1.;
    std::vector< vec3 > results2;
    vec3 answer2( 2. + 0.5 * std::sqrt( 2. ), 2. - 0.5 * std::sqrt( 2. ), 0. );
    verdict(
        Intersection::circle_triangle( trgl0, trgl1, trgl2, O_circle2, N_circle2, r2,
            results2 ) && results2.size() == 1, "Test2 (number of points)" );
    verdict( are_almost_equal( results2[0], answer2 ),
        "Test2 (intersection coordinates)" );

    // The circle is adjacent to the triangle plane on a triangle vertex
    vec3 O_circle3( 1., 1., 4. );
    vec3 N_circle3( 1., 1., 0. );
    double r3 = 4.;
    std::vector< vec3 > results3;
    vec3 answer3 = trgl0;
    verdict(
        Intersection::circle_triangle( trgl0, trgl1, trgl2, O_circle3, N_circle3, r3,
            results3 ) && results3.size() == 1, "Test3 (number of points)" );
    verdict( are_almost_equal( results3[0], answer3 ),
        "Test3 (intersection coordinates)" );

    Logger::out( "TEST", " " );
}

void test_plane_plane_intersection()
{
    Logger::out( "TEST", "Test Plane-Plane intersections" );

    vec3 O_P0( 4., -2., 0. );
    vec3 N_P0( 1., -2., 4. );

    // Two parallel planes
    vec3 O_P11( 6., 0., 1.52 );
    vec3 N_P11( -2., 4., -8. );
    vec3 O_inter_result1;
    vec3 D_inter_result1;
    verdict(
        !Intersection::plane_plane( O_P0, N_P0, O_P11, N_P11, O_inter_result1,
            D_inter_result1 ), "Test parallel planes" );

    // Two times the same plane
    vec3 O_P12( 4., -2., 0. );
    vec3 N_P12( -1., 2., -4. );
    vec3 O_inter_result2;
    vec3 D_inter_result2;
    verdict(
        !Intersection::plane_plane( O_P0, N_P0, O_P12, N_P12, O_inter_result2,
            D_inter_result2 ), "Test same plane" );

    // Two intersecting plane
    vec3 O_P13( -3., -1., 1. );
    vec3 N_P13( 0., 2., -4. );
    vec3 O_inter_result3;
    vec3 D_inter_result3;
    vec3 O_inter_answer3( 2., -3., 0. );
    vec3 D_inter_answer3( 0., 2., 1. );
    verdict(
        Intersection::plane_plane( O_P0, N_P0, O_P13, N_P13, O_inter_result3,
            D_inter_result3 )
            && are_almost_equal( normalize( D_inter_answer3 ),
                normalize( D_inter_result3 ) )
            && are_almost_equal( normalize( D_inter_answer3 ),
                normalize( O_inter_result3 - O_inter_answer3 ) ),
        "Test intersecting planes" );

    Logger::out( "TEST", " " );
}

void test_line_line_intersection()
{
    Logger::out( "TEST", "Test Line-Line intersections" );

    // Two parallel lines
    vec2 O_L0_parallel( 0., 0. );
    vec2 D_L0_parallel( 1.5, 1.5 );
    vec2 O_L1_parallel( 1., 1. );
    vec2 result_parallel;
    verdict(
        !Intersection::line_line( O_L0_parallel, D_L0_parallel, O_L1_parallel,
            D_L0_parallel, result_parallel ), "Test parallel lines" );

    // Two times the same line
    vec2 O_L0_same( 0., 0. );
    vec2 D_L0_same( 1.5, 1.5 );
    vec2 D_L1_same( -2.5, -2.5 );
    vec2 result_same;
    verdict(
        !Intersection::line_line( O_L0_same, D_L0_same, O_L0_same, D_L1_same,
            result_same ), "Test same line" );

    // Two intersecting lines
    vec2 O_L0_inter( 0., 0. );
    vec2 D_L0_inter( 1.5, 1.5 );
    vec2 O_L1_inter( 2., 0. );
    vec2 D_L1_inter( 2.5, -2.5 );
    vec2 result_inter;
    vec2 result_answer( 1., 1. );
    verdict(
        Intersection::line_line( O_L0_inter, D_L0_inter, O_L1_inter, D_L1_inter,
            result_inter ) && are_almost_equal( result_inter, result_answer ),
        "Test intersecting lines" );

    Logger::out( "TEST", " " );
}

void test_segment_segment_intersection()
{
    Logger::out( "TEST", "Test Segment-Segment intersections" );

    // Two non-intersecting segments
    vec2 p0_seg0( 0., 0. );
    vec2 p1_seg0( 1.5, 1.5 );
    vec2 p0_seg1( 2., 2. );
    vec2 p1_seg1( 3., 2. );
    vec2 no_result;
    verdict(
        !Intersection::segment_segment( p0_seg0, p1_seg0, p0_seg1, p1_seg1,
            no_result ), "Test non-intersecting segments" );

    // Two times the same segment
    vec2 p0_seg0_same( 0., 0. );
    vec2 p1_seg0_same( 1.5, 1.5 );
    vec2 result_same;
    verdict(
        !Intersection::segment_segment( p0_seg0_same, p1_seg0_same, p0_seg0_same,
            p1_seg0_same, result_same ), "Test same segment" );

    // Two intersecting segments
    vec2 p0_seg0_inter( 0., 0. );
    vec2 p1_seg0_inter( 1.5, 1.5 );
    vec2 p0_seg1_inter( 2., 0. );
    vec2 p1_seg1_inter( 0., 2. );
    vec2 result_inter;
    vec2 result_answer( 1., 1. );
    verdict(
        Intersection::segment_segment( p0_seg0_inter, p1_seg0_inter, p0_seg1_inter,
            p1_seg1_inter, result_inter )
            && are_almost_equal( result_inter, result_answer ),
        "Test intersecting segments" );

    // Two intersecting segments from same origin
    vec2 p0_seg0_inter2( 0., 0. );
    vec2 p1_seg0_inter2( 1.5, 1.5 );
    vec2 p1_seg1_inter2( 2., 0. );
    vec2 result_inter2;
    vec2 result_answer2( 0., 0. );
    verdict(
        Intersection::segment_segment( p0_seg0_inter2, p1_seg0_inter2,
            p0_seg0_inter2, p1_seg1_inter2, result_inter2 )
            && are_almost_equal( result_inter2, result_answer2 ),
        "Test intersecting segments from same origin" );

    // Two intersecting segments at extremity
    vec2 p0_seg0_inter3( 0., 0. );
    vec2 p1_seg0_inter3( 1., 1. );
    vec2 p0_seg1_inter3( 2., 0. );
    vec2 p1_seg1_inter3( 0., 2. );
    vec2 result_inter3;
    vec2 result_answer3( 1., 1. );
    verdict(
        Intersection::segment_segment( p0_seg0_inter3, p1_seg0_inter3,
            p0_seg1_inter3, p1_seg1_inter3, result_inter3 )
            && are_almost_equal( result_inter3, result_answer3 ),
        "Test intersecting segments at one extremity" );

    Logger::out( "TEST", " " );
}

void test_segment_line_intersection()
{
    Logger::out( "TEST", "Test Segment-Line intersections" );

    // non-intersecting
    vec2 p0_seg( 0., 0. );
    vec2 p1_seg( 1.5, 1.5 );
    vec2 O_line( 2., 2. );
    vec2 D_line( 0., 2. );
    vec2 no_result;
    verdict(
        !Intersection::segment_line( p0_seg, p1_seg, O_line, D_line, no_result ),
        "Test non-intersecting" );

    // Segment is on the line
    vec2 p0_seg0_same( 0., 0. );
    vec2 p1_seg0_same( 1.5, 1.5 );
    vec2 D_line_same( -1., -1. );
    vec2 result_same;
    verdict(
        !Intersection::segment_line( p0_seg0_same, p1_seg0_same, p0_seg0_same,
            D_line_same, result_same ), "Test segment on line" );

    // intersecting
    vec2 p0_seg0_inter( 0., 0. );
    vec2 p1_seg0_inter( 2., 2. );
    vec2 O_line_inter( 2., 0. );
    vec2 D_line_inter( -0.5, 0.5 );
    vec2 result_inter;
    vec2 result_answer( 1., 1. );
    verdict(
        Intersection::segment_line( p0_seg0_inter, p1_seg0_inter, O_line_inter,
            D_line_inter, result_inter )
            && are_almost_equal( result_inter, result_answer ),
        "Test intersecting" );

    // intersecting from same origin
    vec2 p0_seg0_inter2( 0., 0. );
    vec2 p1_seg0_inter2( 1.5, 1.5 );
    vec2 D_line_inter2( 0., 1. );
    vec2 result_inter2;
    vec2 result_answer2( 0., 0. );
    verdict(
        Intersection::segment_line( p0_seg0_inter2, p1_seg0_inter2, p0_seg0_inter2,
            D_line_inter2, result_inter2 )
            && are_almost_equal( result_inter2, result_answer2 ),
        "Test intersecting from same origin" );

    // intersecting segments at extremity
    vec2 p0_seg0_inter3( 0., 0. );
    vec2 p1_seg0_inter3( 1., 1. );
    vec2 p0_seg1_inter3( 0., 2. );
    vec2 D_line_inter3( -0.5, 0.5 );
    vec2 result_inter3;
    vec2 result_answer3( 1., 1. );
    verdict(
        Intersection::segment_line( p0_seg0_inter3, p1_seg0_inter3, p0_seg1_inter3,
            D_line_inter3, result_inter3 )
            && are_almost_equal( result_inter3, result_answer3 ),
        "Test intersecting segments at one extremity" );

    Logger::out( "TEST", " " );
}

int main()
{
    try {
        default_configure();

        Logger::out( "TEST", "Test intersection algorithms" );

        test_line_plane_intersection();
        test_segment_plane_intersection();
        test_segment_triangle_intersection();
        test_circle_plane_intersection();
        test_disk_segment_intersection();
        test_circle_triangle_intersection();
        test_plane_plane_intersection();
        test_line_line_intersection();
        test_segment_segment_intersection();
        test_segment_line_intersection();

    } catch( const RINGMeshException& e ) {
        Logger::err( e.category(), e.what() );
        return 1;
    } catch( const std::exception& e ) {
        Logger::err( "Exception", e.what() );
        return 1;
    }
    Logger::out( "TEST", "SUCCESS" );
    return 0;
}
