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

#include <ringmesh/ringmesh_tests_config.h>

#include <vector>

#include <geogram/basic/vecg.h>

#include <ringmesh/basic/geometry.h>
#include <ringmesh/basic/logger.h>

/*!
 * @author Pierre Anquez
 */

using namespace RINGMesh;

void verdict( bool condition, std::string test_name )
{
    if( !condition )
    {
        throw RINGMeshException( "TEST", test_name, ": KO" );
    }
    else
    {
        Logger::out( "TEST", test_name, ": OK" );
    }
}

template < index_t DIMENSION >
bool are_almost_equal(
    const vecn< DIMENSION >& vec0, const vecn< DIMENSION >& vec1 )
{
    return ( vec0 - vec1 ).length2() < global_epsilon_sq;
}

void test_line_plane_intersection()
{
    Logger::out( "TEST", "Test Line-Plane intersections" );
    Geometry::Plane plane{ { 1., 0., 0. }, { 1., 4., -2 } };

    // Intersection is the origin point of the line
    Geometry::Line3D line1{ { -2., 4., 1. }, { 1., 5., 0. } };
    bool does_line1_intersect_plane;
    vec3 result1;
    std::tie( does_line1_intersect_plane, result1 ) =
        Intersection::line_plane( line1, plane );
    vec3 answer1{ line1.origin };
    verdict( does_line1_intersect_plane && result1 == answer1,
        "True intersection1" );

    // Intersection is a point
    Geometry::Line3D line2{ { -2., 4., 1. }, { -41., 7., -28. } };
    bool does_line2_intersect_plane;
    vec3 result2;
    std::tie( does_line2_intersect_plane, result2 ) =
        Intersection::line_plane( line2, plane );
    vec3 answer2{ 1., -77., -49. };
    verdict( does_line2_intersect_plane && result2 == answer2,
        "True intersection2" );

    // The line is parallel to the plane
    Geometry::Line3D line3{ { 0., 2., 1. }, { 0., 1., 8. } };
    bool does_line3_intersect_plane;
    std::tie( does_line3_intersect_plane, std::ignore ) =
        Intersection::line_plane( line3, plane );
    verdict( !does_line3_intersect_plane, "Line parallel to the plane" );

    // The line is included into the plane
    Geometry::Line3D line4{ { 0., 2., 1. }, { 1., 1., 8. } };
    bool does_line4_intersect_plane;
    std::tie( does_line4_intersect_plane, std::ignore ) =
        Intersection::line_plane( line4, plane );
    verdict( !does_line4_intersect_plane, "Line included into the plane" );

    Logger::out( "TEST", " " );
}

void test_segment_plane_intersection()
{
    Logger::out( "TEST", "Test Segment-Plane intersections" );
    Geometry::Plane plane{ { 1., 1.25, -2 }, { 1., -2., 3. } };

    // Intersection is a point
    Geometry::Segment3D seg1{ { 3., -1., 6.75 }, { -3., -10., -8.25 } };
    bool does_seg1_intersect_plane;
    vec3 result1;
    std::tie( does_seg1_intersect_plane, result1 ) =
        Intersection::segment_plane( seg1, plane );
    vec3 answer1{ 1., -4., 1.75 };
    verdict( does_seg1_intersect_plane && are_almost_equal( result1, answer1 ),
        "True intersection" );

    // Intersection is an extremal point of the segment
    Geometry::Segment3D seg2{ { 3., -1., 6.75 }, { 1., -4., 1.75 } };
    bool does_seg2_intersect_plane;
    vec3 result2;
    std::tie( does_seg2_intersect_plane, result2 ) =
        Intersection::segment_plane( seg2, plane );
    vec3 answer2{ seg2.p1 };
    verdict( does_seg2_intersect_plane && are_almost_equal( result2, answer2 ),
        "Intersection at segment extremity" );

    // Line intersects but not the segment
    Geometry::Segment3D seg3{ { 3., -1., 6.75 }, { 2., -2.5, 4.25 } };
    bool does_seg3_intersect_plane;
    std::tie( does_seg3_intersect_plane, std::ignore ) =
        Intersection::segment_plane( seg3, plane );
    verdict( !does_seg3_intersect_plane,
        "Intersection on line but but into segment" );

    // The segment is parallel to the plane
    Geometry::Segment3D seg4{ { 0., 0., 6.75 }, { 1., 0., 7.25 } };
    bool does_seg4_intersect_plane;
    std::tie( does_seg4_intersect_plane, std::ignore ) =
        Intersection::segment_plane( seg4, plane );
    verdict( !does_seg4_intersect_plane, "Segment parallel to plane" );

    // The segment is included into the plane
    Geometry::Segment3D seg5{ { 5., -2., 5. }, { 1., -2., 3. } };
    bool does_seg5_intersect_plane;
    std::tie( does_seg5_intersect_plane, std::ignore ) =
        Intersection::segment_plane( seg5, plane );
    verdict( !does_seg5_intersect_plane, "Segment included into the plane" );

    Logger::out( "TEST", " " );
}

void test_segment_triangle_intersection()
{
    Logger::out( "TEST", "Test Segment-Triangle intersections" );
    Geometry::Triangle3D triangle{ { 1., 1., 0. }, { 2., 3., 0. },
        { 4., -1., 0. } };

    // True intersection inside the triangle
    Geometry::Segment3D seg1{ { 2., 2., 3. }, { 2., 2., -3. } };
    bool does_seg1_intersect_triangle;
    vec3 result1;
    std::tie( does_seg1_intersect_triangle, result1 ) =
        Intersection::segment_triangle( seg1, triangle );
    vec3 answer1{ 2., 2., 0. };
    verdict(
        does_seg1_intersect_triangle && are_almost_equal( result1, answer1 ),
        "Test1" );

    // No intersection
    Geometry::Segment3D seg2{ { 2., 2., 3. }, { 20., 2., -1. } };
    bool does_seg2_intersect_triangle;
    std::tie( does_seg2_intersect_triangle, std::ignore ) =
        Intersection::segment_triangle( seg2, triangle );
    verdict( !does_seg2_intersect_triangle, "Test2" );

    // Intersection at a triangle vertex
    Geometry::Segment3D seg3{ { 1., 4., 3. }, { 3., 2., -3. } };
    bool does_seg3_intersect_triangle;
    vec3 result3;
    std::tie( does_seg3_intersect_triangle, result3 ) =
        Intersection::segment_triangle( seg3, triangle );
    vec3 answer3{ triangle.p1 };
    verdict(
        does_seg3_intersect_triangle && are_almost_equal( result3, answer3 ),
        "Test3" );

    // Intersection on a triangle edge
    Geometry::Segment3D seg4{ { 1., -4., 3. }, { 4., 4., -3. } };
    bool does_seg4_intersect_triangle;
    vec3 result4;
    std::tie( does_seg4_intersect_triangle, result4 ) =
        Intersection::segment_triangle( seg4, triangle );
    vec3 answer4{ 2.5, 0., 0. };
    verdict(
        does_seg4_intersect_triangle && are_almost_equal( result4, answer4 ),
        "Test4" );

    // Segment included inside the triangle
    Geometry::Segment3D seg5{ { 2., 2., 0. }, { 3., 1., 0. } };
    bool does_seg5_intersect_triangle;
    std::tie( does_seg5_intersect_triangle, std::ignore ) =
        Intersection::segment_triangle( seg5, triangle );
    verdict( !does_seg5_intersect_triangle, "Test5" );

    // Segment is a triangle edge
    Geometry::Segment3D seg6{ triangle.p0, triangle.p1 };
    bool does_seg6_intersect_triangle;
    std::tie( does_seg6_intersect_triangle, std::ignore ) =
        Intersection::segment_triangle( seg6, triangle );
    verdict( !does_seg6_intersect_triangle, "Test6" );

    // Segment in the same plane than triangle, one point inside the other
    // outside
    Geometry::Segment3D seg7{ { 2., 2., 0. }, { 4., 1., -0. } };
    bool does_seg7_intersect_triangle = true;
    std::tie( does_seg7_intersect_triangle, std::ignore ) =
        Intersection::segment_triangle( seg7, triangle );
    verdict( !does_seg7_intersect_triangle, "Test7" );

    Logger::out( "TEST", " " );
}

void test_circle_plane_intersection()
{
    Logger::out( "TEST", "Test Circle-Plane intersections" );
    Geometry::Plane plane{ { 0., 0., -2. }, { 3., 1., -1. } };

    // Circle parallel to the plane
    Geometry::Circle circle1{ { { 0., 0., 1. }, { 2., 3., 4. } }, 4. };
    bool does_circle1_intersect_plane;
    std::tie( does_circle1_intersect_plane, std::ignore ) =
        Intersection::circle_plane( circle1, plane );
    verdict( !does_circle1_intersect_plane, "Test circle parallel to plane" );


    // Circle adjacent to the plane
    Geometry::Circle circle2{ { { -1., 2., 0. }, { 2., 3., 4. } }, 5. };
    bool does_circle2_intersect_plane;
    std::vector< vec3 > results2;
    std::tie( does_circle2_intersect_plane, results2 ) =
        Intersection::circle_plane( circle2, plane );
    vec3 answer2{ 2., 3., -1. };
    verdict( does_circle2_intersect_plane && results2.size() == 1,
        "Test circle adjacent to the plane (number of points)" );
    verdict( are_almost_equal( results2[0], answer2 ),
        "Test circle adjacent to the plane (intersection coordinates)" );

    // Circle crossing the plane
    Geometry::Circle circle3{ { { 1., -1., 0. }, { 2., 3., 0. } }, 2. };
    bool does_circle3_intersect_plane;
    std::vector< vec3 > results3;
    std::tie( does_circle3_intersect_plane, results3 ) =
        Intersection::circle_plane( circle3, plane );
    vec3 answer31{ circle3.plane.origin
                   + vec3{ std::sqrt( 2 ) * std::cos( M_PI / 6 ),
                         std::sqrt( 2 ) * std::cos( M_PI / 6 ), -1. } };
    vec3 answer32{ circle3.plane.origin
                   + vec3{ -std::sqrt( 2 ) * std::cos( M_PI / 6 ),
                         -std::sqrt( 2 ) * std::cos( M_PI / 6 ), -1. } };
    verdict( does_circle3_intersect_plane && results3.size() == 2,
        "Test circle crossing the plane (number of points)" );
    verdict( ( are_almost_equal( results3[0], answer31 )
                 || ( are_almost_equal( results3[0], answer32 ) ) )
                 && ( are_almost_equal( results3[1], answer31 )
                        || ( are_almost_equal( results3[1], answer32 ) ) ),
        "Test circle crossing the plane (intersection coordinates)" );

    // Circle distance from the plane

    Geometry::Circle circle4{ { { -2., -1., 1. }, { 3., 0., 3. } }, 2. };
    bool does_circle4_intersect_plane;
    std::tie( does_circle4_intersect_plane, std::ignore ) =
    	Intersection::circle_plane( circle4, plane );
    verdict( !does_circle4_intersect_plane, "Test circle distance from plane" );

           Logger::out( "TEST", " " );
}

void test_disk_segment_intersection()
{
    Logger::out( "TEST", "Test Disk-Segment intersections" );
    Geometry::Disk disk{ { { 0., 4., 0. }, { 2., 2., 2. } }, 4. };

    // Segment in the disk plane
    Geometry::Segment3D seg1{ { 1., 2., 3. }, { 3., 2., 1. } };
    bool does_seg1_intersect_disk;
    std::tie( does_seg1_intersect_disk, std::ignore ) =
        Intersection::segment_disk( seg1, disk );
    verdict( !does_seg1_intersect_disk, "Test segment inside disk" );

    // Segment adjacent to the disk
    Geometry::Segment3D seg2{ { -2., 0., -2. }, { -2., 4., 4. } };
    bool does_seg2_intersect_disk;
    std::tie( does_seg2_intersect_disk, std::ignore ) =
        Intersection::segment_disk( seg2, disk );
    verdict( !does_seg2_intersect_disk, "Test segment tangent to the disk" );

    // Circle crossing the disk
    Geometry::Segment3D seg3{ { 1., 1., 3. }, { 3., 3., 1. } };
    bool does_seg3_intersect_disk;
    vec3 result3;
    std::tie( does_seg3_intersect_disk, result3 ) =
        Intersection::segment_disk( seg3, disk );
    verdict( does_seg3_intersect_disk,
        "Test circle adjacent to the plane (intersection exists)" );
    vec3 answer3{ 2., 2., 2. };
    verdict( are_almost_equal( result3, answer3 ),
        "Test circle adjacent to the plane (intersection coordinates)" );

    Logger::out( "TEST", " " );
}

void test_circle_triangle_intersection()
{
    Logger::out( "TEST", "Test Circle-Triangle intersections" );
    Geometry::Triangle3D triangle{ { 1., 1., 0. }, { 2., 3., 0. },
        { 4., -1., 0. } };

    // The circle is (almost) adjacent to the triangle plane on a triangle
    // border
    Geometry::Circle circle1{ { { 1., 1., 0. }, { 1.5, 2., 4. } }, 4. };
    bool does_circle1_intersect_triangle;
    std::vector< vec3 > results1;
    std::tie( does_circle1_intersect_triangle, results1 ) =
        Intersection::triangle_circle( triangle, circle1 );
    verdict( !does_circle1_intersect_triangle, "Test1 (number of points)" );

    // One point inside triangle the other one outside
    Geometry::Circle circle2{ { { 1., 1., 0. }, { 2., 2., 0. } }, 1. };
    bool does_circle2_intersect_triangle;
    std::vector< vec3 > results2;
    std::tie( does_circle2_intersect_triangle, results2 ) =
        Intersection::triangle_circle( triangle, circle2 );
    vec3 answer2{ 2. + 0.5 * std::sqrt( 2. ), 2. - 0.5 * std::sqrt( 2. ), 0. };
    verdict( does_circle2_intersect_triangle && results2.size() == 1,
        "Test2 (number of points)" );
    verdict( are_almost_equal( results2[0], answer2 ),
        "Test2 (intersection coordinates)" );

    // The circle is (almost) adjacent to the triangle plane on a triangle
    // vertex
    Geometry::Circle circle3{ { { 1., 1., 0. }, { 1., 1., 4. } }, 4. };
    bool does_circle3_intersect_triangle;
    std::vector< vec3 > results3;
    std::tie( does_circle3_intersect_triangle, results3 ) =
        Intersection::triangle_circle( triangle, circle3 );
    verdict( !does_circle3_intersect_triangle, "Test3 (number of points)" );

    Logger::out( "TEST", " " );
}

void test_plane_plane_intersection()
{
    Logger::out( "TEST", "Test Plane-Plane intersections" );
    Geometry::Plane plane0{ { 1., -2., 4. }, { 4., -2., 0. } };

    // Two parallel planes
    Geometry::Plane plane1{ { -2., 4., -8. }, { 6., 0., 1.52 } };
    bool does_P1_intersect_plane;
    std::tie( does_P1_intersect_plane, std::ignore ) =
        Intersection::plane_plane( plane0, plane1 );
    verdict( !does_P1_intersect_plane, "Test parallel planes" );

    // Two times the same plane
    Geometry::Plane plane2{ { -1., 2., -4. }, { 4., -2., 0. } };
    bool does_P2_intersect_plane;
    std::tie( does_P2_intersect_plane, std::ignore ) =
        Intersection::plane_plane( plane0, plane2 );
    verdict( !does_P2_intersect_plane, "Test same plane" );

    // Two intersecting plane
    Geometry::Plane plane3{ { 0., 2., -4. }, { -3., -1., 1. } };
    bool does_P3_intersect_plane;
    Geometry::Line3D result3;
    std::tie( does_P3_intersect_plane, result3 ) =
        Intersection::plane_plane( plane0, plane3 );
    vec3 O_inter_answer3{ 2., -3., 0. };
    vec3 D_inter_answer3{ 0., 2., 1. };
    verdict( does_P3_intersect_plane
                 && are_almost_equal(
                        normalize( D_inter_answer3 ), result3.direction )
                 && are_almost_equal( normalize( D_inter_answer3 ),
                        normalize( result3.origin - O_inter_answer3 ) ),
        "Test intersecting planes" );
    Logger::out( "TEST", " " );
}

void test_line_line_intersection()
{
    Logger::out( "TEST", "Test Line-Line intersections" );
    Geometry::Line2D line{ { 1.5, 1.5 }, { 0., 0. } };

    // Two parallel lines
    Geometry::Line2D line_parallel{ { 1.5, 1.5 }, { 1., 1. } };
    bool does_parallel_lines_intersect;
    std::tie( does_parallel_lines_intersect, std::ignore ) =
        Intersection::line_line( line, line_parallel );
    verdict( !does_parallel_lines_intersect, "Test parallel lines" );

    // Two times the same line
    Geometry::Line2D Line_same{ { -2.5, -2.5 }, { 0., 0. } };
    bool does_same_lines_intersect;
    std::tie( does_same_lines_intersect, std::ignore ) =
        Intersection::line_line( line, Line_same );
    verdict( !does_same_lines_intersect, "Test same line" );

    // Two intersecting lines
    Geometry::Line2D Line_inter{ { 2.5, -2.5 }, { 2., 0. } };
    bool does_line_intersect_line;
    vec2 result_inter;
    vec2 result_answer{ 1., 1. };
    std::tie( does_line_intersect_line, result_inter ) =
        Intersection::line_line( line, Line_inter );
    verdict( does_line_intersect_line
                 && are_almost_equal( result_inter, result_answer ),
        "Test intersecting lines" );

    Logger::out( "TEST", " " );
}

void test_segment_segment_intersection()
{
    Logger::out( "TEST", "Test Segment-Segment intersections" );
    Geometry::Segment2D segment{ { 0., 0. }, { 1.5, 1.5 } };

    // Two non-intersecting segments
    Geometry::Segment2D segment1{ { 2., 2. }, { 3., 2. } };
    bool do_segments_intersect;
    std::tie( do_segments_intersect, std::ignore ) =
        Intersection::segment_segment( segment, segment1 );
    verdict( !do_segments_intersect, "Test non-intersecting segments" );

    // Two times the same segment
    std::tie( do_segments_intersect, std::ignore ) =
        Intersection::segment_segment( segment, segment );
    verdict( !do_segments_intersect, "Test same segment" );

    // Two intersecting segments
    Geometry::Segment2D segment_inter{ { 2., 0. }, { 0., 2. } };
    vec2 result_inter;
    std::tie( do_segments_intersect, result_inter ) =
        Intersection::segment_segment( segment, segment_inter );
    vec2 result_answer{ 1., 1. };
    verdict( do_segments_intersect
                 && are_almost_equal( result_inter, result_answer ),
        "Test intersecting segments" );

    // Two intersecting segments from same origin
    Geometry::Segment2D segment_inter2{ { 0., 0. }, { 2., 0. } };
    vec2 result_inter2;
    std::tie( do_segments_intersect, result_inter2 ) =
        Intersection::segment_segment( segment, segment_inter2 );
    vec2 result_answer2{ 0., 0. };
    verdict( do_segments_intersect
                 && are_almost_equal( result_inter2, result_answer2 ),
        "Test intersecting segments from same origin" );

    // Two intersecting segments at extremity
    Geometry::Segment2D segment_inter3{ { 0., 0. }, { 1., 1. } };
    vec2 result_inter3;
    std::tie( do_segments_intersect, result_inter3 ) =
        Intersection::segment_segment( segment_inter, segment_inter3 );
    vec2 result_answer3{ 1., 1. };
    verdict( do_segments_intersect
                 && are_almost_equal( result_inter3, result_answer3 ),
        "Test intersecting segments at one extremity" );

    Logger::out( "TEST", " " );
}

void test_segment_line_intersection()
{
    Logger::out( "TEST", "Test Segment-Line intersections" );
    Geometry::Segment2D segment{ { 0., 0. }, { 1.5, 1.5 } };

    // non-intersecting
    Geometry::Line2D line{ { 0., 2. }, { 2., 2. } };
    bool does_segment_intersect_line;
    std::tie( does_segment_intersect_line, std::ignore ) =
        Intersection::segment_line( segment, line );
    verdict( !does_segment_intersect_line, "Test non-intersecting" );

    // Segment is on the line
    Geometry::Line2D line_same{ { -1., -1. }, { 1.5, 1.5 } };
    std::tie( does_segment_intersect_line, std::ignore ) =
        Intersection::segment_line( segment, line_same );
    verdict( !does_segment_intersect_line, "Test segment on line" );

    // intersecting
    Geometry::Segment2D segment_inter{ { 0., 0. }, { 2., 2. } };
    Geometry::Line2D line_inter{ { -0.5, 0.5 }, { 2., 0. } };
    vec2 result_inter;
    std::tie( does_segment_intersect_line, result_inter ) =
        Intersection::segment_line( segment_inter, line_inter );
    vec2 result_answer{ 1., 1. };
    verdict( does_segment_intersect_line
                 && are_almost_equal( result_inter, result_answer ),
        "Test intersecting" );

    // intersecting from same origin
    Geometry::Line2D line_inter2{ { 0., 1. }, { 0., 0. } };
    vec2 result_inter2;
    std::tie( does_segment_intersect_line, result_inter2 ) =
        Intersection::segment_line( segment, line_inter2 );
    vec2 result_answer2{ 0., 0. };
    verdict( does_segment_intersect_line
                 && are_almost_equal( result_inter2, result_answer2 ),
        "Test intersecting from same origin" );

    // intersecting segments at extremity
    Geometry::Segment2D segment_inter3{ { 0., 0. }, { 1., 1. } };
    Geometry::Line2D line_inter3{ { -0.5, 0.5 }, { 0., 2. } };
    vec2 result_inter3;
    std::tie( does_segment_intersect_line, result_inter3 ) =
        Intersection::segment_line( segment_inter3, line_inter3 );
    vec2 result_answer3{ 1., 1. };
    verdict( does_segment_intersect_line
                 && are_almost_equal( result_inter3, result_answer3 ),
        "Test intersecting segments at one extremity" );

    Logger::out( "TEST", " " );
}

void test_line_sphere_intersection()
{
    Logger::out( "TEST", "Test Line-Sphere intersections" );
    Geometry::Sphere3D sphere{ { 2., 2., 2. }, 4. };

    // Line outside the sphere
    Geometry::Line3D line_outside{ { -3., 2., 1. }, { 10., 10., 10. } };
    bool line_outside_intersect_sphere;
    std::tie( line_outside_intersect_sphere, std::ignore ) =
        Intersection::line_sphere( line_outside, sphere );
    verdict( !line_outside_intersect_sphere, "Test line outside sphere" );

    // Line tangent to the sphere
    Geometry::Line3D line_tangent{ { 0., 1., 1. }, { -2., 5., 5. } };
    bool line_tangent_intersect_sphere;
    std::vector< vec3 > tangent_result;
    std::tie( line_tangent_intersect_sphere, tangent_result ) =
        Intersection::line_sphere( line_tangent, sphere );
    verdict( tangent_result.size() == 1, "Test line tangent to the sphere" );
    vec3 tangent_answer{ -2., 2., 2. };
    verdict( are_almost_equal( tangent_result.front(), tangent_answer ),
        "Test line tangent to the sphere (intersection coordinates)" );

    // Line crossing the sphere
    Geometry::Line3D line_cross{ { 3., 0., 0. }, { 0., 2., 2. } };
    bool line_cross_intersect_sphere;
    std::vector< vec3 > cross_result;
    std::tie( line_cross_intersect_sphere, cross_result ) =
        Intersection::line_sphere( line_cross, sphere );
    verdict( cross_result.size() == 2,
        "Test line crossing the sphere (intersection exists)" );
    vec3 answer_cross0{ -2., 2., 2. };
    vec3 answer_cross1{ 6., 2., 2. };
    verdict( are_almost_equal( cross_result.front(), answer_cross0 ),
        "Test line crossing the sphere (intersection coordinates)" );
    verdict( are_almost_equal( cross_result.back(), answer_cross1 ),
        "Test line crossing the sphere (intersection coordinates)" );

    Logger::out( "TEST", " " );
}

void test_segment_sphere_intersection()
{
    Logger::out( "TEST", "Test Segment-Sphere intersections" );
    Geometry::Sphere3D sphere{ { 2., 2., 2. }, 4. };

    // Segment outside the sphere
    Geometry::Segment3D seg_outside{ { 10., 10., 10. }, { 15., 20., 10. } };
    bool segment_outside_intersect_sphere;
    std::tie( segment_outside_intersect_sphere, std::ignore ) =
        Intersection::segment_sphere( seg_outside, sphere );
    verdict( !segment_outside_intersect_sphere, "Test segment outside sphere" );

    // Segment tangent to the sphere
    Geometry::Segment3D seg_tangent{ { -2., 5., 5. }, { -2., -5., -5. } };
    bool segment_tangent_intersect_sphere;
    std::vector< vec3 > tangent_result;
    std::tie( segment_tangent_intersect_sphere, tangent_result ) =
        Intersection::segment_sphere( seg_tangent, sphere );
    verdict( tangent_result.size() == 1, "Test segment tangent to the sphere" );
    vec3 tangent_answer{ -2., 2., 2. };
    verdict( are_almost_equal( tangent_result.front(), tangent_answer ),
        "Test segment tangent to the sphere (intersection coordinates)" );

    // Segment crossing the sphere
    Geometry::Segment3D seg_cross{ { -10., 2., 2. }, { 10., 2., 2. } };
    bool segment_cross_intersect_sphere;
    std::vector< vec3 > cross_result;
    std::tie( segment_cross_intersect_sphere, cross_result ) =
        Intersection::segment_sphere( seg_cross, sphere );
    verdict( cross_result.size() == 2,
        "Test segment crossing the sphere (2 intersections)" );
    vec3 answer_cross0{ -2., 2., 2. };
    vec3 answer_cross1{ 6., 2., 2. };
    verdict( are_almost_equal( cross_result.front(), answer_cross0 ),
        "Test segment crossing the sphere (intersection coordinates)" );
    verdict( are_almost_equal( cross_result.back(), answer_cross1 ),
        "Test segment crossing the sphere (intersection coordinates)" );

    // Segment crossing the sphere
    Geometry::Segment3D seg_cross2{ { 2., 3., 2. }, { 2., 8., 2. } };
    bool segment_cross_intersect_sphere2;
    std::vector< vec3 > cross_result2;
    std::tie( segment_cross_intersect_sphere2, cross_result2 ) =
        Intersection::segment_sphere( seg_cross2, sphere );
    verdict( cross_result2.size() == 1,
        "Test segment crossing the sphere (1 intersection)" );
    vec3 answer_cross3{ 2., 6., 2. };
    verdict( are_almost_equal( cross_result2.front(), answer_cross3 ),
        "Test segment crossing the sphere (intersection coordinates)" );

    Logger::out( "TEST", " " );
}

int main()
{
    try
    {
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
        test_line_sphere_intersection();
        test_segment_sphere_intersection();
    }
    catch( const RINGMeshException& e )
    {
        Logger::err( e.category(), e.what() );
        return 1;
    }
    catch( const std::exception& e )
    {
        Logger::err( "Exception", e.what() );
        return 1;
    }
    Logger::out( "TEST", "SUCCESS" );
    return 0;
}
