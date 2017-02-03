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

#include <ringmesh/basic/geometry.h>

/*!
 * @author Pierre Anquez
 */

using namespace RINGMesh ;

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

}

void test_circle_plane_intersection()
{
    Logger::out( "TEST" ) << "Test Circle-Plane intersections" << std::endl ;

}

void test_disk_segment_intersection()
{
    Logger::out( "TEST" ) << "Test Disk-Segment intersections" << std::endl ;

}

void test_circle_triangle_intersection()
{
    Logger::out( "TEST" ) << "Test Circle-Triangle intersections" << std::endl ;

}

void test_plane_plane_intersection()
{
    Logger::out( "TEST" ) << "Test Plane-Plane intersections" << std::endl ;

}

int main()
{
    try {
        default_configure() ;

        Logger::out( "TEST" ) << "Test intersection algorithms" << std::endl ;

//        test_nn_search_ringmesh() ;
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
