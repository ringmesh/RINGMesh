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

#include <iostream>
#include <ringmesh/basic/box.h>
#include <ringmesh/basic/logger.h>
#include <ringmesh/ringmesh_tests_config.h>
#include <tuple>
/*!
 * @author Marine Hubert
 */

using namespace RINGMesh;

namespace
{
    template < typename T >
    inline T sqr( T x )
    {
        return x * x;
    }
}

void test_boxes()
{
    Logger::out( "TEST", "Test box functions" );
    Box3D B;
    vec3 pt0{ 0, 0, 0 };
    vec3 pt1{ 1, 1, 1 };
    vec3 pt2{ 2, 2, 2 };
    vec3 pt3{ 3, 3, 3 };
    vec3 pt4{ 4, 4, 4 };

    // At initialization, initialized_ = false;
    if( B.initialized() )
    {
        throw RINGMeshException( "TEST", "Error in initialized" );
    }

    // Test addition of a point
    B.add_point( pt1 );
    if( B.min() != pt1 )
    {
        throw RINGMeshException( "TEST", "Error in addition" );
    }

    if( B.max() != pt1 )
    {
        throw RINGMeshException( "TEST", "Error in addition" );
    }

    // Test initialized_ = true;
    if( !B.initialized() )
    {
        throw RINGMeshException( "TEST", "Error in initialized" );
    }

    // Check min
    B.add_point( pt3 );

    if( B.min() != pt1 )
    {
        throw RINGMeshException( "TEST", "Error in minimum" );
    }
    // Check max
    if( B.max() != pt3 )
    {
        throw RINGMeshException( "TEST", "Error in maximum" );
    }

    // Check center
    vec3 center_test = pt2;
    if( B.center() != center_test )
    {
        throw RINGMeshException( "TEST", "Error in center" );
    }

    // Check diagonal
    vec3 diag = pt2;

    if( B.diagonal() != diag )
    {
        throw RINGMeshException( "TEST", "Error in diagonal" );
    }

    // Check add_box
    Box3D A;
    A.add_box( B );
    if( ( B.min() != A.min() ) || ( A.max() != B.max() ) )
    {
        throw RINGMeshException( "TEST", "Error in add_box" );
    }

    // Check clear
    A.clear();
    if( A.initialized() )
    {
        throw RINGMeshException( "TEST", "Error in clear" );
    }

    // Check bboxes_overlap
    Box3D C;

    C.add_point( pt2 );
    C.add_point( pt0 );
    A.add_point( pt3 );
    A.add_point( pt1 );

    if( !C.bboxes_overlap( A ) )
    {
        throw RINGMeshException( "TEST", "Error in overlap maximum" );
    }
    if( !A.bboxes_overlap( C ) )
    {
        throw RINGMeshException( "TEST", "Error in overlap minimum" );
    }

    // Check bbox_union
    Box3D E = A.bbox_union( C );

    if( ( E.min() != C.min() ) || ( E.max() != A.max() ) )
    {
        throw RINGMeshException( "TEST", "Error in box union" );
    }

    // Check bbox intersection
    bool do_boxes_intersect;
    Box3D F;
    std::tie( do_boxes_intersect, F ) = ( C.bbox_intersection( A ) );
    if( !do_boxes_intersect )
    {
        throw RINGMeshException( "TEST", "Error in box intersection " );
    }
    if( ( F.min() != pt1 ) || ( F.max() != pt2 ) )
    {
        throw RINGMeshException( "TEST", "Error in box intersection" );
    }

    // Check contains
    if( !E.contains( pt0 ) )
    {
        throw RINGMeshException( "TEST", "Error in box contains" );
    }
    if( !E.contains( pt2 ) )
    {
        throw RINGMeshException( "TEST", "Error in box contains" );
    }
    if( E.contains( pt4 ) )
    {
        throw RINGMeshException( "TEST", "Error in box contains" );
    }

    // Check signed distance
    if( A.signed_distance( pt2 ) != -1 )
    {
        throw RINGMeshException(
            "TEST", "Error in box signed distance (inside)" );
    }
    if( A.signed_distance( pt1 ) != 0 )
    {
        throw RINGMeshException( "TEST", "Error in box signed distance (min)" );
    }
    if( A.signed_distance( pt4 ) != 3 )
    {
        throw RINGMeshException(
            "TEST", "Error in box signed distance (outside)" );
    }
}

int main()
{
    try
    {
        Logger::out( "TEST", "Test boxes :" );
        test_boxes();
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
