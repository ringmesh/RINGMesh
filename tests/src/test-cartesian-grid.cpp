/*
 * Copyright (c) 2012-2017, Association Scientifique pour la Geologie et ses
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

#include <ringmesh/basic/geometry.h>
#include <ringmesh/basic/logger.h>
#include <ringmesh/mesh/cartesian_grid.h>

/*!
 * @file Tests for the frame structures and the cartesian grids
 * @author Melchior Schuh-Senlis
 */

using namespace RINGMesh;

void test_frames_and_cartesian_grid()
{
    vec3 origin{ 1, 1, 1 };
    vec3 x{ 0.5, 0, 0 };
    vec3 y{ 0, 2, 0 };
    vec3 z{ 0, 0, 4 };

    Frame3D frame{ x, y, z };
    ReferenceFrame3D reference_frame{ origin, frame };

    vec3 test_coordinates{ 5, 2, 1 };
    vec3 frame_test_coordinates =
        ReferenceFrameManipulator3D::coords_from_global_to_frame(
            reference_frame, test_coordinates );
    vec3 test_coordinates2 =
        ReferenceFrameManipulator3D::coords_from_frame_to_global(
            reference_frame, frame_test_coordinates );
    if( test_coordinates != test_coordinates2 )
    {
        throw RINGMeshException(
            "TEST", "Error in coordinate reference change" );
    }
    if( !ReferenceFrameManipulator3D::frame_is_orthogonal( reference_frame ) )
    {
        throw RINGMeshException(
            "TEST", "Error in checking the orthogonality of a frame" );
    }
    vec3 z2{ 1, 1, 2 };
    Frame3D frame2{ x, y, z2 };
    ReferenceFrame3D reference_frame2{ origin, frame2 };
    if( ReferenceFrameManipulator3D::frame_is_orthogonal( reference_frame2 ) )
    {
        throw RINGMeshException(
            "TEST", "Error in checking the orthogonality of a frame" );
    }

    ReferenceFrame3D inverse_frame =
        ReferenceFrameManipulator3D::reference_frame_from_global_to_local(
            reference_frame );
    if( inverse_frame.origin() != vec3{ -2, -0.5, -0.25 }
        || inverse_frame[0] != vec3{ 2, 0, 0 }
        || inverse_frame[1] != vec3{ 0, 0.5, 0 }
        || inverse_frame[2] != vec3{ 0, 0, 0.25 } )
    {
        throw RINGMeshException( "TEST", "Error in reference frame change" );
    }

    ivec3 dimensions_grille{ 10, 8, 9 };
    CartesianGrid< 3 > cartesiangrid{ dimensions_grille, reference_frame };
}

int main()
{
    using namespace RINGMesh;

    try
    {
        default_configure();

        Logger::out( "TEST", "Frames and Cartesian Grid" );

        test_frames_and_cartesian_grid();
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
