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
    vec3 origin{ 100, 200, -4 };
    vec3 x{ 1, 2, -1 };
    vec3 y{ 2, 1, 4 };
    vec3 z{ 3, -2, -1 };

    Frame3D frame{ x, y, z };
    ReferenceFrame3D reference_frame{ origin, frame };

    vec3 test_coordinates{ 20, 10, 30 };
    vec3 frame_test_coordinates =
        ReferenceFrameManipulator3D::coords_from_global_to_frame(
            reference_frame, test_coordinates );
    vec3 operation_test_coordinates =
        ReferenceFrameManipulator3D::coords_from_frame_to_global(
            reference_frame, frame_test_coordinates );
    if( !inexact_equal( test_coordinates, operation_test_coordinates, 1e-7 ) )
    {
        throw RINGMeshException(
            "TEST", "Error in coordinate reference change" );
    }
    if( !ReferenceFrameManipulator3D::is_frame_orthogonal( reference_frame ) )
    {
        throw RINGMeshException(
            "TEST", "Error in checking the orthogonality of a frame" );
    }

    ReferenceFrame3D inverse_frame =
        ReferenceFrameManipulator3D::reference_frame_from_global_to_local(
            reference_frame );
    if( !inexact_equal( inverse_frame,
            ReferenceFrameManipulator3D::
                orthogonal_reference_frame_from_global_to_local(
                    reference_frame ),
            1e-7 ) )
    {
        throw RINGMeshException(
            "TEST", "Error in orthogonal reference frame change" );
    }

    ivec3 grid_dimensions{ 10, 8, 9 };
    CartesianGrid3D cartesiangrid{ grid_dimensions, reference_frame };
    if( cartesiangrid.cell_volume() != 42. )
    {
        throw RINGMeshException( "TEST", "Error in cell volume" );
    }
    if( cartesiangrid.cell_offset_from_global_point( vec3{ 103, 200.5, -3 } )
        != 0 )
    {
        throw RINGMeshException( "TEST", "Error in calculating the offset 1" );
    }
    if( cartesiangrid.cell_offset_from_global_point( vec3{ 107, 201.5, -1 } )
        != 91 )
    {
        throw RINGMeshException( "TEST", "Error in calculating the offset 2" );
    }
    if( cartesiangrid.cell_offset_from_global_point( vec3{ 1, 1, 1 } ) != -1 )
    {
        throw RINGMeshException( "TEST", "Error in calculating the offset 3" );
    }
}

int main()
{
    using namespace RINGMesh;

    try
    {
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
