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

#include <ringmesh/basic/geometry.h>
#include <ringmesh/basic/logger.h>
#include <ringmesh/mesh/cartesian_grid.h>

/*!
 * @file Tests for the frame structures and the cartesian grids
 * @author Melchior Schuh-Senlis
 */

using namespace RINGMesh;

void test_frames()
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
}

void test_cartesian_grids()
{
    vec3 origin{ 100, 200, -4 };
    vec3 x{ 1, 2, -1 };
    vec3 y{ 2, 1, 4 };
    vec3 z{ 3, -2, -1 };

    Frame3D frame{ x, y, z };
    ReferenceFrame3D reference_frame{ origin, frame };
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
    if( cartesiangrid.cell_offset_from_global_point( vec3{ 1, 1, 1 } )
        != NO_ID )
    {
        throw RINGMeshException( "TEST", "Error in calculating the offset 3" );
    }
    if( cartesiangrid.cell_offset( cartesiangrid.local_from_offset( 51 ) )
        != 51 )
    {
        throw RINGMeshException(
            "TEST", "Error in calculating the inverse offset 1" );
    }
    if( cartesiangrid.local_from_offset(
            cartesiangrid.cell_offset( sivec3{ 5, 2, 8 } ) )
        != sivec3{ 5, 2, 8 } )
    {
        throw RINGMeshException(
            "TEST", "Error in calculating the inverse offset 2" );
    }

    std::vector< float > values( cartesiangrid.nb_cells() );
    std::vector< index_t > uint_values( cartesiangrid.nb_cells() );
    for( auto i : range( cartesiangrid.nb_cells() ) )
    {
        values[i] = (float) i;
        uint_values[i] = i;
    }
    std::string name{ "index" };
    cartesiangrid.add_float_attribute( name, values );
    std::string name2{ "uint_index" };
    cartesiangrid.add_index_t_attribute( name2, uint_values );

    for( auto i : range( cartesiangrid.nb_cells_axis( 0 ) ) )
    {
        for( auto j : range( cartesiangrid.nb_cells_axis( 1 ) ) )
        {
            for( auto k : range( cartesiangrid.nb_cells_axis( 2 ) ) )
            {
                sivec3 position{ (signed_index_t) i, (signed_index_t) j,
                    (signed_index_t) k };
                index_t cell_offset = cartesiangrid.cell_offset( position );
                if( cartesiangrid.get_float_attribute_value( name, position )
                    != (float) cell_offset )
                {
                    throw RINGMeshException(
                        "Test", "Error in float attribute value #1." );
                }
                if( cartesiangrid.get_index_t_attribute_value( name2, position )
                    != cell_offset )
                {
                    throw RINGMeshException(
                        "Test", "Error in index_t attribute value #1." );
                }
            }
        }
    }

    GEO::Attribute< float > indices{ cartesiangrid.attributes_manager(), name };
    GEO::Attribute< index_t > index_t_indices{
        cartesiangrid.attributes_manager(), name2
    };
    for( auto i : range( cartesiangrid.nb_cells_axis( 0 ) ) )
    {
        for( auto j : range( cartesiangrid.nb_cells_axis( 1 ) ) )
        {
            for( auto k : range( cartesiangrid.nb_cells_axis( 2 ) ) )
            {
                sivec3 position{ (signed_index_t) i, (signed_index_t) j,
                    (signed_index_t) k };
                index_t cell_offset = cartesiangrid.cell_offset( position );
                if( indices[cell_offset] != (float) cell_offset )
                {
                    throw RINGMeshException(
                        "Test", "Error in float attribute value #2." );
                }
                if( index_t_indices[cell_offset] != cell_offset )
                {
                    throw RINGMeshException(
                        "Test", "Error in index_t attribute value #2." );
                }
            }
        }
    }

    GEO::Attribute< float >* indices_pointer =
        cartesiangrid.get_float_attribute( name );
    GEO::Attribute< index_t >* index_t_indices_pointer =
        cartesiangrid.get_index_t_attribute( name2 );
    for( auto i : range( cartesiangrid.nb_cells_axis( 0 ) ) )
    {
        for( auto j : range( cartesiangrid.nb_cells_axis( 1 ) ) )
        {
            for( auto k : range( cartesiangrid.nb_cells_axis( 2 ) ) )
            {
                sivec3 position{ (signed_index_t) i, (signed_index_t) j,
                    (signed_index_t) k };
                index_t cell_offset = cartesiangrid.cell_offset( position );
                if( indices_pointer->operator[]( cell_offset )
                    != (float) cell_offset )
                {
                    throw RINGMeshException(
                        "Test", "Error in float attribute value #3." );
                }
                if( index_t_indices_pointer->operator[]( cell_offset )
                    != cell_offset )
                {
                    throw RINGMeshException(
                        "Test", "Error in index_t attribute value #3." );
                }
            }
        }
    }
    delete indices_pointer;
    delete index_t_indices_pointer;

    CartesianGridBuilder3D cgbuilder{ cartesiangrid };
    cgbuilder.remove_section_from_cartesian_grid( 2, 4 );
    if( cartesiangrid.nb_cells_axis( 2 ) != 8u )
    {
        throw RINGMeshException( "Test remove section",
            "Wrong number of sections left : ",
            cartesiangrid.nb_cells_axis( 2 ), "instead of 8." );
    }
    for( auto i : range( cartesiangrid.nb_cells_axis( 0 ) ) )
    {
        for( auto j : range( cartesiangrid.nb_cells_axis( 1 ) ) )
        {
            for( auto k : range( 3 ) )
            {
                sivec3 position{ (signed_index_t) i, (signed_index_t) j,
                    (signed_index_t) k };
                index_t cell_offset = cartesiangrid.cell_offset( position );
                if( indices[cell_offset] != (float) cell_offset )
                {
                    throw RINGMeshException(
                        "Test remove section", "Error in attribute value #1." );
                }
            }
        }
    }
    for( auto i : range( cartesiangrid.nb_cells_axis( 0 ) ) )
    {
        for( auto j : range( cartesiangrid.nb_cells_axis( 1 ) ) )
        {
            for( auto k : range( 4, cartesiangrid.nb_cells_axis( 2 ) ) )
            {
                sivec3 position{ (signed_index_t) i, (signed_index_t) j,
                    (signed_index_t) k };
                index_t cell_offset = cartesiangrid.cell_offset( position );
                if( indices[cell_offset] != (float) ( cell_offset + 10 * 8 ) )
                {
                    throw RINGMeshException(
                        "Test remove section", "Error in attribute value #2." );
                }
            }
        }
    }
}

int main()
{
    using namespace RINGMesh;

    try
    {
        Logger::out( "TEST", "Frames and Cartesian Grid" );

        test_frames();
        Logger::out( "TEST", "Frames work" );
        test_cartesian_grids();
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
