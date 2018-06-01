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

    std::vector< float > fvalues( cartesiangrid.nb_cells() );
    std::vector< index_t > uvalues( cartesiangrid.nb_cells() );
    std::vector< int > ivalues( cartesiangrid.nb_cells() );
    std::vector< double > dvalues( cartesiangrid.nb_cells() );
    for( auto i : range( cartesiangrid.nb_cells() ) )
    {
        fvalues[i] = static_cast< float >( i );
        uvalues[i] = i;
        ivalues[i] = static_cast< int >( i );
        dvalues[i] = static_cast< double >( i );
    }
    std::string fname{ "f_index" };
    std::string uname{ "u_index" };
    std::string iname{ "i_index" };
    std::string dname{ "d_index" };
    cartesiangrid.add_attribute< float >( fname, fvalues );
    cartesiangrid.add_attribute< index_t >( uname, uvalues );
    cartesiangrid.add_attribute< int >( iname, ivalues );
    cartesiangrid.add_attribute< double >( dname, dvalues );

    for( auto i : range( cartesiangrid.nb_cells_axis( 0 ) ) )
    {
        for( auto j : range( cartesiangrid.nb_cells_axis( 1 ) ) )
        {
            for( auto k : range( cartesiangrid.nb_cells_axis( 2 ) ) )
            {
                sivec3 position{ (signed_index_t) i, (signed_index_t) j,
                    (signed_index_t) k };
                index_t cell_offset = cartesiangrid.cell_offset( position );
                if( cartesiangrid.get_attribute_value< float >(
                        fname, position )
                        != static_cast< float >( cell_offset )
                    || cartesiangrid.get_attribute_value< index_t >(
                           uname, position )
                           != cell_offset
                    || cartesiangrid.get_attribute_value< int >(
                           iname, position )
                           != static_cast< int >( cell_offset )
                    || cartesiangrid.get_attribute_value< double >(
                           dname, position )
                           != static_cast< double >( cell_offset ) )
                {
                    throw RINGMeshException(
                        "Test", "Error in attribute value #1." );
                }
            }
        }
    }

    GEO::Attribute< float > findices{ cartesiangrid.attributes_manager(),
        fname };
    GEO::Attribute< index_t > uindices{ cartesiangrid.attributes_manager(),
        uname };
    GEO::Attribute< int > iindices{ cartesiangrid.attributes_manager(), iname };
    GEO::Attribute< double > dindices{ cartesiangrid.attributes_manager(),
        dname };
    for( auto i : range( cartesiangrid.nb_cells_axis( 0 ) ) )
    {
        for( auto j : range( cartesiangrid.nb_cells_axis( 1 ) ) )
        {
            for( auto k : range( cartesiangrid.nb_cells_axis( 2 ) ) )
            {
                sivec3 position{ (signed_index_t) i, (signed_index_t) j,
                    (signed_index_t) k };
                index_t cell_offset = cartesiangrid.cell_offset( position );
                if( findices[cell_offset] != static_cast< float >( cell_offset )
                    || uindices[cell_offset] != cell_offset
                    || iindices[cell_offset]
                           != static_cast< int >( cell_offset )
                    || dindices[cell_offset]
                           != static_cast< double >( cell_offset ) )
                {
                    throw RINGMeshException(
                        "Test", "Error in attribute value #2." );
                }
            }
        }
    }

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
                if( uindices[cell_offset] != cell_offset )
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
                if( uindices[cell_offset] != ( cell_offset + 10 * 8 ) )
                {
                    throw RINGMeshException(
                        "Test remove section", "Error in attribute value #2." );
                }
            }
        }
    }
}

void test_cartesian_grid_exception()
{
    bool exceptions_are_thrown{ false };

    vec3 origin{ 100, 200, -4 };
    vec3 x{ 1, 2, -1 };
    vec3 y{ 2, 1, 4 };
    vec3 z{ 3, -2, -1 };

    Frame3D frame{ x, y, z };
    ReferenceFrame3D reference_frame{ origin, frame };
    ivec3 grid_dimensions{ 10, 8, 9 };
    CartesianGrid3D cartesiangrid{ grid_dimensions, reference_frame };

    std::vector< bool > values( cartesiangrid.nb_cells() - 1 );
    std::string name{ "index" };
    try
    {
        cartesiangrid.add_attribute< bool >( name, values );
    }
    catch( const RINGMeshException& e )
    {
        exceptions_are_thrown = true;
        Logger::out( "Test", "Exception \"", e.what(), "\" is well thrown." );
    }
    if( !exceptions_are_thrown )
    {
        throw RINGMeshException(
            "Test", "Exception attribute vector size is not thrown" );
    }
    else
    {
        exceptions_are_thrown = false;
    }
    values.push_back( true );

    try
    {
        cartesiangrid.add_attribute< bool >( name, values );
    }
    catch( const RINGMeshException& e )
    {
        exceptions_are_thrown = true;
        Logger::out( "Test", "Exception \"", e.what(), "\" is well thrown." );
    }
    if( !exceptions_are_thrown )
    {
        throw RINGMeshException(
            "Test", "Exception attribute type is not thrown" );
    }
    else
    {
        exceptions_are_thrown = false;
    }

    std::vector< float > fvalues( cartesiangrid.nb_cells() );
    for( auto i : range( cartesiangrid.nb_cells() ) )
    {
        fvalues[i] = static_cast< float >( i );
    }
    std::string fname{ "f_index" };
    cartesiangrid.add_attribute< float >( fname, fvalues );
    try
    {
        sivec3 position{ -1, 0, 0 };
        cartesiangrid.get_attribute_value< float >( fname, position );
    }
    catch( const RINGMeshException& e )
    {
        exceptions_are_thrown = true;
        Logger::out( "Test", "Exception \"", e.what(), "\" is well thrown." );
    }
    if( !exceptions_are_thrown )
    {
        throw RINGMeshException(
            "Test", "Exception out of frame is not thrown" );
    }
    else
    {
        exceptions_are_thrown = false;
    }
}

int main()
{
    using namespace RINGMesh;

    try
    {
        Logger::out( "TEST", "Frames and Cartesian Grid" );

        test_frames();
        Logger::out( "TEST", "Frames : OK" );
        test_cartesian_grids();
        Logger::out( "TEST", "Cartesian grids : OK" );
        test_cartesian_grid_exception();
        Logger::out( "TEST", "Grid exceptions : OK" );
    }
    catch( const std::exception& e )
    {
        Logger::err( "Exception", e.what() );
        return 1;
    }
    Logger::out( "TEST", "SUCCESS" );
    return 0;
}
