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

#include <vector>

#include <ringmesh/basic/attributes.h>

/*!
 * @author Arnaud Botella
 */

using namespace RINGMesh;

template < typename T >
void test_attributes()
{
    Logger::out( "Test", "Testing type ", typeid( T ).name() );
    AttributesManager manager;
    Attribute< T > attribute0( manager, "attribute" );
    manager.resize( 10 );

    if( attribute0.size() != 10 )
    {
        throw RINGMeshException( "Test", "Failed to resize attribute" );
    }

    Attribute< T > attribute1( manager, "attribute" );
    attribute1.set_value( 5, T( 1 ) );

    for( auto i : range( 10 ) )
    {
        if( i == 5 )
        {
            if( attribute0[i] != T( 1 ) )
            {
                throw RINGMeshException( "Test", "Failed to edit attribute" );
            }
        }
        else
        {
            if( attribute0[i] != T( 0 ) )
            {
                throw RINGMeshException( "Test", "Failed to edit attribute" );
            }
        }
    }
    attribute1.set_constant_value( T( 0 ) );
    if( attribute1.size() != 1 )
    {
        throw RINGMeshException(
            "Test", "Failed to create a constant attribut" );
    }
}

int main()
{
    using namespace RINGMesh;

    try
    {
        default_configure();
        Logger::out( "TEST", "Test attributes" );
        test_attributes< double >();
        test_attributes< bool >();
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
