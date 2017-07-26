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

#include <memory>

#include <ringmesh/basic/factory.h>

/*!
 * @author Arnaud Botella
 */

using namespace RINGMesh;

class ___A___ {

};

class ___B___ {

};

class ___Base___ {
public:
    virtual ~___Base___() = default;
protected:
    ___Base___( ___A___ a, ___B___ b )
        : a_( a ), b_( b )
    {
    }

protected:
    ___A___ a_;
    ___B___ b_;
};

class ___Derived___: public ___Base___ {
public:
    ___Derived___( ___A___& a, ___B___& b )
    : ___Base___( a, b )
    {
    }
};

class ___Derived2___: public ___Base___ {
public:
    ___Derived2___( ___A___ a, ___B___ b )
        : ___Base___( a, b )
    {
    }
};

class ___Derived3___: public ___Base___ {
public:
    ___Derived3___( ___A___& a, ___B___ b )
    : ___Base___( a, b )
    {
    }
};

void verdict( bool is_not_instantiated, std::string name )
{
    if( is_not_instantiated ) {
        throw RINGMeshException( "TEST", "Failed to instantiate the ", name,
            " class" );
    }
}

int main()
{
    using namespace RINGMesh;

    try {
        default_configure();
        Logger::out( "TEST", "Test Factory" );

        Factory2< std::string, ___Base___ > factory;
        factory.register_creator< ___Derived___, ___A___ &, ___B___& >( "Derived" );
        factory.register_creator< ___Derived2___, ___A___, ___B___ >( "Derived2" );
        factory.register_creator< ___Derived3___, ___A___ &, ___B___ >( "Derived3" );

        ___A___ a;
        ___B___ b;
        auto D = factory.create( "Derived", a, b );
        verdict( !D, "Derived" );
        auto D2 = factory.create( "Derived2", ___A___(), ___B___() );
        verdict( !D2, "Derived2" );
        auto D3 = factory.create( "Derived3", a, ___B___() );
        verdict( !D3, "Derived3" );

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
