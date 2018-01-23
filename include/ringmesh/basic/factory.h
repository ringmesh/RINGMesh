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

#pragma once

#include <ringmesh/basic/common.h>

#include <ringmesh/basic/logger.h>
#include <ringmesh/basic/singleton.h>

/*!
 * @file  Factory class
 * @author Arnaud Botella
 */

namespace RINGMesh
{
    /*!
     * Generic factory
     * Example of use with A the base class and B, C derived classes
     *      // Instantiation
     *      using MyFactory = Factory< std::string, A, int, double >;
     *      // Registration
     *      MyFactory::register_creator< B >( "B" );   // B constructor takes an
     * int and a double
     *      MyFactory::register_creator< C >( "C" );   // C constructor takes an
     * int and a double
     *      // Creation
     *      std::unique_ptr< A > c = MyFactory::create( "C", 2, 8.6 );
     */
    template < typename Key, typename BaseClass, typename... Args >
    class basic_api Factory : public Singleton
    {
        static_assert( std::has_virtual_destructor< BaseClass >::value,
            "BaseClass must have a virtual destructor" );

    public:
        template < typename DerivedClass >
        static void register_creator( const Key& key )
        {
            static_assert( std::is_base_of< BaseClass, DerivedClass >::value,
                "DerivedClass is not a subclass of BaseClass" );
            static_assert(
                std::is_constructible< DerivedClass, Args... >::value,
                "DerivedClass is not constructible with Args..." );
            auto& store = get_store();
            if( !store
                     .emplace(
                         key, Creator( create_function_impl< DerivedClass > ) )
                     .second )
            {
                Logger::warn(
                    "Factory", "Trying to register twice the same key" );
            }
        }

        static std::unique_ptr< BaseClass > create(
            const Key& key, const Args&... args )
        {
            auto& store = get_store();
            auto creator = store.find( key );
            if( creator != store.end() )
            {
                return creator->second(
                    std::forward< const Args& >( args )... );
            }
            return {};
        }

        static std::vector< Key > list_creators()
        {
            auto& store = get_store();
            std::vector< Key > creators;
            creators.reserve( store.size() );
            for( const auto& creator : store )
            {
                creators.emplace_back( creator.first );
            }
            return creators;
        }

        static bool has_creator( const Key& key )
        {
            auto& store = get_store();
            return store.find( key ) != store.end();
        }

        using Creator = typename std::add_pointer< std::unique_ptr< BaseClass >(
            const Args&... ) >::type;
        using FactoryStore = std::map< Key, Creator >;

    private:
        template < typename DerivedClass >
        static std::unique_ptr< BaseClass > create_function_impl(
            Args&&... args )
        {
            return std::unique_ptr< BaseClass >{ new DerivedClass{
                std::forward< Args >( args )... } };
        }

        static FactoryStore& get_store()
        {
            return Singleton::instance< Factory >().store_;
        }

    private:
        FactoryStore store_;
    };

} // namespace RINGMesh
