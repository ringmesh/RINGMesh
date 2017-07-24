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

#pragma once

#include <ringmesh/basic/common.h>

#include <typeindex>

/*!
 * @file  Factory class
 * @author Arnaud Botella
 */

namespace RINGMesh {

    /*!
     * Generic factory
     * Example of use with A the base class and B, C Inherited classes
     *      // Instantiate
     *      Factory< std::string, A > factory;
     *      // Registration
     *      factory.register< B >( "B" );                 // B constructor has no argument
     *      factory.register< C, int >( "C" );            // C constructor takes an int
     *      factory.register< C, double, double >( "C" ); // Another C constructor
     *      // Creation
     *      std::unique_ptr< A > c = factory.crete( "C", 2.1, 8.6 );
     */
    template< typename Key, typename BaseClass >
    class Factory {
        static_assert( std::has_virtual_destructor<BaseClass>::value, "BaseClass must have a virtual destructor" );
    public:
        template< typename DerivedClass, typename ... Args >
        void register_creator( const Key& key )
        {
            static_assert( std::is_base_of<BaseClass, DerivedClass>::value, "DerivedClass must be a subclass of BaseClass" );
            static_assert( std::is_constructible<DerivedClass, Args...>::value, "DerivedClass must be constructible with Args..." );
            creators_.emplace(
                CreatorKey { key, create_function_type_index< Args... >() },
                reinterpret_cast< CreateFunc< > >( create_function_impl<
                    DerivedClass, Args... > ) );
        }

        template< typename ... Args >
        std::unique_ptr< BaseClass > create(
            const Key& key,
            Args const&... args ) const
        {
            auto creator = creators_.find(
                { key, create_function_type_index< Args... >() } );
            if( creator != creators_.end() ) {
                return reinterpret_cast< CreateFunc< Args... > >( creator )(
                    std::forward<const Args&>( args )... );
            } else {
                return {};
            }
        }

    private:
        template< typename ... Args >
        static std::type_index create_function_type_index()
        {
            return {typeid( CreateFunc<Args...> )};
        }

        template< typename DerivedClass, typename ... Args >
        static std::unique_ptr< BaseClass > create_function_impl(
            Args const&... args )
        {
            return std::unique_ptr< BaseClass > { new DerivedClass {
                std::forward<const Args&>( args )... } };
        }

        template< typename ... Args >
        using CreateFunc = typename std::add_pointer< std::unique_ptr<BaseClass>( const Args&... ) >::type;
        using CreatorKey = std::pair< Key, std::type_index >;

        std::map< CreatorKey, CreateFunc<> > creators_;
    };
}
