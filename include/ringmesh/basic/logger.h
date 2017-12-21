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

#pragma once

#include <ringmesh/basic/common.h>

#include <mutex>

#include <geogram/basic/logger.h>

/*!
 * @file Logger class declaration
 * @author Arnaud Botella
 */

namespace RINGMesh
{
    class basic_api Logger
    {
    public:
        static void div( const std::string& title )
        {
            std::lock_guard< std::mutex > locking( lock() );
            GEO::Logger::div( title );
        }

        template < typename... Args >
        static void out( const std::string& feature, const Args&... args )
        {
            std::lock_guard< std::mutex > locking( lock() );
            log( GEO::Logger::out( feature ), args... );
        }

        template < typename... Args >
        static void err( const std::string& feature, const Args&... args )
        {
            std::lock_guard< std::mutex > locking( lock() );
            log( GEO::Logger::err( feature ), args... );
        }

        template < typename... Args >
        static void warn( const std::string& feature, const Args&... args )
        {
            std::lock_guard< std::mutex > locking( lock() );
            log( GEO::Logger::warn( feature ), args... );
        }

        static GEO::Logger* instance()
        {
            return GEO::Logger::instance();
        }

    private:
        static void log( std::ostream& os )
        {
            os << std::endl;
        }

        template < class A0, class... Args >
        static void log( std::ostream& os, const A0& a0, const Args&... args )
        {
            os << a0;
            log( os, args... );
        }

        static std::mutex& lock()
        {
            static std::mutex lock;
            return lock;
        }
    };
    class basic_api ThreadSafeConsoleLogger : public GEO::ConsoleLogger
    {
        using base_class = GEO::ConsoleLogger;

    public:
        void div( const std::string& title )
        {
            std::lock_guard< std::mutex > lock( lock_ );
            base_class::div( title );
        }

        void out( const std::string& str )
        {
            std::lock_guard< std::mutex > lock( lock_ );
            base_class::out( str );
        }

        void warn( const std::string& str )
        {
            std::lock_guard< std::mutex > lock( lock_ );
            base_class::warn( str );
        }

        void err( const std::string& str )
        {
            std::lock_guard< std::mutex > lock( lock_ );
            base_class::err( str );
        }

        void status( const std::string& str )
        {
            std::lock_guard< std::mutex > lock( lock_ );
            base_class::status( str );
        }

    private:
        std::mutex lock_{};
    };
} // namespace RINGMesh

#define DEBUG( a ) RINGMesh::Logger::out( "Debug", #a, " = ", a )
