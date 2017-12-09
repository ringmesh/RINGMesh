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

#include <ringmesh/basic/plugin_manager.h>

#include <geogram/basic/file_system.h>

#include <ringmesh/basic/algorithm.h>
#include <ringmesh/basic/logger.h>
#include <ringmesh/basic/pimpl_impl.h>

/*!
 * @file PluginManger class declaration
 * @author Arnaud Botella
 *
 */

namespace
{
    std::vector< std::string > plugins;
} //namespace

#ifdef linux
#include <dlfcn.h>

namespace RINGMesh
{
    class PluginManger::Impl
    {
    public:
        /*!
         * Loads the given library.
         * RTLD_NOW: All undefined symbols in the library are resolved before dlopen() returns.
         *          If this cannot be done, an error is returned.
         * RTLD_GLOBAL: The symbols defined by this library will be made available
         *          for symbol resolution of subsequently loaded libraries.
         */
        void load_library( const std::string& plugin_name ) const
        {
            std::string library{ library_name( plugin_name ) };
            void* plugin_handle = dlopen( library.c_str(),
            RTLD_NOW | RTLD_GLOBAL );
            if( plugin_handle == nullptr )
            {
                throw RINGMeshException( "Plugin", "Could not load ", library,
                    ": ", dlerror() );
            }
        }

    private:
        std::string library_name( const std::string& plugin_name ) const
        {
            return { "lib" + plugin_name + ".so" };
        }
    };
} // namespace RINGMesh
#elif _WIN32
#include <Windows.h>

namespace RINGMesh
{
    class PluginManger::Impl
    {
    public:
        /*!
         * Loads the given library.
         */
        void load_library( const std::string& plugin_path ) const
        {
            void* plugin_handle = LoadLibrary( plugin_path.c_str() );
            if( plugin_handle == nullptr )
            {
                LPTSTR message;
                FormatMessage(
                    FORMAT_MESSAGE_ALLOCATE_BUFFER |
                    FORMAT_MESSAGE_FROM_SYSTEM |
                    FORMAT_MESSAGE_IGNORE_INSERTS,
                    NULL,
                    GetLastError(),
                    MAKELANGID(LANG_NEUTRAL, SUBLANG_DEFAULT), // Default language
                    &message,
                    0,
                    NULL
                );
                throw RINGMeshException( "Plugin", "Could not load ", plugin_path,
                    ": ", message );
            }
        }

    private:

        std::string library_name( const std::string& plugin_name ) const
        {
            return { plugin_name + ".dll" };
        }
    };
} // namespace RINGMesh
#endif

namespace RINGMesh
{
    PImpl< PluginManger::Impl > PluginManger::impl_;

    bool PluginManger::load_module( const std::string& plugin_name )
    {
        try
        {
            Logger::out( "Plugin", "Loading ", plugin_name );
            if( contains( plugins, plugin_name ) )
            {
                throw RINGMeshException( "Plugin", plugin_name, " already loaded" );
            }

            impl_->load_library( plugin_name );
            plugins.push_back( plugin_name );
        }
        catch( RINGMeshException& e )
        {
            Logger::err( e.category(), e.what() );
            return false;
        }
        return true;
    }

} // namespace RINGMesh
