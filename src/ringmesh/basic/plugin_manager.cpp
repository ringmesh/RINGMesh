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
    std::string plugin_directory( const std::string& executable_directory )
    {
        std::string parent_separator;
        for( auto i : RINGMesh::range( 3 ) )
        {
            ringmesh_unused( i );
            std::string cur_directory { executable_directory + parent_separator
                + "/lib/" };
            if( GEO::FileSystem::is_directory( cur_directory ) )
            {
                return cur_directory;
            }
            parent_separator += "/..";
        }
        return executable_directory;
    }

    std::vector< std::string > plugins;
} //namespace

#ifdef linux
#include <dlfcn.h>
#include <limits.h>
#include <unistd.h>

namespace RINGMesh
{

    class PluginManger::Impl
    {
    public:
        Impl()
            : plugin_directory_( plugin_directory( executable_directory() ) )
        {
            DEBUG( plugin_directory_ );
        }
        /*!
         * Loads the given library.
         * RTLD_NOW: All undefined symbols in the library are resolved before dlopen() returns.
         *          If this cannot be done, an error is returned.
         * RTLD_GLOBAL: The symbols defined by this library will be made available
         *          for symbol resolution of subsequently loaded libraries.
         */
        void load_library( const std::string& plugin_path ) const
        {
            void* plugin_handle = dlopen( plugin_path.c_str(),
            RTLD_NOW | RTLD_GLOBAL );
            if( plugin_handle == nullptr )
            {
                throw RINGMeshException( "Plugin", "Could not load ", plugin_path,
                    ": ", dlerror() );
            }
        }

        std::string find_library( const std::string& plugin_name ) const
        {
            std::string library_name { "lib" + plugin_name + ".so" };
            std::string library_path { plugin_directory_ + library_name };
            if( !GEO::FileSystem::is_file( library_path ) )
            {
                throw RINGMeshException( "Plugin", plugin_name, " not found" );
            }
            return library_path;
        }

    private:
        std::string executable_directory() const
        {
            char buff[PATH_MAX];
            ssize_t len = ::readlink( "/proc/self/exe", buff, sizeof( buff ) - 1 );
            if( len == -1 )
            {
                throw RINGMeshException( "PluginManager",
                    "Cannot find the location of the curent executable" );
            }
            buff[len] = '\0';
            return GEO::FileSystem::dir_name( std::string( buff ) );
        }
    private:
        const std::string plugin_directory_;
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

            std::string plugin_path = impl_->find_library( plugin_name );
            impl_->load_library( plugin_path );
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
