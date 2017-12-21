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

#include <geogram/basic/command_line.h>
#include <geogram/basic/file_system.h>
#include <geogram/basic/line_stream.h>

#include <ringmesh/basic/algorithm.h>
#include <ringmesh/basic/logger.h>
#include <ringmesh/basic/pimpl_impl.h>

/*!
 * @file PluginManager class declaration
 * @author Arnaud Botella
 */

namespace
{
    std::vector< std::string > plugins;

    bool read_plugins_configuration_file( const std::string& configuration_file )
    {
        try
        {
            GEO::LineInput file( configuration_file );
            ringmesh_assert( file.OK() );
            while( !file.eof() && file.get_line() )
            {
                file.get_fields();
                ringmesh_assert( file.nb_fields() == 1 );
                if( !RINGMesh::PluginManager::load_plugin( file.field( 0 ) ) )
                {
                    RINGMesh::Logger::err( "Plugin", "Failed to load ", file.field( 0 ) );
                    return false;
                }
            }
        }
        catch( const std::logic_error& ex )
        {
            RINGMesh::Logger::err( "Plugin", "Got an error: ", ex.what() );
            return false;
        }
        return true;
    }

    bool load_plugins_configuration( const std::string& configuration_directory )
    {
        auto config_file = configuration_directory + "/"
            + RINGMesh::PluginManager::configuration_file;
        if( std::ifstream{ config_file.c_str() }.good() )
        {
            return read_plugins_configuration_file( config_file );
        }
        return false;
    }
} //namespace

#if _WIN32
#include <KnownFolders.h>
#include <Shlobj.h>
#include <Windows.h>

namespace RINGMesh
{
    class PluginManager::Impl
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
                std::string message;
                LPTSTR errorText{ nullptr };
                FormatMessage(
                    // use system message tables to retrieve error text
                    FORMAT_MESSAGE_FROM_SYSTEM
                    // allocate buffer on local heap for error text
                    |FORMAT_MESSAGE_ALLOCATE_BUFFER
                    // Important! will fail otherwise, since we're not
                    // (and CANNOT) pass insertion parameters
                    |FORMAT_MESSAGE_IGNORE_INSERTS,
                    NULL,// unused with FORMAT_MESSAGE_FROM_SYSTEM
                    GetLastError(),
                    MAKELANGID(LANG_NEUTRAL, SUBLANG_DEFAULT),
                    (LPTSTR) &errorText,// output
                    0,// minimum size for output buffer
                    NULL);// arguments - see note

                if( errorText != nullptr )
                {
                    message = errorText;
                    LocalFree(errorText);
                }
                throw RINGMeshException( "Plugin", "Could not load ", plugin_path,
                    ": ", message );
            }
        }

        std::string home_directory() const
        {
            PWSTR path{ nullptr };
            HRESULT hr = SHGetKnownFolderPath(FOLDERID_Documents, 0, nullptr, &path);
            std::string result;
            if( SUCCEEDED( hr ) )
            {
                std::stringstream ss;
                ss << path;
                result = ss.str();
            }
            CoTaskMemFree(path);
            return result;
        }

        std::string running_directory() const
        {
            TCHAR path[MAX_PATH];
            GetCurrentDirectory( MAX_PATH, path );
            return path;
        }

    private:
        std::string library_name( const std::string& plugin_name ) const
        {
            return { plugin_name + ".dll" };
        }
    };
} // namespace RINGMesh
#else
#include <dlfcn.h>
#include <pwd.h>

namespace RINGMesh
{
    class PluginManager::Impl
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

        std::string home_directory() const
        {
            std::string homedir{ getenv( "HOME" ) };
            if( homedir.empty() )
            {
                struct passwd pd;
                struct passwd* pwdptr = &pd;
                struct passwd* tempPwdPtr;
                char pwdbuffer[200];
                auto pwdlinelen = sizeof( pwdbuffer );
                getpwuid_r( 22, pwdptr, pwdbuffer, pwdlinelen, &tempPwdPtr );
                homedir = pd.pw_dir;
            }
            return homedir;
        }

        std::string running_directory() const
        {
            char path[PATH_MAX];
            return ::getcwd( path, PATH_MAX );
        }

    private:
        std::string library_name( const std::string& plugin_name ) const
        {
#if __APPLE__
            return { "lib" + plugin_name + ".dylib" };
#else
            return { "lib" + plugin_name + ".so" };
#endif
        }
    };
} // namespace RINGMesh
#endif

namespace RINGMesh
{
    PImpl< PluginManager::Impl > PluginManager::impl_;
    const std::string PluginManager::configuration_file = std::string{ "RINGMesh.ini" };

    bool PluginManager::load_plugin( const std::string& plugin_name )
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

    bool PluginManager::load_plugins()
    {
        auto all_plugin_names = GEO::CmdLine::get_arg( "sys:plugins" );
        if( !all_plugin_names.empty() )
        {
            std::vector< std::string > plugin_names;
            GEO::String::split_string( all_plugin_names, ';', plugin_names );
            for( const auto& plugin_name : plugin_names )
            {
                if( !load_plugin( plugin_name ) )
                {
                    return false;
                }
            }
            return true;
        }
        return load_plugins_configuration( impl_->running_directory() )
            || load_plugins_configuration( impl_->home_directory() );
    }
} // namespace RINGMesh
