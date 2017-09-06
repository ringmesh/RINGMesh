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

#include <ringmesh/io/zip_file.h>

#include <minizip/unzip.h>
#include <minizip/zip.h>

#include <geogram/basic/file_system.h>

#include <ringmesh/basic/pimpl_impl.h>

/*!
 * @file Manage zip files
 * @author Arnaud Botella
 */

namespace RINGMesh {

    class ZipFile::Impl {
    public:
        Impl( const std::string& filename )
            : zip_directory_ { GEO::FileSystem::dir_name( filename ) }
        {
            std::string pwd { set_working_directory() };
            zip_file_ = zipOpen(
                GEO::FileSystem::base_name( filename, false ).c_str(),
                APPEND_STATUS_CREATE );
            ringmesh_assert( zip_file_ != nil );
            unset_working_directory( pwd );
        }

        ~Impl()
        {
            zipClose( zip_file_, NULL );
        }

        void add_file( const std::string& filename )
        {
            std::string pwd { set_working_directory() };
            std::fstream file( filename.c_str(), std::ios::in | std::ios::binary );
            file.seekg( 0, std::ios::end );
            auto size = static_cast< index_t >( file.tellg() );
            file.seekg( 0, std::ios::beg );
            std::vector< char > buffer( size );
            file.read( &buffer[0], size );
            zipOpenNewFileInZip( zip_file_, filename.c_str(), nullptr, nullptr, 0,
                nullptr, 0, nullptr, Z_DEFLATED, Z_DEFAULT_COMPRESSION );
            zipWriteInFileInZip( zip_file_, size == 0 ? "" : &buffer[0], size );
            zipCloseFileInZip( zip_file_ );
            file.close();
            unset_working_directory( pwd );
        }

    private:
        std::string set_working_directory()
        {
            std::string pwd { GEO::FileSystem::get_current_working_directory() };
            bool valid_new_working_directory {
                GEO::FileSystem::set_current_working_directory( zip_directory_ ) };
            if( !valid_new_working_directory ) {
                throw RINGMeshException( "ZipFile",
                    "Output directory does not exist" );
            }
            return pwd;
        }

        void unset_working_directory( const std::string& directory )
        {
            GEO::FileSystem::set_current_working_directory( directory );
        }

    private:
        std::string zip_directory_;
        zipFile zip_file_ { nullptr };
    };

    ZipFile::ZipFile( const std::string& filename )
        : impl_ { filename }
    {
    }

    void ZipFile::add_file( const std::string& filename )
    {
        impl_->add_file( filename );
    }

} // namespace RINGMesh
