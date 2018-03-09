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

#include <ringmesh/io/zip_file.h>

#include <fstream>

#include <minizip/unzip.h>
#include <minizip/zip.h>

#include <geogram/basic/file_system.h>

#include <ringmesh/basic/pimpl_impl.h>

#define MAX_FILENAME 512
#define READ_SIZE 8192

/*!
 * @file Manage zip files
 * @author Arnaud Botella
 */

namespace RINGMesh
{
    class ZipFile::Impl
    {
    public:
        Impl( const std::string& filename )
        {
            zip_file_ = zipOpen( filename.c_str(), APPEND_STATUS_CREATE );
            if( zip_file_ == nullptr )
            {
                throw RINGMeshException(
                    "ZipFile", "Could not read file ", filename );
            }
        }

        ~Impl()
        {
            zipClose( zip_file_, NULL );
        }

        void add_file( const std::string& filename )
        {
            std::fstream file(
                filename.c_str(), std::ios::in | std::ios::binary );
            file.seekg( 0, std::ios::end );
            auto size = static_cast< index_t >( file.tellg() );
            file.seekg( 0, std::ios::beg );
            std::vector< char > buffer( size );
            file.read( &buffer[0], size );
            zipOpenNewFileInZip( zip_file_, filename.c_str(), nullptr, nullptr,
                0, nullptr, 0, nullptr, Z_DEFLATED, Z_DEFAULT_COMPRESSION );
            zipWriteInFileInZip( zip_file_, size == 0 ? "" : &buffer[0], size );
            zipCloseFileInZip( zip_file_ );
            file.close();
        }

    private:
        zipFile zip_file_{ nullptr };
    };

    ZipFile::ZipFile( const std::string& filename ) : impl_{ filename } {}

    ZipFile::~ZipFile() {}

    void ZipFile::add_file( const std::string& filename )
    {
        impl_->add_file( filename );
    }

    class UnZipFile::Impl
    {
    public:
        Impl( const std::string& filename, std::string directory_to_unzip )
            : zip_file_{ unzOpen( filename.c_str() ) },
              directory_to_unzip_( std::move( directory_to_unzip ) )
        {
            if( zip_file_ == nullptr )
            {
                throw RINGMeshException(
                    "UnZipFile", "Could not read file ", filename );
            }
            if( !GEO::FileSystem::is_directory( directory_to_unzip_ ) )
            {
                GEO::FileSystem::create_directory( directory_to_unzip_ );
            }
        }

        ~Impl()
        {
            unzClose( zip_file_ );
        }

        std::string get_file( const std::string& filename )
        {
            unzLocateFile( zip_file_, filename.c_str(), 0 );
            return unzip_current_file( zip_file_, filename );
        }

        void start_extract()
        {
            if( unzGoToFirstFile( zip_file_ ) != UNZ_OK )
            {
                throw RINGMeshException(
                    "UnZipFile", "Unable to uncompress the first file" );
            }
        }

        std::string get_current_file()
        {
            return unzip_current_file(
                zip_file_, get_current_filename().c_str() );
        }

        std::string get_current_filename()
        {
            char char_file_name[MAX_FILENAME];
            if( unzGetCurrentFileInfo64( zip_file_, nullptr, char_file_name,
                    MAX_FILENAME, nullptr, 0, nullptr, 0 )
                != UNZ_OK )
            {
                throw RINGMeshException(
                    "UnZipFile", "Unable to get file name" );
            }
            return { char_file_name };
        }

        bool next_file()
        {
            return unzGoToNextFile( zip_file_ ) == UNZ_OK;
        }

    private:
        std::string unzip_current_file(
            unzFile uz, const std::string& filename )
        {
            char read_buffer[READ_SIZE];
            if( unzOpenCurrentFile( uz ) != UNZ_OK )
            {
                unzClose( uz );
                throw RINGMeshException( "UnZipFile", "Could not open file" );
            }
            const std::string unziped_file{ directory_to_unzip_ + "/"
                                            + filename };
            FILE* out{ fopen( unziped_file.c_str(), "wb" ) };
            if( out == nullptr )
            {
                unzCloseCurrentFile( uz );
                unzClose( uz );
                throw RINGMeshException(
                    "UnZipFile", "Could not open destination file" );
            }
            int error{ UNZ_OK };
            do
            {
                error = unzReadCurrentFile( uz, read_buffer, READ_SIZE );
                if( error < 0 )
                {
                    unzCloseCurrentFile( uz );
                    unzClose( uz );
                    fclose( out );
                    throw RINGMeshException(
                        "UnZipFile", "Invalid error: ", error );
                }
                if( error > 0 )
                {
                    fwrite( read_buffer, static_cast< std::size_t >( error ),
                        std::size_t( 1 ), out );
                }
            } while( error > 0 );
            fclose( out );
            unzCloseCurrentFile( uz );
            return unziped_file;
        }

    private:
        unzFile zip_file_{ nullptr };
        const std::string directory_to_unzip_;
    };

    UnZipFile::UnZipFile(
        const std::string& filename, std::string directory_to_unzip )
        : impl_{ filename, std::move( directory_to_unzip ) }
    {
    }

    UnZipFile::~UnZipFile() {}

    std::string UnZipFile::get_file( const std::string& filename )
    {
        return impl_->get_file( filename );
    }

    void UnZipFile::start_extract()
    {
        impl_->start_extract();
    }

    std::string UnZipFile::get_current_file()
    {
        return impl_->get_current_file();
    }

    std::string UnZipFile::get_current_filename()
    {
        return impl_->get_current_filename();
    }

    bool UnZipFile::next_file()
    {
        return impl_->next_file();
    }

} // namespace RINGMesh
