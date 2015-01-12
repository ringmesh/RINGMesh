/*[
 * Association Scientifique pour la Geologie et ses Applications (ASGA)
 * Copyright (c) 1993-2013 ASGA. All Rights Reserved.
 *
 * This program is a Trade Secret of the ASGA and it is not to be:
 * - reproduced, published, or disclosed to other,
 * - distributed or displayed,
 * - used for purposes or on Sites other than described
 *   in the GOCAD Advancement Agreement,
 * without the prior written authorization of the ASGA. Licencee
 * agrees to attach or embed this Notice on all copies of the program,
 * including partial copies or modified versions thereof.
 ]*/

#include <grgmesh/boundary_model.h>
#include <grgmesh/io.h>

#include <geogram/mesh/mesh_io.h>
#include <geogram/basic/string.h>

#include <third_party/zip/zip.h>
#include <third_party/zip/unzip.h>

#include <iostream>
#include <fstream>
#include <string>

namespace GRGMesh {
    namespace GRGMeshIO {

        bool load_BoundaryModel_from_Model3D(
            const std::string& filename,
            BoundaryModel& model )
        {
            std::ifstream input( filename.c_str() ) ;
            if( !input ) {
                std::cout << "cannot open file:" << filename << std::endl ;
                return false ;
            }

            BoundaryModelBuilder builder( model ) ;
            builder.load_file( input ) ;
            return true ;
        }
        /// Save a \param[in] macro mesh in a .zip file which contains all the mesh file. Type of the export is
        /// determined by the extension given in \param[in] filename
        bool save_macro_mesh( const MacroMesh& mm, const std::string& filename )
        {
            zipFile zf = zipOpen( filename.c_str(), APPEND_STATUS_CREATE ) ;
            zip_fileinfo zfi = { 0 };
            for( index_t i = 0; i < mm.nb_meshes(); i++ ) {
                GEO::MeshIOFlags flags ;
                flags.set_element( GEO::MESH_CELLS ) ;
                const GEO::Mesh& m = mm.mesh( i ) ;
                std::string name = GEO::String::to_string( i ) + ".meshb" ;
                GEO::mesh_save( m, name.c_str(), flags ) ;
                std::fstream file(name.c_str(), std::ios::binary | std::ios::in);
                file.seekg(0, std::ios::beg);
                index_t size = file.tellg();
                std::vector<char> buffer(size);
                zipOpenNewFileInZip( zf, GEO::String::to_string( i ).c_str(), &zfi,
                    NULL, 0, NULL, 0, NULL, Z_DEFLATED, Z_DEFAULT_COMPRESSION) ;
                zipWriteInFileInZip( zf, size == 0 ? "" : &buffer[0], size) ;
                zipCloseFileInZip(file) ;
                file.close() ;
            }
            zipClose(zf, NULL) ;
            return true ;
        }
    }
}
