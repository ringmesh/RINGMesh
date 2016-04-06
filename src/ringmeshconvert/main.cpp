/*
 * Copyright (c) 2012-2016, Association Scientifique pour la Geologie et ses Applications (ASGA)
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

#include <ringmesh/common.h>

#include <geogram/basic/command_line.h>
#include <geogram/basic/stopwatch.h>
#include <geogram/mesh/mesh_io.h>

#include <ringmesh/command_line.h>
#include <ringmesh/geo_model.h>
#include <ringmesh/geo_model_api.h>
#include <ringmesh/geo_model_builder.h>
#include <ringmesh/geo_model_validity.h>
#include <ringmesh/io.h>


/*!
 * @author Arnaud Botella
 */

namespace {
    using namespace RINGMesh ;


    /*!
     * @brief Write in the out stream things to save for CONTACT, INTERFACE and LAYERS
     */
    void save_high_level_bme( std::ofstream& out, const GeoModelElement& E )
    {
        /// First line:  TYPE - ID - NAME - GEOL
        out << E.gme_id() << " " ;
        if( E.has_name() ) {
            out << E.name() << " " ;
        } else {
            out << "no_name " ;
        }
        out << GeoModelElement::geol_name( E.geological_feature() ) << std::endl ;

        /// Second line:  IDS of children
        for( index_t j = 0; j < E.nb_children(); ++j ) {
            out << " " << E.child_id( j ).index ;
        }
        out << std::endl ;
    }

    void save_topology( const GeoModel& M, const std::string& file_name )
    {
        std::ofstream out( file_name.c_str() ) ;
        out.precision( 16 ) ;
        if( out.bad() ) {
            throw RINGMeshException( "I/O",
                "Error when opening the file: " + file_name ) ;
        }

        out << "RINGMESH BOUNDARY MODEL" << std::endl ;
        out << "NAME " << M.name() << std::endl ;

        // Numbers of the different types of elements
        for( index_t i = GME::CORNER; i < GME::NO_TYPE; i++ ) {
            GME::TYPE type = static_cast< GME::TYPE >( i ) ;
            out << "NB_" << GME::type_name( type ) << " " << M.nb_elements( type )
                << std::endl ;
        }
        // Write high-level elements
        for( index_t i = GME::CONTACT; i < GME::NO_TYPE; i++ ) {
            GME::TYPE type = static_cast< GME::TYPE >( i ) ;
            index_t nb = M.nb_elements( type ) ;
            for( index_t j = 0; j < nb; ++j ) {
                save_high_level_bme( out, M.element( GME::gme_t( type, j ) ) ) ;
            }
        }
        // Regions
        for( index_t i = 0; i < M.nb_regions(); ++i ) {
            const Region& E = M.region( i ) ;
            // Save ID - NAME
            out << E.gme_id() << " " ;
            if( E.has_name() ) {
                out << E.name() ;
            } else {
                out << "no_name" ;
            }
            out << std::endl ;
            // Second line Signed ids of boundary surfaces
            for( index_t j = 0; j < E.nb_boundaries(); ++j ) {
                if( E.side( j ) ) {
                    out << "+" ;
                } else {
                    out << "-" ;
                }
                out << E.boundary_gme( j ).index << " " ;
            }
            out << std::endl ;
        }

        // Universe
        out << "UNIVERSE " << std::endl ;
        for( index_t j = 0; j < M.universe().nb_boundaries(); ++j ) {
            if( M.universe().side( j ) ) {
                out << "+" ;
            } else {
                out << "-" ;
            }
            out << M.universe().boundary_gme( j ).index << " " ;
        }
        out << std::endl ;
    }


    /*!
     * @brief Save the GeoModel into a dedicated format bm
     * @todo Write the description of the BM format
     * @todo We need a generic read/write for the attributes !!
     */
    void save_bm_file( const GeoModel& M, const std::string& file_name )
    {
        save_topology( M, file_name ) ;
        std::ofstream out( file_name.c_str(), std::ios::out | std::ios::app ) ;

        out.precision( 16 ) ;

        // Corners
        for( index_t i = 0; i < M.nb_corners(); ++i ) {
            out << M.corner( i ).gme_id() << " " << M.corner( i ).vertex()
                << std::endl ;
        }
        // Lines
        for( index_t i = 0; i < M.nb_lines(); ++i ) {
            const Line& L = M.line( i ) ;
            out << L.gme_id() << std::endl ;
            out << "LINE_VERTICES " << L.nb_vertices() << std::endl ;
            for( index_t j = 0; j < L.nb_vertices(); ++j ) {
                out << L.vertex( j ) << std::endl ;
            }
            out << "IN_BOUNDARY " ;
            for( index_t j = 0; j < L.nb_in_boundary(); ++j ) {
                out << L.in_boundary_gme( j ).index << " " ;
            }
            out << std::endl ;
        }

        // Surfaces
        for( index_t i = 0; i < M.nb_surfaces(); ++i ) {
            const Surface& S = M.surface( i ) ;
            out << S.gme_id() << std::endl ;
            out << "SURFACE_VERTICES " << S.nb_vertices() << std::endl ;
            for( index_t j = 0; j < S.nb_vertices(); ++j ) {
                out << S.vertex( j ) << std::endl ;
            }

            out << "SURFACE_CORNERS " << S.nb_facet_corners() << std::endl ;
            out << "SURFACE_FACETS " << S.nb_cells() << std::endl ;
            for( index_t j = 0; j < S.nb_cells(); ++j ) {
                out << S.nb_vertices_in_facet( j ) << " " ;
                for( index_t v = 0; v < S.nb_vertices_in_facet( j ); ++v ) {
                    out << S.surf_vertex_id( j, v ) << " " ;
                }
                out << std::endl ;
            }
        }
    }

    class BMIOHandler: public GeoModelIOHandler {
    public:
        virtual void load( const std::string& filename, GeoModel& model )
        {
            std::ifstream input( filename.c_str() ) ;
            if( !input ) {
                throw RINGMeshException( "I/O", "Failed to open file " + filename ) ;
            }
            GeoModelBuilderBM builder( model, filename ) ;
            builder.build_model() ;
            GEO::Logger::out( "I/O" ) << " Loaded model " << model.name() << " from "
                << filename << std::endl ;
            print_geomodel( model ) ;
            is_geomodel_valid( model ) ;
        }

        virtual void save( const GeoModel& model, const std::string& filename )
        {
            save_bm_file( model, filename ) ;
        }
    } ;

#ifdef MINIZIP_FIXED
    class MMIOHandler: public GeoModelIOHandler {
    public:
        virtual void load( const std::string& filename, GeoModel& gm )
        {
            unzFile uz = unzOpen( filename.c_str() ) ;
            unz_global_info global_info ;
            if( unzGetGlobalInfo( uz, &global_info ) != UNZ_OK ) {
                unzClose( uz ) ;
                throw RINGMeshException( "ZLIB",
                    "Could not read file global info" ) ;
            }
            for( index_t r = 0; r < gm.nb_regions(); r++ ) {
                char filename[MAX_FILENAME] ;
                unzip_file( uz, filename ) ;
                GEO::MeshIOFlags flags ;
                flags.set_element( GEO::MESH_FACETS ) ;
                flags.set_element( GEO::MESH_CELLS ) ;
                flags.set_element( GEO::MESH_EDGES ) ;
                flags.set_attribute( GEO::MESH_FACET_REGION ) ;
                GEO::Mesh& m = gm.region( r ).mesh() ;
                std::string ext = GEO::FileSystem::extension( filename ) ;
                if( ext == "meshb" ) {
                    GEO::Logger::instance()->set_minimal( true ) ;
                    GEO::mesh_load( GEO::String::to_string( filename ), m, flags ) ;
                    GEO::Logger::instance()->set_minimal( false ) ;
                } else {
                    ringmesh_assert_not_reached;
                }
                GEO::FileSystem::delete_file( filename ) ;  // WHY ?? [Jeanne]

                if( ( r + 1 ) < global_info.number_entry ) {
                    if( unzGoToNextFile( uz ) != UNZ_OK ) {
                        unzClose( uz ) ;
                        throw RINGMeshException( "ZLIB",
                            "Could not read next file" ) ;
                    }
                }
            }
            unzClose( uz ) ;
        }

        /// Save a \param[in] gm macro mesh in a .zip file which contains all the mesh file. Type of the export is
        /// determined by the extension given in \param[in] filename
        virtual void save( const GeoModel& gm, const std::string& filename )
        {
            std::string pwd = GEO::FileSystem::get_current_working_directory() ;
            GEO::FileSystem::set_current_working_directory(
                GEO::FileSystem::dir_name( filename ) ) ;
            zipFile zf = zipOpen( filename.c_str(), APPEND_STATUS_CREATE ) ;
            for( index_t m = 0; m < gm.nb_regions(); m++ ) {
                GEO::MeshIOFlags flags ;
                flags.set_element( GEO::MESH_FACETS ) ;
                flags.set_element( GEO::MESH_CELLS ) ;
                flags.set_element( GEO::MESH_EDGES ) ;
                flags.set_attribute( GEO::MESH_FACET_REGION ) ;

                const GEO::Mesh& cur_mesh = gm.region( m ).mesh() ;
                std::string name_mesh_file = "region_" + GEO::String::to_string( m )
                    + ".meshb" ;

                GEO::Logger::instance()->set_quiet( true ) ;
                GEO::mesh_save( cur_mesh, name_mesh_file, flags ) ;
                GEO::Logger::instance()->set_quiet( false ) ;

                zip_file( zf, name_mesh_file ) ;

                GEO::FileSystem::delete_file( name_mesh_file ) ;

            }
            zipClose( zf, NULL ) ;
            GEO::FileSystem::set_current_working_directory( pwd ) ;
        }
    } ;
#endif

    /************************************************************************/
    /*!
     * Loads a GeoModel from a file
     * @param[in] filename the file to load
     * @param[out] model the model to fill
     */
    void geomodel_surface_load( const std::string& filename, GeoModel& model )
    {
        if( filename.empty() ) {
            throw RINGMeshException( "I/O",
                "No filename provided for structural model, use in:model" ) ;
        }
        std::ifstream input( filename.c_str() ) ;
        if( !input ) {
            throw RINGMeshException( "I/O", "Cannot open file: " + filename ) ;
        }

        GEO::Logger::out( "I/O" ) << "Loading file: " << filename << std::endl ;

        IOHandler_var handler = GeoModelIOHandler::get_handler( filename ) ;
        handler->load( filename, model ) ;
    }

    /*!
     * Saves a GeoModel in a file
     * @param[in] model the model to save
     * @param[in] filename the filename where to save it
     */
    void geomodel_surface_save( const GeoModel& model, const std::string& filename )
    {
        GEO::Logger::out( "I/O" ) << "Saving file " << filename << std::endl ;

        IOHandler_var handler = GeoModelIOHandler::get_handler( filename ) ;
        handler->save( model, filename ) ;
    }


    /*!
     * Loads a GeoModel from a file
     * @param[in] filename the file to load
     * @param][out] model the mesh to fill
     */
    void geomodel_volume_load( const std::string& filename, GeoModel& model )
    {
        GEO::Logger::out( "I/O" ) << "Loading file " << filename << "..."
            << std::endl ;

        IOHandler_var handler = GeoModelIOHandler::get_handler( filename ) ;
        handler->load( filename, model ) ;
    }

    /*!
     * Saves a GeoModel in a file
     * @param[in] model the mesh to save
     * @param[in] filename the file where to save
     */
    void geomodel_volume_save( const GeoModel& model, const std::string& filename )
    {
        GEO::Logger::out( "I/O" ) << "Saving file " << filename << "..."
            << std::endl ;

        IOHandler_var handler = GeoModelIOHandler::get_handler( filename ) ;
        handler->save( model, filename ) ;
    }


}

namespace RINGMesh {

    namespace CmdLine {
        void import_temp_in_out()
        {
            GEO::CmdLine::declare_arg( "in:model", "",
                "Filename of the input structural model" ) ;
            GEO::CmdLine::declare_arg( "in:mesh", "",
                "Filename of the input volumetric mesh" ) ;
            GEO::CmdLine::declare_arg( "out:model", "",
                "Saves the structural model" ) ;
            GEO::CmdLine::declare_arg( "out:mesh", "",
                "Saves the volumetric mesh of the structural model" ) ;
            GEO::CmdLine::declare_arg( "in:old_geomodel", "",
                "Saves the volumetric mesh of the structural model" ) ;
            ringmesh_register_IOHandler_creator( BMIOHandler, "bm" ) ;
#ifdef MINIZIP_FIXED
            ringmesh_register_IOHandler_creator( MMIOHandler, "mm" );
#endif
        }
    }
}
int main( int argc, char** argv )
{
    using namespace RINGMesh ;

    try {

        GEO::initialize() ;
        configure_geogram() ;
        configure_ringmesh() ;

        GEO::Logger::div( "RINGMeshConvert" ) ;
        GEO::Logger::out( "" ) << "Welcome to RINGMeshConvert !" << std::endl ;
        GEO::Logger::out( "" ) << "People working on the project in RING"
            << std::endl ;
        GEO::Logger::out( "" ) << "Arnaud Botella <arnaud.botella@univ-lorraine.fr> "
            << std::endl ;

        CmdLine::import_arg_group( "in" ) ;
        CmdLine::import_arg_group( "out" ) ;
        CmdLine::import_temp_in_out() ;
        if( argc == 1 ) {
            GEO::CmdLine::show_usage() ;
            return 0 ;
        }

        std::vector< std::string > filenames ;
        if( !GEO::CmdLine::parse( argc, argv, filenames ) ) {
            return 1 ;
        }

        GEO::Stopwatch total( "Total time" ) ;

        GeoModel model_in ;

        std::string geomodel_in_name = GEO::CmdLine::get_arg( "in:geomodel" ) ;
        if( geomodel_in_name.empty() ) {
            geomodel_in_name = GEO::CmdLine::get_arg( "in:model" ) ;
            if( geomodel_in_name.empty() ) {
                throw RINGMeshException( "I/O",
                    "Give at least a filename in in:model or in:geomodel" ) ;
            } else {
                geomodel_surface_load( geomodel_in_name, model_in ) ;
            }
        } else {
            geomodel_load( model_in, geomodel_in_name ) ;
        }

        if( GEO::CmdLine::get_arg( "in:geomodel" ).empty() ) {
            std::string mesh_in_name = GEO::CmdLine::get_arg( "in:mesh" ) ;
            if( !mesh_in_name.empty() ) {
                geomodel_volume_load( mesh_in_name, model_in ) ;
            }
        }

        std::string geomodel_out_name = GEO::CmdLine::get_arg( "out:geomodel" ) ;
        if(geomodel_out_name.empty()) {
            std::string model_out_name = GEO::CmdLine::get_arg( "out:model" ) ;
            std::string mesh_out_name = GEO::CmdLine::get_arg( "out:mesh" ) ;
            if( !model_out_name.empty() ) {
                geomodel_surface_save( model_in, model_out_name ) ;
            }
            if( !mesh_out_name.empty() ) {
                 geomodel_volume_save( model_in, mesh_out_name ) ;
             }
        }
        else {
            geomodel_save(model_in,geomodel_out_name) ;
        }

    } catch( const RINGMeshException& e ) {
        GEO::Logger::err( e.category() ) << e.what() << std::endl ;
        return 1 ;
    } catch( const std::exception& e ) {
        GEO::Logger::err( "Exception" ) << e.what() << std::endl ;
        return 1 ;
    }
    return 0 ;
}
