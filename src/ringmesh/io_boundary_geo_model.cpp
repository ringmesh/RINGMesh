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

#include <ringmesh/io.h>

#include <ctime>

#include <geogram/basic/algorithm.h>
#include <geogram/basic/file_system.h>
#include <geogram/basic/line_stream.h>
#include <geogram/basic/logger.h>
#include <geogram/mesh/mesh_io.h>
#include <geogram/basic/attributes.h>
#include <geogram/basic/memory.h>

#include <ringmesh/geo_model.h>
#include <ringmesh/geo_model_api.h>
#include <ringmesh/geo_model_builder.h>
#include <ringmesh/geo_model_element.h>
#include <ringmesh/geo_model_validity.h>
#include <ringmesh/ringmesh_config.h>
#include <ringmesh/utils.h>

/*!
 * @file Implementation of classes to load and save surface GeoModels meshes
 * @author various
 */

namespace {
    using namespace RINGMesh ;

    /*!
     * From a given file name (MyFile.ext), create a MyFile directory
     * in the directory containing that file, or the current working directory.
     * @param[in] filename the filename
     * @param[out] path full path to the created directory
     * @param[out] directory name of the created directory
     */
    void create_directory_from_filename(
        const std::string& filename,
        std::string& path,
        std::string& directory )
    {
        path = GEO::FileSystem::dir_name( filename ) ;
        directory = GEO::FileSystem::base_name( filename, true ) ;
        if( path == "." ) {
            path = GEO::FileSystem::get_current_working_directory() ;
        }
        path += "/" + directory ;
        GEO::FileSystem::create_directory( path ) ;
    }

    /*!
     * @brief Total number of facets in the Surfaces of a BM
     */
    inline index_t nb_facets( const GeoModel& BM )
    {
        index_t result = 0 ;
        for( index_t i = 0; i < BM.nb_surfaces(); ++i ) {
            result += BM.surface( i ).nb_cells() ;
        }
        return result ;
    }

    /*!
     * @brief Write a region information in a stream
     * @details Used by function to save the Model in a .ml file
     *
     * @param[in] count Region index in the file
     * @param[in] region The region to save
     * @param[in,out] out The file output stream
     */
    void save_region( index_t count, const Region& region, std::ostream& out )
    {
        out << "REGION " << count << "  " << region.name() << " " << std::endl ;
        index_t it = 0 ;

        for( index_t i = 0; i < region.nb_boundaries(); ++i ) {
            out << "  " ;
            if( region.side( i ) ) {
                out << "+" ;
            } else {
                out << "-" ;
            }
            out << region.boundary( i ).index() + 1 ;
            it++ ;
            if( it == 5 ) {
                out << std::endl ;
                it = 0 ;
            }
        }
        out << "  0" << std::endl ;
    }

    /*!
     * @brief Write information for on layer in a stream
     * @details Used by function to save the Model in a .ml file
     *
     * @param[in] count Index of the layer in the file
     * @param[in] offset Offset of region indices in the file
     * @param[in] layer The layer to write
     * @param[in,out] out The output file stream
     */
    void save_layer(
        index_t count,
        index_t offset,
        const GeoModelElement& layer,
        std::ostream& out )
    {
        out << "LAYER " << layer.name() << " " << std::endl ;
        index_t it = 0 ;

        for( index_t i = 0; i < layer.nb_children(); ++i ) {
            out << "  " << layer.child_id( i ).index + offset + 1 ;
            it++ ;
            if( it == 5 ) {
                out << std::endl ;
                it = 0 ;
            }
        }
        out << "  0" << std::endl ;
    }

    /*!
     * @brief Write basic header for Gocad coordinate system.
     * @param[in,out] out Output .ml file stream
     */
    void save_coordinate_system( std::ostream& out )
    {
        out << "GOCAD_ORIGINAL_COORDINATE_SYSTEM" << std::endl << "NAME Default"
            << std::endl << "AXIS_NAME \"X\" \"Y\" \"Z\"" << std::endl
            << "AXIS_UNIT \"m\" \"m\" \"m\"" << std::endl << "ZPOSITIVE Elevation"
            << std::endl << "END_ORIGINAL_COORDINATE_SYSTEM" << std::endl ;
    }

    /*!
     * @brief Check if the model can be saved in a skua-gocad .ml file
     * @details It assumes that the model is valid and verifies that:
     *   - all Interfaces have a name and geological feature
     *   - all Surfaces are in an Interface
     *   - all Surfaces are triangulated
     *   - all Regions have a name
     */
    bool check_gocad_validity( const GeoModel& M )
    {
        if( M.nb_interfaces() == 0 ) {
            GEO::Logger::err( "" ) << " The GeoModel " << M.name()
                << " has no Interface" << std::endl ;
            return false ;
        }
        for( index_t i = 0; i < M.nb_interfaces(); ++i ) {
            const GME& E = M.one_interface( i ) ;
            if( !E.has_name() ) {
                GEO::Logger::err( "" ) << E.gme_id() << " has no name" << std::endl ;
                return false ;
            }
            if( !E.has_geological_feature() ) {
                GEO::Logger::err( "" ) << E.gme_id() << " has no geological feature"
                    << std::endl ;
                return false ;
            }
        }
        for( index_t s = 0; s < M.nb_surfaces(); ++s ) {
            const Surface& S = M.surface( s ) ;
            if( !S.has_parent() ) {
                GEO::Logger::err( "" ) << S.gme_id()
                    << " does not belong to any Interface of the model"
                    << std::endl ;
                return false ;
            }
            if( !S.is_simplicial() ) {
                GEO::Logger::err( "" ) << S.gme_id() << " is not triangulated "
                    << std::endl ;
                return false ;
            }
        }
        for( index_t r = 0; r < M.nb_regions(); ++r ) {
            const Region& R = M.region( r ) ;
            if( !R.has_name() ) {
                GEO::Logger::err( "" ) << R.gme_id() << " has no name" << std::endl ;
                return false ;
            }
        }
        return true ;
    }

    /*! Brute force inefficient but I am debugging !!!! */
    bool has_surface_edge( const Surface& S, index_t v0_in, index_t v1_in )
    {
        for( index_t i = 0; i < S.nb_cells(); ++i ) {
            for( index_t j = 0; j < S.nb_vertices_in_facet( i ); ++j ) {
                index_t v0 = S.polytop_vertex_index( i, j ) ;
                index_t v1 = S.polytop_vertex_index( i, S.next_in_facet( i, j ) ) ;
                if( ( v0 == v0_in && v1 == v1_in ) || ( v0 == v1_in && v1 == v0_in ) ) {
                    return true ;
                }
            }
        }
        return false ;
    }

    /*!
     * @brief Save the model in a .ml file if it can
     * @param[in] M the model to save
     * @param[in,out] out Output file stream
     */
    void save_gocad_model3d( const GeoModel& M, std::ostream& out )
    {
        if( !check_gocad_validity( M ) ) {
            throw RINGMeshException( "I/O",
                " The GeoModel " + M.name() + +" cannot be saved in .ml format" ) ;
        }
        out.precision( 16 ) ;

        // Gocad Model3d headers
        out << "GOCAD Model3d 1" << std::endl << "HEADER {" << std::endl << "name: "
            << M.name() << std::endl << "}" << std::endl ;

        save_coordinate_system( out ) ;

        // Gocad::TSurf = RINGMesh::Interface
        for( index_t i = 0; i < M.nb_interfaces(); ++i ) {
            out << "TSURF " << M.one_interface( i ).name() << std::endl ;
        }

        index_t count = 1 ;

        // Gocad::TFace = RINGMesh::Surface
        for( index_t i = 0; i < M.nb_surfaces(); ++i ) {
            const Surface& s = M.surface( i ) ;
            out << "TFACE " << count << "  " ;
            out << GME::geol_name( s.geological_feature() ) ;
            out << " " << s.parent().name() << std::endl ;

            // Print the key facet which is the first three
            // vertices of the first facet
            out << "  " << s.vertex( 0, 0 ) << std::endl ;
            out << "  " << s.vertex( 0, 1 ) << std::endl ;
            out << "  " << s.vertex( 0, 2 ) << std::endl ;

            ++count ;
        }
        // Universe
        index_t offset_layer = count ;
        save_region( count, M.universe(), out ) ;
        ++count ;
        // Regions
        for( index_t i = 0; i < M.nb_regions(); ++i ) {
            save_region( count, M.region( i ), out ) ;
            ++count ;
        }
        // Layers
        for( index_t i = 0; i < M.nb_layers(); ++i ) {
            save_layer( count, offset_layer, M.layer( i ), out ) ;
            ++count ;
        }
        out << "END" << std::endl ;

        // Save the geometry of the Surfaces, Interface per Interface
        for( index_t i = 0; i < M.nb_interfaces(); ++i ) {
            const GME& tsurf = M.one_interface( i ) ;
            // TSurf beginning header
            out << "GOCAD TSurf 1" << std::endl << "HEADER {" << std::endl << "name:"
                << tsurf.name() << std::endl << "name_in_model_list:" << tsurf.name()
                << std::endl << "}" << std::endl ;
            save_coordinate_system( out ) ;

            out << "GEOLOGICAL_FEATURE " << tsurf.name() << std::endl
                << "GEOLOGICAL_TYPE " ;
            out << GME::geol_name( tsurf.geological_feature() ) ;
            out << std::endl ;
            out << "PROPERTY_CLASS_HEADER Z {" << std::endl << "is_z:on" << std::endl
                << "}" << std::endl ;

            index_t vertex_count = 1 ;
            // TFace vertex index = Surface vertex index + offset
            index_t offset = vertex_count ;

            // To collect Corners(BStones) indexes
            // and boundary (Line) first and second vertex indexes
            std::set< index_t > corners ;
            std::set< std::pair< index_t, index_t > > lineindices ;
            for( index_t j = 0; j < tsurf.nb_children(); ++j ) {
                offset = vertex_count ;
                const Surface& S = dynamic_cast< const Surface& >( tsurf.child( j ) ) ;

                out << "TFACE" << std::endl ;
                for( index_t k = 0; k < S.nb_vertices(); ++k ) {
                    out << "VRTX " << vertex_count << " " << S.vertex( k )
                        << std::endl ;
                    vertex_count++ ;
                }
                for( index_t k = 0; k < S.nb_cells(); ++k ) {
                    out << "TRGL " << S.polytop_vertex_index( k, 0 ) + offset << " "
                        << S.polytop_vertex_index( k, 1 ) + offset << " "
                        << S.polytop_vertex_index( k, 2 ) + offset << std::endl ;
                }
                for( index_t k = 0; k < S.nb_boundaries(); ++k ) {
                    const Line& L = dynamic_cast< const Line& >( S.boundary( k ) ) ;
                    index_t v0_model_id = L.model_vertex_id( 0 ) ;
                    index_t v1_model_id = L.model_vertex_id( 1 ) ;

                    std::vector< index_t > v0_surface_ids = S.gme_vertex_indices(
                        v0_model_id ) ;
                    std::vector< index_t > v1_surface_ids = S.gme_vertex_indices(
                        v1_model_id ) ;

                    if( !S.has_inside_border() ) {
                        ringmesh_assert(
                            v0_surface_ids.size() == 1
                                && v1_surface_ids.size() == 1 ) ;
                        index_t v0 = v0_surface_ids[0] ;
                        index_t v1 = v1_surface_ids[0] ;
                        v0 += offset ;
                        v1 += offset ;

                        lineindices.insert(
                            std::pair< index_t, index_t >( v0, v1 ) ) ;
                        corners.insert( v0 ) ;
                    } else {
                        // We need to get the right pair of v0 - v1  (not crossing the inside boundary)
                        // corner and a border
                        int count = 0 ;
                        bool to_break = false ;
                        for( index_t iv0 = 0; iv0 < v0_surface_ids.size(); ++iv0 ) {
                            index_t v0 = v0_surface_ids[iv0] ;
                            for( index_t iv1 = 0; iv1 < v1_surface_ids.size();
                                ++iv1 ) {
                                index_t v1 = v1_surface_ids[iv1] ;
                                if( has_surface_edge( S, v0, v1 ) ) {
                                    lineindices.insert(
                                        std::pair< index_t, index_t >( v0 + offset,
                                            v1 + offset ) ) ;
                                    count++ ;
                                }
                                if( !L.is_inside_border( S ) && count == 1 ) {
                                    to_break = true ;
                                    break ;
                                } else if( count == 2 ) {
                                    to_break = true ;
                                    break ;
                                }
                            }
                            if( to_break ) {
                                corners.insert( v0 + offset ) ;
                                break ;
                            }
                        }
                    }
                    // Set a BSTONE at the line other extremity
                    const Corner& c1 =
                        dynamic_cast< const Corner& >( L.boundary( 1 ) ) ;
                    corners.insert(
                        S.surf_vertex_id( c1.model_vertex_id() ) + offset ) ;
                }
            }
            // Add the remaining bstones that are not already in bstones
            for( std::set< index_t >::iterator it( corners.begin() );
                it != corners.end(); ++it ) {
                out << "BSTONE " << *it << std::endl ;
            }
            for( std::set< std::pair< index_t, index_t > >::iterator it(
                lineindices.begin() ); it != lineindices.end(); ++it ) {
                out << "BORDER " << vertex_count << " " << it->first << " "
                    << it->second << std::endl ;
                vertex_count++ ;
            }
            out << "END" << std::endl ;
        }
    }

    /*! To save the attributes in a Graphite readable file, we need to write the correct
     * keyword for the attribute type - We restrict ourselves to the 3 types
     * int          "integer"
     * double       "real"
     * float        "real"
     * bool         "boolean"
     */
    inline std::string alias_name( const std::string& in )
    {
        if( in == "int" ) {
            return "integer" ;
        } else if( in == "index" ) {
            return "integer" ;
        } else if( in == "double" ) {
            return "real" ;
        } else if( in == "float" ) {
            return "real" ;
        } else if( in == "bool" ) {
            return "boolean" ;
        }
        ringmesh_assert_not_reached;
        return "" ;
    }

    /*!
     * @brief Save the model in smesh format
     * @details No attributes and no boundary marker are transferred
     * @todo Test this function - Create the appropriate handler
     */
    void save_smesh_file( const GeoModel& M, const std::string& file_name )
    {
        std::ofstream out( file_name.c_str() ) ;
        if( out.bad() ) {
            GEO::Logger::err( "I/O" ) << "Error when opening the file: "
                << file_name.c_str() << std::endl ;
            return ;
        }
        out.precision( 16 ) ;

        /// 1. Write the unique vertices
        out << "# Node list" << std::endl ;
        out << "# node count, 3 dim, no attribute, no boundary marker" << std::endl ;
        out << M.mesh.vertices.nb() << " 3 0 0" << std::endl ;
        out << "# node index, node coordinates " << std::endl ;
        for( index_t p = 0; p < M.mesh.vertices.nb(); p++ ) {
            const vec3& V = M.mesh.vertices.vertex( p ) ;
            out << p << " " << " " << V.x << " " << V.y << " " << V.z << std::endl ;
        }

        /// 2. Write the triangles
        out << "# Part 2 - facet list" << std::endl ;
        out << "# facet count, no boundary marker" << std::endl ;
        out << nb_facets( M ) << "  0 " << std::endl ;

        for( index_t i = 0; i < M.nb_surfaces(); ++i ) {
            const Surface& S = M.surface( i ) ;
            for( index_t f = 0; f < S.nb_cells(); f++ ) {
                out << S.nb_vertices_in_facet( f ) << " " ;
                for( index_t v = 0; v < S.nb_vertices_in_facet( f ); v++ ) {
                    out << S.model_vertex_id( f, v ) << " " ;
                }
                out << std::endl ;
            }
        }

        // Do not forget the stupid zeros at the end of the file
        out << std::endl << "0" << std::endl << "0" << std::endl ;
    }

    /************************************************************************/

    class MLIOHandler: public GeoModelIOHandler {
    public:
        /*! Load a .ml (Gocad file)
         * @pre Filename is valid
         */
        virtual void load( const std::string& filename, GeoModel& model )
        {
            std::ifstream input( filename.c_str() ) ;
            if( !input ) {
                throw RINGMeshException( "I/O", "Failed to open file " + filename ) ;
            }
            GeoModelBuilderGocad builder( model, filename ) ;

            time_t start_load, end_load ;
            time( &start_load ) ;

            builder.build_model() ;
            print_geomodel( model ) ;
            // Check validity
            is_geomodel_valid( model ) ;

            time( &end_load ) ;
            GEO::Logger::out( "I/O" ) << " Loaded model " << model.name() << " from "
                << std::endl << filename << " timing: "
                << difftime( end_load, start_load ) << "sec" << std::endl ;
        }

        virtual void save( const GeoModel& model, const std::string& filename )
        {

            std::ofstream out( filename.c_str() ) ;
            save_gocad_model3d( model, out ) ;
        }
    } ;

    class UCDIOHandler: public GeoModelIOHandler {
    public:
        virtual void load( const std::string& filename, GeoModel& model )
        {
            throw RINGMeshException( "I/O",
                "Geological model loading of a from UCD mesh not yet implemented" ) ;
        }

        virtual void save( const GeoModel& model, const std::string& filename )
        {
            std::string full_path ;
            std::string directory ;
            create_directory_from_filename( filename, full_path, directory ) ;

            std::ostringstream oss_cmd ;
            oss_cmd << full_path << "/cmd.lgi" ;
            std::ofstream cmd( oss_cmd.str().c_str() ) ;
            cmd << "cmo/create/3DMesh" << std::endl ;
            for( index_t s = 0; s < model.nb_surfaces(); s++ ) {
                std::ostringstream oss ;
                oss << full_path << "/surface_" << s << ".inp" ;
                std::ofstream out( oss.str().c_str() ) ;
                out.precision( 16 ) ;

                const Surface& surface = model.surface( s ) ;
                out << surface.nb_vertices() << " " << surface.nb_cells() << " 0 0 0"
                    << std::endl ;
                for( index_t v = 0; v < surface.nb_vertices(); v++ ) {
                    out << v << " " << surface.vertex( v ) << std::endl ;
                }

                for( index_t f = 0; f < surface.nb_cells(); f++ ) {
                    out << f << " 0 tri" ;
                    for( index_t v = 0; v < surface.nb_vertices_in_facet( f );
                        v++ ) {
                        out << " " << surface.polytop_vertex_index( f, v ) ;
                    }
                    out << std::endl ;
                }

                cmd << "cmo/create/s_" << s << std::endl ;
                cmd << "read/surface_" << s << ".inp/s_" << s << std::endl ;
                cmd << "surface/surface_" << s << "/" ;
                if( surface.is_on_voi() ) {
                    cmd << "reflect" ;
                } else {
                    cmd << "interface" ;
                }
                cmd << "/sheet/s_" << s << std::endl ;
            }

            for( index_t r = 0; r < model.nb_regions(); r++ ) {
                const Region& region = model.region( r ) ;
                cmd << "region/" << region.name() << "/" ;
                std::string sep = "" ;
                for( index_t s = 0; s < region.nb_boundaries(); s++ ) {
                    cmd << sep << "&" << std::endl ;
                    if( region.side( s ) )
                        cmd << "g" ;
                    else
                        cmd << "l" ;
                    cmd << "e surface_" << region.boundary_gme( s ).index ;
                    sep = " and " ;
                }
                cmd << std::endl ;
                cmd << "mregion/mat" << r << "/" << region.name() << std::endl ;
            }

            cmd << "finish" << std::endl ;
        }
    } ;

}
/************************************************************************/
namespace RINGMesh {
    /*
     * Initializes the possible handlers for IO GeoModel files
     */
    void GeoModelIOHandler::initialize_boundary_geomodel_output()
    {
        ringmesh_register_IOHandler_creator( MLIOHandler, "ml" ) ;
        ringmesh_register_IOHandler_creator( UCDIOHandler, "inp" );
    }
}
