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
 *     * Neither the name of the <organization> nor the
 *       names of its contributors may be used to endorse or promote products
 *       derived from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL <COPYRIGHT HOLDER> BE LIABLE FOR ANY
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

#include <iomanip>
#include <stack>

#include <geogram/basic/command_line.h>
#include <geogram/basic/file_system.h>
#include <geogram/basic/line_stream.h>
#include <geogram/basic/logger.h>

#include <geogram/mesh/mesh_geometry.h>
#include <geogram/mesh/mesh_io.h>

#include <ringmesh/geo_model.h>
#include <ringmesh/geo_model_api.h>
#include <ringmesh/geo_model_builder_so.h>
#include <ringmesh/geo_model_validity.h>
#include <ringmesh/geogram_extension.h>
#include <ringmesh/geometry.h>
#include <ringmesh/well.h>
#include <ringmesh/utils.h>
#include <ringmesh/io.h>

/*!
 * @file Implementation of classes loading volumetric GeoModels
 * @author Arnaud Botella and Antoine Mazuyer
 */

namespace {
    using namespace RINGMesh ;

    static std::string TAB = "\t" ;
    static std::string SPACE = " " ;
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

    /*!
     * @brief Save the connectivity of a GeoModelin in a file
     * @param[in] M the GeoModel
     * @param[in] file_name the file name (lol)
     */
    void save_connectivity( const GeoModel& M, const std::string& file_name )
    {
        std::ofstream out( file_name.c_str() ) ;
        out.precision( 16 ) ;
        if( out.bad() ) {
            throw RINGMeshException( "I/O",
                "Error when opening the file: " + file_name ) ;
        }

        for( index_t t = GME::CORNER; t < GME::REGION; t++ ) {
            GME::TYPE type = static_cast< GME::TYPE >( t ) ;
            for( index_t e = 0; e < M.nb_elements( type ); e++ ) {
                const GeoModelMeshElement& cur_geo_model_element = M.mesh_element(
                    type, e ) ;
                out << "GME" << " " << cur_geo_model_element.type_name( type ) << " "
                    << e << std::endl ;
                for( index_t in_b = 0; in_b < cur_geo_model_element.nb_in_boundary();
                    in_b++ ) {
                    out << cur_geo_model_element.in_boundary_gme( in_b ).index
                        << " " ;
                }
                out << std::endl ;
            }

        }
    }

    /*!
     * @brief Save the topology of a GeoModelin a file
     * @param[in] M the GeoModel
     * @param[in] file_name the file name (lol)
     */
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

    bool save_mesh(
        const GeoModelMeshElement& geo_model_element_mesh,
        const std::string& name,
        const GEO::MeshIOFlags& flags )
    {
        if( geo_model_element_mesh.type() == GME::REGION ) {
            const Region& region = geo_model_element_mesh.model().region(
                geo_model_element_mesh.index() ) ;
            if( !region.is_meshed() ) {
                // a region is not necessary meshed.
                return false ;
            }
        }
        GEO::mesh_save( geo_model_element_mesh.mesh(), name, flags ) ;
        return true ;
    }

    /************************************************************************/
    class AsterIOHandler: public GeoModelIOHandler {
    public:
        virtual void load( const std::string& filename, GeoModel& geomodel )
        {
            throw RINGMeshException( "I/O",
                "Loading of a GeoModel from Code_Aster mesh not implemented yet" ) ;
        }
        virtual void save( const GeoModel& gm, const std::string& filename )
        {
            std::vector< index_t > vertex_exported_id( gm.mesh.vertices.nb(),
                NO_ID ) ;
            std::vector< index_t > atom_exported_id(
                gm.mesh.cells.nb_duplicated_vertices(), NO_ID ) ;
            std::ofstream out( filename.c_str() ) ;
            out.precision( 16 ) ;
            std::vector< bool > vertex_exported( gm.mesh.vertices.nb(), false ) ;
            std::vector< bool > atom_exported(
                gm.mesh.cells.nb_duplicated_vertices(), false ) ;

            index_t nb_vertices_exported = 0 ;
            index_t cur_cell = 0 ;

            const GeoModelMesh& mesh = gm.mesh ;
            /// 1. Write the vertices coordinates (with the duplicate ones)
            out << "COOR_3D" << std::endl ;
            for( index_t c = 0; c < mesh.cells.nb(); c++ ) {
                for( index_t v = 0; v < mesh.cells.nb_vertices( c ); v++ ) {
                    index_t vertex_id, atom_id ;
                    if( mesh.cells.is_corner_duplicated( c, v, atom_id ) ) {
                        vertex_id = mesh.cells.duplicated_vertex( atom_id ) ;
                        if( atom_exported[atom_id] ) continue ;
                        atom_exported[atom_id] = true ;
                        atom_exported_id[atom_id] = nb_vertices_exported ;
                    } else {
                        vertex_id = mesh.cells.vertex( c, v ) ;
                        if( vertex_exported[vertex_id] ) continue ;
                        vertex_exported[vertex_id] = true ;
                        vertex_exported_id[vertex_id] = nb_vertices_exported ;
                    }
                    out << "V" << nb_vertices_exported++ << " "
                        << mesh.vertices.vertex( vertex_id ) << std::endl ;
                }
            }
            out << "FINSF" << std::endl ;

            /// 2. Write tetrahedra
            /// @todo Review: what about other elements ? [AB]
            out << "TETRA4" << std::endl ;
            for( index_t r = 0; r < gm.nb_regions(); r++ ) {
                for( index_t c = 0; c < mesh.cells.nb_tet( r ); c++ ) {
                    out << "C" << cur_cell++ << " " ;
                    for( index_t v = 0; v < mesh.cells.nb_vertices( c ); v++ ) {
                        index_t atom_id ;
                        if( mesh.cells.is_corner_duplicated( c, v, atom_id ) ) {
                            out << "V" << atom_exported_id[atom_id] << " " ;
                        } else {
                            index_t vertex_id = mesh.cells.vertex( c, v ) ;
                            out << "V" << vertex_exported_id[vertex_id] << " " ;
                        }
                    }
                    out << std::endl ;
                }
            }
            out << "FINSF" << std::endl ;

            /// 3. Associate tetrahedra to each region
            cur_cell = 0 ;
            for( index_t r = 0; r < gm.nb_regions(); r++ ) {
                out << "GROUP_MA" << std::endl ;
                out << gm.region( r ).name() << std::endl ;
                for( index_t c = 0; c < mesh.cells.nb_tet( r ); c++ ) {
                    out << "C" << cur_cell++ << std::endl ;
                }
            }
            out << "FINSF" << std::endl ;

            /// 4. Write triangles
            out << "TRIA3" << std::endl ;
            for( index_t i = 0; i < gm.nb_interfaces(); i++ ) {
                const RINGMesh::GeoModelElement& interf = gm.one_interface( i ) ;
                for( index_t s = 0; s < interf.nb_children(); s++ ) {
                    index_t surface_id = interf.child_id( s ).index ;
                    for( index_t f = 0; f < mesh.facets.nb_triangle( surface_id );
                        f++ ) {
                        index_t facet_id = mesh.facets.triangle( surface_id, f ) ;
                        out << "F" << facet_id ;
                        for( index_t v = 0; v < mesh.facets.nb_vertices( facet_id );
                            v++ ) {
                            out << " V"
                                << vertex_exported_id[mesh.facets.vertex( facet_id,
                                    v )] ;
                        }
                        out << std::endl ;
                    }
                }
            }
            out << "FINSF" << std::endl ;

            /// 5. Associate triangles to each surface
            for( index_t i = 0; i < gm.nb_interfaces(); i++ ) {
                const RINGMesh::GeoModelElement& interf = gm.one_interface( i ) ;
                for( index_t s = 0; s < interf.nb_children(); s++ ) {
                    index_t surface_id = interf.child_id( s ).index ;
                    out << "GROUP_MA" << std::endl ;
                    out << interf.name() << std::endl ;
                    for( index_t f = 0; f < mesh.facets.nb_triangle( surface_id );
                        f++ ) {
                        index_t facet_id = mesh.facets.triangle( surface_id, f ) ;
                        out << "F" << facet_id ;
                        out << std::endl ;
                    }
                    out << "FINSF" << std::endl ;
                }
            }

            out << "FIN" << std::endl ;

            out.close() ;
        }
    } ;

    /************************************************************************/

    /*!
     * @brief Save the attribute of one GeomModelElement
     * @param[in] filename path to the text file
     * @param[in] mesh_element the MeshElement of the original GEO::Mesh
     */
    void save_mesh_element_attributes(
        const std::string& filename,
        const GEO::MeshSubElementsStore& mesh_element )
    {

        GEO::vector< std::string > attribute_names ;
        mesh_element.attributes().list_attribute_names( attribute_names ) ;
        for( index_t att = 0; att < mesh_element.attributes().nb(); att++ ) {

            //TODO: this is temps ! Waiting for clean attribute type manager
            //by Bruno, for the moment we only deal with double
            if( !is_attribute_a_double( mesh_element.attributes(),
                attribute_names[att] ) || attribute_names[att] == "point" ) {
                continue ;
            }
            std::ofstream out( filename.c_str(), std::ios::out | std::ios::app ) ;
            out.precision( 16 ) ;
            GEO::Attribute< double > cur_att( mesh_element.attributes(),
                attribute_names[att] ) ;
            index_t attribute_dim = cur_att.dimension() ;
            out << "# " << attribute_names[att] << " double " << attribute_dim
                << std::endl ;

            for( index_t el = 0; el < mesh_element.nb(); el++ ) {
                for( index_t dim = 0; dim < attribute_dim; dim++ ) {
                    out << cur_att[el * attribute_dim + dim] << " " ;
                }
                out << std::endl ;
            }

        }

    }

    /*!
     * @brief Save the attribute of all the GeomModelElement
     * @param[in] root_name it contains the GeoModelElement type and its id
     * @param[in] mesh the mesh on which the attribute are stored
     * @param[in] zf the zip file
     */
    void save_geo_mesh_attributes(
        const std::string& root_name,
        const GEO::Mesh& mesh,
        zipFile& zf )
    {
        std::string vertices_attributes_file_name = root_name
            + "_vertices_attributes.txt" ;
        std::string edges_attributes_file_name = root_name
            + "_edges_attributes.txt" ;
        std::string facets_attributes_file_name = root_name
            + "_facets_attributes.txt" ;
        std::string cells_attributes_file_name = root_name
            + "_cells_attributes.txt" ;

        save_mesh_element_attributes( vertices_attributes_file_name,
            mesh.vertices ) ;
        if( GEO::FileSystem::is_file( vertices_attributes_file_name ) ) {
            zip_file( zf, vertices_attributes_file_name ) ;
            GEO::FileSystem::delete_file( vertices_attributes_file_name ) ;
        }

        save_mesh_element_attributes( edges_attributes_file_name, mesh.edges ) ;
        if( GEO::FileSystem::is_file( edges_attributes_file_name ) ) {
            zip_file( zf, edges_attributes_file_name ) ;
            GEO::FileSystem::delete_file( edges_attributes_file_name ) ;
        }

        save_mesh_element_attributes( facets_attributes_file_name, mesh.facets ) ;
        if( GEO::FileSystem::is_file( facets_attributes_file_name ) ) {
            zip_file( zf, facets_attributes_file_name ) ;
            GEO::FileSystem::delete_file( facets_attributes_file_name ) ;
        }
        save_mesh_element_attributes( cells_attributes_file_name, mesh.cells ) ;
        if( GEO::FileSystem::is_file( cells_attributes_file_name ) ) {
            zip_file( zf, cells_attributes_file_name ) ;
            GEO::FileSystem::delete_file( cells_attributes_file_name ) ;
        }
    }

    /*!
     * @brief Save the GeoModelMeshElement in a meshb file
     * @param[in] geo_model_element_mesh the GeoModelMeshElement you want to save
     * @param[in] zf the zip file
     */
    void save_geo_model_mesh_element(
        const GeoModelMeshElement& geo_model_element_mesh,
        zipFile& zf )
    {
        GME::TYPE type = geo_model_element_mesh.type() ;

        std::string name ;
        build_string_for_geo_model_element_export( type,
            geo_model_element_mesh.index(), name ) ;
        //@TODO this will be removed and integrated in the build_string....
        name += ".geogram" ;
        GEO::Logger* logger = GEO::Logger::instance() ;
        logger->set_quiet( true ) ;
        GEO::MeshIOFlags flags ;
//        flags.set_attribute( GEO::MESH_ALL_ATTRIBUTES ) ;
        bool is_saved = save_mesh( geo_model_element_mesh, name, flags ) ;
        logger->set_quiet( false ) ;

        if( is_saved ) {
            zip_file( zf, name ) ;
            GEO::FileSystem::delete_file( name ) ;
        }

//        std::string root_name ;
//        std::string handler ;
//        char separator = '.' ;
//        GEO::String::split_string( name, separator, root_name, handler ) ;
//        save_geo_mesh_attributes( root_name, geo_model_element_mesh.mesh(), zf ) ;
    }

    class GeoModelHandler: public GeoModelIOHandler {
        virtual void load( const std::string& filename, GeoModel& model )
        {
            std::string pwd = GEO::FileSystem::get_current_working_directory() ;
            GEO::FileSystem::set_current_working_directory(
                GEO::FileSystem::dir_name( filename ) ) ;
            GeoModelBuilderGM builder( model,
                GEO::FileSystem::base_name( filename, false ) ) ;
            builder.build_model() ;
            GEO::Logger::out( "I/O" ) << " Loaded model " << model.name() << " from "
                << filename << std::endl ;
            print_geomodel( model ) ;
            is_geomodel_valid( model ) ;
            GEO::FileSystem::set_current_working_directory( pwd ) ;

        }
        virtual void save( const GeoModel& model, const std::string& filename )
        {
            const std::string pwd =
                GEO::FileSystem::get_current_working_directory() ;
            bool valid_new_working_directory =
                GEO::FileSystem::set_current_working_directory(
                    GEO::FileSystem::dir_name( filename ) ) ;
            if( !valid_new_working_directory ) {
                throw RINGMeshException( "I/O", "Output directory does not exist" ) ;
            }

            zipFile zf = zipOpen(
                GEO::FileSystem::base_name( filename, false ).c_str(),
                APPEND_STATUS_CREATE ) ;
            ringmesh_assert( zf != nil ) ;

            save_topology( model, "topology.txt" ) ;
            zip_file( zf, "topology.txt" ) ;
            GEO::FileSystem::delete_file( "topology.txt" ) ;

            save_connectivity( model, "connectivity.txt" ) ;
            zip_file( zf, "connectivity.txt" ) ;
            GEO::FileSystem::delete_file( "connectivity.txt" ) ;

            for( index_t t = GME::CORNER; t <= GME::REGION; t++ ) {
                GME::TYPE type = static_cast< GME::TYPE >( t ) ;
                for( index_t e = 0; e < model.nb_elements( type ); e++ ) {
                    save_geo_model_mesh_element( model.mesh_element( type, e ),
                        zf ) ;
                }
            }
            zipClose( zf, NULL ) ;
            GEO::FileSystem::set_current_working_directory( pwd ) ;
        }

    } ;

    /************************************************************************/

    class LMIOHandler: public GeoModelIOHandler {
    public:
        virtual void load( const std::string& filename, GeoModel& geomodel )
        {
            throw RINGMeshException( "I/O",
                "Loading of a GeoModel from a mesh not implemented yet" ) ;
        }
        virtual void save( const GeoModel& gm, const std::string& filename )
        {
            gm.mesh.edges.test_and_initialize() ;
            gm.mesh.facets.test_and_initialize() ;
            gm.mesh.cells.test_and_initialize() ;

            GEO::Mesh mesh( 3 ) ;
            gm.mesh.copy_mesh( mesh ) ;

            GEO::Logger::instance()->set_minimal( true ) ;
            GEO::mesh_save( mesh, filename ) ;
            GEO::Logger::instance()->set_minimal( false ) ;
        }
    } ;

    /************************************************************************/
    class TetGenIOHandler: public GeoModelIOHandler {
    public:
        virtual void load( const std::string& filename, GeoModel& geomodel )
        {
            throw RINGMeshException( "I/O",
                "Loading of a GeoModel from TetGen not implemented yet" ) ;
        }
        virtual void save( const GeoModel& gm, const std::string& filename )
        {
            std::string directory = GEO::FileSystem::dir_name( filename ) ;
            std::string file = GEO::FileSystem::base_name( filename ) ;

            std::ostringstream oss_node ;
            oss_node << directory << "/" << file << ".node" ;
            std::ofstream node( oss_node.str().c_str() ) ;
            node.precision( 16 ) ;

            const GeoModelMesh& mesh = gm.mesh ;
            node << mesh.vertices.nb() << " 3 0 0" << std::endl ;
            for( index_t v = 0; v < mesh.vertices.nb(); v++ ) {
                node << v << SPACE << mesh.vertices.vertex( v ) << std::endl ;
            }

            std::ostringstream oss_ele ;
            oss_ele << directory << "/" << file << ".ele" ;
            std::ofstream ele( oss_ele.str().c_str() ) ;
            std::ostringstream oss_neigh ;
            oss_neigh << directory << "/" << file << ".neigh" ;
            std::ofstream neigh( oss_neigh.str().c_str() ) ;

            ele << mesh.cells.nb() << " 4 1" << std::endl ;
            neigh << mesh.cells.nb() << " 4" << std::endl ;
            index_t nb_tet_exported = 0 ;
            for( index_t m = 0; m < gm.nb_regions(); m++ ) {
                for( index_t tet = 0; tet < mesh.cells.nb_tet( m ); tet++ ) {
                    index_t cell = mesh.cells.tet( m, tet ) ;
                    ele << nb_tet_exported << SPACE << mesh.cells.vertex( cell, 0 )
                        << SPACE << mesh.cells.vertex( cell, 1 ) << SPACE
                        << mesh.cells.vertex( cell, 2 ) << SPACE
                        << mesh.cells.vertex( cell, 3 ) << SPACE << m + 1
                        << std::endl ;
                    neigh << nb_tet_exported ;
                    for( index_t f = 0; f < mesh.cells.nb_facets( tet ); f++ ) {
                        neigh << SPACE ;
                        index_t adj = mesh.cells.adjacent( cell, f ) ;
                        if( adj == GEO::NO_CELL ) {
                            neigh << -1 ;
                        } else {
                            neigh << adj ;
                        }
                    }
                    neigh << std::endl ;
                    nb_tet_exported++ ;
                }
            }
        }
    } ;

    /************************************************************************/

    struct RINGMesh2VTK {
        index_t element_type ;
        index_t vertices[8] ;
    } ;

    static RINGMesh2VTK tet_descriptor_vtk = { 10,                  // type
        { 0, 1, 2, 3 }     // vertices
    } ;

    static RINGMesh2VTK hex_descriptor_vtk = { 12,                         // type
        { 0, 4, 5, 1, 2, 6, 7, 3 }     // vertices
    } ;

    static RINGMesh2VTK prism_descriptor_vtk = { 13,                     // type
        { 0, 2, 1, 3, 5, 4 }   // vertices
    } ;

    static RINGMesh2VTK pyramid_descriptor_vtk = { 14,                 // type
        { 0, 1, 2, 3, 4 }  // vertices
    } ;

    static RINGMesh2VTK* cell_type_to_cell_descriptor_vtk[4] = {
        &tet_descriptor_vtk, &hex_descriptor_vtk, &prism_descriptor_vtk,
        &pyramid_descriptor_vtk } ;

    class VTKIOHandler: public GeoModelIOHandler {
    public:
        virtual void load( const std::string& filename, GeoModel& geomodel )
        {
            throw RINGMeshException( "I/O",
                "Loading of a GeoModel from VTK not implemented yet" ) ;
        }
        virtual void save( const GeoModel& gm, const std::string& filename )
        {
            std::ofstream out( filename.c_str() ) ;
            out.precision( 16 ) ;

            out << "# vtk DataFile Version 2.0" << std::endl ;
            out << "Unstructured Grid" << std::endl ;
            out << "ASCII" << std::endl ;
            out << "DATASET UNSTRUCTURED_GRID" << std::endl ;

            const GeoModelMesh& mesh = gm.mesh ;
            out << "POINTS " << mesh.vertices.nb() << " double" << std::endl ;
            for( index_t v = 0; v < mesh.vertices.nb(); v++ ) {
                out << mesh.vertices.vertex( v ) << std::endl ;
            }
            out << std::endl ;

            index_t total_corners = ( 4 + 1 ) * mesh.cells.nb_tet()
                + ( 5 + 1 ) * mesh.cells.nb_pyramid()
                + ( 6 + 1 ) * mesh.cells.nb_prism()
                + ( 8 + 1 ) * mesh.cells.nb_hex() ;
            out << "CELLS " << mesh.cells.nb_cells() << SPACE << total_corners
                << std::endl ;
            for( index_t c = 0; c < mesh.cells.nb(); c++ ) {
                out << mesh.cells.nb_vertices( c ) ;
                const RINGMesh2VTK& descriptor =
                    *cell_type_to_cell_descriptor_vtk[mesh.cells.type( c )] ;
                for( index_t v = 0; v < mesh.cells.nb_vertices( c ); v++ ) {
                    index_t vertex_id = descriptor.vertices[v] ;
                    out << SPACE << mesh.cells.vertex( c, vertex_id ) ;
                }
                out << std::endl ;
            }

            out << "CELL_TYPES " << mesh.cells.nb() << std::endl ;
            for( index_t c = 0; c < mesh.cells.nb(); c++ ) {
                const RINGMesh2VTK& descriptor =
                    *cell_type_to_cell_descriptor_vtk[mesh.cells.type( c )] ;
                out << descriptor.element_type << std::endl ;
            }
            out << std::endl ;

            out << "CELL_DATA " << mesh.cells.nb() << std::endl ;
            out << "SCALARS region int 1" << std::endl ;
            out << "LOOKUP_TABLE default" << std::endl ;
            for( index_t c = 0; c < mesh.cells.nb(); c++ ) {
                out << mesh.cells.region( c ) << std::endl ;
            }
            out << std::endl ;
        }
    } ;

    /************************************************************************/

    static index_t cell_type_mfem[4] = { 4, 5, NO_ID, NO_ID } ;

    static index_t facet_type_mfem[3] = { 2, 3, NO_ID } ;

    class MFEMIOHandler: public GeoModelIOHandler {
    public:
        virtual void load( const std::string& filename, GeoModel& geomodel )
        {
            throw RINGMeshException( "I/O",
                "Loading of a GeoModel from VTK not implemented yet" ) ;
        }
        virtual void save( const GeoModel& gm, const std::string& filename )
        {
            const GeoModelMesh& geomodel_mesh = gm.mesh ;
            index_t nb_cells = geomodel_mesh.cells.nb() ;
            if( geomodel_mesh.cells.nb_tet() != nb_cells
                && geomodel_mesh.cells.nb_hex() != nb_cells ) {
                throw RINGMeshException( "I/O",
                    "Export to MFEM format works only with full tet or full hex format" ) ;
            }
            std::ofstream out( filename.c_str() ) ;
            out.precision( 16 ) ;

            write_header( geomodel_mesh, out ) ;
            write_cells( geomodel_mesh, out ) ;
            write_facets( geomodel_mesh, out ) ;
            write_vertices( geomodel_mesh, out ) ;

        }

        /*!
         * @brief Write the header for the MFEM mesh file
         * @param[in] geomodel_mesh the GeoModelMesh to be saved
         * @param[in] out the ofstream that wrote the MFEM mesh file
         */
        void write_header( const GeoModelMesh& geomodel_mesh, std::ofstream& out )
        {
            // MFEM mesh version
            out << "MFEM mesh v1.0" << std::endl ;
            out << std::endl ;

            // Dimension is always 3 in our case
            out << "dimension" << std::endl ;
            out << "3" << std::endl ;
            out << std::endl ;
        }

        /*!
         * @brief Write the cells for the MFEM mesh file (work only with
         * tetrahedra and hexaedra)
         * @details The structure of the MFEM file for cells is
         * [group_id] [cell_type] [id_vertex_0] [id_vertex_1} .....
         * cell_type is 4 for  tetrahedra and 5 for hexahedra.
         * group_id begin with 1
         * @param[in] geomodel_mesh the GeoModelMesh to be saved
         * @param[in] out the ofstream that wrote the MFEM mesh file
         */
        void write_cells( const GeoModelMesh& geomodel_mesh, std::ofstream& out )
        {
            index_t nb_cells = geomodel_mesh.cells.nb() ;
            out << "elements" << std::endl ;
            out << nb_cells << std::endl ;
            for( index_t c = 0; c < nb_cells; c++ ) {
                out << geomodel_mesh.cells.region( c ) << " " ;
                out << cell_type_mfem[geomodel_mesh.cells.type( c )] << " " ;
                for( index_t v = 0; v < geomodel_mesh.cells.nb_vertices( c ); v++ ) {
                    out << geomodel_mesh.cells.vertex( c, v ) << " " ;
                }
                out << std::endl ;
            }
            out << std::endl ;
        }

        /*!
         * @brief Write the facets for the MFEM mesh file (work only with
         * triangles and quads)
         * @details The structure of the MFEM file for facets is
         * [group_id] [facet_type] [id_vertex_0] [id_vertex_1} .....
         * facet_type is 2 for triangles and 3 for the quads
         * group_id is continuous with the groupd indexes of the cells
         * @param[in] geomodel_mesh the GeoModelMesh to be saved
         * @param[in] out the ofstream that wrote the MFEM mesh file
         */
        void write_facets( const GeoModelMesh& geomodel_mesh, std::ofstream& out )
        {
            index_t offset = geomodel_mesh.model().nb_regions() ;
            out << "boundary" << std::endl ;
            out << geomodel_mesh.facets.nb() << std::endl ;
            for( index_t f = 0; f < geomodel_mesh.facets.nb(); f++ ) {
                index_t not_used = 0 ;
                out << geomodel_mesh.facets.surface( f ) + offset + 1 << " " ;
                out << facet_type_mfem[geomodel_mesh.facets.type( f, not_used )]
                    << " " ;
                for( index_t v = 0; v < geomodel_mesh.facets.nb_vertices( f );
                    v++ ) {
                    out << geomodel_mesh.facets.vertex( f, v ) << " " ;
                }
                out << std::endl ;
            }
            out << std::endl ;
        }

        /*!
         * @brief Write the vertices for the MFEM mesh file
         * @details The structure of the MFEM file for vertices is
         * [x] [y] [z]
         * @param[in] geomodel_mesh the GeoModelMesh to be saved
         * @param[in] out the ofstream that wrote the MFEM mesh file
         */
        void write_vertices( const GeoModelMesh& geomodel_mesh, std::ofstream& out )
        {
            out << "vertices" << std::endl ;
            out << geomodel_mesh.vertices.nb() << std::endl ;
            out << "3" << std::endl ;
            for( index_t v = 0; v < geomodel_mesh.vertices.nb(); v++ ) {
                out << geomodel_mesh.vertices.vertex( v ) << std::endl ;
            }
        }

    } ;
    /************************************************************************/

    class TSolidIOHandler: public GeoModelIOHandler {
    public:
        virtual void load( const std::string& filename, GeoModel& model )
        {
            std::ifstream input( filename.c_str() ) ;
            if( input ) {
                GeoModelBuilderTSolid builder( model, filename ) ;

                time_t start_load, end_load ;
                time( &start_load ) ;

                builder.build_model() ;
                print_geomodel( model ) ;
                // Check boundary model validity
                RINGMesh::is_geomodel_valid( model ) ;

                time( &end_load ) ;

                GEO::Logger::out( "I/O" ) << " Loaded model " << model.name()
                    << " from " << std::endl << filename << " timing: "
                    << difftime( end_load, start_load ) << "sec" << std::endl ;
            } else {
                throw RINGMeshException( "I/O",
                    "Failed loading model from file " + filename ) ;
            }
        }
        virtual void save( const GeoModel& gm, const std::string& filename )
        {
            std::ofstream out( filename.c_str() ) ;
            out.precision( 16 ) ;

            const GeoModel& model = gm ;
            // Print Model3d headers
            out << "GOCAD TSolid 1" << std::endl << "HEADER {" << std::endl
                << "name:" << model.name() << std::endl << "}" << std::endl ;

            out << "GOCAD_ORIGINAL_COORDINATE_SYSTEM" << std::endl << "NAME Default"
                << std::endl << "AXIS_NAME \"X\" \"Y\" \"Z\"" << std::endl
                << "AXIS_UNIT \"m\" \"m\" \"m\"" << std::endl
                << "ZPOSITIVE Elevation" << std::endl
                << "END_ORIGINAL_COORDINATE_SYSTEM" << std::endl ;

            const GeoModelMesh& mesh = gm.mesh ;
            mesh.set_duplicate_mode( GeoModelMeshCells::ALL ) ;

            std::vector< bool > vertex_exported( mesh.vertices.nb(), false ) ;
            std::vector< bool > atom_exported( mesh.cells.nb_duplicated_vertices(),
                false ) ;
            std::vector< index_t > vertex_exported_id( mesh.vertices.nb(), NO_ID ) ;
            std::vector< index_t > atom_exported_id(
                mesh.cells.nb_duplicated_vertices(), NO_ID ) ;
            index_t nb_vertices_exported = 1 ;
            for( index_t r = 0; r < model.nb_regions(); r++ ) {
                const RINGMesh::Region& region = model.region( r ) ;
                out << "TVOLUME " << region.name() << std::endl ;

                // Export not duplicated vertices
                for( index_t c = 0; c < mesh.cells.nb(); c++ ) {
                    index_t cell = mesh.cells.cell( r, c ) ;
                    for( index_t v = 0; v < mesh.cells.nb_vertices( c ); v++ ) {
                        index_t atom_id ;
                        if( !mesh.cells.is_corner_duplicated( cell, v, atom_id ) ) {
                            index_t vertex_id = mesh.cells.vertex( cell, v ) ;
                            if( vertex_exported[vertex_id] ) continue ;
                            vertex_exported[vertex_id] = true ;
                            vertex_exported_id[vertex_id] = nb_vertices_exported ;
                            out << "VRTX " << nb_vertices_exported++ << " "
                                << mesh.vertices.vertex( vertex_id ) << std::endl ;
                        }
                    }
                }

                // Export duplicated vertices
                for( index_t c = 0; c < mesh.cells.nb(); c++ ) {
                    index_t cell = mesh.cells.cell( r, c ) ;
                    for( index_t v = 0; v < mesh.cells.nb_vertices( c ); v++ ) {
                        index_t atom_id ;
                        if( mesh.cells.is_corner_duplicated( cell, v, atom_id ) ) {
                            if( atom_exported[atom_id] ) continue ;
                            atom_exported[atom_id] = true ;
                            atom_exported_id[atom_id] = nb_vertices_exported ;
                            index_t vertex_id = mesh.cells.vertex( cell, v ) ;
                            out << "ATOM " << nb_vertices_exported++ << " "
                                << vertex_exported_id[vertex_id] << std::endl ;
                        }
                    }
                }

                // Mark if a boundary is ending in the region
                std::map< index_t, index_t > sides ;
                for( index_t s = 0; s < region.nb_boundaries(); s++ ) {
                    if( sides.count( region.boundary_gme( s ).index ) > 0 )
                        // a surface is encountered twice, it is ending in the region
                        sides[region.boundary_gme( s ).index] = 2 ;
                    else
                        sides[region.boundary_gme( s ).index] = region.side( s ) ;
                }

                GEO::Attribute< index_t > attribute( mesh.facet_attribute_manager(),
                    surface_att_name ) ;
                for( index_t c = 0; c < mesh.cells.nb(); c++ ) {
                    out << "TETRA" ;
                    for( index_t v = 0; v < mesh.cells.nb_vertices( c ); v++ ) {
                        index_t atom_id ;
                        if( !mesh.cells.is_corner_duplicated( c, v, atom_id ) ) {
                            index_t vertex_id = mesh.cells.vertex( c, v ) ;
                            out << " " << vertex_exported_id[vertex_id] ;
                        } else {
                            out << " " << atom_exported_id[atom_id] ;
                        }
                    }
                    out << std::endl ;
                    out << "# CTETRA " << region.name() ;
                    for( index_t f = 0; f < mesh.cells.nb_facets( c ); f++ ) {
                        out << " " ;
                        index_t facet = NO_ID ;
                        bool side ;
                        if( mesh.cells.is_cell_facet_on_surface( c, f, facet,
                            side ) ) {
                            index_t surface_id = mesh.facets.surface( facet ) ;
                            side ? out << "+" : out << "-" ;
                            out << model.surface( surface_id ).parent().name() ;
                        } else {
                            out << "none" ;
                        }
                    }
                    out << std::endl ;
                }
            }

            out << "MODEL" << std::endl ;
            int tface_count = 1 ;
            for( index_t i = 0; i < model.nb_interfaces(); i++ ) {
                const RINGMesh::GeoModelElement& interf = model.one_interface( i ) ;
                out << "SURFACE " << interf.name() << std::endl ;
                for( index_t s = 0; s < interf.nb_children(); s++ ) {
                    out << "TFACE " << tface_count++ << std::endl ;
                    index_t surface_id = interf.child_id( s ).index ;
                    out << "KEYVERTICES" ;
                    index_t key_facet_id = mesh.facets.facet( surface_id, 0 ) ;
                    for( index_t v = 0; v < mesh.facets.nb_vertices( key_facet_id );
                        v++ ) {
                        out << " "
                            << vertex_exported_id[mesh.facets.vertex( key_facet_id,
                                v )] ;
                    }
                    out << std::endl ;
                    for( index_t f = 0; f < mesh.facets.nb_facets( surface_id );
                        f++ ) {
                        index_t facet_id = mesh.facets.facet( surface_id, f ) ;
                        out << "TRGL" ;
                        for( index_t v = 0; v < mesh.facets.nb_vertices( facet_id );
                            v++ ) {
                            out << " "
                                << vertex_exported_id[mesh.facets.vertex( facet_id,
                                    v )] ;
                        }
                        out << std::endl ;
                    }
                }
            }

            for( index_t r = 0; r < model.nb_regions(); r++ ) {
                const RINGMesh::Region& region = model.region( r ) ;
                out << "MODEL_REGION " << region.name() << " " ;
                region.side( 0 ) ? out << "+" : out << "-" ;
                out << region.boundary_gme( 0 ).index + 1 << std::endl ;
            }

            out << "END" << std::endl ;
        }
    } ;

    /************************************************************************/

    struct RINGMesh2CSMP {
        index_t element_type ;
        index_t nb_vertices ;
        index_t vertices[8] ;
        index_t nb_facets ;
        index_t facet[6] ;
    } ;

    static RINGMesh2CSMP tet_descriptor = { 4,                  // type
        4,                  // nb vertices
        { 0, 1, 2, 3 },     // vertices
        4,                  // nb facets
        { 0, 1, 2, 3 }      // facets
    } ;

    static RINGMesh2CSMP hex_descriptor = { 6,                         // type
        8,                              // nb vertices
        { 0, 4, 5, 1, 2, 6, 7, 3 },     // vertices
        6,                              // nb facets
        { 2, 0, 5, 1, 4, 3 }            // facets
    } ;

    static RINGMesh2CSMP prism_descriptor = { 12,                     // type
        6,                      // nb vertices
        { 0, 1, 2, 3, 4, 5 },   // vertices
        5,                      // nb facets
        { 0, 2, 4, 3, 1 }       // facets
    } ;

    static RINGMesh2CSMP pyramid_descriptor = { 18,                 // type
        5,                  // nb vertices
        { 0, 1, 2, 3, 4 },  // vertices
        5,                  // nb facets
        { 1, 4, 3, 2, 0 }   // facets
    } ;

    static RINGMesh2CSMP* cell_type_to_cell_descriptor[4] = {
        &tet_descriptor, &hex_descriptor, &prism_descriptor, &pyramid_descriptor } ;

    class CSMPIOHandler: public GeoModelIOHandler {
    public:
        CSMPIOHandler()
        {
            clear() ;
        }

        virtual void load( const std::string& filename, GeoModel& geomodel )
        {
            throw RINGMeshException( "I/O",
                "Loading of a GeoModel from CSMP not implemented yet" ) ;
        }
        virtual void save( const GeoModel& gm, const std::string& filename )
        {
            initialize( gm ) ;

            std::string directory = GEO::FileSystem::dir_name( filename ) ;
            std::string file = GEO::FileSystem::base_name( filename ) ;

            std::ostringstream oss_ascii ;
            oss_ascii << directory << "/" << file << ".asc" ;
            std::ofstream ascii( oss_ascii.str().c_str() ) ;
            ascii.precision( 16 ) ;

            ascii << gm.name() << std::endl ;
            ascii << "Model generated from RINGMesh" << std::endl ;

            std::ostringstream oss_data ;
            oss_data << directory << "/" << file << ".dat" ;
            std::ofstream data( oss_data.str().c_str() ) ;
            data.precision( 16 ) ;

            std::ostringstream oss_regions ;
            oss_regions << directory << "/" << file << "-regions.txt" ;
            std::ofstream regions( oss_regions.str().c_str() ) ;
            regions << "'" << oss_regions.str() << std::endl ;
            regions << "no properties" << std::endl ;

            const GeoModelMesh& mesh = gm.mesh ;
            index_t count = 0 ;
            // Conversion from (X,Y,Z) to (X,Z,-Y)
            signed_index_t conversion_sign[3] = { 1, 1, -1 } ;
            index_t conversion_axis[3] = { 0, 2, 1 } ;
            data << mesh.vertices.nb() << " # PX, PY, PZ" << std::endl ;
            for( index_t dim = 0; dim < 3; dim++ ) {
                for( index_t v = 0; v < mesh.vertices.nb(); v++ ) {
                    data << " "
                        << conversion_sign[dim]
                            * mesh.vertices.vertex( v )[conversion_axis[dim]] ;
                    new_line( count, 5, data ) ;
                }
                reset_line( count, data ) ;
            }
            reset_line( count, data ) ;

            index_t nb_families = 0 ;
            std::vector< index_t > nb_triangle_interface( gm.nb_interfaces(), 0 ) ;
            std::vector< index_t > nb_quad_interface( gm.nb_interfaces(), 0 ) ;
            for( index_t i = 0; i < gm.nb_interfaces(); i++ ) {
                const GeoModelElement& interf = gm.one_interface( i ) ;
                for( index_t s = 0; s < interf.nb_children(); s++ ) {
                    index_t s_id = interf.child_id( s ).index ;
                    nb_triangle_interface[i] += mesh.facets.nb_triangle( s_id ) ;
                    nb_quad_interface[i] += mesh.facets.nb_quad( s_id ) ;
                }
                if( nb_triangle_interface[i] > 0 ) nb_families++ ;
                if( nb_quad_interface[i] > 0 ) nb_families++ ;
            }
            for( index_t r = 0; r < gm.nb_regions(); r++ ) {
                if( mesh.cells.nb_tet( r ) > 0 ) nb_families++ ;
                if( mesh.cells.nb_pyramid( r ) > 0 ) nb_families++ ;
                if( mesh.cells.nb_prism( r ) > 0 ) nb_families++ ;
                if( mesh.cells.nb_hex( r ) > 0 ) nb_families++ ;
            }
            if( gm.wells() ) nb_families += gm.wells()->nb_wells() ;

            ascii << nb_families << " # Number of families" << std::endl ;
            ascii << "# Object name" << TAB << "Element type" << TAB << "Material-ID"
                << TAB << "Number of elements" << std::endl ;
            for( index_t r = 0; r < gm.nb_regions(); r++ ) {
                const RINGMesh::GeoModelElement& region = gm.region( r ) ;
                regions << region.name() << std::endl ;
                std::string element_type[4] = {
                    "TETRA_4", "HEXA_8", "PENTA_6", "PYRA_5" } ;
                for( index_t type = GEO::MESH_TET; type < GEO::MESH_CONNECTOR;
                    type++ ) {
                    GEO::MeshCellType T = static_cast< GEO::MeshCellType >( type ) ;
                    if( mesh.cells.nb_cells( r, T ) > 0 ) {
                        ascii << region.name() << TAB << element_type[type] << TAB
                            << 0 << TAB << mesh.cells.nb_cells( r, T ) << std::endl ;
                    }
                }
            }
            for( index_t i = 0; i < gm.nb_interfaces(); i++ ) {
                regions << interface_name( i, gm ) << std::endl ;
                if( nb_triangle_interface[i] > 0 ) {
                    ascii << interface_name( i, gm ) << TAB << "TRI_3" << TAB << 0
                        << TAB << nb_triangle_interface[i] << std::endl ;
                }
                if( nb_quad_interface[i] > 0 ) {
                    ascii << interface_name( i, gm ) << TAB << "QUAD_4" << TAB << 0
                        << TAB << nb_quad_interface[i] << std::endl ;
                }
            }
            if( gm.wells() ) {
                for( index_t w = 0; w < gm.wells()->nb_wells(); w++ ) {
                    const Well& well = gm.wells()->well( w ) ;
                    regions << well.name() << std::endl ;
                    ascii << well.name() << TAB << "BAR_2" << TAB << 0 << TAB
                        << well.nb_edges() << std::endl ;
                }
            }

            data << "# PBFLAGS" << std::endl ;
            for( index_t p = 0; p < mesh.vertices.nb(); p++ ) {
                data << " " << std::setw( 3 ) << point_boundary( p ) ;
                new_line( count, 20, data ) ;
            }
            reset_line( count, data ) ;

            data << "# PBVALS" << std::endl ;
            for( index_t p = 0; p < mesh.vertices.nb(); p++ ) {
                data << " " << std::setw( 3 ) << 0 ;
                new_line( count, 20, data ) ;
            }
            reset_line( count, data ) ;

            index_t nb_total_elements = mesh.cells.nb_cells()
                + mesh.facets.nb_facets() + mesh.edges.nb_edges() ;
            data << nb_total_elements << " # PELEMENT" << std::endl ;
            for( index_t r = 0; r < gm.nb_regions(); r++ ) {
                index_t element_type[4] = { 4, 6, 12, 18 } ;
                for( index_t type = GEO::MESH_TET; type < GEO::MESH_CONNECTOR;
                    type++ ) {
                    GEO::MeshCellType T = static_cast< GEO::MeshCellType >( type ) ;
                    for( index_t el = 0; el < mesh.cells.nb_cells( r, T ); el++ ) {
                        data << " " << std::setw( 3 ) << element_type[type] ;
                        new_line( count, 20, data ) ;
                    }
                }
            }
            for( index_t i = 0; i < gm.nb_interfaces(); i++ ) {
                for( index_t el = 0; el < nb_triangle_interface[i]; el++ ) {
                    data << " " << std::setw( 3 ) << 8 ;
                    new_line( count, 20, data ) ;
                }
                for( index_t el = 0; el < nb_quad_interface[i]; el++ ) {
                    data << " " << std::setw( 3 ) << 14 ;
                    new_line( count, 20, data ) ;
                }
            }
            if( gm.wells() ) {
                for( index_t w = 0; w < gm.wells()->nb_wells(); w++ ) {
                    const Well& well = gm.wells()->well( w ) ;
                    for( index_t e = 0; e < well.nb_edges(); e++ ) {
                        data << " " << std::setw( 3 ) << 2 ;
                        new_line( count, 20, data ) ;
                    }
                }
            }
            reset_line( count, data ) ;

            ascii
                << "# now the elements which make up each object are listed in sequence"
                << std::endl ;
            index_t cur_cell = 0 ;
            for( index_t r = 0; r < gm.nb_regions(); r++ ) {
                const RINGMesh::GeoModelElement& region = gm.region( r ) ;
                std::string element_type[4] = {
                    "TETRA_4", "HEXA_8", "PENTA_6", "PYRA_5" } ;
                for( index_t type = GEO::MESH_TET; type < GEO::MESH_CONNECTOR;
                    type++ ) {
                    GEO::MeshCellType T = static_cast< GEO::MeshCellType >( type ) ;
                    if( mesh.cells.nb_cells( r, T ) > 0 ) {
                        ascii << region.name() << " " << element_type[type] << " "
                            << mesh.cells.nb_cells( r, T ) << std::endl ;
                        for( index_t el = 0; el < mesh.cells.nb_cells( r, T );
                            el++ ) {
                            ascii << cur_cell++ << " " ;
                            new_line( count, 10, ascii ) ;
                        }
                        reset_line( count, ascii ) ;
                    }
                }
            }
            for( index_t i = 0; i < gm.nb_interfaces(); i++ ) {
                if( nb_triangle_interface[i] > 0 ) {
                    ascii << interface_name( i, gm ) << " " << "TRI_3" << " "
                        << nb_triangle_interface[i] << std::endl ;
                    for( index_t el = 0; el < nb_triangle_interface[i]; el++ ) {
                        ascii << cur_cell++ << " " ;
                        new_line( count, 10, ascii ) ;
                    }
                    reset_line( count, ascii ) ;
                }
                if( nb_quad_interface[i] > 0 ) {
                    ascii << interface_name( i, gm ) << " " << "QUAD_4" << " "
                        << nb_quad_interface[i] << std::endl ;
                    for( index_t el = 0; el < nb_quad_interface[i]; el++ ) {
                        ascii << cur_cell++ << " " ;
                        new_line( count, 10, ascii ) ;
                    }
                    reset_line( count, ascii ) ;
                }
            }
            if( gm.wells() ) {
                for( index_t w = 0; w < gm.wells()->nb_wells(); w++ ) {
                    const Well& well = gm.wells()->well( w ) ;
                    ascii << well.name() << " " << "BAR_2" << " " << well.nb_edges()
                        << std::endl ;
                    for( index_t e = 0; e < well.nb_edges(); e++ ) {
                        ascii << cur_cell++ << " " ;
                        new_line( count, 10, ascii ) ;
                    }
                    reset_line( count, ascii ) ;
                }
            }

            index_t nb_plist = 3 * mesh.facets.nb_triangle()
                + 4 * mesh.facets.nb_quad() + 4 * mesh.cells.nb_tet()
                + 5 * mesh.cells.nb_pyramid() + 6 * mesh.cells.nb_prism()
                + 8 * mesh.cells.nb_hex() + 2 * mesh.edges.nb_edges() ;
            data << nb_plist << " # PLIST" << std::endl ;
            for( index_t r = 0; r < gm.nb_regions(); r++ ) {
                for( index_t type = GEO::MESH_TET; type < GEO::MESH_CONNECTOR;
                    type++ ) {
                    GEO::MeshCellType T = static_cast< GEO::MeshCellType >( type ) ;
                    RINGMesh2CSMP& descriptor = *cell_type_to_cell_descriptor[T] ;
                    for( index_t el = 0; el < mesh.cells.nb_cells( r, T ); el++ ) {
                        index_t cell = mesh.cells.cell( r, el, T ) ;
                        for( index_t p = 0; p < descriptor.nb_vertices; p++ ) {
                            index_t csmp_p = descriptor.vertices[p] ;
                            index_t vertex_id = mesh.cells.vertex( cell, csmp_p ) ;
                            data << " " << std::setw( 7 ) << vertex_id ;
                            new_line( count, 10, data ) ;
                        }
                    }
                }
            }
            for( index_t i = 0; i < gm.nb_interfaces(); i++ ) {
                const GeoModelElement& interf = gm.one_interface( i ) ;
                for( index_t s = 0; s < interf.nb_children(); s++ ) {
                    index_t s_id = interf.child_id( s ).index ;
                    for( index_t el = 0; el < mesh.facets.nb_triangle( s_id );
                        el++ ) {
                        index_t tri = mesh.facets.triangle( s_id, el ) ;
                        for( index_t p = 0; p < mesh.facets.nb_vertices( tri );
                            p++ ) {
                            index_t vertex_id = mesh.facets.vertex( tri, p ) ;
                            data << " " << std::setw( 7 ) << vertex_id ;
                            new_line( count, 10, data ) ;
                        }
                    }
                    for( index_t el = 0; el < mesh.facets.nb_quad( s_id ); el++ ) {
                        index_t quad = mesh.facets.quad( s_id, el ) ;
                        for( index_t p = 0; p < mesh.facets.nb_vertices( quad );
                            p++ ) {
                            index_t vertex_id = mesh.facets.vertex( quad, p ) ;
                            data << " " << std::setw( 7 ) << vertex_id ;
                            new_line( count, 10, data ) ;
                        }
                    }
                }
            }
            for( index_t w = 0; w < mesh.edges.nb_wells(); w++ ) {
                for( index_t e = 0; e < mesh.edges.nb_edges( w ); e++ ) {
                    for( index_t v = 0; v < 2; v++ ) {
                        index_t vertex_id = mesh.edges.vertex( w, e, v ) ;
                        data << " " << std::setw( 7 ) << vertex_id ;
                        new_line( count, 10, data ) ;
                    }
                }
            }
            reset_line( count, data ) ;

            index_t nb_facets = 3 * mesh.facets.nb_triangle()
                + 4 * mesh.facets.nb_quad() + 4 * mesh.cells.nb_tet()
                + 5 * mesh.cells.nb_pyramid() + 5 * mesh.cells.nb_prism()
                + 6 * mesh.cells.nb_hex() + 2 * mesh.edges.nb_edges() ;
            data << nb_facets << " # PFVERTS" << std::endl ;
            for( index_t r = 0; r < gm.nb_regions(); r++ ) {
                for( index_t type = GEO::MESH_TET; type < GEO::MESH_CONNECTOR;
                    type++ ) {
                    GEO::MeshCellType T = static_cast< GEO::MeshCellType >( type ) ;
                    RINGMesh2CSMP& descriptor = *cell_type_to_cell_descriptor[T] ;
                    for( index_t el = 0; el < mesh.cells.nb_cells( r, T ); el++ ) {
                        index_t cell = mesh.cells.cell( r, el ) ;
                        for( index_t f = 0; f < descriptor.nb_facets; f++ ) {
                            index_t csmp_f = descriptor.facet[f] ;
                            index_t adj = mesh.cells.adjacent( cell, csmp_f ) ;
                            if( adj == GEO::NO_CELL ) {
                                data << " " << std::setw( 7 ) << -28 ;
                            } else {
                                data << " " << std::setw( 7 ) << adj ;
                            }
                            new_line( count, 10, data ) ;
                        }
                    }
                }
            }
            for( index_t i = 0; i < gm.nb_interfaces(); i++ ) {
                const GeoModelElement& interf = gm.one_interface( i ) ;
                for( index_t s = 0; s < interf.nb_children(); s++ ) {
                    index_t s_id = interf.child_id( s ).index ;
                    for( index_t el = 0; el < mesh.facets.nb_triangle( s_id );
                        el++ ) {
                        index_t tri = mesh.facets.triangle( s_id, el ) ;
                        for( index_t e = 0; e < mesh.facets.nb_vertices( tri );
                            e++ ) {
                            index_t adj = mesh.facets.adjacent( tri, e ) ;
                            if( adj == GEO::NO_FACET ) {
                                data << " " << std::setw( 7 ) << -28 ;
                            } else {
                                data << " " << std::setw( 7 ) << adj ;
                            }
                            new_line( count, 10, data ) ;
                        }
                    }
                    for( index_t el = 0; el < mesh.facets.nb_quad( s_id ); el++ ) {
                        index_t quad = mesh.facets.quad( s_id, el ) ;
                        for( index_t e = 0; e < mesh.facets.nb_vertices( quad );
                            e++ ) {
                            index_t adj = mesh.facets.adjacent( quad, e ) ;
                            if( adj == GEO::NO_FACET ) {
                                data << " " << std::setw( 7 ) << -28 ;
                            } else {
                                data << " " << std::setw( 7 ) << adj ;
                            }
                            new_line( count, 10, data ) ;
                        }
                    }
                }
            }
            index_t edge_offset = mesh.facets.nb() + mesh.cells.nb() ;
            index_t cur_edge = 0 ;
            for( index_t w = 0; w < mesh.edges.nb_wells(); w++ ) {
                data << " " << std::setw( 7 ) << -28 ;
                new_line( count, 10, data ) ;
                if( mesh.edges.nb_edges( w ) > 1 ) {
                    data << " " << std::setw( 7 ) << edge_offset + cur_edge + 1 ;
                    cur_edge++ ;
                    new_line( count, 10, data ) ;
                    for( index_t e = 1; e < mesh.edges.nb_edges( w ) - 1;
                        e++, cur_edge++ ) {
                        data << " " << std::setw( 7 ) << edge_offset + cur_edge - 1 ;
                        new_line( count, 10, data ) ;
                        data << " " << std::setw( 7 ) << edge_offset + cur_edge + 1 ;
                        new_line( count, 10, data ) ;
                    }
                    data << " " << std::setw( 7 ) << edge_offset + cur_edge - 1 ;
                    new_line( count, 10, data ) ;
                }
                data << " " << std::setw( 7 ) << -28 ;
                cur_edge++ ;
                new_line( count, 10, data ) ;
            }
            reset_line( count, data ) ;

            data << nb_total_elements << " # PMATERIAL" << std::endl ;
            for( index_t i = 0; i < nb_total_elements; i++ ) {
                data << " " << std::setw( 3 ) << 0 ;
                new_line( count, 20, data ) ;
            }
        }

    private:
        void new_line(
            index_t& count,
            index_t number_of_counts,
            std::ofstream& out ) const
        {
            count++ ;
            if( count == number_of_counts ) {
                count = 0 ;
                out << std::endl ;
            }
        }
        void reset_line( index_t& count, std::ofstream& out ) const
        {
            if( count != 0 ) {
                count = 0 ;
                out << std::endl ;
            }
        }
        void clear()
        {
            point_boundaries_.clear() ;
            box_model_ = false ;
            back_ = NO_ID ;
            top_ = NO_ID ;
            front_ = NO_ID ;
            bottom_ = NO_ID ;
            left_ = NO_ID ;
            right_ = NO_ID ;
            corner_boundary_flags_.clear() ;
            edge_boundary_flags_.clear() ;
            surface_boundary_flags_.clear() ;
        }
        void initialize( const GeoModel& gm )
        {
            clear() ;

            const GeoModel& model = gm ;
            std::string cmsp_filename = GEO::CmdLine::get_arg( "out:csmp" ) ;
            box_model_ = cmsp_filename != "" ;
            if( box_model_ ) {
                GEO::LineInput parser( cmsp_filename ) ;
                if( !parser.OK() ) {
                    throw RINGMeshException( "I/O",
                        "Cannot open file: " + cmsp_filename ) ;
                }
                parser.get_line() ;
                parser.get_fields() ;
                while( !parser.eof() ) {
                    if( parser.nb_fields() == 0 ) continue ;
                    if( parser.nb_fields() != 3 ) return ;
                    std::string type = parser.field( 1 ) ;
                    index_t interface_id = NO_ID ;
                    if( type == "NAME" ) {
                        std::string name = parser.field( 2 ) ;
                        for( index_t i = 0; i < model.nb_interfaces(); i++ ) {
                            if( model.one_interface( i ).name() == name ) {
                                interface_id = i ;
                                break ;
                            }
                        }
                    } else if( type == "ID" ) {
                        interface_id = parser.field_as_uint( 2 ) ;
                    } else {
                        throw RINGMeshException( "I/O", "Unknown type: " + type ) ;
                    }

                    std::string keyword = parser.field( 0 ) ;
                    if( keyword == "BACK" ) {
                        back_ = interface_id ;
                    } else if( keyword == "TOP" ) {
                        top_ = interface_id ;
                    } else if( keyword == "FRONT" ) {
                        front_ = interface_id ;
                    } else if( keyword == "BOTTOM" ) {
                        bottom_ = interface_id ;
                    } else if( keyword == "LEFT" ) {
                        left_ = interface_id ;
                    } else if( keyword == "RIGHT" ) {
                        right_ = interface_id ;
                    } else {
                        throw RINGMeshException( "I/O",
                            "Unknown keyword: " + keyword ) ;
                    }
                    parser.get_line() ;
                    parser.get_fields() ;
                }

                if( back_ == NO_ID || top_ == NO_ID || front_ == NO_ID
                    || bottom_ == NO_ID || left_ == NO_ID || right_ == NO_ID ) {
                    throw RINGMeshException( "I/O",
                        "Missing box shape information" ) ;
                }

                surface_boundary_flags_[back_] = -7 ;
                surface_boundary_flags_[top_] = -5 ;
                surface_boundary_flags_[front_] = -6 ;
                surface_boundary_flags_[bottom_] = -4 ;
                surface_boundary_flags_[left_] = -2 ;
                surface_boundary_flags_[right_] = -3 ;

                std::set< index_t > back_bottom ;
                back_bottom.insert( back_ ) ;
                back_bottom.insert( bottom_ ) ;
                edge_boundary_flags_[back_bottom] = -16 ;
                std::set< index_t > back_right ;
                back_right.insert( back_ ) ;
                back_right.insert( right_ ) ;
                edge_boundary_flags_[back_right] = -17 ;
                std::set< index_t > back_top ;
                back_top.insert( back_ ) ;
                back_top.insert( top_ ) ;
                edge_boundary_flags_[back_top] = -18 ;
                std::set< index_t > back_left ;
                back_left.insert( back_ ) ;
                back_left.insert( left_ ) ;
                edge_boundary_flags_[back_left] = -19 ;
                std::set< index_t > right_bottom ;
                right_bottom.insert( right_ ) ;
                right_bottom.insert( bottom_ ) ;
                edge_boundary_flags_[right_bottom] = -20 ;
                std::set< index_t > right_top ;
                right_top.insert( right_ ) ;
                right_top.insert( top_ ) ;
                edge_boundary_flags_[right_top] = -21 ;
                std::set< index_t > left_top ;
                left_top.insert( left_ ) ;
                left_top.insert( top_ ) ;
                edge_boundary_flags_[left_top] = -22 ;
                std::set< index_t > left_bottom ;
                left_bottom.insert( left_ ) ;
                left_bottom.insert( bottom_ ) ;
                edge_boundary_flags_[left_bottom] = -23 ;
                std::set< index_t > front_bottom ;
                front_bottom.insert( front_ ) ;
                front_bottom.insert( bottom_ ) ;
                edge_boundary_flags_[front_bottom] = -24 ;
                std::set< index_t > front_right ;
                front_right.insert( front_ ) ;
                front_right.insert( right_ ) ;
                edge_boundary_flags_[front_right] = -25 ;
                std::set< index_t > front_top ;
                front_top.insert( front_ ) ;
                front_top.insert( top_ ) ;
                edge_boundary_flags_[front_top] = -26 ;
                std::set< index_t > front_left ;
                front_left.insert( front_ ) ;
                front_left.insert( left_ ) ;
                edge_boundary_flags_[front_left] = -27 ;

                std::set< index_t > back_top_left ;
                back_top_left.insert( back_ ) ;
                back_top_left.insert( top_ ) ;
                back_top_left.insert( left_ ) ;
                corner_boundary_flags_[back_top_left] = -13 ;
                std::set< index_t > back_top_right ;
                back_top_right.insert( back_ ) ;
                back_top_right.insert( top_ ) ;
                back_top_right.insert( right_ ) ;
                corner_boundary_flags_[back_top_right] = -14 ;
                std::set< index_t > back_bottom_left ;
                back_bottom_left.insert( back_ ) ;
                back_bottom_left.insert( bottom_ ) ;
                back_bottom_left.insert( left_ ) ;
                corner_boundary_flags_[back_bottom_left] = -8 ;
                std::set< index_t > back_bottom_right ;
                back_bottom_right.insert( back_ ) ;
                back_bottom_right.insert( bottom_ ) ;
                back_bottom_right.insert( right_ ) ;
                corner_boundary_flags_[back_bottom_right] = -10 ;
                std::set< index_t > front_top_left ;
                front_top_left.insert( front_ ) ;
                front_top_left.insert( top_ ) ;
                front_top_left.insert( left_ ) ;
                corner_boundary_flags_[front_top_left] = -15 ;
                std::set< index_t > front_top_right ;
                front_top_right.insert( front_ ) ;
                front_top_right.insert( top_ ) ;
                front_top_right.insert( right_ ) ;
                corner_boundary_flags_[front_top_right] = -9 ;
                std::set< index_t > front_bottom_left ;
                front_bottom_left.insert( front_ ) ;
                front_bottom_left.insert( bottom_ ) ;
                front_bottom_left.insert( left_ ) ;
                corner_boundary_flags_[front_bottom_left] = -12 ;
                std::set< index_t > front_bottom_right ;
                front_bottom_right.insert( front_ ) ;
                front_bottom_right.insert( bottom_ ) ;
                front_bottom_right.insert( right_ ) ;
                corner_boundary_flags_[front_bottom_right] = -11 ;
            }

            point_boundaries_.resize( gm.mesh.vertices.nb() ) ;
            for( index_t s = 0; s < model.nb_surfaces(); s++ ) {
                index_t interface_id = model.surface( s ).parent_id().index ;
                for( index_t f = 0; f < gm.mesh.facets.nb_facets( s ); f++ ) {
                    index_t f_id = gm.mesh.facets.facet( s, f ) ;
                    for( index_t v = 0; v < gm.mesh.facets.nb_vertices( f_id );
                        v++ ) {
                        index_t vertex_id = gm.mesh.facets.vertex( f_id, v ) ;
                        point_boundaries_[vertex_id].insert( interface_id ) ;
                    }
                }
            }
        }
        std::string interface_name( index_t i, const GeoModel& gm )
        {
            if( box_model_ ) {
                if( i == back_ ) {
                    return "BACK" ;
                } else if( i == top_ ) {
                    return "TOP" ;
                } else if( i == front_ ) {
                    return "FRONT" ;
                } else if( i == bottom_ ) {
                    return "BOTTOM" ;
                } else if( i == left_ ) {
                    return "LEFT" ;
                } else if( i == right_ ) {
                    return "RIGHT" ;
                }
            }
            return gm.one_interface( i ).name() ;
        }
        signed_index_t point_boundary( index_t p )
        {
            ringmesh_assert( p < point_boundaries_.size() ) ;
            const std::set< unsigned int >& boundaries = point_boundaries_[p] ;
            if( box_model_ ) {
                if( boundaries.size() == 1 ) {
                    std::map< unsigned int, int >::iterator it =
                        surface_boundary_flags_.find( *boundaries.begin() ) ;
                    ringmesh_assert( it != surface_boundary_flags_.end() ) ;
                    return it->second ;
                } else if( boundaries.size() == 2 ) {
                    std::map< std::set< unsigned int >, int >::iterator it =
                        edge_boundary_flags_.find( boundaries ) ;
                    ringmesh_assert( it != edge_boundary_flags_.end() ) ;
                    return it->second ;
                } else if( boundaries.size() == 3 ) {
                    std::map< std::set< unsigned int >, int >::iterator it =
                        corner_boundary_flags_.find( boundaries ) ;
                    ringmesh_assert( it != corner_boundary_flags_.end() ) ;
                    return it->second ;
                } else {
                    return 0 ;
                }
            } else {
                if( boundaries.empty() )
                    return 0 ;
                else
                    return -28 ;
            }
        }
    private:
        std::vector< std::set< index_t > > point_boundaries_ ;

        bool box_model_ ;
        index_t back_ ;
        index_t top_ ;
        index_t front_ ;
        index_t bottom_ ;
        index_t left_ ;
        index_t right_ ;

        std::map< std::set< index_t >, signed_index_t > corner_boundary_flags_ ;
        std::map< std::set< index_t >, signed_index_t > edge_boundary_flags_ ;
        std::map< index_t, signed_index_t > surface_boundary_flags_ ;
    } ;

    /************************************************************************/

    class GPRSIOHandler: public GeoModelIOHandler {
    public:
        struct Pipe {
            Pipe( index_t v0_in, index_t v1_in )
                : v0( v0_in ), v1( v1_in )
            {
            }
            index_t v0 ;
            index_t v1 ;
        } ;
        virtual void load( const std::string& filename, GeoModel& geomodel )
        {
            throw RINGMeshException( "I/O",
                "Loading of a GeoModel from GPRS not implemented yet" ) ;
        }
        virtual void save( const GeoModel& gm, const std::string& filename )
        {
            std::string path = GEO::FileSystem::dir_name( filename ) ;
            std::string directory = GEO::FileSystem::base_name( filename ) ;
            if( path == "." ) {
                path = GEO::FileSystem::get_current_working_directory() ;
            }
            std::ostringstream oss ;
            oss << path << "/" << directory ;
            std::string full_path = oss.str() ;
            GEO::FileSystem::create_directory( full_path ) ;

            std::ostringstream oss_pipes ;
            oss_pipes << full_path << "/pipes.in" ;
            std::ofstream out_pipes( oss_pipes.str().c_str() ) ;

            std::ostringstream oss_vol ;
            oss_vol << full_path << "/vol.in" ;
            std::ofstream out_vol( oss_vol.str().c_str() ) ;
            out_vol.precision( 16 ) ;

            std::ostringstream oss_xyz ;
            oss_xyz << full_path << "/gprs.xyz" ;
            std::ofstream out_xyz( oss_xyz.str().c_str() ) ;
            out_xyz.precision( 16 ) ;

            const GeoModelMesh& mesh = gm.mesh ;
            std::deque< Pipe > pipes ;
            index_t cell_offset = mesh.cells.nb() ;
            for( index_t c = 0; c < mesh.cells.nb(); c++ ) {
                for( index_t f = 0; f < mesh.cells.nb_facets( c ); f++ ) {
                    index_t facet = NO_ID ;
                    bool not_used ;
                    if( mesh.cells.is_cell_facet_on_surface( c, f, facet,
                        not_used ) ) {
                        pipes.push_back( Pipe( c, facet + cell_offset ) ) ;
                    } else {
                        index_t adj = mesh.cells.adjacent( c, f ) ;
                        if( adj != GEO::NO_CELL && adj < c ) {
                            pipes.push_back( Pipe( c, adj ) ) ;
                        }
                    }
                }
            }

            index_t nb_edges = 0 ;
            for( index_t l = 0; l < gm.nb_lines(); l++ ) {
                nb_edges += gm.line( l ).nb_cells() ;
            }
            std::vector< index_t > temp ;
            temp.reserve( 3 ) ;
            std::vector< std::vector< index_t > > edges( nb_edges, temp ) ;
            std::vector< vec3 > edge_vertices( nb_edges ) ;
            index_t count_edge = 0 ;
            for( index_t l = 0; l < gm.nb_lines(); l++ ) {
                const Line& line = gm.line( l ) ;
                for( index_t e = 0; e < line.nb_cells(); e++ ) {
                    edge_vertices[count_edge++ ] = 0.5
                        * ( line.vertex( e ) + line.vertex( e + 1 ) ) ;
                }
            }
            ColocaterANN ann( edge_vertices, false ) ;

            for( index_t f = 0; f < mesh.facets.nb(); f++ ) {
                for( index_t e = 0; e < mesh.facets.nb_vertices( f ); e++ ) {
                    index_t adj = mesh.facets.adjacent( f, e ) ;
                    if( adj != GEO::NO_CELL && adj < f ) {
                        pipes.push_back(
                            Pipe( f + cell_offset, adj + cell_offset ) ) ;
                    } else {
                        const vec3& e0 = mesh.vertices.vertex(
                            mesh.facets.vertex( f, e ) ) ;
                        const vec3& e1 = mesh.vertices.vertex(
                            mesh.facets.vertex( f,
                                ( e + 1 ) % mesh.facets.nb_vertices( f ) ) ) ;
                        vec3 query = 0.5 * ( e0 + e1 ) ;
                        std::vector< index_t > results ;
                        if( ann.get_colocated( query, results ) ) {
                            edges[results[0]].push_back( cell_offset + f ) ;
                        } else {
                            ringmesh_assert_not_reached ;
                        }
                    }
                }
            }

            index_t nb_pipes = pipes.size() ;
            for( index_t e = 0; e < edges.size(); e++ ) {
                nb_pipes += binomial_coef( edges[e].size() ) ;
            }
            out_pipes << nb_pipes << std::endl ;
            for( index_t p = 0; p < pipes.size(); p++ ) {
                const Pipe& pipe = pipes[p] ;
                out_pipes << pipe.v0 << SPACE << pipe.v1 << std::endl ;
            }
            for( index_t e = 0; e < edges.size(); e++ ) {
                const std::vector< index_t > vertices = edges[e] ;
                for( index_t v0 = 0; v0 < vertices.size() - 1; v0++ ) {
                    for( index_t v1 = v0 + 1; v1 < vertices.size(); v1++ ) {
                        out_pipes << vertices[v0] << SPACE << vertices[v1]
                            << std::endl ;
                    }
                }
            }

            out_xyz
                << "Node geometry, not used by GPRS but useful to reconstruct a pipe-network"
                << std::endl ;
            for( index_t c = 0; c < mesh.cells.nb(); c++ ) {
                out_xyz << mesh.cells.center( c ) << std::endl ;
                out_vol << mesh.cells.volume( c ) << std::endl ;
            }
            for( index_t f = 0; f < mesh.facets.nb(); f++ ) {
                out_xyz << mesh.facets.center( f ) << std::endl ;
                out_vol << mesh.facets.area( f ) << std::endl ;
            }
        }
        index_t binomial_coef( index_t n ) const
        {
            switch( n ) {
                case 1:
                    return 0 ;
                case 2:
                    return 1 ;
                case 3:
                    return 3 ;
                case 4:
                    return 6 ;
                case 5:
                    return 10 ;
                case 6:
                    return 15 ;
                case 7:
                    return 21 ;
                case 8:
                    return 28 ;
                case 9:
                    return 36 ;
                case 10:
                    return 45 ;
                default:
                    ringmesh_assert_not_reached ;
                    return 0 ;

            }
        }
    } ;

    /************************************************************************/

//        struct RINGMesh2GMSH {
//                   index_t element_type ;
//                   index_t nb_vertices ;
//                   index_t vertices[8] ;
//                   index_t nb_facets ;
//                   index_t nb_vertices_in_facet[6] ;
//                   index_t facet[6] ;
//                   index_t vertices_in_facet[6][4] ;
//               } ;
//
//               static RINGMesh2GMSH tet_descriptor_gmsh = { 4,                  // type
//                   4,                  // nb vertices
//                   { 0, 1, 2, 3, 5, 8 ,9 ,4,6,7 },     // vertices
//                   4,                  // nb facets
//                   { 3, 3, 3, 3 },     // nb vertices in facet
//                   { 0, 1, 2, 3 },     // facets
//                   { { 1, 3, 2 }, { 0, 2, 3 }, { 3, 1, 0 }, { 0, 1, 2 } } } ;
//
//               static RINGMesh2GMSH hex_descriptor_gmsh = { 6,                         // type
//                   8,                              // nb vertices
//                   { 4, 0, 5, 1, 7, 3, 6, 2 },     // vertices
//                   6,                              // nb facets
//                   { 4, 4, 4, 4, 4, 4 },           // nb vertices in facet
//                   { 4, 2, 1, 3, 0, 5 },           // facets
//                   { { 0, 3, 7, 4 }, { 2, 1, 5, 6 }, { 1, 0, 4, 5 }, { 3, 2, 6, 7 }, {
//                       1, 2, 3, 0 }, { 4, 7, 6, 5 } } } ;
//
//               static RINGMesh2GMSH prism_descriptor_gmsh = { 12,                     // type
//                   6,                      // nb vertices
//                   { 0, 1, 2, 3, 4, 5 },   // vertices
//                   5,                      // nb facets
//                   { 3, 4, 4, 4, 3 },      // nb vertices in facet
//                   { 0, 2, 4, 3, 1 },      // facets
//                   {
//                       { 0, 1, 2 }, { 3, 5, 4 }, { 0, 3, 4, 1 }, { 0, 2, 5, 3 }, {
//                           1, 4, 5, 2 } } } ;
//
//               static RINGMesh2GMSH pyramid_descriptor_gmsh = { 18,                 // type
//                   5,                  // nb vertices
//                   { 0, 1, 2, 3, 4 },  // vertices
//                   5,                  // nb facets
//                   { 3, 3, 3, 3, 4 },  // nb vertices in facet
//                   { 1, 3, 4, 2, 0 },  // facets
//                   { { 0, 1, 2, 3 }, { 0, 4, 1 }, { 0, 3, 4 }, { 2, 4, 3 }, { 2, 1, 4 } } } ;
    class MSHIOHandler: public GeoModelIOHandler {
    public:
        virtual void load( const std::string& filename, GeoModel& geomodel )
        {
            throw RINGMeshException( "I/O",
                "Loading of a GeoModel from GMSH not implemented yet" ) ;
        }
        virtual void save( const GeoModel& gm, const std::string& filename )
        {
            /// @todo after implementing GMMOrder
            throw RINGMeshException( "I/O",
                "Saving of a GeoModel from GMSH not implemented yet" ) ;
//                gm.set_duplicate_mode( FAULT ) ;

            std::ofstream out( filename.c_str() ) ;
            out.precision( 16 ) ;

            out << "$MeshFormat" << std::endl ;
            out << "2.2 0 8" << std::endl ;
            out << "$EndMeshFormat" << std::endl ;

            out << "$Nodes" << std::endl ;
//            out << gm.order.nb_total_vertices() << std::endl ;
//            for( index_t p = 0; p < gm.vertices.nb(); p++ ) {
//
//                const vec3& point = gm.vertices.vertex( p ) ;
//                if( p == 0 ) {
//                    std::cout << "io val " << point.x << std::endl ;
//
//                }
//                out << p + 1 << SPACE << point.x << SPACE << point.y << SPACE
//                    << point.z << std::endl ;
//            }
//            index_t vertex_offset = gm.vertices.nb() ;
//            for( index_t p = 0; p < gm.vertices.nb_duplicated_vertices(); p++ ) {
//                const vec3& point = gm.vertices.duplicated_vertex( p ) ;
//                out << vertex_offset + p + 1 << SPACE << point.x << SPACE << point.y
//                    << SPACE << point.z << std::endl ;
//            }
//            vertex_offset += gm.vertices.nb_duplicated_vertices() ;
//            index_t nb_order_vertices = gm.order.nb() ;
//            for( index_t p = 0; p < nb_order_vertices; p++ ) {
//                out << vertex_offset + p + 1 << SPACE << gm.order.point( p ).x
//                    << SPACE << gm.order.point( p ).y << SPACE
//                    << gm.order.point( p ).z << std::endl ;
//            }
//            out << "$EndNodes" << std::endl ;
//
//            index_t cell_type[4] = { 4, 5, 6, 7 } ;
//            index_t facet_type[5] = { -1, -1, -1, 2, 3 } ;
//            if( gm.get_order() == 2 ) {
//                cell_type[0] = 11 ;
//                cell_type[1] = 17 ;
//                cell_type[2] = 18 ;
//                cell_type[3] = 19 ;
//                facet_type[0] = -1 ;
//                facet_type[1] = -1 ;
//                facet_type[2] = -1 ;
//                facet_type[3] = 9 ;
//                facet_type[4] = 16 ;
//            } else if( gm.get_order() > 2 ) {
//                GEO::Logger::err( "" ) << "The order " << gm.get_order() << " "
//                    << "is not supported"
//                    << " for the gmsh export. The export will take order 1 elements"
//                    << std::endl ;
//            }
//            const GeoModel& model = gm ;
//            index_t offset_region = gm.nb_regions() ;
//            index_t offset_interface = model.nb_interfaces() * 2 ; // one for each side
//            index_t nb_facets = 0 ;
//            std::vector< ColocaterANN* > anns( model.nb_surfaces(), nil ) ;
//            for( index_t s = 0; s < model.nb_surfaces(); s++ ) {
//                if( gm.vertices.is_surface_to_duplicate( s ) )
//                    nb_facets += 2 * gm.facets.nb_facets( s ) ;
//                else
//                    nb_facets += gm.facets.nb_facets( s ) ;
//                anns[s] = new ColocaterANN( model.surface( s ).mesh(),
//                    ColocaterANN::FACETS ) ;
//            }
//            out << "$Elements" << std::endl ;
//            out << gm.cells.nb_cells() + nb_facets << std::endl ;
//            index_t cur_cell = 1 ;
//            for( index_t m = 0; m < gm.nb_regions(); m++ ) {
//                const GEO::Mesh& mesh = gm.mesh( m ) ;
//                GEO::Attribute< index_t > attribute( mesh.facets.attributes(),
//                    surface_att_name ) ;
//                const GeoModelElement& region = model.region( m ) ;
//                std::vector< index_t > surfaces ;
//                surfaces.reserve( region.nb_boundaries() ) ;
//                for( index_t b = 0; b < region.nb_boundaries(); b++ ) {
//                    index_t cur_s_id = region.boundary_gme( b ).index ;
//                    if( !gm.vertices.is_surface_to_duplicate( cur_s_id ) ) continue ;
//                    surfaces.push_back( cur_s_id ) ;
//                }
//                for( index_t c = 0; c < mesh.cells.nb(); c++ ) {
//                    out << cur_cell++ << " " << cell_type[mesh.cells.type( c )]
//                        << " 2 " << m + 1 << SPACE << m ;
//                    for( index_t v = mesh.cells.corners_begin( c );
//                        v < mesh.cells.corners_end( c ); v++ ) {
//                        index_t vertex_id ;
//                        index_t duplicated_vertex_id ;
//                        out << SPACE ;
//                        if( gm.vertices.vertex_id( m, v, vertex_id,
//                            duplicated_vertex_id ) ) {
//                            out << vertex_id + 1 ;
//                        } else {
//                            out << vertex_offset + duplicated_vertex_id + 1 ;
//                        }
//                    }
//                    if( gm.get_order() == 2 ) {
//                        out << SPACE ;
//                        out << gm.order.get_id_on_cell( m, c, 3 ) + 1 ;
//                        out << SPACE ;
//                        out << gm.order.get_id_on_cell( m, c, 0 ) + 1 ;
//                        out << SPACE ;
//                        out << gm.order.get_id_on_cell( m, c, 4 ) + 1 ;
//                        out << SPACE ;
//                        out << gm.order.get_id_on_cell( m, c, 5 ) + 1 ;
//                        out << SPACE ;
//                        out << gm.order.get_id_on_cell( m, c, 1 ) + 1 ;
//                        out << SPACE ;
//                        out << gm.order.get_id_on_cell( m, c, 2 ) + 1 ;
//                    }
//                    out << std::endl ;
//
//                    for( index_t f = 0; f < mesh.cells.nb_facets( c ); f++ ) {
//                        vec3 facet_bary = mesh_cell_facet_center( mesh, c,
//                            f ) ;
//                        vec3 cell_facet_normal = mesh_cell_facet_normal( mesh,
//                            c, f ) ;
//                        for( index_t s = 0; s < surfaces.size(); s++ ) {
//                            index_t surface_id = surfaces[s] ;
//                            std::vector< index_t > result ;
//                            if( anns[surface_id]->get_colocated( facet_bary,
//                                result ) ) {
//                                vec3 facet_normal =
//                                    model.surface( surface_id ).facet_normal(
//                                        result[0] ) ;
//                                bool side = dot( facet_normal, cell_facet_normal )
//                                    > 0 ;
//                                out << cur_cell++ << " "
//                                    << facet_type[mesh.cells.facet_nb_vertices( c,
//                                        f )] << " 2 "
//                                    << offset_region
//                                        + 2
//                                            * model.surface( surface_id ).parent_id().index
//                                        + side + 1 << SPACE
//                                    << offset_region + offset_interface
//                                        + 2 * surface_id + side ;
//                                for( index_t v = 0;
//                                    v < mesh.cells.facet_nb_vertices( c, f ); v++ ) {
//                                    index_t corner_id =
//                                        mesh.cells.corner( c,
//                                            mesh.cells.descriptor( c ).facet_vertex[f][v] ) ;
//                                    index_t vertex_id ;
//                                    index_t duplicated_vertex_id ;
//                                    out << SPACE ;
//                                    if( gm.vertices.vertex_id( m, corner_id,
//                                        vertex_id, duplicated_vertex_id ) ) {
//                                        out << vertex_id + 1 ;
//                                    } else {
//                                        out
//                                            << vertex_offset + duplicated_vertex_id
//                                                + 1 ;
//                                    }
//                                }
//                                out << std::endl ;
//                                break ;
//                            }
//                        }
//                    }
//                }
//            }
//
//            for( index_t i = 0; i < model.nb_interfaces(); i++ ) {
//                const GeoModelElement& interf = model.one_interface( i ) ;
//                for( index_t s = 0; s < interf.nb_children(); s++ ) {
//                    index_t s_id = interf.child_id( s ).index ;
//                    if( gm.vertices.is_surface_to_duplicate( s_id ) ) continue ;
//                    index_t mesh_id = gm.facets.mesh( s_id ) ;
//                    const GEO::Mesh& mesh = gm.mesh( mesh_id ) ;
//                    for( index_t t = 0; t < gm.facets.nb_facets( s_id ); t++ ) {
//                        index_t facet_id = gm.facets.facet( s_id, t ) ;
//                        out << cur_cell++ << SPACE
//                            << facet_type[mesh.facets.nb_vertices( facet_id )]
//                            << " 2 " << offset_region + 2 * i + 1 << SPACE
//                            << offset_region + offset_interface + 2 * s_id ;
//                        for( index_t v = 0; v < mesh.facets.nb_vertices( facet_id );
//                            v++ ) {
//                            index_t v_id = mesh.facets.vertex( facet_id, v ) ;
//                            out << SPACE
//                                << gm.vertices.vertex_id( mesh_id, v_id ) + 1 ;
//                        }
//                        for( index_t v = 0;
//                            v
//                                < mesh.facets.nb_vertices( facet_id )
//                                    * ( gm.get_order() - 1 ); v++ ) {
//                            out << SPACE ;
//                            out << gm.order.get_id_on_facet( s, facet_id, v ) + 1 ;
//                        }
//                        out << std::endl ;
//                    }
//                }
//            }
//            out << "$EndElements" << std::endl ;
//
//            if( GEO::CmdLine::get_arg_bool( "out:kine3d" ) ) {
//                std::string directory = GEO::FileSystem::dir_name( filename ) ;
//                std::string file = GEO::FileSystem::base_name( filename ) ;
//                std::ostringstream oss_kine ;
//                oss_kine << directory << "/" << file << ".gmsh_info" ;
//                std::ofstream kine3d( oss_kine.str().c_str() ) ;
//                for( index_t i = 0; i < model.nb_interfaces(); i++ ) {
//                    const GeoModelElement& interf = model.one_interface( i ) ;
//                    index_t s_id = interf.child_id( 0 ).index ;
//                    kine3d << offset_region + 2 * i + 1 << ":" << interf.name()
//                        << ",1," ;
//                    const RINGMesh::GeoModelElement& E = model.one_interface( i ) ;
//                    if( RINGMesh::GeoModelElement::is_fault(
//                        E.geological_feature() ) ) {
//                        kine3d << "FaultFeatureClass" ;
//                    } else if( RINGMesh::GeoModelElement::is_stratigraphic_limit(
//                        E.geological_feature() ) ) {
//                        kine3d << "HorizonFeatureClass" ;
//                    } else if( E.is_on_voi() ) {
//                        kine3d << "ModelRINGMesh::BoundaryFeatureClass" ;
//                    }
//                    kine3d << std::endl ;
//                    if( gm.vertices.is_surface_to_duplicate( s_id ) ) {
//                        kine3d << offset_region + 2 * i + 1 << ":" << interf.name()
//                            << ",0," ;
//                        const RINGMesh::GeoModelElement& E = model.one_interface(
//                            i ) ;
//                        if( RINGMesh::GeoModelElement::is_fault(
//                            E.geological_feature() ) ) {
//                            kine3d << "FaultFeatureClass" ;
//                        } else if( RINGMesh::GeoModelElement::is_stratigraphic_limit(
//                            E.geological_feature() ) ) {
//                            kine3d << "HorizonFeatureClass" ;
//                        } else if( E.is_on_voi() ) {
//                            kine3d << "ModelRINGMesh::BoundaryFeatureClass" ;
//                        }
//                        kine3d << std::endl ;
//                    }
//                }
//            }
        }
    } ;
}

namespace RINGMesh {

    void geomodel_load( GeoModel& model, const std::string& filename )
    {
        GEO::Logger::out( "I/O" ) << "Loading file " << filename << "..."
            << std::endl ;

        IOHandler_var handler = GeoModelIOHandler::get_handler( filename ) ;
        handler->load( filename, model ) ;
    }

    void geomodel_save( const GeoModel& model, const std::string& filename )
    {
        GEO::Logger::out( "I/O" ) << "Saving file " << filename << "..."
            << std::endl ;

        IOHandler_var handler = GeoModelIOHandler::get_handler( filename ) ;
        handler->save( model, filename ) ;

    }

    /************************************************************************/

    /*
     * Initializes the possible handler for IO files
     */
    void GeoModelIOHandler::initialize_full_geomodel_output()
    {
        ringmesh_register_IOHandler_creator( LMIOHandler, "meshb" ) ;
        ringmesh_register_IOHandler_creator( LMIOHandler, "mesh" );
        ringmesh_register_IOHandler_creator( TetGenIOHandler, "tetgen" );
        ringmesh_register_IOHandler_creator( TSolidIOHandler, "so" );
        ringmesh_register_IOHandler_creator( CSMPIOHandler, "csmp" );
        ringmesh_register_IOHandler_creator( AsterIOHandler, "mail" );
        ringmesh_register_IOHandler_creator( VTKIOHandler, "vtk" );
        ringmesh_register_IOHandler_creator( GPRSIOHandler, "gprs" );
        ringmesh_register_IOHandler_creator( MSHIOHandler, "msh" );
        ringmesh_register_IOHandler_creator( MFEMIOHandler, "mfem" );
        ringmesh_register_IOHandler_creator( GeoModelHandler, "gm" );}

}
