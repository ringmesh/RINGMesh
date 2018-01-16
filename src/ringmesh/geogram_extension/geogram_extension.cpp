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

#include <ringmesh/geogram_extension/geogram_extension.h>

#include <algorithm>

#include <geogram/basic/file_system.h>
#include <geogram/basic/line_stream.h>

#include <geogram/mesh/mesh.h>
#include <geogram/mesh/mesh_geometry.h>
#include <geogram/mesh/mesh_io.h>
#include <geogram/mesh/mesh_repair.h>

#include <ringmesh/basic/logger.h>

/*!
 * @todo Re-orgarnize this mess
 */

namespace
{
    using namespace RINGMesh;

    /*!
     * @brief TSurfMeshIOHandler for importing .ts files into a mesh.
     */
    class TSurfMeshIOHandler final : public GEO::MeshIOHandler
    {
    public:
        TSurfMeshIOHandler()
            : starting_index_( 1 ),
              mesh_dimension_( 3 ),
              nb_vertices_( 0 ),
              nb_triangles_( 0 ),
              z_sign_( 1 )
        {
        }

        /*!
         * @brief Load a TSurf saved in .ts format
         * @warning Assumes there is only one TSurf in the file.
         *          Will undoubtedly crash if it is not the case.
         * @todo Prevent crashing if the file is not as expected
         *
         * @param filename the name of the .ts file to be processed.
         * @param mesh to which vertices and facets will be assigned
         * @param flags Not used for now.
         */
        bool load( const std::string& filename,
            GEO::Mesh& mesh,
            const GEO::MeshIOFlags& flag = GEO::MeshIOFlags() ) final
        {
            ringmesh_unused( flag );
            filename_ = filename;
            if( !is_file_valid() )
            {
                return false;
            }
            else
            {
                read_number_of_vertices_and_triangles();
                read_vertex_property_names();
                read_vertex_property_sizes();
                allocate_vertices();
                allocate_triangles();
                allocate_vertex_properties();
                read_vertices_and_triangles();
                assign_and_repair_mesh( mesh );
                return true;
            }
        }

        /*!
         * @brief Save a Mesh in .ts format
         */
        bool save( const GEO::Mesh& mesh,
            const std::string& filename,
            const GEO::MeshIOFlags& flag = GEO::MeshIOFlags() ) final
        {
            ringmesh_unused( flag );
            if( !mesh.facets.are_simplices() )
            {
                throw RINGMeshException( "I/O",
                    "Cannot save a non triangulated mesh into TSurf format" );
            }

            std::ofstream out( filename.c_str() );
            out.precision( 16 );
            save_header( out, GEO::FileSystem::base_name( filename ) );
            std::vector< std::string > att_v_double_names;
            std::vector< index_t > vertex_attr_dims;
            fill_vertex_attribute_header(
                mesh, out, att_v_double_names, vertex_attr_dims );

            out << "TFACE" << std::endl;
            save_vertices( out, mesh, att_v_double_names, vertex_attr_dims );
            save_triangles( out, mesh );
            out << "END" << std::endl;

            return true;
        }

    private:
        void save_vertices( std::ofstream& out,
            const GEO::Mesh& mesh,
            const std::vector< std::string >& att_v_double_names,
            const std::vector< index_t >& vertex_attr_dims )
        {
            for( auto v : range( mesh.vertices.nb() ) )
            {
                // PVRTX must be used instead of VRTX because
                // properties are not read by Gocad if it is VRTX.
                out << "PVRTX " << v + starting_index_ << " "
                    << mesh.vertices.point( v );
                for( auto attr_dbl_itr : range( att_v_double_names.size() ) )
                {
                    GEO::Attribute< double > cur_attr(
                        mesh.vertices.attributes(),
                        att_v_double_names[attr_dbl_itr] );
                    index_t nb_dimensions = vertex_attr_dims[attr_dbl_itr];
                    for( auto dim_itr : range( nb_dimensions ) )
                    {
                        out << " " << cur_attr[v * nb_dimensions + dim_itr];
                    }
                }
                out << std::endl;
            }
        }
        void save_triangles( std::ofstream& out, const GEO::Mesh& mesh )
        {
            for( auto t : range( mesh.facets.nb() ) )
            {
                out << "TRGL";
                for( auto v : range( 3 ) )
                {
                    out << " " << mesh.facets.vertex( t, v ) + starting_index_;
                }
                out << std::endl;
            }
        }
        void save_header( std::ofstream& out, const std::string& mesh_name )
        {
            out << "GOCAD TSurf" << std::endl;
            out << "HEADER {" << std::endl;
            out << "name: " << mesh_name << std::endl;
            out << "}" << std::endl;
        }
        void fill_vertex_attribute_header( const GEO::Mesh& mesh,
            std::ofstream& out,
            std::vector< std::string >& att_v_double_names,
            std::vector< index_t >& vertex_attr_dims ) const
        {
            GEO::vector< std::string > att_v_names;
            GEO::AttributesManager& mesh_vertex_mgr =
                mesh.vertices.attributes();
            mesh_vertex_mgr.list_attribute_names( att_v_names );
            for( auto att_v : range( mesh_vertex_mgr.nb() ) )
            {
                if( att_v_names[att_v] == "point" )
                {
                    continue;
                }

                if( !GEO::Attribute< double >::is_defined(
                        mesh_vertex_mgr, att_v_names[att_v] ) )
                {
                    continue;
                }
                att_v_double_names.push_back( att_v_names[att_v] );
                index_t cur_dim =
                    mesh_vertex_mgr.find_attribute_store( att_v_names[att_v] )
                        ->dimension();
                vertex_attr_dims.push_back( cur_dim );
            }

            if( !att_v_double_names.empty() )
            {
                index_t nb_attributes =
                    static_cast< index_t >( att_v_double_names.size() );
                out << "PROPERTIES";
                for( auto attr_dbl_itr : range( nb_attributes ) )
                {
                    out << " " << att_v_double_names[attr_dbl_itr];
                }
                out << std::endl;
                out << "PROP_LEGAL_RANGES";
                for( auto attr_dbl_itr : range( nb_attributes ) )
                {
                    ringmesh_unused( attr_dbl_itr );
                    out << " **none**  **none**";
                }
                out << std::endl;
                out << "NO_DATA_VALUES";
                for( auto attr_dbl_itr : range( nb_attributes ) )
                {
                    ringmesh_unused( attr_dbl_itr );
                    out << " -99999";
                }
                out << std::endl;
                out << "READ_ONLY";
                for( auto attr_dbl_itr : range( nb_attributes ) )
                {
                    ringmesh_unused( attr_dbl_itr );
                    out << " 1";
                }
                out << std::endl;
                out << "PROPERTY_CLASSES";
                for( auto attr_dbl_itr : range( nb_attributes ) )
                {
                    out << " " << att_v_double_names[attr_dbl_itr];
                }
                out << std::endl;
                out << "PROPERTY_KINDS";
                for( auto attr_dbl_itr : range( nb_attributes ) )
                {
                    ringmesh_unused( attr_dbl_itr );
                    out << " \"Real Number\"";
                }
                out << std::endl;
                out << "PROPERTY_SUBCLASSES";
                for( auto attr_dbl_itr : range( nb_attributes ) )
                {
                    ringmesh_unused( attr_dbl_itr );
                    out << " QUANTITY Float";
                }
                out << std::endl;
                out << "ESIZES";
                for( auto attr_dbl_itr : range( nb_attributes ) )
                {
                    out << " "
                        << std::to_string( vertex_attr_dims[attr_dbl_itr] );
                }
                out << std::endl;
                out << "UNITS";
                for( auto attr_dbl_itr : range( nb_attributes ) )
                {
                    ringmesh_unused( attr_dbl_itr );
                    out << " unitless";
                }
                out << std::endl;
                for( auto attr_dbl_itr : range( nb_attributes ) )
                {
                    out << "PROPERTY_CLASS_HEADER "
                        << att_v_double_names[attr_dbl_itr] << " {"
                        << std::endl;
                    out << "kind: Real Number" << std::endl;
                    out << "unit: unitless" << std::endl;
                    out << "}" << std::endl;
                }
            }
        }
        // This function read the z_sign too [PA]
        void read_number_of_vertices_and_triangles()
        {
            GEO::LineInput in( filename_ );
            while( !in.eof() && in.get_line() )
            {
                in.get_fields();
                if( in.nb_fields() > 0 )
                {
                    if( in.field_matches( 0, "ZPOSITIVE" ) )
                    {
                        if( in.field_matches( 1, "Elevation" ) )
                        {
                            z_sign_ = 1;
                        }
                        else if( in.field_matches( 1, "Depth" ) )
                        {
                            z_sign_ = -1;
                        }
                    }
                    else if( in.field_matches( 0, "VRTX" )
                             || in.field_matches( 0, "PVRTX" ) )
                    {
                        nb_vertices_++;
                    }
                    else if( in.field_matches( 0, "PATOM" )
                             || in.field_matches( 0, "ATOM" ) )
                    {
                        nb_vertices_++;
                    }
                    else if( in.field_matches( 0, "TRGL" ) )
                    {
                        nb_triangles_++;
                    }
                }
            }
        }

        void read_vertices_and_triangles()
        {
            GEO::LineInput in( filename_ );
            index_t v = 0;
            index_t t = 0;
            while( !in.eof() && in.get_line() )
            {
                in.get_fields();
                if( in.nb_fields() > 0 )
                {
                    if( in.field_matches( 0, "VRTX" )
                        || in.field_matches( 0, "PVRTX" ) )
                    {
                        vertices_[mesh_dimension_ * v] =
                            in.field_as_double( 2 );
                        vertices_[mesh_dimension_ * v + 1] =
                            in.field_as_double( 3 );
                        vertices_[mesh_dimension_ * v + 2] =
                            in.field_as_double( 4 ) * z_sign_;
                        if( in.field_matches( 0, "PVRTX" ) )
                        {
                            index_t offset = 5;
                            for( auto prop_name_itr :
                                range( vertex_property_names_.size() ) )
                            {
                                for( auto v_attr_dim_itr :
                                    range( vertex_attribute_dims_
                                            [prop_name_itr] ) )
                                {
                                    vertex_attributes_[prop_name_itr][v]
                                                      [v_attr_dim_itr] =
                                                          in.field_as_double(
                                                              offset );
                                    ++offset;
                                }
                            }
                        }
                        ++v;
                    }
                    else if( in.field_matches( 0, "PATOM" )
                             || in.field_matches( 0, "ATOM" ) )
                    {
                        index_t v0 = in.field_as_uint( 2 ) - starting_index_;
                        vertices_[mesh_dimension_ * v] =
                            vertices_[mesh_dimension_ * v0];
                        vertices_[mesh_dimension_ * v + 1] =
                            vertices_[mesh_dimension_ * v0 + 1];
                        vertices_[mesh_dimension_ * v + 2] =
                            vertices_[mesh_dimension_ * v0 + 2];
                        ++v;
                    }
                    else if( in.field_matches( 0, "TRGL" ) )
                    {
                        triangles_[3 * t] =
                            index_t( in.field_as_uint( 1 ) - starting_index_ );
                        triangles_[3 * t + 1] =
                            index_t( in.field_as_uint( 2 ) - starting_index_ );
                        triangles_[3 * t + 2] =
                            index_t( in.field_as_uint( 3 ) - starting_index_ );
                        t++;
                    }
                }
            }
        }

        void read_vertex_property_names()
        {
            GEO::LineInput in( filename_ );
            while( !in.eof() && in.get_line() )
            {
                in.get_fields();
                if( in.nb_fields() > 0 )
                {
                    if( !in.field_matches( 0, "PROPERTIES" ) )
                    {
                        continue;
                    }
                    vertex_property_names_.reserve( in.nb_fields() - 1 );
                    for( auto prop_name_itr : range( 1, in.nb_fields() ) )
                    {
                        vertex_property_names_.push_back(
                            in.field( prop_name_itr ) );
                    }
                    return; // No need to continue.
                }
            }
        }

        void read_vertex_property_sizes()
        {
            GEO::LineInput in( filename_ );
            while( !in.eof() && in.get_line() )
            {
                in.get_fields();
                if( in.nb_fields() > 0 )
                {
                    if( !in.field_matches( 0, "ESIZES" ) )
                    {
                        continue;
                    }
                    vertex_property_names_.reserve( in.nb_fields() - 1 );
                    for( auto prop_size_itr : range( 1, in.nb_fields() ) )
                    {
                        vertex_attribute_dims_.push_back(
                            in.field_as_uint( prop_size_itr ) );
                    }
                    return; // No need to continue.
                }
            }
        }

        void assign_and_repair_mesh( GEO::Mesh& mesh )
        {
            GEO::coord_index_t dimension =
                static_cast< GEO::coord_index_t >( mesh_dimension_ );
            mesh.facets.assign_triangle_mesh(
                dimension, vertices_, triangles_, true );

            assign_tsurf_properties_to_geogram_mesh( mesh );

            // Do not use GEO::MESH_REPAIR_DEFAULT because it glues the
            // disconnected edges along internal boundaries
            GEO::mesh_repair( mesh, GEO::MESH_REPAIR_DUP_F );
        }

        void assign_tsurf_properties_to_geogram_mesh( GEO::Mesh& mesh )
        {
            for( auto prop_name_itr : range( vertex_property_names_.size() ) )
            {
                index_t nb_dimensions = vertex_attribute_dims_[prop_name_itr];
                GEO::Attribute< double > attr;
                attr.create_vector_attribute( mesh.vertices.attributes(),
                    vertex_property_names_[prop_name_itr], nb_dimensions );
                for( auto v_itr : range( nb_vertices_ ) )
                {
                    for( auto prop_dim_itr : range( nb_dimensions ) )
                    {
                        attr[v_itr * nb_dimensions + prop_dim_itr] =
                            vertex_attributes_[prop_name_itr][v_itr]
                                              [prop_dim_itr];
                    }
                }
            }
        }

        bool is_file_valid()
        {
            GEO::LineInput in( filename_ );
            if( !in.OK() )
            {
                return false;
            }
            else
            {
                return true;
            }
        }

        void allocate_vertices()
        {
            vertices_.resize( mesh_dimension_ * nb_vertices_ );
        }

        void allocate_triangles()
        {
            triangles_.resize( 3 * nb_triangles_ );
        }
        void allocate_vertex_properties()
        {
            vertex_attributes_.resize( vertex_property_names_.size() );
            for( auto vertex_attributes_itr :
                range( vertex_attributes_.size() ) )
            {
                vertex_attributes_[vertex_attributes_itr].resize(
                    nb_vertices_ );
                for( auto vertex_itr : range( nb_vertices_ ) )
                {
                    vertex_attributes_[vertex_attributes_itr][vertex_itr]
                        .resize(
                            vertex_attribute_dims_[vertex_attributes_itr], 0 );
                }
            }
        }

    private:
        index_t starting_index_;
        index_t mesh_dimension_;
        index_t nb_vertices_;
        index_t nb_triangles_;
        int z_sign_;
        std::string filename_;
        GEO::vector< double > vertices_;
        GEO::vector< index_t > triangles_;
        GEO::vector< std::string > vertex_property_names_;
        GEO::vector< index_t > vertex_attribute_dims_;
        GEO::vector< GEO::vector< GEO::vector< double > > > vertex_attributes_;
    };

    class LINMeshIOHandler final : public GEO::MeshIOHandler
    {
    public:
        bool load( const std::string& filename,
            GEO::Mesh& mesh,
            const GEO::MeshIOFlags& flag = GEO::MeshIOFlags() ) final
        {
            ringmesh_unused( flag );
            GEO::LineInput file( filename );

            while( !file.eof() && file.get_line() )
            {
                file.get_fields();
                if( file.nb_fields() > 0 )
                {
                    if( file.field_matches( 0, "v" ) )
                    {
                        vec3 vertex = load_vertex( file, 1 );
                        mesh.vertices.create_vertex( vertex.data() );
                    }
                    else if( file.field_matches( 0, "s" ) )
                    {
                        mesh.edges.create_edge( file.field_as_uint( 1 ) - 1,
                            file.field_as_uint( 2 ) - 1 );
                    }
                }
            }
            return true;
        }
        bool save( const GEO::Mesh& M,
            const std::string& filename,
            const GEO::MeshIOFlags& ioflags = GEO::MeshIOFlags() ) final
        {
            ringmesh_unused( M );
            ringmesh_unused( filename );
            ringmesh_unused( ioflags );
            throw RINGMeshException(
                "I/O", "Saving a Mesh into .lin format not implemented yet" );
            return false;
        }

    private:
        vec3 load_vertex( GEO::LineInput& file, index_t field ) const
        {
            double x = file.field_as_double( field );
            field++;
            double y = file.field_as_double( field );
            field++;
            double z = file.field_as_double( field );
            return vec3( x, y, z );
        }
    };

} // namespace

namespace RINGMesh
{
    /***********************************************************************/
    /* Loading and saving a GEO::Mesh                                      */

    void ringmesh_geogram_mesh_io_initialize()
    {
        geo_register_MeshIOHandler_creator( TSurfMeshIOHandler, "ts" );
        geo_register_MeshIOHandler_creator( LINMeshIOHandler, "lin" );
    }

    double mesh_cell_signed_volume( const GEO::Mesh& M, index_t c )
    {
        double volume = 0;
        switch( M.cells.type( c ) )
        {
        case GEO::MESH_TET:
            volume = GEO::Geom::tetra_signed_volume(
                GEO::Geom::mesh_vertex( M, M.cells.vertex( c, 0 ) ),
                GEO::Geom::mesh_vertex( M, M.cells.vertex( c, 1 ) ),
                GEO::Geom::mesh_vertex( M, M.cells.vertex( c, 2 ) ),
                GEO::Geom::mesh_vertex( M, M.cells.vertex( c, 3 ) ) );
            break;
        case GEO::MESH_PYRAMID:
            volume = GEO::Geom::tetra_signed_volume(
                GEO::Geom::mesh_vertex( M, M.cells.vertex( c, 0 ) ),
                GEO::Geom::mesh_vertex( M, M.cells.vertex( c, 1 ) ),
                GEO::Geom::mesh_vertex( M, M.cells.vertex( c, 2 ) ),
                GEO::Geom::mesh_vertex( M, M.cells.vertex( c, 4 ) ) );
            volume += GEO::Geom::tetra_signed_volume(
                GEO::Geom::mesh_vertex( M, M.cells.vertex( c, 0 ) ),
                GEO::Geom::mesh_vertex( M, M.cells.vertex( c, 2 ) ),
                GEO::Geom::mesh_vertex( M, M.cells.vertex( c, 3 ) ),
                GEO::Geom::mesh_vertex( M, M.cells.vertex( c, 4 ) ) );
            break;
        case GEO::MESH_PRISM:
        case GEO::MESH_HEX:
        {
            vec3 ori( 0, 0, 0 );
            for( auto f : range( M.cells.nb_facets( c ) ) )
            {
                switch( M.cells.facet_nb_vertices( c, f ) )
                {
                case 3:
                    volume += GEO::Geom::tetra_signed_volume(
                        GEO::Geom::mesh_vertex(
                            M, M.cells.facet_vertex( c, f, 0 ) ),
                        GEO::Geom::mesh_vertex(
                            M, M.cells.facet_vertex( c, f, 1 ) ),
                        GEO::Geom::mesh_vertex(
                            M, M.cells.facet_vertex( c, f, 2 ) ),
                        ori );
                    break;
                case 4:
                    volume += GEO::Geom::tetra_signed_volume(
                        GEO::Geom::mesh_vertex(
                            M, M.cells.facet_vertex( c, f, 0 ) ),
                        GEO::Geom::mesh_vertex(
                            M, M.cells.facet_vertex( c, f, 1 ) ),
                        GEO::Geom::mesh_vertex(
                            M, M.cells.facet_vertex( c, f, 2 ) ),
                        ori );
                    volume += GEO::Geom::tetra_signed_volume(
                        GEO::Geom::mesh_vertex(
                            M, M.cells.facet_vertex( c, f, 0 ) ),
                        GEO::Geom::mesh_vertex(
                            M, M.cells.facet_vertex( c, f, 2 ) ),
                        GEO::Geom::mesh_vertex(
                            M, M.cells.facet_vertex( c, f, 3 ) ),
                        ori );
                    break;
                default:
                    ringmesh_assert_not_reached;
                    return 0;
                }
            }
            break;
        }
        default:
            return 0;
        }
        return volume;
    }

    double mesh_cell_volume( const GEO::Mesh& M, index_t c )
    {
        return std::fabs( mesh_cell_signed_volume( M, c ) );
    }

    vec3 mesh_cell_facet_barycenter(
        const GEO::Mesh& M, index_t cell, index_t f )
    {
        vec3 result( 0., 0., 0. );
        index_t nb_vertices = M.cells.facet_nb_vertices( cell, f );
        for( auto v : range( nb_vertices ) )
        {
            result +=
                GEO::Geom::mesh_vertex( M, M.cells.facet_vertex( cell, f, v ) );
        }
        ringmesh_assert( nb_vertices > 0 );

        return result / static_cast< double >( nb_vertices );
    }

    vec3 mesh_cell_barycenter( const GEO::Mesh& M, index_t cell )
    {
        vec3 result( 0.0, 0.0, 0.0 );
        for( auto v : range( M.cells.nb_vertices( cell ) ) )
        {
            result += GEO::Geom::mesh_vertex( M, M.cells.vertex( cell, v ) );
        }
        return ( 1.0 / M.cells.nb_vertices( cell ) ) * result;
    }

} // namespace RINGMesh
