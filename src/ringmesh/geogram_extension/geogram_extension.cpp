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

#include <ringmesh/geogram_extension/geogram_extension.h>

#include <fstream>

#include <geogram/basic/file_system.h>
#include <geogram/basic/line_stream.h>
#include <geogram/basic/logger.h>

#include <geogram/mesh/mesh.h>
#include <geogram/mesh/mesh_geometry.h>
#include <geogram/mesh/mesh_io.h>
#include <geogram/mesh/mesh_repair.h>

/*!
 * @todo Re-orgarnize this mess
 */

namespace {
    using namespace RINGMesh ;

    /*!
     * @return the cotangent (inverse of the tangent) at point b
     * in the triangle abc.
     */
    double cotangent( const vec3& a, const vec3& b, const vec3& c )
    {
        const vec3 ba = a - b ;
        const vec3 bc = c - b ;
        const vec3 bc_cross_ba = GEO::cross( bc, ba ) ;
        const double bc_dot_ba = GEO::dot( bc, ba ) ;
        ringmesh_assert( bc_cross_ba.length() > global_epsilon ) ;
        return bc_dot_ba / bc_cross_ba.length() ;
    }

    bool is_attribute_a_double(
        GEO::AttributesManager& att_manager,
        const std::string& att_name )
    {
        return GEO::Attribute< double >::is_defined( att_manager, att_name ) ;
    }
}

namespace RINGMesh {

    using GEO::vec3 ;
    using GEO::Mesh ;

    /***********************************************************************/
    /* Loading and saving a GEO::Mesh                                      */

    /*! 
     * @brief TSurfMeshIOHandler for importing .ts files into a mesh.
     */
    class TSurfMeshIOHandler: public GEO::MeshIOHandler {
    public:
        TSurfMeshIOHandler()
            :
                starting_index_( 1 ),
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
        virtual bool load(
            const std::string& filename,
            GEO::Mesh& mesh,
            const GEO::MeshIOFlags& flag = GEO::MeshIOFlags() )
        {
            ringmesh_unused( flag ) ;
            filename_ = filename ;
            if( !is_file_valid() ) {
                return false ;
            } else {
                read_number_of_vertices_and_triangles() ;
                read_vertex_property_names() ;
                read_vertex_property_sizes() ;
                allocate_vertices() ;
                allocate_triangles() ;
                allocate_vertex_properties() ;
                read_vertices_and_triangles() ;
                assign_and_repair_mesh( mesh ) ;
                return true ;
            }
        }

        /*!
         * @brief Save a Mesh in .ts format
         */
        virtual bool save(
            const GEO::Mesh& mesh,
            const std::string& filename,
            const GEO::MeshIOFlags& flag = GEO::MeshIOFlags() )
        {
            ringmesh_unused( flag ) ;
            if( !mesh.facets.are_simplices() ) {
                throw RINGMeshException( "I/O",
                    "Cannot save a non triangulated mesh into TSurf format" ) ;
            }

            std::ofstream out( filename.c_str() ) ;
            out.precision( 16 ) ;
            save_header( out, GEO::FileSystem::base_name( filename ) ) ;
            std::vector< std::string > att_v_double_names ;
            std::vector< index_t > vertex_attr_dims ;
            fill_vertex_attribute_header( mesh, out, att_v_double_names,
                vertex_attr_dims ) ;

            out << "TFACE" << std::endl ;
            save_vertices( out, mesh, att_v_double_names, vertex_attr_dims ) ;
            save_triangles( out, mesh ) ;
            out << "END" << std::endl ;

            return true ;
        }

    private:
        void save_vertices(
            std::ofstream& out,
            const GEO::Mesh& mesh,
            const std::vector< std::string >& att_v_double_names,
            const std::vector< index_t >& vertex_attr_dims )
        {
            for( index_t v = 0; v < mesh.vertices.nb(); v++ ) {
                // PVRTX must be used instead of VRTX because
                // properties are not read by Gocad if it is VRTX.
                out << "PVRTX " << v + starting_index_ << " "
                    << mesh.vertices.point( v ) ;
                for( index_t attr_dbl_itr = 0;
                    attr_dbl_itr < att_v_double_names.size(); ++attr_dbl_itr ) {
                    GEO::ReadOnlyScalarAttributeAdapter cur_attr(
                        mesh.vertices.attributes(),
                        att_v_double_names[attr_dbl_itr] ) ;
                    for( index_t dim_itr = 0;
                        dim_itr < vertex_attr_dims[attr_dbl_itr]; ++dim_itr ) {
                        out << " "
                            << cur_attr[v * vertex_attr_dims[attr_dbl_itr]
                                + dim_itr] ;
                    }
                }
                out << std::endl ;
            }
        }
        void save_triangles( std::ofstream& out, const GEO::Mesh& mesh )
        {
            for( index_t t = 0; t < mesh.facets.nb(); t++ ) {
                out << "TRGL" ;
                for( index_t v = 0; v < 3; v++ ) {
                    out << " " << mesh.facets.vertex( t, v ) + starting_index_ ;
                }
                out << std::endl ;
            }
        }
        void save_header( std::ofstream& out, const std::string& mesh_name )
        {
            out << "GOCAD TSurf" << std::endl ;
            out << "HEADER {" << std::endl ;
            out << "name: " << mesh_name << std::endl ;
            out << "}" << std::endl ;
        }
        void fill_vertex_attribute_header(
            const GEO::Mesh& mesh,
            std::ofstream& out,
            std::vector< std::string >& att_v_double_names,
            std::vector< index_t >& vertex_attr_dims ) const
        {
            GEO::vector< std::string > att_v_names ;
            GEO::AttributesManager& mesh_vertex_mgr = mesh.vertices.attributes() ;
            mesh_vertex_mgr.list_attribute_names( att_v_names ) ;
            for( index_t att_v = 0; att_v < mesh_vertex_mgr.nb();
                att_v++ ) {

                if( att_v_names[att_v] == "point" ) {
                    continue ;
                }

                if( !GEO::ReadOnlyScalarAttributeAdapter::is_defined( mesh_vertex_mgr,
                    att_v_names[att_v] ) ) {
                    continue ;
                }
                att_v_double_names.push_back( att_v_names[att_v] ) ;
                index_t cur_dim =
                    mesh_vertex_mgr.find_attribute_store(
                        att_v_names[att_v] )->dimension() ;
                vertex_attr_dims.push_back( cur_dim ) ;
            }

            if( !att_v_double_names.empty() ) {
                out << "PROPERTIES" ;
                for( index_t attr_dbl_itr = 0;
                    attr_dbl_itr < att_v_double_names.size(); ++attr_dbl_itr ) {
                    out << " " << att_v_double_names[attr_dbl_itr] ;
                }
                out << std::endl ;
                out << "PROP_LEGAL_RANGES" ;
                for( index_t attr_dbl_itr = 0;
                    attr_dbl_itr < att_v_double_names.size(); ++attr_dbl_itr ) {
                    out << " **none**  **none**" ;
                }
                out << std::endl ;
                out << "NO_DATA_VALUES" ;
                for( index_t attr_dbl_itr = 0;
                    attr_dbl_itr < att_v_double_names.size(); ++attr_dbl_itr ) {
                    out << " -99999" ;
                }
                out << std::endl ;
                out << "READ_ONLY" ;
                for( index_t attr_dbl_itr = 0;
                    attr_dbl_itr < att_v_double_names.size(); ++attr_dbl_itr ) {
                    out << " 1" ;
                }
                out << std::endl ;
                out << "PROPERTY_CLASSES" ;
                for( index_t attr_dbl_itr = 0;
                    attr_dbl_itr < att_v_double_names.size(); ++attr_dbl_itr ) {
                    out << " " << att_v_double_names[attr_dbl_itr] ;
                }
                out << std::endl ;
                out << "PROPERTY_KINDS" ;
                for( index_t attr_dbl_itr = 0;
                    attr_dbl_itr < att_v_double_names.size(); ++attr_dbl_itr ) {
                    out << " \"Real Number\"" ;
                }
                out << std::endl ;
                out << "PROPERTY_SUBCLASSES" ;
                for( index_t attr_dbl_itr = 0;
                    attr_dbl_itr < att_v_double_names.size(); ++attr_dbl_itr ) {
                    out << " QUANTITY Float" ;
                }
                out << std::endl ;
                out << "ESIZES" ;
                for( index_t attr_dbl_itr = 0;
                    attr_dbl_itr < att_v_double_names.size(); ++attr_dbl_itr ) {
                    out << " "
                        << GEO::String::to_string( vertex_attr_dims[attr_dbl_itr] ) ;
                }
                out << std::endl ;
                out << "UNITS" ;
                for( index_t attr_dbl_itr = 0;
                    attr_dbl_itr < att_v_double_names.size(); ++attr_dbl_itr ) {
                    out << " unitless" ;
                }
                out << std::endl ;
                for( index_t attr_dbl_itr = 0;
                    attr_dbl_itr < att_v_double_names.size(); ++attr_dbl_itr ) {
                    out << "PROPERTY_CLASS_HEADER "
                        << att_v_double_names[attr_dbl_itr] << " {" << std::endl ;
                    out << "kind: Real Number" << std::endl ;
                    out << "unit: unitless" << std::endl ;
                    out << "}" << std::endl ;
                }
            }
        }
        // This function read the z_sign too [PA]
        void read_number_of_vertices_and_triangles()
        {
            GEO::LineInput in( filename_ ) ;
            while( !in.eof() && in.get_line() ) {
                in.get_fields() ;
                if( in.nb_fields() > 0 ) {
                    if( in.field_matches( 0, "ZPOSITIVE" ) ) {
                        if( in.field_matches( 1, "Elevation" ) ) {
                            z_sign_ = 1 ;
                        } else if( in.field_matches( 1, "Depth" ) ) {
                            z_sign_ = -1 ;
                        }
                    } else if( in.field_matches( 0, "VRTX" )
                        || in.field_matches( 0, "PVRTX" ) ) {
                        nb_vertices_++ ;
                    } else if( in.field_matches( 0, "PATOM" )
                        || in.field_matches( 0, "ATOM" ) ) {
                        nb_vertices_++ ;
                    } else if( in.field_matches( 0, "TRGL" ) ) {
                        nb_triangles_++ ;
                    }
                }
            }
        }

        void read_vertices_and_triangles()
        {
            GEO::LineInput in( filename_ ) ;
            index_t v = 0 ;
            index_t t = 0 ;
            while( !in.eof() && in.get_line() ) {
                in.get_fields() ;
                if( in.nb_fields() > 0 ) {
                    if( in.field_matches( 0, "VRTX" )
                        || in.field_matches( 0, "PVRTX" ) ) {
                        vertices_[mesh_dimension_ * v] = in.field_as_double( 2 ) ;
                        vertices_[mesh_dimension_ * v + 1] = in.field_as_double(
                            3 ) ;
                        vertices_[mesh_dimension_ * v + 2] = in.field_as_double( 4 )
                            * z_sign_ ;
                        if( in.field_matches( 0, "PVRTX" ) )
                        {
                            index_t offset = 5 ;
                            for( index_t prop_name_itr = 0;
                                prop_name_itr < vertex_property_names_.size();
                                ++prop_name_itr ) {
                                for( index_t v_attr_dim_itr = 0;
                                    v_attr_dim_itr
                                        < vertex_attribute_dims_[prop_name_itr];
                                    ++v_attr_dim_itr ) {
                                    vertex_attributes_[prop_name_itr][v][v_attr_dim_itr] = in.field_as_double( offset ) ;
                                    ++offset ;
                                }
                            }
                        }
                        ++v ;
                    } else if( in.field_matches( 0, "PATOM" )
                        || in.field_matches( 0, "ATOM" ) ) {
                        index_t v0 = in.field_as_uint( 2 ) - starting_index_ ;
                        vertices_[mesh_dimension_ * v] = vertices_[mesh_dimension_
                            * v0] ;
                        vertices_[mesh_dimension_ * v + 1] =
                            vertices_[mesh_dimension_ * v0 + 1] ;
                        vertices_[mesh_dimension_ * v + 2] =
                            vertices_[mesh_dimension_ * v0 + 2] ;
                        ++v ;
                    } else if( in.field_matches( 0, "TRGL" ) ) {
                        triangles_[3 * t] = index_t(
                            in.field_as_uint( 1 ) - starting_index_ ) ;
                        triangles_[3 * t + 1] = index_t(
                            in.field_as_uint( 2 ) - starting_index_ ) ;
                        triangles_[3 * t + 2] = index_t(
                            in.field_as_uint( 3 ) - starting_index_ ) ;
                        t++ ;
                    }
                }
            }
        }

        void read_vertex_property_names()
        {
            GEO::LineInput in( filename_ ) ;
            while( !in.eof() && in.get_line() ) {
                in.get_fields() ;
                if( in.nb_fields() > 0 ) {
                    if( !in.field_matches( 0, "PROPERTIES" ) ) {
                        continue ;
                    }
                    vertex_property_names_.reserve( in.nb_fields() -1 ) ;
                    for( index_t prop_name_itr = 1; prop_name_itr < in.nb_fields();
                        ++prop_name_itr ) {
                        vertex_property_names_.push_back(
                            in.field( prop_name_itr ) ) ;
                    }
                    return ; // No need to continue.
                }
            }
        }

        void read_vertex_property_sizes()
        {
            GEO::LineInput in( filename_ ) ;
            while( !in.eof() && in.get_line() ) {
                in.get_fields() ;
                if( in.nb_fields() > 0 ) {
                    if( !in.field_matches( 0, "ESIZES" ) ) {
                        continue ;
                    }
                    vertex_property_names_.reserve( in.nb_fields() - 1 ) ;
                    for( index_t prop_size_itr = 1; prop_size_itr < in.nb_fields();
                        ++prop_size_itr ) {
                        vertex_attribute_dims_.push_back(
                            in.field_as_uint( prop_size_itr ) ) ;
                    }
                    return ; // No need to continue.
                }
            }
        }

        void assign_and_repair_mesh( GEO::Mesh& mesh )
        {
            GEO::coord_index_t dimension =
                static_cast< GEO::coord_index_t >( mesh_dimension_ ) ;
            mesh.facets.assign_triangle_mesh( dimension, vertices_, triangles_,
                true ) ;

            assign_tsurf_properties_to_geogram_mesh( mesh ) ;

            // Do not use GEO::MESH_REPAIR_DEFAULT because it glues the
            // disconnected edges along internal boundaries
            GEO::mesh_repair( mesh, GEO::MESH_REPAIR_DUP_F ) ;
        }

        void assign_tsurf_properties_to_geogram_mesh( GEO::Mesh& mesh )
        {
            for( index_t prop_name_itr = 0;
                prop_name_itr < vertex_property_names_.size(); ++prop_name_itr ) {
                GEO::Attribute< double > attr ;
                attr.create_vector_attribute( mesh.vertices.attributes(),
                    vertex_property_names_[prop_name_itr],
                    vertex_attribute_dims_[prop_name_itr] ) ;
                for( index_t v_itr = 0; v_itr < nb_vertices_; ++v_itr ) {
                    for( index_t prop_dim_itr = 0;
                        prop_dim_itr < vertex_attribute_dims_[prop_name_itr];
                        ++prop_dim_itr ) {
                        attr[v_itr * vertex_attribute_dims_[prop_name_itr]
                            + prop_dim_itr] =
                            vertex_attributes_[prop_name_itr][v_itr][prop_dim_itr] ;
                    }
                }
            }
        }

        bool is_file_valid()
        {
            GEO::LineInput in( filename_ ) ;
            if( !in.OK() ) {
                return false ;
            } else {
                return true ;
            }
        }

        void allocate_vertices()
        {
            vertices_.resize( mesh_dimension_ * nb_vertices_ ) ;
        }

        void allocate_triangles()
        {
            triangles_.resize( 3 * nb_triangles_ ) ;
        }
        void allocate_vertex_properties()
        {
            vertex_attributes_.resize( vertex_property_names_.size() ) ;
            for( index_t vertex_attributes_itr = 0;
                vertex_attributes_itr < vertex_attributes_.size();
                ++vertex_attributes_itr ) {
                vertex_attributes_[vertex_attributes_itr].resize(
                    nb_vertices_ ) ;
                for( index_t vertex_itr = 0;
                    vertex_itr < nb_vertices_;
                    ++vertex_itr ) {
                    vertex_attributes_[vertex_attributes_itr][vertex_itr].resize(
                        vertex_attribute_dims_[vertex_attributes_itr], 0 ) ;
                }
            }
        }

    private:
        index_t starting_index_ ;
        index_t mesh_dimension_ ;
        index_t nb_vertices_ ;
        index_t nb_triangles_ ;
        int z_sign_ ;
        std::string filename_ ;
        GEO::vector< double > vertices_ ;
        GEO::vector< index_t > triangles_ ;
        GEO::vector< std::string > vertex_property_names_ ;
        GEO::vector< index_t > vertex_attribute_dims_ ;
        GEO::vector< GEO::vector< GEO::vector< double > > > vertex_attributes_ ;
    } ;

    class LINMeshIOHandler: public GEO::MeshIOHandler {
    public:
        virtual bool load(
            const std::string& filename,
            GEO::Mesh& mesh,
            const GEO::MeshIOFlags& flag = GEO::MeshIOFlags() )
        {
            ringmesh_unused( flag ) ;
            GEO::LineInput file( filename ) ;

            while( !file.eof() && file.get_line() ) {
                file.get_fields() ;
                if( file.nb_fields() > 0 ) {
                    if( file.field_matches( 0, "v" ) ) {
                        vec3 vertex = load_vertex( file, 1 ) ;
                        mesh.vertices.create_vertex( vertex.data() ) ;
                    } else if( file.field_matches( 0, "s" ) ) {
                        mesh.edges.create_edge( file.field_as_uint( 1 ) - 1,
                            file.field_as_uint( 2 ) - 1 ) ;
                    }
                }
            }
            return true ;

        }
        virtual bool save(
            const GEO::Mesh& M,
            const std::string& filename,
            const GEO::MeshIOFlags& ioflags = GEO::MeshIOFlags() )
        {
            ringmesh_unused( M ) ;
            ringmesh_unused( filename ) ;
            ringmesh_unused( ioflags ) ;
            throw RINGMeshException( "I/O",
                "Saving a Mesh into .lin format not implemented yet" ) ;
            return false ;
        }

    private:
        vec3 load_vertex( GEO::LineInput& file, index_t field ) const
        {
            double x = file.field_as_double( field++ ) ;
            double y = file.field_as_double( field++ ) ;
            double z = file.field_as_double( field++ ) ;
            return vec3( x, y, z ) ;
        }
    } ;

    void ringmesh_mesh_io_initialize()
    {
        geo_register_MeshIOHandler_creator( TSurfMeshIOHandler, "ts" ) ;
        geo_register_MeshIOHandler_creator( LINMeshIOHandler, "lin" ) ;
    }

    /*!
     * Computes the volume of a Mesh cell
     * @param[in] M the mesh
     * @param[in] c the cell index
     * @return the volume of the cell
     */
    double mesh_cell_signed_volume( const GEO::Mesh& M, index_t c )
    {
        double volume = 0 ;
        switch( M.cells.type( c ) ) {
            case GEO::MESH_TET:
                volume = GEO::Geom::tetra_signed_volume(
                    GEO::Geom::mesh_vertex( M, M.cells.vertex( c, 0 ) ),
                    GEO::Geom::mesh_vertex( M, M.cells.vertex( c, 1 ) ),
                    GEO::Geom::mesh_vertex( M, M.cells.vertex( c, 2 ) ),
                    GEO::Geom::mesh_vertex( M, M.cells.vertex( c, 3 ) ) ) ;
                break ;
            case GEO::MESH_PYRAMID:
                volume = GEO::Geom::tetra_signed_volume(
                    GEO::Geom::mesh_vertex( M, M.cells.vertex( c, 0 ) ),
                    GEO::Geom::mesh_vertex( M, M.cells.vertex( c, 1 ) ),
                    GEO::Geom::mesh_vertex( M, M.cells.vertex( c, 2 ) ),
                    GEO::Geom::mesh_vertex( M, M.cells.vertex( c, 4 ) ) ) ;
                volume += GEO::Geom::tetra_signed_volume(
                    GEO::Geom::mesh_vertex( M, M.cells.vertex( c, 0 ) ),
                    GEO::Geom::mesh_vertex( M, M.cells.vertex( c, 2 ) ),
                    GEO::Geom::mesh_vertex( M, M.cells.vertex( c, 3 ) ),
                    GEO::Geom::mesh_vertex( M, M.cells.vertex( c, 4 ) ) ) ;
                break ;
            case GEO::MESH_PRISM:
            case GEO::MESH_HEX: {
                vec3 ori( 0, 0, 0 ) ;
                for( index_t f = 0; f < M.cells.nb_facets( c ); f++ ) {
                    switch( M.cells.facet_nb_vertices( c, f ) ) {
                        case 3:
                            volume += GEO::Geom::tetra_signed_volume(
                                GEO::Geom::mesh_vertex( M,
                                    M.cells.facet_vertex( c, f, 0 ) ),
                                GEO::Geom::mesh_vertex( M,
                                    M.cells.facet_vertex( c, f, 1 ) ),
                                GEO::Geom::mesh_vertex( M,
                                    M.cells.facet_vertex( c, f, 2 ) ), ori ) ;
                            break ;
                        case 4:
                            volume += GEO::Geom::tetra_signed_volume(
                                GEO::Geom::mesh_vertex( M,
                                    M.cells.facet_vertex( c, f, 0 ) ),
                                GEO::Geom::mesh_vertex( M,
                                    M.cells.facet_vertex( c, f, 1 ) ),
                                GEO::Geom::mesh_vertex( M,
                                    M.cells.facet_vertex( c, f, 2 ) ), ori ) ;
                            volume += GEO::Geom::tetra_signed_volume(
                                GEO::Geom::mesh_vertex( M,
                                    M.cells.facet_vertex( c, f, 0 ) ),
                                GEO::Geom::mesh_vertex( M,
                                    M.cells.facet_vertex( c, f, 2 ) ),
                                GEO::Geom::mesh_vertex( M,
                                    M.cells.facet_vertex( c, f, 3 ) ), ori ) ;
                            break ;
                        default:
                            ringmesh_assert_not_reached ;
                            return 0 ;
                    }
                }
                break ;
            }
            default:
                return 0 ;
        }
        return volume ;
    }

    double mesh_cell_volume( const GEO::Mesh& M, index_t c )
    {
        return std::fabs( mesh_cell_signed_volume( M, c ) ) ;
    }
    /*!
     * Computes the Mesh cell facet barycenter
     * @param[in] M the mesh
     * @param[in] cell the cell index
     * @param[in] f the facet index in the cell
     * @return the cell facet center
     */
    vec3 mesh_cell_facet_barycenter( const GEO::Mesh& M, index_t cell, index_t f )
    {
        vec3 result( 0., 0., 0. ) ;
        index_t nb_vertices = M.cells.facet_nb_vertices( cell, f ) ;
        for( index_t v = 0; v < nb_vertices; ++v ) {
            result += GEO::Geom::mesh_vertex( M,
                M.cells.facet_vertex( cell, f, v ) ) ;
        }
        ringmesh_assert( nb_vertices > 0 ) ;

        return result / static_cast< double >( nb_vertices ) ;
    }

    /*!
     * Computes the non weighted barycenter of a volumetric
     * cell of a Mesh
     * @param[in] M the mesh
     * @param[in] cell the cell index
     * @return the cell center
     */
    vec3 mesh_cell_barycenter( const GEO::Mesh& M, index_t cell )
    {
        vec3 result( 0.0, 0.0, 0.0 ) ;
        double count = 0.0 ;
        for( index_t v = 0; v < M.cells.nb_vertices( cell ); ++v ) {
            result += GEO::Geom::mesh_vertex( M, M.cells.vertex( cell, v ) ) ;
            count += 1.0 ;
        }
        return ( 1.0 / count ) * result ;
    }

    /*!
     * Get vertices when an edge is divide into \p nb_parts parts
     * @param[in] mesh the mesh
     * @param[in] edge the edge id in \p mesh
     * @param[in] nb_parts the number of edge division
     * @param[out] points the points which divide the edge
     */
    void divide_edge_in_parts(
        const GEO::Mesh& mesh,
        index_t edge,
        index_t nb_parts,
        std::vector< vec3 >& points )
    {
        points.resize( nb_parts - 1 ) ;
        double pond = 1. / nb_parts ;
        vec3 node0 = GEO::Geom::mesh_vertex( mesh, mesh.edges.vertex( edge, 0 ) ) ;
        vec3 node1 = GEO::Geom::mesh_vertex( mesh, mesh.edges.vertex( edge, 1 ) ) ;
        for( index_t i = 0; i < nb_parts - 1; i++ ) {
            for( index_t j = 0; j < 3; j++ ) {
                points[i][j] = ( i + 1 ) * pond * node1[j]
                    + ( 1. - ( i + 1 ) * pond ) * node0[j] ;
            }
        }

    }

    void divide_edge_in_parts(
        vec3& node0,
        vec3& node1,
        index_t nb_parts,
        std::vector< vec3 >& points )
    {
        if( nb_parts > 0 ) {
            points.resize( nb_parts - 1 ) ;
            double pond = 1. / nb_parts ;
            for( index_t i = 0; i < nb_parts - 1; i++ ) {
                for( index_t j = 0; j < 3; j++ ) {
                    points[i][j] = ( i + 1 ) * pond * node1[j]
                        + ( 1. - ( i + 1 ) * pond ) * node0[j] ;
                }
            }
        }
    }
    /*!
     * Gets the closest local vertex index in a mesh cell of a point
     * @param[in] mesh the mesh
     * @param[in] p the point to test
     * @param[in] t the cell index
     * @return the local vertex index
     */
    index_t get_nearest_vertex_index(
        const GEO::Mesh& mesh,
        const vec3& p,
        index_t t )
    {
        double dist = GEO::Numeric::max_float64() ;
        index_t result = NO_ID ;
        for( index_t v = 0; v < mesh.cells.nb_vertices( t ); v++ ) {
            double distance = length2(
                GEO::Geom::mesh_vertex( mesh, mesh.cells.vertex( t, v ) ) - p ) ;
            if( distance < dist ) {
                result = v ;
            }
        }
        return result ;
    }

    void barycentric_coordinates_point_inside_mesh_facet(
        const vec3& point_inside_facet,
        const GEO::Mesh& mesh,
        index_t facet,
        std::vector< double >& barycentric_coordinates )
    {
        const GEO::MeshVertices& mesh_vertices = mesh.vertices ;
        const GEO::MeshFacets& mesh_facets = mesh.facets ;
        ringmesh_assert( facet < mesh_facets.nb() ) ;
        const index_t nb_vertices = mesh_facets.nb_vertices( facet ) ;
        barycentric_coordinates.resize( nb_vertices, 0. ) ;
        double sum = 0. ;
        for( index_t facet_vertex_itr = 0;
            facet_vertex_itr < mesh_facets.nb_vertices( facet );
            ++facet_vertex_itr ) {
            const index_t cur_vertex_id = mesh_facets.vertex( facet,
                facet_vertex_itr ) ;
            const vec3& cur_vertex_vec = mesh_vertices.point( cur_vertex_id ) ;

            const index_t prev_local_id = mesh_facets.prev_corner_around_facet(
                facet, cur_vertex_id ) ;
            const index_t prev_id = mesh_facets.vertex( facet, prev_local_id ) ;
            const vec3& prev_vec = mesh_vertices.point( prev_id ) ;
            const index_t next_local_id = mesh_facets.next_corner_around_facet(
                facet, cur_vertex_id ) ;
            const index_t next_id = mesh_facets.vertex( facet, next_local_id ) ;
            const vec3& next_vec = mesh_vertices.point( next_id ) ;

            const double numerator = ( cotangent( point_inside_facet, cur_vertex_vec,
                prev_vec )
                + cotangent( point_inside_facet, cur_vertex_vec, next_vec ) ) ;
            const double denominator =
                ( point_inside_facet - cur_vertex_vec ).length2() ;
            ringmesh_assert( denominator > global_epsilon ) ;
            barycentric_coordinates[facet_vertex_itr] = numerator / denominator ;
            sum += barycentric_coordinates[facet_vertex_itr] ;
        }

        ringmesh_assert( sum > global_epsilon ) ;
        for( index_t i = 0; i < nb_vertices; ++i ) {
            barycentric_coordinates[i] /= sum ;
        }
    }

    /*!
     * Rotation of all the vertices of a mesh following
     * a defined rotational matrix.
     *
     * @param mesh[in,out] the mesh to rotate.
     *
     * @param[in] rot_mat matrix which defines the rotation.
     */
    void rotate_mesh( GEO::Mesh& mesh, const GEO::Matrix< 4, double >& rot_mat )
    {
        for( index_t v = 0; v < mesh.vertices.nb(); v++ ) {
            double old_coords[4] ;
            for( index_t i = 0; i < 3; i++ ) {
                old_coords[i] = mesh.vertices.point_ptr( v )[i] ;
            }
            old_coords[3] = 1. ;
            double new_coords[4] ;
            GEO::mult( rot_mat, old_coords, new_coords ) ;

            for( index_t i = 0; i < 3; i++ ) {
                mesh.vertices.point_ptr( v )[i] = new_coords[i] ;
            }
            ringmesh_assert( new_coords[ 3 ] == 1. ) ;
        }
    }

    void print_bounded_attributes( const GEO::Mesh& M )
    {
        {
            GEO::vector< std::string > names ;
            M.vertices.attributes().list_attribute_names( names ) ;
            if( !names.empty() ) {
                std::vector< bool > is_bounded( names.size(), false ) ;
                bool failed = false ;
                for( index_t a = 0; a < names.size(); a++ ) {
                    if( names[a] == "point" ) continue ;
                    is_bounded[a] = M.vertices.attributes().find_attribute_store(
                        names[a] )->has_observers() ;
                    if( is_bounded[a] ) failed = true ;
                }
                if( failed ) {
                    Logger::err( "Attributes" )
                        << "Attributes still bounded on vertices:" ;
                    for( index_t a = 0; a < names.size(); a++ ) {
                        if( is_bounded[a] ) {
                            Logger::err( "Attributes" ) << " " << names[a] ;
                        }
                    }
                    Logger::err( "Attributes" ) << std::endl ;
                }
            }
        }
        {
            GEO::vector< std::string > names ;
            M.edges.attributes().list_attribute_names( names ) ;
            if( !names.empty() ) {
                std::vector< bool > is_bounded( names.size() ) ;
                bool failed = false ;
                for( index_t a = 0; a < names.size(); a++ ) {
                    is_bounded[a] = M.edges.attributes().find_attribute_store(
                        names[a] )->has_observers() ;
                    if( is_bounded[a] ) failed = true ;
                }
                if( failed ) {
                    Logger::err( "Attributes" )
                        << "Attributes still bounded on edges:" ;
                    for( index_t a = 0; a < names.size(); a++ ) {
                        if( is_bounded[a] ) {
                            Logger::err( "Attributes" ) << " " << names[a] ;
                        }
                    }
                    Logger::err( "Attributes" ) << std::endl ;
                }
            }
        }
        {
            GEO::vector< std::string > names ;
            M.facets.attributes().list_attribute_names( names ) ;
            if( !names.empty() ) {
                std::vector< bool > is_bounded( names.size() ) ;
                bool failed = false ;
                for( index_t a = 0; a < names.size(); a++ ) {
                    is_bounded[a] = M.facets.attributes().find_attribute_store(
                        names[a] )->has_observers() ;
                    if( is_bounded[a] ) failed = true ;
                }
                if( failed ) {
                    Logger::err( "Attributes" )
                        << "Attributes still bounded on facets:" ;
                    for( index_t a = 0; a < names.size(); a++ ) {
                        if( is_bounded[a] ) {
                            Logger::err( "Attributes" ) << " " << names[a] ;
                        }
                    }
                    Logger::err( "Attributes" ) << std::endl ;
                }
            }
        }
        {
            GEO::vector< std::string > names ;
            M.facet_corners.attributes().list_attribute_names( names ) ;
            if( !names.empty() ) {
                std::vector< bool > is_bounded( names.size() ) ;
                bool failed = false ;
                for( index_t a = 0; a < names.size(); a++ ) {
                    is_bounded[a] =
                        M.facet_corners.attributes().find_attribute_store( names[a] )->has_observers() ;
                    if( is_bounded[a] ) failed = true ;
                }
                if( failed ) {
                    Logger::err( "Attributes" )
                        << "Attributes still bounded on facet_corners:" ;
                    for( index_t a = 0; a < names.size(); a++ ) {
                        if( is_bounded[a] ) {
                            Logger::err( "Attributes" ) << " " << names[a] ;
                        }
                    }
                    Logger::err( "Attributes" ) << std::endl ;
                }
            }
        }
        {
            GEO::vector< std::string > names ;
            M.cells.attributes().list_attribute_names( names ) ;
            if( !names.empty() ) {
                std::vector< bool > is_bounded( names.size() ) ;
                bool failed = false ;
                for( index_t a = 0; a < names.size(); a++ ) {
                    is_bounded[a] = M.cells.attributes().find_attribute_store(
                        names[a] )->has_observers() ;
                    if( is_bounded[a] ) {
                        failed = true ;
                    }
                }
                if( failed ) {
                    Logger::err( "Attributes" )
                        << "Attributes still bounded on cells:" ;
                    for( index_t a = 0; a < names.size(); a++ ) {
                        if( is_bounded[a] ) {
                            Logger::err( "Attributes" ) << " " << names[a] ;
                        }
                    }
                    Logger::err( "Attributes" ) << std::endl ;
                }
            }
        }
        {
            GEO::vector< std::string > names ;
            M.cell_corners.attributes().list_attribute_names( names ) ;
            if( !names.empty() ) {
                std::vector< bool > is_bounded( names.size() ) ;
                bool failed = false ;
                for( index_t a = 0; a < names.size(); a++ ) {
                    is_bounded[a] = M.cell_corners.attributes().find_attribute_store(
                        names[a] )->has_observers() ;
                    if( is_bounded[a] ) failed = true ;
                }
                if( failed ) {
                    Logger::err( "Attributes" )
                        << "Attributes still bounded on cell_corners:" ;
                    for( index_t a = 0; a < names.size(); a++ ) {
                        if( is_bounded[a] ) {
                            Logger::err( "Attributes" ) << " " << names[a] ;
                        }
                    }
                    Logger::err( "Attributes" ) << std::endl ;
                }
            }
        }
        {
            GEO::vector< std::string > names ;
            M.cell_facets.attributes().list_attribute_names( names ) ;
            if( !names.empty() ) {
                std::vector< bool > is_bounded( names.size() ) ;
                bool failed = false ;
                for( index_t a = 0; a < names.size(); a++ ) {
                    is_bounded[a] = M.cell_facets.attributes().find_attribute_store(
                        names[a] )->has_observers() ;
                    if( is_bounded[a] ) {
                        failed = true ;
                    }
                }
                if( failed ) {
                    Logger::err( "Attributes" )
                        << "Attributes still bounded on cell_facets:" ;
                    for( index_t a = 0; a < names.size(); a++ ) {
                        if( is_bounded[a] ) {
                            Logger::err( "Attributes" ) << " " << names[a] ;
                        }
                    }
                    Logger::err( "Attributes" ) << std::endl ;
                }
            }
        }
    }

}
