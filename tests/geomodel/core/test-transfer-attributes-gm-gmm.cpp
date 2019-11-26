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

#include <ringmesh/ringmesh_tests_config.h>

#include <future>

#include <geogram/basic/attributes.h>

#include <ringmesh/geomodel/core/geomodel.h>
#include <ringmesh/geomodel/core/geomodel_mesh.h>
#include <ringmesh/geomodel/core/geomodel_mesh_entity.h>

#include <ringmesh/io/io.h>

/*!
 * @author Benjamin Chauvin
 */

namespace
{
    using namespace RINGMesh;

    const std::string attribute_names[6] = { "long_int_attr", "bool_attr",
        "double_attr", "vec3_attr", "dim_6_double_attr", "char_attr" };

    std::mutex lock;

    void load_geomodel( GeoModel3D& in, const std::string& filename )
    {
        std::string input_model_file_name( ringmesh_test_data_path );
        input_model_file_name += filename;

        bool loaded_model_is_valid = geomodel_load( in, input_model_file_name );

        if( !loaded_model_is_valid )
        {
            throw RINGMeshException( "RINGMesh Test",
                "Failed when loading model ", in.name(),
                ": the loaded model is not valid." );
        }
    }

    void assign_vertex_attribute_values( index_t vertex_i,
        const vec3& cur_vertex,
        GEO::Attribute< long int >& vertex_long_int_attr,
        GEO::Attribute< bool >& vertex_bool_attr,
        GEO::Attribute< double >& vertex_double_attr,
        GEO::Attribute< vec3 >& vertex_vec3_attr,
        GEO::Attribute< double >& vertex_dim_6_double_attr,
        GEO::Attribute< char >& vertex_char_attr )
    {
        const long int rounded_vertex_xy =
            std::lrint( cur_vertex.x * cur_vertex.y );
        vertex_long_int_attr[vertex_i] = rounded_vertex_xy;
        vertex_bool_attr[vertex_i] = ( rounded_vertex_xy % 2 == 0 );
        vertex_double_attr[vertex_i] = cur_vertex.x;
        vertex_vec3_attr[vertex_i] = cur_vertex;
        vertex_dim_6_double_attr[vertex_i * 6 + 0] = cur_vertex.x;
        vertex_dim_6_double_attr[vertex_i * 6 + 1] = cur_vertex.y;
        vertex_dim_6_double_attr[vertex_i * 6 + 2] = cur_vertex.z;
        vertex_dim_6_double_attr[vertex_i * 6 + 3] = cur_vertex.x;
        vertex_dim_6_double_attr[vertex_i * 6 + 4] = cur_vertex.y;
        vertex_dim_6_double_attr[vertex_i * 6 + 5] = cur_vertex.z;
        vertex_char_attr[vertex_i] = std::to_string( cur_vertex.y ).data()[0];
    }

    void set_vertex_attributes_on_geomodel_regions( const GeoModel3D& geomodel )
    {
        for( index_t reg_i = 0; reg_i < geomodel.nb_regions(); ++reg_i )
        {
            const Region3D& cur_reg = geomodel.region( reg_i );

            GEO::AttributesManager& reg_attr_mgr =
                cur_reg.vertex_attribute_manager();
            GEO::Attribute< long int > vertex_long_int_attr(
                reg_attr_mgr, attribute_names[0] );
            GEO::Attribute< bool > vertex_bool_attr(
                reg_attr_mgr, attribute_names[1] );
            GEO::Attribute< double > vertex_double_attr(
                reg_attr_mgr, attribute_names[2] );
            GEO::Attribute< vec3 > vertex_vec3_attr(
                reg_attr_mgr, attribute_names[3] );
            GEO::Attribute< double > vertex_dim_6_double_attr;
            vertex_dim_6_double_attr.create_vector_attribute(
                reg_attr_mgr, attribute_names[4], 6 );
            GEO::Attribute< char > vertex_char_attr(
                reg_attr_mgr, attribute_names[5] );

            for( index_t vertex_i = 0; vertex_i < cur_reg.nb_vertices();
                 ++vertex_i )
            {
                const vec3& cur_vertex = cur_reg.vertex( vertex_i );
                assign_vertex_attribute_values( vertex_i, cur_vertex,
                    vertex_long_int_attr, vertex_bool_attr, vertex_double_attr,
                    vertex_vec3_attr, vertex_dim_6_double_attr,
                    vertex_char_attr );
            }
        }
    }

    void set_vertex_attributes_on_geomodelmesh( const GeoModel3D& geomodel )
    {
        const GeoModelMesh3D& gmm = geomodel.mesh;
        const GeoModelMeshVertices3D& gmmv = gmm.vertices;
        GEO::AttributesManager& gmmv_attr_mgr = gmmv.attribute_manager();

        GEO::Attribute< long int > vertex_long_int_attr(
            gmmv_attr_mgr, attribute_names[0] );
        GEO::Attribute< bool > vertex_bool_attr(
            gmmv_attr_mgr, attribute_names[1] );
        GEO::Attribute< double > vertex_double_attr(
            gmmv_attr_mgr, attribute_names[2] );
        GEO::Attribute< vec3 > vertex_vec3_attr(
            gmmv_attr_mgr, attribute_names[3] );
        GEO::Attribute< double > vertex_dim_6_double_attr;
        vertex_dim_6_double_attr.create_vector_attribute(
            gmmv_attr_mgr, attribute_names[4], 6 );
        GEO::Attribute< char > vertex_char_attr(
            gmmv_attr_mgr, attribute_names[5] );

        for( index_t v_i = 0; v_i < gmmv.nb(); ++v_i )
        {
            const vec3& cur_vertex = gmmv.vertex( v_i );
            assign_vertex_attribute_values( v_i, cur_vertex,
                vertex_long_int_attr, vertex_bool_attr, vertex_double_attr,
                vertex_vec3_attr, vertex_dim_6_double_attr, vertex_char_attr );
        }
    }

    void assign_cell_attribute_values( index_t cell_i,
        double cell_volume,
        const vec3& cell_barycenter,
        GEO::Attribute< long int >& cell_long_int_attr,
        GEO::Attribute< bool >& cell_bool_attr,
        GEO::Attribute< double >& cell_double_attr,
        GEO::Attribute< vec3 >& cell_vec3_attr,
        GEO::Attribute< double >& cell_dim_6_double_attr,
        GEO::Attribute< char >& cell_char_attr )
    {
        const long int rounded_volume = std::lrint( cell_volume );
        cell_long_int_attr[cell_i] = rounded_volume;
        cell_bool_attr[cell_i] = ( rounded_volume % 2 == 0 );
        cell_double_attr[cell_i] = cell_volume;
        cell_vec3_attr[cell_i] = cell_barycenter;
        cell_dim_6_double_attr[cell_i * 6 + 0] = cell_vec3_attr[cell_i].x;
        cell_dim_6_double_attr[cell_i * 6 + 1] = cell_vec3_attr[cell_i].y;
        cell_dim_6_double_attr[cell_i * 6 + 2] = cell_vec3_attr[cell_i].z;
        cell_dim_6_double_attr[cell_i * 6 + 3] = cell_vec3_attr[cell_i].x;
        cell_dim_6_double_attr[cell_i * 6 + 4] = cell_vec3_attr[cell_i].y;
        cell_dim_6_double_attr[cell_i * 6 + 5] = cell_vec3_attr[cell_i].z;
        cell_char_attr[cell_i] =
            std::to_string( cell_vec3_attr[cell_i].y ).data()[0];
    }

    void set_cell_attributes_on_geomodel_regions( const GeoModel3D& geomodel )
    {
        for( index_t reg_i = 0; reg_i < geomodel.nb_regions(); ++reg_i )
        {
            const Region3D& cur_reg = geomodel.region( reg_i );

            GEO::AttributesManager& reg_attr_mgr =
                cur_reg.cell_attribute_manager();
            GEO::Attribute< long int > cell_long_int_attr(
                reg_attr_mgr, attribute_names[0] );
            GEO::Attribute< bool > cell_bool_attr(
                reg_attr_mgr, attribute_names[1] );
            GEO::Attribute< double > cell_double_attr(
                reg_attr_mgr, attribute_names[2] );
            GEO::Attribute< vec3 > cell_vec3_attr(
                reg_attr_mgr, attribute_names[3] );
            GEO::Attribute< double > cell_dim_6_double_attr;
            cell_dim_6_double_attr.create_vector_attribute(
                reg_attr_mgr, attribute_names[4], 6 );
            GEO::Attribute< char > cell_char_attr(
                reg_attr_mgr, attribute_names[5] );

            for( index_t cell_i = 0; cell_i < cur_reg.nb_mesh_elements();
                 ++cell_i )
            {
                const double cell_volume = cur_reg.mesh_element_size( cell_i );
                assign_cell_attribute_values( cell_i, cell_volume,
                    cur_reg.mesh_element_barycenter( cell_i ),
                    cell_long_int_attr, cell_bool_attr, cell_double_attr,
                    cell_vec3_attr, cell_dim_6_double_attr, cell_char_attr );
            }
        }
    }

    void set_cell_attributes_on_geomodelmesh( const GeoModel3D& geomodel )
    {
        const GeoModelMesh3D& gmm = geomodel.mesh;
        const GeoModelMeshCells3D& gmmc = gmm.cells;
        GEO::AttributesManager& gmmc_attr_mgr = gmmc.attribute_manager();

        GEO::Attribute< long int > cell_long_int_attr(
            gmmc_attr_mgr, attribute_names[0] );
        GEO::Attribute< bool > cell_bool_attr(
            gmmc_attr_mgr, attribute_names[1] );
        GEO::Attribute< double > cell_double_attr(
            gmmc_attr_mgr, attribute_names[2] );
        GEO::Attribute< vec3 > cell_vec3_attr(
            gmmc_attr_mgr, attribute_names[3] );
        GEO::Attribute< double > cell_dim_6_double_attr;
        cell_dim_6_double_attr.create_vector_attribute(
            gmmc_attr_mgr, attribute_names[4], 6 );
        GEO::Attribute< char > cell_char_attr(
            gmmc_attr_mgr, attribute_names[5] );

        for( index_t cell_i = 0; cell_i < gmmc.nb(); ++cell_i )
        {
            const double cell_volume = gmmc.volume( cell_i );
            assign_cell_attribute_values( cell_i, cell_volume,
                gmmc.barycenter( cell_i ), cell_long_int_attr, cell_bool_attr,
                cell_double_attr, cell_vec3_attr, cell_dim_6_double_attr,
                cell_char_attr );
        }
    }

    void check_attribute_exists( GEO::AttributesManager& attr_mgr,
        const std::string& attr_name,
        const std::string& on_mesh,
        const std::string& stored_on )
    {
        if( !attr_mgr.is_defined( attr_name ) )
        {
            throw RINGMeshException( "RINGMesh Test", "Attribute ", attr_name,
                " does not exist on ", on_mesh, " stored on ",
                stored_on + "." );
        }
    }

    void check_attributes_exist_on_mesh( GEO::AttributesManager& attr_mgr,
        const std::string& on_mesh,
        const std::string& stored_on )
    {
        for( index_t i = 0; i < 6; ++i )
            check_attribute_exists(
                attr_mgr, attribute_names[i], on_mesh, stored_on );
    }

    void check_vertex_attr_transfer_from_geomodel_regions_to_geomodelmesh(
        const GeoModel3D& geomodel )
    {
        const GeoModelMesh3D& gmm = geomodel.mesh;
        const GeoModelMeshVertices3D& gmmv = gmm.vertices;
        GEO::AttributesManager& gmmv_attr_mgr = gmmv.attribute_manager();
        check_attributes_exist_on_mesh(
            gmmv_attr_mgr, "geomodelmesh", "vertices" );

        GEO::Attribute< long int > vertex_long_int_attr(
            gmmv_attr_mgr, attribute_names[0] );
        GEO::Attribute< bool > vertex_bool_attr(
            gmmv_attr_mgr, attribute_names[1] );
        GEO::Attribute< double > vertex_double_attr(
            gmmv_attr_mgr, attribute_names[2] );
        GEO::Attribute< vec3 > vertex_vec3_attr(
            gmmv_attr_mgr, attribute_names[3] );
        GEO::Attribute< double > vertex_dim_6_double_attr(
            gmmv_attr_mgr, attribute_names[4] );
        GEO::Attribute< char > vertex_char_attr(
            gmmv_attr_mgr, attribute_names[5] );

        for( index_t vertex_i = 0; vertex_i < gmmv.nb(); ++vertex_i )
        {
            const vec3& cur_vertex = gmmv.vertex( vertex_i );
            const long int rounded_vertex_xy =
                std::lrint( cur_vertex.x * cur_vertex.y );

            if( rounded_vertex_xy != vertex_long_int_attr[vertex_i] )
            {
                throw RINGMeshException( "RINGMesh Test",
                    "Bad vertex transfer from geomodel region "
                    "to geomodelmesh for long int." );
            }

            bool is_pair = ( rounded_vertex_xy % 2 == 0 );
            if( is_pair != vertex_bool_attr[vertex_i] )
            {
                throw RINGMeshException( "RINGMesh Test",
                    "Bad vertex transfer from geomodel region "
                    "to geomodelmesh for bool." );
            }

            if( std::abs( cur_vertex.x - vertex_double_attr[vertex_i] )
                > geomodel.epsilon() )
            {
                throw RINGMeshException( "RINGMesh Test",
                    "Bad vertex transfer from geomodel region "
                    "to geomodelmesh for double." );
            }

            const vec3 diff = cur_vertex - vertex_vec3_attr[vertex_i];
            if( std::abs( diff.x ) > geomodel.epsilon()
                || std::abs( diff.y ) > geomodel.epsilon()
                || std::abs( diff.z ) > geomodel.epsilon() )
            {
                throw RINGMeshException( "RINGMesh Test",
                    "Bad vertex transfer from geomodel region "
                    "to geomodelmesh for vec3." );
            }

            if( std::abs(
                    cur_vertex.x - vertex_dim_6_double_attr[6 * vertex_i + 0] )
                    > geomodel.epsilon()
                || std::abs( cur_vertex.y
                             - vertex_dim_6_double_attr[6 * vertex_i + 1] )
                       > geomodel.epsilon()
                || std::abs( cur_vertex.z
                             - vertex_dim_6_double_attr[6 * vertex_i + 2] )
                       > geomodel.epsilon()
                || std::abs( cur_vertex.x
                             - vertex_dim_6_double_attr[6 * vertex_i + 3] )
                       > geomodel.epsilon()
                || std::abs( cur_vertex.y
                             - vertex_dim_6_double_attr[6 * vertex_i + 4] )
                       > geomodel.epsilon()
                || std::abs( cur_vertex.z
                             - vertex_dim_6_double_attr[6 * vertex_i + 5] )
                       > geomodel.epsilon() )
            {
                throw RINGMeshException( "RINGMesh Test",
                    "Bad vertex transfer from geomodel region "
                    "to geomodelmesh for double dim 6." );
            }

            const char char_vec3_y = std::to_string( cur_vertex.y ).data()[0];
            if( char_vec3_y != vertex_char_attr[vertex_i] )
            {
                throw RINGMeshException( "RINGMesh Test",
                    "Bad vertex transfer from geomodel region "
                    "to geomodelmesh for char." );
            }
        }
    }

    void check_vertex_attr_transfer_from_geomodelmesh_to_geomodel_regions(
        const GeoModel3D& geomodel )
    {
        for( index_t reg_i = 0; reg_i < geomodel.nb_regions(); ++reg_i )
        {
            const Region3D& cur_reg = geomodel.region( reg_i );
            GEO::AttributesManager& reg_v_attr_mgr =
                cur_reg.vertex_attribute_manager();
            check_attributes_exist_on_mesh(
                reg_v_attr_mgr, "geomodel region", "vertices" );

            GEO::Attribute< long int > vertex_long_int_attr(
                reg_v_attr_mgr, attribute_names[0] );
            GEO::Attribute< bool > vertex_bool_attr(
                reg_v_attr_mgr, attribute_names[1] );
            GEO::Attribute< double > vertex_double_attr(
                reg_v_attr_mgr, attribute_names[2] );
            GEO::Attribute< vec3 > vertex_vec3_attr(
                reg_v_attr_mgr, attribute_names[3] );
            GEO::Attribute< double > vertex_dim_6_double_attr(
                reg_v_attr_mgr, attribute_names[4] );
            GEO::Attribute< char > vertex_char_attr(
                reg_v_attr_mgr, attribute_names[5] );

            for( index_t vertex_i = 0; vertex_i < cur_reg.nb_vertices();
                 ++vertex_i )
            {
                const vec3& cur_vertex = cur_reg.vertex( vertex_i );
                const long int rounded_vertex_xy =
                    std::lrint( cur_vertex.x * cur_vertex.y );

                if( rounded_vertex_xy != vertex_long_int_attr[vertex_i] )
                {
                    throw RINGMeshException( "RINGMesh Test",
                        "Bad vertex transfer from geomodelmesh to geomodel "
                        "region for long int." );
                }

                bool is_pair = ( rounded_vertex_xy % 2 == 0 );
                if( is_pair != vertex_bool_attr[vertex_i] )
                {
                    throw RINGMeshException( "RINGMesh Test",
                        "Bad vertex transfer from geomodelmesh to geomodel "
                        "region for bool." );
                }

                if( std::abs( cur_vertex.x - vertex_double_attr[vertex_i] )
                    > geomodel.epsilon() )
                {
                    throw RINGMeshException( "RINGMesh Test",
                        "Bad vertex transfer from geomodelmesh to geomodel "
                        "region for double." );
                }

                const vec3 diff = cur_vertex - vertex_vec3_attr[vertex_i];
                if( std::abs( diff.x ) > geomodel.epsilon()
                    || std::abs( diff.y ) > geomodel.epsilon()
                    || std::abs( diff.z ) > geomodel.epsilon() )
                {
                    throw RINGMeshException( "RINGMesh Test",
                        "Bad vertex transfer from geomodelmesh to geomodel "
                        "region for vec3." );
                }

                if( std::abs( cur_vertex.x
                              - vertex_dim_6_double_attr[6 * vertex_i + 0] )
                        > geomodel.epsilon()
                    || std::abs( cur_vertex.y
                                 - vertex_dim_6_double_attr[6 * vertex_i + 1] )
                           > geomodel.epsilon()
                    || std::abs( cur_vertex.z
                                 - vertex_dim_6_double_attr[6 * vertex_i + 2] )
                           > geomodel.epsilon()
                    || std::abs( cur_vertex.x
                                 - vertex_dim_6_double_attr[6 * vertex_i + 3] )
                           > geomodel.epsilon()
                    || std::abs( cur_vertex.y
                                 - vertex_dim_6_double_attr[6 * vertex_i + 4] )
                           > geomodel.epsilon()
                    || std::abs( cur_vertex.z
                                 - vertex_dim_6_double_attr[6 * vertex_i + 5] )
                           > geomodel.epsilon() )
                {
                    throw RINGMeshException( "RINGMesh Test",
                        "Bad vertex transfer from geomodelmesh to geomodel "
                        "region for double dim 6." );
                }

                const char char_vec3_y =
                    std::to_string( cur_vertex.y ).data()[0];
                if( char_vec3_y != vertex_char_attr[vertex_i] )
                {
                    throw RINGMeshException( "RINGMesh Test",
                        "Bad vertex transfer from geomodelmesh to geomodel "
                        "region for char." );
                }
            }
        }
    }

    void check_cell_attr_transfer_from_geomodel_regions_to_geomodelmesh(
        const GeoModel3D& geomodel )
    {
        const GeoModelMesh3D& gmm = geomodel.mesh;
        const GeoModelMeshCells3D& gmmc = gmm.cells;
        GEO::AttributesManager& gmmc_attr_mgr = gmmc.attribute_manager();
        check_attributes_exist_on_mesh(
            gmmc_attr_mgr, "geomodelmesh", "cells" );

        GEO::Attribute< long int > cell_long_int_attr(
            gmmc_attr_mgr, attribute_names[0] );
        GEO::Attribute< bool > cell_bool_attr(
            gmmc_attr_mgr, attribute_names[1] );
        GEO::Attribute< double > cell_double_attr(
            gmmc_attr_mgr, attribute_names[2] );
        GEO::Attribute< vec3 > cell_vec3_attr(
            gmmc_attr_mgr, attribute_names[3] );
        GEO::Attribute< double > cell_dim_6_double_attr(
            gmmc_attr_mgr, attribute_names[4] );
        GEO::Attribute< char > cell_char_attr(
            gmmc_attr_mgr, attribute_names[5] );

        for( index_t cell_i = 0; cell_i < gmmc.nb_cells(); ++cell_i )
        {
            const double cell_volume = gmmc.volume( cell_i );
            const long int rounded_volume = std::lrint( cell_volume );

            if( rounded_volume != cell_long_int_attr[cell_i] )
            {
                throw RINGMeshException( "RINGMesh Test",
                    "Bad cell transfer from geomodel region "
                    "to geomodelmesh for long int." );
            }

            bool is_pair = ( rounded_volume % 2 == 0 );
            if( is_pair != cell_bool_attr[cell_i] )
            {
                throw RINGMeshException( "RINGMesh Test",
                    "Bad cell transfer from geomodel region "
                    "to geomodelmesh for bool." );
            }

            if( std::abs( cell_volume - cell_double_attr[cell_i] )
                > geomodel.epsilon3() )
            {
                throw RINGMeshException( "RINGMesh Test",
                    "Bad cell transfer from geomodel region "
                    "to geomodelmesh for double." );
            }

            const vec3 cell_barycenter = gmmc.barycenter( cell_i );
            const vec3 diff = cell_barycenter - cell_vec3_attr[cell_i];
            if( std::abs( diff.x ) > geomodel.epsilon()
                || std::abs( diff.y ) > geomodel.epsilon()
                || std::abs( diff.z ) > geomodel.epsilon() )
            {
                throw RINGMeshException( "RINGMesh Test",
                    "Bad cell transfer from geomodel region "
                    "to geomodelmesh for vec3." );
            }

            if( std::abs(
                    cell_barycenter.x - cell_dim_6_double_attr[6 * cell_i + 0] )
                    > geomodel.epsilon()
                || std::abs( cell_barycenter.y
                             - cell_dim_6_double_attr[6 * cell_i + 1] )
                       > geomodel.epsilon()
                || std::abs( cell_barycenter.z
                             - cell_dim_6_double_attr[6 * cell_i + 2] )
                       > geomodel.epsilon()
                || std::abs( cell_barycenter.x
                             - cell_dim_6_double_attr[6 * cell_i + 3] )
                       > geomodel.epsilon()
                || std::abs( cell_barycenter.y
                             - cell_dim_6_double_attr[6 * cell_i + 4] )
                       > geomodel.epsilon()
                || std::abs( cell_barycenter.z
                             - cell_dim_6_double_attr[6 * cell_i + 5] )
                       > geomodel.epsilon() )
            {
                throw RINGMeshException( "RINGMesh Test",
                    "Bad cell transfer from geomodel region "
                    "to geomodelmesh for double dim 6." );
            }

            const char char_vec3_y =
                std::to_string( cell_barycenter.y ).data()[0];
            if( char_vec3_y != cell_char_attr[cell_i] )
            {
                throw RINGMeshException( "RINGMesh Test",
                    "Bad cell transfer from geomodel region "
                    "to geomodelmesh for char." );
            }
        }
    }

    void check_cell_attr_transfer_from_geomodelmesh_to_geomodel_regions(
        const GeoModel3D& geomodel )
    {
        for( index_t reg_i = 0; reg_i < geomodel.nb_regions(); ++reg_i )
        {
            const Region3D& cur_reg = geomodel.region( reg_i );
            GEO::AttributesManager& reg_c_attr_mgr =
                cur_reg.cell_attribute_manager();
            check_attributes_exist_on_mesh(
                reg_c_attr_mgr, "geomodel region", "cells" );

            GEO::Attribute< long int > cell_long_int_attr(
                reg_c_attr_mgr, attribute_names[0] );
            GEO::Attribute< bool > cell_bool_attr(
                reg_c_attr_mgr, attribute_names[1] );
            GEO::Attribute< double > cell_double_attr(
                reg_c_attr_mgr, attribute_names[2] );
            GEO::Attribute< vec3 > cell_vec3_attr(
                reg_c_attr_mgr, attribute_names[3] );
            GEO::Attribute< double > cell_dim_6_double_attr(
                reg_c_attr_mgr, attribute_names[4] );
            GEO::Attribute< char > cell_char_attr(
                reg_c_attr_mgr, attribute_names[5] );

            for( index_t cell_i = 0; cell_i < cur_reg.nb_mesh_elements();
                 ++cell_i )
            {
                const double cell_volume = cur_reg.mesh_element_size( cell_i );
                const long int rounded_volume = std::lrint( cell_volume );

                if( rounded_volume != cell_long_int_attr[cell_i] )
                {
                    throw RINGMeshException( "RINGMesh Test",
                        "Bad cell transfer from geomodelmesh "
                        "to geomodel region for long int." );
                }

                bool is_pair = ( rounded_volume % 2 == 0 );
                if( is_pair != cell_bool_attr[cell_i] )
                {
                    throw RINGMeshException( "RINGMesh Test",
                        "Bad cell transfer from geomodelmesh "
                        "to geomodel region for bool." );
                }

                if( std::abs( cell_volume - cell_double_attr[cell_i] )
                    > geomodel.epsilon3() )
                {
                    throw RINGMeshException( "RINGMesh Test",
                        "Bad cell transfer from geomodelmesh "
                        "to geomodel region for double." );
                }

                const vec3 cell_barycenter =
                    cur_reg.mesh_element_barycenter( cell_i );
                const vec3 diff = cell_barycenter - cell_vec3_attr[cell_i];
                if( std::abs( diff.x ) > geomodel.epsilon()
                    || std::abs( diff.y ) > geomodel.epsilon()
                    || std::abs( diff.z ) > geomodel.epsilon() )
                {
                    throw RINGMeshException( "RINGMesh Test",
                        "Bad cell transfer from geomodelmesh "
                        "to geomodel region for vec3." );
                }

                if( std::abs( cell_barycenter.x
                              - cell_dim_6_double_attr[6 * cell_i + 0] )
                        > geomodel.epsilon()
                    || std::abs( cell_barycenter.y
                                 - cell_dim_6_double_attr[6 * cell_i + 1] )
                           > geomodel.epsilon()
                    || std::abs( cell_barycenter.z
                                 - cell_dim_6_double_attr[6 * cell_i + 2] )
                           > geomodel.epsilon()
                    || std::abs( cell_barycenter.x
                                 - cell_dim_6_double_attr[6 * cell_i + 3] )
                           > geomodel.epsilon()
                    || std::abs( cell_barycenter.y
                                 - cell_dim_6_double_attr[6 * cell_i + 4] )
                           > geomodel.epsilon()
                    || std::abs( cell_barycenter.z
                                 - cell_dim_6_double_attr[6 * cell_i + 5] )
                           > geomodel.epsilon() )
                {
                    throw RINGMeshException( "RINGMesh Test",
                        "Bad cell transfer from geomodelmesh to geomodel "
                        "region for double dim 6." );
                }

                const char char_vec3_y =
                    std::to_string( cell_barycenter.y ).data()[0];
                if( char_vec3_y != cell_char_attr[cell_i] )
                {
                    throw RINGMeshException( "RINGMesh Test",
                        "Bad cell transfer from geomodelmesh "
                        "to geomodel region for char." );
                }
            }
        }
    }

    void check_attr_transfer_from_geomodel_regions_to_geomodelmesh(
        const GeoModel3D& geomodel )
    {
        check_vertex_attr_transfer_from_geomodel_regions_to_geomodelmesh(
            geomodel );
        check_cell_attr_transfer_from_geomodel_regions_to_geomodelmesh(
            geomodel );
    }

    void check_attr_transfer_from_geomodelmesh_to_geomodel_regions(
        const GeoModel3D& geomodel )
    {
        check_vertex_attr_transfer_from_geomodelmesh_to_geomodel_regions(
            geomodel );
        check_cell_attr_transfer_from_geomodelmesh_to_geomodel_regions(
            geomodel );
    }

    void load_file( GeoModel3D& geomodel )
    {
        std::lock_guard< std::mutex > locking( lock );
        load_geomodel( geomodel, "modelA1_volume_meshed.gm" );
    }

    void tests_transfer_from_geomodel_regions_to_geomodelmesh()
    {
        GeoModel3D geomodel;
        load_file( geomodel );

        set_vertex_attributes_on_geomodel_regions( geomodel );
        set_cell_attributes_on_geomodel_regions( geomodel );
        const GeoModelMesh3D& gmm = geomodel.mesh;
        gmm.transfer_attributes_from_gm_regions_to_gmm();
        check_attr_transfer_from_geomodel_regions_to_geomodelmesh( geomodel );
    }

    void tests_transfer_from_geomodelmesh_to_geomodel_regions()
    {
        GeoModel3D geomodel;
        load_file( geomodel );

        set_vertex_attributes_on_geomodelmesh( geomodel );
        set_cell_attributes_on_geomodelmesh( geomodel );
        const GeoModelMesh3D& gmm = geomodel.mesh;
        gmm.transfer_attributes_from_gmm_to_gm_regions();
        check_attr_transfer_from_geomodelmesh_to_geomodel_regions( geomodel );
    }

    void run_tests()
    {
        // long int is not a default attribute type in geogram.
        GEO::geo_register_attribute_type< long int >( "long int" );

        std::vector< std::future< void > > futures;

        futures.emplace_back( std::async( std::launch::async,
            &tests_transfer_from_geomodel_regions_to_geomodelmesh ) );
        futures.emplace_back( std::async( std::launch::async,
            &tests_transfer_from_geomodelmesh_to_geomodel_regions ) );

        for( auto& future : futures )
        {
            future.wait();
        }
    }
} // namespace

int main()
{
    using namespace RINGMesh;

    try
    {
        Logger::out( "TEST", "Tests of attribute transfer between the geomodel "
                             "and the geomodelmesh" );
        run_tests();
    }
    catch( const RINGMeshException& e )
    {
        Logger::err( e.category(), e.what() );
        return 1;
    }
    catch( const std::exception& e )
    {
        Logger::err( "Exception", e.what() );
        return 1;
    }
    Logger::out( "TEST", "SUCCESS" );
    return 0;
}
