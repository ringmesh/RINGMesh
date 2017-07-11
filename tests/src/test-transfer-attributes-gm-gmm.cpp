/*
 * Copyright (c) 2012-2017, Association Scientifique pour la Geologie et ses Applications (ASGA)
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

#include <ringmesh/ringmesh_tests_config.h>

#include <ringmesh/geomodel/geomodel.h>
#include <ringmesh/geomodel/geomodel_mesh.h>
#include <ringmesh/io/io.h>

/*!
 * @author Benjamin Chauvin
 */

namespace {
    using namespace RINGMesh;

    void load_geomodel( GeoModel& in, const std::string& filename )
    {
        std::string input_model_file_name( ringmesh_test_data_path );
        input_model_file_name += filename;

        bool loaded_model_is_valid = geomodel_load( in, input_model_file_name );

        if( !loaded_model_is_valid ) {
            throw RINGMeshException( "RINGMesh Test",
                "Failed when loading model " + in.name()
                    + ": the loaded model is not valid." );
        }
    }

    void set_vertex_attributes_on_geomodel_regions( const GeoModel& geomodel )
    {
        for( index_t reg_i = 0; reg_i < geomodel.nb_regions(); ++reg_i ) {
            const Region& cur_reg = geomodel.region( reg_i );

            GEO::AttributesManager& reg_attr_mgr =
                cur_reg.vertex_attribute_manager();
            GEO::Attribute< long int > vertex_long_int_attr( reg_attr_mgr,
                "vertex_long_int_attr" );
            GEO::Attribute< bool > vertex_bool_attr( reg_attr_mgr,
                "vertex_bool_attr" );
            GEO::Attribute< double > vertex_double_attr( reg_attr_mgr,
                "vertex_double_attr" );
            GEO::Attribute< vec3 > vertex_vec3_attr( reg_attr_mgr,
                "vertex_vec3_attr" );
            GEO::Attribute< double > vertex_dim_6_double_attr;
            vertex_dim_6_double_attr.create_vector_attribute( reg_attr_mgr,
                "vertex_dim_6_double_attr", 6 );
            GEO::Attribute< char > vertex_char_attr( reg_attr_mgr,
                "vertex_char_attr" );

            for( index_t vertex_i = 0; vertex_i < cur_reg.nb_vertices();
                ++vertex_i ) {
                const vec3& cur_vertex = cur_reg.vertex( vertex_i );
                const long int rounded_vertex_xy = std::lrint(
                    cur_vertex.x * cur_vertex.y );
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
                vertex_char_attr[vertex_i] =
                    std::to_string( cur_vertex.y ).data()[0];
            }
        }
    }

    void set_cell_attributes_on_geomodel_regions( const GeoModel& geomodel )
    {
        for( index_t reg_i = 0; reg_i < geomodel.nb_regions(); ++reg_i ) {
            const Region& cur_reg = geomodel.region( reg_i );

            GEO::AttributesManager& reg_attr_mgr = cur_reg.cell_attribute_manager();
            GEO::Attribute< long int > cell_long_int_attr( reg_attr_mgr,
                "cell_long_int_attr" );
            GEO::Attribute< bool > cell_bool_attr( reg_attr_mgr, "cell_bool_attr" );
            GEO::Attribute< double > cell_double_attr( reg_attr_mgr,
                "cell_double_attr" );
            GEO::Attribute< vec3 > cell_vec3_attr( reg_attr_mgr, "cell_vec3_attr" );
            GEO::Attribute< double > cell_dim_6_double_attr;
            cell_dim_6_double_attr.create_vector_attribute( reg_attr_mgr,
                "cell_dim_6_double_attr", 6 );
            GEO::Attribute< char > cell_char_attr( reg_attr_mgr, "cell_char_attr" );

            for( index_t cell_i = 0; cell_i < cur_reg.nb_mesh_elements();
                ++cell_i ) {
                const double cell_volume = cur_reg.mesh_element_size( cell_i );
                const long int rounded_volume = std::lrint( cell_volume );
                cell_long_int_attr[cell_i] = rounded_volume;
                cell_bool_attr[cell_i] = ( rounded_volume % 2 == 0 );
                cell_double_attr[cell_i] = cell_volume;
                cell_vec3_attr[cell_i] = cur_reg.mesh_element_barycenter( cell_i );
                cell_dim_6_double_attr[cell_i * 6 + 0] = cell_vec3_attr[cell_i].x;
                cell_dim_6_double_attr[cell_i * 6 + 1] = cell_vec3_attr[cell_i].y;
                cell_dim_6_double_attr[cell_i * 6 + 2] = cell_vec3_attr[cell_i].z;
                cell_dim_6_double_attr[cell_i * 6 + 3] = cell_vec3_attr[cell_i].x;
                cell_dim_6_double_attr[cell_i * 6 + 4] = cell_vec3_attr[cell_i].y;
                cell_dim_6_double_attr[cell_i * 6 + 5] = cell_vec3_attr[cell_i].z;
                cell_char_attr[cell_i] =
                    std::to_string( cell_vec3_attr[cell_i].y ).data()[0];
            }
        }
    }

    void check_attribute_exists(
        GEO::AttributesManager& attr_mgr,
        const std::string& attr_name,
        const std::string& on_mesh )
    {
        if( !attr_mgr.is_defined( attr_name ) ) {
            throw RINGMeshException( "RINGMesh Test",
                "Attribute " + attr_name + " does not exist on " + on_mesh + "." );
        }
    }

    void check_attributes_exist_on_geomodel( GEO::AttributesManager& gmmv_attr_mgr )
    {
        const std::string on_mesh = "geomodel";
        /// TODO store attribute in string to avoid to write them several times...
        check_attribute_exists( gmmv_attr_mgr, "vertex_long_int_attr", on_mesh );
        check_attribute_exists( gmmv_attr_mgr, "vertex_bool_attr", on_mesh );
        check_attribute_exists( gmmv_attr_mgr, "vertex_double_attr", on_mesh );
        check_attribute_exists( gmmv_attr_mgr, "vertex_vec3_attr", on_mesh );
        check_attribute_exists( gmmv_attr_mgr, "vertex_dim_6_double_attr", on_mesh );
        check_attribute_exists( gmmv_attr_mgr, "vertex_char_attr", on_mesh );
    }

    void check_attributes_exist_on_geomodelmesh(
        GEO::AttributesManager& gmmc_attr_mgr )
    {
        const std::string on_mesh = "geomodelmesh";
        /// TODO store attribute in string to avoid to write them several times...
        check_attribute_exists( gmmc_attr_mgr, "cell_long_int_attr", on_mesh );
        check_attribute_exists( gmmc_attr_mgr, "cell_bool_attr", on_mesh );
        check_attribute_exists( gmmc_attr_mgr, "cell_double_attr", on_mesh );
        check_attribute_exists( gmmc_attr_mgr, "cell_vec3_attr", on_mesh );
        check_attribute_exists( gmmc_attr_mgr, "cell_dim_6_double_attr", on_mesh );
        check_attribute_exists( gmmc_attr_mgr, "cell_char_attr", on_mesh );
    }

    void check_vertex_attr_transfer_from_geomodel_regions_to_geomodelmesh(
        const GeoModel& geomodel )
    {
        const GeoModelMesh& gmm = geomodel.mesh;
        const GeoModelMeshVertices& gmmv = gmm.vertices;
        GEO::AttributesManager& gmmv_attr_mgr = gmmv.attribute_manager();
        check_attributes_exist_on_geomodel( gmmv_attr_mgr );

        GEO::Attribute< long int > vertex_long_int_attr( gmmv_attr_mgr,
            "vertex_long_int_attr" );
        GEO::Attribute< bool > vertex_bool_attr( gmmv_attr_mgr, "vertex_bool_attr" );
        GEO::Attribute< double > vertex_double_attr( gmmv_attr_mgr,
            "vertex_double_attr" );
        GEO::Attribute< vec3 > vertex_vec3_attr( gmmv_attr_mgr, "vertex_vec3_attr" );
        GEO::Attribute< double > vertex_dim_6_double_attr( gmmv_attr_mgr,
            "vertex_dim_6_double_attr" );
        GEO::Attribute< char > vertex_char_attr( gmmv_attr_mgr, "vertex_char_attr" );

        for( index_t vertex_i = 0; vertex_i < gmmv.nb(); ++vertex_i ) {
            const vec3& cur_vertex = gmmv.vertex( vertex_i );
            const long int rounded_vertex_xy = std::lrint(
                cur_vertex.x * cur_vertex.y );

            if( rounded_vertex_xy != vertex_long_int_attr[vertex_i] ) {
                throw RINGMeshException( "RINGMesh Test", "Bad transfer" ); //TODO improve message
            }

            bool is_pair = ( rounded_vertex_xy % 2 == 0 );
            if( is_pair != vertex_bool_attr[vertex_i] ) {
                throw RINGMeshException( "RINGMesh Test", "Bad transfer" ); //TODO improve message
            }

            if( std::abs( cur_vertex.x - vertex_double_attr[vertex_i] )
                > geomodel.epsilon() ) {
                throw RINGMeshException( "RINGMesh Test", "Bad transfer" ); //TODO improve message
            }

            const vec3 diff = cur_vertex - vertex_vec3_attr[vertex_i];
            if( std::abs( diff.x ) > geomodel.epsilon()
                || std::abs( diff.y ) > geomodel.epsilon()
                || std::abs( diff.z ) > geomodel.epsilon() ) {
                throw RINGMeshException( "RINGMesh Test", "Bad transfer" ); //TODO improve message
            }

            if( std::abs( cur_vertex.x - vertex_dim_6_double_attr[6 * vertex_i + 0] )
                > geomodel.epsilon()
                || std::abs(
                    cur_vertex.y - vertex_dim_6_double_attr[6 * vertex_i + 1] )
                    > geomodel.epsilon()
                || std::abs(
                    cur_vertex.z - vertex_dim_6_double_attr[6 * vertex_i + 2] )
                    > geomodel.epsilon()
                || std::abs(
                    cur_vertex.x - vertex_dim_6_double_attr[6 * vertex_i + 3] )
                    > geomodel.epsilon()
                || std::abs(
                    cur_vertex.y - vertex_dim_6_double_attr[6 * vertex_i + 4] )
                    > geomodel.epsilon()
                || std::abs(
                    cur_vertex.z - vertex_dim_6_double_attr[6 * vertex_i + 5] )
                    > geomodel.epsilon() ) {
                throw RINGMeshException( "RINGMesh Test", "Bad transfer" ); //TODO improve message
            }

            const char char_vec3_y = std::to_string( cur_vertex.y ).data()[0];
            if( char_vec3_y != vertex_char_attr[vertex_i] ) {
                throw RINGMeshException( "RINGMesh Test", "Bad transfer" ); //TODO improve message
            }
        }
    }

    void check_cell_attr_transfer_from_geomodel_regions_to_geomodelmesh(
        const GeoModel& geomodel )
    {
        const GeoModelMesh& gmm = geomodel.mesh;
        const GeoModelMeshCells& gmmc = gmm.cells;
        GEO::AttributesManager& gmmc_attr_mgr = gmmc.attribute_manager();
        check_attributes_exist_on_geomodelmesh( gmmc_attr_mgr );

        GEO::Attribute< long int > cell_long_int_attr( gmmc_attr_mgr,
            "cell_long_int_attr" );
        GEO::Attribute< bool > cell_bool_attr( gmmc_attr_mgr, "cell_bool_attr" );
        GEO::Attribute< double > cell_double_attr( gmmc_attr_mgr,
            "cell_double_attr" );
        GEO::Attribute< vec3 > cell_vec3_attr( gmmc_attr_mgr, "cell_vec3_attr" );
        GEO::Attribute< double > cell_dim_6_double_attr( gmmc_attr_mgr,
            "cell_dim_6_double_attr" );
        GEO::Attribute< char > cell_char_attr( gmmc_attr_mgr, "cell_char_attr" );

        for( index_t cell_i = 0; cell_i < gmmc.nb_cells(); ++cell_i ) {
            const double cell_volume = gmmc.volume( cell_i );
            const long int rounded_volume = std::lrint( cell_volume );

            if( rounded_volume != cell_long_int_attr[cell_i] ) {
                throw RINGMeshException( "RINGMesh Test", "Bad transfer" ); //TODO improve message
            }

            bool is_pair = ( rounded_volume % 2 == 0 );
            if( is_pair != cell_bool_attr[cell_i] ) {
                throw RINGMeshException( "RINGMesh Test", "Bad transfer" ); //TODO improve message
            }

            if( std::abs( cell_volume - cell_double_attr[cell_i] )
                > geomodel.epsilon3() ) {
                throw RINGMeshException( "RINGMesh Test", "Bad transfer" ); //TODO improve message
            }

            const vec3 cell_barycenter = gmmc.barycenter( cell_i );
            const vec3 diff = cell_barycenter - cell_vec3_attr[cell_i];
            if( std::abs( diff.x ) > geomodel.epsilon()
                || std::abs( diff.y ) > geomodel.epsilon()
                || std::abs( diff.z ) > geomodel.epsilon() ) {
                throw RINGMeshException( "RINGMesh Test", "Bad transfer" ); //TODO improve message
            }

            if( std::abs(
                cell_barycenter.x - cell_dim_6_double_attr[6 * cell_i + 0] )
                > geomodel.epsilon()
                || std::abs(
                    cell_barycenter.y - cell_dim_6_double_attr[6 * cell_i + 1] )
                    > geomodel.epsilon()
                || std::abs(
                    cell_barycenter.z - cell_dim_6_double_attr[6 * cell_i + 2] )
                    > geomodel.epsilon()
                || std::abs(
                    cell_barycenter.x - cell_dim_6_double_attr[6 * cell_i + 3] )
                    > geomodel.epsilon()
                || std::abs(
                    cell_barycenter.y - cell_dim_6_double_attr[6 * cell_i + 4] )
                    > geomodel.epsilon()
                || std::abs(
                    cell_barycenter.z - cell_dim_6_double_attr[6 * cell_i + 5] )
                    > geomodel.epsilon() ) {
                throw RINGMeshException( "RINGMesh Test", "Bad transfer" ); //TODO improve message
            }

            const char char_vec3_y = std::to_string( cell_barycenter.y ).data()[0];
            if( char_vec3_y != cell_char_attr[cell_i] ) {
                throw RINGMeshException( "RINGMesh Test", "Bad transfer" ); //TODO improve message
            }
        }
    }

    void check_attr_transfer_from_geomodel_regions_to_geomodelmesh(
        const GeoModel& geomodel )
    {
        check_vertex_attr_transfer_from_geomodel_regions_to_geomodelmesh( geomodel );
        check_cell_attr_transfer_from_geomodel_regions_to_geomodelmesh( geomodel );

    }

    void tests_transfer_from_geomodel_regions_to_geomodelmesh()
    {
        GeoModel geomodel;
        load_geomodel( geomodel, "modelA1_volume_meshed.gm" );

        set_vertex_attributes_on_geomodel_regions( geomodel );
        set_cell_attributes_on_geomodel_regions( geomodel );
        const GeoModelMesh& gmm = geomodel.mesh;
        gmm.transfer_attributes_from_gm_regions_to_gmm();
        check_attr_transfer_from_geomodel_regions_to_geomodelmesh( geomodel );
    }

    void tests_transfer_from_geomodelmesh_to_geomodel_regions()
    {

    }

    void run_tests()
    {
        // long int is not a default attribute type in geogram.
        GEO::geo_register_attribute_type< long int >( "long int" );
        tests_transfer_from_geomodel_regions_to_geomodelmesh();
        tests_transfer_from_geomodelmesh_to_geomodel_regions();
    }

}

int main()
{
    using namespace RINGMesh;

    try {
        default_configure();

        Logger::out( "TEST", "Test IO for a GeoModel in .gm" );
        run_tests();

    } catch( const RINGMeshException& e ) {
        Logger::err( e.category(), e.what() );
        return 1;
    } catch( const std::exception& e ) {
        Logger::err( "Exception", e.what() );
        return 1;
    }
    Logger::out( "TEST", "SUCCESS" );
    return 0;
}
