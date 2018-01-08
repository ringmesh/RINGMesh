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

#include <vector>

#include <ringmesh/basic/geometry.h>
#include <ringmesh/geomodel/core/geomodel.h>
#include <ringmesh/geomodel/core/geomodel_mesh_entity.h>

#include <ringmesh/io/io.h>

using namespace RINGMesh;

void error( index_t vertex_id_in_mesh_entity,
    index_t vertex_id_in_geomodel_mesh,
    const gmme_id& mesh_entity_gmme_id )
{
    throw RINGMeshException( "TEST", "Vertex ", vertex_id_in_mesh_entity,
        " in entity ", mesh_entity_gmme_id.type().string(),
        mesh_entity_gmme_id.index(),
        " has not the same coordinates than its equivalent vertex ",
        vertex_id_in_geomodel_mesh, " in the GeoModelMesh" );
}
void test_geomodel_vertices( const GeoModel3D& geomodel )
{
    const GeoModelMeshVertices3D& geomodel_mesh_vertices =
        geomodel.mesh.vertices;
    for( const MeshEntityType& mesh_entity_type :
        geomodel.entity_type_manager().mesh_entity_manager.mesh_entity_types() )
    {
        for( index_t mesh_entity_id :
            range( geomodel.nb_mesh_entities( mesh_entity_type ) ) )
        {
            const GeoModelMeshEntity3D& cur_geomodel_mesh_entity =
                geomodel.mesh_entity(
                    gmme_id( mesh_entity_type, mesh_entity_id ) );
            for( index_t vertex_id_in_mesh_entity :
                range( cur_geomodel_mesh_entity.nb_vertices() ) )
            {
                index_t vertex_id_in_geomodel_mesh =
                    geomodel_mesh_vertices.geomodel_vertex_id(
                        cur_geomodel_mesh_entity.gmme(),
                        vertex_id_in_mesh_entity );
                if( geomodel_mesh_vertices.vertex( vertex_id_in_geomodel_mesh )
                    != cur_geomodel_mesh_entity.vertex(
                           vertex_id_in_mesh_entity ) )
                {
                    error( vertex_id_in_mesh_entity, vertex_id_in_geomodel_mesh,
                        cur_geomodel_mesh_entity.gmme() );
                }
            }
        }
    }
}

void test_GMEVertex( const GeoModel3D& geomodel )
{
    const GeoModelMeshVertices3D& geomodel_mesh_vertices =
        geomodel.mesh.vertices;

    for( index_t vertex_id_in_geomodel_mesh :
        range( geomodel_mesh_vertices.nb() ) )
    {
        std::vector< GMEVertex > vertices_on_geomodel_mesh_entity =
            geomodel_mesh_vertices.gme_vertices( vertex_id_in_geomodel_mesh );
        for( const GMEVertex& cur_vertex_on_geomodel :
            vertices_on_geomodel_mesh_entity )
        {
            const GeoModelMeshEntity3D& cur_geomodel_mesh_entity =
                geomodel.mesh_entity( cur_vertex_on_geomodel.gmme );
            index_t vertex_id_in_mesh_entity = cur_vertex_on_geomodel.v_index;
            if( geomodel_mesh_vertices.vertex( vertex_id_in_geomodel_mesh )
                != cur_geomodel_mesh_entity.vertex( vertex_id_in_mesh_entity ) )
            {
                error( vertex_id_in_mesh_entity, vertex_id_in_geomodel_mesh,
                    cur_geomodel_mesh_entity.gmme() );
            }
        }
    }
}

int main()
{
    using namespace RINGMesh;

    try
    {
        Logger::out( "TEST", "Test GeoModelMeshVertices" );

        std::string input_model_file_name =
            ringmesh_test_data_path + "unit_cube_volume_meshed.gm";

        GeoModel3D in;
        bool loaded_model_is_valid = geomodel_load( in, input_model_file_name );

        if( !loaded_model_is_valid )
        {
            throw RINGMeshException(
                "RINGMesh Test", "Failed when loading model ", in.name() );
        }
        test_geomodel_vertices( in );
        test_GMEVertex( in );
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
