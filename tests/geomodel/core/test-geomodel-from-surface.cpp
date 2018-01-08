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

#include <geogram/basic/command_line.h>

#include <geogram/mesh/mesh_io.h>

#include <ringmesh/geomodel/builder/geomodel_builder_from_mesh.h>
#include <ringmesh/geomodel/core/geomodel.h>
#include <ringmesh/geomodel/core/geomodel_api.h>
#include <ringmesh/geomodel/core/geomodel_mesh_entity.h>
#include <ringmesh/geomodel/tools/geomodel_validity.h>
#include <ringmesh/io/io.h>

#include <ringmesh/mesh/mesh_index.h>

/*!
 * Test the creation of a GeoModel from a conformal surface mesh
 * @todo Test on other datasets: nested spheres.
 * @author Jeanne Pellerin
 */

int main()
{
    using namespace RINGMesh;

    try
    {
        Logger::out( "TEST", "Test GeoModel building from Surface" );

        std::vector< std::future< void > > futures;

        GEO::Mesh in;
        GEO::mesh_load( ringmesh_test_data_path + "modelA6.mesh", in );
        GeoModel3D geomodel;

        GeoModelBuilderSurfaceMesh builder( geomodel, in );
        builder.build_polygonal_surfaces_from_connected_components();
        builder.build_lines_and_corners_from_surfaces();
        builder.build_regions_from_lines_and_surfaces();
        builder.end_geomodel();

        // Checking the validity of loaded model
        ValidityCheckMode validity_mode = ValidityCheckMode::ALL;

#ifdef RINGMESH_DEBUG
        validity_mode =
            validity_mode ^ ValidityCheckMode::POLYGON_INTERSECTIONS;
#endif

        futures.emplace_back(
            std::async( std::launch::async, [&geomodel, &validity_mode] {
                if( !is_geomodel_valid( geomodel, validity_mode ) )
                {
                    throw RINGMeshException( "RINGMesh Test",
                        "Failed when loading model ", geomodel.name(),
                        ": the loaded model is not valid." );
                }
            } ) );

        GEO::Mesh surface_meshes;
        // Compute mesh with duplicated points to compares number
        // of mesh elements and mesh entities
        for( const auto& surface : geomodel.surfaces() )
        {
            index_t vertex_it{ surface_meshes.vertices.create_vertices(
                surface.nb_vertices() ) };
            for( index_t v : range( surface.nb_vertices() ) )
            {
                surface_meshes.vertices.point( vertex_it + v ) =
                    surface.vertex( v );
            }
            index_t facet_it{ surface_meshes.facets.create_triangles(
                surface.nb_mesh_elements() ) };
            for( index_t f : range( surface.nb_mesh_elements() ) )
            {
                for( index_t v :
                    range( surface.nb_mesh_element_vertices( f ) ) )
                {
                    surface_meshes.facets.set_vertex( facet_it + f, v,
                        vertex_it
                            + surface.mesh_element_vertex_index( { f, v } ) );
                }
            }
        }
        surface_meshes.facets.connect();

        futures.emplace_back(
            std::async( std::launch::async, [&surface_meshes] {
                // Save computed mesh
                std::string output_file2( ringmesh_test_output_path );
                output_file2 += "saved_modelA6_dupl_points.mesh";
                GEO::mesh_save( surface_meshes, output_file2 );
            } ) );

        GeoModel3D reloaded_model;
        GeoModelBuilderSurfaceMesh builder2( reloaded_model, surface_meshes );
        builder2.build_polygonal_surfaces_from_connected_components();
        builder2.build_lines_and_corners_from_surfaces();
        builder2.build_regions_from_lines_and_surfaces();
        builder2.end_geomodel();

        futures.emplace_back(
            std::async( std::launch::async, [&reloaded_model, &validity_mode] {
                // Checking if building has been successfully done
                if( !is_geomodel_valid( reloaded_model, validity_mode ) )
                {
                    throw RINGMeshException( "RINGMesh Test",
                        "Failed when reloading model ", reloaded_model.name(),
                        ": the reloaded model is not valid." );
                }
            } ) );

        // Checking number of mesh elements
        if( surface_meshes.vertices.nb() != in.vertices.nb() )
        {
            throw RINGMeshException( "RINGMesh Test",
                "Error when building model: not same number of vertices ",
                "than input mesh." );
        }
        if( surface_meshes.facets.nb() != in.facets.nb() )
        {
            throw RINGMeshException( "RINGMesh Test",
                "Error when building model: not same number of facets ",
                "than input mesh." );
        }
        if( surface_meshes.cells.nb() != in.cells.nb() )
        {
            throw RINGMeshException( "RINGMesh Test",
                "Error when building model: not same number of cells ",
                "than input mesh." );
        }

        // Checking number of GeoModelMeshEntities
        if( reloaded_model.nb_corners() != geomodel.nb_corners() )
        {
            throw RINGMeshException( "RINGMesh Test",
                "Error when reload model: not same number of corners ",
                "between saved model and reload model." );
        }
        if( reloaded_model.nb_lines() != geomodel.nb_lines() )
        {
            throw RINGMeshException( "RINGMesh Test",
                "Error when reload model: not same number of lines ",
                "between saved model and reload model." );
        }
        if( reloaded_model.nb_surfaces() != geomodel.nb_surfaces() )
        {
            throw RINGMeshException( "RINGMesh Test",
                "Error when reload model: not same number of surfaces ",
                "between saved model and reload model." );
        }
        if( reloaded_model.nb_regions() != geomodel.nb_regions() )
        {
            throw RINGMeshException( "RINGMesh Test",
                "Error when reload model: not same number of regions ",
                "between saved model and reload model." );
        }

        for( auto& future : futures )
        {
            future.wait();
        }
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
