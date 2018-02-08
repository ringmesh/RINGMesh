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

#include <array>
#include <iomanip>
#include <iostream>

#include <ringmesh/basic/geometry.h>
#include <ringmesh/geomodel/core/geomodel.h>
#include <ringmesh/geomodel/core/geomodel_api.h>
#include <ringmesh/geomodel/core/geomodel_entity.h>
#include <ringmesh/geomodel/core/geomodel_geological_entity.h>
#include <ringmesh/geomodel/core/geomodel_mesh_entity.h>
#include <ringmesh/mesh/mesh.h>

/*!
 * @file Set of high level API functions
 */

namespace
{
    using namespace RINGMesh;

    template < index_t DIMENSION >
    index_t count_geomodel_edges( const GeoModel< DIMENSION >& geomodel )
    {
        index_t result{ 0 };
        for( const auto& line : geomodel.lines() )
        {
            result += line.nb_mesh_elements();
        }
        return result;
    }

    template < index_t DIMENSION >
    index_t count_geomodel_polygons( const GeoModel< DIMENSION >& geomodel )
    {
        index_t result{ 0 };
        for( const auto& surface : geomodel.surfaces() )
        {
            result += surface.nb_mesh_elements();
        }
        return result;
    }

    template < index_t DIMENSION >
    index_t count_geomodel_cells( const GeoModel< DIMENSION >& geomodel );

    template <>
    index_t count_geomodel_cells( const GeoModel2D& geomodel )
    {
        ringmesh_unused( geomodel );
        return 0;
    }

    template <>
    index_t count_geomodel_cells( const GeoModel3D& geomodel )
    {
        index_t nb_cells{ 0 };
        for( const auto& region : geomodel.regions() )
        {
            nb_cells += region.nb_mesh_elements();
        }
        return nb_cells;
    }

    template < index_t DIMENSION >
    std::tuple< double, double, double, double, double >
        compute_region_volumes_per_cell_type(
            const Region< DIMENSION >& region )
    {
        double tet_volume{ 0 };
        double pyramid_volume{ 0 };
        double prism_volume{ 0 };
        double hex_volume{ 0 };
        double poly_volume{ 0 };
        for( auto c : range( region.nb_mesh_elements() ) )
        {
            index_t nb_vertices{ region.nb_mesh_element_vertices( c ) };
            double volume{ region.mesh().cell_volume( c ) };
            switch( nb_vertices )
            {
            case 4:
                tet_volume += volume;
                break;
            case 5:
                pyramid_volume += volume;
                break;
            case 6:
                prism_volume += volume;
                break;
            case 8:
                hex_volume += volume;
                break;
            default:
                poly_volume += volume;
                break;
            }
        }
        return std::make_tuple(
            tet_volume, pyramid_volume, prism_volume, hex_volume, poly_volume );
    }

    template < index_t DIMENSION >
    std::tuple< double, double, double, double, double, double >
        compute_geomodel_volumes_per_cell_type(
            const GeoModel< DIMENSION >& geomodel )
    {
        double tet_volume{ 0 };
        double pyramid_volume{ 0 };
        double prism_volume{ 0 };
        double hex_volume{ 0 };
        double poly_volume{ 0 };

        for( const auto& region : geomodel.regions() )
        {
            std::tie( tet_volume, pyramid_volume, prism_volume, hex_volume,
                poly_volume ) = compute_region_volumes_per_cell_type( region );
        }
        double total{ tet_volume + pyramid_volume + prism_volume + hex_volume
                      + poly_volume };
        return std::make_tuple( total, tet_volume, pyramid_volume, prism_volume,
            hex_volume, poly_volume );
    }

    void print_one_cell_stat( const index_t& nb_cell,
        const index_t& nb_cell_total,
        const std::string& cell_type )
    {
        Logger::out( "GeoModel", "* ", nb_cell, " ", cell_type, " (",
            nb_cell * 100 / nb_cell_total, "%)" );
    }

    void print_one_cell_volume_stat( const double& cell_volume,
        const double& cell_volume_total,
        const std::string& cell_type )
    {
        Logger::out( "GeoModel", "* ", cell_type, " volume ", cell_volume, " (",
            static_cast< index_t >(
                cell_volume * 100 / cell_volume_total + 0.5 ),
            "%)" );
    }
    template < index_t DIMENSION >
    void print_geomodel_base_mesh_stats( const GeoModel< DIMENSION >& geomodel )
    {
        Logger::out( "GeoModel", "Model ", geomodel.name(), " is made of\n",
            std::setw( 10 ), std::left, geomodel.mesh.vertices.nb(),
            " vertices\n", std::setw( 10 ), std::left,
            count_geomodel_edges( geomodel ), " edges" );

        index_t nb_triangles{ geomodel.mesh.polygons.nb_triangle() };
        index_t nb_quads{ geomodel.mesh.polygons.nb_quad() };
        index_t nb_unclassified_polygons{
            geomodel.mesh.polygons.nb_unclassified_polygon()
        };
        index_t nb_polygons{ geomodel.mesh.polygons.nb_polygons() };
        Logger::out(
            "GeoModel", std::setw( 10 ), std::left, nb_polygons, " polygons" );
        if( nb_triangles > 0 )
        {
            print_one_cell_stat( nb_triangles, nb_polygons, "triangles" );
        }
        if( nb_quads > 0 )
        {
            print_one_cell_stat( nb_quads, nb_polygons, "quads" );
        }
        if( nb_unclassified_polygons > 0 )
        {
            print_one_cell_stat( nb_unclassified_polygons, nb_polygons,
                "unclassified polygons" );
        }
    }
} // namespace

namespace RINGMesh
{
    template < index_t DIMENSION >
    void print_nb_mesh_entities(
        const GeoModel< DIMENSION >& geomodel, const MeshEntityType& type )
    {
        Logger::out( "GeoModel", std::setw( 10 ), std::left,
            geomodel.nb_mesh_entities( type ), " ", type );
    }

    template < index_t DIMENSION >
    void print_nb_geological_entities( const GeoModel< DIMENSION >& geomodel,
        const GeologicalEntityType& type )
    {
        if( geomodel.nb_geological_entities( type ) == 0 )
        {
            return;
        }
        Logger::out( "GeoModel", std::setw( 10 ), std::left,
            geomodel.nb_geological_entities( type ), " ", type );
    }

    template < index_t DIMENSION >
    void print_geomodel( const GeoModel< DIMENSION >& geomodel )
    {
        Logger::out( "GeoModel", "Model ", geomodel.name(), " has\n",
            std::setw( 10 ), std::left, geomodel.mesh.vertices.nb(),
            " vertices\n", std::setw( 10 ), std::left,
            count_geomodel_edges( geomodel ), " edges\n", std::setw( 10 ),
            std::left, count_geomodel_polygons( geomodel ), " polygons" );
        index_t nb_cells{ count_geomodel_cells( geomodel ) };
        if( nb_cells != 0 )
        {
            Logger::out(
                "GeoModel", std::setw( 10 ), std::left, nb_cells, " cells" );
        }
        Logger::out( "GeoModel" );

        const EntityTypeManager< DIMENSION >& manager =
            geomodel.entity_type_manager();
        const std::vector< MeshEntityType >& mesh_entity_types =
            manager.mesh_entity_manager.mesh_entity_types();
        for( const MeshEntityType& type : mesh_entity_types )
        {
            print_nb_mesh_entities( geomodel, type );
        }
        const std::vector< GeologicalEntityType >& geological_entity_types =
            manager.geological_entity_manager.geological_entity_types();
        for( const GeologicalEntityType& type : geological_entity_types )
        {
            print_nb_geological_entities( geomodel, type );
        }
    }

    template <>
    void print_geomodel_mesh_stats( const GeoModel2D& geomodel )
    {
        print_geomodel_base_mesh_stats( geomodel );
    }

    template <>
    void print_geomodel_mesh_stats( const GeoModel3D& geomodel )
    {
        print_geomodel_base_mesh_stats( geomodel );

        index_t nb_tet{ geomodel.mesh.cells.nb_tet() };
        index_t nb_pyramids{ geomodel.mesh.cells.nb_pyramid() };
        index_t nb_prisms{ geomodel.mesh.cells.nb_prism() };
        index_t nb_hex{ geomodel.mesh.cells.nb_hex() };
        index_t nb_poly{ geomodel.mesh.cells.nb_connector() };
        index_t nb_cells{ geomodel.mesh.cells.nb_cells() };
        Logger::out(
            "GeoModel", std::setw( 10 ), std::left, nb_cells, " cells" );
        if( nb_tet > 0 )
        {
            print_one_cell_stat( nb_tet, nb_cells, "tet" );
        }
        if( nb_pyramids > 0 )
        {
            print_one_cell_stat( nb_pyramids, nb_cells, "pyramids" );
        }
        if( nb_prisms > 0 )
        {
            print_one_cell_stat( nb_prisms, nb_cells, "prisms" );
        }
        if( nb_hex > 0 )
        {
            print_one_cell_stat( nb_hex, nb_cells, "hex" );
        }
        if( nb_poly > 0 )
        {
            print_one_cell_stat( nb_poly, nb_cells, "polyhedra" );
        }
        Logger::out( "GeoModel" );
    }

    void print_geomodel_mesh_cell_volumes( const GeoModel3D& geomodel )
    {
        double volume{ 0 };
        double tet_volume{ 0 };
        double pyramid_volume{ 0 };
        double prism_volume{ 0 };
        double hex_volume{ 0 };
        double poly_volume{ 0 };
        std::tie( volume, tet_volume, pyramid_volume, prism_volume, hex_volume,
            poly_volume ) = compute_geomodel_volumes_per_cell_type( geomodel );
        Logger::out( "GeoModel", "Model ", geomodel.name(), " has a volume of ",
            volume );
        if( tet_volume > 0 )
        {
            print_one_cell_volume_stat( tet_volume, volume, "tet" );
        }
        if( pyramid_volume > 0 )
        {
            print_one_cell_volume_stat( pyramid_volume, volume, "pyramid" );
        }
        if( prism_volume > 0 )
        {
            print_one_cell_volume_stat( prism_volume, volume, "prism" );
        }
        if( hex_volume > 0 )
        {
            print_one_cell_volume_stat( hex_volume, volume, "hex" );
        }
        if( poly_volume > 0 )
        {
            print_one_cell_volume_stat( poly_volume, volume, "polyhedron" );
        }
        Logger::out( "GeoModel" );
    }

    template < index_t DIMENSION >
    index_t find_mesh_entity_id_from_name(
        const GeoModel< DIMENSION >& geomodel,
        const MeshEntityType& gmme_type,
        const std::string& name )
    {
        index_t mesh_entity_id{ NO_ID };
        for( auto elt_i : range( geomodel.nb_mesh_entities( gmme_type ) ) )
        {
            const GeoModelMeshEntity< DIMENSION >& cur_gme =
                geomodel.mesh_entity( gmme_type, elt_i );
            if( cur_gme.name() == name )
            {
                if( mesh_entity_id != NO_ID )
                {
                    throw RINGMeshException( "GeoModel",
                        " At least two GeoModelMeshEntity have the same name "
                        "in the GeoModel: ",
                        name );
                }
                mesh_entity_id = cur_gme.index();
            }
        }
        if( mesh_entity_id == NO_ID )
        {
            throw RINGMeshException( "GeoModel", name,
                " does not match with any actual "
                "GeoModelEntity name in the GeoModel" );
        }
        return mesh_entity_id;
    }

    template < index_t DIMENSION >
    index_t find_geological_entity_id_from_name(
        const RINGMesh::GeoModel< DIMENSION >& geomodel,
        const RINGMesh::GeologicalEntityType& gmge_type,
        const std::string& name )
    {
        index_t geological_entity_id{ NO_ID };
        for( auto& cur_gme : geomodel.geol_entities( gmge_type ) )
        {
            if( cur_gme.name() == name )
            {
                if( geological_entity_id != NO_ID )
                {
                    throw RINGMeshException( "GeoModel",
                        "At least two GeoModelGeologicalEntity have the same "
                        "name in the GeoModel : ",
                        name );
                }
                geological_entity_id = cur_gme.index();
            }
        }
        if( geological_entity_id == NO_ID )
        {
            throw RINGMeshException( "GeoModel", name,
                " does not match with any actual "
                "GeoModelEntity name in the GeoModel" );
        }
        return geological_entity_id;
    }

    template void geomodel_core_api print_geomodel( const GeoModel2D& );
    template void geomodel_core_api print_geomodel_mesh_stats(
        const GeoModel2D& );
    template index_t geomodel_core_api find_mesh_entity_id_from_name(
        const GeoModel2D&, const MeshEntityType&, const std::string& );
    template index_t geomodel_core_api find_geological_entity_id_from_name(
        const RINGMesh::GeoModel2D&,
        const RINGMesh::GeologicalEntityType&,
        const std::string& );

    template void geomodel_core_api print_geomodel( const GeoModel3D& );
    template void geomodel_core_api print_geomodel_mesh_stats(
        const GeoModel3D& );
    template index_t geomodel_core_api find_mesh_entity_id_from_name(
        const GeoModel3D&, const MeshEntityType&, const std::string& );
    template index_t geomodel_core_api find_geological_entity_id_from_name(
        const RINGMesh::GeoModel3D&,
        const RINGMesh::GeologicalEntityType&,
        const std::string& );

} // namespace RINGMesh
