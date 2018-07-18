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

#include <ringmesh/basic/common.h>

#include <geogram/basic/command_line.h>
#include <geogram/basic/stopwatch.h>
#include <geogram/mesh/mesh_io.h>

#include <ringmesh/basic/command_line.h>
#include <ringmesh/geomodel/builder/geomodel_builder.h>
#include <ringmesh/geomodel/core/geomodel.h>
#include <ringmesh/geomodel/core/geomodel_api.h>
#include <ringmesh/geomodel/core/geomodel_mesh_entity.h>
#include <ringmesh/geomodel/tools/geomodel_validity.h>
#include <ringmesh/io/io.h>
#include <ringmesh/mesh/mesh_index.h>

#include <metis.h>

/*!
 * @author Antoine Mazuyer
 */

namespace
{
    using namespace RINGMesh;

    template < index_t DIMENSION >
    index_t total_number_of_elements( const GeoModel< DIMENSION >& geomodel );
    template < index_t DIMENSION >
    index_t total_number_of_elements_base(
        const GeoModel< DIMENSION >& geomodel )
    {
        return geomodel.nb_corners() + geomodel.mesh.edges.nb()
               + geomodel.mesh.polygons.nb();
    }
    template <>
    index_t total_number_of_elements( const GeoModel2D& geomodel )
    {
        return total_number_of_elements_base( geomodel );
    }
    template <>
    index_t total_number_of_elements( const GeoModel3D& geomodel )
    {
        // return total_number_of_elements_base( geomodel ) +
        // geomodel.mesh.cells.nb();
        return geomodel.mesh.cells.nb();
    }

    template < index_t DIMENSION >
    index_t total_number_of_corners( const GeoModel< DIMENSION >& geomodel );
    template < index_t DIMENSION >
    index_t total_number_of_corners_base(
        const GeoModel< DIMENSION >& geomodel )
    {
        return 3 * geomodel.mesh.polygons.nb_triangle()
               + 4 * geomodel.mesh.polygons.nb_quad()
               + 2 * geomodel.mesh.edges.nb() + geomodel.nb_corners();
    }
    template <>
    index_t total_number_of_corners( const GeoModel2D& geomodel )
    {
        return total_number_of_corners_base( geomodel );
    }
    template <>
    index_t total_number_of_corners( const GeoModel3D& geomodel )
    {
        /*return total_number_of_corners_base(geomodel) +  */
        return 4 * geomodel.mesh.cells.nb_tet()
               + 5 * geomodel.mesh.cells.nb_pyramid()
               + 6 * geomodel.mesh.cells.nb_prism()
               + 8 * geomodel.mesh.cells.nb_hex();
    }

    template < index_t DIMENSION >
    void write_entity_elements_connectivity(
        const GeoModelMeshEntity< DIMENSION >& geomodel_mesh_entity,
        idx_t* eptr,
        idx_t* eind,
        index_t& offset_eptr,
        index_t& offset_eind )
    {
        const GeoModel< DIMENSION >& geomodel = geomodel_mesh_entity.geomodel();
        for( auto element : range( geomodel_mesh_entity.nb_mesh_elements() ) )
        {
            offset_eptr++;
            auto nb_vertices =
                geomodel_mesh_entity.nb_mesh_element_vertices( element );
            eptr[offset_eptr] = eptr[offset_eptr - 1] + nb_vertices;
            for( auto local_vertex_index : range( nb_vertices ) )
            {
                auto vertex_index = geomodel.mesh.vertices.geomodel_vertex_id(
                    geomodel_mesh_entity.gmme(),
                    geomodel_mesh_entity.mesh_element_vertex_index(
                        { element, local_vertex_index } ) );
                eind[offset_eind++] = vertex_index;
            }
        }
    }

    template < index_t DIMENSION >
    void create_metis_table(
        const GeoModel< DIMENSION >& geomodel, idx_t* eptr, idx_t* eind );

    template <>
    void create_metis_table(
        const GeoModel2D& geomodel, idx_t* eptr, idx_t* eind )
    {
        index_t offset_eptr = 0;
        index_t offset_eind = 0;
        for( auto& corner : geomodel.corners() )
        {
            write_entity_elements_connectivity(
                corner, eptr, eind, offset_eptr, offset_eind );
        }

        for( auto& line : geomodel.lines() )
        {
            write_entity_elements_connectivity(
                line, eptr, eind, offset_eptr, offset_eind );
        }

        for( auto& surface : geomodel.surfaces() )
        {
            write_entity_elements_connectivity(
                surface, eptr, eind, offset_eptr, offset_eind );
        }
    }
    template <>
    void create_metis_table(
        const GeoModel3D& geomodel, idx_t* eptr, idx_t* eind )
    {
        index_t offset_eptr = 0;
        index_t offset_eind = 0;
        /*
        for( auto &corner : geomodel.corners() ) {
            write_entity_elements_connectivity(
                    corner,
                    eptr,
                    eind,
                    offset_eptr,
                    offset_eind
                    );
        }

        for( auto &line : geomodel.lines() ) {
            write_entity_elements_connectivity(
                    line,
                    eptr,
                    eind,
                    offset_eptr,
                    offset_eind);
        }

        for( auto &surface : geomodel.surfaces() ) {
            write_entity_elements_connectivity(
                    surface,
                    eptr,
                    eind,
                    offset_eptr,
                    offset_eind);
        }
        */

        for( auto& region : geomodel.regions() )
        {
            write_entity_elements_connectivity(
                region, eptr, eind, offset_eptr, offset_eind );
        }
    }
    template < index_t DIMENSION >
    void write_attributes(
        const GeoModel< DIMENSION >& geomodel, idx_t* epart, idx_t* npart );

    template < index_t DIMENSION >
    void write_attributes_base(
        const GeoModel< DIMENSION >& geomodel, idx_t* epart, idx_t* npart )
    {
        GEO::AttributesManager& vertex_attribute_manager =
            geomodel.mesh.vertices.attribute_manager();
        GEO::Attribute< index_t > partition_vertex(
            vertex_attribute_manager, "partition" );
        for( auto vertex_index : range( geomodel.mesh.vertices.nb() ) )
        {
            partition_vertex[vertex_index] = npart[vertex_index];
        }
    }
    template <>
    void write_attributes(
        const GeoModel2D& geomodel, idx_t* epart, idx_t* npart )
    {
    }
    template <>
    void write_attributes(
        const GeoModel3D& geomodel, idx_t* epart, idx_t* npart )
    {
        DEBUG( "toto" );
        write_attributes_base( geomodel, epart, npart );
        GEO::AttributesManager& cell_attribute_manager =
            geomodel.mesh.cells.attribute_manager();
        GEO::Attribute< index_t > partition_cells(
            cell_attribute_manager, "partition" );
        for( auto cell_index : range( geomodel.mesh.cells.nb() ) )
        {
            partition_cells[cell_index] = epart[cell_index];
        }
        geomodel.mesh.transfer_attributes_from_gmm_to_gm_regions();
    }

    template < index_t DIMENSION >
    void write_attributes( const GeoModel< DIMENSION >& geomodel, idx_t* part );

    template <>
    void write_attributes( const GeoModel2D& geomodel, idx_t* part )
    {
    }
    template <>
    void write_attributes( const GeoModel3D& geomodel, idx_t* part )
    {
        DEBUG( "toto" );
        GEO::AttributesManager& cell_attribute_manager =
            geomodel.mesh.cells.attribute_manager();
        GEO::Attribute< index_t > partition_cells(
            cell_attribute_manager, "partition" );
        for( auto cell_index : range( geomodel.mesh.cells.nb() ) )
        {
            partition_cells[cell_index] = part[cell_index];
        }
        geomodel.mesh.transfer_attributes_from_gmm_to_gm_regions();
    }
    template < index_t DIMENSION >
    void build_graph(
        const GeoModel< DIMENSION >& geomodel, idx_t* xadj, idx_t* adjncy );

    template <>
    void build_graph( const GeoModel2D& geomodel, idx_t* xadj, idx_t* adjncy )
    {
    }

    template <>
    void build_graph( const GeoModel3D& geomodel, idx_t* xadj, idx_t* adjncy )
    {
        index_t offset = 0;
        xadj[0] = 0;
        for( auto cell_index : range( geomodel.mesh.cells.nb() ) )
        {
            auto nb_facets = geomodel.mesh.cells.nb_facets( cell_index );
            index_t nb_adjency = 0;
            for( auto facet_index : range( nb_facets ) )
            {
                index_t adjacent_cell_index =
                    geomodel.mesh.cells.adjacent( cell_index, facet_index );
                if( adjacent_cell_index != -1 )
                {
                    adjncy[offset++] =
                        static_cast< idx_t >( adjacent_cell_index );
                    nb_adjency++;
                }
                else
                {
                    DEBUG( "HAHAHAHAHAHA" );
                }
            }
            xadj[cell_index + 1] =
                xadj[cell_index] + static_cast< idx_t >( nb_adjency );
            DEBUG( cell_index );
            DEBUG( xadj[cell_index + 1] );
        }
    }

    template < index_t DIMENSION >
    index_t nb_co( const GeoModel< DIMENSION >& geomodel );

    template <>
    index_t nb_co( const GeoModel2D& geomodel )
    {
        return 0;
    }

    template <>
    index_t nb_co( const GeoModel3D& geomodel )
    {
        index_t nb_cos = 0;
        for( auto cell_index : range( geomodel.mesh.cells.nb() ) )
        {
            auto nb_facets = geomodel.mesh.cells.nb_facets( cell_index );
            for( auto facet_index : range( nb_facets ) )
            {
                index_t adjacent_cell_index =
                    geomodel.mesh.cells.adjacent( cell_index, facet_index );
                if( adjacent_cell_index != -1 )
                {
                    nb_cos++;
                }
            }
        }
        return nb_cos;
    }
    template < index_t DIMENSION >
    void partition_geomodel( const std::string& geomodel_in_name )
    {
        GeoModel< DIMENSION > geomodel;
        geomodel_load( geomodel, geomodel_in_name );
        index_t number_of_element = total_number_of_elements( geomodel );
        index_t total_corners = total_number_of_corners( geomodel );
        // eptr vector (see metis documentation)
        idx_t* eptr = new idx_t[static_cast< int >( number_of_element + 1 )];
        eptr[0] = 0;

        // eind vector (see metis documentation)
        idx_t* eind = new idx_t[static_cast< int >( total_corners )];
        DEBUG( number_of_element + 1 );
        DEBUG( total_corners );
        /*
        create_metis_table( geomodel, eptr, eind );
        idx_t * ne = new idx_t(static_cast< idx_t > (number_of_element));
        idx_t * nn = new idx_t(static_cast< idx_t >
        (geomodel.mesh.vertices.nb())); idx_t * nparts = new idx_t(4); idx_t *
        epart = new idx_t[static_cast< idx_t > (number_of_element)]; idx_t *
        npart = new idx_t[static_cast< idx_t > (geomodel.mesh.vertices.nb())] ;
        idx_t * ncommon = new idx_t(3) ;
        idx_t * objval = new idx_t(0);

        DEBUG("c'est passeRRRe");
        METIS_PartMeshNodal( ne, nn, eptr, eind, nullptr, nullptr,
                nparts,nullptr,nullptr,objval,epart,npart);
        */
        DEBUG( "HEIN" );
        idx_t* nvtxs = new idx_t( static_cast< idx_t >( number_of_element ) );
        idx_t* ncon = new idx_t( 1 );
        index_t nb_cos = nb_co( geomodel );
        idx_t* xadj = new idx_t[static_cast< int >( number_of_element + 1 )];
        idx_t* adjncy = new idx_t[static_cast< int >( 2 * nb_cos )];
        idx_t* objval = new idx_t( 0 );
        idx_t* nparts = new idx_t( 10 );
        idx_t* part = new idx_t[static_cast< int >( number_of_element )];
        idx_t options[40];
        METIS_SetDefaultOptions( options );
        options[METIS_OPTION_CONTIG] = 1; // set METIS_OPTION_CONTIG
        build_graph( geomodel, xadj, adjncy );
        DEBUG( "GFINI" );
        METIS_PartGraphKway( nvtxs, ncon, xadj, adjncy, nullptr, nullptr,
            nullptr, nparts, nullptr, nullptr, options, objval, part );

        write_attributes( geomodel, part );
        DEBUG( "HEINNN" );
        std::string geomodel_out_name = GEO::CmdLine::get_arg( "out:geomodel" );
        if( geomodel_out_name.empty() )
        {
            throw RINGMeshException(
                "I/O", "Give the parameter out:geomodel to save the geomodel" );
        }
        geomodel_save( geomodel, geomodel_out_name );
    }

    void show_usage_example()
    {
        Logger::div( "Example" );
        Logger::out( "",
            "ringmesh-partition in:geomodel=path/to/input/geomodel.ext ",
            "out:geomodel=path/to/output/geomodel.vtu" );
    }
}

int main( int argc, char** argv )
{
    using namespace RINGMesh;

    try
    {
        print_header_information();
        Logger::div( "RINGMesh-Convert" );
        Logger::out( "", "Welcome to RINGMesh-Convert !" );

        CmdLine::import_arg_group( "in" );
        CmdLine::import_arg_group( "out" );
        if( argc == 1 )
        {
            GEO::CmdLine::show_usage();
            show_usage_example();
            return 0;
        }

        std::vector< std::string > filenames;
        if( !GEO::CmdLine::parse( argc, argv, filenames ) )
        {
            show_usage_example();
            return 1;
        }

        GEO::Stopwatch total( "Total time" );

        std::string geomodel_in_name = GEO::CmdLine::get_arg( "in:geomodel" );
        if( geomodel_in_name.empty() )
        {
            throw RINGMeshException(
                "I/O", "Give at least a filename in in:geomodel" );
        }

        index_t dimension = find_geomodel_dimension( geomodel_in_name );
        if( dimension == 2 )
        {
            partition_geomodel< 2 >( geomodel_in_name );
        }
        else if( dimension == 3 )
        {
            partition_geomodel< 3 >( geomodel_in_name );
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
    return 0;
}
