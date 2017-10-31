/*
 * Copyright (c) 2012-2017, Association Scientifique pour la Geologie et ses
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

namespace
{
    // The reg phys field in the GMSH format is set to 0 for each element
    static const index_t reg_phys = 0;

    static const index_t adeli_point_type = 15;
    static const index_t adeli_line_type = 1;
    static const index_t adeli_triangle_type = 2;
    static const index_t adeli_tet_type = 4;
    static const index_t adeli_cell_types[4] = { adeli_point_type,
        adeli_line_type, adeli_triangle_type, adeli_tet_type };
    class RegionAndDependentEntities
    {
    public:
        RegionAndDependentEntities(
            const GeoModel3D& geomodel, gmme_id region_gmme )
            : geomodel_( geomodel ), region_gmme_( std::move( region_gmme ) )
        {
            build();
        }
        RegionAndDependentEntities( const GeoModel3D& geomodel,
            gmme_id region_gmme,
            index_t offset_vertices,
            index_t offset_elements,
            index_t offset_corners,
            index_t offset_lines,
            index_t offset_surfaces,
            index_t offset_regions )
            : geomodel_( geomodel ),
              region_gmme_( std::move( region_gmme ) ),
              offset_vertices_( offset_vertices ),
              offset_elements_( offset_elements ),
              offset_corners_( offset_corners ),
              offset_lines_( offset_lines ),
              offset_surfaces_( offset_surfaces ),
              offset_regions_( offset_regions )
        {
            build();
        }

        index_t nb_total_elements() const
        {
            index_t nb_elements =
                geomodel_.mesh_entity( region_gmme_ ).nb_mesh_elements();
            nb_elements += nb_elements_in_entities( surfaces_ );
            nb_elements += nb_elements_in_entities( lines_ );
            nb_elements += static_cast< index_t >( corners_.size() );
            return nb_elements;
        }

        void write_entities( std::ofstream& out )
        {
            write_corners( out );
            write_entities( lines_, offset_lines_, out );
            write_entities( surfaces_, offset_surfaces_, out );
            std::vector< gmme_id > regions( 1 );
            regions[0] = region_gmme_;
            write_entities( regions, offset_regions_, out );
            offset_vertices_ +=
                geomodel_.mesh_entity( region_gmme_ ).nb_vertices();
        }

        index_t offset_vertices() const
        {
            return offset_vertices_;
        }
        index_t offset_elements() const
        {
            return offset_elements_;
        }
        index_t offset_corners() const
        {
            return offset_corners_;
        }
        index_t offset_lines() const
        {
            return offset_lines_;
        }
        index_t offset_surfaces() const
        {
            return offset_surfaces_;
        }
        index_t offset_regions() const
        {
            return offset_regions_;
        }

    private:
        void build()
        {
            get_boundary_entities(
                geomodel_.mesh_entity( region_gmme_ ), surfaces_ );
            sort_unique( surfaces_ );
            sort_unique( lines_ );
            sort_unique( corners_ );
        }
        void write_corners( std::ofstream& out )
        {
            const auto& nn =
                geomodel_.mesh_entity( region_gmme_ ).vertex_nn_search();
            for( const auto& corner_gmme : corners_ )
            {
                std::vector< index_t > element_vertices( 1 );
                element_vertices[0] =
                    nn.get_closest_neighbor(
                        geomodel_.mesh_entity( corner_gmme ).vertex( 0 ) )
                    + offset_vertices_;
                write_mesh_entity_element( adeli_point_type, element_vertices,
                    offset_corners_++, out );
            }
        }
        void write_entities( const std::vector< gmme_id >& entities,
            index_t& offset,
            std::ofstream& out )
        {
            for( const auto& entity_gmme : entities )
            {
                write_entity(
                    geomodel_.mesh_entity( entity_gmme ), offset++, out );
            }
        }

        void write_entity( const GeoModelMeshEntity3D& mesh_entity,
            index_t offset,
            std::ofstream& out )
        {
            for( auto mesh_entity_element :
                range( mesh_entity.nb_mesh_elements() ) )
            {
                std::vector< index_t > element_vertices =
                    get_element_vertices( mesh_entity, mesh_entity_element );
                write_mesh_entity_element(
                    adeli_cell_types[mesh_entity.nb_mesh_element_vertices(
                                         mesh_entity_element )
                                     - 1],
                    element_vertices, offset, out );
            }
        }

        void write_mesh_entity_element( index_t cell_descriptor,
            const std::vector< index_t >& element_vertices,
            index_t offset,
            std::ofstream& out )
        {
            out << offset_elements_++ << " " << cell_descriptor << " "
                << reg_phys << " " << offset << " " << element_vertices.size()
                << " ";
            for( auto element_vertex : element_vertices )
            {
                out << element_vertex << " ";
            }
            out << EOL;
        }

        std::vector< index_t > get_element_vertices(
            const GeoModelMeshEntity3D& mesh_entity,
            index_t element_index ) const
        {
            const auto& nn =
                geomodel_.mesh_entity( region_gmme_ ).vertex_nn_search();
            index_t nb_vertices{ mesh_entity.nb_mesh_element_vertices(
                element_index ) };
            std::vector< index_t > element_vertices( nb_vertices );
            for( auto element_local_vertex : range( nb_vertices ) )
            {
                element_vertices[element_local_vertex] =
                    offset_vertices_
                    + nn.get_closest_neighbor( mesh_entity.vertex(
                          mesh_entity.mesh_element_vertex_index(
                              { element_index, element_local_vertex } ) ) );
            }
            return element_vertices;
        }

        index_t nb_elements_in_entities(
            const std::vector< gmme_id >& entities ) const
        {
            index_t nb_elements{ 0 };
            for( const auto& entity : entities )
            {
                nb_elements +=
                    geomodel_.mesh_entity( entity ).nb_mesh_elements();
            }
            return nb_elements;
        }
        void get_boundary_entities( const GeoModelMeshEntity3D& incident,
            std::vector< gmme_id >& boundaries )
        {
            for( auto boundary_index : range( incident.nb_boundaries() ) )
            {
                const auto& boundary_gmme =
                    incident.boundary_gmme( boundary_index );
                boundaries.emplace_back( boundary_gmme );
                if( boundary_gmme.type() == surface_type_name_static() )
                {
                    get_boundary_entities(
                        geomodel_.mesh_entity( boundary_gmme ), lines_ );
                }
                else if( boundary_gmme.type() == line_type_name_static() )
                {
                    get_boundary_entities(
                        geomodel_.mesh_entity( boundary_gmme ), corners_ );
                }
                else if( boundary_gmme.type() == corner_type_name_static() )
                {
                }
                else
                {
                    ringmesh_assert_not_reached;
                }
            }
        }

    private:
        /// The exported GeoModel
        const GeoModel3D& geomodel_;

        /// The base Region
        const gmme_id region_gmme_;

        /// The Surfaces which are incident to the Region
        std::vector< gmme_id > surfaces_;

        /// The Lines which are incident to the Surface
        std::vector< gmme_id > lines_;

        /// The Corners which are incident to the Lines
        std::vector< gmme_id > corners_;

        index_t offset_vertices_{ 0 };

        index_t offset_elements_{ 0 };

        index_t offset_corners_{ 0 };

        index_t offset_lines_{ 0 };

        index_t offset_surfaces_{ 0 };

        index_t offset_regions_{ 0 };
    };

    // The index begins at 1.
    static const index_t id_offset_adeli = 1;

    /*!
     * @brief export for ADELI
     * https://sif.info-ufr.univ-montp2.fr/?q=content/adeli
     * @details This export is in fact a V1.0 of the .msh file, suitable for
     * running Finite Element Simulation with the ADELI solver.
     * The description of the output file can be found here :
     * http://gmsh.info/doc/texinfo/gmsh.html#MSH-file-format-version-1_002e0
     * First, nodes are written, then the elements. The elements are not
     * separated.
     * Corners are written (with vertex), then Lines (with edges), then Surfaces
     * (with surfaces, then Regions (with tetrahedron)
     */
    class AdeliIOHandler final : public GeoModelOutputHandler3D
    {
    public:
        void save(
            const GeoModel3D& geomodel, const std::string& filename ) final
        {
            std::ofstream out( filename.c_str() );
            out.precision( 16 );
            const RINGMesh::GeoModelMesh3D& geomodel_mesh = geomodel.mesh;
            if( geomodel_mesh.cells.nb() != geomodel_mesh.cells.nb_tet() )
            {
                throw RINGMeshException(
                    "I/O", "Adeli supports only tet meshes" );
            }
            if( geomodel.nb_regions() == 0 )
            {
                throw RINGMeshException(
                    "I/O", "Adeli can't load model with 0 regions" );
            }
            write_regions_vertices( geomodel, out );
            write_regions( geomodel, out );
            out << std::flush;
        }

    private:
        void write_regions( const GeoModel3D& geomodel, std::ofstream& out )
        {
            out << "$ELM" << EOL;
            out << count_regions_and_deps_elements( geomodel ) << EOL;
            std::vector< RegionAndDependentEntities > regions_and_deps;
            regions_and_deps.emplace_back( geomodel,
                geomodel.region( 0 ).gmme(), id_offset_adeli, id_offset_adeli,
                id_offset_adeli, id_offset_adeli, id_offset_adeli,
                id_offset_adeli );
            for( auto region_index : range( 1, geomodel.nb_regions() ) )
            {
                regions_and_deps[region_index - 1].write_entities( out );
                regions_and_deps.emplace_back( geomodel,
                    geomodel.region( region_index ).gmme(),
                    regions_and_deps[region_index - 1].offset_vertices(),
                    regions_and_deps[region_index - 1].offset_elements(),
                    regions_and_deps[region_index - 1].offset_corners(),
                    regions_and_deps[region_index - 1].offset_lines(),
                    regions_and_deps[region_index - 1].offset_surfaces(),
                    regions_and_deps[region_index - 1].offset_regions() );
            }
            regions_and_deps[geomodel.nb_regions() - 1].write_entities( out );
            out << "$ENDELM" << EOL;
        }

        void write_regions_vertices(
            const GeoModel3D& geomodel, std::ofstream& out )
        {
            out << "$NOD" << EOL;
            out << count_regions_vertices( geomodel ) << EOL;
            index_t vertex_index{ id_offset_adeli };
            for( const auto& region : geomodel.regions() )
            {
                for( auto v : range( region.nb_vertices() ) )
                {
                    out << vertex_index++ << " " << region.vertex( v ) << EOL;
                }
            }
            out << "$ENDNOD" << EOL;
        }

        index_t count_regions_vertices( const GeoModel3D& geomodel )
        {
            index_t nb_vertices = 0;
            for( const auto& region : geomodel.regions() )
            {
                nb_vertices += region.nb_vertices();
            }
            return nb_vertices;
        }

        index_t count_regions_and_deps_elements(
            const GeoModel3D& geomodel ) const
        {
            index_t regions_and_deps_elements{ 0 };
            for( const auto& region : geomodel.regions() )
            {
                RegionAndDependentEntities region_and_deps(
                    geomodel, region.gmme() );
                index_t nb_elements_in_region =
                    region_and_deps.nb_total_elements();
                regions_and_deps_elements =
                    regions_and_deps_elements + nb_elements_in_region;
            }
            return regions_and_deps_elements;
        }
    };
}
