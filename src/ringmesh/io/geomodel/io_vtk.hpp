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

#include "common_vtk.h"
namespace
{

    /*!
     * @brief class for Legacy vtk file
     * @details this format is kind of deprecated, please prefer
     * the XML vtu file format
     */
    class VTKIOHandler final : public GeoModelOutputHandler3D
    {
    public:
        void save(
            const GeoModel3D& geomodel, const std::string& filename ) final
        {
            std::ofstream out( filename.c_str() );
            out.precision( 16 );

            write_header( out);

            write_vertices( geomodel, out);

            write_elements( geomodel, out);

            write_elements_types( geomodel, out);

            write_elements_regions( geomodel, out);

            write_vertices_attributes( geomodel, out);

            out << std::flush;
        }

    private:

        /*!
         * @brief write basic informations
         * @param[in,out] out the ofstream being completed
         */
        void write_header( std::ofstream& out) {
            out << "# vtk DataFile Version 2.0" << EOL;
            out << "Unstructured Grid" << EOL;
            out << "ASCII" << EOL;
            out << "DATASET UNSTRUCTURED_GRID" << EOL;
        }

        /*!
         * @brief write the vertices coordinates
         * @details the vertices being written are the ones on the GeoModelMesh
         * @param[in] geomodel the GeoModel being exported
         * @param[in,out] out the ofstream being completed
         */
        void write_vertices(const GeoModel3D& geomodel,
                std::ofstream& out) {
            out << "POINTS " << geomodel.mesh.vertices.nb() << " double" << EOL;
            for( auto v : range( geomodel.mesh.vertices.nb() ) )
            {
                out << geomodel.mesh.vertices.vertex( v ) << EOL;
            }
            out << EOL;
        }

        /*!
         * @brief write the vertices coordinates
         * @details the vertices being written are the ones on the GeoModelMesh
         * @param[in] geomodel the GeoModel being exported
         * @param[in,out] out the ofstream being completed
         */
        void write_vertices_attributes(const GeoModel3D& geomodel,
                std::ofstream& out) {
            auto& attribute_manager = geomodel.mesh.vertices.attribute_manager();
            out << "POINT_DATA " << geomodel.mesh.vertices.nb() << EOL;
            GEO::vector< std::string > attribute_names;
            attribute_manager.list_attribute_names(attribute_names);
            for(auto attribute_index : range(attribute_manager.nb()) ){
                std::string attribute_name = attribute_names[attribute_index];
                if( GEO::AttributeBase< double >::is_defined(attribute_manager,
                            attribute_name) ) {
                    GEO::Attribute< double > cur_attribute(attribute_manager,
                            attribute_name);
                    auto dimension = cur_attribute.dimension();
                    if( dimension == 1 ) {
                        out << "SCALARS " << attribute_name <<  " double 1" << EOL;
                        out << "LOOKUP_TABLE default" << EOL;
                    }
                    if( dimension > 1 ) {
                        out << "VECTORS " << attribute_name <<  " double" << EOL;
                    }
                    for( auto index : range( cur_attribute.size() ) ) {
                        for(auto dim : range(dimension) ) {
                            out << cur_attribute[dimension*index + dim] << SPACE;
                        }
                        out << EOL;
                    }
                }
                else {
                    Logger::warn("I/O", "Attribute ", attribute_name, "is not a double ",
                            "attribute and will not be exported");
                }
            }
            out << EOL;
        }

        /*!
         * @brief write the elements
         * @details in VTK, elements are corner vertices, edges, polygons and
         * cells in the RINGMesh point of view. All these elements are being written.
         * @param[in] geomodel the GeoModel being exported
         * @param[in,out] out the ofstream being completed
         */
        void write_elements(const GeoModel3D& geomodel, std::ofstream& out) {
            const auto& geomodel_mesh = geomodel.mesh;
            index_t total_corners = ( 4 + 1 ) * geomodel_mesh.cells.nb_tet()
                                    + ( 5 + 1 ) * geomodel_mesh.cells.nb_pyramid()
                                    + ( 6 + 1 ) * geomodel_mesh.cells.nb_prism()
                                    + ( 8 + 1 ) * geomodel_mesh.cells.nb_hex()
                                    + ( 3 + 1 ) * geomodel_mesh.polygons.nb_triangle()
                                    + ( 4 + 1 ) * geomodel_mesh.polygons.nb_quad()
                                    + ( 2 + 1 ) * geomodel_mesh.edges.nb()
                                    + geomodel.nb_corners() * 2;
            auto number_of_cells = total_number_of_elements(geomodel);

            // These 2 vectors will be filled during the process of connectivity
            // writing to avoid doing three times the loop over the entity, then the
            // elements....
            elements_types_.resize(number_of_cells) ;
            elements_regions_.resize(number_of_cells) ;

            out << "CELLS " << number_of_cells << SPACE << total_corners
                << EOL;
            index_t offset = 0;
            for( auto &corner : geomodel.corners() ) {
                write_entity_elements_connectivity(
                        corner,
                        ringmesh2vtk_vertex,
                        offset,
                        out);
            }

            for( auto &line : geomodel.lines() ) {
                write_entity_elements_connectivity(
                        line,
                        ringmesh2vtk_edges,
                        offset,
                        out );
            }

            for( auto &surface : geomodel.surfaces() ) {
                write_entity_elements_connectivity(
                        surface,
                        ringmesh2vtk_polygons,
                        offset,
                        out );
            }

            for( auto &region : geomodel.regions() ) {
                write_entity_elements_connectivity(
                        region,
                        ringmesh2vtk_cells,
                        offset,
                        out );
            }
            out << EOL;
        }

        /*!
         * @brief write the connectivity of the elements
         * @details in VTK, elements are corner vertices, edges, polygons and
         * cells in the RINGMesh point of view. All these elements are being written.
         * This function also fill the 2 vectors elements_regions_ and elements_types
         * which are written later.
         * @param[in] geomodel_mesh_entity the GeoModelMeshEntity being writtend
         * @param[in] decriptor the descriptor associated with the elements of the
         * GeoModelMeshEntity. The use is to pass from the RINGMesh indexing to the
         * VTK one
         * @param[in,out] offset this value is used to fill the types vector and the
         * regions vectors with the values.
         * @param[in,out] out the ofstream being completed
         */
        void write_entity_elements_connectivity(
                const GeoModelMeshEntity3D& geomodel_mesh_entity,
                RINGMesh2VTK* descriptor[],
                index_t& offset,
                std::ofstream& out) {
            const GeoModel3D& geomodel = geomodel_mesh_entity.geomodel();
            for( auto element : range( geomodel_mesh_entity.nb_mesh_elements() )) {
                auto nb_vertices =
                    geomodel_mesh_entity.nb_mesh_element_vertices( element );
                out << nb_vertices << SPACE;
                elements_types_[offset] = descriptor[nb_vertices]->entity_type;
                elements_regions_[offset++] = geomodel_mesh_entity.index();
                for( auto local_vertex_index : range( nb_vertices ) ){
                    auto vertex_index_vtk =
                        geomodel.mesh.vertices.geomodel_vertex_id(
                                geomodel_mesh_entity.gmme(),
                                geomodel_mesh_entity.mesh_element_vertex_index(
                                    {element, descriptor[nb_vertices]
                                    ->vertices[local_vertex_index]}));
                    out <<vertex_index_vtk << SPACE;
                }
                out << EOL;
            }
        }

        /*!
         * @brief write the VTK types of the elements
         * @param[in] geomodel the GeoModel being written
         * @param[in,out] out the ofstream being completed
         */
        void write_elements_types(
                const GeoModel3D& geomodel,
                std::ofstream& out
                ) {
            auto number_of_cells = total_number_of_elements(geomodel);
            out << "CELL_TYPES " << number_of_cells << EOL;
            for( auto element_index : range( number_of_cells ) ) {
                out << elements_types_[element_index] << EOL;
            }
            out << EOL;
        }

        /*!
         * @brief write the regions of the elements
         * @param[in] geomodel the GeoModel being written
         * @param[in,out] out the ofstream being completed
         */
        void write_elements_regions(
                const GeoModel3D& geomodel,
                std::ofstream& out
                ) {
            auto number_of_cells = total_number_of_elements(geomodel);
            out << "CELL_DATA " << number_of_cells << EOL;
            out << "SCALARS region int 1" << EOL;
            out << "LOOKUP_TABLE default" << EOL;
            for( auto element_index : range( number_of_cells ) ) {
                out << elements_regions_[element_index] << EOL;
            }
            out << EOL;
        }
    private:
        /// This vector contains the VTK cell types of the corner vertices,
        /// edges, polygons and cells
        std::vector< index_t > elements_types_;
        
        /// This vector contains the VTK cell regions of the corner vertices,
        /// edges, polygons and cells (i.e. the corner index, the line index,
        /// the surface index or the region index in
        std::vector< index_t > elements_regions_;
    };
}
