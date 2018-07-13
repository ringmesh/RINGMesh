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
#include <tinyxml2.h>

namespace
{
    /*!
     * @brief IO for the vtu file format.
     * vtu is the replacement of the legacy vtk format for
     * the unstructured grids. It is written in a XML file.
     */
    class VTUIOHandler final : public GeoModelOutputHandler3D
    {
    public:
        VTUIOHandler() :
            output_(){
                // All nodes of the XML file initialization
                init_all_nodes();
                set_nodes_relations();
            }
        void save(
            const GeoModel3D& geomodel, const std::string& filename ) final
        {
            write_header();
            vtkfile_ = output_.NewElement("VTKFile");

            // Not a vtup (p for parallel format) so one piece
            write_piece( geomodel );

            write_vertices_attributes( geomodel );

            write_vertices( geomodel );

            write_elements_and_regions( geomodel );

            write_xml_file( filename);
        }
    private:
        void init_all_nodes() {
            vtkfile_ = output_.NewElement("VTKFile");
            ugrid_ =output_.NewElement("UnstructuredGrid");
            piece_ =output_.NewElement("Piece");
            pdata_ =output_.NewElement("PointData");
            points_ =output_.NewElement("Points");
            cdata_ =output_.NewElement("CellData");
            cells_ =output_.NewElement("Cells");
        }
        
        void set_nodes_relations() {
            output_.InsertFirstChild(vtkfile_);
            vtkfile_->InsertEndChild(ugrid_);
            ugrid_->InsertEndChild(piece_);
            piece_->InsertEndChild(pdata_);
            piece_->InsertEndChild(points_);
            piece_->InsertEndChild(cdata_);
            piece_->InsertEndChild(cells_);
        }

        /*!
         * @brief set basic informations of the XML file
         */
        void write_header() {
            auto * decl = output_.NewDeclaration();
            output_.InsertFirstChild(decl);
            vtkfile_->SetAttribute("type", "UnstructuredGrid" );
            vtkfile_->SetAttribute("version", "0.1" );
            vtkfile_->SetAttribute("byte_order", "LittleEndian" );
        }

        /*!
         * @brief Write basic information of the mesh
         * @details Write vertices and cells number
         * @warning The number of cells is equal, in a RINGMesh
         * point of view to the number of
         * edges + the number of polygons + the number of cells
         * @param[in] geomodel the GeoModel being exported
         */
        void write_piece(const GeoModel3D& geomodel ) {
            piece_->SetAttribute(
                    "NumberOfPoints",
                    geomodel.mesh.vertices.nb() );
            piece_->SetAttribute(
                    "NumberOfCells",
                    total_number_of_elements(geomodel) );
        }

        /*!
         * @brief Write the vertices attributes
         * @details the attributes being written are the on present on the
         * GeoModelMesh
         * @param[in] geomodel the GeoModel being exported
         */
        void write_vertices_attributes( const GeoModel3D& geomodel ) {
            auto & vertices_attribute_manager =
                geomodel.mesh.vertices.attribute_manager();
            write_attributes( vertices_attribute_manager );
        }

        /*!
         * @brief Write all the attributes in a attribute manager
         * @details This method is used to write the attributes present
         * ont the GeoModelMesh for edges, corners and cells.
         */
        void write_attributes( GEO::AttributesManager& attribute_manager) {
            GEO::vector< std::string > attribute_names;
            attribute_manager.list_attribute_names(attribute_names);
            for(auto attribute_index : range(attribute_manager.nb()) ){
                std::string attribute_name = attribute_names[attribute_index];
                if( GEO::AttributeBase< double >::is_defined(attribute_manager,
                            attribute_name) ) {
                    GEO::Attribute< double > cur_attribute(attribute_manager,
                            attribute_name);
                    auto dimension = cur_attribute.dimension();

                    // Creation of the Data Array
                    auto * data_array = output_.NewElement("DataArray");
                    set_data_array_attributes(data_array, attribute_name, dimension);

                    // Write the Data Array
                    for( auto index : range( cur_attribute.size() ) ) {
                        auto * to_write = output_.NewText("");
                        std::string attribute_string = "\n";
                        for(auto dim : range(dimension) ) {
                            attribute_string += std::to_string(
                                    cur_attribute[index*dimension + dim]);
                            attribute_string += " ";
                        }
                        to_write->SetValue( attribute_string.c_str());
                        data_array->InsertEndChild(to_write);
                    }
                    pdata_->InsertEndChild(data_array);
                }
                else {
                    Logger::warn("I/O", "Attribute ", attribute_name, "is not a double ",
                            "attribute and will not be exported");
                }
            }
        }

        /*!
         * @brief Write the vertices coordinates
         * @details the coordinates are stored in a data array
         * @param[in] geomodel the GeoModel being exported
         */
        void write_vertices( const GeoModel3D& geomodel ){
            auto * vdata_array = output_.NewElement("DataArray");
            set_data_array_attributes( vdata_array, "Points", 3);
            for( auto v : range(geomodel.mesh.vertices.nb() )) {
                auto * to_write = output_.NewText("");
                std::string vertex_string = "\n";
                std::ostringstream vertex_stream;
                vertex_stream << geomodel.mesh.vertices.vertex( v );
                vertex_string += vertex_stream.str();
                to_write->SetValue( vertex_string.c_str());
                vdata_array->InsertEndChild( to_write);
            }
            points_->InsertEndChild( vdata_array);
        }

        /*!
         * @brief Write the elements information
         * @details We write corner vertices, edges, polygons and cells
         * in the RINGMesh point of view. In vtk, all these elements are "elements"
         * An attribute region is also written, corresponding to the corner, the line
         * the surface or the region which belong to the element.
         * @param[in] geomodel the GeoModel being exported
         */
        void write_elements_and_regions( const GeoModel3D& geomodel ) {
            auto * cregion = output_.NewElement("DataArray");
            auto * cconnectivity = output_.NewElement("DataArray");
            auto * coffsets = output_.NewElement("DataArray");
            auto * ctypes = output_.NewElement("DataArray");
            set_data_array_attributes( cregion, "region", "UInt64");
            set_data_array_attributes( cconnectivity, "connectivity", "UInt64");
            set_data_array_attributes( coffsets, "offsets", "UInt64");
            set_data_array_attributes( ctypes, "types", "UInt8");
            index_t offset = 0;
            for( auto &corner : geomodel.corners() ) {
                write_entity_elements(
                        corner,
                        ringmesh2vtk_vertex,
                        cregion, cconnectivity, coffsets, ctypes,
                        offset );
            }

            for( auto &line : geomodel.lines() ) {
                write_entity_elements(
                        line,
                        ringmesh2vtk_edges,
                        cregion, cconnectivity, coffsets, ctypes,
                        offset );
            }

            for( auto &surface : geomodel.surfaces() ) {
                write_entity_elements(
                        surface,
                        ringmesh2vtk_polygons,
                        cregion, cconnectivity, coffsets, ctypes,
                        offset );
            }

            for( auto &region : geomodel.regions() ) {
                write_entity_elements(
                        region,
                        ringmesh2vtk_cells,
                        cregion, cconnectivity, coffsets, ctypes,
                        offset );
            }
            cdata_->InsertEndChild( cregion );
            cells_->InsertEndChild( cconnectivity );
            cells_->InsertEndChild( coffsets );
            cells_->InsertEndChild( ctypes );
        }


        /*!
         * @brief write all the elements of a mesh entity
         * @param[in] geomodel_mesh_entity the entity being written,
         * @param[in] descriptor one structure which contains information
         * on the elements of the entity (i.e. the vtk type and the vtk connectivity
         * @param[in] region_array one element of the XML file containing the region index
         * @param[in] cconnectivity one element of the XML file containing the element
         * connectivity
         * @param[in] coffsets one element of the XML file containing the offsets,
         * i.e. a value which is incremented by the number of vertices of the elements
         * alreay written
         * @param[in] ctypes one element of the XML file containing the types of the elements
         * @param[in,out] offset the offset (i.e. the cumulated number of vertices
         * of elements already written)
         */
        void write_entity_elements( const GeoModelMeshEntity3D& geomodel_mesh_entity,
                RINGMesh2VTK* descriptor[],
                tinyxml2::XMLElement * region_array,
                tinyxml2::XMLElement * cconnectivity,
                tinyxml2::XMLElement * coffsets,
                tinyxml2::XMLElement * ctypes,
                index_t& offset) {
            const GeoModel3D& geomodel = geomodel_mesh_entity.geomodel();
            for( auto element : range( geomodel_mesh_entity.nb_mesh_elements() )) {
                auto * to_write_region = output_.NewText("");
                auto * to_write_connectivity = output_.NewText("");
                auto * to_write_offsets = output_.NewText("");
                auto * to_write_types = output_.NewText("");
                std::string region_string = "\n";
                std::string connectivity_string = "\n";
                std::string offsets_string = "\n";
                std::string types_string = "\n";
                auto nb_vertices =
                    geomodel_mesh_entity.nb_mesh_element_vertices( element );
                for( auto local_vertex_index : range( nb_vertices ) ){
                    auto vertex_index_vtk =
                        geomodel.mesh.vertices.geomodel_vertex_id(
                                geomodel_mesh_entity.gmme(),
                                geomodel_mesh_entity.mesh_element_vertex_index(
                                    {element, descriptor[nb_vertices]
                                    ->vertices[local_vertex_index]}));
                    connectivity_string += std::to_string( vertex_index_vtk );
                    connectivity_string += " ";
                }
                region_string += std::to_string(geomodel_mesh_entity.index());
                types_string += std::to_string(descriptor[nb_vertices]->entity_type);
                offset += nb_vertices;
                offsets_string +=
                    std::to_string(offset);
                to_write_connectivity->SetValue( connectivity_string.c_str());
                to_write_offsets->SetValue( offsets_string.c_str());
                to_write_types->SetValue( types_string.c_str());
                to_write_region->SetValue( region_string.c_str());
                region_array->InsertEndChild(to_write_region);
                cconnectivity->InsertEndChild(to_write_connectivity);
                coffsets->InsertEndChild(to_write_offsets);
                ctypes->InsertEndChild(to_write_types);
            }
        }

        void set_data_array_attributes( 
                tinyxml2::XMLElement * xml_element,
                const std::string& name,
                const std::string& type = "Float64",
                const std::string& format = "ascii") {
            xml_element->SetAttribute(
                    "type", type.c_str());
            xml_element->SetAttribute(
                    "Name", name.c_str());
            xml_element->SetAttribute(
                    "format", format.c_str());
        }

        void set_data_array_attributes( 
                tinyxml2::XMLElement * xml_element,
                const std::string& name,
                index_t dimension,
                const std::string& type = "Float64",
                const std::string& format = "ascii") {
            set_data_array_attributes(xml_element, name, type, format);
            xml_element->SetAttribute(
                    "NumberOfComponents",
                    dimension);
        }

        void write_xml_file(const std::string& filename) {
            output_.SaveFile( filename.c_str() );
        }
    private:
        /// main file of tinyxml2
        tinyxml2::XMLDocument output_;

        /// Pointer to the main node of the XML file
        tinyxml2::XMLElement * vtkfile_;

        /// Pointer to the unstructured grid node
        tinyxml2::XMLElement * ugrid_;
        
        /// Pointer to the only piece of this mesh (see vtup file for more than one pieces)
        tinyxml2::XMLElement* piece_;

        /// Pointer to the node which contains the point data
        tinyxml2::XMLElement* pdata_;

        /// Pointer ot the node which contains the point coordinates
        tinyxml2::XMLElement* points_;

        /// Pointer to the node which contains the cell data
        /// Warning : Here cells are edges, polygons and cells (in RINGMesh sense)
        tinyxml2::XMLElement* cdata_;

        /// Pointer ot the node which contains the point coordinates
        /// Warning : Here cells are edges, polygons and cells (in RINGMesh sense)
        tinyxml2::XMLElement* cells_;
    };
}

