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

#include <tinyxml2.h>
#include<ringmesh/basic/common.h>
#include<ringmesh/geomodel/core/geomodel.h>
#include<ringmesh/io/io.h>
#include<geogram/basic/attributes.h>
using namespace RINGMesh;
namespace
{
    /*!
     * In VTU/VTK, edges, polygons and cells (in the RINGMesh vocbulary)
     * are considered as... cells
     */
    index_t total_number_of_cells( const GeoModel3D& geomodel) {
        /*
        return geomodel.mesh.edges.nb() +
            geomodel.mesh.polygons.nb() +
            geomodel.mesh.cells.nb();
            */
        return geomodel.mesh.cells.nb();
    }
    /*!
     * @brief IO for the vtu file format.
     * vtu is the replacement of the legacy vtk format for
     * the unstructured grids. It is written in a XML file.
     */
    class VTUIOHandler final : public GeoModelOutputHandler3D
    {
    public:
        void save(
            const GeoModel3D& geomodel, const std::string& filename ) final
        {
            // Writer
            tinyxml2::XMLDocument output;
            
            // Header
            tinyxml2::XMLDeclaration * decl = output.NewDeclaration();
            output.InsertFirstChild(decl);
            
            // Base node
            tinyxml2::XMLElement * vtkfile = output.NewElement("VTKFile");
            vtkfile->SetAttribute("type", "UnstructuredGrid" );
            vtkfile->SetAttribute("version", "0.1" );
            vtkfile->SetAttribute("byte_order", "LittleEndian" );
            output.InsertEndChild( vtkfile);

            tinyxml2::XMLElement * ugrid = output.NewElement("UnstructuredGrid");

            // Number of vertices and cells
            tinyxml2::XMLElement * piece = output.NewElement("Piece");
            piece->SetAttribute(
                    "NumberOfPoints",
                    geomodel.mesh.vertices.nb() );
            piece->SetAttribute(
                    "NumberOfCells",
                    total_number_of_cells(geomodel) );
            ugrid->InsertEndChild( piece);

            // Vertices attribute
            tinyxml2::XMLElement * pdata = output.NewElement("PointData");
            GEO::AttributesManager & vertices_attribute_manager =
                geomodel.mesh.vertices.attribute_manager();
            GEO::vector< std::string > vertices_attribute_names;
            vertices_attribute_manager.list_attribute_names(vertices_attribute_names);
            for(auto attribute_index : range(vertices_attribute_manager.nb()) ){
                std::string attribute_name = vertices_attribute_names[attribute_index];
                if( GEO::AttributeBase< double >::is_defined(vertices_attribute_manager,
                            attribute_name) ) {
                    GEO::Attribute< double > cur_attribute( vertices_attribute_manager,
                            attribute_name);
                    index_t dimension = cur_attribute.dimension();

                    // Creation of the Data Array
                    tinyxml2::XMLElement * data_array = output.NewElement("DataArray");
                    set_data_array_attributes(data_array, attribute_name, dimension);

                    // Write the Data Array
                    for( auto vertex_index : range( geomodel.mesh.vertices.nb() ) ) {
                        tinyxml2::XMLText * to_write = output.NewText("");
                        std::string attribute_string = "\n";
                        for(auto dim : range(dimension) ) {
                            attribute_string += std::to_string(
                                    cur_attribute[vertex_index*dimension + dim]);
                            attribute_string += " ";
                        }
                        to_write->SetValue( attribute_string.c_str());
                        data_array->InsertEndChild(to_write);
                    }
                    pdata->InsertEndChild(data_array);
                }
            }
            piece->InsertEndChild( pdata );

            // Write the vertices
            tinyxml2::XMLElement * points = output.NewElement("Points");
            tinyxml2::XMLElement * vdata_array = output.NewElement("DataArray");
            set_data_array_attributes( vdata_array, "Points", 3);
            for( auto v : range(geomodel.mesh.vertices.nb() )) {
                tinyxml2::XMLText * to_write = output.NewText("");
                std::string vertex_string = "\n";
                std::ostringstream vertex_stream;
                vertex_stream << geomodel.mesh.vertices.vertex( v );
                vertex_string += vertex_stream.str();
                to_write->SetValue( vertex_string.c_str());
                vdata_array->InsertEndChild( to_write);
            }
            points->InsertEndChild( vdata_array);
            piece->InsertEndChild( points );

            // Cell attribute
            tinyxml2::XMLElement * cdata = output.NewElement("CellData");
            GEO::AttributesManager & cells_attribute_manager =
                geomodel.mesh.cells.attribute_manager();
            GEO::vector< std::string > cells_attribute_names;
            cells_attribute_manager.list_attribute_names(cells_attribute_names);
            for(auto attribute_index : range(cells_attribute_manager.nb()) ){
                std::cout << "Enter loop attribute" << std::endl;
                std::string attribute_name = cells_attribute_names[attribute_index];
                if( GEO::AttributeBase< double >::is_defined(cells_attribute_manager,
                            attribute_name) ) {
                    GEO::Attribute< double > cur_attribute( cells_attribute_manager,
                            attribute_name);
                    index_t dimension = cur_attribute.dimension();

                    // Creation of the Data Array
                    tinyxml2::XMLElement * data_array = output.NewElement("DataArray");
                    set_data_array_attributes(data_array, attribute_name, dimension);

                    // Write the Data Array
                    for( auto cell_index : range( geomodel.mesh.cells.nb() ) ) {
                        std::cout << "Enter cell loop" << std::endl;
                        tinyxml2::XMLText * to_write = output.NewText("");
                        std::string attribute_string = "\n";
                        for(auto dim : range(dimension) ) {
                            attribute_string += std::to_string(
                                    cur_attribute[cell_index*dimension + dim]);
                            attribute_string += " ";
                        }
                        to_write->SetValue( attribute_string.c_str());
                        data_array->InsertEndChild(to_write);
                    }
                    cdata->InsertEndChild(data_array);
                }
            }
            // Write the regions
            tinyxml2::XMLElement * region_array = output.NewElement("DataArray");
            set_data_array_attributes(region_array, "region", 1);
            for( auto cell_index : range( geomodel.mesh.cells.nb() ) ) {
                tinyxml2::XMLText * to_write = output.NewText("");
                std::string region_string = "\n";
                region_string += std::to_string(geomodel.mesh.cells.region( cell_index));
                to_write->SetValue( region_string.c_str());
                region_array->InsertEndChild(to_write);
            }
            cdata->InsertEndChild(region_array);

            piece->InsertEndChild( cdata );

            // Write the cells
            tinyxml2::XMLElement * cells = output.NewElement("Cells");
            tinyxml2::XMLElement * cconnectivity = output.NewElement("DataArray");
            tinyxml2::XMLElement * coffsets = output.NewElement("DataArray");
            tinyxml2::XMLElement * ctypes = output.NewElement("DataArray");
            set_data_array_attributes( cconnectivity, "connectivity", "UInt64");
            set_data_array_attributes( coffsets, "offsets", "UInt64");
            set_data_array_attributes( ctypes, "types", "UInt8");
            index_t offset = 0;
            for( auto cell_index :range( geomodel.mesh.cells.nb() ) ) {
                tinyxml2::XMLText * to_write_connectivity = output.NewText("");
                tinyxml2::XMLText * to_write_offsets = output.NewText("");
                tinyxml2::XMLText * to_write_types = output.NewText("");
                std::string connectivity_string = "\n";
                std::string offsets_string = "\n";
                std::string types_string = "\n";
                const auto& descriptor =
                    *cell_type_to_cell_descriptor_vtk[to_underlying_type(
                       geomodel. mesh.cells.type( cell_index ) )];
                for( auto local_vertex_index :
                        range( geomodel.mesh.cells.nb_vertices( cell_index ) ) ) {
                    connectivity_string +=
                        std::to_string(
                                geomodel.mesh.cells.vertex(
                                    {cell_index, local_vertex_index }));
                    connectivity_string += " ";
                }
                offset +=geomodel.mesh.cells.nb_vertices( cell_index );
                types_string += std::to_string(descriptor.entity_type);
                offsets_string +=
                    std::to_string(offset);
                to_write_connectivity->SetValue( connectivity_string.c_str());
                to_write_offsets->SetValue( offsets_string.c_str());
                to_write_types->SetValue( types_string.c_str());
                cconnectivity->InsertEndChild(to_write_connectivity);
                coffsets->InsertEndChild(to_write_offsets);
                ctypes->InsertEndChild(to_write_types);
            }
                cells->InsertEndChild(cconnectivity);
                cells->InsertEndChild(coffsets);
                cells->InsertEndChild(ctypes);
            piece->InsertEndChild( cells );

            vtkfile->InsertEndChild(ugrid);


            output.SaveFile( filename.c_str() );
        }
    private:
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

    };
}
