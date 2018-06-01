/*
 * Copyright (c) 2018, Association Scientifique pour la Geologie et ses
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

#include <geogram/basic/attributes.h>

#include <ringmesh/basic/geometry.h>
#include <ringmesh/geomodel/core/geomodel.h>
#include <ringmesh/geomodel/core/geomodel_api.h>

#include <ringmesh/geomodel/core/geomodel_mesh_entity.h>
#include <ringmesh/io/geomodel_builder_resqml.h>

#include <ringmesh/mesh/mesh_index.h>

#include <ringmesh/mesh/surface_mesh.h>
#include <ringmesh/mesh/volume_mesh.h>

#include <fesapi/common/HdfProxy.h>
#include <fesapi/common/EpcDocument.h>
#include <fesapi/resqml2_0_1/LocalDepth3dCrs.h>
#include <fesapi/resqml2_0_1/LocalTime3dCrs.h>
#include <fesapi/resqml2_0_1/UnstructuredGridRepresentation.h>

/*!
 * @brief Implementation of the class to build GeoModel from input
 * RESQML2 .epc file
 * @author Wan-Chiu Li
 */

namespace RINGMesh{

using namespace RESQML2_0_1_NS;


class GeoModelBuilderRESQMLImpl{
public:
    explicit GeoModelBuilderRESQMLImpl(GeoModelBuilderRESQML& builder);
    ~GeoModelBuilderRESQMLImpl() = default;

    bool load_file();
    bool build_fake_geomodel();
    void output_info(const COMMON_NS::EpcDocument& pck);

private:
    GeoModelBuilderRESQML& builder_;
};

GeoModelBuilderRESQMLImpl::GeoModelBuilderRESQMLImpl(
    GeoModelBuilderRESQML& builder
) : builder_(builder){
    build_fake_geomodel();
}

void GeoModelBuilderRESQMLImpl::output_info(const COMMON_NS::EpcDocument& pck){
    std::cout << "EpcDocument name " << pck.getName() << " in " << (pck.getStorageDirectory().empty() ? "working directory." : pck.getStorageDirectory()) << std::endl;

    unsigned int hdfProxyCount = pck.getHdfProxyCount();
    std::cout << "There are " << pck.getHdfProxyCount() << " hdf files associated to this epc document." << std::endl;
    for (unsigned int hdfProxyIndex = 0; hdfProxyIndex < hdfProxyCount; ++hdfProxyIndex) {
        std::cout << "Hdf file relative path : " << pck.getHdfProxy(hdfProxyIndex)->getRelativePath() << std::endl;
    }
    for (size_t warningIndex = 0; warningIndex < pck.getWarnings().size(); ++warningIndex) {
        std::cout << "Warning #" << warningIndex << " : " << pck.getWarnings()[warningIndex] << std::endl;
    }
    std::vector<UnstructuredGridRepresentation*> unstructuredGridRepSet = pck.getUnstructuredGridRepresentationSet();

    std::cout << std::endl << "UNSTRUCTURED GRID REP: " << unstructuredGridRepSet.size()<<std::endl;
}

bool GeoModelBuilderRESQMLImpl::load_file(){

    COMMON_NS::EpcDocument pck(builder_.filename(), COMMON_NS::EpcDocument::READ_ONLY);

    std::string resqmlResult = pck.deserialize();
    if (!resqmlResult.empty()) {
        std::cerr << resqmlResult << std::endl;
        std::cout << "Press enter to continue..." << std::endl;
        std::cin.get();
    }

    output_info(pck);

    std::vector<UnstructuredGridRepresentation*> unstructuredGridRepSet = pck.getUnstructuredGridRepresentationSet();

    for (size_t i = 0; i < unstructuredGridRepSet.size(); ++i)
    {

        if (!unstructuredGridRepSet[i]->isPartial() && unstructuredGridRepSet[i]->hasGeometry())
        {
            unstructuredGridRepSet[i]->loadGeometry();

            std::unique_ptr<double[]> gridPoints( new double[unstructuredGridRepSet[i]->getXyzPointCountOfPatch(0) * 3]);
            unstructuredGridRepSet[i]->getXyzPointsOfAllPatchesInGlobalCrs(&gridPoints[0]);

            const index_t output_region{ 0 };  // always 0 since there is only 1
            const gmme_id cur_region( region_type_name_static(), output_region );

            const ULONG64 nb_vertices = unstructuredGridRepSet[i]->getXyzPointCountOfPatch(0);
            builder_.geometry.create_mesh_entity_vertices( cur_region, nb_vertices );

            for(ULONG64 i=0; i<nb_vertices; ++i){
                 const vec3 vertex(gridPoints[i*3], gridPoints[i*3+1], gridPoints[i*3+2]);
                 builder_.geometry.set_mesh_entity_vertex(cur_region, i,vertex, false  );
            }

            const ULONG64 nb_cells = unstructuredGridRepSet[i]->getCellCount();
            builder_.geometry.create_region_cells( cur_region.index(), CellType::TETRAHEDRON, nb_cells );

            std::unique_ptr<ULONG64[]> faceCountOfCells( new ULONG64[nb_cells] );
            unstructuredGridRepSet[i]->getCumulativeFaceCountPerCell(&faceCountOfCells[0]);

            ULONG64 faceCount = faceCountOfCells[nb_cells - 1];

            std::unique_ptr<ULONG64[]> nodeCountOfFaces( new ULONG64[faceCount] );
            unstructuredGridRepSet[i]->getCumulativeNodeCountPerFace(&nodeCountOfFaces[0]);

            for(ULONG64 cell=0; cell<nb_cells; ++cell){
                ULONG64 end_face = faceCountOfCells[cell];
                ULONG64 start_face = (cell==0)? 0 : faceCountOfCells[cell-1];

                std::vector< index_t > vertices = {
                    unstructuredGridRepSet[i]->getNodeIndicesOfFaceOfCell(cell, 0)[0],
                    unstructuredGridRepSet[i]->getNodeIndicesOfFaceOfCell(cell, 0)[1],
                    unstructuredGridRepSet[i]->getNodeIndicesOfFaceOfCell(cell, 0)[2],
                    0};

                bool found = false;
                for(ULONG64 f=1; f< (end_face-start_face); ++f){
                    ULONG64 nb_nodes = unstructuredGridRepSet[i]->getNodeCountOfFaceOfCell(cell, f);

                    for(ULONG64 node = 0; node<nb_nodes; ++node){
                        const ULONG64 node_index = unstructuredGridRepSet[i]->getNodeIndicesOfFaceOfCell(cell, f)[node];
                        if( node_index != vertices[0] && node_index != vertices[1] && node_index != vertices[2]){
                            vertices[3] = node_index ;
                            found=true;
                            break;
                        }

                        if(found){
                            break;
                        }
                    }
                }
                ringmesh_assert(found);

                builder_.geometry.set_region_element_geometry( cur_region.index(), cell, vertices );
            }
            unstructuredGridRepSet[i]->unloadGeometry();

        }
    }
    return true;
}

bool GeoModelBuilderRESQMLImpl::build_fake_geomodel(){
    //#############################
    // Declaration of the Entities#
    //#############################

    // For the next section, read the documentation to understand
    // the concept of Geological Entity and Mesh Entities
    // Let's to a sum up of the GeoModel we try to build:
    // For the Geological Entities (handle by the class
    // GeoModelGeologicalEntity):
    // 16 Contacts
    index_t nb_contacts = 6;
    // 1 horizons + 6 boundaries = 7 Interfaces
    index_t nb_interfaces = 4;
    // 2 Layers
    index_t nb_layers = 1;

    // For the Meshed Entities, (handle by the class GeoModelMeshEntity)
    // 12 Corners
    index_t nb_corners = 4;
    // 20 Lines
    index_t nb_lines = 6;
    // 11 Surfaces
    index_t nb_surfaces = 4;
    // 2  Regions
    index_t nb_regions = 1;

    // We first create the GeoModelGeoglogicalEntity
    // Create the contacts
    for( index_t contact = 0; contact < nb_contacts; contact++ )
    {
        builder_.geology.create_geological_entity(
            Contact3D::type_name_static() );
        // the static method type_name_static() is available for each
        // GeoModelEntity. It returns an EntityType which is a string
        // corresponding to the Type of the entity.
    }

    // Create the Interfaces
    for( index_t interface_itr = 0; interface_itr < nb_interfaces;
            interface_itr++ )
    {
        builder_.geology.create_geological_entity(
            Interface3D::type_name_static() );
    }

    // Create the Layers
    for( index_t layer = 0; layer < nb_layers; layer++ )
    {
        builder_.geology.create_geological_entity(
            Layer3D::type_name_static() );
    }

    // Then we create the GeoModelMeshEntity
    // Create the Corners
    for( index_t corner = 0; corner < nb_corners; corner++ )
    {
        builder_.topology.create_mesh_entity( Corner3D::type_name_static() );
    }

    // Create the Lines
    for( index_t lines = 0; lines < nb_lines; lines++ )
    {
        builder_.topology.create_mesh_entity( Line3D::type_name_static() );
    }

    // Create the Surfaces
    for( index_t surface = 0; surface < nb_surfaces; surface++ )
    {
        builder_.topology.create_mesh_entity(
            Surface3D::type_name_static() );
    }

    // Create the Regions
    for( index_t region = 0; region < nb_regions; region++ )
    {
        builder_.topology.create_mesh_entity( Region3D::type_name_static() );
    }

    //#############################
    // Setting the Geometry       #
    //#############################

    // We declare the coordinates of the corners. We arrange the corner in a
    // table
    vec3 corners_table[4];
    corners_table[0] = vec3( 0, 0, 300 );
    corners_table[1] = vec3( 700, 0, 350);
    corners_table[2] = vec3( 0, 150, 300 );
    corners_table[3] = vec3( 0, 0, 500 );

    // We associate the coordinates with the corners
    for( index_t corner = 0; corner < nb_corners; corner++ )
    {
        builder_.geometry.set_corner( corner, corners_table[corner] );
    }

    // We associate the coordinates with the lines
    // We create a vector cur_coor_line containing the 2 vertices
    // for each line. Of course, you can have more vertices in a Line
    std::vector< vec3 > cur_coor_line( 2 );
    cur_coor_line[0] = corners_table[0];
    cur_coor_line[1] = corners_table[1];
    builder_.geometry.set_line( 0, cur_coor_line );

    cur_coor_line[0] = corners_table[1];
    cur_coor_line[1] = corners_table[2];
    builder_.geometry.set_line( 1, cur_coor_line );

    cur_coor_line[0] = corners_table[0];
    cur_coor_line[1] = corners_table[2];
    builder_.geometry.set_line( 2, cur_coor_line );

    cur_coor_line[0] = corners_table[0];
    cur_coor_line[1] = corners_table[3];
    builder_.geometry.set_line( 3, cur_coor_line );

    cur_coor_line[0] = corners_table[1];
    cur_coor_line[1] = corners_table[3];
    builder_.geometry.set_line( 4, cur_coor_line );

    cur_coor_line[0] = corners_table[2];
    cur_coor_line[1] = corners_table[3];
    builder_.geometry.set_line( 5, cur_coor_line );

    // We associate the coordinates with the Surfaces
    // We create a vector cur_coor_surface containing 4 vertices.
    // These 4 vertices delimits each surface so each surface
    // will contain one unique quad as a facet.
    // You can defined a more complicated mesh (for example a
    // triangular mesh) with these methods.
    std::vector< index_t > facet( 3, 0 );
    facet[0] = 0;
    facet[1] = 1;
    facet[2] = 2;

    std::vector< index_t > facet_ptr( 2 );
    facet_ptr[0] = 0;
    facet_ptr[1] = 3;
    std::vector< vec3 > cur_coor_surface( 3 );
    cur_coor_surface[0] = corners_table[0];
    cur_coor_surface[1] = corners_table[1];
    cur_coor_surface[2] = corners_table[2];

    builder_.geometry.set_surface_geometry(
        0, cur_coor_surface, facet, facet_ptr );

    cur_coor_surface[0] = corners_table[0];
    cur_coor_surface[1] = corners_table[1];
    cur_coor_surface[2] = corners_table[3];

    builder_.geometry.set_surface_geometry(
        1, cur_coor_surface, facet, facet_ptr );

    cur_coor_surface[0] = corners_table[1];
    cur_coor_surface[1] = corners_table[2];
    cur_coor_surface[2] = corners_table[3];

    builder_.geometry.set_surface_geometry(
        2, cur_coor_surface, facet, facet_ptr );

    cur_coor_surface[0] = corners_table[0];
    cur_coor_surface[1] = corners_table[2];
    cur_coor_surface[2] = corners_table[3];

    builder_.geometry.set_surface_geometry(
        3, cur_coor_surface, facet, facet_ptr );

    //###################################
    // Setting the Boundaries relations #
    //###################################

    // We set the Corners which are incident entities of the lines
    // The add_mesh_entity_boundary_relation method take as first argument
    // the
    // gme_t of the boundary and in second argument
    // the id of the GeoModelMeshentity bounded by the boundary
    // Remember :
    // Lines are bounded by Corners
    // Surfaces are bounded by Lines
    // Region are bounded by Surfaces

    // For corner 0
    // Corner 0 is a boundary of the lines: 0, 3, and 13.
    builder_.topology.add_line_corner_boundary_relation( 0, 0 );
    builder_.topology.add_line_corner_boundary_relation( 2, 0 );
    builder_.topology.add_line_corner_boundary_relation( 3, 0 );

    // For corner 1
    // Corner 1 is a boundary of the lines: 2, 3, and 17.
    builder_.topology.add_line_corner_boundary_relation( 0, 1 );
    builder_.topology.add_line_corner_boundary_relation( 1, 1 );
    builder_.topology.add_line_corner_boundary_relation( 4, 1 );

    // For corner 2
    // Corner 2 is a boundary of the lines: 1, 2, and 19.
    builder_.topology.add_line_corner_boundary_relation( 1, 2 );
    builder_.topology.add_line_corner_boundary_relation( 2, 2 );
    builder_.topology.add_line_corner_boundary_relation( 5, 2 );

    // For corner 3
    // Corner 3 is a boundary of the lines: 0, 1, and 15.
    builder_.topology.add_line_corner_boundary_relation( 3, 3 );
    builder_.topology.add_line_corner_boundary_relation( 4, 3 );
    builder_.topology.add_line_corner_boundary_relation( 5, 3 );


    /////////////////////////////////////////////////////////

    // For line 0
    // Line 0 is a boundary of the surfaces: 0 and 4.
    builder_.topology.add_surface_line_boundary_relation( 0, 0 );
    builder_.topology.add_surface_line_boundary_relation( 1, 0 );

    // For line 1
    // Line 1 is a boundary of the surfaces: 0 and 10.
    builder_.topology.add_surface_line_boundary_relation( 0, 1 );
    builder_.topology.add_surface_line_boundary_relation( 2, 1 );

    // For line 2
    // Line 2 is a boundary of the surfaces: 0 and 6.
    builder_.topology.add_surface_line_boundary_relation( 0, 2 );
    builder_.topology.add_surface_line_boundary_relation( 3, 2 );

    // For line 3
    // Line 3 is a boundary of the surfaces: 0 and 8.
    builder_.topology.add_surface_line_boundary_relation( 1, 3 );
    builder_.topology.add_surface_line_boundary_relation( 3, 3 );

    // For line 4
    // Line 4 is a boundary of the surfaces: 1, 3 and 4.
    builder_.topology.add_surface_line_boundary_relation( 1, 4 );
    builder_.topology.add_surface_line_boundary_relation( 2, 4 );


    // For line 5
    // Line 5 is a boundary of the surfaces: 1, 9 and 10.
    builder_.topology.add_surface_line_boundary_relation( 2, 5 );
    builder_.topology.add_surface_line_boundary_relation( 3, 5 );

    /////////////////////////////////////////////////////////

    // For surface 0
    // Surface 0 is a boundary of the region 0.
    builder_.topology.add_region_surface_boundary_relation(
        0, 0, false ); // TODO side ????

    // For surface 1
    // Surface 1 is a boundary of the region 0.
    builder_.topology.add_region_surface_boundary_relation(
        0, 1, false ); // TODO side ????

    // For surface 2
    // Surface 2 is a boundary of the region 1.
    builder_.topology.add_region_surface_boundary_relation(
        0, 2, false ); // TODO side ????

    // For surface 3
    // Surface 3 is a boundary of the region 1.
    builder_.topology.add_region_surface_boundary_relation(
        0, 3, false ); // TODO side ????


    //#####################################
    // Setting the parent/child relations #
    //#####################################

    // Remember :
    // Child of a Contact is a Line
    // Child of an Interface is a Surface
    // Child of a Layer is a Region

    // We use the method "add_parent_children_relation"
    // First argument is the parent (ie a GeoModelGeologicalEntity)
    // Second argument is the index of the child (ie a GeoModelMeshEntity)

    // For Contact 0
    builder_.geology.add_parent_children_relation(
        gmge_id( Contact3D::type_name_static(), 0 ),
        gmme_id( Line3D::type_name_static(), 0 ) );

    // For Contact 1
    builder_.geology.add_parent_children_relation(
        gmge_id( Contact3D::type_name_static(), 1 ),
        gmme_id( Line3D::type_name_static(), 1 ) );

    // For Contact 2
    builder_.geology.add_parent_children_relation(
        gmge_id( Contact3D::type_name_static(), 2 ),
        gmme_id( Line3D::type_name_static(), 2 ) );

    // For Contact 3
    builder_.geology.add_parent_children_relation(
        gmge_id( Contact3D::type_name_static(), 3 ),
        gmme_id( Line3D::type_name_static(), 3 ) );

    // For Contact 4
    builder_.geology.add_parent_children_relation(
        gmge_id( Contact3D::type_name_static(), 4 ),
        gmme_id( Line3D::type_name_static(), 4 ) );

    // For Contact 5
    builder_.geology.add_parent_children_relation(
        gmge_id( Contact3D::type_name_static(), 5 ),
        gmme_id( Line3D::type_name_static(), 5 ) );


    /////////////////////////////////////////////////

    // For Interface 0
    builder_.geology.add_parent_children_relation(
        gmge_id( Interface3D::type_name_static(), 0 ),
        gmme_id( Surface3D::type_name_static(), 0 ) );

    // For Interface 1
    builder_.geology.add_parent_children_relation(
        gmge_id( Interface3D::type_name_static(), 1 ),
        gmme_id( Surface3D::type_name_static(), 1 ) );

    // For Interface 2
    builder_.geology.add_parent_children_relation(
        gmge_id( Interface3D::type_name_static(), 2 ),
        gmme_id( Surface3D::type_name_static(), 2 ) );

    // For Interface 3
    builder_.geology.add_parent_children_relation(
        gmge_id( Interface3D::type_name_static(), 3 ),
        gmme_id( Surface3D::type_name_static(), 3 ) );


    ///////////////////////////////////////////////////

    // For Layer 0
    builder_.geology.add_parent_children_relation(
        gmge_id( Layer3D::type_name_static(), 0 ),
        gmme_id( Region3D::type_name_static(), 0 ) );

    return true;
}


/*****************************************************************************/

GeoModelBuilderRESQML::GeoModelBuilderRESQML(
    GeoModel3D& geomodel, const std::string& filename
) : GeoModelBuilderFile( geomodel, std::move( filename ) ),
    impl_(new GeoModelBuilderRESQMLImpl(*this))
{

}

GeoModelBuilderRESQML::~GeoModelBuilderRESQML(){
    //needed due to the unique_ptr impl_
}

void GeoModelBuilderRESQML::load_file(){
    ringmesh_assert(impl_->load_file());
}

} // namespace RINGMesh

