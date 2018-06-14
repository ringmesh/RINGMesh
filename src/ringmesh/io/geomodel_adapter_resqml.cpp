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
#include <ringmesh/geomodel/core/geomodel_geological_entity.h>

#include <ringmesh/geomodel/core/geomodel_mesh_entity.h>
#include <ringmesh/io/geomodel_adapter_resqml.h>

#include <ringmesh/mesh/mesh_index.h>

#include <ringmesh/mesh/surface_mesh.h>
#include <ringmesh/mesh/volume_mesh.h>

#include <ringmesh/mesh/mesh_builder.h>

#include <fesapi/common/EpcDocument.h>
#include <fesapi/common/HdfProxy.h>
#include <fesapi/resqml2/AbstractFeature.h>
#include <fesapi/resqml2_0_1/LocalDepth3dCrs.h>
#include <fesapi/resqml2_0_1/LocalTime3dCrs.h>
#include <fesapi/resqml2_0_1/BoundaryFeature.h>
#include <fesapi/resqml2_0_1/FrontierFeature.h>
#include <fesapi/resqml2_0_1/TectonicBoundaryFeature.h>
#include <fesapi/resqml2_0_1/Horizon.h>
#include <fesapi/resqml2_0_1/TriangulatedSetRepresentation.h>
#include <fesapi/resqml2_0_1/UnstructuredGridRepresentation.h>
#include <fesapi/tools/GuidTools.h>

/*!
 * @brief Implementation of the class to build GeoModel from input
 * RESQML2 .epc file
 * @author Wan-Chiu Li
 */

namespace{
//     bool save_volume_rep( VolumeMesh3D& mesh,
//         UnstructuredGridRepresentation& unstructed_grid )
//     {

}


namespace RINGMesh
{
    using namespace RESQML2_0_1_NS;
    using namespace RESQML2_NS;
    using namespace COMMON_NS;
    using GMGE = GeoModelGeologicalEntity<3>;

    class GeoModelAdapterRESQMLImpl
    {
    public:
        GeoModelAdapterRESQMLImpl(
            const GeoModel3D& geomodel,
            const std::string& filename
        );
        ~GeoModelAdapterRESQMLImpl() = default;

        bool init();

        bool save_file();

        void serialize();
        bool save_surfaces();
     //   bool save_volumes();

        AbstractFeature* find_or_create_feature(
            const Interface3D& interface
        );
        AbstractFeatureInterpretation*
        find_or_create_interpretation(
            AbstractFeature& feature
        );
    private:
        const GeoModel3D& geomodel_;
        const std::string& filename_;
        std::unique_ptr<EpcDocument> pck_;
        AbstractHdfProxy* hdf_proxy_;
        LocalDepth3dCrs* local_3d_crs_;
        LocalTime3dCrs* local_time_3d_crs_;

        std::map<gmge_id, AbstractFeature*> geo_entity_2_feature_;

        std::map<AbstractFeature*,AbstractFeatureInterpretation*> 
        feature_2_interp_;
    };

    GeoModelAdapterRESQMLImpl::GeoModelAdapterRESQMLImpl(
        const GeoModel3D& geomodel,
        const std::string& filename
    ) : geomodel_( geomodel ),
        filename_( filename ),
        pck_(nullptr),
        hdf_proxy_(nullptr),
        local_3d_crs_(nullptr),
        local_time_3d_crs_(nullptr)
    {
        ringmesh_assert(init());
    }

    bool GeoModelAdapterRESQMLImpl::init()
    {
        pck_ = std::unique_ptr<EpcDocument>(
                new EpcDocument(filename_, EpcDocument::OVERWRITE));

        hdf_proxy_ = pck_->createHdfProxy("", "Hdf Proxy", 
            pck_->getStorageDirectory(), pck_->getName() + ".h5");

        local_3d_crs_ = pck_->createLocalDepth3dCrs(
            "", "Default local CRS", .0, 
            .0, .0,.0, gsoap_resqml2_0_1::eml20__LengthUom__m, 23031, 
            gsoap_resqml2_0_1::eml20__LengthUom__m, "Unknown", false);

        local_time_3d_crs_ = pck_->createLocalTime3dCrs(
            "", "Default local time CRS", 1.0,     0.1, 15, .0,
            gsoap_resqml2_0_1::eml20__LengthUom__m, 23031,
            gsoap_resqml2_0_1::eml20__TimeUom__s,
            gsoap_resqml2_0_1::eml20__LengthUom__m,
            "Unknown", false);
        return true;
    }

    AbstractFeature* GeoModelAdapterRESQMLImpl::find_or_create_feature(
        const Interface3D& interface
    ){
        AbstractFeature* feature = nullptr;
        auto result = geo_entity_2_feature_.find(interface.gmge());
        if(result == geo_entity_2_feature_.end())
        {
            GMGE::GEOL_FEATURE geo_feat = interface.geological_feature();
            std::string guid(tools::GuidTools::generateUidAsString());
            if(GMGE::is_fault(geo_feat)
            ){
                feature = pck_->createFault(guid, "Fault");
            }else if(GMGE::is_stratigraphic_limit(geo_feat)
            ){
                feature = pck_->createHorizon(guid, "Horizon");
            }else if(GMGE::GEOL_FEATURE::VOI == geo_feat){
                feature = pck_->createFrontier(guid, "VOI");
            }else{
                ringmesh_assert_not_reached;
                return nullptr;
            }
            geo_entity_2_feature_[interface.gmge()] = feature;
        }else{
            feature = result->second;
        }
        return feature;
    }

    AbstractFeatureInterpretation*
    GeoModelAdapterRESQMLImpl::find_or_create_interpretation(
        AbstractFeature& feature
    ){
        AbstractFeatureInterpretation* interp = nullptr;
        auto result = feature_2_interp_.find(&feature);
        if(result == feature_2_interp_.end())
        {
            std::string guid(tools::GuidTools::generateUidAsString());
            if(dynamic_cast<TectonicBoundaryFeature*>(&feature)!=nullptr){
                interp
                    = (AbstractFeatureInterpretation*)
                        pck_->createFaultInterpretation(
                        (TectonicBoundaryFeature*)&feature,
                        guid,
                        "Fault1");
            }else if(dynamic_cast<Horizon*>(&feature)!=nullptr){
                interp
                    = (AbstractFeatureInterpretation*)
                        pck_->createHorizonInterpretation(
                        (Horizon*)&feature,
                        guid,
                        "Horizon1");
            }else if(dynamic_cast<FrontierFeature*>(&feature)!=nullptr){
                interp
                    = (AbstractFeatureInterpretation*)
                        pck_->createGenericFeatureInterpretation(
                        &feature,
                        guid,
                        "VOI1");
            }else{
                ringmesh_assert_not_reached;
                return nullptr;
            }

            feature_2_interp_[&feature] = interp;
        }else{
            interp = result->second;
        }
        return interp;
    }

    void GeoModelAdapterRESQMLImpl::serialize()
    {
        hdf_proxy_->close();

        std::cout << "Start serialization of " << pck_->getName() << " in " << 
            (pck_->getStorageDirectory().empty() ? "working directory." : 
            pck_->getStorageDirectory()) << std::endl;
            pck_->serialize();
    }

    bool GeoModelAdapterRESQMLImpl::save_surfaces()
    {
        for( const auto& cur_interface :
            geomodel_.geol_entities( Interface3D::type_name_static() ) )
        {
            GMGE::GEOL_FEATURE geo_feat = cur_interface.geological_feature();
            if( !cur_interface.has_geological_feature() ||
                !(GMGE::is_fault(geo_feat) ||
                  GMGE::is_stratigraphic_limit(geo_feat)||
                  GMGE::GEOL_FEATURE::VOI == geo_feat)
            )
            {
                Logger::err(
                    "", cur_interface.gmge(),
                    "geological feature is empty or not supported" );
                return false;
            }

            AbstractFeature* feature
                = find_or_create_feature(
                    static_cast<const Interface3D&>(cur_interface));
            AbstractFeatureInterpretation* interp
                = find_or_create_interpretation(*feature);

            std::string guid(tools::GuidTools::generateUidAsString());
            TriangulatedSetRepresentation* rep
                = pck_->createTriangulatedSetRepresentation(
                    interp, local_3d_crs_,
                    guid,
                    "nimp");

            unsigned int interface_vertex_count = 0;
            for( auto i : range( cur_interface.nb_children() ) )
            {
                const Surface3D& surface
                    = static_cast<const Surface3D&>(cur_interface.child( i ));

                std::unique_ptr< double[] > points(
                    new double[surface.nb_vertices() * 3] );

                vec3 p;
                for( auto v : range( surface.nb_vertices() ) )
                {
                    p = surface.vertex(v);
                    points[v*3]   = p[0];
                    points[v*3+1] = p[1];
                    points[v*3+2] = p[2];
                }

                std::unique_ptr< unsigned int[] > node_indices(
                    new unsigned int[surface.nb_mesh_elements() * 3] );
                for( auto t : range( surface.nb_mesh_elements() ) )
                {
                    for( auto v : range( surface.nb_mesh_element_vertices( t )))
                    {
                        node_indices[t*3 + v]
                            = interface_vertex_count
                                + surface.mesh_element_vertex_index({t, v });
                    }
                }

                rep->pushBackTrianglePatch(
                    surface.nb_vertices(), &points[0],
                    surface.nb_mesh_elements(), &node_indices[0],
                    hdf_proxy_);

                interface_vertex_count += surface.nb_vertices();
            }
        }

        return true;
    }

//     bool GeoModelAdapterRESQMLImpl::save_volumes()
//     {
//         std::vector< UnstructuredGridRepresentation* > unstructuredGridRepSet =
//             pck.getUnstructuredGridRepresentationSet();
//         if( unstructuredGridRepSet.empty() )
//         {
//             return true;
//         }
// 
//         std::cout << std::endl
//                   << "UNSTRUCTURED GRID REP: " << unstructuredGridRepSet.size()
//                   << std::endl;
// 
//         for( auto unstructured_grid : unstructuredGridRepSet )
//         {
//             showAllMetadata( unstructured_grid );
// 
//             if( !unstructured_grid->isPartial()
//                 && unstructured_grid->hasGeometry() )
//             {
//                 auto mesh = VolumeMesh3D::create_mesh();
//                 ringmesh_assert( read_volume_rep( *mesh, *unstructured_grid ) );
// 
//                 // the volume mesh from resqml is here, need to find the
//                 // corresponding region of in the GeoModel3D
//                 const auto& nn_search = mesh->cell_facet_nn_search();
// 
//                 int region_index = -1;
//                 for( auto r : range( geomodel_.nb_regions() ) )
//                 {
//                     const Region3D& region = geomodel_.region( r );
// 
//                     bool match = true;
//                     for( auto b : range( region.nb_boundaries() ) )
//                     {
//                         const Surface3D& surface = region.boundary( b );
//                         for( auto p : range( surface.nb_mesh_elements() ) )
//                         {
//                             auto center = surface.mesh_element_barycenter( p );
//                             auto result = nn_search.get_neighbors(
//                                 center, surface.geomodel().epsilon() );
//                             if( result.empty() )
//                             {
//                                 match = false;
//                                 break;
//                             }
//                         }
//                         if( !match )
//                         {
//                             break;
//                         }
//                     }
//                     if( match )
//                     {
//                         region_index = (int) r;
//                         break;
//                     }
//                 }
// 
//                 if( region_index == -1 )
//                 {
//                     return false;
//                 }
// 
//                 // corresponding region found, build its volume mesh
//                 const gmme_id region_id(
//                     region_type_name_static(), (index_t) region_index );
// 
//                 auto mesh_builder = builder_.geometry.create_region_builder(
//                     region_id.index() );
// 
//                 for( auto v : range( mesh->nb_vertices() ) )
//                 {
//                     mesh_builder->create_vertex( mesh->vertex( v ) );
//                 }
// 
//                 mesh_builder->create_cells(
//                     mesh->nb_cells(), CellType::TETRAHEDRON );
// 
//                 for( auto cell : range( mesh->nb_cells() ) )
//                 {read
//                     const index_t nb_vertices = mesh->nb_cell_vertices( cell );
//                     std::vector< index_t > cell_vertices( nb_vertices, 0 );
//                     for( auto v : range( nb_vertices ) )
//                     {
//                         ElementLocalVertex lv( cell, v );
//                         mesh_builder->set_cell_vertex(
//                             lv, mesh->cell_vertex( lv ) );
//                     }
//                 }
// 
//                 mesh_builder->connect_cells();
//             }
//         }
//         return true;
//     }

    bool GeoModelAdapterRESQMLImpl::save_file()
    {
        try {
            save_surfaces(  );

    //         save_volumes( pck );


            serialize();
        }catch (const std::invalid_argument & Exp)
        {
            std::cerr << "Error : " << Exp.what() << ".\n";
        }

        return true;
    }

/*****************************************************************************/

    GeoModelAdapterRESQML::GeoModelAdapterRESQML(
        const GeoModel3D& geomodel, const std::string& filename )
        : impl_( new GeoModelAdapterRESQMLImpl(geomodel, filename))
    {
    }

    GeoModelAdapterRESQML::~GeoModelAdapterRESQML()
    {
        // needed due to the unique_ptr impl_
    }

    void GeoModelAdapterRESQML::save_file()
    {
        ringmesh_assert( impl_->save_file() );
    }

} // namespace RINGMesh
