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
#include <ringmesh/geomodel/core/stratigraphic_column.h>

#include <ringmesh/io/geomodel_adapter_resqml.h>

#include <ringmesh/mesh/mesh_index.h>

#include <ringmesh/mesh/surface_mesh.h>
#include <ringmesh/mesh/volume_mesh.h>

#include <ringmesh/mesh/mesh_builder.h>

#include <fesapi/common/EpcDocument.h>
#include <fesapi/common/HdfProxy.h>
#include <fesapi/resqml2/AbstractFeature.h>
#include <fesapi/resqml2_0_1/BoundaryFeature.h>
#include <fesapi/resqml2_0_1/ContinuousProperty.h>
#include <fesapi/resqml2_0_1/DiscreteProperty.h>
#include <fesapi/resqml2_0_1/FrontierFeature.h>
#include <fesapi/resqml2_0_1/Horizon.h>
#include <fesapi/resqml2_0_1/HorizonInterpretation.h>
#include <fesapi/resqml2_0_1/LocalDepth3dCrs.h>
#include <fesapi/resqml2_0_1/LocalTime3dCrs.h>
#include <fesapi/resqml2_0_1/StratigraphicColumn.h>
#include <fesapi/resqml2_0_1/StratigraphicColumnRankInterpretation.h>
#include <fesapi/resqml2_0_1/StratigraphicUnitFeature.h>
#include <fesapi/resqml2_0_1/StratigraphicUnitInterpretation.h>
#include <fesapi/resqml2_0_1/TectonicBoundaryFeature.h>
#include <fesapi/resqml2_0_1/TriangulatedSetRepresentation.h>
#include <fesapi/resqml2_0_1/UnstructuredGridRepresentation.h>

#include <fesapi/tools/GuidTools.h>

/*!
 * @brief Implementation of the class to write GeoModel to
 * RESQML2 .epc file
 * @author Wan-Chiu Li
 */

namespace RINGMesh
{
    using namespace RESQML2_0_1_NS;
    using namespace RESQML2_NS;
    using namespace COMMON_NS;
    using GMGE = GeoModelGeologicalEntity< 3 >;

    class GeoModelAdapterRESQMLImpl
    {
    public:
        GeoModelAdapterRESQMLImpl(
            const GeoModel3D& geomodel, const std::string& filename );
        ~GeoModelAdapterRESQMLImpl() = default;

        bool init();
        bool save_file();
        void serialize();
        bool write_surfaces();
        bool write_volumes();

        bool write_property( bool first_patch,
            std::vector< AbstractValuesProperty* >& properties,
            GEO::AttributesManager& attri_manager,
            gsoap_resqml2_0_1::resqml2__IndexableElements attachement,
            AbstractRepresentation& rep );
        bool write_stratigraphic_column();

        AbstractFeature* find_or_create_feature( const GMGE& entity );
        AbstractFeatureInterpretation* find_or_create_interpretation(
            AbstractFeature& feature );

    private:
        const GeoModel3D& geomodel_;
        const std::string& filename_;
        std::unique_ptr< EpcDocument > pck_;
        AbstractHdfProxy* hdf_proxy_;
        LocalDepth3dCrs* local_3d_crs_;

        std::map< gmge_id, AbstractFeature* > geo_entity_2_feature_;
        std::map< AbstractFeature*, AbstractFeatureInterpretation* >
            feature_2_interp_;
        std::map< std::string, PropertyKind* > type_2_property_kind_;
    };

    GeoModelAdapterRESQMLImpl::GeoModelAdapterRESQMLImpl(
        const GeoModel3D& geomodel, const std::string& filename )
        : geomodel_( geomodel ),
          filename_( filename ),
          pck_( nullptr ),
          hdf_proxy_( nullptr ),
          local_3d_crs_( nullptr )
    {
        ringmesh_assert( init() );
    }

    bool GeoModelAdapterRESQMLImpl::init()
    {
        pck_ = std::unique_ptr< EpcDocument >(
            new EpcDocument( filename_, EpcDocument::OVERWRITE ) );

        hdf_proxy_ = pck_->createHdfProxy( "", "Hdf Proxy",
            pck_->getStorageDirectory(), pck_->getName() + ".h5" );

        local_3d_crs_ = pck_->createLocalDepth3dCrs( "", "Default local CRS",
            .0, .0, .0, .0, gsoap_resqml2_0_1::eml20__LengthUom__m, 23031,
            gsoap_resqml2_0_1::eml20__LengthUom__m, "Unknown", false );

        return true;
    }

    AbstractFeature* GeoModelAdapterRESQMLImpl::find_or_create_feature(
        const GMGE& entity )
    {
        AbstractFeature* feature = nullptr;
        auto result = geo_entity_2_feature_.find( entity.gmge() );
        if( result == geo_entity_2_feature_.end() )
        {
            GMGE::GEOL_FEATURE geo_feat = entity.geological_feature();
            std::string guid( tools::GuidTools::generateUidAsString() );
            if( GMGE::is_fault( geo_feat ) )
            {
                feature = pck_->createFault( guid, entity.name() );
            }
            else if( GMGE::is_stratigraphic_limit( geo_feat ) )
            {
                feature = pck_->createHorizon( guid, entity.name() );
            }
            else if( GMGE::GEOL_FEATURE::VOI == geo_feat )
            {
                feature = pck_->createFrontier( guid, entity.name() );
            }
            else if( GMGE::GEOL_FEATURE::STRATI_UNIT == geo_feat )
            {
                feature = pck_->createStratigraphicUnit( guid, entity.name() );
            }
            else
            {
                feature = pck_->createBoundaryFeature( guid, entity.name() );
            }
            geo_entity_2_feature_[entity.gmge()] = feature;
        }
        else
        {
            feature = result->second;
        }
        return feature;
    }

    AbstractFeatureInterpretation*
        GeoModelAdapterRESQMLImpl::find_or_create_interpretation(
            AbstractFeature& feature )
    {
        AbstractFeatureInterpretation* interp = nullptr;
        auto result = feature_2_interp_.find( &feature );
        if( result == feature_2_interp_.end() )
        {
            std::string guid( tools::GuidTools::generateUidAsString() );
            if( dynamic_cast< TectonicBoundaryFeature* >( &feature )
                != nullptr )
            {
                interp =
                    (AbstractFeatureInterpretation*)
                        pck_->createFaultInterpretation(
                            (TectonicBoundaryFeature*) &feature, guid, guid );
            }
            else if( dynamic_cast< Horizon* >( &feature ) != nullptr )
            {
                interp = (AbstractFeatureInterpretation*)
                             pck_->createHorizonInterpretation(
                                 (Horizon*) &feature, guid, guid );
            }
            else if( dynamic_cast< FrontierFeature* >( &feature ) != nullptr )
            {
                interp = (AbstractFeatureInterpretation*)
                             pck_->createGenericFeatureInterpretation(
                                 &feature, guid, guid );
            }
            else if( dynamic_cast< StratigraphicUnitFeature* >( &feature )
                     != nullptr )
            {
                interp =
                    (AbstractFeatureInterpretation*)
                        pck_->createStratigraphicUnitInterpretation(
                            (StratigraphicUnitFeature*) &feature, guid, guid );
            }
            else
            {
                interp = (AbstractFeatureInterpretation*)
                             pck_->createBoundaryFeatureInterpretation(
                                 (BoundaryFeature*) &feature, guid, guid );
            }

            feature_2_interp_[&feature] = interp;
        }
        else
        {
            interp = result->second;
        }
        return interp;
    }

    void GeoModelAdapterRESQMLImpl::serialize()
    {
        hdf_proxy_->close();

        std::cout << "Start serialization of " << pck_->getName() << " in "
                  << ( pck_->getStorageDirectory().empty()
                             ? "working directory."
                             : pck_->getStorageDirectory() )
                  << std::endl;

        pck_->serialize();
    }

    bool GeoModelAdapterRESQMLImpl::write_property( bool first_patch,
        std::vector< AbstractValuesProperty* >& properties,
        GEO::AttributesManager& attri_manager,
        gsoap_resqml2_0_1::resqml2__IndexableElements attachment,
        AbstractRepresentation& rep )
    {
        GEO::vector< std::string > attribute_names;
        attri_manager.list_attribute_names( attribute_names );
        if( attribute_names.size() == 0 )
        {
            return true;
        }

        if( first_patch )
        {
            properties.resize( attribute_names.size(), nullptr );

            for( auto i : range( attribute_names.size() ) )
            {
                if( attribute_names[i] == "point" )
                {
                    properties[i] = nullptr;
                    continue;
                }

                GEO::AttributeStore* store =
                    attri_manager.find_attribute_store( attribute_names[i] );
                std::string element_type = GEO::AttributeStore::
                    element_type_name_by_element_typeid_name(
                        store->element_typeid_name() );
                if( !GEO::AttributeStore::element_typeid_name_is_known(
                        element_type ) )
                {
                    Logger::warn( "I/O",
                        "Skipping attribute: ", attribute_names[i],
                        " of type: ", element_type );
                    properties[i] = nullptr;
                    continue;
                }

                std::string guid( tools::GuidTools::generateUidAsString() );
                if( element_type == "double" )
                {
                    properties[i] = pck_->createContinuousProperty( &rep, guid,
                        attribute_names[i], store->dimension(), attachment,
                        "unitless",
                        gsoap_resqml2_0_1::
                            resqml2__ResqmlPropertyKind__continuous );
                }
                else if( element_type == "int" )
                {
                    properties[i] = pck_->createDiscreteProperty( &rep, guid,
                        attribute_names[i], store->dimension(), attachment,
                        gsoap_resqml2_0_1::
                            resqml2__ResqmlPropertyKind__discrete );
                }
            }
        }

        for( auto i : range( attribute_names.size() ) )
        {
            if( properties[i] == nullptr )
            {
                continue;
            }
            GEO::AttributeStore* store =
                attri_manager.find_attribute_store( attribute_names[i] );

            if( dynamic_cast< ContinuousProperty* >( properties[i] )
                != nullptr )
            {
                ContinuousProperty* continuous =
                    static_cast< ContinuousProperty* >( properties[i] );

                continuous->pushBackDoubleHdf5Array1dOfValues(
                    (double*) store->data(), store->size() * store->dimension(),
                    hdf_proxy_ );
            }
            else if( dynamic_cast< DiscreteProperty* >( properties[i] )
                     != nullptr )
            {
                DiscreteProperty* discrete =
                    static_cast< DiscreteProperty* >( properties[i] );

                discrete->pushBackIntHdf5Array1dOfValues( (int*) store->data(),
                    store->size() * store->dimension(), hdf_proxy_, -99999 );
            }
            else
            {
                ringmesh_assert_not_reached;
            }
        }
        return true;
    }

    bool GeoModelAdapterRESQMLImpl::write_surfaces()
    {
        for( const auto& interface :
            geomodel_.geol_entities( Interface3D::type_name_static() ) )
        {
            for( auto i : range( interface.nb_children() ) )
            {
                const Surface3D& surface =
                    static_cast< const Surface3D& >( interface.child( i ) );
                if( !surface.is_simplicial() )
                {
                    Logger::err( "", surface.gmme(),
                        "Non simplicial surfaces are not supported" );
                    return false;
                }
            }
        }

        for( const auto& interface :
            geomodel_.geol_entities( Interface3D::type_name_static() ) )
        {
            GMGE::GEOL_FEATURE geo_feat = interface.geological_feature();
            if( !interface.has_geological_feature()
                || !( GMGE::is_fault( geo_feat )
                       || GMGE::is_stratigraphic_limit( geo_feat )
                       || GMGE::GEOL_FEATURE::VOI == geo_feat ) )
            {
                Logger::warn( "", interface.gmge(),
                    " geological feature is empty or not recognized" );
            }

            AbstractFeature* feature = find_or_create_feature( interface );
            AbstractFeatureInterpretation* interp =
                find_or_create_interpretation( *feature );

            std::string guid( tools::GuidTools::generateUidAsString() );
            TriangulatedSetRepresentation* rep =
                pck_->createTriangulatedSetRepresentation(
                    interp, local_3d_crs_, guid, feature->getTitle() );

            unsigned int interface_vertex_count = 0;
            for( auto i : range( interface.nb_children() ) )
            {
                const Surface3D& surface =
                    static_cast< const Surface3D& >( interface.child( i ) );

                std::unique_ptr< double[] > points(
                    new double[surface.nb_vertices() * 3] );

                vec3 p;
                for( auto v : range( surface.nb_vertices() ) )
                {
                    p = surface.vertex( v );
                    points[v * 3] = p[0];
                    points[v * 3 + 1] = p[1];
                    points[v * 3 + 2] = p[2];
                }

                std::unique_ptr< unsigned int[] > node_indices(
                    new unsigned int[surface.nb_mesh_elements() * 3] );
                for( auto t : range( surface.nb_mesh_elements() ) )
                {
                    for( auto v :
                        range( surface.nb_mesh_element_vertices( t ) ) )
                    {
                        node_indices[t * 3 + v] =
                            interface_vertex_count
                            + surface.mesh_element_vertex_index( { t, v } );
                    }
                }

                rep->pushBackTrianglePatch( surface.nb_vertices(), &points[0],
                    surface.nb_mesh_elements(), &node_indices[0], hdf_proxy_ );

                interface_vertex_count += surface.nb_vertices();
            }

            if( interface_vertex_count == 0 )
            {
                return true;
            }

            std::vector< AbstractValuesProperty* > vertex_properties;
            std::vector< AbstractValuesProperty* > polygon_properties;
            for( auto i : range( interface.nb_children() ) )
            {
                const Surface3D& surface =
                    static_cast< const Surface3D& >( interface.child( i ) );

                write_property( i == 0, vertex_properties,
                    surface.vertex_attribute_manager(),
                    gsoap_resqml2_0_1::resqml2__IndexableElements__nodes,
                    *rep );
                write_property( i == 0, polygon_properties,
                    surface.polygon_attribute_manager(),
                    gsoap_resqml2_0_1::resqml2__IndexableElements__cells,
                    *rep );
            }
        }
        return true;
    }

    bool GeoModelAdapterRESQMLImpl::write_volumes()
    {
        for( const auto& layer :
            geomodel_.geol_entities( Layer3D::type_name_static() ) )
        {
            AbstractFeatureInterpretation* interp = nullptr;
            if( layer.has_geological_feature() )
            {
                AbstractFeature* feature = find_or_create_feature( layer );
                AbstractFeatureInterpretation* interp =
                    find_or_create_interpretation( *feature );
            }
            else
            {
                Logger::warn( "", layer.gmge(),
                    " geological feature is empty or not recognized" );
            }

            // TODO: I do not understand why the interp is not needed...
            ringmesh_unused( interp );

            std::vector< UnstructuredGridRepresentation* > reps;
            for( auto i : range( layer.nb_children() ) )
            {
                const Region3D& region =
                    static_cast< const Region3D& >( layer.child( i ) );

                std::string guid( tools::GuidTools::generateUidAsString() );
                UnstructuredGridRepresentation* rep =
                    pck_->createUnstructuredGridRepresentation( local_3d_crs_,
                        guid, region.name(), region.nb_mesh_elements() );
                reps.push_back( rep );

                std::unique_ptr< double[] > points(
                    new double[region.nb_vertices() * 3] );

                vec3 p;
                for( auto v : range( region.nb_vertices() ) )
                {
                    p = region.vertex( v );
                    points[v * 3] = p[0];
                    points[v * 3 + 1] = p[1];
                    points[v * 3 + 2] = p[2];
                }

                std::unique_ptr< ULONG64[] > cumul_faces_cells(
                    new ULONG64[region.nb_mesh_elements()] );

                std::vector< ULONG64 > face_indices_per_cell;
                std::vector< ULONG64 > cumul_vertices_face;
                std::vector< ULONG64 > node_indices_per_face;
                std::vector< unsigned char > face_righthandness;

                index_t facet_count = 0;
                for( auto t : range( region.nb_mesh_elements() ) )
                {
                    for( auto f : range( region.nb_cell_facets( t ) ) )
                    {
                        for( auto v :
                            range( region.nb_cell_facet_vertices( t, f ) ) )
                        {
                            node_indices_per_face.push_back(
                                region.cell_facet_vertex_index( t, f, v ) );
                        }
                        face_indices_per_cell.push_back( facet_count++ );

                        // TODO compute real face righthandness
                        face_righthandness.push_back( 1 );
                        cumul_vertices_face.push_back(
                            node_indices_per_face.size() );
                    }
                    cumul_faces_cells[t] = facet_count;
                }

                rep->setGeometry( &face_righthandness[0], &points[0],
                    region.nb_vertices(), hdf_proxy_, &face_indices_per_cell[0],
                    &cumul_faces_cells[0], face_righthandness.size(),
                    &node_indices_per_face[0], &cumul_vertices_face[0],
                    gsoap_resqml2_0_1::resqml2__CellShape__polyhedral );
            }

            for( auto i : range( layer.nb_children() ) )
            {
                const Region3D& region =
                    static_cast< const Region3D& >( layer.child( i ) );

                std::vector< AbstractValuesProperty* > vertex_properties;
                std::vector< AbstractValuesProperty* > cell_properties;

                // first_patch is always true because fesapi only support
                // UnstructuredGrid with only one patch
                write_property( true, vertex_properties,
                    region.vertex_attribute_manager(),
                    gsoap_resqml2_0_1::resqml2__IndexableElements__nodes,
                    *reps[i] );
                write_property( true, cell_properties,
                    region.cell_attribute_manager(),
                    gsoap_resqml2_0_1::resqml2__IndexableElements__cells,
                    *reps[i] );
            }
        }

        return true;
    }

    namespace
    {
        gsoap_resqml2_0_1::resqml2__ContactMode get_contact_mode(
            RELATION relation )
        {
            switch( relation )
            {
            case RELATION::BASELAP:
                return gsoap_resqml2_0_1::resqml2__ContactMode__baselap;
            case RELATION::ERODED:
                return gsoap_resqml2_0_1::resqml2__ContactMode__erosion;
            case RELATION::CONFORMABLE:
                return gsoap_resqml2_0_1::resqml2__ContactMode__proportional;
            default:
                ringmesh_assert_not_reached;
                return gsoap_resqml2_0_1::resqml2__ContactMode__proportional;
            }
        }
    } // namespace

    bool GeoModelAdapterRESQMLImpl::write_stratigraphic_column()
    {
        const StratigraphicColumn* column = geomodel_.stratigraphic_column();
        if( column == nullptr )
        {
            // nothing to do
            return true;
        }

        RESQML2_0_1_NS::StratigraphicColumn* stratiColumn =
            pck_->createStratigraphicColumn(
                tools::GuidTools::generateUidAsString(), column->get_name() );
        OrganizationFeature* stratiModelFeature =
            pck_->createStratigraphicModel(
                tools::GuidTools::generateUidAsString(), "stratiModel" );
        StratigraphicColumnRankInterpretation* stratiColumnRank =
            pck_->createStratigraphicColumnRankInterpretationInAge(
                stratiModelFeature, tools::GuidTools::generateUidAsString(),
                "Stratigraphic column rank", 0 );
        stratiColumn->pushBackStratiColumnRank( stratiColumnRank );

        const NestedStratigraphicUnit::RankedUnits units =
            column->get_units_with_rank( 0 );

        for( auto unit_index : range( units.size() ) )
        {
            AbstractFeature* unit_feature =
                find_or_create_feature( *units[unit_index]->get_layer() );
            StratigraphicUnitInterpretation* unit_interp =
                dynamic_cast< StratigraphicUnitInterpretation* >(
                    find_or_create_interpretation( *unit_feature ) );
            ringmesh_assert( unit_interp != nullptr );

            stratiColumnRank->pushBackStratiUnitInterpretation( unit_interp );

            if( unit_index == ( units.size() - 1 ) )
            {
                break;
            }

            AbstractFeature* unit_below_feature =
                find_or_create_feature( *units[unit_index + 1]->get_layer() );
            StratigraphicUnitInterpretation* unit_below_interp =
                dynamic_cast< StratigraphicUnitInterpretation* >(
                    find_or_create_interpretation( *unit_below_feature ) );
            ringmesh_assert( unit_below_interp != nullptr );

            ringmesh_assert(
                units[unit_index]->get_interface_base() != nullptr );
            AbstractFeature* horizon_feature = find_or_create_feature(
                *units[unit_index]->get_interface_base() );
            HorizonInterpretation* horizon_interp =
                dynamic_cast< HorizonInterpretation* >(
                    find_or_create_interpretation( *horizon_feature ) );
            ringmesh_assert( horizon_interp != nullptr );

            stratiColumnRank->pushBackStratigraphicBinaryContact( unit_interp,
                get_contact_mode( units[unit_index]->get_relation_base() ),
                unit_below_interp,
                get_contact_mode( units[unit_index + 1]->get_relation_top() ),
                horizon_interp );
        }
        return true;
    }

    bool GeoModelAdapterRESQMLImpl::save_file()
    {
        try
        {
            write_surfaces();
            write_volumes();
            write_stratigraphic_column();
            serialize();
        }
        catch( const std::invalid_argument& Exp )
        {
            std::cerr << "Error : " << Exp.what() << ".\n";
        }

        return true;
    }

    /*****************************************************************************/

    GeoModelAdapterRESQML::GeoModelAdapterRESQML(
        const GeoModel3D& geomodel, const std::string& filename )
        : impl_( new GeoModelAdapterRESQMLImpl( geomodel, filename ) )
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
