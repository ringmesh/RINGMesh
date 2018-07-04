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

#include <ringmesh/io/geomodel_builder_resqml.h>

#include <ringmesh/mesh/mesh_index.h>

#include <ringmesh/mesh/surface_mesh.h>
#include <ringmesh/mesh/volume_mesh.h>

#include <ringmesh/mesh/mesh_builder.h>

#include <fesapi/common/EpcDocument.h>
#include <fesapi/common/HdfProxy.h>
#include <fesapi/resqml2/AbstractFeatureInterpretation.h>
#include <fesapi/resqml2_0_1/CategoricalProperty.h>
#include <fesapi/resqml2_0_1/ContinuousProperty.h>
#include <fesapi/resqml2_0_1/DiscreteProperty.h>
#include <fesapi/resqml2_0_1/FrontierFeature.h>
#include <fesapi/resqml2_0_1/GeneticBoundaryFeature.h>
#include <fesapi/resqml2_0_1/Horizon.h>
#include <fesapi/resqml2_0_1/LocalDepth3dCrs.h>
#include <fesapi/resqml2_0_1/LocalTime3dCrs.h>
#include <fesapi/resqml2_0_1/PropertyKind.h>
#include <fesapi/resqml2_0_1/PropertyKindMapper.h>
#include <fesapi/resqml2_0_1/StratigraphicColumn.h>
#include <fesapi/resqml2_0_1/StratigraphicColumnRankInterpretation.h>
#include <fesapi/resqml2_0_1/StratigraphicUnitInterpretation.h>
#include <fesapi/resqml2_0_1/TectonicBoundaryFeature.h>
#include <fesapi/resqml2_0_1/TriangulatedSetRepresentation.h>
#include <fesapi/resqml2_0_1/UnstructuredGridRepresentation.h>

/*!
 * @brief Implementation of the class to build GeoModel from input
 * RESQML2 .epc file
 * @author Wan-Chiu Li
 */

namespace RINGMesh
{
    using namespace RESQML2_0_1_NS;
    using namespace RESQML2_NS;
    using namespace COMMON_NS;
    using GMGE = GeoModelGeologicalEntity< 3 >;

    namespace
    {
        void showAllMetadata(
            AbstractObject* obj, const std::string& prefix = "" )
        {
            std::cout << prefix << "Title is : " << obj->getTitle()
                      << std::endl;
            std::cout << prefix << "Guid is : " << obj->getUuid() << std::endl;
            if( !obj->isPartial() )
            {
                for( unsigned int i = 0; i < obj->getAliasCount(); ++i )
                {
                    std::cout
                        << prefix
                        << "Alias is : " << obj->getAliasAuthorityAtIndex( i )
                        << ":" << obj->getAliasTitleAtIndex( i ) << std::endl;
                }
                for( unsigned int i = 0; i < obj->getExtraMetadataCount(); ++i )
                {
                    std::cout << prefix << "Extrametadata is : "
                              << obj->getExtraMetadataKeyAtIndex( i ) << ":"
                              << obj->getExtraMetadataStringValueAtIndex( i )
                              << std::endl;
                }
            }
            else
            {
                std::cout << prefix << "IS PARTIAL!" << std::endl;
            }
            std::cout << prefix
                      << "--------------------------------------------------"
                      << std::endl;
        }

    } // anonymous namespace

    /****************************************************************************/
    struct UnitInfo
    {
        UnitInfo()
            : relation_top_( RELATION::CONFORMABLE ),
              relation_base_( RELATION::CONFORMABLE )
        {
        }

        gmge_id interface_top_;
        gmge_id interface_base_;
        gmge_id layer_;
        RELATION relation_top_;
        RELATION relation_base_;
        std::string name_;
    };

    class GeoModelBuilderRESQMLImpl
    {
    public:
        GeoModelBuilderRESQMLImpl(
            GeoModelBuilderRESQML& builder, GeoModel3D& geomodel );
        ~GeoModelBuilderRESQMLImpl() = default;

        bool load_file();

        void deserialize( EpcDocument& pck );
        void read_strati_column( EpcDocument& pck );
        bool read_surfaces( const EpcDocument& pck );
        bool read_volumes( const EpcDocument& pck );

    protected:
        gmge_id find_layer( const Region3D& region ) const;

    private:
        GeoModelBuilderRESQML& builder_;
        GeoModel3D& geomodel_;

        std::map< gmge_id, AbstractFeatureInterpretation* >
            geo_entity_2_interp_;
        std::map< AbstractFeatureInterpretation*, gmge_id >
            interp_2_geo_entity_;
        std::map< AbstractFeatureInterpretation*, UnitInfo > unit_2_info_;
    };

    GeoModelBuilderRESQMLImpl::GeoModelBuilderRESQMLImpl(
        GeoModelBuilderRESQML& builder, GeoModel3D& geomodel )
        : builder_( builder ), geomodel_( geomodel )
    {
    }

    void GeoModelBuilderRESQMLImpl::deserialize( EpcDocument& pck )
    {
        std::string resqmlResult = pck.deserialize();
        if( !resqmlResult.empty() )
        {
            std::cerr << resqmlResult << std::endl;
            std::cout << "Press enter to continue..." << std::endl;
            std::cin.get();
        }

        std::cout << "EpcDocument name " << pck.getName() << " in "
                  << ( pck.getStorageDirectory().empty()
                             ? "working directory."
                             : pck.getStorageDirectory() )
                  << std::endl;

        unsigned int hdfProxyCount = pck.getHdfProxyCount();
        std::cout << "There are " << pck.getHdfProxyCount()
                  << " hdf files associated to this epc document." << std::endl;
        for( unsigned int hdfProxyIndex = 0; hdfProxyIndex < hdfProxyCount;
             ++hdfProxyIndex )
        {
            std::cout << "Hdf file relative path : "
                      << pck.getHdfProxy( hdfProxyIndex )->getRelativePath()
                      << std::endl;
        }
        for( size_t warningIndex = 0; warningIndex < pck.getWarnings().size();
             ++warningIndex )
        {
            std::cout << "Warning #" << warningIndex << " : "
                      << pck.getWarnings()[warningIndex] << std::endl;
        }
    }

    namespace
    {
        RELATION get_contact_mode(
            gsoap_resqml2_0_1::resqml2__ContactMode mode )
        {
            switch( mode )
            {
            case gsoap_resqml2_0_1::resqml2__ContactMode__baselap:
                return RELATION::BASELAP;
            case gsoap_resqml2_0_1::resqml2__ContactMode__erosion:
                return RELATION::ERODED;
            case gsoap_resqml2_0_1::resqml2__ContactMode__proportional:
                return RELATION::CONFORMABLE;
            default:
                ringmesh_assert_not_reached;
                return RELATION::CONFORMABLE;
            }
        }
    } // namespace

    gmge_id GeoModelBuilderRESQMLImpl::find_layer(
        const Region3D& region ) const
    {
        std::set< const Interface3D* > horizons;
        for( auto b : range( region.nb_boundaries() ) )
        {
            const Surface3D& surface = region.boundary( b );
            const Interface3D& interface = static_cast< const Interface3D& >(
                surface.parent( Interface3D::type_name_static() ) );
            if( GMGE::is_stratigraphic_limit( interface.geological_feature() ) )
            {
                horizons.insert( &interface );
            }
        }

        if( horizons.size() >= 1 && horizons.size() <= 2 )
        {
            std::vector< gmge_id > horizons_list;
            for( auto h : horizons )
            {
                horizons_list.push_back( h->gmge() );
            }
            std::sort( horizons_list.begin(), horizons_list.end() );
            for( auto& unit : unit_2_info_ )
            {
                std::vector< gmge_id > cur_horizons_list;
                if( unit.second.interface_top_.is_defined() )
                {
                    cur_horizons_list.push_back( unit.second.interface_top_ );
                }
                if( unit.second.interface_base_.is_defined() )
                {
                    cur_horizons_list.push_back( unit.second.interface_base_ );
                }
                std::sort( cur_horizons_list.begin(), cur_horizons_list.end() );

                std::vector< gmge_id > result;
                std::set_intersection( horizons_list.begin(),
                    horizons_list.end(), cur_horizons_list.begin(),
                    cur_horizons_list.end(), std::back_inserter( result ) );
                if( result.size() == horizons_list.size() )
                {
                    auto id = interp_2_geo_entity_.find( unit.first );
                    ringmesh_assert( id != interp_2_geo_entity_.end() );
                    return id->second;
                }
            }
        }
    }

    void GeoModelBuilderRESQMLImpl::read_strati_column( EpcDocument& pck )
    {
        const std::vector< RESQML2_0_1_NS::StratigraphicColumn* >
            stratiColumnSet = pck.getStratigraphicColumnSet();
        ringmesh_assert( stratiColumnSet.size() <= 1 );
        for( auto column : stratiColumnSet )
        {
            showAllMetadata( column );
            for( auto rank_interp :
                column->getStratigraphicColumnRankInterpretationSet() )
            {
                std::cout << "\tCOLUMN RANK INTERP" << std::endl;
                showAllMetadata( rank_interp );
                ringmesh_assert( rank_interp->isAChronoStratiRank() );

                RockFeature rock( "rock", ROCKTYPE::NONE );

                for( auto unit :
                    rank_interp->getStratigraphicUnitInterpretationSet() )
                {
                    const gmge_id layer =
                        builder_.geology.create_geological_entity(
                            Layer3D::type_name_static() );

                    builder_.geology.set_geological_entity_geol_feature(
                        layer, GMGE::GEOL_FEATURE::STRATI_UNIT );

                    geo_entity_2_interp_[layer] = unit;
                    interp_2_geo_entity_[unit] = layer;

                    unit_2_info_[unit].layer_ = layer;
                    unit_2_info_[unit].name_ = unit->getTitle();
                }

                for( auto contactIndex :
                    range( rank_interp->getContactCount() ) )
                {
                    unit_2_info_[rank_interp->getSubjectOfContact(
                                     contactIndex )]
                        .relation_base_ = get_contact_mode(
                        rank_interp->getSubjectContactModeOfContact(
                            contactIndex ) );
                    unit_2_info_[rank_interp->getSubjectOfContact(
                                     contactIndex )]
                        .interface_base_ = interp_2_geo_entity_
                        [(AbstractFeatureInterpretation*)
                                rank_interp->getHorizonInterpretationOfContact(
                                    contactIndex )];

                    unit_2_info_[rank_interp->getDirectObjectOfContact(
                                     contactIndex )]
                        .relation_top_ = get_contact_mode(
                        rank_interp->getDirectObjectContactModeOfContact(
                            contactIndex ) );
                    unit_2_info_[rank_interp->getDirectObjectOfContact(
                                     contactIndex )]
                        .interface_top_ = interp_2_geo_entity_
                        [(AbstractFeatureInterpretation*)
                                rank_interp->getHorizonInterpretationOfContact(
                                    contactIndex )];
                }
            }

            // assign regions to layers
            // TODO this is a workaround until the
            // SealedVolumeFrameworkRepresentation
            // is released in fesapi that tells us the layer information of a
            // region
            for( auto& region : geomodel_.regions() )
            {
                gmge_id layer = find_layer( region );
                ringmesh_assert( layer.is_defined() );
                builder_.geology.add_parent_children_relation(
                    layer, region.gmme() );
            }

            // build the stati column
            RockFeature rock( "rock", ROCKTYPE::NONE );
            std::vector< const StratigraphicUnit* > units;

            for( const auto& unit : unit_2_info_ )
            {
                const Interface3D* top =
                    unit.second.interface_top_.is_defined()
                        ? dynamic_cast< const Interface3D* >(
                              &geomodel_.geological_entity(
                                  Interface3D::type_name_static(),
                                  unit.second.interface_top_.index() ) )
                        : nullptr;
                const Interface3D* base =
                    unit.second.interface_base_.is_defined()
                        ? dynamic_cast< const Interface3D* >(
                              &geomodel_.geological_entity(
                                  Interface3D::type_name_static(),
                                  unit.second.interface_base_.index() ) )
                        : nullptr;

                UnsubdividedStratigraphicUnit* sunit =
                    new UnsubdividedStratigraphicUnit( unit.second.name_, top,
                        base,
                        dynamic_cast< const Layer3D& >(
                            geomodel_.geological_entity(
                                Layer3D::type_name_static(), 0 ) ),
                        unit.second.relation_top_, unit.second.relation_base_,
                        rock, 0, 10 );

                units.push_back( sunit );

                std::cout << " relation base: "
                          << (int) unit.second.relation_base_ << std::endl;
                std::cout << " relation top: "
                          << (int) unit.second.relation_top_ << std::endl;
                std::cout << " interface base: "
                          << unit.second.interface_base_.index() << std::endl;
                std::cout << " interface top: "
                          << unit.second.interface_top_.index() << std::endl;
            }

            StratigraphicColumn* strati =
                new StratigraphicColumn( column->getTitle(), units,
                    STRATIGRAPHIC_PARADIGM::CHRONOSTRATIGRAPHIC );
            geomodel_.set_stratigraphic_column( strati );
        }

        if( stratiColumnSet.empty() )
        {
            if( !geomodel_.entity_type_manager()
                     .geological_entity_manager.is_valid_type(
                         Layer3D::type_name_static() ) )
            {
                const gmge_id layer = builder_.geology.create_geological_entity(
                    Layer3D::type_name_static() );

                builder_.geology.set_geological_entity_geol_feature(
                    layer, GMGE::GEOL_FEATURE::STRATI_UNIT );

                for( auto& region : geomodel_.regions() )
                {
                    builder_.geology.add_parent_children_relation(
                        layer, region.gmme() );
                }
            }
        }
    }

    namespace
    {
        bool read_property( AbstractValuesProperty& property,
            const unsigned int patch_index,
            GEO::AttributesManager& attri_manager )
        {
            if( property.getValuesHdfDatatype()
                == AbstractValuesProperty::UNKNOWN )
            {
                throw RINGMeshException( "RESQML2", "Invalid property type: ",
                    AbstractValuesProperty::UNKNOWN );
            }

            showAllMetadata( &property );
            const unsigned int el_per_val = property.getElementCountPerValue();

            const unsigned int val_count =
                property.getValuesCountOfPatch( patch_index );

            if( dynamic_cast< ContinuousProperty* >( &property ) != nullptr )
            {
                ContinuousProperty& continuousProp =
                    static_cast< ContinuousProperty& >( property );

                if( ( property.getValuesHdfDatatype()
                            == AbstractValuesProperty::FLOAT
                        || property.getValuesHdfDatatype()
                               == AbstractValuesProperty::DOUBLE ) )
                {
                    GEO::Attribute< double > attribute;
                    attribute.create_vector_attribute(
                        attri_manager, property.getTitle(), el_per_val );

                    std::unique_ptr< double[] > values( new double[val_count] );
                    continuousProp.getDoubleValuesOfPatch(
                        patch_index, &values[0] );
                    for( auto val : range( val_count ) )
                    {
                        attribute[val] = values[val];
                    }
                }
            }
            else if( dynamic_cast< DiscreteProperty* >( &property ) != nullptr )
            {
                DiscreteProperty* discreteProp =
                    dynamic_cast< DiscreteProperty* >( &property );

                GEO::Attribute< int > attribute;
                attribute.create_vector_attribute(
                    attri_manager, property.getTitle(), el_per_val );

                std::unique_ptr< int[] > values( new int[val_count] );
                discreteProp->getIntValuesOfPatch( patch_index, &values[0] );
                for( auto val : range( val_count ) )
                {
                    attribute[val] = values[val];
                }
            }

            return true;
        }
    } // namespace

    bool GeoModelBuilderRESQMLImpl::read_surfaces( const EpcDocument& pck )
    {
        std::vector< TriangulatedSetRepresentation* > all_tri_set_rep =
            pck.getAllTriangulatedSetRepSet();
        ringmesh_assert( !all_tri_set_rep.empty() );

        std::cout << std::endl
                  << "ALL TRI REP: " << all_tri_set_rep.size() << std::endl;

        for( auto tri_set : all_tri_set_rep )
        {
            showAllMetadata( tri_set );

            const gmge_id interface_id =
                builder_.geology.create_geological_entity(
                    Interface3D::type_name_static() );

            AbstractFeature* feature = nullptr;
            if( tri_set->getInterpretation() != nullptr )
            {
                feature = tri_set->getInterpretation()->getInterpretedFeature();
            }

            if( feature != nullptr )
            {
                if( dynamic_cast< TectonicBoundaryFeature* >( feature )
                    != nullptr )
                {
                    builder_.geology.set_geological_entity_geol_feature(
                        interface_id, GMGE::GEOL_FEATURE::FAULT );
                }
                else if( dynamic_cast< FrontierFeature* >( feature )
                         != nullptr )
                {
                    builder_.geology.set_geological_entity_geol_feature(
                        interface_id, GMGE::GEOL_FEATURE::VOI );
                }
                else if( dynamic_cast< GeneticBoundaryFeature* >( feature )
                         != nullptr )
                {
                    builder_.geology.set_geological_entity_geol_feature(
                        interface_id, GMGE::GEOL_FEATURE::STRATI );
                }

                interp_2_geo_entity_[tri_set->getInterpretation()] =
                    interface_id;
                geo_entity_2_interp_[interface_id] =
                    tri_set->getInterpretation();
            }

            ULONG64 global_point_count = 0;
            for( auto patch : range( tri_set->getPatchCount() ) )
            {
                const ULONG64 pointCount =
                    tri_set->getXyzPointCountOfPatch( patch );

                std::cout << "point Count " << pointCount << std::endl;

                std::cout << "TRI REP GEOMETRY" << std::endl;
                std::unique_ptr< double[] > xyzPoints(
                    new double[pointCount * 3] );
                tri_set->getXyzPointsOfPatch( patch, &xyzPoints[0] );

                std::vector< vec3 > points( pointCount, vec3() );
                for( auto i : range( pointCount ) )
                {
                    points[i] = vec3( xyzPoints[i * 3], xyzPoints[i * 3 + 1],
                        xyzPoints[i * 3 + 2] );
                }

                const unsigned int triangleCount =
                    tri_set->getTriangleCountOfPatch( patch );
                std::cout << "triangle Count " << triangleCount << std::endl;

                std::vector< index_t > trgls( triangleCount * 3, 0 );

                tri_set->getTriangleNodeIndicesOfPatch( patch, &trgls[0] );
                for( auto& node : trgls )
                {
                    node -= (index_t) global_point_count;
                }

                std::vector< index_t > trgls_ptr( triangleCount + 1, 0 );
                for( auto i : range( trgls_ptr.size() ) )
                {
                    trgls_ptr[i] = i * 3;
                }

                const gmme_id children = builder_.topology.create_mesh_entity(
                    Surface3D::type_name_static() );

                builder_.geology.add_parent_children_relation(
                    interface_id, children );

                builder_.geometry.set_surface_geometry(
                    children.index(), points, trgls, trgls_ptr );

                for( auto prop_index :
                    range( tri_set->getValuesPropertySet().size() ) )
                {
                    AbstractValuesProperty* propVal =
                        tri_set->getValuesPropertySet()[prop_index];
                    const Surface3D& cur_surf =
                        geomodel_.surface( children.index() );

                    const gsoap_resqml2_0_1::resqml2__IndexableElements
                        attachment = propVal->getAttachmentKind();
                    if( attachment
                        == gsoap_resqml2_0_1::
                               resqml2__IndexableElements__nodes )
                    {
                        read_property( *propVal, patch,
                            cur_surf.vertex_attribute_manager() );
                    }
                    else if( attachment
                             == gsoap_resqml2_0_1::
                                    resqml2__IndexableElements__cells )
                    {
                        read_property( *propVal, patch,
                            cur_surf.polygon_attribute_manager() );
                    }
                }

                global_point_count += pointCount;
            }
        }

        builder_.build_lines_and_corners_from_surfaces();
        builder_.build_regions_from_lines_and_surfaces();
        builder_.geology.build_contacts();
        return true;
    }

    namespace
    {
        bool read_volume_rep( VolumeMesh3D& mesh,
            UnstructuredGridRepresentation& unstructed_grid )
        {
            auto mesh_builder = VolumeMeshBuilder3D::create_builder( mesh );

            unstructed_grid.loadGeometry();

            std::unique_ptr< double[] > gridPoints(
                new double[unstructed_grid.getXyzPointCountOfPatch( 0 ) * 3] );
            unstructed_grid.getXyzPointsOfAllPatchesInGlobalCrs(
                &gridPoints[0] );

            const ULONG64 nb_vertices =
                unstructed_grid.getXyzPointCountOfPatch( 0 );

            for( auto v : range( nb_vertices ) )
            {
                const vec3 vertex( gridPoints[v * 3], gridPoints[v * 3 + 1],
                    gridPoints[v * 3 + 2] );

                mesh_builder->create_vertex( vertex );
            }

            const ULONG64 nb_cells = unstructed_grid.getCellCount();
            mesh_builder->create_cells(
                (index_t) nb_cells, CellType::TETRAHEDRON );

            std::unique_ptr< ULONG64[] > faceCountOfCells(
                new ULONG64[nb_cells] );
            unstructed_grid.getCumulativeFaceCountPerCell(
                &faceCountOfCells[0] );

            const ULONG64 faceCount = faceCountOfCells[nb_cells - 1];

            std::unique_ptr< ULONG64[] > nodeCountOfFaces(
                new ULONG64[faceCount] );
            unstructed_grid.getCumulativeNodeCountPerFace(
                &nodeCountOfFaces[0] );

            for( auto cell : range( nb_cells ) )
            {
                ULONG64 end_face = faceCountOfCells[cell];
                ULONG64 start_face =
                    ( cell == 0 ) ? 0 : faceCountOfCells[cell - 1];

                std::vector< index_t > vertices = {
                    (index_t) unstructed_grid.getNodeIndicesOfFaceOfCell(
                        cell, 0 )[0],
                    (index_t) unstructed_grid.getNodeIndicesOfFaceOfCell(
                        cell, 0 )[1],
                    (index_t) unstructed_grid.getNodeIndicesOfFaceOfCell(
                        cell, 0 )[2],
                    0
                };

                bool found = false;
                for( unsigned int f = 1; f < ( end_face - start_face ); ++f )
                {
                    ULONG64 nb_nodes =
                        unstructed_grid.getNodeCountOfFaceOfCell( cell, f );

                    for( ULONG64 node = 0; node < nb_nodes; ++node )
                    {
                        const ULONG64 node_index =
                            unstructed_grid.getNodeIndicesOfFaceOfCell(
                                cell, f )[node];
                        if( node_index != vertices[0]
                            && node_index != vertices[1]
                            && node_index != vertices[2] )
                        {
                            vertices[3] = (index_t) node_index;
                            found = true;
                            break;
                        }

                        if( found )
                        {
                            break;
                        }
                    }
                }
                ringmesh_assert( found );

                for( auto v_id : range( 4 ) )
                {
                    mesh_builder->set_cell_vertex(
                        { cell, v_id }, vertices[v_id] );
                }
            }
            unstructed_grid.unloadGeometry();
            mesh_builder->connect_cells();

            return true;
        }
    } // namespace

    bool GeoModelBuilderRESQMLImpl::read_volumes( const EpcDocument& pck )
    {
        std::vector< UnstructuredGridRepresentation* > unstructuredGridRepSet =
            pck.getUnstructuredGridRepresentationSet();
        if( unstructuredGridRepSet.empty() )
        {
            return true;
        }

        std::cout << std::endl
                  << "UNSTRUCTURED GRID REP: " << unstructuredGridRepSet.size()
                  << std::endl;

        for( auto unstructured_grid : unstructuredGridRepSet )
        {
            showAllMetadata( unstructured_grid );

            if( !unstructured_grid->isPartial()
                && unstructured_grid->hasGeometry() )
            {
                auto mesh = VolumeMesh3D::create_mesh();
                bool result = read_volume_rep( *mesh, *unstructured_grid );
                ringmesh_assert( result );

                // the volume mesh from resqml is here, need to find the
                // corresponding region of in the GeoModel3D
                const auto& nn_search = mesh->cell_facet_nn_search();

                int region_index = -1;
                for( auto r : range( geomodel_.nb_regions() ) )
                {
                    const Region3D& region = geomodel_.region( r );

                    bool match = true;
                    for( auto b : range( region.nb_boundaries() ) )
                    {
                        const Surface3D& surface = region.boundary( b );
                        for( auto p : range( surface.nb_mesh_elements() ) )
                        {
                            auto center = surface.mesh_element_barycenter( p );
                            auto result = nn_search.get_neighbors(
                                center, surface.geomodel().epsilon() );
                            if( result.empty() )
                            {
                                match = false;
                                break;
                            }
                        }
                        if( !match )
                        {
                            break;
                        }
                    }
                    if( match )
                    {
                        region_index = (int) r;
                        break;
                    }
                }

                if( region_index == -1 )
                {
                    return false;
                }

                // corresponding region found, build its volume mesh
                const gmme_id region_id(
                    region_type_name_static(), (index_t) region_index );

                auto mesh_builder = builder_.geometry.create_region_builder(
                    region_id.index() );

                for( auto v : range( mesh->nb_vertices() ) )
                {
                    mesh_builder->create_vertex( mesh->vertex( v ) );
                }

                mesh_builder->create_cells(
                    mesh->nb_cells(), CellType::TETRAHEDRON );

                for( auto cell : range( mesh->nb_cells() ) )
                {
                    const index_t nb_vertices = mesh->nb_cell_vertices( cell );
                    std::vector< index_t > cell_vertices( nb_vertices, 0 );
                    for( auto v : range( nb_vertices ) )
                    {
                        ElementLocalVertex lv( cell, v );
                        mesh_builder->set_cell_vertex(
                            lv, mesh->cell_vertex( lv ) );
                    }
                }

                mesh_builder->connect_cells();

                // property
                const Region3D& cur_reg = geomodel_.region( region_id.index() );

                for( auto prop_index :
                    range( unstructured_grid->getValuesPropertySet().size() ) )
                {
                    AbstractValuesProperty* propVal =
                        unstructured_grid->getValuesPropertySet()[prop_index];
                    // TODO there is only one patch in a unstructured_grid

                    const gsoap_resqml2_0_1::resqml2__IndexableElements
                        attachment = propVal->getAttachmentKind();
                    if( attachment
                        == gsoap_resqml2_0_1::
                               resqml2__IndexableElements__nodes )
                    {
                        read_property(
                            *propVal, 0, cur_reg.vertex_attribute_manager() );
                    }
                    else if( attachment
                             == gsoap_resqml2_0_1::
                                    resqml2__IndexableElements__cells )
                    {
                        read_property(
                            *propVal, 0, cur_reg.cell_attribute_manager() );
                    }
                }
            }
        }

        return true;
    }

    bool GeoModelBuilderRESQMLImpl::load_file()
    {
        EpcDocument pck( builder_.filename(), fesapi_resources_path,
            EpcDocument::READ_ONLY );

        deserialize( pck );

        read_surfaces( pck );

        read_volumes( pck );

        read_strati_column( pck );

        return true;
    }

    /****************************************************************************/

    GeoModelBuilderRESQML::GeoModelBuilderRESQML(
        GeoModel3D& geomodel, const std::string& filename )
        : GeoModelBuilderFile( geomodel, std::move( filename ) ),
          impl_( new GeoModelBuilderRESQMLImpl( *this, this->geomodel_ ) )
    {
    }

    GeoModelBuilderRESQML::~GeoModelBuilderRESQML()
    {
        // needed due to the unique_ptr impl_
    }

    void GeoModelBuilderRESQML::load_file()
    {
        bool result = impl_->load_file();
        ringmesh_assert( result );
    }

} // namespace RINGMesh
