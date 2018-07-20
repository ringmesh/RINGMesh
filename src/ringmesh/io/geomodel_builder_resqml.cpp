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
#include <geogram/basic/logger.h>

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
        bool build_fake_geomodel();

    private:
        GeoModelBuilderRESQML& builder_;
        GeoModel3D& geomodel_;

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
            // TODO the following detection is not working if the region doesn't
            // contain a single face of one or even both horizons
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

                if( cur_horizons_list.size() != horizons_list.size() )
                {
                    continue;
                }
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
        return gmge_id();
    }

    void GeoModelBuilderRESQMLImpl::read_strati_column( EpcDocument& pck )
    {
        const std::vector< RESQML2_0_1_NS::StratigraphicColumn* >
            stratiColumnSet = pck.getStratigraphicColumnSet();
        // TODO only one column for the moment
        ringmesh_assert( stratiColumnSet.size() <= 1 );

        for( auto column : stratiColumnSet )
        {
            showAllMetadata( column );

            // TODO now we read only rank 0
            ringmesh_assert(
                column->getStratigraphicColumnRankInterpretationSet().size()
                == 1 );

            for( auto rank_interp :
                column->getStratigraphicColumnRankInterpretationSet() )
            {
                showAllMetadata( rank_interp );
                // TODO only chronostrati
                ringmesh_assert( rank_interp->isAChronoStratiRank() );

                for( auto unit :
                    rank_interp->getStratigraphicUnitInterpretationSet() )
                {
                    const gmge_id layer =
                        builder_.geology.create_geological_entity(
                            Layer3D::type_name_static() );

                    builder_.geology.set_geological_entity_geol_feature(
                        layer, GMGE::GEOL_FEATURE::STRATI_UNIT );

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

                // assign regions to layers
                // TODO this is a workaround until the
                // SealedVolumeFrameworkRepresentation
                // is released in fesapi that tells us the layer information of
                // a region
                for( auto& region : geomodel_.regions() )
                {
                    gmge_id layer = find_layer( region );
                    ringmesh_assert( layer.is_defined() );
                    builder_.geology.add_parent_children_relation(
                        layer, region.gmme() );
                }

                // build the stati column
                RockFeature rock( "rock", ROCKTYPE::NONE );
                std::vector< std::shared_ptr< const StratigraphicUnit > > units;

                for( auto unit_key :
                    rank_interp->getStratigraphicUnitInterpretationSet() )
                {
                    const UnitInfo& unit = unit_2_info_[unit_key];
                    const Interface3D* top =
                        unit.interface_top_.is_defined()
                            ? dynamic_cast< const Interface3D* >(
                                  &geomodel_.geological_entity(
                                      Interface3D::type_name_static(),
                                      unit.interface_top_.index() ) )
                            : nullptr;
                    const Interface3D* base =
                        unit.interface_base_.is_defined()
                            ? dynamic_cast< const Interface3D* >(
                                  &geomodel_.geological_entity(
                                      Interface3D::type_name_static(),
                                      unit.interface_base_.index() ) )
                            : nullptr;

                    std::shared_ptr< const StratigraphicUnit > sunit(
                        new UnsubdividedStratigraphicUnit( unit.name_, base,
                            top,
                            dynamic_cast< const Layer3D& >(
                                geomodel_.geological_entity(
                                    Layer3D::type_name_static(),
                                    unit.layer_.index() ) ),
                            unit.relation_top_, unit.relation_base_, rock, 0,
                            10 ) );

                    units.push_back( sunit );

#ifdef RINGMESH_DEBUG
                    std::cout << " relation base: " << (int) unit.relation_base_
                              << std::endl;
                    std::cout << " relation top: " << (int) unit.relation_top_
                              << std::endl;
                    std::cout
                        << " interface base: " << unit.interface_base_.index()
                        << std::endl;
                    std::cout
                        << " interface top: " << unit.interface_top_.index()
                        << std::endl;
#endif
                }

                StratigraphicColumn* strati =
                    new StratigraphicColumn( column->getTitle(), units,
                        STRATIGRAPHIC_PARADIGM::CHRONOSTRATIGRAPHIC );
                geomodel_.set_stratigraphic_column( strati );
            }
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

    bool GeoModelBuilderRESQMLImpl::build_fake_geomodel()
    {
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
            builder_.topology.create_mesh_entity(
                Corner3D::type_name_static() );
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
            builder_.topology.create_mesh_entity(
                Region3D::type_name_static() );
        }

        //#############################
        // Setting the Geometry       #
        //#############################

        // We declare the coordinates of the corners. We arrange the corner
        // in a table
        vec3 corners_table[4];
        corners_table[0] = vec3( 0, 0, 300 );
        corners_table[1] = vec3( 700, 0, 350 );
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
        // The add_mesh_entity_boundary_relation method take as first
        // argument the gme_t of the boundary and in second argument the id
        // of the GeoModelMeshentity bounded by the boundary Remember :
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
        // Second argument is the index of the child (ie a
        // GeoModelMeshEntity)

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

    bool GeoModelBuilderRESQMLImpl::read_surfaces( const EpcDocument& pck )
    {
        std::vector< TriangulatedSetRepresentation* > all_tri_set_rep =
            pck.getAllTriangulatedSetRepSet();

        if( all_tri_set_rep.empty() )
        {
            return build_fake_geomodel();
        }

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
            }

            ULONG64 global_point_count = 0;
            for( auto patch : range( tri_set->getPatchCount() ) )
            {
                const ULONG64 pointCount =
                    tri_set->getXyzPointCountOfPatch( patch );

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
            for( auto c : range( nb_cells ) )
            {
                const index_t nb_faces =
                    (index_t) unstructed_grid.getFaceCountOfCell( c );
                if( nb_faces == 4 )
                {
                    const index_t cell = mesh_builder->create_cells(
                        (index_t) 1, CellType::TETRAHEDRON );

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
                    for( unsigned int f = 1; f < nb_faces; ++f )
                    {
                        const ULONG64 nb_nodes =
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
                else if( nb_faces == 5 )
                {
                    const index_t cell = mesh_builder->create_cells(
                        (index_t) 1, CellType::PYRAMID );

                    std::vector< index_t > vertices = {
                        (index_t) unstructed_grid.getNodeIndicesOfFaceOfCell(
                            cell, 0 )[1],
                        (index_t) unstructed_grid.getNodeIndicesOfFaceOfCell(
                            cell, 0 )[0],
                        (index_t) unstructed_grid.getNodeIndicesOfFaceOfCell(
                            cell, 0 )[3],
                        (index_t) unstructed_grid.getNodeIndicesOfFaceOfCell(
                            cell, 0 )[2],
                        (index_t) unstructed_grid.getNodeIndicesOfFaceOfCell(
                            cell, 1 )[0]
                    };

                    for( auto v_id : range( 5 ) )
                    {
                        mesh_builder->set_cell_vertex(
                            { cell, v_id }, vertices[v_id] );
                    }
                }
                else if( nb_faces == 6 )
                {
                    const index_t cell = mesh_builder->create_cells(
                        (index_t) 1, CellType::HEXAHEDRON );

                    std::vector< index_t > vertices = {
                        (index_t) unstructed_grid.getNodeIndicesOfFaceOfCell(
                            cell, 0 )[0],
                        (index_t) unstructed_grid.getNodeIndicesOfFaceOfCell(
                            cell, 0 )[1],
                        (index_t) unstructed_grid.getNodeIndicesOfFaceOfCell(
                            cell, 0 )[3],
                        (index_t) unstructed_grid.getNodeIndicesOfFaceOfCell(
                            cell, 0 )[2],
                        (index_t) unstructed_grid.getNodeIndicesOfFaceOfCell(
                            cell, 1 )[0],
                        (index_t) unstructed_grid.getNodeIndicesOfFaceOfCell(
                            cell, 1 )[1],
                        (index_t) unstructed_grid.getNodeIndicesOfFaceOfCell(
                            cell, 1 )[3],
                        (index_t) unstructed_grid.getNodeIndicesOfFaceOfCell(
                            cell, 1 )[2]
                    };

                    for( auto v_id : range( 8 ) )
                    {
                        mesh_builder->set_cell_vertex(
                            { cell, v_id }, vertices[v_id] );
                    }
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

        for( auto unstructured_grid : unstructuredGridRepSet )
        {
            showAllMetadata( unstructured_grid );

            if( unstructured_grid->isPartial()
                || !unstructured_grid->hasGeometry() )
            {
                Logger::err(
                    "Attempted to read partial or empty UnstructuredGrid" );
            }

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
                if( !geomodel_.entity_type_manager()
                         .mesh_entity_manager.is_valid_type(
                             Region3D::type_name_static() ) )
                {
                    return false;
                }
                else
                {
                    region_index = 0;
                }
            }

            // corresponding region found, build its volume mesh
            const gmme_id region_id(
                region_type_name_static(), (index_t) region_index );

            auto mesh_builder =
                builder_.geometry.create_region_builder( region_id.index() );

            for( auto v : range( mesh->nb_vertices() ) )
            {
                mesh_builder->create_vertex( mesh->vertex( v ) );
            }

            for( auto cell : range( mesh->nb_cells() ) )
            {
                mesh_builder->create_cells(
                    (index_t) 1, mesh->cell_type( cell ) );

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
                // TODO there is only one patch in a UnstructuredGrid

                const gsoap_resqml2_0_1::resqml2__IndexableElements attachment =
                    propVal->getAttachmentKind();
                if( attachment
                    == gsoap_resqml2_0_1::resqml2__IndexableElements__nodes )
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
