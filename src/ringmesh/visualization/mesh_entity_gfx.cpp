/*
 * Copyright (c) 2012-2017, Association Scientifique pour la Geologie et ses Applications (ASGA)
 * All rights reserved.
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
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL ASGA BE LIABLE FOR ANY
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 *     http://www.ring-team.org
 *
 *     RING Project
 *     Ecole Nationale Superieure de Geologie - GeoRessources
 *     2 Rue du Doyen Marcel Roubault - TSA 70605
 *     54518 VANDOEUVRE-LES-NANCY
 *     FRANCE
 */

/*! 
 * @file Implementation of visualization of GeoModelEntities
 * @author Benjamin Chauvin and Arnaud Botella
 */

#include <ringmesh/visualization/mesh_entity_gfx.h>

#ifdef RINGMESH_WITH_GRAPHICS

#include <ringmesh/geomodel/geomodel.h>
#include <ringmesh/geomodel/geomodel_entity.h>
#include <ringmesh/geomodel/geomodel_mesh_entity.h>
#include <ringmesh/geomodel/geomodel_geological_entity.h>

#include <ringmesh/visualization/geomodel_gfx.h>
#include <ringmesh/visualization/geogram_gfx.h>

namespace {
    using namespace RINGMesh;

    std::string get_attribute_name_with_coordinate(
        const std::string& name,
        index_t coordinate )
    {
        return name + "[" + GEO::String::to_string( coordinate ) + "]";
    }

    void compute_attribute_range(
        GEO::ReadOnlyScalarAttributeAdapter& attribute,
        double& min,
        double& max )
    {
        if( attribute.is_bound() ) {
            for( index_t i = 0; i < attribute.size(); ++i ) {
                double value = attribute[i];
                min = GEO::geo_min( min, value );
                max = GEO::geo_max( max, value );
            }
        }
    }
}

namespace RINGMesh {

    class CellAttributeGfx: public AttributeGfx {
    public:
        CellAttributeGfx( AttributeGfxManager& manager )
            : AttributeGfx( manager )
        {
        }

        virtual std::string location_name() const override
        {
            return "cells";
        }
        virtual void bind_attribute() override
        {
            std::string attribute_name = get_attribute_name_with_coordinate(
                manager_.name(), manager_.coordinate() );
            manager_.gfx().regions.set_scalar_attribute( GEO::MESH_CELLS,
                attribute_name, manager_.minimum(), manager_.maximum(),
                manager_.colormap() );
        }
        virtual void unbind_attribute() override
        {
            manager_.gfx().regions.unset_scalar_attribute();
        }
        virtual index_t nb_coordinates() override
        {
            const GeoModel< 3 >* geomodel = manager_.gfx().geomodel();
            GEO::AttributeStore* store =
                geomodel->region( 0 ).cell_attribute_manager().find_attribute_store(
                    manager_.name() );

            if( store == nullptr ) return 0;
            return store->dimension();
        }
    private:
        virtual void do_compute_range( double& attribute_min, double& attribute_max ) override
        {
            std::string attribute_name = get_attribute_name_with_coordinate(
                manager_.name(), manager_.coordinate() );
            const GeoModel< 3 >* geomodel = manager_.gfx().geomodel();
            for( index_t r = 0; r < geomodel->nb_regions(); r++ ) {
                GEO::ReadOnlyScalarAttributeAdapter attribute(
                    geomodel->region( r ).cell_attribute_manager(), attribute_name );
                compute_attribute_range( attribute, attribute_min, attribute_max );
            }
        }
    };

    class CellVertexAttributeGfx: public AttributeGfx {
    public:
        CellVertexAttributeGfx( AttributeGfxManager& manager )
            : AttributeGfx( manager )
        {
        }

        virtual std::string location_name() const override
        {
            return "cell_vertices";
        }
        virtual void bind_attribute() override
        {
            std::string attribute_name = get_attribute_name_with_coordinate(
                manager_.name(), manager_.coordinate() );
            manager_.gfx().regions.set_scalar_attribute( GEO::MESH_VERTICES,
                attribute_name, manager_.minimum(), manager_.maximum(),
                manager_.colormap() );
        }
        virtual void unbind_attribute() override
        {
            manager_.gfx().regions.unset_scalar_attribute();
        }
        virtual index_t nb_coordinates() override
        {
            const GeoModel< 3 >* geomodel = manager_.gfx().geomodel();
            GEO::AttributeStore* store =
                geomodel->region( 0 ).vertex_attribute_manager().find_attribute_store(
                    manager_.name() );

            if( store == nullptr ) return 0;
            return store->dimension();
        }
    private:
        virtual void do_compute_range( double& attribute_min, double& attribute_max ) override
        {
            std::string attribute_name = get_attribute_name_with_coordinate(
                manager_.name(), manager_.coordinate() );
            const GeoModel< 3 >* geomodel = manager_.gfx().geomodel();
            for( index_t r = 0; r < geomodel->nb_regions(); r++ ) {
                GEO::ReadOnlyScalarAttributeAdapter attribute(
                    geomodel->region( r ).vertex_attribute_manager(),
                    attribute_name );
                compute_attribute_range( attribute, attribute_min, attribute_max );
            }
        }
    };

    class PolygonAttributeGfx: public AttributeGfx {
    public:
        PolygonAttributeGfx( AttributeGfxManager& manager )
            : AttributeGfx( manager )
        {
        }

        virtual std::string location_name() const override
        {
            return "polygon";
        }
        virtual void bind_attribute() override
        {
            std::string attribute_name = get_attribute_name_with_coordinate(
                manager_.name(), manager_.coordinate() );
            manager_.gfx().surfaces.set_scalar_attribute( GEO::MESH_FACETS,
                attribute_name, manager_.minimum(), manager_.maximum(),
                manager_.colormap() );
        }
        virtual void unbind_attribute() override
        {
            manager_.gfx().surfaces.unset_scalar_attribute();
        }
        virtual index_t nb_coordinates() override
        {
            const GeoModel< 3 >* geomodel = manager_.gfx().geomodel();
            GEO::AttributeStore* store =
                geomodel->surface( 0 ).polygon_attribute_manager().find_attribute_store(
                    manager_.name() );

            if( store == nullptr ) return 0;
            return store->dimension();
        }
    private:
        virtual void do_compute_range( double& attribute_min, double& attribute_max ) override
        {
            std::string attribute_name = get_attribute_name_with_coordinate(
                manager_.name(), manager_.coordinate() );
            const GeoModel< 3 >* geomodel = manager_.gfx().geomodel();
            for( index_t s = 0; s < geomodel->nb_surfaces(); s++ ) {
                GEO::ReadOnlyScalarAttributeAdapter attribute(
                    geomodel->surface( s ).polygon_attribute_manager(),
                    attribute_name );
                compute_attribute_range( attribute, attribute_min, attribute_max );
            }
        }
    };

    class PolygonVertexAttributeGfx: public AttributeGfx {
    public:
        PolygonVertexAttributeGfx( AttributeGfxManager& manager )
            : AttributeGfx( manager )
        {
        }

        virtual std::string location_name() const override
        {
            return "polygon_vertices";
        }
        virtual void bind_attribute() override
        {
            std::string attribute_name = get_attribute_name_with_coordinate(
                manager_.name(), manager_.coordinate() );
            manager_.gfx().surfaces.set_scalar_attribute( GEO::MESH_VERTICES,
                attribute_name, manager_.minimum(), manager_.maximum(),
                manager_.colormap() );
        }
        virtual void unbind_attribute() override
        {
            manager_.gfx().surfaces.unset_scalar_attribute();
        }
        virtual index_t nb_coordinates() override
        {
            const GeoModel< 3 >* geomodel = manager_.gfx().geomodel();
            GEO::AttributeStore* store =
                geomodel->surface( 0 ).vertex_attribute_manager().find_attribute_store(
                    manager_.name() );

            if( store == nullptr ) return 0;
            return store->dimension();
        }
    private:
        virtual void do_compute_range( double& attribute_min, double& attribute_max ) override
        {
            std::string attribute_name = get_attribute_name_with_coordinate(
                manager_.name(), manager_.coordinate() );
            const GeoModel< 3 >* geomodel = manager_.gfx().geomodel();
            for( index_t s = 0; s < geomodel->nb_surfaces(); s++ ) {
                GEO::ReadOnlyScalarAttributeAdapter attribute(
                    geomodel->surface( s ).vertex_attribute_manager(),
                    attribute_name );
                compute_attribute_range( attribute, attribute_min, attribute_max );
            }
        }
    };

    AttributeGfxManager::AttributeGfxManager( GeoModelGfx& gfx )
        :
            gfx_( gfx ),
            location_( nb_locations ),
            coordinate_( 0 ),
            colormap_texture_( 0 ),
            minimum_( 0.0 ),
            maximum_( 0.0 )
    {
        attributes_[polygons].reset( new PolygonAttributeGfx( *this ) );
        attributes_[polygon_vertices].reset(
            new PolygonVertexAttributeGfx( *this ) );
        attributes_[cells].reset( new CellAttributeGfx( *this ) );
        attributes_[cell_vertices].reset( new CellVertexAttributeGfx( *this ) );
    }

    std::string AttributeGfxManager::location_name( Attribute_location location )
    {
        if( location == nb_locations ) {
            return "location";
        } else {
            return attributes_[location]->location_name();
        }
    }

    void AttributeGfxManager::compute_range()
    {
        if( location() < nb_locations ) {
            attributes_[location()]->compute_range();
        }
    }

    void AttributeGfxManager::bind_attribute()
    {
        if( location() < nb_locations ) {
            attributes_[location()]->bind_attribute();
        }
    }

    void AttributeGfxManager::unbind_attribute()
    {
        if( location() < nb_locations ) {
            attributes_[location()]->unbind_attribute();
        }
    }

    index_t AttributeGfxManager::nb_coordinates() const
    {
        if( location() < nb_locations ) {
            return attributes_[location()]->nb_coordinates();
        }
        return 0;
    }

    std::unique_ptr< PointSetMeshGfx > PointSetMeshGfx::create_gfx(
        const PointSetMesh< 3 >& mesh )
    {
        PointSetMeshGfx* gfx = PointSetMeshGfxFactory::create_object(
            mesh.type_name() );
        if( !gfx ) {
            Logger::warn( "PointSetMeshGfx",
                "Could not create mesh data structure: ", mesh.type_name() );
            Logger::warn( "PointSetMeshGfx",
                "Falling back to GeogramPointSetMeshGfx data structure" );

            gfx = new GeogramPointSetMeshGfx;
        }
        gfx->set_mesh( mesh );
        return std::unique_ptr< PointSetMeshGfx >( gfx );
    }

    std::unique_ptr< LineMeshGfx > LineMeshGfx::create_gfx(
        const LineMesh< 3 >& mesh )
    {
        LineMeshGfx* gfx = LineMeshGfxFactory::create_object( mesh.type_name() );
        if( !gfx ) {
            Logger::warn( "LineMeshGfx", "Could not create mesh data structure: ",
                mesh.type_name() );
            Logger::warn( "LineMeshGfx",
                "Falling back to GeogramLineMeshGfx data structure" );

            gfx = new GeogramLineMeshGfx;
        }
        gfx->set_mesh( mesh );
        return std::unique_ptr< LineMeshGfx >( gfx );
    }

    std::unique_ptr< SurfaceMeshGfx > SurfaceMeshGfx::create_gfx(
        const SurfaceMesh< 3 >& mesh )
    {
        SurfaceMeshGfx* gfx = SurfaceMeshGfxFactory::create_object(
            mesh.type_name() );
        if( !gfx ) {
            Logger::warn( "SurfaceMeshGfx", "Could not create mesh data structure: ",
                mesh.type_name() );
            Logger::warn( "SurfaceMeshGfx",
                "Falling back to GeogramSurfaceMeshGfx data structure" );

            gfx = new GeogramSurfaceMeshGfx;
        }
        gfx->set_mesh( mesh );
        return std::unique_ptr< SurfaceMeshGfx >( gfx );
    }

    std::unique_ptr< VolumeMeshGfx > VolumeMeshGfx::create_gfx(
        const VolumeMesh< 3 >& mesh )
    {
        VolumeMeshGfx* gfx = VolumeMeshGfxFactory::create_object( mesh.type_name() );
        if( !gfx ) {
            Logger::warn( "VolumeMeshGfx", "Could not create mesh data structure: ",
                mesh.type_name() );
            Logger::warn( "VolumeMeshGfx",
                "Falling back to GeogramVolumeMeshGfx data structure" );

            gfx = new GeogramVolumeMeshGfx;
        }
        gfx->set_mesh( mesh );
        return std::unique_ptr< VolumeMeshGfx >( gfx );
    }
}

#endif
