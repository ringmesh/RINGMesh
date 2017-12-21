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

/*!
 * @file Implementation of the GeoModel
 * @author Jeanne Pellerin and Arnaud Botella
 */

#include <geogram/basic/command_line.h>

#include <ringmesh/basic/box.h>

#include <ringmesh/geomodel/core/geomodel.h>
#include <ringmesh/geomodel/core/geomodel_geological_entity.h>
#include <ringmesh/geomodel/core/geomodel_mesh_entity.h>

namespace
{
    using namespace RINGMesh;

    template < index_t DIMENSION >
    void compute_mesh_entity_bbox(
        const GeoModelMeshEntity< DIMENSION >& entity, Box< DIMENSION >& bbox )
    {
        for( auto v : range( entity.nb_vertices() ) )
        {
            bbox.add_point( entity.vertex( v ) );
        }
    }

    template < index_t DIMENSION >
    double compute_percentage_bbox_diagonal(
        const GeoModelBase< DIMENSION >& gm )
    {
        Box< DIMENSION > bbox;
        if( gm.nb_surfaces() > 0 )
        {
            for( auto s : range( gm.nb_surfaces() ) )
            {
                compute_mesh_entity_bbox( gm.surface( s ), bbox );
            }
        }
        else if( gm.nb_lines() > 0 )
        {
            for( auto l : range( gm.nb_lines() ) )
            {
                compute_mesh_entity_bbox( gm.line( l ), bbox );
            }
        }
        else
        {
            ringmesh_assert( gm.nb_corners() > 0 );
            for( auto c : range( gm.nb_corners() ) )
            {
                bbox.add_point( gm.corner( c ).vertex( 0 ) );
            }
        }
        return bbox.diagonal().length()
               * GEO::CmdLine::get_arg_double( "epsilon" );
    }
} // namespace

namespace RINGMesh
{
    template < index_t DIMENSION >
    GeoModelBase< DIMENSION >::GeoModelBase( GeoModel< DIMENSION >& geomodel )
        : mesh( geomodel )
    {
    }

    template < index_t DIMENSION >
    GeoModelBase< DIMENSION >::~GeoModelBase()
    {
    }

    template < index_t DIMENSION >
    index_t GeoModelBase< DIMENSION >::nb_mesh_entities(
        const MeshEntityType& type ) const
    {
        const MeshEntityTypeManager< DIMENSION >& manager =
            entity_type_manager().mesh_entity_manager;
        if( manager.is_line( type ) )
        {
            return nb_lines();
        }
        if( manager.is_corner( type ) )
        {
            return nb_corners();
        }
        if( manager.is_surface( type ) )
        {
            return nb_surfaces();
        }
        ringmesh_assert_not_reached;
        return 0;
    }

    template < index_t DIMENSION >
    const GeoModelMeshEntity< DIMENSION >&
        GeoModelBase< DIMENSION >::mesh_entity( const gmme_id& id ) const
    {
        const MeshEntityTypeManager< DIMENSION >& manager =
            entity_type_manager().mesh_entity_manager;
        const MeshEntityType& type = id.type();
        index_t index = id.index();
        if( manager.is_line( type ) )
        {
            return line( index );
        }
        if( manager.is_corner( type ) )
        {
            return corner( index );
        }
        if( manager.is_surface( type ) )
        {
            return surface( index );
        }
        ringmesh_assert_not_reached;
        return surface( 0 );
    }

    template < index_t DIMENSION >
    const std::vector< std::unique_ptr< GeoModelMeshEntity< DIMENSION > > >&
        GeoModelBase< DIMENSION >::mesh_entities(
            const MeshEntityType& type ) const
    {
        const MeshEntityTypeManager< DIMENSION >& manager =
            entity_type_manager().mesh_entity_manager;
        if( manager.is_corner( type ) )
        {
            return corners_;
        }
        if( manager.is_line( type ) )
        {
            return lines_;
        }
        if( manager.is_surface( type ) )
        {
            return surfaces_;
        }
        ringmesh_assert_not_reached;
        return surfaces_;
    }

    template < index_t DIMENSION >
    const Corner< DIMENSION >& GeoModelBase< DIMENSION >::corner(
        index_t index ) const
    {
        ringmesh_assert( index < corners_.size() );
        return *static_cast< const Corner< DIMENSION >* >(
            corners_[index].get() );
    }
    template < index_t DIMENSION >
    const Line< DIMENSION >& GeoModelBase< DIMENSION >::line(
        index_t index ) const
    {
        ringmesh_assert( index < lines_.size() );
        return *static_cast< const Line< DIMENSION >* >( lines_[index].get() );
    }
    template < index_t DIMENSION >
    const Surface< DIMENSION >& GeoModelBase< DIMENSION >::surface(
        index_t index ) const
    {
        ringmesh_assert( index < surfaces_.size() );
        return *static_cast< const Surface< DIMENSION >* >(
            surfaces_[index].get() );
    }

    template < index_t DIMENSION >
    void GeoModelBase< DIMENSION >::set_wells(
        const WellGroup< DIMENSION >* wells )
    {
        wells_ = wells;
    }

    template < index_t DIMENSION >
    double GeoModelBase< DIMENSION >::epsilon() const
    {
        if( epsilon_ == -1 )
        {
            epsilon_ = compute_percentage_bbox_diagonal( *this );
        }
        return epsilon_;
    }

    template < index_t DIMENSION >
    GeoModel< DIMENSION >::GeoModel() : GeoModelBase< DIMENSION >( *this )
    {
    }

    GeoModel< 3 >::GeoModel() : GeoModelBase< 3 >( *this )
    {
    }

    GeoModel< 3 >::~GeoModel()
    {
    }

    GeoModel< 2 >::GeoModel() : GeoModelBase< 2 >( *this )
    {
    }

    GeoModel< 2 >::GeoModel( PlaneReferenceFrame3D plane_reference_frame )
        : GeoModelBase< 2 >( *this ), reference_frame_( plane_reference_frame )
    {
    }

    GeoModel< 2 >::~GeoModel()
    {
    }

    const Region3D& GeoModel< 3 >::region( index_t index ) const
    {
        ringmesh_assert( index < regions_.size() );
        return *static_cast< const Region3D* >( regions_[index].get() );
    }

    const std::vector< std::unique_ptr< GeoModelMeshEntity3D > >&
        GeoModel< 3 >::mesh_entities( const MeshEntityType& type ) const
    {
        if( entity_type_manager().mesh_entity_manager.is_region( type ) )
        {
            return regions_;
        }
        return GeoModelBase3D::mesh_entities( type );
    }

    const GeoModelMeshEntity3D& GeoModel< 3 >::mesh_entity(
        const gmme_id& id ) const
    {
        const MeshEntityType& type = id.type();
        index_t index = id.index();
        if( entity_type_manager().mesh_entity_manager.is_region( type ) )
        {
            return region( index );
        }
        return GeoModelBase3D::mesh_entity( id );
        ringmesh_assert_not_reached;
        return surface( 0 );
    }

    index_t GeoModel< 3 >::nb_mesh_entities( const MeshEntityType& type ) const
    {
        if( entity_type_manager().mesh_entity_manager.is_region( type ) )
        {
            return nb_regions();
        }
        return GeoModelBase3D::nb_mesh_entities( type );
    }

    SurfaceSide GeoModel< 3 >::voi_surfaces() const
    {
        SurfaceSide surface_side;
        auto& voi_surfaces = surface_side.surfaces_;
        voi_surfaces.reserve( nb_surfaces() );
        auto& voi_surface_region_side = surface_side.sides_;
        voi_surface_region_side.reserve( nb_surfaces() );

        for( const auto& cur_surface : surfaces() )
        {
            if( cur_surface.is_on_voi() )
            {
                ringmesh_assert( cur_surface.nb_incident_entities() == 1 );
                voi_surfaces.push_back( cur_surface.index() );
                const auto& incident_region = cur_surface.incident_entity( 0 );

                index_t local_boundary_id{ NO_ID };
                for( auto region_boundary_i :
                    range( incident_region.nb_boundaries() ) )
                {
                    if( incident_region.boundary_gmme( region_boundary_i )
                            .index()
                        == cur_surface.index() )
                    {
                        local_boundary_id = region_boundary_i;
                        break;
                    }
                }
                ringmesh_assert( local_boundary_id != NO_ID );

                voi_surface_region_side.push_back(
                    incident_region.side( local_boundary_id ) );
            }
        }

        return surface_side;
    }

    LineSide GeoModel< 2 >::voi_lines() const
    {
        LineSide line_side;
        auto& voi_lines = line_side.lines_;
        voi_lines.reserve( nb_lines() );
        auto& voi_line_surface_side = line_side.sides_;
        voi_line_surface_side.reserve( nb_lines() );

        for( const auto& cur_line : lines() )
        {
            if( cur_line.is_on_voi() )
            {
                ringmesh_assert( cur_line.nb_incident_entities() == 1 );
                voi_lines.push_back( cur_line.index() );
                const auto& incident_surface = cur_line.incident_entity( 0 );

                index_t local_boundary_id{ NO_ID };
                for( auto surface_boundary_i :
                    range( incident_surface.nb_boundaries() ) )
                {
                    if( incident_surface.boundary_gmme( surface_boundary_i )
                            .index()
                        == cur_line.index() )
                    {
                        local_boundary_id = surface_boundary_i;
                        break;
                    }
                }
                ringmesh_assert( local_boundary_id != NO_ID );

                voi_line_surface_side.push_back(
                    incident_surface.side( local_boundary_id ) );
            }
        }

        return line_side;
    }

    template class geomodel_core_api GeoModel< 2 >;
    template class geomodel_core_api GeoModelBase< 2 >;

    template class geomodel_core_api GeoModelBase< 3 >;

} // namespace RINGMesh
