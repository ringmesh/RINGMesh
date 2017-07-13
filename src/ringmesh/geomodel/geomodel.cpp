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
 * @file Implementation of THE GeoModel
 * @author Jeanne Pellerin and Arnaud Botella 
 */

#include <ringmesh/geomodel/geomodel.h>

#include <geogram/basic/command_line.h>

#include <ringmesh/geomodel/geomodel_mesh_entity.h>
#include <ringmesh/geomodel/geomodel_geological_entity.h>

namespace {
    using namespace RINGMesh;

    void compute_surface_bbox( const GeoModel& gm, index_t surface_id, Box3d& bbox )
    {
        const Surface& surface = gm.surface( surface_id );
        for( index_t v = 0; v < surface.nb_vertices(); v++ ) {
            bbox.add_point( surface.vertex( v ) );
        }
    }

    double compute_percentage_bbox_diagonal( const GeoModel& gm )
    {
        Box3d bbox;
        ringmesh_assert( gm.nb_surfaces() > 0 );
        for( index_t s = 0; s < gm.nb_surfaces(); s++ ) {
            compute_surface_bbox( gm, s, bbox );
        }
        return bbox.diagonal().length() * GEO::CmdLine::get_arg_double( "epsilon" );
    }
}

namespace RINGMesh {

    GeoModel::GeoModel()
        : mesh( *this ), epsilon_( -1 ), wells_( nullptr )
    {
    }

    index_t GeoModel::nb_mesh_entities( const MeshEntityType& type ) const
    {
        if( MeshEntityTypeManager::is_line( type ) ) {
            return nb_lines();
        } else if( MeshEntityTypeManager::is_corner( type ) ) {
            return nb_corners();
        } else if( MeshEntityTypeManager::is_surface( type ) ) {
            return nb_surfaces();
        } else if( MeshEntityTypeManager::is_region( type ) ) {
            return nb_regions();
        } else {
            ringmesh_assert_not_reached;
            return 0;
        }
    }

    const GeoModelMeshEntity& GeoModel::mesh_entity( gmme_id id ) const
    {
        const MeshEntityType& type = id.type();
        index_t index = id.index();
        if( MeshEntityTypeManager::is_line( type ) ) {
            return line( index );
        } else if( MeshEntityTypeManager::is_corner( type ) ) {
            return corner( index );
        } else if( MeshEntityTypeManager::is_surface( type ) ) {
            return surface( index );
        } else if( MeshEntityTypeManager::is_region( type ) ) {
            return region( index );
        }
        ringmesh_assert_not_reached;
        return surface( 0 );
    }

    const std::vector< std::unique_ptr< GeoModelMeshEntity > >& GeoModel::mesh_entities(
        const MeshEntityType& type ) const
    {
        if( MeshEntityTypeManager::is_corner( type ) ) {
            return corners_;
        } else if( MeshEntityTypeManager::is_line( type ) ) {
            return lines_;
        } else if( MeshEntityTypeManager::is_surface( type ) ) {
            return surfaces_;
        } else if( MeshEntityTypeManager::is_region( type ) ) {
            return regions_;
        } else {
            ringmesh_assert_not_reached;
            return surfaces_;
        }
    }

    const Corner& GeoModel::corner( index_t index ) const
    {
        ringmesh_assert( index < corners_.size() );
        return *static_cast< const Corner* >( corners_[index].get() );
    }
    const Line& GeoModel::line( index_t index ) const
    {
        ringmesh_assert( index < lines_.size() );
        return *static_cast< const Line* >( lines_[index].get() );
    }
    const Surface& GeoModel::surface( index_t index ) const
    {
        ringmesh_assert( index < surfaces_.size() );
        return *static_cast< const Surface* >( surfaces_[index].get() );
    }
    const Region& GeoModel::region( index_t index ) const
    {
        ringmesh_assert( index < regions_.size() );
        return *static_cast< const Region* >( regions_[index].get() );
    }

    void GeoModel::set_wells( const WellGroup* wells )
    {
        wells_ = wells;
    }

    double GeoModel::epsilon() const
    {
        if( epsilon_ == -1 ) {
            epsilon_ = compute_percentage_bbox_diagonal( *this );
        }
        return epsilon_;
    }

    GeoModel::SurfaceSide GeoModel::get_voi_surfaces() const
    {
        SurfaceSide surface_side;
        std::vector< index_t >& voi_surfaces = surface_side.surfaces_;
        voi_surfaces.reserve( nb_surfaces() );
        std::vector< bool >& voi_surface_region_side = surface_side.sides_;
        voi_surface_region_side.reserve( nb_surfaces() );

        for( index_t surface_i = 0; surface_i < nb_surfaces(); ++surface_i ) {
            const Surface& cur_surface = surface( surface_i );
            if( cur_surface.is_on_voi() ) {
                ringmesh_assert( cur_surface.nb_incident_entities() == 1 );
                voi_surfaces.push_back( surface_i );
                const Region& incident_region = cur_surface.incident_entity( 0 );

                index_t local_boundary_id = NO_ID;
                for( index_t region_boundary_i = 0;
                    region_boundary_i < incident_region.nb_boundaries();
                    ++region_boundary_i ) {
                    if( incident_region.boundary_gmme( region_boundary_i ).index()
                        == surface_i ) {
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

} // namespace
