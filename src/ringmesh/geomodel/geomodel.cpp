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

    void compute_mesh_entity_bbox( const GeoModelMeshEntity& entity, Box< 3 >& bbox )
    {
        for( index_t v = 0; v < entity.nb_vertices(); v++ ) {
            bbox.add_point( entity.vertex( v ) );
        }
    }

    double compute_percentage_bbox_diagonal( const GeoModel& gm )
    {
        Box< 3 > bbox;
        if( gm.universe().nb_boundaries() > 0 ) {
            const Universe& universe = gm.universe();
            for( index_t s = 0; s < universe.nb_boundaries(); s++ ) {
                compute_mesh_entity_bbox(
                    gm.mesh_entity( universe.boundary_gmme( s ) ), bbox );
            }
        } else {
            if( gm.nb_surfaces() > 0 ) {
                for( index_t s = 0; s < gm.nb_surfaces(); s++ ) {
                    compute_mesh_entity_bbox( gm.surface( s ), bbox );
                }
            } else if( gm.nb_lines() > 0 ) {
                for( index_t l = 0; l < gm.nb_lines(); l++ ) {
                    compute_mesh_entity_bbox( gm.line( l ), bbox );
                }
            } else {
                ringmesh_assert( gm.nb_corners() > 0 );
                for( index_t c = 0; c < gm.nb_corners(); c++ ) {
                    bbox.add_point( gm.corner( c ).vertex( 0 ) );
                }
            }
        }
        return bbox.diagonal().length() * GEO::CmdLine::get_arg_double( "epsilon" );
    }
}

namespace RINGMesh {

    GeoModel::GeoModel()
        : mesh( *this ), epsilon_( -1 ), universe_( *this ), wells_( nullptr )
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

} // namespace
