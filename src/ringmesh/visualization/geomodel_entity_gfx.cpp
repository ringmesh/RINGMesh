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

#include <ringmesh/visualization/geomodel_entity_gfx.h>

#ifdef RINGMESH_WITH_GRAPHICS

#include <ringmesh/geomodel/geomodel.h>
#include <ringmesh/geomodel/geomodel_entity.h>
#include <ringmesh/geomodel/geomodel_mesh_entity.h>
#include <ringmesh/geomodel/geomodel_geological_entity.h>

#include <ringmesh/visualization/geomodel_gfx.h>
#include <ringmesh/visualization/mesh_entity_gfx.h>
#include <ringmesh/visualization/geogram_gfx.h>

namespace RINGMesh {

    GeoModelGfxEntity::GeoModelGfxEntity( GeoModelGfx& gfx )
        : gfx_( gfx )
    {
    }

    void GeoModelGfxEntity::set_scalar_attribute(
        GEO::MeshElementsFlags subelements,
        const std::string& name,
        double attr_min,
        double attr_max,
        GLuint colormap_texture )
    {
        for( std::unique_ptr< MeshEntityGfx >& e : entities_ ) {
            e->set_scalar_attribute( subelements, name, attr_min, attr_max,
                colormap_texture );
        }
    }

    void GeoModelGfxEntity::unset_scalar_attribute()
    {
        for( std::unique_ptr< MeshEntityGfx >& e : entities_ ) {
            e->unset_scalar_attribute();
        }
    }

    void GeoModelGfxEntity::set_vertex_visibility( bool is_visible )
    {
        for( index_t e = 0; e < entities_.size(); e++ ) {
            set_vertex_visibility( e, is_visible );
        }
    }

    void GeoModelGfxEntity::set_vertex_visibility( index_t entity_id, bool is_visible )
    {
        entities_[entity_id]->set_vertex_visible( is_visible );
    }

    void GeoModelGfxEntity::set_vertex_color( float red, float green, float blue )
    {
        for( index_t e = 0; e < entities_.size(); e++ ) {
            set_vertex_color( e, red, green, blue );
        }
    }

    void GeoModelGfxEntity::set_vertex_color( index_t entity_id, float red, float green, float blue )
    {
        entities_[entity_id]->set_vertex_color( red, green, blue );
    }

    void GeoModelGfxEntity::set_vertex_size( index_t s )
    {
        for( index_t e = 0; e < entities_.size(); e++ ) {
            set_vertex_size( e, s );
        }
    }

    void GeoModelGfxEntity::set_vertex_size( index_t entity_id, index_t s )
    {
        entities_[entity_id]->set_vertex_size( s );
    }

    /*****************************************************************/

    CornerGfxEntity::CornerGfxEntity( GeoModelGfx& gfx )
        : GeoModelGfxEntity( gfx )
    {
    }

    PointSetMeshGfx& CornerGfxEntity::corner( index_t corner_id )
    {
        ringmesh_assert( corner_id < entities_.size() );
        return static_cast< PointSetMeshGfx& >( *entities_[corner_id] );
    }

    void CornerGfxEntity::initialize()
    {
        if( entities_.empty() ) {
            entities_.reserve( gfx_.geomodel()->nb_corners() );
            for( index_t e = 0; e < gfx_.geomodel()->nb_corners(); e++ ) {
                entities_.push_back(
                    PointSetMeshGfx::create_gfx(
                        gfx_.geomodel()->corner( e ).low_level_mesh_storage() ) );
            }
        }
    }

    void CornerGfxEntity::draw()
    {
        for( index_t c = 0; c < entities_.size(); c++ ) {
            PointSetMeshGfx& pointset = corner( c );
            if( pointset.get_vertex_visible() ) pointset.draw_vertices();
        }
    }

    /*****************************************************************/

    LineGfxEntity::LineGfxEntity( GeoModelGfx& gfx )
        : GeoModelGfxEntity( gfx )
    {
    }

    LineMeshGfx& LineGfxEntity::line( index_t line_id )
    {
        ringmesh_assert( line_id < entities_.size() );
        return static_cast< LineMeshGfx& >( *entities_[line_id] );
    }

    void LineGfxEntity::initialize()
    {
        if( entities_.empty() ) {
            entities_.reserve( gfx_.geomodel()->nb_lines() );
            for( index_t e = 0; e < gfx_.geomodel()->nb_lines(); e++ ) {
                entities_.push_back(
                    LineMeshGfx::create_gfx(
                        gfx_.geomodel()->line( e ).low_level_mesh_storage() ) );
            }
        }
    }

    void LineGfxEntity::draw()
    {
        for( index_t l = 0; l < entities_.size(); l++ ) {
            LineMeshGfx& line = this->line( l );
            if( line.get_vertex_visible() ) line.draw_vertices();
            if( line.get_edge_visible() ) line.draw_edges();
        }
    }

    void LineGfxEntity::set_line_color( float red, float green, float blue )
    {
        for( index_t l = 0; l < entities_.size(); l++ ) {
            set_line_color( l, red, green, blue );
        }
    }

    void LineGfxEntity::set_line_color( index_t line_id, float red, float green, float blue )
    {
        line( line_id ).set_edge_color( red, green, blue );
    }

    void LineGfxEntity::set_line_visibility( bool is_visible )
    {
        for( index_t l = 0; l < entities_.size(); l++ ) {
            set_line_visibility( l, is_visible );
        }
    }

    void LineGfxEntity::set_line_visibility( index_t line_id, bool is_visible )
    {
        line( line_id ).set_edge_visible( is_visible );
    }

    void LineGfxEntity::set_line_size( index_t size )
    {
        for( index_t l = 0; l < entities_.size(); l++ ) {
            set_line_size( l, size );
        }
    }
    void LineGfxEntity::set_line_size( index_t line_id, index_t size )
    {
        line( line_id ).set_edge_width( size );
    }

    /*****************************************************************/

    SurfaceGfxEntity::SurfaceGfxEntity( GeoModelGfx& gfx )
        : GeoModelGfxEntity( gfx )
    {
    }

    SurfaceMeshGfx& SurfaceGfxEntity::surface( index_t surface_id )
    {
        ringmesh_assert( surface_id < entities_.size() );
        return static_cast< SurfaceMeshGfx& >( *entities_[surface_id] );
    }

    void SurfaceGfxEntity::initialize()
    {
        if( entities_.empty() ) {
            entities_.reserve( gfx_.geomodel()->nb_surfaces() );
            for( index_t e = 0; e < gfx_.geomodel()->nb_surfaces(); e++ ) {
                entities_.push_back(
                    SurfaceMeshGfx::create_gfx(
                        gfx_.geomodel()->surface( e ).low_level_mesh_storage() ) );
            }
        }
    }

    void SurfaceGfxEntity::draw()
    {
        for( index_t s = 0; s < entities_.size(); s++ ) {
            SurfaceMeshGfx& surface = this->surface( s );
            if( surface.get_vertex_visible() ) surface.draw_vertices();
            if( surface.get_surface_visible() ) surface.draw_surface();
        }
    }

    void SurfaceGfxEntity::set_surface_color( float red, float green, float blue )
    {
        for( index_t s = 0; s < entities_.size(); s++ ) {
            set_surface_color( s, red, green, blue );
        }
    }

    void SurfaceGfxEntity::set_surface_color( index_t surface_id, float red, float green, float blue )
    {
        surface( surface_id ).set_surface_color( red, green, blue );
    }

    void SurfaceGfxEntity::set_backface_surface_color( float red, float green, float blue )
    {
        for( index_t s = 0; s < entities_.size(); s++ ) {
            set_backface_surface_color( s, red, green, blue );
        }
    }

    void SurfaceGfxEntity::set_backface_surface_color(
        index_t surface_id,
        float red,
        float green,
        float blue )
    {
        surface( surface_id ).set_backface_surface_color( red, green, blue );
    }

    void SurfaceGfxEntity::set_surface_visibility( bool is_visible )
    {
        for( index_t s = 0; s < entities_.size(); s++ ) {
            set_surface_visibility( s, is_visible );
        }
    }

    void SurfaceGfxEntity::set_surface_visibility( index_t surface_id, bool is_visible )
    {
        surface( surface_id ).set_surface_visible( is_visible );
    }

    void SurfaceGfxEntity::set_mesh_color( float red, float green, float blue )
    {
        for( index_t s = 0; s < entities_.size(); s++ ) {
            set_mesh_color( s, red, green, blue );
        }
    }

    void SurfaceGfxEntity::set_mesh_color( index_t surface_id, float red, float green, float blue )
    {
        surface( surface_id ).set_mesh_color( red, green, blue );
    }

    void SurfaceGfxEntity::set_mesh_visibility( bool is_visible )
    {
        for( index_t s = 0; s < entities_.size(); s++ ) {
            set_mesh_visibility( s, is_visible );
        }
    }

    void SurfaceGfxEntity::set_mesh_visibility( index_t surface_id, bool is_visible )
    {
        surface( surface_id ).set_mesh_visibility( is_visible );
    }

    void SurfaceGfxEntity::set_mesh_size( index_t size )
    {
        for( index_t s = 0; s < entities_.size(); s++ ) {
            set_mesh_size( s, size );
        }
    }

    void SurfaceGfxEntity::set_mesh_size( index_t surface_id, index_t size )
    {
        surface( surface_id ).set_mesh_width( size );
    }

    /*****************************************************************/

    RegionGfxEntity::RegionGfxEntity( GeoModelGfx& gfx )
        : GeoModelGfxEntity( gfx )
    {
    }

    VolumeMeshGfx& RegionGfxEntity::region( index_t region_id )
    {
        ringmesh_assert( region_id < entities_.size() );
        return static_cast< VolumeMeshGfx& >( *entities_[region_id] );
    }

    void RegionGfxEntity::initialize()
    {
        if( entities_.empty() ) {
            entities_.reserve( gfx_.geomodel()->nb_regions() );
            for( index_t e = 0; e < gfx_.geomodel()->nb_regions(); e++ ) {
                entities_.push_back(
                    VolumeMeshGfx::create_gfx(
                        gfx_.geomodel()->region( e ).low_level_mesh_storage() ) );
            }
        }
    }

    void RegionGfxEntity::draw()
    {
        for( index_t r = 0; r < entities_.size(); r++ ) {
            VolumeMeshGfx& region = this->region( r );
            if( region.get_vertex_visible() ) region.draw_vertices();
            if( region.get_region_visible() ) region.draw_volume();
        }
    }

    void RegionGfxEntity::set_mesh_color( float red, float green, float blue )
    {
        for( index_t r = 0; r < entities_.size(); r++ ) {
            set_mesh_color( r, red, green, blue );
        }
    }

    void RegionGfxEntity::set_mesh_color( index_t region_id, float red, float green, float blue )
    {
        region( region_id ).set_mesh_color( red, green, blue );
    }

    void RegionGfxEntity::set_cell_colors_by_type()
    {
        for( index_t r = 0; r < entities_.size(); r++ ) {
            set_cell_colors_by_type( r );
        }
    }

    void RegionGfxEntity::set_draw_cells( CellType type, bool x )
    {
        for( index_t r = 0; r < entities_.size(); r++ ) {
            set_draw_cells( r, type, x );
        }
    }

    void RegionGfxEntity::set_draw_cells( index_t region_id, CellType type, bool x )
    {
        region( region_id ).set_draw_cells( type, x );
    }

    void RegionGfxEntity::set_cell_colors_by_type( index_t region_id )
    {
        region( region_id ).set_cell_colors_by_type();
    }

    void RegionGfxEntity::set_mesh_visibility( bool is_visible )
    {
        for( index_t r = 0; r < entities_.size(); r++ ) {
            set_mesh_visibility( r, is_visible );
        }
    }

    void RegionGfxEntity::set_mesh_visibility( index_t region_id, bool is_visible )
    {
        region( region_id ).set_mesh_visibility( is_visible );
    }

    void RegionGfxEntity::set_mesh_size( index_t size )
    {
        for( index_t r = 0; r < entities_.size(); r++ ) {
            set_mesh_size( r, size );
        }
    }

    void RegionGfxEntity::set_mesh_size( index_t region_id, index_t size )
    {
        region( region_id ).set_mesh_width( size );
    }

    void RegionGfxEntity::set_region_color( float red, float green, float blue )
    {
        for( index_t r = 0; r < entities_.size(); r++ ) {
            set_region_color( r, red, green, blue );
        }
    }

    void RegionGfxEntity::set_region_color( index_t region_id, float red, float green, float blue )
    {
        region( region_id ).set_cells_color( red, green, blue );
    }

    void RegionGfxEntity::set_region_visibility( bool is_visible )
    {
        for( index_t r = 0; r < entities_.size(); r++ ) {
            set_region_visibility( r, is_visible );
        }
    }

    void RegionGfxEntity::set_region_visibility( index_t region_id, bool is_visible )
    {
        region( region_id ).set_region_visible( is_visible );
    }

    void RegionGfxEntity::set_cell_type_visibility( CellType t, bool is_visible )
    {
        for( index_t r = 0; r < entities_.size(); r++ ) {
            set_cell_type_visibility( r, t, is_visible );
        }
    }

    void RegionGfxEntity::set_cell_type_visibility(
        index_t region_id,
        CellType t,
        bool is_visible )
    {
        region( region_id ).set_draw_cells( t, is_visible );
    }

    void RegionGfxEntity::set_shrink( double shrink )
    {
        for( index_t r = 0; r < entities_.size(); r++ ) {
            set_shrink( r, shrink );
        }
    }

    void RegionGfxEntity::set_shrink( index_t region_id, double shrink )
    {
        region( region_id ).set_shrink( shrink );
    }
}

#endif
