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
 * @author Benjamin Chaunvin and Arnaud Botella
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

    void GeoModelGfxEntity::set_vertex_visibility( bool b )
    {
        for( index_t e = 0; e < entities_.size(); e++ ) {
            set_vertex_visibility( e, b );
        }
    }

    void GeoModelGfxEntity::set_vertex_visibility( index_t e, bool b )
    {
        entities_[e]->set_vertices_visible( b );
    }

    void GeoModelGfxEntity::set_vertex_color( float r, float g, float b )
    {
        for( index_t e = 0; e < entities_.size(); e++ ) {
            set_vertex_color( e, r, g, b );
        }
    }

    void GeoModelGfxEntity::set_vertex_color( index_t e, float r, float g, float b )
    {
        entities_[e]->set_vertices_color( r, g, b );
    }

    void GeoModelGfxEntity::set_vertex_size( index_t s )
    {
        for( index_t e = 0; e < entities_.size(); e++ ) {
            set_vertex_size( e, s );
        }
    }

    void GeoModelGfxEntity::set_vertex_size( index_t e, index_t s )
    {
        entities_[e]->set_vertices_size( s );
    }

    void GeoModelGfxEntity::set_mesh_element_color( float r, float g, float b )
    {
        for( index_t k = 0; k < entities_.size(); k++ ) {
            set_mesh_element_color( k, r, g, b );
        }
    }

    void GeoModelGfxEntity::set_mesh_element_visibility( bool b )
    {
        for( index_t l = 0; l < entities_.size(); l++ ) {
            set_mesh_element_visibility( l, b );
        }
    }

    void GeoModelGfxEntity::set_mesh_element_size( index_t s )
    {
        for( index_t l = 0; l < entities_.size(); l++ ) {
            set_mesh_element_size( l, s );
        }
    }

    void GeoModelGfxEntity::set_mesh_element_color(
        index_t e,
        float r,
        float g,
        float b )
    {
        ringmesh_unused( e );
        ringmesh_unused( r );
        ringmesh_unused( g );
        ringmesh_unused( b );
    }
    void GeoModelGfxEntity::set_mesh_element_visibility( index_t e, bool b )
    {
        ringmesh_unused( e );
        ringmesh_unused( b );
    }
    void GeoModelGfxEntity::set_mesh_element_size( index_t e, index_t s )
    {
        ringmesh_unused( e );
        ringmesh_unused( s );
    }

    /*****************************************************************/

    CornerGfxEnity::CornerGfxEnity( GeoModelGfx& gfx )
        : GeoModelGfxEntity( gfx )
    {
    }

    PointSetGfx& CornerGfxEnity::corner( index_t c )
    {
        ringmesh_assert( c < entities_.size() );
        return static_cast< PointSetGfx& >( *entities_[c] );
    }

    void CornerGfxEnity::initialize()
    {
        if( entities_.empty() ) {
            entities_.reserve( gfx_.geomodel()->nb_corners() );
            for( index_t e = 0; e < gfx_.geomodel()->nb_corners(); e++ ) {
                entities_.emplace_back(
                    new GeogramPointSetGfx( gfx_, gfx_.geomodel()->corner( e ) ) );
            }
        }
    }

    void CornerGfxEnity::set_mesh_element_color(
        index_t e,
        float r,
        float g,
        float b )
    {
        set_vertex_color( e, r, g, b );
    }
    void CornerGfxEnity::set_mesh_element_visibility( index_t e, bool b )
    {
        set_vertex_visibility( e, b );
    }
    void CornerGfxEnity::set_mesh_element_size( index_t e, index_t s )
    {
        set_vertex_size( e, s );
    }

    void CornerGfxEnity::draw()
    {
        for( index_t c = 0; c < entities_.size(); c++ ) {
            if( entities_[c]->get_vertices_visible() ) entities_[c]->draw_vertices();
        }
    }

    /*****************************************************************/

    LineGfxEntity::LineGfxEntity( GeoModelGfx& gfx )
        : GeoModelGfxEntity( gfx )
    {
    }

    LineGfx& LineGfxEntity::line( index_t l )
    {
        ringmesh_assert( l < entities_.size() );
        return static_cast< LineGfx& >( *entities_[l] );
    }

    void LineGfxEntity::initialize()
    {
        if( entities_.empty() ) {
            entities_.reserve( gfx_.geomodel()->nb_lines() );
            for( index_t e = 0; e < gfx_.geomodel()->nb_lines(); e++ ) {
                entities_.emplace_back(
                    new GeogramLineGfx( gfx_, gfx_.geomodel()->line( e ) ) );
            }
        }
    }

    void LineGfxEntity::draw()
    {
        for( index_t l = 0; l < entities_.size(); l++ ) {
            LineGfx& line = this->line( l );
            if( line.get_vertices_visible() ) line.draw_vertices();
            if( line.get_edges_visible() ) line.draw_edges();
        }
    }

    void LineGfxEntity::set_mesh_element_color(
        index_t l,
        float r,
        float g,
        float b )
    {
        line( l ).set_edges_color( r, g, b );
    }

    void LineGfxEntity::set_mesh_element_visibility( index_t l, bool b )
    {
        line( l ).set_edges_visible( b );
    }

    void LineGfxEntity::set_mesh_element_size( index_t l, index_t s )
    {
        line( l ).set_edges_width( s );
    }

    /*****************************************************************/

    SurfaceGfxEntity::SurfaceGfxEntity( GeoModelGfx& gfx )
        : GeoModelGfxEntity( gfx )
    {
    }

    SurfaceGfx& SurfaceGfxEntity::surface( index_t s )
    {
        ringmesh_assert( s < entities_.size() );
        return static_cast< SurfaceGfx& >( *entities_[s] );
    }

    void SurfaceGfxEntity::initialize()
    {
        if( entities_.empty() ) {
            entities_.reserve( gfx_.geomodel()->nb_surfaces() );
            for( index_t e = 0; e < gfx_.geomodel()->nb_surfaces(); e++ ) {
                entities_.emplace_back(
                    new GeogramSurfaceGfx( gfx_, gfx_.geomodel()->surface( e ) ) );
            }
        }
    }

    void SurfaceGfxEntity::draw()
    {
        for( index_t s = 0; s < entities_.size(); s++ ) {
            SurfaceGfx& surface = this->surface( s );
            if( surface.get_vertices_visible() ) surface.draw_vertices();
            if( surface.get_surface_visible() ) surface.draw_surface();
        }
    }

    void SurfaceGfxEntity::set_mesh_element_color(
        index_t s,
        float r,
        float g,
        float b )
    {
        surface( s ).set_surface_color( r, g, b );
    }

    void SurfaceGfxEntity::set_backface_surfaces_color( float r, float g, float b )
    {
        for( index_t s = 0; s < entities_.size(); s++ ) {
            set_backface_surface_color( s, r, g, b );
        }
    }

    void SurfaceGfxEntity::set_backface_surface_color(
        index_t s,
        float r,
        float g,
        float b )
    {
        surface( s ).set_backface_surface_color( r, g, b );
    }

    void SurfaceGfxEntity::set_mesh_element_visibility( index_t s, bool b )
    {
        surface( s ).set_surface_visible( b );
    }

    void SurfaceGfxEntity::set_mesh_color( float r, float g, float b )
    {
        for( index_t s = 0; s < entities_.size(); s++ ) {
            set_mesh_color( s, r, g, b );
        }
    }

    void SurfaceGfxEntity::set_mesh_color( index_t s, float r, float g, float b )
    {
        surface( s ).set_mesh_color( r, g, b );
    }

    void SurfaceGfxEntity::set_mesh_visibility( bool b )
    {
        for( index_t s = 0; s < entities_.size(); s++ ) {
            set_mesh_visibility( s, b );
        }
    }

    void SurfaceGfxEntity::set_mesh_visibility( index_t s, bool b )
    {
        surface( s ).set_mesh_visibility( b );
    }

    void SurfaceGfxEntity::set_mesh_size( index_t size )
    {
        for( index_t s = 0; s < entities_.size(); s++ ) {
            set_mesh_size( s, size );
        }
    }

    void SurfaceGfxEntity::set_mesh_size( index_t s, index_t size )
    {
        surface( s ).set_mesh_width( size );
    }

    /*****************************************************************/

    RegionGfxEntity::RegionGfxEntity( GeoModelGfx& gfx )
        : GeoModelGfxEntity( gfx )
    {
    }

    VolumeGfx& RegionGfxEntity::region( index_t r )
    {
        ringmesh_assert( r < entities_.size() );
        return static_cast< VolumeGfx& >( *entities_[r] );
    }

    void RegionGfxEntity::initialize()
    {
        if( entities_.empty() ) {
            entities_.reserve( gfx_.geomodel()->nb_regions() );
            for( index_t e = 0; e < gfx_.geomodel()->nb_regions(); e++ ) {
                entities_.emplace_back(
                    new GeogramVolumeGfx( gfx_, gfx_.geomodel()->region( e ) ) );
            }
        }
    }

    void RegionGfxEntity::draw()
    {
        for( index_t r = 0; r < entities_.size(); r++ ) {
            VolumeGfx& region = this->region( r );
            if( region.get_vertices_visible() ) region.draw_vertices();
            if( region.get_region_visible() ) region.draw_volume();
        }
    }

    void RegionGfxEntity::set_mesh_color( float r, float g, float b )
    {
        for( index_t m = 0; m < entities_.size(); m++ ) {
            set_mesh_color( m, r, g, b );
        }
    }

    void RegionGfxEntity::set_mesh_color( index_t m, float r, float g, float b )
    {
        region( m ).set_mesh_color( r, g, b );
    }

    void RegionGfxEntity::set_cells_colors_by_type()
    {
        for( index_t m = 0; m < entities_.size(); m++ ) {
            set_cells_colors_by_type( m );
        }
    }

    void RegionGfxEntity::set_draw_cells( GEO::MeshCellType type, bool x )
    {
        for( index_t m = 0; m < entities_.size(); m++ ) {
            set_draw_cells( m, type, x );
        }
    }

    void RegionGfxEntity::set_draw_cells( index_t m, GEO::MeshCellType type, bool x )
    {
        region( m ).set_draw_cells( type, x );
    }

    void RegionGfxEntity::set_cells_colors_by_type( index_t m )
    {
        region( m ).set_cells_colors_by_type();
    }

    void RegionGfxEntity::set_mesh_visibility( bool b )
    {
        for( index_t m = 0; m < entities_.size(); m++ ) {
            set_mesh_visibility( m, b );
        }
    }

    void RegionGfxEntity::set_mesh_visibility( index_t m, bool b )
    {
        region( m ).set_mesh_visibility( b );
    }

    void RegionGfxEntity::set_mesh_size( index_t s )
    {
        for( index_t m = 0; m < entities_.size(); m++ ) {
            set_mesh_size( m, s );
        }
    }

    void RegionGfxEntity::set_mesh_size( index_t m, index_t s )
    {
        region( m ).set_mesh_width( s );
    }

    void RegionGfxEntity::set_mesh_element_color(
        index_t m,
        float r,
        float g,
        float b )
    {
        region( m ).set_cells_color( r, g, b );
    }

    void RegionGfxEntity::set_mesh_element_visibility( index_t m, bool b )
    {
        region( m ).set_region_visible( b );
    }
    void RegionGfxEntity::set_cell_type_visibility( GEO::MeshCellType t, bool b )
    {
        for( index_t m = 0; m < entities_.size(); m++ ) {
            set_cell_type_visibility( m, t, b );
        }
    }
    void RegionGfxEntity::set_cell_type_visibility(
        index_t m,
        GEO::MeshCellType t,
        bool b )
    {
        region( m ).set_draw_cells( t, b );
    }

    void RegionGfxEntity::set_shrink( double s )
    {
        for( index_t m = 0; m < entities_.size(); m++ ) {
            set_shrink( m, s );
        }
    }

    void RegionGfxEntity::set_shrink( index_t m, double s )
    {
        region( m ).set_shrink( s );
    }
}

#endif
