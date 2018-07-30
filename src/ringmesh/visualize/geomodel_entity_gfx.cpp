/*
 * Copyright (c) 2012-2018, Association Scientifique pour la Geologie et ses
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
 * @file Implementation of visualization of GeoModelEntities
 * @author Benjamin Chauvin and Arnaud Botella
 */

#include <ringmesh/visualize/geomodel_entity_gfx.h>

#ifdef RINGMESH_WITH_GRAPHICS

#include <ringmesh/geomodel/core/geomodel.h>
#include <ringmesh/geomodel/core/geomodel_entity.h>
#include <ringmesh/geomodel/core/geomodel_geological_entity.h>
#include <ringmesh/geomodel/core/geomodel_mesh_entity.h>

#include <ringmesh/visualize/geogram_gfx.h>
#include <ringmesh/visualize/geomodel_gfx.h>
#include <ringmesh/visualize/mesh_entity_gfx.h>

namespace RINGMesh
{
    template < index_t DIMENSION >
    GeoModelGfxEntity< DIMENSION >::GeoModelGfxEntity(
        GeoModelGfx< DIMENSION >& gfx )
        : gfx_( gfx )
    {
    }

    template < index_t DIMENSION >
    void GeoModelGfxEntity< DIMENSION >::set_scalar_attribute(
        GEO::MeshElementsFlags subelements,
        const std::string& name,
        double attr_min,
        double attr_max,
        GLuint colormap_texture )
    {
        for( std::unique_ptr< MeshEntityGfx< DIMENSION > >& e : entities_ )
        {
            e->set_scalar_attribute(
                subelements, name, attr_min, attr_max, colormap_texture );
        }
    }

    template < index_t DIMENSION >
    void GeoModelGfxEntity< DIMENSION >::unset_scalar_attribute()
    {
        for( std::unique_ptr< MeshEntityGfx< DIMENSION > >& e : entities_ )
        {
            e->unset_scalar_attribute();
        }
    }

    template < index_t DIMENSION >
    void GeoModelGfxEntity< DIMENSION >::set_vertex_visibility(
        bool is_visible )
    {
        for( auto e : range( entities_.size() ) )
        {
            set_vertex_visibility( e, is_visible );
        }
    }

    template < index_t DIMENSION >
    void GeoModelGfxEntity< DIMENSION >::set_vertex_visibility(
        index_t entity_id, bool is_visible )
    {
        entities_[entity_id]->set_vertex_visible( is_visible );
    }

    template < index_t DIMENSION >
    void GeoModelGfxEntity< DIMENSION >::set_vertex_color(
        float red, float green, float blue )
    {
        for( auto e : range( entities_.size() ) )
        {
            set_vertex_color( e, red, green, blue );
        }
    }

    template < index_t DIMENSION >
    void GeoModelGfxEntity< DIMENSION >::set_vertex_color(
        index_t entity_id, float red, float green, float blue )
    {
        entities_[entity_id]->set_vertex_color( red, green, blue );
    }

    template < index_t DIMENSION >
    void GeoModelGfxEntity< DIMENSION >::set_vertex_size( index_t s )
    {
        for( auto e : range( entities_.size() ) )
        {
            set_vertex_size( e, s );
        }
    }

    template < index_t DIMENSION >
    void GeoModelGfxEntity< DIMENSION >::set_vertex_size(
        index_t entity_id, index_t s )
    {
        entities_[entity_id]->set_vertex_size( s );
    }

    /*****************************************************************/

    template < index_t DIMENSION >
    CornerGfxEntity< DIMENSION >::CornerGfxEntity(
        GeoModelGfx< DIMENSION >& gfx )
        : GeoModelGfxEntity< DIMENSION >( gfx )
    {
    }

    template < index_t DIMENSION >
    PointSetMeshGfx< DIMENSION >& CornerGfxEntity< DIMENSION >::corner(
        index_t corner_id )
    {
        ringmesh_assert( corner_id < this->entities_.size() );
        return static_cast< PointSetMeshGfx< DIMENSION >& >(
            *this->entities_[corner_id] );
    }

    template < index_t DIMENSION >
    void CornerGfxEntity< DIMENSION >::initialize()
    {
        if( this->entities_.empty() )
        {
            this->entities_.reserve( this->gfx_.geomodel()->nb_corners() );
            for( const auto& corner : this->gfx_.geomodel()->corners() )
            {
                this->entities_.push_back(
                    PointSetMeshGfx< DIMENSION >::create_gfx( corner.mesh() ) );
            }
        }
    }

    template < index_t DIMENSION >
    void CornerGfxEntity< DIMENSION >::draw()
    {
        for( auto c : range( this->entities_.size() ) )
        {
            PointSetMeshGfx< DIMENSION >& pointset = corner( c );
            if( pointset.get_vertex_visible() )
                pointset.draw_vertices();
        }
    }

    /*****************************************************************/

    template < index_t DIMENSION >
    LineGfxEntity< DIMENSION >::LineGfxEntity( GeoModelGfx< DIMENSION >& gfx )
        : GeoModelGfxEntity< DIMENSION >( gfx )
    {
    }

    template < index_t DIMENSION >
    LineMeshGfx< DIMENSION >& LineGfxEntity< DIMENSION >::line(
        index_t line_id )
    {
        ringmesh_assert( line_id < this->entities_.size() );
        return static_cast< LineMeshGfx< DIMENSION >& >(
            *this->entities_[line_id] );
    }

    template < index_t DIMENSION >
    void LineGfxEntity< DIMENSION >::initialize()
    {
        if( this->entities_.empty() )
        {
            this->entities_.reserve( this->gfx_.geomodel()->nb_lines() );
            for( const auto& line : this->gfx_.geomodel()->lines() )
            {
                this->entities_.push_back(
                    LineMeshGfx< DIMENSION >::create_gfx( line.mesh() ) );
            }
        }
    }

    template < index_t DIMENSION >
    void LineGfxEntity< DIMENSION >::draw()
    {
        for( auto l : range( this->entities_.size() ) )
        {
            LineMeshGfx< DIMENSION >& line = this->line( l );
            if( line.get_vertex_visible() )
                line.draw_vertices();
            if( line.get_edge_visible() )
                line.draw_edges();
        }
    }

    template < index_t DIMENSION >
    void LineGfxEntity< DIMENSION >::set_line_color(
        float red, float green, float blue )
    {
        for( auto l : range( this->entities_.size() ) )
        {
            set_line_color( l, red, green, blue );
        }
    }

    template < index_t DIMENSION >
    void LineGfxEntity< DIMENSION >::set_line_color(
        index_t line_id, float red, float green, float blue )
    {
        line( line_id ).set_edge_color( red, green, blue );
    }

    template < index_t DIMENSION >
    void LineGfxEntity< DIMENSION >::set_line_visibility( bool is_visible )
    {
        for( auto l : range( this->entities_.size() ) )
        {
            set_line_visibility( l, is_visible );
        }
    }

    template < index_t DIMENSION >
    void LineGfxEntity< DIMENSION >::set_line_visibility(
        index_t line_id, bool is_visible )
    {
        line( line_id ).set_edge_visible( is_visible );
    }

    template < index_t DIMENSION >
    void LineGfxEntity< DIMENSION >::set_line_size( index_t size )
    {
        for( auto l : range( this->entities_.size() ) )
        {
            set_line_size( l, size );
        }
    }

    template < index_t DIMENSION >
    void LineGfxEntity< DIMENSION >::set_line_size(
        index_t line_id, index_t size )
    {
        line( line_id ).set_edge_width( size );
    }

    /*****************************************************************/

    template < index_t DIMENSION >
    SurfaceGfxEntity< DIMENSION >::SurfaceGfxEntity(
        GeoModelGfx< DIMENSION >& gfx )
        : GeoModelGfxEntity< DIMENSION >( gfx )
    {
    }

    template < index_t DIMENSION >
    SurfaceMeshGfx< DIMENSION >& SurfaceGfxEntity< DIMENSION >::surface(
        index_t surface_id )
    {
        ringmesh_assert( surface_id < this->entities_.size() );
        return static_cast< SurfaceMeshGfx< DIMENSION >& >(
            *this->entities_[surface_id] );
    }

    template < index_t DIMENSION >
    void SurfaceGfxEntity< DIMENSION >::initialize()
    {
        if( this->entities_.empty() )
        {
            this->entities_.reserve( this->gfx_.geomodel()->nb_surfaces() );
            for( const auto& surface : this->gfx_.geomodel()->surfaces() )
            {
                this->entities_.push_back(
                    SurfaceMeshGfx< DIMENSION >::create_gfx( surface.mesh() ) );
            }
        }
    }

    template < index_t DIMENSION >
    void SurfaceGfxEntity< DIMENSION >::draw()
    {
        for( auto s : range( this->entities_.size() ) )
        {
            SurfaceMeshGfx< DIMENSION >& surface = this->surface( s );
            if( surface.get_vertex_visible() )
                surface.draw_vertices();
            if( surface.get_surface_visible() )
                surface.draw_surface();
        }
    }

    template < index_t DIMENSION >
    void SurfaceGfxEntity< DIMENSION >::set_surface_color(
        float red, float green, float blue )
    {
        for( auto s : range( this->entities_.size() ) )
        {
            set_surface_color( s, red, green, blue );
        }
    }

    template < index_t DIMENSION >
    void SurfaceGfxEntity< DIMENSION >::set_surface_color(
        index_t surface_id, float red, float green, float blue )
    {
        surface( surface_id ).set_surface_color( red, green, blue );
    }

    template < index_t DIMENSION >
    void SurfaceGfxEntity< DIMENSION >::set_backface_surface_color(
        float red, float green, float blue )
    {
        for( auto s : range( this->entities_.size() ) )
        {
            set_backface_surface_color( s, red, green, blue );
        }
    }

    template < index_t DIMENSION >
    void SurfaceGfxEntity< DIMENSION >::set_backface_surface_color(
        index_t surface_id, float red, float green, float blue )
    {
        surface( surface_id ).set_backface_surface_color( red, green, blue );
    }

    template < index_t DIMENSION >
    void SurfaceGfxEntity< DIMENSION >::set_surface_visibility(
        bool is_visible )
    {
        for( auto s : range( this->entities_.size() ) )
        {
            set_surface_visibility( s, is_visible );
        }
    }

    template < index_t DIMENSION >
    void SurfaceGfxEntity< DIMENSION >::set_surface_visibility(
        index_t surface_id, bool is_visible )
    {
        surface( surface_id ).set_surface_visible( is_visible );
    }

    template < index_t DIMENSION >
    void SurfaceGfxEntity< DIMENSION >::set_mesh_color(
        float red, float green, float blue )
    {
        for( auto s : range( this->entities_.size() ) )
        {
            set_mesh_color( s, red, green, blue );
        }
    }

    template < index_t DIMENSION >
    void SurfaceGfxEntity< DIMENSION >::set_mesh_color(
        index_t surface_id, float red, float green, float blue )
    {
        surface( surface_id ).set_mesh_color( red, green, blue );
    }

    template < index_t DIMENSION >
    void SurfaceGfxEntity< DIMENSION >::set_mesh_visibility( bool is_visible )
    {
        for( auto s : range( this->entities_.size() ) )
        {
            set_mesh_visibility( s, is_visible );
        }
    }

    template < index_t DIMENSION >
    void SurfaceGfxEntity< DIMENSION >::set_mesh_visibility(
        index_t surface_id, bool is_visible )
    {
        surface( surface_id ).set_mesh_visibility( is_visible );
    }

    template < index_t DIMENSION >
    void SurfaceGfxEntity< DIMENSION >::set_mesh_size( index_t size )
    {
        for( auto s : range( this->entities_.size() ) )
        {
            set_mesh_size( s, size );
        }
    }

    template < index_t DIMENSION >
    void SurfaceGfxEntity< DIMENSION >::set_mesh_size(
        index_t surface_id, index_t size )
    {
        surface( surface_id ).set_mesh_width( size );
    }

    /*****************************************************************/

    template < index_t DIMENSION >
    RegionGfxEntity< DIMENSION >::RegionGfxEntity( GeoModelGfx3D& gfx )
        : GeoModelGfxEntity< DIMENSION >( gfx )
    {
    }

    template < index_t DIMENSION >
    VolumeMeshGfx< DIMENSION >& RegionGfxEntity< DIMENSION >::region(
        index_t region_id )
    {
        ringmesh_assert( region_id < entities_.size() );
        return static_cast< VolumeMeshGfx< DIMENSION >& >(
            *entities_[region_id] );
    }

    template < index_t DIMENSION >
    void RegionGfxEntity< DIMENSION >::initialize()
    {
        if( entities_.empty() )
        {
            entities_.reserve( gfx_.geomodel()->nb_regions() );
            for( const auto& region : gfx_.geomodel()->regions() )
            {
                entities_.push_back(
                    VolumeMeshGfx< DIMENSION >::create_gfx( region.mesh() ) );
            }
        }
    }

    template < index_t DIMENSION >
    void RegionGfxEntity< DIMENSION >::draw()
    {
        for( auto r : range( entities_.size() ) )
        {
            VolumeMeshGfx< DIMENSION >& region = this->region( r );
            if( region.get_vertex_visible() )
                region.draw_vertices();
            if( region.get_region_visible() )
                region.draw_volume();
        }
    }

    template < index_t DIMENSION >
    void RegionGfxEntity< DIMENSION >::set_mesh_color(
        float red, float green, float blue )
    {
        for( auto r : range( entities_.size() ) )
        {
            set_mesh_color( r, red, green, blue );
        }
    }

    template < index_t DIMENSION >
    void RegionGfxEntity< DIMENSION >::set_mesh_color(
        index_t region_id, float red, float green, float blue )
    {
        region( region_id ).set_mesh_color( red, green, blue );
    }

    template < index_t DIMENSION >
    void RegionGfxEntity< DIMENSION >::set_cell_colors_by_type()
    {
        for( auto r : range( entities_.size() ) )
        {
            set_cell_colors_by_type( r );
        }
    }

    template < index_t DIMENSION >
    void RegionGfxEntity< DIMENSION >::set_draw_cells( CellType type, bool x )
    {
        for( auto r : range( entities_.size() ) )
        {
            set_draw_cells( r, type, x );
        }
    }

    template < index_t DIMENSION >
    void RegionGfxEntity< DIMENSION >::set_draw_cells(
        index_t region_id, CellType type, bool x )
    {
        region( region_id ).set_draw_cells( type, x );
    }

    template < index_t DIMENSION >
    void RegionGfxEntity< DIMENSION >::set_cell_colors_by_type(
        index_t region_id )
    {
        region( region_id ).set_cell_colors_by_type();
    }

    template < index_t DIMENSION >
    void RegionGfxEntity< DIMENSION >::set_mesh_visibility( bool is_visible )
    {
        for( auto r : range( entities_.size() ) )
        {
            set_mesh_visibility( r, is_visible );
        }
    }

    template < index_t DIMENSION >
    void RegionGfxEntity< DIMENSION >::set_mesh_visibility(
        index_t region_id, bool is_visible )
    {
        region( region_id ).set_mesh_visibility( is_visible );
    }

    template < index_t DIMENSION >
    void RegionGfxEntity< DIMENSION >::set_mesh_size( index_t size )
    {
        for( auto r : range( entities_.size() ) )
        {
            set_mesh_size( r, size );
        }
    }

    template < index_t DIMENSION >
    void RegionGfxEntity< DIMENSION >::set_mesh_size(
        index_t region_id, index_t size )
    {
        region( region_id ).set_mesh_width( size );
    }

    template < index_t DIMENSION >
    void RegionGfxEntity< DIMENSION >::set_region_color(
        float red, float green, float blue )
    {
        for( auto r : range( entities_.size() ) )
        {
            set_region_color( r, red, green, blue );
        }
    }

    template < index_t DIMENSION >
    void RegionGfxEntity< DIMENSION >::set_region_color(
        index_t region_id, float red, float green, float blue )
    {
        region( region_id ).set_cells_color( red, green, blue );
    }

    template < index_t DIMENSION >
    void RegionGfxEntity< DIMENSION >::set_region_visibility( bool is_visible )
    {
        for( auto r : range( entities_.size() ) )
        {
            set_region_visibility( r, is_visible );
        }
    }

    template < index_t DIMENSION >
    void RegionGfxEntity< DIMENSION >::set_region_visibility(
        index_t region_id, bool is_visible )
    {
        region( region_id ).set_region_visible( is_visible );
    }

    template < index_t DIMENSION >
    void RegionGfxEntity< DIMENSION >::set_cell_type_visibility(
        CellType t, bool is_visible )
    {
        for( auto r : range( entities_.size() ) )
        {
            set_cell_type_visibility( r, t, is_visible );
        }
    }

    template < index_t DIMENSION >
    void RegionGfxEntity< DIMENSION >::set_cell_type_visibility(
        index_t region_id, CellType t, bool is_visible )
    {
        region( region_id ).set_draw_cells( t, is_visible );
    }

    template < index_t DIMENSION >
    void RegionGfxEntity< DIMENSION >::set_shrink( double shrink )
    {
        for( auto r : range( entities_.size() ) )
        {
            set_shrink( r, shrink );
        }
    }

    template < index_t DIMENSION >
    void RegionGfxEntity< DIMENSION >::set_shrink(
        index_t region_id, double shrink )
    {
        region( region_id ).set_shrink( shrink );
    }

    template class visualize_api GeoModelGfxEntity< 2 >;
    template class visualize_api CornerGfxEntity< 2 >;
    template class visualize_api LineGfxEntity< 2 >;
    template class visualize_api SurfaceGfxEntity< 2 >;

    template class visualize_api GeoModelGfxEntity< 3 >;
    template class visualize_api CornerGfxEntity< 3 >;
    template class visualize_api LineGfxEntity< 3 >;
    template class visualize_api SurfaceGfxEntity< 3 >;
    template class visualize_api RegionGfxEntity< 3 >;
} // namespace RINGMesh

#endif
