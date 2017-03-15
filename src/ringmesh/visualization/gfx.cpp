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

#include <ringmesh/visualization/gfx.h>

#ifdef RINGMESH_WITH_GRAPHICS

#include <ringmesh/geomodel/geomodel.h>
#include <ringmesh/geomodel/geomodel_entity.h>
#include <ringmesh/geomodel/geomodel_mesh_entity.h>
#include <ringmesh/geomodel/geomodel_geological_entity.h>

namespace {
    using namespace RINGMesh ;

    std::string get_attribute_name_with_coordinate(
        const std::string& name,
        index_t coordinate )
    {
        return name + "[" + GEO::String::to_string( coordinate ) + "]" ;
    }

    void compute_attribute_range(
        GEO::ReadOnlyScalarAttributeAdapter& attribute,
        double& min,
        double& max )
    {
        if( attribute.is_bound() ) {
            for( index_t i = 0; i < attribute.size(); ++i ) {
                double value = attribute[i] ;
                min = GEO::geo_min( min, value ) ;
                max = GEO::geo_max( max, value ) ;
            }
        }
    }
}
namespace RINGMesh {

    class MeshEntityGfx: public GEO::MeshGfx {
    ringmesh_disable_copy( MeshEntityGfx ) ;
    public:
        MeshEntityGfx(
            const GeoModelGfx& gfx,
            const GEO::Mesh& mesh,
            bool vertice_visible )
            : vertices_visible_( vertice_visible ), gfx_( gfx )
        {
            set_mesh( &mesh ) ;
        }
        virtual ~MeshEntityGfx()
        {
        }
        void need_to_update()
        {
            buffer_objects_dirty_ = true ;
            attributes_buffer_objects_dirty_ = true ;
        }

        void draw_vertices()
        {
            GEO::MeshGfx::draw_vertices() ;
        }
        void draw_edges()
        {
            GEO::MeshGfx::draw_edges() ;
        }

        void set_vertices_visible( bool b )
        {
            vertices_visible_ = b ;
        }
        bool get_vertices_visible() const
        {
            return vertices_visible_ ;
        }

    protected:
        bool vertices_visible_ ;

        const GeoModelGfx& gfx_ ;

    } ;

    class CornerGfx: public MeshEntityGfx {
    public:
        CornerGfx( const GeoModelGfx& gfx, const Corner& corner )
            : MeshEntityGfx( gfx, corner.gfx_mesh(), true )
        {
            set_points_color( 1, 0, 0 ) ;
        }
    } ;

    class LineGfx: public MeshEntityGfx {
    public:
        LineGfx( const GeoModelGfx& gfx, const Line& line )
            : MeshEntityGfx( gfx, line.gfx_mesh(), false ), edges_visible_( true )
        {
            set_points_color( 1, 1, 1 ) ;
            set_mesh_color( 1, 1, 1 ) ;
        }
        void set_edges_visible( bool b )
        {
            edges_visible_ = b ;
        }
        bool get_edges_visible() const
        {
            return edges_visible_ ;
        }

    private:
        bool edges_visible_ ;

    } ;

    class SurfaceGfx: public MeshEntityGfx {
    public:
        SurfaceGfx( const GeoModelGfx& gfx, const Surface& surface )
            :
                MeshEntityGfx( gfx, surface.gfx_mesh(), false ),
                surface_visible_( true )
        {
        }

        void draw_surface()
        {
            GEO::MeshGfx::draw_surface() ;
        }

        void set_surface_visible( bool b )
        {
            surface_visible_ = b ;
        }
        bool get_surface_visible() const
        {
            return surface_visible_ ;
        }
    private:
        bool surface_visible_ ;

    } ;

    class RegionGfx: public MeshEntityGfx {
    public:
        RegionGfx( const GeoModelGfx& gfx, const Region& region )
            :
                MeshEntityGfx( gfx, region.gfx_mesh(), false ),
                region_visible_( true ),
                surface_visible_( false ),
                edges_visible_( false )
        {
            set_points_color( 0.0, 0.0, 0.0 ) ;
        }
        void set_edges_visible( bool b )
        {
            edges_visible_ = b ;
        }
        bool get_edges_visible() const
        {
            return edges_visible_ ;
        }
        void set_surface_visible( bool b )
        {
            surface_visible_ = b ;
        }
        bool get_surface_visible() const
        {
            return surface_visible_ ;
        }
        void set_region_visible( bool b )
        {
            region_visible_ = b ;
        }
        bool get_region_visible() const
        {
            return region_visible_ ;
        }

    private:
        bool region_visible_ ;
        bool surface_visible_ ;
        bool edges_visible_ ;

    } ;

    /*****************************************************************/

    GeoModelGfxManager::GeoModelGfxManager( GeoModelGfx& gfx )
        : gfx_( gfx )
    {
    }

    GeoModelGfxManager::~GeoModelGfxManager()
    {
        for( index_t e = 0; e < entities_.size(); e++ ) {
            delete entities_[e] ;
        }
    }

    void GeoModelGfxManager::need_to_update()
    {
        for( index_t e = 0; e < entities_.size(); e++ ) {
            entities_[e]->need_to_update() ;
        }
    }

    void GeoModelGfxManager::set_scalar_attribute(
        GEO::MeshElementsFlags subelements,
        const std::string& name,
        double attr_min,
        double attr_max,
        GLuint colormap_texture )
    {
        for( index_t e = 0; e < entities_.size(); e++ ) {
            entities_[e]->set_scalar_attribute( subelements, name, attr_min,
                attr_max, colormap_texture ) ;
        }
    }

    void GeoModelGfxManager::unset_scalar_attribute()
    {
        for( index_t e = 0; e < entities_.size(); e++ ) {
            entities_[e]->unset_scalar_attribute() ;
        }
    }

    /*!
     * Sets the vertex entity visibility to all the entities
     * @param[in] b the visibility
     */
    void GeoModelGfxManager::set_vertex_visibility( bool b )
    {
        for( index_t e = 0; e < entities_.size(); e++ ) {
            set_vertex_visibility( e, b ) ;
        }
    }
    /*!
     * Sets the vertex entity visibility
     * @param[in] e the entity index
     * @param[in] b the visibility
     */
    void GeoModelGfxManager::set_vertex_visibility( index_t e, bool b )
    {
        entities_[e]->set_vertices_visible( b ) ;
    }

    /*!
     * Sets the entity color to all the entities
     * @param[in] r the red component of the color in [0.0, 1.0]
     * @param[in] g the green component of the color in [0.0, 1.0]
     * @param[in] b the blue component of the color in [0.0, 1.0]
     */
    void GeoModelGfxManager::set_vertex_color( float r, float g, float b )
    {
        for( index_t e = 0; e < entities_.size(); e++ ) {
            set_vertex_color( e, r, g, b ) ;
        }
    }
    /*!
     * Sets the vertex color
     * @param[in] e the entity index
     * @param[in] r the red component of the color in [0.0, 1.0]
     * @param[in] g the green component of the color in [0.0, 1.0]
     * @param[in] b the blue component of the color in [0.0, 1.0]
     */
    void GeoModelGfxManager::set_vertex_color( index_t e, float r, float g, float b )
    {
        entities_[e]->set_points_color( r, g, b ) ;
    }
    /*!
     * Sets the vertex size to all the elements
     * @param[in] s the size
     */
    void GeoModelGfxManager::set_vertex_size( index_t s )
    {
        for( index_t e = 0; e < entities_.size(); e++ ) {
            set_vertex_size( e, s ) ;
        }
    }
    /*!
     * Sets the vertex size to all the elements
     * @param[in] e the entity index
     * @param[in] s the size
     */
    void GeoModelGfxManager::set_vertex_size( index_t e, index_t s )
    {
        entities_[e]->set_points_size( static_cast< float >( s ) ) ;
    }

    /*!
     * Sets the entity color to all the elements
     * @param[in] r the red component of the color in [0.0, 1.0]
     * @param[in] g the green component of the color in [0.0, 1.0]
     * @param[in] b the blue component of the color in [0.0, 1.0]
     */
    void GeoModelGfxManager::set_mesh_element_color( float r, float g, float b )
    {
        for( index_t k = 0; k < entities_.size(); k++ ) {
            set_mesh_element_color( k, r, g, b ) ;
        }
    }
    /*!
     * Sets the entity visibility to all the elements
     * @param[in] b the visibility
     */
    void GeoModelGfxManager::set_mesh_element_visibility( bool b )
    {
        for( index_t l = 0; l < entities_.size(); l++ ) {
            set_mesh_element_visibility( l, b ) ;
        }
    }
    /*!
     * Sets the mesh_element entity size to all the elements
     * @param[in] s the size
     */
    void GeoModelGfxManager::set_mesh_element_size( index_t s )
    {
        for( index_t l = 0; l < entities_.size(); l++ ) {
            set_mesh_element_size( l, s ) ;
        }
    }

    void GeoModelGfxManager::set_mesh_element_color(
        index_t e,
        float r,
        float g,
        float b )
    {
        ringmesh_unused( e ) ;
        ringmesh_unused( r ) ;
        ringmesh_unused( g ) ;
        ringmesh_unused( b ) ;
    }
    void GeoModelGfxManager::set_mesh_element_visibility( index_t e, bool b )
    {
        ringmesh_unused( e ) ;
        ringmesh_unused( b ) ;
    }
    void GeoModelGfxManager::set_mesh_element_size( index_t e, index_t s )
    {
        ringmesh_unused( e ) ;
        ringmesh_unused( s ) ;
    }

    /*****************************************************************/

    CornerGfxManager::CornerGfxManager( GeoModelGfx& gfx )
        : GeoModelGfxManager( gfx )
    {
    }

    void CornerGfxManager::initialize()
    {
        if( entities_.empty() ) {
            entities_.resize( gfx_.geomodel()->nb_corners(), nil ) ;
            for( index_t e = 0; e < entities_.size(); e++ ) {
                entities_[e] = new CornerGfx( gfx_, gfx_.geomodel()->corner( e ) ) ;
            }
        }
    }

    void CornerGfxManager::set_mesh_element_color(
        index_t e,
        float r,
        float g,
        float b )
    {
        set_vertex_color( e, r, g, b ) ;
    }
    void CornerGfxManager::set_mesh_element_visibility( index_t e, bool b )
    {
        set_vertex_visibility( e, b ) ;
    }
    void CornerGfxManager::set_mesh_element_size( index_t e, index_t s )
    {
        set_vertex_size( e, s ) ;
    }

    /*!
     * Draws the corners
     */
    void CornerGfxManager::draw()
    {
        for( index_t c = 0; c < entities_.size(); c++ ) {
            if( entities_[c]->get_vertices_visible() )
                entities_[c]->draw_vertices() ;
        }
    }

    /*****************************************************************/

    LineGfxManager::LineGfxManager( GeoModelGfx& gfx )
        : GeoModelGfxManager( gfx )
    {
    }

    void LineGfxManager::initialize()
    {
        if( entities_.empty() ) {
            entities_.resize( gfx_.geomodel()->nb_lines(), nil ) ;
            for( index_t e = 0; e < entities_.size(); e++ ) {
                entities_[e] = new LineGfx( gfx_, gfx_.geomodel()->line( e ) ) ;
            }
        }
    }

    /*!
     * Draws the lines
     */
    void LineGfxManager::draw()
    {
        for( index_t l = 0; l < entities_.size(); l++ ) {
            LineGfx* line = dynamic_cast< LineGfx* >( entities_[l] ) ;
            if( line->get_vertices_visible() ) line->draw_vertices() ;
            if( line->get_edges_visible() ) line->draw_edges() ;
        }
    }
    /*!
     * Sets the line color
     * @param[in] l the line index
     * @param[in] r the red component of the color in [0.0, 1.0]
     * @param[in] g the green component of the color in [0.0, 1.0]
     * @param[in] b the blue component of the color in [0.0, 1.0]
     */
    void LineGfxManager::set_mesh_element_color(
        index_t l,
        float r,
        float g,
        float b )
    {
        entities_[l]->set_mesh_color( r, g, b ) ;
    }
    /*!
     * Sets the line visibility
     * @param[in] l the line index
     * @param[in] b the visibility
     */
    void LineGfxManager::set_mesh_element_visibility( index_t l, bool b )
    {
        dynamic_cast< LineGfx* >( entities_[l] )->set_edges_visible( b ) ;
    }
    /*!
     * Sets the mesh_element line size
     * @param[in] l the line index
     * @param[in] s the size
     */
    void LineGfxManager::set_mesh_element_size( index_t l, index_t s )
    {
        entities_[l]->set_mesh_width( s ) ;
    }

    /*****************************************************************/

    SurfaceGfxManager::SurfaceGfxManager( GeoModelGfx& gfx )
        : GeoModelGfxManager( gfx )
    {
    }

    void SurfaceGfxManager::initialize()
    {
        if( entities_.empty() ) {
            entities_.resize( gfx_.geomodel()->nb_surfaces(), nil ) ;
            for( index_t e = 0; e < entities_.size(); e++ ) {
                entities_[e] = new SurfaceGfx( gfx_,
                    gfx_.geomodel()->surface( e ) ) ;
            }
        }
    }

    /*!
     * Draws the surfaces
     */
    void SurfaceGfxManager::draw()
    {
        for( index_t s = 0; s < entities_.size(); s++ ) {
            SurfaceGfx* surface = dynamic_cast< SurfaceGfx* >( entities_[s] ) ;
            if( surface->get_vertices_visible() ) surface->draw_vertices() ;
            if( surface->get_surface_visible() ) surface->draw_surface() ;
        }
    }
    /*!
     * Sets the surface color
     * @param[in] e the surface index
     * @param[in] r the red component of the color in [0.0, 1.0]
     * @param[in] g the green component of the color in [0.0, 1.0]
     * @param[in] b the blue component of the color in [0.0, 1.0]
     */
    void SurfaceGfxManager::set_mesh_element_color(
        index_t e,
        float r,
        float g,
        float b )
    {
        entities_[e]->set_surface_color( r, g, b ) ;
    }
    /*!
     * Sets the backface surface color to all the elements
     * @param[in] r the red component of the color in [0.0, 1.0]
     * @param[in] g the green component of the color in [0.0, 1.0]
     * @param[in] b the blue component of the color in [0.0, 1.0]
     */
    void SurfaceGfxManager::set_backface_surfaces_color( float r, float g, float b )
    {
        for( index_t s = 0; s < entities_.size(); s++ ) {
            set_backface_surface_color( s, r, g, b ) ;
        }
    }
    /*!
     * Sets the backsurface surface color
     * @param[in] s the surface index
     * @param[in] r the red component of the color in [0.0, 1.0]
     * @param[in] g the green component of the color in [0.0, 1.0]
     * @param[in] b the blue component of the color in [0.0, 1.0]
     */
    void SurfaceGfxManager::set_backface_surface_color(
        index_t s,
        float r,
        float g,
        float b )
    {
        entities_[s]->set_backface_surface_color( r, g, b ) ;
    }
    /*!
     * Sets the surface visibility
     * @param[in] s the surface index
     * @param[in] b the visibility
     */
    void SurfaceGfxManager::set_mesh_element_visibility( index_t s, bool b )
    {
        dynamic_cast< SurfaceGfx* >( entities_[s] )->set_surface_visible( b ) ;
    }
    /*!
     * Sets the mesh surface color to all the elements
     * @param[in] r the red component of the color in [0.0, 1.0]
     * @param[in] g the green component of the color in [0.0, 1.0]
     * @param[in] b the blue component of the color in [0.0, 1.0]
     */
    void SurfaceGfxManager::set_mesh_color( float r, float g, float b )
    {
        for( index_t s = 0; s < entities_.size(); s++ ) {
            set_mesh_color( s, r, g, b ) ;
        }
    }
    /*!
     * Sets the mesh surface color
     * @param[in] s the surface index
     * @param[in] r the red component of the color in [0.0, 1.0]
     * @param[in] g the green component of the color in [0.0, 1.0]
     * @param[in] b the blue component of the color in [0.0, 1.0]
     */
    void SurfaceGfxManager::set_mesh_color( index_t s, float r, float g, float b )
    {
        entities_[s]->set_mesh_color( r, g, b ) ;
    }
    /*!
     * Sets the mesh surface visibility to all the elements
     * @param[in] b the visibility
     */
    void SurfaceGfxManager::set_mesh_visibility( bool b )
    {
        for( index_t s = 0; s < entities_.size(); s++ ) {
            set_mesh_visibility( s, b ) ;
        }
    }
    /*!
     * Sets the mesh surface visibility
     * @param[in] s the surface index
     * @param[in] b the visibility
     */
    void SurfaceGfxManager::set_mesh_visibility( index_t s, bool b )
    {
        entities_[s]->set_show_mesh( b ) ;
    }
    /*!
     * Sets the mesh surface size to all the elements
     * @param[in] size the size
     */
    void SurfaceGfxManager::set_mesh_size( index_t size )
    {
        for( index_t s = 0; s < entities_.size(); s++ ) {
            set_mesh_size( s, size ) ;
        }
    }
    /*!
     * Sets the mesh surface size
     * @param[in] s the surface index
     * @param[in] size the size
     */
    void SurfaceGfxManager::set_mesh_size( index_t s, index_t size )
    {
        entities_[s]->set_mesh_width( size ) ;
    }

    /*****************************************************************/

    RegionGfxManager::RegionGfxManager( GeoModelGfx& gfx )
        : GeoModelGfxManager( gfx )
    {
    }

    void RegionGfxManager::initialize()
    {
        if( entities_.empty() ) {
            entities_.resize( gfx_.geomodel()->nb_regions(), nil ) ;
            for( index_t e = 0; e < entities_.size(); e++ ) {
                entities_[e] = new RegionGfx( gfx_, gfx_.geomodel()->region( e ) ) ;
            }
        }
    }

    /*!
     * Draws the Regions
     */
    void RegionGfxManager::draw()
    {
        for( index_t m = 0; m < entities_.size(); m++ ) {
            RegionGfx* region = dynamic_cast< RegionGfx* >( entities_[m] ) ;
            if( region->get_vertices_visible() ) region->draw_vertices() ;
            if( region->get_edges_visible() ) region->draw_edges() ;
            if( region->get_surface_visible() ) region->draw_surface() ;
            if( region->get_region_visible() ) region->draw_volume() ;
        }
    }

    /*!
     * Sets the edge color to all the meshes
     * @param[in] r the red component of the color in [0.0, 1.0]
     * @param[in] g the green component of the color in [0.0, 1.0]
     * @param[in] b the blue component of the color in [0.0, 1.0]
     */
    void RegionGfxManager::set_edge_color( float r, float g, float b )
    {
        for( index_t k = 0; k < entities_.size(); k++ ) {
            set_edge_color( k, r, g, b ) ;
        }
    }
    /*!
     * Sets the edge color
     * @param[in] m the mesh index
     * @param[in] r the red component of the color in [0.0, 1.0]
     * @param[in] g the green component of the color in [0.0, 1.0]
     * @param[in] b the blue component of the color in [0.0, 1.0]
     */
    void RegionGfxManager::set_edge_color( index_t m, float r, float g, float b )
    {
        ringmesh_assert( m < entities_.size() ) ;
        entities_[m]->set_mesh_color( r, g, b ) ; //TODO function not good?
    }
    /*!
     * Sets the edge visibility to all the meshes
     * @param[in] b the visibility
     */
    void RegionGfxManager::set_edge_visibility( bool b )
    {
        for( index_t m = 0; m < entities_.size(); m++ ) {
            set_edge_visibility( m, b ) ;
        }
    }
    /*!
     * Sets the edge visibility
     * @param[in] m the mesh index
     * @param[in] b the visibility
     */
    void RegionGfxManager::set_edge_visibility( index_t m, bool b )
    {
        ringmesh_assert( m < entities_.size() ) ;
        dynamic_cast< RegionGfx* >( entities_[m] )->set_edges_visible( b ) ;
    }
    /*!
     * Sets the edge line size to all the lines
     * @param[in] s the size
     */
    void RegionGfxManager::set_edge_size( index_t s )
    {
        for( index_t m = 0; m < entities_.size(); m++ ) {
            set_edge_size( m, s ) ;
        }
    }
    /*!
     * Sets the edge line size
     * @param[in] l the line index
     * @param[in] s the size
     */
    void RegionGfxManager::set_edge_size( index_t l, index_t s )
    {
        ringmesh_assert( l < entities_.size() ) ;
        entities_[l]->set_mesh_width( s ) ;
    }

    /*!
     * Sets the mesh region color to all the regions
     * @param[in] r the red component of the color in [0.0, 1.0]
     * @param[in] g the green component of the color in [0.0, 1.0]
     * @param[in] b the blue component of the color in [0.0, 1.0]
     */
    void RegionGfxManager::set_mesh_color( float r, float g, float b )
    {
        for( index_t m = 0; m < entities_.size(); m++ ) {
            set_mesh_color( m, r, g, b ) ;
        }
    }

    /*!
     * Sets the mesh region color
     * @param[in] m the region index
     * @param[in] r the red component of the color in [0.0, 1.0]
     * @param[in] g the green component of the color in [0.0, 1.0]
     * @param[in] b the blue component of the color in [0.0, 1.0]
     */
    void RegionGfxManager::set_mesh_color( index_t m, float r, float g, float b )
    {
        ringmesh_assert( m < entities_.size() ) ;
        entities_[m]->set_mesh_color( r, g, b ) ;
    }

    /*!
     * Toggles the cell region color per cell type to all the regions
     */
    void RegionGfxManager::set_color_cell_type()
    {
        for( index_t m = 0; m < entities_.size(); m++ ) {
            set_color_cell_type( m ) ;
        }
    }

    /*!
     * Toggles the cell type to all the regions
     */
    void RegionGfxManager::set_draw_cells( GEO::MeshCellType type, bool x )
    {
        for( index_t m = 0; m < entities_.size(); m++ ) {
            set_draw_cells( m, type, x ) ;
        }
    }

    /*!
     * Toggles the cell type display
     */
    void RegionGfxManager::set_draw_cells(
        index_t m,
        GEO::MeshCellType type,
        bool x )
    {
        ringmesh_assert( m < entities_.size() ) ;
        entities_[m]->set_draw_cells( type, x ) ;
    }

    /*!
     * Toggles the cell region color per cell type
     * @param[in] m the region index
     */
    void RegionGfxManager::set_color_cell_type( index_t m )
    {
        ringmesh_assert( m < entities_.size() ) ;
        entities_[m]->set_cells_colors_by_type() ;
    }

    /*!
     * Sets the mesh region visibility to all the regions
     * @param[in] b the visibility
     */
    void RegionGfxManager::set_mesh_visibility( bool b )
    {
        for( index_t m = 0; m < entities_.size(); m++ ) {
            set_mesh_visibility( m, b ) ;
        }
    }

    /*!
     * Sets the mesh region visibility
     * @param[in] m the region index
     * @param[in] b the visibility
     */
    void RegionGfxManager::set_mesh_visibility( index_t m, bool b )
    {
        ringmesh_assert( m < entities_.size() ) ;
        entities_[m]->set_show_mesh( b ) ;
    }

    /*!
     * Sets the mesh region size to all the regions
     * @param[in] s the size
     */
    void RegionGfxManager::set_mesh_size( index_t s )
    {
        for( index_t m = 0; m < entities_.size(); m++ ) {
            set_mesh_size( m, s ) ;
        }
    }

    /*!
     * Sets the mesh region size
     * @param[in] m the region index
     * @param[in] s the size
     */
    void RegionGfxManager::set_mesh_size( index_t m, index_t s )
    {
        ringmesh_assert( m < entities_.size() ) ;
        entities_[m]->set_mesh_width( s ) ;
    }

    /*!
     * Sets the cell region color
     * @param[in] m the region index
     * @param[in] r the red component of the color in [0.0, 1.0]
     * @param[in] g the green component of the color in [0.0, 1.0]
     * @param[in] b the blue component of the color in [0.0, 1.0]
     */
    void RegionGfxManager::set_mesh_element_color(
        index_t m,
        float r,
        float g,
        float b )
    {
        ringmesh_assert( m < entities_.size() ) ;
        entities_[m]->set_cells_color( r, g, b ) ;
    }

    /*!
     * Sets the cell region visibility to all the regions
     * @param[in] m the region index
     * @param[in] b the visibility
     */
    void RegionGfxManager::set_mesh_element_visibility( index_t m, bool b )
    {
        ringmesh_assert( m < entities_.size() ) ;
        dynamic_cast< RegionGfx* >( entities_[m] )->set_region_visible( b ) ;
    }
    void RegionGfxManager::set_cell_type_visibility( GEO::MeshCellType t, bool b )
    {
        for( index_t m = 0; m < entities_.size(); m++ ) {
            set_cell_type_visibility( m, t, b ) ;
        }
    }
    void RegionGfxManager::set_cell_type_visibility(
        index_t m,
        GEO::MeshCellType t,
        bool b )
    {
        ringmesh_assert( m < entities_.size() ) ;
        entities_[m]->set_draw_cells( t, b ) ;
    }

    /*!
     * Sets the cell region shrink to all the regions
     * @param[in] s the shrink
     */
    void RegionGfxManager::set_shrink( double s )
    {
        for( index_t m = 0; m < entities_.size(); m++ ) {
            set_shrink( m, s ) ;
        }
    }

    /*!
     * Sets the cell region shrink
     * @param[in] m the region index
     * @param[in] s the shrink
     */
    void RegionGfxManager::set_shrink( index_t m, double s )
    {
        ringmesh_assert( m < entities_.size() ) ;
        entities_[m]->set_shrink( s ) ;
    }

    /*****************************************************************/

    class AttributeGfx {
    public:
        AttributeGfx( AttributeGfxManager& manager )
            : manager_( manager )
        {
        }
        virtual ~AttributeGfx()
        {
        }

        virtual std::string location_name() const = 0 ;
        void compute_range()
        {
            double attribute_min = max_float64() ;
            double attribute_max = min_float64() ;
            do_compute_range( attribute_min, attribute_max ) ;
            manager_.set_minimum( attribute_min ) ;
            manager_.set_maximum( attribute_max ) ;
        }
        virtual void bind_attribute() = 0 ;
        virtual void unbind_attribute() = 0 ;
        virtual index_t nb_coordinates() = 0 ;

    private:
        virtual void do_compute_range(
            double& attribute_min,
            double& attribute_max ) = 0 ;

    protected:
        AttributeGfxManager& manager_ ;
    } ;

    class CellAttributeGfx: public AttributeGfx {
    public:
        CellAttributeGfx( AttributeGfxManager& manager )
            : AttributeGfx( manager )
        {
        }

        virtual std::string location_name() const override
        {
            return "cells" ;
        }
        virtual void bind_attribute() override
        {
            std::string attribute_name = get_attribute_name_with_coordinate(
                manager_.name(), manager_.coordinate() ) ;
            manager_.gfx().regions.set_scalar_attribute( GEO::MESH_CELLS,
                attribute_name, manager_.minimum(), manager_.maximum(),
                manager_.colormap() ) ;
        }
        virtual void unbind_attribute() override
        {
            manager_.gfx().regions.unset_scalar_attribute() ;
        }
        virtual index_t nb_coordinates() override
        {
            const GeoModel* geomodel = manager_.gfx().geomodel() ;
            GEO::AttributeStore* store =
                geomodel->region( 0 ).cell_attribute_manager().find_attribute_store(
                    manager_.name() ) ;

            if( store == nil ) return 0 ;
            return store->dimension() ;
        }
    private:
        virtual void do_compute_range( double& attribute_min, double& attribute_max ) override
        {
            std::string attribute_name = get_attribute_name_with_coordinate(
                manager_.name(), manager_.coordinate() ) ;
            const GeoModel* geomodel = manager_.gfx().geomodel() ;
            for( index_t r = 0; r < geomodel->nb_regions(); r++ ) {
                GEO::ReadOnlyScalarAttributeAdapter attribute(
                    geomodel->region( r ).cell_attribute_manager(),
                    attribute_name ) ;
                compute_attribute_range( attribute, attribute_min, attribute_max ) ;
            }
        }
    } ;

    class CellVertexAttributeGfx: public AttributeGfx {
    public:
        CellVertexAttributeGfx( AttributeGfxManager& manager )
            : AttributeGfx( manager )
        {
        }

        virtual std::string location_name() const override
        {
            return "cell_vertices" ;
        }
        virtual void bind_attribute() override
        {
            std::string attribute_name = get_attribute_name_with_coordinate(
                manager_.name(), manager_.coordinate() ) ;
            manager_.gfx().regions.set_scalar_attribute( GEO::MESH_VERTICES,
                attribute_name, manager_.minimum(), manager_.maximum(),
                manager_.colormap() ) ;
        }
        virtual void unbind_attribute() override
        {
            manager_.gfx().regions.unset_scalar_attribute() ;
        }
        virtual index_t nb_coordinates() override
        {
            const GeoModel* geomodel = manager_.gfx().geomodel() ;
            GEO::AttributeStore* store =
                geomodel->region( 0 ).vertex_attribute_manager().find_attribute_store(
                    manager_.name() ) ;

            if( store == nil ) return 0 ;
            return store->dimension() ;
        }
    private:
        virtual void do_compute_range( double& attribute_min, double& attribute_max ) override
        {
            std::string attribute_name = get_attribute_name_with_coordinate(
                manager_.name(), manager_.coordinate() ) ;
            const GeoModel* geomodel = manager_.gfx().geomodel() ;
            for( index_t r = 0; r < geomodel->nb_regions(); r++ ) {
                GEO::ReadOnlyScalarAttributeAdapter attribute(
                    geomodel->region( r ).vertex_attribute_manager(),
                    attribute_name ) ;
                compute_attribute_range( attribute, attribute_min, attribute_max ) ;
            }
        }
    } ;

    class FacetAttributeGfx: public AttributeGfx {
    public:
        FacetAttributeGfx( AttributeGfxManager& manager )
            : AttributeGfx( manager )
        {
        }

        virtual std::string location_name() const override
        {
            return "facet" ;
        }
        virtual void bind_attribute() override
        {
            std::string attribute_name = get_attribute_name_with_coordinate(
                manager_.name(), manager_.coordinate() ) ;
            manager_.gfx().surfaces.set_scalar_attribute( GEO::MESH_FACETS,
                attribute_name, manager_.minimum(), manager_.maximum(),
                manager_.colormap() ) ;
        }
        virtual void unbind_attribute() override
        {
            manager_.gfx().surfaces.unset_scalar_attribute() ;
        }
        virtual index_t nb_coordinates() override
        {
            const GeoModel* geomodel = manager_.gfx().geomodel() ;
            GEO::AttributeStore* store =
                geomodel->surface( 0 ).facet_attribute_manager().find_attribute_store(
                    manager_.name() ) ;

            if( store == nil ) return 0 ;
            return store->dimension() ;
        }
    private:
        virtual void do_compute_range( double& attribute_min, double& attribute_max ) override
        {
            std::string attribute_name = get_attribute_name_with_coordinate(
                manager_.name(), manager_.coordinate() ) ;
            const GeoModel* geomodel = manager_.gfx().geomodel() ;
            for( index_t s = 0; s < geomodel->nb_surfaces(); s++ ) {
                GEO::ReadOnlyScalarAttributeAdapter attribute(
                    geomodel->surface( s ).facet_attribute_manager(),
                    attribute_name ) ;
                compute_attribute_range( attribute, attribute_min, attribute_max ) ;
            }
        }
    } ;

    class FacetVertexAttributeGfx: public AttributeGfx {
    public:
        FacetVertexAttributeGfx( AttributeGfxManager& manager )
            : AttributeGfx( manager )
        {
        }

        virtual std::string location_name() const override
        {
            return "facet_vertices" ;
        }
        virtual void bind_attribute() override
        {
            std::string attribute_name = get_attribute_name_with_coordinate(
                manager_.name(), manager_.coordinate() ) ;
            manager_.gfx().surfaces.set_scalar_attribute( GEO::MESH_VERTICES,
                attribute_name, manager_.minimum(), manager_.maximum(),
                manager_.colormap() ) ;
        }
        virtual void unbind_attribute() override
        {
            manager_.gfx().surfaces.unset_scalar_attribute() ;
        }
        virtual index_t nb_coordinates() override
        {
            const GeoModel* geomodel = manager_.gfx().geomodel() ;
            GEO::AttributeStore* store =
                geomodel->surface( 0 ).vertex_attribute_manager().find_attribute_store(
                    manager_.name() ) ;

            if( store == nil ) return 0 ;
            return store->dimension() ;
        }
    private:
        virtual void do_compute_range( double& attribute_min, double& attribute_max ) override
        {
            std::string attribute_name = get_attribute_name_with_coordinate(
                manager_.name(), manager_.coordinate() ) ;
            const GeoModel* geomodel = manager_.gfx().geomodel() ;
            for( index_t s = 0; s < geomodel->nb_surfaces(); s++ ) {
                GEO::ReadOnlyScalarAttributeAdapter attribute(
                    geomodel->surface( s ).vertex_attribute_manager(),
                    attribute_name ) ;
                compute_attribute_range( attribute, attribute_min, attribute_max ) ;
            }
        }
    } ;

    AttributeGfxManager::AttributeGfxManager( GeoModelGfx& gfx )
        :
            gfx_( gfx ),
            location_( nb_locations ),
            coordinate_( 0 ),
            colormap_texture_( 0 ),
            minimum_( 0.0 ),
            maximum_( 0.0 )
    {
        attributes_[facets] = new FacetAttributeGfx( *this ) ;
        attributes_[facet_vertices] = new FacetVertexAttributeGfx( *this ) ;
        attributes_[cells] = new CellAttributeGfx( *this ) ;
        attributes_[cell_vertices] = new CellVertexAttributeGfx( *this ) ;
    }

    AttributeGfxManager::~AttributeGfxManager()
    {
        for( index_t i = 0; i < nb_locations; i++ ) {
            delete attributes_[i] ;
        }
    }

    std::string AttributeGfxManager::location_name( Attribute_location location )
    {
        if( location == nb_locations )
            return "location" ;
        else
            return attributes_[location]->location_name() ;
    }

    void AttributeGfxManager::compute_range()
    {
        if( location() < nb_locations ) {
            attributes_[location()]->compute_range() ;
        }
    }

    void AttributeGfxManager::bind_attribute()
    {
        if( location() < nb_locations ) {
            attributes_[location()]->bind_attribute() ;
        }
    }

    void AttributeGfxManager::unbind_attribute()
    {
        if( location() < nb_locations ) {
            attributes_[location()]->unbind_attribute() ;
        }
    }

    index_t AttributeGfxManager::nb_coordinates() const
    {
        if( location() < nb_locations ) {
            return attributes_[location()]->nb_coordinates() ;
        }
        return 0 ;
    }

    /*****************************************************************/

    GeoModelGfx::GeoModelGfx()
        :
            geomodel_( nil ),
            corners( *this ),
            lines( *this ),
            surfaces( *this ),
            regions( *this ),
            attribute( *this )
    {
    }

    GeoModelGfx::~GeoModelGfx()
    {
    }

    /*!
     * Sets the GeoModel associated to the graphics
     * @param[in] geomodel the GeoModel
     */
    void GeoModelGfx::set_geomodel( const GeoModel& geomodel )
    {
        geomodel_ = &geomodel ;
        initialize() ;
    }

    /*!
     * Gets the GeoModel associated to the graphics
     * @return the GeoModel
     */
    const GeoModel* GeoModelGfx::geomodel() const
    {
        return geomodel_ ;
    }

    /*!
     * Initializes the database according the GeoModel dimensions
     */
    void GeoModelGfx::initialize()
    {
        ringmesh_assert( geomodel_ ) ;
        corners.initialize() ;
        lines.initialize() ;
        surfaces.initialize() ;
        regions.initialize() ;
    }

    void GeoModelGfx::need_to_update()
    {
        corners.need_to_update() ;
        lines.need_to_update() ;
        surfaces.need_to_update() ;
        regions.need_to_update() ;
    }
}

#endif
