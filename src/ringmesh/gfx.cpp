/*
 * Copyright (c) 2012-2016, Association Scientifique pour la Geologie et ses Applications (ASGA)
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

#include <ringmesh/gfx.h>

#ifdef RINGMESH_WITH_GRAPHICS

#include <geogram/basic/file_system.h>
#include <geogram/basic/logger.h>
#include <geogram/mesh/mesh_geometry.h>

#include <geogram_gfx/basic/GL.h>
#include <geogram_gfx/glup_viewer/glup_viewer.h>
#include <geogram_gfx/third_party/quicktext/glQuickText.h>

#include <ringmesh/command_line.h>
#include <ringmesh/geo_model_entity.h>
#include <ringmesh/io.h>

#define define_color( name, r, g, b )\
    class name: public GetColor {\
    public:\
        virtual Color get_color() {\
            return Color( r, g, b ) ;\
        }\
    }; \
    ringmesh_register_color_creator( name, #name ) \

namespace {
    using namespace RINGMesh ;

    std::string path_to_label(
        const std::string& viewer_path,
        const std::string& path )
    {
        if( GEO::String::string_starts_with( path, viewer_path ) ) {
            return path.substr( viewer_path.length(),
                path.length() - viewer_path.length() ) ;
        }
        return path ;
    }

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
        virtual void draw_edges()
        {
            index_t w = get_mesh_width() ;
            set_mesh_width( w + 1 ) ;
            GEO::MeshGfx::draw_edges() ;
            set_mesh_width( w ) ;
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

        virtual void draw_surface()
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
        entities_[e]->set_points_size( s ) ;
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
            entities_.resize( gfx_.geo_model()->nb_corners(), nil ) ;
            for( index_t e = 0; e < entities_.size(); e++ ) {
                entities_[e] = new CornerGfx( gfx_, gfx_.geo_model()->corner( e ) ) ;
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
            entities_.resize( gfx_.geo_model()->nb_lines(), nil ) ;
            for( index_t e = 0; e < entities_.size(); e++ ) {
                entities_[e] = new LineGfx( gfx_, gfx_.geo_model()->line( e ) ) ;
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
            entities_.resize( gfx_.geo_model()->nb_surfaces(), nil ) ;
            for( index_t e = 0; e < entities_.size(); e++ ) {
                entities_[e] = new SurfaceGfx( gfx_,
                    gfx_.geo_model()->surface( e ) ) ;
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
            entities_.resize( gfx_.geo_model()->nb_regions(), nil ) ;
            for( index_t e = 0; e < entities_.size(); e++ ) {
                entities_[e] = new RegionGfx( gfx_, gfx_.geo_model()->region( e ) ) ;
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

        virtual std::string location_name() const
        {
            return "cells" ;
        }
        virtual void bind_attribute()
        {
            std::string attribute_name = get_attribute_name_with_coordinate(
                manager_.name(), manager_.coordinate() ) ;
            const GeoModel* model = manager_.gfx().geo_model() ;
            for( index_t r = 0; r < model->nb_regions(); r++ ) {
                manager_.gfx().regions.set_scalar_attribute( GEO::MESH_CELLS,
                    attribute_name, manager_.minimum(), manager_.maximum(),
                    manager_.colormap() ) ;
            }
        }
        virtual void unbind_attribute()
        {
            const GeoModel* model = manager_.gfx().geo_model() ;
            for( index_t r = 0; r < model->nb_regions(); r++ ) {
                manager_.gfx().regions.unset_scalar_attribute() ;
            }
        }
        virtual index_t nb_coordinates()
        {
            const GeoModel* model = manager_.gfx().geo_model() ;
            GEO::AttributeStore* store =
                model->region( 0 ).cell_attribute_manager().find_attribute_store(
                    manager_.name() ) ;

            if( store == nil ) return 0 ;
            return store->dimension() ;
        }
    private:
        virtual void do_compute_range( double& attribute_min, double& attribute_max )
        {
            std::string attribute_name = get_attribute_name_with_coordinate(
                manager_.name(), manager_.coordinate() ) ;
            const GeoModel* model = manager_.gfx().geo_model() ;
            for( index_t r = 0; r < model->nb_regions(); r++ ) {
                GEO::ReadOnlyScalarAttributeAdapter attribute(
                    model->region( r ).cell_attribute_manager(), attribute_name ) ;
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

        virtual std::string location_name() const
        {
            return "cell_vertices" ;
        }
        virtual void bind_attribute()
        {
            std::string attribute_name = get_attribute_name_with_coordinate(
                manager_.name(), manager_.coordinate() ) ;
            const GeoModel* model = manager_.gfx().geo_model() ;
            for( index_t r = 0; r < model->nb_regions(); r++ ) {
                manager_.gfx().regions.set_scalar_attribute( GEO::MESH_VERTICES,
                    attribute_name, manager_.minimum(), manager_.maximum(),
                    manager_.colormap() ) ;
            }
        }
        virtual void unbind_attribute()
        {
            const GeoModel* model = manager_.gfx().geo_model() ;
            for( index_t r = 0; r < model->nb_regions(); r++ ) {
                manager_.gfx().regions.unset_scalar_attribute() ;
            }
        }
        virtual index_t nb_coordinates()
        {
            const GeoModel* model = manager_.gfx().geo_model() ;
            GEO::AttributeStore* store =
                model->region( 0 ).vertex_attribute_manager().find_attribute_store(
                    manager_.name() ) ;

            if( store == nil ) return 0 ;
            return store->dimension() ;
        }
    private:
        virtual void do_compute_range( double& attribute_min, double& attribute_max )
        {
            std::string attribute_name = get_attribute_name_with_coordinate(
                manager_.name(), manager_.coordinate() ) ;
            const GeoModel* model = manager_.gfx().geo_model() ;
            for( index_t r = 0; r < model->nb_regions(); r++ ) {
                GEO::ReadOnlyScalarAttributeAdapter attribute(
                    model->region( r ).vertex_attribute_manager(), attribute_name ) ;
                compute_attribute_range( attribute, attribute_min, attribute_max ) ;
            }
        }
    } ;

    AttributeGfxManager::AttributeGfxManager( GeoModelGfx& gfx )
        :
            gfx_( gfx ),
            location_( nb_locations ),
            coordinate_( 0 ),
            minimum_( 0.0 ),
            maximum_( 0.0 )
    {
        attributes_[cells] = new CellAttributeGfx( *this ) ;
        attributes_[cell_vertices] = new CellVertexAttributeGfx( *this ) ;
    }

    AttributeGfxManager::~AttributeGfxManager()
    {
        delete attributes_[cells] ;
        delete attributes_[cell_vertices] ;
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
            model_( nil ),
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
     * @param[in] model the GeoModel
     */
    void GeoModelGfx::set_geo_model( const GeoModel& model )
    {
        model_ = &model ;
        initialize() ;
    }

    /*!
     * Gets the GeoModel associated to the graphics
     * @return the GeoModel
     */
    const GeoModel* GeoModelGfx::geo_model() const
    {
        return model_ ;
    }

    /*!
     * Initializes the database according the GeoModel dimensions
     */
    void GeoModelGfx::initialize()
    {
        ringmesh_assert( model_ ) ;
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

    /*****************************************************************/

    RINGMeshApplication::GeoModelViewer::GeoModelViewer(
        RINGMeshApplication& app,
        const std::string& filename )
        : app_( app )
    {
        is_visible_ = true ;

        show_corners_ = true ;
        show_lines_ = true ;
        show_surface_ = true ;
        show_volume_ = false ;
        colored_cells_ = false ;
        show_voi_ = false ;
        show_colored_regions_ = false ;
        show_colored_layers_ = false ;
        show_colormap_ = false ;

        show_hex_ = true ;
        show_prism_ = true ;
        show_pyramid_ = true ;
        show_tetra_ = true ;

        shrink_ = 0.0 ;
        mesh_visible_ = true ;
        meshed_regions_ = false ;

        show_attributes_ = false ;
        attribute_min_ = 0.0f ;
        attribute_max_ = 0.0f ;
        reset_attribute_name() ;

        geomodel_load( GM_, filename ) ;
        for( GEO::index_t s = 0; s < GM_.nb_surfaces(); s++ ) {
            const RINGMesh::Surface& S = GM_.surface( s ) ;
            for( GEO::index_t v = 0; v < S.nb_vertices(); ++v ) {
                const vec3& p = S.vertex( v ) ;
                bbox_.add_point( p ) ;
            }
        }

        meshed_regions_ = GM_.region( 0 ).is_meshed() ;
        if( meshed_regions_ ) {
            show_volume_ = true ;
        }
        GM_gfx_.set_geo_model( GM_ ) ;
        if( !app.colormaps_.empty() ) {
            GM_gfx_.attribute.set_colormap( app.colormaps_[0].texture ) ;
        }
    }

    void RINGMeshApplication::GeoModelViewer::reset_attribute_name()
    {
        GM_gfx_.attribute.set_name( "name" ) ;
    }

    void RINGMeshApplication::GeoModelViewer::toggle_colored_cells()
    {
        show_colored_regions_.new_status = false ;
        show_colored_layers_.new_status = false ;
        GM_gfx_.regions.set_color_cell_type() ;
    }
    void RINGMeshApplication::GeoModelViewer::toggle_colored_regions()
    {
        colored_cells_.new_status = false ;
        show_colored_layers_.new_status = false ;
        for( GEO::index_t r = 0; r < GM_.nb_regions(); r++ ) {
            GM_gfx_.regions.set_mesh_element_color( r,
                std::fmod( GEO::Numeric::random_float32(), 1 ),
                std::fmod( GEO::Numeric::random_float32(), 1 ),
                std::fmod( GEO::Numeric::random_float32(), 1 ) ) ;
        }
    }

    void RINGMeshApplication::GeoModelViewer::toggle_colored_layers()
    {
        colored_cells_.new_status = false ;
        show_colored_regions_.new_status = false ;
        for( GEO::index_t l = 0; l < GM_.nb_layers(); l++ ) {
            float red = std::fmod( GEO::Numeric::random_float32(), 1 ) ;
            float green = std::fmod( GEO::Numeric::random_float32(), 1 ) ;
            float blue = std::fmod( GEO::Numeric::random_float32(), 1 ) ;
            const GeoModelEntity& cur_layer = GM_.layer( l ) ;
            for( index_t r = 0; r < cur_layer.nb_children(); ++r )
                GM_gfx_.regions.set_mesh_element_color( cur_layer.child( r ).index(),
                    red, green, blue ) ;
        }
    }

    void RINGMeshApplication::GeoModelViewer::draw_scene()
    {
        GM_gfx_.surfaces.set_mesh_visibility( mesh_visible_ ) ;
        GM_gfx_.regions.set_mesh_visibility( mesh_visible_ ) ;

        if( show_attributes_ ) {
            GM_gfx_.attribute.bind_attribute() ;
        } else {
            GM_gfx_.attribute.unbind_attribute() ;
        }

        if( show_corners_ ) {
            GM_gfx_.corners.draw() ;
        }

        if( show_lines_ ) {
            if( app_.white_bg_ ) {
                GM_gfx_.lines.GeoModelGfxManager::set_mesh_element_color( 0.0f, 0.0f,
                    0.0f ) ;
            } else {
                GM_gfx_.lines.GeoModelGfxManager::set_mesh_element_color( 1.0f, 1.0f,
                    1.0f ) ;
            }
            GM_gfx_.lines.draw() ;
        }

        if( show_surface_ ) {
            if( app_.white_bg_ ) {
                GM_gfx_.surfaces.GeoModelGfxManager::set_mesh_element_color( 0.9f,
                    0.9f, 0.9f ) ;
                GM_gfx_.surfaces.set_mesh_color( 0.0f, 0.0f, 0.0f ) ;
            } else {
                GM_gfx_.surfaces.GeoModelGfxManager::set_mesh_element_color( 0.1f,
                    0.1f, 0.1f ) ;
                GM_gfx_.surfaces.set_mesh_color( 1.0f, 1.0f, 1.0f ) ;
            }
            for( GEO::index_t s = 0; s < GM_.nb_surfaces(); s++ ) {
                if( GM_.surface( s ).is_on_voi() ) {
                    GM_gfx_.surfaces.set_mesh_element_visibility( s, show_voi_ ) ;
                }
            }
            GM_gfx_.surfaces.draw() ;
        }

        if( show_volume_ && meshed_regions_ ) {
            if( colored_cells_.need_to_update() ) {
                colored_cells_.update() ;
                if( colored_cells_.new_status ) {
                    toggle_colored_cells() ;
                }
            } else if( show_colored_regions_.need_to_update() ) {
                show_colored_regions_.update() ;
                if( show_colored_regions_.new_status ) {
                    toggle_colored_regions() ;
                }
            } else if( show_colored_layers_.need_to_update() ) {
                show_colored_layers_.update() ;
                if( show_colored_layers_.new_status ) {
                    toggle_colored_layers() ;
                }
            }
            if( !colored_cells_.new_status && !show_colored_regions_.new_status
                && !show_colored_layers_.new_status ) {
                colored_cells_.update() ;
                show_colored_regions_.update() ;
                show_colored_layers_.update() ;
                if( app_.white_bg_ ) {
                    GM_gfx_.regions.GeoModelGfxManager::set_mesh_element_color( 0.9f,
                        0.9f, 0.9f ) ;
                    GM_gfx_.regions.set_mesh_color( 0.0f, 0.0f, 0.0f ) ;
                } else {
                    GM_gfx_.regions.GeoModelGfxManager::set_mesh_element_color( 0.1f,
                        0.1f, 0.1f ) ;
                    GM_gfx_.regions.set_mesh_color( 1.0f, 1.0f, 1.0f ) ;
                }
            }
            GM_gfx_.regions.set_draw_cells( GEO::MESH_HEX, show_hex_ ) ;
            GM_gfx_.regions.set_draw_cells( GEO::MESH_PRISM, show_prism_ ) ;
            GM_gfx_.regions.set_draw_cells( GEO::MESH_PYRAMID, show_pyramid_ ) ;
            GM_gfx_.regions.set_draw_cells( GEO::MESH_TET, show_tetra_ ) ;
            GM_gfx_.regions.set_shrink( shrink_ ) ;
            GM_gfx_.regions.draw() ;
        }
    }

    void RINGMeshApplication::GeoModelViewer::set_attribute_names(
        const GEO::AttributesManager& attributes )
    {
        GEO::vector< std::string > attribute_names ;
        attributes.list_attribute_names( attribute_names ) ;
        for( index_t i = 0; i < attribute_names.size(); ++i ) {
            const GEO::AttributeStore* store = attributes.find_attribute_store(
                attribute_names[i] ) ;
            if( GEO::ReadOnlyScalarAttributeAdapter::can_be_bound_to( store ) ) {
                if( ImGui::Button( attribute_names[i].c_str() ) ) {
                    GM_gfx_.attribute.set_name( attribute_names[i] ) ;
                    GM_gfx_.attribute.set_coordinate( 0 ) ;
                    autorange() ;
                    ImGui::CloseCurrentPopup() ;
                }
            }
        }
    }

    void RINGMeshApplication::GeoModelViewer::autorange()
    {
        GM_gfx_.attribute.compute_range() ;
        attribute_max_ = GM_gfx_.attribute.maximum() ;
        attribute_min_ = GM_gfx_.attribute.minimum() ;
    }

    void RINGMeshApplication::GeoModelViewer::draw_object_properties()
    {
        ImGui::Checkbox( "Attributes", &show_attributes_ ) ;
        if( show_attributes_ ) {
            if( ImGui::Button(
                GM_gfx_.attribute.location_name( GM_gfx_.attribute.location() ).c_str(),
                ImVec2( -1, 0 ) ) ) {
                ImGui::OpenPopup( "##Locations" ) ;
            }
            if( ImGui::BeginPopup( "##Locations" ) ) {
                for( index_t i = 0; i < AttributeGfxManager::nb_locations; i++ ) {
                    AttributeGfxManager::Attribute_location l =
                        static_cast< AttributeGfxManager::Attribute_location >( i ) ;
                    if( ImGui::Button(
                        GM_gfx_.attribute.location_name( l ).c_str() ) ) {
                        GM_gfx_.attribute.set_location( l ) ;
                        reset_attribute_name() ;
                        ImGui::CloseCurrentPopup() ;
                    }
                }
                ImGui::EndPopup() ;
            }

            if( ImGui::Button( GM_gfx_.attribute.name().c_str(),
                ImVec2( -1, 0 ) ) ) {
                ImGui::OpenPopup( "##Attributes" ) ;
            }
            if( ImGui::BeginPopup( "##Attributes" ) ) {
                switch( GM_gfx_.attribute.location() ) {
                    case AttributeGfxManager::cells:
                        set_attribute_names(
                            GM_.region( 0 ).cell_attribute_manager() ) ;
                        break ;
                    case AttributeGfxManager::cell_vertices:
                        set_attribute_names(
                            GM_.region( 0 ).vertex_attribute_manager() ) ;
                        break ;
                }
                ImGui::EndPopup() ;
            }
            if( GM_gfx_.attribute.location() != AttributeGfxManager::nb_locations
                && GM_gfx_.attribute.nb_coordinates() > 1 ) {
                if( ImGui::Button(
                    GEO::String::to_string( GM_gfx_.attribute.coordinate() ).c_str(),
                    ImVec2( -1, 0 ) ) ) {
                    ImGui::OpenPopup( "##Coordinates" ) ;
                }
                if( ImGui::BeginPopup( "##Coordinates" ) ) {
                    for( index_t i = 0; i < GM_gfx_.attribute.nb_coordinates();
                        i++ ) {
                        if( ImGui::Button( GEO::String::to_string( i ).c_str() ) ) {
                            GM_gfx_.attribute.set_coordinate( i ) ;
                            autorange() ;
                            ImGui::CloseCurrentPopup() ;
                        }
                    }
                    ImGui::EndPopup() ;
                }
            }
            if( ImGui::InputFloat( "min", &attribute_min_ ) ) {
                GM_gfx_.attribute.set_minimum( attribute_min_ ) ;
            }
            if( ImGui::InputFloat( "max", &attribute_max_ ) ) {
                GM_gfx_.attribute.set_maximum( attribute_max_ ) ;
            }
            if( ImGui::Button( "autorange", ImVec2( -1, 0 ) ) ) {
                autorange() ;
            }
            if( ImGui::ImageButton(
                app_.convert_to_ImTextureID( GM_gfx_.attribute.colormap() ),
                ImVec2( 115, 8 ) ) ) {
                ImGui::OpenPopup( "##Colormap" ) ;
            }
            if( ImGui::BeginPopup( "##Colormap" ) ) {
                for( index_t i = 0; i < app_.colormaps_.size(); ++i ) {
                    if( ImGui::ImageButton(
                        app_.convert_to_ImTextureID( app_.colormaps_[i].texture ),
                        ImVec2( 100, 8 ) ) ) {
                        GM_gfx_.attribute.set_colormap(
                            app_.colormaps_[i].texture ) ;
                        ImGui::CloseCurrentPopup() ;
                    }
                }
                ImGui::EndPopup() ;
            }
            ImGui::Checkbox( "Colormap [M]", &show_colormap_ ) ;
        }

        ImGui::Separator() ;
        ImGui::Checkbox( "VOI [V]", &show_voi_ ) ;
        ImGui::Checkbox( "Mesh [m]", &mesh_visible_ ) ;

        ImGui::Separator() ;
        ImGui::Checkbox( "Corner [c]", &show_corners_ ) ;
        ImGui::Checkbox( "Line [e]", &show_lines_ ) ;
        ImGui::Checkbox( "Surface [s]", &show_surface_ ) ;

        if( meshed_regions_ ) {
            ImGui::Separator() ;
            ImGui::Checkbox( "Region [v]", &show_volume_ ) ;
            if( show_volume_ ) {
                ImGui::Checkbox( "Col. cells [C]", &colored_cells_.new_status ) ;
                ImGui::Checkbox( "Col. regions [r]",
                    &show_colored_regions_.new_status ) ;
                ImGui::Checkbox( "Col. layers [R]",
                    &show_colored_layers_.new_status ) ;
                ImGui::SliderFloat( "Shrk.", &shrink_, 0.0f, 1.0f, "%.1f" ) ;
                ImGui::Checkbox( "Hex", &show_hex_ ) ;
                ImGui::Checkbox( "Prism", &show_prism_ ) ;
                ImGui::Checkbox( "Pyramid", &show_pyramid_ ) ;
                ImGui::Checkbox( "Tetra", &show_tetra_ ) ;
            }
        }
    }

    void RINGMeshApplication::GeoModelViewer::draw_colormap()
    {
        GLUPboolean clipping_save = glupIsEnabled( GLUP_CLIPPING ) ;
        glupDisable( GLUP_CLIPPING ) ;

        glupMatrixMode( GLUP_TEXTURE_MATRIX ) ;
        glupLoadIdentity() ;

        glupMatrixMode( GLUP_PROJECTION_MATRIX ) ;
        glupPushMatrix() ;
        glupLoadIdentity() ;

        glupMatrixMode( GLUP_MODELVIEW_MATRIX ) ;
        glupPushMatrix() ;
        glupLoadIdentity() ;

        const float z = -1.0f ;
        const float w = 0.3 ;
        const float h = 0.1 ;
        const float x1 = 0. ;
        const float y1 = -0.9 ;
        const float tmin = float( GM_gfx_.attribute.minimum() ) ;
        const float tmax = float( GM_gfx_.attribute.maximum() ) ;
        GEO::glupMapTexCoords1d( tmin, tmax, 1. ) ;

        glupColor3f( 1.0f, 1.0f, 1.0f ) ;
        glupDisable( GLUP_LIGHTING ) ;
        glupEnable( GLUP_TEXTURING ) ;
        glupTextureMode( GLUP_TEXTURE_REPLACE ) ;
        glupTextureType( GLUP_TEXTURE_1D ) ;
        glupEnable( GLUP_DRAW_MESH ) ;
        glupSetColor3f( GLUP_MESH_COLOR, 0.0f, 0.0f, 0.0f ) ;
        glupSetMeshWidth( 2 ) ;
        glupSetCellsShrink( 0.0f ) ;

        glupBegin( GLUP_QUADS ) ;
        glupTexCoord1f( tmin ) ;
        glupVertex3f( x1 - w, y1, z ) ;
        glupTexCoord1f( tmax ) ;
        glupVertex3f( x1 + w, y1, z ) ;
        glupTexCoord1f( tmax ) ;
        glupVertex3f( x1 + w, y1 + h, z ) ;
        glupTexCoord1f( tmin ) ;
        glupVertex3f( x1 - w, y1 + h, z ) ;
        glupEnd() ;

        glupTextureType( GLUP_TEXTURE_2D ) ;
        glupMatrixMode( GLUP_TEXTURE_MATRIX ) ;
        glupLoadIdentity() ;
        glupMatrixMode( GLUP_MODELVIEW_MATRIX ) ;

        glupSetColor4f( GLUP_FRONT_AND_BACK_COLOR, 0.0f, 0.0f, 0.0f, 1.0f ) ;

        const float font_sz = 0.003f ;
        const float font_height = 0.4f
            * float( glQuickText::getFontHeight( font_sz ) ) ;

        glQuickText::printfAt( x1 - w, y1 - font_height, z, font_sz,
            GEO::String::to_string( GM_gfx_.attribute.minimum() ).c_str() ) ;

        glQuickText::printfAt( x1 + w, y1 - font_height, z, font_sz,
            GEO::String::to_string( GM_gfx_.attribute.maximum() ).c_str() ) ;

        glupMatrixMode( GLUP_PROJECTION_MATRIX ) ;
        glupPopMatrix() ;

        glupMatrixMode( GLUP_MODELVIEW_MATRIX ) ;
        glupPopMatrix() ;

        if( clipping_save ) {
            glupEnable( GLUP_CLIPPING ) ;
        }
    }

    /*****************************************************************/

    RINGMeshApplication::MeshViewer::MeshViewer(
        RINGMeshApplication& app,
        const std::string& filename )
        : app_( app )
    {
        is_visible_ = true ;

        show_vertices_ = false ;
        vertices_size_ = 1.0f ;

        show_surface_ = true ;
        show_surface_colors_ = true ;
        show_mesh_ = true ;
        show_surface_borders_ = false ;

        show_volume_ = false ;
        cells_shrink_ = 0.0f ;
        show_colored_cells_ = false ;
        show_hexes_ = true ;

        show_attributes_ = false ;
        current_colormap_texture_ = 0 ;
        attribute_min_ = 0.0f ;
        attribute_max_ = 0.0f ;
        attribute_ = "vertices.point_fp32[0]" ;
        attribute_name_ = "point_fp32[0]" ;
        attribute_subelements_ = GEO::MESH_VERTICES ;

        GEO::mesh_load( filename, mesh_ ) ;
        name_ = GEO::FileSystem::base_name( filename, true ) ;
        mesh_gfx_.set_mesh( &mesh_ ) ;

        for( index_t v = 0; v < mesh_.vertices.nb(); v++ ) {
            bbox_.add_point( mesh_.vertices.point_ptr( v ) ) ;
        }
    }

    void RINGMeshApplication::MeshViewer::draw_object_properties()
    {
        ImGui::Checkbox( "attributes", &show_attributes_ ) ;
        if( show_attributes_ ) {
            if( attribute_min_ == 0.0f && attribute_max_ == 0.0f ) {
                autorange() ;
            }
            if( ImGui::Button( ( attribute_ + "##Attribute" ).c_str(),
                ImVec2( -1, 0 ) ) ) {
                ImGui::OpenPopup( "##Attributes" ) ;
            }
            if( ImGui::BeginPopup( "##Attributes" ) ) {
                std::vector< std::string > attributes ;
                GEO::String::split_string( attribute_names(), ';', attributes ) ;
                for( index_t i = 0; i < attributes.size(); ++i ) {
                    if( ImGui::Button( attributes[i].c_str() ) ) {
                        set_attribute( attributes[i] ) ;
                        ImGui::CloseCurrentPopup() ;
                    }
                }
                ImGui::EndPopup() ;
            }
            ImGui::InputFloat( "min", &attribute_min_ ) ;
            ImGui::InputFloat( "max", &attribute_max_ ) ;
            if( ImGui::Button( "autorange", ImVec2( -1, 0 ) ) ) {
                autorange() ;
            }
            if( ImGui::ImageButton(
                app_.convert_to_ImTextureID( current_colormap_texture_ ),
                ImVec2( 115, 8 ) ) ) {
                ImGui::OpenPopup( "##Colormap" ) ;
            }
            if( ImGui::BeginPopup( "##Colormap" ) ) {
                for( index_t i = 0; i < app_.colormaps_.size(); ++i ) {
                    if( ImGui::ImageButton(
                        app_.convert_to_ImTextureID( app_.colormaps_[i].texture ),
                        ImVec2( 100, 8 ) ) ) {
                        current_colormap_texture_ = app_.colormaps_[i].texture ;
                        ImGui::CloseCurrentPopup() ;
                    }
                }
                ImGui::EndPopup() ;
            }
        }

        ImGui::Separator() ;
        ImGui::Checkbox( "Vertices [p]", &show_vertices_ ) ;
        if( show_vertices_ ) {
            ImGui::SliderFloat( "sz.", &vertices_size_, 0.1f, 5.0f, "%.1f" ) ;
        }

        if( mesh_.facets.nb() != 0 ) {
            ImGui::Separator() ;
            ImGui::Checkbox( "Surface [S]", &show_surface_ ) ;
            if( show_surface_ ) {
                ImGui::Checkbox( "colors [c]", &show_surface_colors_ ) ;
                ImGui::Checkbox( "mesh [m]", &show_mesh_ ) ;
                ImGui::Checkbox( "borders [B]", &show_surface_borders_ ) ;
            }
        }

        if( mesh_.cells.nb() != 0 ) {
            ImGui::Separator() ;
            ImGui::Checkbox( "Volume [V]", &show_volume_ ) ;
            if( show_volume_ ) {
                ImGui::SliderFloat( "shrk.", &cells_shrink_, 0.0f, 1.0f, "%.2f" ) ;
                if( !mesh_.cells.are_simplices() ) {
                    ImGui::Checkbox( "colored cells [C]", &show_colored_cells_ ) ;
                    ImGui::Checkbox( "hexes [j]", &show_hexes_ ) ;
                }
            }
        }
    }

    void RINGMeshApplication::MeshViewer::draw_scene()
    {
        mesh_gfx_.set_lighting( app_.lighting_ ) ;

        if( show_attributes_ ) {
            mesh_gfx_.set_scalar_attribute( attribute_subelements_, attribute_name_,
                double( attribute_min_ ), double( attribute_max_ ),
                current_colormap_texture_, 1 ) ;
        } else {
            mesh_gfx_.unset_scalar_attribute() ;
        }

        if( show_vertices_ ) {
            mesh_gfx_.set_points_size( vertices_size_ ) ;
            mesh_gfx_.draw_vertices() ;
        }

        if( app_.white_bg_ ) {
            mesh_gfx_.set_mesh_color( 0.0, 0.0, 0.0 ) ;
        } else {
            mesh_gfx_.set_mesh_color( 1.0, 1.0, 1.0 ) ;
        }

        if( show_surface_colors_ ) {
            if( mesh_.cells.nb() == 0 ) {
                mesh_gfx_.set_surface_color( 0.5f, 0.75f, 1.0f ) ;
                mesh_gfx_.set_backface_surface_color( 1.0f, 0.0f, 0.0f ) ;
            } else {
                mesh_gfx_.set_surface_color( 0.7f, 0.0f, 0.0f ) ;
                mesh_gfx_.set_backface_surface_color( 1.0f, 1.0f, 0.0f ) ;
            }
        } else {
            if( app_.white_bg_ ) {
                mesh_gfx_.set_surface_color( 0.9f, 0.9f, 0.9f ) ;
            } else {
                mesh_gfx_.set_surface_color( 0.1f, 0.1f, 0.1f ) ;
            }
        }

        mesh_gfx_.set_show_mesh( show_mesh_ ) ;

        if( show_surface_ ) {
            mesh_gfx_.draw_surface() ;
        }

        if( show_surface_borders_ ) {
            mesh_gfx_.draw_surface_borders() ;
        }

        if( show_mesh_ ) {
            mesh_gfx_.draw_edges() ;
        }

        if( show_volume_ ) {
            if( glupIsEnabled( GLUP_CLIPPING )
                && glupGetClipMode() == GLUP_CLIP_SLICE_CELLS ) {
                mesh_gfx_.set_lighting( false ) ;
            }

            mesh_gfx_.set_shrink( double( cells_shrink_ ) ) ;
            mesh_gfx_.set_draw_cells( GEO::MESH_HEX, show_hexes_ ) ;
            if( show_colored_cells_ ) {
                mesh_gfx_.set_cells_colors_by_type() ;
            } else {
                mesh_gfx_.set_cells_color( 0.9f, 0.9f, 0.9f ) ;
            }
            mesh_gfx_.draw_volume() ;

            mesh_gfx_.set_lighting( app_.lighting_ ) ;
        }
    }

    void RINGMeshApplication::MeshViewer::autorange()
    {
        if( attribute_subelements_ != GEO::MESH_NONE ) {
            attribute_min_ = 0.0 ;
            attribute_max_ = 0.0 ;
            const GEO::MeshSubElementsStore& subelements = mesh_.get_subelements_by_type(
                attribute_subelements_ ) ;
            GEO::ReadOnlyScalarAttributeAdapter attribute( subelements.attributes(),
                attribute_name_ ) ;
            if( attribute.is_bound() ) {
                attribute_min_ = GEO::Numeric::max_float32() ;
                attribute_max_ = GEO::Numeric::min_float32() ;
                for( index_t i = 0; i < subelements.nb(); ++i ) {
                    attribute_min_ = GEO::geo_min( attribute_min_,
                        float( attribute[i] ) ) ;
                    attribute_max_ = GEO::geo_max( attribute_max_,
                        float( attribute[i] ) ) ;
                }
            }
        }
    }

    std::string RINGMeshApplication::MeshViewer::attribute_names()
    {
        return mesh_.get_scalar_attributes() ;
    }

    void RINGMeshApplication::MeshViewer::set_attribute(
        const std::string& attribute )
    {
        attribute_ = attribute ;
        std::string subelements_name ;
        GEO::String::split_string( attribute_, '.', subelements_name, attribute_name_ ) ;
        attribute_subelements_ = mesh_.name_to_subelements_type( subelements_name ) ;
        if( attribute_min_ == 0.0f && attribute_max_ == 0.0f ) {
            autorange() ;
        }
    }

    /*****************************************************************/

    RINGMeshApplication::RINGMeshApplication( int argc, char** argv )
        :
            GEO::Application( argc, argv, "<filename>" ),
            current_viewer_( NO_ID ),
            current_viewer_type_( NONE )
    {
        GEO::CmdLine::declare_arg( "attributes", true, "load mesh attributes" ) ;
        GEO::CmdLine::declare_arg( "single_precision", false,
            "use single precision vertices (FP32)" ) ;
        configure_ringmesh() ;

        std::vector< std::string > ringmesh_extensions ;
        GeoModelIOHandlerFactory::list_creators( ringmesh_extensions ) ;
        ringmesh_file_extensions_ = GEO::String::join_strings( ringmesh_extensions,
            ';' ) ;

        std::vector< std::string > geogram_extensions ;
        GEO::MeshIOHandlerFactory::list_creators( geogram_extensions ) ;
        geogram_file_extensions_ = GEO::String::join_strings( geogram_extensions,
            ';' ) ;

        GEO::Logger::div( "RINGMeshView" ) ;
        GEO::Logger::out( "" ) << "Welcome to RINGMeshView !" << std::endl ;
        GEO::Logger::out( "" ) << "People working on the project in RING"
            << std::endl ;
        GEO::Logger::out( "" ) << "Arnaud Botella <arnaud.botella@univ-lorraine.fr> "
            << std::endl ;
        GEO::Logger::out( "" )
            << "Benjamin Chauvin <benjamin.chauvin@univ-lorraine.fr> " << std::endl ;
        GEO::Logger::out( "" )
            << "Antoine Mazuyer <antoine.mazuyer@univ-lorraine.fr> " << std::endl ;
    }

    RINGMeshApplication::~RINGMeshApplication()
    {
        for( index_t i = 0; i < models_.size(); i++ ) {
            delete models_[i] ;
        }
        for( index_t i = 0; i < meshes_.size(); i++ ) {
            delete meshes_[i] ;
        }
    }

    RINGMeshApplication* RINGMeshApplication::instance()
    {
        RINGMeshApplication* result =
            dynamic_cast< RINGMeshApplication* >( GEO::Application::instance() ) ;
        ringmesh_assert( result != nil ) ;
        return result ;
    }

    void RINGMeshApplication::browse_geogram( const std::string& path )
    {
        std::vector< std::string > files ;
        GEO::FileSystem::get_directory_entries( path, files ) ;
        for( GEO::index_t i = 0; i < files.size(); ++i ) {
            if( GEO::FileSystem::is_directory( files[i] ) ) {
                if( ImGui::BeginMenu( path_to_label( path_, files[i] ).c_str() ) ) {
                    browse_geogram( files[i] ) ;
                    ImGui::EndMenu() ;
                }
            } else {
                if( can_load_geogram( files[i] ) ) {
                    if( ImGui::MenuItem(
                        path_to_label( path_, files[i] ).c_str() ) ) {
                        load_geogram( files[i] ) ;
                    }
                }
            }
        }
    }

    bool RINGMeshApplication::can_load_geogram(const std::string& filename) {
        std::string extensions_str = supported_geogram_read_file_extensions() ;
        if( extensions_str == "" ) {
            return false ;
        }
        if( extensions_str == "*" ) {
            return true ;
        }
        std::string extension = GEO::FileSystem::extension( filename ) ;
        std::vector< std::string > extensions ;
        GEO::String::split_string( extensions_str, ';', extensions ) ;
        for( index_t i = 0; i < extensions.size(); ++i ) {
            if( extensions[i] == extension ) {
                return true ;
            }
        }
        return false ;
    }

    bool RINGMeshApplication::load_geogram( const std::string& filename )
    {
        if( !filename.empty() ) {
            meshes_.push_back( new MeshViewer( *this, filename ) ) ;
            current_viewer_ = meshes_.size() - 1 ;
            current_viewer_type_ = MESH ;
        }

        update_region_of_interest() ;
        return true ;
    }

    void RINGMeshApplication::draw_load_menu()
    {
        if( ImGui::BeginMenu( "Load..." ) ) {
            ImGui::Selectable( ".." ) ;
            if( ImGui::IsItemClicked() ){
                path_ += "/.." ;
            }
            browse( path_ ) ;
            ImGui::EndMenu() ;
        }
    }

    void RINGMeshApplication::draw_application_menus()
    {
        if( ImGui::BeginMenu( "Debug" ) ) {
            if( ImGui::BeginMenu( "Load..." ) ) {
                ImGui::Selectable( ".." ) ;
                if( ImGui::IsItemClicked() ){
                    path_ += "/.." ;
                }
                browse_geogram(path_);
                ImGui::EndMenu();
            }
            ImGui::EndMenu();
        }
    }

    void RINGMeshApplication::init_graphics()
    {
        GEO::Application::init_graphics() ;
        glup_viewer_add_key_func( 'c', show_corners, "corners" ) ;
        glup_viewer_add_key_func( 'e', show_lines, "lines" ) ;
        glup_viewer_add_key_func( 's', show_surface, "surface" ) ;
        glup_viewer_add_key_func( 'v', show_volume, "toggle volume" ) ;
        glup_viewer_add_key_func( 'V', show_voi, "toggle VOI" ) ;
        glup_viewer_add_key_func( 'm', mesh_visible, "mesh" ) ;
        glup_viewer_add_key_func( 'M', show_colormap, "colormap" ) ;
        glup_viewer_add_key_func( 'x', increment_shrink, "shrink cells" ) ;
        glup_viewer_add_key_func( 'X', decrement_shrink, "unshrink cells" ) ;
        glup_viewer_add_key_func( 'C', colored_cells, "toggle colored cells" ) ;
        glup_viewer_add_key_func( 'r', show_colored_regions,
            "toggle colored regions" ) ;
        glup_viewer_add_key_func( 'R', show_colored_layers,
            "toggle colored layers" ) ;

        init_colormaps() ;
        for( index_t i = 0; i < models_.size(); i++ ) {
            models_[i]->GM_gfx_.attribute.set_colormap( colormaps_[0].texture ) ;
        }
        glup_viewer_disable( GLUP_VIEWER_BACKGROUND ) ;
    }

    void RINGMeshApplication::show_corners()
    {
        if( instance()->current_viewer_type_ != GEOMODEL ) return ;
        GeoModelViewer& viewver = *instance()->models_[instance()->current_viewer_] ;
        viewver.show_corners_ = !viewver.show_corners_ ;
    }
    void RINGMeshApplication::show_lines()
    {
        if( instance()->current_viewer_type_ != GEOMODEL ) return ;
        GeoModelViewer& viewver = *instance()->models_[instance()->current_viewer_] ;
        viewver.show_lines_ = !viewver.show_lines_ ;
    }
    void RINGMeshApplication::show_surface()
    {
        if( instance()->current_viewer_type_ != GEOMODEL ) return ;
        GeoModelViewer& viewver = *instance()->models_[instance()->current_viewer_] ;
        viewver.show_surface_ = !viewver.show_surface_ ;
    }
    void RINGMeshApplication::show_volume()
    {
        if( instance()->current_viewer_type_ != GEOMODEL ) return ;
        GeoModelViewer& viewver = *instance()->models_[instance()->current_viewer_] ;
        viewver.show_volume_ = !viewver.show_volume_ ;
    }
    void RINGMeshApplication::show_voi()
    {
        if( instance()->current_viewer_type_ != GEOMODEL ) return ;
        GeoModelViewer& viewver = *instance()->models_[instance()->current_viewer_] ;
        viewver.show_voi_ = !viewver.show_voi_ ;
    }
    void RINGMeshApplication::mesh_visible()
    {
        if( instance()->current_viewer_type_ != GEOMODEL ) return ;
        GeoModelViewer& viewver = *instance()->models_[instance()->current_viewer_] ;
        viewver.mesh_visible_ = !viewver.mesh_visible_ ;
    }
    void RINGMeshApplication::show_colormap()
    {
        if( instance()->current_viewer_type_ != GEOMODEL ) return ;
        GeoModelViewer& viewver = *instance()->models_[instance()->current_viewer_] ;
        viewver.show_colormap_ = !viewver.show_colormap_ ;
    }
    void RINGMeshApplication::colored_cells()
    {
        if( instance()->current_viewer_type_ != GEOMODEL ) return ;
        GeoModelViewer& viewver = *instance()->models_[instance()->current_viewer_] ;
        viewver.colored_cells_.new_status = !viewver.colored_cells_.new_status ;
    }
    void RINGMeshApplication::show_colored_regions()
    {
        if( instance()->current_viewer_type_ != GEOMODEL ) return ;
        GeoModelViewer& viewver = *instance()->models_[instance()->current_viewer_] ;
        viewver.show_colored_regions_.new_status =
            !viewver.show_colored_regions_.new_status ;
    }
    void RINGMeshApplication::show_colored_layers()
    {
        if( instance()->current_viewer_type_ != GEOMODEL ) return ;
        GeoModelViewer& viewver = *instance()->models_[instance()->current_viewer_] ;
        viewver.show_colored_layers_.new_status =
            !viewver.show_colored_layers_.new_status ;
    }
    void RINGMeshApplication::increment_shrink()
    {
        if( instance()->current_viewer_type_ != GEOMODEL ) return ;
        GeoModelViewer& viewver = *instance()->models_[instance()->current_viewer_] ;
        viewver.shrink_ = std::min( viewver.shrink_ + 0.1f, 1.f ) ;
    }
    void RINGMeshApplication::decrement_shrink()
    {
        if( instance()->current_viewer_type_ != GEOMODEL ) return ;
        GeoModelViewer& viewver = *instance()->models_[instance()->current_viewer_] ;
        viewver.shrink_ = std::max( viewver.shrink_ - 0.1f, 0.f ) ;
    }

    bool RINGMeshApplication::load( const std::string& filename )
    {
        if( !filename.empty() ) {
            models_.push_back( new GeoModelViewer( *this, filename ) ) ;
            current_viewer_ = models_.size() - 1 ;
            current_viewer_type_ = GEOMODEL ;
        }

        update_region_of_interest() ;
        return true ;
    }

    void RINGMeshApplication::update_region_of_interest()
    {
        Box3d bbox ;

        for( index_t i = 0; i < models_.size(); i++ ) {
            if( models_[i]->is_visible_ ) {
                bbox.add_box( models_[i]->bbox_ ) ;
            }
        }
        for( index_t i = 0; i < meshes_.size(); i++ ) {
            if( meshes_[i]->is_visible_ ) {
                bbox.add_box( meshes_[i]->bbox_ ) ;
            }
        }

        glup_viewer_set_region_of_interest( float( bbox.min()[0] ),
            float( bbox.min()[1] ), float( bbox.min()[2] ), float( bbox.max()[0] ),
            float( bbox.max()[1] ), float( bbox.max()[2] ) ) ;
    }

    void RINGMeshApplication::draw_scene()
    {
        if( current_viewer_ == NO_ID ) return ;

        for( index_t i = 0; i < meshes_.size(); i++ ) {
            if( meshes_[i]->is_visible_ ) meshes_[i]->draw_scene() ;
        }
        for( index_t i = 0; i < models_.size(); i++ ) {
            if( models_[i]->is_visible_ ) models_[i]->draw_scene() ;
        }

        if( current_viewer_type_ == GEOMODEL ) {
            GeoModelViewer& viewer = *models_[current_viewer_] ;
            if( viewer.show_colormap_ ) {
                viewer.draw_colormap() ;
            }
        }
    }

    std::string RINGMeshApplication::supported_read_file_extensions()
    {
        return ringmesh_file_extensions_ ;
    }
    std::string RINGMeshApplication::supported_geogram_read_file_extensions()
    {
        return geogram_file_extensions_ ;
    }

    void RINGMeshApplication::draw_viewer_properties()
    {
        GEO::Application::draw_viewer_properties() ;

        if( !models_.empty() ) {
            ImGui::Separator() ;
            ImGui::Text( "GeoModel" ) ;
            for( index_t i = 0; i < models_.size(); i++ ) {
                GeoModelViewer& viewer = *models_[i] ;
                if( ImGui::Checkbox( viewer.GM_.name().c_str(),
                    &viewer.is_visible_ ) ) {
                    current_viewer_ = i ;
                    current_viewer_type_ = GEOMODEL ;
                    update_region_of_interest() ;
                }
            }
        }

        if( !meshes_.empty() ) {
            ImGui::Separator() ;
            ImGui::Text( "Mesh" ) ;
            for( index_t i = 0; i < meshes_.size(); i++ ) {
                MeshViewer& viewer = *meshes_[i] ;
                if( ImGui::Checkbox( viewer.name_.c_str(),
                    &viewer.is_visible_ ) ) {
                    current_viewer_ = i ;
                    current_viewer_type_ = MESH ;
                    update_region_of_interest() ;
                }
            }
        }
    }

    void RINGMeshApplication::draw_object_properties()
    {
        if( current_viewer_ == NO_ID ) return ;
        switch( current_viewer_type_ ) {
            case GEOMODEL :
                ringmesh_assert( current_viewer_ < models_.size() ) ;
                models_[current_viewer_]->draw_object_properties() ;
                return ;
            case MESH:
                ringmesh_assert( current_viewer_ < meshes_.size() ) ;
                meshes_[current_viewer_]->draw_object_properties() ;
                return;
        }
    }
}

#endif
