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

#include <geogram/mesh/mesh_geometry.h>
#include <geogram/basic/logger.h>

#include <geogram_gfx/glup_viewer/glup_viewer.h>

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

    GeoModelGfx::GeoModelGfx()
        :
            model_( nil ),
            corners_(),
            lines_(),
            surfaces_(),
            regions_(),
            attribute_location_( nb_locations ),
            attribute_coordinate_( 0 ),
            attribute_min_( 0.0 ),
            attribute_max_( 0.0 )
    {
    }

    GeoModelGfx::~GeoModelGfx()
    {
        clear_storage() ;
    }

    void GeoModelGfx::clear_storage()
    {
        for( index_t c = 0; c < corners_.size(); c++ ) {
            delete corners_[c] ;
        }
        for( index_t l = 0; l < lines_.size(); l++ ) {
            delete lines_[l] ;
        }
        for( index_t s = 0; s < surfaces_.size(); s++ ) {
            delete surfaces_[s] ;
        }
        for( index_t r = 0; r < regions_.size(); r++ ) {
            delete regions_[r] ;
        }
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
        clear_storage() ;
        if( corners_.empty() && lines_.empty() && surfaces_.empty() ) {
            corners_.resize( model_->nb_corners(), nil ) ;
            lines_.resize( model_->nb_lines(), nil ) ;
            surfaces_.resize( model_->nb_surfaces(), nil ) ;
            regions_.resize( model_->nb_regions(), nil ) ;

            for( index_t c = 0; c < corners_.size(); c++ ) {
                corners_[c] = new CornerGfx( *this, model_->corner( c ) ) ;
            }
            for( index_t l = 0; l < lines_.size(); l++ ) {
                lines_[l] = new LineGfx( *this, model_->line( l ) ) ;
            }
            for( index_t s = 0; s < surfaces_.size(); s++ ) {
                surfaces_[s] = new SurfaceGfx( *this, model_->surface( s ) ) ;
            }
            for( index_t r = 0; r < model_->nb_regions(); r++ ) {
                regions_[r] = new RegionGfx( *this, model_->region( r ) ) ;
            }
        }
    }

    void GeoModelGfx::need_to_update()
    {
        for( index_t c = 0; c < corners_.size(); c++ ) {
            corners_[c]->need_to_update() ;
        }
        for( index_t l = 0; l < lines_.size(); l++ ) {
            lines_[l]->need_to_update() ;
        }
        for( index_t s = 0; s < surfaces_.size(); s++ ) {
            surfaces_[s]->need_to_update() ;
        }
        for( index_t r = 0; r < model_->nb_regions(); r++ ) {
            regions_[r]->need_to_update() ;
        }
    }

    void GeoModelGfx::compute_cell_vertex_attribute_range(
        index_t coordinate,
        const std::string& name )
    {
        attribute_min_ = max_float64() ;
        attribute_max_ = min_float64() ;
        std::string attribute_name = get_attribute_name_with_coordinate( name,
            coordinate ) ;
        for( index_t r = 0; r < regions_.size(); r++ ) {
            GEO::ReadOnlyScalarAttributeAdapter attribute(
                model_->region( r ).vertex_attribute_manager(), attribute_name ) ;
            compute_attribute_range( attribute, attribute_min_, attribute_max_ ) ;
        }
    }

    void GeoModelGfx::compute_cell_attribute_range(
        index_t coordinate,
        const std::string& name )
    {
        attribute_min_ = max_float64() ;
        attribute_max_ = min_float64() ;
        std::string attribute_name = get_attribute_name_with_coordinate( name,
            coordinate ) ;
        for( index_t r = 0; r < regions_.size(); r++ ) {
            GEO::ReadOnlyScalarAttributeAdapter attribute(
                model_->region( r ).cell_attribute_manager(), attribute_name ) ;
            compute_attribute_range( attribute, attribute_min_, attribute_max_ ) ;
        }
    }

    void GeoModelGfx::set_cell_vertex_attribute(
        const std::string& name,
        index_t coordinate,
        GLuint colormap_texture )
    {
        std::string attribute_name = get_attribute_name_with_coordinate( name,
            coordinate ) ;
        for( index_t r = 0; r < regions_.size(); r++ ) {
            regions_[r]->set_scalar_attribute( GEO::MESH_VERTICES, attribute_name,
                attribute_min_, attribute_max_, colormap_texture ) ;
        }

    }

    void GeoModelGfx::set_cell_attribute(
        const std::string& name,
        index_t coordinate,
        GLuint colormap_texture )
    {
        std::string attribute_name = get_attribute_name_with_coordinate( name,
            coordinate ) ;
        for( index_t r = 0; r < regions_.size(); r++ ) {
            regions_[r]->set_scalar_attribute( GEO::MESH_CELLS, attribute_name,
                attribute_min_, attribute_max_, colormap_texture ) ;
        }

    }

    /*!
     * Draws the corners
     */
    void GeoModelGfx::draw_corners()
    {
        for( index_t c = 0; c < corners_.size(); c++ ) {
            if( corners_[c]->get_vertices_visible() ) corners_[c]->draw_vertices() ;
        }
    }
    /*!
     * Sets the corner color to all the corners
     * @param[in] r the red component of the color in [0.0, 1.0]
     * @param[in] g the green component of the color in [0.0, 1.0]
     * @param[in] b the blue component of the color in [0.0, 1.0]
     */
    void GeoModelGfx::set_corners_color( float r, float g, float b )
    {
        for( index_t c = 0; c < corners_.size(); c++ ) {
            set_corner_color( c, r, g, b ) ;
        }
    }
    /*!
     * Sets the corner color
     * @param[in] c the corner index
     * @param[in] r the red component of the color in [0.0, 1.0]
     * @param[in] g the green component of the color in [0.0, 1.0]
     * @param[in] b the blue component of the color in [0.0, 1.0]
     */
    void GeoModelGfx::set_corner_color( index_t c, float r, float g, float b )
    {
        corners_[c]->set_points_color( r, g, b ) ;
    }
    /*!
     * Sets the corner visibility to all the corners
     * @param[in] b the visibility
     */
    void GeoModelGfx::set_corners_visibility( bool b )
    {
        for( index_t c = 0; c < corners_.size(); c++ ) {
            set_corner_visibility( c, b ) ;
        }
    }
    /*!
     * Sets the corner visibility
     * @param[in] c the corner index
     * @param[in] b the visibility
     */
    void GeoModelGfx::set_corner_visibility( index_t c, bool b )
    {
        corners_[c]->set_vertices_visible( b ) ;
    }
    /*!
     * Sets the corner size to all the corners
     * @param[in] s the size
     */
    void GeoModelGfx::set_corners_size( index_t s )
    {
        for( index_t c = 0; c < corners_.size(); c++ ) {
            set_corner_size( c, s ) ;
        }
    }
    /*!
     * Sets the corner size
     * @param[in] c the corner index
     * @param[in] s the size
     */
    void GeoModelGfx::set_corner_size( index_t c, index_t s )
    {
        corners_[c]->set_points_size( float( s ) ) ;
    }

    /*!
     * Draws the lines
     */
    void GeoModelGfx::draw_lines()
    {
        for( index_t l = 0; l < lines_.size(); l++ ) {
            if( lines_[l]->get_vertices_visible() ) lines_[l]->draw_vertices() ;
            if( lines_[l]->get_edges_visible() ) lines_[l]->draw_edges() ;
        }
    }
    /*!
     * Sets the line color to all the lines
     * @param[in] r the red component of the color in [0.0, 1.0]
     * @param[in] g the green component of the color in [0.0, 1.0]
     * @param[in] b the blue component of the color in [0.0, 1.0]
     */
    void GeoModelGfx::set_edge_lines_color( float r, float g, float b )
    {
        for( index_t k = 0; k < lines_.size(); k++ ) {
            set_edge_line_color( k, r, g, b ) ;
        }
    }
    /*!
     * Sets the line color
     * @param[in] l the line index
     * @param[in] r the red component of the color in [0.0, 1.0]
     * @param[in] g the green component of the color in [0.0, 1.0]
     * @param[in] b the blue component of the color in [0.0, 1.0]
     */
    void GeoModelGfx::set_edge_line_color( index_t l, float r, float g, float b )
    {
        lines_[l]->set_mesh_color( r, g, b ) ;
    }
    /*!
     * Sets the line visibility to all the lines
     * @param[in] b the visibility
     */
    void GeoModelGfx::set_edge_lines_visibility( bool b )
    {
        for( index_t l = 0; l < lines_.size(); l++ ) {
            set_edge_line_visibility( l, b ) ;
        }
    }
    /*!
     * Sets the line visibility
     * @param[in] l the line index
     * @param[in] b the visibility
     */
    void GeoModelGfx::set_edge_line_visibility( index_t l, bool b )
    {
        lines_[l]->set_edges_visible( b ) ;
    }
    /*!
     * Sets the edge line size to all the lines
     * @param[in] s the size
     */
    void GeoModelGfx::set_edge_lines_size( index_t s )
    {
        for( index_t l = 0; l < lines_.size(); l++ ) {
            set_edge_line_size( l, s ) ;
        }
    }
    /*!
     * Sets the edge line size
     * @param[in] l the line index
     * @param[in] s the size
     */
    void GeoModelGfx::set_edge_line_size( index_t l, index_t s )
    {
        lines_[l]->set_mesh_width( s ) ;
    }
    /*!
     * Sets the vertex line color to all the lines
     * @param[in] r the red component of the color in [0.0, 1.0]
     * @param[in] g the green component of the color in [0.0, 1.0]
     * @param[in] b the blue component of the color in [0.0, 1.0]
     */
    void GeoModelGfx::set_vertex_lines_color( float r, float g, float b )
    {
        for( index_t k = 0; k < lines_.size(); k++ ) {
            set_vertex_line_color( k, r, g, b ) ;
        }
    }
    /*!
     * Sets the vertex line color
     * @param[in] l the line index
     * @param[in] r the red component of the color in [0.0, 1.0]
     * @param[in] g the green component of the color in [0.0, 1.0]
     * @param[in] b the blue component of the color in [0.0, 1.0]
     */
    void GeoModelGfx::set_vertex_line_color( index_t l, float r, float g, float b )
    {
        lines_[l]->set_points_color( r, g, b ) ;
    }
    /*!
     * Sets the vertex line visibility to all the lines
     * @param[in] b the visibility
     */
    void GeoModelGfx::set_vertex_lines_visibility( bool b )
    {
        for( index_t l = 0; l < lines_.size(); l++ ) {
            set_vertex_line_visibility( l, b ) ;
        }
    }
    /*!
     * Sets the vertex line visibility
     * @param[in] l the line index
     * @param[in] b the visibility
     */
    void GeoModelGfx::set_vertex_line_visibility( index_t l, bool b )
    {
        lines_[l]->set_vertices_visible( b ) ;
    }
    /*!
     * Sets the vertex line size to all the lines
     * @param[in] s the size
     */
    void GeoModelGfx::set_vertex_lines_size( index_t s )
    {
        for( index_t l = 0; l < lines_.size(); l++ ) {
            set_vertex_line_size( l, s ) ;
        }
    }
    /*!
     * Sets the vertex line size
     * @param[in] l the line index
     * @param[in] s the size
     */
    void GeoModelGfx::set_vertex_line_size( index_t l, index_t s )
    {
        lines_[l]->set_points_size( float( s ) ) ;
    }

    /*!
     * Draws the surfaces
     */
    void GeoModelGfx::draw_surfaces()
    {
        for( index_t s = 0; s < surfaces_.size(); s++ ) {
            if( surfaces_[s]->get_vertices_visible() )
                surfaces_[s]->draw_vertices() ;
            if( surfaces_[s]->get_surface_visible() ) surfaces_[s]->draw_surface() ;
        }
    }
    void GeoModelGfx::set_surfaces_lighting( bool value )
    {
        for( index_t s = 0; s < surfaces_.size(); s++ ) {
            surfaces_[s]->set_lighting( value ) ;
        }
    }
    /*!
     * Sets the surface color to all the surfaces
     * @param[in] r the red component of the color in [0.0, 1.0]
     * @param[in] g the green component of the color in [0.0, 1.0]
     * @param[in] b the blue component of the color in [0.0, 1.0]
     */
    void GeoModelGfx::set_surfaces_color( float r, float g, float b )
    {
        for( index_t s = 0; s < surfaces_.size(); s++ ) {
            set_surface_color( s, r, g, b ) ;
        }
    }
    /*!
     * Sets the surface color
     * @param[in] s the surface index
     * @param[in] r the red component of the color in [0.0, 1.0]
     * @param[in] g the green component of the color in [0.0, 1.0]
     * @param[in] b the blue component of the color in [0.0, 1.0]
     */
    void GeoModelGfx::set_surface_color( index_t s, float r, float g, float b )
    {
        surfaces_[s]->set_surface_color( r, g, b ) ;
    }
    /*!
     * Sets the backface surface color to all the surfaces
     * @param[in] r the red component of the color in [0.0, 1.0]
     * @param[in] g the green component of the color in [0.0, 1.0]
     * @param[in] b the blue component of the color in [0.0, 1.0]
     */
    void GeoModelGfx::set_backface_surfaces_color( float r, float g, float b )
    {
        for( index_t s = 0; s < surfaces_.size(); s++ ) {
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
    void GeoModelGfx::set_backface_surface_color(
        index_t s,
        float r,
        float g,
        float b )
    {
        surfaces_[s]->set_backface_surface_color( r, g, b ) ;
    }
    /*!
     * Sets the surface visibility to all the surfaces
     * @param[in] b the visibility
     */
    void GeoModelGfx::set_surfaces_visibility( bool b )
    {
        for( index_t s = 0; s < surfaces_.size(); s++ ) {
            set_surface_visibility( s, b ) ;
        }
    }
    /*!
     * Sets the surface visibility
     * @param[in] s the surface index
     * @param[in] b the visibility
     */
    void GeoModelGfx::set_surface_visibility( index_t s, bool b )
    {
        surfaces_[s]->set_surface_visible( b ) ;
    }
    /*!
     * Sets the mesh surface color to all the surfaces
     * @param[in] r the red component of the color in [0.0, 1.0]
     * @param[in] g the green component of the color in [0.0, 1.0]
     * @param[in] b the blue component of the color in [0.0, 1.0]
     */
    void GeoModelGfx::set_mesh_surfaces_color( float r, float g, float b )
    {
        for( index_t s = 0; s < surfaces_.size(); s++ ) {
            set_mesh_surface_color( s, r, g, b ) ;
        }
    }
    /*!
     * Sets the mesh surface color
     * @param[in] s the surface index
     * @param[in] r the red component of the color in [0.0, 1.0]
     * @param[in] g the green component of the color in [0.0, 1.0]
     * @param[in] b the blue component of the color in [0.0, 1.0]
     */
    void GeoModelGfx::set_mesh_surface_color( index_t s, float r, float g, float b )
    {
        surfaces_[s]->set_mesh_color( r, g, b ) ;
    }
    /*!
     * Sets the mesh surface visibility to all the surfaces
     * @param[in] b the visibility
     */
    void GeoModelGfx::set_mesh_surfaces_visibility( bool b )
    {
        for( index_t s = 0; s < surfaces_.size(); s++ ) {
            set_mesh_surface_visibility( s, b ) ;
        }
    }
    /*!
     * Sets the mesh surface visibility
     * @param[in] s the surface index
     * @param[in] b the visibility
     */
    void GeoModelGfx::set_mesh_surface_visibility( index_t s, bool b )
    {
        surfaces_[s]->set_show_mesh( b ) ;
    }
    /*!
     * Sets the mesh surface size to all the surfaces
     * @param[in] size the size
     */
    void GeoModelGfx::set_mesh_surfaces_size( index_t size )
    {
        for( index_t s = 0; s < surfaces_.size(); s++ ) {
            set_mesh_surface_size( s, size ) ;
        }
    }
    /*!
     * Sets the mesh surface size
     * @param[in] s the surface index
     * @param[in] size the size
     */
    void GeoModelGfx::set_mesh_surface_size( index_t s, index_t size )
    {
        surfaces_[s]->set_mesh_width( size ) ;
    }
    /*!
     * Sets the vertex surface color to all the surfaces
     * @param[in] r the red component of the color in [0.0, 1.0]
     * @param[in] g the green component of the color in [0.0, 1.0]
     * @param[in] b the blue component of the color in [0.0, 1.0]
     */
    void GeoModelGfx::set_vertex_surfaces_color( float r, float g, float b )
    {
        for( index_t s = 0; s < surfaces_.size(); s++ ) {
            set_vertex_surface_color( s, r, g, b ) ;
        }
    }
    /*!
     * Sets the vertex surface color
     * @param[in] s the surface index
     * @param[in] r the red component of the color in [0.0, 1.0]
     * @param[in] g the green component of the color in [0.0, 1.0]
     * @param[in] b the blue component of the color in [0.0, 1.0]
     */
    void GeoModelGfx::set_vertex_surface_color(
        index_t s,
        float r,
        float g,
        float b )
    {
        surfaces_[s]->set_points_color( r, g, b ) ;
    }
    /*!
     * Sets the vertex surface visibility to all the surfaces
     * @param[in] b the visibility
     */
    void GeoModelGfx::set_vertex_surfaces_visibility( bool b )
    {
        for( index_t s = 0; s < surfaces_.size(); s++ ) {
            set_vertex_surface_visibility( s, b ) ;
        }
    }
    /*!
     * Sets the vertex surface visibility
     * @param[in] s the surface index
     * @param[in] b the visibility
     */
    void GeoModelGfx::set_vertex_surface_visibility( index_t s, bool b )
    {
        surfaces_[s]->set_vertices_visible( b ) ;
    }
    /*!
     * Sets the vertex surface size to all the surfaces
     * @param[in] size the size
     */
    void GeoModelGfx::set_vertex_surfaces_size( index_t size )
    {
        for( index_t s = 0; s < surfaces_.size(); s++ ) {
            set_vertex_surface_size( s, size ) ;
        }
    }
    /*!
     * Sets the vertex surface size
     * @param[in] s the surface index
     * @param[in] size the size
     */
    void GeoModelGfx::set_vertex_surface_size( index_t s, index_t size )
    {
        surfaces_[s]->set_points_size( float( size ) ) ;
    }

    /*!
     * Draws the MacroMesh
     */
    void GeoModelGfx::draw_regions()
    {
        for( index_t m = 0; m < regions_.size(); m++ ) {
            if( regions_[m]->get_vertices_visible() ) regions_[m]->draw_vertices() ;
            if( regions_[m]->get_edges_visible() ) regions_[m]->draw_edges() ;
            if( regions_[m]->get_surface_visible() ) regions_[m]->draw_surface() ;
            if( regions_[m]->get_region_visible() ) regions_[m]->draw_volume() ;
        }
    }

    /*!
     * Sets the vertex region color to all the regions
     * @param[in] r the red component of the color in [0.0, 1.0]
     * @param[in] g the green component of the color in [0.0, 1.0]
     * @param[in] b the blue component of the color in [0.0, 1.0]
     */
    void GeoModelGfx::set_vertex_regions_color( float r, float g, float b )
    {
        for( index_t m = 0; m < regions_.size(); m++ ) {
            set_vertex_region_color( m, r, g, b ) ;
        }
    }

    /*!
     * Sets the vertex region color
     * @param[in] m the region index
     * @param[in] r the red component of the color in [0.0, 1.0]
     * @param[in] g the green component of the color in [0.0, 1.0]
     * @param[in] b the blue component of the color in [0.0, 1.0]
     */
    void GeoModelGfx::set_vertex_region_color( index_t m, float r, float g, float b )
    {
        ringmesh_assert( m < regions_.size() ) ;
        regions_[m]->set_points_color( r, g, b ) ;
    }

    /*!
     * Sets the vertex region visibility to all the regions
     * @param[in] b the visibility
     */
    void GeoModelGfx::set_vertex_regions_visibility( bool b )
    {
        for( index_t m = 0; m < regions_.size(); m++ ) {
            set_vertex_region_visibility( m, b ) ;
        }
    }

    /*!
     * Sets the vertex region visibility
     * @param[in] m the region index
     * @param[in] b the visibility
     */
    void GeoModelGfx::set_vertex_region_visibility( index_t m, bool b )
    {
        ringmesh_assert( m < regions_.size() ) ;
        regions_[m]->set_vertices_visible( b ) ;
    }

    /*!
     * Sets the vertex region size to all the regions
     * @param[in] s the size
     */
    void GeoModelGfx::set_vertex_regions_size( index_t s )
    {
        for( index_t m = 0; m < regions_.size(); m++ ) {
            set_vertex_region_size( m, s ) ;
        }
    }

    /*!
     * Sets the vertex region size
     * @param[in] m the region index
     * @param[in] s the size
     */
    void GeoModelGfx::set_vertex_region_size( index_t m, index_t s )
    {
        ringmesh_assert( m < regions_.size() ) ;
        regions_[m]->set_points_size( float( s ) ) ;
    }

    /*!
     * Sets the edge color to all the meshes
     * @param[in] r the red component of the color in [0.0, 1.0]
     * @param[in] g the green component of the color in [0.0, 1.0]
     * @param[in] b the blue component of the color in [0.0, 1.0]
     */
    void GeoModelGfx::set_edge_regions_color( float r, float g, float b )
    {
        for( index_t k = 0; k < regions_.size(); k++ ) {
            set_edge_region_color( k, r, g, b ) ;
        }
    }
    /*!
     * Sets the edge color
     * @param[in] m the mesh index
     * @param[in] r the red component of the color in [0.0, 1.0]
     * @param[in] g the green component of the color in [0.0, 1.0]
     * @param[in] b the blue component of the color in [0.0, 1.0]
     */
    void GeoModelGfx::set_edge_region_color( index_t m, float r, float g, float b )
    {
        ringmesh_assert( m < regions_.size() ) ;
        regions_[m]->set_mesh_color( r, g, b ) ; //TODO function not good?
    }
    /*!
     * Sets the edge visibility to all the meshes
     * @param[in] b the visibility
     */
    void GeoModelGfx::set_edge_regions_visibility( bool b )
    {
        for( index_t m = 0; m < regions_.size(); m++ ) {
            set_edge_region_visibility( m, b ) ;
        }
    }
    /*!
     * Sets the edge visibility
     * @param[in] m the mesh index
     * @param[in] b the visibility
     */
    void GeoModelGfx::set_edge_region_visibility( index_t m, bool b )
    {
        ringmesh_assert( m < regions_.size() ) ;
        regions_[m]->set_edges_visible( b ) ;
    }
    /*!
     * Sets the edge line size to all the lines
     * @param[in] s the size
     */
    void GeoModelGfx::set_edge_regions_size( index_t s )
    {
        for( index_t m = 0; m < regions_.size(); m++ ) {
            set_edge_region_size( m, s ) ;
        }
    }
    /*!
     * Sets the edge line size
     * @param[in] l the line index
     * @param[in] s the size
     */
    void GeoModelGfx::set_edge_region_size( index_t l, index_t s )
    {
        ringmesh_assert( l < regions_.size() ) ;
        regions_[l]->set_mesh_width( s ) ;
    }

    void GeoModelGfx::set_cell_regions_lighting( bool value )
    {
        for( index_t r = 0; r < regions_.size(); r++ ) {
            regions_[r]->set_lighting( value ) ;
        }
    }

    /*!
     * Sets the mesh region color to all the regions
     * @param[in] r the red component of the color in [0.0, 1.0]
     * @param[in] g the green component of the color in [0.0, 1.0]
     * @param[in] b the blue component of the color in [0.0, 1.0]
     */
    void GeoModelGfx::set_cell_mesh_regions_color( float r, float g, float b )
    {
        for( index_t m = 0; m < regions_.size(); m++ ) {
            set_cell_mesh_region_color( m, r, g, b ) ;
        }
    }

    /*!
     * Sets the mesh region color
     * @param[in] m the region index
     * @param[in] r the red component of the color in [0.0, 1.0]
     * @param[in] g the green component of the color in [0.0, 1.0]
     * @param[in] b the blue component of the color in [0.0, 1.0]
     */
    void GeoModelGfx::set_cell_mesh_region_color(
        index_t m,
        float r,
        float g,
        float b )
    {
        ringmesh_assert( m < regions_.size() ) ;
        regions_[m]->set_mesh_color( r, g, b ) ;
    }

    /*!
     * Toggles the cell region color per cell type to all the regions
     */
    void GeoModelGfx::set_cell_regions_color_type()
    {
        for( index_t m = 0; m < regions_.size(); m++ ) {
            set_cell_region_color_type( m ) ;
        }
    }

    /*!
     * Toggles the cell region color per cell type
     * @param[in] m the region index
     */
    void GeoModelGfx::set_cell_region_color_type( index_t m )
    {
        ringmesh_assert( m < regions_.size() ) ;
        regions_[m]->set_cells_colors_by_type() ;
    }

    /*!
     * Sets the mesh region visibility to all the regions
     * @param[in] b the visibility
     */
    void GeoModelGfx::set_cell_mesh_regions_visibility( bool b )
    {
        for( index_t m = 0; m < regions_.size(); m++ ) {
            set_cell_mesh_region_visibility( m, b ) ;
        }
    }

    /*!
     * Sets the mesh region visibility
     * @param[in] m the region index
     * @param[in] b the visibility
     */
    void GeoModelGfx::set_cell_mesh_region_visibility( index_t m, bool b )
    {
        ringmesh_assert( m < regions_.size() ) ;
        regions_[m]->set_show_mesh( b ) ;
    }

    /*!
     * Sets the mesh region size to all the regions
     * @param[in] s the size
     */
    void GeoModelGfx::set_cell_mesh_regions_size( index_t s )
    {
        for( index_t m = 0; m < regions_.size(); m++ ) {
            set_cell_mesh_region_size( m, s ) ;
        }
    }

    /*!
     * Sets the mesh region size
     * @param[in] m the region index
     * @param[in] s the size
     */
    void GeoModelGfx::set_cell_mesh_region_size( index_t m, index_t s )
    {
        ringmesh_assert( m < regions_.size() ) ;
        regions_[m]->set_mesh_width( s ) ;
    }

    /*!
     * Sets the cell region color to all the regions
     * @param[in] r the red component of the color in [0.0, 1.0]
     * @param[in] g the green component of the color in [0.0, 1.0]
     * @param[in] b the blue component of the color in [0.0, 1.0]
     */
    void GeoModelGfx::set_cell_regions_color( float r, float g, float b )
    {
        for( index_t m = 0; m < regions_.size(); m++ ) {
            set_cell_region_color( m, r, g, b ) ;
        }
    }

    /*!
     * Sets the cell region color
     * @param[in] m the region index
     * @param[in] r the red component of the color in [0.0, 1.0]
     * @param[in] g the green component of the color in [0.0, 1.0]
     * @param[in] b the blue component of the color in [0.0, 1.0]
     */
    void GeoModelGfx::set_cell_region_color( index_t m, float r, float g, float b )
    {
        ringmesh_assert( m < regions_.size() ) ;
        regions_[m]->set_cells_color( r, g, b ) ;
    }

    /*!
     * Sets the cell region visibility to all the regions
     * @param[in] b the visibility
     */
    void GeoModelGfx::set_cell_regions_visibility( bool b )
    {
        for( index_t m = 0; m < regions_.size(); m++ ) {
            set_cell_region_visibility( m, b ) ;
        }
    }

    /*!
     * Sets the cell region visibility to all the regions
     * @param[in] m the region index
     * @param[in] b the visibility
     */
    void GeoModelGfx::set_cell_region_visibility( index_t m, bool b )
    {
        ringmesh_assert( m < regions_.size() ) ;
        regions_[m]->set_region_visible( b ) ;
    }
    void GeoModelGfx::set_cell_regions_type_visibility( GEO::MeshCellType t, bool b )
    {
        for( index_t m = 0; m < regions_.size(); m++ ) {
            set_cell_region_type_visibility( m, t, b ) ;
        }
    }
    void GeoModelGfx::set_cell_region_type_visibility(
        index_t m,
        GEO::MeshCellType t,
        bool b )
    {
        ringmesh_assert( m < regions_.size() ) ;
        regions_[m]->set_draw_cells( t, b ) ;
    }

    /*!
     * Sets the cell region shrink to all the regions
     * @param[in] s the shrink
     */
    void GeoModelGfx::set_cell_regions_shrink( double s )
    {
        for( index_t m = 0; m < regions_.size(); m++ ) {
            set_cell_region_shrink( m, s ) ;
        }
    }

    /*!
     * Sets the cell region shrink
     * @param[in] m the region index
     * @param[in] s the shrink
     */
    void GeoModelGfx::set_cell_region_shrink( index_t m, double s )
    {
        ringmesh_assert( m < regions_.size() ) ;
        regions_[m]->set_shrink( s ) ;
    }

    /*****************************************************************/

    RINGMeshApplication::RINGMeshApplication( int argc, char** argv )
        : GEO::Application( argc, argv, "<filename>" )
    {
        configure_ringmesh() ;

        GM_ = nil ;

        show_corners_ = true ;
        show_lines_ = true ;
        show_surface_ = true ;
        show_volume_ = false ;
        colored_cells_ = false ;
        show_voi_ = true ;
        show_colored_regions_ = false ;
        show_colored_layers_ = false ;

        shrink_ = 0.0 ;
        mesh_visible_ = true ;
        meshed_regions_ = false ;

        show_attributes_ = false ;
        attribute_min_ = 0.0f ;
        attribute_max_ = 0.0f ;
        reset_attribute_name() ;

        std::vector< std::string > extensions ;
        GeoModelIOHandlerFactory::list_creators( extensions ) ;
        file_extensions_ = GEO::String::join_strings( extensions, ';' ) ;

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
        if( GM_ ) {
            delete GM_ ;
        }
    }

    void RINGMeshApplication::reset_attribute_name()
    {
        GM_gfx_.set_attribute_name( "name" ) ;
    }

    RINGMeshApplication* RINGMeshApplication::instance()
    {
        RINGMeshApplication* result =
            dynamic_cast< RINGMeshApplication* >( GEO::Application::instance() ) ;
        ringmesh_assert( result != nil ) ;
        return result ;
    }

    void RINGMeshApplication::init_graphics()
    {
        GEO::Application::init_graphics() ;
        glup_viewer_add_toggle( 'c', &show_corners_, "corners" ) ;
        glup_viewer_add_toggle( 'e', &show_lines_, "lines" ) ;
        glup_viewer_add_toggle( 's', &show_surface_, "surface" ) ;
        glup_viewer_add_toggle( 'v', &show_volume_, "toggle volume" ) ;
        glup_viewer_add_toggle( 'V', &show_voi_, "toggle VOI" ) ;
        glup_viewer_add_toggle( 'm', &mesh_visible_, "mesh" ) ;
        glup_viewer_add_key_func( 'x', increment_shrink, "shrink cells" ) ;
        glup_viewer_add_key_func( 'X', decrement_shrink, "unshrink cells" ) ;
        glup_viewer_add_key_func( 'C', toggle_colored_cells,
            "toggle colored cells" ) ;
        glup_viewer_add_key_func( 'r', toggle_colored_regions,
            "toggle colored regions" ) ;
        glup_viewer_add_key_func( 'R', &toggle_colored_layers,
            "toggle colored layers" ) ;

        init_colormaps();
        GM_gfx_.set_colormap( colormaps_[0].texture ) ;
    }

    void RINGMeshApplication::increment_shrink()
    {
        instance()->shrink_ = std::min( instance()->shrink_ + 0.1f, 1.f ) ;
    }
    void RINGMeshApplication::decrement_shrink()
    {
        instance()->shrink_ = std::max( instance()->shrink_ - 0.1f, 0.f ) ;
    }
    void RINGMeshApplication::toggle_colored_cells()
    {
        instance()->colored_cells_ = !instance()->colored_cells_ ;
        if( instance()->colored_cells_ ) {
            instance()->show_colored_regions_ = false ;
            instance()->show_colored_layers_ = false ;
            instance()->GM_gfx_.set_cell_regions_color_type() ;
        }
    }
    void RINGMeshApplication::toggle_colored_regions()
    {
        instance()->show_colored_regions_ = !instance()->show_colored_regions_ ;
        if( instance()->show_colored_regions_ ) {
            instance()->colored_cells_ = false ;
            instance()->show_colored_layers_ = false ;
            for( GEO::index_t r = 0; r < instance()->GM_->nb_regions(); r++ ) {
                instance()->GM_gfx_.set_cell_region_color( r,
                    std::fmod( GEO::Numeric::random_float32(), 1 ),
                    std::fmod( GEO::Numeric::random_float32(), 1 ),
                    std::fmod( GEO::Numeric::random_float32(), 1 ) ) ;
            }
        }
    }

    void RINGMeshApplication::toggle_colored_layers()
    {
        instance()->show_colored_layers_ = !instance()->show_colored_layers_ ;
        if( instance()->show_colored_layers_ ) {
            instance()->colored_cells_ = false ;
            instance()->show_colored_regions_ = false ;
            for( GEO::index_t l = 0; l < instance()->GM_->nb_layers(); l++ ) {
                float red = std::fmod( GEO::Numeric::random_float32(), 1 ) ;
                float green = std::fmod( GEO::Numeric::random_float32(), 1 ) ;
                float blue = std::fmod( GEO::Numeric::random_float32(), 1 ) ;
                const GeoModelEntity& cur_layer = instance()->GM_->layer( l ) ;
                for( index_t r = 0; r < cur_layer.nb_children(); ++r )
                    instance()->GM_gfx_.set_cell_region_color( cur_layer.child( r ).index(),
                        red, green, blue ) ;
            }
        }
    }

    bool RINGMeshApplication::load( const std::string& filename )
    {
        double xyzmin[3] ;
        double xyzmax[3] ;
        for( GEO::index_t c = 0; c < 3; c++ ) {
            xyzmin[c] = GEO::Numeric::max_float64() ;
            xyzmax[c] = GEO::Numeric::min_float64() ;
        }

        if( !filename.empty() ) {
            try {
                if( GM_ ) {
                    delete GM_ ;
                }
                GM_ = new GeoModel ;
                geomodel_load( *GM_, filename ) ;
                meshed_regions_ = GM_->region( 0 ).is_meshed() ;
                if( meshed_regions_ ) {
                    show_volume_ = true ;
                }

            } catch( const RINGMeshException& e ) {
                GEO::Logger::err( e.category() ) << e.what() << std::endl ;
                return false ;
            }
        } else {
            GEO::Logger::err( "I/O" ) << "Give at least a filename in geomodel"
                << std::endl ;
            return false ;
        }
        GM_gfx_.set_geo_model( *GM_ ) ;

        for( GEO::index_t s = 0; s < GM_->nb_surfaces(); s++ ) {
            const RINGMesh::Surface& S = GM_->surface( s ) ;
            for( GEO::index_t v = 0; v < S.nb_vertices(); ++v ) {
                const vec3& p = S.vertex( v ) ;
                for( GEO::coord_index_t c = 0; c < 3; ++c ) {
                    xyzmin[c] = GEO::geo_min( xyzmin[c], p[c] ) ;
                    xyzmax[c] = GEO::geo_max( xyzmax[c], p[c] ) ;
                }
            }
        }

        glup_viewer_set_region_of_interest( float( xyzmin[0] ), float( xyzmin[1] ),
            float( xyzmin[2] ), float( xyzmax[0] ), float( xyzmax[1] ),
            float( xyzmax[2] ) ) ;

        return true ;
    }

    void RINGMeshApplication::draw_scene()
    {
        if( white_bg_ ) {
            GM_gfx_.set_surfaces_color( 0.9f, 0.9f, 0.9f ) ;

            GM_gfx_.set_edge_lines_color( 0.0f, 0.0f, 0.0f ) ;
            GM_gfx_.set_mesh_surfaces_color( 0.0f, 0.0f, 0.0f ) ;
            GM_gfx_.set_cell_mesh_regions_color( 0.0f, 0.0f, 0.0f ) ;
            GM_gfx_.set_edge_regions_color( 0.0f, 0.0f, 0.0f ) ;

            GM_gfx_.set_vertex_regions_color( 0.0f, 0.0f, 0.0f ) ;
        } else {
            GM_gfx_.set_surfaces_color( 0.1f, 0.1f, 0.1f ) ;

            GM_gfx_.set_edge_lines_color( 1.0f, 1.0f, 1.0f ) ;
            GM_gfx_.set_mesh_surfaces_color( 1.0f, 1.0f, 1.0f ) ;
            GM_gfx_.set_cell_mesh_regions_color( 1.0f, 1.0f, 1.0f ) ;
            GM_gfx_.set_edge_regions_color( 1.0f, 1.0f, 1.0f ) ;

            GM_gfx_.set_vertex_regions_color( 1.0f, 1.0f, 1.0f ) ;
        }

        GM_gfx_.set_mesh_surfaces_visibility( mesh_visible_ ) ;
        GM_gfx_.set_cell_mesh_regions_visibility( mesh_visible_ ) ;

        if( show_attributes_ ) {
            switch( GM_gfx_.attribute_location() ) {
                case GeoModelGfx::cells:
                    GM_gfx_.set_cell_attribute( GM_gfx_.attribute_name(),
                        GM_gfx_.attribute_coordinate(), GM_gfx_.colormap() ) ;
                    break ;
                case GeoModelGfx::cell_vertices:
                    GM_gfx_.set_cell_vertex_attribute( GM_gfx_.attribute_name(),
                        GM_gfx_.attribute_coordinate(), GM_gfx_.colormap() ) ;
                    break ;
            }
        } else {
            switch( GM_gfx_.attribute_location() ) {
                case GeoModelGfx::cells:
                    GM_gfx_.set_cell_attribute( "", 0, GM_gfx_.colormap() ) ;
                    break ;
                case GeoModelGfx::cell_vertices:
                    GM_gfx_.set_cell_vertex_attribute( "", 0, GM_gfx_.colormap() ) ;
                    break ;
            }
        }

        if( show_corners_ ) {
            GM_gfx_.draw_corners() ;
        }

        if( show_lines_ ) {
            GM_gfx_.draw_lines() ;
        }

        if( show_surface_ ) {
            for( GEO::index_t s = 0; s < GM_->nb_surfaces(); s++ ) {
                if( GM_->surface( s ).is_on_voi() ) {
                    GM_gfx_.set_surface_visibility( s, show_voi_ ) ;
                }
            }
            GM_gfx_.draw_surfaces() ;
        }

        if( show_volume_ && meshed_regions_ ) {
            if( !colored_cells_ && !show_colored_regions_ && !show_colored_layers_ ) {
                GM_gfx_.set_cell_regions_color( 0.9f, 0.9f, 0.9f ) ;
            }
            GM_gfx_.set_cell_regions_shrink( shrink_ ) ;
            GM_gfx_.draw_regions() ;
        }
    }

    std::string RINGMeshApplication::supported_read_file_extensions()
    {
        return file_extensions_ ;
    }

    void RINGMeshApplication::set_attribute_names(
        const GEO::AttributesManager& attributes )
    {
        GEO::vector< std::string > attribute_names ;
        attributes.list_attribute_names( attribute_names ) ;
        for( index_t i = 0; i < attribute_names.size(); ++i ) {
            const GEO::AttributeStore* store = attributes.find_attribute_store(
                attribute_names[i] ) ;
            if( GEO::ReadOnlyScalarAttributeAdapter::can_be_bound_to( store ) ) {
                if( ImGui::Button( attribute_names[i].c_str() ) ) {
                    GM_gfx_.set_attribute_name( attribute_names[i] ) ;
                    GM_gfx_.set_attribute_coordinate( 0 ) ;
                    autorange() ;
                    ImGui::CloseCurrentPopup() ;
                }
            }
        }
    }
    void RINGMeshApplication::autorange()
    {
        switch( GM_gfx_.attribute_location() ) {
            case GeoModelGfx::cells:
                GM_gfx_.compute_cell_attribute_range( GM_gfx_.attribute_coordinate(),
                    GM_gfx_.attribute_name() ) ;
                break ;
            case GeoModelGfx::cell_vertices:
                GM_gfx_.compute_cell_vertex_attribute_range(
                    GM_gfx_.attribute_coordinate(), GM_gfx_.attribute_name() ) ;
                break ;
        }
        attribute_max_ = GM_gfx_.attribute_max() ;
        attribute_min_ = GM_gfx_.attribute_min() ;
    }

    index_t RINGMeshApplication::nb_coordinates() const
    {
        GEO::AttributeStore* store = nil ;
        switch( GM_gfx_.attribute_location() ) {
            case GeoModelGfx::cells:
                store =
                    GM_->region( 0 ).cell_attribute_manager().find_attribute_store(
                        GM_gfx_.attribute_name() ) ;
                break ;
            case GeoModelGfx::cell_vertices:
                store =
                    GM_->region( 0 ).vertex_attribute_manager().find_attribute_store(
                        GM_gfx_.attribute_name() ) ;
                break ;
            default:
                return false ;
        }
        if( store == nil ) return 0 ;
        return store->dimension() ;
    }

    void RINGMeshApplication::draw_object_properties()
    {
        ImGui::Checkbox( "Attributes", &show_attributes_ ) ;
        if( show_attributes_ ) {
            if( ImGui::Button(
                GM_gfx_.attribute_location_name( GM_gfx_.attribute_location() ).c_str(),
                ImVec2( -1, 0 ) ) ) {
                ImGui::OpenPopup( "##Locations" ) ;
            }
            if( ImGui::BeginPopup( "##Locations" ) ) {
                for( index_t i = 0; i < GeoModelGfx::nb_locations; i++ ) {
                    GeoModelGfx::GeoModelAttribute_location l =
                        static_cast< GeoModelGfx::GeoModelAttribute_location >( i ) ;
                    if( ImGui::Button(
                        GeoModelGfx::attribute_location_name( l ).c_str() ) ) {
                        GM_gfx_.set_attribute_location( l ) ;
                        reset_attribute_name() ;
                        ImGui::CloseCurrentPopup() ;
                    }
                }
                ImGui::EndPopup() ;
            }

            if( ImGui::Button( GM_gfx_.attribute_name().c_str(), ImVec2( -1, 0 ) ) ) {
                ImGui::OpenPopup( "##Attributes" ) ;
            }
            if( ImGui::BeginPopup( "##Attributes" ) ) {
                switch( GM_gfx_.attribute_location() ) {
                    case GeoModelGfx::cells:
                        set_attribute_names( GM_->region( 0 ).cell_attribute_manager() ) ;
                        break ;
                    case GeoModelGfx::cell_vertices :
                        set_attribute_names( GM_->region( 0 ).vertex_attribute_manager() ) ;
                        break ;
                }
                ImGui::EndPopup() ;
            }
            if( nb_coordinates() > 1 ) {
                if( ImGui::Button(
                    GEO::String::to_string( GM_gfx_.attribute_coordinate() ).c_str(),
                    ImVec2( -1, 0 ) ) ) {
                    ImGui::OpenPopup( "##Coordinates" ) ;
                }
                if( ImGui::BeginPopup( "##Coordinates" ) ) {
                    for( index_t i = 0; i < nb_coordinates(); i++ ) {
                        if( ImGui::Button( GEO::String::to_string( i ).c_str() ) ) {
                            GM_gfx_.set_attribute_coordinate( i ) ;
                            autorange() ;
                            ImGui::CloseCurrentPopup() ;
                        }
                    }
                    ImGui::EndPopup() ;
                }
            }
            if( ImGui::InputFloat( "min", &attribute_min_ ) ) {
                GM_gfx_.set_attribute_min( attribute_min_ ) ;
            }
            if( ImGui::InputFloat( "max", &attribute_max_ ) ) {
                GM_gfx_.set_attribute_max( attribute_max_ ) ;
            }
            if( ImGui::Button( "autorange", ImVec2( -1, 0 ) ) ) {
                autorange() ;
            }
            if( ImGui::ImageButton(
                convert_to_ImTextureID( GM_gfx_.colormap() ),
                ImVec2( 115, 8 ) ) ) {
                ImGui::OpenPopup( "##Colormap" ) ;
            }
            if( ImGui::BeginPopup( "##Colormap" ) ) {
                for( index_t i = 0; i < colormaps_.size(); ++i ) {
                    if( ImGui::ImageButton(
                        convert_to_ImTextureID( colormaps_[i].texture ),
                        ImVec2( 100, 8 ) ) ) {
                        GM_gfx_.set_colormap( colormaps_[i].texture ) ;
                        ImGui::CloseCurrentPopup() ;
                    }
                }
                ImGui::EndPopup() ;
            }
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
                ImGui::Checkbox( "Col. cells [C]", &colored_cells_ ) ;
                ImGui::Checkbox( "Col. regions [r]", &show_colored_regions_ ) ;
                ImGui::Checkbox( "Col. layers [R]", &show_colored_layers_ ) ;
                ImGui::SliderFloat( "Shrk.", &shrink_, 0.0f, 1.0f, "%.1f" ) ;
            }
        }
    }

}

#endif
