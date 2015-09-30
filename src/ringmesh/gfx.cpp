/*
 * Copyright (c) 2012-2015, Association Scientifique pour la Geologie et ses Applications (ASGA)
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *     * Redistributions in binary form must reproduce the above copyright
 *       notice, this list of conditions and the following disclaimer in the
 *       documentation and/or other materials provided with the distribution.
 *     * Neither the name of the <organization> nor the
 *       names of its contributors may be used to endorse or promote products
 *       derived from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL <COPYRIGHT HOLDER> BE LIABLE FOR ANY
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 *  Contacts:
 *     Arnaud.Botella@univ-lorraine.fr
 *     Antoine.Mazuyer@univ-lorraine.fr
 *     Jeanne.Pellerin@wias-berlin.de
 *
 *     http://www.ring-team.org
 *
 *     RING Project
 *     Ecole Nationale Supérieure de Géologie - Georessources
 *     2 Rue du Doyen Marcel Roubault - TSA 70605
 *     54518 VANDOEUVRE-LES-NANCY
 *     FRANCE
 */

/*! \author Jeanne Pellerin and Arnaud Botella */

#ifdef RINGMESH_WITH_GRAPHICS

#include <ringmesh/gfx.h>
#include <ringmesh/geo_model.h>
#include <ringmesh/geo_model_element.h>
#include <ringmesh/macro_mesh.h>

#include <geogram/basic/logger.h>

namespace RINGMesh {

    class MeshElementGfx {
    ringmesh_disable_copy( MeshElementGfx ) ;
    public:
        MeshElementGfx( const GEO::Mesh& mesh, bool vertice_visible )
            : vertices_visible_( vertice_visible )
        {
            gfx_.set_mesh( &mesh ) ;
        }
        virtual ~MeshElementGfx(){}

        GEO::MeshGfx& gfx()
        {
            return gfx_ ;
        }

        virtual void draw_vertices() {
            gfx_.draw_vertices() ;
        }
        virtual void draw_edges() {
            index_t w = gfx_.get_mesh_width() ;
            gfx_.set_mesh_width( w+1 ) ;
            gfx_.draw_edges() ;
            gfx_.set_mesh_width( w ) ;
        }
        virtual void draw_surface() {
            gfx_.draw_surface() ;
        }
        virtual void draw_volume() {
            gfx_.draw_volume() ;
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
        GEO::MeshGfx gfx_ ;

        bool vertices_visible_ ;
    } ;

    class CornerGfx: public MeshElementGfx {
    public:
        CornerGfx( const Corner& corner )
            : MeshElementGfx( corner.mesh(), true )
        {
            gfx_.set_points_color( 1, 0, 0 ) ;
        }
    } ;

    class LineGfx: public MeshElementGfx {
    public:
        LineGfx( const Line& line )
            : MeshElementGfx( line.mesh(), false ), edges_visible_( true )
        {
            gfx_.set_points_color( 1, 1, 1 ) ;
            gfx_.set_mesh_color( 1, 1, 1 ) ;
        }
        LineGfx( const GEO::Mesh& mesh )
            : MeshElementGfx( mesh, false ), edges_visible_( true )
        {
            gfx_.set_points_color( 1, 1, 1 ) ;
            gfx_.set_mesh_color( 1, 1, 1 ) ;
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

    class SurfaceGfx: public MeshElementGfx {
    public:
        SurfaceGfx( const Surface& surface )
            : MeshElementGfx( surface.mesh(), false ), surface_visible_( true )
        {
        }
        SurfaceGfx( const GEO::Mesh& mesh )
            : MeshElementGfx( mesh, false ), surface_visible_( true )
        {
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

    class RegionGfx: public MeshElementGfx {
    public:
        RegionGfx( const GEO::Mesh& mesh )
            : MeshElementGfx( mesh, false ), region_visible_( true )
        {
        }

//        void set_edges_visible( bool b )
//        {
//            edges_visible_ = b ;
//        }
//        bool get_edges_visible() const
//        {
//            return edges_visible_ ;
//        }
//        void set_surface_visible( bool b )
//        {
//            surface_visible_ = b ;
//        }
//        bool get_surface_visible() const
//        {
//            return surface_visible_ ;
//        }
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

    } ;

    GeoModelGfx::GeoModelGfx()
        : model_( nil ), corners_(), lines_(), surfaces_()
    {
    }

    GeoModelGfx::~GeoModelGfx()
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
    const GeoModel* GeoModelGfx::geo_model()
    {
        return model_ ;
    }

    /*!
     * Initializes the database according the GeoModel dimensions
     */
    void GeoModelGfx::initialize()
    {
        ringmesh_debug_assert( model_ ) ;
        if( corners_.empty() && lines_.empty() && surfaces_.empty() ) {
            corners_.resize( model_->nb_corners(), nil ) ;
            lines_.resize( model_->nb_lines(), nil ) ;
            surfaces_.resize( model_->nb_surfaces(), nil ) ;

            for( index_t c = 0; c < corners_.size(); c++ ) {
                corners_[c] = new CornerGfx( model_->corner( c ) ) ;
            }
            for( index_t l = 0; l < lines_.size(); l++ ) {
                lines_[l] = new LineGfx( model_->line( l ) ) ;
            }
            for( index_t s = 0; s < surfaces_.size(); s++ ) {
                surfaces_[s] = new SurfaceGfx( model_->surface( s ) ) ;
            }
        }
    }

    /*!
     * Draws the corners
     */
    void GeoModelGfx::draw_corners()
    {
        GEO::Logger::instance()->set_quiet( true ) ;
        for( index_t c = 0; c < corners_.size(); c++ ) {
            if( corners_[c]->get_vertices_visible() )
                corners_[c]->draw_vertices() ;
        }
        GEO::Logger::instance()->set_quiet( false ) ;
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
        corners_[c]->gfx().set_points_color( r, g, b ) ;
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
        corners_[c]->gfx().set_points_size( s ) ;
    }

    /*!
     * Draws the lines
     */
    void GeoModelGfx::draw_lines()
    {
        GEO::Logger::instance()->set_quiet( true ) ;
        for( index_t l = 0; l < lines_.size(); l++ ) {
            if( lines_[l]->get_vertices_visible() )
                lines_[l]->draw_vertices() ;
            if( lines_[l]->get_edges_visible() ) lines_[l]->draw_edges() ;
        }
        GEO::Logger::instance()->set_quiet( false ) ;
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
    void GeoModelGfx::set_edge_line_color(
        index_t l,
        float r,
        float g,
        float b )
    {
        lines_[l]->gfx().set_mesh_color( r, g, b ) ;
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
        lines_[l]->gfx().set_mesh_width( s ) ;
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
    void GeoModelGfx::set_vertex_line_color(
        index_t l,
        float r,
        float g,
        float b )
    {
        lines_[l]->gfx().set_points_color( r, g, b ) ;
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
        lines_[l]->gfx().set_points_size( s ) ;
    }

    /*!
     * Draws the surfaces
     */
    void GeoModelGfx::draw_surfaces()
    {
        GEO::Logger::instance()->set_quiet( true ) ;
        for( index_t s = 0; s < surfaces_.size(); s++ ) {
            if( surfaces_[s]->get_vertices_visible() )
                surfaces_[s]->draw_vertices() ;
            if( surfaces_[s]->get_surface_visible() )
                surfaces_[s]->draw_surface() ;
        }
        GEO::Logger::instance()->set_quiet( false ) ;
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
        surfaces_[s]->gfx().set_surface_color( r, g, b ) ;
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
        surfaces_[s]->gfx().set_backface_surface_color( r, g, b ) ;
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
    void GeoModelGfx::set_mesh_surface_color(
        index_t s,
        float r,
        float g,
        float b )
    {
        surfaces_[s]->gfx().set_mesh_color( r, g, b ) ;
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
        surfaces_[s]->gfx().set_show_mesh( b ) ;
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
        surfaces_[s]->gfx().set_mesh_width( size ) ;
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
        surfaces_[s]->gfx().set_points_color( r, g, b ) ;
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
        surfaces_[s]->gfx().set_points_size( size ) ;
    }

    /*!
     * Constructor. Does nothing.
     */
    MacroMeshGfx::MacroMeshGfx()
        : mm_( nil ), meshes_(), surfaces_(), edges_()
    {
    }

    /*!
     * Destructor. Deletes the RegionGfx, SurfaceGfx and LineGfx
     * inside meshes_, surfaces_ and edges_.
     */
    MacroMeshGfx::~MacroMeshGfx()
    {
        for( index_t m = 0; m < meshes_.size(); m++ ) {
            delete meshes_[m] ;
            delete surfaces_[m] ;
            delete edges_[m] ;
        }
    }

    /*!
     * Sets the MacroMesh associated to the graphics
     * @param[in] mm the MacroMesh
     */
    void MacroMeshGfx::set_macro_mesh( const MacroMesh& mm )
    {
        mm_ = &mm ;
        initialize() ;
    }

    /*!
     * Initializes the database according the MacroMesh dimensions
     */
    void MacroMeshGfx::initialize()
    {
        ringmesh_debug_assert( mm_ ) ;
        if( meshes_.empty() ) {
            meshes_.resize( mm_->nb_meshes(), nil ) ;
            surfaces_.resize( mm_->nb_meshes(), nil ) ;
            edges_.resize( mm_->nb_meshes(), nil ) ;

            for( index_t m = 0; m < meshes_.size(); m++ ) {
                meshes_[m] = new RegionGfx( mm_->mesh( m ) ) ;
                surfaces_[m] = new SurfaceGfx( mm_->mesh( m ) ) ;
                edges_[m] = new LineGfx( mm_->mesh( m ) ) ;
            }
        }
        set_surface_regions_visibility( false ) ;
    }

    /*!
     * Draws the MacroMesh
     */
    void MacroMeshGfx::draw()
    {
        for( index_t m = 0; m < meshes_.size(); m++ ) {
            if( meshes_[m]->get_vertices_visible() )
                meshes_[m]->draw_vertices() ;
            if( edges_[m]->get_edges_visible() ) edges_[m]->draw_edges() ;
            if( surfaces_[m]->get_surface_visible() )
                surfaces_[m]->draw_surface() ;
            if( meshes_[m]->get_region_visible() ) meshes_[m]->draw_volume() ;
        }
    }

    /*!
     * Sets the vertex region color to all the regions
     * @param[in] r the red component of the color in [0.0, 1.0]
     * @param[in] g the green component of the color in [0.0, 1.0]
     * @param[in] b the blue component of the color in [0.0, 1.0]
     */
    void MacroMeshGfx::set_vertex_regions_color( float r, float g, float b )
    {
        for( index_t m = 0; m < meshes_.size(); m++ ) {
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
    void MacroMeshGfx::set_vertex_region_color(
        index_t m,
        float r,
        float g,
        float b )
    {
        ringmesh_debug_assert( m < meshes_.size() ) ;
        meshes_[m]->gfx().set_points_color( r, g, b ) ;
    }

    /*!
     * Sets the vertex region visibility to all the regions
     * @param[in] b the visibility
     */
    void MacroMeshGfx::set_vertex_regions_visibility( bool b )
    {
        for( index_t m = 0; m < meshes_.size(); m++ ) {
            set_vertex_region_visibility( m, b ) ;
        }
    }

    /*!
     * Sets the vertex region visibility
     * @param[in] m the region index
     * @param[in] b the visibility
     */
    void MacroMeshGfx::set_vertex_region_visibility( index_t m, bool b )
    {
        ringmesh_debug_assert( m < meshes_.size() ) ;
        meshes_[m]->set_vertices_visible( b ) ;
    }

    /*!
     * Sets the vertex region size to all the regions
     * @param[in] s the size
     */
    void MacroMeshGfx::set_vertex_regions_size( index_t s )
    {
        for( index_t m = 0; m < meshes_.size(); m++ ) {
            set_vertex_region_size( m, s ) ;
        }
    }

    /*!
     * Sets the vertex region size
     * @param[in] m the region index
     * @param[in] s the size
     */
    void MacroMeshGfx::set_vertex_region_size( index_t m, index_t s )
    {
        ringmesh_debug_assert( m < meshes_.size() ) ;
        meshes_[m]->gfx().set_points_size( s ) ;
    }

    /*!
     * Sets the edge color to all the meshes
     * @param[in] r the red component of the color in [0.0, 1.0]
     * @param[in] g the green component of the color in [0.0, 1.0]
     * @param[in] b the blue component of the color in [0.0, 1.0]
     */
    void MacroMeshGfx::set_edge_regions_color( float r, float g, float b )
    {
        for( index_t k = 0; k < edges_.size(); k++ ) {
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
    void MacroMeshGfx::set_edge_region_color( index_t m, float r, float g, float b )
    {
        ringmesh_debug_assert( m < meshes_.size() ) ;
        edges_[m]->gfx().set_mesh_color( r, g, b ) ; //TODO function not good?
    }
    /*!
     * Sets the edge visibility to all the meshes
     * @param[in] b the visibility
     */
    void MacroMeshGfx::set_edge_regions_visibility( bool b )
    {
        for( index_t m = 0; m < edges_.size(); m++ ) {
            set_edge_region_visibility( m, b ) ;
        }
    }
    /*!
     * Sets the edge visibility
     * @param[in] m the mesh index
     * @param[in] b the visibility
     */
    void MacroMeshGfx::set_edge_region_visibility( index_t m, bool b )
    {
        ringmesh_debug_assert( m < meshes_.size() ) ;
        edges_[m]->set_edges_visible( b ) ;
    }
    /*!
     * Sets the edge line size to all the lines
     * @param[in] s the size
     */
    void MacroMeshGfx::set_edge_regions_size( index_t s )
    {
        for( index_t m = 0; m < edges_.size(); m++ ) {
            set_edge_region_size( m, s ) ;
        }
    }
    /*!
     * Sets the edge line size
     * @param[in] l the line index
     * @param[in] s the size
     */
    void MacroMeshGfx::set_edge_region_size( index_t l, index_t s )
    {
        ringmesh_debug_assert( l < meshes_.size() ) ;
        edges_[l]->gfx().set_mesh_width( s ) ;
    }

    /*!
     * Sets the surface color to all the surfaces
     * @param[in] r the red component of the color in [0.0, 1.0]
     * @param[in] g the green component of the color in [0.0, 1.0]
     * @param[in] b the blue component of the color in [0.0, 1.0]
     */
    void MacroMeshGfx::set_surface_regions_color( float r, float g, float b )
    {
        for( index_t s = 0; s < surfaces_.size(); s++ ) {
            set_surface_region_color( s, r, g, b ) ;
        }
    }
    /*!
     * Sets the surface color
     * @param[in] s the surface index
     * @param[in] r the red component of the color in [0.0, 1.0]
     * @param[in] g the green component of the color in [0.0, 1.0]
     * @param[in] b the blue component of the color in [0.0, 1.0]
     */
    void MacroMeshGfx::set_surface_region_color(
        index_t s,
        float r,
        float g,
        float b )
    {
        ringmesh_debug_assert( s < meshes_.size() ) ;
        surfaces_[s]->gfx().set_surface_color( r, g, b ) ;
    }
    /*!
     * Sets the backface surface color to all the surfaces
     * @param[in] r the red component of the color in [0.0, 1.0]
     * @param[in] g the green component of the color in [0.0, 1.0]
     * @param[in] b the blue component of the color in [0.0, 1.0]
     */
    void MacroMeshGfx::set_backface_surface_regions_color(
        float r,
        float g,
        float b )
    {
        for( index_t s = 0; s < surfaces_.size(); s++ ) {
            set_backface_surface_region_color( s, r, g, b ) ;
        }
    }
    /*!
     * Sets the backsurface surface color
     * @param[in] s the surface index
     * @param[in] r the red component of the color in [0.0, 1.0]
     * @param[in] g the green component of the color in [0.0, 1.0]
     * @param[in] b the blue component of the color in [0.0, 1.0]
     */
    void MacroMeshGfx::set_backface_surface_region_color(
        index_t s,
        float r,
        float g,
        float b )
    {
        ringmesh_debug_assert( s < meshes_.size() ) ;
        surfaces_[s]->gfx().set_backface_surface_color( r, g, b ) ;
    }
    /*!
     * Sets the surface visibility to all the surfaces
     * @param[in] b the visibility
     */
    void MacroMeshGfx::set_surface_regions_visibility( bool b )
    {
        for( index_t s = 0; s < surfaces_.size(); s++ ) {
            set_surface_region_visibility( s, b ) ;
        }
    }
    /*!
     * Sets the surface visibility
     * @param[in] s the surface index
     * @param[in] b the visibility
     */
    void MacroMeshGfx::set_surface_region_visibility( index_t s, bool b )
    {
        ringmesh_debug_assert( s < meshes_.size() ) ;
        surfaces_[s]->set_surface_visible( b ) ;
    }
    /*!
     * Sets the mesh surface color to all the surfaces
     * @param[in] r the red component of the color in [0.0, 1.0]
     * @param[in] g the green component of the color in [0.0, 1.0]
     * @param[in] b the blue component of the color in [0.0, 1.0]
     */
    void MacroMeshGfx::set_mesh_surface_regions_color( float r, float g, float b )
    {
        for( index_t s = 0; s < surfaces_.size(); s++ ) {
            set_mesh_surface_region_color( s, r, g, b ) ;
        }
    }
    /*!
     * Sets the mesh surface color
     * @param[in] s the surface index
     * @param[in] r the red component of the color in [0.0, 1.0]
     * @param[in] g the green component of the color in [0.0, 1.0]
     * @param[in] b the blue component of the color in [0.0, 1.0]
     */
    void MacroMeshGfx::set_mesh_surface_region_color(
        index_t s,
        float r,
        float g,
        float b )
    {
        ringmesh_debug_assert( s < meshes_.size() ) ;
        surfaces_[s]->gfx().set_mesh_color( r, g, b ) ;
    }
    /*!
     * Sets the mesh surface visibility to all the surfaces
     * @param[in] b the visibility
     */
    void MacroMeshGfx::set_mesh_surface_regions_visibility( bool b )
    {
        for( index_t s = 0; s < surfaces_.size(); s++ ) {
            set_mesh_surface_region_visibility( s, b ) ;
        }
    }
    /*!
     * Sets the mesh surface visibility
     * @param[in] s the surface index
     * @param[in] b the visibility
     */
    void MacroMeshGfx::set_mesh_surface_region_visibility( index_t s, bool b )
    {
        ringmesh_debug_assert( s < meshes_.size() ) ;
        surfaces_[s]->gfx().set_show_mesh( b ) ;
    }
    /*!
     * Sets the mesh surface size to all the surfaces
     * @param[in] size the size
     */
    void MacroMeshGfx::set_mesh_surface_regions_size( index_t size )
    {
        for( index_t s = 0; s < surfaces_.size(); s++ ) {
            set_mesh_surface_region_size( s, size ) ;
        }
    }
    /*!
     * Sets the mesh surface size
     * @param[in] s the surface index
     * @param[in] size the size
     */
    void MacroMeshGfx::set_mesh_surface_region_size( index_t s, index_t size )
    {
        ringmesh_debug_assert( s < meshes_.size() ) ;
        surfaces_[s]->gfx().set_mesh_width( size ) ;
    }

    /*!
     * Sets the mesh region color to all the regions
     * @param[in] r the red component of the color in [0.0, 1.0]
     * @param[in] g the green component of the color in [0.0, 1.0]
     * @param[in] b the blue component of the color in [0.0, 1.0]
     */
    void MacroMeshGfx::set_cell_mesh_regions_color( float r, float g, float b )
    {
        for( index_t m = 0; m < meshes_.size(); m++ ) {
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
    void MacroMeshGfx::set_cell_mesh_region_color(
        index_t m,
        float r,
        float g,
        float b )
    {
        ringmesh_debug_assert( m < meshes_.size() ) ;
        meshes_[m]->gfx().set_mesh_color( r, g, b ) ;
    }

    /*!
     * Toggles the cell region color per cell type to all the regions
     */
    void MacroMeshGfx::set_cell_regions_color_type()
    {
        for( index_t m = 0; m < meshes_.size(); m++ ) {
            set_cell_region_color_type( m ) ;
        }
    }

    /*!
     * Toggles the cell region color per cell type
     * @param[in] m the region index
     */
    void MacroMeshGfx::set_cell_region_color_type( index_t m )
    {
        ringmesh_debug_assert( m < meshes_.size() ) ;
        meshes_[m]->gfx().set_cells_colors_by_type() ;
    }

    /*!
     * Sets the mesh region visibility to all the regions
     * @param[in] b the visibility
     */
    void MacroMeshGfx::set_cell_mesh_regions_visibility( bool b )
    {
        for( index_t m = 0; m < meshes_.size(); m++ ) {
            set_cell_mesh_region_visibility( m, b ) ;
        }
    }

    /*!
     * Sets the mesh region visibility
     * @param[in] m the region index
     * @param[in] b the visibility
     */
    void MacroMeshGfx::set_cell_mesh_region_visibility( index_t m, bool b )
    {
        ringmesh_debug_assert( m < meshes_.size() ) ;
        meshes_[m]->gfx().set_show_mesh( b ) ;
    }

    /*!
     * Sets the mesh region size to all the regions
     * @param[in] s the size
     */
    void MacroMeshGfx::set_cell_mesh_regions_size( index_t s )
    {
        for( index_t m = 0; m < meshes_.size(); m++ ) {
            set_cell_mesh_region_size( m, s ) ;
        }
    }

    /*!
     * Sets the mesh region size
     * @param[in] m the region index
     * @param[in] s the size
     */
    void MacroMeshGfx::set_cell_mesh_region_size( index_t m, index_t s )
    {
        ringmesh_debug_assert( m < meshes_.size() ) ;
        meshes_[m]->gfx().set_mesh_width( s ) ;
    }

    /*!
     * Sets the cell region color to all the regions
     * @param[in] r the red component of the color in [0.0, 1.0]
     * @param[in] g the green component of the color in [0.0, 1.0]
     * @param[in] b the blue component of the color in [0.0, 1.0]
     */
    void MacroMeshGfx::set_cell_regions_color( float r, float g, float b )
    {
        for( index_t m = 0; m < meshes_.size(); m++ ) {
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
    void MacroMeshGfx::set_cell_region_color( index_t m, float r, float g, float b )
    {
        ringmesh_debug_assert( m < meshes_.size() ) ;
        meshes_[m]->gfx().set_cells_color( r, g, b ) ;
    }

    /*!
     * Sets the cell region visibility to all the regions
     * @param[in] b the visibility
     */
    void MacroMeshGfx::set_cell_regions_visibility( bool b )
    {
        for( index_t m = 0; m < meshes_.size(); m++ ) {
            set_cell_region_visibility( m, b ) ;
        }
    }

    /*!
     * Sets the cell region visibility to all the regions
     * @param[in] m the region index
     * @param[in] b the visibility
     */
    void MacroMeshGfx::set_cell_region_visibility( index_t m, bool b )
    {
        ringmesh_debug_assert( m < meshes_.size() ) ;
        meshes_[m]->set_region_visible( b ) ;
    }
    void MacroMeshGfx::set_cell_regions_type_visibility(
        GEO::MeshCellType t,
        bool b )
    {
        for( index_t m = 0; m < meshes_.size(); m++ ) {
            set_cell_region_type_visibility( m, t, b ) ;
        }
    }
    void MacroMeshGfx::set_cell_region_type_visibility(
        index_t m,
        GEO::MeshCellType t,
        bool b )
    {
        ringmesh_debug_assert( m < meshes_.size() ) ;
        meshes_[m]->gfx().set_draw_cells( t, b ) ;
    }

    /*!
     * Sets the cell region shrink to all the regions
     * @param[in] s the shrink
     */
    void MacroMeshGfx::set_cell_regions_shrink( double s )
    {
        for( index_t m = 0; m < meshes_.size(); m++ ) {
            set_cell_region_shrink( m, s ) ;
        }
    }

    /*!
     * Sets the cell region shrink
     * @param[in] m the region index
     * @param[in] s the shrink
     */
    void MacroMeshGfx::set_cell_region_shrink( index_t m, double s )
    {
        ringmesh_debug_assert( m < meshes_.size() ) ;
        meshes_[m]->gfx().set_shrink( s ) ;
    }

} // namespace

#endif
