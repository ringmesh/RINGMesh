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
 *     http://www.gocad.org
 *
 *     GOCAD Project
 *     Ecole Nationale Supérieure de Géologie - Georessources
 *     2 Rue du Doyen Marcel Roubault - TSA 70605
 *     54518 VANDOEUVRE-LES-NANCY
 *     FRANCE
 */

/*! \author Jeanne Pellerin and Arnaud Botella */

#ifdef RINGMESH_WITH_GRAPHICS

#include <ringmesh/gfx.h>
#include <ringmesh/boundary_model.h>
#include <ringmesh/boundary_model_element.h>

namespace RINGMesh {

    class RINGMESH_API BoundaryModelMeshElementGfx {
    ringmesh_disable_copy( BoundaryModelMeshElementGfx ) ;
    public:
        BoundaryModelMeshElementGfx( const BoundaryModelMeshElement& mesh )
            : vertices_visible_( false )
        {
            gfx_.set_mesh( &mesh.mesh() ) ;
        }

        GEO::MeshGfx& gfx()
        {
            return gfx_ ;
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

    class RINGMESH_API CornerGfx: public BoundaryModelMeshElementGfx {
    public:
        CornerGfx( const Corner& corner )
            : BoundaryModelMeshElementGfx( corner )
        {
            vertices_visible_ = true ;
            gfx_.set_points_color( 1, 0, 0 ) ;
        }
    } ;

    class RINGMESH_API LineGfx: public BoundaryModelMeshElementGfx {
    public:
        LineGfx( const Line& line )
            : BoundaryModelMeshElementGfx( line )
        {
            vertices_visible_ = false ;
            gfx_.set_points_color( 1, 1, 1 ) ;
            edges_visible_ = true ;
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

    class RINGMESH_API SurfaceGfx: public BoundaryModelMeshElementGfx {
    public:
        SurfaceGfx( const Surface& surface )
            : BoundaryModelMeshElementGfx( surface )
        {
            vertices_visible_ = false ;
            surface_visible_ = true ;
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


    BoundaryModelGfx::BoundaryModelGfx()
        : model_( nil )
    {
    }

    void BoundaryModelGfx::set_boundary_model( const BoundaryModel& model )
    {
        model_ = &model ;
        initialize() ;
    }

    void BoundaryModelGfx::initialize()
    {
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

    void BoundaryModelGfx::draw_corners()
    {
        for( index_t c = 0; c < corners_.size(); c++ ) {
            if( corners_[c]->get_vertices_visible() ) corners_[c]->gfx().draw_vertices() ;
        }
    }
    void BoundaryModelGfx::set_corners_color( float r, float g, float b )
    {
        for( index_t c = 0; c < corners_.size(); c++ ) {
            set_corner_color( c, r, g, b ) ;
        }
    }
    void BoundaryModelGfx::set_corner_color( index_t c, float r, float g, float b )
    {
        corners_[c]->gfx().set_points_color( r, g, b ) ;
    }
    void BoundaryModelGfx::set_corners_visibility( bool b )
    {
        for( index_t c = 0; c < corners_.size(); c++ ) {
            set_corner_visibility( c, b ) ;
        }
    }
    void BoundaryModelGfx::set_corner_visibility( index_t c, bool b )
    {
        corners_[c]->set_vertices_visible( b ) ;
    }
    void BoundaryModelGfx::set_corners_size( index_t s )
    {
        for( index_t c = 0; c < corners_.size(); c++ ) {
            set_corner_size( c, s ) ;
        }
    }
    void BoundaryModelGfx::set_corner_size( index_t c, index_t s )
    {
        corners_[c]->gfx().set_points_size( s ) ;
    }

    void BoundaryModelGfx::draw_lines()
    {
        for( index_t l = 0; l < lines_.size(); l++ ) {
            if( lines_[l]->get_vertices_visible() ) lines_[l]->gfx().draw_vertices() ;
            if( lines_[l]->get_edges_visible() ) lines_[l]->gfx().draw_edges() ;
        }
    }
    void BoundaryModelGfx::set_edge_lines_color( float r, float g, float b )
    {
        for( index_t k = 0; k < lines_.size(); k++ ) {
            set_edge_line_color( k, r, g, b ) ;
        }
    }
    void BoundaryModelGfx::set_edge_line_color( index_t l, float r, float g, float b )
    {
        lines_[l]->gfx().set_mesh_color( r, g, b ) ;
    }
    void BoundaryModelGfx::set_edge_lines_visibility( bool b )
    {
        for( index_t l = 0; l < lines_.size(); l++ ) {
            set_edge_line_visibility( l, b ) ;
        }
    }
    void BoundaryModelGfx::set_edge_line_visibility( index_t l, bool b )
    {
        lines_[l]->set_edges_visible( b ) ;
    }
    void BoundaryModelGfx::set_edge_lines_size( index_t s )
    {
        for( index_t l = 0; l < lines_.size(); l++ ) {
            set_edge_line_size( l, s ) ;
        }
    }
    void BoundaryModelGfx::set_edge_line_size( index_t l, index_t s )
    {
        lines_[l]->gfx().set_mesh_width( s ) ;
    }

    void BoundaryModelGfx::set_vertex_lines_color( float r, float g, float b )
    {
        for( index_t k = 0; k < lines_.size(); k++ ) {
            set_vertex_line_color( k, r, g, b ) ;
        }
    }
    void BoundaryModelGfx::set_vertex_line_color( index_t l, float r, float g, float b )
    {
        lines_[l]->gfx().set_points_color( r, g, b ) ;
    }
    void BoundaryModelGfx::set_vertex_lines_visibility( bool b )
    {
        for( index_t l = 0; l < lines_.size(); l++ ) {
            set_vertex_line_visibility( l, b ) ;
        }
    }
    void BoundaryModelGfx::set_vertex_line_visibility( index_t l, bool b )
    {
        lines_[l]->set_vertices_visible( b ) ;
    }
    void BoundaryModelGfx::set_vertex_lines_size( index_t s )
    {
        for( index_t l = 0; l < lines_.size(); l++ ) {
            set_vertex_line_size( l, s ) ;
        }
    }
    void BoundaryModelGfx::set_vertex_line_size( index_t l, index_t s )
    {
        lines_[l]->gfx().set_points_size( s ) ;
    }

    void BoundaryModelGfx::draw_surfaces()
    {
        for( index_t s = 0; s < surfaces_.size(); s++ ) {
            if( surfaces_[s]->get_vertices_visible() ) surfaces_[s]->gfx().draw_vertices() ;
            if( surfaces_[s]->get_surface_visible() ) surfaces_[s]->gfx().draw_surface() ;
        }
    }
    void BoundaryModelGfx::set_surfaces_color( float r, float g, float b )
    {
        for( index_t s = 0; s < surfaces_.size(); s++ ) {
            set_surface_color( s, r, g, b ) ;
        }
    }
    void BoundaryModelGfx::set_surface_color( index_t s, float r, float g, float b )
    {
        surfaces_[s]->gfx().set_surface_color( r, g, b ) ;
    }
    void BoundaryModelGfx::set_surfaces_visibility( bool b )
    {
        for( index_t s = 0; s < surfaces_.size(); s++ ) {
            set_surface_visibility( s, b ) ;
        }
    }
    void BoundaryModelGfx::set_surface_visibility( index_t s, bool b )
    {
        surfaces_[s]->set_surface_visible( b ) ;
    }
    void BoundaryModelGfx::set_mesh_surfaces_color( float r, float g, float b )
    {
        for( index_t s = 0; s < surfaces_.size(); s++ ) {
            set_mesh_surface_color( s, r, g, b ) ;
        }
    }
    void BoundaryModelGfx::set_mesh_surface_color( index_t s, float r, float g, float b )
    {
        surfaces_[s]->gfx().set_mesh_color( r, g, b ) ;
    }
    void BoundaryModelGfx::set_mesh_surfaces_visibility( bool b )
    {
        for( index_t s = 0; s < surfaces_.size(); s++ ) {
            set_mesh_surface_visibility( s, b ) ;
        }
    }
    void BoundaryModelGfx::set_mesh_surface_visibility( index_t s, bool b )
    {
        surfaces_[s]->gfx().set_show_mesh( b ) ;
    }
    void BoundaryModelGfx::set_mesh_surfaces_size( index_t size )
    {
        for( index_t s = 0; s < surfaces_.size(); s++ ) {
            set_mesh_surface_size( s, size ) ;
        }
    }
    void BoundaryModelGfx::set_mesh_surface_size( index_t s, index_t size )
    {
        surfaces_[s]->gfx().set_mesh_width( size ) ;
    }

    void BoundaryModelGfx::set_vertex_surfaces_color( float r, float g, float b )
    {
        for( index_t s = 0; s < surfaces_.size(); s++ ) {
            set_vertex_surface_color( s, r, g, b ) ;
        }
    }
    void BoundaryModelGfx::set_vertex_surface_color( index_t s, float r, float g, float b )
    {
        surfaces_[s]->gfx().set_points_color( r, g, b ) ;
    }
    void BoundaryModelGfx::set_vertex_surfaces_visibility( bool b )
    {
        for( index_t s = 0; s < surfaces_.size(); s++ ) {
            set_vertex_surface_visibility( s, b ) ;
        }
    }
    void BoundaryModelGfx::set_vertex_surface_visibility( index_t s, bool b )
    {
        surfaces_[s]->set_vertices_visible( b ) ;
    }
    void BoundaryModelGfx::set_vertex_surfaces_size( index_t size )
    {
        for( index_t s = 0; s < surfaces_.size(); s++ ) {
            set_vertex_surface_size( s, size ) ;
        }
    }
    void BoundaryModelGfx::set_vertex_surface_size( index_t s, index_t size )
    {
        surfaces_[s]->gfx().set_points_size( size ) ;
    }

} // namespace

#endif
