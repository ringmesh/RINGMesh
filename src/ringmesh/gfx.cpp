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

    BoundaryModelMeshElementGfx::BoundaryModelMeshElementGfx(
        const BoundaryModelMeshElement& mesh )
        : vertices_visible_( false )
    {
        gfx_.set_mesh( &mesh.mesh() ) ;
    }

    CornerGfx::CornerGfx( const Corner& corner )
        : BoundaryModelMeshElementGfx( corner )
    {
        vertices_visible_ = true ;
    }

    LineGfx::LineGfx( const Line& line )
        : BoundaryModelMeshElementGfx( line )
    {
        vertices_visible_ = false ;
    }

    SurfaceGfx::SurfaceGfx( const Surface& surface )
        : BoundaryModelMeshElementGfx( surface )
    {
        vertices_visible_ = false ;
    }

    BoundaryModelGfx::BoundaryModelGfx()
        : model_( nil )
    {
    }

    void BoundaryModelGfx::set_boundary_model( const BoundaryModel& model )
    {
        model_ = &model ;
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

    void BoundaryModelGfx::draw_corners()
    {
        for( index_t c = 0; c < corners_.size(); c++ ) {
            if( corners_[c]->get_vertices_visible() ) corners_[c]->gfx().draw_vertices() ;
        }
    }

} // namespace

#endif
