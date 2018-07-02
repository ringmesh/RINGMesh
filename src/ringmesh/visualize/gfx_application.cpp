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

#include <ringmesh/visualize/gfx_application.h>

#ifdef RINGMESH_WITH_GRAPHICS

#include <algorithm>

#include <geogram/basic/command_line.h>
#include <geogram/basic/file_system.h>

#include <geogram/mesh/mesh_io.h>

#include <geogram_gfx/basic/GL.h>

#include <geogram_gfx/glup_viewer/glup_viewer.h>

#include <geogram_gfx/third_party/quicktext/glQuickText.h>

#include <ringmesh/basic/command_line.h>

#include <ringmesh/geomodel/builder/geomodel_builder.h>
#include <ringmesh/geomodel/core/geomodel_entity.h>
#include <ringmesh/geomodel/core/geomodel_geological_entity.h>
#include <ringmesh/geomodel/core/geomodel_mesh_entity.h>

#include <ringmesh/io/io.h>

// RINGMesh colormaps
#include <ringmesh/visualize/colormaps/bgyr.xpm>
#include <ringmesh/visualize/colormaps/jet.xpm>
#include <ringmesh/visualize/colormaps/lgbt.xpm>
#include <ringmesh/visualize/colormaps/minmax_transparent.xpm>

namespace
{
    using namespace RINGMesh;

    typedef std::vector< std::vector< ImColor > > ColorTable;

    ImColor black( 0, 0, 0 );
    ImColor dark_grey( 128, 128, 128 );
    ImColor grey( 192, 192, 192 );
    ImColor white( 255, 255, 255 );

    ImColor violet( 71, 61, 139 );
    ImColor blue( 0, 0, 255 );
    ImColor other_blue( 100, 151, 237 );
    ImColor light_blue( 136, 207, 235 );

    ImColor grass_green( 85, 107, 47 );
    ImColor green( 50, 205, 50 );
    ImColor light_green( 175, 255, 47 );
    ImColor brown( 160, 81, 45 );

    ImColor red( 255, 0, 0 );
    ImColor orange( 255, 162, 0 );
    ImColor yellow( 255, 255, 0 );
    ImColor pink( 255, 0, 255 );

    std::string path_to_label(
        const std::string& viewer_path, const std::string& path )
    {
        if( GEO::String::string_starts_with( path, viewer_path ) )
        {
            return path.substr(
                viewer_path.length(), path.length() - viewer_path.length() );
        }
        return path;
    }

    bool GetChar( void* data, int idx, const char** out_text )
    {
        *out_text = static_cast< const std::vector< std::string >* >( data )
                        ->at( static_cast< long unsigned int >( idx ) )
                        .c_str();
        return true;
    }

    template < index_t DIMENSION >
    void compute_mesh_entity_bbox(
        const GeoModelMeshEntity< DIMENSION >& entity, Box< DIMENSION >& bbox )
    {
        for( auto v : range( entity.nb_vertices() ) )
        {
            bbox.add_point( entity.vertex( v ) );
        }
    }

    void home_rotation()
    {
        float axis[3] = { 1, 0.5, 1 };
        float angle{ 80 };
        glup_viewer_set_scene_rotation( axis, angle );
    }
}
namespace RINGMesh
{
    template < index_t DIMENSION >
    RINGMeshApplication::GeoModelViewerBase< DIMENSION >::GeoModelViewerBase(
        RINGMeshApplication& app, const std::string& filename )
        : app_( app )
    {
        corner_style_.color_ = red;
        corner_style_.size_ = 1;
        corner_style_.visible_vertices_ = false;
        corner_style_.vertex_color_ = pink;
        corner_style_.vertex_size_ = 0;

        line_style_.color_ = black;
        line_style_.size_ = 1;
        line_style_.visible_vertices_ = false;
        line_style_.vertex_color_ = orange;
        line_style_.vertex_size_ = 3;

        surface_style_.color_ = grey;
        surface_style_.size_ = 1;
        surface_style_.visible_vertices_ = false;
        surface_style_.vertex_color_ = light_blue;
        surface_style_.vertex_size_ = 3;

        mesh_color_ = black;
        reset_attribute_name();

        geomodel_load( GM_, filename );
        // Computation of the BBox is set with surface vertices
        // or with those of lines and corners if the model has no surface
        if( GM_.nb_surfaces() > 0 )
        {
            for( const auto& surface : GM_.surfaces() )
            {
                compute_mesh_entity_bbox( surface, bbox_ );
            }
        }
        else if( GM_.nb_lines() > 0 )
        {
            for( const auto& line : GM_.lines() )
            {
                compute_mesh_entity_bbox( line, bbox_ );
            }
        }
        else
        {
            for( const auto& corner : GM_.corners() )
            {
                compute_mesh_entity_bbox( corner, bbox_ );
            }
        }

        const std::vector< MeshEntityType >& types =
            GM_.entity_type_manager().mesh_entity_manager.mesh_entity_types();
        entity_types_.reserve( types.size() + 1 );
        entity_types_.emplace_back( "All" );
        for( const MeshEntityType& type : types )
        {
            entity_types_.emplace_back( type.string() );
        }
        for( auto i : range( GM_.nb_geological_entity_types() ) )
        {
            entity_types_.emplace_back(
                GM_.geological_entity_type( i ).string() );
        }
        GM_gfx_.set_geomodel( GM_ );
        if( !app.colormaps_.empty() )
        {
            GM_gfx_.attribute.set_colormap( app.colormaps_[0].texture );
        }
    }

    template < index_t DIMENSION >
    void RINGMeshApplication::GeoModelViewerBase<
        DIMENSION >::reset_attribute_name()
    {
        GM_gfx_.attribute.set_name( "name" );
    }

    template < index_t DIMENSION >
    void RINGMeshApplication::GeoModelViewerBase< DIMENSION >::draw_scene()
    {
        if( selected_entity_type_ != 0 )
        {
            index_t selected_entity_type_casted =
                static_cast< index_t >( selected_entity_type_ );
            const std::string& type =
                entity_types_[selected_entity_type_casted];

            if( selected_entity_type_casted
                < GM_.entity_type_manager()
                          .mesh_entity_manager.nb_mesh_entity_types()
                      + 1 )
            {
                MeshEntityType mesh_type{ type };
                selected_entity_id_ = std::min(
                    static_cast< int >( GM_.nb_mesh_entities( mesh_type ) - 1 ),
                    selected_entity_id_ );
                gmme_id entity_id(
                    mesh_type, static_cast< index_t >( selected_entity_id_ ) );
                toggle_mesh_entity_and_boundaries_visibility( entity_id );
            }
            else
            {
                GeologicalEntityType geol_type{ type };
                selected_entity_id_ =
                    std::min( static_cast< int >(
                                  GM_.nb_geological_entities( geol_type ) - 1 ),
                        selected_entity_id_ );
                gmge_id entity_id(
                    geol_type, static_cast< index_t >( selected_entity_id_ ) );
                toggle_geological_entity_visibility( entity_id );
            }
        }

        if( show_attributes_ )
        {
            GM_gfx_.attribute.bind_attribute();
        }
        else
        {
            GM_gfx_.attribute.unbind_attribute();
        }

        if( show_corners_ )
        {
            GM_gfx_.corners.set_vertex_color( corner_style_.color_.x,
                corner_style_.color_.y, corner_style_.color_.z );
            GM_gfx_.corners.set_vertex_size(
                static_cast< index_t >( corner_style_.size_ ) );
            GM_gfx_.corners.draw();
        }

        if( show_lines_ )
        {
            GM_gfx_.lines.set_line_color( line_style_.color_.x,
                line_style_.color_.y, line_style_.color_.z );
            GM_gfx_.lines.set_line_size(
                static_cast< index_t >( line_style_.size_ ) );
            if( selected_entity_type_ == 0 )
            {
                GM_gfx_.lines.set_vertex_visibility(
                    line_style_.visible_vertices_ );
            }
            if( line_style_.visible_vertices_ )
            {
                GM_gfx_.lines.set_vertex_size(
                    static_cast< index_t >( line_style_.vertex_size_ ) );
                GM_gfx_.lines.set_vertex_color( line_style_.vertex_color_.x,
                    line_style_.vertex_color_.y, line_style_.vertex_color_.z );
            }
            GM_gfx_.lines.draw();
        }

        if( show_surface_ )
        {
            GM_gfx_.surfaces.set_mesh_visibility( mesh_visible_ );
            GM_gfx_.surfaces.set_mesh_color(
                mesh_color_.x, mesh_color_.y, mesh_color_.z );
            GM_gfx_.surfaces.set_surface_color( surface_style_.color_.x,
                surface_style_.color_.y, surface_style_.color_.z );
            GM_gfx_.surfaces.set_mesh_size(
                static_cast< index_t >( surface_style_.size_ ) );
            if( selected_entity_type_ == 0 )
            {
                GM_gfx_.surfaces.set_vertex_visibility(
                    surface_style_.visible_vertices_ );
            }
            if( surface_style_.visible_vertices_ )
            {
                GM_gfx_.surfaces.set_vertex_size(
                    static_cast< index_t >( surface_style_.vertex_size_ ) );
                GM_gfx_.surfaces.set_vertex_color(
                    surface_style_.vertex_color_.x,
                    surface_style_.vertex_color_.y,
                    surface_style_.vertex_color_.z );
            }
            if( selected_entity_type_ == 0 )
            {
                for( const auto& surface : GM_.surfaces() )
                {
                    if( surface.is_on_voi() )
                    {
                        GM_gfx_.surfaces.set_surface_visibility(
                            surface.index(), show_voi_ );
                    }
                }
            }
            GM_gfx_.surfaces.draw();
        }
    }

    template < index_t DIMENSION >
    void RINGMeshApplication::GeoModelViewerBase< DIMENSION >::
        set_attribute_names( const std::vector< std::string >& names )
    {
        for( const std::string& name : names )
        {
            if( ImGui::Button( name.c_str() ) )
            {
                GM_gfx_.attribute.set_name( name );
                GM_gfx_.attribute.set_coordinate( 0 );
                autorange();
                ImGui::CloseCurrentPopup();
            }
        }
    }

    template < index_t DIMENSION >
    void RINGMeshApplication::GeoModelViewerBase< DIMENSION >::autorange()
    {
        GM_gfx_.attribute.compute_range();
        attribute_max_ = static_cast< float >( GM_gfx_.attribute.maximum() );
        attribute_min_ = static_cast< float >( GM_gfx_.attribute.minimum() );
    }

    template < index_t DIMENSION >
    void RINGMeshApplication::GeoModelViewerBase<
        DIMENSION >::update_all_entity_visibility( bool value )
    {
        GM_gfx_.corners.set_vertex_visibility( value );
        GM_gfx_.lines.set_line_visibility( value );
        GM_gfx_.surfaces.set_surface_visibility( value );
        if( !value || line_style_.visible_vertices_ )
        {
            GM_gfx_.lines.set_vertex_visibility( value );
        }
        if( !value || surface_style_.visible_vertices_ )
        {
            GM_gfx_.surfaces.set_vertex_visibility( value );
        }
    }

    template < index_t DIMENSION >
    void RINGMeshApplication::GeoModelViewerBase<
        DIMENSION >::update_entity_visibility()
    {
        index_t selected_entity_type_casted =
            static_cast< index_t >( selected_entity_type_ );
        const std::string& type = entity_types_[selected_entity_type_casted];
        if( selected_entity_type_ == 0 )
        {
            update_all_entity_visibility( true );
        }
        else
        {
            update_all_entity_visibility( false );
            if( selected_entity_type_casted
                < GM_.entity_type_manager()
                          .mesh_entity_manager.nb_mesh_entity_types()
                      + 1 )
            {
                MeshEntityType mesh_type{ type };
                selected_entity_id_ = std::min(
                    static_cast< int >( GM_.nb_mesh_entities( mesh_type ) - 1 ),
                    selected_entity_id_ );
                gmme_id entity_id(
                    mesh_type, static_cast< index_t >( selected_entity_id_ ) );
                toggle_mesh_entity_and_boundaries_visibility( entity_id );
            }
            else
            {
                GeologicalEntityType geol_type{ type };
                selected_entity_id_ =
                    std::min( static_cast< int >(
                                  GM_.nb_geological_entities( geol_type ) - 1 ),
                        selected_entity_id_ );
                gmge_id entity_id(
                    geol_type, static_cast< index_t >( selected_entity_id_ ) );
                toggle_geological_entity_visibility( entity_id );
            }
        }
    }

    template < index_t DIMENSION >
    void RINGMeshApplication::GeoModelViewerBase< DIMENSION >::
        toggle_mesh_entity_and_boundaries_visibility( const gmme_id& entity_id )
    {
        const MeshEntityTypeManager< DIMENSION >& manager =
            GM_.entity_type_manager().mesh_entity_manager;
        if( manager.is_corner( entity_id.type() ) )
        {
            toggle_corner_visibility( entity_id.index() );
        }
        else if( manager.is_line( entity_id.type() ) )
        {
            toggle_line_and_boundaries_visibility( entity_id.index() );
        }
        else if( manager.is_surface( entity_id.type() ) )
        {
            toggle_surface_and_boundaries_visibility( entity_id.index() );
        }
        else
        {
            ringmesh_assert_not_reached;
        }
    }

    template < index_t DIMENSION >
    void RINGMeshApplication::GeoModelViewerBase<
        DIMENSION >::toggle_corner_visibility( index_t corner_id )
    {
        GM_gfx_.corners.set_vertex_visibility( corner_id, true );
    }

    template < index_t DIMENSION >
    void RINGMeshApplication::GeoModelViewerBase<
        DIMENSION >::toggle_line_and_boundaries_visibility( index_t line_id )
    {
        GM_gfx_.lines.set_line_visibility( line_id, true );
        GM_gfx_.lines.set_vertex_visibility(
            line_id, line_style_.visible_vertices_ );
        const Line< DIMENSION >& line = GM_.line( line_id );
        for( auto i : range( line.nb_boundaries() ) )
        {
            toggle_corner_visibility( line.boundary_gmme( i ).index() );
        }
    }

    template < index_t DIMENSION >
    void RINGMeshApplication::GeoModelViewerBase< DIMENSION >::
        toggle_surface_and_boundaries_visibility( index_t surface_id )
    {
        GM_gfx_.surfaces.set_surface_visibility( surface_id, true );
        GM_gfx_.surfaces.set_vertex_visibility(
            surface_id, surface_style_.visible_vertices_ );
        const Surface< DIMENSION >& surface = GM_.surface( surface_id );
        for( auto i : range( surface.nb_boundaries() ) )
        {
            toggle_line_and_boundaries_visibility(
                surface.boundary_gmme( i ).index() );
        }
    }

    template < index_t DIMENSION >
    void RINGMeshApplication::GeoModelViewerBase< DIMENSION >::
        toggle_geological_entity_visibility( const gmge_id& entity_id )
    {
        const GeoModelGeologicalEntity< DIMENSION >& entity =
            GM_.geological_entity( entity_id );
        for( auto i : range( entity.nb_children() ) )
        {
            const gmme_id& child_id = entity.child_gmme( i );
            toggle_mesh_entity_and_boundaries_visibility( child_id );
        }
    }

    template < index_t DIMENSION >
    void RINGMeshApplication::GeoModelViewerBase<
        DIMENSION >::draw_object_properties()
    {
        if( ImGui::Combo( "Type", &selected_entity_type_, GetChar,
                static_cast< void* >( &entity_types_ ),
                static_cast< int >( entity_types_.size() ) ) )
        {
            update_entity_visibility();
        }
        if( selected_entity_type_ > 0 )
        {
            if( ImGui::InputInt( "Id", &selected_entity_id_, 1 ) )
            {
                selected_entity_id_ = std::max( 0, selected_entity_id_ );
                update_entity_visibility();
            }
        }
        ImGui::Separator();
        ImGui::Checkbox( "Attributes", &show_attributes_ );
        if( show_attributes_ )
        {
            if( ImGui::Button( GM_gfx_.attribute.location_name().c_str(),
                    ImVec2( -1, 0 ) ) )
            {
                ImGui::OpenPopup( "##Locations" );
            }
            if( ImGui::BeginPopup( "##Locations" ) )
            {
                std::vector< std::string > locations =
                    GM_gfx_.attribute.registered_locations();
                for( const std::string& location : locations )
                {
                    if( ImGui::Button( location.c_str() ) )
                    {
                        GM_gfx_.attribute.set_location( location );
                        reset_attribute_name();
                        ImGui::CloseCurrentPopup();
                    }
                }
                ImGui::EndPopup();
            }

            if( ImGui::Button(
                    GM_gfx_.attribute.name().c_str(), ImVec2( -1, 0 ) ) )
            {
                ImGui::OpenPopup( "##Attributes" );
            }
            if( ImGui::BeginPopup( "##Attributes" ) )
            {
                set_attribute_names( GM_gfx_.attribute.get_attribute_names() );
                ImGui::EndPopup();
            }
            if( GM_gfx_.attribute.location_name() != "location"
                && GM_gfx_.attribute.nb_coordinates() > 1 )
            {
                if( ImGui::Button(
                        std::to_string( GM_gfx_.attribute.coordinate() )
                            .c_str(),
                        ImVec2( -1, 0 ) ) )
                {
                    ImGui::OpenPopup( "##Coordinates" );
                }
                if( ImGui::BeginPopup( "##Coordinates" ) )
                {
                    for( auto i : range( GM_gfx_.attribute.nb_coordinates() ) )
                    {
                        if( ImGui::Button( std::to_string( i ).c_str() ) )
                        {
                            GM_gfx_.attribute.set_coordinate( i );
                            autorange();
                            ImGui::CloseCurrentPopup();
                        }
                    }
                    ImGui::EndPopup();
                }
            }
            if( ImGui::InputFloat( "min", &attribute_min_ ) )
            {
                GM_gfx_.attribute.set_minimum(
                    static_cast< double >( attribute_min_ ) );
            }
            if( ImGui::InputFloat( "max", &attribute_max_ ) )
            {
                GM_gfx_.attribute.set_maximum(
                    static_cast< double >( attribute_max_ ) );
            }
            if( ImGui::Button( "autorange", ImVec2( -1, 0 ) ) )
            {
                autorange();
            }
            if( ImGui::ImageButton(
                    app_.convert_to_ImTextureID( GM_gfx_.attribute.colormap() ),
                    ImVec2( 115, 8 ) ) )
            {
                ImGui::OpenPopup( "##Colormap" );
            }
            if( ImGui::BeginPopup( "##Colormap" ) )
            {
                for( const auto& colormap : app_.colormaps_ )
                {
                    if( ImGui::ImageButton(
                            app_.convert_to_ImTextureID( colormap.texture ),
                            ImVec2( 100, 8 ) ) )
                    {
                        GM_gfx_.attribute.set_colormap( colormap.texture );
                        ImGui::CloseCurrentPopup();
                    }
                }
                ImGui::EndPopup();
            }
            ImGui::Checkbox( "Colormap [M]", &show_colormap_ );
        }

        ImGui::Separator();
        ImGui::Checkbox( "VOI [V]", &show_voi_ );
        ImGui::Checkbox( "Mesh [m]", &mesh_visible_ );
        ImGui::SameLine();
        ImGui::ColorEdit3WithPalette( "Mesh color", &mesh_color_.x );

        ImGui::Separator();
        ImGui::Checkbox( "Corner [c]", &show_corners_ );
        draw_entity_style_editor( "Corner color", corner_style_ );

        ImGui::Separator();
        ImGui::Checkbox( "Line [e]", &show_lines_ );
        draw_entity_style_editor( "Line color", line_style_ );
        ImGui::Checkbox( "Vertices##Line", &line_style_.visible_vertices_ );
        if( line_style_.visible_vertices_ )
        {
            draw_entity_vertex_style_editor( "Line vertex color", line_style_ );
        }

        ImGui::Separator();
        ImGui::Checkbox( "Surface [s]", &show_surface_ );
        draw_entity_style_editor( "Surface color", surface_style_ );
        ImGui::Checkbox(
            "Vertices##Surface", &surface_style_.visible_vertices_ );
        if( surface_style_.visible_vertices_ )
        {
            draw_entity_vertex_style_editor(
                "Surface vertex color", surface_style_ );
        }
    }

    template < index_t DIMENSION >
    void RINGMeshApplication::GeoModelViewerBase<
        DIMENSION >::draw_entity_style_editor( const std::string& label,
        EntityStyle& style )
    {
        ImGui::SameLine();
        ImGui::ColorEdit3WithPalette( label.c_str(), &style.color_.x );
        ImGui::InputInt( "", &style.size_, 1 );
        style.size_ = std::max( style.size_, 0 );
    }

    template < index_t DIMENSION >
    void RINGMeshApplication::GeoModelViewerBase<
        DIMENSION >::draw_entity_vertex_style_editor( const std::string& label,
        EntityStyle& style )
    {
        ImGui::SameLine();
        ImGui::ColorEdit3WithPalette( label.c_str(), &style.vertex_color_.x );
        ImGui::InputInt( "", &style.vertex_size_, 1 );
        style.vertex_size_ = std::max( style.vertex_size_, 0 );
        style.vertex_size_ = std::min( style.vertex_size_, 50 );
    }

    template < index_t DIMENSION >
    void RINGMeshApplication::GeoModelViewerBase< DIMENSION >::draw_colormap()
    {
        GLUPboolean clipping_save = glupIsEnabled( GLUP_CLIPPING );
        glupDisable( GLUP_CLIPPING );

        glupMatrixMode( GLUP_TEXTURE_MATRIX );
        glupLoadIdentity();

        glupMatrixMode( GLUP_PROJECTION_MATRIX );
        glupPushMatrix();
        glupLoadIdentity();

        glupMatrixMode( GLUP_MODELVIEW_MATRIX );
        glupPushMatrix();
        glupLoadIdentity();

        const double z = -1.0;
        const double w = 0.3;
        const double h = 0.1;
        const double x1 = 0.;
        const double y1 = -0.9;
        const double tmin = GM_gfx_.attribute.minimum();
        const double tmax = GM_gfx_.attribute.maximum();
        GEO::glupMapTexCoords1d( tmin, tmax, 1. );

        glupColor3d( 1.0, 1.0, 1.0 );
        glupDisable( GLUP_LIGHTING );
        glupEnable( GLUP_TEXTURING );
        glupTextureMode( GLUP_TEXTURE_REPLACE );
        glupTextureType( GLUP_TEXTURE_1D );
        glupEnable( GLUP_DRAW_MESH );
        glupSetColor3d( GLUP_MESH_COLOR, 0.0, 0.0, 0.0 );
        glupSetMeshWidth( 2 );
        glupSetCellsShrink( 0.0f );

        glupBegin( GLUP_QUADS );
        glupTexCoord1d( tmin );
        glupVertex3d( x1 - w, y1, z );
        glupTexCoord1d( tmax );
        glupVertex3d( x1 + w, y1, z );
        glupTexCoord1d( tmax );
        glupVertex3d( x1 + w, y1 + h, z );
        glupTexCoord1d( tmin );
        glupVertex3d( x1 - w, y1 + h, z );
        glupEnd();

        glupTextureType( GLUP_TEXTURE_2D );
        glupMatrixMode( GLUP_TEXTURE_MATRIX );
        glupLoadIdentity();
        glupMatrixMode( GLUP_MODELVIEW_MATRIX );

        glupSetColor4d( GLUP_FRONT_AND_BACK_COLOR, 0.0, 0.0, 0.0, 1.0 );

        const double font_sz = 0.003;
        const double font_height = 0.4 * glQuickText::getFontHeight( font_sz );

        std::string min_value = std::to_string( GM_gfx_.attribute.minimum() );
        const double nb_min_letter = static_cast< double >( min_value.size() );
        glQuickText::printfAt( x1 - w - font_height * nb_min_letter * 0.3,
            y1 - font_height, z, font_sz, min_value.c_str() );

        std::string max_value = std::to_string( GM_gfx_.attribute.maximum() );
        const double nb_max_letter = static_cast< double >( max_value.size() );
        glQuickText::printfAt( x1 + w - font_height * nb_max_letter * 0.3,
            y1 - font_height, z, font_sz, max_value.c_str() );

        glupMatrixMode( GLUP_PROJECTION_MATRIX );
        glupPopMatrix();

        glupMatrixMode( GLUP_MODELVIEW_MATRIX );
        glupPopMatrix();

        if( clipping_save )
        {
            glupEnable( GLUP_CLIPPING );
        }
    }

    RINGMeshApplication::GeoModelViewer< 2 >::GeoModelViewer(
        RINGMeshApplication& app, const std::string& filename )
        : RINGMeshApplication::GeoModelViewerBase< 2 >( app, filename )
    {
    }

    RINGMeshApplication::GeoModelViewer< 3 >::GeoModelViewer(
        RINGMeshApplication& app, const std::string& filename )
        : RINGMeshApplication::GeoModelViewerBase< 3 >( app, filename )
    {
        home_rotation();
        volume_style_.color_ = grey;
        volume_style_.size_ = 1;
        volume_style_.visible_vertices_ = false;
        volume_style_.vertex_color_ = light_green;
        volume_style_.vertex_size_ = 3;

        if( GM_.nb_regions() > 0 )
        {
            meshed_regions_ = GM_.region( 0 ).is_meshed();
        }
        if( meshed_regions_ )
        {
            show_volume_ = true;
        }
    }

    void RINGMeshApplication::GeoModelViewer< 3 >::toggle_colored_cells()
    {
        show_colored_regions_.new_status = false;
        show_colored_layers_.new_status = false;
        GM_gfx_.regions.set_cell_colors_by_type();
    }

    void RINGMeshApplication::GeoModelViewer< 3 >::toggle_colored_regions()
    {
        colored_cells_.new_status = false;
        show_colored_layers_.new_status = false;
        for( auto r : range( GM_.nb_regions() ) )
        {
            GM_gfx_.regions.set_region_color( r,
                std::fmod( GEO::Numeric::random_float32(), 1.f ),
                std::fmod( GEO::Numeric::random_float32(), 1.f ),
                std::fmod( GEO::Numeric::random_float32(), 1.f ) );
        }
    }

    void RINGMeshApplication::GeoModelViewer< 3 >::toggle_colored_layers()
    {
        // To disable the key 'R'. If no layer within the model, layer is not
        // a valid type.
        if( !GM_.entity_type_manager().geological_entity_manager.is_valid_type(
                Layer3D::type_name_static() ) )
        {
            show_colored_layers_.new_status = false;
            return;
        }
        colored_cells_.new_status = false;
        show_colored_regions_.new_status = false;
        for( auto l :
            range( GM_.nb_geological_entities( Layer3D::type_name_static() ) ) )
        {
            float red = std::fmod( GEO::Numeric::random_float32(), 1.f );
            float green = std::fmod( GEO::Numeric::random_float32(), 1.f );
            float blue = std::fmod( GEO::Numeric::random_float32(), 1.f );
            const GeoModelGeologicalEntity3D& cur_layer =
                GM_.geological_entity( Layer3D::type_name_static(), l );
            for( auto r : range( cur_layer.nb_children() ) )
                GM_gfx_.regions.set_region_color(
                    cur_layer.child( r ).index(), red, green, blue );
        }
    }

    void RINGMeshApplication::GeoModelViewer< 3 >::draw_scene()
    {
        GeoModelViewerBase3D::draw_scene();

        if( show_volume_ && meshed_regions_ )
        {
            GM_gfx_.regions.set_mesh_visibility( mesh_visible_ );
            if( colored_cells_.need_to_update() )
            {
                colored_cells_.update();
                if( colored_cells_.new_status )
                {
                    toggle_colored_cells();
                }
            }
            else if( show_colored_regions_.need_to_update() )
            {
                show_colored_regions_.update();
                if( show_colored_regions_.new_status )
                {
                    toggle_colored_regions();
                }
            }
            else if( show_colored_layers_.need_to_update() )
            {
                show_colored_layers_.update();
                if( show_colored_layers_.new_status )
                {
                    toggle_colored_layers();
                }
            }
            if( !colored_cells_.new_status && !show_colored_regions_.new_status
                && !show_colored_layers_.new_status )
            {
                colored_cells_.update();
                show_colored_regions_.update();
                show_colored_layers_.update();
                GM_gfx_.regions.set_mesh_color(
                    mesh_color_.x, mesh_color_.y, mesh_color_.z );
                GM_gfx_.regions.set_region_color( volume_style_.color_.x,
                    volume_style_.color_.y, volume_style_.color_.z );
            }
            GM_gfx_.regions.set_mesh_size(
                static_cast< index_t >( volume_style_.size_ ) );
            if( selected_entity_type_ == 0 )
            {
                GM_gfx_.regions.set_vertex_visibility(
                    volume_style_.visible_vertices_ );
            }
            if( volume_style_.visible_vertices_ )
            {
                GM_gfx_.regions.set_vertex_size(
                    static_cast< index_t >( volume_style_.vertex_size_ ) );
                GM_gfx_.regions.set_vertex_color( volume_style_.vertex_color_.x,
                    volume_style_.vertex_color_.y,
                    volume_style_.vertex_color_.z );
            }
            GM_gfx_.regions.set_draw_cells( CellType::HEXAHEDRON, show_hex_ );
            GM_gfx_.regions.set_draw_cells( CellType::PRISM, show_prism_ );
            GM_gfx_.regions.set_draw_cells( CellType::PYRAMID, show_pyramid_ );
            GM_gfx_.regions.set_draw_cells(
                CellType::TETRAHEDRON, show_tetra_ );
            GM_gfx_.regions.set_shrink( static_cast< double >( shrink_ ) );
            GM_gfx_.regions.draw();
        }
    }

    void RINGMeshApplication::GeoModelViewer< 3 >::update_all_entity_visibility(
        bool value )
    {
        GeoModelViewerBase3D::update_all_entity_visibility( value );
        GM_gfx_.regions.set_region_visibility( value );
        if( volume_style_.visible_vertices_ )
        {
            GM_gfx_.regions.set_vertex_visibility( value );
        }
    }

    void RINGMeshApplication::GeoModelViewer<
        3 >::toggle_region_and_boundaries_visibility( index_t region_id )
    {
        GM_gfx_.regions.set_region_visibility( region_id, true );
        GM_gfx_.regions.set_vertex_visibility(
            region_id, volume_style_.visible_vertices_ );
        const Region3D& region = GM_.region( region_id );
        for( auto i : range( region.nb_boundaries() ) )
        {
            toggle_surface_and_boundaries_visibility(
                region.boundary_gmme( i ).index() );
        }
    }

    void RINGMeshApplication::GeoModelViewer< 3 >::
        toggle_mesh_entity_and_boundaries_visibility( const gmme_id& entity_id )
    {
        const MeshEntityTypeManager3D& manager =
            GM_.entity_type_manager().mesh_entity_manager;
        if( manager.is_region( entity_id.type() ) )
        {
            toggle_region_and_boundaries_visibility( entity_id.index() );
        }
        else
        {
            GeoModelViewerBase3D::toggle_mesh_entity_and_boundaries_visibility(
                entity_id );
        }
    }

    void RINGMeshApplication::GeoModelViewer< 3 >::draw_object_properties()
    {
        GeoModelViewerBase3D::draw_object_properties();

        if( meshed_regions_ )
        {
            ImGui::Separator();
            ImGui::Checkbox( "Region [v]", &show_volume_ );
            draw_entity_style_editor( "Volume color", volume_style_ );
            ImGui::Checkbox(
                "Vertices##Region", &volume_style_.visible_vertices_ );
            if( volume_style_.visible_vertices_ )
            {
                draw_entity_vertex_style_editor(
                    "Volume vertex color", volume_style_ );
            }
            if( show_volume_ )
            {
                ImGui::Checkbox( "Col. cells [C]", &colored_cells_.new_status );
                ImGui::Checkbox(
                    "Col. regions [r]", &show_colored_regions_.new_status );
                if( GM_.entity_type_manager()
                        .geological_entity_manager.is_valid_type(
                            Layer3D::type_name_static() ) )
                {
                    ImGui::Checkbox(
                        "Col. layers [R]", &show_colored_layers_.new_status );
                }
                ImGui::SliderFloat( "Shrk.", &shrink_, 0.0f, 1.0f, "%.1f" );
                ImGui::Checkbox( "Hex", &show_hex_ );
                ImGui::Checkbox( "Prism", &show_prism_ );
                ImGui::Checkbox( "Pyramid", &show_pyramid_ );
                ImGui::Checkbox( "Tetra", &show_tetra_ );
            }
        }
    }
    /*****************************************************************/

    RINGMeshApplication::MeshViewer::MeshViewer(
        RINGMeshApplication& app, const std::string& filename )
        : app_( app )
    {
        vertices_color_ = green;

        if( !filename.empty() )
        {
            GEO::mesh_load( filename, mesh_ );
            name_ = GEO::FileSystem::base_name( filename, true );
        }
        mesh_gfx_.set_mesh( &mesh_ );

        for( auto v : range( mesh_.vertices.nb() ) )
        {
            bbox_.add_point( mesh_.vertices.point( v ) );
        }
    }

    void RINGMeshApplication::MeshViewer::draw_object_properties()
    {
        ImGui::Checkbox( "attributes", &show_attributes_ );
        if( show_attributes_ )
        {
            if( attribute_min_ == 0.0f && attribute_max_ == 0.0f )
            {
                autorange();
            }
            if( ImGui::Button(
                    ( attribute_ + "##Attribute" ).c_str(), ImVec2( -1, 0 ) ) )
            {
                ImGui::OpenPopup( "##Attributes" );
            }
            if( ImGui::BeginPopup( "##Attributes" ) )
            {
                std::vector< std::string > attributes;
                GEO::String::split_string( attribute_names(), ';', attributes );
                for( const std::string& att : attributes )
                {
                    if( ImGui::Button( att.c_str() ) )
                    {
                        set_attribute( att );
                        ImGui::CloseCurrentPopup();
                    }
                }
                ImGui::EndPopup();
            }
            ImGui::InputFloat( "min", &attribute_min_ );
            ImGui::InputFloat( "max", &attribute_max_ );
            if( ImGui::Button( "autorange", ImVec2( -1, 0 ) ) )
            {
                autorange();
            }
            if( ImGui::ImageButton(
                    app_.convert_to_ImTextureID( current_colormap_texture_ ),
                    ImVec2( 115, 8 ) ) )
            {
                ImGui::OpenPopup( "##Colormap" );
            }
            if( ImGui::BeginPopup( "##Colormap" ) )
            {
                for( const auto& colormap : app_.colormaps_ )
                {
                    if( ImGui::ImageButton(
                            app_.convert_to_ImTextureID( colormap.texture ),
                            ImVec2( 100, 8 ) ) )
                    {
                        current_colormap_texture_ = colormap.texture;
                        ImGui::CloseCurrentPopup();
                    }
                }
                ImGui::EndPopup();
            }
        }

        ImGui::Separator();
        ImGui::Checkbox( "Vertices [p]", &show_vertices_ );
        if( show_vertices_ )
        {
            ImGui::SliderFloat( "size", &vertices_size_, 0.1f, 5.0f, "%.1f" );
            ImGui::ColorEdit3WithPalette( "Vertex color", &vertices_color_.x );
        }

        if( mesh_.edges.nb() != 0 )
        {
            ImGui::Separator();
            ImGui::Checkbox( "Edges", &show_edges_ );
            if( show_edges_ )
            {
                ImGui::InputInt( "", &edges_size_, 1 );
                edges_size_ = std::max( edges_size_, 0 );
                ImGui::ColorEdit3WithPalette( "Edge color", &edges_color_.x );
            }
        }

        if( mesh_.facets.nb() != 0 )
        {
            ImGui::Separator();
            ImGui::Checkbox( "Surface [S]", &show_surface_ );
            if( show_surface_ )
            {
                ImGui::Checkbox( "colors [c]", &show_surface_colors_ );
                ImGui::Checkbox( "mesh [m]", &show_mesh_ );
                ImGui::Checkbox( "borders [B]", &show_surface_borders_ );
            }
        }

        if( mesh_.cells.nb() != 0 )
        {
            ImGui::Separator();
            ImGui::Checkbox( "Volume [V]", &show_volume_ );
            if( show_volume_ )
            {
                ImGui::SliderFloat(
                    "shrk.", &cells_shrink_, 0.0f, 1.0f, "%.2f" );
                if( !mesh_.cells.are_simplices() )
                {
                    ImGui::Checkbox(
                        "colored cells [C]", &show_colored_cells_ );
                    ImGui::Checkbox( "hexes [j]", &show_hexes_ );
                }
            }
        }
    }

    void RINGMeshApplication::MeshViewer::draw_scene()
    {
        mesh_gfx_.set_lighting( app_.lighting_ );

        if( show_attributes_ )
        {
            mesh_gfx_.set_scalar_attribute( attribute_subelements_,
                attribute_name_, double( attribute_min_ ),
                double( attribute_max_ ), current_colormap_texture_, 1 );
        }
        else
        {
            mesh_gfx_.unset_scalar_attribute();
        }

        if( show_vertices_ )
        {
            mesh_gfx_.set_points_size( vertices_size_ );
            mesh_gfx_.set_points_color(
                vertices_color_.x, vertices_color_.y, vertices_color_.z );
            mesh_gfx_.draw_vertices();
        }

        if( show_edges_ )
        {
            mesh_gfx_.set_mesh_width( static_cast< index_t >( edges_size_ ) );
            mesh_gfx_.set_mesh_color(
                edges_color_.x, edges_color_.y, edges_color_.z );
            mesh_gfx_.draw_edges();
        }

        mesh_gfx_.set_mesh_color( 0.0, 0.0, 0.0 );

        if( show_surface_colors_ )
        {
            if( mesh_.cells.nb() == 0 )
            {
                mesh_gfx_.set_surface_color( 0.5f, 0.75f, 1.0f );
                mesh_gfx_.set_backface_surface_color( 1.0f, 0.0f, 0.0f );
            }
            else
            {
                mesh_gfx_.set_surface_color( 0.7f, 0.0f, 0.0f );
                mesh_gfx_.set_backface_surface_color( 1.0f, 1.0f, 0.0f );
            }
        }
        else
        {
            mesh_gfx_.set_surface_color( 0.9f, 0.9f, 0.9f );
        }

        mesh_gfx_.set_show_mesh( show_mesh_ );

        if( show_surface_ )
        {
            mesh_gfx_.draw_surface();
        }

        if( show_surface_borders_ )
        {
            mesh_gfx_.draw_surface_borders();
        }

        if( show_volume_ )
        {
            if( glupIsEnabled( GLUP_CLIPPING )
                && glupGetClipMode() == GLUP_CLIP_SLICE_CELLS )
            {
                mesh_gfx_.set_lighting( false );
            }

            mesh_gfx_.set_shrink( double( cells_shrink_ ) );
            mesh_gfx_.set_draw_cells( GEO::MESH_HEX, show_hexes_ );
            if( show_colored_cells_ )
            {
                mesh_gfx_.set_cells_colors_by_type();
            }
            else
            {
                mesh_gfx_.set_cells_color( 0.9f, 0.9f, 0.9f );
            }
            mesh_gfx_.draw_volume();

            mesh_gfx_.set_lighting( app_.lighting_ );
        }
    }

    void RINGMeshApplication::MeshViewer::autorange()
    {
        if( attribute_subelements_ != GEO::MESH_NONE )
        {
            attribute_min_ = 0.0;
            attribute_max_ = 0.0;
            const GEO::MeshSubElementsStore& subelements =
                mesh_.get_subelements_by_type( attribute_subelements_ );
            GEO::ReadOnlyScalarAttributeAdapter attribute(
                subelements.attributes(), attribute_name_ );
            if( attribute.is_bound() )
            {
                attribute_min_ = GEO::Numeric::max_float32();
                attribute_max_ = GEO::Numeric::min_float32();
                for( auto i : range( subelements.nb() ) )
                {
                    attribute_min_ =
                        GEO::geo_min( attribute_min_, float( attribute[i] ) );
                    attribute_max_ =
                        GEO::geo_max( attribute_max_, float( attribute[i] ) );
                }
            }
        }
    }

    std::string RINGMeshApplication::MeshViewer::attribute_names()
    {
        return mesh_.get_scalar_attributes();
    }

    void RINGMeshApplication::MeshViewer::set_attribute(
        const std::string& attribute )
    {
        attribute_ = attribute;
        std::string subelements_name;
        GEO::String::split_string(
            attribute_, '.', subelements_name, attribute_name_ );
        attribute_subelements_ =
            mesh_.name_to_subelements_type( subelements_name );
        if( attribute_min_ == 0.0f && attribute_max_ == 0.0f )
        {
            autorange();
        }
    }

    /*****************************************************************/

    RINGMeshApplication::RINGMeshApplication( int argc, char** argv )
        : GEO::Application( argc, argv, "<filename>" )
    {
        GEO::CmdLine::declare_arg( "attributes", true, "load mesh attributes" );
        GEO::CmdLine::declare_arg(
            "single_precision", false, "use single precision vertices (FP32)" );

        auto ringmesh_2d_extensions =
            GeoModelInputHandlerFactory2D::list_creators();
        auto ringmesh_3d_extensions =
            GeoModelInputHandlerFactory3D::list_creators();
        ringmesh_file_extensions_ =
            GEO::String::join_strings( ringmesh_2d_extensions, ';' )
            + GEO::String::join_strings( ringmesh_3d_extensions, ';' );

        std::vector< std::string > geogram_extensions;
        GEO::MeshIOHandlerFactory::list_creators( geogram_extensions );
        geogram_file_extensions_ =
            GEO::String::join_strings( geogram_extensions, ';' );

        print_header_information();

        Logger::div( "RINGMesh-View" );
        Logger::out( "", "Welcome to RINGMesh-View !" );
    }

    void RINGMeshApplication::quit()
    {
        glup_viewer_exit_main_loop();
        Logger::instance()->unregister_client( console_ );
    }

    RINGMeshApplication* RINGMeshApplication::instance()
    {
        RINGMeshApplication* result = dynamic_cast< RINGMeshApplication* >(
            GEO::Application::instance() );
        ringmesh_assert( result != nullptr );
        return result;
    }

    void RINGMeshApplication::browse_geogram( const std::string& path )
    {
        std::vector< std::string > files;
        GEO::FileSystem::get_directory_entries( path, files );
        std::sort( files.begin(), files.end() );
        for( const std::string& file : files )
        {
            if( GEO::FileSystem::is_directory( file ) )
            {
                if( ImGui::BeginMenu( path_to_label( path_, file ).c_str() ) )
                {
                    browse_geogram( file );
                    ImGui::EndMenu();
                }
            }
            else
            {
                if( can_load_geogram( file ) )
                {
                    if( ImGui::MenuItem(
                            path_to_label( path_, file ).c_str() ) )
                    {
                        load_geogram( file );
                    }
                }
            }
        }
    }

    bool RINGMeshApplication::can_load_geogram( const std::string& filename )
    {
        std::string extensions_str = supported_geogram_read_file_extensions();
        if( extensions_str == "" )
        {
            return false;
        }
        if( extensions_str == "*" )
        {
            return true;
        }
        std::string extension = GEO::FileSystem::extension( filename );
        std::vector< std::string > extensions;
        GEO::String::split_string( extensions_str, ';', extensions );
        for( const std::string& ext : extensions )
        {
            if( ext == extension )
            {
                return true;
            }
        }
        return false;
    }

    bool RINGMeshApplication::load_geogram( const std::string& filename )
    {
        if( !filename.empty() )
        {
            meshes_.emplace_back( new MeshViewer( *this, filename ) );
            current_viewer_ = static_cast< index_t >( meshes_.size() - 1 );
            current_viewer_type_ = ViewerType::MESH;
        }

        update_region_of_interest();
        return true;
    }

    void RINGMeshApplication::draw_application_menus()
    {
        if( ImGui::BeginMenu( "Debug" ) )
        {
            if( ImGui::BeginMenu( "Load..." ) )
            {
                ImGui::Selectable( ".." );
                if( ImGui::IsItemClicked() )
                {
                    path_ += "/..";
                }
                browse_geogram( path_ );
                ImGui::EndMenu();
            }
            ImGui::EndMenu();
        }
        if( ImGui::BeginMenu( "Create..." ) )
        {
            if( ImGui::MenuItem( "point" ) )
            {
                GEO::Command::set_current(
                    "create_point(std::string name=\"debug\","
                    " double x=0, double y=0, double z=0)",
                    this, &RINGMeshApplication::create_point );
            }
            if( ImGui::MenuItem( "AABB" ) )
            {
                GEO::Command::set_current(
                    "create_aabbox(std::string name=\"box\","
                    " double xmin=0, double ymin=0, double zmin=0,"
                    " double xmax=1, double ymax=1, double zmax=1)",
                    this, &RINGMeshApplication::create_aabbox );
            }
            ImGui::EndMenu();
        }
    }

    void RINGMeshApplication::create_point(
        std::string name, double x, double y, double z )
    {
        MeshViewer* viewer{ nullptr };
        for( auto& i : meshes_ )
        {
            if( i->name_ == name )
            {
                viewer = i.get();
                break;
            }
        }
        if( !viewer )
        {
            meshes_.emplace_back( new MeshViewer( *this, "" ) );
            viewer = meshes_.back().get();
        }
        vec3 point{ x, y, z };
        viewer->mesh_.vertices.create_vertex( point.data() );
        viewer->mesh_gfx_.set_mesh( &viewer->mesh_ );
        viewer->bbox_.add_point( point );
        viewer->name_ = name;
        viewer->show_vertices_ = true;
        current_viewer_ = static_cast< index_t >( meshes_.size() - 1 );
        current_viewer_type_ = ViewerType::MESH;
        update_region_of_interest();
    }

    void RINGMeshApplication::create_aabbox( std::string name,
        double xmin,
        double ymin,
        double zmin,
        double xmax,
        double ymax,
        double zmax )
    {
        vec3 min{ xmin, ymin, zmin };
        vec3 max{ xmax, ymax, zmax };
        MeshViewer* viewer{ nullptr };
        for( auto& i : meshes_ )
        {
            if( i->name_ == name )
            {
                viewer = i.get();
                break;
            }
        }
        if( !viewer )
        {
            meshes_.emplace_back( new MeshViewer( *this, "" ) );
            viewer = meshes_.back().get();
        }
        const index_t prev_nbv{ viewer->mesh_.vertices.nb() };
        vec3 box_other_vertex1{ min[0], min[1], max[2] };
        vec3 box_other_vertex2{ min[0], max[1], max[2] };
        vec3 box_other_vertex3{ max[0], min[1], max[2] };
        vec3 box_other_vertex4{ max[0], max[1], min[2] };
        vec3 box_other_vertex5{ max[0], min[1], min[2] };
        vec3 box_other_vertex6{ min[0], max[1], min[2] };
        viewer->mesh_.vertices.create_vertex( min.data() );
        viewer->mesh_.vertices.create_vertex( box_other_vertex1.data() );
        viewer->mesh_.vertices.create_vertex( box_other_vertex2.data() );
        viewer->mesh_.vertices.create_vertex( box_other_vertex3.data() );
        viewer->mesh_.vertices.create_vertex( box_other_vertex4.data() );
        viewer->mesh_.vertices.create_vertex( box_other_vertex5.data() );
        viewer->mesh_.vertices.create_vertex( box_other_vertex6.data() );
        viewer->mesh_.vertices.create_vertex( max.data() );
        viewer->mesh_.edges.create_edge( prev_nbv + 0, prev_nbv + 1 );
        viewer->mesh_.edges.create_edge( prev_nbv + 0, prev_nbv + 5 );
        viewer->mesh_.edges.create_edge( prev_nbv + 0, prev_nbv + 6 );
        viewer->mesh_.edges.create_edge( prev_nbv + 1, prev_nbv + 2 );
        viewer->mesh_.edges.create_edge( prev_nbv + 1, prev_nbv + 3 );
        viewer->mesh_.edges.create_edge( prev_nbv + 2, prev_nbv + 6 );
        viewer->mesh_.edges.create_edge( prev_nbv + 2, prev_nbv + 7 );
        viewer->mesh_.edges.create_edge( prev_nbv + 3, prev_nbv + 5 );
        viewer->mesh_.edges.create_edge( prev_nbv + 3, prev_nbv + 7 );
        viewer->mesh_.edges.create_edge( prev_nbv + 4, prev_nbv + 5 );
        viewer->mesh_.edges.create_edge( prev_nbv + 4, prev_nbv + 6 );
        viewer->mesh_.edges.create_edge( prev_nbv + 4, prev_nbv + 7 );
        viewer->mesh_.facets.create_quad(
            prev_nbv + 0, prev_nbv + 6, prev_nbv + 4, prev_nbv + 5 );
        viewer->mesh_.facets.create_quad(
            prev_nbv + 0, prev_nbv + 1, prev_nbv + 2, prev_nbv + 6 );
        viewer->mesh_.facets.create_quad(
            prev_nbv + 0, prev_nbv + 5, prev_nbv + 3, prev_nbv + 1 );
        viewer->mesh_.facets.create_quad(
            prev_nbv + 7, prev_nbv + 2, prev_nbv + 1, prev_nbv + 3 );
        viewer->mesh_.facets.create_quad(
            prev_nbv + 7, prev_nbv + 3, prev_nbv + 5, prev_nbv + 4 );
        viewer->mesh_.facets.create_quad(
            prev_nbv + 7, prev_nbv + 4, prev_nbv + 6, prev_nbv + 2 );
        viewer->mesh_gfx_.set_mesh( &viewer->mesh_ );
        viewer->bbox_.add_point( min );
        viewer->bbox_.add_point( max );
        viewer->name_ = name;
        viewer->show_vertices_ = false;
        viewer->show_mesh_ = true;
        viewer->show_surface_borders_ = false;
        viewer->show_surface_colors_ = false;
        current_viewer_ = static_cast< index_t >( meshes_.size() - 1 );
        current_viewer_type_ = ViewerType::MESH;
        update_region_of_interest();
    }

    void RINGMeshApplication::init_graphics()
    {
        GEO::Application::init_graphics();

        init_ringmesh_colormaps();
        glupEnable( GLUP_ALPHA_DISCARD );
        glup_viewer_disable( GLUP_VIEWER_BACKGROUND );
    }

    bool RINGMeshApplication::load( const std::string& filename )
    {
        if( !filename.empty() && GEO::FileSystem::is_file( filename ) )
        {
            index_t dimension = find_geomodel_dimension( filename );
            if( dimension == 2 )
            {
                geomodels2d_.emplace_back(
                    new GeoModelViewer2D( *this, filename ) );
                current_viewer_ =
                    static_cast< index_t >( geomodels2d_.size() - 1 );
                current_viewer_type_ = ViewerType::GEOMODEL2D;
            }
            else if( dimension == 3 )
            {
                geomodels3d_.emplace_back(
                    new GeoModelViewer3D( *this, filename ) );
                current_viewer_ =
                    static_cast< index_t >( geomodels3d_.size() - 1 );
                current_viewer_type_ = ViewerType::GEOMODEL3D;
            }
            else
            {
                ringmesh_assert_not_reached;
            }
        }

        update_region_of_interest();
        return true;
    }

    bool RINGMeshApplication::save( const std::string& filename )
    {
        if( !filename.empty() )
        {
            if( current_viewer_type_ == ViewerType::GEOMODEL2D )
            {
                GeoModelViewerBase2D& viewer = *geomodels2d_[current_viewer_];
                geomodel_save( viewer.GM_, filename );
            }
            if( current_viewer_type_ == ViewerType::GEOMODEL3D )
            {
                GeoModelViewerBase3D& viewer = *geomodels3d_[current_viewer_];
                geomodel_save( viewer.GM_, filename );
            }
        }
        return true;
    }

    void RINGMeshApplication::update_region_of_interest()
    {
        Box3D bbox;
        for( std::unique_ptr< GeoModelViewer2D >& geomodel : geomodels2d_ )
        {
            if( geomodel->is_visible_ )
            {
                vec2 min = geomodel->bbox_.min();
                vec2 max = geomodel->bbox_.max();
                bbox.add_point( vec3( min.x, min.y, 0. ) );
                bbox.add_point( vec3( max.x, max.y, 0. ) );
            }
        }
        for( std::unique_ptr< GeoModelViewer3D >& geomodel : geomodels3d_ )
        {
            if( geomodel->is_visible_ )
            {
                bbox.add_box( geomodel->bbox_ );
            }
        }
        for( std::unique_ptr< MeshViewer >& mesh : meshes_ )
        {
            if( mesh->is_visible_ )
            {
                bbox.add_box( mesh->bbox_ );
            }
        }

        if( bbox.initialized() )
        {
            glup_viewer_set_region_of_interest( float( bbox.min()[0] ),
                float( bbox.min()[1] ), float( bbox.min()[2] ),
                float( bbox.max()[0] ), float( bbox.max()[1] ),
                float( bbox.max()[2] ) );
        }
    }

    void RINGMeshApplication::init_ringmesh_colormaps()
    {
        init_colormaps();
        init_colormap( "bgyr", bgyr_xpm );
        init_colormap( "jet", jet_xpm );
        init_colormap( "lgbt", lgbt_xpm );
        init_colormap( "minmax_transparent", minmax_transparent_xpm );
    }

    void RINGMeshApplication::draw_scene()
    {
        if( current_viewer_ == NO_ID )
            return;

        for( std::unique_ptr< MeshViewer >& mesh : meshes_ )
        {
            if( mesh->is_visible_ )
            {
                mesh->draw_scene();
            }
        }
        for( std::unique_ptr< GeoModelViewer2D >& geomodel : geomodels2d_ )
        {
            if( geomodel->is_visible_ )
            {
                geomodel->draw_scene();
            }
        }
        for( std::unique_ptr< GeoModelViewer3D >& geomodel : geomodels3d_ )
        {
            if( geomodel->is_visible_ )
            {
                geomodel->draw_scene();
            }
        }

        if( current_viewer_type_ == ViewerType::GEOMODEL2D )
        {
            GeoModelViewerBase2D& viewer = *geomodels2d_[current_viewer_];
            if( viewer.show_colormap_ )
            {
                viewer.draw_colormap();
            }
        }
        if( current_viewer_type_ == ViewerType::GEOMODEL3D )
        {
            GeoModelViewerBase3D& viewer = *geomodels3d_[current_viewer_];
            if( viewer.show_colormap_ )
            {
                viewer.draw_colormap();
            }
        }
    }

    std::string RINGMeshApplication::supported_read_file_extensions()
    {
        return ringmesh_file_extensions_;
    }
    std::string RINGMeshApplication::supported_geogram_read_file_extensions()
    {
        return geogram_file_extensions_;
    }

    std::string RINGMeshApplication::supported_write_file_extensions()
    {
        return ringmesh_file_extensions_;
    }

    template < index_t DIMENSION >
    void RINGMeshApplication::draw_geomodel_viewer_properties(
        std::vector< std::unique_ptr< GeoModelViewer< DIMENSION > > >&
            geomodels,
        int& id )
    {
        if( !geomodels.empty() )
        {
            ImGui::Separator();
            std::ostringstream oss;
            oss << "GeoModel" << DIMENSION << "D";
            ImGui::Text( "%s", oss.str().c_str() );
            for( auto i : range( geomodels.size() ) )
            {
                GeoModelViewer< DIMENSION >& viewer = *geomodels[i];
                ImGui::PushID( id++ );
                if( ImGui::Checkbox(
                        viewer.GM_.name().c_str(), &viewer.is_visible_ ) )
                {
                    current_viewer_ = i;
                    current_viewer_type_ = viewer.type();
                    update_region_of_interest();
                }
                ImGui::SameLine( ImGui::GetWindowWidth() - 30 );
                if( ImGui::Button( "X" ) )
                {
                    geomodels.erase( geomodels.begin() + i );
                    if( current_viewer_type_ == viewer.type()
                        && current_viewer_ >= i )
                    {
                        current_viewer_--;
                    }
                    if( geomodels.empty() )
                    {
                        current_viewer_type_ = ViewerType::NONE;
                    }
                    break;
                }
                ImGui::PopID();
            }
        }
    }

    void RINGMeshApplication::draw_viewer_properties()
    {
        if( ImGui::Button( "home [H]", ImVec2( -1, 0 ) ) )
        {
            glup_viewer_home();
            bool geomodel_2d_displayed{ false };
            for( auto& geomodel2d : geomodels2d_ )
            {
                if( geomodel2d->is_visible_ )
                {
                    geomodel_2d_displayed = true;
                    break;
                }
            }
            if( !geomodel_2d_displayed )
            {
                home_rotation();
            }
        }

        ImGui::Separator();
        ImGui::Checkbox( "Lighting [L]", &lighting_ );
        if( lighting_ )
        {
            ImGui::Checkbox(
                "edit light [l]", (bool*) glup_viewer_is_enabled_ptr(
                                      GLUP_VIEWER_ROTATE_LIGHT ) );
        }

        ImGui::Separator();
        ImGui::Checkbox( "Clipping [F1]",
            (bool*) glup_viewer_is_enabled_ptr( GLUP_VIEWER_CLIP ) );
        if( glup_viewer_is_enabled( GLUP_VIEWER_CLIP ) )
        {
            ImGui::Combo( "mode", (int*) &clip_mode_,
                "std. GL\0cells\0straddle\0slice\0\0" );
            ImGui::Checkbox( "edit clip [F2]",
                (bool*) glup_viewer_is_enabled_ptr( GLUP_VIEWER_EDIT_CLIP ) );
            ImGui::Checkbox( "fixed clip [F3]",
                (bool*) glup_viewer_is_enabled_ptr( GLUP_VIEWER_FIXED_CLIP ) );
        }

        ImGui::Separator();
        ImGui::Text( "Background" );
        if( ImGui::ColorEdit3WithPalette(
                "Color", background_color_1_.data() ) )
        {
            glup_viewer_enable( GLUP_VIEWER_BACKGROUND );
            background_color_2_ = background_color_1_;
        }

        int id{ 0 };
        draw_geomodel_viewer_properties( geomodels2d_, id );
        draw_geomodel_viewer_properties( geomodels3d_, id );

        if( !meshes_.empty() )
        {
            ImGui::Separator();
            ImGui::Text( "Mesh" );
            for( auto i : range( meshes_.size() ) )
            {
                MeshViewer& viewer = *meshes_[i];
                ImGui::PushID( id++ );
                if( ImGui::Checkbox(
                        viewer.name_.c_str(), &viewer.is_visible_ ) )
                {
                    current_viewer_ = i;
                    current_viewer_type_ = ViewerType::MESH;
                    update_region_of_interest();
                }
                ImGui::SameLine( ImGui::GetWindowWidth() - 30 );
                if( ImGui::Button( "X" ) )
                {
                    meshes_.erase( meshes_.begin() + i );
                    if( current_viewer_type_ == ViewerType::MESH
                        && current_viewer_ >= i )
                    {
                        current_viewer_--;
                    }
                    if( meshes_.empty() )
                    {
                        current_viewer_type_ = ViewerType::NONE;
                    }
                    break;
                }
                ImGui::PopID();
            }
        }
    }

    void RINGMeshApplication::draw_object_properties()
    {
        if( current_viewer_ == NO_ID )
            return;
        switch( current_viewer_type_ )
        {
        case ViewerType::GEOMODEL2D:
            ringmesh_assert( current_viewer_ < geomodels2d_.size() );
            geomodels2d_[current_viewer_]->draw_object_properties();
            return;
        case ViewerType::GEOMODEL3D:
            ringmesh_assert( current_viewer_ < geomodels3d_.size() );
            geomodels3d_[current_viewer_]->draw_object_properties();
            return;
        case ViewerType::MESH:
            ringmesh_assert( current_viewer_ < meshes_.size() );
            meshes_[current_viewer_]->draw_object_properties();
            return;
        default:
            return;
        }
    }
} // namespace RINGMesh

#endif
