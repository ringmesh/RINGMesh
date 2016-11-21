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

#include <ringmesh/visualization/gfx_application.h>

#ifdef RINGMESH_WITH_GRAPHICS

#include <geogram/basic/file_system.h>
#include <geogram/basic/logger.h>
#include <geogram/mesh/mesh_geometry.h>

#include <geogram_gfx/basic/GL.h>
#include <geogram_gfx/glup_viewer/glup_viewer.h>
#include <geogram_gfx/third_party/quicktext/glQuickText.h>

#include <ringmesh/basic/command_line.h>
#include <ringmesh/geomodel/geo_model_entity.h>
#include <ringmesh/geomodel/geo_model_mesh_entity.h>
#include <ringmesh/geomodel/geo_model_geological_entity.h>
#include <ringmesh/io/io.h>

namespace {
    using namespace RINGMesh ;

    typedef std::vector< std::vector< ImColor > > ColorTable ;

    ImColor black( 0, 0, 0 ) ;
    ImColor dark_grey( 128, 128, 128 ) ;
    ImColor grey( 192, 192, 192 ) ;
    ImColor white( 255, 255, 255 ) ;

    ImColor violet( 71, 61, 139 ) ;
    ImColor blue( 0, 0, 255 ) ;
    ImColor other_blue( 100, 151, 237 ) ;
    ImColor light_blue( 136, 207, 235 ) ;

    ImColor grass_green( 85, 107, 47 ) ;
    ImColor green( 50, 205, 50 ) ;
    ImColor light_green( 175, 255, 47 ) ;
    ImColor brown( 160, 81, 45 ) ;

    ImColor red( 255, 0, 0 ) ;
    ImColor orange( 255, 162, 0 ) ;
    ImColor yellow( 255, 255, 0 ) ;
    ImColor pink( 255, 0, 255 ) ;

    ColorTable create_color_table()
    {
        ColorTable color_table_init( 4 ) ;
        color_table_init[0].push_back( black ) ;
        color_table_init[0].push_back( dark_grey ) ;
        color_table_init[0].push_back( grey ) ;
        color_table_init[0].push_back( white ) ;

        color_table_init[1].push_back( violet ) ;
        color_table_init[1].push_back( blue ) ;
        color_table_init[1].push_back( other_blue ) ;
        color_table_init[1].push_back( light_blue ) ;

        color_table_init[2].push_back( grass_green ) ;
        color_table_init[2].push_back( green ) ;
        color_table_init[2].push_back( light_green ) ;
        color_table_init[2].push_back( brown ) ;

        color_table_init[3].push_back( red ) ;
        color_table_init[3].push_back( orange ) ;
        color_table_init[3].push_back( yellow ) ;
        color_table_init[3].push_back( pink ) ;
        return color_table_init ;
    }

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

    bool GetChar( void* data, int idx, const char** out_text )
    {
        *out_text = static_cast< const std::vector< std::string >* >( data )->at(
            static_cast< long unsigned int >( idx ) ).c_str() ;
        return true ;
    }

}
namespace RINGMesh {

    ColorTable RINGMeshApplication::GeoModelViewer::color_table_ =
        create_color_table() ;

    RINGMeshApplication::GeoModelViewer::GeoModelViewer(
        RINGMeshApplication& app,
        const std::string& filename )
        : app_( app )
    {
        is_visible_ = true ;

        show_corners_ = true ;
        corner_style_.color_ = red ;
        corner_style_.size_ = 1 ;
        show_lines_ = true ;
        line_style_.color_ = black ;
        line_style_.size_ = 1 ;
        show_surface_ = true ;
        surface_style_.color_ = grey ;
        surface_style_.size_ = 1 ;
        show_volume_ = false ;
        volume_style_.color_ = grey ;
        volume_style_.size_ = 1 ;
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
        mesh_color_ = black ;
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
        selected_entity_type_ = 0 ;
        selected_entity_id_ = 0 ;
        entity_types_.push_back( "All" ) ;
        entity_types_.push_back( Corner::type_name_static() ) ;
        entity_types_.push_back( Line::type_name_static() ) ;
        entity_types_.push_back( Surface::type_name_static() ) ;
        entity_types_.push_back( Region::type_name_static() ) ;
        for( index_t i = 0; i < GM_.nb_geological_entity_types(); i++ ) {
            entity_types_.push_back( GM_.geological_entity_type( i ) ) ;
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
        // To disable the key 'R'.
        if( GM_.nb_geological_entities( Layer::type_name_static() ) == 0 ) {
            show_colored_layers_.new_status = false ;
            return ;
        }
        colored_cells_.new_status = false ;
        show_colored_regions_.new_status = false ;
        for( GEO::index_t l = 0;
            l < GM_.nb_geological_entities( Layer::type_name_static() ); l++ ) {
            float red = std::fmod( GEO::Numeric::random_float32(), 1 ) ;
            float green = std::fmod( GEO::Numeric::random_float32(), 1 ) ;
            float blue = std::fmod( GEO::Numeric::random_float32(), 1 ) ;
            const GeoModelGeologicalEntity& cur_layer = GM_.geological_entity(
                Layer::type_name_static(), l ) ;
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
            GM_gfx_.corners.GeoModelGfxManager::set_mesh_element_color(
                corner_style_.color_.Value.x, corner_style_.color_.Value.y,
                corner_style_.color_.Value.z ) ;
            GM_gfx_.corners.GeoModelGfxManager::set_mesh_element_size(
                static_cast< index_t >( corner_style_.size_ ) ) ;
            GM_gfx_.corners.draw() ;
        }

        if( show_lines_ ) {
            GM_gfx_.lines.GeoModelGfxManager::set_mesh_element_color(
                line_style_.color_.Value.x, line_style_.color_.Value.y,
                line_style_.color_.Value.z ) ;
            GM_gfx_.lines.GeoModelGfxManager::set_mesh_element_size(
                static_cast< index_t >( line_style_.size_ ) ) ;
            GM_gfx_.lines.draw() ;
        }

        if( show_surface_ ) {
            GM_gfx_.surfaces.set_mesh_color( mesh_color_.Value.x,
                mesh_color_.Value.y, mesh_color_.Value.z ) ;
            GM_gfx_.surfaces.GeoModelGfxManager::set_mesh_element_color(
                surface_style_.color_.Value.x, surface_style_.color_.Value.y,
                surface_style_.color_.Value.z ) ;
            GM_gfx_.surfaces.set_mesh_size(
                static_cast< index_t >( surface_style_.size_ ) ) ;
            if( selected_entity_type_ == 0 ) {
                for( GEO::index_t s = 0; s < GM_.nb_surfaces(); s++ ) {
                    if( GM_.surface( s ).is_on_voi() ) {
                        GM_gfx_.surfaces.set_mesh_element_visibility( s,
                            show_voi_ ) ;
                    }
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
                GM_gfx_.regions.set_mesh_color( mesh_color_.Value.x,
                    mesh_color_.Value.y, mesh_color_.Value.z ) ;
                GM_gfx_.regions.GeoModelGfxManager::set_mesh_element_color(
                    volume_style_.color_.Value.x, volume_style_.color_.Value.y,
                    volume_style_.color_.Value.z ) ;
            }
            GM_gfx_.regions.set_mesh_size(
                static_cast< index_t >( volume_style_.size_ ) ) ;
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
        attribute_max_ = static_cast< float >( GM_gfx_.attribute.maximum() ) ;
        attribute_min_ = static_cast< float >( GM_gfx_.attribute.minimum() ) ;
    }

    void RINGMeshApplication::GeoModelViewer::update_entity_visibility()
    {
        index_t selected_entity_type_casted =
            static_cast< index_t >( selected_entity_type_ ) ;
        const std::string& type = entity_types_[selected_entity_type_casted] ;
        if( selected_entity_type_ == 0 ) {
            GM_gfx_.corners.GeoModelGfxManager::set_mesh_element_visibility( true ) ;
            GM_gfx_.lines.GeoModelGfxManager::set_mesh_element_visibility( true ) ;
            GM_gfx_.surfaces.GeoModelGfxManager::set_mesh_element_visibility(
                true ) ;
            GM_gfx_.regions.GeoModelGfxManager::set_mesh_element_visibility( true ) ;
        } else {
            GM_gfx_.corners.GeoModelGfxManager::set_mesh_element_visibility(
                false ) ;
            GM_gfx_.lines.GeoModelGfxManager::set_mesh_element_visibility( false ) ;
            GM_gfx_.surfaces.GeoModelGfxManager::set_mesh_element_visibility(
                false ) ;
            GM_gfx_.regions.GeoModelGfxManager::set_mesh_element_visibility(
                false ) ;
            if( selected_entity_type_casted
                < EntityTypeManager::nb_mesh_entity_types() + 1 ) {
                selected_entity_id_ = std::min(
                    static_cast< int >( GM_.nb_mesh_entities( type ) - 1 ),
                    selected_entity_id_ ) ;
                gme_t entity_id( type,
                    static_cast< index_t >( selected_entity_id_ ) ) ;
                toggle_mesh_entity_and_boundaries_visibility( entity_id ) ;
            } else {
                selected_entity_id_ = std::min(
                    static_cast< int >( GM_.nb_geological_entities( type ) - 1 ),
                    selected_entity_id_ ) ;
                gme_t entity_id( type,
                    static_cast< index_t >( selected_entity_id_ ) ) ;
                toggle_geological_entity_visibility( entity_id ) ;
            }
        }
    }

    void RINGMeshApplication::GeoModelViewer::toggle_mesh_entity_and_boundaries_visibility(
        const gme_t& entity_id )
    {
        if( EntityTypeManager::is_corner( entity_id.type ) ) {
            toggle_corner_visibility( entity_id.index ) ;
        } else if( EntityTypeManager::is_line( entity_id.type ) ) {
            toggle_line_and_boundaries_visibility( entity_id.index ) ;
        } else if( EntityTypeManager::is_surface( entity_id.type ) ) {
            toggle_surface_and_boundaries_visibility( entity_id.index ) ;
        } else if( EntityTypeManager::is_region( entity_id.type ) ) {
            toggle_region_and_boundaries_visibility( entity_id.index ) ;
        } else {
            ringmesh_assert_not_reached ;
        }
    }

    void RINGMeshApplication::GeoModelViewer::toggle_corner_visibility(
        index_t corner_id )
    {
        GM_gfx_.corners.set_mesh_element_visibility( corner_id, true ) ;
    }

    void RINGMeshApplication::GeoModelViewer::toggle_line_and_boundaries_visibility(
        index_t line_id )
    {
        GM_gfx_.lines.set_mesh_element_visibility( line_id, true ) ;
        const Line& line = GM_.line( line_id ) ;
        toggle_corner_visibility( line.boundary_gme( 0 ).index ) ;
        toggle_corner_visibility( line.boundary_gme( 1 ).index ) ;
    }

    void RINGMeshApplication::GeoModelViewer::toggle_surface_and_boundaries_visibility(
        index_t surface_id )
    {
        GM_gfx_.surfaces.set_mesh_element_visibility( surface_id, true ) ;
        const Surface& surface = GM_.surface( surface_id ) ;
        for( index_t i = 0; i < surface.nb_boundaries(); i++ ) {
            toggle_line_and_boundaries_visibility(
                surface.boundary_gme( i ).index ) ;
        }
    }

    void RINGMeshApplication::GeoModelViewer::toggle_region_and_boundaries_visibility(
        index_t region_id )
    {
        GM_gfx_.regions.set_mesh_element_visibility( region_id, true ) ;
        const Region& region = GM_.region( region_id ) ;
        for( index_t i = 0; i < region.nb_boundaries(); i++ ) {
            toggle_surface_and_boundaries_visibility(
                region.boundary_gme( i ).index ) ;
        }
    }

    void RINGMeshApplication::GeoModelViewer::toggle_geological_entity_visibility(
        const gme_t& entity_id )
    {
        const GeoModelGeologicalEntity& entity = GM_.geological_entity( entity_id ) ;
        for( index_t i = 0; i < entity.nb_children(); i++ ) {
            const gme_t& child_id = entity.child_gme( i ) ;
            toggle_mesh_entity_and_boundaries_visibility( child_id ) ;
        }
    }

    void RINGMeshApplication::GeoModelViewer::draw_object_properties()
    {
        if( ImGui::Combo( "Type", &selected_entity_type_, GetChar,
            static_cast< void* >( &entity_types_ ),
            static_cast< int >( entity_types_.size() ) ) ) {
            update_entity_visibility() ;
        }
        if( selected_entity_type_ > 0 ) {
            if( ImGui::InputInt( "Id", &selected_entity_id_, 1 ) ) {
                selected_entity_id_ = std::max( 0, selected_entity_id_ ) ;
                update_entity_visibility() ;
            }
        }
        ImGui::Separator() ;
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
                    case AttributeGfxManager::facets:
                        set_attribute_names(
                            GM_.surface( 0 ).facet_attribute_manager() ) ;
                        break ;
                    case AttributeGfxManager::facet_vertices:
                        set_attribute_names(
                            GM_.surface( 0 ).vertex_attribute_manager() ) ;
                        break ;
                    case AttributeGfxManager::cells:
                        set_attribute_names(
                            GM_.region( 0 ).cell_attribute_manager() ) ;
                        break ;
                    case AttributeGfxManager::cell_vertices:
                        set_attribute_names(
                            GM_.region( 0 ).vertex_attribute_manager() ) ;
                        break ;
                    default:
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
        ImGui::SameLine() ;
        ImGui::PushStyleColor( ImGuiCol_Button, mesh_color_ ) ;
        if( ImGui::Button( "  ##MeshColor" ) ) {
            ImGui::OpenPopup( "##MeshColorTable" ) ;
        }
        ImGui::PopStyleColor() ;
        if( ImGui::BeginPopup( "##MeshColorTable" ) ) {
            show_color_table_popup( mesh_color_ ) ;
        }

        ImGui::Separator() ;
        ImGui::Checkbox( "Corner [c]", &show_corners_ ) ;
        draw_entity_style_editor( "##CornerColor", corner_style_ ) ;

        ImGui::Checkbox( "Line [e]", &show_lines_ ) ;
        draw_entity_style_editor( "##LineColor", line_style_ ) ;

        ImGui::Checkbox( "Surface [s]", &show_surface_ ) ;
        draw_entity_style_editor( "##SurfaceColor", surface_style_ ) ;

        if( meshed_regions_ ) {
            ImGui::Separator() ;
            ImGui::Checkbox( "Region [v]", &show_volume_ ) ;
            draw_entity_style_editor( "##VolumeColor", volume_style_ ) ;
            if( show_volume_ ) {
                ImGui::Checkbox( "Col. cells [C]", &colored_cells_.new_status ) ;
                ImGui::Checkbox( "Col. regions [r]",
                    &show_colored_regions_.new_status ) ;
                if( GM_.nb_geological_entities( Layer::type_name_static() ) != 0 ) {
                    ImGui::Checkbox( "Col. layers [R]",
                        &show_colored_layers_.new_status ) ;
                }
                ImGui::SliderFloat( "Shrk.", &shrink_, 0.0f, 1.0f, "%.1f" ) ;
                ImGui::Checkbox( "Hex", &show_hex_ ) ;
                ImGui::Checkbox( "Prism", &show_prism_ ) ;
                ImGui::Checkbox( "Pyramid", &show_pyramid_ ) ;
                ImGui::Checkbox( "Tetra", &show_tetra_ ) ;
            }
        }
    }

    void RINGMeshApplication::GeoModelViewer::draw_entity_style_editor(
        const std::string& label,
        EntityStyle& style )
    {
        ImGui::PushStyleColor( ImGuiCol_Button, style.color_ ) ;
        if( ImGui::Button( ( "  " + label ).c_str() ) ) {
            ImGui::OpenPopup( label.c_str() ) ;
        }
        ImGui::PopStyleColor() ;
        if( ImGui::BeginPopup( label.c_str() ) ) {
            show_color_table_popup( style.color_ ) ;
        }
        ImGui::SameLine() ;
        ImGui::InputInt( "", &style.size_, 1 ) ;
        style.size_ = std::max( style.size_, 0 ) ;
    }

    void RINGMeshApplication::GeoModelViewer::show_color_table_popup(
        ImColor& color )
    {
        int id = 0 ;
        for( index_t i = 0; i < color_table_.size(); i++ ) {
            for( index_t j = 0; j < color_table_[i].size(); j++ ) {
                if( j > 0 ) {
                    ImGui::SameLine() ;
                }
                ImGui::PushID( id++ ) ;
                ImGui::PushStyleColor( ImGuiCol_Button, color_table_[i][j] ) ;
                if( ImGui::Button( "  " ) ) {
                    color = color_table_[i][j] ;
                    ImGui::CloseCurrentPopup() ;
                }
                ImGui::PopStyleColor() ;
                ImGui::PopID() ;
            }
        }
        ImGui::EndPopup() ;
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
        const float w = 0.3f ;
        const float h = 0.1f ;
        const float x1 = 0.f ;
        const float y1 = -0.9f ;
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

        std::string min_value = GEO::String::to_string(
            GM_gfx_.attribute.minimum() ) ;
        float nb_min_letter = static_cast< float >( min_value.size() ) ;
        glQuickText::printfAt( x1 - w - font_height * nb_min_letter * 0.3,
            y1 - font_height, z, font_sz, min_value.c_str() ) ;

        std::string max_value = GEO::String::to_string(
            GM_gfx_.attribute.maximum() ) ;
        float nb_max_letter = static_cast< float >( max_value.size() ) ;
        glQuickText::printfAt( x1 + w - font_height * nb_max_letter * 0.3,
            y1 - font_height, z, font_sz, max_value.c_str() ) ;

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

        if( !filename.empty() ) {
            GEO::mesh_load( filename, mesh_ ) ;
            name_ = GEO::FileSystem::base_name( filename, true ) ;
        }
        mesh_gfx_.set_mesh( &mesh_ ) ;

        for( index_t v = 0; v < mesh_.vertices.nb(); v++ ) {
            bbox_.add_point( mesh_.vertices.point( v ) ) ;
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
            const GEO::MeshSubElementsStore& subelements =
                mesh_.get_subelements_by_type( attribute_subelements_ ) ;
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
        GEO::String::split_string( attribute_, '.', subelements_name,
            attribute_name_ ) ;
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

        Logger::div( "RINGMeshView" ) ;
        Logger::out( "" ) << "Welcome to RINGMeshView !" << std::endl ;
        Logger::out( "" ) << "This project is developped by the RINGMesh"
            << " developpers team:" << std::endl ;
        Logger::out( "" ) << "RINGMesh-dev <georessources-ringmesh-dev@univ-lorraine.fr> "
            << std::endl ;
        Logger::out( "" ) << "You can have access to the full code through "
            << "its Bitbucket repository: "	<< std::endl ;
        Logger::out( "" ) << "https://bitbucket.org/ring_team/ringmesh" 
            << std::endl ; 
        Logger::out( "" ) << "More information on this project and other " 
            << "projects of the team: " << std::endl ;
        Logger::out( "" ) << "http://www.ring-team.org" << std::endl ; 
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
        std::sort( files.begin(), files.end() ) ;
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

    bool RINGMeshApplication::can_load_geogram( const std::string& filename )
    {
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
            current_viewer_ = static_cast< index_t >( meshes_.size() - 1 ) ;
            current_viewer_type_ = MESH ;
        }

        update_region_of_interest() ;
        return true ;
    }

    void RINGMeshApplication::draw_application_menus()
    {
        if( ImGui::BeginMenu( "Debug" ) ) {
            if( ImGui::BeginMenu( "Load..." ) ) {
                ImGui::Selectable( ".." ) ;
                if( ImGui::IsItemClicked() ) {
                    path_ += "/.." ;
                }
                browse_geogram( path_ ) ;
                ImGui::EndMenu() ;
            }
            if( ImGui::MenuItem( "Create point" ) ) {
                GEO::Command::set_current( "create_point(std::string name=\"debug\","
                    " double x=0, double y=0, double z=0)", this,
                    &RINGMeshApplication::create_point ) ;
            }
            ImGui::EndMenu() ;
        }
    }

    void RINGMeshApplication::create_point(
        std::string name,
        double x,
        double y,
        double z )
    {
        MeshViewer* viewer = nil ;
        for( index_t i = 0; i < meshes_.size(); i++ ) {
            if( meshes_[i]->name_ == name ) {
                viewer = meshes_[i] ;
                break ;
            }
        }
        if( !viewer ) {
            meshes_.push_back( new MeshViewer( *this, "" ) ) ;
            viewer = meshes_.back() ;
        }
        vec3 point( x, y, z ) ;
        viewer->mesh_.vertices.create_vertex( point.data() ) ;
        viewer->mesh_gfx_.set_mesh( &viewer->mesh_ ) ;
        viewer->bbox_.add_point( point ) ;
        viewer->name_ = name ;
        current_viewer_ = static_cast< index_t >( meshes_.size() - 1 ) ;
        current_viewer_type_ = MESH ;
        update_region_of_interest() ;
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
            current_viewer_ = static_cast< index_t >( models_.size() - 1 ) ;
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
                ImGui::PushID( static_cast< int >( i ) ) ;
                if( ImGui::Checkbox( viewer.GM_.name().c_str(),
                    &viewer.is_visible_ ) ) {
                    current_viewer_ = i ;
                    current_viewer_type_ = GEOMODEL ;
                    update_region_of_interest() ;
                }
                ImGui::PopID() ;
            }
        }

        if( !meshes_.empty() ) {
            ImGui::Separator() ;
            ImGui::Text( "Mesh" ) ;
            for( index_t i = 0; i < meshes_.size(); i++ ) {
                MeshViewer& viewer = *meshes_[i] ;
                ImGui::PushID( static_cast< int >( i ) ) ;
                if( ImGui::Checkbox( viewer.name_.c_str(), &viewer.is_visible_ ) ) {
                    current_viewer_ = i ;
                    current_viewer_type_ = MESH ;
                    update_region_of_interest() ;
                }
                ImGui::PopID() ;
            }
        }
    }

    void RINGMeshApplication::draw_object_properties()
    {
        if( current_viewer_ == NO_ID ) return ;
        switch( current_viewer_type_ ) {
            case GEOMODEL:
                ringmesh_assert( current_viewer_ < models_.size() ) ;
                models_[current_viewer_]->draw_object_properties() ;
                return ;
            case MESH:
                ringmesh_assert( current_viewer_ < meshes_.size() ) ;
                meshes_[current_viewer_]->draw_object_properties() ;
                return ;
            default:
                return ;
        }
    }
}

#endif
