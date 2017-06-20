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

#pragma once

#include <ringmesh/basic/common.h>

#ifdef RINGMESH_WITH_GRAPHICS

#include <geogram_gfx/glup_viewer/glup_viewer_gui.h>

#include <ringmesh/basic/box.h>
#include <ringmesh/geomodel/geomodel.h>

#include <ringmesh/visualization/geomodel_gfx.h>

/*!
 * @file Classes for GeoModel visualization
 * @author Benjamin Chauvin and Arnaud Botella
 */

namespace RINGMesh {

    class RINGMESH_API RINGMeshApplication: public GEO::Application {
    ringmesh_disable_copy( RINGMeshApplication );
    public:
        enum ViewerType {
            GEOMODEL, MESH, NONE
        };

        RINGMeshApplication( int argc, char** argv );
        virtual ~RINGMeshApplication() = default;

        virtual void quit();

    private:
        static RINGMeshApplication* instance();

        virtual std::string supported_read_file_extensions() override;
        std::string supported_geogram_read_file_extensions();
        virtual void init_graphics() override;
        virtual bool load( const std::string& filename ) override;
        virtual void draw_scene() override;
        virtual void draw_object_properties() override;
        virtual void draw_viewer_properties() override;
        virtual void draw_application_menus() override;

        bool load_geogram( const std::string& filename );
        bool can_load_geogram( const std::string& filename );
        void browse_geogram( const std::string& path );
        void update_region_of_interest();

        static void increment_shrink();
        static void decrement_shrink();
        static void show_corners();
        static void show_lines();
        static void show_surface();
        static void show_volume();
        static void show_voi();
        static void mesh_visible();
        static void show_colormap();
        static void colored_cells();
        static void show_colored_regions();
        static void show_colored_layers();

        static void show_color_table_popup( ImColor& color );

        void create_point(
            std::string name = "debug",
            double x = 0,
            double y = 0,
            double z = 0 );

    private:
        class GeoModelViewer {
        ringmesh_disable_copy( GeoModelViewer );
        public:
            struct OldNewStatus {
                void operator=( bool value )
                {
                    old_status = value;
                    new_status = value;
                }
                bool need_to_update() const
                {
                    return old_status != new_status;
                }
                void update()
                {
                    old_status = new_status;
                }
                bool old_status;
                bool new_status;
            };
            struct EntityStyle {
                ImColor color_;
                int size_;
                bool visible_vertices_;
                ImColor vertex_color_;
                int vertex_size_;
            };

        public:
            GeoModelViewer( RINGMeshApplication& app, const std::string& filename );

            void draw_scene();
            void draw_object_properties();
            void draw_viewer_properties();

            void draw_colormap();
            void toggle_colored_cells();
            void toggle_colored_regions();
            void toggle_colored_layers();

            void reset_attribute_name();
            void set_attribute_names( const GEO::AttributesManager& attributes );
            void autorange();
            void draw_entity_style_editor(
                const std::string& label,
                EntityStyle& style );
            void draw_entity_vertex_style_editor(
                const std::string& label,
                EntityStyle& style );
            void update_entity_visibility();

            void toggle_corner_visibility( index_t corner_id );
            void toggle_line_and_boundaries_visibility( index_t line_id );
            void toggle_surface_and_boundaries_visibility( index_t surface_id );
            void toggle_region_and_boundaries_visibility( index_t region_id );
            void toggle_geological_entity_visibility( const gmge_id& entity_id );
            void toggle_mesh_entity_and_boundaries_visibility(
                const gmme_id& entity_id );

        public:
            RINGMeshApplication& app_;
            bool is_visible_;
            GeoModel< 3 > GM_;
            GeoModelGfx GM_gfx_;
            Box< 3 > bbox_;
            std::vector< std::string > entity_types_;
            int selected_entity_type_;
            int selected_entity_id_;

            bool show_corners_;
            EntityStyle corner_style_;
            bool show_lines_;
            EntityStyle line_style_;
            bool show_surface_;
            EntityStyle surface_style_;
            bool show_volume_;
            EntityStyle volume_style_;
            bool show_voi_;
            OldNewStatus colored_cells_;
            OldNewStatus show_colored_regions_;
            OldNewStatus show_colored_layers_;
            bool show_colormap_;

            bool show_hex_;
            bool show_prism_;
            bool show_pyramid_;
            bool show_tetra_;

            float shrink_;
            bool mesh_visible_;
            ImColor mesh_color_;
            bool meshed_regions_;

            bool show_attributes_;
            float attribute_min_;
            float attribute_max_;
        };

        class MeshViewer {
        ringmesh_disable_copy( MeshViewer );
        public:
            MeshViewer( RINGMeshApplication& app, const std::string& filename );

            void draw_object_properties();
            void draw_scene();

            void autorange();
            std::string attribute_names();
            void set_attribute( const std::string& attribute );

        public:
            RINGMeshApplication& app_;
            bool is_visible_;
            GEO::Mesh mesh_;
            GEO::MeshGfx mesh_gfx_;
            Box< 3 > bbox_;
            std::string name_;

            bool show_vertices_;
            float vertices_size_;
            ImColor vertices_color_;

            bool show_surface_;
            bool show_surface_colors_;
            bool show_mesh_;
            bool show_surface_borders_;

            bool show_volume_;
            float cells_shrink_;
            bool show_colored_cells_;
            bool show_hexes_;

            bool show_attributes_;
            GLuint current_colormap_texture_;
            std::string attribute_;
            GEO::MeshElementsFlags attribute_subelements_;
            std::string attribute_name_;
            float attribute_min_;
            float attribute_max_;
        };

    private:
        std::vector< std::unique_ptr< GeoModelViewer > > geomodels_;
        std::vector< std::unique_ptr< MeshViewer > > meshes_;
        std::string ringmesh_file_extensions_;
        std::string geogram_file_extensions_;
        index_t current_viewer_;
        ViewerType current_viewer_type_;

        static std::vector< std::vector< ImColor > > color_table_;

    };
}

#endif
