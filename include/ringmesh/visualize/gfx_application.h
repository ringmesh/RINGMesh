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

#pragma once

#include <ringmesh/visualize/common.h>

#ifdef RINGMESH_WITH_GRAPHICS

#include <geogram_gfx/glup_viewer/glup_viewer_gui.h>

#include <ringmesh/basic/box.h>
#include <ringmesh/geomodel/builder/geomodel_builder.h>
#include <ringmesh/geomodel/core/geomodel.h>

#include <ringmesh/visualize/geomodel_gfx.h>

/*!
 * @file Classes for GeoModel visualization
 * @author Benjamin Chauvin and Arnaud Botella
 */

namespace RINGMesh
{
    class visualize_api RINGMeshApplication : public GEO::Application
    {
        ringmesh_disable_copy_and_move( RINGMeshApplication );

    public:
        enum struct ViewerType
        {
            GEOMODEL2D,
            GEOMODEL3D,
            MESH,
            NONE
        };

        RINGMeshApplication( int argc, char** argv );
        ~RINGMeshApplication() = default;

        virtual void quit();

    protected:
        static RINGMeshApplication* instance();

        std::string supported_read_file_extensions() override;
        std::string supported_geogram_read_file_extensions();
        void init_graphics() override;
        bool load( const std::string& filename ) override;
        bool save( const std::string& filename ) override;
        void draw_scene() override;
        void draw_object_properties() override;
        void draw_viewer_properties() override;
        void draw_application_menus() override;

        bool load_geogram( const std::string& filename );
        bool can_load_geogram( const std::string& filename );
        void browse_geogram( const std::string& path );
        void update_region_of_interest();
        void init_ringmesh_colormaps();

    private:
        void create_point( std::string name = "debug",
            double x = 0,
            double y = 0,
            double z = 0 );

        void create_aabbox( std::string name = "box",
            double xmin = 0,
            double ymin = 0,
            double zmin = 0,
            double xmax = 1,
            double ymax = 1,
            double zmax = 1 );

    private:
        template < index_t DIMENSION >
        class GeoModelViewerBase
        {
            ringmesh_disable_copy_and_move( GeoModelViewerBase );

        public:
            struct OldNewStatus
            {
                OldNewStatus& operator=( bool value )
                {
                    old_status = value;
                    new_status = value;
                    return *this;
                }
                bool need_to_update() const
                {
                    return old_status != new_status;
                }
                void update()
                {
                    old_status = new_status;
                }
                bool old_status{ false };
                bool new_status{ false };
            };
            struct EntityStyle
            {
                ImVec4 color_;
                int size_{ 1 };
                bool visible_vertices_{ false };
                ImVec4 vertex_color_;
                int vertex_size_{ 1 };
            };

        public:
            GeoModelViewerBase(
                RINGMeshApplication& app, const std::string& filename );
            virtual ~GeoModelViewerBase() = default;

            virtual void draw_scene();
            virtual void draw_object_properties();
            void draw_viewer_properties();
            void draw_colormap();

            void reset_attribute_name();
            void set_attribute_names( const std::vector< std::string >& names );
            void autorange();
            void draw_entity_style_editor(
                const std::string& label, EntityStyle& style );
            void draw_entity_vertex_style_editor(
                const std::string& label, EntityStyle& style );
            void update_entity_visibility();
            virtual void update_all_entity_visibility( bool value );

            void toggle_corner_visibility( index_t corner_id );
            void toggle_line_and_boundaries_visibility( index_t line_id );
            void toggle_surface_and_boundaries_visibility( index_t surface_id );
            void toggle_geological_entity_visibility(
                const gmge_id& entity_id );
            virtual void toggle_mesh_entity_and_boundaries_visibility(
                const gmme_id& entity_id );

            virtual ViewerType type() = 0;

        public:
            /*this struc exist only to enable the use of an uncompressed
            std::vector<Bool> by using a std::vector<SpeBool>*/
            struct SpeBool
            {
                bool value_{ true };
            };

            RINGMeshApplication& app_;
            bool is_visible_{ true };
            GeoModel< DIMENSION > GM_;
            GeoModelBuilder< DIMENSION > GM_builder_;
            GeoModelGfx< DIMENSION > GM_gfx_;
            Box< DIMENSION > bbox_;
            std::vector< std::string > entity_types_;
            int selected_entity_type_{ 0 };
            int selected_entity_id_{ 0 };

            // add variable here

            bool show_corners_{ true };
            EntityStyle corner_style_;
            bool show_lines_{ true };
            EntityStyle line_style_;
            bool show_surface_{ true };
            EntityStyle surface_style_;
            bool show_voi_{ false };
            bool show_colormap_{ false };

            std::vector< SpeBool > surface_is_visible_; // determine if each
                                                        // surface is visible or
                                                        // not
            std::vector< SpeBool > line_is_visible_; // eq for lines
            std::vector< SpeBool > corner_is_visible_; // eq for corners

            std::vector< char* >
                temp_surface_name_; // use to modify surface names
            std::vector< char* > temp_line_name_; // eq for line
            std::vector< char* > temp_corner_name_; // eq for corners
            char* temp_name_; // eq for GM.name_

            bool mesh_visible_{ true };
            ImVec4 mesh_color_;
            bool show_attributes_{ false };
            float attribute_min_{ 0 };
            float attribute_max_{ 0 };
        };

        ALIAS_2D_AND_3D( GeoModelViewerBase );

        template < index_t DIMENSION >
        class visualize_api GeoModelViewer final
            : public GeoModelViewerBase< DIMENSION >
        {
        };

        ALIAS_2D_AND_3D( GeoModelViewer );

        class MeshViewer
        {
            ringmesh_disable_copy_and_move( MeshViewer );

        public:
            MeshViewer( RINGMeshApplication& app, const std::string& filename );
            ~MeshViewer() = default;

            void draw_object_properties();
            void draw_scene();

            void autorange();
            std::string attribute_names();
            void set_attribute( const std::string& attribute );

        public:
            RINGMeshApplication& app_;
            bool is_visible_{ true };
            GEO::Mesh mesh_;
            GEO::MeshGfx mesh_gfx_;
            Box3D bbox_;
            std::string name_;
            char* temp_name_; // use to modify name

            bool show_vertices_{ false };
            float vertices_size_{ 1 };
            ImVec4 vertices_color_;

            bool show_edges_{ true };
            int edges_size_{ 1 };
            ImVec4 edges_color_{ 0, 0, 0, 0 };

            bool show_surface_{ true };
            bool show_surface_colors_{ true };
            bool show_mesh_{ true };
            bool show_surface_borders_{ false };

            bool show_volume_{ false };
            float cells_shrink_{ 0 };
            bool show_colored_cells_{ false };
            bool show_hexes_{ true };

            bool show_attributes_{ false };
            GLuint current_colormap_texture_{ 0 };
            std::string attribute_ = std::string{ "vertices.point_fp32[0]" };
            GEO::MeshElementsFlags attribute_subelements_{ GEO::MESH_VERTICES };
            std::string attribute_name_;
            float attribute_min_{ 0 };
            float attribute_max_{ 0 };
        };

        template < index_t DIMENSION >
        void draw_geomodel_viewer_properties(
            std::vector< std::unique_ptr< GeoModelViewer< DIMENSION > > >&
                geomodels,
            int& id );

    protected:
        std::vector< std::unique_ptr< GeoModelViewer3D > > geomodels3d_;
        std::vector< std::unique_ptr< GeoModelViewer2D > > geomodels2d_;
        std::vector< std::unique_ptr< MeshViewer > > meshes_;
        std::string ringmesh_file_extensions_;
        std::string geogram_file_extensions_;
        index_t current_viewer_{ NO_ID };
        ViewerType current_viewer_type_{ ViewerType::NONE };

        ImVec4 backgound_color_{ 0., 0., 0., 0. };
    };

    template <>
    class RINGMeshApplication::GeoModelViewer< 2 > final
        : public GeoModelViewerBase< 2 >
    {
    public:
        GeoModelViewer( RINGMeshApplication& app, const std::string& filename );
        ViewerType type() override
        {
            return ViewerType::GEOMODEL2D;
        }
    };

    template <>
    class RINGMeshApplication::GeoModelViewer< 3 > final
        : public GeoModelViewerBase< 3 >
    {
    public:
        GeoModelViewer( RINGMeshApplication& app, const std::string& filename );
        ViewerType type() override
        {
            return ViewerType::GEOMODEL3D;
        }

        void draw_scene() override;
        void draw_object_properties() override;
        void update_all_entity_visibility( bool value ) override;

        void toggle_colored_cells();
        void toggle_colored_regions();
        void toggle_colored_layers();
        void toggle_region_and_boundaries_visibility( index_t region_id );
        void toggle_mesh_entity_and_boundaries_visibility(
            const gmme_id& entity_id ) override;

    public:
        bool show_hex_{ true };
        bool show_prism_{ true };
        bool show_pyramid_{ true };
        bool show_tetra_{ true };

        float shrink_{ 0 };
        bool meshed_regions_{ false };
        bool show_volume_{ false };
        EntityStyle volume_style_;

        OldNewStatus colored_cells_;
        OldNewStatus show_colored_regions_;
        OldNewStatus show_colored_layers_;
    };
} // namespace RINGMesh

#endif
