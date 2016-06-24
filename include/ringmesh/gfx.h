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

#ifndef __RINGMESH_GFX__
#define __RINGMESH_GFX__

#include <ringmesh/common.h>

#ifdef RINGMESH_WITH_GRAPHICS

#include <geogram/basic/factory.h>
#include <geogram/basic/command_line.h>

#include <geogram_gfx/glup_viewer/glup_viewer_gui.h>
#include <geogram_gfx/mesh/mesh_gfx.h>

#include <ringmesh/geo_model.h>

#include <ringmesh/geo_model.h>
#include <ringmesh/geo_model.h>
/*!
 * @file Classes for GeoModel visualization
 * @author Benjamin Chauvin and Arnaud Botella
 */

namespace RINGMesh {
    class GeoModel ;
    class CornerGfx ;
    class LineGfx ;
    class SurfaceGfx ;
    class RegionGfx ;
    class MeshEntityGfx ;
}

namespace RINGMesh {

    class RINGMESH_API GeoModelGfx {
    ringmesh_disable_copy( GeoModelGfx ) ;
    public:
        GeoModelGfx() ;
        ~GeoModelGfx() ;

        void set_geo_model( const GeoModel& model ) ;
        const GeoModel* geo_model() const ;
        void initialize() ;
        void need_to_update() ;

        void set_cell_vertex_attribute(
            const std::string& name,
            index_t coordinate,
            GLuint colormap_texture) ;
        void set_cell_attribute(
            const std::string& name,
            index_t coordinate,
            GLuint colormap_texture ) ;

        void compute_cell_vertex_attribute_range( index_t coordinate,
            const std::string& name ) ;
        void compute_cell_attribute_range(
            index_t coordinate,
            const std::string& name ) ;

        void set_attribute_max( double max ) {
            attribute_max_ = max ;
        }
        double get_attribute_max() const {
            return attribute_max_ ;
        }
        void set_attribute_min( double min ) {
            attribute_min_ = min ;
        }
        double get_attribute_min() const {
            return attribute_min_ ;
        }

        void draw_corners() ;
        void draw_lines() ;
        void draw_surfaces() ;
        void draw_regions() ;

        // Settings for the corners
        void set_corners_color( float r, float g, float b ) ;
        void set_corner_color( index_t c, float r, float g, float b ) ;
        void set_corners_visibility( bool b ) ;
        void set_corner_visibility( index_t c, bool b ) ;
        void set_corners_size( index_t s ) ;
        void set_corner_size( index_t c, index_t s ) ;

        // Settings for the lines
        void set_edge_lines_color( float r, float g, float b ) ;
        void set_edge_line_color( index_t c, float r, float g, float b ) ;
        void set_edge_lines_visibility( bool b ) ;
        void set_edge_line_visibility( index_t c, bool b ) ;
        void set_edge_lines_size( index_t s ) ;
        void set_edge_line_size( index_t c, index_t s ) ;

        void set_vertex_lines_color( float r, float g, float b ) ;
        void set_vertex_line_color( index_t c, float r, float g, float b ) ;
        void set_vertex_lines_visibility( bool b ) ;
        void set_vertex_line_visibility( index_t c, bool b ) ;
        void set_vertex_lines_size( index_t s ) ;
        void set_vertex_line_size( index_t c, index_t s ) ;

        // Settings for the surfaces
        void set_surfaces_lighting( bool value ) ;
        void set_surfaces_color( float r, float g, float b ) ;
        void set_surface_color( index_t c, float r, float g, float b ) ;
        void set_backface_surfaces_color( float r, float g, float b ) ;
        void set_backface_surface_color( index_t c, float r, float g, float b ) ;
        void set_surfaces_visibility( bool b ) ;
        void set_surface_visibility( index_t c, bool b ) ;

        void set_mesh_surfaces_color( float r, float g, float b ) ;
        void set_mesh_surface_color( index_t c, float r, float g, float b ) ;
        void set_mesh_surfaces_visibility( bool b ) ;
        void set_mesh_surface_visibility( index_t c, bool b ) ;
        void set_mesh_surfaces_size( index_t s ) ;
        void set_mesh_surface_size( index_t c, index_t s ) ;

        void set_vertex_surfaces_color( float r, float g, float b ) ;
        void set_vertex_surface_color( index_t c, float r, float g, float b ) ;
        void set_vertex_surfaces_visibility( bool b ) ;
        void set_vertex_surface_visibility( index_t c, bool b ) ;
        void set_vertex_surfaces_size( index_t s ) ;
        void set_vertex_surface_size( index_t c, index_t s ) ;

        // Settings for the regions
        void set_vertex_regions_color( float m, float g, float b ) ;
        void set_vertex_region_color( index_t m, float r, float g, float b ) ;
        void set_vertex_regions_visibility( bool b ) ;
        void set_vertex_region_visibility( index_t m, bool b ) ;
        void set_vertex_regions_size( index_t s ) ;
        void set_vertex_region_size( index_t m, index_t s ) ;

        void set_edge_regions_color( float r, float g, float b ) ;
        void set_edge_region_color( index_t m, float r, float g, float b ) ;
        void set_edge_regions_visibility( bool b ) ;
        void set_edge_region_visibility( index_t m, bool b ) ;
        void set_edge_regions_size( index_t s ) ;
        void set_edge_region_size( index_t l, index_t s ) ;

        void set_cell_mesh_regions_color( float r, float g, float b ) ;
        void set_cell_mesh_region_color( index_t m, float r, float g, float b ) ;
        void set_cell_mesh_regions_visibility( bool b ) ;
        void set_cell_mesh_region_visibility( index_t m, bool b ) ;
        void set_cell_mesh_regions_size( index_t s ) ;
        void set_cell_mesh_region_size( index_t m, index_t s ) ;

        void set_cell_regions_lighting( bool value )  ;
        void set_cell_regions_color( float m, float g, float b ) ;
        void set_cell_region_color( index_t m, float r, float g, float b ) ;
        void set_cell_regions_color_type() ;
        void set_cell_region_color_type( index_t m ) ;
        void set_cell_regions_visibility( bool b ) ;
        void set_cell_region_visibility( index_t m, bool b ) ;
        void set_cell_regions_type_visibility( GEO::MeshCellType t, bool b ) ;
        void set_cell_region_type_visibility(
            index_t m,
            GEO::MeshCellType t,
            bool b ) ;
        void set_cell_regions_shrink( double s ) ;
        void set_cell_region_shrink( index_t m, double s ) ;

    private:

    private:
        /// The GeoModel associated to the graphics
        const GeoModel* model_ ;

        /// The graphics associated to each Corner
        std::vector< CornerGfx* > corners_ ;
        /// The graphics associated to each Line
        std::vector< LineGfx* > lines_ ;
        /// The graphics associated to each Surface
        std::vector< SurfaceGfx* > surfaces_ ;
        /// The graphics associated to each Region
        std::vector< RegionGfx* > regions_ ;
        double attribute_max_ ;
        double attribute_min_ ;
    } ;

    class RINGMESH_API RINGMeshApplication: public GEO::Application {
    public:
        RINGMeshApplication( int argc, char** argv ) ;
        ~RINGMeshApplication() ;
    private:
        static RINGMeshApplication* instance() ;

        virtual std::string supported_read_file_extensions();
        virtual void init_graphics() ;
        virtual bool load( const std::string& filename ) ;
        virtual void draw_scene() ;

        static void toggle_points() ;
        static void toggle_volume() ;
        static void toggle_voi() ;
        static void toggle_mesh() ;
        static void increment_shrink() ;
        static void decrement_shrink() ;
        static void toggle_colored_cells() ;
        static void toggle_colored_regions() ;
        static void toggle_colored_layers() ;

    private:
        GeoModel* GM_ ;
        GeoModelGfx GM_gfx_ ;
        std::string file_extensions_;

        bool show_borders_ ;
        bool show_corners_ ;
        bool show_lines_ ;
        bool show_surface_ ;
        bool show_volume_ ;
        bool colored_cells_ ;
        bool show_voi_ ;
        bool show_colored_regions_ ;
        bool show_colored_layers_ ;
        bool show_points_ ;

        float shrink_ ;
        bool mesh_visible_ ;
        bool meshed_regions_ ;

    } ;
}

#endif
#endif
