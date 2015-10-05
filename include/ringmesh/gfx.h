///*
// * Copyright (c) 2012-2015, Association Scientifique pour la Geologie et ses Applications (ASGA)
// * All rights reserved.
// *
// * Redistribution and use in source and binary forms, with or without
// * modification, are permitted provided that the following conditions are met:
// *     * Redistributions of source code must retain the above copyright
// *       notice, this list of conditions and the following disclaimer.
// *     * Redistributions in binary form must reproduce the above copyright
// *       notice, this list of conditions and the following disclaimer in the
// *       documentation and/or other materials provided with the distribution.
// *     * Neither the name of the <organization> nor the
// *       names of its contributors may be used to endorse or promote products
// *       derived from this software without specific prior written permission.
// *
// * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
// * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// * DISCLAIMED. IN NO EVENT SHALL <COPYRIGHT HOLDER> BE LIABLE FOR ANY
// * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
// * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
// *
// *  Contacts:
// *     Arnaud.Botella@univ-lorraine.fr
// *     Antoine.Mazuyer@univ-lorraine.fr
// *     Jeanne.Pellerin@wias-berlin.de
// *
// *     http://www.ring-team.org
// *
// *     RING Project
// *     Ecole Nationale Superieure de Geologie - Georessources
// *     2 Rue du Doyen Marcel Roubault - TSA 70605
// *     54518 VANDOEUVRE-LES-NANCY
// *     FRANCE
// */
//
//#ifndef __RINGMESH_GFX__
//#define __RINGMESH_GFX__
//
//#ifdef RINGMESH_WITH_GRAPHICS
//
//#include <ringmesh/common.h>
//
//#include <geogram_gfx/mesh/mesh_gfx.h>
//
//namespace RINGMesh {
//    class GeoModel ;
//    class MacroMesh ;
//    class CornerGfx ;
//    class LineGfx ;
//    class SurfaceGfx ;
//    class RegionGfx ;
//}
//
//namespace RINGMesh {
//
//    class RINGMESH_API GeoModelGfx {
//    ringmesh_disable_copy( GeoModelGfx ) ;
//    public:
//        GeoModelGfx() ;
//        ~GeoModelGfx() ;
//
//        void set_geo_model( const GeoModel& model ) ;
//        const GeoModel* geo_model() ;
//        void initialize() ;
//
//        void draw_corners() ;
//        void draw_lines() ;
//        void draw_surfaces() ;
//
//        void set_corners_color( float r, float g, float b ) ;
//        void set_corner_color( index_t c, float r, float g, float b ) ;
//        void set_corners_visibility( bool b ) ;
//        void set_corner_visibility( index_t c, bool b ) ;
//        void set_corners_size( index_t s ) ;
//        void set_corner_size( index_t c, index_t s ) ;
//
//        void set_edge_lines_color( float r, float g, float b ) ;
//        void set_edge_line_color( index_t c, float r, float g, float b ) ;
//        void set_edge_lines_visibility( bool b ) ;
//        void set_edge_line_visibility( index_t c, bool b ) ;
//        void set_edge_lines_size( index_t s ) ;
//        void set_edge_line_size( index_t c, index_t s ) ;
//
//        void set_vertex_lines_color( float r, float g, float b ) ;
//        void set_vertex_line_color( index_t c, float r, float g, float b ) ;
//        void set_vertex_lines_visibility( bool b ) ;
//        void set_vertex_line_visibility( index_t c, bool b ) ;
//        void set_vertex_lines_size( index_t s ) ;
//        void set_vertex_line_size( index_t c, index_t s ) ;
//
//        void set_surfaces_color( float r, float g, float b ) ;
//        void set_surface_color( index_t c, float r, float g, float b ) ;
//        void set_backface_surfaces_color( float r, float g, float b ) ;
//        void set_backface_surface_color( index_t c, float r, float g, float b ) ;
//        void set_surfaces_visibility( bool b ) ;
//        void set_surface_visibility( index_t c, bool b ) ;
//
//        void set_mesh_surfaces_color( float r, float g, float b ) ;
//        void set_mesh_surface_color( index_t c, float r, float g, float b ) ;
//        void set_mesh_surfaces_visibility( bool b ) ;
//        void set_mesh_surface_visibility( index_t c, bool b ) ;
//        void set_mesh_surfaces_size( index_t s ) ;
//        void set_mesh_surface_size( index_t c, index_t s ) ;
//
//        void set_vertex_surfaces_color( float r, float g, float b ) ;
//        void set_vertex_surface_color( index_t c, float r, float g, float b ) ;
//        void set_vertex_surfaces_visibility( bool b ) ;
//        void set_vertex_surface_visibility( index_t c, bool b ) ;
//        void set_vertex_surfaces_size( index_t s ) ;
//        void set_vertex_surface_size( index_t c, index_t s ) ;
//
//    private:
//        /// The GeoModel associated to the graphics
//        const GeoModel* model_ ;
//
//        /// The graphics associated to each Corner
//        std::vector< CornerGfx* > corners_ ;
//        /// The graphics associated to each Line
//        std::vector< LineGfx* > lines_ ;
//        /// The graphics associated to each Surface
//        std::vector< SurfaceGfx* > surfaces_ ;
//
//    } ;
//
//    class RINGMESH_API MacroMeshGfx {
//    ringmesh_disable_copy( MacroMeshGfx ) ;
//    public:
//        MacroMeshGfx() ;
//        ~MacroMeshGfx() ;
//
//        const MacroMesh* macro_mesh() const {
//            return mm_ ;
//        }
//        void set_macro_mesh( const MacroMesh& model ) ;
//        void initialize() ;
//
//        void draw() ;
//
//        void set_vertex_regions_color( float m, float g, float b ) ;
//        void set_vertex_region_color( index_t m, float r, float g, float b ) ;
//        void set_vertex_regions_visibility( bool b ) ;
//        void set_vertex_region_visibility( index_t m, bool b ) ;
//        void set_vertex_regions_size( index_t s ) ;
//        void set_vertex_region_size( index_t m, index_t s ) ;
//
//        void set_edge_regions_color( float r, float g, float b ) ;
//        void set_edge_region_color( index_t m, float r, float g, float b ) ;
//        void set_edge_regions_visibility( bool b ) ;
//        void set_edge_region_visibility( index_t m, bool b ) ;
//        void set_edge_regions_size( index_t s ) ;
//        void set_edge_region_size( index_t l, index_t s ) ;
//
//        void set_surface_regions_color( float r, float g, float b ) ;
//        void set_surface_region_color( index_t m, float r, float g, float b ) ;
//        void set_backface_surface_regions_color( float r, float g, float b ) ;
//        void set_backface_surface_region_color( index_t m, float r, float g, float b ) ;
//        void set_surface_regions_visibility( bool b ) ;
//        void set_surface_region_visibility( index_t m, bool b ) ;
//        void set_mesh_surface_regions_color( float r, float g, float b ) ;
//        void set_mesh_surface_region_color( index_t m, float r, float g, float b ) ;
//        void set_mesh_surface_regions_visibility( bool b ) ;
//        void set_mesh_surface_region_visibility( index_t m, bool b ) ;
//        void set_mesh_surface_regions_size( index_t s ) ;
//        void set_mesh_surface_region_size( index_t m, index_t s ) ;
//
//        void set_cell_mesh_regions_color( float m, float g, float b ) ;
//        void set_cell_mesh_region_color( index_t m, float r, float g, float b ) ;
//        void set_cell_mesh_regions_visibility( bool b ) ;
//        void set_cell_mesh_region_visibility( index_t m, bool b ) ;
//        void set_cell_mesh_regions_size( index_t s ) ;
//        void set_cell_mesh_region_size( index_t m, index_t s ) ;
//
//        void set_cell_regions_color( float m, float g, float b ) ;
//        void set_cell_region_color( index_t m, float r, float g, float b ) ;
//        void set_cell_regions_color_type() ;
//        void set_cell_region_color_type( index_t m ) ;
//        void set_cell_regions_visibility( bool b ) ;
//        void set_cell_region_visibility( index_t m, bool b ) ;
//        void set_cell_regions_type_visibility( GEO::MeshCellType t, bool b ) ;
//        void set_cell_region_type_visibility(
//            index_t m,
//            GEO::MeshCellType t,
//            bool b ) ;
//        void set_cell_regions_shrink( double s ) ;
//        void set_cell_region_shrink( index_t m, double s ) ;
//
//    private:
//        /// The MacroMesh associated to the graphics
//        const MacroMesh* mm_ ;
//
//        /// The graphics associated to each Mesh
//        std::vector< RegionGfx* > meshes_ ;
//        std::vector< SurfaceGfx* > surfaces_ ;
//        std::vector< LineGfx* > edges_ ;
//    } ;
//}
//
//#endif
//#endif
