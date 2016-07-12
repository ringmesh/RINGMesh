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

#include <geogram_gfx/glup_viewer/glup_viewer_gui.h>

/*!
 * @file Classes for GeoModel visualization
 * @author Benjamin Chauvin and Arnaud Botella
 */

namespace RINGMesh {
    class GeoModel ;
    class GeoModelGfx ;
    class AttributeGfx ;
    class CornerGfx ;
    class LineGfx ;
    class SurfaceGfx ;
    class RegionGfx ;
    class MeshEntityGfx ;
}

namespace RINGMesh {

    class RINGMESH_API GeoModelGfxManager {
    ringmesh_disable_copy( GeoModelGfxManager ) ;
    public:
        GeoModelGfxManager( GeoModelGfx& gfx ) ;
        virtual ~GeoModelGfxManager() ;

        virtual void draw() = 0 ;
        virtual void initialize() = 0 ;
        void need_to_update() ;
        void set_scalar_attribute(
            GEO::MeshElementsFlags subelements,
            const std::string& name,
            double attr_min,
            double attr_max,
            GLuint colormap_texture ) ;
        void unset_scalar_attribute() ;

        void set_vertex_color( float r, float g, float b ) ;
        void set_vertex_color( index_t e, float r, float g, float b ) ;
        void set_vertex_visibility( bool b ) ;
        void set_vertex_visibility( index_t e, bool b ) ;
        void set_vertex_size( index_t s ) ;
        void set_vertex_size( index_t e, index_t s ) ;

        void set_mesh_element_color( float r, float g, float b ) ;
        void set_mesh_element_visibility( bool b ) ;
        void set_mesh_element_size( index_t s ) ;
        virtual void set_mesh_element_color( index_t e, float r, float g, float b ) ;
        virtual void set_mesh_element_visibility( index_t e, bool b ) ;
        virtual void set_mesh_element_size( index_t e, index_t s ) ;

    protected:
        GeoModelGfx& gfx_ ;
        std::vector< MeshEntityGfx* > entities_ ;

    } ;

    class RINGMESH_API CornerGfxManager: public GeoModelGfxManager {
    public:
        CornerGfxManager( GeoModelGfx& gfx ) ;

        virtual void draw() ;
        virtual void initialize() ;

        virtual void set_mesh_element_color( index_t e, float r, float g, float b ) ;
        virtual void set_mesh_element_visibility( index_t e, bool b ) ;
        virtual void set_mesh_element_size( index_t e, index_t s ) ;

    } ;

    class RINGMESH_API LineGfxManager: public GeoModelGfxManager {
    public:
        LineGfxManager( GeoModelGfx& gfx ) ;

        virtual void draw() ;
        virtual void initialize() ;

        virtual void set_mesh_element_color( index_t e, float r, float g, float b ) ;
        virtual void set_mesh_element_visibility( index_t e, bool b ) ;
        virtual void set_mesh_element_size( index_t e, index_t s ) ;

    } ;



    class RINGMESH_API SurfaceGfxManager: public GeoModelGfxManager {
    public:
        SurfaceGfxManager( GeoModelGfx& gfx ) ;

        virtual void draw() ;
        virtual void initialize() ;

        virtual void set_mesh_element_color( index_t e, float r, float g, float b ) ;
        virtual void set_mesh_element_visibility( index_t e, bool b ) ;

        void set_backface_surfaces_color( float r, float g, float b ) ;
        void set_backface_surface_color( index_t c, float r, float g, float b ) ;

        void set_mesh_color( float r, float g, float b ) ;
        void set_mesh_color( index_t c, float r, float g, float b ) ;
        void set_mesh_visibility( bool b ) ;
        void set_mesh_visibility( index_t c, bool b ) ;
        void set_mesh_size( index_t s ) ;
        void set_mesh_size( index_t c, index_t s ) ;
    } ;

    class RINGMESH_API RegionGfxManager: public GeoModelGfxManager {
    public:
        RegionGfxManager( GeoModelGfx& gfx ) ;

        virtual void draw() ;
        virtual void initialize() ;

        virtual void set_mesh_element_color( index_t e, float r, float g, float b ) ;
        virtual void set_mesh_element_visibility( index_t e, bool b ) ;

        void set_edge_color( float r, float g, float b ) ;
        void set_edge_color( index_t m, float r, float g, float b ) ;
        void set_edge_visibility( bool b ) ;
        void set_edge_visibility( index_t m, bool b ) ;
        void set_edge_size( index_t s ) ;
        void set_edge_size( index_t l, index_t s ) ;

        void set_mesh_color( float r, float g, float b ) ;
        void set_mesh_color( index_t m, float r, float g, float b ) ;
        void set_mesh_visibility( bool b ) ;
        void set_mesh_visibility( index_t m, bool b ) ;
        void set_mesh_size( index_t s ) ;
        void set_mesh_size( index_t m, index_t s ) ;

        void set_draw_cells( GEO::MeshCellType type, bool x ) ;
        void set_draw_cells( index_t m, GEO::MeshCellType type, bool x ) ;
        void set_color_cell_type() ;
        void set_color_cell_type( index_t m ) ;
        void set_cell_type_visibility( GEO::MeshCellType t, bool b ) ;
        void set_cell_type_visibility( index_t m, GEO::MeshCellType t, bool b ) ;
        void set_shrink( double s ) ;
        void set_shrink( index_t m, double s ) ;

    } ;


    class RINGMESH_API AttributeGfxManager {
    ringmesh_disable_copy( AttributeGfxManager ) ;
    public:
        enum Attribute_location {
            cells, cell_vertices, nb_locations
        } ;
        AttributeGfxManager( GeoModelGfx& gfx ) ;
        ~AttributeGfxManager() ;

        GeoModelGfx& gfx() {
            return gfx_ ;
        }

        void compute_range() ;

        void unbind_attribute() ;
        void bind_attribute() ;

        std::string location_name( Attribute_location location ) ;

        void set_maximum( double max ) {
            maximum_ = max ;
        }
        double maximum() const {
            return maximum_ ;
        }
        void set_minimum( double min ) {
            minimum_ = min ;
        }
        double minimum() const {
            return minimum_ ;
        }
        void set_location( Attribute_location location ) {
            location_ = location ;
        }
        Attribute_location location() const {
            return location_ ;
        }
        index_t nb_coordinates() const ;
        void set_coordinate( const index_t& coordinate ) {
            coordinate_ = coordinate ;
        }
        const index_t& coordinate() const {
            return coordinate_ ;
        }
        void set_name( const std::string& name ) {
            name_ = name ;
        }
        const std::string& name() const {
            return name_ ;
        }
        void set_colormap( GLuint colormap ) {
            colormap_texture_ = colormap ;
        }
        GLuint colormap() const {
            return colormap_texture_ ;
        }

    private:
        GeoModelGfx& gfx_ ;

        std::string name_ ;
        Attribute_location location_ ;
        index_t coordinate_ ;
        GLuint colormap_texture_;
        double minimum_ ;
        double maximum_ ;

        AttributeGfx* attributes_[nb_locations] ;

    };

    class RINGMESH_API GeoModelGfx {
    ringmesh_disable_copy( GeoModelGfx ) ;
    public:

        GeoModelGfx() ;
        ~GeoModelGfx() ;

        void set_geo_model( const GeoModel& model ) ;
        const GeoModel* geo_model() const ;
        void initialize() ;
        void need_to_update() ;

    private:
        /// The GeoModel associated to the graphics
        const GeoModel* model_ ;

    public:
        CornerGfxManager corners ;
        LineGfxManager lines ;
        SurfaceGfxManager surfaces ;
        RegionGfxManager regions ;
        AttributeGfxManager attribute ;
    } ;

    /*****************************************************************/

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
        virtual void draw_object_properties() ;
        virtual void draw_viewer_properties() ;
        void draw_colormap() ;
        void toggle_colored_cells() ;
        void toggle_colored_regions() ;
        void toggle_colored_layers() ;

        static void increment_shrink() ;
        static void decrement_shrink() ;

        void reset_attribute_name() ;
        void set_attribute_names( const GEO::AttributesManager& attributes ) ;
        void autorange() ;

        struct OldNewStatus {
            void operator=( bool value )
            {
                old_status = value ;
                new_status = value ;
            }
            bool need_to_update() const
            {
                return old_status != new_status ;
            }
            void update() {
                old_status = new_status ;
            }
            bool old_status ;
            bool new_status ;
        };

    private:
        GeoModel* GM_ ;
        GeoModelGfx GM_gfx_ ;
        std::string file_extensions_;

        bool show_corners_ ;
        bool show_lines_ ;
        bool show_surface_ ;
        bool show_volume_ ;
        bool show_voi_ ;
        OldNewStatus colored_cells_ ;
        OldNewStatus show_colored_regions_ ;
        OldNewStatus show_colored_layers_ ;
        bool show_colormap_ ;

        bool show_hex_ ;
        bool show_prism_ ;
        bool show_pyramid_ ;
        bool show_tetra_ ;

        float shrink_ ;
        bool mesh_visible_ ;
        bool meshed_regions_ ;

        bool show_attributes_;
        float attribute_min_ ;
        float attribute_max_ ;

    } ;
}

#endif
#endif
