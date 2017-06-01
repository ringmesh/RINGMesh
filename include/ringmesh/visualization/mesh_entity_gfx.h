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

#include <memory>

#include <geogram/basic/factory.h>

#include <geogram_gfx/glup_viewer/glup_viewer_gui.h>

/*!
 * @file Classes for mesh entity visualization
 * @author Benjamin Chauvin and Arnaud Botella
 */

namespace RINGMesh {
    class GeoModelGfx;
    class AttributeGfx;

    template< index_t DIMENSION > class PointSetMesh;
    template< index_t DIMENSION > class LineMesh;
    template< index_t DIMENSION > class SurfaceMesh;
    template< index_t DIMENSION > class VolumeMesh;
}

namespace RINGMesh {

    class RINGMESH_API MeshEntityGfx {
    ringmesh_disable_copy( MeshEntityGfx );
    public:
        virtual ~MeshEntityGfx() = default;

        virtual void draw_vertices() = 0;
        virtual void set_vertex_color( float red, float green, float blue ) = 0;
        virtual void set_vertex_size( index_t s ) = 0;

        virtual void set_scalar_attribute(
            GEO::MeshElementsFlags subelements,
            const std::string& name,
            double attr_min,
            double attr_max,
            GLuint colormap_texture ) = 0;
        virtual void unset_scalar_attribute() = 0;

        void set_vertex_visible( bool is_visible )
        {
            vertex_visible_ = is_visible;
        }

        bool get_vertex_visible() const
        {
            return vertex_visible_;
        }

    protected:
        MeshEntityGfx()
        {
            set_vertex_visible( false );
        }

    private:
        bool vertex_visible_;
    };

    class RINGMESH_API PointSetMeshGfx: public MeshEntityGfx {
    public:
        virtual void set_mesh( const PointSetMesh< 3 >& mesh ) = 0;

        static std::unique_ptr< PointSetMeshGfx > create_gfx(
            const PointSetMesh< 3 >& mesh );

    protected:
        PointSetMeshGfx()
        {
            set_vertex_visible( true );
        }

    };
    using PointSetMeshGfxFactory = GEO::Factory0< PointSetMeshGfx >;
#define ringmesh_register_point_set_gfx(type) \
    geo_register_creator(RINGMesh::PointSetMeshGfxFactory, type ## Gfx, type::type_name_static())

    class RINGMESH_API LineMeshGfx: public MeshEntityGfx {
    public:
        virtual void set_mesh( const LineMesh< 3 >& mesh ) = 0;

        static std::unique_ptr< LineMeshGfx > create_gfx(
            const LineMesh< 3 >& mesh );

        void set_edge_visible( bool is_visible )
        {
            edge_visible_ = is_visible;
        }
        bool get_edge_visible() const
        {
            return edge_visible_;
        }

        virtual void draw_edges() = 0;
        virtual void set_edge_color( float red, float green, float blue ) = 0;
        virtual void set_edge_width( index_t s ) = 0;

    protected:
        LineMeshGfx()
        {
            set_vertex_visible( false );
            set_edge_visible( true );
        }

    private:
        bool edge_visible_;
    };
    using LineMeshGfxFactory = GEO::Factory0< LineMeshGfx >;
#define ringmesh_register_line_gfx(type) \
    geo_register_creator(RINGMesh::LineMeshGfxFactory, type ## Gfx, type::type_name_static())

    class RINGMESH_API SurfaceMeshGfx: public MeshEntityGfx {
    public:
        virtual void set_mesh( const SurfaceMesh< 3 >& mesh ) = 0;

        static std::unique_ptr< SurfaceMeshGfx > create_gfx(
            const SurfaceMesh< 3 >& mesh );

        virtual void draw_surface() = 0;
        virtual void set_surface_color( float red, float green, float blue ) = 0;
        virtual void set_backface_surface_color(
            float red,
            float green,
            float blue ) = 0;
        virtual void set_mesh_color( float red, float green, float blue ) = 0;
        virtual void set_mesh_visibility( bool is_visible ) = 0;
        virtual void set_mesh_width( index_t s ) = 0;

        void set_surface_visible( bool is_visible )
        {
            surface_visible_ = is_visible;
        }
        bool get_surface_visible() const
        {
            return surface_visible_;
        }

    protected:
        SurfaceMeshGfx()
        {
            set_vertex_visible( false );
            set_surface_visible( true );
        }
    private:
        bool surface_visible_;
    };
    using SurfaceMeshGfxFactory = GEO::Factory0< SurfaceMeshGfx >;
#define ringmesh_register_surface_gfx(type) \
    geo_register_creator(RINGMesh::SurfaceMeshGfxFactory, type ## Gfx, type::type_name_static())

    class RINGMESH_API VolumeMeshGfx: public MeshEntityGfx {
    public:
        virtual void set_mesh( const VolumeMesh< 3 >& mesh ) = 0;

        static std::unique_ptr< VolumeMeshGfx > create_gfx(
            const VolumeMesh< 3 >& mesh );

        void set_region_visible( bool is_visible )
        {
            region_visible_ = is_visible;
        }
        bool get_region_visible() const
        {
            return region_visible_;
        }

        virtual void draw_volume() = 0;
        virtual void set_draw_cells( GEO::MeshCellType type, bool x ) = 0;
        virtual void set_cell_colors_by_type() = 0;
        virtual void set_cells_color( float red, float green, float blue ) = 0;
        virtual void set_mesh_color( float red, float green, float blue ) = 0;
        virtual void set_mesh_visibility( bool b ) = 0;
        virtual void set_mesh_width( index_t s ) = 0;
        virtual void set_shrink( double s ) = 0;

    protected:
        VolumeMeshGfx()
        {
            set_vertex_visible( false );
            set_region_visible( true );
        }

    private:
        bool region_visible_;
    };
    using VolumeMeshGfxFactory = GEO::Factory0< VolumeMeshGfx >;
#define ringmesh_register_volume_gfx(type) \
    geo_register_creator(RINGMesh::VolumeMeshGfxFactory, type ## Gfx, type::type_name_static())

    class RINGMESH_API AttributeGfxManager {
    ringmesh_disable_copy( AttributeGfxManager );
    public:
        enum Attribute_location {
            polygons, polygon_vertices, cells, cell_vertices, nb_locations
        };
        AttributeGfxManager( GeoModelGfx& gfx );

        GeoModelGfx& gfx()
        {
            return gfx_;
        }

        void compute_range();

        void unbind_attribute();
        void bind_attribute();

        std::string location_name( Attribute_location location );

        void set_maximum( double max )
        {
            maximum_ = max;
        }
        double maximum() const
        {
            return maximum_;
        }
        void set_minimum( double min )
        {
            minimum_ = min;
        }
        double minimum() const
        {
            return minimum_;
        }
        void set_location( Attribute_location location )
        {
            location_ = location;
        }
        Attribute_location location() const
        {
            return location_;
        }
        index_t nb_coordinates() const;
        void set_coordinate( const index_t& coordinate )
        {
            coordinate_ = coordinate;
        }
        const index_t& coordinate() const
        {
            return coordinate_;
        }
        void set_name( const std::string& name )
        {
            name_ = name;
        }
        const std::string& name() const
        {
            return name_;
        }
        void set_colormap( GLuint colormap )
        {
            colormap_texture_ = colormap;
        }
        GLuint colormap() const
        {
            return colormap_texture_;
        }

    private:
        GeoModelGfx& gfx_;

        std::string name_;
        Attribute_location location_;
        index_t coordinate_;
        GLuint colormap_texture_;
        double minimum_;
        double maximum_;

        std::unique_ptr< AttributeGfx > attributes_[nb_locations];

    };

    class RINGMESH_API AttributeGfx {
    public:
        AttributeGfx( AttributeGfxManager& manager )
            : manager_( manager )
        {
        }
        virtual ~AttributeGfx() = default;

        virtual std::string location_name() const = 0;
        void compute_range()
        {
            double attribute_min = max_float64();
            double attribute_max = min_float64();
            do_compute_range( attribute_min, attribute_max );
            manager_.set_minimum( attribute_min );
            manager_.set_maximum( attribute_max );
        }
        virtual void bind_attribute() = 0;
        virtual void unbind_attribute() = 0;
        virtual index_t nb_coordinates() = 0;

    private:
        virtual void do_compute_range(
            double& attribute_min,
            double& attribute_max ) = 0;

    protected:
        AttributeGfxManager& manager_;
    };
}

#endif
