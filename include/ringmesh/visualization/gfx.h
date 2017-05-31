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

#include <geogram_gfx/glup_viewer/glup_viewer_gui.h>

/*!
 * @file Classes for GeoModel visualization
 * @author Benjamin Chauvin and Arnaud Botella
 */

namespace RINGMesh {
    class GeoModel;
    class GeoModelGfx;
    class AttributeGfx;
    class CornerGfx;
    class LineGfx;
    class SurfaceGfx;
    class RegionGfx;
    class MeshEntityGfx;
}

namespace RINGMesh {

    class RINGMESH_API GeoModelGfxManager {
    ringmesh_disable_copy( GeoModelGfxManager );
    public:
        GeoModelGfxManager( GeoModelGfx& gfx );
        virtual ~GeoModelGfxManager() = default;

        virtual void draw() = 0;
        virtual void initialize() = 0;
        void need_to_update();
        void set_scalar_attribute(
            GEO::MeshElementsFlags subelements,
            const std::string& name,
            double attr_min,
            double attr_max,
            GLuint colormap_texture );
        void unset_scalar_attribute();

        /*!
         * Sets the entity color to all the entities
         * @param[in] r the red component of the color in [0.0, 1.0]
         * @param[in] g the green component of the color in [0.0, 1.0]
         * @param[in] b the blue component of the color in [0.0, 1.0]
         */
        void set_vertex_color( float r, float g, float b ); /*!
         * Sets the vertex color
         * @param[in] e the entity index
         * @param[in] r the red component of the color in [0.0, 1.0]
         * @param[in] g the green component of the color in [0.0, 1.0]
         * @param[in] b the blue component of the color in [0.0, 1.0]
         */
        void set_vertex_color( index_t e, float r, float g, float b );
        /*!
         * Sets the vertex entity visibility to all the entities
         * @param[in] b the visibility
         */
        void set_vertex_visibility( bool b );
        /*!
         * Sets the vertex entity visibility
         * @param[in] e the entity index
         * @param[in] b the visibility
         */
        void set_vertex_visibility( index_t e, bool b );
        /*!
         * Sets the vertex size to all the elements
         * @param[in] s the size
         */
        void set_vertex_size( index_t s );
        /*!
         * Sets the vertex size to all the elements
         * @param[in] e the entity index
         * @param[in] s the size
         */
        void set_vertex_size( index_t e, index_t s );

        /*!
         * Sets the entity color to all the elements
         * @param[in] r the red component of the color in [0.0, 1.0]
         * @param[in] g the green component of the color in [0.0, 1.0]
         * @param[in] b the blue component of the color in [0.0, 1.0]
         */
        void set_mesh_element_color( float r, float g, float b );
        /*!
         * Sets the entity visibility to all the elements
         * @param[in] b the visibility
         */
        void set_mesh_element_visibility( bool b );
        /*!
         * Sets the mesh_element entity size to all the elements
         * @param[in] s the size
         */
        void set_mesh_element_size( index_t s );
        virtual void set_mesh_element_color( index_t e, float r, float g, float b );
        virtual void set_mesh_element_visibility( index_t e, bool b );
        virtual void set_mesh_element_size( index_t e, index_t s );

    protected:
        GeoModelGfx& gfx_;
        std::vector< std::unique_ptr< MeshEntityGfx > > entities_;

    };

    class MeshEntityGfx: public GEO::MeshGfx {
    ringmesh_disable_copy( MeshEntityGfx );
    public:
        MeshEntityGfx(
            const GeoModelGfx& gfx,
            const GEO::Mesh& mesh,
            bool vertice_visible )
            : vertices_visible_( vertice_visible ), gfx_( gfx )
        {
            set_mesh( &mesh );
        }
        virtual ~MeshEntityGfx() = default;

        void need_to_update()
        {
            buffer_objects_dirty_ = true;
            attributes_buffer_objects_dirty_ = true;
        }

        void draw_vertices()
        {
            GEO::MeshGfx::draw_vertices();
        }
        void draw_edges()
        {
            GEO::MeshGfx::draw_edges();
        }

        void set_vertices_visible( bool b )
        {
            vertices_visible_ = b;
        }
        bool get_vertices_visible() const
        {
            return vertices_visible_;
        }

    protected:
        bool vertices_visible_;

        const GeoModelGfx& gfx_;

    };

    class RINGMESH_API CornerGfxManager: public GeoModelGfxManager {
    public:
        CornerGfxManager( GeoModelGfx& gfx );

        /*!
         * Draws the corners
         */
        virtual void draw() override;
        virtual void initialize() override;

        virtual void set_mesh_element_color( index_t c, float r, float g, float b )
            override;
        virtual void set_mesh_element_visibility( index_t c, bool b ) override;
        virtual void set_mesh_element_size( index_t c, index_t s ) override;

    };

    class RINGMESH_API LineGfxManager: public GeoModelGfxManager {
    public:
        LineGfxManager( GeoModelGfx& gfx );

        /*!
         * Draws the lines
         */
        virtual void draw() override;
        virtual void initialize() override;
        /*!
         * Sets the line color
         * @param[in] l the line index
         * @param[in] r the red component of the color in [0.0, 1.0]
         * @param[in] g the green component of the color in [0.0, 1.0]
         * @param[in] b the blue component of the color in [0.0, 1.0]
         */
        virtual void set_mesh_element_color( index_t l, float r, float g, float b )
            override;
        /*!
         * Sets the line visibility
         * @param[in] l the line index
         * @param[in] b the visibility
         */
        virtual void set_mesh_element_visibility( index_t l, bool b ) override;
        /*!
         * Sets the mesh_element line size
         * @param[in] l the line index
         * @param[in] s the size
         */
        virtual void set_mesh_element_size( index_t l, index_t s ) override;

    };

    class RINGMESH_API SurfaceGfxManager: public GeoModelGfxManager {
    public:
        SurfaceGfxManager( GeoModelGfx& gfx );

        /*!
         * Draws the surfaces
         */
        virtual void draw() override;
        virtual void initialize() override;
        /*!
         * Sets the surface color
         * @param[in] s the surface index
         * @param[in] r the red component of the color in [0.0, 1.0]
         * @param[in] g the green component of the color in [0.0, 1.0]
         * @param[in] b the blue component of the color in [0.0, 1.0]
         */
        virtual void set_mesh_element_color( index_t s, float r, float g, float b )
            override;
        /*!
         * Sets the surface visibility
         * @param[in] s the surface index
         * @param[in] b the visibility
         */
        virtual void set_mesh_element_visibility( index_t s, bool b ) override;
        /*!
         * Sets the backface surface color to all the elements
         * @param[in] r the red component of the color in [0.0, 1.0]
         * @param[in] g the green component of the color in [0.0, 1.0]
         * @param[in] b the blue component of the color in [0.0, 1.0]
         */
        void set_backface_surfaces_color( float r, float g, float b );
        /*!
         * Sets the backsurface surface color
         * @param[in] s the surface index
         * @param[in] r the red component of the color in [0.0, 1.0]
         * @param[in] g the green component of the color in [0.0, 1.0]
         * @param[in] b the blue component of the color in [0.0, 1.0]
         */
        void set_backface_surface_color( index_t s, float r, float g, float b );
        /*!
         * Sets the mesh surface color to all the elements
         * @param[in] r the red component of the color in [0.0, 1.0]
         * @param[in] g the green component of the color in [0.0, 1.0]
         * @param[in] b the blue component of the color in [0.0, 1.0]
         */
        void set_mesh_color( float r, float g, float b );
        /*!
         * Sets the mesh surface color
         * @param[in] s the surface index
         * @param[in] r the red component of the color in [0.0, 1.0]
         * @param[in] g the green component of the color in [0.0, 1.0]
         * @param[in] b the blue component of the color in [0.0, 1.0]
         */
        void set_mesh_color( index_t s, float r, float g, float b );
        /*!
         * Sets the mesh surface visibility to all the elements
         * @param[in] b the visibility
         */
        void set_mesh_visibility( bool b );
        /*!
         * Sets the mesh surface visibility
         * @param[in] s the surface index
         * @param[in] b the visibility
         */
        void set_mesh_visibility( index_t s, bool b );
        /*!
         * Sets the mesh surface size to all the elements
         * @param[in] s the size
         */
        void set_mesh_size( index_t s );
        /*!
         * Sets the mesh surface size
         * @param[in] s the surface index
         * @param[in] size the size
         */
        void set_mesh_size( index_t s, index_t size );
    };

    class RINGMESH_API RegionGfxManager: public GeoModelGfxManager {
    public:
        RegionGfxManager( GeoModelGfx& gfx );

        /*!
         * Draws the Regions
         */
        virtual void draw() override;
        virtual void initialize() override;

        /*!
         * Sets the cell region color
         * @param[in] m the region index
         * @param[in] r the red component of the color in [0.0, 1.0]
         * @param[in] g the green component of the color in [0.0, 1.0]
         * @param[in] b the blue component of the color in [0.0, 1.0]
         */
        virtual void set_mesh_element_color( index_t e, float r, float g, float b )
            override;
        /*!
         * Sets the cell region visibility to all the regions
         * @param[in] m the region index
         * @param[in] b the visibility
         */
        virtual void set_mesh_element_visibility( index_t e, bool b ) override;

        /*!
         * Sets the edge color to all the meshes
         * @param[in] r the red component of the color in [0.0, 1.0]
         * @param[in] g the green component of the color in [0.0, 1.0]
         * @param[in] b the blue component of the color in [0.0, 1.0]
         */
        void set_edge_color( float r, float g, float b );
        /*!
         * Sets the edge color
         * @param[in] m the mesh index
         * @param[in] r the red component of the color in [0.0, 1.0]
         * @param[in] g the green component of the color in [0.0, 1.0]
         * @param[in] b the blue component of the color in [0.0, 1.0]
         */
        void set_edge_color( index_t m, float r, float g, float b );
        /*!
         * Sets the edge visibility to all the meshes
         * @param[in] b the visibility
         */
        void set_edge_visibility( bool b );
        /*!
         * Sets the edge visibility
         * @param[in] m the mesh index
         * @param[in] b the visibility
         */
        void set_edge_visibility( index_t m, bool b );
        /*!
         * Sets the edge line size to all the lines
         * @param[in] s the size
         */
        void set_edge_size( index_t s );
        /*!
         * Sets the edge line size
         * @param[in] l the line index
         * @param[in] s the size
         */
        void set_edge_size( index_t l, index_t s );

        /*!
         * Sets the mesh region color to all the regions
         * @param[in] r the red component of the color in [0.0, 1.0]
         * @param[in] g the green component of the color in [0.0, 1.0]
         * @param[in] b the blue component of the color in [0.0, 1.0]
         */
        void set_mesh_color( float r, float g, float b );
        /*!
         * Sets the mesh region color
         * @param[in] m the region index
         * @param[in] r the red component of the color in [0.0, 1.0]
         * @param[in] g the green component of the color in [0.0, 1.0]
         * @param[in] b the blue component of the color in [0.0, 1.0]
         */
        void set_mesh_color( index_t m, float r, float g, float b );
        /*!
         * Sets the mesh region visibility to all the regions
         * @param[in] b the visibility
         */
        void set_mesh_visibility( bool b );
        /*!
         * Sets the mesh region visibility
         * @param[in] m the region index
         * @param[in] b the visibility
         */
        void set_mesh_visibility( index_t m, bool b );
        /*!
         * Sets the mesh region size to all the regions
         * @param[in] s the size
         */
        void set_mesh_size( index_t s );
        /*!
         * Sets the mesh region size
         * @param[in] m the region index
         * @param[in] s the size
         */
        void set_mesh_size( index_t m, index_t s );

        /*!
         * Toggles the cell type to all the regions
         */
        void set_draw_cells( GEO::MeshCellType type, bool x );
        /*!
         * Toggles the cell type display
         */
        void set_draw_cells( index_t m, GEO::MeshCellType type, bool x );
        /*!
         * Toggles the cell region color per cell type to all the regions
         */
        void set_color_cell_type();
        /*!
         * Toggles the cell region color per cell type
         * @param[in] m the region index
         */
        void set_color_cell_type( index_t m );
        void set_cell_type_visibility( GEO::MeshCellType t, bool b );
        void set_cell_type_visibility( index_t m, GEO::MeshCellType t, bool b );
        /*!
         * Sets the cell region shrink to all the regions
         * @param[in] s the shrink
         */
        void set_shrink( double s );
        /*!
         * Sets the cell region shrink
         * @param[in] m the region index
         * @param[in] s the shrink
         */
        void set_shrink( index_t m, double s );

    };

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

    class AttributeGfx {
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

    class RINGMESH_API GeoModelGfx {
    ringmesh_disable_copy( GeoModelGfx );
    public:

        GeoModelGfx();
        ~GeoModelGfx();

        /*!
         * Sets the GeoModel associated to the graphics
         * @param[in] geomodel the GeoModel
         */
        void set_geomodel( const GeoModel& geomodel );
        /*!
         * Gets the GeoModel associated to the graphics
         * @return the GeoModel
         */
        const GeoModel* geomodel() const;
        /*!
         * Initializes the database according the GeoModel dimensions
         */
        void initialize();
        void need_to_update();

    private:
        /// The GeoModel associated to the graphics
        const GeoModel* geomodel_;

    public:
        CornerGfxManager corners;
        LineGfxManager lines;
        SurfaceGfxManager surfaces;
        RegionGfxManager regions;
        AttributeGfxManager attribute;
    };
}

#endif
