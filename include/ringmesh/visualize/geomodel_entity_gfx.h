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

#include <memory>

#include <geogram_gfx/glup_viewer/glup_viewer_gui.h>

/*!
 * @file Classes for GeoModel entity visualization
 * @author Benjamin Chauvin and Arnaud Botella
 */

namespace RINGMesh
{
    FORWARD_DECLARATION_DIMENSION_CLASS( PointSetMeshGfx );
    FORWARD_DECLARATION_DIMENSION_CLASS( LineMeshGfx );
    FORWARD_DECLARATION_DIMENSION_CLASS( SurfaceMeshGfx );
    FORWARD_DECLARATION_DIMENSION_CLASS( VolumeMeshGfx );
    FORWARD_DECLARATION_DIMENSION_CLASS( MeshEntityGfx );
    FORWARD_DECLARATION_DIMENSION_CLASS( GeoModelGfx );

    ALIAS_3D( GeoModelGfx );
} // namespace RINGMesh

namespace RINGMesh
{
    template < index_t DIMENSION >
    class visualize_api GeoModelGfxEntity
    {
        ringmesh_disable_copy_and_move( GeoModelGfxEntity );
        ringmesh_template_assert_2d_or_3d( DIMENSION );

    public:
        explicit GeoModelGfxEntity( GeoModelGfx< DIMENSION >& gfx );
        virtual ~GeoModelGfxEntity() = default;

        virtual void draw() = 0;
        virtual void initialize() = 0;
        void set_scalar_attribute( GEO::MeshElementsFlags subelements,
            const std::string& name,
            double attr_min,
            double attr_max,
            GLuint colormap_texture );
        void unset_scalar_attribute();

        /*!
         * Sets the entity color to all the entities
         * @param[in] red the red component of the color in [0.0, 1.0]
         * @param[in] green the green component of the color in [0.0, 1.0]
         * @param[in] blue the blue component of the color in [0.0, 1.0]
         */
        void set_vertex_color( float red, float green, float blue );
        /*!
         * Sets the vertex color
         * @param[in] entity_id the entity index
         * @param[in] red the red component of the color in [0.0, 1.0]
         * @param[in] green the green component of the color in [0.0, 1.0]
         * @param[in] blue the blue component of the color in [0.0, 1.0]
         */
        void set_vertex_color(
            index_t entity_id, float red, float green, float blue );
        /*!
         * Sets the vertex entity visibility to all the entities
         * @param[in] is_visible the visibility
         */
        void set_vertex_visibility( bool is_visible );
        /*!
         * Sets the vertex entity visibility
         * @param[in] entity_id the entity index
         * @param[in] is_visible the visibility
         */
        void set_vertex_visibility( index_t entity_id, bool is_visible );
        /*!
         * Sets the vertex size to all the elements
         * @param[in] size the size
         */
        void set_vertex_size( index_t size );
        /*!
         * Sets the vertex size to all the elements
         * @param[in] entity_id the entity index
         * @param[in] size the size
         */
        void set_vertex_size( index_t entity_id, index_t size );

    protected:
        GeoModelGfx< DIMENSION >& gfx_;
        std::vector< std::unique_ptr< MeshEntityGfx< DIMENSION > > > entities_;
    };

    template < index_t DIMENSION >
    class visualize_api CornerGfxEntity final
        : public GeoModelGfxEntity< DIMENSION >
    {
    public:
        explicit CornerGfxEntity( GeoModelGfx< DIMENSION >& gfx );

        PointSetMeshGfx< DIMENSION >& corner( index_t corner_id );

        /*!
         * Draws the corners
         */
        void draw() override;
        void initialize() override;
    };

    template < index_t DIMENSION >
    class visualize_api LineGfxEntity final
        : public GeoModelGfxEntity< DIMENSION >
    {
    public:
        explicit LineGfxEntity( GeoModelGfx< DIMENSION >& gfx );

        LineMeshGfx< DIMENSION >& line( index_t line_id );

        /*!
         * Draws the lines
         */
        void draw() override;
        void initialize() override;

        /*!
         * Sets the line color
         * @param[in] red the red component of the color in [0.0, 1.0]
         * @param[in] green the green component of the color in [0.0, 1.0]
         * @param[in] blue the blue component of the color in [0.0, 1.0]
         */
        void set_line_color( float red, float green, float blue );
        /*!
         * Sets the line color
         * @param[in] line_id the line index
         * @param[in] red the red component of the color in [0.0, 1.0]
         * @param[in] green the green component of the color in [0.0, 1.0]
         * @param[in] blue the blue component of the color in [0.0, 1.0]
         */
        void set_line_color(
            index_t line_id, float red, float green, float blue );
        /*!
         * Sets the line visibility
         * @param[in] is_visible the visibility
         */
        void set_line_visibility( bool is_visible );
        /*!
         * Sets the line visibility
         * @param[in] line_id the line index
         * @param[in] is_visible the visibility
         */
        void set_line_visibility( index_t line_id, bool is_visible );
        /*!
         * Sets the line size
         * @param[in] size the size
         */
        void set_line_size( index_t size );
        /*!
         * Sets the line size
         * @param[in] line_id the line index
         * @param[in] size the size
         */
        void set_line_size( index_t line_id, index_t size );
    };

    template < index_t DIMENSION >
    class visualize_api SurfaceGfxEntity final
        : public GeoModelGfxEntity< DIMENSION >
    {
    public:
        explicit SurfaceGfxEntity( GeoModelGfx< DIMENSION >& gfx );

        SurfaceMeshGfx< DIMENSION >& surface( index_t surface_id );

        /*!
         * Draws the surfaces
         */
        void draw() override;
        void initialize() override;
        /*!
         * Sets the surfaces color
         * @param[in] red the red component of the color in [0.0, 1.0]
         * @param[in] green the green component of the color in [0.0, 1.0]
         * @param[in] blue the blue component of the color in [0.0, 1.0]
         */
        void set_surface_color( float red, float green, float blue );
        /*!
         * Sets the surface color
         * @param[in] s the surface index
         * @param[in] red the red component of the color in [0.0, 1.0]
         * @param[in] green the green component of the color in [0.0, 1.0]
         * @param[in] blue the blue component of the color in [0.0, 1.0]
         */
        void set_surface_color(
            index_t surface_id, float red, float green, float blue );
        /*!
         * Sets the surface visibility
         * @param[in] is_visible the visibility
         */
        void set_surface_visibility( bool is_visible );
        /*!
         * Sets the surface visibility
         * @param[in] surface_id the surface index
         * @param[in] is_visible the visibility
         */
        void set_surface_visibility( index_t surface_id, bool is_visible );
        /*!
         * Sets the backface surface color to all the elements
         * @param[in] red the red component of the color in [0.0, 1.0]
         * @param[in] green the green component of the color in [0.0, 1.0]
         * @param[in] blue the blue component of the color in [0.0, 1.0]
         */
        void set_backface_surface_color( float red, float green, float blue );
        /*!
         * Sets the backsurface surface color
         * @param[in] surface_id the surface index
         * @param[in] red the red component of the color in [0.0, 1.0]
         * @param[in] green the green component of the color in [0.0, 1.0]
         * @param[in] blue the blue component of the color in [0.0, 1.0]
         */
        void set_backface_surface_color(
            index_t surface_id, float red, float green, float blue );
        /*!
         * Sets the mesh surface color to all the elements
         * @param[in] red the red component of the color in [0.0, 1.0]
         * @param[in] green the green component of the color in [0.0, 1.0]
         * @param[in] blue the blue component of the color in [0.0, 1.0]
         */
        void set_mesh_color( float red, float green, float blue );
        /*!
         * Sets the mesh surface color
         * @param[in] surface_id the surface index
         * @param[in] red the red component of the color in [0.0, 1.0]
         * @param[in] green the green component of the color in [0.0, 1.0]
         * @param[in] blue the blue component of the color in [0.0, 1.0]
         */
        void set_mesh_color(
            index_t surface_id, float red, float green, float blue );
        /*!
         * Sets the mesh surface visibility to all the elements
         * @param[in] is_visible the visibility
         */
        void set_mesh_visibility( bool is_visible );
        /*!
         * Sets the mesh surface visibility
         * @param[in] surface_id the surface index
         * @param[in] is_visible the visibility
         */
        void set_mesh_visibility( index_t surface_id, bool is_visible );
        /*!
         * Sets the mesh surface size to all the elements
         * @param[in] size the size
         */
        void set_mesh_size( index_t size );
        /*!
         * Sets the mesh surface size
         * @param[in] surface_id the surface index
         * @param[in] size the size
         */
        void set_mesh_size( index_t surface_id, index_t size );
    };

    template < index_t DIMENSION >
    class visualize_api RegionGfxEntity final : public GeoModelGfxEntity< 3 >
    {
    public:
        explicit RegionGfxEntity( GeoModelGfx3D& gfx );

        VolumeMeshGfx< DIMENSION >& region( index_t region_id );

        /*!
         * Draws the Regions
         */
        void draw() override;
        void initialize() override;

        /*!
         * Sets the region color
         * @param[in] red the red component of the color in [0.0, 1.0]
         * @param[in] green the green component of the color in [0.0, 1.0]
         * @param[in] blue the blue component of the color in [0.0, 1.0]
         */
        void set_region_color( float red, float green, float blue );
        /*!
         * Sets the region color
         * @param[in] region_id the region index
         * @param[in] red the red component of the color in [0.0, 1.0]
         * @param[in] green the green component of the color in [0.0, 1.0]
         * @param[in] blue the blue component of the color in [0.0, 1.0]
         */
        void set_region_color(
            index_t region_id, float red, float green, float blue );
        /*!
         * Sets the region visibility to all the regions
         * @param[in] is_visible the visibility
         */
        void set_region_visibility( bool is_visible );
        /*!
         * Sets the region visibility to all the regions
         * @param[in] region_id the region index
         * @param[in] is_visible the visibility
         */
        void set_region_visibility( index_t region_id, bool is_visible );

        /*!
         * Sets the mesh region color to all the regions
         * @param[in] red the red component of the color in [0.0, 1.0]
         * @param[in] green the green component of the color in [0.0, 1.0]
         * @param[in] blue the blue component of the color in [0.0, 1.0]
         */
        void set_mesh_color( float red, float green, float blue );
        /*!
         * Sets the mesh region color
         * @param[in] region_id the region index
         * @param[in] red the red component of the color in [0.0, 1.0]
         * @param[in] green the green component of the color in [0.0, 1.0]
         * @param[in] blue the blue component of the color in [0.0, 1.0]
         */
        void set_mesh_color(
            index_t region_id, float red, float green, float blue );
        /*!
         * Sets the mesh region visibility to all the regions
         * @param[in] is_visible the visibility
         */
        void set_mesh_visibility( bool is_visible );
        /*!
         * Sets the mesh region visibility
         * @param[in] region_id the region index
         * @param[in] is_visible the visibility
         */
        void set_mesh_visibility( index_t region_id, bool is_visible );
        /*!
         * Sets the mesh region size to all the regions
         * @param[in] size the size
         */
        void set_mesh_size( index_t size );
        /*!
         * Sets the mesh region size
         * @param[in] region_id the region index
         * @param[in] size the size
         */
        void set_mesh_size( index_t region_id, index_t size );

        /*!
         * Toggles the cell type to all the regions
         */
        void set_draw_cells( CellType type, bool x );
        /*!
         * Toggles the cell type display
         */
        void set_draw_cells( index_t region_id, CellType type, bool x );
        /*!
         * Toggles the cell region color per cell type to all the regions
         */
        void set_cell_colors_by_type();
        /*!
         * Toggles the cell region color per cell type
         * @param[in] region_id the region index
         */
        void set_cell_colors_by_type( index_t region_id );
        void set_cell_type_visibility( CellType t, bool is_visible );
        void set_cell_type_visibility(
            index_t region_id, CellType t, bool is_visible );
        /*!
         * Sets the cell region shrink to all the regions
         * @param[in] shrink the shrink
         */
        void set_shrink( double shrink );
        /*!
         * Sets the cell region shrink
         * @param[in] region_id the region index
         * @param[in] shrink the shrink
         */
        void set_shrink( index_t region_id, double shrink );
    };

    ALIAS_2D_AND_3D( RegionGfxEntity );
} // namespace RINGMesh

#endif
