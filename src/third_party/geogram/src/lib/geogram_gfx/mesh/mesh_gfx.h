/*
 *  Copyright (c) 2012-2014, Bruno Levy
 *  All rights reserved.
 *
 *  Redistribution and use in source and binary forms, with or without
 *  modification, are permitted provided that the following conditions are met:
 *
 *  * Redistributions of source code must retain the above copyright notice,
 *  this list of conditions and the following disclaimer.
 *  * Redistributions in binary form must reproduce the above copyright notice,
 *  this list of conditions and the following disclaimer in the documentation
 *  and/or other materials provided with the distribution.
 *  * Neither the name of the ALICE Project-Team nor the names of its
 *  contributors may be used to endorse or promote products derived from this
 *  software without specific prior written permission.
 * 
 *  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 *  AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 *  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 *  ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 *  LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 *  CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 *  SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 *  INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 *  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 *  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 *  POSSIBILITY OF SUCH DAMAGE.
 *
 *  If you modify this software, you should include a notice giving the
 *  name of the person performing the modification, the date of modification,
 *  and the reason for such modification.
 *
 *  Contact: Bruno Levy
 *
 *     Bruno.Levy@inria.fr
 *     http://www.loria.fr/~levy
 *
 *     ALICE Project
 *     LORIA, INRIA Lorraine, 
 *     Campus Scientifique, BP 239
 *     54506 VANDOEUVRE LES NANCY CEDEX 
 *     FRANCE
 *
 */

#ifndef __GEOGRAM_GFX_MESH_GFX__
#define __GEOGRAM_GFX_MESH_GFX__

#include <geogram_gfx/basic/common.h>
#include <geogram/mesh/mesh.h>

/**
 * \file geogram_gfx/mesh/mesh_gfx.h
 * \brief A class to display a mesh using OpenGL.
 */

namespace GEO {

    class MeshGfxImpl;
    
    /**
     * \brief Draws a mesh using OpenGL.
     */
    class GEOGRAM_GFX_API MeshGfx {
    public:

        /**
         * \brief MeshGfx constructor.
         */
        MeshGfx();

        /**
         * \brief MeshGfx destructor.
         */
        ~MeshGfx();

        /**
         * \brief Draws the vertices of the mesh.
         */
        void draw_vertices();

        /**
         * \brief Draws the edges of the mesh.
         */
        void draw_edges();
        
        /**
         * \brief Draws the surfacic part of the mesh.
         */
        void draw_surface();

        /**
         * \brief Draws the borders of the surfacic 
         *  part of the mesh.
         */
        void draw_surface_borders();
        
        /**
         * \brief Draws the volumetric part of the mesh.
         */
        void draw_volume();


        /**
         * \brief Gets the mesh visibility flag.
         * \details The mesh visibility flags specifies
         *  whether mesh edges should be drawn. The used
         *  color can be specified by set_mesh_color()
         * \retval true if mesh edges should be displayed
         * \retval false otherwise
         */
        bool get_show_mesh() const;

        /**
         * \brief Sets the mesh visibility flag.
         * \param[in] x the new value of the mesh visibility flag.
         * \details The mesh visibility flags specifies
         *  whether mesh edges should be drawn. The used
         *  color can be specified by set_mesh_color()
         * \note For now, regions are only implemented for
         *   triangulated meshes and tetrahedralized meshes
         *   (not implemented yet for hybrid surfacic and 
         *    volumetric meshes).
         */
        void set_show_mesh(bool x);

        /**
         * \brief Gets the mesh width
         * \details The mesh width is taken into account 
         *   when the mesh visibility flag is set 
         *   (by set_show_mesh()), when drawing facets
         *   and cells.
         * \return the mesh width
         */
        index_t get_mesh_width() const;

        /**
         * \brief Sets the mesh width
         * \details The mesh width is taken into account 
         *   when the mesh visibility flag is set 
         *   (by set_show_mesh()), when drawing facets
         *   and cells.
         * \param[in] x the mesh width (minimum is 1)
         */
        void set_mesh_width(index_t x);

        /**
         * \brief Gets the mesh border width
         * \details The mesh border width is the one used
         *   by draw_surface_borders()
         * \return the mesh border width
         */
        index_t get_mesh_border_width() const;

        /**
         * \brief Sets the mesh border width
         * \details The mesh border width is the one used
         *   by draw_surface_borders()
         * \param[in] x the mesh width (minimum is 1)
         */
        void set_mesh_border_width(index_t x);
        
        /**
         * \brief Gets the cells shrink coefficient.
         * \details The cells shrink coefficient is used 
         *  to display cells slighly smaller than what they are.
         *  Cells shrinking is only supported in GLSL mode.
         * \return the cells shrink coefficient, betwe 0.0 (no shrink) 
         *  and 1.0 (full shrink)
         */
        double get_shrink() const;

        /**
         * \brief Sets the cells shrink coefficient.
         * \details The cells shrink coefficient is used 
         *  to display cells slighly smaller than what they are.
         *  Cells shrinking is only supported in GLSL mode.
         * \param[in] x the cells shrink coefficient, betwe 0.0 (no shrink) 
         *  and 1.0 (full shrink)
         */
        void set_shrink(double x);


        /**
         * \brief Gets the animate flag
         * \details When animate mode is activated and the
         *  mesh has 6d vertices, then an animation is displayed.
         *  The first three coordinates correspond to the vertex
         *  position at initial time (t=0). The last three coordinates 
         *  correspond to the vertex position at final time (t=1).
         * \retval true if animation is used
         * \retval false otherwise
         * \see get_time(), set_time()
         */
        bool get_animate() const;

        /**
         * \brief Gets the animate flag
         * \details When animate mode is activated and the
         *  mesh has 6d vertices, then an animation is displayed.
         *  The first three coordinates correspond to the vertex
         *  position at initial time (t=0). The last three coordinates 
         *  correspond to the vertex position at final time (t=1).
         * \param[in] x true if animation should be used, false otherwise
         * \see get_time(), set_time()
         */
        void set_animate(bool x);
        
        /**
         * \brief Gets the time of the animation.
         * \details Used if animate mode is set.
         * \return the time of the animation, betwe 0.0 (initial) 
         *  and 1.0 (final)
         * \see get_animate(), set_animate()
         */
        double get_time() const;

        /**
         * \brief Gets the time of the animation.
         * \details Used if animate mode is set.
         * \param[in] x the time of the animation, betwe 0.0 (initial) 
         *  and 1.0 (final)
         * \see get_animate(), set_animate()
         */
        void set_time(double x);
        
        /**
         * \brief Gets the cell visibility flag.
         * \details It is possible to specify cell visibility
         *  flags for each individual cell type.
         * \param[in] type one of MESH_TET, MESH_HEX, MESH_PRISM, MESH_PYRAMID
         * \retval true if the cells of \p type should be displayed
         * \retval false otherwise
         */
        bool get_draw_cells(MeshCellType type) const;


        /**
         * \brief Sets the cell visibility flag.
         * \details It is possible to specify cell visibility
         *  flags for each individual cell type.
         * \param[in] type one of MESH_TET, MESH_HEX, MESH_PRISM, MESH_PYRAMID
         * \param[in] x true if mesh cells of type \p type should be displayed,
         *  false otherwise.
         */
        void set_draw_cells(MeshCellType type, bool x);

        /**
         * \brief Sets the points color
         * \details Specifies the color used to display points
         * \param[in] r,g,b the components of the points color,
         *  in (0.0 .. 1.0)
         * \see draw_points()
         */
        void set_points_color(float r, float g, float b);

        /**
         * \brief Gets the points color
         * \param[out] r,g,b the components of the points color,
         *  in (0.0 .. 1.0)
         * \see draw_points()
         */
        void get_points_color(float& r, float& g, float& b) const;

        /**
         * \brief Sets the point size
         * \param[in] x the point size (minimum 1)
         * \see draw_points()
         */
        void set_points_size(float x);

        /**
         * \brief Gets the point size
         * \return the point size
         * \see draw_points()
         */
        float get_points_size() const;
        
        /**
         * \brief Sets the mesh color
         * \details Specifies the mesh color to be used if 
         *  mesh edges should be displayed.
         * \param[in] r,g,b the components of the mesh color,
         *  in (0.0 .. 1.0)
         * \see set_show_mesh(), draw_surface(), draw_volume()
         */
        void set_mesh_color(float r, float g, float b);

        /**
         * \brief Gets the mesh color
         * \param[out] r,g,b the components of the mesh color,
         *  in (0.0 .. 1.0)
         * \see set_show_mesh(), draw_surface(), draw_volume()
         */
        void get_mesh_color(float& r, float& g, float& b) const;
        
        /**
         * \brief Sets the surface color
         * \details Specifies the color used to display the
         *  surfacic part of the mesh. It specifies the color 
         *  of both frontfacing and backfacing faces.
         * \param[in] r,g,b the components of the surface color,
         *  in (0.0 .. 1.0)
         * \see draw_surface(), set_backface_surface_color()
         */
        void set_surface_color(float r, float g, float b);

        /**
         * \brief Gets the surface color
         * \param[out] r,g,b the components of the surface color,
         *  in (0.0 .. 1.0)
         * \see draw_surface()
         */
        void get_surface_color(float& r, float& g, float& b) const;
        
        /**
         * \brief Sets the surface color for backfacing faces.
         * \details Specifies the color used to display the
         *  backfaces of the surfacic part of the mesh. 
         * \param[in] r,g,b the components of the surface color,
         *  in (0.0 .. 1.0)
         * \see set_show_mesh(), draw_surface(), draw_volume()
         */
        void set_backface_surface_color(float r, float g, float b);

        /**
         * \brief Sets the color used to display mesh cells.
         * \param[in] r,g,b the components of the cells color,
         *  in (0.0 .. 1.0)
         * \see set_cells_colors_by_type(), draw_volume()
         */
        void set_cells_color(float r, float g, float b);

        /**
         * \brief Gets the cells color
         * \param[out] r,g,b the components of the cells color,
         *  in (0.0 .. 1.0)
         * \see set_cells_colors_by_type(), draw_volume()
         */
        void get_cells_color(float& r, float& g, float& b) const;

        /**
         * \brief Sets a different color for each mesh cell type
         * \details it uses the following colors:
         *  - tets: red
         *  - hexes: white
         *  - prisms: green
         *  - pyramids: blue
         * \see set_cells_color(), draw_volume()
         */
        void set_cells_colors_by_type();
        
        /**
         * \brief Gets the lighing flag
         * \retval true if lighting should be used
         * \retval false otherwise
         */
        bool get_lighting() const;

        /**
         * \brief Sets the lighting flag
         * \param[in] x true if lighting should be used, false
         *  otherwise.
         */
        void set_lighting(bool x);

        /**
         * \brief Sets the mesh
         * \param[in] M a pointer to the mesh that should be
         *  displayed.
         */
        void set_mesh(const Mesh* M);

        /**
         * \brief Gets the mesh
         * \return a pointer to the mesh that will be displayed.
         */
        const Mesh* mesh() const;

    protected:

        /**
         * \brief Creates or updates all the used OpenGL 
         *  objects/program/textures.
         * \details May throw an exception if OpenGL driver / hardware 
         *  does not support the used features.
         */
        void setup();

        /**
         * \brief Replaces the current implementation with
         *  another one.
         * \details Colors, attributes, drawing modes are
         *  copied from the current implementation to the
         *  new one. On exit, ownership of new_impl is 
         *  transfered to this MeshGfx.
         */
        void replace_implementation(MeshGfxImpl* new_impl);

        /**
         * \brief Gets the supported GLSL version.
         * \details The supported GLSL version is determined
         *  from hardware/driver capabilities and user-defined
         *  parameters.
         */
        double supported_GLSL_version();
        
        MeshGfxImpl* impl_;
        bool initialized_;
        
        /**
         * \brief Forbids MeshGfx copy..
         */
        MeshGfx(const MeshGfx& rhs);
        
        /**
         * \brief Forbids MeshGfx copy..
         */
        MeshGfx& operator=(const MeshGfx& rhs);
    };

}

#endif
