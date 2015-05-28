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

#include <geogram_gfx/basic/common.h>
#include <geogram/mesh/mesh.h>
#include <geogram/mesh/mesh_io.h>

/**
 * \file geogram_gfx/mesh/mesh_gfx.h
 * \brief A class to display a mesh using OpenGL.
 */

namespace GEO {

    class Mesh;

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
        bool get_show_mesh() const {
            return show_mesh_;
        }

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
        void set_show_mesh(bool x) {
            show_mesh_ = x;
        }


        /**
         * \brief Gets the mesh width
         * \details The mesh width is taken into account 
         *   when the mesh visibility flag is set 
         *   (by set_show_mesh()), when drawing facets
         *   and cells.
         * \return the mesh width
         */
        index_t get_mesh_width() const {
            return mesh_width_;
        }

        /**
         * \brief Sets the mesh width
         * \details The mesh width is taken into account 
         *   when the mesh visibility flag is set 
         *   (by set_show_mesh()), when drawing facets
         *   and cells.
         * \param[in] x the mesh width (minimum is 1)
         */
        void set_mesh_width(index_t x) {
            mesh_width_ = x;
        }

        /**
         * \brief Gets the mesh border width
         * \details The mesh border width is the one used
         *   by draw_surface_borders()
         * \return the mesh border width
         */
        index_t get_mesh_border_width() const {
            return mesh_border_width_;
        }

        /**
         * \brief Sets the mesh border width
         * \details The mesh border width is the one used
         *   by draw_surface_borders()
         * \param[in] x the mesh width (minimum is 1)
         */
        void set_mesh_border_width(index_t x) {
            mesh_border_width_ = x;
        }
        
        /**
         * \brief Gets the region visibility flag.
         * \return the value of the region visibility flag.
         * \details If activated, cell and facet regions are displayed
         *  using a colormap.
         * \note For now, regions are only implemented for
         *   triangulated meshes and tetrahedralized meshes
         *   (not implemented yet for hybrid surfacic and 
         *    volumetric meshes).
         */
        bool get_show_regions() const {
            return show_regions_;
        }

        /**
         * \brief Gets the region visibility flag.
         * \param[in] x the new value of the region visibility flag.
         * \details If activated, cell and facet regions are displayed
         *  using a colormap.
         */
        void set_show_regions(bool x) {
            show_regions_ = x;
        }
        
        /**
         * \brief Gets the cells shrink coefficient.
         * \details The cells shrink coefficient is used 
         *  to display cells slighly smaller than what they are.
         *  Cells shrinking is only supported in GLSL mode.
         * \return the cells shrink coefficient, betwe 0.0 (no shrink) 
         *  and 1.0 (full shrink)
         */
        double get_shrink() const {
            return shrink_;
        }

        /**
         * \brief Sets the cells shrink coefficient.
         * \details The cells shrink coefficient is used 
         *  to display cells slighly smaller than what they are.
         *  Cells shrinking is only supported in GLSL mode.
         * \param[in] x the cells shrink coefficient, betwe 0.0 (no shrink) 
         *  and 1.0 (full shrink)
         */
        void set_shrink(double x) {
            shrink_ = x;
            geo_clamp(shrink_, 0.0, 1.0);
        }


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
        bool get_animate() const {
            return animate_;
        }

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
        void set_animate(bool x) {
            animate_ = x;
        }
        
        /**
         * \brief Gets the time of the animation.
         * \details Used if animate mode is set.
         * \return the time of the animation, betwe 0.0 (initial) 
         *  and 1.0 (final)
         * \see get_animate(), set_animate()
         */
        double get_time() const {
            return time_;
        }

        /**
         * \brief Gets the time of the animation.
         * \details Used if animate mode is set.
         * \param[in] x the time of the animation, betwe 0.0 (initial) 
         *  and 1.0 (final)
         * \see get_animate(), set_animate()
         */
        void set_time(double x) {
            time_ = x;
            geo_clamp(time_, 0.0, 1.0);
        }
        
        /**
         * \brief Gets the cell visibility flag.
         * \details It is possible to specify cell visibility
         *  flags for each individual cell type.
         * \param[in] type one of MESH_TET, MESH_HEX, MESH_PRISM, MESH_PYRAMID
         * \retval true if the cells of \p type should be displayed
         * \retval false otherwise
         */
        bool get_draw_cells(MeshCellType type) const {
            geo_assert(type < MESH_NB_CELL_TYPES);
            return draw_cells_[type];
        }


        /**
         * \brief Sets the cell visibility flag.
         * \details It is possible to specify cell visibility
         *  flags for each individual cell type.
         * \param[in] type one of MESH_TET, MESH_HEX, MESH_PRISM, MESH_PYRAMID
         * \param[in] x true if mesh cells of type \p type should be displayed,
         *  false otherwise.
         */
        void set_draw_cells(MeshCellType type, bool x) {
            geo_assert(type < MESH_NB_CELL_TYPES);
            draw_cells_[type] = x;
        }


        /**
         * \brief Sets the points color
         * \details Specifies the color used to display points
         * \param[in] r,g,b the components of the points color,
         *  in (0.0 .. 1.0)
         * \see draw_points()
         */
        void set_points_color(float r, float g, float b) {
            set_color(PRG_POINTS, r, g, b);
        }

        /**
         * \brief Gets the points color
         * \param[out] r,g,b the components of the points color,
         *  in (0.0 .. 1.0)
         * \see draw_points()
         */
        void get_points_color(float& r, float& g, float& b) const {
            get_color(PRG_POINTS, r, g, b);
        }


        /**
         * \brief Sets the point size
         * \param[in] x the point size (minimum 1)
         * \see draw_points()
         */
        void set_points_size(float x) {
            points_size_ = x;
        }

        /**
         * \brief Gets the point size
         * \return the point size
         * \see draw_points()
         */
        float get_points_size() const {
            return points_size_;
        }
        
        /**
         * \brief Sets the mesh color
         * \details Specifies the mesh color to be used if 
         *  mesh edges should be displayed.
         * \param[in] r,g,b the components of the mesh color,
         *  in (0.0 .. 1.0)
         * \see set_show_mesh(), draw_surface(), draw_volume()
         */
        void set_mesh_color(float r, float g, float b) {
            mesh_color_[0] = r;
            mesh_color_[1] = g;
            mesh_color_[2] = b;
            mesh_color_[3] = 1.0f;                        
        }


        /**
         * \brief Gets the mesh color
         * \param[out] r,g,b the components of the mesh color,
         *  in (0.0 .. 1.0)
         * \see set_show_mesh(), draw_surface(), draw_volume()
         */
        void get_mesh_color(float& r, float& g, float& b) const {
            r = mesh_color_[0];
            g = mesh_color_[1];
            b = mesh_color_[2];
        }
        
        /**
         * \brief Sets the surface color
         * \details Specifies the color used to display the
         *  surfacic part of the mesh. It specifies the color 
         *  of both frontfacing and backfacing faces.
         * \param[in] r,g,b the components of the surface color,
         *  in (0.0 .. 1.0)
         * \see draw_surface(), set_backface_surface_color()
         */
        void set_surface_color(float r, float g, float b) {
            set_color(PRG_TRI,  r, g, b);
            set_color(PRG_QUAD, r, g, b);            
        }

        /**
         * \brief Gets the surface color
         * \param[out] r,g,b the components of the surface color,
         *  in (0.0 .. 1.0)
         * \see draw_surface()
         */
        void get_surface_color(float& r, float& g, float& b) const {
            get_color(PRG_TRI, r, g, b);
        }
        
        /**
         * \brief Sets the surface color for backfacing faces.
         * \details Specifies the color used to display the
         *  backfaces of the surfacic part of the mesh. 
         * \param[in] r,g,b the components of the surface color,
         *  in (0.0 .. 1.0)
         * \see set_show_mesh(), draw_surface(), draw_volume()
         */
        void set_backface_surface_color(float r, float g, float b) {
            set_back_color(PRG_TRI,  r, g, b);
            set_back_color(PRG_QUAD, r, g, b);            
        }

        /**
         * \brief Sets the color used to display mesh cells.
         * \param[in] r,g,b the components of the cells color,
         *  in (0.0 .. 1.0)
         * \see set_cells_colors_by_type(), draw_volume()
         */
        void set_cells_color(float r, float g, float b) {
            set_color(PRG_TET,     r, g, b);
            set_color(PRG_HEX,     r, g, b);
            set_color(PRG_PRISM,   r, g, b);
            set_color(PRG_PYRAMID, r, g, b);
        }

        /**
         * \brief Gets the cells color
         * \param[out] r,g,b the components of the cells color,
         *  in (0.0 .. 1.0)
         * \see set_cells_colors_by_type(), draw_volume()
         */
        void get_cells_color(float& r, float& g, float& b) const {
            get_color(PRG_TET, r, g, b);
        }

        
        /**
         * \brief Sets a different color for each mesh cell type
         * \details it uses the following colors:
         *  - tets: red
         *  - hexes: white
         *  - prisms: green
         *  - pyramids: blue
         * \see set_cells_color(), draw_volume()
         */
        void set_cells_colors_by_type() {
            set_color(PRG_TET,     1.0f, 0.0f, 0.0f);
            set_color(PRG_HEX,     0.9f, 0.9f, 0.9f);
            set_color(PRG_PRISM,   0.0f, 1.0f, 0.0f);
            set_color(PRG_PYRAMID, 0.0f, 0.0f, 1.0f);            
        }

        /**
         * \brief Gets the lighing flag
         * \retval true if lighting should be used
         * \retval false otherwise
         */
        bool get_lighting() const {
            return lighting_;
        }

        /**
         * \brief Sets the lighting flag
         * \param[in] x true if lighting should be used, false
         *  otherwise.
         */
        void set_lighting(bool x) {
            lighting_ = x;
        }


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
        const Mesh* mesh() const {
            return mesh_;
        }
        
    protected:

        /**
         * \brief Draws a surfacic mesh solely composed 
         *  of triangles
         */
        void draw_triangles();

        /**
         * \brief Draws a surfacic mesh solely composed 
         *   of triangles and quads.
         */
        void draw_triangles_and_quads();

        /**
         * \brief Draws a surfacic mesh composed of 
         *  arbitrary polygons.
         */
        void draw_polygons();


        /**
         * \brief Draws a surfacic mesh animation.
         * \details The mesh is in 6d, the first three coordinates
         *  correspond to vertex location at time t=0; and the last
         *  three coordinates correspond to vertex location at time 
         *  t=1. Positions are interpolated with time_.
         *  \see set_time(), set_animate()
         */
        void draw_triangles_animation();
        

        /**
         * \brief Draws a mesh solely composed of tetrahedra.
         */
        void draw_tets();


        /**
         * \brief Draws a surfacic mesh animation.
         * \details The mesh is in 6d, the first three coordinates
         *  correspond to vertex location at time t=0; and the last
         *  three coordinates correspond to vertex location at time 
         *  t=1. Positions are interpolated with time_.
         *  \see set_time(), set_animate()
         */
        void draw_tets_animation();

        /**
         * \brief Draws a mesh composed of arbitrary cells.
         */
        void draw_cells();

        /**
         * \brief Draws a mesh composed of arbitrary cells, 
         *  without using any GLSL shader.
         */
        void draw_cells_no_shader();

        /**
         * \brief Sends all the cell of a given type to
         *  OpenGL using glDrawElements() call.
         * \details A cache (cell draw cache) memorizes the
         *  first index and length of each continuous chunk
         *  of cells with the same type. 
         * \param[in] cell_type the type of the cells to 
         *  be sent to OpenGL (one of MESH_TET, MESH_HEX,
         *  MESH_PRISM, MESH_PYRAMID)
         * \param[in] mode the OpenGL primitive to be used
         * \see clear_cell_draw_cache()
         */
        void draw_mesh_cells_as_opengl_elements(
            MeshCellType cell_type, GLenum mode
        );


        /**
         * \brief Clears the cell draw cache.
         * \details The cell draw cache memorizes the
         *  first index and length of each continuous chunk
         *  of cells with the same type. 
         * \see draw_mesh_cells_as_opengl_elements()
         */
        void clear_cell_draw_cache();

        /**
         * \brief Sends all the cell of a given type to
         *  OpenGL as points with attributes.
         * \details Each cell generates an OpenGL vertex.
         *  The vertices of the cell are generic attributes
         *  of the generated OpenGL vertices.
         * \param[in] cell_type the type of the cells to 
         *  be sent to OpenGL (one of MESH_TET, MESH_HEX,
         *  MESH_PRISM, MESH_PYRAMID)
         */
        void draw_mesh_cells_as_opengl_points(
            MeshCellType cell_type
        );
        
        
        /**
         * \brief Creates shaders and Vertex Buffer Arrays,
         *  and binds Vertex Buffer Objects depending on
         *  what should be drawn after.
         * \param[in] what specifies what should be drawn,
         *   one of MESH_VERTICES, MESH_FACETS, MESH_CELLS
         */
        void begin_draw(MeshElementsFlags what);

        /**
         * \brief Unbinds Vertex Buffer Objects.
         */
        void end_draw();


        /**
         * \brief Symbolic constants referring to a GPU program in the
         *  array programs[].
         */
        enum ShaderName {
            PRG_POINTS        =0,
            PRG_TRI           =1,
            PRG_QUAD          =2,
            PRG_TET           =3,
            PRG_HEX           =4,
            PRG_PRISM         =5,
            PRG_PYRAMID       =6,
            PRG_NB            =7
        } ;

        /**
         * \brief Starts using the shader of the specified
         *  name.
         * \details In addition, the front and back colors attached
         *  to the shader are sent to OpenGL using set_colors()
         * \param name one of (PRG_POINTS, PRG_TRI, PRG_QUAD,
         *  PRG_TET, PRG_HEX, PRG_PRISM, PRG_PYRAMID)
         */
        void begin_shader(ShaderName name);

        /**
         * \brief Stops using the shader of the specified name
         */
        void end_shader();

        /**
         * \brief Sends to OpenGL the front and back colors 
         *  attached to a given shader.
         */
        void set_colors(ShaderName name);
        
        /**
         * \brief Creates all the shaders.
         */
        void setup_shaders();

        /**
         * \brief Deletes all the shaders.
         */
        void delete_shaders();
        
        /**
         * \brief Creates the vertex buffer objects.
         */
        void setup_VBOs();

        /**
         * \brief Deletes all the vertex buffer objects.
         */
        void delete_VBOs();

        
    protected:
        const Mesh* mesh_;
        
        GLuint     vertices_VBO_;
        GLuint     edge_indices_VBO_;
        GLuint     facet_indices_VBO_;
        GLuint     cell_indices_VBO_;
        GLuint     facet_region_VBO_;
        GLuint     cell_region_VBO_;

        GLuint     colormap_TEX_;

        bool VBO_dirty_;
        
        bool show_mesh_;
        bool show_regions_;
        double shrink_;
        bool animate_;
        double time_;
        bool lighting_;

        /** 
         * \brief true if the surface has only triangles and quads.
         */
        bool triangles_and_quads_;

        /**
         * \brief true if shaders are used.
         */
        bool GLSL_mode_;

        /**
         * \brief true if tesselation shaders are used.
         */
        bool GLSL_tesselation_;

        /**
         * \brief GLSL version supported by the OpenGL driver.
         */
        double GLSL_version_;
        
        /**
         * \brief true if shaders are already initialized.
         */
        bool shaders_init_;
        
        /**
         * \brief GPU programs are grouped in an array, so that
         *  setting the same uniform variable in all of them can
         *  be done easily with a for() loop.
         */
        GLuint programs_[PRG_NB];
        
        /**
         * \brief Default frontfacing color to be used 
         *  with a given program.
         */
        GLfloat colors_[PRG_NB][4];

        /**
         * \brief Default backfacing color to be used 
         *  with a given program.
         */
        GLfloat back_colors_[PRG_NB][4];
        
        /**
         * \brief Mesh color
         */
        GLfloat mesh_color_[4];
        
        /**
         * \brief Toggles cell drawing by type.
         */
        bool draw_cells_[MESH_NB_CELL_TYPES];
        
        /**
         * \brief Drawing instructions for cells can be 'cached'. 
         * \details They are cached for each cell type independantly.
         *   This array indicates for each cell type whether the cache
         *   is up to date.
         */
        bool cell_draw_cached_[MESH_NB_CELL_TYPES];

        /**
         * \brief Drawing instructions for cells can be 'cached'. 
         * \details They are cached for each cell type independantly.
         *   This array indicates for each drawing call of cell type 
         *  how many vertices should be issued in the call.
         */
        vector<GLsizei> cell_draw_nb_vertices_[GEO::MESH_NB_CELL_TYPES];

        /**
         * \brief Drawing instructions for cells can be 'cached'. 
         * \details They are cached for each cell type independantly.
         *   This array indicates for each drawing call of cell type 
         *  the first index of the call.
         */
        vector<void*>   cell_draw_start_index_[GEO::MESH_NB_CELL_TYPES];

        float points_size_;
        index_t mesh_width_;
        index_t mesh_border_width_;
        
        /**
         * \brief Defines the default color for one of the programs.
         * \param[in] index index of the program, in 0..PRG_NB - 1
         * \param[in] r the red component, in 0.0f..1.0f
         * \param[in] g the green component, in 0.0f..1.0f
         * \param[in] b the blue component, in 0.0f..1.0f
         */
        inline void set_front_color(GLuint index, float r, float g, float b) {
            colors_[index][0] = r;
            colors_[index][1] = g;
            colors_[index][2] = b;
            colors_[index][3] = 1.0f;
        }

        /**
         * \brief Defines the backfacing default color for one 
         *  of the programs.
         * \param[in] index index of the program, in 0..PRG_NB - 1
         * \param[in] r the red component, in 0.0f..1.0f
         * \param[in] g the green component, in 0.0f..1.0f
         * \param[in] b the blue component, in 0.0f..1.0f
         */
        inline void set_back_color(GLuint index, float r, float g, float b) {
            back_colors_[index][0] = r;
            back_colors_[index][1] = g;
            back_colors_[index][2] = b;
            back_colors_[index][3] = 1.0f;
        }

        /**
         * \brief Sets both frontfacing and backfacing colors
         *  for a given program.
         * \param[in] index index of the program, in 0..PRG_NB - 1
         * \param[in] r the red component, in 0.0f..1.0f
         * \param[in] g the green component, in 0.0f..1.0f
         * \param[in] b the blue component, in 0.0f..1.0f
         */
        inline void set_color(GLuint index, float r, float g, float b) {
            set_front_color(index, r, g, b);
            set_back_color(index, r, g, b);
        }

        /**
         * \brief Gets the color used by a given program.
         * \param[in] index index of the program, in 0..PRG_NB - 1
         * \param[out] r the red component, in 0.0f..1.0f
         * \param[out] g the green component, in 0.0f..1.0f
         * \param[out] b the blue component, in 0.0f..1.0f
         */
        inline void get_color(
            GLuint index, float& r, float& g, float& b
        ) const {
            r = colors_[index][0];
            g = colors_[index][1];
            b = colors_[index][2];
        }
    };

    
}
