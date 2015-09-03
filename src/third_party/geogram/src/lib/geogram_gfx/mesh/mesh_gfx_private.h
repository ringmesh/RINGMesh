
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

#ifndef __GEOGRAM_GFX_MESH_GFX_PRIVATE__
#define __GEOGRAM_GFX_MESH_GFX_PRIVATE__

/**
 * \file geogram_gfx/mesh/mesh_gfx_private.h
 * \brief Internal classes used by the implementation of MeshGfx
 */

#include <geogram_gfx/basic/common.h>
#include <geogram_gfx/basic/GLSL.h>
#include <geogram/mesh/mesh.h>
#include <geogram/mesh/mesh_io.h>


namespace GEO {


    /**
     * \brief Baseclass for internal Implementations of MeshGfx.
     */
    class GEOGRAM_GFX_API MeshGfxImpl {
    public:

        /**
         * \brief MeshGfxImpl constructor.
         */
        MeshGfxImpl();
        
        /**
         * \brief MeshGfxImpl destructor.
         */
        virtual ~MeshGfxImpl() ;


        /**
         * \brief Draws the vertices of the mesh.
         */
        virtual void draw_vertices() = 0;

        /**
         * \brief Draws the edges of the mesh.
         */
        virtual void draw_edges() = 0;
        
        /**
         * \brief Draws the surfacic part of the mesh.
         */
        virtual void draw_surface() = 0;

        /**
         * \brief Draws the borders of the surfacic 
         *  part of the mesh.
         */
        virtual void draw_surface_borders() = 0;
        
        /**
         * \brief Draws the volumetric part of the mesh.
         */
        virtual void draw_volume() = 0;


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
            set_color(PRG_LINES, r, g, b);
        }


        /**
         * \brief Gets the mesh color
         * \param[out] r,g,b the components of the mesh color,
         *  in (0.0 .. 1.0)
         * \see set_show_mesh(), draw_surface(), draw_volume()
         */
        void get_mesh_color(float& r, float& g, float& b) const {
            get_color(PRG_LINES, r, g, b);
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
        virtual void set_mesh(const Mesh* M);

        /**
         * \brief Gets the mesh
         * \return a pointer to the mesh that will be displayed.
         */
        const Mesh* mesh() const {
            return mesh_;
        }

        /**
         * \brief Creates OpenGL buffers and shaders
         *  if need be.
         * \details May throw an exception if some OpenGL functionalities
         *   are not supported by the hardware/driver. 
         *   Must be called before using any draw_xxx() function.
         */
        virtual void setup();


        /**
         * \brief Copies all the drawing attributes from
         *  another MeshGfxImpl.
         * \param[in] rhs a const reference to the MeshGfxImpl
         *  from which attributes should be copied
         */
        void copy_drawing_attributes(const MeshGfxImpl& rhs);
        
    protected:

        /**
         * \brief Defines the default color for one of the programs.
         * \param[in] index index of the program, in 0..PRG_NB - 1
         * \param[in] r the red component, in 0.0f..1.0f
         * \param[in] g the green component, in 0.0f..1.0f
         * \param[in] b the blue component, in 0.0f..1.0f
         */
        void set_front_color(GLuint index, float r, float g, float b) {
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
        void set_back_color(GLuint index, float r, float g, float b) {
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
        void set_color(GLuint index, float r, float g, float b) {
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
        void get_color(
            GLuint index, float& r, float& g, float& b
        ) const {
            r = colors_[index][0];
            g = colors_[index][1];
            b = colors_[index][2];
        }

        /**
         * \brief Gets the back color used by a given program.
         * \param[in] index index of the program, in 0..PRG_NB - 1
         * \param[out] r the red component, in 0.0f..1.0f
         * \param[out] g the green component, in 0.0f..1.0f
         * \param[out] b the blue component, in 0.0f..1.0f
         */
        void get_back_color(
            GLuint index, float& r, float& g, float& b
        ) const {
            r = back_colors_[index][0];
            g = back_colors_[index][1];
            b = back_colors_[index][2];
        }
        
        /**
         * \brief Creates shaders and Vertex Buffer Arrays,
         *  and binds Vertex Buffer Objects depending on
         *  what should be drawn after.
         * \param[in] what specifies what should be drawn,
         *   one of MESH_VERTICES, MESH_FACETS, MESH_CELLS
         */
        virtual void begin_draw(MeshElementsFlags what);

        /**
         * \brief Unbinds Vertex Buffer Objects.
         */
        virtual void end_draw();

        /**
         * \brief Symbolic constants referring to a GPU program in the
         *  array programs[].
         */
        enum ShaderName {
            PRG_POINTS  = 0,
            PRG_LINES   = 1,
            PRG_TRI     = 2,
            PRG_QUAD    = 3,
            PRG_TET     = 4,
            PRG_HEX     = 5,
            PRG_PRISM   = 6,
            PRG_PYRAMID = 7,
            PRG_NB      = 8
        } ;

        /**
         * \brief Starts using the shader of the specified
         *  name.
         * \details In addition, the front and back colors attached
         *  to the shader are sent to OpenGL using set_colors()
         * \param name one of (PRG_POINTS, PRG_LINES, PRG_TRI, PRG_QUAD,
         *  PRG_TET, PRG_HEX, PRG_PRISM, PRG_PYRAMID)
         */
        virtual void begin_shader(ShaderName name);

        /**
         * \brief Stops using the shader of the specified name
         */
        virtual void end_shader();

        /**
         * \brief Sends to OpenGL the front and back colors 
         *  attached to a given shader.
         */
        void set_colors(ShaderName name);
        
        /**
         * \brief Creates all the shaders.
         */
        virtual void setup_shaders();

        /**
         * \brief Deletes all the shaders.
         */
        virtual void delete_shaders();
        
        /**
         * \brief Creates the vertex buffer objects.
         */
        virtual void setup_VBOs();

        /**
         * \brief Deletes all the vertex buffer objects.
         */
        virtual void delete_VBOs();


        /**
         * \brief Updates the content of an OpenGL buffer object, 
         *   and resizes it if need be.
         * \param[in,out] buffer_id OpenGL opaque id of the buffer object. 
         *   0 means uninitialized.
         *   may be changed on exit if the buffer needed to be resized.
         * \param[in] target buffer object target 
         *   (GL_ARRAY_BUFFER, GL_INDEX_BUFFER ...)
         * \param[in] new size of the buffer data, in bytes
         * \param[in] data pointer to the data to be copied into the buffer, 
         *  of length new_size
         */
        static void update_buffer_object(
            GLuint& buffer_id, GLenum target, size_t new_size, const void* data
        );


        /**
         * \brief Tests whether current OpenGL polygon mode is filled.
         * \retval true if OpenGL polygon mode is filled
         * \retval false if OpenGL polygon mode is wireframe
         */
        static bool glFillsPolygons() {
            GLint polymode[2];
            glGetIntegerv(GL_POLYGON_MODE, polymode);
            return (polymode[0] == GL_FILL);
        }

        /**
         * \brief Tests whether an OpenGL program is currently used.
         * \retval true if an OpenGL program is currently used
         * \retval false otherwise
         */
        static bool glUsesProgram() {
            GLint program;
            glGetIntegerv(GL_CURRENT_PROGRAM, &program);
            return (program != 0);
        }

        
        /**
         * \brief Sends a vertex of a mesh to OpenGL
         * \param[in] M a const reference to the mesh
         * \param[in] v the index of the vertex in \p M
         */
        void glMeshVertex(index_t v) {
            // TODO: test glDrawElements() with a single
            // element instead.
            if(mesh_->vertices.single_precision()) {
                glVertex3fv(
                    mesh_->vertices.single_precision_point_ptr(v)
                );
            } else {
                glVertex3dv(
                    mesh_->vertices.point_ptr(v)
                );
            }
        }
        

        /**
         * \brief Sends a triangle normal to OpenGL
         * \param[in] pi,pj,pk pointers to the single-precision coordinates
         *  of the free vertices.
         */
        static inline void glTriangleNormal(
            const float* pi, const float* pj, const float* pk
        ) {
            float N[3];
            N[0] = -(pi[1]-pj[1])*(pk[2]-pj[2]) + (pi[2]-pj[2])*(pk[1]-pj[1]);
            N[1] = -(pi[2]-pj[2])*(pk[0]-pj[0]) + (pi[0]-pj[0])*(pk[2]-pj[2]);
            N[2] = -(pi[0]-pj[0])*(pk[1]-pj[1]) + (pi[1]-pj[1])*(pk[0]-pj[0]);
            glNormal3fv(N);
        }

        /**
         * \brief Sends a triangle normal to OpenGL
         * \param[in] pi,pj,pk pointers to the double-precision coordinates
         *  of the free vertices.
         */
        static inline void glTriangleNormal(
            const double* pi, const double* pj, const double* pk
        ) {
            double N[3];
            N[0] = -(pi[1]-pj[1])*(pk[2]-pj[2]) + (pi[2]-pj[2])*(pk[1]-pj[1]);
            N[1] = -(pi[2]-pj[2])*(pk[0]-pj[0]) + (pi[0]-pj[0])*(pk[2]-pj[2]);
            N[2] = -(pi[0]-pj[0])*(pk[1]-pj[1]) + (pi[1]-pj[1])*(pk[0]-pj[0]);
            glNormal3dv(N);
        }


        /**
         * \brief Sends a mesh triangle normal to OpenGL
         * \param[in] i,j,k indices of the vertices of the triangle
         */
        inline void glMeshTriangleNormal(index_t i, index_t j, index_t k) {
            if(mesh_->vertices.single_precision()) {
                const float* pi = mesh_->vertices.single_precision_point_ptr(i);
                const float* pj = mesh_->vertices.single_precision_point_ptr(j);
                const float* pk = mesh_->vertices.single_precision_point_ptr(k);
                glTriangleNormal(pi,pj,pk);
            } else {
                const double* pi = mesh_->vertices.point_ptr(i);
                const double* pj = mesh_->vertices.point_ptr(j);
                const double* pk = mesh_->vertices.point_ptr(k);
                glTriangleNormal(pi,pj,pk);
            }
        }

        /**
         * \brief Sends a mesh surface facet normal to OpenGL
         * \param[in] f index of the facet 
         */
        void glMeshFacetNormal(index_t f);

        /**
         * \brief Sends a mesh cell facet normal to OpenGL
         * \param[in] f index of the facet 
         */
        inline void glMeshCellFacetNormal(index_t c, index_t lf) {
            glMeshTriangleNormal(
                mesh_->cells.facet_vertex(c,lf,0),
                mesh_->cells.facet_vertex(c,lf,1),
                mesh_->cells.facet_vertex(c,lf,2)
            );
        }

        void draw_triangles();
        void draw_triangles_and_quads();
        void draw_polygons();        
        void draw_triangles_animation();
        void draw_tets_as_lines_adjacency();
        void draw_tets_animation_as_lines_adjacency();
        void draw_mesh_cells_as_opengl_elements(
            MeshCellType cell_type, GLenum mode
        );
        void draw_mesh_cells_as_opengl_points(
            MeshCellType cell_type
        );

        void clear_cell_draw_cache();
        
    protected:
        const Mesh* mesh_;
        
        GLuint     vertices_VBO_;
        GLuint     edge_indices_VBO_;
        GLuint     facet_indices_VBO_;
        GLuint     cell_indices_VBO_;

        bool VBO_dirty_;
        
        bool show_mesh_;
        double shrink_;
        bool animate_;
        double time_;
        bool lighting_;

        /** 
         * \brief true if the surface has only triangles and quads.
         */
        bool triangles_and_quads_;

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
         * \brief Toggles cell drawing by type.
         */
        bool draw_cells_[MESH_NB_CELL_TYPES];
        
        float points_size_;
        index_t mesh_width_;
        index_t mesh_border_width_;

        bool shaders_init_;

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
    };

    /************************************************************************************/

    /**
     * \brief Internal implementation of MeshGfx that
     *  does not use GLSL shaders (plain OpenGL)
     */
    class MeshGfxImplNoShader : public MeshGfxImpl {
    public:

        /**
         * \brief MeshGfxImplNoShader constructor.
         */
        MeshGfxImplNoShader();

        /**
         * \brief MeshGfxImplNoShader destructor.
         */
        virtual ~MeshGfxImplNoShader();
        
        virtual void draw_vertices();
        virtual void draw_edges();
        virtual void draw_surface();
        virtual void draw_surface_borders();
        virtual void draw_volume();
    };
    
    /************************************************************************************/    

    /**
     * \brief Internal implementation of MeshGfx that
     *   uses GLSL 1.5
     */
    class MeshGfxImplGLSL150 : public MeshGfxImpl {
    public:

        /**
         * \brief MeshGfxImplGLSL150 constructor.
         */
        MeshGfxImplGLSL150();

        /**
         * \brief MeshGfxImplGLSL150 destructor.
         */
        virtual ~MeshGfxImplGLSL150();
        
        virtual void draw_vertices();
        virtual void draw_edges();
        virtual void draw_surface();
        virtual void draw_surface_borders();
        virtual void draw_volume();

    protected:
        virtual void setup_shaders();
        virtual void begin_shader(ShaderName name);        
    };
    
    /************************************************************************************/    

    /**
     * \brief Internal implementation of MeshGfx that
     *   uses GLSL 4.4
     */
    class MeshGfxImplGLSL440 : public MeshGfxImpl {
    public:

        /**
         * \brief MeshGfxImplGLSL440 constructor.
         */
        MeshGfxImplGLSL440();

        /**
         * \brief MeshGfxImplGLSL150 destructor.
         */
        virtual ~MeshGfxImplGLSL440();
        
        virtual void draw_vertices();
        virtual void draw_edges();
        virtual void draw_surface();
        virtual void draw_surface_borders();
        virtual void draw_volume();

    protected:
        virtual void setup_shaders();
        virtual void begin_shader(ShaderName name);        
    };

    /************************************************************************************/        
}

#endif
