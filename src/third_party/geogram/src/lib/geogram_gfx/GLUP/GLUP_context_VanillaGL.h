/*
 *  Copyright (c) 2012-2016, Bruno Levy
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

#ifndef __GEOGRAM_GFX_GLUP_GLUP_CONTEXT_VANILLAGL___
#define __GEOGRAM_GFX_GLUP_GLUP_CONTEXT_VANILLAGL___

#include <geogram_gfx/basic/common.h>
#include <geogram_gfx/GLUP/GLUP.h>
#include <geogram_gfx/GLUP/GLUP_context.h>

/**
 * \file geogram_gfx/GLUP/GLUP_context_VanillaGL.h
 * \brief Internal implementation of GLUP using plain old OpenGL.
 */

namespace GLUP {
    using namespace GEO;

    /*********************************************************************/

    /**
     * \brief Implementation of GLUP using Vanilla (old-style) OpenGL.
     * \details This implementation does not use any shader. It is used
     *  as a fallback when the initialization of the other ones fails.
     *  Some primitive may be not implemented, degraded or of very low
     *  performance.
     */
    class Context_VanillaGL : public Context {
    public:
        /**
         * \brief Context_VanillaGL constructor.
         */
        Context_VanillaGL();

        /**
         * \copydoc Context::profile_name()
         */
        virtual const char* profile_name() const;
        
        /**
         * \copydoc Context::primitive_supports_array_mode()
         */
        virtual bool primitive_supports_array_mode(GLUPprimitive prim) const;

        /**
         * \copydoc Context::begin()
         */
        virtual void begin(GLUPprimitive primitive);

        /**
         * \copydoc Context::end()
         */
        virtual void end();
        
    protected:

        /**
         * \brief Configures texturing-related OpenGL state
         *  variables according to the GLUP state variables.
         * \details This function is called by begin()
         */
        void configure_OpenGL_texturing();

        /**
         * \brief Configures lighting-related OpenGL state
         *  variables according to the GLUP state variables.
         * \details This function is called by begin()
         */
        void configure_OpenGL_lighting();

        /**
         * \brief Configures clipping-related OpenGL state
         *  variables according to the GLUP state variables.
         * \details This function is called by begin(). 
         */
        void configure_OpenGL_clipping();
        
        /**
         * \brief Configures lighting-related OpenGL state
         *  variables according to the GLUP state variables.
         * \details This function is called by begin(). It needs
         *  to be called after configure_OpenGL_texturing() and
         *  configure_OpenGL_lighting() since it overrides texturing and
         *  lighting settings.
         */
        void configure_OpenGL_picking();        

        
        /**
         * \copydoc Context::setup()
         */
        virtual void setup();
        
        /**
         * \copydoc Context::flush_immediate_buffers()
         */
        virtual void flush_immediate_buffers();


        /**
         * \brief Flushes the immediate buffer with the
         *  current drawing modes. 
         * \details This function is separated from 
         *  flush_immediate_buffers(), since we need to
         *  flush the buffer twice when mesh drawing is
         *  enabled.
         */
        virtual void flush_immediate_buffers_once();
        
        /**
         * \copydoc Context::setup_immediate_buffers()
         */
        virtual void setup_immediate_buffers();


        /**
         * \copydoc Context::setup_primitives()
         */
        virtual void setup_primitives();
        
        /**
         * \copydoc Context::get_state_variable_address()
         */
        Memory::pointer get_state_variable_address(const char* name);


        /**
         * \brief Shrinks the cells in the immediate buffer.
         * \details Applies the shrinking factor (state variable
         *   "cells_shrink") to all the cells stored in the current
         *   immediate buffer. Since there is no function to query
         *   the content of the current buffer, modidying it is 
         *   acceptable.
         */
        void shrink_cells_in_immediate_buffers();

        /**
         * \brief Updates v_is_visible_[] according to
         *  current clipping plane.
         */
        void classify_vertices_in_immediate_buffers();

        /**
         * \brief Tests whether the cell starting at a given vertex
         *  in the immediate buffer is clipped, according to current
         *  clipping mode and current primitive type.
         * \param[in] first_v index of the first vertex of the cell in
         *  the immediate buffer
         * \retval true if the cell starting at \p first_v in the 
         *  immediate buffer is clipped-out
         * \retval false otherwise
         */
        bool cell_is_clipped(index_t first_v);

        /**
         * \brief Tests whether cells should be sliced.
         * \retval true if cells should be sliced
         * \retval false otherwise
         */
        bool clip_slice_cells() const {
            return (
                uniform_state_.toggle[GLUP_CLIPPING].get() &&
                uniform_state_.clipping_mode.get() == GLUP_CLIP_SLICE_CELLS &&
                immediate_state_.primitive() >= GLUP_TETRAHEDRA
            ) ;
        }
        
        /**
         * \brief Sends a vertex and its optional attributes to OpenGL.
         * \param[in] v the index of the vertex from the immediate buffer.
         */
        void output_vertex(index_t v) {
            if(immediate_state_.buffer[GLUP_COLOR_ATTRIBUTE].is_enabled()) {
                glColor4fv(
                    immediate_state_.buffer[GLUP_COLOR_ATTRIBUTE].element_ptr(v)
                );
            }
            if(immediate_state_.buffer[GLUP_TEX_COORD_ATTRIBUTE].is_enabled()) {
                glTexCoord4fv(
                    immediate_state_.buffer[GLUP_TEX_COORD_ATTRIBUTE].
                    element_ptr(v)
                );
            }
            glVertex4fv(
                immediate_state_.buffer[GLUP_VERTEX_ATTRIBUTE].element_ptr(v)
            );
        }

        /**
         * \brief Sends a triangle normal to OpenGL
         * \param[in] v1,v2,v3 the indices of the three vertices from
         *  the immediate buffer.
         */
        void output_normal(index_t v1, index_t v2, index_t v3);

        /**
         * \brief Sends a quad normal to OpenGL
         * \param[in] v1,v2,v3,v4 the indices of the four vertices from
         *  the immediate buffer.
         */
        void output_normal(index_t v1, index_t v2, index_t v3, index_t v4);


        /**
         * \brief Sends a picking id to OpenGL and encodes it as a color.
         * \details The current base picking id is added to the id.
         *  If picking is deactivated or constant by object, 
         *  it does nothing.
         */
        void output_picking_id(index_t id) {
            if(pick_primitives_) {
                glPickingIdAsColor(
                    index_t(uniform_state_.base_picking_id.get()) + id
                );
            }
        }
        
        /**
         * \brief Sends a flat-shaded triangle to OpenGL
         * \param[in] v1,v2,v3 the indices of the three vertices from the
         *  immediate buffer.
         */
        void flat_shaded_triangle(index_t v1, index_t v2, index_t v3) {
            if(uniform_state_.toggle[GLUP_LIGHTING].get()) {
                output_normal(v1,v2,v3);
            }
            output_vertex(v1);
            output_vertex(v2);
            output_vertex(v3);
        }

        /**
         * \brief Sends a flat-shaded quad to OpenGL
         * \param[in] v1,v2,v3,v4 the indices of the three vertices from the
         *  immediate buffer.
         */
        void flat_shaded_quad(index_t v1, index_t v2, index_t v3, index_t v4) {
            if(uniform_state_.toggle[GLUP_LIGHTING].get()) {
                output_normal(v1,v2,v3,v4);
            }
            output_vertex(v1);
            output_vertex(v2);
            output_vertex(v4);
            output_vertex(v3);            
        }

        /**
         * \brief Sends the contents of the immediate buffers to 
         *  OpenGL, as point primitives.
         */
        void draw_immediate_buffer_GLUP_POINTS();

        /**
         * \brief Sends the contents of the immediate buffers to 
         *  OpenGL, as line primitives.
         */
        void draw_immediate_buffer_GLUP_LINES();

        /**
         * \brief Sends the contents of the immediate buffers to 
         *  OpenGL, as triangle primitives.
         */
        void draw_immediate_buffer_GLUP_TRIANGLES();

        /**
         * \brief Sends the contents of the immediate buffers to 
         *  OpenGL, as quad primitives.
         */
        void draw_immediate_buffer_GLUP_QUADS();

        /**
         * \brief Sends the contents of the immediate buffers to 
         *  OpenGL, as tetrahedra primitives.
         */
        void draw_immediate_buffer_GLUP_TETRAHEDRA();

        /**
         * \brief Sends the contents of the immediate buffers to 
         *  OpenGL, as hexahedra primitives.
         */
        void draw_immediate_buffer_GLUP_HEXAHEDRA();

        /**
         * \brief Sends the contents of the immediate buffers to 
         *  OpenGL, as prism primitives.
         */
        void draw_immediate_buffer_GLUP_PRISMS();

        /**
         * \brief Sends the contents of the immediate buffers to 
         *  OpenGL, as pyramid primitives.
         */
        void draw_immediate_buffer_GLUP_PYRAMIDS();

        /**
         * \brief Sends the contents of the immediate buffers to 
         *  OpenGL, as connectors primitives.
         */
        void draw_immediate_buffer_GLUP_CONNECTORS();
        
        /**
         * \brief Computes the intersection between the clipping plane and
         *  a segment.
         * \param[in] v1 index of the first extremity of the segment in the
         *  immediate buffer
         * \param[in] v2 index of the second extremity of the segment in the
         *  immediate buffer
         * \param[in] v1 index of where to wrote the intersection in the 
         *  isect_xxx arrays
         */
        void compute_intersection(index_t v1, index_t v2, index_t vi) {
            const GLUPfloat* eqn = world_clip_plane_;
            const GLUPfloat* p1 = immediate_state_.buffer[0].element_ptr(v1);
            const GLUPfloat* p2 = immediate_state_.buffer[0].element_ptr(v2);
            
            GLUPfloat t = -eqn[3] -(
                eqn[0]*p1[0] +
                eqn[1]*p1[1] +
                eqn[2]*p1[2]
            );

            GLUPfloat d =
                eqn[0]*(p2[0]-p1[0]) +
                eqn[1]*(p2[1]-p1[1]) +
                eqn[2]*(p2[2]-p1[2]) ;
            
            if(fabs(d) < 1e-6) {
                t = 0.5f;
            } else {
                t /= d;
            }

            GLUPfloat s = 1.0f - t;
            
            isect_point_[4*vi+0] = s*p1[0] + t*p2[0];
            isect_point_[4*vi+1] = s*p1[1] + t*p2[1];
            isect_point_[4*vi+2] = s*p1[2] + t*p2[2];
            isect_point_[4*vi+3] = 1.0f;
            
            if(immediate_state_.buffer[1].is_enabled()) {
                const GLUPfloat* c1 =
                    immediate_state_.buffer[1].element_ptr(v1);
                const GLUPfloat* c2 =
                    immediate_state_.buffer[1].element_ptr(v2);
                isect_color_[4*vi+0] = s*c1[0] + t*c2[0];
                isect_color_[4*vi+1] = s*c1[1] + t*c2[1];
                isect_color_[4*vi+2] = s*c1[2] + t*c2[2];
                isect_color_[4*vi+3] = s*c1[3] + t*c2[3];                
            }
            
            if(immediate_state_.buffer[2].is_enabled()) {
                const GLUPfloat* tex1 =
                    immediate_state_.buffer[2].element_ptr(v1);
                const GLUPfloat* tex2 =
                    immediate_state_.buffer[2].element_ptr(v2);
                
                isect_tex_coord_[4*vi+0] = s*tex1[0] + t*tex2[0];
                isect_tex_coord_[4*vi+1] = s*tex1[1] + t*tex2[1];
                isect_tex_coord_[4*vi+2] = s*tex1[2] + t*tex2[2];
                isect_tex_coord_[4*vi+3] = s*tex1[3] + t*tex2[3];
            }
        }

        /**
         * \brief Assemble the configuration code of a primitive
         *  relative to the clipping plane.
         * \param[in] first_v index of the first vertex of the 
         *  primitive in the immediate buffer
         * \param[in] nb_v number of vertices of the primitive
         * \return an integer with the i-th bit set if vertex i
         *  is visible, and unset if it is clipped.
         */
        index_t get_config(index_t first_v, index_t nb_v) {
            index_t result = 0;
            for(index_t lv=0; lv<nb_v; ++lv) {
                if(v_is_visible_[first_v+lv]) {
                    result = result | (1u << lv);
                }
            }
            return result;
        }

        /**
         * \brief Draws all the primitives from the immediate buffer using
         *  the marching cells algorithm.
         * \details This function is used when clipping is enabled and when
         *  clippping mode is GLUP_CLIP_SLICE_CELLS
         */
        void draw_immediate_buffer_with_marching_cells(
            const MarchingCell& cell
        );
        
    private:
        std::map<std::string, GLsizei> variable_to_offset_;

        /**
         * \brief Indicates for a given vertex whether it is clipped or
         *  is visible, according to the current clipping plane.
         */
        bool v_is_visible_[IMMEDIATE_BUFFER_SIZE];

        /**
         * \brief Indicates whether a picking id should be send to 
         *  OpenGL for each primitive.
         */
        bool pick_primitives_;

        /**
         * \brief computed intersections.
         * \details Used when clipping mode is GLUP_CLIP_SLICE_CELLS.
         */
        GLUPfloat isect_point_[12*4];

        /**
         * \brief computed colors of intersections.
         * \details Used when clipping mode is GLUP_CLIP_SLICE_CELLS.
         */
        GLUPfloat isect_color_[12*4];

        /**
         * \brief computed texture coordinates of intersections.
         * \details Used when clipping mode is GLUP_CLIP_SLICE_CELLS.
         */
        GLUPfloat isect_tex_coord_[12*4];


        /**
         * \brief Cached pointer to uniform state variable.
         */
        GLUPfloat* world_clip_plane_;
    };

    /*********************************************************************/

    
}

#endif
