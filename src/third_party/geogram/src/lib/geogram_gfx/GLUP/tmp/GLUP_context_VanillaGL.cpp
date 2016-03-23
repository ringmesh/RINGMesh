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

#include <geogram_gfx/GLUP/GLUP_context_VanillaGL.h>
#include <geogram_gfx/basic/GLSL.h>
#include <geogram/basic/command_line.h>
#include <geogram/basic/progress.h>
#include <geogram/basic/logger.h>

namespace GLUP {
    using namespace GEO;

    /***********************************************************************/

    Context_VanillaGL::Context_VanillaGL() {
    }

    const char* Context_VanillaGL::profile_name() const {
        return "VanillaGL";
    }
    
    void Context_VanillaGL::begin(GLUPprimitive primitive) {

        if(!primitive_info_[primitive].implemented) {
            Logger::warn("GLUP")
                << "glupBegin(): "
                << glup_primitive_name(primitive)
                << " not implemented in this profile" << std::endl;
        }

        update_uniform_buffer();

        if(uniform_state_.toggle[GLUP_VERTEX_COLORS].get()) {
            immediate_state_.buffer[GLUP_COLOR_ATTRIBUTE].enable();
        } else {
            immediate_state_.buffer[GLUP_COLOR_ATTRIBUTE].disable();
        }
        
        if(uniform_state_.toggle[GLUP_TEXTURING].get()) {
            immediate_state_.buffer[GLUP_TEX_COORD_ATTRIBUTE].enable();
        } else {
            immediate_state_.buffer[GLUP_TEX_COORD_ATTRIBUTE].disable();
        }
        
        immediate_state_.begin(primitive);
        prepare_to_draw(primitive);
        configure_OpenGL_texturing();
        configure_OpenGL_lighting();
        configure_OpenGL_clipping();        
        configure_OpenGL_picking();
    }

    void Context_VanillaGL::configure_OpenGL_texturing() {
        glDisable(GL_TEXTURE_1D);
        glDisable(GL_TEXTURE_2D);
        glDisable(GL_TEXTURE_3D);            
        if( !uniform_state_.toggle[GLUP_PICKING].get() &&
            uniform_state_.toggle[GLUP_TEXTURING].get()
        ) {
            switch(uniform_state_.texture_type.get()) {
            case GLUP_TEXTURE_1D:
                glEnable(GL_TEXTURE_1D);
                break;
            case GLUP_TEXTURE_2D:
                glEnable(GL_TEXTURE_2D);
                break;
            case GLUP_TEXTURE_3D:
                glEnable(GL_TEXTURE_3D);
                break;
            }
            glMatrixMode(GL_TEXTURE);
            glLoadMatrixf(get_matrix(GLUP_TEXTURE_MATRIX));
            glMatrixMode(GL_MODELVIEW);

            switch(uniform_state_.texture_mode.get()) {
            case GLUP_TEXTURE_REPLACE:
                // Yes, it's GL_MODULATE also for texture replace,
                // because GL_TEXTURE_REPLACE removes the shading !
                if(uniform_state_.toggle[GLUP_LIGHTING].get()) {
                    glTexEnvi(
                        GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE
                    );
                } else {
                    glTexEnvi(
                        GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE
                    );
                }
                break;
            case GLUP_TEXTURE_MODULATE:
                glTexEnvi(
                    GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE
                );
                break;
            case GLUP_TEXTURE_ADD:
                glTexEnvi(
                    GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_ADD
                );
                break;
            }
        }
    }
   
    void Context_VanillaGL::configure_OpenGL_lighting() {
        GLUPprimitive primitive = immediate_state_.primitive();
        switch(primitive) {
        case GLUP_POINTS:
        case GLUP_LINES:
            glDisable(GL_LIGHTING);
            break;
        default:
            if(uniform_state_.toggle[GLUP_LIGHTING].get()) {
                glEnable(GL_LIGHTING);
                glEnable(GL_NORMALIZE);
            } else {
                glDisable(GL_LIGHTING);                
            }
            break;
        }
        if(
            uniform_state_.toggle[GLUP_VERTEX_COLORS].get() ||
            primitive == GLUP_POINTS ||
            primitive == GLUP_LINES
        ) {
            glEnable(GL_COLOR_MATERIAL);
        } else {
            glEnable(GL_COLOR_MATERIAL);
            glColor3f(1.0f, 1.0f, 1.0f);
            glDisable(GL_COLOR_MATERIAL);            
        }

        
        glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);
        glColor4fv(uniform_state_.color[GLUP_FRONT_COLOR].get_pointer());
        
        glMaterialfv(
            GL_FRONT, GL_DIFFUSE,
            uniform_state_.color[GLUP_FRONT_COLOR].get_pointer()
        );        
        glMaterialfv(
            GL_BACK, GL_DIFFUSE,
            uniform_state_.color[GLUP_BACK_COLOR].get_pointer()
        );
        static GLfloat specular[4] = { 1.0f, 1.0f, 1.0f, 1.0f };
        static GLfloat ambient[4]  = { 0.2f, 0.2f, 0.2f, 1.0f };
        glMaterialfv(
            GL_FRONT_AND_BACK, GL_AMBIENT, ambient
        );
        glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, specular);        
        glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, 30);
        glLightModeli(GL_LIGHT_MODEL_COLOR_CONTROL, GL_SEPARATE_SPECULAR_COLOR);
    }

    void Context_VanillaGL::configure_OpenGL_picking() {
        if(uniform_state_.toggle[GLUP_PICKING].get()) {
            
            glEnable(GL_COLOR_MATERIAL);
            glDisable(GL_LIGHTING);
            
            // Disable colors and texture coordinates
            immediate_state_.buffer[1].disable();
            immediate_state_.buffer[2].disable();

            //   Disable buffers for interpolated clipping
            // attributes.
            glDisableClientState(GL_COLOR_ARRAY);
            glDisableClientState(GL_TEXTURE_COORD_ARRAY);
            
            uniform_state_.base_picking_id.set(0);
            
            switch(uniform_state_.picking_mode.get()) {
            case GLUP_PICK_PRIMITIVE: {
                pick_primitives_ = true;
            } break;
            case GLUP_PICK_CONSTANT: {
                pick_primitives_ = false;                
                glPickingIdAsColor(index_t(uniform_state_.picking_id.get()));
            } break;
            }
        } else {
            pick_primitives_ = false;
        }
    }

    void Context_VanillaGL::configure_OpenGL_clipping() {
        if(clip_slice_cells()) {
            glNormal3f(
                -world_clip_plane_[0],
                -world_clip_plane_[1],
                -world_clip_plane_[2]            
            );
            glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_FALSE);
            
            // Bind the intersection buffers as vertex array, colors and texture
            // coordinates (if enabled for the last two ones).
            glEnableClientState(GL_VERTEX_ARRAY);
        
            glVertexPointer(
                4,
                GL_FLOAT,
                0,  // stride
                isect_point_ 
            );
        
            if(immediate_state_.buffer[1].is_enabled()) {
                glEnableClientState(GL_COLOR_ARRAY);
                glColorPointer(
                    4,
                    GL_FLOAT,
                    0,  // stride
                    isect_color_ 
                );
            }
        
            if(immediate_state_.buffer[2].is_enabled()) {
                glEnableClientState(GL_TEXTURE_COORD_ARRAY);
                glTexCoordPointer(
                    4,
                    GL_FLOAT,
                    0,  // stride
                    isect_tex_coord_ 
                );
            }
        }
    }
    
    void Context_VanillaGL::end() {
        flush_immediate_buffers();
        if(clip_slice_cells()) {
            glDisableClientState(GL_VERTEX_ARRAY);
            glDisableClientState(GL_COLOR_ARRAY);
            glDisableClientState(GL_TEXTURE_COORD_ARRAY);
        }
    }
    
    void Context_VanillaGL::setup() {
        uniform_buffer_dirty_=true;
        matrices_dirty_=true;        
        uniform_binding_point_=1;

        std::map<std::string, GLsizei> type_to_size;
        type_to_size["bool"] = sizeof(int);
        type_to_size["vec4"] = 4*sizeof(GLfloat);
        type_to_size["vec3"] = 3*sizeof(GLfloat);
        type_to_size["mat4"] = 4*4*sizeof(GLfloat);
        type_to_size["mat3"] = 4*3*sizeof(GLfloat); // yes, 4, there is padding
        type_to_size["float"] = sizeof(GLfloat);
        type_to_size["int"] = sizeof(GLint);

        // Parse uniform state declaration in order to "emulate" it...
        uniform_buffer_size_ = 0;
        std::istringstream input(uniform_state_declaration());
        std::string line;
        while(std::getline(input,line)) {
            std::vector<std::string> words;
            GEO::String::split_string(line, ' ', words);
            if(
                (words.size() == 2 && words[1][words[1].length()-1] == ';') ||
                (words.size() == 3 && words[2] == ";")
            ) {
                std::string vartype = words[0];
                std::string varname = words[1];
                if(varname[varname.length()-1] == ';') {
                    varname = varname.substr(0, varname.length()-1);
                }
                if(type_to_size.find(vartype) != type_to_size.end()) {
                    variable_to_offset_[varname] = uniform_buffer_size_;
                    uniform_buffer_size_ += type_to_size[vartype];
                }
            }
        }
        
        uniform_buffer_data_ = new Memory::byte[uniform_buffer_size_];
        Memory::clear(uniform_buffer_data_, size_t(uniform_buffer_size_));

        setup_state_variables();
        setup_immediate_buffers();
        setup_primitives();

        world_clip_plane_ = uniform_state_.world_clip_plane.get_pointer();
    }

    void Context_VanillaGL::shrink_cells_in_immediate_buffers() {
        if(
            uniform_state_.cells_shrink.get() == 0.0f   ||
            immediate_state_.primitive() == GLUP_POINTS ||
            immediate_state_.primitive() == GLUP_LINES  ||
            (uniform_state_.clipping_mode.get() == GLUP_CLIP_SLICE_CELLS &&
             uniform_state_.toggle[GLUP_CLIPPING].get())
        ) {
            return;
        }
        
        GLfloat s = uniform_state_.cells_shrink.get();
        GLfloat g[3];
        index_t nb_v = nb_vertices_per_primitive[immediate_state_.primitive()];
        index_t v=0;
        while(v < immediate_state_.nb_vertices()) {
            g[0] = 0.0f;
            g[1] = 0.0f;
            g[2] = 0.0f;
            for(index_t lv=0; lv<nb_v; ++lv) {
                GLfloat* p = immediate_state_.buffer[0].element_ptr(v+lv);
                g[0] += p[0];
                g[1] += p[1];
                g[2] += p[2];
            }
            g[0] /= GLfloat(nb_v);
            g[1] /= GLfloat(nb_v);
            g[2] /= GLfloat(nb_v);
            for(index_t lv=0; lv<nb_v; ++lv) {
                GLfloat* p = immediate_state_.buffer[0].element_ptr(v+lv);
                p[0] = s*g[0] + (1.0f - s)*p[0];
                p[1] = s*g[1] + (1.0f - s)*p[1];
                p[2] = s*g[2] + (1.0f - s)*p[2];                    
            }
            v += nb_v;
        }
    }

    void Context_VanillaGL::classify_vertices_in_immediate_buffers() {
        if(!uniform_state_.toggle[GLUP_CLIPPING].get()) {
            return;
        }
        if(uniform_state_.clipping_mode.get() == GLUP_CLIP_STANDARD) {
            return;
        }
        if(
            immediate_state_.primitive() == GLUP_POINTS ||
            immediate_state_.primitive() == GLUP_LINES  ||
            immediate_state_.primitive() == GLUP_TRIANGLES ||
            immediate_state_.primitive() == GLUP_QUADS
        ) {
            return;
        }

        for(index_t v=0; v<immediate_state_.nb_vertices(); ++v) {
            float* p = immediate_state_.buffer[0].element_ptr(v);
            float s = 0.0;
            for(index_t i=0; i<4; ++i) {
                s += world_clip_plane_[i]*p[i];
            }
            v_is_visible_[v] = (s >= 0);
        }
    }

    bool Context_VanillaGL::cell_is_clipped(index_t first_v) {
        if(!uniform_state_.toggle[GLUP_CLIPPING].get()) {
            return false;
        }
        if(uniform_state_.clipping_mode.get() == GLUP_CLIP_STANDARD) {
            return false;
        }
        index_t nb_visible=0;
        index_t nb_in_cell =
            nb_vertices_per_primitive[immediate_state_.primitive()];
        for(index_t lv=0; lv<nb_in_cell; ++lv) {
            nb_visible += (v_is_visible_[first_v + lv]);
        }
        switch(uniform_state_.clipping_mode.get()) {
        case GLUP_CLIP_STRADDLING_CELLS:
            return (nb_visible == 0 || nb_visible == nb_in_cell);
            break;
        case GLUP_CLIP_WHOLE_CELLS:
            return (nb_visible == 0);
            break;
        case GLUP_CLIP_SLICE_CELLS:
            return false;
            break;
        }
        return false;
    }

    void Context_VanillaGL::flush_immediate_buffers() {
        shrink_cells_in_immediate_buffers();
        classify_vertices_in_immediate_buffers();        

        glPushAttrib(
            GL_ENABLE_BIT | GL_LIGHTING_BIT | GL_POLYGON_BIT | GL_TEXTURE_BIT
        );        
        if(
            uniform_state_.clipping_mode.get() != GLUP_CLIP_STANDARD &&
            immediate_state_.primitive() != GLUP_POINTS &&
            immediate_state_.primitive() != GLUP_LINES &&
            immediate_state_.primitive() != GLUP_TRIANGLES &&
            immediate_state_.primitive() != GLUP_QUADS                       
        ) {
            glDisable(GL_CLIP_PLANE0);
        }

        flush_immediate_buffers_once();

        if(uniform_state_.toggle[GLUP_PICKING].get()) {
            if(pick_primitives_) {
                uniform_state_.base_picking_id.set(
                    uniform_state_.base_picking_id.get() +
                    int(immediate_state_.nb_primitives())
                );
            }
        } else if(
            uniform_state_.toggle[GLUP_DRAW_MESH].get() &&
            immediate_state_.primitive() != GLUP_POINTS &&
            immediate_state_.primitive() != GLUP_LINES
        ) {
            // Do it one more time for the mesh
            
            glDisable(GL_LIGHTING);
            glDisable(GL_TEXTURE_1D);
            glDisable(GL_TEXTURE_2D);
            glDisable(GL_TEXTURE_3D);
            glDisable(GL_COLOR_MATERIAL);
            
            // Disable vertex attributes.
            bool va1_enabled = immediate_state_.buffer[1].is_enabled();
            bool va2_enabled = immediate_state_.buffer[2].is_enabled();
            immediate_state_.buffer[1].disable();
            immediate_state_.buffer[2].disable();

            if(clip_slice_cells()) {
                glDisableClientState(GL_COLOR_ARRAY);
                glDisableClientState(GL_TEXTURE_COORD_ARRAY);
            }
            
            glColor3fv(uniform_state_.color[GLUP_MESH_COLOR].get_pointer());
            glLineWidth(uniform_state_.mesh_width.get());
            glPolygonMode(GL_FRONT_AND_BACK,GL_LINE);

            flush_immediate_buffers_once();

            // Restore previous state for vertex attributes.
            if(va1_enabled) {
                immediate_state_.buffer[1].enable();
                if(clip_slice_cells()) {
                    glEnableClientState(GL_COLOR_ARRAY);                    
                }
            }
            if(va2_enabled) {
                immediate_state_.buffer[2].enable();
                if(clip_slice_cells()) {
                    glEnableClientState(GL_TEXTURE_COORD_ARRAY);
                }
            }
        }
        glPopAttrib();        
        immediate_state_.reset();        
    }
    
    void Context_VanillaGL::flush_immediate_buffers_once() {
        switch(immediate_state_.primitive()) {
        case GLUP_POINTS:
            draw_immediate_buffer_GLUP_POINTS();
            break;
        case GLUP_LINES:
            draw_immediate_buffer_GLUP_LINES();
            break;
        case GLUP_TRIANGLES: 
            draw_immediate_buffer_GLUP_TRIANGLES();
            break;
        case GLUP_QUADS:
            draw_immediate_buffer_GLUP_QUADS();
            break;
        case GLUP_TETRAHEDRA:
            draw_immediate_buffer_GLUP_TETRAHEDRA();
            break;
        case GLUP_HEXAHEDRA:
            draw_immediate_buffer_GLUP_HEXAHEDRA();
            break;
        case GLUP_PRISMS:
            draw_immediate_buffer_GLUP_PRISMS();
            break;
        case GLUP_PYRAMIDS:
            draw_immediate_buffer_GLUP_PYRAMIDS();
            break;
        case GLUP_CONNECTORS:
            draw_immediate_buffer_GLUP_CONNECTORS();
            break;
        default:
            Logger::warn("GLUP VanillaGL")
                << glup_primitive_name(immediate_state_.primitive())
                << ":not implemented yet"
                << std::endl;
            break;
        }
    }

    void Context_VanillaGL::setup_immediate_buffers() {
    }

    void Context_VanillaGL::setup_primitives() {
        primitive_info_.resize(GLUP_NB_PRIMITIVES);
        primitive_info_[GLUP_POINTS].implemented = true;
        primitive_info_[GLUP_LINES].implemented = true;
        primitive_info_[GLUP_TRIANGLES].implemented = true;
        primitive_info_[GLUP_QUADS].implemented = true;
        primitive_info_[GLUP_TETRAHEDRA].implemented = true;
        primitive_info_[GLUP_HEXAHEDRA].implemented = true;
        primitive_info_[GLUP_PYRAMIDS].implemented = true;
        primitive_info_[GLUP_PRISMS].implemented = true;
        primitive_info_[GLUP_CONNECTORS].implemented = true;
    }
    
    Memory::pointer Context_VanillaGL::get_state_variable_address(
        const char* name
    ) {
        geo_assert(variable_to_offset_.find(name) != variable_to_offset_.end());
        return uniform_buffer_data_ + variable_to_offset_[name];
    }
    
    bool Context_VanillaGL::primitive_supports_array_mode(
        GLUPprimitive prim
    ) const {
        geo_argused(prim);
        return false;
    }

    void Context_VanillaGL::draw_immediate_buffer_GLUP_POINTS() {
        glDisable(GL_LIGHTING);
        glEnable(GL_POINT_SMOOTH);
        glBegin(GL_POINTS);
        for(index_t v=0; v<immediate_state_.nb_vertices(); ++v) {
            output_picking_id(v);
            output_vertex(v);
        }
        glEnd();
    }

    void Context_VanillaGL::draw_immediate_buffer_GLUP_LINES() {
        glBegin(GL_LINES);
        index_t v = 0;
        while(v < immediate_state_.nb_vertices()) {
            output_picking_id(v/2);
            output_vertex(v);
            output_vertex(v+1);
            v += 2;
        }
        glEnd();
    }

    void Context_VanillaGL::draw_immediate_buffer_GLUP_TRIANGLES() {
        glBegin(GL_TRIANGLES);
        index_t v = 0;
        while(v < immediate_state_.nb_vertices()) {
            output_picking_id(v/3);
            flat_shaded_triangle(v,v+1,v+2);
            v += 3;
        }
        glEnd();
    }

    void Context_VanillaGL::draw_immediate_buffer_GLUP_QUADS() {
        glBegin(GL_QUADS);
        index_t v = 0;
        while(v < immediate_state_.nb_vertices()) {
            output_picking_id(v/4);
            flat_shaded_quad(v,v+1,v+3,v+2);
            v += 4;
        }
        glEnd();
    }


    void Context_VanillaGL::draw_immediate_buffer_with_marching_cells(
        const MarchingCell& cell
    ) {
        
        index_t v0=0;
        while(v0 < immediate_state_.nb_vertices()) {
            index_t config = get_config(v0, cell.nb_vertices());

            //   Compute all the intersection vertices (plus their
            // attributes if enabled).
            for(index_t i=0; i<cell.config_size(config); ++i) {
                index_t e = cell.config_edges(config)[i];
                compute_intersection(
                    v0+cell.edge_vertex(e,0), v0+cell.edge_vertex(e,1), e
                );
            }            


            if(cell.config_size(config) != 0) {
                output_picking_id(v0/cell.nb_vertices());
            
                //   With the bound intersection buffer,
                // we can draw the intersection polygon with
                // a single OpenGL call ! The marching cells table directly
                // refers to the vertices in the intersection buffer.
                glDrawElements(
                    GL_POLYGON,
                    GLsizei(cell.config_size(config)),
                    GL_UNSIGNED_INT,
                    cell.config_edges(config)
                );
            }
            
            v0 += cell.nb_vertices();
        }

    }
    
    void Context_VanillaGL::draw_immediate_buffer_GLUP_TETRAHEDRA() {
        if(clip_slice_cells()) {
            draw_immediate_buffer_with_marching_cells(marching_tet_);
        } else {
            glBegin(GL_TRIANGLES);
            index_t v0 = 0;
            while(v0 < immediate_state_.nb_vertices()) {            
                if(!cell_is_clipped(v0)) {
                    output_picking_id(v0/4);                
                    index_t v1 = v0+1;
                    index_t v2 = v0+2;
                    index_t v3 = v0+3;
                    flat_shaded_triangle(v0,v1,v2);                           
                    flat_shaded_triangle(v1,v0,v3);                             
                    flat_shaded_triangle(v0,v2,v3);                             
                    flat_shaded_triangle(v2,v1,v3);
                }
                v0 += 4;
            }
            glEnd();
        }
    }

    void Context_VanillaGL::draw_immediate_buffer_GLUP_HEXAHEDRA() {
        if(clip_slice_cells()) {
            draw_immediate_buffer_with_marching_cells(marching_hex_);
        } else {
            glBegin(GL_QUADS);
            index_t v0 = 0;
            while(v0 < immediate_state_.nb_vertices()) {
                if(!cell_is_clipped(v0)) {
                    output_picking_id(v0/8);                
                    index_t v1 = v0+1;
                    index_t v2 = v0+2;
                    index_t v3 = v0+3;
                    index_t v4 = v0+4;
                    index_t v5 = v0+5;
                    index_t v6 = v0+6;
                    index_t v7 = v0+7;
                    flat_shaded_quad(v0,v2,v4,v6);
                    flat_shaded_quad(v3,v1,v7,v5);
                    flat_shaded_quad(v1,v0,v5,v4);
                    flat_shaded_quad(v2,v3,v6,v7);
                    flat_shaded_quad(v1,v3,v0,v2);
                    flat_shaded_quad(v4,v6,v5,v7);
                }
                v0 += 8;
            }
            glEnd();
        }
    }

    void Context_VanillaGL::draw_immediate_buffer_GLUP_PRISMS() {
        if(clip_slice_cells()) {
            draw_immediate_buffer_with_marching_cells(marching_prism_);
        } else {
            glBegin(GL_QUADS);
            index_t v0 = 0;
            while(v0 < immediate_state_.nb_vertices()) {
                if(!cell_is_clipped(v0)) {
                    output_picking_id(v0/6);                
                    index_t v1 = v0+1;
                    index_t v2 = v0+2;
                    index_t v3 = v0+3;
                    index_t v4 = v0+4;
                    index_t v5 = v0+5;
                    flat_shaded_quad(v0,v3,v1,v4);
                    flat_shaded_quad(v0,v2,v3,v5);
                    flat_shaded_quad(v1,v4,v2,v5);
                }
                v0 += 6;
            }
            glEnd();
            glBegin(GL_TRIANGLES);
            v0 = 0;
            while(v0 < immediate_state_.nb_vertices()) {
                if(!cell_is_clipped(v0)) {
                    output_picking_id(v0/6);                                
                    index_t v1 = v0+1;
                    index_t v2 = v0+2;
                    index_t v3 = v0+3;
                    index_t v4 = v0+4;
                    index_t v5 = v0+5;
                    flat_shaded_triangle(v0,v1,v2);
                    flat_shaded_triangle(v5,v4,v3);
                }
                v0 += 6;
            }
            glEnd();
        }
    }

    void Context_VanillaGL::draw_immediate_buffer_GLUP_PYRAMIDS() {
        if(clip_slice_cells()) {
            draw_immediate_buffer_with_marching_cells(marching_pyramid_);
        } else {
            glBegin(GL_QUADS);
            index_t v0 = 0;
            while(v0 < immediate_state_.nb_vertices()) {
                if(!cell_is_clipped(v0)) {
                    output_picking_id(v0/5);                                
                    index_t v1 = v0+1;
                    index_t v2 = v0+2;
                    index_t v3 = v0+3;
                    flat_shaded_quad(v0,v1,v3,v2);
                }
                v0 += 5;
            }
            glEnd();
            glBegin(GL_TRIANGLES);
            v0 = 0;
            while(v0 < immediate_state_.nb_vertices()) {
                if(!cell_is_clipped(v0)) {
                    output_picking_id(v0/5);                
                    index_t v1 = v0+1;
                    index_t v2 = v0+2;
                    index_t v3 = v0+3;
                    index_t v4 = v0+4;
                    flat_shaded_triangle(v0,v4,v1);
                    flat_shaded_triangle(v0,v3,v4);
                    flat_shaded_triangle(v2,v4,v3);
                    flat_shaded_triangle(v2,v1,v4);
                }
                v0 += 5;
            }
            glEnd();
        }
    }

    void Context_VanillaGL::draw_immediate_buffer_GLUP_CONNECTORS() {
        if(clip_slice_cells()) {
            draw_immediate_buffer_with_marching_cells(marching_connector_);
        } else {
            glBegin(GL_QUADS);
            index_t v0 = 0;
            while(v0 < immediate_state_.nb_vertices()) {
                if(!cell_is_clipped(v0)) {
                    output_picking_id(v0/4);                                
                    index_t v1 = v0+1;
                    index_t v2 = v0+2;
                    index_t v3 = v0+3;
                    flat_shaded_quad(v0,v1,v3,v2);
                }
                v0 += 4;
            }
            glEnd();
            glBegin(GL_TRIANGLES);
            v0 = 0;
            while(v0 < immediate_state_.nb_vertices()) {
                if(!cell_is_clipped(v0)) {
                    output_picking_id(v0/4);                
                    index_t v1 = v0+1;
                    index_t v2 = v0+2;
                    index_t v3 = v0+3;
                    flat_shaded_triangle(v2,v1,v0);
                    flat_shaded_triangle(v3,v2,v0);
                }
                v0 += 4;
            }
            glEnd();
        }
    }
    
    void Context_VanillaGL::output_normal(index_t v1, index_t v2, index_t v3) {
        GLfloat* p1 = immediate_state_.buffer[0].element_ptr(v1);
        GLfloat* p2 = immediate_state_.buffer[0].element_ptr(v2);
        GLfloat* p3 = immediate_state_.buffer[0].element_ptr(v3);

        // scale vector components, else it can generate floating
        // point exceptions when manipulating very small vectors
        // (e.g. at the poles of the sphere generated in Graphite).
        const float s = 100.0f;
        
        GLfloat U[3];
        U[0] = s*(p2[0] - p1[0]);
        U[1] = s*(p2[1] - p1[1]);
        U[2] = s*(p2[2] - p1[2]);
        
        GLfloat V[3];            
        V[0] = s*(p3[0] - p1[0]);
        V[1] = s*(p3[1] - p1[1]);
        V[2] = s*(p3[2] - p1[2]);
        
        glNormal3f(
            U[1]*V[2] - U[2]*V[1],
            U[2]*V[0] - U[0]*V[2],
            U[0]*V[1] - U[1]*V[0]                
        );
    }


    void Context_VanillaGL::output_normal(
        index_t v1, index_t v2, index_t v3, index_t v4
    ) {
        GLfloat* p1 = immediate_state_.buffer[0].element_ptr(v1);
        GLfloat* p2 = immediate_state_.buffer[0].element_ptr(v2);
        GLfloat* p3 = immediate_state_.buffer[0].element_ptr(v3);
        GLfloat* p4 = immediate_state_.buffer[0].element_ptr(v4);        

        // scale vector components, else it can generate floating
        // point exceptions when manipulating very small vectors
        // (e.g. at the poles of the sphere generated in Graphite).
        const float s = 100.0f;
        
        GLfloat U1[3];
        U1[0] = s*(p2[0] - p1[0]);
        U1[1] = s*(p2[1] - p1[1]);
        U1[2] = s*(p2[2] - p1[2]);
        
        GLfloat V1[3];            
        V1[0] = s*(p4[0] - p1[0]);
        V1[1] = s*(p4[1] - p1[1]);
        V1[2] = s*(p4[2] - p1[2]);

        GLfloat U2[3];
        U2[0] = s*(p4[0] - p3[0]);
        U2[1] = s*(p4[1] - p3[1]);
        U2[2] = s*(p4[2] - p3[2]);
        
        GLfloat V2[3];            
        V2[0] = s*(p2[0] - p3[0]);
        V2[1] = s*(p2[1] - p3[1]);
        V2[2] = s*(p2[2] - p3[2]);
        
        glNormal3f(
            (U1[1]*V1[2]-U1[2]*V1[1]) - (U2[1]*V2[2]-U2[2]*V2[1]),
            (U1[2]*V1[0]-U1[0]*V1[2]) - (U2[2]*V2[0]-U2[0]*V2[2]),
            (U1[0]*V1[1]-U1[1]*V1[0]) - (U2[0]*V2[1]-U2[1]*V2[0])               
        );
    }
}

