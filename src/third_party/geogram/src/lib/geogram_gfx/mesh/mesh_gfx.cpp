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

#include <geogram_gfx/mesh/mesh_gfx.h>
#include <geogram/basic/logger.h>
#include <geogram/basic/command_line.h>
#include <geogram_gfx/basic/GLSL.h>

// TODO: implement attribute display for cell facets.
// TODO: use vertex arrays for attribute display for
//   vertex attributes whenever possible.

namespace {
    using namespace GEO;
}

namespace GEO {

    MeshGfx::MeshGfx() {
        show_mesh_ = true;
        mesh_width_ = 1;
        mesh_border_width_ = 2;
        shrink_ = 0.0;
        animate_ = false;
        time_ = 0.0f;
        draw_cells_[MESH_TET] = true;
        draw_cells_[MESH_HEX] = true;
        draw_cells_[MESH_PRISM] = true;
        draw_cells_[MESH_PYRAMID] = true;
        draw_cells_[MESH_CONNECTOR] = true;
        points_size_ = 1.0f;
        set_points_color(0.0f, 1.0f, 0.0f);
        set_mesh_color(0.0f, 0.0f, 0.0f);
        set_surface_color(0.0f, 0.5f, 1.0f);
        set_backface_surface_color(1.0f, 0.0f, 0.0f);
        set_cells_color(0.9f, 0.9f, 0.9f);
        cells_colors_by_type_ = false;
        lighting_ = true;
        picking_mode_ = MESH_NONE;
        object_picking_id_ = index_t(-1);
        mesh_ = nil;
        triangles_and_quads_ = true;

        buffer_objects_dirty_ = false;
        attributes_buffer_objects_dirty_ = false;
        
        vertices_VAO_ = 0;
        edges_VAO_ = 0;
        facets_VAO_ = 0;
        cells_VAO_ = 0;
        
        vertices_VBO_  = 0;
        edge_indices_VBO_ = 0;
        facet_indices_VBO_ = 0;
        cell_indices_VBO_ = 0;
        vertices_attribute_VBO_ = 0;
        
        do_animation_ = false;

        attribute_subelements_ = MESH_NONE;
        attribute_min_ = 0.0;
        attribute_max_ = 0.0;
        attribute_colormap_texture_ = 0;
        attribute_repeat_ = 1;
    }

    MeshGfx::~MeshGfx() {
        if(vertices_VAO_ != 0) {
            glDeleteVertexArrays(1,&vertices_VAO_);
            vertices_VAO_ = 0;
        }

        if(edges_VAO_ != 0) {
            glDeleteVertexArrays(1,&edges_VAO_);
            edges_VAO_ = 0;
        }

        if(facets_VAO_ != 0) {
            glDeleteVertexArrays(1,&facets_VAO_);
            facets_VAO_ = 0;
        }

        if(cells_VAO_ != 0) {
            glDeleteVertexArrays(1,&cells_VAO_);
            cells_VAO_ = 0;
        }
        
        if(vertices_VBO_ != 0) {
            glDeleteBuffers(1,&vertices_VBO_);
            vertices_VBO_ = 0;
        }

        if(edge_indices_VBO_ != 0) {
            glDeleteBuffers(1,&edge_indices_VBO_);
            edge_indices_VBO_ = 0;
        }
        
        if(facet_indices_VBO_ != 0) {
            glDeleteBuffers(1,&facet_indices_VBO_);
            facet_indices_VBO_ = 0;
        }

        if(cell_indices_VBO_ != 0) {
            glDeleteBuffers(1,&cell_indices_VBO_);
            cell_indices_VBO_ = 0;
        }

        if(vertices_attribute_VBO_ != 0) {
            glDeleteBuffers(1,&vertices_attribute_VBO_);
            vertices_attribute_VBO_ = 0;
        }
    }

    bool MeshGfx::can_use_array_mode(GLUPprimitive prim) const {
        if(do_animation_) {
            return false;
        }
        if(!glupPrimitiveSupportsArrayMode(prim)) {
            return false;
        }
        if(attribute_subelements_ != MESH_NONE) {
            return false;
        }
        return true;
    }
    
    void MeshGfx::draw_vertices() {
        set_GLUP_parameters();
        set_GLUP_picking(MESH_VERTICES);
        update_buffer_objects_if_needed();
        
        glupEnable(GLUP_LIGHTING);
        glupSetColor3fv(GLUP_FRONT_COLOR, points_color_);
        glPointSize(points_size_ * 5.0f);

        if(vertices_selection_ == "") {

            if(
                vertices_VAO_ != 0 &&
                glupPrimitiveSupportsArrayMode(GLUP_POINTS)
                // can_use_array_mode(GLUP_POINTS)
            ) {
                glBindVertexArray(vertices_VAO_);
                if(attribute_subelements_ == MESH_VERTICES) {
                    begin_attributes();
                }
                glupDrawArrays(GLUP_POINTS, 0, GLUPsizei(mesh_->vertices.nb()));
                glBindVertexArray(0);
                if(attribute_subelements_ == MESH_VERTICES) {
                    end_attributes();
                }
            } else {
                if(attribute_subelements_ == MESH_VERTICES) {
                    begin_attributes();
                    glupBegin(GLUP_POINTS);
                    for(index_t v=0; v<mesh_->vertices.nb(); ++v) {
                        draw_vertex_with_attribute(v);
                    }
                    glupEnd();
                    end_attributes();
                } else {
                    glupBegin(GLUP_POINTS);
                    for(index_t v=0; v<mesh_->vertices.nb(); ++v) {
                        draw_vertex(v);
                    }
                    glupEnd();
                }
            }
            
        }  else if(picking_mode_ == MESH_NONE) {
            Attribute<bool> v_selection;
            v_selection.bind_if_is_defined(
                mesh_->vertices.attributes(), vertices_selection_
            );
            if(v_selection.is_bound()) {
                glupBegin(GLUP_POINTS);            
                for(index_t v=0; v<mesh_->vertices.nb(); ++v) {
                    if(v_selection[v]) {
                        draw_vertex(v);
                    }
                }
                glupEnd();
            }
        }
        
        glupDisable(GLUP_PICKING);        
    }

    void MeshGfx::draw_edges() {
        set_GLUP_parameters();
        set_GLUP_picking(MESH_EDGES);
        update_buffer_objects_if_needed();
        
        glupSetColor3fv(GLUP_FRONT_COLOR, mesh_color_);
        glLineWidth(GLfloat(mesh_width_));
        if(can_use_array_mode(GLUP_LINES) && edges_VAO_ != 0) {
            glBindVertexArray(edges_VAO_);
            glupDrawElements(
                GLUP_LINES,
                GLUPsizei(mesh_->edges.nb()*2),
                GL_UNSIGNED_INT,
                0
            );
            glBindVertexArray(0);
        } else {
            if(attribute_subelements_ == MESH_VERTICES) {
                begin_attributes();
                glupBegin(GLUP_LINES);
                for(index_t e=0; e<mesh_->edges.nb(); ++e) {
                    index_t v1 = mesh_->edges.vertex(e,0);
                    index_t v2 = mesh_->edges.vertex(e,1);
                    draw_vertex_with_attribute(v1);
                    draw_vertex_with_attribute(v2);
                }
                glupEnd();
                end_attributes();            
            } else if(attribute_subelements_ == MESH_EDGES) {
                begin_attributes();
                glupBegin(GLUP_LINES);
                for(index_t e=0; e<mesh_->edges.nb(); ++e) {
                    index_t v1 = mesh_->edges.vertex(e,0);
                    index_t v2 = mesh_->edges.vertex(e,1);
                    if(picking_mode_ == MESH_NONE) {
                        glupTexCoord1d(attribute_[e]);
                    }
                    draw_vertex(v1);
                    draw_vertex(v2);
                }
                glupEnd();
                end_attributes();            
            } else {
                glupBegin(GLUP_LINES);
                for(index_t e=0; e<mesh_->edges.nb(); ++e) {
                    index_t v1 = mesh_->edges.vertex(e,0);
                    index_t v2 = mesh_->edges.vertex(e,1);
                    draw_vertex(v1);
                    draw_vertex(v2);
                }
                glupEnd();
            }
        }
    }

    void MeshGfx::draw_surface() {
        set_GLUP_parameters();
        set_GLUP_picking(MESH_FACETS);
        update_buffer_objects_if_needed();

        glupSetCellsShrink(0.0f);

        if(attribute_subelements_ != MESH_NONE) {
            glupSetColor3f(GLUP_FRONT_AND_BACK_COLOR, 1.0f, 1.0f, 1.0f);
        } else {
            glupSetColor3fv(GLUP_FRONT_COLOR, surface_color_);
            glupSetColor3fv(GLUP_BACK_COLOR, backface_surface_color_);
        }

        if(
            (picking_mode_ == MESH_NONE) && (
            attribute_subelements_ == MESH_VERTICES ||
            attribute_subelements_ == MESH_FACETS ||
            attribute_subelements_ == MESH_FACET_CORNERS) 
        ) {
            if(mesh_->facets.are_simplices()) {
                begin_attributes();
                glupBegin(GLUP_TRIANGLES);
                for(index_t f=0; f<mesh_->facets.nb(); ++f) {
                    for(
                        index_t c=mesh_->facets.corners_begin(f);
                        c<mesh_->facets.corners_end(f); ++c) {
                        index_t v=mesh_->facet_corners.vertex(c);
                        draw_surface_vertex_with_attribute(v,f,c);
                    }
                }
                glupEnd();
                end_attributes();
            } else if(triangles_and_quads_) {
                begin_attributes();
                glupBegin(GLUP_TRIANGLES);
                for(index_t f=0; f<mesh_->facets.nb(); ++f) {
                    if(mesh_->facets.nb_vertices(f) != 3) {
                        continue;
                    }
                    for(
                        index_t c=mesh_->facets.corners_begin(f);
                        c<mesh_->facets.corners_end(f); ++c) {
                        index_t v=mesh_->facet_corners.vertex(c);
                        draw_surface_vertex_with_attribute(v,f,c);
                    }
                }
                glupEnd();
                glupBegin(GLUP_QUADS);
                for(index_t f=0; f<mesh_->facets.nb(); ++f) {
                    if(mesh_->facets.nb_vertices(f) != 4) {
                        continue;
                    }
                    for(
                        index_t c=mesh_->facets.corners_begin(f);
                        c<mesh_->facets.corners_end(f); ++c) {
                        index_t v=mesh_->facet_corners.vertex(c);
                        draw_surface_vertex_with_attribute(v,f,c);
                    }
                }
                glupEnd();
                end_attributes();
            } else {
                begin_attributes();
                glupDisable(GLUP_DRAW_MESH);
                glupBegin(GLUP_TRIANGLES);
                for(index_t f=0; f<mesh_->facets.nb(); ++f) {
                    index_t c1 = mesh_->facets.corners_begin(f);
                    index_t v1 = mesh_->facet_corners.vertex(c1);
                    for(
                        index_t c2 = c1+1;
                        c2+1<mesh_->facets.corners_end(f); ++c2
                    ) {
                        index_t c3=c2+1;
                        index_t v2=mesh_->facet_corners.vertex(c2);
                        index_t v3=mesh_->facet_corners.vertex(c3);
                        draw_surface_vertex_with_attribute(v1,f,c1);
                        draw_surface_vertex_with_attribute(v2,f,c2);
                        draw_surface_vertex_with_attribute(v3,f,c3);
                    }
                }
                glupEnd();
                end_attributes();
                if(show_mesh_ && (picking_mode_ == MESH_NONE)) {
                    glupDisable(GLUP_VERTEX_COLORS);                            
                    draw_surface_mesh_with_lines();
                }
            }
        } else if(mesh_->facets.are_simplices() && !do_animation_) {
            if(
                facets_VAO_ != 0 &&
                can_use_array_mode(GLUP_TRIANGLES)
            ) {
                glBindVertexArray(facets_VAO_);
                glupDrawElements(
                    GLUP_TRIANGLES,
                    GLUPsizei(mesh_->facet_corners.nb()),
                    GL_UNSIGNED_INT,
                    0
                );
                glBindVertexArray(0);
            } else {
                glupBegin(GLUP_TRIANGLES);
                for(index_t t=0; t<mesh_->facets.nb(); ++t) {
                    draw_vertex(mesh_->facets.vertex(t,0));
                    draw_vertex(mesh_->facets.vertex(t,1));
                    draw_vertex(mesh_->facets.vertex(t,2));
                }
                glupEnd();
            }
        } else if(triangles_and_quads_ && (picking_mode_ == MESH_NONE)) {
            glupBegin(GLUP_TRIANGLES);
            for(index_t f=0; f<mesh_->facets.nb(); ++f) {
                if(mesh_->facets.nb_vertices(f) == 3) {
                    draw_vertex(mesh_->facets.vertex(f,0));
                    draw_vertex(mesh_->facets.vertex(f,1));
                    draw_vertex(mesh_->facets.vertex(f,2));
                }
            }
            glupEnd();
            glupBegin(GLUP_QUADS);
            for(index_t f=0; f<mesh_->facets.nb(); ++f) {
                if(mesh_->facets.nb_vertices(f) == 4) {
                    draw_vertex(mesh_->facets.vertex(f,0));
                    draw_vertex(mesh_->facets.vertex(f,1));
                    draw_vertex(mesh_->facets.vertex(f,2));
                    draw_vertex(mesh_->facets.vertex(f,3));
                }
            }
            glupEnd();
            
        } else {
            glupDisable(GLUP_DRAW_MESH);

            // Using vertex colors to do the picking.
            if(picking_mode_ != MESH_NONE) {
                glupDisable(GLUP_PICKING);
                glupDisable(GLUP_LIGHTING);
                glupEnable(GLUP_VERTEX_COLORS);
            }
            
            glupBegin(GLUP_TRIANGLES);
            bool picking_vertex_colors = false;
            if(picking_mode_ != MESH_NONE) {
                picking_vertex_colors = (
                    (picking_mode_ & MESH_FACETS) != 0 &&
                    object_picking_id_ == index_t(-1)
                );
                set_GLUP_vertex_color_from_picking_id(object_picking_id_);
            }
            for(index_t f=0; f<mesh_->facets.nb(); ++f) {
                if(picking_vertex_colors) {
                    set_GLUP_vertex_color_from_picking_id(f);      
                }
                index_t v1 = mesh_->facets.vertex(f,0);
                for(index_t lv=1; lv+1<mesh_->facets.nb_vertices(f); ++lv) {
                    index_t v2 = mesh_->facets.vertex(f,lv);
                    index_t v3 = mesh_->facets.vertex(f,lv+1);
                    draw_vertex(v1);
                    draw_vertex(v2);
                    draw_vertex(v3);
                }
            }
            glupEnd();
            if(show_mesh_ && (picking_mode_ == MESH_NONE)) {
                glupDisable(GLUP_VERTEX_COLORS);                            
                draw_surface_mesh_with_lines();
            }
        }
    }

    void MeshGfx::draw_surface_mesh_with_lines() {
        glLineWidth(GLfloat(mesh_width_));
        glupSetColor3fv(GLUP_FRONT_AND_BACK_COLOR, mesh_color_);
        glupBegin(GLUP_LINES);
        for(index_t f=0; f<mesh_->facets.nb(); ++f) {
            for(index_t c1=mesh_->facets.corners_begin(f);
                c1 < mesh_->facets.corners_end(f); ++c1
                ) {
                index_t c2 =
                    mesh_->facets.next_corner_around_facet(f,c1);
                index_t v1 = mesh_->facet_corners.vertex(c1);
                index_t v2 = mesh_->facet_corners.vertex(c2);
                draw_vertex(v1);
                draw_vertex(v2);
            }
        }
        glupEnd();
    }
    
    void MeshGfx::draw_surface_borders() {
        if(picking_mode_ != MESH_NONE) {
            return;
        }
        set_GLUP_parameters();
        glupSetColor3fv(GLUP_FRONT_COLOR, mesh_color_);
        glLineWidth(GLfloat(mesh_border_width_));
        glupBegin(GLUP_LINES);
        for(index_t f=0; f<mesh_->facets.nb(); ++f) {
            for(
                index_t c1=mesh_->facets.corners_begin(f);
                c1<mesh_->facets.corners_end(f); ++c1
            ) {
                if(mesh_->facet_corners.adjacent_facet(c1) == NO_FACET) {
                    index_t v1 = mesh_->facet_corners.vertex(c1);
                    index_t c2 = mesh_->facets.next_corner_around_facet(f,c1);
                    index_t v2 = mesh_->facet_corners.vertex(c2);
                    draw_vertex(v1);
                    draw_vertex(v2);                    
                }
            }
        }
        glupEnd();
    }

    void MeshGfx::draw_volume() {
        if(mesh_->cells.nb() == 0) {
            return;
        }
        set_GLUP_parameters();
        glupSetCellsShrink(GLUPfloat(shrink_));
        if(mesh_->cells.are_simplices() && !do_animation_ &&
           can_use_array_mode(GLUP_TETRAHEDRA)
        ) {
            if(!draw_cells_[MESH_TET]) {
                return;
            }

            glupSetColor3fv(GLUP_FRONT_AND_BACK_COLOR, cells_color_[MESH_TET]);

            update_buffer_objects_if_needed();
            if(
                cells_VAO_ != 0 &&
                can_use_array_mode(GLUP_TETRAHEDRA)
            ) {
                glBindVertexArray(cells_VAO_);
                glupDrawElements(
                    GLUP_TETRAHEDRA,
                    GLUPsizei(mesh_->cells.nb()*4),
                    GL_UNSIGNED_INT,
                    0
                );
                glBindVertexArray(0);
            } else {
                if( (picking_mode_ == MESH_NONE) && (
                        attribute_subelements_ == MESH_VERTICES ||
                        attribute_subelements_ == MESH_CELLS ||
                        attribute_subelements_ == MESH_CELL_CORNERS
                    )
                ) {
                    begin_attributes();
                    glupBegin(GLUP_TETRAHEDRA);
                    for(index_t t=0; t<mesh_->cells.nb(); ++t) {
                        index_t v0 = mesh_->cells.vertex(t,0);
                        index_t v1 = mesh_->cells.vertex(t,1);
                        index_t v2 = mesh_->cells.vertex(t,2);
                        index_t v3 = mesh_->cells.vertex(t,3);
                        index_t c0 = 4*t;
                        draw_volume_vertex_with_attribute(v0, t, c0);
                        draw_volume_vertex_with_attribute(v1, t, c0+1);
                        draw_volume_vertex_with_attribute(v2, t, c0+2);
                        draw_volume_vertex_with_attribute(v3, t, c0+3);
                    }
                    glupEnd();
                    end_attributes();
                } else {
                    glupBegin(GLUP_TETRAHEDRA);
                    for(index_t t=0; t<mesh_->cells.nb(); ++t) {
                        draw_vertex(mesh_->cells.vertex(t,0));
                        draw_vertex(mesh_->cells.vertex(t,1));
                        draw_vertex(mesh_->cells.vertex(t,2));
                        draw_vertex(mesh_->cells.vertex(t,3));                
                    }
                    glupEnd();
                }
            }
        } else {
            static GLUPprimitive geogram_to_glup[MESH_NB_CELL_TYPES] = {
                GLUP_TETRAHEDRA,
                GLUP_HEXAHEDRA,
                GLUP_PRISMS,
                GLUP_PYRAMIDS,
                GLUP_CONNECTORS
            };
            if( (picking_mode_ == MESH_NONE) && (
                    attribute_subelements_ == MESH_VERTICES ||
                    attribute_subelements_ == MESH_CELLS ||
                    attribute_subelements_ == MESH_CELL_CORNERS
                 )
            ) {
                begin_attributes();
                for(index_t type=MESH_TET; type<MESH_NB_CELL_TYPES; ++type) {
                    if(!draw_cells_[type]) {
                        continue;
                    }
                    glupBegin(geogram_to_glup[type]);
                    for(index_t cell=0; cell<mesh_->cells.nb(); ++cell) {
                        if(index_t(mesh_->cells.type(cell)) != type) {
                            continue;
                        }
                        index_t c0 = mesh_->cells.corners_begin(cell);
                        for(index_t lv=0;
                            lv<mesh_->cells.nb_vertices(cell); ++lv
                        ) {
                            draw_volume_vertex_with_attribute(
                                mesh_->cells.vertex(cell,lv),
                                cell,
                                c0+lv
                            );
                        }
                    }
                    glupEnd();
                }
                end_attributes();                        
            } else {
                bool has_cells[MESH_NB_CELL_TYPES];
                for(index_t type=0; type<MESH_NB_CELL_TYPES; ++type) {
                    has_cells[type] = false;
                }
                has_cells[MESH_TET] = true;
                for(index_t type=MESH_TET; type < MESH_NB_CELL_TYPES; ++type) {
                    if(!draw_cells_[type] || !has_cells[type]) {
                        continue;
                    }
                    glupSetColor3fv(
                        GLUP_FRONT_AND_BACK_COLOR, cells_color_[type]
                    );
                    glupBegin(geogram_to_glup[type]);
                    for(index_t cell=0; cell<mesh_->cells.nb(); ++cell) {
                        index_t this_cell_type =
                            index_t(mesh_->cells.type(cell));
                        if(this_cell_type != type) {
                            has_cells[mesh_->cells.type(cell)] = true;
                            continue;
                        }
                        for(index_t lv=0;
                            lv<mesh_->cells.nb_vertices(cell); ++lv
                        ) {
                            draw_vertex(mesh_->cells.vertex(cell,lv));
                        }
                    }
                    glupEnd();
                }
            }
        }
    }

    void MeshGfx::set_mesh(const Mesh* mesh) {
        mesh_ = mesh;
        triangles_and_quads_ = true;
        if(mesh_ != nil) {
            for(index_t f = 0; f<mesh_->facets.nb(); ++f) {
                index_t nb = mesh_->facets.nb_vertices(f);
                if(nb != 3 && nb != 4) {
                    triangles_and_quads_ = false;
                    break;
                }
            }
        }
        buffer_objects_dirty_ = true;
        attributes_buffer_objects_dirty_ = true;
    }

    void MeshGfx::set_GLUP_parameters() {
        if(glupCurrentContext() == nil) {
            glupMakeCurrent(glupCreateContext());
        }
        glupCopyFromGLState(GLUP_ALL_ATTRIBUTES);
        if(show_mesh_) {
            glupEnable(GLUP_DRAW_MESH);
        } else {
            glupDisable(GLUP_DRAW_MESH);
        }
        glupSetColor3fv(GLUP_MESH_COLOR, mesh_color_);
        glupSetMeshWidth(GLUPint(mesh_width_));
        glPointSize(points_size_);
        if(lighting_) {
            glupEnable(GLUP_LIGHTING);
        } else {
            glupDisable(GLUP_LIGHTING);
        }

        do_animation_ =
            (animate_ && mesh_->vertices.dimension() >= 6);
    }

    void MeshGfx::set_GLUP_picking(MeshElementsFlags what) {
        if(picking_mode_ == MESH_NONE && object_picking_id_ == index_t(-1)) {
            glupDisable(GLUP_PICKING);
        } else {
            glupEnable(GLUP_PICKING);
            if(
                (object_picking_id_ == index_t(-1)) &&
                ((picking_mode_ & what) != 0)
            ) {
                glupPickingMode(GLUP_PICK_PRIMITIVE);                    
            } else {
                glupPickingMode(GLUP_PICK_CONSTANT);
                glupPickingId(object_picking_id_);
            }
        }
    }

    void MeshGfx::set_GLUP_vertex_color_from_picking_id(index_t id) {
        GLubyte r = GLubyte( id        & 255);
        GLubyte g = GLubyte((id >> 8)  & 255);
        GLubyte b = GLubyte((id >> 16) & 255);
        GLubyte a = GLubyte((id >> 24) & 255);
        glupColor4f(
            GLfloat(r)/255.0f,
            GLfloat(g)/255.0f,
            GLfloat(b)/255.0f,
            GLfloat(a)/255.0f
        );
    }

    void MeshGfx::bind_vertices_VBO() {
        glBindBuffer(GL_ARRAY_BUFFER, vertices_VBO_);
        glEnableVertexAttribArray(0);

        GLint dim = GLint(geo_min(3u, mesh_->vertices.dimension()));
        
        if(mesh_->vertices.single_precision()) {
            GLsizei stride = GLsizei(
                mesh_->vertices.dimension() * sizeof(float)
            );
            glVertexAttribPointer(
                0,        // Attribute 0
                dim,      // nb coordinates per vertex
                GL_FLOAT, // input coordinates representation
                GL_FALSE, // do not normalize
                stride,   // offset between two consecutive vertices
                0         // addr. relative to bound VBO 
            );
        } else {
            GLsizei stride = GLsizei(
                mesh_->vertices.dimension() * sizeof(double)
            );
            glVertexAttribPointer(
                0,         // Attribute 0
                dim,         // nb coordinates per vertex
                GL_DOUBLE, // input coordinates representation
                GL_FALSE,  // do not normalize
                stride,    // offset between two consecutive vertices
                0          // addr. relative to bound VBO 
            );
        }
    }

    
    void MeshGfx::update_buffer_objects_if_needed() {
        if(mesh_->vertices.nb() == 0) {
            return;
        }

        if(!buffer_objects_dirty_) {
            update_attribute_buffer_objects_if_needed();
            return;
        }

        if(!strcmp(glupCurrentProfileName(),"VanillaGL")) {
            return;
        }
        
        if(mesh_->vertices.single_precision()) {
            size_t size = mesh_->vertices.nb() *
                mesh_->vertices.dimension() * sizeof(float);
            update_or_check_buffer_object(
                vertices_VBO_, GL_ARRAY_BUFFER,
                size, mesh_->vertices.single_precision_point_ptr(0),
                buffer_objects_dirty_
            );
        } else {
            size_t size = mesh_->vertices.nb() *
                mesh_->vertices.dimension() * sizeof(double);
            
            update_or_check_buffer_object(
                vertices_VBO_, GL_ARRAY_BUFFER,
                size, mesh_->vertices.point_ptr(0),
                buffer_objects_dirty_
            );
        }

        if(vertices_VAO_ == 0) {
            glGenVertexArrays(1, &vertices_VAO_);
        }
        glBindVertexArray(vertices_VAO_);
        bind_vertices_VBO();
        glBindVertexArray(0);

        if(mesh_->edges.nb()) {
            update_or_check_buffer_object(
                edge_indices_VBO_, GL_ELEMENT_ARRAY_BUFFER,
                mesh_->edges.nb() * 2 * sizeof(int),
                mesh_->edges.vertex_index_ptr(0),
                buffer_objects_dirty_
            );
            if(edges_VAO_ == 0) {
                glGenVertexArrays(1, &edges_VAO_);
            }
            glBindVertexArray(edges_VAO_);
            bind_vertices_VBO();
            glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, edge_indices_VBO_);
            glBindVertexArray(0);
        }
        
        if(mesh_->facets.nb() != 0 && mesh_->facets.are_simplices()) {
            update_or_check_buffer_object(
                facet_indices_VBO_, GL_ELEMENT_ARRAY_BUFFER,
                mesh_->facet_corners.nb() * sizeof(int),
                mesh_->facet_corners.vertex_index_ptr(0),
                buffer_objects_dirty_
            );
            if(facets_VAO_ == 0) {
                glGenVertexArrays(1, &facets_VAO_);
            }
            glBindVertexArray(facets_VAO_);
            bind_vertices_VBO();
            glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, facet_indices_VBO_);
            glBindVertexArray(0);
        }
        
        if(mesh_->cells.nb() != 0 && mesh_->cells.are_simplices()) {
            update_or_check_buffer_object(
                cell_indices_VBO_, GL_ELEMENT_ARRAY_BUFFER,
                mesh_->cell_corners.nb() * sizeof(int),
                mesh_->cell_corners.vertex_index_ptr(0),
                buffer_objects_dirty_
            );
            if(cells_VAO_ == 0) {
                glGenVertexArrays(1, &cells_VAO_);
            }
            glBindVertexArray(cells_VAO_);
            bind_vertices_VBO();
            glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, cell_indices_VBO_);
            glBindVertexArray(0);
        }
        
        buffer_objects_dirty_ = false;
        update_attribute_buffer_objects_if_needed();
    }

    void MeshGfx::set_scalar_attribute(
        MeshElementsFlags subelements,
        const std::string& name,
        double attr_min, double attr_max,
        GLuint colormap_texture,
        index_t repeat
    ) {
        if(
            subelements != attribute_subelements_ ||
            attribute_name_ != name
        ) {
            attributes_buffer_objects_dirty_ = true;
        }
        attribute_subelements_ = subelements;
        attribute_name_ = name;
        attribute_min_ = attr_min;
        attribute_max_ = attr_max;
        attribute_repeat_ = repeat;
        attribute_colormap_texture_ = colormap_texture;

        const MeshSubElementsStore& mesh_subelements =
            mesh_->get_subelements_by_type(attribute_subelements_);

        if(!ReadOnlyScalarAttributeAdapter::is_defined(
               mesh_subelements.attributes(), attribute_name_
           )
        ) {
            attribute_subelements_ = MESH_NONE;
        }
        
    }
    
    void MeshGfx::update_attribute_buffer_objects_if_needed() {
        if(mesh_->vertices.nb() == 0) {
            return;
        }

        if(!attributes_buffer_objects_dirty_) {
            return;
        }

        if(!strcmp(glupCurrentProfileName(),"VanillaGL")) {
            return;
        }

        if(attribute_subelements_ == MESH_VERTICES) {
            attribute_.bind_if_is_defined(
                mesh_->vertices.attributes(), attribute_name_
            );
            if(attribute_.attribute_store()->dimension() > 4) {
                attribute_.unbind();
            }
        }
        
        if(attribute_.is_bound()) {
            size_t element_size = attribute_.attribute_store()->element_size();
            GLint dimension = GLint(attribute_.attribute_store()->dimension());
            index_t nb_items = attribute_.size();
            GLsizei stride = GLsizei(element_size) * dimension;
            const void* data = attribute_.attribute_store()->data();
            const GLvoid* offset = (const GLvoid*)(
                attribute_.attribute_store()->element_size() *
                index_t(attribute_.element_index())
            );

            update_or_check_buffer_object(
                vertices_attribute_VBO_, GL_ARRAY_BUFFER,
                element_size*index_t(dimension)*nb_items,
                data,
                attributes_buffer_objects_dirty_
            );

            glBindVertexArray(vertices_VAO_);

            glBindBuffer(GL_ARRAY_BUFFER, vertices_attribute_VBO_);
            
            glEnableVertexAttribArray(2); // 2 = tex coords

            switch(attribute_.element_type()) {
            case ReadOnlyScalarAttributeAdapter::ET_UINT8:
                glVertexAttribPointer(
                    2, dimension, GL_UNSIGNED_BYTE, GL_FALSE, stride, offset
                );
            case ReadOnlyScalarAttributeAdapter::ET_INT8:
                glVertexAttribPointer(
                    2, dimension, GL_BYTE, GL_FALSE, stride, offset
                );
                break;
            case ReadOnlyScalarAttributeAdapter::ET_UINT32:
                glVertexAttribPointer(
                    2, dimension, GL_UNSIGNED_INT, GL_FALSE, stride, offset
                );
                break;
            case ReadOnlyScalarAttributeAdapter::ET_INT32:
                glVertexAttribPointer(
                    2, dimension, GL_INT, GL_FALSE, stride, offset
                );
                break;
            case ReadOnlyScalarAttributeAdapter::ET_FLOAT32:
                glVertexAttribPointer(
                    2, dimension, GL_FLOAT, GL_FALSE, stride, offset
                );
                break;
            case ReadOnlyScalarAttributeAdapter::ET_FLOAT64:
                glVertexAttribPointer(
                    2, dimension, GL_DOUBLE, GL_FALSE, stride, offset
                );
                break;
            case ReadOnlyScalarAttributeAdapter::ET_NONE:
                geo_assert_not_reached;
                break;
            }
            glBindVertexArray(0);
            glBindBuffer(GL_ARRAY_BUFFER, 0);
            
            attribute_.unbind();            
        } else {
            glBindVertexArray(vertices_VAO_);
            glDisableVertexAttribArray(2); // 2 = tex coords
            glBindVertexArray(0);
        }
        
        attributes_buffer_objects_dirty_ = false;
    }
    
    void MeshGfx::begin_attributes() {
        if(picking_mode_ != MESH_NONE) {
            return;
        }
        if(attribute_subelements_ == MESH_NONE) {
            return;
        }
        const MeshSubElementsStore& subelements =
            mesh_->get_subelements_by_type(attribute_subelements_);
        attribute_.bind_if_is_defined(
            subelements.attributes(), attribute_name_
        );
        if(!attribute_.is_bound()) {
            return;
        }
        
        glupEnable(GLUP_TEXTURING);
        glupTextureType(GLUP_TEXTURE_2D); // TODO: 1D
        glupTextureMode(GLUP_TEXTURE_REPLACE);

        glActiveTexture(GL_TEXTURE0);
        glBindTexture(GL_TEXTURE_2D, attribute_colormap_texture_);

        // Setup a texture matrix that rescales attribute range
        // from [attribute_min_,attribute_max_] to [0,1]
        glupMatrixMode(GLUP_TEXTURE_MATRIX);
        GLUPdouble M[16];
        Memory::clear(M, sizeof(GLUPdouble)*16);
        double d = attribute_max_ - attribute_min_;
        if(::fabs(d) < 1e-10) {
            d = 0.0;
        } else {
            d = 1.0 / d;
        }
        d *= double(geo_max(attribute_repeat_, 1u));
        M[0] =  d;
        M[12] = -d*attribute_min_;
        M[15] = 1.0;
        glupLoadMatrixd(M);
        glupMatrixMode(GLUP_MODELVIEW_MATRIX);

        glupSetColor3f(GLUP_FRONT_AND_BACK_COLOR, 1.0f, 1.0f, 1.0f);
    }

    void MeshGfx::end_attributes() {
        if(attribute_.is_bound()) {
            glupDisable(GLUP_TEXTURING);
            attribute_.unbind();
        }
    }
}

