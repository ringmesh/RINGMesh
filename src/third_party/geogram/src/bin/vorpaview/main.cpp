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

#include <geogram_gfx/basic/GLSL.h>
#include <geogram_gfx/glup_viewer/glup_viewer.h>
#include <geogram/mesh/mesh.h>
#include <geogram/mesh/mesh_io.h>
#include <geogram_gfx/mesh/mesh_gfx.h>
#include <geogram/basic/command_line.h>
#include <geogram/basic/command_line_args.h>
#include <geogram/basic/file_system.h>
#include <geogram/basic/logger.h>
#include <geogram/basic/stopwatch.h>

namespace {

    GEO::Mesh M(3);
    GEO::MeshGfx M_gfx;

    bool show_colors   = true;
    bool white_bg      = true;
    bool lighting      = true;
    
    bool show_vertices = false;
    GLfloat vertices_size = 1.0f;

    bool show_surface  = true;
    bool show_mesh = true;
    bool show_surface_borders = false;
    
    bool show_volume   = false;
    GLfloat cells_shrink = 0.0f;    
    bool show_hexes = true;    
    bool show_colored_cells = false;

    GLfloat OTM_time = 0.0f;
    GLfloat OTM_speed = 1.0f;

    /**
     * \brief Increments the time of the Optimal Transport animation.
     */
    void increment_time() {
        OTM_time += 0.05f;
        if(OTM_time > 1.0f) {
            OTM_time = 1.0f;
        }
    }

    /**
     * \brief Decrements the time of the Optimal Transport animation.
     */
    void decrement_time() {
        OTM_time -= 0.05f;
        if(OTM_time < 0.0f) {
            OTM_time = 0.0f;
        }
    }


    /**
     * \brief Increments cells shrinkage (shrinks more).
     */
    void increment_shrink() {
        cells_shrink += 0.05f;
        if(cells_shrink > 1.0f) {
            cells_shrink = 1.0f;
        }
    }

    /**
     * \brief Decrements cells shrinkage (shrinks less).
     */
    void decrement_shrink() {
        cells_shrink -= 0.05f;
        if(cells_shrink < 0.0f) {
            cells_shrink = 0.0f;
        }
    }

    void load_mesh(const std::string& filename);
    
    /**
     * \brief Initializes OpenGL objects.
     * \details Specifed as glup_viewer_set_init_func() callback.
     */
    void init() {
        GEO::Graphics::initialize();

        glup_viewer_disable(GLUP_VIEWER_BACKGROUND);        
        glup_viewer_set_background_color(1.0, 1.0, 1.0);
        glup_viewer_add_toggle('b', &white_bg, "white background");
        glup_viewer_add_toggle('c', &show_colors, "colors");

        glup_viewer_add_key_func('r', decrement_time, "Decrement time");
        glup_viewer_add_key_func('t', increment_time, "Increment time");

        glup_viewer_add_key_func('x', decrement_shrink, "Decrement shrink");
        glup_viewer_add_key_func('w', increment_shrink, "Increment shrink");
        
        glup_viewer_add_toggle('p', &show_vertices, "vertices");
        
        glup_viewer_add_toggle('S', &show_surface, "surface");
        glup_viewer_add_toggle('B', &show_surface_borders, "borders");
        glup_viewer_add_toggle('m', &show_mesh, "mesh");

        glup_viewer_add_toggle('V', &show_volume, "volume");
        glup_viewer_add_toggle('j',  &show_hexes, "hexes");
        glup_viewer_add_toggle('C', &show_colored_cells, "colored cells");

#ifdef GEO_OS_EMSCRIPTEN        
        if(GEO::FileSystem::is_file("morph.tet6")) {
            load_mesh("morph.tet6");
            show_vertices = false;
            glup_viewer_enable(GLUP_VIEWER_IDLE_REDRAW);
        }
#endif        
    }

    /**
     * \brief Draws the current mesh according to current rendering mode.
     * \details Specifed as glup_viewer_set_display_func() callback.
     */
    void display() {

        if(glup_viewer_is_enabled(GLUP_VIEWER_IDLE_REDRAW)) {
            OTM_time = float(
                         sin(double(OTM_speed) * GEO::SystemStopwatch::now())
                       );
            OTM_time = 0.5f * (OTM_time + 1.0f);
        }
        
        M_gfx.set_lighting(lighting);
        M_gfx.set_time(double(OTM_time));
        
        if(M_gfx.mesh() != &M) {
            M_gfx.set_mesh(&M);
        }
        
        if(show_vertices) {
            M_gfx.set_points_size(vertices_size);
            M_gfx.draw_vertices();
        }

        if(white_bg) {
            glup_viewer_set_background_color(1.0, 1.0, 1.0);
            M_gfx.set_mesh_color(0.0, 0.0, 0.0);
        } else {
            glup_viewer_set_background_color(0.0, 0.0, 0.0);
            M_gfx.set_mesh_color(1.0, 1.0, 1.0);
        }
        
        if(show_colors) {
            if(M.cells.nb() == 0) {
                M_gfx.set_surface_color(0.5f, 0.75f, 1.0f);
                M_gfx.set_backface_surface_color(1.0f, 0.0f, 0.0f);
            } else {
                M_gfx.set_surface_color(0.7f, 0.0f, 0.0f);
                M_gfx.set_backface_surface_color(1.0f, 1.0f, 0.0f);
            }
        } else {
            if(white_bg) {
                M_gfx.set_surface_color(0.9f, 0.9f, 0.9f);
            } else {
                M_gfx.set_surface_color(0.1f, 0.1f, 0.1f);
            }
        }

        M_gfx.set_show_mesh(show_mesh);

        if(show_surface) {
            M_gfx.draw_surface();
        }
        
        if(show_surface_borders) {
            M_gfx.draw_surface_borders();
        }

        if(show_mesh) {
            M_gfx.draw_edges();
        }

        if(show_volume) {
            M_gfx.set_shrink(double(cells_shrink));
            M_gfx.set_draw_cells(GEO::MESH_HEX, show_hexes);
            if(show_colored_cells) {
                M_gfx.set_cells_colors_by_type();
            } else {
                M_gfx.set_cells_color(0.9f, 0.9f, 0.9f);
            }
            M_gfx.draw_volume();
        }
    }

    /**
     * \brief Gets the bounding box of a mesh animation.
     * \details In animated mode, the mesh animation is stored as a mesh with
     *  6d coordinates, that correspond to the geometric location
     *  at the vertices at time 0 and at time 1.
     * \param[in] M_in the mesh
     * \param[out] xyzmin a pointer to the three minimum coordinates
     * \param[out] xyzmax a pointer to the three maximum coordinates
     * \param[in] animate true if displaying a mesh animation
     */
    void get_bbox(
        const GEO::Mesh& M_in, double* xyzmin, double* xyzmax,
        bool animate
    ) {
        geo_assert(M_in.vertices.dimension() >= GEO::index_t(animate ? 6 : 3));
        for(GEO::index_t c = 0; c < 3; c++) {
            xyzmin[c] = GEO::Numeric::max_float64();
            xyzmax[c] = GEO::Numeric::min_float64();
        }

        for(GEO::index_t v = 0; v < M_in.vertices.nb(); ++v) {
            if(M_in.vertices.single_precision()) {
                const float* p = M_in.vertices.single_precision_point_ptr(v);
                for(GEO::coord_index_t c = 0; c < 3; ++c) {
                    xyzmin[c] = GEO::geo_min(xyzmin[c], double(p[c]));
                    xyzmax[c] = GEO::geo_max(xyzmax[c], double(p[c]));
                    if(animate) {
                        xyzmin[c] = GEO::geo_min(xyzmin[c], double(p[c+3]));
                        xyzmax[c] = GEO::geo_max(xyzmax[c], double(p[c+3]));
                    }
                }
            } else {
                const double* p = M_in.vertices.point_ptr(v);
                for(GEO::coord_index_t c = 0; c < 3; ++c) {
                    xyzmin[c] = GEO::geo_min(xyzmin[c], p[c]);
                    xyzmax[c] = GEO::geo_max(xyzmax[c], p[c]);
                    if(animate) {
                        xyzmin[c] = GEO::geo_min(xyzmin[c], p[c+3]);
                        xyzmax[c] = GEO::geo_max(xyzmax[c], p[c+3]);
                    }
                }
            }
        }
    }

    /**
     * \brief Inverts the normals of a mesh.
     * \details In color mode, this swaps the red and the blue sides.
     */
    void invert_normals() {
        for(GEO::index_t f = 0; f < M.facets.nb(); ++f) {
            M.facets.flip(f);
        }
        M_gfx.set_mesh(&M);
    }

    /**
     * \brief Loads a mesh from a file.
     */
    void load_mesh(const std::string& filename) {
        M_gfx.set_mesh(nil);
        
        GEO::MeshIOFlags flags;
        if(GEO::CmdLine::get_arg_bool("attributes")) {
            flags.set_attribute(GEO::MESH_FACET_REGION);
            flags.set_attribute(GEO::MESH_CELL_REGION);            
        } 
        if(!GEO::mesh_load(filename, M, flags)) {
            return;
        }

        if(
            GEO::FileSystem::extension(filename) == "obj6" ||
            GEO::FileSystem::extension(filename) == "tet6"
        ) {
            GEO::Logger::out("Vorpaview")
                << "Displaying optimal transport" << std::endl;

            M_gfx.set_animate(true);
            double xyzmin[3];
            double xyzmax[3];
            get_bbox(M, xyzmin, xyzmax, true);
            glup_viewer_set_region_of_interest(
                float(xyzmin[0]), float(xyzmin[1]), float(xyzmin[2]),
                float(xyzmax[0]), float(xyzmax[1]), float(xyzmax[2])
            );
        } else {
            M_gfx.set_animate(false);            
            M.vertices.set_dimension(3);
            double xyzmin[3];
            double xyzmax[3];
            get_bbox(M, xyzmin, xyzmax, false);
            glup_viewer_set_region_of_interest(
                float(xyzmin[0]), float(xyzmin[1]), float(xyzmin[2]),
                float(xyzmax[0]), float(xyzmax[1]), float(xyzmax[2])
            );
        }
    }

    /**
     * \brief Loads a mesh from a file icon dropped into the window.
     * \details Specifed as glup_viewer_set_drag_drop_func() callback.
     */
    void dropped_file_cb(char* filename) {
        load_mesh(std::string(filename));
    }
}

/**
 * \brief Drawns and manages the graphic user interface.
 */
static void overlay() {
    ImGui::SetNextWindowPos(ImVec2(20, 20), ImGuiSetCond_Once);
    ImGui::SetNextWindowSize(ImVec2(140, 400), ImGuiSetCond_Once);
    
    ImGui::Begin("Tweaks [T]");

    if(M.vertices.dimension() >= 6) {
        ImGui::Separator();
        ImGui::Checkbox(
            "Animate [a]",
            (bool*)glup_viewer_is_enabled_ptr(GLUP_VIEWER_IDLE_REDRAW)
        );
        ImGui::SliderFloat("spd.", &OTM_speed, 1.0f, 10.0f, "%.1f");
        ImGui::SliderFloat("t.", &OTM_time, 0.0f, 1.0f, "%.2f");
    }
    
    ImGui::Separator();    
    ImGui::Checkbox("Vertices [p]", &show_vertices);
    ImGui::SliderFloat("sz.", &vertices_size, 0.1f, 5.0f, "%.1f");

    if(M.facets.nb() != 0) {
        ImGui::Separator();
        ImGui::Checkbox("Surface [S]", &show_surface);
        if(show_surface) {
            ImGui::Checkbox("mesh [m]", &show_mesh);
            ImGui::Checkbox("borders [B]", &show_surface_borders);
        }
    }

    if(M.cells.nb() != 0) {
        ImGui::Separator();
        ImGui::Checkbox("Volume [V]", &show_volume);
        if(show_volume) {
            ImGui::Checkbox("mesh [m]", &show_mesh);            
            ImGui::SliderFloat(
                "shrk.", &cells_shrink, 0.0f, 1.0f, "%.2f"
            );        
            if(!M.cells.are_simplices()) {
                ImGui::Checkbox("colored cells [C]", &show_colored_cells);  
                ImGui::Checkbox("hexes [j]", &show_hexes);
            }
        }
    }

    ImGui::Separator();
    ImGui::Checkbox(
        "Lighting [L]", &lighting
    );
    if(lighting) {
        ImGui::Checkbox(
            "edit [l]",
            (bool*)glup_viewer_is_enabled_ptr(GLUP_VIEWER_ROTATE_LIGHT)
        );
    }
    
    ImGui::Separator();
    ImGui::Checkbox(
        "Clipping [F1]", (bool*)glup_viewer_is_enabled_ptr(GLUP_VIEWER_CLIP)
    );
    if(glup_viewer_is_enabled(GLUP_VIEWER_CLIP)) {
        ImGui::Checkbox(
            "edit [F2]",
            (bool*)glup_viewer_is_enabled_ptr(GLUP_VIEWER_EDIT_CLIP)
        );
        ImGui::Checkbox(
            "fixed [F3]",
            (bool*)glup_viewer_is_enabled_ptr(GLUP_VIEWER_FIXED_CLIP)
        );
    }

    /*
    ImGui::Separator();
    ImGui::Text("Colors");
    ImGui::Checkbox("colors [c]", &show_colors);
    ImGui::Checkbox("white bkgnd [b]", &white_bg);
    */
    
    ImGui::End();
}

int main(int argc, char** argv) {

    GEO::initialize();
    GEO::Logger::instance()->set_quiet(false);

    GEO::CmdLine::import_arg_group("standard");
    GEO::CmdLine::import_arg_group("algo");
    GEO::CmdLine::import_arg_group("gfx");    

    // Default value for activating GLSL is determined by the
    // name of the executable, so that Windows users can
    // determine defaut behavior simply by changing the name
    // of the executable.
    const std::string& program_name = GEO::FileSystem::base_name(argv[0]);
    if(program_name == "vorpaview0") {
        GEO::CmdLine::set_arg("gfx:GLUP_profile","VanillaGL");
    }
    
    GEO::CmdLine::declare_arg(
        "attributes", true, "load mesh attributes"
    );

    GEO::CmdLine::declare_arg(
        "single_precision", true, "use single precision vertices (FP32)"
    );
    
    std::vector<std::string> filenames;
    if(!GEO::CmdLine::parse(argc, argv, filenames, "<filename>")) {
        return 1;
    }

    if(GEO::CmdLine::get_arg_bool("single_precision")) {
        M.vertices.set_single_precision();
    }
    
    if(filenames.size() == 1) {
        load_mesh(filenames[0]);
    } 

    if(M.cells.nb() != 0) {
        show_volume = true;
    } else if(M.facets.nb() == 0) {
        show_vertices = true;
    }

    glup_viewer_set_window_title(
        (char*) "||||||(G)||E||(O)|(G)||R||/A\\|M|||||||"
    );
    glup_viewer_set_init_func(init);
    glup_viewer_set_display_func(display);
    glup_viewer_set_overlay_func(overlay);
    glup_viewer_set_drag_drop_func(dropped_file_cb);
    
    glup_viewer_add_toggle('L', &lighting, "toggle lighting");
    glup_viewer_add_toggle(
        'a', glup_viewer_is_enabled_ptr(GLUP_VIEWER_IDLE_REDRAW), "animate"
    );
    glup_viewer_add_key_func('n', invert_normals, "invert normals");

    glup_viewer_add_toggle(
        'T', glup_viewer_is_enabled_ptr(GLUP_VIEWER_TWEAKBARS), "tweakbars"
    );

    if(GEO::CmdLine::get_arg_bool("gfx:full_screen")) {
       glup_viewer_enable(GLUP_VIEWER_FULL_SCREEN);
    }

    M_gfx.set_points_color(0.0, 1.0, 0.0);
    
    glup_viewer_main_loop(argc, argv);

    return 0;
}
