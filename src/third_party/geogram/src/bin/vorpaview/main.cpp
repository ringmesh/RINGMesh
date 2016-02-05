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
#include <geogram_gfx/glut_viewer/glut_viewer.h>
#include <geogram/mesh/mesh.h>
#include <geogram/mesh/mesh_io.h>
#include <geogram_gfx/mesh/mesh_gfx.h>
#include <geogram/basic/command_line.h>
#include <geogram/basic/command_line_args.h>
#include <geogram/basic/file_system.h>
#include <geogram/basic/logger.h>

namespace {

    GEO::Mesh M(3);
    GEO::MeshGfx M_gfx;

    bool hexes = true;
    bool show_colors   = true;
    bool show_borders  = false;
    bool white_bg      = true;
    bool show_surface  = true;
    bool show_volume   = false;
    bool show_vertices = false;
    bool colored_cells = false;

    /**
     * \brief Zooms in.
     * \details Zooming factor is 1.1x.
     */
    void zoom_in() {
        *glut_viewer_float_ptr(GLUT_VIEWER_ZOOM) *= 1.1f;
    }

    /**
     * \brief Zooms out.
     * \details De-zooming factor is (1/1.1)x.
     */
    void zoom_out() {
        *glut_viewer_float_ptr(GLUT_VIEWER_ZOOM) /= 1.1f;
    }

    /**
     * \brief Toggles black or white background color.
     */
    void toggle_background() {
        white_bg = !white_bg;
        if(white_bg) {
            glut_viewer_set_background_color(1.0, 1.0, 1.0);
            M_gfx.set_mesh_color(0.0, 0.0, 0.0);
        } else {
            glut_viewer_set_background_color(0.0, 0.0, 0.0);
            M_gfx.set_mesh_color(1.0, 1.0, 1.0);
        }
    }

    /**
     * \brief Toggles mesh display.
     */
    void toggle_mesh() {
        M_gfx.set_show_mesh(!M_gfx.get_show_mesh());
    }

    /**
     * \brief Increases cell shrinking factor.
     */
    void inc_shrink() {
        M_gfx.set_shrink(M_gfx.get_shrink() + 0.1);
    }

    /**
     * \brief Decreases cell shrinking factor.
     */
    void dec_shrink() {
        M_gfx.set_shrink(M_gfx.get_shrink() - 0.1);
    }
    
    /**
     * \brief Toggles color / BW display.
     */
    void toggle_colors() {
        show_colors = !show_colors;
    }

    /**
     * \brief Toggles lighting / constant color mode.
     */
    void toggle_lighting() {
        M_gfx.set_lighting(!M_gfx.get_lighting());
    }

    /**
     * \brief Cycles between three possible points size.
     */
    void cycle_points_size() {
        if(M_gfx.get_points_size() == 2.0f) {
            M_gfx.set_points_size(0.1f);
        } else if(M_gfx.get_points_size() == 0.1f) {
            M_gfx.set_points_size(1.0f);            
        } else {
            M_gfx.set_points_size(2.0f);            
        }
    }
    
    /**
     * \brief Increments the time of the Optimal Transport animation.
     */
    void increment_time() {
        M_gfx.set_time(M_gfx.get_time() + 0.05);
    }

    /**
     * \brief Decrements the time of the Optimal Transport animation.
     */
    void decrement_time() {
        M_gfx.set_time(M_gfx.get_time() - 0.05);        
    }
    
    /**
     * \brief Initializes OpenGL objects.
     * \details Specifed as glut_viewer_set_init_func() callback.
     */
    void init() {
        GEO::Graphics::initialize();
        
        glut_viewer_set_background_color(1.0, 1.0, 1.0);
        glut_viewer_add_toggle(
            'T', glut_viewer_is_enabled_ptr(GLUT_VIEWER_TWEAKBARS),
            "Toggle tweakbars"
        );
        glut_viewer_add_key_func('b', toggle_background, "Toggle background");
        glut_viewer_add_key_func('c', toggle_colors, "colors");
        glut_viewer_add_toggle('B', &show_borders, "borders");
        glut_viewer_add_key_func('z', zoom_in, "Zoom in");
        glut_viewer_add_key_func('Z', zoom_out, "Zoom out");
        glut_viewer_add_key_func('r', decrement_time, "Decrement time");
        glut_viewer_add_key_func('t', increment_time, "Increment time");
        glut_viewer_disable(GLUT_VIEWER_TWEAKBARS);
        glut_viewer_disable(GLUT_VIEWER_BACKGROUND);
        glut_viewer_add_key_func('m', toggle_mesh, "mesh");
        glut_viewer_add_toggle('j',  &hexes, "hexes");
    }

    /**
     * \brief Draws the current mesh according to current rendering mode.
     * \details Specifed as glut_viewer_set_display_func() callback.
     */
    void display() {

        if(M_gfx.mesh() != &M) {
            M_gfx.set_mesh(&M);
        }
        
        if(show_borders) {
            M_gfx.draw_surface_borders();
        }

        M_gfx.set_draw_cells(GEO::MESH_HEX, hexes);
        
        GLfloat shininess = 20.0f;
        glMaterialfv(GL_FRONT_AND_BACK, GL_SHININESS, &shininess);
        static float spec[4] = {0.6f, 0.6f, 0.6f, 1.0f};
        glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, spec);
        
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

        if(show_vertices) {
            M_gfx.draw_vertices();
        }
        
        if(show_surface) {
            M_gfx.draw_surface();
        }

        if(show_volume) {
            M_gfx.draw_volume();
        }

    }

    /**
     * \brief Gets the bounding box of a mesh animation.
     * \details In animated mode, the mesh animation is stored as a mesh with
     *  6d coordinates, that correspond to the geometric location
     *  at the vertices at time 0 and at time 1.
     * \param[in] M the mesh
     * \param[out] xyzmin a pointer to the three minimum coordinates
     * \param[out] xyzmax a pointer to the three maximum coordinates
     * \param[in] animate true if displaying a mesh animation
     */
    void get_bbox(
        const GEO::Mesh& M, double* xyzmin, double* xyzmax,
        bool animate
    ) {
        geo_assert(M.vertices.dimension() >= GEO::index_t(animate ? 6 : 3));
        for(GEO::index_t c = 0; c < 3; c++) {
            xyzmin[c] = GEO::Numeric::max_float64();
            xyzmax[c] = GEO::Numeric::min_float64();
        }

        for(GEO::index_t v = 0; v < M.vertices.nb(); ++v) {
            if(M.vertices.single_precision()) {
                const float* p = M.vertices.single_precision_point_ptr(v);
                for(GEO::coord_index_t c = 0; c < 3; ++c) {
                    xyzmin[c] = GEO::geo_min(xyzmin[c], double(p[c]));
                    xyzmax[c] = GEO::geo_max(xyzmax[c], double(p[c]));
                    if(animate) {
                        xyzmin[c] = GEO::geo_min(xyzmin[c], double(p[c+3]));
                        xyzmax[c] = GEO::geo_max(xyzmax[c], double(p[c+3]));
                    }
                }
            } else {
                const double* p = M.vertices.point_ptr(v);
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
            glut_viewer_set_region_of_interest(
                float(xyzmin[0]), float(xyzmin[1]), float(xyzmin[2]),
                float(xyzmax[0]), float(xyzmax[1]), float(xyzmax[2])
            );
        } else {
            M.vertices.set_dimension(3);
            double xyzmin[3];
            double xyzmax[3];
            get_bbox(M, xyzmin, xyzmax, false);
            glut_viewer_set_region_of_interest(
                float(xyzmin[0]), float(xyzmin[1]), float(xyzmin[2]),
                float(xyzmax[0]), float(xyzmax[1]), float(xyzmax[2])
            );
        }
    }

    /**
     * \brief Loads a mesh from a file icon dropped into the window.
     * \details Specifed as glut_viewer_set_drag_drop_func() callback.
     */
    void dropped_file_cb(char* filename) {
        load_mesh(std::string(filename));
    }


    void toggle_colored_cells() {
        colored_cells = !colored_cells;
        if(colored_cells) {
            M_gfx.set_cells_colors_by_type();
        } else {
            M_gfx.set_cells_color(0.9f, 0.9f, 0.9f);
        }
    }

    void toggle_regions() {
//       M_gfx.set_show_regions(!M_gfx.get_show_regions());
    }

    GEO::index_t cur_RVD=0;
    void load_RVD() {
        load_mesh("RVD_" + GEO::String::to_string(cur_RVD) + ".meshb");
        M_gfx.set_mesh(&M);
    }
    
    void load_RVD_inc() {
        ++cur_RVD;
        load_RVD();
    }

    void load_RVD_dec() {
        --cur_RVD;
        load_RVD();
    }
    
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
        GEO::CmdLine::set_arg("gfx:GLSL",false);
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

    glut_viewer_set_window_title(
        (char*) "||||||(G)||E||(O)|(G)||R||/A\\|M|||||||"
    );
    glut_viewer_set_init_func(init);
    glut_viewer_set_display_func(display);
    glut_viewer_set_drag_drop_func(dropped_file_cb);
    glut_viewer_add_toggle('V', &show_volume, "volume");
    glut_viewer_add_toggle('p', &show_vertices, "vertices");
    glut_viewer_add_key_func('P', &cycle_points_size, "change points size");    
    glut_viewer_add_toggle('S', &show_surface, "surface");
    glut_viewer_add_key_func('L', toggle_lighting, "toggle lighting");
    glut_viewer_add_key_func('n', invert_normals, "invert normals");
    glut_viewer_add_key_func('x', dec_shrink, "unshrink cells");
    glut_viewer_add_key_func('w', inc_shrink, "shrink cells");
    glut_viewer_add_key_func('C', toggle_colored_cells, "toggle colored cells");
    glut_viewer_add_key_func('R', toggle_regions, "toggle regions");    

    glut_viewer_add_key_func('D', load_RVD_inc, "load next RVD");
    glut_viewer_add_key_func('F', load_RVD_dec, "load prev RVD");    
    
    if(GEO::CmdLine::get_arg_bool("gfx:full_screen")) {
       glut_viewer_enable(GLUT_VIEWER_FULL_SCREEN);
    }

    M_gfx.set_points_color(0.0, 1.0, 0.0);
    
    glut_viewer_main_loop(argc, argv);

    // Note: when 'q' is pressed, exit() is called
    // because there is no simple way of exiting from
    // glut's event loop, therefore this line is not
    // reached and memory is not properly freed on exit.
    // TODO: add a function in freeglut to exit event loop.

    return 0;
}
