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
 *     http:www.loria.fr/~levy
 *
 *     ALICE Project
 *     LORIA, INRIA Lorraine,
 *     Campus Scientifique, BP 239
 *     54506 VANDOEUVRE LES NANCY CEDEX
 *     FRANCE
 *
 *
 * MODYFIED BY:
 *
 * Copyright (c) 2012-2015, Association Scientifique pour la Geologie et ses Applications (ASGA)
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *     * Redistributions in binary form must reproduce the above copyright
 *       notice, this list of conditions and the following disclaimer in the
 *       documentation and/or other materials provided with the distribution.
 *     * Neither the name of the <organization> nor the
 *       names of its contributors may be used to endorse or promote products
 *       derived from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL <COPYRIGHT HOLDER> BE LIABLE FOR ANY
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 *  Contacts:
 *     Arnaud.Botella@univ-lorraine.fr
 *     Antoine.Mazuyer@univ-lorraine.fr
 *     Jeanne.Pellerin@wias-berlin.de
 *
 *     http:www.ring-team.org
 *
 *     RING Project
 *     Ecole Nationale Superieure de Geologie - Georessources
 *     2 Rue du Doyen Marcel Roubault - TSA 70605
 *     54518 VANDOEUVRE-LES-NANCY
 *     FRANCE
 */
#include <ringmesh/common.h>

#ifdef RINGMESH_WITH_GRAPHICS

#include <ringmesh/geo_model.h>
#include <ringmesh/gfx.h>
#include <ringmesh/io.h>

#include <geogram_gfx/basic/GLSL.h>
#include <geogram_gfx/glut_viewer/glut_viewer.h>
#include <geogram/basic/command_line.h>
#include <geogram/basic/command_line_args.h>
#include <geogram/basic/file_system.h>
#include <geogram/basic/logger.h>

#include <algorithm>

namespace {

    RINGMesh::GeoModel GM ;
    RINGMesh::GeoModelGfx GM_gfx ;

    bool show_borders = false ;
    bool white_bg = true ;
    bool show_corners = true ;
    bool show_lines = true ;
    bool show_surface = true ;
    bool show_volume = false ;
    bool show_wells = false ;
    bool colored_cells = false ;
    bool show_voi = true ;
    bool show_colored_regions = false ;

    double shrink = 0.0 ;
    bool mesh_visible = true ;
    bool meshed_regions = false ;

    /**
     * \brief Toggles the well display
     */
    void toggle_wells()
    {
        show_wells = !show_wells ;
        GM_gfx.set_edge_regions_visibility( show_wells ) ;
    }

    /**
     * \brief Toggles the volume display
     */
    void toggle_volume()
    {
        show_volume = !show_volume ;
        GM_gfx.set_cell_regions_visibility( show_volume ) ;
    }

    /**
     * \brief Toggles the VOI display
     */
    void toggle_voi()
    {
        show_voi = !show_voi ;
        for( GEO::index_t s = 0; s < GM.nb_surfaces(); s++ ) {
            if( GM.surface( s ).is_on_voi() ) {
                GM_gfx.set_surface_visibility( s, show_voi ) ;
            }
        }
    }

    /**
     * \brief Zooms in.
     * \details Zooming factor is 1.1x.
     */
    void zoom_in()
    {
        *glut_viewer_float_ptr( GLUT_VIEWER_ZOOM ) *= 1.1f ;
    }

    /**
     * \brief Zooms out.
     * \details De-zooming factor is (1/1.1)x.
     */
    void zoom_out()
    {
        *glut_viewer_float_ptr( GLUT_VIEWER_ZOOM ) /= 1.1f ;
    }

    /**
     * \brief Toggles black or white background color.
     */
    void toggle_background()
    {
        white_bg = !white_bg ;
        if( white_bg ) {
            glut_viewer_set_background_color( 1.0, 1.0, 1.0 ) ;
            GM_gfx.set_edge_lines_color( 0.0, 0.0, 0.0 ) ;
            GM_gfx.set_mesh_surfaces_color( 0.0, 0.0, 0.0 ) ;
            GM_gfx.set_cell_mesh_regions_color( 0.0, 0.0, 0.0 ) ;
            GM_gfx.set_edge_regions_color( 0.0, 0.0, 0.0 ) ;
        } else {
            glut_viewer_set_background_color( 0.0, 0.0, 0.0 ) ;
            GM_gfx.set_edge_lines_color( 1.0, 1.0, 1.0 ) ;
            GM_gfx.set_mesh_surfaces_color( 1.0, 1.0, 1.0 ) ;
            GM_gfx.set_cell_mesh_regions_color( 1.0, 1.0, 1.0 ) ;
            GM_gfx.set_edge_regions_color( 1.0, 1.0, 1.0 ) ;
        }
    }

    /**
     * \brief Toggles mesh display.
     */
    void toggle_mesh()
    {
        mesh_visible = !mesh_visible ;
        GM_gfx.set_mesh_surfaces_visibility( mesh_visible ) ;
        GM_gfx.set_cell_mesh_regions_visibility( mesh_visible ) ;
    }

    /**
     * \brief Increases cell shrinking factor.
     */
    void inc_shrink()
    {
        shrink = std::min( shrink + 0.1, 1. ) ;
        GM_gfx.set_cell_regions_shrink( shrink ) ;
    }

    /**
     * \brief Decreases cell shrinking factor.
     */
    void dec_shrink()
    {
        shrink = std::max( shrink - 0.1, 0. ) ;
        GM_gfx.set_cell_regions_shrink( shrink ) ;
    }

    /**
     * \brief Toggles lighting / constant color mode.
     */
    void toggle_lighting()
    {
        //GM_gfx.set_lighting( !GM_gfx.get_lighting() ) ;
    }

    /**
     * \brief Initializes OpenGL objects.
     * \details Specifed as glut_viewer_set_init_func() callback.
     */
    void init()
    {
        GEO::Graphics::initialize() ;

        glut_viewer_set_background_color( 1.0, 1.0, 1.0 ) ;
        glut_viewer_add_toggle( 'T',
            glut_viewer_is_enabled_ptr( GLUT_VIEWER_TWEAKBARS ),
            "Toggle tweakbars" ) ;
        glut_viewer_add_key_func( 'b', toggle_background, "Toggle background" ) ;
        glut_viewer_add_toggle( 'B', &show_borders, "borders" ) ;
        glut_viewer_add_key_func( 'z', zoom_in, "Zoom in" ) ;
        glut_viewer_add_key_func( 'Z', zoom_out, "Zoom out" ) ;
        glut_viewer_disable( GLUT_VIEWER_TWEAKBARS ) ;
        glut_viewer_disable( GLUT_VIEWER_BACKGROUND ) ;
        glut_viewer_add_key_func( 'm', toggle_mesh, "mesh" ) ;
    }

    /**
     * \brief Draws the current mesh according to current rendering mode.
     * \details Specifed as glut_viewer_set_display_func() callback.
     */
    void display()
    {

        if( GM_gfx.geo_model() != &GM ) {
            GM_gfx.set_geo_model( GM ) ;
        }

        GLfloat shininess = 20.0f ;
        glMaterialfv( GL_FRONT_AND_BACK, GL_SHININESS, &shininess ) ;
        static float spec[4] = { 0.6f, 0.6f, 0.6f, 1.0f } ;
        glMaterialfv( GL_FRONT_AND_BACK, GL_SPECULAR, spec ) ;

        if( white_bg ) {
            GM_gfx.set_surfaces_color( 0.9f, 0.9f, 0.9f ) ;

            GM_gfx.set_edge_lines_color( 0.0, 0.0, 0.0 ) ;
            GM_gfx.set_mesh_surfaces_color( 0.0, 0.0, 0.0 ) ;
            GM_gfx.set_cell_mesh_regions_color( 0.0, 0.0, 0.0 ) ;
            GM_gfx.set_edge_regions_color( 0.0, 0.0, 0.0 ) ;
        } else {
            GM_gfx.set_surfaces_color( 0.1f, 0.1f, 0.1f ) ;

            GM_gfx.set_edge_lines_color( 1.0, 1.0, 1.0 ) ;
            GM_gfx.set_mesh_surfaces_color( 1.0, 1.0, 1.0 ) ;
            GM_gfx.set_cell_mesh_regions_color( 1.0, 1.0, 1.0 ) ;
            GM_gfx.set_edge_regions_color( 1.0, 1.0, 1.0 ) ;
        }

        if( show_corners ) {
            GM_gfx.draw_corners() ;
        }

        if( show_lines ) {
            GM_gfx.draw_lines() ;
        }

        if( show_surface ) {
            GM_gfx.draw_surfaces() ;
        }

        if( show_volume || show_wells ) {
            GM_gfx.draw_regions() ;
        }

    }

    /**
     * \brief Gets the bounding box of a mesh animation.
     * \details In animated mode, the mesh animation is stored as a mesh with
     *  6d coordinates, that correspond to the geometric location
     *  at the vertices at time 0 and at time 1.
     * \param[in] GM the GeoModel mesh.
     * \param[out] xyzmin a pointer to the three minimum coordinates
     * \param[out] xyzmax a pointer to the three maximum coordinates
     * \param[in] animate true if displaying a mesh animation
     */
    void get_bbox( const RINGMesh::GeoModel& GM, double* xyzmin, double* xyzmax )
    {
        for( GEO::index_t s = 0; s < GM.nb_surfaces(); s++ ) {
            GEO::Mesh& M = GM.surface( s ).mesh() ;
            for( GEO::index_t v = 0; v < M.vertices.nb(); ++v ) {
                const double* p = M.vertices.point_ptr( v ) ;
                for( GEO::coord_index_t c = 0; c < 3; ++c ) {
                    xyzmin[c] = GEO::geo_min( xyzmin[c], p[c] ) ;
                    xyzmax[c] = GEO::geo_max( xyzmax[c], p[c] ) ;
                }
            }
        }
    }

    /**
     * \brief Inverts the normals of a mesh.
     * \details In color mode, this swaps the red and the blue sides.
     */
    void invert_normals()
    {
        for( GEO::index_t s = 0; s < GM.nb_surfaces(); s++ ) {
            GEO::Mesh& M = GM.surface( s ).mesh() ;
            for( GEO::index_t f = 0; f < M.facets.nb(); ++f ) {
                M.facets.flip( f ) ;
            }
        }
    }

    /**
     * \brief Loads a mesh from a file.
     */
    void load_mesh()
    {
        if( !RINGMesh::model_load( GEO::CmdLine::get_arg( "model" ), GM ) ) {
            return ;
        }

        double xyzmin[3] ;
        double xyzmax[3] ;
        for( GEO::index_t c = 0; c < 3; c++ ) {
            xyzmin[c] = GEO::Numeric::max_float64() ;
            xyzmax[c] = GEO::Numeric::min_float64() ;
        }

        if( GEO::CmdLine::get_arg( "mesh" ) != "" ) {
            if( !RINGMesh::mesh_load( GEO::CmdLine::get_arg( "mesh" ), GM ) ) {
                return ;
            }
            meshed_regions = true ;
        }
        get_bbox( GM, xyzmin, xyzmax ) ;

        glut_viewer_set_region_of_interest( float( xyzmin[0] ), float( xyzmin[1] ),
            float( xyzmin[2] ), float( xyzmax[0] ), float( xyzmax[1] ),
            float( xyzmax[2] ) ) ;
    }

    void toggle_colored_cells()
    {
        colored_cells = !colored_cells ;
        if( colored_cells ) {
            GM_gfx.set_cell_regions_color_type() ;
        } else {
            GM_gfx.set_cell_regions_color( 0.9f, 0.9f, 0.9f ) ;
        }
    }

    void toggle_colored_regions()
    {
        show_colored_regions = !show_colored_regions ;

        if( show_colored_regions && meshed_regions ) {
            for( GEO::index_t r = 0; r < GM.nb_regions(); r++ ) {
                GM_gfx.set_cell_region_color( r,
                    std::fmod( GEO::Numeric::random_float32(), 1 ),
                    std::fmod( GEO::Numeric::random_float32(), 1 ),
                    std::fmod( GEO::Numeric::random_float32(), 1 ) ) ;
            }
        } else {
            GM_gfx.set_cell_regions_color( 0.9f, 0.9f, 0.9f ) ;
        }
    }

}

int main( int argc, char** argv )
{
    GEO::Logger::div( "RINGMeshView" ) ;
    GEO::Logger::out( "" ) << "Welcome to RINGMeshView !" << std::endl ;
    GEO::Logger::out( "" ) << "People working on the project in RING" << std::endl ;
    GEO::Logger::out( "" ) << "Arnaud Botella <arnaud.botella@univ-lorraine.fr> "
        << std::endl ;
    GEO::Logger::out( "" ) << "Benjamin Chauvin <benjamin.chauvin@univ-lorraine.fr> "
        << std::endl ;
    GEO::Logger::out( "" ) << "Antoine Mazuyer <antoine.mazuyer@univ-lorraine.fr> "
        << std::endl ;

    GEO::CmdLine::declare_arg( "model", "", "filename of the structural model" ) ;
    GEO::CmdLine::declare_arg( "mesh", "", "filename of the volumetric mesh" ) ;

    if( argc == 1 ) {
        GEO::CmdLine::show_usage() ;
        return 0 ;
    }

    std::vector< std::string > filenames ;
    if( !GEO::CmdLine::parse( argc, argv, filenames ) ) {
        return 1 ;
    }

    load_mesh() ;

    if( meshed_regions ) {
        toggle_volume() ;
    }

    glut_viewer_set_window_title( (char*) "RINGMeshView" ) ;
    glut_viewer_set_init_func( init ) ;
    glut_viewer_set_display_func( display ) ;
    glut_viewer_add_toggle( 'c', &show_corners, "corners" ) ;
    glut_viewer_add_toggle( 'e', &show_lines, "lines" ) ;
    glut_viewer_add_toggle( 's', &show_surface, "surface" ) ;
    glut_viewer_add_key_func( 'r', &toggle_colored_regions,
        "toggle colored regions" ) ;
    glut_viewer_add_key_func( 'v', &toggle_volume, "toggle volume" ) ;
    glut_viewer_add_key_func( 'w', &toggle_wells, "toggle wells" ) ;
    glut_viewer_add_key_func( 'V', toggle_voi, "toggle VOI" ) ;
    glut_viewer_add_key_func( 'L', toggle_lighting, "toggle lighting" ) ;
    glut_viewer_add_key_func( 'n', invert_normals, "invert normals" ) ;
    glut_viewer_add_key_func( 'X', dec_shrink, "unshrink cells" ) ;
    glut_viewer_add_key_func( 'x', inc_shrink, "shrink cells" ) ;
    glut_viewer_add_key_func( 'C', toggle_colored_cells, "toggle colored cells" ) ;

    if( GEO::CmdLine::get_arg_bool( "gfx:full_screen" ) ) {
        glut_viewer_enable( GLUT_VIEWER_FULL_SCREEN ) ;
    }

//    GM_gfx.set_points_color( 0.0, 1.0, 0.0 ) ;

    glut_viewer_main_loop( argc, argv ) ;

    // Note: when 'q' is pressed, exit() is called
    // because there is no simple way of exiting from
    // glut's event loop, therefore this line is not
    // reached and memory is not properly freed on exit.
    // TODO: add a function in freeglut to exit event loop.

    return 0 ;
}

#else
#include <geogram/basic/logger.h>
int main() {
    GEO::Logger::out("RINGMeshView") 
        << "To compile RINGMesh viewer you need to configure " 
        << "the project with the RINGMESH_WITH_GRAPHICS option ON"
        << std::endl ;
    return 0 ;
}
#endif
