/*
 *  OGF/Graphite: Geometry and Graphics Programming Library + Utilities
 *  Copyright (C) 2000-2015 INRIA - Project ALICE
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 *
 *  If you modify this software, you should include a notice giving the
 *  name of the person performing the modification, the date of modification,
 *  and the reason for such modification.
 *
 *  Contact for Graphite: Bruno Levy - Bruno.Levy@inria.fr
 *  Contact for this Plugin: Nicolas Ray - nicolas.ray@inria.fr
 *
 *     Project ALICE
 *     LORIA, INRIA Lorraine,
 *     Campus Scientifique, BP 239
 *     54506 VANDOEUVRE LES NANCY CEDEX
 *     FRANCE
 *
 *  Note that the GNU General Public License does not permit incorporating
 *  the Software into proprietary programs.
 *
 * As an exception to the GPL, Graphite can be linked with the following
 * (non-GPL) libraries:
 *     Qt, tetgen, SuperLU, WildMagic and CGAL
 */


#ifndef H_HEXDOM_ALGO_PGP_OPT_H
#define H_HEXDOM_ALGO_PGP_OPT_H

#include <exploragram/basic/common.h>
#include <exploragram/hexdom/basic.h> 
#include <geogram/mesh/mesh.h>

namespace GEO {
    
    class EXPLORAGRAM_API PGPopt {
    public:

        // what we need to create the parametric space U
        PGPopt(Mesh* p_m);

        void optimize_corr(double max_corr_prop = .35);
        void optimize_PGP();


        vec3 wish_angle_edge_geom(index_t e, bool inv);
        vec3 wish_angle_corr(index_t e, bool inv);
        vec3 wish_angle(index_t e, bool inv);
        bool is_PGP_singular(index_t c, index_t lf);

        bool face_is_resp(index_t c, index_t lf);

        // export functions
        vec3 change_tet_basis(vec3 in, vec3 P[4], vec3 P_img[4]);


        bool in_tet(vec3 test, vec3 P[4], double eps = .00001);

        bool tet_is_PGP_singular_fct(index_t t);


        void make_compatible(index_t  ref, index_t v);

        /**
         * All informations required to produce a hex-dom mesh are exported
         * into a 2D non-manifold mesh with uv coordinates
         *
         * This new mesh is the only input of the following steps in the
         * hex dom generation pipeline
         */
        void export_boundary_with_uv(Mesh* hex, Attribute<vec2>& uv, Attribute<index_t> &singtri);


        void get_grid_vertices(index_t t, std::vector<vec3>& psetX, std::vector<vec3>& psetU, bool dual);
        void export_hexes(Mesh* hex);

        void export_points(Mesh* pts);


        //// cubecover
        index_t get_the_only_non_zero_lc(index_t c, index_t cf, Attribute<index_t>& CCedgeid);
        index_t number_edges_on_cut_graph(Attribute<index_t>& CCedgeid);
        
        
        //bool constraint_boundary is set to false only to compute the correction one form
        void cubcover(bool compute_only_corr =false);




        index_t edge_from_vertices(index_t v0, index_t v1, bool &inv) {
            inv = v0 > v1;
            if (inv) std::swap(v0, v1);

            index_t start = v2e[v0];
            index_t end = v0 < m->vertices.nb() - 1 ? v2e[v0 + 1] : m->edges.nb();
            for (index_t e = start; e < end; e++) {
                if (m->edges.vertex(e, 1) == v1) return e;
            }
            geo_assert_not_reached;
         
        }


        void snap_U_to_round(double eps = 0.05) {
            FOR(v, m->vertices.nb()) FOR(d, 3) {
                if (std::abs(U[v][d] - round(U[v][d])) < eps) {
                    U[v][d] = round(U[v][d]);
                }
            }
        }



        Mesh* m;
        Attribute<vec3> U;
        Attribute<mat3> B;

        // A PGP solution is not only a mesh + attrib: it also requires these datas
       
        Attribute<vec3> corr;
        Attribute<vec3i> tij;

        vector<index_t> v2e;
        vector<vector<index_t> > v2eopp;
    };



}
#endif
