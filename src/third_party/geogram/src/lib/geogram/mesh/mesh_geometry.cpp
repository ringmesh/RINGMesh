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

#include <geogram/mesh/mesh_geometry.h>
#include <geogram/delaunay/LFS.h>
#include <geogram/voronoi/CVT.h>
#include <geogram/basic/geometry.h>

namespace {

    using namespace GEO;

    /**
     * \brief Computes a sizing field using local feature size
     * \details The sizing field is stored into the vertex weights
     *  of the mesh
     * \param[in] M the mesh
     * \param[in] LFS the local feature size
     * \param[in] gradation power to be applied to the sizing field
     */
    void compute_sizing_field_lfs(
        Mesh& M, const LocalFeatureSize& LFS, double gradation
    ) {
        Logger::out("LFS") << "Computing sizing field" << std::endl;
        for(index_t v = 0; v < M.nb_vertices(); v++) {
            double lfs2 = LFS.squared_lfs(M.vertex_ptr(v));
            double weight = pow(lfs2, -2.0 * gradation);
            M.set_weight(v, weight);
        }
    }
}

/****************************************************************************/

namespace GEO {

    namespace Geom {

        double mesh_dihedral_angle(const Mesh& M, index_t c) {
            geo_debug_assert(M.is_triangulated());
            index_t f1 = c/3;
            signed_index_t sf2 = M.corner_adjacent_facet(c);
            geo_debug_assert(sf2 >= 0);
            index_t f2 = index_t(sf2);
            vec3 n1 = normalize(mesh_facet_normal(M,f1));
            vec3 n2 = normalize(mesh_facet_normal(M,f2));
            double sign = 1.0;
            if(dot(cross(n1,n2),mesh_corner_vector(M,c)) > 0.0) {
                sign = -1.0;
            }
            return sign*acos(dot(n1, n2));
        }

        double mesh_area(const Mesh& M, coord_index_t dim) {
            double result = 0.0;
            for(index_t f = 0; f < M.nb_facets(); f++) {
                result += M.facet_area(f, dim);
            }
            return result;
        }
    }

    void compute_normals(Mesh& M) {
        if(M.dimension() < 6) {
            M.set_dimension(6);
        } else {
            for(index_t i = 0; i < M.nb_vertices(); i++) {
                Geom::mesh_vertex_normal_ref(M, i) = vec3(0.0, 0.0, 0.0);
            }
        }
        for(index_t f = 0; f < M.nb_facets(); f++) {
            vec3 N = Geom::mesh_facet_normal(M, f);
            for(
                index_t corner = M.facet_begin(f);
                corner < M.facet_end(f); corner++
            ) {
                index_t v = M.corner_vertex_index(corner);
                Geom::mesh_vertex_normal_ref(M, v) += N;
            }
        }
        for(index_t i = 0; i < M.nb_vertices(); i++) {
            Geom::mesh_vertex_normal_ref(M, i) = normalize(
                Geom::mesh_vertex_normal(M, i)
            );
        }
    }

    void simple_Laplacian_smooth(Mesh& M, index_t nb_iter, bool normals_only) {
        geo_assert(M.dimension() >= 6);
        std::vector<vec3> p(M.nb_vertices());
        std::vector<double> c(M.nb_vertices());

        for(index_t k = 0; k < nb_iter; k++) {
            p.assign(M.nb_vertices(), vec3(0.0, 0.0, 0.0));
            c.assign(M.nb_vertices(), 0);
            for(index_t f = 0; f < M.nb_facets(); f++) {
                index_t b = M.facet_begin(f);
                index_t e = M.facet_end(f);
                for(index_t c1 = b; c1 != e; c1++) {
                    index_t c2 = (c1 == e - 1) ? b : c1 + 1;
                    index_t v1 = M.corner_vertex_index(c1);
                    index_t v2 = M.corner_vertex_index(c2);
                    if(v1 < v2) {
                        double a = 1.0;
                        c[v1] += a;
                        c[v2] += a;
                        if(normals_only) {
                            p[v1] += a * Geom::mesh_vertex_normal(M, v2);
                            p[v2] += a * Geom::mesh_vertex_normal(M, v1);
                        } else {
                            p[v1] += a * Geom::mesh_vertex(M, v2);
                            p[v2] += a * Geom::mesh_vertex(M, v1);
                        }
                    }
                }
            }
            for(index_t v = 0; v < M.nb_vertices(); v++) {
                if(normals_only) {
                    double l = length(p[v]);
                    if(l > 1e-30) {
                        Geom::mesh_vertex_normal_ref(M, v) = (1.0 / l) * p[v];
                    }
                } else {
                    Geom::mesh_vertex_ref(M, v) = 1.0 / c[v] * p[v];
                }
            }
        }
        if(!normals_only) {
            compute_normals(M);
        }
    }

    void get_bbox(const Mesh& M, double* xyzmin, double* xyzmax) {
        geo_assert(M.dimension() >= 3);
        for(index_t c = 0; c < 3; c++) {
            xyzmin[c] = Numeric::max_float64();
            xyzmax[c] = Numeric::min_float64();
        }
        for(index_t v = 0; v < M.nb_vertices(); v++) {
            const double* p = M.vertex_ptr(v);
            for(index_t c = 0; c < 3; c++) {
                xyzmin[c] = geo_min(xyzmin[c], p[c]);
                xyzmax[c] = geo_max(xyzmax[c], p[c]);
            }
        }
    }

    void get_bbox(
        const SinglePrecisionMesh& M, double* xyzmin, double* xyzmax
    ) {
        geo_assert(M.dimension() >= 3);
        for(index_t c = 0; c < 3; c++) {
            xyzmin[c] = Numeric::max_float64();
            xyzmax[c] = Numeric::min_float64();
        }
        for(index_t v = 0; v < M.nb_vertices(); v++) {
            const float* p = M.vertex_ptr(v);
            for(index_t c = 0; c < 3; c++) {
                xyzmin[c] = geo_min(xyzmin[c], double(p[c]));
                xyzmax[c] = geo_max(xyzmax[c], double(p[c]));
            }
        }
    }

    double bbox_diagonal(const Mesh& M) {
        geo_assert(M.dimension() >= 3);
        double xyzmin[3];
        double xyzmax[3];
        get_bbox(M, xyzmin, xyzmax);
        return ::sqrt(
            geo_sqr(xyzmax[0] - xyzmin[0]) +
            geo_sqr(xyzmax[1] - xyzmin[1]) +
            geo_sqr(xyzmax[2] - xyzmin[2])
        );
    }

    void set_anisotropy(Mesh& M, double s) {
        if(M.dimension() < 6) {
            compute_normals(M);
        }
        if(s == 0.0) {
            unset_anisotropy(M);
            return;
        }
        s *= bbox_diagonal(M);
        for(index_t i = 0; i < M.nb_vertices(); i++) {
            Geom::mesh_vertex_normal_ref(M, i) =
                s * normalize(Geom::mesh_vertex_normal(M, i));
        }
    }

    void unset_anisotropy(Mesh& M) {
        for(index_t i = 0; i < M.nb_vertices(); i++) {
            Geom::mesh_vertex_normal_ref(M, i) = normalize(
                Geom::mesh_vertex_normal(M, i)
            );
        }
    }

    void compute_sizing_field(
        Mesh& M, double gradation, index_t nb_lfs_samples
    ) {
        if(nb_lfs_samples != 0) {
            Logger::out("LFS") << "Sampling surface" << std::endl;
            CentroidalVoronoiTesselation CVT(&M, 3);
            CVT.compute_initial_sampling(nb_lfs_samples);
            CVT.Lloyd_iterations(5);
            CVT.Newton_iterations(10);
            Logger::out("LFS") << "Computing medial axis" << std::endl;
            LocalFeatureSize LFS(CVT.nb_points(), CVT.embedding(0));
            compute_sizing_field_lfs(M, LFS, gradation);
        } else {
            if(M.dimension() == 3) {
                LocalFeatureSize LFS(M.nb_vertices(), M.vertex_ptr(0));
                compute_sizing_field_lfs(M, LFS, gradation);
            } else {
                std::vector<double> pts;
                pts.reserve(M.nb_vertices() * 3);
                for(index_t v = 0; v < M.nb_vertices(); v++) {
                    pts.push_back(M.vertex_ptr(v)[0]);
                    pts.push_back(M.vertex_ptr(v)[1]);
                    pts.push_back(M.vertex_ptr(v)[2]);
                }
                LocalFeatureSize LFS(M.nb_vertices(), &pts[0]);
                compute_sizing_field_lfs(M, LFS, gradation);
            }
        }
    }

    void normalize_embedding_area(Mesh& M) {
        if(M.dimension() == 3) {
            M.remove_vertices_weights();
            return;
        }
        std::vector<double> area3d(M.nb_vertices(), 0.0);
        std::vector<double> areaNd(M.nb_vertices(), 0.0);
        for(index_t f = 0; f < M.nb_facets(); f++) {
            double A3d = M.facet_area(f, 3);
            double ANd = M.facet_area(f);
            for(index_t c = M.facet_begin(f); c < M.facet_end(f); c++) {
                index_t v = M.corner_vertex_index(c);
                area3d[v] += A3d;
                areaNd[v] += ANd;
            }
        }
        for(index_t v = 0; v < M.nb_vertices(); v++) {
            double A3d = area3d[v];
            double ANd = areaNd[v];
            ANd = geo_max(ANd, 1e-6);
            double w = ::pow(A3d / ANd, 2.0);
            M.set_weight(v, w);
        }
    }

    double mesh_cell_volume(
        const Mesh& M, index_t c
    ) {
        // Only implemented for tetrahedra
        // TODO: other cell types
        geo_assert(M.cell_type(c) == MESH_TET);
        geo_assert(M.dimension() >= 3);
        const double* p0 = M.vertex_ptr(M.tet_vertex_index(c,0));
        const double* p1 = M.vertex_ptr(M.tet_vertex_index(c,1));
        const double* p2 = M.vertex_ptr(M.tet_vertex_index(c,2));
        const double* p3 = M.vertex_ptr(M.tet_vertex_index(c,3));
        return ::fabs(Geom::tetra_signed_volume(p0,p1,p2,p3));
    }

    
    double mesh_cells_volume(const Mesh& M) {
        double result = 0.0;
        for(index_t c=0; c<M.nb_cells(); ++c) {
            result += mesh_cell_volume(M,c);
        }
        return result;
    }
    
}

