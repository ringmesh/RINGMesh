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

#include <geogram/voronoi/generic_RVD_cell.h>
#include <geogram/numerics/predicates.h>

namespace GEOGen {

    index_t ConvexCell::plus1mod3_[3] = {1, 2, 0};
    index_t ConvexCell::minus1mod3_[3] = {2, 0, 1};

    std::ostream& ConvexCell::show_stats(std::ostream& os) const {
        unsigned int nb_free = 0;
        unsigned int nb_conflict = 0;
        unsigned int nb_used = 0;
        for(unsigned int t = 0; t < max_t(); t++) {
            switch(triangles_[t].status_) {
                case TRI_IS_FREE:
                    nb_free++;
                    break;
                case TRI_IS_USED:
                    nb_used++;
                    break;
                case TRI_IS_CONFLICT:
                    nb_conflict++;
                    break;
            }
        }
        return os << "Nb tot = " << max_t()
            << " free=" << nb_free
            << " used=" << nb_used
            << " conflict=" << nb_conflict
            << std::endl;
    }

    Sign ConvexCell::side_exact(
        const Mesh* mesh, const Delaunay* delaunay,
        const GEOGen::Vertex& q,
        const double* pi, const double* pj,
        coord_index_t dim,
        bool symbolic_is_surface
    ) const {

        switch(q.sym().nb_boundary_facets()) {
            case 0:
            {
                // The point q is the intersection between
                //   three bisectors [pi b0], [pi b1] and [pi b2]
                // (and a tet [q0 q1 q2 q3])

                unsigned int b0 = q.sym().bisector(0);
                unsigned int b1 = q.sym().bisector(1);
                unsigned int b2 = q.sym().bisector(2);

                if(dim == 3) {
                    // 3d is a special case for side4()
                    //   (intrinsic dim == ambient dim)
                    // therefore embedding tet q0,q1,q2,q3 is not needed.
                    return GEO::PCK::side4_3d_SOS(
                        pi,
                        delaunay->vertex_ptr(b0),
                        delaunay->vertex_ptr(b1),
                        delaunay->vertex_ptr(b2),
                        pj
                    );
                } else {
                    geo_debug_assert(cell_id() >= 0);
                    index_t t = index_t(cell_id());
                    return GEO::PCK::side4_SOS(
                        pi,
                        delaunay->vertex_ptr(b0),
                        delaunay->vertex_ptr(b1),
                        delaunay->vertex_ptr(b2),
                        pj,
                        mesh->vertices.point_ptr(mesh->cells.tet_vertex(t, 0)),
                        mesh->vertices.point_ptr(mesh->cells.tet_vertex(t, 1)),
                        mesh->vertices.point_ptr(mesh->cells.tet_vertex(t, 2)),
                        mesh->vertices.point_ptr(mesh->cells.tet_vertex(t, 3)),
                        dim
                    );
                }
            } break;

            case 1:
            {
                // The point q is the intersection between
                //   a facet (f0,f1,f2) of the surface and two
                //   bisectors [pi b0] and [pi b1].

                unsigned int b0 = q.sym().bisector(0);
                unsigned int b1 = q.sym().bisector(1);
                unsigned int f = q.sym().boundary_facet(0);

                if(symbolic_is_surface) {
                    index_t c = mesh->facets.corners_begin(f);
                    const double* q0 = mesh->vertices.point_ptr(
                        mesh->facet_corners.vertex(c)
                    );
                    const double* q1 = mesh->vertices.point_ptr(
                        mesh->facet_corners.vertex(c+1)
                    );
                    const double* q2 = mesh->vertices.point_ptr(
                        mesh->facet_corners.vertex(c+2)
                    );
                                                                                
                    return GEO::PCK::side3_SOS(
                        pi,
                        delaunay->vertex_ptr(b0),
                        delaunay->vertex_ptr(b1),
                        pj,
                        q0, q1, q2, dim
                    );
                    
                } else {
                    index_t t = f / 4;
                    index_t lf = f % 4;
                    index_t j0 = mesh->cells.tet_vertex(
                        t, GEO::MeshCells::local_tet_facet_vertex_index(lf, 0)
                    );
                    index_t j1 = mesh->cells.tet_vertex(
                        t, GEO::MeshCells::local_tet_facet_vertex_index(lf, 1)
                    );
                    index_t j2 = mesh->cells.tet_vertex(
                        t, GEO::MeshCells::local_tet_facet_vertex_index(lf, 2)
                    );
                    
                    return GEO::PCK::side3_SOS(
                        pi,
                        delaunay->vertex_ptr(b0),
                        delaunay->vertex_ptr(b1),
                        pj,
                        mesh->vertices.point_ptr(j0),
                        mesh->vertices.point_ptr(j1),
                        mesh->vertices.point_ptr(j2),
                        dim
                    );
                }
            } break;

            case 2:
            {
                // The point q is the intersection between
                //   two facets of the surface (i.e. an edge [e0 e1])
                //   and one bisector [pi b0].
                // i.e. it's a vertex of the surface.
                index_t b0 = q.sym().bisector(0);
                index_t e0, e1;
                q.sym().get_boundary_edge(e0, e1);
                return GEO::PCK::side2_SOS(
                    pi, delaunay->vertex_ptr(b0), pj,
                    mesh->vertices.point_ptr(e0),
                    mesh->vertices.point_ptr(e1),
                    dim
                );
            } break;

            case 3:
            {
                // The point q is the intersection between
                //   three facets of the surface
                //   (i.e. a vertex v0 of the surface).
                unsigned int v0 = q.sym().get_boundary_vertex();
                return GEO::PCK::side1_SOS(
                    pi, pj, mesh->vertices.point_ptr(v0), dim
                );
            } break;
        }
        geo_assert_not_reached;
        return Sign(0);
    }

    void ConvexCell::initialize_from_mesh_tetrahedron(
        const Mesh* mesh, index_t t, bool symbolic
    ) {
        clear();

        index_t v0 = mesh->cells.tet_vertex(t, 0);
        index_t v1 = mesh->cells.tet_vertex(t, 1);
        index_t v2 = mesh->cells.tet_vertex(t, 2);
        index_t v3 = mesh->cells.tet_vertex(t, 3);

        signed_index_t t0 = signed_index_t(mesh->cells.tet_adjacent(t, 0));
        signed_index_t t1 = signed_index_t(mesh->cells.tet_adjacent(t, 1));
        signed_index_t t2 = signed_index_t(mesh->cells.tet_adjacent(t, 2));
        signed_index_t t3 = signed_index_t(mesh->cells.tet_adjacent(t, 3));

        create_vertex();
        create_vertex();
        create_vertex();
        create_vertex();

        set_cell_id(signed_index_t(t));

        set_vertex_id(0, (t0 == signed_index_t(GEO::NO_CELL)) ? 0 : -t0 - 1);
        set_vertex_id(1, (t1 == signed_index_t(GEO::NO_CELL)) ? 0 : -t1 - 1);
        set_vertex_id(2, (t2 == signed_index_t(GEO::NO_CELL)) ? 0 : -t2 - 1);
        set_vertex_id(3, (t3 == signed_index_t(GEO::NO_CELL)) ? 0 : -t3 - 1);

        create_triangle(mesh->vertices.point_ptr(v0), 2, 1, 3, 2, 1, 3);
        create_triangle(mesh->vertices.point_ptr(v1), 3, 0, 2, 3, 0, 2);
        create_triangle(mesh->vertices.point_ptr(v2), 0, 3, 1, 0, 3, 1);
        create_triangle(mesh->vertices.point_ptr(v3), 2, 0, 1, 2, 0, 1);

        if(symbolic) {

            index_t f0 = global_facet_id(mesh, t, 0);
            index_t f1 = global_facet_id(mesh, t, 1);
            index_t f2 = global_facet_id(mesh, t, 2);
            index_t f3 = global_facet_id(mesh, t, 3);

            triangle_dual(0).sym().set_boundary_vertex(v0);
            triangle_dual(0).sym().add_boundary_facet(f1);
            triangle_dual(0).sym().add_boundary_facet(f2);
            triangle_dual(0).sym().add_boundary_facet(f3);

            triangle_dual(1).sym().set_boundary_vertex(v1);
            triangle_dual(1).sym().add_boundary_facet(f2);
            triangle_dual(1).sym().add_boundary_facet(f3);
            triangle_dual(1).sym().add_boundary_facet(f0);

            triangle_dual(2).sym().set_boundary_vertex(v2);
            triangle_dual(2).sym().add_boundary_facet(f3);
            triangle_dual(2).sym().add_boundary_facet(f0);
            triangle_dual(2).sym().add_boundary_facet(f1);

            triangle_dual(3).sym().set_boundary_vertex(v3);
            triangle_dual(3).sym().add_boundary_facet(f0);
            triangle_dual(3).sym().add_boundary_facet(f1);
            triangle_dual(3).sym().add_boundary_facet(f2);
        }
    }
}

