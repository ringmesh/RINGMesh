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

#include <geogram/basic/common.h>
#include <geogram/basic/logger.h>
#include <geogram/basic/command_line.h>
#include <geogram/basic/command_line_args.h>
#include <geogram/basic/line_stream.h>
#include <geogram/basic/stopwatch.h>
#include <geogram/basic/geometry.h>
#include <geogram/voronoi/generic_RVD_cell.h>
#include <geogram/mesh/mesh_halfedges.h>
#include <geogram/mesh/mesh_io.h>
#include <geogram/mesh/mesh_topology.h>

namespace {

    using namespace GEO;
    using GEOGen::ConvexCell;

    /**
     * \brief Creates a mesh from a cube.
     * \details Creates a mesh that represents a cube centered on the
     * origin, with a edge length of 2.
     * \param[out] M the resulting mesh
     */
    void initialize_mesh_with_box(Mesh& M) {
        M.clear();;;
        M.vertices.set_dimension(3);

        M.vertices.create_vertex(vec3(-1, -1, -1).data());
        M.vertices.create_vertex(vec3(-1, -1, 1).data());
        M.vertices.create_vertex(vec3(-1, 1, -1).data());
        M.vertices.create_vertex(vec3(-1, 1, 1).data());
        M.vertices.create_vertex(vec3(1, -1, -1).data());
        M.vertices.create_vertex(vec3(1, -1, 1).data());
        M.vertices.create_vertex(vec3(1, 1, -1).data());
        M.vertices.create_vertex(vec3(1, 1, 1).data());

        M.facets.create_quad(7,6,2,3);
        M.facets.create_quad(1,3,2,0);
        M.facets.create_quad(5,7,3,1);
        M.facets.create_quad(4,6,7,5);
        M.facets.create_quad(4,5,1,0);
        M.facets.create_quad(6,4,0,2);

        M.facets.connect();
    }

    /**
     * \brief Copies a Mesh into a ConvexCell
     * \details The surface mesh in \p M represents the boundary
     *  of the ConvexCell \p M.
     * \param[out] C a reference to the ConvexCell
     * \param[in] M a reference to the input Mesh
     */
    void initialize_convex_cell_with_mesh(
        ConvexCell& C, Mesh& M
    ) {
        for(index_t f = 0; f < M.facets.nb(); ++f) {
            C.create_vertex();
        }
        std::vector<MeshHalfedges::Halfedge> v2h(M.vertices.nb());

        MeshHalfedges MH(M);
        for(index_t f = 0; f < M.facets.nb(); ++f) {
            for(index_t c = M.facets.corners_begin(f);
                c < M.facets.corners_end(f); ++c
            ) {
                index_t v = M.facet_corners.vertex(c);
                v2h[v] = MeshHalfedges::Halfedge(f, c);
            }
        }

        for(index_t v = 0; v < M.vertices.nb(); ++v) {
            index_t fi[3];
            index_t va[3];
            index_t cur = 0;
            MeshHalfedges::Halfedge H = v2h[v];
            do {
                geo_assert(cur < 3);
                fi[cur] = H.facet;
                index_t ca = M.facets.next_corner_around_facet(
                    H.facet, H.corner
                );
                va[cur] = M.facet_corners.vertex(ca);
                bool ok = MH.move_to_prev_around_vertex(H);
                geo_assert(ok);
                ++cur;
            } while(H != v2h[v]);
            // Note: va[] order is different, because of
            //   Mesh numbering -> Triangulation numering
            // conversion !
            C.create_triangle(
                M.vertices.point_ptr(v), 1.0,
                fi[0], fi[1], fi[2], va[2], va[0], va[1]
            );
            C.triangle_dual(v).sym().add_boundary_facet(fi[0]);
            C.triangle_dual(v).sym().add_boundary_facet(fi[1]);
            C.triangle_dual(v).sym().add_boundary_facet(fi[2]);
            C.triangle_dual(v).sym().set_boundary_vertex(v);
        }
        C.set_symbolic_is_surface(true);
    }

    /**
     * \brief Converts a ConvexCell into a Mesh.
     * \details The constructed Mesh \p M is surfacic and
     *  represents the boundary of the ConvexCell \p C.
     * \param[in] C a const reference to the input ConvexCell
     * \param[in] M a reference to the constructed Mesh
     */
    void convex_cell_to_mesh(const ConvexCell& C, Mesh& M) {
        std::vector<index_t> tri_to_v(C.max_t());
        M.clear();
        M.vertices.set_dimension(3);

        index_t cur_v = 0;
        for(index_t t = 0; t < C.max_t(); ++t) {
            if(C.triangle_is_valid(t)) {
                M.vertices.create_vertex(C.triangle_dual(t).point());
                tri_to_v[t] = cur_v;
                ++cur_v;
            }
        }
        vector<index_t> facet_vertices;
        for(index_t v = 0; v < C.max_v(); v++) {
            facet_vertices.resize(0);
            signed_index_t t = C.vertex_triangle(v);
            if(t != -1) {
                

                ConvexCell::Corner first_c(
                    index_t(t), C.find_triangle_vertex(index_t(t), v)
                );
                ConvexCell::Corner c = first_c;
                do {
                    facet_vertices.push_back(tri_to_v[c.t]);
                    C.move_to_next_around_vertex(c);
                } while(c != first_c);

                index_t f = M.facets.create_polygon(facet_vertices.size());
                for(index_t lv=0; lv<facet_vertices.size(); ++lv) {
                    M.facets.set_vertex(f,lv,facet_vertices[lv]);
                }
                
            }
        }
        M.facets.connect();
    }

    /**
     * \brief Saves a Delaunay triangulation to a file
     * \details Saves the vertices of the Delaunay triangulation \p delaunay
     * to file \p filename;
     * \param[in] delaunay a Delaunay triangulation
     * \param[in] filename path to the file to save
     * \retval true if the file could be successfully saved
     * \retval false otherwise
     */
    bool save_points(const Delaunay* delaunay, const std::string& filename) {
        std::ofstream out(filename.c_str());
        if(!out) {
            Logger::out("I/O")
                << "Could not save points to file " << filename << "..."
                << std::endl;
            return false;
        }

        Logger::out("I/O")
            << "Saving points to file " << filename << "..."
            << std::endl;

        out << delaunay->nb_vertices() << std::endl;
        for(index_t i = 0; i < delaunay->nb_vertices(); ++i) {
            for(coord_index_t c = 0; c < delaunay->dimension(); ++c) {
                out << delaunay->vertex_ptr(i)[c] << ' ';
            }
            out << std::endl;
        }
        return true;
    }

    /**
     * \brief Loads a XYZ point file
     * \details Load points from file \p filename in XYZ format and stores the
     * points coordinates to output vector \p points.
     * \param[out] points output vector of points coordinates
     * \param[in] filename path to the file to load
     * \retval true if the file could be successfully loaded
     * \retval false otherwise
     */
    bool load_points(
        vector<double>& points,
        const std::string& filename
    ) {
        try {
            LineInput in(filename);
            if(!in.OK()) {
                return false;
            }
            index_t nb_points = 0;
            index_t cur = 0;
            while(!in.eof() && in.get_line()) {
                in.get_fields();
                switch(in.nb_fields()) {
                    case 1:
                    {
                        nb_points = in.field_as_uint(0);
                        points.resize(3 * nb_points);
                    } break;
                    case 3:
                    {
                        if(cur >= nb_points) {
                            Logger::err("I/O")
                                << "too many points in .xyz file"
                                << std::endl;
                            return false;
                        }
                        points[3 * cur] = in.field_as_double(0);
                        points[3 * cur + 1] = in.field_as_double(1);
                        points[3 * cur + 2] = in.field_as_double(2);
                        ++cur;
                    } break;
                    default:
                    {
                        Logger::err("I/O")
                            << "invalid number of fields in .xyz file"
                            << std::endl;
                        return false;
                    }
                }
            }
        }
        catch(const std::exception& ex) {
            Logger::err("I/O") << ex.what() << std::endl;
            return false;
        }
        return true;
    }

    /**
     * \brief Creates a random point-set on the surface of a sphere
     * \details Generates \p nb_vertices vertices randomly distributed on the
     * surface of a sphere centered on [0.5, 0.5, 0.5] with radius 1. The
     * generated vertices are stored to output vector \p vertices.
     * \param[out] vertices output vector of vertices coordinates
     * \param[in] nb_vertices number of vertices to generate
     */
    void init_sphere(vector<double>& vertices, index_t nb_vertices) {
        vertices.resize((nb_vertices + 1) * 3);
        vertices[0] = 0.0;
        vertices[1] = 0.0;
        vertices[2] = 0.0;
        for(index_t i = 1; i < nb_vertices + 1; ++i) {
            vec3 v(
                Numeric::random_float64(),
                Numeric::random_float64(),
                Numeric::random_float64()
            );
            v = normalize(v - vec3(0.5, 0.5, 0.5));
            vertices[3 * i] = v[0];
            vertices[3 * i + 1] = v[1];
            vertices[3 * i + 2] = v[2];
        }
    }

    /**
     * \brief Creates a random point-set on the surface of a cone
     * \details Generates \p nb_vertices vertices randomly distributed on the
     * surface of a cone. The generated vertices are stored to output vector
     * \p vertices.
     * \param[out] points output vector of vertices coordinates
     * \param[in] nb_vertices number of vertices to generate
     * \param[in] use_random_vertices if true, random angular sampling is used
     */
     void init_cone(
        vector<double>& vertices, index_t nb_vertices, bool use_random_vertices=true
     ) {
        vertices.resize((nb_vertices + 1) * 3);
        vertices[0] = 0.0;
        vertices[1] = 0.0;
        vertices[2] = 0.0;
        for(index_t i = 1; i < nb_vertices + 1; ++i) {
            vec3 N;
            if(use_random_vertices) {
                N = vec3(
                    Numeric::random_float64(), Numeric::random_float64(), 0.0
                );
                N = 2.0 * N - vec3(1.0, 1.0, 0.0);
            } else {
                double s = ::sin(double(i) * 2.0 * M_PI / double(nb_vertices - 1));
                double c = ::cos(double(i) * 2.0 * M_PI / double(nb_vertices - 1));
                N = vec3(s, c, 0.0);
            }
            N = normalize(N);
            N = N - vec3(0.0, 0.0, 1.0);
            vertices[3 * i] = N[0];
            vertices[3 * i + 1] = N[1];
            vertices[3 * i + 2] = N[2];
        }
    }
}

int main(int argc, char** argv) {

    using namespace GEO;
    using GEOGen::ConvexCell;

    GEO::initialize();
    int result = 0;

    try {

        Stopwatch W("Total time");

        std::vector<std::string> filenames;
        std::string output_filename = "out.eobj";
        std::string points_filename;

        CmdLine::import_arg_group("standard");
        CmdLine::import_arg_group("algo");

        CmdLine::declare_arg("nb_clip", 10, "number of clipping planes");
        CmdLine::declare_arg(
            "nb_clip_times", 1, "number of times clipping is done"
        );
        CmdLine::declare_arg(
            "shape", "sphere", "one of sphere,cone,file"
        );
        CmdLine::declare_arg(
            "save_points", false, "save generated points into points.xyz"
        );

        if(
            !CmdLine::parse(
                argc, argv, filenames, "<pointsfile> <outputfile>"
            )
        ) {
            return 1;
        }

        Delaunay_var delaunay = Delaunay::create(3);
        index_t nb_vertices = CmdLine::get_arg_uint("nb_clip");
        index_t nb_times = CmdLine::get_arg_uint("nb_clip_times");
        std::string shape = CmdLine::get_arg("shape");
        bool exact = (CmdLine::get_arg("algo:predicates") == "exact");

        if(exact) {
            Logger::out("Predicates")
                << "Using exact predicates" << std::endl;
        }

        vector<double> vertices;
        if(shape == "sphere") {
            Logger::out("Program")
                << "Using random sphere shape" << std::endl;
            init_sphere(vertices, nb_vertices);
        } else if(shape == "cone") {
            Logger::out("Program")
                << "Using random cone shape" << std::endl;
            init_cone(vertices, nb_vertices);
        } else if(shape == "file") {
            Logger::out("Program")
                << "Taking shape from file" << std::endl;

            if(filenames.size() == 0) {
                Logger::err("Program")
                    << "Missing input shape file argument" << std::endl;
                return 1;
            }

            points_filename = filenames[0];
            filenames.erase(filenames.begin());

            Logger::out("I/O")
                << "Loading shape from file " << points_filename
                << std::endl;

            if(!load_points(vertices, points_filename)) {
                Logger::err("I/O")
                    << "Could not load file: " << points_filename
                    << std::endl;
                return 1;
            }
            nb_vertices = (vertices.size() / 3) - 1;
        } else {
            Logger::err("Program")
                << shape << ": invalid shape" << std::endl;
            return 1;
        }

        if(filenames.size() >= 1) {
            output_filename = filenames[0];
        }

        if(filenames.size() > 1) {
            Logger::warn("Program")
                << "Extraneous files ignored" << std::endl;
        }

        delaunay->set_vertices(nb_vertices + 1, &(vertices[0]));

        if(CmdLine::get_arg_bool("save_points")) {
            if(!save_points(delaunay, "points.xyz")) {
                return 1;
            }
        }

        CmdLine::set_arg("nb_clip", delaunay->nb_vertices() - 1);

        Mesh M;
        initialize_mesh_with_box(M);
        ConvexCell C(3);
        initialize_convex_cell_with_mesh(C, M);

        for(index_t k = 0; k < nb_times; ++k) {
            for(index_t i = 1; i < nb_vertices; ++i) {
                C.clip_by_plane<3>(&M, delaunay, 0, i, exact, exact);
            }
        }

        Mesh C_mesh;
        convex_cell_to_mesh(C, C_mesh);

        Logger::out("I/O")
            << "Saving mesh to file " << output_filename
            << std::endl;

        if(!mesh_save(C_mesh, output_filename)) {
            return 1;
        }

        if(!meshes_have_same_topology(M, C_mesh, true)) {
            return 1;
        }
    }
    catch(const std::exception& e) {
        std::cerr << "Received an exception: " << e.what() << std::endl;
        return 1;
    }

    return result;
}

