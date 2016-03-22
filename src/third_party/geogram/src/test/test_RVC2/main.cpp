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

#include <geogram/mesh/mesh.h>
#include <geogram/mesh/mesh_io.h>
#include <geogram/mesh/mesh_geometry.h>
#include <geogram/mesh/mesh_repair.h>
#include <geogram/delaunay/delaunay.h>
#include <geogram/delaunay/delaunay_nn.h>
#include <geogram/points/nn_search.h>
#include <geogram/voronoi/RVD.h>
#include <geogram/basic/matrix.h>
#include <geogram/basic/logger.h>
#include <geogram/basic/command_line.h>
#include <geogram/basic/command_line_args.h>
#include <geogram/basic/file_system.h>
#include <geogram/basic/progress.h>
#include <stdarg.h>
#include "least_square_normal.h"


namespace {

    using namespace GEO;

   /**
    * \brief Creates a mesh from a cube.
    * \param[out] M the resulting mesh
    */
    void initialize_mesh_with_box(Mesh& M) {
        M.clear();
        M.vertices.set_dimension(3);
        const double d = 1.0;
        M.vertices.create_vertex(vec3(-d, -d, -d).data());
        M.vertices.create_vertex(vec3(-d, -d, d).data());
        M.vertices.create_vertex(vec3(-d, d, -d).data());
        M.vertices.create_vertex(vec3(-d, d, d).data());
        M.vertices.create_vertex(vec3(d, -d, -d).data());
        M.vertices.create_vertex(vec3(d, -d, d).data());
        M.vertices.create_vertex(vec3(d, d, -d).data());
        M.vertices.create_vertex(vec3(d, d, d).data());
        M.facets.create_quad(7,6,2,3);
        M.facets.create_quad(1,3,2,0);
        M.facets.create_quad(5,7,3,1);
        M.facets.create_quad(4,6,7,5);
        M.facets.create_quad(4,5,1,0);
        M.facets.create_quad(6,4,0,2);
        M.facets.connect();
    }

        // \param[in,out] M is scaled then centered at center
    void center_scale_mesh(Mesh& M, vec3 center, double radius) {
        double xyz_min[3];
        double xyz_max[3];
        get_bbox(M, xyz_min, xyz_max);
        double scale = 2.0*radius / Geom::distance(xyz_min, xyz_max, 3);
        vec3 g = 0.5*(vec3(xyz_min) + vec3(xyz_max));
        for(index_t i=0; i<M.vertices.nb(); ++i) {
            M.vertices.point(i) -= g;
            M.vertices.point(i) *= scale;
            M.vertices.point(i) += center;
        }
    }

        // \param[in] M_ref is used as a reference for calculating the transformation
        // \param[in,out] submeshes are scaled then centered at center
    void center_scale_submeshes(Mesh& M_ref, vec3 center, double radius, vector<Mesh*>& submeshes) {
        double xyz_min[3];
        double xyz_max[3];
        get_bbox(M_ref, xyz_min, xyz_max);
        double scale = 2.0*radius / Geom::distance(xyz_min, xyz_max, 3);
        vec3 g = 0.5*(vec3(xyz_min) + vec3(xyz_max));
        for(index_t sub = 0; sub != submeshes.size(); ++sub){
                Mesh& M = *submeshes[sub];
                        for(index_t i=0; i<M.vertices.nb(); ++i) {
                                M.vertices.point(i) -= g;
                                M.vertices.point(i) *= scale;
                                M.vertices.point(i) += center;
                        }
        }
    }

    vec3 M_x_V(Matrix<double, 3>& M, const vec3& V){
        vec3 ret;
        for(index_t i = 0; i != 3; ++i){
                ret.data()[i] = 0.;
                for(index_t j = 0; j != 3; ++j){
                                ret.data()[i] += M(i,j)*V.data()[j];
                }
        }
        return ret;
    }


        // \param[in,out] M is scaled then rotated then centered at center
    void center_rotate_scale_mesh(Mesh& M, vec3 center, Matrix<double,3>& Rot, double radius){
        double xyz_min[3];
        double xyz_max[3];
        get_bbox(M, xyz_min, xyz_max);
        double scale = 2.0*radius / Geom::distance(xyz_min, xyz_max, 3);
        vec3 g = 0.5*(vec3(xyz_min) + vec3(xyz_max));
        for(index_t i=0; i<M.vertices.nb(); ++i) {
            M.vertices.point(i) -= g;
            M.vertices.point(i) *= scale;
            M.vertices.point(i) = M_x_V(Rot, M.vertices.point(i));
            M.vertices.point(i) += center;
        }
    }

        // \param[in] M_ref is used as a reference for calculating the transformation
        // \param[in,out] submeshes are scaled then then rotated then centered at center
    void center_rotate_scale_submeshes(Mesh& M_ref, vec3 center, 
                Matrix<double,3>& Rot, double radius, vector<Mesh*>& submeshes){
                double xyz_min[3];
                double xyz_max[3];
                get_bbox(M_ref, xyz_min, xyz_max);
                double scale = 2.0*radius / Geom::distance(xyz_min, xyz_max, 3);
                vec3 g = 0.5*(vec3(xyz_min) + vec3(xyz_max));
                for(index_t sub = 0; sub != submeshes.size(); ++sub){
                        Mesh& M = *submeshes[sub];
                        for(index_t i=0; i<M.vertices.nb(); ++i) {
                                M.vertices.point(i) -= g;
                                M.vertices.point(i) *= scale;
                                M.vertices.point(i) = M_x_V(Rot, M.vertices.point(i));
                                M.vertices.point(i) += center;
                        }
                }
    }

        // Computes the rotation matrix from vector A to vector B
        // can be later used to align a shape to a vector
    void rotation_matrix_A_to_B(const vec3& Ain, const vec3& Bin, Matrix<double,3>& Result){
        geo_debug_assert(Ain.length() != 0. && Bin.length() != 0.);
        vec3 A;
        vec3 B;
                if(Ain.length() != 1.){ A = Ain / Ain.length(); } else { A = Ain; }
                if(Bin.length() != 1.){ B = Bin / Bin.length(); } else { B = Bin; }

        // Cross product
        vec3 V = cross(A,B);
        
        // Sine of angle between A and B
        double s = V.length();
        if(std::fabs(s) < 1e-8){ // A and B already aligned
                Result.load_identity();
                if(A.x == -B.x && A.y == -B.y && A.z == -B.z)
                        Result *= -1.;
                return;
        }

        // Skew-symmetric cross-product matrix of a^b
        Matrix<double, 3> ssm;
        ssm.load_zero();
        ssm(0,1) = -V.z; ssm(0,2) = V.y;
        ssm(1,0) = V.z; ssm(1,2) = -V.x;
        ssm(2,0) = -V.y; ssm(2,1) = V.x;
        // Cosine of angle between A and B
        double c = A.x * B.x + A.y * B.y + A.z * B.z;
        Matrix<double, 3> sq_ssm = ssm * ssm;
        sq_ssm *= (1.-c)/(s*s);


                // Building result
                // source : http://math.stackexchange.com/a/476311
        Result.load_identity();
        Result += ssm;
        Result += sq_ssm;
    }

        // Will create an independant surface mesh in vector out for
        // each cell in volumetric mesh in
        // Only works if cells are tetras for the moment
        // Warning : this method use new to create the meshes in out, the user
        // is responsible for the destruction of the output meshes
    void volumetric_mesh_to_surface_meshes(Mesh& in, vector<Mesh*>& out){
        out.clear();
        out.reserve(in.cells.nb());
                for(index_t c = 0; c != in.cells.nb(); ++c){
                        geo_assert(in.cells.type(c) == GEO::MESH_TET);
                        Mesh* m = new Mesh(3);  
                        // Copy the vertices of the cell in the new mesh
                        m->vertices.create_vertices(in.cells.nb_vertices(c));
                        for(index_t lv = 0; lv != in.cells.nb_vertices(c); ++lv){
                                m->vertices.point(lv) = in.vertices.point(in.cells.vertex(c, lv));
                        }

                        // Create facets in the new mesh for each facet of the input cell
                        for(index_t lf = 0; lf != in.cells.nb_facets(c); ++lf){
                                index_t v1 = in.cells.local_tet_facet_vertex_index(lf, 0);
                                index_t v2 = in.cells.local_tet_facet_vertex_index(lf, 1);
                                index_t v3 = in.cells.local_tet_facet_vertex_index(lf, 2);
                                m->facets.create_triangle(v1, v2, v3);
                        }
                        m->facets.connect();
                        out.push_back(m);
                } 
    }

    /**
     * \brief Shrinks a mesh.
     * \param[in,out] M the mesh to be shrunk
     * \param[in] factor the shrinking factor (1.0 means
     *   no shrinking, 0.5 means average shrinking).
     */
    void shrink_mesh(Mesh& M, double factor) {
        double xyz_min[3];
        double xyz_max[3];
        get_bbox(M, xyz_min, xyz_max);
        vec3 g = 0.5*(vec3(xyz_min) + vec3(xyz_max));
        for(index_t i=0; i<M.vertices.nb(); ++i) {
            M.vertices.point(i) -= g;
            M.vertices.point(i) *= factor;
            M.vertices.point(i) += g;
        }
    }

        // This shrink uses the vector center as a reference center
        // for the scaling
    void shrink_mesh(Mesh& M, double factor, vec3 center) {
        for(index_t i=0; i<M.vertices.nb(); ++i) {
            M.vertices.point(i) -= center;
            M.vertices.point(i) *= factor;
            M.vertices.point(i) += center;
        }
    }


        // Estimate the normal of a pointset at point using the nb_neigh
        // nearest neighbors
        vec3 estimate_normal(NearestNeighborSearch& nns, index_t point,
                index_t nb_neigh){
        if(nb_neigh > nns.nb_points()) nb_neigh = nns.nb_points();
        // Neighbors searsh, copy-pasted from 
        // Delaunay_NearestNeighbors::get_neighbors_internal()
                index_t* closest_pt_ix = (index_t*) alloca(sizeof(index_t) * nb_neigh);
                double* closest_pt_dist = (double*) alloca(sizeof(double) * nb_neigh);
                nns.get_nearest_neighbors(
                                nb_neigh, point, closest_pt_ix, closest_pt_dist);
                vector<index_t> neighbors(nb_neigh);
                index_t nb_neigh_result = 0;
                for(index_t j = 0; j < nb_neigh; j++) {
                        geo_debug_assert(signed_index_t(closest_pt_ix[j]) >= 0);
                        if(closest_pt_ix[j] != point) {
                                // Check for duplicated points
                                if(closest_pt_dist[j] == 0.0) {
                                        // If i is not the first one (in the
                                        // duplicated points), then we 'disconnect' it
                                        // (no neighbor !)
                                        geo_debug_assert(signed_index_t(closest_pt_ix[j]) >= 0);
                                        if(closest_pt_ix[j] < point) {
                                                return vec3(0,0,0);
                                        }
                                        // Else, i is the first one, and we simply
                                        // skip (do not store) the connection with
                                        // closest_pt_ix[j].
                                } else {
                                        neighbors[nb_neigh_result] = closest_pt_ix[j];
                                        nb_neigh_result++;
                                        if(nb_neigh_result == nb_neigh - 1) {
                                                break;
                                        }
                                }
                        }
                }
                if(nb_neigh_result != nb_neigh) neighbors.resize(nb_neigh_result);

                // Normal estimation
                LeastSquaresNormal least_squares_normal_;
                least_squares_normal_.begin();
                for (index_t jj = 0; jj < neighbors.size(); jj++) {
                        const double* nei_coord = nns.point_ptr(neighbors[jj]);
                        vec3 nei(nei_coord[0], nei_coord[1], nei_coord[2]);
                        least_squares_normal_.add_point(nei);
                }
                least_squares_normal_.end();
                return least_squares_normal_.get_normal();
    }
}


int main(int argc, char** argv) {

    GEO::initialize();
    GEO::Logger::instance()->set_quiet(false);
    // GEO::CmdLine::import_arg_group("standard");
    GEO::CmdLine::import_arg_group("algo");
    GEO::CmdLine::declare_arg_percent("size", 10.0, "elements size, in bbox diagonal percent");
    GEO::CmdLine::declare_arg("shrink", 0.9, "cells shrink for visualization (shrink=1 does nothing)");
    GEO::CmdLine::declare_arg("border_only", false, "output only RVC facets on the border");
    GEO::CmdLine::declare_arg("interior_only", false, "output only RVC facets not on the border");
    GEO::CmdLine::declare_arg("use_normal", true, "z-axis of shapes used for the intersection are aligned to surface normals");
    GEO::CmdLine::declare_arg("nb_neighbors", 20, "number of neighbor points used to estimate normals of the surface");
    
    std::vector<std::string> filenames;
    if(!GEO::CmdLine::parse(argc, argv, filenames, "points_filename <thickening_filename>")) {
        return 1;
    }
    
    if(filenames.size() != 1 && filenames.size() != 2) {
        return 1;
    }

    GEO::Mesh points;
    GEO::MeshIOFlags flags;
    flags.reset_element(GEO::MESH_FACETS);
    flags.reset_element(GEO::MESH_CELLS);    
    GEO::mesh_load(filenames[0], points, flags);
    GEO::mesh_repair(points);

    double diag = GEO::bbox_diagonal(points);
    double size = GEO::CmdLine::get_arg_percent("size",diag);
    double shrink = GEO::CmdLine::get_arg_double("shrink");
    bool border_only = GEO::CmdLine::get_arg_bool("border_only");
    bool interior_only = GEO::CmdLine::get_arg_bool("interior_only");
    bool align_to_normal = GEO::CmdLine::get_arg_bool("use_normal");
    index_t nb_neighbors = GEO::CmdLine::get_arg_uint("nb_neighbors");
    
    //   Since we compute restricted Voronoi cells one cell at a
    // time, the mesh argument of the restricted Voronoi diagram
    // is not used.
    GEO::Mesh dummy_mesh;
    
    //  Create a Delaunay API that encapsulates a Kd-tree
    GEO::Delaunay_var delaunay = Delaunay::create(3,"NN");
    delaunay->set_vertices(points.vertices.nb(), points.vertices.point_ptr(0));
    
    // The Delaunay_NearestNeighbors* is needed to compute an estimation of
    // the normal, this cast will only work if the delaunay is of type "NN"
    // the Delaunay_NearestNeighbors cannot be created directly because his
    // destructor is private
    Delaunay_NearestNeighbors* dnn = dynamic_cast<Delaunay_NearestNeighbors*>(&*delaunay);
    NearestNeighborSearch* nns = dnn->nn_search();
    
    GEO::RestrictedVoronoiDiagram_var RVD = 
        GEO::RestrictedVoronoiDiagram::create(delaunay, &dummy_mesh);
    
    GEO::Mesh cell_base;
    GEO::Mesh cell;
    GEO::Mesh clipped;
    GEO::Attribute<signed_index_t> facet_id;
    if(border_only || interior_only) {
        facet_id.bind(clipped.facets.attributes(),"id");
    }
    
    if(filenames.size() == 2) {
        mesh_load(filenames[1],cell_base);
    } else {
        initialize_mesh_with_box(cell_base);
    }
    
    index_t progress_divider =
        (points.vertices.nb() > 10000) ? 100 : 1;
    
    GEO::ProgressTask task("RVC",points.vertices.nb()/progress_divider);
    // For each point, create a shape on the point
    // and clip it with the Voronoi cell of the point.

        // If the mesh is a surface with node connectivity of 3, 
        // it can be used directly. Else, it need to be decomposed in
        // submeshes that will be clipped one by one
        // The decomposition will only work if the input mesh is composed
        // of tetras
        vector<Mesh*> submeshes_base;
        vector<Mesh*> submeshes;
        if(cell_base.cells.nb() == 0){
                // TODO : verify the vertices connectivity are 3
                submeshes_base.push_back(&cell_base);
        } else {
                // The following method use new to create the submeshes, 
                // do not forget to delete them at the end
                volumetric_mesh_to_surface_meshes(cell_base, submeshes_base);
        }
        submeshes.resize(submeshes_base.size());
        for(index_t s = 0; s != submeshes.size(); ++s){
                submeshes[s] = new Mesh(3);
        }

        Mesh output(3);
        // The loop over points 
        vec3 z_axis(0.,0.,1.);
        Matrix<double,3> Rotation;
    for(GEO::index_t i=0; i<points.vertices.nb(); ++i) {

        if(!(i%progress_divider)) {
            task.progress(i/progress_divider);
        }

        // Import submeshes aligned to the z-axis
        // use new to create submeshes, do not forget to delete
        cell.copy(cell_base); // used to compute the transformation information
        for(index_t s = 0; s != submeshes_base.size(); ++s){
                submeshes[s]->copy(*submeshes_base[s]);
        }
        
        if(align_to_normal){
                vec3 normal = estimate_normal(*nns, i, nb_neighbors);
                        rotation_matrix_A_to_B(z_axis, normal, Rotation);
                        center_rotate_scale_submeshes(cell, points.vertices.point(i), Rotation, size, submeshes);
        } else {
                        center_scale_submeshes(cell, points.vertices.point(i), size, submeshes);
        }
                for(index_t m = 0; m != submeshes.size(); ++m){
                        RVD->compute_RVC(i,*(submeshes[m]),clipped,facet_id.is_bound());

                        if(shrink != 1.0) {
                                shrink_mesh(clipped, shrink, points.vertices.point(i));
                        }

                        //  Append the generated mesh to the output mesh.
                        index_t new_v = output.vertices.create_vertices(clipped.vertices.nb());
                        for(index_t j=0; j<clipped.vertices.nb(); ++j) {
                                output.vertices.point(new_v+j) = clipped.vertices.point(j);
                        }
                        for(index_t f=0; f<clipped.facets.nb(); ++f) {
                                if(border_only && facet_id[f] >= 0) {
                                        continue;
                                }
                                if(interior_only && facet_id[f] <= 0){
                                        continue;
                                }
                                index_t new_f = output.facets.create_polygon(clipped.facets.nb_vertices(f));
                                index_t lv = 0;
                                for(index_t c=clipped.facets.corners_begin(f); c<clipped.facets.corners_end(f); ++c) {
                                        index_t v = clipped.facet_corners.vertex(c) + new_v;
                                        output.facets.set_vertex(new_f, lv, v);
                                        ++lv;
                                }
                        }
                }
    }
        // Output of the program
        mesh_save(output, "RVC_output.obj");

    // Free memory allocated to submeshes
        for(index_t m = 0; m < submeshes.size(); ++m){
                delete submeshes[m];
        }
        if(submeshes_base.size() >= 2){
                for(index_t m = 0; m < submeshes_base.size(); ++m){
                        delete submeshes_base[m];
                }
        }
    
    return 0;
}
