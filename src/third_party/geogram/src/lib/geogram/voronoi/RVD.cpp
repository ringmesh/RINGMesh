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

#include <geogram/voronoi/RVD.h>
#include <geogram/voronoi/generic_RVD.h>
#include <geogram/voronoi/RVD_mesh_builder.h>
#include <geogram/voronoi/integration_simplex.h>
#include <geogram/mesh/mesh_partition.h>
#include <geogram/mesh/mesh_sampling.h>
#include <geogram/mesh/mesh_repair.h>
#include <geogram/mesh/mesh_AABB.h>
#include <geogram/delaunay/delaunay.h>
#include <geogram/basic/geometry_nd.h>
#include <geogram/basic/process.h>
#include <geogram/basic/command_line.h>
#include <geogram/basic/argused.h>

/*
 * There are three levels of implementation:
 * Level 1: RestrictedVoronoiDiagram is the abstract API seen from client code
 * Level 2: RVD_Nd_Impl<DIM> implements RestrictedVoronoiDiagram
 * Level 3: RVD_Nd_Impl<DIM>::GenRestrictedVoronoiDiagram is
 *  an instanciation of GEOGen::RestrictedVoronoiDiagram (from generic_RVD.h)
 *
 * Warning: there are approx. 1000 lines of boring code ahead.
 */

namespace {

    using namespace GEO;

    /**
     * \brief Generic implementation of RestrictedVoronoiDiagram.
     * \tparam DIM dimension
     */
    template <int DIM>
    class RVD_Nd_Impl : public GEO::RestrictedVoronoiDiagram {

        /** \brief This class type */
        typedef RVD_Nd_Impl<DIM> thisclass;

        /** \brief The base class of this class */
        typedef RestrictedVoronoiDiagram baseclass;

    public:
        /** \brief Implementation based on the generic version. */
        typedef GEOGen::RestrictedVoronoiDiagram<DIM>
            GenRestrictedVoronoiDiagram;

        /** \brief Representation of points. */
        typedef vecng<DIM, double> Point;

        /** \brief Representation of vectors. */
        typedef vecng<DIM, double> Vector;

        /** \brief Represents a point and its symbolic information. */
        typedef typename GenRestrictedVoronoiDiagram::Vertex Vertex;

        /**
         * \brief Specifies the computation done by the threads.
         */
        enum ThreadMode {
            MT_NONE,            /**< uninitialized                          */
            MT_LLOYD,           /**< Lloyd iteration                        */
            MT_NEWTON,          /**< Newton optimization                    */
            MT_INT_SMPLX        /**< Newton with integration simplex        */
        };

        /**
         * \brief Gets a mesh vertex from its index.
         * \param[in] v index of the vertex
         * \return a const reference to a Point
         */
        const Point& mesh_vertex(index_t v) {
            return *(const Point*) mesh_->vertex_ptr(v);
        }

        /**
         * \brief Gets a mesh vertex from a corner index.
         * \param[in] c index of the corner
         * \return a const reference to a Point
         */
        const Point& mesh_corner_vertex(index_t c) {
            index_t v = mesh_->corner_vertex_index(c);
            return mesh_vertex(v);
        }

        /**
         * \brief Gets a Delaunay vertex from its index.
         * \param[in] v index of the Delaunay vertex
         * \return a const reference to a Point
         */
        const Point& seed(index_t v) {
            return *(const Point*) delaunay_->vertex_ptr(v);
        }

        /**
         * \brief Creates a RVD_Nd_Impl.
         *
         * \details The dimension is determined by \p mesh->dimension().
         * \param[in] delaunay the Delaunay triangulation
         * \param[in] mesh the input mesh
         * \param[in] R3_embedding gives for each vertex
         *  its mapping in 3D space.
         * \param[in] R3_embedding_stride gives the stride between
         *  two consecutive vertices in R3_embedding
         */
        RVD_Nd_Impl(
            Delaunay* delaunay, Mesh* mesh,
            const double* R3_embedding, index_t R3_embedding_stride
        ) :
            RestrictedVoronoiDiagram(
                delaunay, mesh, R3_embedding, R3_embedding_stride
            ),
            RVD_(delaunay, mesh) {
            use_exact_projection_ = false;
            is_slave_ = false;
            master_ = nil;
            has_weights_ = false;
            for(index_t v = 0; v < mesh->nb_vertices(); v++) {
                if(mesh->weight(v) != 1.0) {
                    has_weights_ = true;
                    break;
                }
            }
            parts_ = nil;
            nb_parts_ = 0;
            funcval_ = 0.0;
            simplex_func_ = nil;
            arg_vectors_ = nil;
            arg_scalars_ = nil;
            thread_mode_ = MT_NONE;
            nb_triangles_ = 0;
        }

        /**
         * \brief Constructor for parts, used in multithreading mode.
         */
        RVD_Nd_Impl() :
            RestrictedVoronoiDiagram(nil, nil, nil, 0),
            RVD_(nil, nil) {
            use_exact_projection_ = false;
            is_slave_ = true;
            master_ = nil;
            mesh_ = nil;
            parts_ = nil;
            nb_parts_ = 0;
            facets_begin_ = -1;
            facets_end_ = -1;
            funcval_ = 0.0;
            simplex_func_ = nil;
            arg_vectors_ = nil;
            arg_scalars_ = nil;
            thread_mode_ = MT_NONE;
            nb_triangles_ = 0;
        }

        virtual void set_delaunay(Delaunay* delaunay) {
            baseclass::set_delaunay(delaunay);
            RVD_.set_delaunay(delaunay);
            for(index_t p = 0; p < nb_parts_; ++p) {
                parts_[p].set_delaunay(delaunay);
            }
        }

        virtual void set_check_SR(bool x) {
            RVD_.set_check_SR(x);
            for(index_t p = 0; p < nb_parts_; ++p) {
                parts_[p].set_check_SR(x);
            }
        }

        virtual void set_exact_predicates(bool x) {
            RVD_.set_exact_predicates(x);
            for(index_t p = 0; p < nb_parts_; ++p) {
                parts_[p].set_exact_predicates(x);
            }
        }

        virtual bool exact_predicates() const {
            return RVD_.exact_predicates();
        }

        /********************************************************************/

        /**
         * \brief Place holder, "no locking" policy.
         * \details NoLocks is used by algorithms templated
         *  by locking policy, for the single-threaded instances
         *  that do not need synchronization. The multi-threaded
         *  instances are parameterized by SpinLockArray.
         */
        class NoLocks {
        public:
            /**
             * \brief Acquires a spinlock.
             * \details Does nothing in this version
             * \param[in] i index of the spinlock to acquire
             */
            void acquire_spinlock(index_t i) {
                geo_argused(i);
            }

            /**
             * \brief Releases a spinlock.
             * \details Does nothing in this version
             * \param[in] i index of the spinlock to release
             */
            void release_spinlock(index_t i) {
                geo_argused(i);
            }
        };

        // ____________________________________________________________________

        /**
         * \brief Implementation class for surfacic Lloyd relaxation.
         * \details To be used as a template argument
         *    to RVD::for_each_triangle().
         * This version ignores the weights.
         *
         * Computes for each RVD cell:
         * - mg[v] (v's Voronoi cell's total area times centroid)
         * - m[v]  (v's total area)
         * \tparam LOCKS locking policy
         *   (can be one of Process::SpinLockArray, NoLocks)
         */
        template <class LOCKS>
        class ComputeCentroids {
        public:
            /**
             * \brief Constructs a ComputeCentroids.
             * \param[out] mg where to store the centroids
             * \param[out] m where to store the masses
             * \param[in] locks the array of locks 
             *  (or NoLocks in single thread mode)
             */
            ComputeCentroids(
                double* mg,
                double* m,
                LOCKS& locks
            ) :
                mg_(mg),
                m_(m),
                locks_(locks) {
            }

            /**
             * \brief The callback called for each integration simplex.
             * \param[in] v index of current center vertex
             * \param[in] p1 first vertex of current integration simplex
             * \param[in] p2 second vertex of current integration simplex
             * \param[in] p3 third vertex of current integration simplex
             */
            void operator() (
                index_t v,
                const double* p1,
                const double* p2,
                const double* p3
            ) const {
                double cur_m = Geom::triangle_area(p1, p2, p3, DIM);
                double s = cur_m / 3.0;
                locks_.acquire_spinlock(v);
                m_[v] += cur_m;
                double* cur_mg_out = mg_ + v * DIM;
                for(coord_index_t coord = 0; coord < DIM; coord++) {
                    cur_mg_out[coord] += s * (p1[coord] + p2[coord] + p3[coord]);
                }
                locks_.release_spinlock(v);
            }

        private:
            double* mg_;
            double* m_;
            LOCKS& locks_;
        };

        /**
         * \brief Implementation class for surfacic Lloyd relaxation.
         * \details To be used as a template
         *    argument to RVD::for_each_triangle().
         * This version takes the weights into account.
         *
         * Computes for each RVD cell:
         * - mg[v] (v's Voronoi cell's total area times centroid)
         * - m[v]  (v's total area)
         * \tparam LOCKS locking policy
         *   (can be one of Process::SpinLockArray, NoLocks)
         */
        template <class LOCKS>
        class ComputeCentroidsWeighted {
        public:
            /**
             * \brief Constructs a ComputeCentroidsWeighted.
             * \param[out] mg where to store the centroids
             * \param[out] m where to store the masses
             * \param[in] locks the array of locks 
             *  (or NoLocks in single thread mode)
             */
            ComputeCentroidsWeighted(
                double* mg,
                double* m,
                LOCKS& locks
            ) :
                mg_(mg),
                m_(m),
                locks_(locks) {
            }

            /**
             * \brief The callback called for each integration simplex.
             * \param[in] v index of current center vertex
             * \param[in] v1 first vertex of current integration simplex
             * \param[in] v2 second vertex of current integration simplex
             * \param[in] v3 third vertex of current integration simplex
             */
            void operator() (
                index_t v,
                const Vertex& v1,
                const Vertex& v2,
                const Vertex& v3
            ) const {
                double cur_m;
                double cur_Vg[DIM];
                Geom::triangle_centroid(
                    v1.point(), v2.point(), v3.point(),
                    v1.weight(), v2.weight(), v3.weight(),
                    cur_Vg, cur_m, DIM
                );
                locks_.acquire_spinlock(v);
                m_[v] += cur_m;
                double* cur_mg_out = mg_ + v * DIM;
                for(coord_index_t coord = 0; coord < DIM; coord++) {
                    cur_mg_out[coord] += cur_Vg[coord];
                }
                locks_.release_spinlock(v);
            }

        private:
            double* mg_;
            double* m_;
            LOCKS& locks_;
        };

        virtual void compute_centroids_on_surface(double* mg, double* m) {
            create_threads();
            if(nb_parts() == 0) {
                if(master_ != nil) {
                    if(has_weights_) {
                        RVD_.for_each_triangle(
                            ComputeCentroidsWeighted<Process::SpinLockArray>(
                                mg, m, master_->spinlocks_
                            )
                        );
                    } else {
                        RVD_.for_each_triangle(
                            ComputeCentroids<Process::SpinLockArray>(
                                mg, m, master_->spinlocks_
                            )
                        );
                    }
                } else {
                    NoLocks nolocks;
                    if(has_weights_) {
                        RVD_.for_each_triangle(
                            ComputeCentroidsWeighted<NoLocks>(
                                mg, m, nolocks
                            )
                        );
                    } else {
                        RVD_.for_each_triangle(
                            ComputeCentroids<NoLocks>(mg, m, nolocks)
                        );
                    }
                }
            } else {
                thread_mode_ = MT_LLOYD;
                arg_vectors_ = mg;
                arg_scalars_ = m;
                spinlocks_.resize(delaunay_->nb_vertices());
                parallel_for(
                    parallel_for_member_callback(this, &thisclass::run_thread),
                    0, nb_parts()
                );
            }
        }

        /********************************************************************/

        /**
         * \brief Implementation class for surfacic Lloyd relaxation.
         * \details To be used as a template argument
         *    to RVD::for_each_volumetric_integration_simplex().
         * This version ignores the weights.
         *
         * Computes for each RVD cell:
         * - mg[v] (v's Voronoi cell's total area times centroid)
         * - m[v]  (v's total area)
         * \tparam LOCKS locking policy
         *   (can be one of Process::SpinLockArray, NoLocks)
         */
        template <class LOCKS>
        class ComputeCentroidsVolumetric {
        public:
            /**
             * \brief Constructs a ComputeCentroidsVolumetric.
             * \param[out] mg where to store the centroids
             * \param[out] m where to store the masses
             * \param[in] delaunay the Delaunay triangulation
             * \param[in] locks the array of locks
             *  (or NoLocks in single thread mode)
             */
            ComputeCentroidsVolumetric(
                double* mg,
                double* m,
                const Delaunay* delaunay,
                LOCKS& locks
            ) :
                mg_(mg),
                m_(m),
                delaunay_(delaunay),
                locks_(locks) {
            }

            /**
             * \brief The callback called for each integration simplex.
             * \param[in] v index of current center vertex
             * \param[in] v_adj (unused here) is the index of the Voronoi cell
             *  adjacent to t accross facet (\p v1, \p v2, \p v3) or
             *  -1 if it does not exists
             *  \param[in] t (unused here) is the index of the current
             *   tetrahedron
             *  \param[in] t_adj (unused here) is the index of the
             *   tetrahedron adjacent to t accross facet (\p v1, \p v2, \p v3)
             *   or -1 if it does not exists
             * \param[in] p0 first vertex of current integration simplex
             * \param[in] p1 second vertex of current integration simplex
             * \param[in] p2 third vertex of current integration simplex
             * \param[in] p3 fourth vertex of current integration simplex
             */
            void operator() (
                index_t v, signed_index_t v_adj,
                index_t t, signed_index_t t_adj,
                const double* p0,
                const double* p1,
                const double* p2,
                const double* p3
            ) const {
                geo_argused(v_adj);
                geo_argused(t);
                geo_argused(t_adj);
                double cur_m = Geom::tetra_volume<DIM>(
                    p0, p1, p2, p3
                );
                double s = cur_m / 4.0;
                locks_.acquire_spinlock(v);
                m_[v] += cur_m;
                double* cur_mg_out = mg_ + v * DIM;
                for(coord_index_t coord = 0; coord < DIM; coord++) {
                    cur_mg_out[coord] += s * (
                        p0[coord] + p1[coord] + p2[coord] + p3[coord]
                    );
                }
                locks_.release_spinlock(v);
            }

        private:
            double* mg_;
            double* m_;
            const Delaunay* delaunay_;
            LOCKS& locks_;
        };

        virtual void compute_centroids_in_volume(double* mg, double* m) {
            create_threads();
            if(nb_parts() == 0) {
                if(master_ != nil) {
                    RVD_.for_each_tetrahedron(
                        ComputeCentroidsVolumetric<Process::SpinLockArray>(
                            mg, m, RVD_.delaunay(), master_->spinlocks_
                        )
                    );
                } else {
                    NoLocks nolocks;
                    RVD_.for_each_tetrahedron(
                        ComputeCentroidsVolumetric<NoLocks>(
                            mg, m, RVD_.delaunay(), nolocks
                        )
                    );
                }
            } else {
                thread_mode_ = MT_LLOYD;
                arg_vectors_ = mg;
                arg_scalars_ = m;
                spinlocks_.resize(delaunay_->nb_vertices());
                parallel_for(
                    parallel_for_member_callback(this, &thisclass::run_thread),
                    0, nb_parts()
                );
            }
        }

        /********************************************************************/

        /**
         * \brief Implementation class for Newton-based restricted CVT.
         * \details To be used as a template argument
         *    to RVD::for_each_triangle().
         * This version ignores the weights.
         *
         * Computes for each RVD cell:
         * - g (gradient)
         * - f (CVT energy)
         * \tparam LOCKS locking policy
         *    (can be one of Process::SpinLockArray, NoLocks)
         */
        template <class LOCKS>
        class ComputeCVTFuncGrad {
        public:
            /**
             * \brief Constructs a ComputeCVTFuncGrad.
             * \param[in] RVD the restricted Voronoi diagram
             * \param[out] f the computed function value
             * \param[out] g the computed gradient of f, 
             *   allocated by caller, and managed
             *  by caller
             * \param[in] locks the array of locks 
             *  (or NoLocks in single thread mode)
             */
            ComputeCVTFuncGrad(
                const GenRestrictedVoronoiDiagram& RVD,
                double& f,
                double* g,
                LOCKS& locks
            ) :
                f_(f),
                g_(g),
                locks_(locks),
                RVD_(RVD) {
            }

            /**
             * \brief The callback called for each integration simplex.
             * \param[in] v index of current center vertex
             * \param[in] p1 first vertex of current integration simplex
             * \param[in] p2 second vertex of current integration simplex
             * \param[in] p3 third vertex of current integration simplex
             */
            void operator() (
                index_t v,
                const double* p1,
                const double* p2,
                const double* p3
            ) const {

                const double* p0 = RVD_.delaunay()->vertex_ptr(v);

                double t_area = Geom::triangle_area(p1, p2, p3, DIM);

                double cur_f = 0.0;
                for(index_t c = 0; c < DIM; c++) {
                    double u0 = p0[c] - p1[c];
                    double u1 = p0[c] - p2[c];
                    double u2 = p0[c] - p3[c];
                    cur_f += u0 * u0;
                    cur_f += u1 * (u0 + u1);
                    cur_f += u2 * (u0 + u1 + u2);
                }

                f_ += t_area * cur_f / 6.0;

                locks_.acquire_spinlock(v);
                for(index_t c = 0; c < DIM; c++) {
                    double Gc = (1.0 / 3.0) * (p1[c] + p2[c] + p3[c]);
                    g_[DIM * v + c] += (2.0 * t_area) * (p0[c] - Gc);
                }
                locks_.release_spinlock(v);
            }

            double& f_;
            double* g_;
            LOCKS& locks_;
            const GenRestrictedVoronoiDiagram& RVD_;
        };

        /**
         * \brief Implementation class for Newton-based restricted CVT.
         * \details To be used as a template argument
         *    to RVD::for_each_triangle().
         * This version takes the weights into account.
         *
         * Computes for each RVD cell:
         * - g (gradient)
         * - f (CVT energy)
         * \tparam LOCKS locking policy
         *   (can be one of Process::SpinLockArray, NoLocks)
         */
        template <class LOCKS>
        class ComputeCVTFuncGradWeighted {
        public:
            /**
             * \brief Constructs a ComputeCVTFuncGradWeighted.
             * \param[in] RVD the restricted Voronoi diagram
             * \param[out] f the computed function value
             * \param[out] g the computed gradient of f, 
             *  allocated by caller, and managed by caller
             * \param[in] locks the array of locks 
             *  (or NoLocks in single thread mode)
             */
            ComputeCVTFuncGradWeighted(
                const GenRestrictedVoronoiDiagram& RVD,
                double& f,
                double* g,
                LOCKS& locks
            ) :
                f_(f),
                g_(g),
                locks_(locks),
                RVD_(RVD) {
            }

            /**
             * \brief The callback called for each integration simplex.
             * \param[in] v index of current center vertex
             * \param[in] v1 first vertex of current integration simplex
             * \param[in] v2 second vertex of current integration simplex
             * \param[in] v3 third vertex of current integration simplex
             */
            void operator() (
                index_t v,
                const Vertex& v1,
                const Vertex& v2,
                const Vertex& v3
            ) const {

                const double* p0 = RVD_.delaunay()->vertex_ptr(v);

                const double* p1 = v1.point();
                const double* p2 = v2.point();
                const double* p3 = v3.point();

                double t_area = Geom::triangle_area(p1, p2, p3, DIM);

                double Sp = v1.weight() + v2.weight() + v3.weight();
                double rho[3], alpha[3];
                rho[0] = v1.weight();
                rho[1] = v2.weight();
                rho[2] = v3.weight();
                alpha[0] = Sp + rho[0];
                alpha[1] = Sp + rho[1];
                alpha[2] = Sp + rho[2];

                double dotprod_00 = 0.0;
                double dotprod_10 = 0.0;
                double dotprod_11 = 0.0;
                double dotprod_20 = 0.0;
                double dotprod_21 = 0.0;
                double dotprod_22 = 0.0;
                for(unsigned int c = 0; c < DIM; c++) {
                    double sp0 = p0[c] - p1[c];
                    double sp1 = p0[c] - p2[c];
                    double sp2 = p0[c] - p3[c];
                    dotprod_00 += sp0 * sp0;
                    dotprod_10 += sp1 * sp0;
                    dotprod_11 += sp1 * sp1;
                    dotprod_20 += sp2 * sp0;
                    dotprod_21 += sp2 * sp1;
                    dotprod_22 += sp2 * sp2;
                }

                double cur_f = 0.0;
                cur_f += (alpha[0] + rho[0]) * dotprod_00;  // 0 0
                cur_f += (alpha[1] + rho[0]) * dotprod_10;  // 1 0
                cur_f += (alpha[1] + rho[1]) * dotprod_11;  // 1 1
                cur_f += (alpha[2] + rho[0]) * dotprod_20;  // 2 0
                cur_f += (alpha[2] + rho[1]) * dotprod_21;  // 2 1
                cur_f += (alpha[2] + rho[2]) * dotprod_22;  // 2 2

                f_ += t_area * cur_f / 30.0;
                double* g_out = g_ + v * DIM;
                locks_.acquire_spinlock(v);
                for(index_t c = 0; c < DIM; c++) {
                    g_out[c] += (t_area / 6.0) * (
                        4.0 * Sp * p0[c] - (
                            alpha[0] * p1[c] +
                            alpha[1] * p2[c] +
                            alpha[2] * p3[c]
                        )
                    );
                }
                locks_.release_spinlock(v);
            }

            double& f_;
            double* g_;
            LOCKS& locks_;
            const GenRestrictedVoronoiDiagram& RVD_;
        };

        virtual void compute_CVT_func_grad_on_surface(double& f, double* g) {
            create_threads();
            if(nb_parts() == 0) {
                if(master_ != nil) {
                    if(has_weights_) {
                        RVD_.for_each_triangle(
                            ComputeCVTFuncGradWeighted<Process::SpinLockArray>(
                                RVD_, f, g, master_->spinlocks_
                            )
                        );
                    } else {
                        RVD_.for_each_triangle(
                            ComputeCVTFuncGrad<Process::SpinLockArray>(
                                RVD_, f, g, master_->spinlocks_
                            )
                        );
                    }
                } else {
                    NoLocks nolocks;
                    if(has_weights_) {
                        RVD_.for_each_triangle(
                            ComputeCVTFuncGradWeighted<NoLocks>(
                                RVD_, f, g, nolocks
                            )
                        );
                    } else {
                        RVD_.for_each_triangle(
                            ComputeCVTFuncGrad<NoLocks>(
                                RVD_, f, g, nolocks
                            )
                        );
                    }
                }
            } else {
                thread_mode_ = MT_NEWTON;
                arg_vectors_ = g;
                spinlocks_.resize(delaunay_->nb_vertices());
                for(index_t t = 0; t < nb_parts(); t++) {
                    part(t).funcval_ = 0.0;
                }
                parallel_for(
                    parallel_for_member_callback(this, &thisclass::run_thread),
                    0, nb_parts()
                );
                for(index_t t = 0; t < nb_parts(); t++) {
                    f += part(t).funcval_;
                }
            }
        }

        /********************************************************************/

        /**
         * \brief Implementation class for Newton-based restricted CVT
         *  in volume.
         * \details To be used as a template argument
         *    to RVD::for_each_volumetric_integration_simplex().
         * This version ignores the weights.
         *
         * Computes for each RVD cell:
         * - g (gradient)
         * - f (CVT energy)
         * \tparam LOCKS locking policy
         *    (can be one of Process::SpinLockArray, NoLocks)
         */
        template <class LOCKS>
        class ComputeCVTFuncGradVolumetric {
        public:
            /**
             * \brief Constructs a ComputeCentroidsFuncGradVolumetric.
             * \param[in] RVD the restricted Voronoi diagram
             * \param[out] f the computed function value
             * \param[out] g the computed gradient of f,
             *  allocated by caller, and managed
             *  by caller
             * \param[in] locks the array of locks
             *  (or NoLocks in single thread mode)
             */
            ComputeCVTFuncGradVolumetric(
                const GenRestrictedVoronoiDiagram& RVD,
                double& f,
                double* g,
                LOCKS& locks
            ) :
                f_(f),
                g_(g),
                locks_(locks),
                RVD_(RVD) {
            }

            /**
             * \brief The callback called for each integration simplex.
             * \param[in] v index of current center vertex
             * \param[in] v_adj (unused here) is the index of the Voronoi cell
             *  adjacent to t accross facet (\p v1, \p v2, \p v3) or
             *  -1 if it does not exists
             *  \param[in] t (unused here) is the index of the current
             *   tetrahedron
             *  \param[in] t_adj (unused here) is the index of the
             *   tetrahedron adjacent to t accross facet (\p v1, \p v2, \p v3)
             *   or -1 if it does not exists
             * \param[in] p1 first vertex of current integration simplex
             * \param[in] p2 second vertex of current integration simplex
             * \param[in] p3 third vertex of current integration simplex
             */
            void operator() (
                index_t v,
                signed_index_t v_adj,
                index_t t,
                signed_index_t t_adj,
                const double* p1,
                const double* p2,
                const double* p3
            ) const {
                geo_argused(v_adj);
                geo_argused(t);
                geo_argused(t_adj);
                const double* p0 = RVD_.delaunay()->vertex_ptr(v);

                double mi = Geom::tetra_volume<DIM>(p0, p1, p2, p3);

                // fi = (mi/10)*(U.U + V.V + W.W + U.V + V.W + W.U)
                // where: U = p1-p0 ; V = p2-p0 and W=p3-p0
                double fi = 0.0;
                for(coord_index_t c = 0; c < DIM; ++c) {
                    double Uc = p1[c] - p0[c];
                    double Vc = p2[c] - p0[c];
                    double Wc = p3[c] - p0[c];
                    fi += geo_sqr(Uc) + geo_sqr(Vc) + geo_sqr(Wc);
                    fi += (Uc * Vc + Vc * Wc + Wc * Uc);
                }
                fi *= (mi / 10.0);
                f_ += fi;

                // gi = 2*mi(p0 - 1/4(p0 + p1 + p2 + p3))
                double* g_out = g_ + v * DIM;
                locks_.acquire_spinlock(v);
                for(coord_index_t c = 0; c < DIM; ++c) {
                    g_out[c] += 2.0 * mi * (
                        0.75 * p0[c] 
                        - 0.25 * p1[c] - 0.25 * p2[c] - 0.25 * p3[c]
                    );
                }
                locks_.release_spinlock(v);
            }

            double& f_;
            double* g_;
            LOCKS& locks_;
            const GenRestrictedVoronoiDiagram& RVD_;
        };

        virtual void compute_CVT_func_grad_in_volume(double& f, double* g) {
            create_threads();
            if(nb_parts() == 0) {
                if(master_ != nil) {
                    RVD_.for_each_volumetric_integration_simplex(
                        ComputeCVTFuncGradVolumetric<Process::SpinLockArray>(
                            RVD_, f, g, master_->spinlocks_
                        )
                    );
                } else {
                    NoLocks nolocks;
                    RVD_.for_each_volumetric_integration_simplex(
                        ComputeCVTFuncGradVolumetric<NoLocks>(
                            RVD_, f, g, nolocks
                        )
                    );
                }
            } else {
                thread_mode_ = MT_NEWTON;
                arg_vectors_ = g;
                spinlocks_.resize(delaunay_->nb_vertices());
                for(index_t t = 0; t < nb_parts(); t++) {
                    part(t).funcval_ = 0.0;
                }
                parallel_for(
                    parallel_for_member_callback(this, &thisclass::run_thread),
                    0, nb_parts()
                );
                for(index_t t = 0; t < nb_parts(); t++) {
                    f += part(t).funcval_;
                }
            }
        }

        /********************************************************************/

        /**
         * \brief Implementation class for computing function integrals
         *  and gradients over integration simplices.
         * \details To be used as a template argument
         *    to RVD::for_each_triangle() and 
         *    RVD::for_each_volumetric_integration_simplex()
         */
        class ComputeCVTFuncGradIntegrationSimplex {
        public:
            /**
             * \brief Constructs a ComputeCVTFuncGradIntegrationSimplex.
             * \param[in] RVD the Restricted Voronoi Diagram
             * \param[i] F the IntegrationSimplex
             */
            ComputeCVTFuncGradIntegrationSimplex(
                const GenRestrictedVoronoiDiagram& RVD,
                IntegrationSimplex* F
            ) : 
                f_(0.0),
                RVD_(RVD), 
                simplex_func_(F) {
            }

            /**
             * \brief The callback called for each surfacic integration simplex.
             * \param[in] i index of current center vertex
             * \param[in] v1 first vertex of current integration simplex
             * \param[in] v2 second vertex of current integration simplex
             * \param[in] v3 third vertex of current integration simplex
             */
            void operator() (
                index_t i,
                const Vertex& v1,
                const Vertex& v2,
                const Vertex& v3
            ) {
                f_ += simplex_func_->eval(
                    i,v1,v2,v3,RVD_.current_facet()
                );
            }

            /**
             * \brief The callback called for each volumetric
             *   integration simplex.
             * \param[in] v index of current center vertex
             * \param[in] v_adj index of the Voronoi cell adjacent to t accross
             *    facet (\p v1, \p v2, \p v3) or -1 if it does not exists
             * \param[in] t index of the current tetrahedron
             * \param[in] t_adj index of the tetrahedron adjacent to t accross
             *    facet (\p v1, \p v2, \p v3) or -1 if it does not exists
             * \param[in] v1 first vertex of current integration simplex
             * \param[in] v2 second vertex of current integration simplex
             * \param[in] v3 third vertex of current integration simplex
             */
            void operator() (
                index_t v,
                signed_index_t v_adj,
                index_t t,
                signed_index_t t_adj,
                const Vertex& v1,
                const Vertex& v2,
                const Vertex& v3
            ) {
                geo_argused(v_adj);
                geo_argused(t_adj);
                f_ += simplex_func_->eval(
                    v,v1,v2,v3,t
                );
            }

            /**
             * \brief Gets the function value.
             * \return The function value accumulated so far.
             */
            double f() const {
                return f_;
            }

        private:
            double f_;
            const GenRestrictedVoronoiDiagram& RVD_;
            IntegrationSimplex* simplex_func_;
        };

        virtual void compute_integration_simplex_func_grad(
            double& f, double* g, IntegrationSimplex* F
        ) {
            create_threads();
            if(nb_parts() == 0) {
                if(master_ == nil) {
                    F->set_points_and_gradient(
                        delaunay()->dimension(),
                        delaunay()->nb_vertices(), 
                        delaunay()->vertex_ptr(0),
                        g
                    );
                }
                ComputeCVTFuncGradIntegrationSimplex C(RVD_,F);
                bool sym = RVD_.symbolic();
                RVD_.set_symbolic(true);
                if(F->volumetric()) {
                    RVD_.for_each_volumetric_integration_simplex(
                        C,
                        false,  /* Visit inner tets   */
                        false   /* Coherent triangles */
                    );
                } else {
                    RVD_.for_each_triangle(C);
                }
                RVD_.set_symbolic(sym);
                funcval_ = C.f();
                f = C.f();
            } else {
                spinlocks_.resize(delaunay_->nb_vertices());
                F->set_points_and_gradient(
                    delaunay()->dimension(),
                    delaunay()->nb_vertices(), 
                    delaunay()->vertex_ptr(0),
                    g, 
                    &spinlocks_
                );
                thread_mode_ = MT_INT_SMPLX;
                arg_vectors_ = g;
                simplex_func_ = F;
                funcval_ = 0.0;
                for(index_t t = 0; t < nb_parts(); t++) {
                    part(t).arg_vectors_ = g;
                    part(t).simplex_func_ = F;
                    part(t).funcval_ = 0.0;
                }

                parallel_for(
                    parallel_for_member_callback(this, &thisclass::run_thread),
                    0, nb_parts()
                );

                f = 0.0;
                for(index_t t = 0; t < nb_parts(); t++) {
                    f += part(t).funcval_;
                }
            }
        } 

        /********************************************************************/

        /**
         * \brief Implementation class for explicitly constructing
         *    a surfacic mesh that corresponds to the surfacic
         *    restricted Voronoi diagram.
         * \details To be used as a template argument
         *    to RVD::for_each_polygon(). The current Vornoi cell is
         *    reported in facet region.
         * \tparam BUILDER a class that implements iterative mesh building,
         *    e.g., MeshBuilder.
         */
        template <class BUILDER>
        class BuildRVD {
        public:
            /**
             * \brief Constructs a new BuildRVD.
             * \param[in] RVD_in the restricted Voronoi diagram
             * \param[in] builder the lesh builder
             */
            BuildRVD(
                const GenRestrictedVoronoiDiagram& RVD_in,
                BUILDER& builder
            ) :
                RVD(RVD_in),
                builder_(builder),
                current_facet_(-1) {
                builder_.begin_surface();
            }

            /**
             * \brief The destructor
             * \details Terminates the current facet
             *    and the current surface.
             */
            ~BuildRVD() {
                if(current_facet_ != -1) {
                    builder_.end_reference_facet();
                }
                builder_.end_surface();
            }

            /**
             * \brief The callback called for each restricted Voronoi cell.
             * \param[in] v index of current center vertex
             * \param[in] P current restricted Voronoi cell
             */
            void operator() (
                index_t v,
                const typename GenRestrictedVoronoiDiagram::Polygon& P
            ) {
                index_t f = RVD.current_facet();
                if(signed_index_t(f) != current_facet_) {
                    if(current_facet_ != -1) {
                        builder_.end_reference_facet();
                    }
                    current_facet_ = signed_index_t(f);
                    builder_.begin_reference_facet(f);
                }
                builder_.begin_facet(v);
                for(index_t i = 0; i < P.nb_vertices(); i++) {
                    const Vertex& ve = P.vertex(i);
                    builder_.add_vertex_to_facet(ve.point(), ve.sym());
                }
                builder_.end_facet();
            }

        private:
            const GenRestrictedVoronoiDiagram& RVD;
            BUILDER& builder_;
            signed_index_t current_facet_;
        };

        /**
         * \brief Implementation class for explicitly constructing
         *    a volumetric mesh that corresponds to the volumetric
         *    restricted Voronoi diagram.
         * \details To be used as a template argument
         *    to RVD::for_each_volumetric_integration_simplex(). The current Voronoi cell
         *    is reported in tetrahedron region.
         * \note For the moment, vertices are duplicated (will be fixed
         *  in a future version).
         */
        class BuildVolumetricRVD {
        public:
            /**
             * Constructs a new BuildVolumetricRVD.
             * \param[in] RVD the volumetric restricted Voronoi diagram
             * \param[in] dim dimension of the points (can be smaller
             *  than actual dimension of the RVD).
             * \param[out] vertices coordinates of the generated vertices
             * \param[out] triangle_vertex_indices generated triangles, as
             *  vertex indices triplets
             * \param[out] tet_vertex_indices generated tetrahedra, as
             *  vertex indices 4-uples
             * \param[out] triangle_regions each generated triangle has
             *  a region index, that corresponds to the index of the
             *  Voronoi cell the triangle belongs to
             * \param[out] tet_regions each generated tetrahedron has
             *  a region index, that corresponds to the index of the
             *  Voronoi cell the tetrahedron belongs to
             * \param[in] cell_borders_only if true, only the surfacic
             *  borders of the volumetric cells are saved in the mesh, else
             *  volumetric cells are tetrahedralized.
             * \pre dim <= delaunay->dimension()
             */
            BuildVolumetricRVD(
                GenRestrictedVoronoiDiagram& RVD,
                coord_index_t dim,
                vector<double>& vertices,
                vector<index_t>& triangle_vertex_indices,
                vector<index_t>& tet_vertex_indices,
                vector<signed_index_t>& triangle_regions,
                vector<signed_index_t>& tet_regions,
                bool cell_borders_only
            ) :
                delaunay_(RVD.delaunay()),
                mesh_(RVD.mesh()),
                dim_(dim),
                vertices_(vertices),
                triangle_vertex_indices_(triangle_vertex_indices),
                tet_vertex_indices_(tet_vertex_indices),
                triangle_regions_(triangle_regions),
                tet_regions_(tet_regions),
                cell_borders_only_(cell_borders_only)
            {
                vertices_.clear();
                triangle_vertex_indices_.clear();
                tet_vertex_indices_.clear();
                triangle_regions_.clear();
                tet_regions_.clear();

                // The first vertices are copied from Delaunay,
                // the other ones will be created during the traversal
                nb_vertices_ = delaunay_->nb_vertices();
                vertices_.resize(nb_vertices_ * dim);
                for(index_t v = 0; v < delaunay_->nb_vertices(); ++v) {
                    for(coord_index_t c = 0; c < dim; ++c) {
                        vertices_[v * dim + c] = delaunay_->vertex_ptr(v)[c];
                    }
                }
                vertex_map_.set_first_vertex_index(nb_vertices_);
            }

            /**
             * \brief The callback called for each integration simplex.
             * \param[in] v index of current center vertex
             * \param[in] v_adj (unused here) is the index of the Voronoi cell
             *  adjacent to t accross facet (\p v1, \p v2, \p v3) or
             *  -1 if it does not exists
             *  \param[in] t (unused here) is the index of the current
             *   tetrahedron
             *  \param[in] t_adj (unused here) is the index of the
             *   tetrahedron adjacent to t accross facet (\p v1, \p v2, \p v3)
             *   or -1 if it does not exists
             * \param[in] v1 first vertex of current integration simplex
             * \param[in] v2 second vertex of current integration simplex
             * \param[in] v3 third vertex of current integration simplex
             */
            void operator() (
                index_t v, signed_index_t v_adj,
                index_t t, signed_index_t t_adj,
                const Vertex& v1, const Vertex& v2, const Vertex& v3
            ) {
                geo_argused(v_adj);
                geo_argused(t);

                if(cell_borders_only_) {
                    if(signed_index_t(v) > v_adj) {
                        index_t iv2 = find_or_create_vertex(v, v1);
                        index_t iv3 = find_or_create_vertex(v, v2);
                        index_t iv4 = find_or_create_vertex(v, v3);
                        triangle_vertex_indices_.push_back(iv4);
                        triangle_vertex_indices_.push_back(iv3);
                        triangle_vertex_indices_.push_back(iv2);
                        triangle_regions_.push_back(signed_index_t(v));
                    }
                } else {
                    index_t iv1 = v;
                    index_t iv2 = find_or_create_vertex(v, v1);
                    index_t iv3 = find_or_create_vertex(v, v2);
                    index_t iv4 = find_or_create_vertex(v, v3);

                    // Triangle v1,v2,v3 is on border if there is
                    // no adjacent seed and no adjacent tet.
                    if(v_adj == -1 && t_adj == -1) {
                        triangle_vertex_indices_.push_back(iv4);
                        triangle_vertex_indices_.push_back(iv3);
                        triangle_vertex_indices_.push_back(iv2);
                        triangle_regions_.push_back(signed_index_t(v));
                    }

                    tet_vertex_indices_.push_back(iv1);
                    tet_vertex_indices_.push_back(iv2);
                    tet_vertex_indices_.push_back(iv3);
                    tet_vertex_indices_.push_back(iv4);
                    tet_regions_.push_back(signed_index_t(v));
                }
            }

            /**
             * \brief The callback called for each tetrahedron
             * \param[in] v index of current center vertex
             * \param[in] v_adj (unused here) is the index of the Voronoi cell
             *  adjacent to t accross facet (\p v1, \p v2, \p v3) or
             *  -1 if it does not exists
             *  \param[in] t (unused here) is the index of the current
             *   tetrahedron
             *  \param[in] t_adj (unused here) is the index of the
             *   tetrahedron adjacent to t accross facet (\p v1, \p v2, \p v3)
             *   or -1 if it does not exists
             * \param[in] v1 first vertex of current tetrahedron
             * \param[in] v2 second vertex of current tetrahedron
             * \param[in] v3 third vertex of current tetrahedron
             * \param[in] v4 fourth vertex of current tetrahedron
             */
            void operator() (
                index_t v, signed_index_t v_adj,
                index_t t, signed_index_t t_adj,
                const Vertex& v1, const Vertex& v2,
                const Vertex& v3, const Vertex& v4
            ) {
                geo_argused(v_adj);
                geo_argused(t);
                geo_argused(t_adj);
                index_t iv1 = vertices_.size() / dim_;
                for(index_t c = 0; c < dim_; ++c) {
                    vertices_.push_back(v1.point()[c]);
                }
                index_t iv2 = vertices_.size() / dim_;
                for(index_t c = 0; c < dim_; ++c) {
                    vertices_.push_back(v2.point()[c]);
                }
                index_t iv3 = vertices_.size() / dim_;
                for(index_t c = 0; c < dim_; ++c) {
                    vertices_.push_back(v3.point()[c]);
                }
                index_t iv4 = vertices_.size() / dim_;
                for(index_t c = 0; c < dim_; ++c) {
                    vertices_.push_back(v4.point()[c]);
                }
                tet_vertex_indices_.push_back(iv1);
                tet_vertex_indices_.push_back(iv2);
                tet_vertex_indices_.push_back(iv3);
                tet_vertex_indices_.push_back(iv4);
                tet_regions_.push_back(signed_index_t(v));
            }

        protected:
            /**
             * \brief Retreives the index of a vertex given its symbolic
             * representation.
             * \param[in] center_vertex_id index of current Voronoi seed
             * \param[in] v symbolic and geometric representation of the vertex
             * \return the index of the vertex
             */
            index_t find_or_create_vertex(
                index_t center_vertex_id, const Vertex& v
            ) {
                index_t result = vertex_map_.find_or_create_vertex(
                    center_vertex_id, v.sym()
                );
                if(result >= nb_vertices_) {
                    geo_assert(result == nb_vertices_);
                    nb_vertices_ = result + 1;
                    for(coord_index_t c = 0; c < dim_; ++c) {
                        vertices_.push_back(v.point()[c]);
                    }
                }
                return result;
            }

        private:
            const Delaunay* delaunay_;
            const Mesh* mesh_;
            coord_index_t dim_;
            vector<double>& vertices_;
            vector<index_t>& triangle_vertex_indices_;
            vector<index_t>& tet_vertex_indices_;
            vector<signed_index_t>& triangle_regions_;
            vector<signed_index_t>& tet_regions_;
            RVDVertexMap vertex_map_;
            index_t nb_vertices_;
            bool cell_borders_only_;
        };

        virtual void compute_RVD(
            Mesh& M, coord_index_t dim, bool cell_borders_only,
            bool integration_simplices
        ) {
            bool sym = RVD_.symbolic();
            RVD_.set_symbolic(true);
            if(volumetric_) {
                if(dim == 0) {
                    dim = dimension();
                }
                vector<double> vertices;
                vector<index_t> triangle_vertices;
                vector<index_t> tet_vertices;
                vector<signed_index_t> triangle_regions;
                vector<signed_index_t> tet_regions;
                if(cell_borders_only) {
                    RVD_.for_each_volumetric_integration_simplex(
                        BuildVolumetricRVD(
                            RVD_, dim,
                            vertices,
                            triangle_vertices,
                            tet_vertices,
                            triangle_regions,
                            tet_regions,
                            cell_borders_only
                        ),
                        false, // Do not visit inner tetrahedra.
                        true   // Ensure that polygonal facets are triangulated
                               // coherently.
                    );
                } else {
                    if(integration_simplices) {
                        RVD_.for_each_volumetric_integration_simplex(
                            BuildVolumetricRVD(
                                RVD_, dim,
                                vertices,
                                triangle_vertices,
                                tet_vertices,
                                triangle_regions,
                                tet_regions,
                                cell_borders_only
                            ),
                            false, // Do not visit inner tetrahedra.
                            true // Ensure that polygonal facets are triangulated
                                 // coherently.
                        );
                    } else {
                        RVD_.for_each_tetrahedron(
                            BuildVolumetricRVD(
                                RVD_, dim,
                                vertices,
                                triangle_vertices,
                                tet_vertices,
                                triangle_regions,
                                tet_regions,
                                cell_borders_only
                            )
                        );
                    }
                }
                M.clear();
                GEO::MeshMutator::set_dimension(M, dim);
                GEO::MeshMutator::vertices(M).swap(vertices);
                GEO::MeshMutator::set_nb_facets(M, triangle_vertices.size() / 3);
                GEO::MeshMutator::corner_adjacent_facets(M).assign(
                    triangle_vertices.size(), -1
                );
                GEO::MeshMutator::corner_vertices(M).swap(triangle_vertices);
                GEO::MeshMutator::set_triangulated(M, true);
                GEO::MeshMutator::set_nb_cells(M, tet_vertices.size() / 4);
                GEO::MeshMutator::tet_vertices(M).swap(tet_vertices);
                GEO::MeshMutator::tet_adjacents(M).assign(
                    tet_vertices.size(), -1
                );
                GEO::MeshMutator::facet_regions(M).swap(triangle_regions);
                GEO::MeshMutator::tet_regions(M).swap(tet_regions);
                GEO::MeshMutator::set_attributes(
                    M, MeshAttributes(
                        M.attributes() | MESH_FACET_REGION | MESH_CELL_REGION
                    )
                );
                M.update_cached_variables();
                if(cell_borders_only) {
                    mesh_repair(M, MESH_REPAIR_TOPOLOGY);
                } else {
                    M.connect_facets();
                }
            } else {
                RVDMeshBuilder builder(
                    &M, mesh_, delaunay_
                );
                if(dim != 0) {
                    builder.set_dimension(dim);
                }
                RVD_.for_each_polygon(
                    BuildRVD<RVDMeshBuilder>(RVD_, builder)
                );
            }
            RVD_.set_symbolic(sym);
            M.show_stats("RVD");
        }

        /********************************************************************/

        /**
         * \brief Does the actual computation for a specific part
         *    in multithread mode.
         * \param[in] t the index of the part.
         * \pre \p t < nb_parts()
         */
        void run_thread(index_t t) {
            geo_assert(t < nb_parts());
            thisclass& T = part(t);
            switch(thread_mode_) {
                case MT_LLOYD:
                {
                    T.compute_centroids(arg_vectors_, arg_scalars_);
                } break;
                case MT_NEWTON:
                {
                    T.compute_CVT_func_grad(T.funcval_, arg_vectors_);
                } break;
                case MT_INT_SMPLX:
                {
                    T.compute_integration_simplex_func_grad(
                        T.funcval_, arg_vectors_, simplex_func_
                    );
                } break;
                default:
                {
                    geo_assert_not_reached;
                } break;
            }
        }

        /**
         * \brief Removes the duplicated items in a vector.
         * \details In addition, this function sorts the vector.
         * \param[in] V the vector from which duplicates should be removed
         */
        template <class T>
        inline void remove_duplicates(vector<T>& V) {
            std::sort(V.begin(), V.end());
            V.erase(std::unique(V.begin(), V.end()), V.end());
        }

        virtual bool compute_initial_sampling_on_surface(
            double* p, index_t nb_points
        ) {
            geo_assert(mesh_->is_triangulated());

            // We do that here, since this triggers partitioning,
            // that improves data locality. Then data locality is
            // inherited by the generated points.
            create_threads();

            if(facets_begin_ == -1 && facets_end_ == -1) {
                Logger::out("RVD")
                    << "Computing initial sampling on surface, using dimension="
                    << index_t(dimension_) << std::endl;
            }

            return mesh_generate_random_samples_on_surface<DIM>(
                *mesh_, p, nb_points, true, facets_begin_, facets_end_
            );
        }

        virtual bool compute_initial_sampling_in_volume(
            double* p, index_t nb_points
        ) {
            geo_assert(mesh_->nb_tets() != 0);

            // We do that here, since this triggers partitioning,
            // that improves data locality. Then data locality is
            // inherited by the generated points.
            create_threads();

            if(tets_begin_ == -1 && tets_end_ == -1) {
                Logger::out("RVD")
                    << "Computing initial sampling in volume, using dimension="
                    << index_t(dimension_) << std::endl;
            }

            return mesh_generate_random_samples_in_volume<DIM>(
                *mesh_, p, nb_points, true, tets_begin_, tets_end_
            );
        }

        /**
         * \brief Creates the data structures for fast projection.
         * \details It decomposes the surface into triangles, and
         * stores the vertices in a KdTree.
         */
        void prepare_projection() {
            if(!mesh_vertices_.is_nil()) {
                return;
            }

            // Step 1: get triangles
            for(index_t f = 0; f < mesh_->nb_facets(); f++) {
                index_t i = mesh_->facet_begin(f);
                for(index_t j = i + 1; j + 1 < mesh_->facet_end(f); j++) {
                    triangles_.push_back(mesh_->corner_vertex_index(i));
                    triangles_.push_back(mesh_->corner_vertex_index(j));
                    triangles_.push_back(mesh_->corner_vertex_index(j + 1));
                }
            }
            nb_triangles_ = index_t(triangles_.size() / 3);

            // Step 2: get vertices stars
            //    Step 2.1: get one-ring neighborhood
            vector<vector<index_t> > stars2(mesh_->nb_vertices());
            for(index_t t = 0; t < nb_triangles_; t++) {
                stars2[triangles_[3 * t]].push_back(t);
                stars2[triangles_[3 * t + 1]].push_back(t);
                stars2[triangles_[3 * t + 2]].push_back(t);
            }

            //   Step 2.2: get two-ring neighborhood
            stars_.resize(mesh_->nb_vertices());
            for(index_t i = 0; i < stars2.size(); i++) {
                vector<index_t> Ni;
                for(index_t j = 0; j < stars2[i].size(); j++) {
                    index_t t = stars2[i][j];
                    for(index_t iv = 0; iv < 3; iv++) {
                        index_t v = triangles_[3 * t + iv];
                        if(v != i) {
                            Ni.push_back(v);
                        }
                    }
                }
                remove_duplicates(Ni);
                for(index_t j = 0; j < Ni.size(); j++) {
                    index_t k = Ni[j];
                    stars_[i].insert(
                        stars_[i].end(), stars2[k].begin(), stars2[k].end()
                    );
                }
                remove_duplicates(stars_[i]);
            }

            // Step 3: create search structure
            mesh_vertices_ = Delaunay::create(dimension_, "NN");
            index_t nb_vertices = mesh_->nb_vertices();
            vector<double> mesh_vertices(nb_vertices * dimension_);
            for(index_t i = 0; i < nb_vertices; i++) {
                for(index_t coord = 0; coord < dimension_; coord++) {
                    mesh_vertices[i * dimension_ + coord] =
                        mesh_->vertex_ptr(i)[coord];
                }
            }
            mesh_vertices_->set_vertices(
                nb_vertices, &mesh_vertices[0]
            );
        }

        virtual void project_points_on_surface(
            index_t nb_points, double* points, vec3* nearest, bool do_project
        ) {

            prepare_projection();

            if(use_exact_projection_) {
                for(index_t p = 0; p < nb_points; p++) {
                    Point P(points + p * dimension_);
                    double d2 = Numeric::max_float64();
                    for(index_t t = 0; t < nb_triangles_; t++) {
                        const Point& p1 = mesh_vertex(triangles_[3 * t]);
                        const Point& p2 = mesh_vertex(triangles_[3 * t + 1]);
                        const Point& p3 = mesh_vertex(triangles_[3 * t + 2]);

                        double l1, l2, l3;
                        Point nearestP;

                        double cur_d2 = Geom::point_triangle_squared_distance(
                            P, p1, p2, p3, nearestP, l1, l2, l3
                        );

                        if(cur_d2 < d2) {
                            d2 = cur_d2;
                            const vec3& p1_R3 = R3_embedding(triangles_[3 * t]);
                            const vec3& p2_R3 = R3_embedding(triangles_[3 * t + 1]);
                            const vec3& p3_R3 = R3_embedding(triangles_[3 * t + 2]);
                            nearest[p] = l1 * p1_R3 + l2 * p2_R3 + l3 * p3_R3;
                            if(do_project) {
                                for(coord_index_t
                                    coord = 0; coord < dimension_; coord++) {
                                    (points + p * dimension_)[coord] =
                                        nearestP[coord];
                                }
                            }
                        }
                    }
                }
                return;
            }

            // find nearest point on surface in star of nearest vertex
            for(index_t p = 0; p < nb_points; p++) {
                Point P(points + p * dimension_);
                index_t v = mesh_vertices_->nearest_vertex(
                    points + p * dimension_
                );
                double d2 = Numeric::max_float64();
                nearest[p] = R3_embedding(v);
                for(index_t i = 0; i < stars_[v].size(); i++) {
                    index_t t = stars_[v][i];
                    const Point& p1 = mesh_vertex(triangles_[3 * t]);
                    const Point& p2 = mesh_vertex(triangles_[3 * t + 1]);
                    const Point& p3 = mesh_vertex(triangles_[3 * t + 2]);

                    double l1, l2, l3;
                    Point nearestP;
                    double cur_d2 = Geom::point_triangle_squared_distance(
                        P, p1, p2, p3, nearestP, l1, l2, l3
                    );
                    if(cur_d2 < d2) {
                        d2 = cur_d2;
                        const vec3& p1_R3 = R3_embedding(triangles_[3 * t]);
                        const vec3& p2_R3 = R3_embedding(triangles_[3 * t + 1]);
                        const vec3& p3_R3 = R3_embedding(triangles_[3 * t + 2]);
                        nearest[p] = l1 * p1_R3 + l2 * p2_R3 + l3 * p3_R3;
                        if(do_project) {
                            for(coord_index_t coord = 0;
                                coord < dimension_; coord++
                            ) {
                                (points + p * dimension_)[coord] =
                                    nearestP[coord];
                            }
                        }
                    }
                }
            }
        }

        /********************************************************************/

        /**
         * \brief Implementation class for computing the restricted Delaunay
         *    triangulation.
         * \details To be used as a template argument
         *    to RVD::for_each_primal_triangle().
         */
        class GetPrimalTriangles {
        public:
            /**
             * \brief Creates a new GetPrimalTriangles.
             * \param[out] triangles where to store the triangles
             */
            GetPrimalTriangles(
                vector<index_t>& triangles
            ) :
                triangles_(triangles) {
            }

            /**
             * \brief The callback called for each primal triangle.
             * \param[in] v1 index of the first vertex
             * \param[in] v2 index of the second vertex
             * \param[in] v3 index of the third vertex
             */
            void operator() (index_t v1, index_t v2, index_t v3) {
                triangles_.push_back(v1);
                triangles_.push_back(v2);
                triangles_.push_back(v3);
            }

        private:
            vector<index_t>& triangles_;
        };

        // TODO: re-implement 'project' mode

        /**
         * \brief Implementation class for computing the restricted Delaunay
         * triangulation of the connected components.
         *
         * \details The difference with GetPrimalTriangles is that when a
         * restricted Voronoi cell has multiple connected components,
         * more triangles are generated to account for the topology.
         * To be used as a template argument to RVD::for_each_polygon().
         * The RestrictedVoronoiDiagram needs to be in connected-components
         *  priority mode.
         */
        class GetConnectedComponentsPrimalTriangles {
        public:
            /** Internal representation of the polygons. */
            typedef typename GenRestrictedVoronoiDiagram::Polygon Polygon;

            /** Internal representation of the vertices. */
            typedef typename GenRestrictedVoronoiDiagram::Vertex Vertex;

            static const index_t UNINITIALIZED = index_t(-1);
            static const index_t MULTI_COMP    = index_t(-2);

            /**
             * \brief Constructs a new GetConnectedComponentsPrimalTriangles.
             * \param[in] RVD the restricted Voronoi diagram
             * \param[out] triangles where to store the triangles
             * \param[out] vertices where to store the vertices
             * \param[in] dimension dimension of the restricted Voronoi diagram
             * \param[in] project if true, projects onto the input mesh
             * \param[in] use_RVC_centroids if true, replaces vertices with
             *    the centroids of the restricted Voronoi cells
             * \param[in] seed_is_locked indicates for each seed whether it
             *    is locked. Locked seeds are not moved, projected or
             *    replaced by centroids. If left empty, then no seed is locked.
             * \param[out] vertex_to_seed a pointer to a vector that stores for
             *    each generated vertex the index of the seed it corresponds to,
             *    or nil if unused
             * \param[in] keep_original_seeds if true, whenever a restricted
             *   Voronoi cell has a single connected component, the seed is
             *   used instead of the centroid of the component
             */
            GetConnectedComponentsPrimalTriangles(
                const GenRestrictedVoronoiDiagram& RVD,
                vector<index_t>& triangles,
                vector<double>& vertices,
                coord_index_t dimension,
                RDTMode mode,
                const std::vector<bool>& seed_is_locked
            ) :
                RVD_(RVD),
                dimension_(dimension),
                triangles_(triangles),
                vertices_(vertices),
                m_(0.0),
                cur_seed_(-1),
                cur_vertex_(0),
                use_RVC_centroids_((mode & RDT_RVC_CENTROIDS) != 0),
                seed_is_locked_(seed_is_locked),
                prefer_seeds_((mode & RDT_PREFER_SEEDS) != 0) {
                if(prefer_seeds_) {
                    seed_to_vertex_.assign(
                        RVD.delaunay()->nb_vertices(),index_t(UNINITIALIZED)
                    );
                }
            }
            
            /**
             * \brief The callback called for each restricted Voronoi cell.
             * \param[in] s1 index of current center vertex
             * \param[in] P current restricted Voronoi cell
             */
            void operator() (index_t s1, const Polygon& P) {
                if(RVD_.connected_component_changed()) {
                    if(cur_seed_ != -1) {
                        end_connected_component();
                    }
                    begin_connected_component(s1);
                }

                // Accumulate mass and barycenter
                index_t vbase = cur_vertex_ * dimension_;
                for(index_t i = 1; i + 1 < P.nb_vertices(); ++i) {
                    double cur_m = Geom::triangle_area(
                        P.vertex(0).point(),
                        P.vertex(i).point(),
                        P.vertex(i + 1).point(), dimension_
                    );
                    for(coord_index_t c = 0; c < dimension_; ++c) {
                        vertices_[vbase + c] += cur_m / 3.0 * (
                            P.vertex(0).point()[c] +
                            P.vertex(i).point()[c] +
                            P.vertex(i + 1).point()[c]
                        );
                    }
                    m_ += cur_m;
                }

                // Detect Voronoi vertices and generate
                // triangles.
                // Note: they can be generated several times,
                // since we cannot know in advance whether
                // the other instances of the Voronoi vertex
                // were finalized or not (i.e. have their
                // three vertices ready).
                // Duplicate triangles are then filtered-out
                // by client code.
                for(index_t i = 0; i < P.nb_vertices(); ++i) {
                    const Vertex& V = P.vertex(i);
                    if(V.sym().nb_bisectors() == 2) {
                        index_t s2 = V.sym().bisector(0);
                        index_t s3 = V.sym().bisector(1);
                        index_t f = V.sym().boundary_facet(0);

                        index_t v1 = RVD_.current_connected_component();
                        signed_index_t v2 =
                            RVD_.get_facet_seed_connected_component(f,s2);
                        signed_index_t v3 =
                            RVD_.get_facet_seed_connected_component(f,s3);

                        if(v2 >= 0 && v3 >= 0) {
                            triangles_.push_back(v1);
                            triangles_.push_back(index_t(v2));
                            triangles_.push_back(index_t(v3));
                        }
                    }
                }
            }

            /**
             * \brief The destructor
             * \details Terminates the current connected component.
             */
            ~GetConnectedComponentsPrimalTriangles() {
                if(cur_seed_ != -1) {
                    end_connected_component();
                }
                if(prefer_seeds_) {
                    
                    // Construct an axis-aligned bounding box tree,
                    // do not reorder the mesh (needs to be pre-reordered)
                    MeshFacetsAABB AABB(*const_cast<Mesh*>(RVD_.mesh()),false);
                    
                    for(index_t s=0; s<seed_to_vertex_.size(); ++s) {
                        if(
                           seed_to_vertex_[s] != MULTI_COMP && 
                           seed_to_vertex_[s] != UNINITIALIZED
                        ) {
                            index_t vbase = seed_to_vertex_[s] * dimension_;

                            const double* seed_ptr = 
                                RVD_.delaunay()->vertex_ptr(s);

                            const double* vertex_ptr = &(vertices_[vbase]);

                            //  If the seed is nearer to the surface than the
                            // centroid of the connected component of the
                            // restricted Voronoi cell, then use the seed.
                            if(
                                AABB.squared_distance(vec3(seed_ptr)) <
                                AABB.squared_distance(vec3(vertex_ptr))
                            ) {
                                for(coord_index_t c = 0; c < dimension_; ++c) {
                                    vertices_[vbase + c] = seed_ptr[c];
                                }
                            }
                        }
                    }
                }
            }

        protected:
            /**
             * \brief Tests whether a given seed is locked.
             */
            bool seed_is_locked(index_t s) {
                return
                    seed_is_locked_.size() > 0 &&
                    seed_is_locked_[s]
                ;
            }

            /**
             * \brief Starts a new connected component.
             * \param[in] s the seed the connected component
             *    is associated with.
             */
            void begin_connected_component(index_t s) {
                cur_seed_ = signed_index_t(s);
                for(coord_index_t c = 0; c < dimension_; ++c) {
                    vertices_.push_back(0.0);
                }
                m_ = 0.0;
            }

            /**
             * \brief Terminates the current connected component.
             */
            void end_connected_component() {
                if(
                    !use_RVC_centroids_ ||
                    seed_is_locked(index_t(cur_seed_))
                ) {
                    // Copy seed
                    index_t vbase = cur_vertex_ * dimension_;
                    const double* seed_ptr =
                        RVD_.delaunay()->vertex_ptr(index_t(cur_seed_));
                    for(coord_index_t c = 0; c < dimension_; ++c) {
                        vertices_[vbase + c] = seed_ptr[c];
                    }
                } else {
                    // Use restricted Voronoi
                    // cell component's centroid.
                    double scal = (m_ < 1e-30 ? 0.0 : 1.0 / m_);
                    index_t vbase = cur_vertex_ * dimension_;
                    for(coord_index_t c = 0; c < dimension_; ++c) {
                        vertices_[vbase + c] *= scal;
                    }
                }
                if(prefer_seeds_) {
                    if(seed_to_vertex_[index_t(cur_seed_)] == UNINITIALIZED) {
                        seed_to_vertex_[index_t(cur_seed_)] = cur_vertex_;
                    } else {
                        seed_to_vertex_[index_t(cur_seed_)] = MULTI_COMP;
                    }
                }
                ++cur_vertex_;
            }

        private:
            const GenRestrictedVoronoiDiagram& RVD_;
            coord_index_t dimension_;
            vector<index_t>& triangles_;
            vector<double>& vertices_;
            double m_;
            signed_index_t cur_seed_;
            index_t cur_vertex_;
            bool project_;   // TODO: not implemented yet.
            bool use_RVC_centroids_;
            const std::vector<bool>& seed_is_locked_;
            bool prefer_seeds_;
            vector<index_t> seed_to_vertex_;
        };

        /**
         * \brief Implementation class for computing the restricted Delaunay
         *    triangulation in volume mode.
         * \details To be used as a template argument
         *    to RVD::for_each_primal_tetrahedron().
         */
        class GetPrimalTetrahedra {
        public:
            /**
             * \brief Creates a new GetPrimalTetrahedra.
             * \param[out] tetrahedra where to store the tetrahedra
             */
            GetPrimalTetrahedra(
                vector<index_t>& tetrahedra
            ) :
                tetrahedra_(tetrahedra) {
            }

            /**
             * \brief The callback called for each primal tetrahedron.
             * \param[in] v1 index of the first vertex
             * \param[in] v2 index of the second vertex
             * \param[in] v3 index of the third vertex
             * \param[in] v4 index of the fourth vertex
             */
            void operator() (index_t v1, index_t v2, index_t v3, index_t v4) {
                tetrahedra_.push_back(v1);
                tetrahedra_.push_back(v2);
                tetrahedra_.push_back(v3);
                tetrahedra_.push_back(v4);
            }

        private:
            vector<index_t>& tetrahedra_;
        };

        virtual void compute_RDT(
            vector<index_t>& simplices,
            vector<double>& embedding,
            RDTMode mode,
            const vector<bool>& seed_is_locked
        ) {
            if(volumetric_) {
                // For the moment, only simple mode is supported
                simplices.clear();
                RVD_.for_each_primal_tetrahedron(
                    GetPrimalTetrahedra(simplices)
                );
                // Reorient the tetrahedra
                index_t nb_tetrahedra = simplices.size() / 4;
                for(index_t t = 0; t < nb_tetrahedra; ++t) {
                    // const-cast is required because orient3d() is 
                    //  declared in C and does
                    // not know about constness...
                    double* p1 = const_cast<double*>(
                        delaunay()->vertex_ptr(simplices[4 * t])
                    );
                    double* p2 = const_cast<double*>(
                        delaunay()->vertex_ptr(simplices[4 * t + 1])
                    );
                    double* p3 = const_cast<double*>(
                        delaunay()->vertex_ptr(simplices[4 * t + 2])
                    );
                    double* p4 = const_cast<double*>(
                        delaunay()->vertex_ptr(simplices[4 * t + 3])
                    );
                    if(orient3d(p1, p2, p3, p4) > 0) {
                        geo_swap(simplices[4 * t], simplices[4 * t + 1]);
                    }
                }
                embedding.clear();
                embedding.reserve(dimension_ * delaunay_->nb_vertices());
                for(index_t i = 0; i < delaunay_->nb_vertices(); i++) {
                    for(coord_index_t coord = 0; coord < dimension_; coord++) {
                        embedding.push_back(delaunay_->vertex_ptr(i)[coord]);
                    }
                }
            } else {
                if((mode & RDT_MULTINERVE) != 0) {
                    simplices.clear();
                    embedding.clear();
                    bool sym = RVD_.symbolic();
                    RVD_.set_symbolic(true);
                    RVD_.set_connected_components_priority(true);
                    RVD_.for_each_polygon(
                        GetConnectedComponentsPrimalTriangles(
                            RVD_, simplices, embedding, RVD_.dimension(),
                            mode, seed_is_locked
                        )
                    );
                    RVD_.set_symbolic(sym);
                    RVD_.set_connected_components_priority(false);
                } else {
                    // Simple mode: compute RDT, without any post-processing
                    simplices.clear();
                    RVD_.for_each_primal_triangle(
                        GetPrimalTriangles(simplices)
                    );
                    embedding.clear();
                    embedding.reserve(dimension_ * delaunay_->nb_vertices());
                    for(index_t i = 0; i < delaunay_->nb_vertices(); i++) {
                        for(
                            coord_index_t coord = 0; 
                            coord < dimension_; ++coord
                        ){
                             embedding.push_back(
                                 delaunay_->vertex_ptr(i)[coord]
                             );
                        }
                    }
                }
            }
        }

        virtual void create_threads() {
            // TODO: check if number of facets is not smaller than
            // number of threads
            // TODO: create parts even if facets range is specified
            // (and subdivide facets range)
            if(is_slave_ || facets_begin_ != -1 || facets_end_ != -1) {
                return;
            }
            index_t nb_parts_in = Process::maximum_concurrent_threads();
            if(nb_parts() != nb_parts_in) {
                if(nb_parts_in == 1) {
                    delete_threads();
                } else {
                    vector<index_t> facet_ptr;
                    vector<index_t> tet_ptr;
                    Logger::out("RVD")
                        << "Creating "
                        << nb_parts_in << " parts" << std::endl;
                    mesh_partition(
                        *mesh_, MESH_PARTITION_HILBERT,
                        facet_ptr, tet_ptr, nb_parts_in
                    );
                    delete_threads();
                    parts_ = new thisclass[nb_parts_in];
                    nb_parts_ = nb_parts_in;
                    for(index_t i = 0; i < nb_parts(); ++i) {
                        part(i).mesh_ = mesh_;
                        part(i).set_delaunay(delaunay_);
                        part(i).R3_embedding_base_ = R3_embedding_base_;
                        part(i).R3_embedding_stride_ = R3_embedding_stride_;
                        part(i).has_weights_ = has_weights_;
                        part(i).master_ = this;
                        part(i).RVD_.set_mesh(mesh_);
                        part(i).set_facets_range(
                            facet_ptr[i], facet_ptr[i + 1]
                        );
                        part(i).set_exact_predicates(RVD_.exact_predicates());
                        part(i).set_volumetric(volumetric());
                    }
                    if(mesh_->nb_tets() != 0) {
                        for(index_t i = 0; i < nb_parts(); ++i) {
                            part(i).set_tetrahedra_range(
                                tet_ptr[i], tet_ptr[i + 1]
                            );
                        }
                    }
                    geo_assert(!Process::is_running_threads());
                }
            }
        }

        virtual void set_volumetric(bool x) {
            volumetric_ = x;
            for(index_t i = 0; i < nb_parts(); ++i) {
                part(i).set_volumetric(x);
            }
        }

        virtual void set_facets_range(
            index_t facets_begin, index_t facets_end
        ) {
            RVD_.set_facets_range(facets_begin, facets_end);
            facets_begin_ = signed_index_t(facets_begin);
            facets_end_ = signed_index_t(facets_end);
        }

        virtual void set_tetrahedra_range(
            index_t tets_begin, index_t tets_end
        ) {
            RVD_.set_tetrahedra_range(tets_begin, tets_end);
            tets_begin_ = signed_index_t(tets_begin);
            tets_end_ = signed_index_t(tets_end);
        }

        virtual void delete_threads() {
            delete[] parts_;
            parts_ = nil;
            nb_parts_ = 0;
        }

        /**
         * \brief Gets the number of parts (or number of threads).
         */
        index_t nb_parts() const {
            return nb_parts_;
        }

        /**
         * \brief Gets a given part from its index.
         * \param[in] i index of the part
         * \pre \p i < nb_parts()
         */
        thisclass& part(index_t i) {
            geo_debug_assert(i < nb_parts());
            return parts_[i];
        }

    protected:

        GenRestrictedVoronoiDiagram RVD_;

        // For projection
        bool use_exact_projection_;
        index_t nb_triangles_;
        vector<index_t> triangles_;
        vector<vector<index_t> > stars_;
        Delaunay_var mesh_vertices_;

        // One of MT_NONE, MT_LLOYD, MT_NEWTON
        ThreadMode thread_mode_;

        bool is_slave_;

        // Variables for 'master' in multithreading mode
        thisclass* parts_;
        index_t nb_parts_;
        Process::SpinLockArray spinlocks_;

        // Newton mode with int. simplex
        IntegrationSimplex* simplex_func_; 

        // master stores argument for compute_centroids() and
        // compute_CVT_func_grad() to pass it to the parts.
        double* arg_vectors_;
        double* arg_scalars_;

        // Variables for 'slaves' in multithreading mode
        thisclass* master_;
        double funcval_;  // Newton mode: function value

    protected:
        /**
         * \brief Destructor
         */
        virtual ~RVD_Nd_Impl() {
            delete_threads();
        }

    private:
        /** \brief Forbids construction by copy. */
        RVD_Nd_Impl(const thisclass&);

        /** \brief Forbids assignment. */
        thisclass& operator= (const thisclass&);
    };
}

/****************************************************************************/

namespace GEO {

    RestrictedVoronoiDiagram* RestrictedVoronoiDiagram::create(
        Delaunay* delaunay, Mesh* mesh,
        const double* R3_embedding, index_t R3_embedding_stride
    ) {
        RestrictedVoronoiDiagram* result = nil;
        geo_assert(delaunay != nil);
        coord_index_t dim = delaunay->dimension();
        switch(dim) {
            case 2:
                result = new RVD_Nd_Impl<2>(
                    delaunay, mesh, R3_embedding, R3_embedding_stride
                );
                break;
            case 3:
                result = new RVD_Nd_Impl<3>(
                    delaunay, mesh, R3_embedding, R3_embedding_stride
                );
                break;
            case 4:
                result = new RVD_Nd_Impl<4>(
                    delaunay, mesh, R3_embedding, R3_embedding_stride
                );
                break;
            case 6:
                result = new RVD_Nd_Impl<6>(
                    delaunay, mesh, R3_embedding, R3_embedding_stride
                );
                break;
            default:
                geo_assert_not_reached;
                return nil;
                break;
        }
        if(CmdLine::get_arg("algo:predicates") == "exact") {
            Logger::out("RVD") << "Using exact predicates" << std::endl;
            result->set_exact_predicates(true);
        } else {
            Logger::out("RVD") << "Using fast predicates" << std::endl;
        }
        return result;
    }

    void RestrictedVoronoiDiagram::set_delaunay(Delaunay* delaunay) {
        delaunay_ = delaunay;
        if(delaunay_ != nil) {
            dimension_ = delaunay->dimension();
        } else {
            dimension_ = 0;
        }
    }

    RestrictedVoronoiDiagram::~RestrictedVoronoiDiagram() {
    }

    RestrictedVoronoiDiagram::RestrictedVoronoiDiagram(
        Delaunay* delaunay, Mesh* mesh,
        const double* R3_embedding, index_t R3_embedding_stride
    ) :
        dimension_(0),
        mesh_(mesh),
        R3_embedding_base_(R3_embedding),
        R3_embedding_stride_(R3_embedding_stride) {
        set_delaunay(delaunay);
        has_weights_ = false;
        facets_begin_ = -1;
        facets_end_ = -1;
        tets_begin_ = -1;
        tets_end_ = -1;
        volumetric_ = false;
    }


    void RestrictedVoronoiDiagram::compute_RDT(
        Mesh& RDT,
        RDTMode mode,
        const vector<bool>& seed_is_locked
    ) {
        vector<index_t> simplices;
        vector<double> embedding;
        compute_RDT(
            simplices, embedding, 
            mode, seed_is_locked
        );
        if(volumetric()) {
            RDT.assign_tet_mesh(dimension(),embedding,simplices,true);
        } else {
            RDT.assign_triangle_mesh(dimension(),embedding,simplices,true);
            mesh_repair(RDT); // Needed to reorient triangles
        }
    }
    
}

