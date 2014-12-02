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

#ifndef __GEOGRAM_MESH_MESH_BUILDER__
#define __GEOGRAM_MESH_MESH_BUILDER__

#include <geogram/basic/common.h>
#include <geogram/mesh/mesh.h>

/**
 * \file geogram/mesh/mesh_builder.h
 * \brief The incremental builder for Mesh
 */

namespace GEOGen {

    /**
     * \brief Incremental builder for Mesh.
     */
    template <class MESH>
    class MeshBuilder {

        /** \brief This class type */
        typedef MeshBuilder<MESH> thisclass;

    public:
        /**
         * \brief The type of the mesh coordinates.
         * \details Double (default) or float (to be used
         *  with SinglePrecisionMesh).
         */
        typedef typename MESH::coord_t coord_t;

        /**
         * \brief Constructs a new MeshBuilder.
         * \param[in] target the target Mesh to be built
         */
        MeshBuilder(
            MESH* target = nil
        ) :
            target_(target),
            state_(NONE) {
        }

        /**
         * \brief Sets the target mesh to be built.
         * \param[in] m a pointer to the mesh to be built
         */
        void set_target(MESH* m) {
            target_ = m;
        }

        /**
         * \brief Gets the target mesh.
         * \return a pointer to the target mesh
         */
        MESH* target() {
            return target_;
        }

        /**
         * \brief Reserves memory in advance for the facets.
         * \details If the number of facets is known in advance,
         *  pre-reserving the memory improves perfomances by
         *  avoiding reallocations.
         * \param[in] nb_f the number of facets that will be created
         */
        void reserve_facets(index_t nb_f) {
            target_->reserve_facets(nb_f);
        }

        /**
         * \brief Starts a new mesh.
         * \details This clears the target.
         */
        void begin_mesh() {
            target_->clear();
            geo_assert(state_ == NONE);
            state_ = IN_SURFACE;
        }

        /**
         * \brief Creates all the vertices of the mesh.
         * \details This overwrites the vertices of the target
         *  mesh. Alternatively, vertices can be created
         *  one by one using add_vertex() or add_vertex_by_ptr().
         * \param[in] nb_vertices number of vertices to create
         * \param[in] dimension dimension of the target mesh
         */
        void create_vertices(index_t nb_vertices, coord_index_t dimension) {
            target_->create_vertices(nb_vertices, dimension);
        }

        /**
         * \brief Creates a new vertex.
         * \details The dimension of the mesh needs to be known. Alternatively
         *  the other form of add_vertex_by_ptr() that takes the vertex
         *  dimension as an argument can be used instead.
         * \param[in] p pointer to the coordinates of the vertex
         * \param[in] w weight associated with the vertex
         */
        void add_vertex_by_ptr(const double* p, double w = 1.0) {
            geo_assert(state_ == IN_SURFACE || state_ == IN_FACET);
            for(index_t i = 0; i < target_->dimension(); i++) {
                target_->vertices_.push_back(coord_t(p[i]));
            }
            if(w != 1.0 && !target_->has_weights()) {
                target_->create_vertices_weights();
            }
            if(target_->has_weights()) {
                target_->weights_.push_back(coord_t(w));
            }
            target_->nb_vertices_ =
                target_->vertices_.size() / target_->dimension();
        }

        /**
         * \brief Creates a new vertex.
         * \param[in] p pointer to the coordinates of the vertex
         * \param[in] w weight associated with the vertex
         * \param[in] dim dimension of the vertex. The first created
         *  vertex determines the dimension of the mesh.
         */
        void add_vertex_by_ptr(const double* p, double w, coord_index_t dim) {
            // The first added point determines the dimension.
            if(target_->dimension_ == 0) {
                target_->dimension_ = dim;
            }
            geo_debug_assert(dim == target_->dimension());
            add_vertex_by_ptr(p, w);
        }

        /**
         * \brief Creates a new vertex.
         * \param[in] P the coordinates of the vertex.
         *  The first added point determines the dimension of the mesh.
         * \param[in] w weight associated with the vertex
         */
        template <class POINT>
        void add_vertex(const POINT& P, double w = 1.0) {

            if(target_->dimension_ == 0) {
                GEO::Logger::out("Mesh")
                    << "Set dimension to " << P.dimension() << std::endl;
                target_->dimension_ = coord_index_t(P.dimension());
            }
            geo_debug_assert(P.dimension() == target_->dimension());
            add_vertex_by_ptr(P.data(), w);
        }

        /**
         * \brief Begins a new facet.
         * \pre this function needs to be called between a begin_surface() /
         *  end_surface() pair.
         */
        void begin_facet() {
            geo_assert(state_ == IN_SURFACE);
            state_ = IN_FACET;
            target_->begin_facet();
        }

        /**
         * \brief Adds a vertex to the current facet.
         * \param[in] v index of the vertex to be added to the current facet
         * \pre this function needs to be called between a begin_facet() /
         *  end_facet() pair.
         */
        void add_vertex_to_facet(index_t v) {
            geo_assert(state_ == IN_FACET);
            geo_debug_assert(v < target_->nb_vertices());
            target_->add_corner(v);
        }

        /**
         * \brief Terminates the current facet.
         * \param[in] facet_region optional region index associated with
         *  the facet. It is stored if the target has the facet region
         *  attribute.
         */
        void end_facet(signed_index_t facet_region = 0) {
            geo_assert(state_ == IN_FACET);
            state_ = IN_SURFACE;
            target_->end_facet(facet_region);
        }

        /**
         * \brief Terminates the current mesh.
         * \param[in] verbose if true, statistics are displayed.
         * \param[in] connect if true, the adjacencies between
         *  mesh elements (facets and tets) are computed.
         */
        void end_mesh(bool verbose = true, bool connect = true) {
            geo_debug_assert(state_ == IN_SURFACE);
            state_ = NONE;
            if(connect) {
                target_->connect_facets();
                target_->connect_cells();
                target_->remove_degree2_vertices();
            }
            target_->update_cached_variables();
            if(verbose) {
                target_->show_stats();
            }
        }

    private:
        /**
         * \brief Encodes the current state of the mesh builder.
         * \details Used for sanity checks (assertions that the
         *  correct calling sequence is respected).
         */
        enum State {
            NONE,
            IN_SURFACE,
            IN_FACET
        };

        MESH* target_;
        State state_;
    };
}

namespace GEO {

    /**
     * The MeshBuilder for the default Mesh class (with coordinates as doubles).
     */
    typedef GEOGen::MeshBuilder<Mesh> MeshBuilder;
}

#endif

