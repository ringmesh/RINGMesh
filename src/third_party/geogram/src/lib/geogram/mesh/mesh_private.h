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

#ifndef __GEOGRAM_MESH_MESH_PRIVATE__
#define __GEOGRAM_MESH_MESH_PRIVATE__

#include <geogram/basic/common.h>
#include <geogram/mesh/mesh.h>

/**
 * \file geogram/mesh/mesh_private.h
 * \brief Direct access to Mesh members
 */

namespace GEOGen {

    /**
     * \brief Used internally by algorithms that modify a mesh.
     * \note Provides low-level access, without any checks, for
     *  internal use in Vorpaline only (or use at your own risk !).
     */
    template <class MESH>
    class MeshMutator {
    public:
        /**
         * \brief The type of the mesh coordinates.
         * \details Double (default) or float (to be used
         *  with SinglePrecisionMesh).
         */
        typedef typename MESH::coord_t coord_t;

        /**
         * \brief Flips the order of the vertices in a facet.
         * \param[in,out] M the mesh
         * \param[in] f the index of the facet
         */
        static void flip_facet(MESH& M, index_t f) {
            index_t d = M.facet_size(f);

            // Allocated on the stack (more multithread-friendly
            // and no need to free)
            index_t* corner_vertex_index =
                (index_t*) alloca(sizeof(index_t) * d);

            signed_index_t* corner_adjacent_facet =
                (signed_index_t*) alloca(sizeof(signed_index_t) * d);

            index_t c0 = M.facet_begin(f);
            for(index_t i = 0; i < d; i++) {
                corner_vertex_index[i] = M.corner_vertex_index(c0 + i);
                corner_adjacent_facet[i] = M.corner_adjacent_facet(c0 + i);
            }
            for(index_t i = 0; i < d; i++) {
                index_t i_v = d - 1 - i;
                index_t i_f = (i_v == 0) ? d - 1 : i_v - 1;
                M.set_corner_vertex_index(c0 + i, corner_vertex_index[i_v]);
                M.set_corner_adjacent_facet(c0 + i, corner_adjacent_facet[i_f]);
            }
        }

        /**
         * \brief Gets the internal corner vertex index vector of a mesh.
         * \param[in] M the mesh
         * \return a reference to the internal corner vertex index vector
         *  of the mesh \p M
         */
        static GEO::vector<index_t>& corner_vertices(MESH& M) {
            return M.corner_vertices_;
        }

        /**
         * \brief Gets the internal corner adjacent facets vector of a mesh.
         * \param[in] M the mesh
         * \return a reference to the internal ajdancet facets vector of
         *  the mesh \p M
         */
        static GEO::vector<signed_index_t>& corner_adjacent_facets(
            MESH& M
        ) {
            return M.corner_adjacent_facets_;
        }

        /**
         * \brief Gets the internal facet pointer vector of a mesh.
         * \param[in] M the mesh
         * \return a reference to the internal facet pointer vector
         *  of the mesh
         */
        static GEO::vector<index_t>& facet_ptr(MESH& M) {
            return M.facet_ptr_;
        }

        /**
         * \brief Gets the internal vertices coordinates vector of a
         *  mesh.
         * \param[in] M the mesh
         * \return a reference to the internal vector of vertices coordinates
         *  of mesh \p M
         */
        static GEO::vector<coord_t>& vertices(MESH& M) {
            return M.vertices_;
        }

        /**
         * \brief Gets the internal vertices coordinates vector of a
         *  mesh.
         * \param[in] M the mesh
         * \return a const reference reference to the internal vector 
         *  of vertices coordinates of mesh \p M
         */
        static const GEO::vector<coord_t>& vertices(const MESH& M) {
            return M.vertices_;
        }
        
        /**
         * \brief Gets the internal weights vector of a mesh.
         * \param[in] M the mesh
         * \return a reference to the internal weights vector of mesh \p M
         */
        static GEO::vector<coord_t>& weights(MESH& M) {
            return M.weights_;
        }

        /**
         * \brief Gets the internal facets region vector of a mesh.
         * \param[in] M the mesh
         * \return a reference to the internal facets region vector of
         *  mesh \p M
         */
        static GEO::vector<signed_index_t>& facet_regions(MESH& M) {
            return M.facet_regions_;
        }

        /**
         * \brief Gets the internal tetrahedra vertices vector of a mesh.
         * \param[in] M the mesh
         * \return a reference to the internal tetrahedra vertices vector of
         *  mesh \p M
         */
        static GEO::vector<index_t>& tet_vertices(MESH& M) {
            return M.cell_vertices_;
        }

        /**
         * \brief Gets the internal cell vertices vector of a mesh.
         * \param[in] M the mesh
         * \return a reference to the internal cell vertices vector of
         *  mesh \p M
         */
        static GEO::vector<index_t>& cell_vertices(MESH& M) {
            return M.cell_vertices_;
        }
        
        /**
         * \brief Gets the internal tetrahedra adjacency vector of a mesh.
         * \param[in] M the mesh
         * \return a reference to the internal tetrahedra adjacency vector
         *  of mesh \p M
         */
        static GEO::vector<signed_index_t>& tet_adjacents(MESH& M) {
            return M.cell_adjacents_;
        }

        /**
         * \brief Gets the internal tetrahedra region vector of a mesh.
         * \param[in] M the mesh
         * \return a reference to the internal tetrahedra regions vector
         *  of mesh \p M
         */
        static GEO::vector<signed_index_t>& tet_regions(MESH& M) {
            return M.cell_regions_;
        }


        /**
         * \brief Gets the internal cell pointers vector of a mesh.
         * \param[in] M the mesh
         * \return a reference to the internal cell pointers vector
         *  of mesh \p M
         */
        static GEO::vector<index_t>& cell_ptr(MESH& M) {
            return M.cell_ptr_;
        }


        /**
         * \brief Gets the internal cell types vector of a mesh.
         * \param[in] M the mesh
         * \return a reference to the internal cell types vector
         *  of mesh \p M
         */
        static GEO::vector<GEO::Numeric::int8>& cell_types(MESH& M) {
            return M.cell_types_;
        }
        
        /**
         * \brief Sets the number of vertices of a mesh
         * \param[out] M the mesh
         * \param[in] x the number of vertices
         */
        static void set_nb_vertices(MESH& M, index_t x) {
            M.nb_vertices_ = x;
        }

        /**
         * \brief Sets the number of facets of a mesh
         * \param[out] M the mesh
         * \param[in] x the number of facets
         */
        static void set_nb_facets(MESH& M, index_t x) {
            M.nb_facets_ = x;
        }

        /**
         * \brief Sets the number of cells of a mesh
         * \param[out] M the mesh
         * \param[in] x the number of cells
         */
        static void set_nb_cells(MESH& M, index_t x) {
            M.nb_cells_ = x;
        }
        
        /**
         * \brief Specifies whether a mesh is triangulated.
         * \details Only changes the triangulated_ flag, does not touch
         *  anything else in the mesh.
         * \param[out] M the mesh
         * \param[in] x true if the mesh is triangulated, false otherwise
         */
        static void set_triangulated(MESH& M, bool x) {
            M.triangulated_ = x;
        }

        /**
         * \brief Sets the attributes of a mesh.
         * \details Only changes the attributes_ flags, does not touch
         *  anything else in the mesh, does not perform any attributes
         *  allocation/deallocation.
         * \param[out] M the mesh
         * \param[in] x a set of MeshAttributes combined with bitwise or
         */
        static void set_attributes(MESH& M, MeshAttributes x) {
            M.attributes_ = x;
        }

        /**
         * \brief Sets the dimension of a mesh.
         * \details Only changes the dimension_ field in the mesh, does
         *  not perform any points allocation/deallocation.
         * \param[out] M the mesh
         * \param[in] dim the dimension
         */
        static void set_dimension(MESH& M, coord_index_t dim) {
            M.dimension_ = dim;
        }
    };
}

namespace GEO {

    /**
     * \brief The MeshMutator for the default Vorpaline Mesh (with point
     * coordinates as doubles).
     */
    typedef GEOGen::MeshMutator<Mesh> MeshMutator;
}

#endif

