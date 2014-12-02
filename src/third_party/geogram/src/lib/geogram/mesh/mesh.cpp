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
#include <geogram/basic/string.h>
#include <algorithm>

namespace GEO {

    /**
     * \brief Saves a cell to a stream, in Alias|Wavefront format.
     * \details Used for debugging purposes (debugging "a la petite cuilliere").
     * \param[in] M a const reference to the mesh
     * \param[in] c the index of the cell to be saved
     * \param[in,out] index_base the first index to be used for the vertices. 
     *  Incremented on exit so that multiple cells may be output to the same
     *  stream.
     * \param[in] vertices a pointer to the continuous array of coordinates.
     * \param[out] out the stream
     */
    void save_cell(
        const GEOGen::MeshBase& M, index_t c, index_t& index_base,
        const float* vertices, std::ostream& out
    ) {
        for(index_t lv=0; lv<M.cell_nb_vertices(c); ++lv) {
            index_t v = M.cell_vertex_index(c,lv);
            const float* pv = vertices + v*3;
            out << "v " << pv[0] << " " << pv[1] << " " << pv[2] << std::endl;
        }
        for(index_t lf=0; lf<M.cell_nb_facets(c); ++lf) {
            out << " f ";
            for(index_t lv=0; lv<M.cell_facet_nb_vertices(c,lf); ++lv) {
                index_t v = M.cell_facet_vertex_index(c,lf,lv);
                // It's a bit ridiculous, we could directly get
                // it from cell descriptor, but it is protected...
                signed_index_t vv = M.find_cell_vertex(c,v);
                geo_assert(vv != -1);
                out << index_t(vv) + index_base << " ";
            }
            out << std::endl;
        }
        index_base += M.cell_nb_vertices(c);
    }
}


namespace GEOGen {

    MeshBase::CellDescriptor MeshBase::tet_descriptor_ = {
        4,         // nb_vertices
        4,         // nb_facets
        {3,3,3,3}, // nb_vertices in facet
        {          // facets
            {1,3,2},
            {0,2,3},
            {3,1,0},
            {0,1,2}
        },
        6,         // nb_edges
        {          // edges
            {1,2}, {2,3}, {3,1}, {0,1}, {0,2}, {0,3}
        }         
    };


    MeshBase::CellDescriptor MeshBase::hex_descriptor_ = {
        8,             // nb_vertices
        6,             // nb_facets
        {4,4,4,4,4,4}, // nb_vertices in facet
        {              // facets
            {0,2,6,4},
            {3,1,5,7},
            {1,0,4,5},
            {2,3,7,6},
            {1,3,2,0},
            {4,6,7,5}
        },
        12,            // nb_edges
        {              // edges
            {0,1},{1,3},{3,2},{2,0},{4,5},{5,7},
            {7,6},{6,4},{0,4},{1,5},{3,7},{2,6}
        }         
    };

    MeshBase::CellDescriptor MeshBase::prism_descriptor_ = {
        6,             // nb_vertices
        5,             // nb_facets
        {3,3,4,4,4},   // nb_vertices in facet
        {              // facets
            {0,1,2},
            {3,5,4},
            {0,3,4,1},
            {0,2,5,3},
            {1,4,5,2}
        },
        9,             // nb_edges
        {              // edges
            {0,1},{1,2},{2,3},{3,4},{4,5},{5,3},{0,3},{1,4},{2,5}
        }         
    };


    MeshBase::CellDescriptor MeshBase::pyramid_descriptor_ = {
        5,             // nb_vertices
        5,             // nb_facets
        {4,3,3,3,3},   // nb_vertices in facet
        {              // facets
            {0,1,2,3},
            {0,4,1},
            {0,3,4},
            {2,4,3},
            {2,1,4}
        },
        8,             // nb_edges
        {              // edges
            {0,1},{1,2},{2,3},{3,0},{0,4},{1,4},{2,4},{3,4}
        }         
    };

    MeshBase::CellDescriptor MeshBase::connector_descriptor_ = {
        4,             // nb_vertices
        3,             // nb_facets
        {4,3,3},   // nb_vertices in facet
        {              // facets
            {0,1,2,3},
            {2,1,0},
            {3,2,0}
        },
        5,             // nb_edges
        {              // edges
            {0,1},{1,2},{2,3},{3,0},{0,2}
        }         
    };

    MeshBase::CellDescriptor* MeshBase::cell_type_to_cell_descriptor_[5] = { 
        &tet_descriptor_, 
        &hex_descriptor_, 
        &prism_descriptor_, 
        &pyramid_descriptor_, 
        &connector_descriptor_
    };

    MeshBase::MeshBase(coord_index_t dim) :
        nb_facets_(0),
        nb_cells_(0),
        triangulated_(true),
        tetrahedralized_(true),
        nb_vertices_(0),
        dimension_(dim),
        in_facet_(false),
        facets_reserve_(0),
        attributes_(GEO::MESH_NO_ATTRIBUTES)
    {
        // Following two lines for debugging purposes
        debug_coords_float_ = nil; // TODO remove it
        debug_coords_double_ = nil; // TODO remove it        
    }

    void MeshBase::clear(bool keep_memory) {
        geo_assert(!in_facet_);
        
        if(keep_memory) {
            corner_vertices_.resize(0);
            corner_adjacent_facets_.resize(0);
            facet_ptr_.resize(0);
            facet_regions_.resize(0);
            cell_vertices_.resize(0);
            cell_adjacents_.resize(0);
            cell_regions_.resize(0);
            cell_ptr_.resize(0);
            cell_types_.resize(0);
        } else {
            corner_vertices_.clear();
            corner_adjacent_facets_.clear();
            facet_ptr_.clear();
            facet_regions_.clear();
            cell_vertices_.clear();
            cell_adjacents_.clear();
            cell_regions_.clear();
            cell_ptr_.clear();
            cell_types_.clear();
        }
        nb_vertices_ = 0;
        nb_facets_ = 0;
        nb_cells_ = 0;
        triangulated_ = true;
        tetrahedralized_ = true;
        facets_reserve_ = 0;
        attributes_ = GEO::MESH_NO_ATTRIBUTES;
    }

    void MeshBase::reserve_facets_and_corners(
        index_t nb_facets, index_t nb_corners
    ) {
        corner_vertices_.reserve(nb_corners);
        corner_adjacent_facets_.reserve(nb_corners);
        if(!triangulated_) {
            facet_ptr_.reserve(nb_facets + 1);
        }
        if(has_attribute(GEO::MESH_FACET_REGION)) {
            facet_regions_.reserve(nb_facets);
        }
        facets_reserve_ = nb_facets;
    }


    void MeshBase::create_tets(index_t nb_tets) {
        cell_vertices_.resize(nb_tets * 4);
        cell_adjacents_.resize(nb_tets * 4);
        tetrahedralized_ = true;
        nb_cells_ = nb_tets;
        cell_ptr_.clear();
        cell_types_.clear();
    }

    void MeshBase::end_facet(signed_index_t facet_region) {
        geo_assert(in_facet_);

        if(triangulated_) {
            index_t facet_size = corner_vertices_.size() - 3 * nb_facets_;
            if(facet_size != 3) {
                triangulated_ = false;
                facet_ptr_.reserve(facets_reserve_ + 1);
                facet_ptr_.resize(nb_facets_ + 1);
                for(index_t i = 0; i <= nb_facets_; i++) {
                    facet_ptr_[i] = 3 * i;
                }
            }
        }
        // Note: this is not 'else', because triangulated_ may
        // change inside the previous block !!!
        if(!triangulated_) {
            if(*(facet_ptr_.rbegin()) == corner_vertices_.size()) {
                GEO::Logger::warn("Mesh")
                    << "Created empty facet" << std::endl;
            }
            facet_ptr_.push_back((index_t) corner_vertices_.size());
        }
        if(has_attribute(GEO::MESH_FACET_REGION)) {
            facet_regions_.push_back(facet_region);
        }
        in_facet_ = false;
        nb_facets_++;
    }

    void MeshBase::dissociate_facet(index_t f) {
        geo_debug_assert(f < nb_facets());
        for(index_t i = facet_begin(f); i < facet_end(f); i++) {
            set_corner_adjacent_facet(i, -1);
        }
    }

    void MeshBase::remove_degree2_vertices() {
        // Sanity check: detect degree 2 vertices
        std::set<index_t> to_dissociate;
        for(index_t f = 0; f < nb_facets(); f++) {
            for(index_t i1 = facet_begin(f); i1 < facet_end(f); i1++) {
                index_t i2 = i1 + 1;
                if(i2 == facet_end(f)) {
                    i2 = facet_begin(f);
                }
                signed_index_t f1 = corner_adjacent_facet(i1);
                signed_index_t f2 = corner_adjacent_facet(i2);
                if(f1 != -1 && f1 == f2) {
                    to_dissociate.insert(f);
                    to_dissociate.insert(index_t(f1));
                }
            }
        }
        if(!to_dissociate.empty()) {
            GEO::Logger::warn("Mesh")
                << to_dissociate.size()
                << " facets with degree 2 vertices (fixed)"
                << std::endl;
        }
        for(std::set<index_t>::iterator it = to_dissociate.begin();
            it != to_dissociate.end(); ++it) {
            dissociate_facet(*it);
        }
    }

    bool MeshBase::has_edge(
        index_t f,
        index_t v1, index_t v2
    ) const {
        index_t b = facet_begin(f);
        index_t n = facet_end(f) - b;
        for(index_t i = 0; i < n; i++) {
            index_t w1 = corner_vertex_index(b + i);
            index_t w2 = corner_vertex_index(b + ((i + 1) % n));
            if(w1 == v1 && w2 == v2) {
                return true;
            }
        }
        return false;
    }


    bool MeshBase::cell_facets_match(
        index_t c1, index_t f1, index_t c2, index_t f2
    ) const {
        index_t nbv = cell_facet_nb_vertices(c1,f1);
        if(cell_facet_nb_vertices(c2,f2) != nbv) {
            return false;
        }
        for(index_t offset=0; offset<nbv; ++offset) {
            bool match=true;
            for(index_t v1=0; v1<nbv; ++v1) {
                index_t v2 = (nbv-v1+offset)%nbv;
                if(
                    cell_facet_vertex_index(c1,f1,v1) !=
                    cell_facet_vertex_index(c2,f2,v2)
                ) {
                    match=false;
                    break;
                }
            }
            if(match) {
                return true;
            }
        }
        return false;
    }

    index_t MeshBase::add_cell(
        GEO::MeshCellType type, index_t* vertices, signed_index_t* adjacents,
        signed_index_t region
    ) {
        if(type == GEO::MESH_TET && is_tetrahedralized()) {
            for(index_t i=0; i<4; ++i) {
                cell_vertices_.push_back(vertices[i]);
            }
            for(index_t i=0; i<4; ++i) {                    
                cell_adjacents_.push_back(
                    (adjacents == nil) ? -1 : adjacents[i]);
            }
        } else {
            if(is_tetrahedralized()) {
                tetrahedralized_ = false;
                cell_types_.assign(
                    nb_cells(), GEO::Numeric::int8(GEO::MESH_TET)
                );
                cell_ptr_.resize(0);
                for(index_t c=0; c<nb_cells(); ++c) {
                    cell_ptr_.push_back(4*c);
                }
                cell_ptr_.push_back(4*nb_cells());
            }
            cell_types_.push_back(GEO::Numeric::int8(type));
            const CellDescriptor& desc = cell_type_to_cell_descriptor(type);
            index_t cell_size = GEO::geo_max(desc.nb_vertices, desc.nb_facets);
            for(index_t i=0; i<cell_size; ++i) {
                cell_vertices_.push_back(
                    i<desc.nb_vertices ? vertices[i] : index_t(-1)
                );
                cell_adjacents_.push_back(
                    (adjacents != nil && i<desc.nb_facets) ? adjacents[i] : -1
                );
            }
            cell_ptr_.push_back(*(cell_ptr_.rbegin()) + cell_size);
            geo_debug_assert(*cell_ptr_.rbegin() == cell_vertices_.size());
            geo_debug_assert(cell_vertices_.size() == cell_adjacents_.size());
        }
        if(has_attribute(GEO::MESH_CELL_REGION)) {
            cell_regions_.push_back(region);
        }
        ++nb_cells_;
        return index_t(nb_cells()-1);
    }

    index_t MeshBase::add_tet(
        index_t v1, index_t v2, index_t v3, index_t v4,
        signed_index_t adj1, signed_index_t adj2,
        signed_index_t adj3, signed_index_t adj4,
        signed_index_t region
    ) {
        index_t vertices[4];
        vertices[0]=v1;
        vertices[1]=v2;
        vertices[2]=v3;
        vertices[3]=v4;
        signed_index_t adjacent[4];
        adjacent[0]=adj1;
        adjacent[1]=adj2;
        adjacent[2]=adj3;
        adjacent[3]=adj4;
        return add_cell(GEO::MESH_TET, vertices, adjacent, region);
    }

    index_t MeshBase::add_hex(
        index_t v1, index_t v2, index_t v3, index_t v4,
        index_t v5, index_t v6, index_t v7, index_t v8,             
        signed_index_t adj1, signed_index_t adj2, signed_index_t adj3,
        signed_index_t adj4, signed_index_t adj5, signed_index_t adj6,
        signed_index_t region
    ) {
        index_t vertices[8];
        vertices[0]=v1;
        vertices[1]=v2;
        vertices[2]=v3;
        vertices[3]=v4;
        vertices[4]=v5;
        vertices[5]=v6;
        vertices[6]=v7;
        vertices[7]=v8;
        signed_index_t adjacent[6];
        adjacent[0]=adj1;
        adjacent[1]=adj2;
        adjacent[2]=adj3;
        adjacent[3]=adj4;
        adjacent[2]=adj5;
        adjacent[3]=adj6;
        return add_cell(GEO::MESH_HEX, vertices, adjacent, region);
    }

    index_t MeshBase::add_prism(
        index_t v1, index_t v2, index_t v3, 
        index_t v4, index_t v5, index_t v6, 
        signed_index_t adj1, signed_index_t adj2, signed_index_t adj3,
        signed_index_t adj4, signed_index_t adj5,
        signed_index_t region
    ) {
        index_t vertices[6];
        vertices[0]=v1;
        vertices[1]=v2;
        vertices[2]=v3;
        vertices[3]=v4;
        vertices[4]=v5;
        vertices[5]=v6;
        signed_index_t adjacent[5];
        adjacent[0]=adj1;
        adjacent[1]=adj2;
        adjacent[2]=adj3;
        adjacent[3]=adj4;
        adjacent[2]=adj5;
        return add_cell(GEO::MESH_PRISM, vertices, adjacent, region);
    }

    index_t MeshBase::add_pyramid(
        index_t v1, index_t v2, index_t v3, 
        index_t v4, index_t v5, 
        signed_index_t adj1, signed_index_t adj2, signed_index_t adj3,
        signed_index_t adj4, signed_index_t adj5,
        signed_index_t region
    ) {
        index_t vertices[5];
        vertices[0]=v1;
        vertices[1]=v2;
        vertices[2]=v3;
        vertices[3]=v4;
        vertices[4]=v5;
        signed_index_t adjacent[5];
        adjacent[0]=adj1;
        adjacent[1]=adj2;
        adjacent[2]=adj3;
        adjacent[3]=adj4;
        adjacent[2]=adj5;
        return add_cell(GEO::MESH_PYRAMID, vertices, adjacent, region);
    }

    index_t MeshBase::add_connector(
        index_t v1, index_t v2, index_t v3, index_t v4, 
        signed_index_t adj1, signed_index_t adj2, signed_index_t adj3,
        signed_index_t region
    ) {
        index_t vertices[4];
        vertices[0]=v1;
        vertices[1]=v2;
        vertices[2]=v3;
        vertices[3]=v4;
        signed_index_t adjacent[3];
        adjacent[0]=adj1;
        adjacent[1]=adj2;
        adjacent[2]=adj3;
        return add_cell(GEO::MESH_CONNECTOR, vertices, adjacent, region);
    }


    void MeshBase::connect_facets() {
            
        static const index_t NO_CORNER = index_t(-1);
        static const index_t NO_FACET  = index_t(-1);
        
        // Chains the corners around each vertex.
        GEO::vector<index_t>
            next_corner_around_vertex(nb_corners(), NO_CORNER);
        
        // Gives for each vertex a corner incident to it.
        GEO::vector<index_t> v2c(nb_vertices(), NO_CORNER);
        
        // Gives for each corner the facet incident to it
        // (or use c/3 if the surface is triangulated).
        GEO::vector<index_t> c2f;
        if(!is_triangulated()) {
            c2f.assign(nb_corners(), NO_FACET);
        }
        
        // Step 1: chain corners around vertices and compute v2c
        for(index_t f = 0; f < nb_facets(); f++) {
            for(index_t c = facet_begin(f); c < facet_end(f); c++) {
                index_t v = corner_vertex_index(c);
                next_corner_around_vertex[c] = v2c[v];
                v2c[v] = c;
            }
        }
        
        // compute c2f is needed
        if(!is_triangulated()) {
            for(index_t f = 0; f < nb_facets(); f++) {
                for(index_t c = facet_begin(f); c < facet_end(f); c++) {
                    c2f[c] = f;
                }
            }
        }
        
        // Step 2: connect
        for(index_t f1 = 0; f1 < nb_facets(); f1++) {
            for(index_t c1 = facet_begin(f1); c1 < facet_end(f1); c1++) {
                if(corner_adjacent_facet(c1) == signed_index_t(NO_FACET)) {
                    index_t v2 = corner_vertex_index(
                        next_around_facet(f1, c1)
                    );
                    //   Traverse all the corners c2 incident to v1, and
                    // find among them the one that is opposite to c1
                    for(
                        index_t c2 = next_corner_around_vertex[c1];
                        c2 != NO_CORNER; c2 = next_corner_around_vertex[c2]
                    ) {
                        if(c2 != c1) {
                            index_t f2 = is_triangulated() ? c2/3 : c2f[c2];
                            index_t c3 = prev_around_facet(f2, c2);
                            index_t v3 = corner_vertex_index(c3);
                            if(v3 == v2) {
                                set_corner_adjacent_facet(c1, f2);
                                set_corner_adjacent_facet(c3, f1);
                                break;
                            }
                        }
                    }
                }
            }
        }
    }

    void MeshBase::connect_tets() {
        geo_assert(is_tetrahedralized());
        if(nb_tets() == 0) {
            return;
        }

        cell_adjacents_.assign(nb_tets() * 4, -1);
        
        GEO::vector<signed_index_t> next_tet_around_vertex(
            nb_tets() * 4, -1
        );
        GEO::vector<signed_index_t> v2t(nb_vertices(), -1);
        
        // Step 1: chain tets around vertices and compute v2t
        for(index_t t = 0; t < nb_tets(); ++t) {
            for(index_t lv = 0; lv < 4; ++lv) {
                index_t v = tet_vertex_index(t, lv);
                next_tet_around_vertex[4 * t + lv] = v2t[v];
                v2t[v] = signed_index_t(t);
            }
        }
        
        // Step 2: connect tets
        for(index_t t1 = 0; t1 < nb_tets(); ++t1) {
            for(index_t lf1 = 0; lf1 < 4; ++lf1) {
                if(tet_adjacent(t1, lf1) == -1) {
                    index_t v1 = tet_facet_vertex_index(t1, lf1, 0);
                    index_t v2 = tet_facet_vertex_index(t1, lf1, 1);
                    index_t v3 = tet_facet_vertex_index(t1, lf1, 2);
                    for(
                        signed_index_t t2 = v2t[v1]; t2 != -1;
                        t2 = next_tet_around_vertex[
                            4 * t2 + find_tet_vertex(index_t(t2), v1)
                        ]
                    ) {
                        signed_index_t lf2 =
                            find_tet_facet(index_t(t2), v3, v2, v1);
                        if(lf2 != -1) {
                            set_tet_adjacent(t1, lf1, t2);
                            set_tet_adjacent(index_t(t2), index_t(lf2), t1);
                            break;
                        }
                    }
                }
            }
        }
    }
    

    /**
     * \brief Tests whether two indices triplets match
     *  up to a circular permutation.
     * \param[in] v1 index of the first vertex of the first triangle
     * \param[in] v2 index of the second vertex of the first triangle
     * \param[in] v3 index of the third vertex of the first triangle
     * \param[in] w1 index of the first vertex of the second triangle
     * \param[in] w2 index of the second vertex of the second triangle
     * \param[in] w3 index of the third vertex of the second triangle
     * \retval true if (\p v1, \p v2, \p v3) = (\p w1, \p w2, \p w3)
     *  up to a circular permutation
     * \retval false otherwise
     */
    inline bool triangles_equal(
        index_t v1, index_t v2, index_t v3,
        index_t w1, index_t w2, index_t w3
    ) {
        return (
            (v1 == w1 && v2 == w2 && v3 == w3) ||
            (v1 == w2 && v2 == w3 && v3 == w1) ||
            (v1 == w3 && v2 == w1 && v3 == w2) 
        );
    }
    
    bool MeshBase::triangular_facet_matches_quad_facet(
        index_t c1, index_t lf1,
        index_t c2, index_t lf2
    ) const {
        geo_debug_assert(cell_facet_nb_vertices(c1,lf1) == 3);
        geo_debug_assert(cell_facet_nb_vertices(c2,lf2) == 4);
        
        index_t v1 = cell_facet_vertex_index(c1,lf1,0);
        index_t v2 = cell_facet_vertex_index(c1,lf1,1);
        index_t v3 = cell_facet_vertex_index(c1,lf1,2);
        index_t w1 = cell_facet_vertex_index(c2,lf2,0);
        index_t w2 = cell_facet_vertex_index(c2,lf2,1);
        index_t w3 = cell_facet_vertex_index(c2,lf2,2);
        index_t w4 = cell_facet_vertex_index(c2,lf2,3);

        // Note: subtriangles in (w1,w2,w3,w4) are
        // in reverse order since two facets can be
        // connected only if they have opposite
        // orientations.
        return (
            triangles_equal(v1,v2,v3,w4,w3,w2) ||
            triangles_equal(v1,v2,v3,w3,w2,w1) ||
            triangles_equal(v1,v2,v3,w2,w1,w4) ||            
            triangles_equal(v1,v2,v3,w1,w4,w3)
        ) ;
    }

    bool MeshBase::triangular_facets_have_common_edge(
        index_t c1, index_t f1,
        index_t c2, index_t f2,
        index_t& e1, index_t& e2
    ) const {
        geo_debug_assert(cell_facet_nb_vertices(c1,f1) == 3);
        geo_debug_assert(cell_facet_nb_vertices(c2,f2) == 3);
        e1 = index_t(-1);
        e2 = index_t(-1);
        for(e1=0; e1<3; ++e1) {
            for(e2=0; e2<3; ++e2) {
                if(
                    cell_facet_vertex_index(c1, f1, (e1+1)%3) ==
                    cell_facet_vertex_index(c2, f2, (e2+2)%3)  &&
                    cell_facet_vertex_index(c1, f1, (e1+2)%3) ==
                    cell_facet_vertex_index(c2, f2, (e2+1)%3) 
                ) {
                    return true;
                }
            }
        }
        return false;
    }

    /**
     * \brief Creates a connector between a quadrandular facet and two 
     *  triangular facets.
     * \details This function is used by the implementation of 
     *  MeshBase::connect_cells()
     * \param[in] M a reference to the MeshBase
     * \param[in] c1 index of the cell that has the quadrangular facet
     * \param[in] lf1 index of the quadrangular facet in \p c1
     * \param[in] matches a const reference to a vector of 
     *  (cell index, facet index) pairs that are candidate triangles to be 
     *  connected to the quadrangular facet. Each of them
     *  has three vertices in common with the quadrangular facet. 
     *  It may contain more than two (cell,facet) index pairs. 
     *  In this case, among them we select the pair of triangular facets
     *  that have an edge in common. 
     * \retval true if a connector was created. A connector is created if 
     *  among the candidate triangular facets there are exactly two facets 
     *  on the border with an edge in common.
     * \retval false otherwise
     */
    static bool create_connector(
        MeshBase& M, index_t c1, index_t lf1,
        const std::vector< std::pair<index_t, index_t> >& matches
    ) {
        if(matches.size() == 0) {
            return false;
        }
        
        if(matches.size() == 1) {
            GEO::Logger::warn("Mesh")
                << "Found only one triangular facet adjacent to a quad facet"
                << std::endl;
            return false;
        }

        // Find among the matches two facets that have an edge in common
        // Yes, there can be more than two candidate facets with three
        //  vertices in common with the quad facet ! But among them,
        //  only two of them have an edge in common.
        index_t adj_c1 = index_t(-1);
        index_t adj_lf1 = index_t(-1);
        index_t adj_c2 = index_t(-1);
        index_t adj_lf2 = index_t(-1);
        index_t e1 = 4;
        index_t e2 = 4;

        index_t nb_found=0;
        for(index_t i=0; i<index_t(matches.size()); ++i) {
            for(index_t j=i+1; j<index_t(matches.size()); ++j) {
                if(M.triangular_facets_have_common_edge(
                       matches[i].first, matches[i].second,
                       matches[j].first, matches[j].second,
                       e1, e2
                )) {
                    adj_c1 = matches[i].first;
                    adj_lf1 = matches[i].second;
                    adj_c2 = matches[j].first;
                    adj_lf2 = matches[j].second;
                    ++nb_found;
                }
            }
        }

        // Sanity check: make sure that we only found a single pair
        // of triangular facets with a common edge that matches the quad.
        if(nb_found > 1) {
            GEO::Logger::warn("Mesh")
                << "Found more than two triangular facets adjacent to a quad"
                << std::endl;
            return false;
        }

        if(nb_found == 0) {
            GEO::Logger::warn("Mesh")
                << "Triangular facets adjacent to a quad have no common edge"
                << std::endl;
            return false;
        }

        // Sanity check: make sure the triangular facets
        // are on the border.
        if(
            M.cell_adjacent(adj_c1, adj_lf1) != -1 ||
            M.cell_adjacent(adj_c2, adj_lf2) != -1
        ) {
            GEO::Logger::warn("Mesh")
                << "Matching tet facets are not on border"
                << std::endl;
            return false;
        }

        // v1 and v2 are on the common edge
        index_t v1 = M.cell_facet_vertex_index(
            adj_c1, adj_lf1, (e1+1)%3
        );
        
        index_t v2 = M.cell_facet_vertex_index(
            adj_c1, adj_lf1, (e1+2)%3
        );
                    
        // w1 and w2 are the opposite vertices
        index_t w1 = M.cell_facet_vertex_index(adj_c1, adj_lf1, e1);
        index_t w2 = M.cell_facet_vertex_index(adj_c2, adj_lf2, e2);
        
        // Create the connector
        index_t conn = M.add_connector(
            v1, w2, v2, w1,
            signed_index_t(c1),
            signed_index_t(adj_c1),
            signed_index_t(adj_c2)
        );

        // Connect the cells with the connector
        M.set_cell_adjacent(c1, lf1, signed_index_t(conn));
        M.set_cell_adjacent(adj_c1, adj_lf1, signed_index_t(conn));
        M.set_cell_adjacent(adj_c2, adj_lf2, signed_index_t(conn));

        return true;
    }
    
    void MeshBase::connect_cells() {
        if(is_tetrahedralized()) {
            connect_tets();
            return;
        }
        cell_adjacents_.assign(cell_adjacents_.size(), -1);
        GEO::vector<signed_index_t> next_cell_around_vertex(
            cell_vertices_.size(), -1
        );
        GEO::vector<signed_index_t> v2cell(nb_vertices(), -1);
            
        // Step 1: chain cells around vertices and compute v2cell
        for(index_t c = 0; c < nb_cells(); ++c) {
            for(index_t lv = 0; lv < cell_nb_vertices(c); ++lv) {
                index_t v = cell_vertex_index(c, lv);
                next_cell_around_vertex[cell_vertices_begin(c) + lv] =
                    v2cell[v];
                v2cell[v] = signed_index_t(c);
            }
        }

        // Step 2: connect cells
        // (c1,lf1) traverse all the cell facets
        for(index_t c1 = 0; c1 < nb_cells(); ++c1) {
            for(index_t lf1 = 0; lf1 < cell_nb_facets(c1); ++lf1) {

                // If (c1,lf1) is on the border, try to connect it
                if(cell_adjacent(c1, lf1) == -1) {

                    // v1 is one of the vertices of (c1,lf1)
                    index_t v1 = cell_facet_vertex_index(c1,lf1,0);

                    // c2 traverses all the cells incident to v1
                    for(
                        index_t c2 = index_t(v2cell[v1]); c2 != index_t(-1);
                        c2 = index_t(next_cell_around_vertex[
                                cell_vertices_begin(c2) +
                                index_t(find_cell_vertex(c2,v1))
                        ])
                    ) {

                        // If we find a cell facet lf2 compatible with (c1,lf1)
                        // in c2, then connect (c1,lf1) to c2 and
                        // (c2,lf2) to c1.
                        signed_index_t lf2 =
                            find_cell_facet(index_t(c2), c1, lf1);
                        if(lf2 != -1) {
                            set_cell_adjacent(c1, lf1, signed_index_t(c2));
                            set_cell_adjacent(
                                index_t(c2), index_t(lf2), signed_index_t(c1)
                            );
                            break;
                        }
                    }
                }
            }
        }
        
        // Step 3: Create connectors, i.e. artificial cells that represent
        // non-conformal connections between two triangular facets and
        // a quadrangular facet.

        // Backup nb_cells since we are creating new cells (connectors)
        // during this loop.
        index_t nb_cells0 = nb_cells();

        index_t weird=0;

        // For each quadrangular face, we compute the list of candidate
        // triangular faces to be connected with it (a vector of
        // (cell index, facet index) pairs).
        std::vector< std::pair<index_t, index_t> > matches;
                
        // (c1,f1) traverse all quadrangular cell facets on the border
        for(index_t c1=0; c1 < nb_cells0; ++c1) {
            if(cell_type(c1) == GEO::MESH_TET) {
                continue;
            }
            for(index_t lf1=0; lf1<cell_nb_facets(c1); ++lf1) {
                if(
                    cell_facet_nb_vertices(c1,lf1) != 4 ||
                    cell_adjacent(c1,lf1) != -1
                ) {
                    continue;
                }

                // Now c2 traverses all the cells incident to one of
                // the vertex of (c1,lf1)
                matches.resize(0);
                for(index_t lv1=0; lv1<cell_facet_nb_vertices(c1,lf1); ++lv1) {
                    index_t v1 = cell_facet_vertex_index(c1,lf1,lv1);
                    for(
                        index_t c2 = index_t(v2cell[v1]); c2 != index_t(-1);
                        c2 = index_t(next_cell_around_vertex[
                                cell_vertices_begin(c2) +
                                index_t(find_cell_vertex(c2,v1))
                        ])
                    ) {
                        geo_debug_assert(find_cell_vertex(c2,v1) != -1);
                        if(c2 == c1 || cell_type(c2) == GEO::MESH_HEX) {
                            continue;
                        }

                        // Among all the triangular facets of c2, find the ones
                        // that can be connected to (c1,lf1)
                        for(index_t lf2=0; lf2<cell_nb_facets(c2); ++lf2) {
                            if(cell_facet_nb_vertices(c2,lf2) != 3) {
                                continue;
                            }
                            if(triangular_facet_matches_quad_facet(
                                   c2,lf2,c1,lf1
                            )) {
                                matches.push_back(std::make_pair(c2,lf2));
                            } 
                        }
                    }
                }

                // Make sure we get each match once only
                std::sort(matches.begin(), matches.end());
                matches.erase(
                    std::unique(matches.begin(),matches.end()),
                    matches.end()
                );

                if(
                    matches.size() != 0 &&
                    !create_connector(*this,c1,lf1,matches)
                ) {
                    if(debug_coords_float_ != nil) {
                        std::string filename =
                            "debug_cell_" +
                            GEO::String::to_string(weird) + ".obj";
                        std::ofstream out(filename.c_str());
                        index_t index_base = 1;
                        GEO::save_cell(
                            *this, c1, index_base, debug_coords_float_, out
                        );
                        for(index_t i=0; i<matches.size(); ++i) {
                            GEO::save_cell(
                                *this, matches[i].first,
                                index_base, debug_coords_float_, out
                            );
                        }
                    }
                    ++weird;
                }
            }
        }
        if(weird != 0) {
            GEO::Logger::warn("Mesh") << "Encountered "
                                      << weird
                                      << " invalid connector configurations"
                                      << std::endl;
        }
    }
    
    void MeshBase::compute_tets_boundaries() {
        geo_debug_assert(is_tetrahedralized());
        triangulated_ = true;
        facet_ptr_.clear();
        corner_vertices_.clear();
        for(index_t t = 0; t < nb_tets(); ++t) {
            for(index_t lf = 0; lf < 4; ++lf) {
                if(tet_adjacent(t, lf) == -1) {
                    index_t v1 = local_tet_facet_vertex_index(lf, 0);
                    index_t v2 = local_tet_facet_vertex_index(lf, 1);
                    index_t v3 = local_tet_facet_vertex_index(lf, 2);
                    corner_vertices_.push_back(tet_vertex_index(t, v3));
                    corner_vertices_.push_back(tet_vertex_index(t, v2));
                    corner_vertices_.push_back(tet_vertex_index(t, v1));
                }
            }
        }
        corner_adjacent_facets_.assign(corner_vertices_.size(), -1);
        nb_facets_ = corner_vertices_.size() / 3;
        connect_facets();
    }

    void MeshBase::compute_cells_boundaries() {
        if(is_tetrahedralized()) {
            compute_tets_boundaries();
            return;
        }
        triangulated_ = false;
        facet_ptr_.clear();
        corner_vertices_.clear();
        facet_ptr_.push_back(0);
        for(index_t c=0; c<nb_cells(); ++c) {
            for(index_t lf=0; lf<cell_nb_facets(c); ++lf) {
                if(cell_adjacent(c,lf) == -1) {
                    for(index_t lv=0; lv<cell_facet_nb_vertices(c,lf); ++lv) {
                        corner_vertices_.push_back(
                            cell_facet_vertex_index(c,lf,lv)
                        );
                    }
                    facet_ptr_.push_back(corner_vertices_.size());
                }
            }
        }
        nb_facets_ = facet_ptr_.size()-1;
        corner_adjacent_facets_.assign(corner_vertices_.size(),-1);
        connect_facets();
    }

    void MeshBase::show_stats(const std::string& tag) const {
        index_t nb_borders = 0;
        for(index_t i = 0; i < nb_corners(); i++) {
            if(corner_adjacent_facet(i) < 0) {
                nb_borders++;
            }
        }
        GEO::Logger::out(tag)
            << "nb_v:" << nb_vertices()
            << " nb_f:" << nb_facets()
            << " nb_b:" << nb_borders
            << " tri:" << is_triangulated()
            << " dim:" << index_t(dimension())
            << std::endl;

        // TODO: display has_weights flag (but it is
        // only available in derived class)
        
        if(nb_cells() != 0) {
            if(is_tetrahedralized()) {
                GEO::Logger::out(tag) << " nb_tets:"
                                      << nb_tets() << std::endl;
            } else {
                
                index_t nb_cells_by_type[GEO::MESH_NB_CELL_TYPES];
                for(index_t i=0; i<GEO::MESH_NB_CELL_TYPES; ++i) {
                    nb_cells_by_type[i] = 0;
                }
                
                for(index_t c=0; c<nb_cells(); ++c) {
                    geo_debug_assert(cell_type(c) < GEO::MESH_NB_CELL_TYPES);
                    ++nb_cells_by_type[cell_type(c)];
                }
                
                GEO::Logger::out(tag) << " Hybrid - nb_cells:"
                                      << nb_cells() << " "
                                      << " Tet:" << nb_cells_by_type[0]
                                      << " Hex:" << nb_cells_by_type[1]
                                      << " Psm:" << nb_cells_by_type[2]
                                      << " Pmd:" << nb_cells_by_type[3]
                                      << " Cnx:" << nb_cells_by_type[4]
                                      << std::endl;
            }
        }
    }
    
}

