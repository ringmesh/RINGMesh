/*
 * Copyright (c) 2012-2017, Association Scientifique pour la Geologie et ses
 * Applications (ASGA). All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *     * Redistributions in binary form must reproduce the above copyright
 *       notice, this list of conditions and the following disclaimer in the
 *       documentation and/or other materials provided with the distribution.
 *     * Neither the name of ASGA nor the
 *       names of its contributors may be used to endorse or promote products
 *       derived from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
 * THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 * PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL ASGA BE LIABLE FOR ANY DIRECT,
 * INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF
 * THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 *     http://www.ring-team.org
 *
 *     RING Project
 *     Ecole Nationale Superieure de Geologie - GeoRessources
 *     2 Rue du Doyen Marcel Roubault - TSA 70605
 *     54518 VANDOEUVRE-LES-NANCY
 *     FRANCE
 */

#pragma once

#include <ringmesh/basic/common.h>

#include <memory>

#include <geogram/mesh/mesh.h>
#include <geogram/mesh/mesh_io.h>

#include <ringmesh/geogram_extension/geogram_extension.h>

#include <ringmesh/mesh/mesh.h>
#include <ringmesh/mesh/mesh_index.h>

namespace RINGMesh
{
    FORWARD_DECLARATION_DIMENSION_CLASS( CartesianGridBuilder );
}

/*!
 * @file Generic class for creating cartesian grids
 * @author Melchior Schuh-Senlis
 */

namespace RINGMesh
{
    // This is an array template of doubles
    template < index_t DIMENSION >
    using intvecn = GEO::vecng< DIMENSION, signed_index_t >;
    // This is an array of 3 integers
    using intvec3 = intvecn< 3 >;
    // This is an array of 2 integers
    using intvec2 = intvecn< 2 >;

    template < index_t DIMENSION >
    class RINGMESH_API CartesianGrid
    {
        friend class CartesianGridBuilder< DIMENSION >;

    public:
        CartesianGrid( intvecn< DIMENSION > nb_cells,
            index_t nb_total_cells,
            ReferenceFrame< DIMENSION > vec_cartesian_axis )
            : nb_cells_( nb_cells ),
              nb_total_cells_( nb_total_cells ),
              vec_cartesian_axis_( vec_cartesian_axis )
        {
            //		vec_reference_axis_ = ;
        }

        //		CartesianGrid( index_t nb_cell_i, index_t nb_cell_j, index_t
        // nb_cell_k,
        //				vec3 vec_axis_i, vec3 vec_axis_j, vec3 vec_axis_k,
        //				vec3 origin_position) : nb_cell_i_(nb_cell_i),
        // nb_cell_j_(nb_cell_j), nb_cell_k_(nb_cell_k),
        //				nb_cells_(nb_cell_i*nb_cell_j*nb_cell_k),
        // vec_axis_i_(vec_axis_i), vec_axis_j_(vec_axis_j),
        //				vec_axis_k_(vec_axis_k),
        // origin_position_(origin_position)
        //		{
        //			vec_x_base_ijk_[0] = ;
        //		}

        static MeshType type_name_static()
        {
            return 'CartesianGrid';
        }
        void save_mesh( const std::string& filename ) const override
        {
            // TODO
        }
        void print_mesh_bounded_attributes() const override
        {
            // TODO
        }
        //		GEO::AttributesManager& vertex_attribute_manager() const
        // override
        //		{
        //			return mesh_->vertices.attributes();
        //		}
        MeshType type_name() const override
        {
            return type_name_static();
        }
        //		static std::string default_extension_static()
        //		{
        //			return "geogram";
        //		}
        //		std::string default_extension() const override
        //		{
        //			return default_extension_static();
        //		}
        //		index_t nb_vertices() const override
        //		{
        //			return mesh_->vertices.nb();
        //		}

    private:
        /**
         * \brief Forbids copy.
         * \details This is to make sure that client code does
         *   not unintentionally copy a Mesh (for
         *   instance by passing it by-value to a function).
         *   Use copy() instead.
         */
        CartesianGrid( const CartesianGrid& rhs );

        /**
         * \brief Forbids copy.
         * \details This is to make sure that client code does
         *   not unintentionally copies a Mesh (for
         *   instance by passing it by-value to a function).
         *   Use copy() instead.
         */
        const CartesianGrid& operator=( const CartesianGrid& rhs );

    protected:
        //		index_t nb_cell_i_;
        //		index_t nb_cell_j_;
        //		index_t nb_cell_k_;
        intvecn< DIMENSION > nb_cells_;
        index_t nb_total_cells_;

        //		vec3 vec_axis_i_;
        //		vec3 vec_axis_j_;
        //		vec3 vec_axis_k_;
        ReferenceFrame< DIMENSION > vec_cartesian_axis_;

        //		vec3 vec_x_base_ijk_;
        //		vec3 vec_y_base_ijk_;
        //		vec3 vec_z_base_ijk_;
        ReferenceFrame< DIMENSION > vec_reference_axis_;

    public:
        vecn< DIMENSION > cell_center_reference_coords(
            intvecn< DIMENSION > cartesian_coords ) const
        {
            vecn< DIMENSION > coords;
            for( index_t i = 0; i < DIMENSION; i++ )
            {
            }
            //			coords[0] =
            // origin_position_[0]+i*vec_axis_i_[0]+j*vec_axis_j_[0]+k*vec_axis_k_[0];
            //			coords[1] =
            // origin_position_[1]+i*vec_axis_i_[1]+j*vec_axis_j_[1]+k*vec_axis_k_[1];
            //			coords[2] =
            // origin_position_[2]+i*vec_axis_i_[2]+j*vec_axis_j_[2]+k*vec_axis_k_[2];
            return coords;
        }

        intvec3 containing_cell( vec3 vertex ) const
        {
            intvec3 coord;
            coord[0] = ( vertex[0] - origin_position_[0] ) / ;
            return coord;
        }

        index_t cell_offset( intvec3 coords ) const
        {
            return coords[0]
                   + nb_cells_i_ * ( coords[1] + nb_cells_j_ * coords[2] );
        }

        intvec3 ijk_from_offset( index_t offset ) const
        {
            intvec3 coords;
            coords[2] = offset / ( nb_cell_i_ * nb_cell_j_ );
            coords[1] = ( offset - coords[2] * nb_cells_i_ * nb_cells_j_ )
                        / nb_cells_i_;
            coords[0] = offset - coords[2] * nb_cells_i_ * nb_cells_j_
                        - coords[1] * nb_cells_i_;
            return coords;
        }
    };
}
