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

//#include <memory>

//#include <geogram/mesh/mesh.h>
//#include <geogram/mesh/mesh_io.h>

//#include <ringmesh/geogram_extension/geogram_extension.h>

#include <ringmesh/basic/frame.h>
//#include <ringmesh/mesh/mesh.h>
//#include <ringmesh/mesh/mesh_index.h>

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
    template < index_t DIMENSION >
    class RINGMESH_API CartesianGrid
    {
        friend class CartesianGridBuilder< DIMENSION >;

    public:
        CartesianGrid( intvecn< DIMENSION > nb_cells,
            index_t nb_total_cells,
            ReferenceFrame< DIMENSION > vec_cartesian_axis )
            : nb_cells_( nb_cells ), cartesian_frame_( vec_cartesian_axis )
        {
            nb_total_cells_ = 1;
            for( auto i : RINGMesh::range( DIMENSION ) )
            {
                nb_total_cells_ *= nb_cells;
            }
        }

        static MeshType type_name_static()
        {
            return 'CartesianGrid';
        }
        void save_mesh( const std::string& filename ) const
        {
            // TODO
        }
        void print_mesh_bounded_attributes() const
        {
            // TODO
        }
        //		GEO::AttributesManager& vertex_attribute_manager() const
        //override
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
        intvecn< DIMENSION > nb_cells_;
        index_t nb_total_cells_;

        ReferenceFrame< DIMENSION > cartesian_frame_;

    public:
        vecn< DIMENSION > cell_center_reference_coords(
            const intvecn< DIMENSION > cartesian_coords ) const
        {
            vecn< DIMENSION > cart_double_coords;
            for( auto i : RINGMesh::range( DIMENSION ) )
            {
                cart_double_coords[i] =
                    static_cast< double >( cartesian_coords[i] );
            }
            return cartesian_frame_.coords_to_base( cart_double_coords );
        }

        intvecn< DIMENSION > containing_cell(
            const vecn< DIMENSION > vertex ) const
        {
            intvecn< DIMENSION > coord;
            for( auto i : RINGMesh::range( DIMENSION ) )
            {
                coord[i] = std::floor( vertex[i] + 0.5 );
            }
            return coord;
        }

        index_t cell_offset( const intvecn< DIMENSION > coords ) const
        {
            index_t offset = 0;
            index_t mult = 1;
            for( auto i : RINGMesh::range( DIMENSION ) )
            {
                offset += coords[i] * mult;
                mult *= nb_cells_[i];
            }
            return offset;
        }

        intvecn< DIMENSION > ijk_from_offset( const index_t offset ) const
        {
            intvecn< DIMENSION > coords;
            index_t off = 0;
            index_t div = nb_total_cells_ / nb_cells_[DIMENSION - 1];
            for( auto i : RINGMesh::range( DIMENSION ) )
            {
                coords[DIMENSION - 1 - i] = ( offset - off ) / div;
                off += coords[DIMENSION - 1 - i] * div;
                div /= nb_cells_[DIMENSION - 1 - i];
            }
            return coords;
        }
    };
}
