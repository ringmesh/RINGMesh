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

#include <ringmesh/basic/frame.h>
#include <ringmesh/mesh/mesh.h>

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
        ringmesh_disable_copy_and_move( CartesianGrid );
        friend class CartesianGridBuilder< DIMENSION >;

    public:
        CartesianGrid( intvecn< DIMENSION > nb_cells_in_each_direction,
            ReferenceFrame< DIMENSION > vec_cartesian_axis )
            : nb_cells_in_each_direction_(
                  std::move( nb_cells_in_each_direction ) ),
              cartesian_frame_( std::move( vec_cartesian_axis ) )
        {
            nb_total_cells_ = 1;
            for( auto i : range( DIMENSION ) )
            {
                nb_total_cells_ *= nb_cells_in_each_direction[i];
            }
            inverse_cartesian_frame_ =
                FrameManipulator::reference_frame_from_global_to_local(
                    cartesian_frame_ );
        }

        static MeshType type_name_static()
        {
            return "CartesianGrid";
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

        vecn< DIMENSION >& cell_center_global_coords(
            const intvecn< DIMENSION >& cartesian_coords ) const
        {
            vecn< DIMENSION > cartesian_double_coords;
            for( auto i : range( DIMENSION ) )
            {
                cartesian_double_coords[i] =
                    static_cast< double >( cartesian_coords[i] );
            }
            return FrameManipulator::coords_from_frame_to_global(
                cartesian_frame_, cartesian_double_coords );
        }

        intvecn< DIMENSION >& containing_cell_from_global_vertex(
            const vecn< DIMENSION >& reference_vertex ) const
        {
            vecn< DIMENSION > frame_vertex =
                FrameManipulator::coords_from_frame_to_global(
                    inverse_cartesian_frame_, reference_vertex );
            return this->containing_cell_from_local_vertex( frame_vertex );
        }

        intvecn< DIMENSION > containing_cell_from_local_vertex(
            const vecn< DIMENSION >& vertex ) const
        {
            intvecn< DIMENSION > coord;
            for( auto i : range( DIMENSION ) )
            {
                coord[i] = std::floor( vertex[i] + 0.5 );
            }
            return coord;
        }

        index_t cell_offset( const intvecn< DIMENSION >& coords ) const
        {
            index_t offset{ 0 };
            index_t mult{ 1 };
            for( auto i : range( DIMENSION ) )
            {
                offset += coords[i] * mult;
                mult *= nb_cells_in_each_direction_[i];
            }
            return offset;
        }

        intvecn< DIMENSION > ijk_from_offset( const index_t offset ) const
        {
            intvecn< DIMENSION > coords;
            index_t off{ 0 };
            index_t div{ nb_total_cells_
                         / nb_cells_in_each_direction_[DIMENSION - 1] };
            for( auto i : range( DIMENSION ) )
            {
                coords[DIMENSION - 1 - i] = ( offset - off ) / div;
                off += coords[DIMENSION - 1 - i] * div;
                div /= nb_cells_in_each_direction_[DIMENSION - 1 - i];
            }
            return coords;
        }

    protected:
        intvecn< DIMENSION > nb_cells_in_each_direction_;
        index_t nb_total_cells_;

        ReferenceFrame< DIMENSION > cartesian_frame_;
        ReferenceFrame< DIMENSION > inverse_cartesian_frame_;
    };
    ALIAS_2D_AND_3D( CartesianGrid );

    //    class RINGMESH_API CartesianGridVolumeMesh : public VolumeMesh
    //	{
    //
    //	};
}
