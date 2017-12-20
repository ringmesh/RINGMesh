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
#include <ringmesh/basic/geometry.h>

#include <ringmesh/mesh/mesh.h>

#include <geogram/basic/attributes.h>

namespace GEO
{
    class AttributesManager;
} // namespace GEO
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
    /**
     * Template class for Cartesian grids of different dimensions
     * Each value of the grid represents the cell for which the point
     * is the corner with the smallest coordinates.
     */
    template < index_t DIMENSION >
    class CartesianGrid
    {
        ringmesh_disable_copy_and_move( CartesianGrid );
        friend class CartesianGridBuilder< DIMENSION >;

    public:
        /*!
         * Constructor for the cartesian grid
         * \param[in] nb_cells_in_each_direction number of cells in each
         * direction of the grid.
         * \param[in] vec_cartesian_axis used to define
         * the cartesian grid cells : the origin of the frame is the position of
         * the cell (0,0,0) in the cartesian grid, the vectors are the
         * directions and length of the grid cells.
         */
        CartesianGrid( ivecn< DIMENSION > nb_cells_in_each_direction,
            ReferenceFrame< DIMENSION > vec_cartesian_axis )
        {
            check_and_update_number_of_cells( nb_cells_in_each_direction );
            check_and_update_frame( vec_cartesian_axis );
            inverse_cartesian_frame_ = ReferenceFrameManipulator< DIMENSION >::
                orthogonal_reference_frame_from_global_to_local(
                    cartesian_frame_ );
            attributes_manager_.resize( nb_total_cells_ );
        }

        void save_mesh( const std::string& filename ) const
        {
            // TODO
        }

        //        void print_mesh_bounded_attributes( const std::string&
        //        output_location ) const
        //        {
        //        	print_bounded_attributes( attributes_manager_,
        //        output_location );
        //        }

        vecn< DIMENSION > cell_center_global_coords(
            const sivecn< DIMENSION >& cartesian_coords ) const
        {
            vecn< DIMENSION > cartesian_double_coords;
            for( auto i : range( DIMENSION ) )
            {
                cartesian_double_coords[i] =
                    static_cast< double >( cartesian_coords[i] );
            }
            return ReferenceFrameManipulator<
                DIMENSION >::coords_from_frame_to_global( cartesian_frame_,
                cartesian_double_coords );
        }

        sivecn< DIMENSION > containing_cell_from_global_point(
            const vecn< DIMENSION >& reference_vertex ) const
        {
            // Since coords_from_frame_to_global is faster than
            // coords_from_global_to_frame, we use it with the inverse matrix
            return this->containing_cell_from_local_point(
                ReferenceFrameManipulator< DIMENSION >::
                    coords_from_frame_to_global(
                        inverse_cartesian_frame_, reference_vertex ) );
//        	return this->containing_cell_from_local_point(
//        			ReferenceFrameManipulator< DIMENSION >::coords_from_global_to_frame(
//        			cartesian_frame_, reference_vertex ) );

        }

        sivecn< DIMENSION > containing_cell_from_local_point(
            const vecn< DIMENSION >& vertex ) const
        {
            sivecn< DIMENSION > coord;
            for( auto i : range( DIMENSION ) )
            {
                coord[i] = std::floor( vertex[i] );
            }
            return coord;
        }

        index_t cell_offset( const sivecn< DIMENSION >& coords ) const
        {
            index_t offset{ 0 };
            index_t mult{ 1 };
            for( auto i : range( DIMENSION ) )
            {
                if( coords[i] >= 0
                    && coords[i] < nb_cells_in_each_direction_[i] )
                {
                    offset += coords[i] * mult;
                    mult *= nb_cells_in_each_direction_[i];
                }
                else
                {
                    Logger::warn( "Point ", coords,
                        " has indexes outside of the cartesian grid limits." );
                    return -1;
                }
            }
            return offset;
        }

        index_t cell_offset_from_global_point(
            const vecn< DIMENSION > coords ) const
        {
            return cell_offset( containing_cell_from_global_point( coords ) );
        }

        sivecn< DIMENSION > local_from_offset( const index_t offset ) const
        {
            sivecn< DIMENSION > coords;
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

        index_t nb_cells() const
        {
            return nb_total_cells_;
        }

        index_t nb_cells_axis( index_t i ) const
        {
            return nb_cells_in_each_direction_[i];
        }

        ReferenceFrame< DIMENSION > grid_vectors() const
		{
        	return cartesian_frame_;
		}

        ReferenceFrame< DIMENSION > inverse_grid_vectors() const
		{
        	return inverse_cartesian_frame_;
		}

        GEO::AttributesManager& attributes_manager() const
        {
            return attributes_manager_;
        }

        double cell_volume() const
        {
            double cell_volume{ 1. };
            for( auto i : range( DIMENSION ) )
            {
                cell_volume *= cartesian_frame_[i].length2();
            }
            return sqrt( cell_volume );
        }

    private:
        void check_and_update_number_of_cells(
            ivecn< DIMENSION >& nb_cells_in_each_direction )
        {
            nb_total_cells_ = 1;
            for( auto i : range( DIMENSION ) )
            {
                if( nb_cells_in_each_direction[i] < 1 )
                {
                    throw RINGMeshException( "RINGMesh Test",
                        "Error: You are trying to create a Cartesian Grid "
                        "with no cell in direction ",
                        i,
                        ", and Cartesian Grids must have at least one cell "
                        "in each direction." );
                }
                nb_total_cells_ *= nb_cells_in_each_direction[i];
            }
            nb_cells_in_each_direction_ =
                std::move( nb_cells_in_each_direction );
        }

        void check_and_update_frame(
            ReferenceFrame< DIMENSION >& vec_cartesian_axis )
        {
            if( !ReferenceFrameManipulator< DIMENSION >::is_frame_orthogonal(
                    vec_cartesian_axis ) )
            {
                throw RINGMeshException( "RINGMesh Test",
                    "Error: the frame you are giving for the "
                    "Cartesian Grid is not orthogonal. " );
            }
            cartesian_frame_ = std::move( vec_cartesian_axis );
        }

    protected:
        CartesianGrid() = default;

        void change_frame( ReferenceFrame< DIMENSION >& vec_cartesian_axis )
        {
            check_and_update_frame( vec_cartesian_axis );
            Logger::warn( "You are currently changing the frame of "
                          "the Cartesian grid, this will affect where the "
                          "values of the attributes in the grid are stored in "
                          "space." );
            inverse_cartesian_frame_ = ReferenceFrameManipulator< DIMENSION >::
                reference_frame_from_global_to_local( cartesian_frame_ );
        }

        void change_attribute_manager(
            GEO::AttributesManager attributes_manager )
        {
            attributes_manager_ = std::move( attributes_manager );
        }

    protected:
        ivecn< DIMENSION > nb_cells_in_each_direction_;
        index_t nb_total_cells_;

        ReferenceFrame< DIMENSION > cartesian_frame_;
        ReferenceFrame< DIMENSION > inverse_cartesian_frame_;

        GEO::AttributesManager attributes_manager_;
    };
    ALIAS_2D_AND_3D( CartesianGrid );

    template < index_t DIMENSION >
    class CartesianGridBuilder
    {
    public:
        static std::unique_ptr< CartesianGridBuilder< DIMENSION > >
            create_builder( CartesianGrid< DIMENSION >& cartesian_grid );

        explicit CartesianGridBuilder(
            CartesianGrid< DIMENSION >& cartesian_grid )
            : cartesian_grid_( dynamic_cast< CartesianGrid< DIMENSION >& >(
                  cartesian_grid ) )
        {
        }

        void resize_vec_axis( index_t axis_id, double new_size )
        {
            if( new_size == 0 )
            {
                throw RINGMeshException( "RINGMesh Test",
                    "Error : you are trying to resize a cell to length 0 in "
                    "direction ",
                    axis_id );
            }
            else
            {
                cartesian_grid_.cartesian_frame_[axis_id] *=
                    ( new_size
                        / cartesian_grid_.cartesian_frame_[axis_id].length() );
            }
        }

        void change_frame( ReferenceFrame< DIMENSION >& vec_cartesian_axis )
        {
            cartesian_grid_.change_frame( vec_cartesian_axis );
        }

        void change_attribute_manager(
            GEO::AttributesManager attributes_manager )
        {
            cartesian_grid_.change_attribute_manager( attributes_manager );
        }

        void remove_section_from_cartesian_grid(
            index_t axis_id, index_t section_position )
        {
            // TODO
            //        	if( cartesian_grid_.nb_cells_in_each_direction[axis_id]
            //        < 2 )
            //			{
            //				throw RINGMeshException( "RINGMesh Test",
            //					"Error: You are trying to remove a section in
            // direction",
            //					axis_id, ", but it would reduce the number of
            // cells
            // in
            // this directions below 1." );
            //			}
            //        	cartesian_grid_.nb_total_cells_ -=
            //        cartesian_grid_.nb_total_cells_ /
            //        cartesian_grid_.nb_cells_in_each_direction_[ axis_id ];
            //        	cartesian_grid_.nb_cells_in_each_direction_[ axis_id ]
            //        -= 1;
            //        	cartesian_grid_.attributes_manager_.resize(
            //        cartesian_grid_.nb_total_cells_ );
        }

    protected:
        CartesianGrid< DIMENSION >& cartesian_grid_;
    };
    ALIAS_2D_AND_3D( CartesianGridBuilder );

    //    class RINGMESH_API CartesianGridVolumeMesh : public VolumeMesh3D
    //    {
    //    public:
    //        CartesianGridVolumeMesh() = default;
    //
    //        static MeshType type_name_static()
    //        {
    //            return "CartesianGrid";
    //        }
    //
    //        MeshType type_name() const override
    //        {
    //            return type_name_static();
    //        }
    //
    //        index_t nb_cells() const override
    //        {
    //            return cartesian_grid_.nb_cells();
    //        }
    //
    //        //        GEO::AttributesManager& cell_attribute_manager() const
    //        //        override
    //        //        {
    //        //            return cartesian_grid_.attributes_manager();
    //        //        }
    //
    //        CellType cell_type( index_t cell_id ) const override
    //        {
    //            return static_cast< CellType >( 1 );
    //        }
    //
    //        double cell_volume( index_t cell_id ) const override
    //        {
    //            return cartesian_grid_.cell_volume();
    //        }
    //
    //    private:
    //        CartesianGrid3D cartesian_grid_;
    //    };
}
