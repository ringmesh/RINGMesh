/*
 * Copyright (c) 2012-2018, Association Scientifique pour la Geologie et ses
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

#include <type_traits>

#include <ringmesh/basic/common.h>
#include <ringmesh/basic/frame.h>
#include <ringmesh/basic/geometry.h>
#include <ringmesh/basic/logger.h>

#include <ringmesh/mesh/common.h>

#include <geogram/basic/attributes.h>

namespace GEO
{
    class AttributesManager;
} // namespace GEO
namespace RINGMesh
{
    FORWARD_DECLARATION_DIMENSION_CLASS( CartesianGridBaseBuilder );
    FORWARD_DECLARATION_DIMENSION_CLASS( CartesianGridBuilder );
}

/*!
 * @file Generic class for creating cartesian grids
 * @author Melchior Schuh-Senlis
 */

namespace RINGMesh
{
    /**
     * Template base class for Cartesian grids of different dimensions
     * Each value of the grid represents the cell for which the point
     * is the corner with the smallest coordinates.
     * After the grid is created, it can only be modified through its
     * associated builder.
     */
    template < index_t DIMENSION >
    class CartesianGridBase
    {
        ringmesh_disable_copy_and_move( CartesianGridBase );
        friend class CartesianGridBaseBuilder< DIMENSION >;

    public:
        CartesianGridBase() = default;

        /*!
         * Constructor for the cartesian grid
         * \param[in] nb_cells_in_each_direction number of cells in each
         * direction of the grid.
         * \param[in] vec_cartesian_axis used to define
         * the cartesian grid cells : the origin of the frame is the position of
         * the corner (0,0,0) of the cell (0,0,0) in the cartesian grid, the
         * vectors are the directions and size in these directions of the grid
         * cells.
         */
        CartesianGridBase( ivecn< DIMENSION > nb_cells_in_each_direction,
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
            throw RINGMeshException( "SAVE",
                "The save_mesh function is not implemented"
                "for cartesian grids yet." );
        }

        /*!
         * @return the coordinates in the global frame of the given grid point.
         */
        vecn< DIMENSION > cell_corner_global_coords(
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

        /*!
         * @return the grid coordinates of the given point from its global
         * coordinates.
         */
        vecn< DIMENSION > local_coords_from_global_point(
            const vecn< DIMENSION >& reference_vertex ) const
        {
            // Since coords_from_frame_to_global is faster than
            // coords_from_global_to_frame, we use it with the inverse matrix
            return ReferenceFrameManipulator< DIMENSION >::
                coords_from_frame_to_global(
                    inverse_cartesian_frame_, reference_vertex );
        }

        /*!
         * @return the grid cell containing the given point from its global
         * coordinates.
         */
        sivecn< DIMENSION > containing_cell_from_global_point(
            const vecn< DIMENSION >& reference_vertex ) const
        {
            // Since coords_from_frame_to_global is faster than
            // coords_from_global_to_frame, we use it with the inverse matrix
            return this->containing_cell_from_local_point(
                ReferenceFrameManipulator< DIMENSION >::
                    coords_from_frame_to_global(
                        inverse_cartesian_frame_, reference_vertex ) );
        }

        /*!
         * @return the grid cell containing the given point from its coordinates
         * in the cartesian grid.
         */
        sivecn< DIMENSION > containing_cell_from_local_point(
            const vecn< DIMENSION >& vertex ) const
        {
            sivecn< DIMENSION > coord;
            for( auto i : range( DIMENSION ) )
            {
                coord[i] =
                    static_cast< signed_index_t >( std::floor( vertex[i] ) );
            }
            return coord;
        }

        /*!
         * @return the index of the cell if we consider the grid as a single
         * vector of cells in memory.
         */
        index_t cell_offset( const sivecn< DIMENSION >& coords ) const
        {
            index_t offset{ 0 };
            index_t mult{ 1 };
            for( auto i : range( DIMENSION ) )
            {
                if( coords[i] >= 0
                    && coords[i] < static_cast< signed_index_t >(
                                       nb_cells_in_each_direction_[i] ) )
                {
                    offset += static_cast< index_t >( coords[i] ) * mult;
                    mult *= nb_cells_in_each_direction_[i];
                }
                else
                {
                    Logger::warn( "Point ", coords,
                        " has indices outside of the cartesian grid limits." );
                    return NO_ID;
                }
            }
            return offset;
        }

        /*!
         * @return the index of the cell containing the given point, from its
         * global coordinates, in the vector containing all the cells of the
         * grid.
         */
        index_t cell_offset_from_global_point(
            const vecn< DIMENSION > coords ) const
        {
            return cell_offset( containing_cell_from_global_point( coords ) );
        }

        /*!
         * @return the coordinates of the cell from its index in the vector
         * containing all the cells of the grid.
         */
        sivecn< DIMENSION > local_from_offset( const index_t offset ) const
        {
            sivecn< DIMENSION > coords;
            index_t off{ 0 };
            index_t div{ nb_total_cells_ };
            for( auto i : range( DIMENSION ) )
            {
                div /= nb_cells_in_each_direction_[DIMENSION - 1 - i];
                index_t coordi{ ( offset - off ) / div };
                off += coordi * div;
                coords[DIMENSION - 1 - i] =
                    static_cast< signed_index_t >( coordi );
            }
            return coords;
        }

        /*!
         * @return the total number of cells in grid.
         */
        index_t nb_cells() const
        {
            return nb_total_cells_;
        }

        /*!
         * @return the number of cells in the grid in the given direction.
         */
        index_t nb_cells_axis( index_t i ) const
        {
            return nb_cells_in_each_direction_[i];
        }

        /*!
         * @return the vector containing the number of cells of the grid
         * in each direction.
         */
        const ivecn< DIMENSION >& nb_cells_vector() const
        {
            return nb_cells_in_each_direction_;
        }

        /*!
         * @return the reference frame of the grid.
         */
        const ReferenceFrame< DIMENSION >& grid_vectors() const
        {
            return cartesian_frame_;
        }

        /*!
         * @return the size of the cells in the given direction \i
         * (length of the vector of the grid reference frame in direction \i)
         */
        double grid_vector_size( index_t i ) const
        {
            return cartesian_frame_[i].length();
        }

        /*!
         * @return the reference frame containing the grid coordinates of the
         * global frame.
         */
        const ReferenceFrame< DIMENSION >& inverse_grid_vectors() const
        {
            return inverse_cartesian_frame_;
        }

        /*!
         * @return the attribute manager containing the grid point attributes
         * values.
         */
        GEO::AttributesManager& attributes_manager()
        {
            return attributes_manager_;
        }

        /*!
         * @return the attribute manager containing the grid point attributes
         * values.
         */
        const GEO::AttributesManager& attributes_manager() const
        {
            return attributes_manager_;
        }

        /*!
         * Returns a std::string containing the name of the typename T.
         */
        template < typename T >
        static std::string type_to_string()
        {
            if( std::is_same< T, index_t >::value )
            {
                return "index_t";
            }
            else if( std::is_same< T, int >::value )
            {
                return "int";
            }
            else if( std::is_same< T, float >::value )
            {
                return "float";
            }
            else if( std::is_same< T, double >::value )
            {
                return "double";
            }
            else
            {
                throw RINGMeshException( "CartesianGrid",
                		"This type cannot be implemented as an attribute." );
            }
        }

        /*!
         * Creates and adds a new attribute to the Cartesian grid from its
         * name and values on the grid.
         * \param[in] attribute_name the name of the new attribute to add
         * \param[in] values a vector of the values of the attribute on all
         * the Cartesian grid cells.
         */
        template < typename T >
        void add_attribute(
            const std::string& attribute_name, const std::vector< T >& values )
        {
            if( values.size() != nb_total_cells_ )
            {
                throw RINGMeshException( "CartesianGrid",
                    "The attribute you're trying to add doesn't have the same"
                    " number of values than the grid." );
            }
            attributes_manager_.bind_attribute_store( attribute_name,
                GEO::AttributeStore::
                    create_attribute_store_by_element_type_name(
                        type_to_string< T >(), 1 ) );
            GEO::Attribute< T > attribute{ attributes_manager_,
                attribute_name };
            for( auto i : range( nb_total_cells_ ) )
            {
                attribute[i] = values[i];
            }
        }

        /*!
         * Creates a new attribute with the name \attribute_name and
         * the same value \single_value on every grid point.
         */
        template < typename T >
        void add_attribute( const std::string& attribute_name, T single_value = 0 )
        {
            attributes_manager_.bind_attribute_store( attribute_name,
                GEO::AttributeStore::
                    create_attribute_store_by_element_type_name(
                        type_to_string< T >(), 1 ) );
            GEO::Attribute< T > attribute{ attributes_manager_,
                attribute_name };
            attribute.fill( single_value );
        }

        /*!
         * Returns the value of the attribute \attribute_name at the grid
         * position \position.
         */
        template < typename T >
        T& get_attribute_value(
            const std::string& attribute_name, sivecn< DIMENSION > position )
        {
            auto index = cell_offset( position );
            if( index == NO_ID )
            {
                throw RINGMeshException( "CartesianGrid",
                    "Points outside of the grid have no attribute value, give "
                    "a point position inside the grid range." );
            }
            GEO::Attribute< T > attribute{ attributes_manager_,
                attribute_name };
            return attribute[index];
        }

        /*!
         * @return the volume of a grid cell.
         */
        double cell_volume() const
        {
            double cell_volume{ 1. };
            for( auto i : range( DIMENSION ) )
            {
                cell_volume *= cartesian_frame_[i].length2();
            }
            return std::sqrt( cell_volume );
        }

    private:
        /*!
         * Function to call to change the number of cells of the grid.
         */
        void check_and_update_number_of_cells(
            ivecn< DIMENSION >& nb_cells_in_each_direction )
        {
            nb_total_cells_ = 1;
            for( auto i : range( DIMENSION ) )
            {
                if( nb_cells_in_each_direction[i] < 1 )
                {
                    throw RINGMeshException( "CartesianGrid",
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
                throw RINGMeshException( "CartesianGrid",
                    "Error: the frame you are giving for the "
                    "Cartesian Grid is not orthogonal. " );
            }
            cartesian_frame_ = std::move( vec_cartesian_axis );
        }

    protected:
        /*!
         * Function to call to change the reference frame of the grid.
         */
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
    ALIAS_2D_AND_3D( CartesianGridBase );

    template < index_t DIMENSION >
    class CartesianGrid final : public CartesianGridBase< DIMENSION >
    {
        friend class CartesianGridBuilder< DIMENSION >;
    };
    ALIAS_2D_AND_3D( CartesianGrid );

    /*!
     * Implementation of the Cartesian grid in 3D, allowing to include
     * the 6 planes that surround it as a class attribute.
     */
    template <>
    class CartesianGrid< 3 > final : public CartesianGridBase< 3 >
    {
        friend class CartesianGridBuilder< 3 >;

    public:
        CartesianGrid() = default;

        CartesianGrid( ivec3 nb_cells_in_each_direction,
            ReferenceFrame3D vec_cartesian_axis )
            : CartesianGridBase(
                  nb_cells_in_each_direction, vec_cartesian_axis ),
              grid_cage_()
        {
            create_grid_cage();
        }

        const std::vector< Geometry::Plane >& grid_cage() const
        {
            return grid_cage_;
        }

    protected:
        void create_grid_cage()
        {
            vec3 highest_coordinates_point{
                cartesian_frame_.origin()
                + cartesian_frame_[0] * nb_cells_axis( 0 )
                + cartesian_frame_[1] * nb_cells_axis( 1 )
                + cartesian_frame_[2] * nb_cells_axis( 2 )
            };
            grid_cage_.reserve( 6 );
            for( auto i : range( 3 ) )
            {
                grid_cage_.emplace_back(
                    cartesian_frame_[i], cartesian_frame_.origin() );
                grid_cage_.emplace_back(
                    cartesian_frame_[i], highest_coordinates_point );
            }
        }

        /// The 6 planes of the grid cage are ordered in this way :
        /// First the 2 with a normal to the first axis of the grid,
        /// then the 2 with a normal to the second axis of the grid,
        /// and finally the 2 with a normal to the third axis of the grid.
        /// Each time, the first of the 2 planes in question is the one which
        /// contains the origin of the grid, and the second is the one which
        /// contains the point with the highest local coordinates in the grid.
        std::vector< Geometry::Plane > grid_cage_;
    };

    /*!
     * Implementation of the cartesian grid in 2D, allowing to include
     * the 4 segments that surround it as a class attribute.
     */
    template <>
    class CartesianGrid< 2 > final : public CartesianGridBase< 2 >
    {
        friend class CartesianGridBuilder< 2 >;

    public:
        CartesianGrid() = default;

        CartesianGrid( ivec2 nb_cells_in_each_direction,
            ReferenceFrame2D vec_cartesian_axis )
            : CartesianGridBase(
                  nb_cells_in_each_direction, vec_cartesian_axis )
        {
        }

        /*!
         * Returns the 4 segments that represent the grid cage.
         * They are ordered in this way :
         * first the two vectors in the direction of the first grid vector.
         * then the two vectors in the direction of the second grid vector.
         */
        const std::vector< Geometry::Segment2D > grid_cage() const
        {
            std::vector< Geometry::Segment2D > cage;
            cage.reserve( 4 );
            cage.emplace_back( cartesian_frame_.origin(),
                cartesian_frame_.origin()
                    + cartesian_frame_[0] * nb_cells_axis( 0 ) );
            cage.emplace_back( cartesian_frame_.origin() + cartesian_frame_[1],
                cartesian_frame_.origin()
                    + cartesian_frame_[0] * nb_cells_axis( 0 )
                    + cartesian_frame_[1] * nb_cells_axis( 1 ) );
            cage.emplace_back( cartesian_frame_.origin(),
                cartesian_frame_.origin()
                    + cartesian_frame_[1] * nb_cells_axis( 1 ) );
            cage.emplace_back( cartesian_frame_.origin()
                                   + cartesian_frame_[0] * nb_cells_axis( 0 ),
                cartesian_frame_.origin()
                    + cartesian_frame_[0] * nb_cells_axis( 0 )
                    + cartesian_frame_[1] * nb_cells_axis( 1 ) );

            return cage;
        }
    };

    /*!
     * Builder class associated to CartesianGridBase.
     * It can be instantiated to modify a Cartesian grid base.
     */
    template < index_t DIMENSION >
    class CartesianGridBaseBuilder
    {
        ringmesh_disable_copy_and_move( CartesianGridBaseBuilder );

    public:
        virtual ~CartesianGridBaseBuilder() = default;

        explicit CartesianGridBaseBuilder(
            CartesianGridBase< DIMENSION >& cartesian_grid_base )
            : cartesian_grid_base_(
                  dynamic_cast< CartesianGridBase< DIMENSION >& >(
                      cartesian_grid_base ) )
        {
        }

        void initialise_grid( ivecn< DIMENSION > nb_cells_in_each_direction,
            ReferenceFrame< DIMENSION > vec_cartesian_axis )
        {
            cartesian_grid_base_.check_and_update_number_of_cells(
                nb_cells_in_each_direction );
            cartesian_grid_base_.check_and_update_frame( vec_cartesian_axis );
            cartesian_grid_base_.inverse_cartesian_frame_ =
                ReferenceFrameManipulator< DIMENSION >::
                    orthogonal_reference_frame_from_global_to_local(
                        cartesian_grid_base_.cartesian_frame_ );
            cartesian_grid_base_.attributes_manager_.resize(
                cartesian_grid_base_.nb_total_cells_ );
        }

        /*!
         * Changes the length of vector \axis_id of the reference frame of the
         * grid associated to this builder, and sets it to \new_size.
         */
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
                cartesian_grid_base_.cartesian_frame_[axis_id] *=
                    ( new_size
                        / cartesian_grid_base_.cartesian_frame_[axis_id]
                              .length() );
            }
        }

        /*!
         * Replaces the reference frame of the grid associated to this builder
         * by input reference frame.
         */
        void change_frame( ReferenceFrame< DIMENSION >& vec_cartesian_axis )
        {
            cartesian_grid_base_.change_frame( vec_cartesian_axis );
        }

        /*!
         * Replaces the attribute manager of the grid associated to this builder
         * by input attribute manager.
         */
        void change_attribute_manager(
            GEO::AttributesManager attributes_manager )
        {
            cartesian_grid_base_.change_attribute_manager( attributes_manager );
        }

        /*!
         * Removes a section of cells of the grid, normal to the vector \axis_id
         * of its reference frame, and with coordinate \section_position on this
         * axis.
         * \param[in] axis_id axis number between 0 and DIMENSION
         */
        void remove_section_from_cartesian_grid(
            index_t axis_id, index_t section_position )
        {
            if( axis_id > DIMENSION - 1 )
            {
                throw RINGMeshException( "CartesianGrid",
                    "Error: Give an axis_id between 0 and the "
                    "dimension of the grid -1." );
            }
            if( cartesian_grid_base_.nb_cells_axis( axis_id ) < 2 )
            {
                throw RINGMeshException( "CartesianGrid",
                    "Error: You are trying to remove a section in direction",
                    axis_id,
                    ", but it would reduce the number of cells in this "
                    "directions below 1." );
            }
            if( section_position
                > cartesian_grid_base_.nb_cells_axis( axis_id ) )
            {
                throw RINGMeshException( "CartesianGrid",
                    "Error: Give a correct position for the "
                    "section you wish to remove." );
            }

            cartesian_grid_base_.nb_total_cells_ -=
                cartesian_grid_base_.nb_total_cells_
                / cartesian_grid_base_.nb_cells_in_each_direction_[axis_id];
            cartesian_grid_base_.nb_cells_in_each_direction_[axis_id] -= 1;

            GEO::vector< index_t > compression_vec =
                compression_vector( axis_id, section_position );
            cartesian_grid_base_.attributes_manager_.compress(
                compression_vec );
            cartesian_grid_base_.attributes_manager_.resize(
                cartesian_grid_base_.nb_total_cells_ );
        }

    private:
        /*!
         * Returns the compression vector needed to remove the values of the
         * removed section from the attributes.
         */
        GEO::vector< index_t > compression_vector(
            index_t axis_id, index_t section_position )
        {
            index_t vec_size{ cartesian_grid_base_.attributes_manager_.size() };
            GEO::vector< index_t > compression_vec{ vec_size, 0 };
            index_t iterator{ 0 };
            for( auto i : range( vec_size ) )
            {
                if( cartesian_grid_base_.local_from_offset( i )[axis_id]
                    != static_cast< signed_index_t >( section_position ) )
                {
                    compression_vec[i] = iterator;
                    iterator++;
                }
                else
                {
                    compression_vec[i] = index_t( -1 );
                }
            }
            return compression_vec;
        }

    protected:
        CartesianGridBase< DIMENSION >& cartesian_grid_base_;
    };
    ALIAS_2D_AND_3D( CartesianGridBaseBuilder );

    template < index_t DIMENSION >
    class CartesianGridBuilder final
        : public CartesianGridBaseBuilder< DIMENSION >
    {
    };
    ALIAS_2D_AND_3D( CartesianGridBuilder );

    /*!
     * 3D instantiation of the builder for Cartesian grids, in order to update
     * the grid cage when applying modifications to it.
     */
    template <>
    class CartesianGridBuilder< 3 > final : public CartesianGridBaseBuilder< 3 >
    {
    public:
        virtual ~CartesianGridBuilder() = default;

        explicit CartesianGridBuilder( CartesianGrid< 3 >& cartesian_grid )
            : CartesianGridBaseBuilder< 3 >( cartesian_grid ),
              cartesian_grid_( cartesian_grid )
        {
        }

        void resize_vec_axis( index_t axis_id, double new_size )
        {
            CartesianGridBaseBuilder< 3 >::resize_vec_axis( axis_id, new_size );
            cartesian_grid_.grid_cage_[axis_id * 2 + 1].origin =
                cartesian_grid_.cartesian_frame_.origin()
                + cartesian_grid_.cartesian_frame_[axis_id]
                      * cartesian_grid_.nb_cells_axis( axis_id );
        }

        void change_frame( ReferenceFrame< 3 >& vec_cartesian_axis )
        {
            cartesian_grid_.change_frame( vec_cartesian_axis );
            cartesian_grid_.create_grid_cage();
        }

        void remove_section_from_cartesian_grid(
            index_t axis_id, index_t section_position )
        {
            CartesianGridBaseBuilder< 3 >::remove_section_from_cartesian_grid(
                axis_id, section_position );
            cartesian_grid_.grid_cage_[axis_id * 2 + 1].origin -=
                cartesian_grid_.cartesian_frame_[axis_id];
        }

    private:
        CartesianGrid< 3 >& cartesian_grid_;
    };
}
