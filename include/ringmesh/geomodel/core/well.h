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

#include <ringmesh/geomodel/core/common.h>

#include <array>
#include <memory>

#include <geogram/basic/attributes.h>

#include <ringmesh/basic/geometry.h>

/*!
 * @file Well related classe declarations
 * @author Arnaud Botella
 */

namespace RINGMesh
{
    FORWARD_DECLARATION_DIMENSION_CLASS( GeoModel );
    FORWARD_DECLARATION_DIMENSION_CLASS( Well );
    FORWARD_DECLARATION_DIMENSION_CLASS( NNSearch );
    FORWARD_DECLARATION_DIMENSION_CLASS( PointSetMesh );
    FORWARD_DECLARATION_DIMENSION_CLASS( LineMesh );
    struct ElementLocalVertex;
} // namespace RINGMesh

namespace RINGMesh
{
    template < index_t DIMENSION >
    class geomodel_core_api WellEntity
    {
        ringmesh_disable_copy_and_move( WellEntity );
        ringmesh_template_assert_2d_or_3d( DIMENSION );

    protected:
        explicit WellEntity( const Well< DIMENSION >* well );
        virtual ~WellEntity() = default;

    public:
        /*!
         * Gets the associated well
         */
        const Well< DIMENSION >& well() const
        {
            return *well_;
        }

    protected:
        /// Pointer to the Well owning this entity
        const Well< DIMENSION >* well_;
    };

    ALIAS_2D_AND_3D( WellEntity );

    // --------------------------------------------------------------------------

    template < index_t DIMENSION >
    class geomodel_core_api WellCorner : public WellEntity< DIMENSION >
    {
    public:
        WellCorner( const Well< DIMENSION >* well,
            const vecn< DIMENSION >& point,
            bool is_on_surface,
            index_t id );

        ~WellCorner();

        const vecn< DIMENSION >& point() const;

        bool is_on_surface() const
        {
            return is_on_surface_;
        }

        bool id() const
        {
            return id_;
        }

        GEO::AttributesManager& vertex_attribute_manager() const;

    private:
        /// True is the corner is on a surface, false if is in a region
        bool is_on_surface_;
        /// The id of the corresponding surface or region
        index_t id_;
        std::unique_ptr< PointSetMesh< DIMENSION > > mesh_;
    };

    ALIAS_2D_AND_3D( WellCorner );

    // --------------------------------------------------------------------------

    template < index_t DIMENSION >
    class geomodel_core_api WellPart : public WellEntity< DIMENSION >
    {
    public:
        /*!
         * Create a WellPart
         * @param[in] well the associated well
         * @param[in] id the position in the parts_ vector of the associated
         * well
         */
        WellPart( const Well< DIMENSION >* well, index_t id );

        ~WellPart();

        /*!
         * Sets the corber id
         * @param[in] c the corner id (0 or 1)
         * @param[in] id the corner id in the corners_ vector the the well
         */
        void set_corner( index_t c, index_t id )
        {
            corners_[c] = id;
        }
        /*!
         * Gets the id of a corner
         * @param[in] c the corner id (0 or 1)
         * @return the corresponding id
         */
        index_t corner( index_t c ) const
        {
            ringmesh_assert( c < 2 );
            return corners_[c];
        }

        /*!
         * Create the associated Mesh of the part
         * @param[in] points the points of the mesh
         * @pre the points should be oriented in the order of the well path
         */
        void set_points( const std::vector< vecn< DIMENSION > >& points );

        /*!
         * Gets the number of edges
         */
        index_t nb_edges() const;
        index_t nb_vertices() const;

        /*!
         * Gets the length of the part
         */
        double length() const;

        /*!
         * Sets the id of the part corresponding to the position in the parts_
         * vector of the well
         * @param[in] id the id to set
         */
        void set_id( index_t id )
        {
            id_ = id;
        }
        /*!
         * Gets the id of the part
         */
        index_t id() const
        {
            return id_;
        }
        const vecn< DIMENSION >& vertex( index_t v ) const;
        const vecn< DIMENSION >& edge_vertex(
            const ElementLocalVertex& well_edge_local_vertex ) const;

        GEO::AttributesManager& vertex_attribute_manager() const;
        GEO::AttributesManager& edge_attribute_manager() const;

        const NNSearch< DIMENSION >& vertices_nn_search() const;

    private:
        /// id of the part corresponding to the position in the parts_ vector of
        /// the well
        index_t id_;
        /// id in the corners_ vector the the well
        std::array< index_t, 2 > corners_;
        std::unique_ptr< LineMesh< DIMENSION > > mesh_;
    };

    ALIAS_2D_AND_3D( WellPart );

    // --------------------------------------------------------------------------

    template < index_t DIMENSION >
    class Edge
    {
    public:
        Edge( const vecn< DIMENSION >& v0, const vecn< DIMENSION >& v1 )
        {
            vertices_[0] = v0;
            vertices_[1] = v1;
        }

        const vecn< DIMENSION >& vertex( index_t i ) const
        {
            return vertices_[i];
        }

        vecn< DIMENSION > barycenter() const
        {
            return ( vertices_[0] + vertices_[1] ) * 0.5;
        }

    private:
        std::array< vecn< DIMENSION >, 2 > vertices_;
    };

    ALIAS_2D_AND_3D( Edge );

    // --------------------------------------------------------------------------

    template < index_t DIMENSION >
    class geomodel_core_api Well
    {
        ringmesh_disable_copy_and_move( Well );

    public:
        Well() = default;
        ~Well() = default;

        /*!
         * Copies information and resize the number of parts and corners
         * @param[in,out] well the current well information will be copied into
         * this one
         */
        void copy_corners_and_informations( Well< DIMENSION >& well ) const;

        /*!
         * Gets the edges of a part
         * @param[in] part_id the part id
         * @param[in,out] edges the edges of the part
         */
        void get_part_edges(
            index_t part_id, std::vector< Edge< DIMENSION > >& edges ) const;
        /*!
         * Gets all the edges of a corresponding region
         * @param[in] region the region id
         * @param[in,out] edges the corresponding edges
         */
        void get_region_edges(
            index_t region, std::vector< Edge< DIMENSION > >& edges ) const;

        /*!
         * Creates a new corner
         * @param[in] vertex the geometric position of the corner
         * @param[in] corner_info the corner_info_t corresponding to the corner
         * to create
         * @return the id of the created corner
         */
        index_t create_corner(
            const vecn< DIMENSION >& vertex, bool is_on_surface, index_t id )
        {
            index_t corner_id = static_cast< index_t >( corners_.size() );
            corners_.emplace_back( new WellCorner< DIMENSION >(
                this, vertex, is_on_surface, id ) );
            return corner_id;
        }
        /*!
         * Finds if a corner at a given geometric position exist
         * @param[in] vertex the geometric position to test
         * @param[in] epsilon the numerical precision used to compare the
         * vertices
         * @return the id of the corner or NO_ID if not found any corresponding
         * to \p p
         */
        index_t find_corner(
            const vecn< DIMENSION >& vertex, double epsilon ) const;
        /*!
         * Gets a corner
         * @param[in] c the id of the corner
         */
        const WellCorner< DIMENSION >& corner( index_t c ) const
        {
            ringmesh_assert( c < corners_.size() );
            return *corners_[c];
        }

        /*!
         * Creates a new part
         * @param[in] region the region id corresponding to the new part
         * @return the id of the created part
         */
        index_t create_part( index_t region )
        {
            index_t part_id = static_cast< index_t >( parts_.size() );
            parts_.emplace_back( new WellPart< DIMENSION >( this, part_id ) );
            part_region_id_.push_back( region );
            return part_id;
        }
        /*!
         * Gets a part
         * @param[in] part_id the part id
         */
        const WellPart< DIMENSION >& part( index_t part_id ) const
        {
            ringmesh_assert( part_id < parts_.size() );
            return *parts_[part_id];
        }
        /*!
         * Gets a part
         * @param[in] part_id the part id
         */
        WellPart< DIMENSION >& part( index_t part_id )
        {
            ringmesh_assert( part_id < parts_.size() );
            return *parts_[part_id];
        }
        /*!
         * Gets the region id of a part
         * @param[in] part_id the part id
         * @return the region id of the part
         */
        index_t part_region_id( index_t part_id ) const
        {
            ringmesh_assert( part_id < nb_parts() );
            return part_region_id_[part_id];
        }

        /*!
         * Gets the number of corners
         */
        index_t nb_corners() const
        {
            return static_cast< index_t >( corners_.size() );
        }
        /*!
         * Gets the number of parts
         */
        index_t nb_parts() const
        {
            return static_cast< index_t >( parts_.size() );
        }
        /*!
         * Gets the number of edges of the well
         */
        index_t nb_edges() const;
        /*!
         * Sets the well name
         * @param[in] name the name to set
         */
        void set_name( const std::string& name )
        {
            name_ = name;
        }
        /*!
         * Gets the well name
         */
        const std::string& name() const
        {
            return name_;
        }

    private:
        /// Vector of the corners of the well
        std::vector< std::unique_ptr< WellCorner< DIMENSION > > > corners_;
        /// Vector of the parts of the well
        std::vector< std::unique_ptr< WellPart< DIMENSION > > > parts_;
        /// Vector of the region id of the parts
        std::vector< index_t > part_region_id_;
        /// Name of the well
        std::string name_;
        /// Number of edges in the well
        mutable index_t nb_edges_{ NO_ID };
    };

    ALIAS_2D_AND_3D( Well );

    // --------------------------------------------------------------------------

    /*!
     * Set of wells associated to a GeoModel
     */
    template < index_t DIMENSION >
    class geomodel_core_api WellGroup
    {
        ringmesh_disable_copy_and_move( WellGroup );

    public:
        WellGroup();
        virtual ~WellGroup() = default;

        /*!
         * Gets all the edges contained in a region
         * @param[in] region the region id
         * @param[in,out] edges the edges of the region
         */
        void get_region_edges(
            index_t region, std::vector< Edge< DIMENSION > >& edges ) const;

        /*!
         * Gets all the edges contained in a region
         * @param[in] region the region id
         * @param[in,out] edges the edges of the region, one vector per well
         */
        void get_region_edges( index_t region,
            std::vector< std::vector< Edge< DIMENSION > > >& edges ) const;

        /*!
         * Gets the associated GeoModel
         */
        const GeoModel< DIMENSION >* geomodel() const
        {
            return geomodel_;
        }
        /*!
         * Sets the associated GeoModel
         */
        void set_geomodel( GeoModel< DIMENSION >* geomodel )
        {
            geomodel_ = geomodel;
        }

        /*!
         * Finds if a well with the same name already exist
         * @param[in] name the name to test
         * @return the id of the well or NO_ID
         */
        index_t find_well( const std::string& name ) const;

        /*!
         * Creates new wells
         * @param[in] nb_wells the number of wells to create
         */
        void create_wells( index_t nb_wells );

        /*!
         * Add a well from its mesh and makes it conformal to the associated
         * GeoModel
         * @param[in] mesh the mesh of the well
         * @param[in] name the name of the well
         */
        void add_well(
            const LineMesh< DIMENSION >& mesh, const std::string& name );

        /*!
         * Gets the number of wells
         */
        index_t nb_wells() const
        {
            return static_cast< index_t >( wells_.size() );
        }
        /*!
         * Gets the well
         * @param[in] w the well id
         * @return the corresponding well
         */
        const Well< DIMENSION >& well( index_t w ) const
        {
            return *wells_[w];
        }

    private:
        void compute_conformal_mesh(
            const LineMesh< DIMENSION >& in, LineMesh< DIMENSION >& out );

    private:
        /// Vector of the wells
        std::vector< Well< DIMENSION >* > wells_;
        /// Associated GeoModel
        GeoModel< DIMENSION >* geomodel_;
    };

    ALIAS_2D_AND_3D( WellGroup );
} // namespace RINGMesh
