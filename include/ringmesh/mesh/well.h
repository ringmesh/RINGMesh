/*
 * Copyright (c) 2012-2017, Association Scientifique pour la Geologie et ses Applications (ASGA)
 * All rights reserved.
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
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL ASGA BE LIABLE FOR ANY
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
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

/*!
 * @file Well related classe declarations 
 * @author Arnaud Botella
 */

namespace RINGMesh {
    class GeoModel ;
    class Well ;
    class Mesh0D ;
    class Mesh1D ;
}

namespace RINGMesh {

    class RINGMESH_API WellEntity {
    ringmesh_disable_copy( WellEntity ) ;
    protected:
        WellEntity( const Well* well ) ;
        virtual ~WellEntity()
        {
        }

    public:
        /*!
         * Gets the associated well
         */
        const Well& well() const
        {
            return *well_ ;
        }

    protected:
        /// Pointer to the Well owning this entity
        const Well* well_ ;
    } ;

// --------------------------------------------------------------------------

    class RINGMESH_API WellCorner: public WellEntity {
    public:
        WellCorner(
            const Well* well,
            const vec3& point,
            bool is_on_surface,
            index_t id ) ;
        virtual ~WellCorner() ;

        const vec3& point() const ;

        bool is_on_surface() const
        {
            return is_on_surface_ ;
        }

        bool id() const
        {
            return id_ ;
        }

    private:
        /// True is the corner is on a surface, false if is in a region
        bool is_on_surface_ ;
        /// The id of the corresponding surface or region
        index_t id_ ;
        Mesh0D* mesh_ ;
    } ;

// --------------------------------------------------------------------------

    class RINGMESH_API WellPart: public WellEntity {
    public:
        WellPart( const Well* well, index_t id ) ;
        ~WellPart() ;

        /*!
         * Sets the corber id
         * @param[in] c the corner id (0 or 1)
         * @param[in] id the corner id in the corners_ vector the the well
         */
        void set_corner( index_t c, index_t id )
        {
            corners_[c] = id ;
        }
        /*!
         * Gets the id of a corner
         * @param[in] c the corner id (0 or 1)
         * @return the corresponding id
         */
        index_t corner( index_t c ) const
        {
            ringmesh_assert( c < 2 ) ;
            return corners_[c] ;
        }

        void set_points( const std::vector< vec3 >& points ) ;

        index_t nb_edges() const ;
        double length() const ;

        /*!
         * Sets the id of the part corresponding to the position in the parts_ vector of the well
         * @param[in] id the id to set
         */
        void set_id( index_t id )
        {
            id_ = id ;
        }
        /*!
         * Gets the id of the part
         */
        index_t id() const
        {
            return id_ ;
        }
        const vec3& vertex( index_t v ) const ;
        const vec3& edge_vertex( index_t edge, index_t v ) const ;

    private:
        /// id of the part corresponding to the position in the parts_ vector of the well
        index_t id_ ;
        /// id in the corners_ vector the the well
        index_t corners_[2] ;
        Mesh1D* mesh_ ;
    } ;

// --------------------------------------------------------------------------

    class RINGMESH_API Edge {
    public:
        Edge( const vec3& v0, const vec3& v1 )
        {
            vertices_[0] = v0 ;
            vertices_[1] = v1 ;
        }

        const vec3& vertex( index_t i ) const
        {
            return vertices_[i] ;
        }

        vec3 barycenter() const
        {
            return ( vertices_[0] + vertices_[1] ) * 0.5 ;
        }
    private:
        vec3 vertices_[2] ;
    } ;

// --------------------------------------------------------------------------

    class RINGMESH_API Well {
    ringmesh_disable_copy( Well ) ;
    public:
        Well() ;
        ~Well() ;

        void copy_corners_and_informations( Well& well ) const ;

        void get_part_edges( index_t part_id, std::vector< Edge >& edges ) const ;
        void get_region_edges( index_t part_id, std::vector< Edge >& edges ) const ;

        /*!
         * Creates a new corner
         * @param[in] vertex the geometric position of the corner
         * @param[in] corner_info the corner_info_t corresponding to the corner to create
         * @return the id of the created corner
         */
        index_t create_corner( const vec3& vertex, bool is_on_surface, index_t id )
        {
            index_t corner_id = static_cast< index_t >( corners_.size() ) ;
            corners_.push_back( new WellCorner( this, vertex, is_on_surface, id ) ) ;
            return corner_id ;
        }
        index_t find_corner( const vec3& vertex, double epsilon ) const ;
        /*!
         * Gets a corner
         * @param[in] c the id of the corner
         */
        const WellCorner& corner( index_t c ) const
        {
            ringmesh_assert( c < corners_.size() ) ;
            return *corners_[c] ;
        }

        /*!
         * Creates a new part
         * @param[in] region the region id corresponding to the new part
         * @return the id of the created part
         */
        index_t create_part( index_t region )
        {
            index_t part_id = static_cast< index_t >( parts_.size() ) ;
            parts_.push_back( new WellPart( this, part_id ) ) ;
            part_region_id_.push_back( region ) ;
            return part_id ;
        }
        /*!
         * Gets a part
         * @param[in] part_id the part id
         */
        const WellPart& part( index_t part_id ) const
        {
            ringmesh_assert( part_id < parts_.size() ) ;
            return *parts_[part_id] ;
        }
        /*!
         * Gets a part
         * @param[in] part_id the part id
         */
        WellPart& part( index_t part_id )
        {
            ringmesh_assert( part_id < parts_.size() ) ;
            return *parts_[part_id] ;
        }
        /*!
         * Gets the region id of a part
         * @param[in] part_id the part id
         * @return the region id of the part
         */
        index_t part_region_id( index_t part_id ) const
        {
            ringmesh_assert( part_id < nb_parts() ) ;
            return part_region_id_[part_id] ;
        }

        /*!
         * Gets the number of corners
         */
        index_t nb_corners() const
        {
            return static_cast< index_t >( corners_.size() ) ;
        }
        /*!
         * Gets the number of parts
         */
        index_t nb_parts() const
        {
            return static_cast< index_t >( parts_.size() ) ;
        }
        index_t nb_edges() const ;
        /*!
         * Sets the well name
         * @param[in] name the name to set
         */
        void set_name( const std::string& name )
        {
            name_ = name ;
        }
        /*!
         * Gets the well name
         */
        const std::string& name() const
        {
            return name_ ;
        }

    private:
        /// Vector of the corners of the well
        std::vector< WellCorner* > corners_ ;
        /// Vector of the parts of the well
        std::vector< WellPart* > parts_ ;
        /// Vector of the region id of the parts
        std::vector< index_t > part_region_id_ ;
        /// Name of the well
        std::string name_ ;
        /// Number of edges in the well
        index_t nb_edges_ ;
    } ;

// --------------------------------------------------------------------------

    /*! 
     * @todo Comment
     */
    class RINGMESH_API WellGroup {
    ringmesh_disable_copy( WellGroup ) ;
    public:
        WellGroup() ;
        virtual ~WellGroup() ;

        void get_region_edges( index_t region, std::vector< Edge >& edges ) const ;

        void get_region_edges(
            index_t region,
            std::vector< std::vector< Edge > >& edges ) const ;

        /*!
         * Gets the associated GeoModel
         */
        const GeoModel* geomodel() const
        {
            return geomodel_ ;
        }
        /*!
         * Sets the associated GeoModel
         */
        void set_geomodel( RINGMesh::GeoModel* geomodel )
        {
            geomodel_ = geomodel ;
        }
        index_t find_well( const std::string& name ) const ;

        void create_wells( index_t nb_wells ) ;
        void add_well( const Mesh1D& mesh, const std::string& name ) ;

        /*!
         * Gets the number of wells
         */
        index_t nb_wells() const
        {
            return static_cast< index_t >( wells_.size() ) ;
        }
        /*!
         * Gets the well
         * @param[in] w the well id
         * @return the corresponding well
         */
        const Well& well( index_t w ) const
        {
            return *wells_[w] ;
        }

    private:
        void compute_conformal_mesh( const Mesh1D& in, Mesh1D& out ) ;

    protected:
        /// Vector of the wells
        std::vector< Well* > wells_ ;
        /// Associated GeoModel
        GeoModel* geomodel_ ;
    } ;
}
