/*
 * Copyright (c) 2012-2016, Association Scientifique pour la Geologie et ses Applications (ASGA)
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

#include <ringmesh/geo_model_builder_from_mesh.h>

#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <set>
#include <stack>

#include <geogram/basic/line_stream.h>
#include <geogram/basic/logger.h>
#include <geogram/mesh/mesh_io.h>
#include <geogram/mesh/mesh_repair.h>
#include <geogram/points/colocate.h>

#include <ringmesh/io.h>
#include <ringmesh/algorithm.h>
#include <ringmesh/geo_model.h>
#include <ringmesh/geo_model_api.h>
#include <ringmesh/geo_model_validity.h>
#include <ringmesh/geometry.h>
#include <ringmesh/geogram_mesh_repair.h>
#include <ringmesh/utils.h>

/*!
 * @file ringmesh/geo_model_builder.cpp
 * @brief Implementation of the classes to build GeoModel from various inputs
 * @author Jeanne Pellerin
 */

namespace RINGMesh {

    bool is_surface_mesh( const GEO::Mesh& mesh )
    {
        return mesh.facets.nb() != 0 ;
    }

    bool is_volume_mesh( const GEO::Mesh& mesh )
    {
        return mesh.cells.nb() != 0 ;
    }


    /*************************************************************************/

    /*!
     * @brief Implementation detail: abstract base class to create a GeoModelEntities 
     *        from SIMPLICIAL meshes
     * @details Manages the correspondence between a Mesh entities and
     *          GeoModelEntities only known by indices.
     *
     * @warning Implemented only for TRIANGULATED Surface and TETRAHEDRALIZED Region.
     * @note Used by GeoModelBuilderMesh.
     */
    class GeoModelEntityFromMesh {
    public:
        typedef std::pair< index_t, index_t > index_pair ;
        typedef std::map< index_t, index_t > index_map ;

        virtual ~GeoModelEntityFromMesh()
        {
        }

        /*! Check that the attribute is defined.
         * If not, returns false otherwise bind it.
         */
        void initialize()
        {
            if( is_gme_attribute_defined() ) {
                bind_gme_attribute() ;
            }
        }

        bool is_valid()
        {
            bool attribute_is_bounded = gme_attribute_.is_bound() ;
            return attribute_is_bounded ;
        }

        index_t count_attribute_values_and_simplexes()
        {
            if( !is_valid() ) {
                return 0 ;
            }
            index_t nb = nb_mesh_simplexes() ;
            for( index_t i = 0; i != nb; ++i ) {
                index_t value = gme_attribute_[i] ;
                if( !is_attribute_value( value ) ) {
                    nb_simplexes_per_attribute_value_[value] = 0 ;
                }
                ++nb_simplexes_per_attribute_value_[value] ;
            }
            return nb_attribute_values() ;
        }

        /*! Number of different values for the attribute */
        index_t nb_attribute_values()
        {
            return index_t( nb_simplexes_per_attribute_value_.size() ) ;
        }

        /*! Sets a mapping from the attribute values on the Mesh and 
         * the indices of the GeoModelEntities to fill
         */
        void set_gme_id_attribute_mapping(
            const std::vector< index_t >& gme_id_to_attribute_in )
        {
            gme_id_to_attribute_value_ = gme_id_to_attribute_in ;

            for( index_t i = 0; i != nb_gme(); ++i ) {
                index_t value = attribute_value_from_gme( i ) ;
                if( is_attribute_value( value ) ) {
                    attribute_value_to_gme_id_[value] = i ;
                } else {
                    GEO::Logger::err( "Debug" )
                        << "Invalid mapping between Mesh attribute and GeoModelEntity"
                        << std::endl ;
                    gme_id_to_attribute_value_[i] = NO_ID ;
                }
            }
        }

        /*! Default mapping: first GeoModelEntity index corresponds to first attribute
         *  value, second to second, etc.
         */
        void set_default_gme_id_attribute_mapping( index_t nb_geomodel_entities )
        {
            std::vector< index_t > default_mapping( nb_geomodel_entities, NO_ID ) ;
            index_t count = 0 ;
            index_map::const_iterator it(
                nb_simplexes_per_attribute_value_.begin() ) ;
            while( it != nb_simplexes_per_attribute_value_.end()
                && count < nb_geomodel_entities ) {
                default_mapping[count] = it->first ;
                ++count ;
                ++it ;
            }
            set_gme_id_attribute_mapping( default_mapping ) ;
        }

        /*!
         * Computes the simplex vertex indices for each GeoModelEntity
         * as well as the mapping from the mesh simplices to the GME mesh simplices
         * for eventual attribute copying
         */
        void compute_gme_simplexes()
        {
            if( !is_valid() ) {
                return ;
            }
            allocate_mesh_simplex_to_gme() ;
            allocate_gme_vertices() ;

            index_t nb = nb_mesh_simplexes() ;
            std::vector< index_t > gme_simplex_counter_( nb_gme(), 0 ) ;
            for( index_t mesh_simplex = 0; mesh_simplex < nb; ++mesh_simplex ) {
                index_t attribute_value = gme_attribute_[mesh_simplex] ;
                if( attribute_value_has_gme_id( attribute_value ) ) {
                    index_t gme_id = attribute_value_to_gme_id_[attribute_value] ;
                    index_t gme_simplex_id = gme_simplex_counter_[gme_id] ;

                    assign_one_gme_simplex_vertices( mesh_simplex, gme_id,
                        gme_simplex_id ) ;
                    assign_mesh_simplex_to_gme_simplex( mesh_simplex, gme_id,
                        gme_simplex_id ) ;
                    ++gme_simplex_counter_[gme_id] ;
                } else {
                    assign_mesh_simplex_to_no_gme_simplex( mesh_simplex ) ;
                }
            }
        }

        /*!
         * Computes adjacencies between simplices for each GeoModelEntity
         * compute_gme_simplexes SHOULD have been called before, nothing done if it wasn't
         */
        void compute_gme_adjacencies()
        {
            if( gme_simplex_vertices_.empty() ) {
                return ;
            }
            allocate_gme_corner_adjacent_gme_simplex() ;

            index_t nb = nb_mesh_simplexes() ;
            for( index_t mesh_simplex = 0; mesh_simplex < nb; ++mesh_simplex ) {
                GMESimplex gme = mesh_simplex_to_gme_simplex_[mesh_simplex] ;
                index_t gme_id = gme.gme_id ;
                index_t gme_simplex = gme.gme_simplex_id ;

                if( gme_id == NO_ID || gme_simplex == NO_ID ) {
                    continue ;
                }

                for( index_t v = 0; v < nb_vertices_per_simplex(); ++v ) {
                    index_t gme_vertex = nb_vertices_per_simplex() * gme_simplex
                        + v ;
                    index_t mesh_vertex = nb_vertices_per_simplex() * mesh_simplex
                        + v ;
                    // Number of vertices = number of edges / number of facets
                    index_t adjacent_simplex = adjacent_simplex_index(
                        mesh_vertex ) ;

                    GMESimplex adjacent_gme_simplex ;
                    if( adjacent_simplex != NO_ID ) {
                        adjacent_gme_simplex =
                            mesh_simplex_to_gme_simplex_[adjacent_simplex] ;
                    }
                    if( adjacent_gme_simplex.gme_id == gme_id ) {
                        gme_corner_adjacent_gme_simplex_[gme_id][gme_vertex] =
                            adjacent_gme_simplex.gme_simplex_id ;
                    }
                }
            }
        }

        /*! Simplex vertex indices in the Mesh for one GeoModelEntity
         */
        const std::vector< index_t >& gme_simplices( index_t gme_id ) const
        {
            ringmesh_assert( gme_id < nb_gme() ) ;
            return gme_simplex_vertices_[gme_id] ;
        }

        const std::vector< index_t >& adjacent_gme_simplices( index_t gme_id ) const
        {
            ringmesh_assert( gme_id < nb_gme() ) ;
            return gme_corner_adjacent_gme_simplex_[gme_id] ;
        }

        template< typename T >
        void copy_simplex_attribute_from_mesh_to_geomodel(
            GEO::Attribute< T >& mesh_attribute,
            AttributeVector< T >& model_attributes ) const
        {
            for( index_t i = 0; i < nb_mesh_simplexes(); ++i ) {
                const GMESimplex& copy_to = mesh_simplex_to_gme_simplex_[i] ;
                model_attributes[copy_to.gme_id][copy_to.gme_simplex_id] =
                    mesh_attribute[i] ;
            }
        }

    protected:
        // A simplex in a GeoModelEntity
        struct GMESimplex {
            GMESimplex()
                : gme_id( NO_ID ), gme_simplex_id( NO_ID )
            {
            }
            GMESimplex( index_t gme_in, index_t simplex_in )
                : gme_id( gme_in ), gme_simplex_id( simplex_in )
            {
            }
            index_t gme_id ;
            index_t gme_simplex_id ;
        } ;

    protected:
        GeoModelEntityFromMesh(
            const GEO::Mesh& M,
            const std::string& attribute_name )
            : mesh_( M ), gme_attribute_name_( attribute_name )
        {
        }

        /*! Number of GeoModelEntities of the considered type */
        index_t nb_gme() const
        {
            return gme_id_to_attribute_value_.size() ;
        }

        /*! Is the value indeed taken by the attribute on the mesh */
        bool is_attribute_value( index_t value )
        {
            return nb_simplexes_per_attribute_value_.count( value ) == 1 ;
        }

        bool attribute_value_has_gme_id( index_t attribute_value ) const
        {
            return attribute_value_to_gme_id_.count( attribute_value ) == 1 ;
        }

        index_t attribute_value_from_gme( index_t gme_id ) const
        {
            return gme_id_to_attribute_value_[gme_id] ;
        }

        bool is_gme_attribute_defined()
        {
            GEO::AttributesManager& manager = mesh_simplex_attribute_manager() ;
            return manager.is_defined( gme_attribute_name_ ) ;
        }

        void bind_gme_attribute()
        {
            GEO::AttributesManager& manager = mesh_simplex_attribute_manager() ;
            gme_attribute_.bind( manager, gme_attribute_name_ ) ;
        }

        void allocate_gme_vertices()
        {
            gme_simplex_vertices_.resize( nb_gme() ) ;
            for( index_t i = 0; i < nb_gme(); ++i ) {
                index_t value = attribute_value_from_gme( i ) ;
                if( is_attribute_value( value ) ) {
                    index_t nb_simplexes = nb_simplexes_per_attribute_value_[value] ;
                    gme_simplex_vertices_[i].resize(
                        nb_vertices_per_simplex() * nb_simplexes ) ;
                }
            }
        }

        void allocate_gme_corner_adjacent_gme_simplex()
        {
            gme_corner_adjacent_gme_simplex_.resize( nb_gme() ) ;
            for( index_t i = 0; i < nb_gme(); ++i ) {
                gme_corner_adjacent_gme_simplex_[i].resize(
                    gme_simplex_vertices_[i].size(), NO_ID ) ;
            }
        }

        void allocate_mesh_simplex_to_gme()
        {
            mesh_simplex_to_gme_simplex_.resize( nb_mesh_simplexes() ) ;
        }

        virtual void assign_one_gme_simplex_vertices(
            index_t mesh_simplex_id,
            index_t gme_id,
            index_t gme_simplex_id )
        {
            index_t from = gme_simplex_id * nb_vertices_per_simplex() ;
            for( index_t v = 0; v != nb_vertices_per_simplex(); ++v ) {
                gme_simplex_vertices_[gme_id][from + v] = mesh_vertex_index(
                    mesh_simplex_id, v ) ;
            }
        }

        void assign_mesh_simplex_to_no_gme_simplex( index_t mesh_simplex_id )
        {
            assign_mesh_simplex_to_gme_simplex( mesh_simplex_id, NO_ID, NO_ID ) ;
        }

        void assign_mesh_simplex_to_gme_simplex(
            index_t mesh_simplex_id,
            index_t gme_id,
            index_t gme_simplex_id )
        {
            mesh_simplex_to_gme_simplex_[mesh_simplex_id] = GMESimplex( gme_id,
                gme_simplex_id ) ;
        }

        virtual index_t nb_mesh_simplexes() const = 0 ;
        virtual GEO::AttributesManager& mesh_simplex_attribute_manager() = 0 ;
        virtual index_t nb_vertices_per_simplex() const = 0 ;
        virtual index_t mesh_vertex_index(
            index_t simplex_id,
            index_t lv ) const = 0 ;
        virtual index_t adjacent_simplex_index(
            index_t facet_or_edge_id ) const = 0 ;

    protected:
        // THE Mesh
        const GEO::Mesh& mesh_ ;
        // Name of the attribute on the Mesh simplices identifying the GMEs
        std::string gme_attribute_name_ ;
        // Attribute giving the GeoModelEntity index on Mesh simplices
        GEO::Attribute< index_t > gme_attribute_ ;
        // Number of simplices of the Mesh with a given attribute value
        std::map< index_t, index_t > nb_simplexes_per_attribute_value_ ;
        // Mapping from the attribute value on the Mesh to a GME index in a GeoModel
        std::map< index_t, index_t > attribute_value_to_gme_id_ ;
        // Mapping from a GME index in a GeoModel to the attribute value on the Mesh  
        std::vector< index_t > gme_id_to_attribute_value_ ;
        // Vertex indices (in the Mesh) of the corners of the simplices per GME
        std::vector< std::vector< index_t > > gme_simplex_vertices_ ;
        // Mapping from a mesh simplex index to a simplex index in a GME
        std::vector< GMESimplex > mesh_simplex_to_gme_simplex_ ;
        // For each GME, store the adjacent simplex (in the GME) for each simplex corner
        std::vector< std::vector< index_t > > gme_corner_adjacent_gme_simplex_ ;
    } ;

    class GeoModelSurfaceFromMesh: public GeoModelEntityFromMesh {
    public:
        GeoModelSurfaceFromMesh(
            const GEO::Mesh& M,
            const std::string& attribute_name )
            : GeoModelEntityFromMesh( M, attribute_name )
        {
        }

        GEO::AttributesManager& mesh_simplex_attribute_manager()
        {
            return mesh_.facets.attributes() ;
        }

        index_t nb_mesh_simplexes() const
        {
            return mesh_.facets.nb() ;
        }

        index_t nb_vertices_per_simplex() const
        {
            return 3 ;
        }

        index_t mesh_vertex_index( index_t simplex_id, index_t vertex ) const
        {
            return mesh_.facets.vertex( simplex_id, vertex ) ;
        }
        index_t adjacent_simplex_index( index_t corner_id ) const
        {
            return mesh_.facet_corners.adjacent_facet( corner_id ) ;
        }
    } ;

    class GeoModelRegionFromMesh: public GeoModelEntityFromMesh {
    public:
        GeoModelRegionFromMesh(
            const GEO::Mesh& M,
            const std::string& attribute_name )
            : GeoModelEntityFromMesh( M, attribute_name )
        {
        }

        GEO::AttributesManager& mesh_simplex_attribute_manager()
        {
            return mesh_.cells.attributes() ;
        }

        index_t nb_mesh_simplexes() const
        {
            return mesh_.cells.nb() ;
        }

        virtual index_t nb_vertices_per_simplex() const
        {
            return 4 ;
        }

        index_t mesh_vertex_index( index_t simplex_id, index_t vertex ) const
        {
            return mesh_.cells.vertex( simplex_id, vertex ) ;
        }

        index_t adjacent_simplex_index( index_t facet_id ) const
        {
            return mesh_.cell_facets.adjacent_cell( facet_id ) ;
        }
    } ;

    /*************************************************************************/

    GeoModelBuilderMesh::~GeoModelBuilderMesh()
    {
        delete surface_builder_ ;
        surface_builder_ = nil ;
        delete region_builder_ ;
        region_builder_ = nil ;
    }

    bool GeoModelBuilderMesh::is_mesh_valid_for_surface_building() const
    {
        bool valid = true ;
        if( !is_surface_mesh( mesh_ ) ) {
            Logger::warn( "GMBuilder" ) << "The given mesh is not a surface mesh "
                << std::endl ;
            valid = false ;
        }
        if( !mesh_.facets.are_simplices() ) {
            Logger::warn( "GMBuilder" ) << "The given mesh is not triangulated"
                << std::endl ;
            valid = false ;
        }
        if( !mesh_.facets.attributes().is_defined( surface_attribute_name_ ) ) {
            Logger::warn( "GMBuilder" ) << "The attribute "
                << surface_attribute_name_ << " is not defined on the Mesh facets"
                << std::endl ;
            valid = false ;
        }
        if( RINGMesh::has_mesh_colocate_vertices( mesh_, epsilon ) ) {
            Logger::warn( "GMBuilder" )
                << " The Mesh has colocated vertices. Repair it beforehand "
                << std::endl ;
            valid = false ;
        }
        return valid ;
    }

    bool GeoModelBuilderMesh::is_mesh_valid_for_region_building() const
    {
        bool valid = true ;
        if( !is_volume_mesh( mesh_ ) ) {
            Logger::warn( "GMBuilder" ) << "The given mesh is not a volumetric mesh "
                << std::endl ;
            valid = false ;
        }
        if( !mesh_.cells.are_simplices() ) {
            Logger::warn( "GMBuilder" ) << "The given mesh is not tetrahedralized"
                << std::endl ;
            valid = false ;
        }
        if( !mesh_.cells.attributes().is_defined( region_attribute_name_ ) ) {
            Logger::warn( "GMBuilder" ) << "The attribute " << region_attribute_name_
                << " is not defined on the Mesh cells " << std::endl ;
            valid = false ;
        }
        if( RINGMesh::has_mesh_colocate_vertices( mesh_, epsilon ) ) {
            Logger::warn( "GMBuilder" )
                << " The Mesh has colocated vertices. Repair it beforehand. "
                << std::endl ;
            valid = false ;
        }
        return valid ;
    }

    void create_and_fill_connected_component_attribute(
        GEO::Mesh& mesh,
        const std::string& connected_component_attribute )
    {
        GEO::Attribute< index_t > connected_component ;
        GEO::AttributesManager& manager = mesh.facets.attributes() ;
        connected_component.bind( manager, connected_component_attribute ) ;

        index_t nb_facets = mesh.facets.nb() ;
        std::vector< bool > visited( nb_facets, false ) ;

        ///@todo This algorithm is implemented over and over again in RINGMesh
        /// and Geogram. Couldn't we do better ?
        index_t nb_connected_components = 0 ;
        for( index_t f = 0; f < nb_facets; f++ ) {
            if( !visited[f] ) {
                nb_connected_components++ ;

                std::stack< index_t > facet_stack ;
                facet_stack.push( f ) ;

                while( !facet_stack.empty() ) {
                    index_t f_from_stack = facet_stack.top() ;
                    facet_stack.pop() ;
                    visited[f_from_stack] = true ;
                    connected_component[f_from_stack] = nb_connected_components ;

                    for( index_t v = 0; v < 3; ++v ) {
                        index_t neighbor_facet = mesh.facets.adjacent( f_from_stack,
                            v ) ;
                        if( neighbor_facet != NO_ID && !visited[neighbor_facet] ) {
                            visited[neighbor_facet] = true ;
                            facet_stack.push( neighbor_facet ) ;
                        }
                    }
                }
            }
        }
    }

    void GeoModelBuilderMesh::prepare_surface_mesh_from_connected_components(
        GEO::Mesh& mesh,
        const std::string& created_facet_attribute )
    {
        // Remove duplicate facets, and triangulate the mesh
        // Side effects: fixes facet orientation and split non-manifold vertices
        // AND empty all attributes.
        GEO::mesh_repair( mesh,
            GEO::MeshRepairMode(
                GEO::MESH_REPAIR_DUP_F | GEO::MESH_REPAIR_TRIANGULATE ) ) ;

        create_and_fill_connected_component_attribute( mesh,
            created_facet_attribute ) ;

        // Remove colocated vertices
        repair_colocate_vertices( mesh, epsilon ) ;
    }

    /*! @details Adds separately each connected component of the mesh
     *          as a Surface of the model under construction.
     *          All the facets of the input mesh are visited and added to a
     *          Surface of the GeoModel.
     *          Connected components of the mesh are determined with a
     *          propagation (or "coloriage" algorithm) using the adjacent_facet
     *          information provided on the input GEO::Mesh.
     *
     * @todo Old code - old building - to delimit connected components
     * vertices are duplicated in the input mesh
     *
     */
    void GeoModelBuilderSurfaceMesh::build_polygonal_surfaces_from_connected_components()
    {
        std::vector< index_t > global_vertex_id_to_id_in_cc( mesh_.vertices.nb(),
            NO_ID ) ;

        std::vector< bool > visited( mesh_.facets.nb(), false ) ;
        for( index_t i = 0; i < mesh_.facets.nb(); i++ ) {
            if( !visited[i] ) {
                std::vector< index_t > cc_corners ;
                std::vector< index_t > cc_facets_ptr ;
                std::vector< vec3 > cc_vertices ;

                /// @todo Review : This should not be necessary as each vertex should
                /// be in one and only one connected component. To test. [JP]
                std::fill( global_vertex_id_to_id_in_cc.begin(),
                    global_vertex_id_to_id_in_cc.end(), NO_ID ) ;

                // First facet begin at corner 0
                cc_facets_ptr.push_back( 0 ) ;

                // Propagate from facet #i 
                std::stack< index_t > S ;
                S.push( i ) ;
                while( !S.empty() ) {
                    index_t f = S.top() ;
                    S.pop() ;
                    visited[f] = true ;

                    for( index_t c = mesh_.facets.corners_begin( f );
                        c < mesh_.facets.corners_end( f ); ++c ) {
                        index_t v = mesh_.facet_corners.vertex( c ) ;
                        if( global_vertex_id_to_id_in_cc[v] == NO_ID ) {
                            global_vertex_id_to_id_in_cc[v] = cc_vertices.size() ;
                            cc_vertices.push_back( mesh_.vertices.point( v ) ) ;
                        }
                        cc_corners.push_back( global_vertex_id_to_id_in_cc[v] ) ;

                        index_t n = mesh_.facet_corners.adjacent_facet( c ) ;
                        if( n != NO_ID && !visited[n] ) {
                            visited[n] = true ;
                            S.push( n ) ;
                        }
                    }
                    cc_facets_ptr.push_back( cc_corners.size() ) ;
                }

                gme_t surface_gme = create_entity( GME::SURFACE ) ;
                set_surface_geometry( surface_gme.index, cc_vertices, cc_corners,
                    cc_facets_ptr ) ;
            }
        }
    }

    void GeoModelBuilderMesh::create_and_build_surfaces()
    {
        create_geomodel_entities( GME::SURFACE, nb_surface_attribute_values_ ) ;
        build_surfaces() ;
    }

    void GeoModelBuilderMesh::build_surfaces()
    {
        if( !is_mesh_valid_for_surface_building() ) {
            return ;
        }

        index_t nb_surfaces = model().nb_surfaces() ;
        surface_builder_->set_default_gme_id_attribute_mapping( nb_surfaces ) ;
        surface_builder_->compute_gme_simplexes() ;
        surface_builder_->compute_gme_adjacencies() ;
        for( index_t i = 0; i != nb_surfaces; ++i ) {
            const std::vector< index_t >& triangle_vertices =
                surface_builder_->gme_simplices( i ) ;
            const std::vector< index_t >& adjacent_triangles =
                surface_builder_->adjacent_gme_simplices( i ) ;
            // Set the Surface facets
            // Set the adjacencies so that internal boundaries are not lost 
            // when connecting the surface facets.
            set_surface_geometry_with_adjacencies( i, triangle_vertices,
                adjacent_triangles ) ;
        }
    }

    void GeoModelBuilderMesh::create_and_build_regions()
    {
        create_geomodel_entities( GME::REGION, nb_region_attribute_values_ ) ;
        build_regions() ;
    }

    void GeoModelBuilderMesh::build_regions()
    {
        if( !is_mesh_valid_for_region_building() ) {
            return ;
        }

        index_t nb_regions = model().nb_regions() ;
        region_builder_->set_default_gme_id_attribute_mapping( nb_regions ) ;

        region_builder_->compute_gme_simplexes() ;
        for( index_t i = 0; i != nb_regions; ++i ) {
            const std::vector< index_t >& tet_vertices =
                region_builder_->gme_simplices( i ) ;
            set_region_geometry( i, tet_vertices ) ;
        }
    }

    void GeoModelBuilderMesh::add_mesh_vertices_to_model()
    {
        index_t nb_vertices = mesh_.vertices.nb() ;
        for( index_t i = 0; i < nb_vertices; ++i ) {
            model().mesh.vertices.add_vertex( mesh_.vertices.point( i ) ) ;
        }
    }

    void GeoModelBuilderMesh::initialize_surface_builder()
    {
        surface_builder_ = new GeoModelSurfaceFromMesh( mesh_,
            surface_attribute_name_ ) ;
        surface_builder_->initialize() ;
        nb_surface_attribute_values_ =
            surface_builder_->count_attribute_values_and_simplexes() ;
    }

    void GeoModelBuilderMesh::initialize_region_builder()
    {
        region_builder_ = new GeoModelRegionFromMesh( mesh_,
            region_attribute_name_ ) ;
        region_builder_->initialize() ;
        nb_region_attribute_values_ =
            region_builder_->count_attribute_values_and_simplexes() ;
    }

    void GeoModelBuilderMesh::copy_facet_attribute_from_mesh(
        const std::string& attribute_name )
    {
        if( !is_facet_attribute_defined< index_t >( mesh_, attribute_name ) ) {
            GEO::Logger::warn( "GMBuilder" ) << "No INDEX_T attribute named "
                << attribute_name << " on mesh facets to copy " << std::endl ;
            return ;
        }
        GEO::Attribute< index_t > attribute( mesh_.facets.attributes(),
            attribute_name ) ;
        AttributeVector< index_t > attributes ;
        create_attributes_on_geomodel_entity_facets< index_t >( model(),
            GME::SURFACE, attribute_name, attributes ) ;
        surface_builder_->copy_simplex_attribute_from_mesh_to_geomodel< index_t >(
            attribute, attributes ) ;
    }

    void GeoModelBuilderMesh::copy_cell_attribute_from_mesh(
        const std::string& attribute_name )
    {
        if( !is_cell_attribute_defined< index_t >( mesh_, attribute_name ) ) {
            GEO::Logger::warn( "GMBuilder" ) << "No INDEX_T attribute named "
                << attribute_name << " on mesh cells to copy " << std::endl ;
            return ;
        }
        GEO::Attribute< index_t > attribute( mesh_.cells.attributes(),
            attribute_name ) ;
        AttributeVector< index_t > attributes ;
        create_attributes_on_geomodel_entity_cells< index_t >( model(), GME::REGION,
            attribute_name, attributes ) ;
        region_builder_->copy_simplex_attribute_from_mesh_to_geomodel< index_t >(
            attribute, attributes ) ;
    }

} // namespace
