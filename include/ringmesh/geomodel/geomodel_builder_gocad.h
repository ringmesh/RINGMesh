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

#include <memory>

#include <geogram/basic/factory.h>
#include <geogram/basic/line_stream.h>

#include <ringmesh/geomodel/geomodel_builder_file.h>

namespace RINGMesh {
    class GeoModelBuilderTSolid;
    class GeoModelBuilderTSolidImpl;
    class GeoModelBuilderML;
    struct VertexMap;
    struct TSolidLoadingStorage;

    CLASS_DIMENSION_ALIASES( Box );
}

namespace RINGMesh {

    void RINGMESH_API initialize_gocad_import_factories();

    class RINGMESH_API GeoModelBuilderGocad: public GeoModelBuilderFile< 3 > {
    public:
        GeoModelBuilderGocad( GeoModel3D& geomodel, std::string filename )
            :
                GeoModelBuilderFile( geomodel, std::move( filename ) ),
				file_line_(filename_), file_name_(filename_)
        {
            /*! @todo Review: A constructor is not supposed to throw, the object is left in an
             * undefined state [JP] */
            if( !file_line_.OK() ) {
                throw RINGMeshException( "I/O", "Failed to open file ", filename_ );
            }
        }

        /*!
         * @brief Parses the file and loads the GeoModel
         * @details The GeoModel loaded by this function is not valid because
         * some computation are still not done (i.e., surface internal borders,
         * lines and corners computation, boundary links between region and
         * surface, contacts)
         */
        void read_file();

    protected:
        virtual void read_line() = 0;

    protected:
		std::string file_name_;
        GEO::LineInput file_line_;
    };

    class GocadBaseParser{
    ringmesh_disable_copy(GocadBaseParser);
    protected:
        GocadBaseParser() = default;
        virtual ~GocadBaseParser() = default;

        GeoModelBuilderGocad& builder()
        {
            ringmesh_assert( builder_ != nullptr );
            return *builder_;
        }

        GeoModel3D& geomodel()
        {
            ringmesh_assert( geomodel_ != nullptr );
            return *geomodel_;
        }

        void set_builder( GeoModelBuilderGocad& builder )
        {
            builder_ = &builder;
        }

        void set_geomodel( GeoModel3D& geomodel )
        {
            geomodel_ = &geomodel;
        }

    private:
        GeoModelBuilderGocad* builder_ { nullptr };
        GeoModel3D* geomodel_ { nullptr };
    };

    struct GocadLoadingStorage {
        GocadLoadingStorage();

        /*!
         * @brief Ends a polygon (by adding the size of list of polygon corners at the
         * end of the vector)
         */
        void end_polygon()
        {
            index_t nb_polygon_corners =
                static_cast< index_t >( cur_surf_polygon_corners_gocad_id_.size() );
            cur_surf_polygon_ptr_.push_back( nb_polygon_corners );
        }

        // The orientation of positive Z
        int z_sign_{ 1 };

        std::vector< vec3 > vertices_;

        std::vector< std::vector< double > > attributes_;

        // The vertices and the atoms
        index_t nb_attribute_fields_ = 0;

        // Current interface index
        index_t cur_interface_ { NO_ID };

        // Current surface index
        index_t cur_surface_ { NO_ID };

        // List of polygon corners for the current surface (gocad indices)
        std::vector< index_t > cur_surf_polygon_corners_gocad_id_;

        // Starting indices (in cur_surf_polygons_corner_gocad_id_) of each
        // polygon of the current surface
        std::vector< index_t > cur_surf_polygon_ptr_;
    };

    class GocadLineParser: public GocadBaseParser {
    ringmesh_disable_copy(GocadLineParser);
    public:
        static std::unique_ptr< GocadLineParser > create(
            const std::string& keyword,
            GeoModelBuilderGocad& gm_builder,
            GeoModel3D& geomodel );
        virtual void execute(
            GEO::LineInput& line,
            GocadLoadingStorage& load_storage ) = 0;
    protected:
        GocadLineParser() = default;
    };

    using GocadLineParserFactory = GEO::Factory0< GocadLineParser >;
#define ringmesh_register_GocadLineParser_creator(type, name) \
     geo_register_creator(GocadLineParserFactory, type, name)

    /*!
     * @brief Structure which maps the vertex indices in Gocad::TSolid to the
     * pair (region, index in region) in the RINGMesh::GeoModel
     */
    struct VertexMap {
        index_t gocad_vertex_id( index_t id ) const
        {
            ringmesh_assert( id < vertices_gocad_id_.size() );
            return vertices_gocad_id_[id];
        }

        index_t region_id( index_t id ) const
        {
            ringmesh_assert( id < vertices_region_id_.size() );
            return vertices_region_id_[id];
        }

        void add_new_region( index_t region_id,
            std::string region_name )
        {
            region_ids_.push_back( region_id );
            region_names_.push_back( region_name );
        }

        void record_vertex_with_its_region( index_t gocad_vertex_id, index_t region_id )
        {
            vertices_gocad_id_.push_back( gocad_vertex_id );
            vertices_region_id_.push_back( region_id );
        }

        index_t nb_regions() const
        {
            ringmesh_assert( region_names_.size() == region_ids_.size() );
            return static_cast<index_t>( region_ids_.size() );
        }

        bool find_region_id_from_name( const std::string& region_name,
            index_t& region_id )
        {
            for( size_t i = 0; i < region_names_.size(); i++ ){
                if( region_name.compare( region_names_[i] ) == 0 ){
                    region_id = region_ids_[i];
                    ringmesh_assert( region_id != NO_ID );
                    return true;
                }
            }
            region_id = NO_ID;
            return false;
        }

        const std::vector< index_t >& get_regions()
        {
            return region_ids_;
        }

        void get_tetra_corners_with_this_region_id( index_t region_id,
            std::vector< index_t >& region_tetra_corners_local ){

            unsigned int counter = 0;
            for( index_t tetra_region_id : vertices_region_id_ ) {
                if( tetra_region_id == region_id ){
                    index_t local_i = local_id( vertices_gocad_id_[counter] );
                    region_tetra_corners_local.push_back( local_i );
                }
                counter++;
            }
        }

        void get_vertices_attributes_list_from_gocad_ids(
            const std::vector< std::vector< double > >& stored_attributes,
            index_t region_id,
            const std::map< index_t, index_t >& lighttsolid_atom_map,
            std::vector< std::vector< double > >& region_tetra_attributes ) const {

            index_t gocad_id = 0;
            for( std::vector< double > attrib : stored_attributes ){
                if( gocad_ids2region_ids_[gocad_id] == region_id ){
                    if( lighttsolid_atom_map.find( gocad_id ) == lighttsolid_atom_map.end() ){
                        region_tetra_attributes.push_back( stored_attributes[gocad_id] );
                    } else {
                        index_t corresponding_gocad_id = lighttsolid_atom_map.find( gocad_id )->second;
                        if( region( corresponding_gocad_id ) != region( gocad_id ) ){
                            region_tetra_attributes.push_back( stored_attributes[corresponding_gocad_id] );
                        }
                    }
                }
                gocad_id++;
            }
        }

        void get_vertices_list_and_local_ids_from_gocad_ids(
            const std::vector< vec3 >& stored_vertices,
            index_t region_id,
            const std::map< index_t, index_t >& lighttsolid_atom_map,
            std::vector< vec3 >& region_tetra_vertices,
            std::vector< index_t >& local_ids ){

            local_ids.clear();

            index_t gocad_id = 0;
            for( vec3 vertex : stored_vertices ){
                if( gocad_ids2region_ids_[gocad_id] == region_id ){
                    if( lighttsolid_atom_map.find( gocad_id ) == lighttsolid_atom_map.end() ){
                        region_tetra_vertices.push_back( stored_vertices[gocad_id] );
                        local_ids.push_back( gocad_id );
                    } else {
                        index_t corresponding_gocad_id = lighttsolid_atom_map.find( gocad_id )->second;
                        if( region( corresponding_gocad_id ) != region( gocad_id ) ){
                            region_tetra_vertices.push_back( stored_vertices[corresponding_gocad_id] );
                            local_ids.push_back( gocad_id );
                        }
                    }
                }
                gocad_id++;
            }
        }

        index_t local_id( index_t gocad_vertex_id ) const
        {
            ringmesh_assert( gocad_vertex_id < gocad_ids2local_ids_.size() );
            return gocad_ids2local_ids_[gocad_vertex_id];
        }

        index_t region( index_t gocad_vertex_id ) const
        {
            ringmesh_assert( gocad_vertex_id < gocad_ids2region_ids_.size() );
            return gocad_ids2region_ids_[gocad_vertex_id];
        }

        void add_vertex( index_t local_vertex_id, index_t region_id )
        {
            gocad_ids2local_ids_.push_back( local_vertex_id );
            gocad_ids2region_ids_.push_back( region_id );
        }

        void reserve( index_t capacity )
        {
            gocad_ids2local_ids_.reserve( capacity );
            gocad_ids2region_ids_.reserve( capacity );
        }

        void reserve_nb_vertices( index_t capacity )
        {
            vertices_gocad_id_.reserve( capacity );
            vertices_region_id_.reserve( capacity );
        }

        void fill_with_lighttsolid_region_ids()
        {
            ringmesh_assert( vertices_gocad_id_.size() == vertices_region_id_.size() );
            size_t lighttsolid_vertices_nb = vertices_gocad_id_.size();
            size_t lighttsolid_region_nb = nb_regions();

            // For every region ...
            for( index_t rgion_id = 0; rgion_id < lighttsolid_region_nb; rgion_id++ ) {
                // we want to record the region ids of the LightTSolid vertices ...
                for( index_t lighttsolid_vertices_id = 0;
                    lighttsolid_vertices_id < lighttsolid_vertices_nb; lighttsolid_vertices_id++ ) {
                    // that are in this region.
                    if( region_id( lighttsolid_vertices_id ) == rgion_id ) {
                        index_t gocad_vertex_i = gocad_vertex_id( lighttsolid_vertices_id );

                        gocad_ids2region_ids_[gocad_vertex_i]
                            = region_id( lighttsolid_vertices_id );
                    }
                }
            }
        }

        void fill_with_lighttsolid_local_ids()
        {
            ringmesh_assert( vertices_gocad_id_.size() == vertices_region_id_.size() );
            size_t lighttsolid_vertices_nb = vertices_gocad_id_.size();
            size_t lighttsolid_region_nb = nb_regions();

            // For every region ...
            for( index_t rgion_id = 0; rgion_id < lighttsolid_region_nb; rgion_id++ ) {
                // we want to record the local ids of the LightTSolid vertices ...
                std::vector< index_t > local_ids = local_ids_[rgion_id];
                for( index_t lighttsolid_vertices_id = 0;
                    lighttsolid_vertices_id < lighttsolid_vertices_nb; lighttsolid_vertices_id++ ) {
                    // that are in this region.
                    if( region_id( lighttsolid_vertices_id ) == rgion_id ) {

                        index_t gocad_vertex_i = gocad_vertex_id( lighttsolid_vertices_id );
                        for( index_t local_id = 0; local_id < local_ids.size(); local_id++ ){
                            index_t gocad_id = local_ids[local_id];
                            if( gocad_id == gocad_vertex_i ){
                                gocad_ids2local_ids_[gocad_vertex_i]
                                    = local_id;
                            }
                        }
                    }
                }
            }
        }

        void deal_with_same_region_atoms( std::map< index_t, index_t > lighttsolid_atom_map )
        {
            for( std::pair< index_t, index_t > pair : lighttsolid_atom_map ){
                if( region( pair.first ) == region( pair.second ) ){
                    gocad_ids2local_ids_[pair.first]
                        = gocad_ids2local_ids_[pair.second];
                }
            }
        }

    public:
        /*!
         * LightTSolids
         * One std::vector< index_t > per region. Regions vector are recorded in the same order as region_ids_.
         * Mapping the indices of vertices from Gocad .so file to the local (in region) indices of vertices.
         * No duplicates: Only the Vertex and the Atoms linking to a different region Vertex are recorded here.
         */
        std::vector< std::vector< index_t > > local_ids_;

    private:
        /*!
         * LightTSolids
         * The gocad IDs of the vertices of the LightTSolid as read in the Tetra lines of the file.
         * Duplicates: One vertex is recorded multiple times due to its presence in multiple tetras.
         */
        std::vector< index_t > vertices_gocad_id_;
        /*!
         * LightTSolids
         * The region IDs of the vertices of the LightTSolid as read in the # CTETRA lines of the file.
         * Duplicates: One vertex is recorded multiple times due to its presence in multiple tetras.
         */
        std::vector< index_t > vertices_region_id_;

        /*!
         * LightTSolids
         * All the region IDs are stored in order here. Ex: 0, 1, 2, 3.
         * No duplicates.
         */
        std::vector< index_t > region_ids_;
        /*!
         * LightTSolids
         * The names of the regions matching the regions IDs.
         * No duplicates.
         */
        std::vector< std::string > region_names_;

    private:
        /*!
         * TSolids & LightTSolids
         * Mapping the indices of vertices from Gocad .so file
         * to the local (in region) indices of vertices
         */
        std::vector< index_t > gocad_ids2local_ids_;
        /*!
         * TSolids & LightTSolids
         * Mapping the indices of vertices from Gocad .so file
         * to the region containing them
         */
        std::vector< index_t > gocad_ids2region_ids_;
    };

    /*!
     * @brief Structure used to load a GeoModel by GeoModelBuilderTSolid
     */
    struct TSolidLoadingStorage: public GocadLoadingStorage {
        // Current region index
        index_t cur_region_ { NO_ID };

        // Map between gocad and GeoModel vertex indices
        VertexMap vertex_map_;

        // Region tetrahedron corners
        std::vector< index_t > tetra_corners_;

        // Names of the attributes for the TSolid
        std::vector< std::string > vertex_attribute_names_;

        // Dimensions of the attributes for the TSolid
        std::vector< index_t > vertex_attribute_dims_;

        //// Current lighttsolid gocad vertex index 1
        index_t cur_gocad_vrtx_id1_{ NO_ID };
        index_t cur_gocad_vrtx_id2_{ NO_ID };
        index_t cur_gocad_vrtx_id3_{ NO_ID };
        index_t cur_gocad_vrtx_id4_{ NO_ID };

        //// LightTSolid map between atoms and vertex
        std::map< index_t, index_t > lighttsolid_atom_map_;

        // The vertices and the atoms
        index_t nb_vertices_;
    };

    class TSolidLineParser: public GocadBaseParser {
    ringmesh_disable_copy(TSolidLineParser);
    public:
        TSolidLineParser() = default;

        static std::unique_ptr< TSolidLineParser > create(
            const std::string& keyword,
            GeoModelBuilderTSolid& gm_builder,
            GeoModel3D& geomodel );
        virtual void execute(
            GEO::LineInput& line,
            TSolidLoadingStorage& load_storage ) = 0;
        virtual void execute_light(
            GEO::LineInput& line,
            TSolidLoadingStorage& load_storage ) = 0;

        /*!
        * @brief Creates an empty entity of type GeoModelEntity::REGION and sets
        * its name from .so file
        * @param[in] region_name Name of the new region
        * @param[in] geomodel_builder Builder of the geomodel
        * @return The index of the initialized region
        */
        index_t initialize_region(
            const std::string& region_name,
            GeoModelBuilderGocad& geomodel_builder )
        {
            gmme_id cur_region =
                geomodel_builder.topology.create_mesh_entity< Region >();
            geomodel_builder.info.set_mesh_entity_name( cur_region, region_name );
            return cur_region.index();
        }

        void assign_attributes_to_mesh(
            index_t region_id,
            const Region< 3 >& region, TSolidLoadingStorage& load_storage,
            const std::vector< std::vector< double > >& region_attributes ){

            for( index_t attrib_itr : range( load_storage.vertex_attribute_names_.size() ) ) {
                std::string name = load_storage.vertex_attribute_names_[attrib_itr];

                if( region.vertex_attribute_manager().is_defined( name ) ) {
                    Logger::warn( "Transfer attribute", "The attribute ", name,
                        " already exists on the ", region.gmme() );
                    continue;
                }
                GEO::Attribute< double > attr;
                index_t nb_dimensions = load_storage.vertex_attribute_dims_[attrib_itr];
                attr.create_vector_attribute( region.vertex_attribute_manager(),
                    load_storage.vertex_attribute_names_[attrib_itr], nb_dimensions );
                // Does it resize all the past attributes to the size of the current attribute? 
                // Problematic, isn't it?
                region.vertex_attribute_manager().resize(
                    region_attributes.size() * nb_dimensions + nb_dimensions );
                for( index_t v_itr : range( region_attributes.size() ) ) {
                    for( index_t attrib_dim_itr : range( nb_dimensions ) ) {
                        attr[v_itr * nb_dimensions + attrib_dim_itr] =
                            region_attributes[v_itr][attrib_dim_itr];
                    }
                }
            }
        }
    };

    using TSolidLineParserFactory = GEO::Factory0< TSolidLineParser >;
#define ringmesh_register_TSolidLineParser_creator(type, name) \
     geo_register_creator(TSolidLineParserFactory, type, name)

    /*!
     * @brief Builds a meshed GeoModel from a Gocad TSolid (file.so)
     */
    class RINGMESH_API GeoModelBuilderTSolid final : public GeoModelBuilderGocad{
    public:
        static const index_t NB_TYPE = 2;
        GeoModelBuilderTSolid( GeoModel3D& geomodel, std::string filename );
        virtual ~GeoModelBuilderTSolid() = default;

    private:
        void read_number_of_vertices();
        void read_type();
        void load_file() final;
        void read_file();

        /*!
         * @brief Reads the first word of the current line (keyword)
         * and executes the good action with the information of the line
         * @details Uses the TsolidLineParser factory
         */
        void read_line() final;

        /*!
         * @brief Computes internal borders of a given surface
         * @details A surface polygon edge is an internal border if it is shared
         * by at least two surfaces. Adjacency of such a polygon edge is set to
         * GEO::NO_FACET.
         * @param[in] geomodel GeoModel to consider
         * @param[in] surface_id Index of the surface
         * @param[in] surface_nns Unique pointers to the NNSearchs of surfaces
         */
        void compute_surface_internal_borders(
            index_t surface_id,
            const std::vector< std::unique_ptr< NNSearch3D > >& surface_nns,
            const std::vector< Box3D >& surface_boxes );

        /*!
         * @brief Computes the NNSearchs of the centers of polygon edges for
         * each surface and their Box3d
         * @param[in] geomodel GeoModel to consider
         * @param[out] surface_nns Unique pointers to the NNSearchs of surfaces
         * @param[out] surface_boxes Bounding Box of surfaces
         */
        void compute_polygon_edge_centers_nn_and_surface_boxes(
            std::vector< std::unique_ptr< NNSearch3D > >& surface_nns,
            std::vector< Box3D >& surface_boxes ) const;

        /*!
         * @brief Computes internal borders of the geomodel surfaces
         * @details An surface polygon edge is an internal border if it is shared
         * by at least two surfaces. Adjacency of such a polygon edge is set to
         * GEO::NO_FACET.
         * @param[in] geomodel GeoModel to consider
         */
        void compute_surfaces_internal_borders();

    protected:
        TSolidLoadingStorage tsolid_load_storage_;

    private:
        index_t file_type_{ 0 };
        std::unique_ptr< GeoModelBuilderTSolidImpl > type_impl_[NB_TYPE];

        friend class RINGMesh::GocadLineParser;
    };

    class GeoModelBuilderTSolidImpl {
    public:
        GeoModelBuilderTSolidImpl(
            GeoModelBuilderTSolid& builder,
            GeoModel< 3 >& geomodel,
            GEO::LineInput& file_line,
            TSolidLoadingStorage& tsolid_load_storage )
            : builder_( builder ), geomodel_( geomodel ),
            file_line_( file_line ), tsolid_load_storage_( tsolid_load_storage )
        {
        }
        virtual ~GeoModelBuilderTSolidImpl() = default;

        virtual void read_line() = 0;

    protected:
        TSolidLoadingStorage& tsolid_load_storage_;
        GEO::LineInput& file_line_;
        GeoModelBuilderTSolid& builder_;
        GeoModel< 3 >& geomodel_;
    };

    class GeoModelBuilderTSolidImpl_TSolid final: public GeoModelBuilderTSolidImpl {
    public:
        GeoModelBuilderTSolidImpl_TSolid(
            GeoModelBuilderTSolid& builder,
            GeoModel< 3 >& geomodel,
            GEO::LineInput& file_line,
            TSolidLoadingStorage& tsolid_load_storage )
            : GeoModelBuilderTSolidImpl( builder, geomodel, file_line, tsolid_load_storage )
        {
        }
        virtual ~GeoModelBuilderTSolidImpl_TSolid() = default;

        virtual void read_line() override
        {
            std::string keyword = file_line_.field( 0 );
            std::unique_ptr< TSolidLineParser > tsolid_parser = TSolidLineParser::create(
                keyword, this->builder_, geomodel_ );
            if( tsolid_parser ) {
                tsolid_parser->execute( file_line_, tsolid_load_storage_ );
            } else {
                std::unique_ptr< GocadLineParser > gocad_parser =
                    GocadLineParser::create( keyword, this->builder_, geomodel_ );
                if( gocad_parser ) {
                    gocad_parser->execute( file_line_, tsolid_load_storage_ );
                }
            }
        }
    };

    class GeoModelBuilderTSolidImpl_LightTSolid final: public GeoModelBuilderTSolidImpl {
    public:
        GeoModelBuilderTSolidImpl_LightTSolid(
            GeoModelBuilderTSolid& builder,
            GeoModel< 3 >& geomodel,
            GEO::LineInput& file_line,
            TSolidLoadingStorage& tsolid_load_storage )
            : GeoModelBuilderTSolidImpl( builder, geomodel, file_line, tsolid_load_storage )
        {
        }
        virtual ~GeoModelBuilderTSolidImpl_LightTSolid() = default;

        virtual void read_line() override
        {
            std::string keyword = file_line_.field( 0 );
            std::unique_ptr< TSolidLineParser > tsolid_parser = TSolidLineParser::create(
                keyword, this->builder_, geomodel_ );
            if( tsolid_parser ) {
                tsolid_parser->execute_light( file_line_, tsolid_load_storage_ );
            } else {
                std::unique_ptr< GocadLineParser > gocad_parser =
                    GocadLineParser::create( keyword, this->builder_, geomodel_ );
                if( gocad_parser ) {
                    gocad_parser->execute( file_line_, tsolid_load_storage_ );
                }
            }
        }
    };

    struct MLLoadingStorage: public GocadLoadingStorage {
        MLLoadingStorage();

        bool is_header_read_ { false };

        /// Offset to read in the tface vertices in the tsurf vertices
        index_t tface_vertex_ptr_ { 0 };
    };
    class MLLineParser: public GocadBaseParser {
    ringmesh_disable_copy(MLLineParser);
    public:
        MLLineParser() = default;

        static std::unique_ptr< MLLineParser > create(
            const std::string& keyword,
            GeoModelBuilderML& gm_builder,
            GeoModel3D& geomodel );
        virtual void execute(
            GEO::LineInput& line,
            MLLoadingStorage& load_storage ) = 0;
    };

    using MLLineParserFactory = GEO::Factory0< MLLineParser >;
#define ringmesh_register_MLLineParser_creator(type, name) \
     geo_register_creator(MLLineParserFactory, type, name)

    /*!
     * @brief Build a GeoModel from a Gocad Model3D (file_model.ml)
     */
    class RINGMESH_API GeoModelBuilderML final : public GeoModelBuilderGocad {
    public:
        GeoModelBuilderML( GeoModel3D& geomodel, std::string filename )
            : GeoModelBuilderGocad( geomodel, std::move( filename ) )
        {
        }
        virtual ~GeoModelBuilderML() = default;

    private:
        /*!
         * @brief Loads and builds a GeoModel from a Gocad .ml file
         * @warning Pretty unstable. Crashes if the file is not exactly what is expected.
         * @details Correspondance between Gocad::Model3D entities
         * and GeoModel entities is :
         *  - Gocad TSurf  <-> GeoModel Interface
         *  - Gocad TFace  <-> GeoModel Surface
         *  - Gocad Region <-> GeoModel Region
         *  - Gocad Layer  <-> GeoModel Layer
         * @param[in] ml_file_name Input .ml file stream
         * @param[in] ignore_file_borders If true, BORDER and BSTONE entries in the files
         * are ignored and the Lines and Corners of the GeoModel are deduced from the
         * connectivity of its Surfaces. By default set to false.
         */
        void load_file() final;

        /*!
         * @brief Reads the first word of the current line (keyword)
         * and executes the good action with the information of the line
         * @details Uses the MLLineParser factory
         */
        void read_line() final;

    private:
        MLLoadingStorage ml_load_storage_;
    };

}
