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
    class GeoModelBuilderML;
    template< index_t DIMENSION > class Box;
    struct VertexMap;
    struct TSolidLoadingStorage;
}

namespace RINGMesh {

    void RINGMESH_API initialize_gocad_import_factories();

    class RINGMESH_API GeoModelBuilderGocad: public GeoModelBuilderFile< 3 > {
    public:
        GeoModelBuilderGocad( GeoModel< 3 >& geomodel, std::string filename )
            :
                GeoModelBuilderFile( geomodel, std::move( filename ) ),
                file_line_( filename_ )
        {
            /*! @todo Review: A constructor is not supposed to throw, the object is left in an
             * undefined state [JP] */
            if( !file_line_.OK() ) {
                throw RINGMeshException( "I/O", "Failed to open file " + filename_ );
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
        GEO::LineInput file_line_;
    };

    class GocadBaseParser: public GEO::Counted {
    ringmesh_disable_copy(GocadBaseParser);
    protected:
        GocadBaseParser()
            : Counted(), builder_( nullptr ), geomodel_( nullptr )
        {
        }
        virtual ~GocadBaseParser() = default;

        GeoModelBuilderGocad& builder()
        {
            ringmesh_assert( builder_ != nullptr );
            return *builder_;
        }

        GeoModel< 3 >& geomodel()
        {
            ringmesh_assert( geomodel_ != nullptr );
            return *geomodel_;
        }

        void set_builder( GeoModelBuilderGocad& builder )
        {
            builder_ = &builder;
        }

        void set_geomodel( GeoModel< 3 >& geomodel )
        {
            geomodel_ = &geomodel;
        }

    private:
        GeoModelBuilderGocad* builder_;
        GeoModel< 3 >* geomodel_;
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
        int z_sign_;

        std::vector< vec3 > vertices_;

        // Current interface index
        index_t cur_interface_;

        // Current surface index
        index_t cur_surface_;

        // List of polygon corners for the current surface (gocad indices)
        std::vector< index_t > cur_surf_polygon_corners_gocad_id_;

        // Starting indices (in cur_surf_polygons_corner_gocad_id_) of each
        // polygon of the current surface
        std::vector< index_t > cur_surf_polygon_ptr_;
    };

    class GocadLineParser: public GocadBaseParser {
    ringmesh_disable_copy(GocadLineParser);
    public:
        GocadLineParser() = default;
        static std::unique_ptr< GocadLineParser > create(
            const std::string& keyword,
            GeoModelBuilderGocad& gm_builder,
            GeoModel< 3 >& geomodel );
        virtual void execute(
            GEO::LineInput& line,
            GocadLoadingStorage& load_storage ) = 0;
    };

    using GocadLineParserFactory = GEO::Factory0< GocadLineParser >;
#define ringmesh_register_GocadLineParser_creator(type, name) \
     geo_register_creator(GocadLineParserFactory, type, name)

    /*!
     * @brief Structure which maps the vertex indices in Gocad::TSolid to the
     * pair (region, index in region) in the RINGMesh::GeoModel
     */
    struct VertexMap {
        VertexMap() = default;

        index_t local_id( index_t gocad_vertex_id ) const
        {
            return gocad_vertices2region_vertices_[gocad_vertex_id];
        }

        index_t region( index_t gocad_vertex_id ) const
        {
            return gocad_vertices2region_id_[gocad_vertex_id];
        }

        void add_vertex( index_t local_vertex_id, index_t region_id )
        {
            gocad_vertices2region_vertices_.push_back( local_vertex_id );
            gocad_vertices2region_id_.push_back( region_id );
        }

        index_t nb_vertex() const
        {
            ringmesh_assert(
                gocad_vertices2region_vertices_.size()
                    == gocad_vertices2region_id_.size() );
            return static_cast< index_t >( gocad_vertices2region_vertices_.size() );
        }

        void reserve( index_t capacity )
        {
            gocad_vertices2region_vertices_.reserve( capacity );
            gocad_vertices2region_id_.reserve( capacity );
        }

    private:
        /*!
         * Mapping the indices of vertices from Gocad .so file
         * to the local (in region) indices of vertices
         */
        std::vector< index_t > gocad_vertices2region_vertices_;
        /*!
         * Mapping the indices of vertices from Gocad .so file
         * to the region containing them
         */
        std::vector< index_t > gocad_vertices2region_id_;
    };

    /*!
     * @brief Structure used to load a GeoModel by GeoModelBuilderTSolid
     */
    struct TSolidLoadingStorage: public GocadLoadingStorage {
        TSolidLoadingStorage();

        // Current region index
        index_t cur_region_;

        // Map between gocad and GeoModel vertex indices
        VertexMap vertex_map_;

        // Region tetrahedron corners
        std::vector< index_t > tetra_corners_;

    };
    class TSolidLineParser: public GocadBaseParser {
    ringmesh_disable_copy(TSolidLineParser);
    public:
        TSolidLineParser() = default;

        static std::unique_ptr< TSolidLineParser > create(
            const std::string& keyword,
            GeoModelBuilderTSolid& gm_builder,
            GeoModel< 3 >& geomodel );
        virtual void execute(
            GEO::LineInput& line,
            TSolidLoadingStorage& load_storage ) = 0;
    };

    using TSolidLineParserFactory = GEO::Factory0< TSolidLineParser >;
#define ringmesh_register_TSolidLineParser_creator(type, name) \
     geo_register_creator(TSolidLineParserFactory, type, name)

    /*!
     * @brief Builds a meshed GeoModel from a Gocad TSolid (file.so)
     */
    class RINGMESH_API GeoModelBuilderTSolid final : public GeoModelBuilderGocad {
    public:
        GeoModelBuilderTSolid( GeoModel< 3 >& geomodel, std::string filename )
            : GeoModelBuilderGocad( geomodel, std::move( filename ) )
        {
        }
        virtual ~GeoModelBuilderTSolid() = default;

    private:
        virtual void load_file() final;

        /*!
         * @brief Reads the first word of the current line (keyword)
         * and executes the good action with the information of the line
         * @details Uses the TsolidLineParser factory
         */
        virtual void read_line() final;

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
            const std::vector< std::unique_ptr< NNSearch< 3 > > >& surface_nns,
            const std::vector< Box< 3 > >& surface_boxes );

        /*!
         * @brief Computes the NNSearchs of the centers of polygon edges for
         * each surface and their Box3d
         * @param[in] geomodel GeoModel to consider
         * @param[out] surface_nns Unique pointers to the NNSearchs of surfaces
         * @param[out] surface_boxes Bounding Box of surfaces
         */
        void compute_polygon_edge_centers_nn_and_surface_boxes(
            std::vector< std::unique_ptr< NNSearch< 3 > > >& surface_nns,
            std::vector< Box< 3 > >& surface_boxes );

        /*!
         * @brief Computes internal borders of the geomodel surfaces
         * @details An surface polygon edge is an internal border if it is shared
         * by at least two surfaces. Adjacency of such a polygon edge is set to
         * GEO::NO_FACET.
         * @param[in] geomodel GeoModel to consider
         */
        void compute_surfaces_internal_borders();

    private:
        TSolidLoadingStorage tsolid_load_storage_;
        friend class RINGMesh::GocadLineParser;
    };

    struct MLLoadingStorage: public GocadLoadingStorage {
        MLLoadingStorage();

        bool is_header_read_;

        /// Offset to read in the tface vertices in the tsurf vertices
        index_t tface_vertex_ptr_;
    };
    class MLLineParser: public GocadBaseParser {
    ringmesh_disable_copy(MLLineParser);
    public:
        MLLineParser() = default;

        static std::unique_ptr< MLLineParser > create(
            const std::string& keyword,
            GeoModelBuilderML& gm_builder,
            GeoModel< 3 >& geomodel );
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
        GeoModelBuilderML( GeoModel< 3 >& geomodel, std::string filename )
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
        virtual void load_file() final;

        /*!
         * @brief Reads the first word of the current line (keyword)
         * and executes the good action with the information of the line
         * @details Uses the MLLineParser factory
         */
        virtual void read_line() final;

    private:
        MLLoadingStorage ml_load_storage_;
    };

}
