/*
 * Copyright (c) 2012-2015, Association Scientifique pour la Geologie et ses Applications (ASGA)
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *     * Redistributions in binary form must reproduce the above copyright
 *       notice, this list of conditions and the following disclaimer in the
 *       documentation and/or other materials provided with the distribution.
 *     * Neither the name of the <organization> nor the
 *       names of its contributors may be used to endorse or promote products
 *       derived from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL <COPYRIGHT HOLDER> BE LIABLE FOR ANY
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 *  Contacts:
 *     Arnaud.Botella@univ-lorraine.fr
 *     Antoine.Mazuyer@univ-lorraine.fr
 *     Jeanne.Pellerin@wias-berlin.de
 *
 *     http://www.ring-team.org
 *
 *     RING Project
 *     Ecole Nationale Superieure de Geologie - GeoRessources
 *     2 Rue du Doyen Marcel Roubault - TSA 70605
 *     54518 VANDOEUVRE-LES-NANCY
 *     FRANCE
 */

/*! \author Jeanne Pellerin */

#ifndef __RINGMESH_GEO_MODEL_BUILDER__
#define __RINGMESH_GEO_MODEL_BUILDER__

#include <ringmesh/common.h>
#include <ringmesh/geo_model.h>
#include <ringmesh/geo_model_editor.h>

#include <vector>
#include <string>
#include <stack>

/*!
* @file ringmesh/geo_model_builder.h
* @brief Classes to build GeoModel from various inputs
* @author Jeanne Pellerin
*/

namespace GEO {
    class Mesh ;
}

namespace RINGMesh {
    /*!
     * @brief Base class for all classes building a GeoModel.
     * @details Derive from this class to build or modify a GeoModel. 
     * @note NON Geometry related functions are in GeoModelEditor class.
     */
    class RINGMESH_API GeoModelBuilder : public GeoModelEditor {
    public:
        GeoModelBuilder( GeoModel& model )
            : GeoModelEditor( model )
        {
        }
        virtual ~GeoModelBuilder()
        {
        }

        /*!
        * @brief Finish up model building, complete missing information
        * and check model correctness.
        */
        bool end_model() ;       
    
        /*!
         * \name Set element geometry from geometrical positions   
         * @{
         */
        void set_element_vertex(
            GME::gme_t t,
            index_t v,
            const vec3& point,
            bool update ) ;

        void set_element_vertices(
            const GME::gme_t& id,
            const std::vector< vec3 >& points,
            bool clear ) ;

        void set_corner(
            const GME::gme_t& corner_id,
            const vec3& point ) ;

        void set_line(
            const GME::gme_t& id,
            const std::vector< vec3 >& vertices ) ;

        void set_surface_geometry(
            const GME::gme_t& surface_id,
            const std::vector< vec3 >& surface_vertices,
            const std::vector< index_t >& surface_facets,
            const std::vector< index_t >& surface_facet_ptr ) ;

        /*! @}
         * \name Set element geometry using GeoModel vertices
         * @{
         */
        void set_element_vertex( const GME::gme_t& id, index_t v, index_t model_vertex ) ;

        void set_element_vertices(
            const GME::gme_t& id,
            const std::vector< index_t >& model_vertices,
            bool clear ) ;

        index_t add_unique_vertex( const vec3& p ) ;

        void set_corner(
            const GME::gme_t& corner_id,
            index_t unique_vertex ) ;

        void set_line(
            const GME::gme_t& id,
            const std::vector< index_t >& unique_vertices ) ;

        void set_surface_geometry(
            const GME::gme_t& surface_id,
            const std::vector< index_t >& surface_vertices,
            const std::vector< index_t >& surface_facets,
            const std::vector< index_t >& surface_facet_ptr ) ;

        void set_surface_geometry(
            const GME::gme_t& surface_id,
            const std::vector< index_t >& corners,
            const std::vector< index_t >& facet_ptr ) ;

        /*!
         * @}
         */   

        index_t find_or_create_duplicate_vertex(
            GeoModelMeshElement& S,
            index_t model_vertex_id,
            index_t surface_vertex_id ) ;

        void cut_surface_by_line( Surface& S, const Line& L ) ;

        void compute_surface_adjacencies( const GME::gme_t& surface_id ) ;


    private:
        void create_surface_geometry(
            const GME::gme_t& surface_id,
            const std::vector< index_t >& facets,
            const std::vector< index_t >& facet_ptr ) ;
    } ;

    
    // Forward declaration of a class used to 
    // implement GeoModelBuilderSurface
    class RegionBuildingInformation ;

    /*!
     * @brief Builder of a GeoModel from its Surfaces without 
     *        information on its Regions or Lines 
     * @note Manifold surface connected components are supposed disjoints
     *       It would be possible if needed to implement this for 
     *       non disjoint Surfaces using non-manifold edge identification.
     */
    class RINGMESH_API GeoModelBuilderSurface: public GeoModelBuilder {
    public:
        GeoModelBuilderSurface( GeoModel& model, bool build_regions = true )
            : GeoModelBuilder( model ),
            build_regions_( build_regions )
        {
        }
        virtual ~GeoModelBuilderSurface()
        {
            for( index_t i = 0; i < regions_info_.size(); ++i ) {
                delete regions_info_[ i ] ;
            }
        }
        /*!
        * @brief Create the model Surfaces from the connected components
        *       of the input surface mesh
        * @pre The input mesh is a surface mesh. Facet adjacencies are available.
        */
        bool set_surfaces( const GEO::Mesh& mesh ) ;
       
        /*!
        * @brief From a GeoModel in which only Surface are defined, create
        * corners, contacts and optionally regions.   
        * @return True if a valid model has been built, else returns false.
        * @pre The GeoModel should have at least one Surface. Nothing is done if not.
        */
        bool build_model() ;

    protected:
        /*!
        * @brief From the topology of the Surface of the GeoModel, build
        * its Lines and Corners
        */
        bool build_lines_and_corners() ;

        /*!
        * @brief Build the regions of the GeoModel from information collected
        * at Line building step.
        */
        bool build_regions() ;

    protected:
        /*! Build or not the Regions of the GeoModel from the the Surfaces 
         * and Lines */ 
        bool build_regions_ ;
        /*! Internal information filled at Line building step */
        std::vector< RegionBuildingInformation* > regions_info_ ;
    } ;


    /*!
    * @brief Build a GeoModel from a Gocad Model3D (file_model.ml)
    */
    class RINGMESH_API GeoModelBuilderGocad : public GeoModelBuilderSurface {
    public:
        GeoModelBuilderGocad( GeoModel& model )
            : GeoModelBuilderSurface( model, false )
        {}
        virtual ~GeoModelBuilderGocad()
        {}

        /*!
        * @brief Load and build a GeoModel from a Gocad .ml file
        * @warning Pretty unstable. Crashes if the file is not exactly as expected.
        * @details The correspondance between Gocad::Model3D elements 
        * and GeoModel elements is :
        *  - Gocad TSurf  <-> GeoModel Interface
        *  - Gocad TFace  <-> GeoModel Surface
        *  - Gocad Region <-> GeoModel Region
        *  - Gocad Layer  <-> GeoModel Layer
        *
        * @param[in] ml_file_name Input .ml file stream
        * @param[in] ignore_file_borders If true, BORDER and BSTONE entries in the files
        * are ignored and the Lines and Corners of the GeoModel are deduced from the 
        * connectivity of its Surfaces. By default set to false.
        */
        bool load_ml_file(
            const std::string& ml_file_name, 
            bool ignore_file_borders = false ) ;

    protected:
        GME::gme_t determine_line_vertices(
            const Surface& S,
            index_t id0,
            index_t id1,
            std::vector< index_t >& border_vertex_model_ids ) const ;

        GME::gme_t determine_line_vertices(
            const Surface& S,
            index_t first_vertex,
            index_t second_vertex,
            std::vector< vec3 >& border_vertex_model_vertices ) const ;

    private:
        void build_contacts() ;

        void create_surface(
            const std::string& interface_name,
            const std::string& type,
            const vec3& p0,
            const vec3& p1,
            const vec3& p2 ) ;

        /*!
        * @brief Triangle that set the orientation of a TFACE
        *        in a .ml file
        */
        struct KeyFacet {
            KeyFacet( const vec3& p0, const vec3& p1, const vec3& p2 )
                : p0_( p0 ), p1_( p1 ), p2_( p2 )
            {}
            vec3 p0_ ;
            vec3 p1_ ;
            vec3 p2_ ;
        } ;

        /*!
        * @brief Check if the surface triangle orientations match the one of the key facet
        */
        bool check_key_facet_orientation( index_t surface ) const ;

        index_t find_key_facet(
            index_t surface_id,
            const vec3& p0,
            const vec3& p1,
            const vec3& p2,
            bool& same_orientation ) const ;

    private:
        std::vector< KeyFacet > key_facets_ ;
    } ;

    /*!
    * @brief Build a GeoModel from a file_model.bm
    */
    class RINGMESH_API GeoModelBuilderBM : public GeoModelBuilder {
    public:
        GeoModelBuilderBM( GeoModel& model )
            : GeoModelBuilder( model )
        {}
        virtual ~GeoModelBuilderBM()
        {}

        bool load_file( const std::string& bm_file_name ) ;

    private:
        static GME::TYPE match_nb_elements( const char* s ) ;

        static GME::TYPE match_type( const char* s ) ;

        static bool match_high_level_type( const char* s )
        {
            return GME::child_allowed( match_type( s ) ) ;
        }
    } ;
}

#endif
