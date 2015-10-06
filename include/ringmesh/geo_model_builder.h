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
 *     Ecole Nationale Superieure de Geologie - Georessources
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

namespace GEO {
    class Mesh ;
}

namespace RINGMesh {
    /*!
     * @brief Base class for all classes building a GeoModel.
     * @details Derive from this class to build or modify a GeoModel
     * 
     * NON Geometry related function have been moved in the GeoModelEditor
     * 
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
      
        bool end_model() ;

        bool complete_element_connectivity() ;
    
        /*! @}
         * \name Set element geometry from geometrical positions   
         * @{
         */
        void set_element_vertex(
            GME::gme_t t,
            index_t v,
            const vec3& point )
        {
            mesh_element( t ).set_vertex( v, point, false ) ;
        }

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

        void set_surface_adjacencies( const GME::gme_t& surface_id ) ;

        /*!
         * @}
         */       
    } ;

    /*!
     * @brief Build a GeoModel from a Gocad Model3D (file_model.ml)
     */
    class RINGMESH_API GeoModelBuilderGocad: public GeoModelBuilder {
    public:
        GeoModelBuilderGocad( GeoModel& model )
            : GeoModelBuilder( model )
        {
        }
        virtual ~GeoModelBuilderGocad()
        {
        }

        bool load_ml_file( const std::string& ml_file_name ) ;

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
            {
            }

        public:
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
    class RINGMESH_API GeoModelBuilderBM: public GeoModelBuilder {
    public:
        GeoModelBuilderBM( GeoModel& model )
            : GeoModelBuilder( model )
        {
        }
        virtual ~GeoModelBuilderBM()
        {
        }

        bool load_file( const std::string& bm_file_name ) ;

    private:
        static GME::TYPE match_nb_elements( const char* s ) ;

        static GME::TYPE match_type( const char* s ) ;

        static bool match_high_level_type( const char* s )
        {
            return GME::child_allowed( match_type( s ) ) ;
        }
    } ;

    /*!
     * @brief Builder of a GeoModel from a surface mesh
     *        in which the manifold surface connected components are disjoints
     */
    class RINGMESH_API GeoModelBuilderSurface: public GeoModelBuilder {
    public:
        GeoModelBuilderSurface( GeoModel& model )
            : GeoModelBuilder( model )
        {
        }
        virtual ~GeoModelBuilderSurface()
        {
        }

        void set_surfaces( const GEO::Mesh& mesh ) ;

        bool build_model( bool build_regions = true ) ;
    } ;

}

#endif
