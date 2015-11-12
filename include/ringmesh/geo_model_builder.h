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


#ifndef __RINGMESH_GEO_MODEL_BUILDER__
#define __RINGMESH_GEO_MODEL_BUILDER__

#include <vector>
#include <string>
#include <stack>

#include <geogram/basic/line_stream.h>

#include <ringmesh/common.h>
#include <ringmesh/geo_model.h>
#include <ringmesh/geo_model_editor.h>



/*!
* @file ringmesh/geo_model_builder.h
* @brief Classes to build GeoModel from various inputs
* @author Jeanne Pellerin
*/

namespace GEO {
    class Mesh ;
}

namespace RINGMesh {

    struct GeoModelIOFlags {
        bool compute_corners ;
        bool compute_lines ;
        bool compute_surfaces ;
        bool compute_regions_brep ;
        bool compute_regions_mesh ;
        // What about other elements ? As they are not mandatory
        // We will see later [JP]
    };

    /*! @todo We need to keep track of the status of the GeoModel
     *  and to know what elements are built, is their topolgy set,
     *  are their mesh assigned?
     */


    // Internal implementation class
    class GeoModelRegionFromSurfaces ;

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
        virtual ~GeoModelBuilder() ;

        bool set_options( GeoModelIOFlags options )
        {
            /*! @todo Check that the options are consistent */
            options_ = options ;
            return true ;
        }

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

        void set_surface_geometry(
            const GME::gme_t& surface_id,
            const std::vector< index_t >& triangle_corners ) ;

        index_t find_or_create_duplicate_vertex(
            GeoModelMeshElement& S,
            index_t model_vertex_id,
            index_t surface_vertex_id ) ;

        void cut_surface_by_line( Surface& S, const Line& L ) ;

        void compute_surface_adjacencies( const GME::gme_t& surface_id ) ;

        /*!
         * @}
         * \name Model building functions
         */   

        /*!
        * @brief From the Surfaces of the GeoModel, build its Lines and Corners
        */
        bool build_lines_and_corners_from_surfaces() ;

        /*!
        * @brief Build the regions of the GeoModel from the Surfaces 
        * @details Call build_lines_and_corners_from_surfaces first
        */
        bool build_brep_regions_from_surfaces() ;
    
        /* 
         * @brief From a GeoModel in which only Surface are defined, create
         * corners, contacts. If the compute_regions_brep is true, then regions
         * too are computed.
         * @return True if a model has been built.
         * @note Valdity is not checked
         * @pre The GeoModel should have at least one Surface. Nothing is done if not.
         */
        bool build_model_from_surfaces() ;

        /*!
        * @brief Finish up model building and complete missing information.
        * @return True except if the model has no Surface
        */
        bool end_model() ;

    protected:
        /*! Elements to compute from the available elements */
        GeoModelIOFlags options_ ; 

        /*! Internal information filled at Line building step */
        std::vector< GeoModelRegionFromSurfaces* > regions_info_ ;

    private:
        void assign_surface_mesh_facets(
            const GME::gme_t& surface_id,
            const std::vector< index_t >& facets,
            const std::vector< index_t >& facet_ptr ) ;
    } ;

    // Implementation
    class GeoModelSurfaceFromMesh ;
    class GeoModelRegionFromMesh ;
    
    /*!
     * @brief Builder of a GeoModel from a Mesh 
     * @details It either fills existing Regions and Surfaces
     * or creates and fills them.
     * @warning Volumetric version is only implemented for a tetrahedral Mesh
     */
    class RINGMESH_API GeoModelBuilderMesh: public GeoModelBuilder {
    public:
        GeoModelBuilderMesh( GeoModel& model, const GEO::Mesh& mesh )
            : GeoModelBuilder(model), 
              mesh_( mesh ),
              surface_builder_(nil),
              region_builder_(nil)
        {};
    
        virtual ~GeoModelBuilderMesh() ;

       /*!
        * @brief Create and fill the model Surfaces 
        * from the surface connected components of the Mesh
        */
        bool build_surfaces_from_connected_components() ;

        /*! 
         * @brief Create and fill the Surfaces 
         *  from a Integer attribute on the Mesh facets
         */
        bool build_surfaces_from_attribute_value( 
            const std::string& facet_attribute_name ) ;

        /*! 
         * @brief Fill the Surface meshes
         * from an Integer facet attribute giving the Surface index
         * of each Mesh facet 
         */
        bool fill_surface_meshes_from_attribute_value(
            const std::string& facet_attribute_name ) ;

        
        bool build_regions_from_connected_components() ;

        bool build_regions_from_attribute_value(
            const std::string& region_attribute_name ) ;
        
        /*! 
        * The given mesh is volumetric to fill the region of the Model
        * There is an attribute "region" that flag the tets region per region
        */
        bool fill_region_meshes_from_attribute_value( 
            const std::string& region_attribute_name ) ;

        /*! @}
        * \name Copy attributes from the Mesh to the GeoModel
        * @{
        */
        template< class T >
        void copy_facet_attribute_from_mesh( 
            const std::string& attribute_name ) ;
        
        template< class T >
        void copy_cell_attribute_from_mesh(
            const std::string& attribute_name ) ;

    protected:
        /*!
         * @brief Check that the mesh is valid to create a GeoModel
         * @details It shall have no duplicated vertices, and no duplicated facets.
         */
        bool has_mesh_colocated_vertices() ;

        /*!
         * @brief All vertices of the Mesh are added to the GeoModelMeshVertices 
         */
        void add_mesh_vertices_to_model() ;

        void initialize_surface_builder( const std::string& attribute_name ) ;
        void initialize_region_builder( const std::string& attribute_name ) ;

    protected:
        const GEO::Mesh& mesh_ ;
        GeoModelSurfaceFromMesh* surface_builder_ ;
        GeoModelRegionFromMesh* region_builder_ ;
    } ;


    class RINGMESH_API GeoModelBuilderFile : public GeoModelBuilder {
    public:
        GeoModelBuilderFile( GeoModel& model, const std::string& filename ) ;
        
        virtual ~GeoModelBuilderFile()
        {
        }
        virtual bool load_file() = 0 ;
        virtual bool build_model()
        {
            if( load_file() ) {
                return end_model() ;
            } else { 
                return false ; 
            }
        }
        /*! @todo Implement function to read the lines of the 
         *        file and wrap the GEO::LineInput which is not that easy to use 
         */
    protected:
        GEO::LineInput in_ ;
    };

    /*!
    * @brief Build a GeoModel from a Gocad Model3D (file_model.ml)
    */
    class RINGMESH_API GeoModelBuilderGocad : public GeoModelBuilderFile {
    public:
        GeoModelBuilderGocad( GeoModel& model, const std::string& filename )
            : GeoModelBuilderFile( model, filename )
        {}
        virtual ~GeoModelBuilderGocad()
        {}
        bool load_file() ;

    private:
        void build_contacts() ;

        GME::gme_t determine_line_vertices(
            const Surface& S,
            index_t id0,
            index_t id1,
            std::vector< vec3 >& border_vertex_model_ids ) const ;

        void create_surface(
            const std::string& interface_name,
            const std::string& type,
            const vec3& p0,
            const vec3& p1,
            const vec3& p2 ) ;
        
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

        std::vector< KeyFacet > key_facets_ ;
    } ;

    /*!
    * @brief Build a GeoModel from a file_model.bm
    */
    class RINGMESH_API GeoModelBuilderBM : public GeoModelBuilderFile {
    public:
        GeoModelBuilderBM( GeoModel& model, const std::string& filename )
            : GeoModelBuilderFile( model, filename )
        {}
        virtual ~GeoModelBuilderBM()
        {}

        bool load_file() ;

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
