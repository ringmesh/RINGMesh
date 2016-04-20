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

#ifndef __RINGMESH_GEO_MODEL_BUILDER__
#define __RINGMESH_GEO_MODEL_BUILDER__

#include <ringmesh/common.h>
#include <third_party/zlib/unzip.h>

#include <vector>
#include <string>
#include <stack>

#include <geogram/basic/line_stream.h>

#include <ringmesh/geo_model_editor.h>
#define MAX_FILENAME 512
#define READ_SIZE 8192
/*!
 * @file ringmesh/geo_model_builder.h
 * @brief Classes to build GeoModel from various inputs
 * @author Jeanne Pellerin
 */

namespace GEO {
    class Mesh ;
}

namespace RINGMesh {
    class GeoModelRegionFromSurfaces ;
    class GeoModelElementFromMesh ;
}

namespace RINGMesh {
    /*!
     * @brief First draft of flags to build a GeoModel
     * @todo Implements functions to set, access the values, depending on what ?
     * To check the consistency of the options. What do we do about the other elements ? [JP] 
     * @todo We need to keep track of the status of the GeoModel when building it:
     * same flags or some others ?    
     */
    class GeoModelBuildingFlags {
    public:
        GeoModelBuildingFlags()
        {
            compute_corners = false ;
            compute_lines = false ;
            compute_surfaces = false ;
            compute_regions_brep = false ;
            compute_regions_mesh = false ;
        }
        bool compute_corners ;
        bool compute_lines ;
        bool compute_surfaces ;
        bool compute_regions_brep ;
        bool compute_regions_mesh ;
    } ;

    /*!
     * @brief Base class for all classes building a GeoModel.
     * @details Derive from this class to build or modify a GeoModel. 
     * @note NON Geometry related modifications are in GeoModelEditor class.
     * @todo To refactor and rename. We need a GeoModelTopologyEditor 
     * and a GeoModelGeometryEditor
     */
    class RINGMESH_API GeoModelBuilder: public GeoModelEditor {
    public:
        GeoModelBuilder( GeoModel& model )
            : GeoModelEditor( model ), options_()
        {
        }
        virtual ~GeoModelBuilder() ;

        /*!
         * @todo Implements sot that it returns true if the input options are consistent
         */
        void set_options( const GeoModelBuildingFlags& options )
        {
            options_ = options ;
        }

        /*!
         * @brief Copy all element meshes from the input geomodel
         * @pre The model under construction has exaclty the same number of elements
         * than the input geomodel.
         */
        void copy_meshes( const GeoModel& from ) ;

        void copy_meshes( const GeoModel& from, GME::TYPE element_type ) ;

        void assign_mesh_to_element( const Mesh& mesh, GME::gme_t to ) ;

        /*!
         * \name Set element geometry from geometrical positions
         * @{
         */
        void set_element_vertex(
            const GME::gme_t& t,
            index_t v,
            const vec3& point,
            bool update ) ;

        void set_element_vertices(
            const GME::gme_t& element_id,
            const std::vector< vec3 >& points,
            bool clear ) ;

        void set_corner( index_t corner_id, const vec3& point ) ;

        void set_line( index_t id, const std::vector< vec3 >& vertices ) ;

        void set_surface_geometry(
            index_t surface_id,
            const std::vector< vec3 >& surface_vertices,
            const std::vector< index_t >& surface_facets,
            const std::vector< index_t >& surface_facet_ptr ) ;

        void set_region_geometry(
            index_t region_id,
            const std::vector< vec3 >& points,
            const std::vector< index_t >& tetras ) ;

        /*! @}
         * \name Set element geometry using global GeoModel vertices
         * @{
         */
        void set_element_vertex(
            const GME::gme_t& id,
            index_t v,
            index_t model_vertex ) ;

        void set_element_vertices(
            const GME::gme_t& element_id,
            const std::vector< index_t >& model_vertices,
            bool clear ) ;

        index_t add_unique_vertex( const vec3& p ) ;

        void set_corner( index_t corner_id, index_t unique_vertex ) ;

        void set_line( index_t id, const std::vector< index_t >& unique_vertices ) ;

        void set_surface_geometry(
            index_t surface_id,
            const std::vector< index_t >& surface_vertices,
            const std::vector< index_t >& surface_facets,
            const std::vector< index_t >& surface_facet_ptr ) ;

        void set_surface_geometry(
            index_t surface_id,
            const std::vector< index_t >& corners,
            const std::vector< index_t >& facet_ptr ) ;

        void set_surface_geometry(
            index_t surface_id,
            const std::vector< index_t >& triangle_corners ) ;

        void set_surface_geometry_with_adjacencies(
            index_t surface_id,
            const std::vector< index_t >& triangle_corners,
            const std::vector< index_t >& adjacent_triangles ) ;

        void set_region_geometry(
            index_t region_id,
            const std::vector< index_t >& tet_corners ) ;

        /*! @}
         * \name Misc
         * @{
         */
        index_t find_or_create_duplicate_vertex(
            GeoModelMeshElement& S,
            index_t model_vertex_id,
            index_t surface_vertex_id ) ;

        void cut_surface_by_line( Surface& S, const Line& L ) ;

        void compute_surface_adjacencies( index_t surface_id ) ;

        GME::gme_t find_or_create_corner( const vec3& point ) ;
        GME::gme_t find_or_create_corner( index_t model_point_id ) ;
        GME::gme_t find_or_create_line( const std::vector< vec3 >& vertices ) ;
        GME::gme_t find_or_create_line(
            const std::vector< index_t >& incident_surfaces,
            GME::gme_t first_corner,
            GME::gme_t second_corner ) ;

        void recompute_geomodel_mesh() ;

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
         * @pre Function build_lines_and_corners_from_surfaces must have been called before
         */
        bool build_brep_regions_from_surfaces() ;

        /*
         * @brief From a GeoModel in which only Surface are defined, create corners, contacts
         * and regions depending on the building flags
         * @note Valdity is not checked
         */
        void build_model_from_surfaces() ;

        /*!
         * @brief Finish up model building and complete missing information.
         */
        void end_model() ;

    protected:
        void build_contacts() ;
        void set_surface_facet_adjacencies(
                index_t surface_id,
                const std::vector< index_t >& facets_id,
                const std::vector< index_t >& edges_id,
                const std::vector< index_t >& adjacent_triangles );

    protected:
        /*! Elements to compute from the available elements */
        GeoModelBuildingFlags options_ ;

        /*! Internal information */
        std::vector< GeoModelRegionFromSurfaces* > regions_info_ ;

    private:
        void assign_surface_mesh_facets(
            index_t surface_id,
            const std::vector< index_t >& facets,
            const std::vector< index_t >& facet_ptr ) ;

        void assign_surface_triangle_mesh(
            index_t surface_id,
            const std::vector< index_t >& triangle_vertices ) ;
        void update_facet_corner(
            Surface& S,
            const std::vector< index_t >& facets,
            index_t old,
            index_t neu ) ;
        void assign_surface_triangle_mesh(
            index_t surface_id,
            const std::vector< index_t >& triangle_vertices,
            const std::vector< index_t >& adjacent_triangles ) ;

        void assign_region_tet_mesh(
            index_t region_id,
            const std::vector< index_t >& tet_vertices ) const ;

        void duplicate_surface_vertices_along_line( Surface& S, const Line& L ) ;
        void disconnect_surface_facets_along_line_edges(
            Surface& S,
            const Line& L ) ;
    } ;

    /*!
     * @brief To build a GeoModel from a set of disconnected polygonal surfaces
     */
    class RINGMESH_API GeoModelBuilderSurfaceMesh: public GeoModelBuilder {
    public:
        GeoModelBuilderSurfaceMesh( GeoModel& model, const GEO::Mesh& mesh )
            : GeoModelBuilder( model ), mesh_( mesh )
        {
            options_.compute_lines = true ;
            options_.compute_corners = true ;
            options_.compute_regions_brep = true ;
        }
        void build_polygonal_surfaces_from_connected_components() ;

    private:
        const GEO::Mesh& mesh_ ;
    } ;

    /*!
     * @brief Builder of a GeoModel from a simplicial surface or volumetric Mesh 
     * @details Regions and Surfaces are identified with an attribute of type index_t
     * on the mesh cells or facet 
     */
    class RINGMESH_API GeoModelBuilderMesh: public GeoModelBuilder {
    public:
        GeoModelBuilderMesh(
            GeoModel& model,
            const GEO::Mesh& mesh,
            const std::string& surface_attribute_name,
            const std::string& region_attribute_name )
            :
                GeoModelBuilder( model ),
                mesh_( mesh ),
                surface_builder_( nil ),
                region_builder_( nil ),
                surface_attribute_name_( surface_attribute_name ),
                region_attribute_name_( region_attribute_name )
        {
            initialize_surface_builder() ;
            initialize_region_builder() ;
            add_mesh_vertices_to_model() ;
        }

        virtual ~GeoModelBuilderMesh() ;

        /*!
         * @brief Prepare a Mesh so that it can be used to build one GeoModel Surfaces
         * @details Repairs the mesh, triangulates it, computes a connected component 
         * attribute of type index_t on the mesh facets and removes colocated vertices. 
         */
        static void prepare_surface_mesh_from_connected_components(
            GEO::Mesh& mesh,
            const std::string& created_surface_attribute ) ;

        void create_and_build_surfaces() ;
        void build_surfaces() ;

        void create_and_build_regions() ;
        void build_regions() ;

        void copy_facet_attribute_from_mesh( const std::string& attribute_name ) ;
        void copy_cell_attribute_from_mesh( const std::string& attribute_name ) ;

    protected:
        /*!
         * @brief Set the unique vertices used to build the GeoModel
         * @details They are cleared when end_model() is called
         */
        void add_mesh_vertices_to_model() ;

        void initialize_surface_builder() ;
        void initialize_region_builder() ;

    private:
        bool is_mesh_valid_for_surface_building() const ;
        bool is_mesh_valid_for_region_building() const ;

    protected:
        const GEO::Mesh& mesh_ ;
        GeoModelElementFromMesh* surface_builder_ ;
        GeoModelElementFromMesh* region_builder_ ;
        std::string surface_attribute_name_ ;
        std::string region_attribute_name_ ;
        index_t nb_surface_attribute_values_ ;
        index_t nb_region_attribute_values_ ;
    } ;

    /*!
     * @brief Abstract class to load and build GeoModels from files 
     */
    class RINGMESH_API GeoModelBuilderFile: public GeoModelBuilder {
    public:
        GeoModelBuilderFile( GeoModel& model, const std::string& filename ) ;

        virtual ~GeoModelBuilderFile()
        {
        }
        void build_model()
        {
            load_file() ;
            end_model() ;
        }

        ///TODO these are temporary protected here. after they will be only in GeoModelBuilderGM
    protected:
        static GME::TYPE match_nb_elements( const char* s ) ;
        static GME::TYPE match_type( const char* s ) ;
        static bool match_high_level_type( const char* s )
        {
            return GME::child_allowed( match_type( s ) ) ;
        }

    private:
        virtual void load_file() = 0 ;
        /*! @todo Implement function to read the lines of the 
         *        file and wrap the GEO::LineInput which is not that easy to use 
         */
    protected:
        std::string filename_ ;
    } ;

    /*!
     * @brief Build a GeoModel from a Gocad Model3D (file_model.ml)
     */
    class RINGMESH_API GeoModelBuilderGocad: public GeoModelBuilderFile {
    public:
        GeoModelBuilderGocad( GeoModel& model, const std::string& filename )
            : GeoModelBuilderFile( model, filename ), file_line_( filename )
        {
            if( !file_line_.OK() ) {
                throw RINGMeshException( "I/O", "Failed to open file " + filename ) ;
            }
        }
        virtual ~GeoModelBuilderGocad()
        {
        }

    private:
        void load_file() ;

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

        /*!
         * Read the coordinates system information of files exported from Gocad.
         * @param[in] in The orientation of z-axis in Gocad. "Elevation" for
         * increasing z toward top and "Depth" for increasing z toward bottom.
         * @return Return 1 if Elevation direction, -1 if Depth direction.
         */
        int read_gocad_coordinates_system( const std::string& in ) ;

    private:
        GEO::LineInput file_line_ ;

        /*!
         * @brief Triangle that set the orientation of a TFACE
         *        in a .ml file
         */
        struct KeyFacet {
            KeyFacet( const vec3& p0, const vec3& p1, const vec3& p2 )
                : p0_( p0 ), p1_( p1 ), p2_( p2 )
            {
            }
            vec3 p0_ ;
            vec3 p1_ ;
            vec3 p2_ ;
        } ;
        std::vector< KeyFacet > key_facets_ ;
    } ;

    /*!
     * @brief Build a GeoModel from a file_model.bm
     * @TODO this class gonna disapear soon
     */
    class RINGMESH_API GeoModelBuilderBM: public GeoModelBuilderFile {
    public:
        GeoModelBuilderBM( GeoModel& model, const std::string& filename )
            : GeoModelBuilderFile( model, filename ), file_line_( filename )
        {
            if( !file_line_.OK() ) {
                throw RINGMeshException( "I/O", "Failed to open file " + filename ) ;
            }
        }
        virtual ~GeoModelBuilderBM()
        {
        }

    private:
        void load_file() ;

    private:
        GEO::LineInput file_line_ ;

    } ;

    class RINGMESH_API GeoModelBuilderGM: public GeoModelBuilderFile {
    public:
        GeoModelBuilderGM( GeoModel& model, const std::string& filename )
            : GeoModelBuilderFile( model, filename )
        {
        }
        virtual ~GeoModelBuilderGM()
        {
        }

    private:
        /*!
         * @brief Load the connectivities. These are how corners are
         * connected to lines, lines connected to surfaces and surfaces
         * connected to regions
         */
        void load_connectivities( GEO::LineInput& file_line ) ;
        /*!
         * @brief Load elements of one type from a zip file
         * @param[in] gme_t the GeoModelElement type
         * @param[in] uz the zip file
         */
        void load_elements( GME::TYPE gme_t, unzFile& uz ) ;

        void load_file() ;
        /*!
         * @brief Unzip a file in a zip file and set it to the current unZIP file
         */
        void unzip_one_file( unzFile& uz, const char filename[MAX_FILENAME] ) ;

        /*!
         * @brief Load the topology. Topology is how corners, lines, surfaces and
         * regions are organized into contacts, interfaces and layers. It also contains
         * basics information on the GeoModel.
         */
        void load_topology( GEO::LineInput& file_line ) ;
    } ;
}

#endif
