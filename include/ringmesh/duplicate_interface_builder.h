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

#ifndef __RINGMESH_DUPLICATE_INTERFACE_BUILDER__
#define __RINGMESH_DUPLICATE_INTERFACE_BUILDER__

#include <ringmesh/common.h>
#include <ringmesh/geo_model_builder.h>

/*!
 * @file ringmesh/duplicate_interface_builder.h
 * @brief Class to duplicate GeoModel Interface to
 * enable sliding along them (faults) and unconformal
 * mesh generation.
 * @author Benjamin Chauvin
 */

namespace RINGMesh {
    class GeoModel ;
}

namespace GEO {
    class Mesh ;
}

namespace RINGMesh {

    class RINGMESH_API DuplicateInterfaceBuilder: public GeoModelBuilder {
    ringmesh_disable_copy(DuplicateInterfaceBuilder) ;
    public:
        DuplicateInterfaceBuilder( GeoModel& model ) ;
        virtual ~DuplicateInterfaceBuilder() ;
        void duplicate_fault_network() ;
    private:
        void homogenize_normal_orientation_surface_all_interfaces() ;
        void homogenize_normal_orientation_surface_one_interface(
            const GeoModelElement& fault_interface,
            std::vector< index_t >& surfaces_to_inverse_normals ) ;
        void get_new_surfaces(
            const GeoModelElement& interface_to_duplicate,
            std::vector< std::vector< index_t > >& to_erase_by_type ) ;
        void build_merged_surfaces(
            const std::map< index_t, std::vector< index_t > >& surfaces_boundary_regions,
            const std::string& side_name,
            std::vector< std::vector< index_t > >& to_erase_by_type,
            const GME::gme_t& sided_interface_gme_t,
            const GeoModelElement& interface_to_duplicate ) ;
        void compute_translation_vectors_duplicated_fault_network_surfaces_and_regions(
            index_t first_new_interface_index,
            const std::vector< std::vector< index_t > >& to_erase_by_type ) ;
        void compute_translation_vectors_duplicated_fault_network(
            index_t first_new_interface_index,
            const std::vector< std::vector< index_t > >& to_erase_by_type ) ;
        void translate_duplicated_fault_network(
            const std::vector< std::vector< index_t > >& to_erase_by_type ) ;
        void save_normals_on_one_new_interface(
            const std::vector< std::vector< index_t > >& to_erase_by_type,
            const GeoModelElement& interface_gme ) const ;
        vec3 get_local_translation_normal(
            const Surface& surface,
            index_t vertex_id_in_surface ) const ;
        vec3 get_local_translation_vector( const vec3& normal ) const ;
        void store_displacement_in_gme(
            const GeoModelMeshElement& gmme,
            index_t vertex_id_in_gmme,
            const vec3& translation ) const ;
        void initialize_translation_attributes(
            const std::vector< std::vector< index_t > >& to_erase_by_type ) ;
        bool is_region_on_right_side_of_sided_interface(
            index_t region_to_check_id,
            const vec3& normal_on_vertex_interface,
            index_t vertex_id_in_region,
            const vec3& vertex_pos ) const ;
        bool is_surface_on_right_side_of_sided_interface(
            index_t surface_to_check_id,
            const vec3& normal_on_vertex_interface,
            index_t vertex_id_in_surface,
            const vec3& vertex_pos ) const ;
        bool is_surface_or_region_on_the_right_side_of_the_fault(
            const GME::gme_t& cur_gme_t,
            const vec3& normal_on_vertex_interface,
            index_t vertex_id_in_gmme,
            const vec3& vertex_pos,
            const GeoModelElement& interface_gme,
            const GeoModelElement& other_side_interface_gme ) const ;
        void set_no_displacement_on_fault_real_extension(
            const std::vector< std::vector< index_t > >& to_erase_by_type ) ;
        bool does_surface_belong_to_interface(
            const Surface& surface,
            const GeoModelElement& interface ) const ;
        void save_normals_on_one_old_interface(
            const GeoModelElement& interface_gme ) const ;
        vec3 get_normal_on_surface_vertex(
            const Surface& surface,
            index_t vertex_id_on_surface ) const ;
        void homogenize_surfaces_around_surface(
            const GeoModelElement& fault_interface,
            const Surface& first_child,
            std::vector< bool >& already_seen,
            std::vector< index_t >& surfaces_to_inverse_normals ) ;
        void inverse_normal_attribute_one_surface( const Surface& surface ) const ;
        void update_region_polarity( const Surface& surface ) ;
        void mutual_cut_between_new_merged_surfaces(
            index_t first_new_interface_index,
            const std::vector< std::vector< index_t > >& to_erase_by_type ) ; /// @todo not a good name
        void add_fake_internal_boudnary_lines_to_merged_surface(
            const std::map< index_t, index_t >& all_surface_lines,
            const std::string& side_name,
            const GME::gme_t& sided_interface_gme_t,
            const GeoModelElement& interface_to_duplicate,
            const GME::gme_t& new_surface_gme_t,
            std::vector< std::vector< index_t > >& to_erase_by_type,
            index_t region_index ) ;
        void save_normal_on_one_surface( const Surface& surface ) const ;
        void split_merged_surface(
            const GME::gme_t& new_surface_gme_t,
            const std::string& side_name,
            const GME::gme_t& sided_interface_gme_t,
            std::vector< std::vector< index_t > >& to_erase_by_type,
            index_t region_index ) ;
        void define_global_motion_relation(
            const std::vector< std::vector< index_t > >& to_erase_by_type ) ;
        void initialize_gme_vertices_links(
            const std::vector< std::vector< index_t > >& to_erase_by_type ) ;

    private:
        class GMEVertexLink {
        ringmesh_disable_copy(GMEVertexLink) ;
        public:
            GMEVertexLink(
                const GMEVertex& gme_vertex,
                const GeoModel& model,
                const std::vector< GMEVertexLink* >& gme_vertices_links ) ;
            ~GMEVertexLink()
            {
            }
            void displace( const vec3& displacement_vector ) ;
            void add_linked_gme_vertex( index_t new_linked_gme_vertex )
            {
                linked_gme_vertices_.push_back( new_linked_gme_vertex ) ;
            }

        private:
            bool has_moved_ ;
            const GeoModel& model_ ;
            const GMEVertex gme_vertex_ ;
            std::vector< index_t > linked_gme_vertices_ ;
            const std::vector< GMEVertexLink* >& gme_vertices_links_ ;
        } ;
    private:
        bool all_meshed_ ;
        std::vector< GMEVertexLink* > gme_vertices_links_ ;
    } ;
}

#endif
