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

#ifndef __RINGMESH_DUPLICATE_FNTK_BUILDER__
#define __RINGMESH_DUPLICATE_FNTK_BUILDER__

#include <ringmesh/basic/common.h>
#include <ringmesh/geomodel/geomodel_builder.h>

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

namespace RINGMesh {

    class RINGMESH_API DuplicateInterfaceBuilder: public GeoModelBuilder {
    ringmesh_disable_copy(DuplicateInterfaceBuilder) ;
    public:
        DuplicateInterfaceBuilder( GeoModel& model) ;
        virtual ~DuplicateInterfaceBuilder() ;
        void duplicate_fault_network( bool gap ) ;
    private:
        const GeoModelGeologicalEntity& interface( index_t interface_id ) const ;
        void check_geomodel_validity_for_duplication() ;
        void build_new_fault_surfaces(
            std::vector< std::vector< index_t > >& to_erase_by_type ) ;
        void flag_corners_lines_contacts_to_be_deleted(
            std::vector< std::vector< index_t > >& to_erase_by_type ) ;
        void add_hole_between_faults(
            std::vector< std::vector< index_t > >& to_erase_by_type,
            index_t nb_initial_interfaces ) ;
        void delete_old_entities(
            std::vector< std::vector< index_t > >& to_erase_by_type ) ;
        void rebuild_valid_geomodel() ;

        void homogenize_normal_orientation_surface_all_interfaces() ;
        void homogenize_normal_orientation_surface_one_interface(
            index_t fault_interface_id,
            std::vector< index_t >& surfaces_to_inverse_normals ) ;
        void get_new_surfaces(
            index_t interface_to_duplicate_id,
            std::vector< std::vector< index_t > >& to_erase_by_type ) ;
        void build_merged_surfaces(
            const std::map< index_t, std::vector< index_t > >& surfaces_boundary_regions,
            const std::string& side_name,
            std::vector< std::vector< index_t > >& to_erase_by_type,
            index_t sided_interface_id,
            index_t interface_to_duplicate_id ) ;
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
            index_t interface_id ) const ;
        void invert_normals_of_surface_list(
            std::vector< index_t >& surfaces_to_inverse_normals ) ;
        vec3 get_local_translation_normal(
            const Surface& surface,
            index_t vertex_id_in_surface ) const ;
        vec3 get_local_translation_vector( const vec3& normal ) const ;
        void store_displacement_in_gme(
            const GeoModelMeshEntity& gmme,
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
            const gme_t& cur_gme_t,
            const vec3& normal_on_vertex_interface,
            index_t vertex_id_in_gmme,
            const vec3& vertex_pos,
            index_t interface_id,
            index_t other_side_interface_id ) const ;
        void set_no_displacement_on_fault_real_extension(
            const std::vector< std::vector< index_t > >& to_erase_by_type ) ;
        bool does_surface_belong_to_interface(
            index_t surface_id,
            index_t interface_id ) const ;
        void save_normals_on_one_old_interface( index_t interface_id ) const ;
        vec3 get_normal_on_surface_vertex(
            const Surface& surface,
            index_t vertex_id_on_surface ) const ;
        void homogenize_surfaces_around_surface(
            index_t fault_interface_id,
            const Surface& first_child,
            std::vector< bool >& already_seen,
            std::vector< index_t >& surfaces_to_inverse_normals ) ;
        void inverse_normal_attribute_one_surface( const Surface& surface ) const ;
        void update_region_polarity( const Surface& surface ) ;
        void add_fake_internal_boudnary_lines_to_merged_surface(
            const std::map< index_t, index_t >& all_surface_lines,
            const std::string& side_name,
            index_t sided_interface_id,
            index_t interface_to_duplicate_id,
            index_t new_surface_id,
            std::vector< std::vector< index_t > >& to_erase_by_type,
            index_t region_index ) ;
        void save_normal_on_one_surface( const Surface& surface ) const ;
        void split_merged_surface(
            index_t new_surface_id,
            const std::string& side_name,
            index_t sided_interface_id,
            std::vector< std::vector< index_t > >& to_erase_by_type,
            index_t region_index ) ;
        void define_global_motion_relation(
            const std::vector< std::vector< index_t > >& to_erase_by_type ) ;
        void do_define_motion_relation(
            const std::vector< std::vector< index_t > >& to_erase_by_type ) ;
        void do_define_motion_relation_on_not_voi_surface(
            const Surface& cur_surface,
            index_t surf_facet_itr,
            index_t v_id_in_surf,
            const Region& reg1,
            GEO::Attribute< index_t >& id_in_link_vector_surf,
            GEO::Attribute< index_t >& id_in_link_vector_reg1 ) ;
        void do_define_motion_relation_on_voi_surface(
            const Surface& cur_surface,
            index_t surf_facet_itr,
            index_t v_id_in_surf,
            const Region& reg1,
            GEO::Attribute< index_t >& id_in_link_vector_surf,
            GEO::Attribute< index_t >& id_in_link_vector_reg1 ) ;
        void initialize_gme_vertices_links(
            const std::vector< std::vector< index_t > >& to_erase_by_type ) ;
        void link_surf_vertex_id_to_reg_vertex_id(
            index_t link_id_surf,
            index_t link_id_reg ) ;
        void fill_reg_nn_searches() ;
        index_t find_region_vertex_id_from_surface_facet_among_colocated_points_in_one_region(
            const std::vector< GMEVertex >& gme_vertices,
            std::vector< index_t > found_gmev_reg,
            const Surface& cur_surface,
            index_t surf_facet_itr,
            const Region& reg,
            index_t surf_v_id_in_gmm ) const ;
        index_t find_reg_vertex_id_in_facet_reg_matching_surf_vertex_id_in_gmm(
            const Region& reg,
            index_t reg_facet_id,
            index_t surf_v_id_in_gmm ) const ;
        void set_no_displacement_on_gme_sharing_vertex(
            index_t vertex_id_in_gmm,
            const std::vector< std::vector< index_t > >& to_erase_by_type ) ;
        bool displace_corner(
            const Corner& corner,
            const Line& line_one_in_boundary ) const ;

        /// TODO copy paste from removal of remove entity, and it is in private in GeoModelEditor
        index_t entity_type_to_index( const EntityType& type ) const
        {
            return find( all_entity_types_, type ) ;
        }
        /// TODO copy paste from removal of remove entity, and it is in private in GeoModelEditor
        void fill_entity_type_to_index_map()
        {
            const EntityTypeManager& manager = entity_type_manager() ;
            all_entity_types_.insert( all_entity_types_.end(),
                manager.mesh_entity_types().begin(),
                manager.mesh_entity_types().end() ) ;

            all_entity_types_.insert( all_entity_types_.end(),
                manager.geological_entity_types().begin(),
                manager.geological_entity_types().end() ) ;
        }
        /// TODO copy paste from removal of remove entity, and it is in private in GeoModelEditor
        const EntityType& index_to_entity_type( index_t index ) const
        {
            return all_entity_types_.at( index ) ;
        }

        void remove_gap() ;

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
            const GeoModel& geomodel_ ;
            const GMEVertex gme_vertex_ ;
            std::vector< index_t > linked_gme_vertices_ ;
            const std::vector< GMEVertexLink* >& gme_vertices_links_ ;
        private:
            static const std::string id_in_link_vector_attribute_name_ ;
        } ;
    private:
        bool all_meshed_ ;
        std::vector< GMEVertexLink* > gme_vertices_links_ ;
        std::vector< const NNSearch* > reg_nn_searches_ ;

        /// TODO copy paste from removal of remove entity
        std::vector< EntityType > all_entity_types_ ;
    private:
        static const std::string translation_attribute_name_ ;
        static const std::string normal_attribute_name_ ;
    } ;
}

#endif
