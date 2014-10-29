/*[
* Association Scientifique pour la Geologie et ses Applications (ASGA)
* Copyright (c) 1993-2013 ASGA. All Rights Reserved.
*
* This program is a Trade Secret of the ASGA and it is not to be:
* - reproduced, published, or disclosed to other,
* - distributed or displayed,
* - used for purposes or on Sites other than described
*   in the GOCAD Advancement Agreement,
* without the prior written authorization of the ASGA. Licencee
* agrees to attach or embed this Notice on all copies of the program,
* including partial copies or modified versions thereof.
]*/

#ifndef __GRGMESH_TETRA_GEN__
#define __GRGMESH_TETRA_GEN__

#include <grgmeshlib/common.h>
#include <grgmeshlib/smart_pointer.h>
#include <grgmeshlib/vecn.h>
#include <grgmeshlib/utils.h>
#include <grgmeshlib/mixed_mesh.h>

#include <third_party/tetgen.h>

#include <vector>

#ifdef USE_MG_TETRA
    extern "C" {
        #include <meshgems/meshgems.h>
        #include <meshgems/tetra.h>
    }
#endif

namespace GRGMesh {

    class BoundaryModelElement ;
    class MixedMesh ;
    class TetraGen ;

    typedef SmartPointer< TetraGen > TetraGen_var ;

    static const std::vector< vec3 > vector_vec3 ;
    static const std::vector< std::vector< Edge > > vector_edge ;
    class GRGMESH_API TetraGen : public Counted {
    public:
        virtual ~TetraGen() {} ;
        static TetraGen_var instantiate(
            const Tetra_method& method,
            MixedMesh& tetmesh,
            const BoundaryModelElement* region,
            bool add_steiner_points = true,
            const std::vector< vec3 >& internal_vertices = vector_vec3,
            const std::vector< std::vector< Edge > >& well_vertices = vector_edge,
            MixedMesh* background = nil ) ;

        virtual bool tetrahedralize() = 0 ;

        uint32 nb_points() const { return points_.size() ; }
        uint32 nb_internal_points() const { return internal_points_.size() ; }
        uint32 nb_total_points() const { return nb_points() + nb_internal_points() ; }
        uint32 nb_triangles() const { return triangles_.size() / 3 ; }
        uint32 point_index( uint32 f, uint32 v ) const { return triangles_[3*f+v] ; }
        const vec3& point( uint32 f, uint32 v ) const { return points_[triangles_[3*f+v]] ; }
        const vec3& point( uint32 v ) const { return points_[v] ; }
        int32 surface_id( uint32 f ) const {
            for( uint32 i = 1; i < surface_id_.size(); i++ ) {
                if( f < surface_ptr_[i] ) return surface_id_[i-1] ;
            }
            return  surface_id_.back() ;
        }
        int32* surface_id_ptr( uint32 f ) {
            for( uint32 i = 1; i < surface_id_.size(); i++ ) {
                if( f < surface_ptr_[i] ) return &surface_id_[i-1] ;
            }
            return  &surface_id_.back() ;
        }
        int32 well_id( uint32 f ) const {
            for( uint32 i = 1; i < well_ptr_.size(); i++ ) {
                if( f < well_ptr_[i] ) return i-1 ;
            }
            return  well_ptr_.size()-1 ;
        }

    protected:
        TetraGen(
            MixedMesh& tetmesh,
            const BoundaryModelElement* region,
            const std::vector< vec3 >& internal_vertices,
            const std::vector< std::vector< Edge > >& well_edges,
            MixedMesh* background ) ;

        void initialize_sorage( uint32 nb_points, uint32 nb_tets ) ;
        void set_point( uint32 index, double* point ) ;
        void set_tetra( uint32 index, int* tet ) ;
        void set_tetra_adjacent( uint32 index, uint32 face, int32 adj ) ;
        void set_face_marker(
            uint32 tet1,
            uint32 tet2,
            uint32 marker ) ;
        void set_tetra_face_marker(
            uint32 tet,
            uint32 adj,
            uint32 marker ) ;

        void store_edge_attrib() const ;

    protected:
        std::vector< vec3 > points_ ;
        std::vector< vec3 > internal_points_ ;
        std::vector< Edge > well_edges_ ;
        std::vector< uint32 > well_ptr_ ;
        std::vector< int32 > well_indices_ ;
        std::vector< int32 > triangles_ ;
        std::vector< int32 > surface_id_ ;
        std::vector< uint32 > surface_ptr_ ;
        MixedMesh& tetmesh_ ;
        MixedMeshMutator tetmesh_mutator_ ;
        double resolution_ ;
        MixedMesh* background_ ;
    } ;


    class GRGMESH_API TetraGen_TetGen: public TetraGen {
    public:
        TetraGen_TetGen(
            MixedMesh& tetmesh,
            const BoundaryModelElement* region,
            bool add_steiner_points,
            const std::vector< vec3 >& internal_vertices,
            const std::vector< std::vector< Edge > >& well_vertices,
            MixedMesh* background ) ;
        virtual ~TetraGen_TetGen() {} ;

        virtual bool tetrahedralize() ;

    private:
        tetgenio tetgen_input_ ;
        tetgenio tetgen_output_ ;
        tetgenio tetgen_background_ ;
        tetgenbehavior tetgen_args_ ;
    } ;

#ifdef USE_MG_TETRA
    class GRGMESH_API TetraGen_MG_Tetra: public TetraGen {
    public:
        TetraGen_MG_Tetra(
            MixedMesh& tetmesh,
            const BoundaryModelElement* region,
            bool add_steiner_points,
            const std::vector< vec3 >& internal_vertices,
            const std::vector< std::vector< Edge > >& well_vertices,
            MixedMesh* background ) ;
        virtual ~TetraGen_MG_Tetra() ;

        virtual bool tetrahedralize() ;

        static status_t my_message_cb( message_t * msg, void *user_data ) ;
        static status_t get_size_value(
            meshgems_integer i,
            meshgems_real* size,
            void *user_data ) ;

    private:
        double get_resolution_value( int32 i ) ;

    private:
        bool add_steiner_points_ ;
        context_t* context_ ;
        mesh_t* mesh_input_ ;
        mesh_t* mesh_output_ ;
        mesh_t* mesh_background_ ;
        sizemap_t* sizemap_ ;
        tetra_session_t* tms_ ;
    } ;
#endif

}

#endif
