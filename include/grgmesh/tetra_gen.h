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

#include <grgmesh/common.h>
#include <grgmesh/utils.h>

#include <geogram/mesh/mesh.h>
#include <geogram/mesh/mesh_builder.h>
#include <geogram/basic/counted.h>
#include <geogram/basic/smart_pointer.h>
#include <geogram/third_party/tetgen/tetgen.h>

#include <vector>

#ifdef USE_MG_TETRA
    extern "C" {
        #include <meshgems/meshgems.h>
        #include <meshgems/tetra.h>
    }
#endif

namespace GRGMesh {

    class BoundaryModelElement ;
    class TetraGen ;

    typedef GEO::SmartPointer< TetraGen > TetraGen_var ;

    static const std::vector< vec3 > vector_vec3 ;
    static const std::vector< std::vector< Edge > > vector_edge ;
    class GRGMESH_API TetraGen : public GEO::Counted {
    public:
        virtual ~TetraGen() ;
        static TetraGen_var instantiate(
            const TetraMethod& method,
            GEO::Mesh& tetmesh,
            const BoundaryModelElement* region,
            bool add_steiner_points = true,
            const std::vector< vec3 >& internal_vertices = vector_vec3,
            const std::vector< std::vector< Edge > >& well_vertices = vector_edge,
            GEO::Mesh* background = nil ) ;

        virtual bool tetrahedralize() = 0 ;

        index_t nb_points() const { return points_.size() ; }
        index_t nb_internal_points() const { return internal_points_.size() ; }
        index_t nb_total_points() const { return nb_points() + nb_internal_points() ; }
        index_t nb_triangles() const { return triangles_.size() / 3 ; }
        index_t point_index( index_t f, index_t v ) const { return triangles_[3*f+v] ; }
        const vec3& point( index_t f, index_t v ) const { return points_[triangles_[3*f+v]] ; }
        const vec3& point( index_t v ) const { return points_[v] ; }
        signed_index_t surface_id( index_t f ) const {
            for( index_t i = 1; i < surface_id_.size(); i++ ) {
                if( f < surface_ptr_[i] ) return surface_id_[i-1] ;
            }
            return  surface_id_.back() ;
        }
        signed_index_t* surface_id_ptr( index_t f ) {
            for( index_t i = 1; i < surface_id_.size(); i++ ) {
                if( f < surface_ptr_[i] ) return &surface_id_[i-1] ;
            }
            return  &surface_id_.back() ;
        }
        signed_index_t well_id( index_t f ) const {
            for( index_t i = 1; i < well_ptr_.size(); i++ ) {
                if( f < well_ptr_[i] ) return i-1 ;
            }
            return  well_ptr_.size()-1 ;
        }

    protected:
        TetraGen(
            GEO::Mesh& tetmesh,
            const BoundaryModelElement* region,
            const std::vector< vec3 >& internal_vertices,
            const std::vector< std::vector< Edge > >& well_edges,
            GEO::Mesh* background ) ;

        void initialize_storage( index_t nb_points, index_t nb_tets, index_t nb_triangles, index_t nb_lines ) ;
        void set_point( index_t index, double* point ) ;
        void set_tetra( index_t index, int* tet, index_t nb_lines, index_t nb_triangles ) ;
        void set_triangle( index_t index, int * triangle, index_t nb_lines ) ;
        void set_line( index_t index, int * line ) ;
        void set_tetra_adjacent( index_t index, index_t face, signed_index_t adj ) ;
        void set_face_marker(
            index_t tri,
            index_t marker ) ;
        void set_tetra_face_marker(
            index_t tet,
            index_t adj,
            index_t marker ) ;

        void store_edge_attrib() const ;

    protected:
        std::vector< vec3 > points_ ;
        std::vector< vec3 > internal_points_ ;
        std::vector< Edge > well_edges_ ;
        std::vector< index_t > well_ptr_ ;
        std::vector< signed_index_t > well_indices_ ;
        std::vector< signed_index_t > triangles_ ;
        std::vector< signed_index_t > surface_id_ ;
        std::vector< index_t > surface_ptr_ ;
        GEO::Mesh& tetmesh_ ;
        GEO::MeshBuilder tetmesh_builder_ ;
        double resolution_ ;
        GEO::Mesh* background_ ;
        const BoundaryModelElement* region_ ;
    } ;


    class GRGMESH_API TetraGen_TetGen: public TetraGen {
    public:
        TetraGen_TetGen(
            GEO::Mesh& tetmesh,
            const BoundaryModelElement* region,
            bool add_steiner_points,
            const std::vector< vec3 >& internal_vertices,
            const std::vector< std::vector< Edge > >& well_vertices,
            GEO::Mesh* background ) ;
        virtual ~TetraGen_TetGen() {} ;

        virtual bool tetrahedralize() ;

    private:
        GEO_3rdParty::tetgenio tetgen_input_ ;
        GEO_3rdParty::tetgenio tetgen_output_ ;
        GEO_3rdParty::tetgenio tetgen_background_ ;
        GEO_3rdParty::tetgenbehavior tetgen_args_ ;
    } ;

#ifdef USE_MG_TETRA
    class GRGMESH_API TetraGen_MG_Tetra: public TetraGen {
    public:
        TetraGen_MG_Tetra(
            GEO::Mesh& tetmesh,
            const BoundaryModelElement* region,
            bool add_steiner_points,
            const std::vector< vec3 >& internal_vertices,
            const std::vector< std::vector< Edge > >& well_vertices,
            GEO::Mesh* background ) ;
        virtual ~TetraGen_MG_Tetra() ;

        virtual bool tetrahedralize() ;

        static status_t my_message_cb( message_t * msg, void *user_data ) ;
        static status_t get_size_value(
            meshgems_integer i,
            meshgems_real* size,
            void *user_data ) ;

    private:
        double get_resolution_value( signed_index_t i ) ;

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
