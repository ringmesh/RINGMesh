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

#ifndef __GRGMESH_MIXED_MESH__
#define __GRGMESH_MIXED_MESH__

#include <grgmeshlib/common.h>
#include <grgmeshlib/mesh.h>

#include <vector> 

namespace GRGMesh {

    /*!
     * Mesh which can handle different type of elements.
     */
    class GRGMESH_API MixedMesh: public Mesh {
        friend class MixedMeshMutator ;
        typedef AttributeManager< LINE > LineAttributeManager ;
        typedef AttributeManager< TRGL > TriangleAttributeManager ;
        typedef AttributeManager< QUAD > QuadAttributeManager ;
        typedef AttributeManager< TETRA > TetraAttributeManager ;
        typedef AttributeManager< PYRAMID > PyramidAttributeManager ;
        typedef AttributeManager< PRISM > PrismAttributeManager ;
        typedef AttributeManager< HEXA > HexaAttributeManager ;

    public:
        MixedMesh()
        {
        }
        virtual ~MixedMesh()
        {
        }
        MixedMesh( const MixedMesh& rhs ) { copy( rhs ) ; }
        MixedMesh& operator=( const MixedMesh& rhs ) { copy( rhs ) ; return *this ; }

        virtual uint32 nb_cells() const
        {
            return nb_lines() + nb_triangles() + nb_quad() + nb_tetra()
                + nb_pyramids() + nb_prisms() + nb_hexa() ;
        }
        uint32 nb_lines() const { return std::max( uint64(0), cells_[LINE].size()-1 ) ; }
        uint32 nb_triangles() const { return std::max( uint64(0), cells_[TRGL].size()-1 ) ; }
        uint32 nb_quad() const { return std::max( uint64(0), cells_[QUAD].size()-1 ) ; }
        uint32 nb_tetra() const { return std::max( uint64(0), cells_[TETRA].size()-1 ) ; }
        uint32 nb_pyramids() const { return std::max( uint64(0), cells_[PYRAMID].size()-1 ) ; }
        uint32 nb_prisms() const { return std::max( uint64(0), cells_[PRISM].size()-1 ) ; }
        uint32 nb_hexa() const { return std::max( uint64(0), cells_[HEXA].size()-1 ) ; }
        virtual uint8 nb_vertices_in_cell( uint32 c ) const
        {
            grgmesh_debug_assert( c < nb_cells() ) ;
            return cell_descriptor( c )->nb_vertices ;
        }
        virtual uint8 nb_facets_in_cell( uint32 c ) const
        {
            grgmesh_debug_assert( c < nb_cells() ) ;
            return cell_descriptor( c )->nb_facets ;
        }
        virtual uint8 nb_vertices_in_cell_facet( uint32 c, uint8 f ) const
        {
            grgmesh_debug_assert( c < nb_cells() ) ;
            return cell_descriptor( c )->nb_vertices_in_facet[f] ;
        }
        virtual uint32 cell_begin( uint32 c ) const {
            uint32 real_index ;
            return cells_[cell_type( c, real_index)][real_index] ;
        }
        virtual uint32 cell_end( uint32 c ) const {
            uint32 real_index ;
            return cells_[cell_type( c+1, real_index)][real_index] ;
        }

        ElementType cell_type( uint32 c, uint32& c_index = dummy_uint32 ) const ;

        const CellDescriptor* cell_descriptor( uint32 c ) const
        {
            const CellDescriptor* result = cell_descriptor_[cell_type( c )] ;
            grgmesh_debug_assert( result != 0 ) ;
            return result ;
        }

        uint32 cell_facet_vertex( uint32 c, uint8 f, uint8 v ) const
        {
            const CellDescriptor* desc = cell_descriptor( c ) ;
            grgmesh_debug_assert( f < desc->nb_facets ) ; grgmesh_debug_assert( v < desc->nb_vertices_in_facet[f] ) ;
            return cell_vertex_index( c, desc->facet[f][v] ) ;
        }

        uint32 tetra_vertex_index( uint32 c, uint8 f, uint8 v ) const {
            return vertex_index( cells_[TETRA][4*c + cell_descriptor_[TETRA]->facet[f][v]] ) ;
        }

        uint32 line_vertex_index( uint32 l, uint8 v ) const {
            return vertex_index( cells_[LINE][l] + v ) ;
        }
        uint32 triangle_vertex_index( uint32 t, uint8 v ) const {
            return vertex_index( cells_[TRGL][t] + v ) ;
        }
        uint32 quad_vertex_index( uint32 q, uint8 v ) const {
            return vertex_index( cells_[QUAD][q] + v ) ;
        }
        uint32 tetra_vertex_index( uint32 t, uint8 v ) const {
            return vertex_index( cells_[TETRA][t] + v ) ;
        }
        uint32 pyramid_vertex_index( uint32 p, uint8 v ) const {
            return vertex_index( cells_[PYRAMID][p] + v ) ;
        }
        uint32 prism_vertex_index( uint32 p, uint8 v ) const {
            return vertex_index( cells_[PRISM][p] + v ) ;
        }
        uint32 hexa_vertex_index( uint32 h, uint8 v ) const {
            return vertex_index( cells_[HEXA][h] + v ) ;
        }

        const vec3& line_vertex( uint32 l, uint8 v ) const {
            return vertex( line_vertex_index( l, v ) ) ;
        }
        const vec3& triangle_vertex( uint32 t, uint8 v ) const {
            return vertex( triangle_vertex_index( t, v ) ) ;
        }
        const vec3& quad_vertex( uint32 q, uint8 v ) const {
            return vertex( quad_vertex_index( q, v ) ) ;
        }
        const vec3& tetra_vertex( uint32 t, uint8 v ) const {
            return vertex( tetra_vertex_index( t, v ) ) ;
        }
        const vec3& pyramid_vertex( uint32 p, uint8 v ) const {
            return vertex( pyramid_vertex_index( p, v ) ) ;
        }
        const vec3& prism_vertex( uint32 p, uint8 v ) const {
            return vertex( prism_vertex_index( p, v ) ) ;
        }
        const vec3& hexa_vertex( uint32 h, uint8 v ) const {
            return vertex( hexa_vertex_index( h, v ) ) ;
        }

        float64 get_nearest_point_in_cell( const vec3& p, uint32 c, vec3& nearest_p ) const ;
        float64 get_nearest_point_in_line( const vec3& p, uint32 c, vec3& nearest_p ) const ;
        float64 get_nearest_point_in_triangle( const vec3& p, uint32 c, vec3& nearest_p ) const ;
        float64 get_nearest_point_in_quad( const vec3& p, uint32 c, vec3& nearest_p ) const ;
        float64 get_nearest_point_in_tetra( const vec3& p, uint32 c, vec3& nearest_p ) const ;
        float64 get_nearest_point_in_pyramid( const vec3& p, uint32 c, vec3& nearest_p ) const ;
        float64 get_nearest_point_in_prism( const vec3& p, uint32 c, vec3& nearest_p ) const ;
        float64 get_nearest_point_in_hexa( const vec3& p, uint32 c, vec3& nearest_p ) const ;

        LineAttributeManager* line_attribute_manager() const
        {
            return const_cast< LineAttributeManager* >( &line_attribute_manager_ ) ;
        }
        TriangleAttributeManager* triangle_attribute_manager() const
        {
            return const_cast< TriangleAttributeManager* >( &triangle_attribute_manager_ ) ;
        }
        QuadAttributeManager* quad_attribute_manager() const
        {
            return const_cast< QuadAttributeManager* >( &quad_attribute_manager_ ) ;
        }
        TetraAttributeManager* tetra_attribute_manager() const
        {
            return const_cast< TetraAttributeManager* >( &tetra_attribute_manager_ ) ;
        }
        PyramidAttributeManager* pyramid_attribute_manager() const
        {
            return const_cast< PyramidAttributeManager* >( &pyramid_attribute_manager_ ) ;
        }
        PrismAttributeManager* prism_attribute_manager() const
        {
            return const_cast< PrismAttributeManager* >( &prism_attribute_manager_ ) ;
        }
        HexaAttributeManager* hexa_attribute_manager() const
        {
            return const_cast< HexaAttributeManager* >( &hexa_attribute_manager_ ) ;
        }
    private:
        void copy( const MixedMesh& rhs )
        {
            for( uint8 i = 0; i < 7; i++) cells_[i] = rhs.cells_[i] ;
        }

    private:
        std::vector< uint32 > cells_[7] ;
        static const CellDescriptor* cell_descriptor_[7] ;


        LineAttributeManager line_attribute_manager_ ;
        TriangleAttributeManager triangle_attribute_manager_ ;
        QuadAttributeManager quad_attribute_manager_ ;
        TetraAttributeManager tetra_attribute_manager_ ;
        PyramidAttributeManager pyramid_attribute_manager_ ;
        PrismAttributeManager prism_attribute_manager_ ;
        HexaAttributeManager hexa_attribute_manager_ ;
    } ;

    template< class ATTRIBUTE >
    class LineAttribute: public Attribute< LINE, ATTRIBUTE > {
    public:
        typedef Attribute< LINE, ATTRIBUTE > superclass ;

        void bind( MixedMesh* mesh, const std::string& name )
        {
            superclass::bind( mesh->line_attribute_manager(), mesh->nb_lines(),
                name ) ;
        }

        void bind( MixedMesh* mesh )
        {
            superclass::bind( mesh->line_attribute_manager(),
                mesh->nb_lines() ) ;
        }

        LineAttribute()
        {
        }

        LineAttribute( MixedMesh* mesh )
        {
            bind( mesh ) ;
        }

        LineAttribute( MixedMesh* mesh, const std::string& name )
        {
            bind( mesh, name ) ;
        }

        static bool is_defined( MixedMesh* mesh, const std::string& name )
        {
            return superclass::is_defined( mesh->line_attribute_manager(), name ) ;
        }
    } ;

    template< class ATTRIBUTE >
    class TriangleAttribute: public Attribute< TRGL, ATTRIBUTE > {
    public:
        typedef Attribute< TRGL, ATTRIBUTE > superclass ;

        void bind( MixedMesh* mesh, const std::string& name )
        {
            superclass::bind( mesh->triangle_attribute_manager(), mesh->nb_triangles(),
                name ) ;
        }

        void bind( MixedMesh* mesh )
        {
            superclass::bind( mesh->triangle_attribute_manager(),
                mesh->nb_triangles() ) ;
        }

        TriangleAttribute()
        {
        }

        TriangleAttribute( MixedMesh* mesh )
        {
            bind( mesh ) ;
        }

        TriangleAttribute( MixedMesh* mesh, const std::string& name )
        {
            bind( mesh, name ) ;
        }

        static bool is_defined( MixedMesh* mesh, const std::string& name )
        {
            return superclass::is_defined( mesh->triangle_attribute_manager(), name ) ;
        }
    } ;

    template< class ATTRIBUTE >
    class QuadAttribute: public Attribute< QUAD, ATTRIBUTE > {
    public:
        typedef Attribute< QUAD, ATTRIBUTE > superclass ;

        void bind( MixedMesh* mesh, const std::string& name )
        {
            superclass::bind( mesh->quad_attribute_manager(), mesh->nb_quad(),
                name ) ;
        }

        void bind( MixedMesh* mesh )
        {
            superclass::bind( mesh->quad_attribute_manager(),
                mesh->nb_quad() ) ;
        }

        QuadAttribute()
        {
        }

        QuadAttribute( MixedMesh* mesh )
        {
            bind( mesh ) ;
        }

        QuadAttribute( MixedMesh* mesh, const std::string& name )
        {
            bind( mesh, name ) ;
        }

        static bool is_defined( MixedMesh* mesh, const std::string& name )
        {
            return superclass::is_defined( mesh->quad_attribute_manager(), name ) ;
        }
    } ;

    template< class ATTRIBUTE >
    class TetraAttribute: public Attribute< TETRA, ATTRIBUTE > {
    public:
        typedef Attribute< TETRA, ATTRIBUTE > superclass ;

        void bind( MixedMesh* mesh, const std::string& name )
        {
            superclass::bind( mesh->tetra_attribute_manager(), mesh->nb_tetra(),
                name ) ;
        }

        void bind( MixedMesh* mesh )
        {
            superclass::bind( mesh->tetra_attribute_manager(),
                mesh->nb_tetra() ) ;
        }

        TetraAttribute()
        {
        }

        TetraAttribute( MixedMesh* mesh )
        {
        	bind( mesh ) ;
        }

        TetraAttribute( MixedMesh* mesh, const std::string& name )
        {
        	bind( mesh, name ) ;
        }


        static bool is_defined( MixedMesh* mesh, const std::string& name )
        {
            return superclass::is_defined( mesh->tetra_attribute_manager(), name ) ;
        }
    } ;

    template< class ATTRIBUTE >
    class PyramidAttribute: public Attribute< PYRAMID, ATTRIBUTE > {
    public:
        typedef Attribute< PYRAMID, ATTRIBUTE > superclass ;

        void bind( MixedMesh* mesh, const std::string& name )
        {
            superclass::bind( mesh->pyramid_attribute_manager(), mesh->nb_pyramids(),
                name ) ;
        }

        void bind( MixedMesh* mesh )
        {
            superclass::bind( mesh->pyramid_attribute_manager(),
                mesh->nb_pyramids() ) ;
        }

        PyramidAttribute()
        {
        }

        PyramidAttribute( MixedMesh* mesh )
        {
            bind( mesh ) ;
        }

        PyramidAttribute( MixedMesh* mesh, const std::string& name )
        {
            bind( mesh, name ) ;
        }

        static bool is_defined( MixedMesh* mesh, const std::string& name )
        {
            return superclass::is_defined( mesh->pyramid_attribute_manager(), name ) ;
        }
    } ;

    template< class ATTRIBUTE >
    class PrismAttribute: public Attribute< PRISM, ATTRIBUTE > {
    public:
        typedef Attribute< PRISM, ATTRIBUTE > superclass ;

        void bind( MixedMesh* mesh, const std::string& name )
        {
            superclass::bind( mesh->prism_attribute_manager(), mesh->nb_prisms(),
                name ) ;
        }

        void bind( MixedMesh* mesh )
        {
            superclass::bind( mesh->prism_attribute_manager(),
                mesh->nb_prisms() ) ;
        }

        PrismAttribute()
        {
        }

        PrismAttribute( MixedMesh* mesh )
        {
            bind( mesh ) ;
        }

        PrismAttribute( MixedMesh* mesh, const std::string& name )
        {
            bind( mesh, name ) ;
        }

        static bool is_defined( MixedMesh* mesh, const std::string& name )
        {
            return superclass::is_defined( mesh->prism_attribute_manager(), name ) ;
        }
    } ;

    template< class ATTRIBUTE >
    class HexaAttribute: public Attribute< HEXA, ATTRIBUTE > {
    public:
        typedef Attribute< HEXA, ATTRIBUTE > superclass ;

        void bind( MixedMesh* mesh, const std::string& name )
        {
            superclass::bind( mesh->hexa_attribute_manager(), mesh->nb_hexa(),
                name ) ;
        }

        void bind( MixedMesh* mesh )
        {
            superclass::bind( mesh->hexa_attribute_manager(),
                mesh->nb_hexa() ) ;
        }

        HexaAttribute()
        {
        }

        HexaAttribute( MixedMesh* mesh )
        {
            bind( mesh  ) ;
        }

        HexaAttribute( MixedMesh* mesh, const std::string& name )
        {
            bind( mesh, name ) ;
        }

        static bool is_defined( MixedMesh* mesh, const std::string& name )
        {
            return superclass::is_defined( mesh->hexa_attribute_manager(), name ) ;
        }
    } ;

    class GRGMESH_API MixedMeshMutator: public MeshMutator {
    public:
        MixedMeshMutator( MixedMesh& mixed_mesh )
            : MeshMutator( mixed_mesh ), mixed_mesh_( mixed_mesh )
        {
        }
        MixedMeshMutator( const MixedMesh& mixed_mesh )
            :
                MeshMutator( const_cast< MixedMesh& >( mixed_mesh ) ),
                mixed_mesh_( const_cast< MixedMesh& >( mixed_mesh ) )
        {
        }
        virtual ~MixedMeshMutator() {}

        std::vector< uint32 >* cells() { return mixed_mesh_.cells_ ; }
        std::vector< uint32 >& lines() { return mixed_mesh_.cells_[LINE] ; }
        std::vector< uint32 >& triangles() { return mixed_mesh_.cells_[TRGL] ; }
        std::vector< uint32 >& quad() { return mixed_mesh_.cells_[QUAD] ; }
        std::vector< uint32 >& tetra() { return mixed_mesh_.cells_[TETRA] ; }
        std::vector< uint32 >& pyramids() { return mixed_mesh_.cells_[PYRAMID] ; }
        std::vector< uint32 >& prisms() { return mixed_mesh_.cells_[PRISM] ; }
        std::vector< uint32 >& hexa() { return mixed_mesh_.cells_[HEXA] ; }

    private:
        MixedMesh& mixed_mesh_ ;
    } ;

}

#endif
