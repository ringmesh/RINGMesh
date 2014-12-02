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

#include <grgmesh/common.h>
#include <grgmesh/matrix.h>
#include <grgmesh/mesh.h>
#include <vector> 

namespace GRGMesh {

    struct AdjInfo {

    } ;

    /*!
     * Mesh which can handle different type of elements.
     */
    class GRGMESH_API MixedMesh: public Mesh {
        friend class MixedMeshMutator ;
        friend class MixedMeshBuilder ;
        typedef AttributeManager< LINE > LineAttributeManager ;
        typedef AttributeManager< TRGL > TriangleAttributeManager ;
        typedef AttributeManager< QUAD > QuadAttributeManager ;
        typedef AttributeManager< TETRA > TetraAttributeManager ;
        typedef AttributeManager< PYRAMID > PyramidAttributeManager ;
        typedef AttributeManager< PRISM > PrismAttributeManager ;
        typedef AttributeManager< HEXA > HexaAttributeManager ;
        typedef SparseMatrix< AdjInfo > AdjMatrix ;

    public:
        MixedMesh()
        {
        }
        virtual ~MixedMesh()
        {
        }
        MixedMesh( const MixedMesh& rhs ) { copy( rhs ) ; }
        MixedMesh& operator=( const MixedMesh& rhs ) { copy( rhs ) ; return *this ; }

        //    _  _            _                    __
        //   | \| |_  _ _ __ | |__  ___ _ _   ___ / _|
        //   | .` | || | '  \| '_ \/ -_) '_| / _ \  _|
        //   |_|\_|\_,_|_|_|_|_.__/\___|_|   \___/_|
        //
        virtual uint32 nb_cells() const
        {
            return nb_lines() + nb_triangles() + nb_quad() + nb_tetra()
                + nb_pyramids() + nb_prisms() + nb_hexa() ;
        } 
        // Sans les cast �a passe pas a la compile dans VS
        uint32 nb_lines() const { return std::max( uint64(0), (uint64) cells_[LINE].size() ) ; }
        uint32 nb_triangles() const { return std::max( uint64(0),  (uint64) cells_[TRGL].size() ) ; }
        uint32 nb_quad() const { return std::max( uint64(0),  (uint64) cells_[QUAD].size() ) ; }
        uint32 nb_tetra() const { return std::max( uint64(0),  (uint64) cells_[TETRA].size() ) ; }
        uint32 nb_pyramids() const { return std::max( uint64(0),  (uint64) cells_[PYRAMID].size() ) ; }
        uint32 nb_prisms() const { return std::max( uint64(0),  (uint64) cells_[PRISM].size() ) ; }
        uint32 nb_hexa() const { return std::max( uint64(0),  (uint64) cells_[HEXA].size() ) ; }
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

        //     ___     _ _   _
        //    / __|___| | | | |_ _  _ _ __  ___
        //   | (__/ -_) | | |  _| || | '_ \/ -_)
        //    \___\___|_|_|  \__|\_, | .__/\___|
        //                       |__/|_|
        const CellDescriptor* cell_descriptor( uint32 c ) const {
            const CellDescriptor* result = cell_descriptor_[cell_type( c )] ;
            return result ;
        }
        ElementType cell_type( uint32 c, uint32& c_index = dummy_uint32 ) const ;
        uint32 global_index( const ElementType& type, const uint32 index ) const ;

        //     ___     _ _    __             _                 _             _         _
        //    / __|___| | |  / _|__ _ __ ___| |_  __ _____ _ _| |_ _____ __ (_)_ _  __| |_____ __
        //   | (__/ -_) | | |  _/ _` / _/ -_)  _| \ V / -_) '_|  _/ -_) \ / | | ' \/ _` / -_) \ /
        //    \___\___|_|_| |_| \__,_\__\___|\__|  \_/\___|_|  \__\___/_\_\ |_|_||_\__,_\___/_\_\
        //
        uint32 cell_vertex_index( uint32 c, uint8 f, uint8 v ) const {
            uint32 c_index ;
            const ElementType& type = cell_type( c, c_index ) ;
            return vertex_index( cells_[type][c_index] + cell_descriptor_[type]->facet[f][v] ) ;
        }
        uint32 line_vertex_index( uint32 c, uint8 f, uint8 v ) const {
            return vertex_index( cells_[LINE][c] + cell_descriptor_[LINE]->facet[f][v] ) ;
        }
        uint32 triangle_vertex_index( uint32 c, uint8 f, uint8 v ) const {
            return vertex_index( cells_[TRGL][c] + cell_descriptor_[TRGL]->facet[f][v] ) ;
        }
        uint32 quad_vertex_index( uint32 c, uint8 f, uint8 v ) const {
            return vertex_index( cells_[QUAD][c] + cell_descriptor_[QUAD]->facet[f][v] ) ;
        }
        uint32 tetra_vertex_index( uint32 c, uint8 f, uint8 v ) const {
            return vertex_index( cells_[TETRA][c] + cell_descriptor_[TETRA]->facet[f][v] ) ;
        }
        uint32 pyramid_vertex_index( uint32 c, uint8 f, uint8 v ) const {
            return vertex_index( cells_[PYRAMID][c] + cell_descriptor_[PYRAMID]->facet[f][v] ) ;
        }
        uint32 prism_vertex_index( uint32 c, uint8 f, uint8 v ) const {
            return vertex_index( cells_[PRISM][c] + cell_descriptor_[PRISM]->facet[f][v] ) ;
        }
        uint32 hexa_vertex_index( uint32 c, uint8 f, uint8 v ) const {
            return vertex_index( cells_[HEXA][c] + cell_descriptor_[HEXA]->facet[f][v] ) ;
        }
        void flip_vertex( uint32 c, uint8 v1, uint8 v2 ) {
            uint32 c_index ;
            const ElementType& type = cell_type( c, c_index ) ;
            uint32 new_v1 = vertex_indices_[ cells_[TRGL][c] + v2 ] ;
            uint32 new_v2 = vertex_indices_[ cells_[TRGL][c] + v1 ] ;
            vertex_indices_[ cells_[TRGL][c] + v1 ] = new_v1 ;
            vertex_indices_[ cells_[TRGL][c] + v2 ] = new_v2 ;

        }

        //     ___     _ _    __             _                 _
        //    / __|___| | |  / _|__ _ __ ___| |_  __ _____ _ _| |_ _____ __
        //   | (__/ -_) | | |  _/ _` / _/ -_)  _| \ V / -_) '_|  _/ -_) \ /
        //    \___\___|_|_| |_| \__,_\__\___|\__|  \_/\___|_|  \__\___/_\_\
        //
        const vec3& cell_vertex( uint32 c, uint8 f, uint8 v ) const {
            return vertex( cell_vertex_index( c, f, v ) ) ;
        }
        const vec3& line_vertex( uint32 l, uint8 f, uint8 v ) const {
            return vertex( line_vertex_index( l, f, v ) ) ;
        }
        const vec3& triangle_vertex( uint32 t, uint8 f, uint8 v ) const {
            return vertex( triangle_vertex_index( t, f, v ) ) ;
        }
        const vec3& quad_vertex( uint32 q, uint8 f, uint8 v ) const {
            return vertex( quad_vertex_index( q, f, v ) ) ;
        }
        const vec3& tetra_vertex( uint32 t, uint8 f, uint8 v ) const {
            return vertex( tetra_vertex_index( t, f, v ) ) ;
        }
        const vec3& pyramid_vertex( uint32 p, uint8 f, uint8 v ) const {
            return vertex( pyramid_vertex_index( p, f, v ) ) ;
        }
        const vec3& prism_vertex( uint32 p, uint8 f, uint8 v ) const {
            return vertex( prism_vertex_index( p, f, v ) ) ;
        }
        const vec3& hexa_vertex( uint32 h, uint8 f, uint8 v ) const {
            return vertex( hexa_vertex_index( h, f, v ) ) ;
        }

        //     ___     _ _               _             _         _
        //    / __|___| | | __ _____ _ _| |_ _____ __ (_)_ _  __| |_____ __
        //   | (__/ -_) | | \ V / -_) '_|  _/ -_) \ / | | ' \/ _` / -_) \ /
        //    \___\___|_|_|  \_/\___|_|  \__\___/_\_\ |_|_||_\__,_\___/_\_\
        //
        uint32 cell_vertex_index( uint32 c, uint8 v ) const
        {
            uint32 c_index ;
            const ElementType& type = cell_type( c, c_index ) ;
            return vertex_index( cells_[type][c_index] + v ) ;
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

        //     ___     _ _               _
        //    / __|___| | | __ _____ _ _| |_ _____ __
        //   | (__/ -_) | | \ V / -_) '_|  _/ -_) \ /
        //    \___\___|_|_|  \_/\___|_|  \__\___/_\_\
        //
        const vec3& cell_vertex( uint32 c, uint8 v ) const {
            return vertex( cell_vertex_index( c, v ) ) ;
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

        //     ___     _ _    __             _     _                              _
        //    / __|___| | |  / _|__ _ __ ___| |_  | |__  __ _ _ _ _  _ __ ___ _ _| |_ ___ _ _
        //   | (__/ -_) | | |  _/ _` / _/ -_)  _| | '_ \/ _` | '_| || / _/ -_) ' \  _/ -_) '_|
        //    \___\___|_|_| |_| \__,_\__\___|\__| |_.__/\__,_|_|  \_, \__\___|_||_\__\___|_|
        //                                                        |__/
        vec3 cell_facet_barycenter( uint32 c, uint8 f ) const ;
        vec3 triangle_facet_barycenter( uint32 c, uint8 f ) const ;
        vec3 quad_facet_barycenter( uint32 c, uint8 f ) const ;
        vec3 tetra_facet_barycenter( uint32 c, uint8 f ) const ;
        vec3 pyramid_facet_barycenter( uint32 c, uint8 f ) const ;
        vec3 prism_facet_barycenter( uint32 c, uint8 f ) const ;
        vec3 hexa_facet_barycenter( uint32 c, uint8 f ) const ;

        //     ___     _ _   _                              _
        //    / __|___| | | | |__  __ _ _ _ _  _ __ ___ _ _| |_ ___ _ _
        //   | (__/ -_) | | | '_ \/ _` | '_| || / _/ -_) ' \  _/ -_) '_|
        //    \___\___|_|_| |_.__/\__,_|_|  \_, \__\___|_||_\__\___|_|
        //                                  |__/
        vec3 cell_barycenter( uint32 c ) const ;
        vec3 triangle_barycenter( uint32 c ) const ;
        vec3 quad_barycenter( uint32 c ) const ;
        vec3 tetra_barycenter( uint32 c ) const ;
        vec3 pyramid_barycenter( uint32 c ) const ;
        vec3 prism_barycenter( uint32 c ) const ;
        vec3 hexa_barycenter( uint32 c ) const ;

        //    _  _                     _               _     _
        //   | \| |___ __ _ _ _ ___ __| |_   _ __  ___(_)_ _| |_
        //   | .` / -_) _` | '_/ -_|_-<  _| | '_ \/ _ \ | ' \  _|
        //   |_|\_\___\__,_|_| \___/__/\__| | .__/\___/_|_||_\__|
        //                                  |_|
        float64 get_nearest_point_in_cell( const vec3& p, uint32 c, vec3& nearest_p ) const ;
        float64 get_nearest_point_in_line( const vec3& p, uint32 c, vec3& nearest_p ) const ;
        float64 get_nearest_point_in_triangle( const vec3& p, uint32 c, vec3& nearest_p ) const ;
        float64 get_nearest_point_in_quad( const vec3& p, uint32 c, vec3& nearest_p ) const ;
        float64 get_nearest_point_in_tetra( const vec3& p, uint32 c, vec3& nearest_p ) const ;
        float64 get_nearest_point_in_pyramid( const vec3& p, uint32 c, vec3& nearest_p ) const ;
        float64 get_nearest_point_in_prism( const vec3& p, uint32 c, vec3& nearest_p ) const ;
        float64 get_nearest_point_in_hexa( const vec3& p, uint32 c, vec3& nearest_p ) const ;

        //      _  _   _       _ _         _
        //     /_\| |_| |_ _ _(_) |__ _  _| |_ ___   _ __  __ _ _ _  __ _ __ _ ___ _ _
        //    / _ \  _|  _| '_| | '_ \ || |  _/ -_) | '  \/ _` | ' \/ _` / _` / -_) '_|
        //   /_/ \_\__|\__|_| |_|_.__/\_,_|\__\___| |_|_|_\__,_|_||_\__,_\__, \___|_|
        //                                                               |___/
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

        //      _      _  _
        //     /_\  __| |(_)__ _ __ ___ _ _  __ _  _
        //    / _ \/ _` || / _` / _/ -_) ' \/ _| || |
        //   /_/ \_\__,_|/ \__,_\__\___|_||_\__|\_, |
        //             |__/                     |__/
        void compute_adjacencies() {
            //todo
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

        AdjMatrix adjacencies_ ;
    } ;

    template< class ATTRIBUTE >
    class GRGMESH_API LineAttribute: public Attribute< LINE, ATTRIBUTE > {
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
    class GRGMESH_API TriangleAttribute: public Attribute< TRGL, ATTRIBUTE > {
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
    class GRGMESH_API QuadAttribute: public Attribute< QUAD, ATTRIBUTE > {
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
    class GRGMESH_API TetraAttribute: public Attribute< TETRA, ATTRIBUTE > {
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
    class GRGMESH_API PyramidAttribute: public Attribute< PYRAMID, ATTRIBUTE > {
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
    class GRGMESH_API PrismAttribute: public Attribute< PRISM, ATTRIBUTE > {
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
    class GRGMESH_API HexaAttribute: public Attribute< HEXA, ATTRIBUTE > {
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

        //    ___      _   _
        //   / __| ___| |_| |_ ___ _ _ ___
        //   \__ \/ -_)  _|  _/ -_) '_(_-<
        //   |___/\___|\__|\__\___|_| /__/
        //
        void set_cell( const ElementType& type, uint32 id, uint32 index ) {
            mixed_mesh_.cells_[type][id] = index ;
        }
        void set_line( uint32 id, uint32 index ) {
            mixed_mesh_.cells_[LINE][id] = index ;
        }
        void set_triangle( uint32 id, uint32 index ) {
            mixed_mesh_.cells_[TRGL][id] = index ;
        }
        void set_quad( uint32 id, uint32 index ) {
            mixed_mesh_.cells_[QUAD][id] = index ;
        }
        void set_tetra( uint32 id, uint32 index ) {
            mixed_mesh_.cells_[TETRA][id] = index ;
        }
        void set_pyramid( uint32 id, uint32 index ) {
            mixed_mesh_.cells_[PYRAMID][id] = index ;
        }
        void set_prism( uint32 id, uint32 index ) {
            mixed_mesh_.cells_[PRISM][id] = index ;
        }
        void set_hexa( uint32 id, uint32 index ) {
            mixed_mesh_.cells_[HEXA][id] = index ;
        }


    private:
        MixedMesh& mixed_mesh_ ;
    } ;

    class GRGMESH_API MixedMeshBuilder: public MeshBuilder {
    public:
        MixedMeshBuilder( MixedMesh& mixed_mesh )
            : MeshBuilder( mixed_mesh ), mixed_mesh_( mixed_mesh )
        {
        }
        MixedMeshBuilder( const MixedMesh& mixed_mesh )
            :
                MeshBuilder( mixed_mesh ),
                mixed_mesh_( const_cast< MixedMesh& >( mixed_mesh ) )
        {
        }
        virtual ~MixedMeshBuilder() {}

        //    ___
        //   | _ \___ ______ _ ___ _____
        //   |   / -_|_-< -_) '_\ V / -_)
        //   |_|_\___/__|___|_|  \_/\___|
        //
        void reserve_cells( const ElementType& type, uint32 nb ) {
            mixed_mesh_.cells_[type].reserve( nb ) ;
        }
        void reserve_lines( uint32 nb ) {
            mixed_mesh_.cells_[LINE].reserve( nb ) ;
        }
        void reserve_triangles( uint32 nb ) {
            mixed_mesh_.cells_[TRGL].reserve( nb ) ;
        }
        void reserve_quad( uint32 nb ) {
            mixed_mesh_.cells_[QUAD].reserve( nb ) ;
        }
        void reserve_tetra( uint32 nb ) {
            mixed_mesh_.cells_[TETRA].reserve( nb ) ;
        }
        void reserve_pyramids( uint32 nb ) {
            mixed_mesh_.cells_[PYRAMID].reserve( nb ) ;
        }
        void reserve_prisms( uint32 nb ) {
            mixed_mesh_.cells_[PRISM].reserve( nb ) ;
        }
        void reserve_hexa( uint32 nb ) {
            mixed_mesh_.cells_[HEXA].reserve( nb ) ;
        }

        //    ___        _
        //   | _ \___ __(_)______
        //   |   / -_|_-< |_ / -_)
        //   |_|_\___/__/_/__\___|
        //
        void resize_cells( const ElementType& type, uint32 nb, uint32 id = dummy_uint32 ) {
            mixed_mesh_.cells_[type].resize( nb, id ) ;
        }
        void resize_lines( uint32 nb, uint32 id = dummy_uint32 ) {
            mixed_mesh_.cells_[LINE].resize( nb, id ) ;
        }
        void resize_triangles( uint32 nb, uint32 id = dummy_uint32 ) {
            mixed_mesh_.cells_[TRGL].resize( nb ) ;
        }
        void resize_quad( uint32 nb, uint32 id = dummy_uint32 ) {
            mixed_mesh_.cells_[QUAD].resize( nb, id ) ;
        }
        void resize_tetra( uint32 nb, uint32 id = dummy_uint32 ) {
            mixed_mesh_.cells_[TETRA].resize( nb, id ) ;
        }
        void resize_pyramids( uint32 nb, uint32 id = dummy_uint32 ) {
            mixed_mesh_.cells_[PYRAMID].resize( nb, id ) ;
        }
        void resize_prisms( uint32 nb, uint32 id = dummy_uint32 ) {
            mixed_mesh_.cells_[PRISM].resize( nb, id ) ;
        }
        void resize_hexa( uint32 nb, uint32 id = dummy_uint32 ) {
            mixed_mesh_.cells_[HEXA].resize( nb, id ) ;
        }

        //      _      _    _          _ _   _         _
        //     /_\  __| |__| |  __ ___| | | (_)_ _  __| |_____ __
        //    / _ \/ _` / _` | / _/ -_) | | | | ' \/ _` / -_) \ /
        //   /_/ \_\__,_\__,_| \__\___|_|_| |_|_||_\__,_\___/_\_\
        //
        void add_cell( const ElementType& type, uint32 index ) {
            mixed_mesh_.cells_[type].push_back( index ) ;
        }
        void add_line( uint32 index ) {
            mixed_mesh_.cells_[LINE].push_back( index ) ;
        }
        void add_triangle( uint32 index ) {
            mixed_mesh_.cells_[TRGL].push_back( index ) ;
        }
        void add_quad( uint32 index ) {
            mixed_mesh_.cells_[QUAD].push_back( index ) ;
        }
        void add_tetra( uint32 index ) {
            mixed_mesh_.cells_[TETRA].push_back( index ) ;
        }
        void add_pyramid( uint32 index ) {
            mixed_mesh_.cells_[PYRAMID].push_back( index ) ;
        }
        void add_prism( uint32 index ) {
            mixed_mesh_.cells_[PRISM].push_back( index ) ;
        }
        void add_hexa( uint32 index ) {
            mixed_mesh_.cells_[HEXA].push_back( index ) ;
        }

        //      _      _    _          _ _
        //     /_\  __| |__| |  __ ___| | |
        //    / _ \/ _` / _` | / _/ -_) | |
        //   /_/ \_\__,_\__,_| \__\___|_|_|
        //
        uint32 add_cell( const ElementType& type ) {
            uint32 index = mixed_mesh_.vertex_indices_.size() ;
            mixed_mesh_.cells_[type].push_back( index ) ;
            return index ;
        }
        uint32 add_line() {
            uint32 index = mixed_mesh_.vertex_indices_.size() ;
            mixed_mesh_.cells_[LINE].push_back( index ) ;
            return index ;
        }
        uint32 add_triangle() {
            uint32 index = mixed_mesh_.vertex_indices_.size() ;
            mixed_mesh_.cells_[TRGL].push_back( index ) ;
            return index ;
        }
        uint32 add_quad() {
            uint32 index = mixed_mesh_.vertex_indices_.size() ;
            mixed_mesh_.cells_[QUAD].push_back( index ) ;
            return index ;
        }
        uint32 add_tetra() {
            uint32 index = mixed_mesh_.vertex_indices_.size() ;
            mixed_mesh_.cells_[TETRA].push_back( index ) ;
            return index ;
        }
        uint32 add_pyramid() {
            uint32 index = mixed_mesh_.vertex_indices_.size() ;
            mixed_mesh_.cells_[PYRAMID].push_back( index ) ;
            return index ;
        }
        uint32 add_prism() {
            uint32 index = mixed_mesh_.vertex_indices_.size() ;
            mixed_mesh_.cells_[PRISM].push_back( index ) ;
            return index ;
        }
        uint32 add_hexa() {
            uint32 index = mixed_mesh_.vertex_indices_.size() ;
            mixed_mesh_.cells_[HEXA].push_back( index ) ;
            return index ;
        }

        //    ___      _            _ _   _         _
        //   / __| ___| |_   __ ___| | | (_)_ _  __| |_____ __
        //   \__ \/ -_)  _| / _/ -_) | | | | ' \/ _` / -_) \ /
        //   |___/\___|\__| \__\___|_|_| |_|_||_\__,_\___/_\_\
        //
        void add_cell( const ElementType& type, uint32 id, uint32 index ) {
            mixed_mesh_.cells_[type][id] = index ;
        }
        void add_line( uint32 id, uint32 index ) {
            mixed_mesh_.cells_[LINE][id] = index ;
        }
        void add_triangle( uint32 id, uint32 index ) {
            mixed_mesh_.cells_[TRGL][id] = index ;
        }
        void add_quad( uint32 id, uint32 index ) {
            mixed_mesh_.cells_[QUAD][id] = index ;
        }
        void add_tetra( uint32 id, uint32 index ) {
            mixed_mesh_.cells_[TETRA][id] = index ;
        }
        void add_pyramid( uint32 id, uint32 index ) {
            mixed_mesh_.cells_[PYRAMID][id] = index ;
        }
        void add_prism( uint32 id, uint32 index ) {
            mixed_mesh_.cells_[PRISM][id] = index ;
        }
        void add_hexa( uint32 id, uint32 index ) {
            mixed_mesh_.cells_[HEXA][id] = index ;
        }

    private:
        MixedMesh& mixed_mesh_ ;
    } ;

}

#endif
