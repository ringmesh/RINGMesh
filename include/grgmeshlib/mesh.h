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

#ifndef __GRGMESH_MESH__
#define __GRGMESH_MESH__

#include <grgmeshlib/common.h>
#include <grgmeshlib/vecn.h>
#include <grgmeshlib/attribute.h>

#include <vector>
#include <map>

namespace GRGMesh {

    class GRGMESH_API Mesh {
        friend class MeshMutator ;
        friend class MeshBuilder ;
        typedef AttributeManager< VERTEX > VertexAttributeManager ;

    public:
        void clear()
        {
            vertices_.clear() ;
            vertex_indices_.clear() ;
        }

        uint32 nb_vertices() const
        {
            return vertices_.size() ;
        }
        virtual uint32 nb_cells() const = 0 ;

        virtual uint32 cell_begin( uint32 c ) const = 0 ;
        virtual uint32 cell_end( uint32 c ) const = 0 ;
        virtual uint8 nb_vertices_in_cell( uint32 c ) const = 0 ;
        virtual uint8 nb_facets_in_cell( uint32 c ) const = 0 ;
        virtual uint8 nb_vertices_in_cell_facet( uint32 c, uint8 f ) const = 0 ;

        uint32 cell_vertex_index( uint32 c, uint32 v ) const
        {
            grgmesh_debug_assert( v < nb_vertices_in_cell( c ) ) ;
            return vertex_indices_[cell_begin( c ) + v] ;
        }
        const vec3& cell_vertex( uint32 c, uint32 v ) const
        {
            grgmesh_debug_assert( v < nb_vertices_in_cell( c ) ) ;
            return vertices_[vertex_indices_[cell_begin( c ) + v]] ;
        }
        const vec3& vertex( uint32 v ) const { return vertices_[v] ; }
        uint32 vertex_index( uint32 i ) const { return vertex_indices_[i] ; }

        virtual ElementType cell_type( uint32 c, uint32& c_index = dummy_uint32  ) const = 0 ;
        virtual const CellDescriptor* cell_descriptor( uint32 c ) const = 0 ;
        uint32 cell_facet_vertex( uint32 c, uint32 f, uint32 v ) const
        {
            const CellDescriptor* desc = cell_descriptor( c ) ;
            grgmesh_debug_assert( f < desc->nb_facets ) ;
            grgmesh_debug_assert( v < desc->nb_vertices_in_facet[f] ) ;
            return cell_vertex_index( c, desc->facet[f][v] ) ;
        }

        vec3 cell_centroid( uint32 c ) const
        {
            vec3 result ;
            for( uint32 i = cell_begin( c ); i < cell_end( c ); i++ ) {
                result += cell_vertex( c, i ) ;
            }
            return result / nb_vertices_in_cell( c );
        }

        VertexAttributeManager* vertex_attribute_manager() const
        {
            return const_cast< VertexAttributeManager* >( &vertex_attribute_manager_ ) ;
        }

    protected:
        Mesh()
        {
        }
        virtual ~Mesh()
        {
        }
        void copy( const Mesh& rhs )
        {
        }

    protected:
        ///List of the vertices
        std::vector< vec3 > vertices_ ;
        ///Mapping between the list of the vertices in a cell and the actual vertices
        std::vector< uint32 > vertex_indices_ ;

    private:
        VertexAttributeManager vertex_attribute_manager_ ;

    } ;

    template< class ATTRIBUTE >
    class VertexAttribute: public Attribute< VERTEX, ATTRIBUTE > {
    public:
        typedef Attribute< VERTEX, ATTRIBUTE > superclass ;

        void bind( Mesh* mesh, const std::string& name )
        {
            superclass::bind( mesh->vertex_attribute_manager(), mesh->nb_vertices(),
                name ) ;
        }

        void bind( Mesh* mesh )
        {
            superclass::bind( mesh->vertex_attribute_manager(),
                mesh->nb_vertices() ) ;
        }

        VertexAttribute()
        {
        }

        VertexAttribute( Mesh* mesh )
        {
            bind( mesh ) ;
        }

        VertexAttribute( Mesh* mesh, const std::string& name )
        {
            bind( mesh, name ) ;
        }

        static bool is_defined( Mesh* mesh, const std::string& name )
        {
            return superclass::is_defined( mesh->vertex_attribute_manager(), name ) ;
        }
    } ;



    class GRGMESH_API MeshMutator {
    public:
        MeshMutator( Mesh& mesh )
            : mesh_( mesh )
        {
        }
        MeshMutator( const Mesh& mesh )
            : mesh_( const_cast< Mesh& >( mesh ) )
        {
        }
        virtual ~MeshMutator() {}

        void set_vertex( uint32 id, const vec3& v )
        {
            mesh_.vertices_[id] = v ;
        }
        void set_vertex_index( uint32 id, uint32 v )
        {
            mesh_.vertex_indices_[id] = v ;
        }

    protected:
        Mesh& mesh_ ;
    } ;

    class GRGMESH_API MeshBuilder {
    public:
        MeshBuilder( Mesh& mesh )
            : mesh_( mesh )
        {
        }
        MeshBuilder( const Mesh& mesh )
            : mesh_( const_cast< Mesh& >( mesh ) )
        {
        }
        virtual ~MeshBuilder()
        {
        }

        void reserve_vertices( uint32 nb )
        {
            mesh_.vertices_.reserve( nb ) ;
        }
        void resize_vertices( uint32 nb, const vec3& v = dummy_vec3 )
        {
            mesh_.vertices_.resize( nb, v ) ;
        }
        void add_vertex( const vec3& v )
        {
            mesh_.vertices_.push_back( v ) ;
        }
        void add_vertex( uint32 id, const vec3& v )
        {
            mesh_.vertices_[id] = v ;
        }

        void reserve_vertex_indices( uint32 nb )
        {
            mesh_.vertex_indices_.reserve( nb ) ;
        }
        void resize_vertex_indices( uint32 nb, uint32 id = dummy_uint32 )
        {
            mesh_.vertex_indices_.resize( nb, id ) ;
        }
        void add_vertex_index( uint32 v )
        {
            mesh_.vertex_indices_.push_back( v ) ;
        }
        void add_vertex_index( uint32 id, uint32 v )
        {
            mesh_.vertex_indices_[id] = v ;
        }

    protected:
        Mesh& mesh_ ;
    } ;

}

#endif
