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

#include <ringmesh/geo_model_api.h>
#include <ringmesh/geo_model.h>
#include <ringmesh/geometry.h>
#include <ringmesh/geogram_extension.h>
#include <ringmesh/tetra_gen.h>

#include <geogram/basic/logger.h>
#include <geogram/basic/progress.h>
#include <geogram/basic/geometry_nd.h>
#include <geogram/mesh/mesh_geometry.h>

#include <iostream>
#include <iomanip>

namespace RINGMesh {
    /*!
    * @brief Total number of facets in the geomodel Surfaces
    */
    index_t nb_facets( const GeoModel& geomodel )
    {
        index_t result = 0 ;
        for( index_t i = 0; i < geomodel.nb_surfaces(); ++i ) {
            result += geomodel.surface( i ).nb_cells() ;
        }
        return result ;
    }

    void print_model( const GeoModel& geomodel )
    {
        GEO::Logger::out( "GeoModel" ) << "Model " << geomodel.name() << " has "
            << std::endl
            << std::setw( 10 ) << std::left
            << geomodel.mesh.vertices.nb() << " vertices "
            << std::endl
            << std::setw( 10 ) << std::left
            << nb_facets( geomodel ) << " facets "
            << std::endl << std::endl ;

        for( index_t t = GME::CORNER; t < GME::NO_TYPE; ++t ) {
            GME::TYPE T = static_cast<GME::TYPE>( t ) ;
            GEO::Logger::out( "GeoModel" ) << std::setw( 10 ) << std::left
                << geomodel.nb_elements( T ) << " " << GME::type_name( T )
                << std::endl ;
        }
    }

	///@todo A class encapsulating the copy from a GeoModel to a Mesh ?
    /// See what has been done in GeoModelMeshBuilder	

    void add_geomodel_vertices_to_mesh( const GeoModel& geomodel, GEO::Mesh& M )
    {
        index_t nbv = geomodel.mesh.vertices.nb() ;
        M.vertices.create_vertices( nbv ) ;

        // We need to copy the point one after another since we do not have access
        // to the storage of the geomodel.vertices
        // I do not want to provide this access [JP]
        for( index_t v = 0; v < nbv; ++v ) {
            M.vertices.point( v ) = geomodel.mesh.vertices.vertex( v ) ;
        }

		GEO::Attribute<index_t> corner_attribute(M.vertices.attributes(), "region");
		corner_attribute.fill(NO_ID);
		for (index_t i = 0; i < geomodel.nb_corners(); ++i){
			index_t vertex_index = geomodel.corner(i).model_vertex_id();
			corner_attribute[vertex_index] = i;
		}
		corner_attribute.unbind();
    }

	void add_line_edges_to_mesh(const Line& line, GEO::Mesh& M)
	{
		index_t from = M.edges.create_edges(line.nb_cells());
		for (index_t i = 0; i < line.nb_cells(); ++i){
			index_t v0 = line.model_vertex_id(i, 0);
			index_t v1 = line.model_vertex_id(i, 1);
			M.edges.set_vertex(from + i, 0, v0);
			M.edges.set_vertex(from + i, 1, v1);
		}
	}

    void create_and_fill_line_index_attribute(const GeoModel& geomodel,
        const std::string& attribute_name,
        GEO::Mesh& M)
    {
        GEO::Attribute<index_t> line_attribute(M.edges.attributes(), "region");
        line_attribute.fill(NO_ID);
        index_t edge_counter = 0;
        for (index_t i = 0; i < geomodel.nb_lines(); ++i){
            index_t nb_line_edges = geomodel.line(i).nb_cells();
            index_t line_edges_start = edge_counter;
            index_t line_edges_end = edge_counter + nb_line_edges;

            for (index_t e = line_edges_start; e != line_edges_end; ++e) {
                line_attribute[e] = i;
            }
            edge_counter += nb_line_edges;
        }
        line_attribute.unbind();
    }

	void add_geomodel_line_edges_to_mesh(const GeoModel& geomodel, GEO::Mesh& M) 
	{
		for (index_t i = 0; i < geomodel.nb_lines(); ++i) {
			const Line& line(geomodel.line(i)) ;
			add_line_edges_to_mesh(line, M);
		}
		
        create_and_fill_line_index_attribute(geomodel, "region", M);
	}

    void add_surface_facets_to_mesh( const Surface& surface, GEO::Mesh& M )
    {
        for( index_t j = 0; j < surface.nb_cells(); ++j ) {
            index_t nbv = surface.nb_vertices_in_facet( j ) ;
            GEO::vector< index_t > ids( nbv ) ;
            for( index_t v = 0; v < nbv; ++v ) {
                ids[ v ] = surface.model_vertex_id( j, v ) ;
            }
            M.facets.create_polygon( ids ) ;
        }    
    }

    void add_surface_triangles_to_mesh( const Surface& surface, GEO::Mesh& M )
    {
        GEO::vector< index_t > triangles( 3*surface.nb_cells() ) ;
        for( index_t j = 0; j < surface.nb_cells(); ++j ) {      
            for( index_t v = 0; v < 3; ++v ) {
                triangles[ 3*j+v ] = surface.model_vertex_id( j, v ) ;
            }           
        }    
        M.facets.assign_triangle_mesh( triangles, true ) ;
    }

	void create_and_fill_surface_index_attribute(const GeoModel& geomodel, 
												 const std::string& attribute_name,
												 GEO::Mesh& M )
	{
		GEO::Attribute<index_t> surface_attribute(M.facets.attributes(), attribute_name);
		surface_attribute.fill(NO_ID);
		index_t facet_counter = 0;
		for (index_t i = 0; i < geomodel.nb_surfaces(); ++i){
			index_t nb_surface_facets = geomodel.surface(i).nb_cells();
			index_t surface_facet_start = facet_counter;
			index_t surface_facet_end = facet_counter + nb_surface_facets;

			for (index_t f = surface_facet_start; f != surface_facet_end; ++f) {
				surface_attribute[f] = i;
			}
            facet_counter += nb_surface_facets;
		}
		surface_attribute.unbind();
	}

    void add_geomodel_surface_facets_to_mesh( const GeoModel& geomodel, GEO::Mesh& M )
    {
        for( index_t i = 0; i < geomodel.nb_surfaces(); ++i ) {
            const Surface& S = geomodel.surface( i ) ;
            if( S.is_simplicial() ) {
                add_surface_triangles_to_mesh( S, M ) ;
            } else {
                add_surface_facets_to_mesh( S, M ) ;
            }
        }
        create_and_fill_surface_index_attribute(geomodel, "region", M);
    }

    void add_region_tets_to_mesh( const Region& region, GEO::Mesh& M )
    {
        GEO::vector< index_t > tets( region.nb_cells()*4 ) ;
        for( index_t j = 0; j < region.nb_cells(); ++j ) {
            for( index_t v = 0; v < 4; ++v ) {
                tets[ 4*j+v ] = region.model_vertex_id( j, v ) ;
            }
        }
        M.cells.assign_tet_mesh( tets, true ) ;
    }

    void create_and_fill_region_index_attribute(const GeoModel& geomodel,
        const std::string& attribute_name,
        GEO::Mesh& M)
    {
        GEO::Attribute<index_t> region_attribute(M.cells.attributes(), attribute_name);
        region_attribute.fill(NO_ID);
        index_t cell_counter = 0;
        for (index_t i = 0; i < geomodel.nb_regions(); ++i){
            index_t nb_region_cells = geomodel.region(i).nb_cells();
            index_t region_cells_start = cell_counter;
            index_t region_cells_end = cell_counter + nb_region_cells;

            for (index_t j = region_cells_start; j != region_cells_end; ++j) {
                region_attribute[j] = i;
            }
            cell_counter += nb_region_cells;
        }
        region_attribute.unbind();
    }

    /*! @todo to implement for other types of cells
     */
    void add_geomodel_region_tets_to_mesh( const GeoModel& geomodel, GEO::Mesh& M )
    {
        for( index_t i = 0; i < geomodel.nb_regions(); ++i ) {
            const Region& region = geomodel.region( i ) ;            
            if( region.is_simplicial() ) {
                add_region_tets_to_mesh( region, M ) ;
            } else {
                ringmesh_assert( false ) ;
            }           
        }
        create_and_fill_region_index_attribute(geomodel, "region", M);
    }


    void build_mesh_from_geomodel( const GeoModel& geomodel, GEO::Mesh& M ) 
    {
        // Keep the attributes when clearing the mesh, otherwise we crash
        M.clear( true ) ;

		add_geomodel_vertices_to_mesh(geomodel, M);
		add_geomodel_line_edges_to_mesh(geomodel, M);
        add_geomodel_surface_facets_to_mesh(geomodel, M);
        add_geomodel_region_tets_to_mesh(geomodel, M);
    }


    double model_element_size( const GeoModelElement& E )
    {
        double result = 0. ;
        if( E.nb_children() ) {
            // Sum up the size of children elements
            for( index_t i = 0; i < E.nb_children(); ++i ) {
                result += model_element_size( E.child( i ) ) ;
            }
            return result ;
        } else {
            switch( E.type() ) {
                case GeoModelElement::REGION: {
                    const Region& R = dynamic_cast< const Region& >( E ) ;
                    // Compute the volume if this is a region
                    for( index_t i = 0; i < R.nb_boundaries(); i++ ) {
                        const Surface& surface =
                            dynamic_cast< const Surface& >( R.boundary( i ) ) ;

                        for( index_t t = 0; t < surface.nb_cells(); t++ ) {
                            const vec3& p0 = surface.vertex( t, 0 ) ;
                            for( index_t v = 1;
                                v + 1 < surface.nb_vertices_in_facet( t ); ++v ) {
                                double cur_volume = ( dot( p0,
                                    cross( surface.vertex( t, v ),
                                        surface.vertex( t, v + 1 ) ) ) )
                                    / static_cast< double >( 6 ) ;
                                R.side( i ) ? result -= cur_volume : result +=
                                                  cur_volume ;
                            }
                        }
                    }
                    return fabs( result ) ;
                }

                case GeoModelElement::SURFACE: {
                    const Surface& S = dynamic_cast< const Surface& >( E ) ;
                    const GEO::Mesh& mesh = S.mesh() ;
                    for( index_t i = 0; i < S.nb_cells(); i++ ) {
                        result += GEO::Geom::mesh_facet_area( mesh, i ) ;
                    }
                    return result ;
                }
                case GeoModelElement::LINE: {
                    const Line& L = dynamic_cast< const Line& >( E ) ;
                    for( index_t i = 1; i < L.nb_vertices(); ++i ) {
                        result += GEO::Geom::distance( L.vertex( i ),
                                                       L.vertex( i - 1 ) ) ;
                    }
                    return result ;
                }
                case GeoModelElement::CORNER: {
                    return 0 ;
                }
            }
            ringmesh_assert_not_reached;
            return result ;
        }
    }

    double model_element_cell_size( const GeoModelElement& E, index_t c )
    {
        double result = 0. ;

        switch( E.type() ) {
            case GeoModelElement::REGION: {
                const Region& R = dynamic_cast< const Region& >( E ) ;
                const GEO::Mesh& mesh = R.mesh() ;
                return RINGMesh::mesh_cell_volume( mesh, c ) ;
            }
            case GeoModelElement::SURFACE: {
                const Surface& S = dynamic_cast< const Surface& >( E ) ;
                const GEO::Mesh& mesh = S.mesh() ;
                return GEO::Geom::mesh_facet_area( mesh, c ) ;
            }
            case GeoModelElement::LINE: {
                const Line& L = dynamic_cast< const Line& >( E ) ;
                const GEO::Mesh& mesh = L.mesh() ;
                return GEO::Geom::distance( L.vertex( c, 0 ), L.vertex( c, 1 ) ) ;
            }
        }
        ringmesh_assert_not_reached ;
        return result ;
    }

    vec3 model_element_center( const GeoModelElement& E )
    {
        vec3 result( 0., 0., 0. ) ;
        index_t nb_vertices = 0 ;

        if( GeoModelElement::has_mesh( E.type() ) ) {
            // @todo Improve efficiency, overload the functions to avoid
            // casting each time
            const GeoModelMeshElement& M =
                dynamic_cast< const GeoModelMeshElement& >( E ) ;
            for( index_t v = 0; v < M.nb_vertices(); v++ ) {
                result += M.vertex( v ) ;
            }
            return result / static_cast< double >( M.nb_vertices() ) ;
        }
        else if( E.nb_children() > 0 ) {
            for( index_t i = 0; i < E.nb_children(); ++i ) {
                const GeoModelMeshElement& F =
                    dynamic_cast< const GeoModelMeshElement& >( E.child( i ) ) ;
                nb_vertices += F.nb_vertices() ;
                result += model_element_center( F ) * F.nb_vertices() ;
            }
            return result / static_cast< double >( nb_vertices ) ;
        } else {
            return result ;
        }
    }

    vec3 model_element_cell_center( const GeoModelMeshElement& E, index_t c )
    {
        vec3 result( 0., 0., 0. ) ;       
        const GEO::Mesh& mesh = E.mesh() ;
        switch( E.type() ) {
            case GeoModelElement::REGION: {                
                return RINGMesh::mesh_cell_center( mesh, c ) ;
            }
            case GeoModelElement::SURFACE: {             
                return GEO::Geom::mesh_facet_center( mesh, c ) ;
            }
            case GeoModelElement::LINE: {
                index_t v0 = mesh.edges.vertex( c, 0 ) ;
                index_t v1 = mesh.edges.vertex( c, 1 ) ;
                return 0.5 * ( mesh.vertices.point( v0 ) + mesh.vertices.point( v1 ) ) ;
            }
            case GeoModelElement::CORNER: {
                return mesh.vertices.point(0) ;
            }
        }
        ringmesh_assert_not_reached ;
        return result ;
    }

    void translate( GeoModel& M, const vec3& translation_vector )
    {
        for( index_t v = 0; v < M.mesh.vertices.nb(); ++v ) {
            // Coordinates are not directly modified to 
            // update the matching vertices in geomodel entities
            vec3 p = M.mesh.vertices.vertex( v ) ;
            for( index_t i = 0; i < 3; i++ ) {
                p[i] += translation_vector[i] ;
            }
            M.mesh.vertices.update_point( v, p ) ;
        }
    }

    void rotate( GeoModel& M, const vec3& origin, const vec3& axis, double theta, bool degrees )
    {
        if( length( axis ) < epsilon ) {
            GEO::Logger::err( "GeoModel" )
                << "Rotation around an epsilon length axis is impossible"
                << std::endl ;
            return ;
        }
        GEO::Matrix< double, 4 > rot_mat ;
        rotation_matrix_about_arbitrary_axis( origin, axis, theta, degrees, rot_mat ) ;

        for( index_t v = 0; v < M.mesh.vertices.nb(); ++v ) {
            const vec3& p = M.mesh.vertices.vertex( v ) ;
            double old[4] = { p[0], p[1], p[2], 1. } ;
            double new_p[4] = { 0, 0, 0, 1. } ;
            GEO::mult( rot_mat, old, new_p ) ;
            /*! @todo You need an epsilon tolerance here [JP] */
            ringmesh_debug_assert( new_p[3] == 1. ) ;

            M.mesh.vertices.update_point( v, vec3( new_p[0], new_p[1], new_p[2] ) ) ;
        }
    }

    void tetrahedralize( GeoModel& M, const std::string& method, index_t region_id, bool add_steiner_points )
    {
        /* @todo Review: Maybe rethink these functions
         *       to have a function that can mesh a region of a geomodel
         *       taking only one vector of points [JP]
         */
        std::vector< std::vector< vec3 > > internal_vertices( M.nb_regions() ) ;
        tetrahedralize( M, method, region_id, add_steiner_points,
            internal_vertices ) ;
    }

    void tetrahedralize( GeoModel& M, const std::string& method, index_t region_id, bool add_steiner_points,
                         const std::vector< std::vector< vec3 > >& internal_vertices )
    {
        if( region_id == NO_ID ) {
            GEO::Logger::out( "Info" ) << "Using " << method << std::endl ;
            GEO::ProgressTask progress( "Compute", M.nb_regions() ) ;
            for( index_t i = 0; i < M.nb_regions(); i++ ) {
                tetrahedralize( M, method, i, add_steiner_points, internal_vertices ) ;                
                progress.next() ;
            }
        } else {
            TetraGen_var tetragen = TetraGen::create( M.region( region_id ).mesh(),
                method ) ;
            tetragen->set_boundaries( M.region( region_id ), M.wells() ) ;
            tetragen->set_internal_points( internal_vertices[region_id] ) ;
            GEO::Logger::instance()->set_quiet( true ) ;
            tetragen->tetrahedralize( add_steiner_points ) ;
            GEO::Logger::instance()->set_quiet( false ) ;
        }

        // The GeoModelMesh should be updated, just erase everything
        // and it will be re-computed during its next access.
        M.mesh.vertices.clear() ;
    }

}
