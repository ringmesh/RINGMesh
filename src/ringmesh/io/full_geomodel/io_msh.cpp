/*
 * Copyright (c) 2012-2017, Association Scientifique pour la Geologie et ses Applications (ASGA)
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
 *     http://www.ring-team.org
 *
 *     RING Project
 *     Ecole Nationale Superieure de Geologie - GeoRessources
 *     2 Rue du Doyen Marcel Roubault - TSA 70605
 *     54518 VANDOEUVRE-LES-NANCY
 *     FRANCE
 */

namespace {
//        struct RINGMesh2GMSH {
//                   index_t entity_type ;
//                   index_t nb_vertices ;
//                   index_t vertices[8] ;
//                   index_t nb_facets ;
//                   index_t nb_vertices_in_facet[6] ;
//                   index_t facet[6] ;
//                   index_t vertices_in_facet[6][4] ;
//               } ;
//
//               static RINGMesh2GMSH tet_descriptor_gmsh = { 4,                  // type
//                   4,                  // nb vertices
//                   { 0, 1, 2, 3, 5, 8 ,9 ,4,6,7 },     // vertices
//                   4,                  // nb facets
//                   { 3, 3, 3, 3 },     // nb vertices in facet
//                   { 0, 1, 2, 3 },     // facets
//                   { { 1, 3, 2 }, { 0, 2, 3 }, { 3, 1, 0 }, { 0, 1, 2 } } } ;
//
//               static RINGMesh2GMSH hex_descriptor_gmsh = { 6,                         // type
//                   8,                              // nb vertices
//                   { 4, 0, 5, 1, 7, 3, 6, 2 },     // vertices
//                   6,                              // nb facets
//                   { 4, 4, 4, 4, 4, 4 },           // nb vertices in facet
//                   { 4, 2, 1, 3, 0, 5 },           // facets
//                   { { 0, 3, 7, 4 }, { 2, 1, 5, 6 }, { 1, 0, 4, 5 }, { 3, 2, 6, 7 }, {
//                       1, 2, 3, 0 }, { 4, 7, 6, 5 } } } ;
//
//               static RINGMesh2GMSH prism_descriptor_gmsh = { 12,                     // type
//                   6,                      // nb vertices
//                   { 0, 1, 2, 3, 4, 5 },   // vertices
//                   5,                      // nb facets
//                   { 3, 4, 4, 4, 3 },      // nb vertices in facet
//                   { 0, 2, 4, 3, 1 },      // facets
//                   {
//                       { 0, 1, 2 }, { 3, 5, 4 }, { 0, 3, 4, 1 }, { 0, 2, 5, 3 }, {
//                           1, 4, 5, 2 } } } ;
//
//               static RINGMesh2GMSH pyramid_descriptor_gmsh = { 18,                 // type
//                   5,                  // nb vertices
//                   { 0, 1, 2, 3, 4 },  // vertices
//                   5,                  // nb facets
//                   { 3, 3, 3, 3, 4 },  // nb vertices in facet
//                   { 1, 3, 4, 2, 0 },  // facets
//                   { { 0, 1, 2, 3 }, { 0, 4, 1 }, { 0, 3, 4 }, { 2, 4, 3 }, { 2, 1, 4 } } } ;
    class MSHIOHandler: public GeoModelIOHandler {
    public:
        virtual bool load( const std::string& filename, GeoModel& geomodel )
        {
            throw RINGMeshException( "I/O",
                "Loading of a GeoModel from GMSH not implemented yet" ) ;
            return false ;
        }
        virtual void save( const GeoModel& geomodel, const std::string& filename )
        {
            /// @todo after implementing GMMOrder
            throw RINGMeshException( "I/O",
                "Saving of a GeoModel from GMSH not implemented yet" ) ;
//                geomodel.set_duplicate_mode( FAULT ) ;

            std::ofstream out( filename.c_str() ) ;
            out.precision( 16 ) ;

            out << "$MeshFormat" << std::endl ;
            out << "2.2 0 8" << std::endl ;
            out << "$EndMeshFormat" << std::endl ;

            out << "$Nodes" << std::endl ;
//            out << gm.order.nb_total_vertices() << std::endl ;
//            for( index_t p = 0; p < gm.vertices.nb(); p++ ) {
//
//                const vec3& point = gm.vertices.vertex( p ) ;
//                if( p == 0 ) {
//                    std::cout << "io val " << point.x << std::endl ;
//
//                }
//                out << p + 1 << SPACE << point.x << SPACE << point.y << SPACE
//                    << point.z << std::endl ;
//            }
//            index_t vertex_offset = gm.vertices.nb() ;
//            for( index_t p = 0; p < gm.vertices.nb_duplicated_vertices(); p++ ) {
//                const vec3& point = gm.vertices.duplicated_vertex( p ) ;
//                out << vertex_offset + p + 1 << SPACE << point.x << SPACE << point.y
//                    << SPACE << point.z << std::endl ;
//            }
//            vertex_offset += gm.vertices.nb_duplicated_vertices() ;
//            index_t nb_order_vertices = gm.order.nb() ;
//            for( index_t p = 0; p < nb_order_vertices; p++ ) {
//                out << vertex_offset + p + 1 << SPACE << gm.order.point( p ).x
//                    << SPACE << gm.order.point( p ).y << SPACE
//                    << gm.order.point( p ).z << std::endl ;
//            }
//            out << "$EndNodes" << std::endl ;
//
//            index_t cell_type[4] = { 4, 5, 6, 7 } ;
//            index_t facet_type[5] = { -1, -1, -1, 2, 3 } ;
//            if( gm.get_order() == 2 ) {
//                cell_type[0] = 11 ;
//                cell_type[1] = 17 ;
//                cell_type[2] = 18 ;
//                cell_type[3] = 19 ;
//                facet_type[0] = -1 ;
//                facet_type[1] = -1 ;
//                facet_type[2] = -1 ;
//                facet_type[3] = 9 ;
//                facet_type[4] = 16 ;
//            } else if( gm.get_order() > 2 ) {
//                Logger::err( "" ) << "The order " << gm.get_order() << " "
//                    << "is not supported"
//                    << " for the gmsh export. The export will take order 1 entities"
//                    << std::endl ;
//            }
//            const GeoModel& geomodel = gm ;
//            index_t offset_region = gm.nb_regions() ;
//            index_t offset_interface = geomodel.nb_interfaces() * 2 ; // one for each side
//            index_t nb_facets = 0 ;
//            std::vector< ColocaterANN* > anns( geomodel.nb_surfaces(), nil ) ;
//            for( index_t s = 0; s < geomodel.nb_surfaces(); s++ ) {
//                if( gm.vertices.is_surface_to_duplicate( s ) )
//                    nb_facets += 2 * gm.facets.nb_facets( s ) ;
//                else
//                    nb_facets += gm.facets.nb_facets( s ) ;
//                anns[s] = new ColocaterANN( geomodel.surface( s ).mesh(),
//                    ColocaterANN::FACETS ) ;
//            }
//            out << "$Entities" << std::endl ;
//            out << gm.cells.nb_cells() + nb_facets << std::endl ;
//            index_t cur_cell = 1 ;
//            for( index_t m = 0; m < gm.nb_regions(); m++ ) {
//                const GEO::Mesh& mesh = gm.mesh( m ) ;
//                GEO::Attribute< index_t > attribute( mesh.facets.attributes(),
//                    surface_att_name ) ;
//                const GeoModelEntity& region = geomodel.region( m ) ;
//                std::vector< index_t > surfaces ;
//                surfaces.reserve( region.nb_boundaries() ) ;
//                for( index_t b = 0; b < region.nb_boundaries(); b++ ) {
//                    index_t cur_s_id = region.boundary_gme( b ).index ;
//                    if( !gm.vertices.is_surface_to_duplicate( cur_s_id ) ) continue ;
//                    surfaces.push_back( cur_s_id ) ;
//                }
//                for( index_t c = 0; c < mesh.cells.nb(); c++ ) {
//                    out << cur_cell++ << " " << cell_type[mesh.cells.type( c )]
//                        << " 2 " << m + 1 << SPACE << m ;
//                    for( index_t v = mesh.cells.corners_begin( c );
//                        v < mesh.cells.corners_end( c ); v++ ) {
//                        index_t vertex_id ;
//                        index_t duplicated_vertex_id ;
//                        out << SPACE ;
//                        if( gm.vertices.vertex_id( m, v, vertex_id,
//                            duplicated_vertex_id ) ) {
//                            out << vertex_id + 1 ;
//                        } else {
//                            out << vertex_offset + duplicated_vertex_id + 1 ;
//                        }
//                    }
//                    if( gm.get_order() == 2 ) {
//                        out << SPACE ;
//                        out << gm.order.get_id_on_cell( m, c, 3 ) + 1 ;
//                        out << SPACE ;
//                        out << gm.order.get_id_on_cell( m, c, 0 ) + 1 ;
//                        out << SPACE ;
//                        out << gm.order.get_id_on_cell( m, c, 4 ) + 1 ;
//                        out << SPACE ;
//                        out << gm.order.get_id_on_cell( m, c, 5 ) + 1 ;
//                        out << SPACE ;
//                        out << gm.order.get_id_on_cell( m, c, 1 ) + 1 ;
//                        out << SPACE ;
//                        out << gm.order.get_id_on_cell( m, c, 2 ) + 1 ;
//                    }
//                    out << std::endl ;
//
//                    for( index_t f = 0; f < mesh.cells.nb_facets( c ); f++ ) {
//                        vec3 facet_bary = mesh_cell_facet_center( mesh, c,
//                            f ) ;
//                        vec3 cell_facet_normal = mesh_cell_facet_normal( mesh,
//                            c, f ) ;
//                        for( index_t s = 0; s < surfaces.size(); s++ ) {
//                            index_t surface_id = surfaces[s] ;
//                            std::vector< index_t > result ;
//                            if( anns[surface_id]->get_colocated( facet_bary,
//                                result ) ) {
//                                vec3 facet_normal =
//                                    geomodel.surface( surface_id ).facet_normal(
//                                        result[0] ) ;
//                                bool side = dot( facet_normal, cell_facet_normal )
//                                    > 0 ;
//                                out << cur_cell++ << " "
//                                    << facet_type[mesh.cells.facet_nb_vertices( c,
//                                        f )] << " 2 "
//                                    << offset_region
//                                        + 2
//                                            * geomodel.surface( surface_id ).parent_id().index
//                                        + side + 1 << SPACE
//                                    << offset_region + offset_interface
//                                        + 2 * surface_id + side ;
//                                for( index_t v = 0;
//                                    v < mesh.cells.facet_nb_vertices( c, f ); v++ ) {
//                                    index_t corner_id =
//                                        mesh.cells.corner( c,
//                                            mesh.cells.descriptor( c ).facet_vertex[f][v] ) ;
//                                    index_t vertex_id ;
//                                    index_t duplicated_vertex_id ;
//                                    out << SPACE ;
//                                    if( gm.vertices.vertex_id( m, corner_id,
//                                        vertex_id, duplicated_vertex_id ) ) {
//                                        out << vertex_id + 1 ;
//                                    } else {
//                                        out
//                                            << vertex_offset + duplicated_vertex_id
//                                                + 1 ;
//                                    }
//                                }
//                                out << std::endl ;
//                                break ;
//                            }
//                        }
//                    }
//                }
//            }
//
//            for( index_t i = 0; i < geomodel.nb_interfaces(); i++ ) {
//                const GeoModelEntity& interf = geomodel.one_interface( i ) ;
//                for( index_t s = 0; s < interf.nb_children(); s++ ) {
//                    index_t s_id = interf.child_id( s ).index ;
//                    if( gm.vertices.is_surface_to_duplicate( s_id ) ) continue ;
//                    index_t mesh_id = gm.facets.mesh( s_id ) ;
//                    const GEO::Mesh& mesh = gm.mesh( mesh_id ) ;
//                    for( index_t t = 0; t < gm.facets.nb_facets( s_id ); t++ ) {
//                        index_t facet_id = gm.facets.facet( s_id, t ) ;
//                        out << cur_cell++ << SPACE
//                            << facet_type[mesh.facets.nb_vertices( facet_id )]
//                            << " 2 " << offset_region + 2 * i + 1 << SPACE
//                            << offset_region + offset_interface + 2 * s_id ;
//                        for( index_t v = 0; v < mesh.facets.nb_vertices( facet_id );
//                            v++ ) {
//                            index_t v_id = mesh.facets.vertex( facet_id, v ) ;
//                            out << SPACE
//                                << gm.vertices.vertex_id( mesh_id, v_id ) + 1 ;
//                        }
//                        for( index_t v = 0;
//                            v
//                                < mesh.facets.nb_vertices( facet_id )
//                                    * ( gm.get_order() - 1 ); v++ ) {
//                            out << SPACE ;
//                            out << gm.order.get_id_on_facet( s, facet_id, v ) + 1 ;
//                        }
//                        out << std::endl ;
//                    }
//                }
//            }
//            out << "$EndEntities" << std::endl ;
//
//            if( GEO::CmdLine::get_arg_bool( "out:kine3d" ) ) {
//                std::string directory = GEO::FileSystem::dir_name( filename ) ;
//                std::string file = GEO::FileSystem::base_name( filename ) ;
//                std::ostringstream oss_kine ;
//                oss_kine << directory << "/" << file << ".gmsh_info" ;
//                std::ofstream kine3d( oss_kine.str().c_str() ) ;
//                for( index_t i = 0; i < geomodel.nb_interfaces(); i++ ) {
//                    const GeoModelEntity& interf = geomodel.one_interface( i ) ;
//                    index_t s_id = interf.child_id( 0 ).index ;
//                    kine3d << offset_region + 2 * i + 1 << ":" << interf.name()
//                        << ",1," ;
//                    const RINGMesh::GeoModelEntity& E = geomodel.one_interface( i ) ;
//                    if( RINGMesh::GeoModelEntity::is_fault(
//                        E.geological_feature() ) ) {
//                        kine3d << "FaultFeatureClass" ;
//                    } else if( RINGMesh::GeoModelEntity::is_stratigraphic_limit(
//                        E.geological_feature() ) ) {
//                        kine3d << "HorizonFeatureClass" ;
//                    } else if( E.is_on_voi() ) {
//                        kine3d << "ModelRINGMesh::BoundaryFeatureClass" ;
//                    }
//                    kine3d << std::endl ;
//                    if( gm.vertices.is_surface_to_duplicate( s_id ) ) {
//                        kine3d << offset_region + 2 * i + 1 << ":" << interf.name()
//                            << ",0," ;
//                        const RINGMesh::GeoModelEntity& E = geomodel.one_interface(
//                            i ) ;
//                        if( RINGMesh::GeoModelEntity::is_fault(
//                            E.geological_feature() ) ) {
//                            kine3d << "FaultFeatureClass" ;
//                        } else if( RINGMesh::GeoModelEntity::is_stratigraphic_limit(
//                            E.geological_feature() ) ) {
//                            kine3d << "HorizonFeatureClass" ;
//                        } else if( E.is_on_voi() ) {
//                            kine3d << "ModelRINGMesh::BoundaryFeatureClass" ;
//                        }
//                        kine3d << std::endl ;
//                    }
//                }
//            }
        }
    } ;

}
