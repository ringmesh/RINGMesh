/*
 * Copyright (c) 2012-2018, Association Scientifique pour la Geologie et ses
 * Applications (ASGA). All rights reserved.
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
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
 * THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 * PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL ASGA BE LIABLE FOR ANY DIRECT,
 * INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF
 * THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 *     http://www.ring-team.org
 *
 *     RING Project
 *     Ecole Nationale Superieure de Geologie - GeoRessources
 *     2 Rue du Doyen Marcel Roubault - TSA 70605
 *     54518 VANDOEUVRE-LES-NANCY
 *     FRANCE
 */

namespace
{
    /*!
     * @brief Total number of polygons in the Surfaces of a geomodel
     */
    index_t nb_polygons( const GeoModel3D& geomodel )
    {
        index_t result{ 0 };
        for( const auto& surface : geomodel.surfaces() )
        {
            result += surface.nb_mesh_elements();
        }
        return result;
    }

    /*!
     * @brief Write a region information in a stream
     * @details Used by function to save the Model in a .ml file
     *
     * @param[in] count Region index in the file
     * @param[in] region The region to save
     * @param[in,out] out The file output stream
     */
    void save_region( index_t count, const Region3D& region, std::ostream& out )
    {
        out << "REGION " << count << "  " << region.name() << " " << EOL;
        index_t it{ 0 };

        for( auto i : range( region.nb_boundaries() ) )
        {
            out << "  ";
            if( region.side( i ) )
            {
                out << "+";
            }
            else
            {
                out << "-";
            }
            out << region.boundary( i ).index() + 1;
            it++;
            if( it == 5 )
            {
                out << EOL;
                it = 0;
            }
        }
        out << "  0" << EOL;
    }

    void save_universe(
        index_t count, const GeoModel3D& geomodel, std::ostream& out )
    {
        const SurfaceSide surface_region_sides = geomodel.voi_surfaces();

        out << "REGION " << count << "  Universe " << EOL;
        index_t it{ 0 };

        for( auto i : range( surface_region_sides.surfaces_.size() ) )
        {
            out << "  ";
            if( surface_region_sides.sides_[i] )
            {
                out << "+";
            }
            else
            {
                out << "-";
            }
            out << surface_region_sides.surfaces_[i] + 1;
            it++;
            if( it == 5 )
            {
                out << EOL;
                it = 0;
            }
        }
        out << "  0" << EOL;
    }

    /*!
     * @brief Write information for on layer in a stream
     * @details Used by function to save the Model in a .ml file
     *
     * @param[in] count Index of the layer in the file
     * @param[in] offset Offset of region indices in the file
     * @param[in] layer The layer to write
     * @param[in,out] out The output file stream
     */
    void save_layer( index_t offset,
        const GeoModelGeologicalEntity3D& layer,
        std::ostream& out )
    {
        out << "LAYER " << layer.name() << " " << EOL;
        index_t it{ 0 };

        for( auto i : range( layer.nb_children() ) )
        {
            out << "  " << layer.child_gmme( i ).index() + offset + 1;
            it++;
            if( it == 5 )
            {
                out << EOL;
                it = 0;
            }
        }
        out << "  0" << EOL;
    }

    /*!
     * @brief Write basic header for Gocad coordinate system.
     * @param[in,out] out Output .ml file stream
     */
    void save_coordinate_system( std::ostream& out )
    {
        out << "GOCAD_ORIGINAL_COORDINATE_SYSTEM" << EOL << "NAME Default"
            << EOL << "AXIS_NAME \"X\" \"Y\" \"Z\"" << EOL
            << "AXIS_UNIT \"m\" \"m\" \"m\"" << EOL << "ZPOSITIVE Elevation"
            << EOL << "END_ORIGINAL_COORDINATE_SYSTEM" << EOL;
    }

    /*!
     * @brief Check if the geomodel can be saved in a skua-gocad .ml file
     * @details It assumes that the geomodel is valid and verifies that:
     *   - all Interfaces have a name and geological feature
     *   - all Surfaces are in an Interface
     *   - all Surfaces are triangulated
     *   - all Regions have a name
     */
    bool check_gocad_validity( const GeoModel3D& geomodel )
    {
        auto nb_interfaces =
            geomodel.nb_geological_entities( Interface3D::type_name_static() );
        if( nb_interfaces == 0 )
        {
            Logger::err(
                "", " The GeoModel ", geomodel.name(), " has no Interface" );
            return false;
        }
        for( auto& cur_interface :
            geomodel.geol_entities( Interface3D::type_name_static() ) )
        {
            if( !cur_interface.has_geological_feature() )
            {
                Logger::err(
                    "", cur_interface.gmge(), " has no geological feature" );
                return false;
            }
        }
        for( const auto& surface : geomodel.surfaces() )
        {
            if( !surface.has_parent() )
            {
                Logger::err( "", surface.gmme(),
                    " does not belong to any Interface of the geomodel" );
                return false;
            }
            if( !surface.is_simplicial() )
            {
                Logger::err( "", surface.gmme(), " is not triangulated " );
                return false;
            }
        }
        return true;
    }

    /*! Brute force inefficient but I am debugging !!!! */
    bool has_surface_edge(
        const Surface3D& surface, index_t v0_in, index_t v1_in )
    {
        for( auto i : range( surface.nb_mesh_elements() ) )
        {
            for( auto j : range( surface.nb_mesh_element_vertices( i ) ) )
            {
                auto v0 = surface.mesh_element_vertex_index( { i, j } );
                auto v1 = surface.mesh_element_vertex_index(
                    surface.mesh().next_polygon_vertex( { i, j } ) );
                if( ( v0 == v0_in && v1 == v1_in )
                    || ( v0 == v1_in && v1 == v0_in ) )
                {
                    return true;
                }
            }
        }
        return false;
    }

    /*!
     * @brief Save the geomodel in a .ml file if it can
     * @param[in] geomodel the geomodel to save
     * @param[in,out] out Output file stream
     */
    void save_gocad_model3d( const GeoModel3D& geomodel, std::ostream& out )
    {
        if( !check_gocad_validity( geomodel ) )
        {
            throw RINGMeshException( "I/O", " The GeoModel ", geomodel.name(),
                " cannot be saved in .ml format" );
        }
        out.precision( 16 );

        // Gocad Model3d headers
        out << "GOCAD Model3d 1" << EOL << "HEADER {" << EOL
            << "name: " << geomodel.name() << EOL << "}" << EOL;

        save_coordinate_system( out );

        // Gocad::TSurf = RINGMesh::Interface
        for( auto& tsurf :
            geomodel.geol_entities( Interface3D::type_name_static() ) )
        {
            out << "TSURF " << tsurf.name() << EOL;
        }

        index_t count{ 1 };

        // Gocad::TFace = RINGMesh::Surface
        for( const auto& surface : geomodel.surfaces() )
        {
            const gmge_id& parent_interface =
                surface.parent_gmge( Interface3D::type_name_static() );
            if( !parent_interface.is_defined() )
            {
                throw RINGMeshException( "I/O", "Failed to save GeoModel",
                    " in .ml Gocad format because Surface ", surface.index(),
                    " has no Interface parent)" );
            }
            const auto& cur_geol_feature =
                geomodel.geological_entity( parent_interface )
                    .geological_feature();

            out << "TFACE " << count << "  ";
            out << GeoModelGeologicalEntity3D::geol_name( cur_geol_feature );
            out << " "
                << surface.parent( Interface3D::type_name_static() ).name()
                << EOL;

            // Print the key polygon which is the first three
            // vertices of the first polygon
            out << "  " << surface.mesh_element_vertex( { 0, 0 } ) << EOL;
            out << "  " << surface.mesh_element_vertex( { 0, 1 } ) << EOL;
            out << "  " << surface.mesh_element_vertex( { 0, 2 } ) << EOL;

            ++count;
        }
        // Universe
        auto offset_layer = count;
        save_universe( count, geomodel, out );
        ++count;
        // Regions
        for( const auto& region : geomodel.regions() )
        {
            save_region( count, region, out );
            ++count;
        }
        // Layers
        if( geomodel.entity_type_manager()
                .geological_entity_manager.is_valid_type(
                    Layer3D::type_name_static() ) )
        {
            for( auto& layer :
                geomodel.geol_entities( Layer3D::type_name_static() ) )
            {
                save_layer( offset_layer, layer, out );
                ++count;
            }
        }
        out << "END" << EOL;

        const auto& geomodel_vertices = geomodel.mesh.vertices;
        // Save the geometry of the Surfaces, Interface per Interface
        for( auto& tsurf :
            geomodel.geol_entities( Interface3D::type_name_static() ) )
        {
            // TSurf beginning header
            out << "GOCAD TSurf 1" << EOL << "HEADER {" << EOL
                << "name:" << tsurf.name() << EOL
                << "name_in_model_list:" << tsurf.name() << EOL << "}" << EOL;
            save_coordinate_system( out );

            out << "GEOLOGICAL_FEATURE " << tsurf.name() << EOL
                << "GEOLOGICAL_TYPE ";
            out << GeoModelGeologicalEntity< 3 >::geol_name(
                tsurf.geological_feature() );
            out << EOL;
            out << "PROPERTY_CLASS_HEADER Z {" << EOL << "is_z:on" << EOL << "}"
                << EOL;

            index_t vertex_count{ 1 };
            // TFace vertex index = Surface vertex index + offset
            auto offset = vertex_count;

            // To collect Corners(BStones) indexes
            // and boundary (Line) first and second vertex indexes
            std::set< index_t > corners;
            std::set< std::pair< index_t, index_t > > lineindices;
            for( auto j : range( tsurf.nb_children() ) )
            {
                offset = vertex_count;
                const auto& surface =
                    dynamic_cast< const Surface3D& >( tsurf.child( j ) );

                out << "TFACE" << EOL;
                for( auto k : range( surface.nb_vertices() ) )
                {
                    out << "VRTX " << vertex_count << " " << surface.vertex( k )
                        << EOL;
                    vertex_count++;
                }
                for( auto k : range( surface.nb_mesh_elements() ) )
                {
                    out << "TRGL "
                        << surface.mesh_element_vertex_index( { k, 0 } )
                               + offset
                        << " "
                        << surface.mesh_element_vertex_index( { k, 1 } )
                               + offset
                        << " "
                        << surface.mesh_element_vertex_index( { k, 2 } )
                               + offset
                        << EOL;
                }
                for( auto k : range( surface.nb_boundaries() ) )
                {
                    const auto& line = surface.boundary( k );
                    auto v0_model_id =
                        geomodel_vertices.geomodel_vertex_id( line.gmme(), 0 );
                    auto v1_model_id =
                        geomodel_vertices.geomodel_vertex_id( line.gmme(), 1 );

                    auto v0_surface_ids =
                        geomodel_vertices.mesh_entity_vertex_id(
                            surface.gmme(), v0_model_id );
                    auto v1_surface_ids =
                        geomodel_vertices.mesh_entity_vertex_id(
                            surface.gmme(), v1_model_id );

                    if( !surface.has_inside_border() )
                    {
                        auto v0 = v0_surface_ids[0];
                        auto v1 = v1_surface_ids[0];
                        v0 += offset;
                        v1 += offset;

                        lineindices.insert(
                            std::pair< index_t, index_t >( v0, v1 ) );
                        corners.insert( v0 );
                    }
                    else
                    {
                        // We need to get the right pair of v0 - v1  (not
                        // crossing the inside boundary)
                        // corner and a border
                        int count = 0;
                        bool to_break = false;
                        for( auto v0 : v0_surface_ids )
                        {
                            for( auto v1 : v1_surface_ids )
                            {
                                if( has_surface_edge( surface, v0, v1 ) )
                                {
                                    lineindices.insert(
                                        std::pair< index_t, index_t >(
                                            v0 + offset, v1 + offset ) );
                                    count++;
                                }
                                if( !line.is_inside_border( surface )
                                    && count == 1 )
                                {
                                    to_break = true;
                                    break;
                                }
                                else if( count == 2 )
                                {
                                    to_break = true;
                                    break;
                                }
                            }
                            if( to_break )
                            {
                                corners.insert( v0 + offset );
                                break;
                            }
                        }
                    }

                    // Set a BSTONE at the line other extremity
                    const auto& c1_id = line.boundary_gmme( 1 );
                    auto gme_vertices =
                        geomodel_vertices.mesh_entity_vertex_id( surface.gmme(),
                            geomodel_vertices.geomodel_vertex_id( c1_id ) );
                    corners.insert( gme_vertices.front() + offset );
                }
            }
            // Add the remaining bstones that are not already in bstones
            for( auto it( corners.begin() ); it != corners.end(); ++it )
            {
                out << "BSTONE " << *it << EOL;
            }
            for( auto it( lineindices.begin() ); it != lineindices.end(); ++it )
            {
                out << "BORDER " << vertex_count << " " << it->first << " "
                    << it->second << EOL;
                vertex_count++;
            }
            out << "END" << EOL;
        }
    }

    class MLIOHandler final : public GeoModelInputHandler3D,
                              public GeoModelOutputHandler3D
    {
    public:
        /*! Load a .ml (Gocad file)
         * @pre Filename is valid
         */
        void load( const std::string& filename, GeoModel3D& geomodel ) final
        {
            std::ifstream input( filename.c_str() );
            if( !input )
            {
                throw RINGMeshException(
                    "I/O", "Failed to open file ", filename );
            }
            GeoModelBuilderML builder( geomodel, filename );
            builder.build_geomodel();
        }

        void save(
            const GeoModel3D& geomodel, const std::string& filename ) final
        {
            std::ofstream out( filename.c_str() );
            save_gocad_model3d( geomodel, out );
            out << std::flush;
        }
    };
}
