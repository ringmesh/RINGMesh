/*
 * Copyright (c) 2012-2017, Association Scientifique pour la Geologie et ses
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

/*!
 * @brief Classes to load and build a GeoModel2D from a .svg
 * @author Arnaud Botella
 */

namespace
{
    class GeoModelBuilderSVG final : public GeoModelBuilderFile< 2 >
    {
    public:
        GeoModelBuilderSVG( GeoModel2D& geomodel, const std::string& filename )
            : GeoModelBuilderFile< 2 >( geomodel, filename ), height_( 0. )
        {
        }
        virtual ~GeoModelBuilderSVG() = default;

        void load_file()
        {
            tinyxml2::XMLDocument svg_file;
            if( svg_file.LoadFile( filename().c_str() )
                != tinyxml2::XML_SUCCESS )
            {
                throw RINGMeshException(
                    "I/O", "Error while loading svg file." );
            }
            tinyxml2::XMLElement* svg_node =
                svg_file.FirstChildElement( "svg" );
            if( svg_node == nullptr )
            {
                throw RINGMeshException(
                    "I/O", "Error while getting root of svg file." );
            }
            height_ = svg_node->DoubleAttribute( "height" );
            create_lines( svg_node );
            create_corners();
            build_surfaces_from_corners_and_lines();
        }

    private:
        void create_lines( tinyxml2::XMLElement* svg_node )
        {
            tinyxml2::XMLElement* group = svg_node->FirstChildElement( "g" );
            while( group != nullptr )
            {
                vec2 group_translation = get_transform( group );
                tinyxml2::XMLElement* path = group->FirstChildElement( "path" );
                while( path != nullptr )
                {
                    vec2 translation =
                        group_translation + get_transform( path );
                    std::string data = path->Attribute( "d" );
                    std::vector< vec2 > vertices =
                        get_path_vertices( data, translation );
                    index_t line_id =
                        topology
                            .create_mesh_entity( Line2D::type_name_static() )
                            .index();
                    geometry.set_line( line_id, vertices );
                    path = path->NextSiblingElement( "path" );
                }
                group = group->NextSiblingElement( "g" );
            }
        }

        vec2 get_transform( tinyxml2::XMLElement* group )
        {
            const char* attribute = group->Attribute( "transform" );
            if( attribute == nullptr )
            {
                return vec2();
            }
            std::string action, parameters;
            GEO::String::split_string( attribute, '(', action, parameters );
            if( action != "translate" )
            {
                throw RINGMeshException( "I/O", "Forbidden transformation ",
                    action, " found in the svg file. " );
            }
            parameters.pop_back();
            std::vector< std::string > coordinates;
            GEO::String::split_string( parameters, ',', coordinates );
            vec2 translation;
            translation.x = GEO::String::to_double( coordinates.front() );
            translation.y = GEO::String::to_double( coordinates.back() );
            return translation;
        }

        void create_corners()
        {
            std::vector< vec2 > point_extremities;
            point_extremities.reserve( geomodel_.nb_lines() * 2 );
            for( const auto& line : line_range< 2 >( geomodel_ ) )
            {
                point_extremities.push_back( line.vertex( 0 ) );
                point_extremities.push_back(
                    line.vertex( line.nb_vertices() - 1 ) );
            }

            NNSearch2D nn_search( point_extremities );
            std::vector< index_t > index_map;
            std::vector< vec2 > unique_points;
            std::tie( std::ignore, index_map, unique_points ) =
                nn_search.get_colocated_index_mapping_and_unique_points(
                    geomodel_.epsilon() );

            topology.create_mesh_entities( Corner2D::type_name_static(),
                static_cast< index_t >( unique_points.size() ) );
            for( index_t c : range( geomodel_.nb_corners() ) )
            {
                geometry.set_corner( c, unique_points[c] );
            }
            index_t index = 0;
            for( const auto& line : line_range< 2 >( geomodel_ ) )
            {
                const auto line_id = line.gmme();
                const index_t corner0 = index_map[index++];
                const index_t corner1 = index_map[index++];
                topology.add_line_corner_boundary_relation(
                    line_id.index(), corner0 );
                topology.add_line_corner_boundary_relation(
                    line_id.index(), corner1 );

                // Update line vertex extremities with corner coordinates
                geometry.set_mesh_entity_vertex(
                    line_id, 0, unique_points[corner0], false );
                geometry.set_mesh_entity_vertex( line_id,
                    line.nb_vertices() - 1, unique_points[corner1], false );
            }
        }

        std::vector< vec2 > get_path_vertices(
            std::string& data, const vec2& translation ) const
        {
            std::vector< std::string > tokens = get_tokens( data );
            std::vector< vec2 > vertices;
            vertices.reserve( tokens.size() );
            bool is_absolute = true;
            for( index_t i = 0; i < tokens.size(); i++ )
            {
                std::string& token = tokens[i];
                if( std::isalpha( token.front() ) )
                {
                    ringmesh_assert( token.size() == 1 );
                    is_absolute = is_command_absolute( token );
                }
                else
                {
                    vec2 vertex;
                    GEO::String::from_string( token, vertex.x );
                    GEO::String::from_string( tokens[++i], vertex.y );
                    if( is_absolute || vertices.empty() )
                    {
                        vertex += translation;
                        vertex.y = height_ - vertex.y;
                        vertices.push_back( vertex );
                    }
                    else
                    {
                        vertex.y = -vertex.y;
                        vertices.push_back( vertices.back() + vertex );
                    }
                }
            }
            return vertices;
        }

        void throw_if_forbidden_command( const std::string& data ) const
        {
            std::string forbidden_letters( "HhVvCcSs" );
            for( char letter : forbidden_letters )
            {
                if( data.find( letter ) != std::string::npos )
                {
                    throw RINGMeshException( "I/O", "Forbidden command ",
                        letter, " found in the svg file. ",
                        "Flatten your paths before importing in RINGMesh." );
                }
            }
        }

        std::vector< std::string > get_tokens( std::string& data ) const
        {
            throw_if_forbidden_command( data );
            std::replace( data.begin(), data.end(), ',', ' ' );
            std::string letters( "MmLl" );
            for( char letter : letters )
            {
                std::string string_letter = GEO::String::to_string( letter );
                replace_all_string_occurens(
                    data, string_letter, " " + string_letter + " " );
            }
            std::vector< std::string > tokens;
            GEO::String::split_string( data, ' ', tokens );
            return tokens;
        }

        void replace_all_string_occurens( std::string& context,
            const std::string& to_replace,
            const std::string& replacement ) const
        {
            std::size_t look_here = 0;
            std::size_t found_here;
            while( ( found_here = context.find( to_replace, look_here ) )
                   != std::string::npos )
            {
                context.replace( found_here, to_replace.size(), replacement );
                look_here = found_here + replacement.size();
            }
        }

        bool is_command_absolute( const std::string& letter ) const
        {
            return letter == "M" || letter == "L";
        }

    private:
        double height_;
    };

    class SVGIOHandler final : public GeoModelInputHandler2D
    {
    public:
        void load( const std::string& filename, GeoModel2D& geomodel ) final
        {
            GeoModelBuilderSVG builder( geomodel, filename );
            builder.build_geomodel();
        }
    };
}
