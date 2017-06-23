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

    class GeoModelBuilderSVG final : public GeoModelBuilderFile< 2 > {
    public:
        GeoModelBuilderSVG( GeoModel< 2 >& geomodel, const std::string& filename )
            : GeoModelBuilderFile< 2 >( geomodel, filename )
        {
        }
        virtual ~GeoModelBuilderSVG() = default;

        void load_file()
        {
            tinyxml2::XMLDocument svg_file;
            if( svg_file.LoadFile( filename_.c_str() ) != tinyxml2::XML_SUCCESS ) {
                throw RINGMeshException( "I/O", "Error while loading svg file." );
            }
            tinyxml2::XMLNode* svg_node = svg_file.FirstChild();
            if( svg_node == nullptr ) {
                throw RINGMeshException( "I/O",
                    "Error while getting root of svg file." );
            }

            tinyxml2::XMLElement* group = svg_node->FirstChildElement( "g" );
            while( group != nullptr ) {
                DEBUG( group->GetText() );
                group = group->NextSiblingElement();
            }
        }
    };

    class SVGIOHandler final: public GeoModelIOHandler< 2 > {
    public:
        virtual ~SVGIOHandler() = default;

        virtual bool load( const std::string& filename, GeoModel< 2 >& geomodel ) final
        {
            GeoModelBuilderSVG builder( geomodel, filename );
            builder.build_geomodel();
            return false;
        }
        virtual void save(
            const GeoModel< 2 >& geomodel,
            const std::string& filename ) final
        {
            throw RINGMeshException( "I/O",
                "Saving a GeoModel in svg not implemented yet" );
        }
    };

}
