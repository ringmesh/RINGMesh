/*
 * Copyright (c) 2012-2016, Association Scientifique pour la Geologie et ses Applications (ASGA)
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

#include <ringmesh/io/io.h>

#include <iomanip>
#include <fstream>

#include <geogram/basic/file_system.h>

#include <ringmesh/geomodel/geo_model.h>
#include <ringmesh/geomodel/geo_model_geological_entity.h>


namespace {
    using namespace RINGMesh ;


    /************************************************************************/
    class TOTALIOHandler: public GeoModelIOHandler {
    public:
        virtual void load( const std::string& filename, GeoModel& geomodel )
        {
            throw RINGMeshException( "I/O",
                "Loading of a GeoModel from TOTAL mesh not implemented yet" ) ;
        }
        virtual void save( const GeoModel& gm, const std::string& filename )
        {
        	std::ofstream out( filename.c_str() ) ;
        	out.precision( 16 ) ;

        	write_vertices( out, gm ) ;
        	write_facets( out, gm ) ;
        	write_cells( out, gm ) ;
        }

    private:
        void write_vertices( std::ofstream& out, const GeoModel& gm )
        {
        	out << "## VERTICES ##" << std::endl ;
        	const GeoModelMeshVertices& vertices = gm.mesh.vertices ;
        	for( index_t v = 0; v < vertices.nb(); v++ ) {
        		out << vertices.vertex( v ) << std::endl ;
        	}
    		out << std::endl ;
        }


        void write_cells( std::ofstream& out, const GeoModel& gm )
        {
        	out << "## CELLS ##" << std::endl ;
        	const GeoModelMeshCells& cells = gm.mesh.cells ;
        	for( index_t c = 0; c < cells.nb(); c++ ) {
        		for( index_t v = 0; v < cells.nb_vertices( c ); v++ ) {
        			out << cells.vertex( c, v ) ;
        		}
        		out << std::endl ;
        	}
    		out << std::endl ;
        }

        void write_facets( std::ofstream& out, const GeoModel& gm )
        {
        	out << "## FACETS ##" << std::endl ;
        	const GeoModelMeshFacets& facets = gm.mesh.facets ;
        	for( index_t i = 0; i < gm.nb_geological_entities( Interface::type_name_static() ); i++ ) {
        		const GeoModelGeologicalEntity& interf = gm.geological_entity( Interface::type_name_static(), i ) ;
        		out << "# " << interf.name() << " #" << std::endl ;
        		for( index_t s = 0; s < interf.nb_children(); s++ ) {
        			index_t surface = interf.child_gme( s ).index ;
        			write_facet_type( out, facets, surface, GeoModelMeshFacets::TRIANGLE ) ;
        			write_facet_type( out, facets, surface, GeoModelMeshFacets::QUAD ) ;
        		}
        		out << std::endl ;
        	}
    		out << std::endl ;
        }

        void write_facet_type( std::ofstream& out, const GeoModelMeshFacets& facets, index_t surface, GeoModelMeshFacets::FacetType type )
        {
        	for( index_t f = 0; f < facets.nb_facets( surface, type ); f++ ) {
        		index_t facet = facets.facet( surface, f, type ) ;
        		for( index_t v = 0; v < facets.nb_vertices( facet ); v++ ) {
        			out << facets.vertex( facet, v ) ;
        		}
        		out << std::endl ;
        	}
        }
    } ;

}
