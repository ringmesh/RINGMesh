/*
 *  Copyright (c) 2012-2014, Bruno Levy
 *  All rights reserved.
 *
 *  Redistribution and use in source and binary forms, with or without
 *  modification, are permitted provided that the following conditions are met:
 *
 *  * Redistributions of source code must retain the above copyright notice,
 *  this list of conditions and the following disclaimer.
 *  * Redistributions in binary form must reproduce the above copyright notice,
 *  this list of conditions and the following disclaimer in the documentation
 *  and/or other materials provided with the distribution.
 *  * Neither the name of the ALICE Project-Team nor the names of its
 *  contributors may be used to endorse or promote products derived from this
 *  software without specific prior written permission.
 * 
 *  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 *  AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 *  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 *  ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 *  LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 *  CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 *  SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 *  INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 *  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 *  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 *  POSSIBILITY OF SUCH DAMAGE.
 *
 *  If you modify this software, you should include a notice giving the
 *  name of the person performing the modification, the date of modification,
 *  and the reason for such modification.
 *
 *  Contact: Bruno Levy
 *
 *     Bruno.Levy@inria.fr
 *     http://www.loria.fr/~levy
 *
 *     ALICE Project
 *     LORIA, INRIA Lorraine, 
 *     Campus Scientifique, BP 239
 *     54506 VANDOEUVRE LES NANCY CEDEX 
 *     FRANCE
 *
 */

#include <geogram/mesh/mesh_io.h>
#include <geogram/mesh/mesh.h>
#include <geogram/mesh/mesh_builder.h>
#include <geogram/mesh/mesh_private.h>
#include <geogram/mesh/mesh_geometry.h>
#include <geogram/basic/line_stream.h>
#include <geogram/basic/b_stream.h>
#include <geogram/basic/file_system.h>
#include <geogram/basic/command_line.h>
#include <geogram/basic/argused.h>
#include <fstream>

extern "C" {
#include <geogram/third_party/LM6/libmesh6.h>
}

#include <geogram/third_party/rply/rply.h>

// TODO: take ioflags into account.

namespace {

    using namespace GEO;

    /**
     * \brief Creates mesh attributes.
     * \param[in,out] M the mesh
     * \param[in] attributes the list of mesh attributes to create
     * \note For now, only MESH_FACET_REGION is implemented
     */
    template <class MESH>
    void create_attributes(
        MESH& M, MeshAttributes attributes
    ) {
        if(attributes & MESH_FACET_REGION) {
            GEOGen::MeshMutator<MESH>::set_attributes(
                M, MeshAttributes(M.attributes() | MESH_FACET_REGION)
            );
            GEOGen::MeshMutator<MESH>::facet_regions(M).resize(M.nb_facets());
        }
        geo_assert(!(attributes & MESH_VERTEX_REGION));
        geo_assert(!(attributes & MESH_VERTEX_NORMAL));
        geo_assert(!(attributes & MESH_VERTEX_COLOR));
        geo_assert(!(attributes & MESH_CELL_REGION));
    }

    /************************************************************************/

    /**
     * \brief IO handler for AliasWavefront OBJ format.
     * \see http://en.wikipedia.org/wiki/Wavefront_.obj_file
     */
    class GEOGRAM_API OBJIOHandler : public MeshIOHandler {
    public:
        /**
         * \brief Creates a OBJ IO handler.
         * \param[in] dim dimension of the vertices (3 for regular 3d mesh)
         */
        OBJIOHandler(coord_index_t dimension = 3) :
            dimension_(dimension) {
        }

        /**
         * \brief Loads a mesh from a file in OBJ format.
         * \param[in] filename name of the file
         * \param[out] M the loaded mesh
         * \param[in] ioflags specifies which attributes and 
         *  elements should be read
         * \param[in] dim dimension of the vertices (3 for regular 3d meshes)
         * \tparam MESH the mesh class
         * \return true on success, false otherwise
         */
        template <class MESH>
        bool load_generic(
            const std::string& filename, MESH& M,
            const MeshIOFlags& ioflags
        ) {
            vector<double> P(dimension_);
            if(M.dimension() != dimension_) {
                M.set_dimension(dimension_);
            }
            GEOGen::MeshBuilder<MESH> builder(&M);

            LineInput in(filename);
            if(!in.OK()) {
                return false;
            }

            bool read_facet_regions = false;
            bool first_facet_attribute = true;
            builder.begin_mesh();
            while(!in.eof() && in.get_line()) {
                in.get_fields();
                if(in.nb_fields() >= 1) {
                    if(in.field_matches(0, "v")) {
                        for(coord_index_t c = 0; c < dimension_; c++) {
                            if(index_t(c + 1) < in.nb_fields()) {
                                P[c] = in.field_as_double(index_t(c + 1));
                            } else {
                                P[c] = 0.0;
                            }
                        }
                        builder.add_vertex_by_ptr(&P[0]);
                    } else if(in.field_matches(0, "f")) {
                        if(in.nb_fields() < 4) {
                            Logger::err("I/O")
                                << "Line " << in.line_number()
                                << ": facet only has " << in.nb_fields()
                                << " corners (at least 3 required)"
                                << std::endl;
                            return false;
                        }
                        builder.begin_facet();
                        for(index_t i = 1; i < in.nb_fields(); i++) {
                            for(char* ptr = in.field(i); *ptr != '\0'; ptr++) {
                                if(*ptr == '/') {
                                    *ptr = '\0';
                                    break;
                                }
                            }
                            index_t vertex_index = in.field_as_uint(i);
                            if(
                                (vertex_index < 1) ||
                                (index_t(vertex_index) > M.nb_vertices())
                            ) {
                                Logger::err("I/O")
                                    << "Line " << in.line_number()
                                    << ": facet corner #" << i
                                    << " references an invalid vertex: "
                                    << vertex_index
                                    << std::endl;
                                return false;
                            }
                            builder.add_vertex_to_facet(vertex_index - 1);
                        }
                        builder.end_facet();
                    } else if(
                        ioflags.has_attribute(MESH_FACET_REGION) &&
                        in.field_matches(0, "#")
                    ) {
                        if(
                            in.nb_fields() >= 5 &&
                            in.field_matches(1, "attribute") &&
                            in.field_matches(3, "facet")
                        ) {
                            if(
                                first_facet_attribute &&
                                in.field_matches(2, "chart") &&
                                in.field_matches(4, "integer")
                            ) {
                                read_facet_regions = true;
                            } else {
                                first_facet_attribute = false;
                            }
                        } else if(
                            read_facet_regions &&
                            in.nb_fields() >= 5 &&
                            in.field_matches(1, "attrs") &&
                            in.field_matches(2, "f")
                        ) {
                            index_t facet_index = in.field_as_uint(3);
                            signed_index_t facet_region = in.field_as_int(4);

                            if(
                                (facet_index < 1) ||
                                (facet_index > M.nb_vertices())
                            ) {
                                Logger::err("I/O")
                                    << "Line " << in.line_number()
                                    << ": facet attributes "
                                    << "reference an invalid facet: "
                                    << facet_index
                                    << std::endl;
                                return false;
                            }
                            if(!M.has_attribute(MESH_FACET_REGION)) {
                                create_attributes(M, MESH_FACET_REGION);
                            }
                            GEOGen::MeshMutator<MESH>::facet_regions(M)[
                                facet_index - 1
                            ] = facet_region;
                        }
                    }
                }
            }
            builder.end_mesh(false);
            return true;
        }

        virtual bool save(
            const Mesh& M, const std::string& filename,
            const MeshIOFlags& ioflags
        ) {
            geo_assert(M.dimension() >= dimension_);
            std::ofstream out(filename.c_str());
            if(!out) {
                Logger::err("I/O")
                    << "Could not create file \'" 
                    << filename << "\'" << std::endl;
                return false;
            }
            
            std::vector<std::string> args;
            CmdLine::get_args(args);
            for(index_t i = 0; i < args.size(); i++) {
                out << "# vorpaline " << args[i] << std::endl;
            }

            for(index_t v = 0; v < M.nb_vertices(); v++) {
                const double* p = M.vertex_ptr(v);
                out << "v ";
                for(index_t c = 0; c < dimension_; c++) {
                    out << p[c] << ' ';
                }
                out << std::endl;
            }
            for(index_t f = 0; f < M.nb_facets(); f++) {
                out << "f ";
                for(index_t i = M.facet_begin(f); i < M.facet_end(f); i++) {
                    out << M.corner_vertex_index(i) + 1 << " ";
                }
                out << std::endl;
            }

            if(
                ioflags.has_attribute(MESH_FACET_REGION) &&
                M.has_attribute(MESH_FACET_REGION)
            ) {
                out << "# attribute chart facet integer" << std::endl;
                for(index_t f = 0; f < M.nb_facets(); f++) {
                    out << "# attrs f "
                        << f + 1 << " "
                        << M.facet_region(f) << std::endl;
                }
            }

            return true;
        }

        virtual bool load(
            const std::string& filename, Mesh& M,
            const MeshIOFlags& ioflags = MeshIOFlags()
        ) {
            return load_generic(filename, M, ioflags);
        }

        virtual bool load(
            const std::string& filename, SinglePrecisionMesh& M,
            const MeshIOFlags& ioflags = MeshIOFlags()
        ) {
            return load_generic(filename, M, ioflags);
        }

    protected:
        virtual ~OBJIOHandler() {
        }

    private:
        coord_index_t dimension_;
    };

    /************************************************************************/

    /**
     * \brief IO handler for the OBJ6 file format
     * \see OBJIOHandler
     */
    class GEOGRAM_API OBJ6IOHandler : public OBJIOHandler {
    public:
        OBJ6IOHandler() :
            OBJIOHandler(6) {
        }

    protected:
        virtual ~OBJ6IOHandler() {
        }
    };

    /************************************************************************/

    /**
     * \brief IO handler for LM5/LM6/Gamma mesh file format
     * \see http://www-roc.inria.fr/gamma/gamma/Membres/CIPD/Loic.Marechal/Research/LM5.html
     */
    class GEOGRAM_API LMIOHandler : public MeshIOHandler {
    public:
        /**
         * \brief Loads a mesh in LM5/LM6/Gamma mesh file format
         * \param[in] filename name of the file
         * \param[out] M the loaded mesh
         * \param[in] ioflags specifies which attributes and elements should be read
         * \tparam MESH the mesh class
         * \return true on success, false otherwise
         */
        template <class MESH>
        bool load_generic(
            const std::string& filename, MESH& M, const MeshIOFlags& ioflags
        ) {
            typedef typename GEOGen::MeshBuilder<MESH>::coord_t coord_t;
            GEOGen::MeshBuilder<MESH> builder(&M);

            int ver, dim;
            int InpMsh = GmfOpenMesh(
                filename.c_str(), GmfRead, &ver, &dim
            );
            if(!InpMsh) {
                Logger::err("I/O") << "Could not open file: "
                    << filename << std::endl;
                return false;
            }

            if(ver != GmfFloat && ver != GmfDouble) {
                Logger::err("I/O") << "Invalid version: " << ver << std::endl;
                GmfCloseMesh(InpMsh);
                return false;
            }

            bool use_doubles = (ver == GmfDouble);

            if(dim != 3) {
                Logger::err("I/O") << "Invalid dimension: " << dim << std::endl;
                GmfCloseMesh(InpMsh);
                return false;
            }

            long NmbVer = GmfStatKwd(InpMsh, GmfVertices);
            long NmbTris = GmfStatKwd(InpMsh, GmfTriangles);
            long NmbQuad = GmfStatKwd(InpMsh, GmfQuadrilaterals);
            long NmbTet = GmfStatKwd(InpMsh, GmfTetrahedra);
            long NmbHex = GmfStatKwd(InpMsh, GmfHexahedra);
            long NmbPrisms = GmfStatKwd(InpMsh, GmfPrisms);
            long NmbPyramids = GmfStatKwd(InpMsh, GmfPyramids);            
            long NmbCells = NmbTet+NmbHex+NmbPrisms+NmbPyramids;
            
            if(NmbVer == 0) {
                Logger::err("I/O") << "File contains no vertices" << std::endl;
                GmfCloseMesh(InpMsh);
                return false;
            }

            if(NmbTris == 0 && NmbQuad == 0 && ioflags.has_element(MESH_FACETS)) {
                Logger::warn("I/O") << "File contains no faces" << std::endl;
            }

            if((NmbCells > 0) && !ioflags.has_element(MESH_CELLS)) {
                Logger::warn("I/O") << "File seems to contain a volume"
                    << std::endl;
            }

            builder.begin_mesh();

            builder.create_vertices(index_t(NmbVer), 3);
            
            // Read vertices
            if(!GmfGotoKwd(InpMsh, GmfVertices)) {
                Logger::err("I/O") << "Failed to access Vertices section"
                    << std::endl;
                GmfCloseMesh(InpMsh);
                return false;
            }

            if(use_doubles) {
                for(index_t i = 0; i < index_t(NmbVer); i++) {
                    double x, y, z;
                    int ref;
                    if(!GmfGetLin(InpMsh, GmfVertices, &x, &y, &z, &ref)) {
                        Logger::err("I/O") << "Failed to read vertex #" << i
                            << std::endl;
                        GmfCloseMesh(InpMsh);
                        return false;
                    }
                    M.vertex_ptr(i)[0] = coord_t(x);
                    M.vertex_ptr(i)[1] = coord_t(y);
                    M.vertex_ptr(i)[2] = coord_t(z);
                }
            } else {
                for(index_t i = 0; i < index_t(NmbVer); i++) {
                    float x, y, z;
                    int ref;
                    if(!GmfGetLin(InpMsh, GmfVertices, &x, &y, &z, &ref)) {
                        Logger::err("I/O") << "Failed to read vertex #" << i
                            << std::endl;
                        GmfCloseMesh(InpMsh);
                        return false;
                    }
                    M.vertex_ptr(i)[0] = coord_t(x);
                    M.vertex_ptr(i)[1] = coord_t(y);
                    M.vertex_ptr(i)[2] = coord_t(z);
                }
            }

            if(ioflags.has_element(MESH_FACETS)) {
                builder.reserve_facets(index_t(NmbTris + NmbQuad));
                // Read triangles
                if(NmbTris > 0) {
                    if(!GmfGotoKwd(InpMsh, GmfTriangles)) {
                        Logger::err("I/O")
                            << "Failed to access Triangles section"
                            << std::endl;
                        GmfCloseMesh(InpMsh);
                        return false;
                    }

                    for(index_t i = 0; i < index_t(NmbTris); i++) {
                        int v[3], ref;
                        int res = GmfGetLin(
                            InpMsh, GmfTriangles, &v[0], &v[1], &v[2], &ref
                        );
                        if(!res) {
                            Logger::err("I/O")
                                << "Failed to read triangle facet #" << i
                                << std::endl;
                            GmfCloseMesh(InpMsh);
                            return false;
                        }
                        if(
                            ref != 0 &&
                            ioflags.has_attribute(MESH_FACET_REGION)
                        ) {
                            create_attributes(M, MESH_FACET_REGION);
                        }
                        builder.begin_facet();
                        for(index_t j = 0; j < 3; ++j) {
                            if(
                                v[j] < 1 ||
                                index_t(v[j]) > builder.target()->nb_vertices()
                            ) {
                                Logger::err("I/O")
                                    << "Error: triangle facet #" << i
                                    << " references an invalid vertex: " << v[j]
                                    << std::endl;
                                GmfCloseMesh(InpMsh);
                                return false;
                            }
                            builder.add_vertex_to_facet(index_t(v[j] - 1));
                        }
                        builder.end_facet(ref);
                    }
                }

                // Read quads
                if(NmbQuad > 0) {
                    if(!GmfGotoKwd(InpMsh, GmfQuadrilaterals)) {
                        Logger::err("I/O")
                            << "Failed to access Quadrilaterals section"
                            << std::endl;
                        GmfCloseMesh(InpMsh);
                        return false;
                    }

                    for(index_t i = 0; i < index_t(NmbQuad); i++) {
                        int v[4], ref;
                        int res = GmfGetLin(
                            InpMsh, GmfQuadrilaterals,
                            &v[0], &v[1], &v[2], &v[3], &ref
                        );
                        if(!res) {
                            Logger::err("I/O")
                                << "Failed to read quadrilateral facet #" << i
                                << std::endl;
                            GmfCloseMesh(InpMsh);
                            return false;
                        }
                        if(
                            ref != 0 &&
                            ioflags.has_attribute(MESH_FACET_REGION)
                        ) {
                            create_attributes(M, MESH_FACET_REGION);
                        }
                        builder.begin_facet();
                        for(index_t j = 0; j < 4; ++j) {
                            if(
                                v[j] < 1 ||
                                index_t(v[j]) > builder.target()->nb_vertices()
                            ) {
                                Logger::err("I/O")
                                    << "Error: quadrilateral facet #" << i
                                    << " references an invalid vertex: " << v[j]
                                    << std::endl;
                                GmfCloseMesh(InpMsh);
                                return false;
                            }
                            builder.add_vertex_to_facet(index_t(v[j] - 1));
                        }
                        builder.end_facet(ref);
                    }
                }
            }

            builder.end_mesh(false);
            
            if(ioflags.has_element(MESH_CELLS)) {
                if(NmbTet > 0) {
                    if(!GmfGotoKwd(InpMsh, GmfTetrahedra)) {
                        Logger::err("I/O")
                            << "Failed to access Tetrahedra section"
                            << std::endl;
                        GmfCloseMesh(InpMsh);
                        return false;
                    }

                    builder.target()->create_tets(index_t(NmbTet));
                    for(index_t i = 0; i < index_t(NmbTet); i++) {
                        int v[4], ref;
                        int res = GmfGetLin(
                            InpMsh, GmfTetrahedra,
                            &v[0], &v[1], &v[2], &v[3], &ref
                        );
                        if(!res) {
                            Logger::err("I/O") << "Failed to read tetrahedron #" << i
                                << std::endl;
                            GmfCloseMesh(InpMsh);
                            return false;
                        }
                        for(index_t j = 0; j < 4; ++j) {
                            if(
                                v[j] < 1 ||
                                index_t(v[j]) > builder.target()->nb_vertices()
                            ) {
                                Logger::err("I/O")
                                    << "Error: tetrahedron  #" << i
                                    << " references an invalid vertex: " << v[j]
                                    << std::endl;
                                GmfCloseMesh(InpMsh);
                                return false;
                            }
                            builder.target()->set_tet_vertex_index(i, j, index_t(v[j] - 1));
                        }
                        // TODO: tet regions.
                    }
                }

                if(NmbHex > 0) {
                    if(!GmfGotoKwd(InpMsh, GmfHexahedra)) {
                        Logger::err("I/O")
                            << "Failed to access Hexahedra section"
                            << std::endl;
                        GmfCloseMesh(InpMsh);
                        return false;
                    }
                    
                    for(index_t i = 0; i < index_t(NmbHex); i++) {
                        int v[8], ref;

                        // Swapping vertices 1<->0 and 4<->5 to
                        // account for differences in the indexing
                        // convetions in .mesh/.meshb files w.r.t.
                        // geogram internal conventions.
                        int res = GmfGetLin(
                            InpMsh, GmfHexahedra,
                            &v[1], &v[0], &v[2], &v[3],
                            &v[5], &v[4], &v[6], &v[7],                            
                            &ref
                        );
                        if(!res) {
                            Logger::err("I/O") << "Failed to read hexahedron #" << i
                                << std::endl;
                            GmfCloseMesh(InpMsh);
                            return false;
                        }
                        for(index_t j = 0; j < 8; ++j) {
                            if(
                                v[j] < 1 ||
                                index_t(v[j]) > builder.target()->nb_vertices()
                            ) {
                                Logger::err("I/O")
                                    << "Error: hexahedron  #" << i
                                    << " references an invalid vertex: " << v[j]
                                    << std::endl;
                                GmfCloseMesh(InpMsh);
                                return false;
                            }
                        }
                        builder.target()->add_hex(
                            index_t(v[0]-1), index_t(v[1]-1), index_t(v[2]-1), index_t(v[3]-1),
                            index_t(v[4]-1), index_t(v[5]-1), index_t(v[6]-1), index_t(v[7]-1),
                            -1,-1,-1,-1,-1,-1, // faces adjacencies
                            ref
                       );
                    }
                }
                
                if(NmbPrisms > 0) {
                    if(!GmfGotoKwd(InpMsh, GmfPrisms)) {
                        Logger::err("I/O")
                            << "Failed to access Prisms section"
                            << std::endl;
                        GmfCloseMesh(InpMsh);
                        return false;
                    }
                    
                    for(index_t i = 0; i < index_t(NmbPrisms); i++) {
                        int v[6], ref;
                        int res = GmfGetLin(
                            InpMsh, GmfPrisms,
                            &v[0], &v[1], &v[2],
                            &v[3], &v[4], &v[5], 
                            &ref
                        );
                        if(!res) {
                            Logger::err("I/O") << "Failed to read prism #" << i
                                << std::endl;
                            GmfCloseMesh(InpMsh);
                            return false;
                        }
                        for(index_t j = 0; j < 6; ++j) {
                            if(
                                v[j] < 1 ||
                                index_t(v[j]) > builder.target()->nb_vertices()
                            ) {
                                Logger::err("I/O")
                                    << "Error: prism  #" << i
                                    << " references an invalid vertex: " << v[j]
                                    << std::endl;
                                GmfCloseMesh(InpMsh);
                                return false;
                            }
                        }
                        builder.target()->add_prism(
                            index_t(v[0]-1), index_t(v[1]-1), index_t(v[2]-1), 
                            index_t(v[3]-1), index_t(v[4]-1), index_t(v[5]-1),
                            -1,-1,-1,-1,-1, // faces adjacencies                            
                            ref
                       );
                    }
                }

                if(NmbPyramids > 0) {
                    if(!GmfGotoKwd(InpMsh, GmfPyramids)) {
                        Logger::err("I/O")
                            << "Failed to access Pyramids section"
                            << std::endl;
                        GmfCloseMesh(InpMsh);
                        return false;
                    }
                    
                    for(index_t i = 0; i < index_t(NmbPyramids); i++) {
                        int v[5], ref;
                        int res = GmfGetLin(
                            InpMsh, GmfPyramids,
                            &v[0], &v[1], &v[2], &v[3], &v[4],
                            &ref
                        );
                        if(!res) {
                            Logger::err("I/O") << "Failed to read pyramid #" << i
                                << std::endl;
                            GmfCloseMesh(InpMsh);
                            return false;
                        }
                        for(index_t j = 0; j < 5; ++j) {
                            if(
                                v[j] < 1 ||
                                index_t(v[j]) > builder.target()->nb_vertices()
                            ) {
                                Logger::err("I/O")
                                    << "Error: prism  #" << i
                                    << " references an invalid vertex: " << v[j]
                                    << std::endl;
                                GmfCloseMesh(InpMsh);
                                return false;
                            }
                        }
                        builder.target()->add_pyramid(
                            index_t(v[0]-1), index_t(v[1]-1), index_t(v[2]-1), 
                            index_t(v[3]-1), index_t(v[4]-1),
                            -1,-1,-1,-1,-1, // faces adjacencies  
                            ref
                       );
                    }
                }

                builder.target()->connect_cells(); 
            } 
            GmfCloseMesh(InpMsh);
            return true;
        }

        virtual bool save(
            const Mesh& M, const std::string& filename,
            const MeshIOFlags& ioflags
        ) {

            bool use_doubles = CmdLine::get_arg_bool("sys:use_doubles");
            int OutMsh = GmfOpenMesh(
                filename.c_str(), GmfWrite,
                (use_doubles ? GmfDouble : GmfFloat), 3
            );

            if(OutMsh == 0) {
                Logger::err("I/O")
                    << "Could not create file \'" << filename << "\'"
                    << std::endl;
                return false;
            }

            {
                index_t nb_vertices = M.nb_vertices();
                GmfSetKwd(OutMsh, GmfVertices, nb_vertices);
                for(index_t i = 0; i < nb_vertices; i++) {
                    const double* p = M.vertex_ptr(i);
                    // Note: GmfSetLin uses variadic arguments and thus expects
                    // double arguments.
                    GmfSetLin(OutMsh, GmfVertices, p[0], p[1], p[2], 0);
                }
            }

            int nb_tris = 0;
            int nb_quads = 0;
            int nb_other = 0;
            for(index_t f = 0; f < M.nb_facets(); f++) {
                switch(M.facet_size(f)) {
                    case 3:
                        nb_tris++;
                        break;
                    case 4:
                        nb_quads++;
                        break;
                    default:
                        nb_other++;
                        break;
                }
            }

            bool save_facet_regions = ioflags.has_attribute(MESH_FACET_REGION);

            if(nb_tris > 0) {
                GmfSetKwd(OutMsh, GmfTriangles, nb_tris);
                for(index_t f = 0; f < M.nb_facets(); f++) {
                    if(M.facet_size(f) == 3) {
                        index_t v1 = M.corner_vertex_index(M.facet_begin(f));
                        index_t v2 = M.corner_vertex_index(M.facet_begin(f)+1);
                        index_t v3 = M.corner_vertex_index(M.facet_begin(f)+2);
                        GmfSetLin(
                            OutMsh, GmfTriangles,
                            int(v1 + 1), int(v2 + 1), int(v3 + 1),
                            save_facet_regions ? M.facet_region(f) : 0
                        );
                    }
                }
            }

            if(nb_quads > 0) {
                GmfSetKwd(OutMsh, GmfQuadrilaterals, nb_quads);
                for(index_t f = 0; f < M.nb_facets(); f++) {
                    if(M.facet_size(f) == 4) {
                        index_t v1 = M.corner_vertex_index(M.facet_begin(f));
                        index_t v2 = M.corner_vertex_index(M.facet_begin(f)+1);
                        index_t v3 = M.corner_vertex_index(M.facet_begin(f)+2);
                        index_t v4 = M.corner_vertex_index(M.facet_begin(f)+3);
                        GmfSetLin(
                            OutMsh, GmfQuadrilaterals,
                            int(v1 + 1), int(v2 + 1), int(v3 + 1), int(v4 + 1),
                            save_facet_regions ? M.facet_region(f) : 0
                        );
                    }
                }
            }

            if(nb_other > 0) {
                Logger::warn("I/O")
                    << "Encountered " << nb_other 
                    << " non-tri / non-quad facets"
                    << std::endl;
            }

            if(M.nb_cells() > 0 && ioflags.has_element(MESH_CELLS)) {
                
                bool save_cell_regions =
                    ioflags.has_attribute(MESH_CELL_REGION) &&
                    M.has_attribute(MESH_CELL_REGION);
                
                if(M.is_tetrahedralized()) {
                    GmfSetKwd(OutMsh, GmfTetrahedra, M.nb_tets());
                    for(index_t t = 0; t < M.nb_tets(); ++t) {
                        GmfSetLin(
                            OutMsh, GmfTetrahedra,
                            M.tet_vertex_index(t, 0) + 1,
                            M.tet_vertex_index(t, 1) + 1,
                            M.tet_vertex_index(t, 2) + 1,
                            M.tet_vertex_index(t, 3) + 1,
                            save_cell_regions ? M.tet_region(t) : 0
                        );
                    }
                } else {
                    index_t nb_tets=0;
                    index_t nb_hexes=0;
                    index_t nb_prisms=0;
                    index_t nb_pyramids=0;
                    for(index_t c=0; c<M.nb_cells(); ++c) {
                        switch(M.cell_type(c)) {
                        case MESH_TET:
                            ++nb_tets;
                            break;
                        case MESH_HEX:
                            ++nb_hexes;
                            break;
                        case MESH_PRISM:
                            ++nb_prisms;
                            break;
                        case MESH_PYRAMID:
                            ++nb_pyramids;
                            break;
                        default:
                            break;
                        }
                    }
                    if(nb_tets > 0) {
                        GmfSetKwd(OutMsh, GmfTetrahedra, nb_tets);
                        for(index_t c=0; c<M.nb_cells(); ++c) {
                            if(M.cell_type(c) == MESH_TET) {
                                GmfSetLin(
                                    OutMsh, GmfTetrahedra,
                                    M.cell_vertex_index(c, 0) + 1,
                                    M.cell_vertex_index(c, 1) + 1,
                                    M.cell_vertex_index(c, 2) + 1,
                                    M.cell_vertex_index(c, 3) + 1,
                                    save_cell_regions ? M.cell_region(c) : 0
                                );
                            }
                        }
                    }
                    if(nb_hexes > 0) {
                        GmfSetKwd(OutMsh, GmfHexahedra, nb_hexes);
                        for(index_t c=0; c<M.nb_cells(); ++c) {
                            if(M.cell_type(c) == MESH_HEX) {
                                // Swapping vertices 1<->0 and 4<->5 to
                                // account for differences in the indexing
                                // convetions in .mesh/.meshb files w.r.t.
                                // geogram internal conventions.
                                GmfSetLin(
                                    OutMsh, GmfHexahedra,
                                    M.cell_vertex_index(c, 1) + 1,
                                    M.cell_vertex_index(c, 0) + 1,
                                    M.cell_vertex_index(c, 2) + 1,
                                    M.cell_vertex_index(c, 3) + 1,
                                    M.cell_vertex_index(c, 5) + 1,
                                    M.cell_vertex_index(c, 4) + 1,
                                    M.cell_vertex_index(c, 6) + 1,
                                    M.cell_vertex_index(c, 7) + 1,
                                    save_cell_regions ? M.cell_region(c) : 0
                                );
                            }
                        }
                    }
                    if(nb_prisms > 0) {
                        GmfSetKwd(OutMsh, GmfPrisms, nb_prisms);
                        for(index_t c=0; c<M.nb_cells(); ++c) {
                            if(M.cell_type(c) == MESH_PRISM) {
                                GmfSetLin(
                                    OutMsh, GmfPrisms,
                                    M.cell_vertex_index(c, 0) + 1,
                                    M.cell_vertex_index(c, 1) + 1,
                                    M.cell_vertex_index(c, 2) + 1,
                                    M.cell_vertex_index(c, 3) + 1,
                                    M.cell_vertex_index(c, 4) + 1,
                                    M.cell_vertex_index(c, 5) + 1,
                                    save_cell_regions ? M.cell_region(c) : 0
                                );
                            }
                        }
                    }
                    if(nb_pyramids > 0) {
                        GmfSetKwd(OutMsh, GmfPyramids, nb_pyramids);
                        for(index_t c=0; c<M.nb_cells(); ++c) {
                            if(M.cell_type(c) == MESH_PYRAMID) {
                                GmfSetLin(
                                    OutMsh, GmfPyramids,
                                    M.cell_vertex_index(c, 0) + 1,
                                    M.cell_vertex_index(c, 1) + 1,
                                    M.cell_vertex_index(c, 2) + 1,
                                    M.cell_vertex_index(c, 3) + 1,
                                    M.cell_vertex_index(c, 4) + 1,
                                    save_cell_regions ? M.cell_region(c) : 0
                                );
                            }
                        }
                    }
                }
            }
            GmfCloseMesh(OutMsh);

            // If file is in ASCII, append parameters as comments
            // at the end of the file.
            if(FileSystem::extension(filename) == "mesh") {
                FILE* f = fopen(filename.c_str(), "a");
                std::vector<std::string> args;
                CmdLine::get_args(args);
                for(index_t i = 0; i < args.size(); i++) {
                    fprintf(f, "# vorpaline %s\n", args[i].c_str());
                }
                fclose(f);
            }

            return true;
        }

        virtual bool load(
            const std::string& filename, Mesh& M,
            const MeshIOFlags& ioflags = MeshIOFlags()
        ) {
            return load_generic(filename, M, ioflags);
        }

        virtual bool load(
            const std::string& filename, SinglePrecisionMesh& M,
            const MeshIOFlags& ioflags = MeshIOFlags()
        ) {
            return load_generic(filename, M, ioflags);
        }

    protected:
        virtual ~LMIOHandler() {
        }
    };

    /************************************************************************/

    /**
     * \brief IO handler for the PLY file format
     * \details ASCII and binary, single and double precision are supported
     * \see http://en.wikipedia.org/wiki/PLY_(file_format)
     */
    class GEOGRAM_API PLYIOHandler : public MeshIOHandler {
    public:
        /**
         * \brief Helper class to read files in PLY format
         * \tparam MESH the mesh class
         */
        template <class MESH>
        class PlyLoader {
        public:
            /** The type for vertices coordinates. */
            typedef typename MESH::coord_t coord_t;

            /**
             * \brief Construts a new PlyLoader.
             * \param[in] filename name of the file to be loaded
             * \param[out] M the loaded mesh
             */
            PlyLoader(const std::string& filename, MESH& M) :
                builder_(&M),
                filename_(filename),
                current_vertex_(max_index_t()),
                current_vertex_ptr_(nil),
                has_colors_(false),
                color_mult_(1.0),
                current_color_(max_index_t()),
                current_color_ptr_(nil),
                tristrip_index_(0),
                load_colors_(false) {
            }

            /**
             * \brief Specifies whether vertices colors should be loaded.
             */
            void set_load_colors(bool x) {
                load_colors_ = x;
            }

            /**
             * \brief Loads the file.
             * \return true on success, false otherwise
             */
            bool load() {
                p_ply ply = ply_open(filename_.c_str(), nil, 0, nil);

                if(ply == nil) {
                    Logger::err("I/O")
                        << "Could not open file: " << filename_
                        << std::endl;
                    return false;
                }

                if(!ply_read_header(ply)) {
                    Logger::err("I/O")
                        << "Invalid PLY header"
                        << std::endl;
                    ply_close(ply);
                    return false;
                }

                current_vertex_ = 0;
                current_color_ = 0;
                if(load_colors_) {
                    check_for_colors(ply);
                } else {
                    has_colors_ = false;
                }

                long nvertices = ply_set_read_cb(
                    ply, "vertex", "x", PlyLoader::vertex_cb, this, 0
                );
                ply_set_read_cb(
                    ply, "vertex", "y", PlyLoader::vertex_cb, this, 1
                );
                ply_set_read_cb(
                    ply, "vertex", "z", PlyLoader::vertex_cb, this, 2
                );

                long nfaces = ply_set_read_cb(
                    ply, "face", "vertex_indices", PlyLoader::face_cb, this, 0
                );
                nfaces += ply_set_read_cb(
                    ply, "face", "vertex_index", PlyLoader::face_cb, this, 0
                );
                long ntstrips = ply_set_read_cb(
                    ply, "tristrips", "vertex_indices",
                    PlyLoader::tristrip_cb, this, 0
                );
                ntstrips += ply_set_read_cb(
                    ply, "tristrips", "vertex_index",
                    PlyLoader::tristrip_cb, this, 0
                );

                if(nvertices == 0) {
                    Logger::err("I/O") 
                        << "File contains no vertices" << std::endl;
                    ply_close(ply);
                    return false;
                }

                builder_.begin_mesh();
                builder_.create_vertices(
                    index_t(nvertices), has_colors_ ? 9 : 3
                );
                builder_.reserve_facets(index_t(nfaces));

                if(!ply_read(ply)) {
                    Logger::err("I/O")
                        << "Problem occured while parsing PLY file"
                        << std::endl;
                    ply_close(ply);
                    builder_.end_mesh(false);
                    return false;
                }

                ply_close(ply);
                builder_.end_mesh(false);

                if(current_vertex_ != builder_.target()->nb_vertices()) {
                    Logger::err("I/O")
                        << "File does not contain enough vertex data"
                        << std::endl;
                    return false;
                }

                if(
                    has_colors_ &&
                    current_color_ != builder_.target()->nb_vertices()
                ) {
                    Logger::err("I/O")
                        << "File does not contain enough color data"
                        << std::endl;
                    return false;
                }

                return true;
            }

        protected:
            /**
             * \brief Detects whether the input file has colors
             * \param[in] ply the p_ply handle to the file
             */
            void check_for_colors(p_ply ply) {
                p_ply_element element = nil;

                bool has_r = false;
                bool has_g = false;
                bool has_b = false;

                bool has_red = false;
                bool has_green = false;
                bool has_blue = false;

                for(;;) {
                    element = ply_get_next_element(ply, element);
                    if(element == nil) {
                        break;
                    }
                    const char* elt_name = nil;
                    ply_get_element_info(element, &elt_name, nil);

                    if(!strcmp(elt_name, "vertex")) {
                        p_ply_property property = nil;
                        for(;;) {
                            property = ply_get_next_property(element, property);
                            if(property == nil) {
                                break;
                            }
                            const char* prop_name = nil;
                            ply_get_property_info(
                                property, &prop_name, nil, nil, nil
                            );
                            has_r = has_r || !strcmp(prop_name, "r");
                            has_g = has_g || !strcmp(prop_name, "g");
                            has_b = has_b || !strcmp(prop_name, "b");
                            has_red = has_red || !strcmp(prop_name, "red");
                            has_green = has_green || !strcmp(prop_name, "green");
                            has_blue = has_blue || !strcmp(prop_name, "blue");
                        }
                    }
                }

                if(has_r && has_g && has_b) {
                    has_colors_ = true;
                    color_mult_ = 1.0;
                    ply_set_read_cb(
                        ply, "vertex", "r", PlyLoader::color_cb, this, 0
                    );
                    ply_set_read_cb(
                        ply, "vertex", "g", PlyLoader::color_cb, this, 1
                    );
                    ply_set_read_cb(
                        ply, "vertex", "b", PlyLoader::color_cb, this, 2
                    );
                } else if(has_red && has_green && has_blue) {
                    has_colors_ = true;
                    color_mult_ = 1.0 / 255.0;
                    ply_set_read_cb(
                        ply, "vertex", "red", PlyLoader::color_cb, this, 0
                    );
                    ply_set_read_cb(
                        ply, "vertex", "green", PlyLoader::color_cb, this, 1
                    );
                    ply_set_read_cb(
                        ply, "vertex", "blue", PlyLoader::color_cb, this, 2
                    );
                } else {
                    has_colors_ = false;
                }
            }

            /**
             * \brief Gets the PlyLoader associated with an opaque p_ply_argument.
             * \details Used to pass a Plyloader through libply callbacks
             * \param[in] argument the opaque p_ply_argument
             * \return the PlyLoader associated with \p argument
             */
            static PlyLoader* loader(p_ply_argument argument) {
                PlyLoader* result = nil;
                ply_get_argument_user_data(argument, (void**) (&result), nil);
                geo_debug_assert(result != nil);
                return result;
            }

            /**
             * \brief The vertex callback, called for each vertex of the input file.
             * \param[in] argument the generic opaque argument (from which this PlyLoader is retreived).
             * \return callback status code, zero for errors, non-zero for success.
             */
            static int vertex_cb(p_ply_argument argument) {
                return loader(argument)->add_vertex_data(argument);
            }

            /**
             * \brief The facet callback, called for each facet of the input file.
             * \param[in] argument the generic opaque argument (from which this PlyLoader is retreived).
             * \return callback status code, zero for errors, non-zero for success.
             */
            static int face_cb(p_ply_argument argument) {
                return loader(argument)->add_face_data(argument);
            }

            /**
             * \brief The triangle strip callback, called for each triangle strip of the input file.
             * \param[in] argument the generic opaque argument (from which this PlyLoader is retreived).
             * \return callback status code, zero for errors, non-zero for success.
             */
            static int tristrip_cb(p_ply_argument argument) {
                return loader(argument)->add_tristrip_data(argument);
            }

            /**
             * \brief The color callback, called for each color data of the input file.
             * \param[in] argument the generic opaque argument (from which this PlyLoader is retreived).
             * \return callback status code, zero for errors, non-zero for success.
             */
            static int color_cb(p_ply_argument argument) {
                return loader(argument)->add_color_data(argument);
            }

            /**
             * \brief Decodes vertex data from a generic callback argument.
             * \param[in] argument the generic callback argument.
             * \return callback status code, zero for errors, non-zero for success.
             */
            int add_vertex_data(p_ply_argument argument) {
                long coord;
                ply_get_argument_user_data(argument, nil, &coord);
                switch(coord) {
                    case 0:
                        geo_debug_assert(builder_.target()->dimension() >= 3);
                        if(current_vertex_ >= builder_.target()->nb_vertices()) {
                            Logger::err("I/O")
                                << "File contains extraneous vertex data"
                                << std::endl;
                            return 0;
                        }
                        current_vertex_ptr_ = builder_.target()->vertex_ptr(current_vertex_++);
                    // PASS THROUGH
                    case 1:
                    case 2:
                        current_vertex_ptr_[coord] = coord_t(ply_get_argument_value(argument));
                        return 1;
                }

                Logger::err("I/O")
                    << "In vertex #" << current_vertex_
                    << ": invalid coordinate index: " << coord
                    << std::endl;
                return 0;
            }

            /**
             * \brief Decodes facet data from a generic callback argument.
             * \param[in] argument the generic callback argument.
             * \return callback status code, zero for errors, non-zero for success.
             */
            int add_face_data(p_ply_argument argument) {
                long length, value_index;
                ply_get_argument_property(argument, nil, &length, &value_index);
                if(value_index < 0) {
                    // Ignore negative values - this is not an error! (facet markers)
                    return 1;
                }

                index_t vertex_index = index_t(ply_get_argument_value(argument));
                if(index_t(vertex_index) >= builder_.target()->nb_vertices()) {
                    Logger::err("I/O")
                        << "Facet corner #" << builder_.target()->nb_facets()
                        << " references an invalid vertex: "
                        << vertex_index
                        << std::endl;
                    return 0;
                }

                if(value_index == 0) {
                    builder_.begin_facet();
                }
                builder_.add_vertex_to_facet(vertex_index);
                if(value_index == length - 1) {
                    builder_.end_facet();
                }
                return 1;
            }

            /**
             * \brief Decodes triangle strip data from a generic callback argument.
             * \param[in] argument the generic callback argument.
             * \return callback status code, zero for errors, non-zero for success.
             */
            int add_tristrip_data(p_ply_argument argument) {
                long length, value_index;
                ply_get_argument_property(argument, nil, &length, &value_index);
                if(value_index < 0) {
                    // Ignore negative values - this is not an error!
                    return 1;
                }

                // NOTE: negative vertex_index values have a special meaning here:
                // they tell the loader to start a new strip

                signed_index_t vertex_index = signed_index_t(ply_get_argument_value(argument));
                if(vertex_index >= signed_index_t(builder_.target()->nb_vertices())) {
                    Logger::err("I/O")
                        << "Invalid vertex reference in tristrip: "
                        << vertex_index
                        << std::endl;
                    return 0;
                }

                if(value_index == 0) {
                    begin_tristrip();
                }
                if(vertex_index >= 0) {
                    add_to_tristrip(index_t(vertex_index));
                } else {
                    end_tristrip();
                    begin_tristrip();
                }
                if(value_index == length - 1) {
                    end_tristrip();
                }
                return 1;
            }

            /**
             * \brief Starts a new triangle strip.
             */
            void begin_tristrip() {
                tristrip_index_ = 0;
            }

            /**
             * \brief Terminates the current triangle strip.
             */
            void end_tristrip() {
            }

            /**
             * \brief Adds a vertex to the current triangle strip.
             * \param[in] vertex_index the index of the vertex
             */
            void add_to_tristrip(index_t vertex_index) {
                if(tristrip_index_ >= 2) {
                    builder_.begin_facet();
                    builder_.add_vertex_to_facet(tristrip_points_[0]);
                    builder_.add_vertex_to_facet(tristrip_points_[1]);
                    builder_.add_vertex_to_facet(vertex_index);
                    builder_.end_facet();
                }

                tristrip_points_[tristrip_index_ & 1] = vertex_index;
                tristrip_index_++;
            }

            /**
             * \brief Adds color data to the current vertex
             * \param[in] argument the generic callback argument
             * \return callback status code, zero for errors, non-zero for success.
             */
            int add_color_data(p_ply_argument argument) {
                long coord;
                ply_get_argument_user_data(argument, nil, &coord);
                double value;
                switch(coord) {
                    case 0:
                        geo_debug_assert(builder_.target()->dimension() >= 9);
                        if(current_color_ >= builder_.target()->nb_vertices()) {
                            Logger::err("I/O")
                                << "File contains extraneous color data"
                                << std::endl;
                            return 0;
                        }
                        current_color_ptr_ = builder_.target()->vertex_ptr(current_color_++) + 6;
                    // PASS THROUGH
                    case 1:
                    case 2:
                        value = double(ply_get_argument_value(argument)) * color_mult_;
                        current_color_ptr_[coord] = coord_t(value);
                        return 1;
                }

                Logger::err("I/O")
                    << "In vertex #" << current_color_
                    << ": invalid color index: " << coord
                    << std::endl;
                return 0;
            }

        protected:
            GEOGen::MeshBuilder<MESH> builder_;
            std::string filename_;

            index_t current_vertex_;
            coord_t* current_vertex_ptr_;

            bool has_colors_;
            double color_mult_;
            index_t current_color_;
            coord_t* current_color_ptr_;

            index_t tristrip_points_[2];
            index_t tristrip_index_;

            bool load_colors_;
        };

        virtual bool save(
            const Mesh& M, const std::string& filename, const MeshIOFlags& ioflags
        ) {
            geo_argused(ioflags);
            p_ply oply = ply_create(
                filename.c_str(), PLY_LITTLE_ENDIAN, NULL, 0, NULL
            );

            if(oply == nil) {
                return false;
            }

            // Create element and properties for vertices
            e_ply_type coord_type = PLY_FLOAT;
            ply_add_element(oply, "vertex", long(M.nb_vertices()));
            ply_add_property(oply, "x", coord_type, coord_type, coord_type);
            ply_add_property(oply, "y", coord_type, coord_type, coord_type);
            ply_add_property(oply, "z", coord_type, coord_type, coord_type);

            // Create element and properties for facets (determine best index types)
            index_t max_facet_size = 0;
            for(index_t f = 0; f < M.nb_facets(); f++) {
                max_facet_size = geo_max(max_facet_size, M.facet_size(f));
            }
            e_ply_type facet_len_type = PLY_UCHAR;
            if(max_facet_size > 65535) {
                facet_len_type = PLY_UINT;
            } else if(max_facet_size > 255) {
                facet_len_type = PLY_USHORT;
            }
            e_ply_type facet_idx_type = PLY_INT;
            ply_add_element(oply, "face", long(M.nb_facets()));
            ply_add_property(
                oply, "vertex_indices", PLY_LIST, facet_len_type, facet_idx_type
            );

            std::vector<std::string> args;
            CmdLine::get_args(args);
            for(index_t i = 0; i < args.size(); i++) {
                ply_add_comment(oply, ("vorpaline " + args[i]).c_str());
            }

            // Write header
            if(!ply_write_header(oply)) {
                ply_close(oply);
                return false;
            }

            // Write vertices
            for(index_t v = 0; v < M.nb_vertices(); v++) {
                for(index_t coord = 0; coord < 3; coord++) {
                    ply_write(oply, M.vertex_ptr(v)[coord]);
                }
            }

            // Write facets
            for(index_t f = 0; f < M.nb_facets(); f++) {
                ply_write(oply, double(M.facet_size(f)));
                for(index_t c = M.facet_begin(f); c < M.facet_end(f); c++) {
                    ply_write(oply, double(M.corner_vertex_index(c)));
                }
            }

            ply_close(oply);
            return true;
        }

        virtual bool load(
            const std::string& filename, Mesh& M,
            const MeshIOFlags& ioflags = MeshIOFlags()
        ) {
            geo_argused(ioflags);
            PlyLoader<Mesh> loader(filename, M);
            return loader.load();
        }

        virtual bool load(
            const std::string& filename, SinglePrecisionMesh& M,
            const MeshIOFlags& ioflags = MeshIOFlags()
        ) {
            geo_argused(ioflags);
            PlyLoader<SinglePrecisionMesh> loader(filename, M);
            return loader.load();
        }

    protected:
        virtual ~PLYIOHandler() {
        }
    };

    /************************************************************************/

    /**
     * \brief IO handler for the OFF file format
     * \see http://www.geomview.org/docs/html/OFF.html
     */
    class GEOGRAM_API OFFIOHandler : public MeshIOHandler {
    public:
        /**
         * \brief Loads a mesh from a file in OFF format.
         * \param[in] filename name of the file
         * \param[out] M the loaded mesh
         * \param[in] ioflags specifies which attributes and elements should be read
         * \tparam MESH the mesh class
         * \return true on success, false otherwise
         */
        template <class MESH>
        bool load_generic(
            const std::string& filename, MESH& M, const MeshIOFlags& ioflags
        ) {
            geo_argused(ioflags);

            // Note: Vertices indexes start by 0 in off format.

            LineInput in(filename);
            if(!in.OK()) {
                return false;
            }

            if(!in.get_line()) {
                Logger::err("I/O")
                    << "Unexpected end of file"
                    << std::endl;
                return false;
            }
            in.get_fields();
            if(in.nb_fields() == 0 || !in.field_matches(0, "OFF")) {
                Logger::err("I/O")
                    << "Line " << in.line_number()
                    << ": unrecognized header"
                    << std::endl;
                return false;
            }

            if(!in.get_line()) {
                Logger::err("I/O")
                    << "Line " << in.line_number()
                    << ": unexpected end of file"
                    << std::endl;
                return false;
            }
            in.get_fields();
            if(in.nb_fields() != 3) {
                Logger::err("I/O")
                    << "Line " << in.line_number()
                    << ": unrecognized header"
                    << std::endl;
                return false;
            }

            index_t nb_vertices = in.field_as_uint(0);
            // second field is nb edges (unused)
            index_t nb_facets = in.field_as_uint(2);

            typedef typename MESH::coord_t coord_t;
            GEOGen::MeshBuilder<MESH> builder(&M);
            builder.begin_mesh();
            builder.create_vertices(nb_vertices, 3);
            builder.reserve_facets(nb_facets);

            for(index_t i = 0; i < nb_vertices; i++) {
                if(!in.get_line()) {
                    Logger::err("I/O")
                        << "Line " << in.line_number()
                        << ": unexpected end of file"
                        << std::endl;
                    return false;
                }
                in.get_fields();
                if(in.nb_fields() != 3) {
                    Logger::err("I/O")
                        << "Line " << in.line_number()
                        << ": invalid number of fields: expected 3 coordinates, got "
                        << in.nb_fields()
                        << std::endl;
                    return false;
                }
                double x = in.field_as_double(0);
                double y = in.field_as_double(1);
                double z = in.field_as_double(2);

                coord_t* p = builder.target()->vertex_ptr(i);
                p[0] = coord_t(x);
                p[1] = coord_t(y);
                p[2] = coord_t(z);
            }

            while(!in.eof() && in.get_line()) {
                in.get_fields();
                if(in.nb_fields() < 4) {
                    Logger::err("I/O")
                        << "Line " << in.line_number()
                        << ": facet line only has " << in.nb_fields()
                        << " fields (expected 1 count + at least 3 corner fields)"
                        << std::endl;
                    return false;
                }
                index_t nb_facet_vertices = in.field_as_uint(0);
                if(in.nb_fields() != nb_facet_vertices + 1) {
                    Logger::err("I/O")
                        << "Line " << in.line_number()
                        << ": facet has " << in.nb_fields() - 1
                        << " actual vertices ("
                        << nb_facet_vertices << " expected)"
                        << std::endl;
                    return false;
                }
                builder.begin_facet();
                for(index_t j = 0; j < nb_facet_vertices; j++) {
                    index_t vertex_index = in.field_as_uint(j + 1);
                    if(vertex_index >= builder.target()->nb_vertices()) {
                        Logger::err("I/O")
                            << "Line " << in.line_number()
                            << ": facet corner #" << j
                            << " references an invalid vertex: "
                            << vertex_index
                            << std::endl;
                        return false;
                    }
                    builder.add_vertex_to_facet(vertex_index);
                }
                builder.end_facet();
            }
            builder.end_mesh(false);
            return true;
        }

        /**
         * \brief Saves a mesh into a file in OFF format.
         * \param[in] M The mesh to save
         * \param[in] filename name of the file
         * \param[in] ioflags specifies which attributes and elements should be saved
         * \return true on success, false otherwise
         */
        virtual bool save(
            const Mesh& M, const std::string& filename, const MeshIOFlags& ioflags
        ) {
            geo_argused(ioflags);

            std::ofstream output(filename.c_str());
            if(!output) {
                return false;
            }
            output << "OFF" << std::endl;

            output << M.nb_vertices() << " "
                << M.nb_facets() << " "
                << M.nb_corners() / 2
                << std::endl;

            // Output Vertices
            for(index_t v = 0; v < M.nb_vertices(); v++) {
                const double* p = M.vertex_ptr(v);
                output << p[0] << " " << p[1] << " " << p[2] << std::endl;
            }

            // Output facets
            for(index_t f = 0; f < M.nb_facets(); f++) {
                output << M.facet_size(f) << " ";
                for(index_t c = M.facet_begin(f); c < M.facet_end(f); c++) {
                    output << M.corner_vertex_index(c) << " ";
                }
                output << std::endl;
            }
            return true;
        }

        virtual bool load(
            const std::string& filename, Mesh& M,
            const MeshIOFlags& ioflags = MeshIOFlags()
        ) {
            return load_generic(filename, M, ioflags);
        }

        virtual bool load(
            const std::string& filename, SinglePrecisionMesh& M,
            const MeshIOFlags& ioflags = MeshIOFlags()
        ) {
            return load_generic(filename, M, ioflags);
        }

    protected:
        virtual ~OFFIOHandler() {
        }
    };

    /************************************************************************/

    /**
     * \brief IO handler for the STL file format (ascii and binary)
     * \see http://en.wikipedia.org/wiki/STL_(file_format)
     */
    class GEOGRAM_API STLIOHandler : public MeshIOHandler {
    public:
        /**
         * \brief Loads a mesh from a file in STL format (ascii version).
         * \param[in] filename name of the file
         * \param[out] M the loaded mesh
         * \param[in] ioflags specifies which attributes and elements should be read
         * \tparam MESH the mesh class
         * \return true on success, false otherwise
         */
        template <class MESH>
        bool load_ascii(
            const std::string& filename, MESH& M, const MeshIOFlags& ioflags
        ) {
            geo_argused(ioflags);

            LineInput in(filename);
            if(!in.OK()) {
                return false;
            }

            int current_chart = -1;
            index_t current = 0;
            bool facet_opened = false;
            GEOGen::MeshBuilder<MESH> builder(&M);
            builder.begin_mesh();
            while(!in.eof() && in.get_line()) {
                in.get_fields();
                if(in.field_matches(0, "outer")) {
                    builder.begin_facet();
                    facet_opened = true;
                } else if(in.field_matches(0, "endloop")) {
                    builder.end_facet(geo_max(current_chart, 0));
                    facet_opened = false;
                    // TODO: store current_chart in attribute
                } else if(in.field_matches(0, "vertex")) {
                    if(in.nb_fields() < 4) {
                        Logger::err("I/O")
                            << "Line " << in.line_number()
                            << ": vertex line has " << in.nb_fields() - 1
                            << " fields (at least 3 required)"
                            << std::endl;
                        return false;
                    }
                    vec3 p;
                    p.x = in.field_as_double(1);
                    p.y = in.field_as_double(2);
                    p.z = in.field_as_double(3);
                    builder.add_vertex(p);
                    builder.add_vertex_to_facet(current);
                    current++;
                } else if(in.field_matches(0, "solid")) {
                    // TODO: create attribute if needed
                    current_chart++;
                }
            }

            if(facet_opened) {
                Logger::err("I/O")
                    << "Line " << in.line_number()
                    << ": current facet is not closed"
                    << std::endl;
                return false;
            }

            builder.end_mesh(false);
            return true;
        }

        /**
         * \brief Loads a mesh from a file in STL format (binary version).
         * \param[in] filename name of the file
         * \param[out] M the loaded mesh
         * \param[in] ioflags specifies which attributes and elements should be read
         * \tparam MESH the mesh class
         * \return true on success, false otherwise
         */
        template <class MESH>
        bool load_binary(
            const std::string& filename, MESH& M, const MeshIOFlags& ioflags
        ) {
            geo_argused(ioflags);

            GEOGen::MeshBuilder<MESH> builder(&M);
            BinaryInputStream in(filename, BinaryStream::GEO_LITTLE_ENDIAN);
            char header[80];
            in.read_opaque_data(header, 80);
            if(!in.OK()) {
                throw "failed to read header";
            }
            Numeric::uint32 nb_triangles;
            in >> nb_triangles;
            if(!in.OK()) {
                throw "failed to read number of triangles";
            }
            index_t cur = 0;
            builder.begin_mesh();
            for(index_t t = 0; t < nb_triangles; t++) {
                Numeric::float32 N[3];
                Numeric::float32 XYZ[9];
                in >> N[0] >> N[1] >> N[2];
                for(index_t i = 0; i < 9; i++) {
                    in >> XYZ[i];
                }
                if(!in.OK()) {
                    throw "failed to read triangle";
                }
                Numeric::uint16 attrib;
                in >> attrib;
                builder.add_vertex(vec3(XYZ[0], XYZ[1], XYZ[2]));
                builder.add_vertex(vec3(XYZ[3], XYZ[4], XYZ[5]));
                builder.add_vertex(vec3(XYZ[6], XYZ[7], XYZ[8]));
                builder.begin_facet();
                builder.add_vertex_to_facet(cur);
                builder.add_vertex_to_facet(cur + 1);
                builder.add_vertex_to_facet(cur + 2);
                builder.end_facet(attrib);
                cur += 3;
            }
            builder.end_mesh(false);
            return true;
        }

        /**
         * \brief Loads a mesh from a file in STL format.
         * \details Supports both ascii and binary STL.
         * \param[in] filename name of the file
         * \param[out] M the loaded mesh
         * \param[in] ioflags specifies which attributes and elements should be read
         * \tparam MESH the mesh class
         * \return true on success, false otherwise
         */
        template <class MESH>
        bool load_generic(
            const std::string& filename, MESH& M, const MeshIOFlags& ioflags
        ) {
            FILE* F = fopen(filename.c_str(), "rb");
            if(F == nil) {
                return false;
            }

            // The safe way of checking whether an STL file is
            // binary is to check whether the size of the file
            // matches the size deduced from the number of triangles
            // (many binary STL files start with SOLID although it
            //  is supposed to be only for ASCII STL files)
            fseek(F, 80, SEEK_SET);
            Numeric::uint32 nb_triangles;
            if(fread(&nb_triangles, sizeof(nb_triangles), 1, F) != 1) {
                Logger::err("I/O")
                    << "Cannot deduce the format of STL file"
                    << std::endl;
                fclose(F);
                return false;
            }
            fseek(F, 0, SEEK_END);
            long file_size = ftell(F);
            fclose(F);
            bool result;
            if(file_size == long(nb_triangles * 50 + 84)) {
                result = load_binary(filename, M, ioflags);
            } else {
                result = load_ascii(filename, M, ioflags);
            }
            return result;
        }

        /**
         * \brief Writes a point to a binary stream.
         * \param[in] out the binary stream
         * \param[in] V the vector to write
         */
        inline void write_stl_vector(BinaryOutputStream& out, const vec3& V) {
            out << Numeric::float32(V.x);
            out << Numeric::float32(V.y);
            out << Numeric::float32(V.z);
        }

        /**
         * \brief Saves a mesh into a file in STL binary format.
         * \param[in] M The mesh to save
         * \param[in] filename name of the file
         * \param[in] ioflags specifies which attributes and elements should be saved
         * \return true on success, false otherwise
         */
        virtual bool save(
            const Mesh& M, const std::string& filename, const MeshIOFlags& ioflags
        ) {
            geo_argused(ioflags);

            BinaryOutputStream out(filename, BinaryStream::GEO_LITTLE_ENDIAN);
            char header[80];
            Memory::clear(header, 80);
            strcpy(header, "generated with \\V(O|R|P/A|L|I|N|E");
            out.write_opaque_data(header, 80);
            Numeric::uint32 nb_triangles = 0;
            for(index_t f = 0; f < M.nb_facets(); f++) {
                nb_triangles += (M.facet_size(f) - 2);
            }
            out << nb_triangles;
            for(index_t f = 0; f < M.nb_facets(); f++) {
                vec3 N = Geom::mesh_facet_normal(M, f);
                index_t c1 = M.facet_begin(f);
                const vec3& p1 = Geom::mesh_vertex(M, M.corner_vertex_index(c1));
                for(index_t c2 = M.facet_begin(f) + 1;
                    c2 + 1 < M.facet_end(f); c2++
                ) {
                    const vec3& p2 = Geom::mesh_vertex(
                        M, M.corner_vertex_index(c2)
                    );
                    const vec3& p3 = Geom::mesh_vertex(
                        M, M.corner_vertex_index(c2 + 1)
                    );
                    Numeric::uint16 attribute = Numeric::uint16(
                        M.has_attribute(MESH_FACET_REGION) ? M.facet_region(f) : 0
                    );
                    write_stl_vector(out, N);
                    write_stl_vector(out, p1);
                    write_stl_vector(out, p2);
                    write_stl_vector(out, p3);
                    out << attribute;
                }
            }
            return true;
        }

        virtual bool load(
            const std::string& filename, Mesh& M,
            const MeshIOFlags& ioflags = MeshIOFlags()
        ) {
            return load_generic(filename, M, ioflags);
        }

        virtual bool load(
            const std::string& filename, SinglePrecisionMesh& M,
            const MeshIOFlags& ioflags = MeshIOFlags()
        ) {
            return load_generic(filename, M, ioflags);
        }

    protected:
        virtual ~STLIOHandler() {
        }
    };

    /************************************************************************/

    /**
     * \brief IO handler for the XYZ file format
     * \details Currtently only loading is supported
     */
    class GEOGRAM_API XYZIOHandler : public MeshIOHandler {
    public:
        /**
         * \brief Loads a pointset from a file in XYZ format.
         * \param[in] filename name of the file
         * \param[out] M the mesh where to store the points
         * \param[in] ioflags specifies which attributes and elements should be read
         * \tparam MESH the mesh class
         * \return true on success, false otherwise
         */
        template <class MESH>
        bool load_generic(
            const std::string& filename, MESH& M, const MeshIOFlags& ioflags
        ) {
            geo_argused(ioflags);
            typedef typename MESH::coord_t coord_t;

            M.clear();
            LineInput in(filename);
            if(!in.OK()) {
                return false;
            }
            index_t nb_points = 0;
            coord_index_t dimension = 0;
            vector<coord_t>& points =
                GEOGen::MeshMutator<MESH>::vertices(M);
            while(!in.eof() && in.get_line()) {
                in.get_fields();
                switch(in.nb_fields()) {
                    case 1:
                        nb_points = in.field_as_uint(0);
                        break;
                    case 3:
                    case 4:
                    case 6:
                    {
                        if(dimension == 0) {
                            // dimension = in.nb_fields() ;
                            dimension = 3;
                            M.set_dimension(dimension);
                            if(nb_points != 0) {
                                points.reserve(M.dimension() * nb_points);
                            }
                        }
                        for(coord_index_t c = 0; c < M.dimension(); c++) {
                            double x =
                                (c < in.nb_fields()) ? in.field_as_double(c) : 0.0;
                            points.push_back(coord_t(x));
                        }
                    }
                    break;
                    default:
                        Logger::err("I/O")
                            << "Line " << in.line_number()
                            << ": wrong number of fields"
                            << std::endl;
                        return false;
                        break;
                }
            }
            M.update_cached_variables();
            return true;
        }

        virtual bool save(
            const Mesh& M, const std::string& filename, const MeshIOFlags& ioflags
        ) {
            geo_argused(M);
            geo_argused(filename);
            geo_argused(ioflags);

            Logger::err("I/O")
                << "Saving a mesh in XYZ format is not supported"
                << std::endl;

            return false;
        }

        virtual bool load(
            const std::string& filename, Mesh& M,
            const MeshIOFlags& ioflags = MeshIOFlags()
        ) {
            return load_generic(filename, M, ioflags);
        }

        virtual bool load(
            const std::string& filename, SinglePrecisionMesh& M,
            const MeshIOFlags& ioflags = MeshIOFlags()
        ) {
            return load_generic(filename, M, ioflags);
        }

    protected:
        virtual ~XYZIOHandler() {
        }
    };

    /************************************************************************/

    /**
     * \brief IO handler for the PTS file format
     * \details Currtently only loading is supported
     */
    class GEOGRAM_API PTSIOHandler : public MeshIOHandler {
    public:
        /**
         * \brief Loads a pointset from a file in PTS format.
         * \param[in] filename name of the file
         * \param[out] M the mesh where to store the points
         * \param[in] ioflags specifies which attributes and elements should be read
         * \tparam MESH the mesh class
         * \return true on success, false otherwise
         */
        template <class MESH>
        bool load_generic(
            const std::string& filename, MESH& M, const MeshIOFlags& ioflags
        ) {
            geo_argused(ioflags);
            typedef typename MESH::coord_t coord_t;

            M.clear();
            LineInput in(filename);
            if(!in.OK()) {
                return false;
            }
            M.set_dimension(3);
            vector<coord_t>& points =
                GEOGen::MeshMutator<MESH>::vertices(M);
            while(!in.eof() && in.get_line()) {
                in.get_fields();
                if(in.nb_fields() == 4 && in.field_matches(0,"v")) {
                    points.push_back(coord_t(in.field_as_double(1)));
                    points.push_back(coord_t(in.field_as_double(2)));
                    points.push_back(coord_t(in.field_as_double(3)));
                } else {
                    Logger::err("I/O")
                        << "Line " << in.line_number()
                        << ": wrong number of fields"
                        << std::endl;
                    return false;
                }
            }
            M.update_cached_variables();
            return true;
        }

        virtual bool save(
            const Mesh& M, const std::string& filename, const MeshIOFlags& ioflags
        ) {
            geo_argused(M);
            geo_argused(filename);
            geo_argused(ioflags);

            Logger::err("I/O")
                << "Saving a mesh in XYZ format is not supported"
                << std::endl;

            return false;
        }

        virtual bool load(
            const std::string& filename, Mesh& M,
            const MeshIOFlags& ioflags = MeshIOFlags()
        ) {
            return load_generic(filename, M, ioflags);
        }

        virtual bool load(
            const std::string& filename, SinglePrecisionMesh& M,
            const MeshIOFlags& ioflags = MeshIOFlags()
        ) {
            return load_generic(filename, M, ioflags);
        }

    protected:
        virtual ~PTSIOHandler() {
        }
    };

    /************************************************************************/


    /**
     * \brief IO handler for the TET file format
     */
    class GEOGRAM_API TETIOHandler : public MeshIOHandler {
    public:
        /**
         * \brief Creates a TET IO handler.
         * \param[in] dim dimension of the vertices (3 for regular 3d mesh)
         */
        TETIOHandler(coord_index_t dimension = 3) :
            dimension_(dimension) {
        }

        /**
         * \brief Loads a mesh from a file in TET format.
         * \details Only tetrahedral cells are supported for now.
         * \param[in] filename name of the file
         * \param[out] M the loaded mesh
         * \param[in] ioflags specifies which attributes and elements 
         *  should be read
         * \tparam MESH the mesh class
         * \return true on success, false otherwise
         */
        template <class MESH>
        bool load_generic(
            const std::string& filename, MESH& M, const MeshIOFlags& ioflags
        ) {
            M.clear();
            LineInput in(filename);
            if(!in.OK()) {
                return false;
            }
            if(!in.get_line()) {
                Logger::err("I/O")
                    << "Unexpected end of file"
                    << std::endl;
                return false;
            }
            in.get_fields();

            index_t nb_vertices = 0;
            index_t nb_cells = 0;
            bool has_arbitrary_cells = false;
            
            if(
                in.nb_fields() == 4 &&
                in.field_matches(1, "vertices") &&
                (
                    in.field_matches(3, "tets") ||
                    in.field_matches(3, "cells")
                )
            ) {
                nb_vertices = in.field_as_uint(0);
                nb_cells = in.field_as_uint(2);
                has_arbitrary_cells = in.field_matches(3, "cells");
            } else {
                if(in.nb_fields() != 2 || !in.field_matches(1, "vertices")) {
                    Logger::err("I/O")
                        << "Line " << in.line_number()
                        << " expected <number_of_vertices> vertices"
                        << std::endl;
                    return false;
                }
                nb_vertices = in.field_as_uint(0);
            
                if(!in.get_line()) {
                    Logger::err("I/O")
                        << "Unexpected end of file"
                        << std::endl;
                    return false;
                }
                in.get_fields();
                if(
                    in.nb_fields() != 2 || (
                        !in.field_matches(1, "tets") &&
                        !in.field_matches(1, "cells")
                    )
                ) {
                    Logger::err("I/O")
                        << "Line " << in.line_number()
                        << " expected <number_of_tets> tets"
                        << std::endl;
                    return false;
                }
                nb_cells = in.field_as_uint(0);
                has_arbitrary_cells = in.field_matches(1, "cells");                
            } 

            GEO::vector<typename MESH::coord_t> vertices;
            vertices.reserve(nb_vertices * dimension_);
            for(index_t v = 0; v < nb_vertices; ++v) {
                if(!in.get_line()) {
                    Logger::err("I/O")
                        << "Unexpected end of file"
                        << std::endl;
                    return false;
                }
                in.get_fields();
                if(in.nb_fields() != index_t(dimension_)) {
                    Logger::err("I/O")
                        << "Line " << in.line_number()
                        << " expected " << dimension_ << " point coordinates"
                        << std::endl;
                }
                for(coord_index_t c = 0; c < dimension_; ++c) {
                    vertices.push_back(
                        (typename MESH::coord_t)(in.field_as_double(c))
                    );
                }
            }
            GEOGen::MeshMutator<MESH>::set_dimension(M, dimension_);
            GEOGen::MeshMutator<MESH>::vertices(M).swap(vertices);
            GEOGen::MeshMutator<MESH>::set_nb_vertices(M, nb_vertices);
            M.update_cached_variables();

            if(ioflags.has_element(MESH_CELLS)) {
                if(has_arbitrary_cells) {
                    for(index_t t = 0; t < nb_cells; ++t) {
                        if(!in.get_line()) {
                            Logger::err("I/O")
                                << "Unexpected end of file"
                                << std::endl;
                            return false;
                        }
                        in.get_fields();
                        if(in.nb_fields() >= 2 && in.field_matches(0,"#") && in.field_matches(1,"C")) {
                            if(in.nb_fields() != 6) {
                                Logger::err("I/O")                                
                                    << "Line " << in.line_number()
                                    << " expected # C v1 v2 v3 v4"
                                    << std::endl;
                                return false;
                            }
                            M.add_connector(
                                in.field_as_uint(2), in.field_as_uint(3), in.field_as_uint(4), in.field_as_uint(5)
                            );
                        } else {
                            if(in.nb_fields() > 0) {
                                index_t nb_vertices_in_cell = in.field_as_uint(0);
                                switch(nb_vertices_in_cell) {
                                case 4: {
                                    M.add_tet(
                                        in.field_as_uint(1), in.field_as_uint(2),
                                        in.field_as_uint(3), in.field_as_uint(4)
                                    );                                    
                                } break;
                                case 8: {
                                    M.add_hex(
                                        in.field_as_uint(1), in.field_as_uint(2),
                                        in.field_as_uint(3), in.field_as_uint(4),
                                        in.field_as_uint(5), in.field_as_uint(6),
                                        in.field_as_uint(7), in.field_as_uint(8)
                                    );
                                } break;
                                case 6: {
                                    M.add_prism(
                                        in.field_as_uint(1), in.field_as_uint(2), in.field_as_uint(3),
                                        in.field_as_uint(4), in.field_as_uint(5), in.field_as_uint(6)
                                    );
                                } break;
                                case 5: {
                                    M.add_pyramid(
                                        in.field_as_uint(1), in.field_as_uint(2), in.field_as_uint(3),
                                        in.field_as_uint(4), in.field_as_uint(5)
                                    );
                                } break;
                                default: {
                                Logger::err("I/O")                                
                                    << "Line " << in.line_number()
                                    << " unexpected number of vertices in cell:"
                                    << nb_vertices_in_cell
                                    << std::endl;
                                return false;
                                }
                                }
                            }
                        }
                    }
                    M.connect_cells(); 
                } else {
                    GEO::vector<index_t> tets;
                    tets.reserve(4 * nb_cells);
                    for(index_t t = 0; t < nb_cells; ++t) {
                        if(!in.get_line()) {
                            Logger::err("I/O")
                                << "Unexpected end of file"
                                << std::endl;
                            return false;
                        }
                        in.get_fields();
                        if(in.nb_fields() != 5 || in.field_as_int(0) != 4) {
                            Logger::err("I/O")
                                << "Line " << in.line_number()
                                << " expected 4 v1 v2 v3 v4"
                                << std::endl;
                        } else {
                            for(index_t i = 0; i < 4; ++i) {
                                index_t v = in.field_as_uint(i + 1);
                                if(i >= nb_vertices) {
                                    Logger::err("I/O")
                                        << "Line " << in.line_number()
                                        << "invalid vertex index"
                                        << std::endl;
                                    return false;
                                }
                                tets.push_back(v);
                            }
                        }
                    }
                    GEOGen::MeshMutator<MESH>::tet_adjacents(M).assign(
                        tets.size(), -1
                    );
                    GEOGen::MeshMutator<MESH>::set_nb_cells(M,tets.size()/4);                    
                    GEOGen::MeshMutator<MESH>::tet_vertices(M).swap(tets);
                    M.connect_tets();
                }
            }
            return true;
        }

        /**
         * \brief Saves a mesh into a file in TET format.
         * \param[in] M The mesh to save
         * \param[in] filename name of the file
         * \param[in] ioflags specifies which attributes and elements 
         * should be saved
         * \return true on success, false otherwise
         */
        virtual bool save(
            const Mesh& M, const std::string& filename, 
            const MeshIOFlags& ioflags
        ) {
            if(M.dimension() < dimension_) {
                return false;
            }
            if(!ioflags.has_element(MESH_CELLS)) {
                return false;
            }
            std::ofstream out(filename.c_str());
            if(!out) {
                return false;
            }
            if(M.is_tetrahedralized()) {
                out << M.nb_vertices() << " vertices" << std::endl;
                out << M.nb_tets() << " tets" << std::endl;
                for(index_t v = 0; v < M.nb_vertices(); ++v) {
                    for(coord_index_t c = 0; c < dimension_; ++c) {
                        out << M.vertex_ptr(v)[c] << " ";
                    }
                    out << std::endl;
                }
                for(index_t t = 0; t < M.nb_tets(); ++t) {
                    out << "4";
                    for(index_t i = 0; i < 4; ++i) {
                        out << " " << M.tet_vertex_index(t, i);
                    }
                    out << std::endl;
                }
            } else {
                out << M.nb_vertices() << " vertices" << std::endl;
                out << M.nb_cells() << " cells" << std::endl;
                for(index_t v = 0; v < M.nb_vertices(); ++v) {
                    for(coord_index_t c = 0; c < dimension_; ++c) {
                        out << M.vertex_ptr(v)[c] << " ";
                    }
                    out << std::endl;
                }
                bool has_connectors = false;
                for(index_t c=0; c<M.nb_cells(); ++c) {
                    switch(M.cell_type(c)) {
                    case MESH_TET: 
                    case MESH_HEX: 
                    case MESH_PRISM: 
                    case MESH_PYRAMID: {
                        out << M.cell_nb_vertices(c) << " ";
                        for(index_t lv=0; lv<M.cell_nb_vertices(c); ++lv) {
                            out << M.cell_vertex_index(c,lv) << " ";
                        }
                        out << std::endl;
                    } break;
                    case MESH_CONNECTOR: {
                        has_connectors=true;
                    } break;
                    default: {
                        geo_assert_not_reached;
                    } break;
                    }
                }
                if(has_connectors) {
                    for(index_t c=0; c<M.nb_cells(); ++c) {
                        if(M.cell_type(c) == MESH_CONNECTOR) {
                            out << "# C"
                                << " " << M.cell_vertex_index(c,0)
                                << " " << M.cell_vertex_index(c,1)
                                << " " << M.cell_vertex_index(c,2)
                                << " " << M.cell_vertex_index(c,3)
                                << std::endl;
                        }
                    }
                }
            }
            return true;
        }

        virtual bool load(
            const std::string& filename, Mesh& M,
            const MeshIOFlags& ioflags = MeshIOFlags()
        ) {
            return load_generic(filename, M, ioflags);
        }

        virtual bool load(
            const std::string& filename, SinglePrecisionMesh& M,
            const MeshIOFlags& ioflags = MeshIOFlags()
        ) {
            return load_generic(filename, M, ioflags);
        }

    protected:
        virtual ~TETIOHandler() {
        }

    private:
        coord_index_t dimension_;
    };

    /************************************************************************/

    /**
     * \brief IO handler for the TET6 file format
     * \see TETIOHandler
     */
    class GEOGRAM_API TET6IOHandler : public TETIOHandler {
    public:
        TET6IOHandler() :
            TETIOHandler(6) {
        }

    protected:
        virtual ~TET6IOHandler() {
        }
    };

    /************************************************************************/

    /**
     * \brief Loads a mesh from a file.
     * \details File format is determined from extension.
     * \param[in] filename name of the file
     * \param[out] M the loaded mesh
     * \param[in] ioflags specifies which attributes and elements should be read
     * \tparam MESH the mesh class
     * \return true on success, false otherwise
     */
    template <class MESH>
    bool mesh_load_generic(
        const std::string& filename, MESH& M, const MeshIOFlags& ioflags
    ) {
        Logger::out("I/O")
            << "Loading input file " << filename << "..."
            << std::endl;

        M.clear();

        // TODO: maybe remove this line, and do that on a case-by-case basis
        GEOGen::MeshMutator<MESH>::set_attributes(
            M, ioflags.attributes()
        );

        bool result = false;
        MeshIOHandler_var handler = MeshIOHandler::get_handler(filename);
        if(handler != nil) {
            try {
                result = handler->load(filename, M, ioflags);
            }
            catch(const std::exception& ex) {
                Logger::err("I/O") << ex.what() << std::endl;
                result = false;
            }
        }

        if(!result) {
            Logger::err("I/O")
                << "Could not load file: " << filename
                << std::endl;
            return false;
        }

        M.show_stats("I/O");

        if(!M.is_valid()) {
            Logger::warn("I/O") << "Loaded mesh is invalid" << std::endl;
        }

        if(M.nb_vertices() != 0) {
            typedef typename MESH::coord_t coord_t;
            index_t nb = M.nb_vertices() * M.dimension();
            coord_t* p = M.vertex_ptr(0);
            bool has_nan = false;
            for(index_t i = 0; i < nb; i++) {
                if(Numeric::is_nan(*p)) {
                    has_nan = true;
                    *p = 0.0;
                }
                p++;
            }
            if(has_nan) {
                Logger::warn("I/O") << "Found NaNs in input file" << std::endl;
            }
        }

        if(M.nb_facets() == 0 && M.is_tetrahedralized() && M.nb_tets() != 0) {
            M.compute_tets_boundaries();
        }

        return true;
    }
}

/****************************************************************************/

namespace GEO {

    MeshIOFlags::MeshIOFlags() {
        dimension_ = 3;
        attributes_ = MeshAttributes(MESH_NO_ATTRIBUTES);
        elements_ = MeshElements(MESH_VERTICES | MESH_FACETS);
    }

    /************************************************************************/

    bool mesh_load(
        const std::string& filename, Mesh& M, const MeshIOFlags& ioflags
    ) {
        return mesh_load_generic(filename, M, ioflags);
    }

    bool mesh_load(
        const std::string& filename, SinglePrecisionMesh& M,
        const MeshIOFlags& ioflags
    ) {
        return mesh_load_generic(filename, M, ioflags);
    }

    bool mesh_save(
        const Mesh& M, const std::string& filename, const MeshIOFlags& ioflags
    ) {
        Logger::out("I/O")
            << "Saving file " << filename << "..."
            << std::endl;

        MeshIOHandler_var handler = MeshIOHandler::get_handler(filename);
        if(handler != nil && handler->save(M, filename, ioflags)) {
            return true;
        }

        Logger::err("I/O")
            << "Could not save file: " << filename
            << std::endl;
        return false;
    }

    /************************************************************************/

    MeshIOHandler* MeshIOHandler::create(const std::string& format) {
        geo_register_MeshIOHandler_creator(LMIOHandler, "mesh");
        geo_register_MeshIOHandler_creator(LMIOHandler, "meshb");
        geo_register_MeshIOHandler_creator(OBJIOHandler, "obj");
        geo_register_MeshIOHandler_creator(OBJIOHandler, "eobj");
        geo_register_MeshIOHandler_creator(OBJ6IOHandler, "obj6");
        geo_register_MeshIOHandler_creator(PLYIOHandler, "ply");
        geo_register_MeshIOHandler_creator(OFFIOHandler, "off");
        geo_register_MeshIOHandler_creator(STLIOHandler, "stl");
        geo_register_MeshIOHandler_creator(XYZIOHandler, "xyz");
        geo_register_MeshIOHandler_creator(PTSIOHandler, "pts");
        geo_register_MeshIOHandler_creator(TETIOHandler, "tet");
        geo_register_MeshIOHandler_creator(TET6IOHandler, "tet6");

        MeshIOHandler* handler = MeshIOHandlerFactory::create_object(format);
        if(handler != nil) {
            return handler;
        }

        Logger::err("I/O")
            << "Unsupported file format: " << format
            << std::endl;
        return nil;
    }

    MeshIOHandler* MeshIOHandler::get_handler(const std::string& filename) {
        std::string ext = FileSystem::extension(filename);
        return create(ext);
    }

    MeshIOHandler::~MeshIOHandler() {
    }
}

