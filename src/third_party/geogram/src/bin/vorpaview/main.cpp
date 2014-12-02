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

#include <geogram_gfx/third_party/glew/glew.h>
#include <geogram_gfx/glut_viewer/glut_viewer.h>
#include <geogram/mesh/mesh.h>
#include <geogram/mesh/mesh_io.h>
#include <geogram/mesh/mesh_geometry.h>
#include <geogram/mesh/mesh_private.h>
#include <geogram/basic/command_line.h>
#include <geogram/basic/command_line_args.h>

namespace {

    GEO::SinglePrecisionMesh M;

    // Vertex buffer objects
    GLuint vertices_VBO = 0;
    GLuint facet_indices_VBO = 0;
    GLuint tet_indices_VBO = 0;

    GLboolean show_mesh = GL_FALSE;
    GLboolean show_colors = GL_TRUE;
    GLboolean show_borders = GL_FALSE;
    GLboolean white_bg = GL_TRUE;
    GLboolean show_surface = GL_TRUE;
    GLboolean show_volume = GL_FALSE;
    GLboolean show_vertices = GL_FALSE;
    GLboolean lighting = GL_TRUE;

    GLuint program_points = 0;
    GLuint program_tri = 0;
    GLuint program_tet = 0;

    /************************************************************************/
    //
    // Shaders for automatic generation of per-facet normals and
    //  fragment-level mesh rendering.

    // Note: GLSL version 150 is the lowest version that supports
    //  geometry shaders (we use it for largest compatibility).
    // Latest GLSL version is 430 (to be used if we use SSBOs to implement
    //  irregular mesh rendering)

    /**
     * \brief The vertex shader.
     * \details It does nothing (pass-through), all the
     *   work is done by the geometry shader.
     */
    const char* vshader_source =
        "#version 150 compatibility                                         \n"
        " void main(void) {                                                 \n"
        "    gl_Position = gl_Vertex ;                                      \n"
        " }                                                                 \n";

    /**
     * \brief The geometry shader for triangles.
     * \details Fetches the triangles from the mesh,
     * does the clipping, transform and lighting
     * computations, and initializes interpolated barycentric
     * coordinates (used by the fragment shader to display
     * the mesh).
     */
    const char* gshader_tri_source =
        "#version 150 compatibility                                         \n"
        "layout(triangles) in;                                              \n"
        "layout(triangle_strip, max_vertices = 3) out;                      \n"
        "flat out float diffuse;                                            \n"
        "flat out float specular;                                           \n"
        "uniform bool lighting=true;                                        \n"
        "out vec2 bary;                                                     \n"
        "vec4 project(vec4 V) {                                             \n"
        "   return gl_ModelViewProjectionMatrix * V ;                       \n"
        "}                                                                  \n"
        "float clip(vec4 V) {                                               \n"
        "   return dot(gl_ModelViewMatrix*V,gl_ClipPlane[0]);               \n"
        "}                                                                  \n"
        "float cosangle(vec3 N, vec3 L) {                                   \n"
        "   float s = inversesqrt(dot(N,N)*dot(L,L)) ;                      \n"
        "   return s*dot(N,L) ;                                             \n"
        "}                                                                  \n"
        "void main() {                                                      \n"
        " if(lighting) {                                                    \n"
        "    vec3 p1 = gl_in[0].gl_Position.xyz ;                           \n"
        "    vec3 p2 = gl_in[1].gl_Position.xyz ;                           \n"
        "    vec3 p3 = gl_in[2].gl_Position.xyz ;                           \n"
        "    vec3 N = gl_NormalMatrix * cross(p2-p1,p3-p1) ;                \n"
        "    vec3 L = gl_LightSource[0].position.xyz ;                      \n"
        "    diffuse = cosangle(N,L) ;                                      \n"
        "    float NdotH = cosangle(N,gl_LightSource[0].halfVector.xyz) ;   \n"
        "    specular = pow(abs(NdotH),gl_FrontMaterial.shininess);         \n"
        " } else {                                                          \n"
        "    diffuse = -1.0 ; specular = 0.0 ;                              \n"
        " }                                                                 \n"
        " gl_Position=project(gl_in[0].gl_Position);                        \n"
        "   gl_ClipDistance[0] = clip(gl_in[0].gl_Position);                \n"
        "   bary = vec2(0.0,0.0) ; EmitVertex();                            \n"
        " gl_Position=project(gl_in[1].gl_Position);                        \n"
        "   gl_ClipDistance[0] = clip(gl_in[1].gl_Position);                \n"
        "   bary = vec2(1.0,0.0) ; EmitVertex();                            \n"
        " gl_Position=project(gl_in[2].gl_Position);                        \n"
        "   gl_ClipDistance[0] = clip(gl_in[2].gl_Position);                \n"
        "   bary = vec2(0.0,1.0) ; EmitVertex();                            \n"
        " EndPrimitive();                                                   \n"
        "}                                                                  \n";

    /**
     * \brief The geometry shader for tetrahedra.
     * \details Fetches the tetrahedra from the mesh,
     * does the clipping, transform and lighting
     * computations, and initializes interpolated barycentric
     * coordinates (used by the fragment shader to display
     * the mesh).
     */
    const char* gshader_tet_source =
        "#version 150 compatibility                                         \n"
        "layout(lines_adjacency) in;                                        \n"
        "layout(triangle_strip, max_vertices = 12) out;                     \n"
        "flat out float diffuse;                                            \n"
        "flat out float specular;                                           \n"
        "uniform bool lighting=true;                                        \n"
        "uniform bool clipping=true;                                        \n"
        "out vec2 bary;                                                     \n"
        "vec4 project(vec4 V) {                                             \n"
        "   return gl_ModelViewProjectionMatrix * V ;                       \n"
        "}                                                                  \n"
        "float clip(vec4 V) {                                               \n"
        "   return dot(gl_ModelViewMatrix*V,gl_ClipPlane[0]);               \n"
        "}                                                                  \n"
        "float cosangle(vec3 N, vec3 L) {                                   \n"
        "   float s = inversesqrt(dot(N,N)*dot(L,L)) ;                      \n"
        "   return s*dot(N,L) ;                                             \n"
        "}                                                                  \n"
        "void flat_shaded_triangle(                                         \n"
        "     vec4 p1,  vec4 p2,  vec4 p3,                                  \n"
        "     vec4 pp1, vec4 pp2, vec4 pp3                                  \n"
        "  ) {                                                              \n"
        "   if(lighting) {                                                  \n"
        "      vec3 N = gl_NormalMatrix * cross((p2-p1).xyz,(p3-p1).xyz) ;  \n"
        "      vec3 L = gl_LightSource[0].position.xyz ;                    \n"
        "      diffuse = cosangle(N,L) ;                                    \n"
        "      float NdotH = cosangle(N,gl_LightSource[0].halfVector.xyz) ; \n"
        "      specular = pow(abs(NdotH),gl_FrontMaterial.shininess);       \n"
        "   } else {                                                        \n"
        "       diffuse = -1.0 ; specular = 0.0 ;                           \n"
        "   }                                                               \n"
        "   gl_Position=pp1; bary = vec2(0.0,0.0) ; EmitVertex();           \n"
        "   gl_Position=pp2; bary = vec2(1.0,0.0) ; EmitVertex();           \n"
        "   gl_Position=pp3; bary = vec2(0.0,1.0) ; EmitVertex();           \n"
        "   EndPrimitive();                                                 \n"
        "}                                                                  \n"
        "void main() {                                                      \n"
        "    vec4 pp0 = project(gl_in[0].gl_Position);                      \n"
        "    vec4 pp1 = project(gl_in[1].gl_Position);                      \n"
        "    vec4 pp2 = project(gl_in[2].gl_Position);                      \n"
        "    vec4 pp3 = project(gl_in[3].gl_Position);                      \n"
        "    if(clipping) {                                                 \n"
        "       if(                                                         \n"
        "           clip(gl_in[0].gl_Position) < 0.0 &&                     \n"
        "           clip(gl_in[1].gl_Position) < 0.0 &&                     \n"
        "           clip(gl_in[2].gl_Position) < 0.0 &&                     \n"
        "           clip(gl_in[3].gl_Position) < 0.0                        \n"
        "       ) {                                                         \n"
        "            return ;                                               \n"
        "       }                                                           \n"
        "    }                                                              \n"
        "    flat_shaded_triangle(                                          \n"
        "      gl_in[0].gl_Position,                                        \n"
        "      gl_in[1].gl_Position,                                        \n"
        "      gl_in[2].gl_Position,                                        \n"
        "      pp0, pp1, pp2                                                \n"
        "    );                                                             \n"
        "    flat_shaded_triangle(                                          \n"
        "      gl_in[1].gl_Position,                                        \n"
        "      gl_in[0].gl_Position,                                        \n"
        "      gl_in[3].gl_Position,                                        \n"
        "      pp1, pp0, pp3                                                \n"
        "    );                                                             \n"
        "    flat_shaded_triangle(                                          \n"
        "      gl_in[0].gl_Position,                                        \n"
        "      gl_in[2].gl_Position,                                        \n"
        "      gl_in[3].gl_Position,                                        \n"
        "      pp0, pp2, pp3                                                \n"
        "    );                                                             \n"
        "    flat_shaded_triangle(                                          \n"
        "      gl_in[2].gl_Position,                                        \n"
        "      gl_in[1].gl_Position,                                        \n"
        "      gl_in[3].gl_Position,                                        \n"
        "      pp2, pp1, pp3                                                \n"
        "    );                                                             \n"
        "}                                                                  \n"
    ;

    /**
     * \brief The fragment shader for polygons.
     * \details Does front/back shading in different colors
     *  and fragment-level mesh display.
     */
    const char* fshader_source =
        "#version 150 compatibility                                         \n"
        "flat in float diffuse;                                             \n"
        "flat in float specular;                                            \n"
        "in vec2 bary;                                                      \n"
        "uniform bool mesh = false ;                                        \n"
        "uniform vec3 mesh_color = vec3(0.0, 0.0, 0.0) ;                    \n"
        "out vec4 frag_color ;                                              \n"
        "float edge_factor(){                                               \n"
        "    vec3 bary3 = vec3(bary.x, bary.y, 1.0-bary.x-bary.y) ;         \n"
        "    vec3 d = fwidth(bary3);                                        \n"
        "    vec3 a3 = smoothstep(vec3(0.0,0.0,0.0), d*1.3, bary3);         \n"
        "    return min(min(a3.x, a3.y), a3.z);                             \n"
        "}                                                                  \n"
        "void main() {                                                      \n"
        "    float s = gl_FrontFacing ? 1.0 : -1.0 ;                        \n"
        "    vec4 Kdiff = gl_FrontFacing ?                                  \n"
        "         gl_FrontMaterial.diffuse : gl_BackMaterial.diffuse ;      \n"
        "    float sdiffuse = s * diffuse ;                                 \n"
        "    vec4 result = vec4(0.1, 0.1, 0.1, 1.0);                        \n"
        "    if(sdiffuse > 0.0) {                                           \n"
        "       result += sdiffuse*Kdiff +                                  \n"
        "                 specular*gl_FrontMaterial.specular;               \n"
        "    }                                                              \n"
        "    frag_color = mesh ?                                            \n"
        "                  mix(vec4(mesh_color,1.0),result,edge_factor()) : \n"
        "                  result ;                                         \n"
        "}                                                                  \n";

    /**
     * \brief The fragment shader for points.
     * \details Makes the points appear as small spheres.
     */
    const char* points_fshader_source =
        "#version 150 compatibility                                         \n"
        "out vec4 frag_color ;                                              \n"
        "void main() {                                                      \n"
        "   vec2 V = gl_TexCoord[0].xy - vec2(0.5, 0.5);                    \n"
        "   float d = 1.0-4.0*dot(V,V);                                     \n"
        "   if(d < 0.0) {                                                   \n"
        "      discard;                                                     \n"
        "   }                                                               \n"
        "   frag_color = d*gl_Color;                                        \n"
        "}                                                                  \n";

    /************************************************************************/

    /**
     * \brief Compiles a shader for a specific target.
     * \details Errors are detected and displayed to std::err.
     *
     * \return the OpenGL opaque Id of the created shader object.
     */
    GLuint compile_shader(const char* source, GLenum target) {
        GLuint s_handle = glCreateShader(target);
        if(s_handle == 0) {
            GEO::Logger::err("GLSL") << "Could not create shader"
                << std::endl;
            exit(1);
        }
        glShaderSource(s_handle, 1, &source, 0);
        glCompileShader(s_handle);
        GLint compile_status;
        glGetShaderiv(s_handle, GL_COMPILE_STATUS, &compile_status);
        if(!compile_status) {
            GLchar compiler_message[4096];
            glGetShaderInfoLog(
                s_handle, sizeof(compiler_message), 0, compiler_message
            );
            std::cerr << "GLSL compiler status :"
                << compile_status << std::endl;
            std::cerr << "GLSL compiler message:"
                << compiler_message << std::endl;
        }
        return s_handle;
    }

    /**
     * \brief Creates a program from the vertex, geometry and fragment shaders.
     * \details Errors are detected and displayed to std::cerr.
     * \param[in] vshader id of the vertex shader
     * \param[in] gshader id of the geometry shader
     * \param[in] fshader id of the fragment shader
     * \return the id of the program
     */
    GLuint setup_program(GLuint vshader, GLuint gshader, GLuint fshader) {
        GLuint program = glCreateProgram();
        if(vshader != 0) {
            glAttachShader(program, vshader);
        }
        if(gshader != 0) {
            glAttachShader(program, gshader);
        }
        if(fshader != 0) {
            glAttachShader(program, fshader);
        }
        glLinkProgram(program);

        GLint link_status;
        glGetProgramiv(program, GL_LINK_STATUS, &link_status);
        if(!link_status) {
            GLchar linker_message[4096];
            glGetProgramInfoLog(
                program, sizeof(linker_message), 0, linker_message
            );
            std::cerr << "GLSL linker status :" << link_status << std::endl;
            std::cerr << "GLSL linker message:" << linker_message << std::endl;
        }
        return program;
    }

    /**
     * \brief Creates all the shaders used by vorpaview.
     * \details Compiles the three shaders (vertex, geometry and fragment),
     *  links them into a program, and displays the potential error
     *  messages to std::cerr. There is one program for displaying
     *  triangles and one program for displaying tetrahedra.
     */
    void setup_shaders() {
        GLuint vshader = compile_shader(vshader_source, GL_VERTEX_SHADER);
        GLuint gshader_tri =
            compile_shader(gshader_tri_source, GL_GEOMETRY_SHADER);
        GLuint gshader_tet =
            compile_shader(gshader_tet_source, GL_GEOMETRY_SHADER);
        GLuint fshader = compile_shader(fshader_source, GL_FRAGMENT_SHADER);
        GLuint points_fshader = compile_shader(
            points_fshader_source, GL_FRAGMENT_SHADER
        );

        program_points = setup_program(0, 0, points_fshader);
        program_tri = setup_program(vshader, gshader_tri, fshader);
        program_tet = setup_program(vshader, gshader_tet, fshader);
    }

    /**
     * \brief Zooms in.
     * \details Zooming factor is 1.1x.
     */
    void zoom_in() {
        *glut_viewer_float_ptr(GLUT_VIEWER_ZOOM) *= 1.1f;
    }

    /**
     * \brief Zooms out.
     * \details De-zooming factor is (1/1.1)x.
     */
    void zoom_out() {
        *glut_viewer_float_ptr(GLUT_VIEWER_ZOOM) /= 1.1f;
    }

    /**
     * \brief Toggles black or white background color.
     */
    void toggle_background() {
        white_bg = !white_bg;
        if(white_bg) {
            glut_viewer_set_background_color(1.0, 1.0, 1.0);
        } else {
            glut_viewer_set_background_color(0.0, 0.0, 0.0);
        }
        if(program_tri != 0) {
            GLint loc = glGetUniformLocation(program_tri, "mesh_color");
            glUseProgram(program_tri);
            if(white_bg) {
                glUniform3f(loc, 0.0, 0.0, 0.0);
            } else {
                glUniform3f(loc, 1.0, 1.0, 1.0);
            }
            glUseProgram(0);
        }
        if(program_tet != 0) {
            GLint loc = glGetUniformLocation(program_tet, "mesh_color");
            glUseProgram(program_tet);
            if(white_bg) {
                glUniform3f(loc, 0.0, 0.0, 0.0);
            } else {
                glUniform3f(loc, 1.0, 1.0, 1.0);
            }
            glUseProgram(0);
        }
    }

    /**
     * \brief Toggles mesh display.
     */
    void toggle_mesh() {
        show_mesh = !show_mesh;
        if(program_tri != 0 && M.is_triangulated()) {
            glUseProgram(program_tri);
            GLint loc = glGetUniformLocation(program_tri, "mesh");
            glUniform1i(loc, show_mesh);
            glUseProgram(0);
        }
        if(program_tet != 0) {
            glUseProgram(program_tet);
            GLint loc = glGetUniformLocation(program_tet, "mesh");
            glUniform1i(loc, show_mesh);
            glUseProgram(0);
        }
    }

    /**
     * \brief Toggles color / BW display.
     */
    void toggle_colors() {
        show_colors = !show_colors;
    }

    /**
     * \brief Toggles lighting / constant color mode.
     */
    void toggle_lighting() {
        lighting = !lighting;
        if(program_tri != 0) {
            glUseProgram(program_tri);
            GLint loc = glGetUniformLocation(program_tri, "lighting");
            glUniform1i(loc, lighting);
            glUseProgram(0);
        }
        if(program_tet != 0) {
            glUseProgram(program_tet);
            GLint loc = glGetUniformLocation(program_tet, "lighting");
            glUniform1i(loc, lighting);
            glUseProgram(0);
        }
    }

    /**
     * \brief Initializes OpenGL objects.
     * \details Specifed as glut_viewer_set_init_func() callback.
     */
    void init() {
        glewInit();
        setup_shaders();
        glut_viewer_set_background_color(1.0, 1.0, 1.0);
        glut_viewer_add_toggle(
            'T', glut_viewer_is_enabled_ptr(GLUT_VIEWER_TWEAKBARS),
            "Toggle tweakbars"
        );
        glut_viewer_add_key_func('b', toggle_background, "Toggle background");
        glut_viewer_add_key_func('c', toggle_colors, "colors");
        glut_viewer_add_toggle('B', &show_borders, "borders");
        glut_viewer_add_key_func('z', zoom_in, "Zoom in");
        glut_viewer_add_key_func('Z', zoom_out, "Zoom out");
        glut_viewer_disable(GLUT_VIEWER_TWEAKBARS);
        glut_viewer_disable(GLUT_VIEWER_BACKGROUND);
        glut_viewer_add_key_func('m', toggle_mesh, "mesh");

        if(GLEW_EXT_geometry_shader4) {
            GEO::Logger::out("OpenGL")
                << "Geometry shaders are supported, good!"
                << std::endl;
        } else {
            GEO::Logger::warn("OpenGL")
                << "Geometry shaders are not supported"
                << std::endl;
        }
    }

    /**
     * \brief Creates OpenGL-side arrays (VertexBufferObjects).
     * \details The first call creates the VertexBufferObjects
     *  that store the point coordinates and the triangle
     *  corner indices in the graphic board.
     */
    void create_VBOs_if_needed() {
        if(vertices_VBO == 0 && M.nb_vertices() != 0) {
            glGenBuffers(1, &vertices_VBO);
            glBindBuffer(GL_ARRAY_BUFFER, vertices_VBO);
            glBufferData(
                GL_ARRAY_BUFFER,
                GLsizeiptr(M.nb_vertices() * M.dimension() * sizeof(float)),
                M.vertex_ptr(0), GL_STATIC_DRAW
            );
            glEnableClientState(GL_VERTEX_ARRAY);
            unsigned int stride =
                (unsigned int) (M.dimension() * sizeof(float));
            geo_assert(M.dimension() == 3);
            glVertexPointer(3, GL_FLOAT, GLsizei(stride), 0);
        }
        if(facet_indices_VBO == 0 && M.nb_facets() != 0) {
            glGenBuffers(1, &facet_indices_VBO);
            glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, facet_indices_VBO);
            glBufferData(
                GL_ELEMENT_ARRAY_BUFFER,
                GLsizeiptr(M.nb_corners() * sizeof(int)),
                M.corner_vertex_index_ptr(0), GL_STATIC_DRAW
            );
        }
        if(tet_indices_VBO == 0 && M.is_tetrahedralized() && M.nb_tets() != 0) {
            glGenBuffers(1, &tet_indices_VBO);
            glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, tet_indices_VBO);
            glBufferData(
                GL_ELEMENT_ARRAY_BUFFER,
                GLsizeiptr(M.nb_tets() * 4 * sizeof(int)),
                M.tet_vertex_index_ptr(0, 0), GL_STATIC_DRAW
            );
        }
    }

    /**
     * \brief Binds the triangles array as current Vertex Buffer Object.
     * \details Creates the triangles Vertex Buffer Object if need be.
     */
    void bind_triangles_VBO() {
        create_VBOs_if_needed();
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, facet_indices_VBO);
    }

    /**
     * \brief Binds the tetrahedra array as current Vertex Buffer Object.
     * \details Creates the tetrahedra Vertex Buffer Object if need be.
     */
    void bind_tets_VBO() {
        create_VBOs_if_needed();
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, tet_indices_VBO);
    }

    /**
     * \brief Deletes all the previously created Vertex Buffer Objects.
     */
    void delete_VBOs_if_needed() {
        if(vertices_VBO != 0) {
            glDeleteBuffers(1, &vertices_VBO);
            vertices_VBO = 0;
        }
        if(facet_indices_VBO != 0) {
            glDeleteBuffers(1, &facet_indices_VBO);
            facet_indices_VBO = 0;
        }
        if(tet_indices_VBO != 0) {
            glDeleteBuffers(1, &tet_indices_VBO);
            tet_indices_VBO = 0;
        }
    }

    /**
     * \brief Draws a triangle mesh using VertexBufferObjects.
     * \details The implementation uses a Geometry Shader (does
     *   not work if the graphic board is too old).
     */
    void draw_triangles_VBOs() {
        if(M.nb_facets() == 0) {
            return;
        }
        bind_triangles_VBO();
        // Note: the fourth argument (0) corresponds to the bound VBO.
        glDrawElements(
            GL_TRIANGLES, GLsizei(M.nb_corners()), GL_UNSIGNED_INT, 0
        );
    }

    /**
     * \brief Draws a tetrahedral mesh using VertexBufferObjects.
     * \details The implementation uses a Geometry Shader (does
     *  not work if the graphic board is too old).
     */
    void draw_tets_VBOs() {
        if(M.nb_tets() == 0) {
            return;
        }
        bind_tets_VBO();
        // Note: the fourth argument (0) corresponds to the bound VBO.
        glDrawElements(
            GL_LINES_ADJACENCY, GLsizei(M.nb_tets() * 4), GL_UNSIGNED_INT, 0
        );
    }

    /**
     * \brief Sends the normal to a triangle to OpenGL
     * \param[in] M a const reference to the mesh
     * \param[in] i index of the first vertex of the triangle
     * \param[in] j index of the second vertex of the triangle
     * \param[in] k index of the third vertex of the triangle
     */

    void glMeshTriangleNormal(
        const GEO::SinglePrecisionMesh& M,
        GEO::index_t i, GEO::index_t j, GEO::index_t k
    ) {
        const float* pi = M.vertex_ptr(i);
        const float* pj = M.vertex_ptr(j);
        const float* pk = M.vertex_ptr(k);
        float N[3];
        N[0] = (pi[1]-pj[1])*(pk[2]-pj[2]) - (pi[2]-pj[2])*(pk[1]-pj[1]);
        N[1] = (pi[2]-pj[2])*(pk[0]-pj[0]) - (pi[0]-pj[0])*(pk[2]-pj[2]);
        N[2] = (pi[0]-pj[0])*(pk[1]-pj[1]) - (pi[1]-pj[1])*(pk[0]-pj[0]);
        glNormal3fv(N);
    }
    
    /**
     * \brief Draws the cells of a hybrid mesh.
     */
    void draw_cells() {

        static bool cell_type_has_triangles[GEO::MESH_NB_CELL_TYPES] = {
            true,  // tets
            false, // hexes
            true,  // prisms
            true,  // pyramids
            true   // connectors
        };

        static bool cell_type_has_quads[GEO::MESH_NB_CELL_TYPES] = {
            false,  // tets
            true, // hexes
            true,  // prisms
            true,  // pyramids
            true   // connectors
        };

        glCullFace(GL_FRONT);
        glEnable(GL_CULL_FACE);
        
        // Try using loaded vertex array with glDrawElements() 
        // Try using a shader, and activating by-cells clipping
        //   (we may need a shader per cell type, or a shader per
        //    facet type)

        glBegin(GL_TRIANGLES);
        for(GEO::index_t c=0; c<M.nb_cells(); ++c) {
            GEO::MeshCellType type = M.cell_type(c);
            if(!cell_type_has_triangles[type]) {
                continue;
            }
            for(GEO::index_t f=0; f<M.cell_nb_facets(c); ++f) {
                if(M.cell_facet_nb_vertices(c,f) == 3) {
                    GEO::index_t i = M.cell_facet_vertex_index(c,f,0);
                    GEO::index_t j = M.cell_facet_vertex_index(c,f,1);
                    GEO::index_t k = M.cell_facet_vertex_index(c,f,2);
                    glMeshTriangleNormal(M,i,j,k);
                    glVertex3fv(M.vertex_ptr(i));
                    glVertex3fv(M.vertex_ptr(j));
                    glVertex3fv(M.vertex_ptr(k));                    
                }
            }
        }
        glEnd();
        
        glBegin(GL_QUADS);
        for(GEO::index_t c=0; c<M.nb_cells(); ++c) {
            GEO::MeshCellType type = M.cell_type(c);            
            if(!cell_type_has_quads[type]) {
                continue;
            }
            for(GEO::index_t f=0; f<M.cell_nb_facets(c); ++f) {
                if(M.cell_facet_nb_vertices(c,f) == 4) {
                    GEO::index_t i = M.cell_facet_vertex_index(c,f,0);
                    GEO::index_t j = M.cell_facet_vertex_index(c,f,1);
                    GEO::index_t k = M.cell_facet_vertex_index(c,f,2);
                    GEO::index_t l = M.cell_facet_vertex_index(c,f,3);
                    glMeshTriangleNormal(M,i,j,k);
                    glVertex3fv(M.vertex_ptr(i));
                    glVertex3fv(M.vertex_ptr(j));
                    glVertex3fv(M.vertex_ptr(k));
                    glVertex3fv(M.vertex_ptr(l));
                }
            }
        }
        glEnd();

        glDisable(GL_CULL_FACE);        
    }
    
    /**
     * \brief Draws a polygon mesh using VertexBufferObjects.
     *
     * \note This function needs optimization (for now, it issues
     * one OpenGL call per polygon, there is probably a means
     * of avoiding that...).
     */
    void draw_polygons_VBOs() {
        if(M.nb_facets() == 0) {
            return;
        }
        create_VBOs_if_needed();
        for(unsigned int f = 0; f < M.nb_facets(); f++) {
            unsigned int b = M.facet_begin(f);
            unsigned int e = M.facet_end(f);
            // Note: the fourth argument (void*)(b*sizeof(int))
            // is relative to the address of the bound
            // VBO.
            glDrawElements(
                GL_POLYGON, GLsizei(e - b),
                GL_UNSIGNED_INT, (void*) (b * sizeof(int))
            );
        }
    }

    /**
     * \brief Draws a surfacic mesh.
     * \details This function is optimized for triangle meshes, and
     *   much much slower for polygon meshes.
     */
    void draw_surface() {
        if(M.nb_facets() == 0) {
            return;
        }
        if(M.is_triangulated()) {
            draw_triangles_VBOs();
        } else {
            draw_polygons_VBOs();
        }
    }

    /**
     * \brief Draws the borders of a surface, and the borders of
     *   the facet regions if present.
     * \details This function is optimized for triangle meshes, and
     *   much much slower for polygon meshes.
     */
    void draw_surface_borders() {
        if(M.nb_facets() == 0) {
            return;
        }
        glDisable(GL_LIGHTING);
        glLineWidth(3);
        if(white_bg) {
            glColor3f(0.0, 0.0, 0.0);
        } else {
            glColor3f(1.0, 1.0, 1.0);
        }
        glBegin(GL_LINES);
        for(GEO::index_t f = 0; f < M.nb_facets(); f++) {
            for(GEO::index_t c1 = M.facet_begin(f);
                c1 < M.facet_end(f); ++c1
            ) {
                if(M.corner_adjacent_facet(c1) == -1) {
                    unsigned int c2 = M.next_around_facet(f, c1);
                    glVertex3fv(M.vertex_ptr(M.corner_vertex_index(c1)));
                    glVertex3fv(M.vertex_ptr(M.corner_vertex_index(c2)));
                }
            }
        }
        glEnd();

        if(!M.has_attribute(GEO::MESH_FACET_REGION)) {
            return;
        }

        glLineWidth(2);
        glBegin(GL_LINES);
        for(GEO::index_t f = 0; f < M.nb_facets(); f++) {
            GEO::signed_index_t f_rgn = M.facet_region(f);
            for(GEO::index_t c1 = M.facet_begin(f);
                c1 < M.facet_end(f); ++c1
            ) {
                GEO::signed_index_t adj_f = M.corner_adjacent_facet(c1);
                if(adj_f != -1 && M.facet_region(GEO::index_t(adj_f)) != f_rgn) {
                    unsigned int c2 = M.next_around_facet(f, c1);
                    glVertex3fv(M.vertex_ptr(M.corner_vertex_index(c1)));
                    glVertex3fv(M.vertex_ptr(M.corner_vertex_index(c2)));
                }
            }
        }
        glEnd();
    }

    /**
     * \brief Draws a volumetric mesh (with tetrahedra).
     */
    void draw_volume() {
        if(M.nb_cells() == 0) {
            return;
        }
        glCullFace(GL_FRONT);
        glEnable(GL_CULL_FACE);
        if(M.is_tetrahedralized()) {
            draw_tets_VBOs();
        }
        glDisable(GL_CULL_FACE);
    }

    /**
     * \brief Draws all the vertices of the mesh from the Vertex Buffer Object.
     */
    void draw_points_VBO() {
        create_VBOs_if_needed();
        glDrawArrays(GL_POINTS, 0, GLsizei(M.nb_vertices()));
    }

    /**
     * \brief Draws all the vertices of the mesh.
     */
    void draw_points() {
        glDisable(GL_LIGHTING);
        glPointSize(15);
        glEnable(GL_POINT_SPRITE);
        glTexEnvi(GL_POINT_SPRITE, GL_COORD_REPLACE, GL_TRUE);
        glColor3f(0.0f, 1.0f, 0.0f);
        glUseProgram(program_points);
        draw_points_VBO();
        glUseProgram(0);
    }

    /**
     * \brief Draws the current mesh according to current rendering mode.
     * \details Specifed as glut_viewer_set_display_func() callback.
     */
    void display() {
        if(program_tri != 0) {
            glUseProgram(program_tri);
            GLint loc = glGetUniformLocation(program_tri, "mesh");
            glUniform1i(loc, M.is_triangulated() && show_mesh);
            glUseProgram(0);
        }
        if(program_tet != 0) {
            glUseProgram(program_tet);
            GLint loc = glGetUniformLocation(program_tet, "mesh");
            glUniform1i(loc, show_mesh);
            glUseProgram(0);
        }

        if(show_borders) {
            draw_surface_borders();
        }

        GLfloat shininess = 20.0f;
        glMaterialfv(GL_FRONT_AND_BACK, GL_SHININESS, &shininess);
        if(show_colors) {
            if(M.nb_cells() == 0) {
                static float diff_front[4] = {0.5f, 0.75f, 1.0f, 1.0f};
                static float diff_back[4] = {1.0f, 0.0f, 0.0f, 1.0f};
                glMaterialfv(GL_FRONT, GL_DIFFUSE, diff_front);
                glMaterialfv(GL_BACK, GL_DIFFUSE, diff_back);
                static float spec[4] = {0.6f, 0.6f, 0.6f, 1.0f};
                glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, spec);
            } else {
                static float diff_back[4] = {1.0f, 1.0f, 0.0f, 0.7f};
                static float diff_front[4] = {0.7f, 0.0f, 0.0f, 1.0f};
                glMaterialfv(GL_FRONT, GL_DIFFUSE, diff_front);
                glMaterialfv(GL_BACK, GL_DIFFUSE, diff_back);
                static float spec[4] = {0.6f, 0.6f, 0.6f, 1.0f};
                glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, spec);
            }
        } else {
            static float spec[4] = {0.6f, 0.6f, 0.6f, 1.0f};
            glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, spec);
            if(white_bg) {
                static float diff[4] = {0.9f, 0.9f, 0.9f, 1.0f};
                glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, diff);
            } else {
                static float diff[4] = {0.1f, 0.1f, 0.1f, 1.0f};
                glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, diff);
            }
        }

        if(show_surface) {
            glUseProgram(program_tri);
            glColor3f(0.0f, 0.5f, 1.0f);
            draw_surface();
            glUseProgram(0);
        }

        if(show_volume) {
            static float tet_color[4] = {1.0f, 1.0f, 1.0f, 1.0f};
            glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, tet_color);

            if(M.is_tetrahedralized()) {
                glUseProgram(program_tet);
                draw_volume();
                glUseProgram(0);
            } else {
                glEnable(GL_LIGHTING);
                draw_cells();
                glDisable(GL_LIGHTING);                
            }
        }

        // If the surface is triangulated, then the mesh
        //  is drawn by the fragment shader (that
        //  changes the color of the fragments near
        //  the edges of the triangles),
        // If the surface has polygons, then the mesh
        //  is drawn "the standard way" below:
        if(show_mesh && !M.is_triangulated()) {
            glDisable(GL_LIGHTING);
            glLineWidth(1);
            if(white_bg) {
                glColor3f(0.0, 0.0, 0.0);
            } else {
                glColor3f(1.0, 1.0, 1.0);
            }
            glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
            draw_surface();
            glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
        }


        if(show_mesh && show_volume && !M.is_tetrahedralized()) {
            glDisable(GL_LIGHTING);
            glLineWidth(1);
            if(white_bg) {
                glColor3f(0.0, 0.0, 0.0);
            } else {
                glColor3f(1.0, 1.0, 1.0);
            }
            glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
            draw_cells();
            glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
        }
        
        if(show_vertices) {
            draw_points();
        }
    }

    /**
     * \brief Loads a mesh from a file.
     */
    void load_mesh(const std::string& filename) {

        delete_VBOs_if_needed();
        GEO::MeshIOFlags flags;
        flags.set_attribute(GEO::MESH_FACET_REGION);
        flags.set_element(GEO::MESH_CELLS);
        if(!GEO::mesh_load(filename, M, flags)) {
            // We need to specify dimension, else
            // we got an invalid mesh that triggers
            // an assertion fail later.
            M.set_dimension(3);
            return;
        }
        if(M.nb_facets() == 0 && M.nb_cells() > 0) {
            M.compute_cells_boundaries();
        }
        M.set_dimension(3);
        double xyzmin[3];
        double xyzmax[3];
        GEO::get_bbox(M, xyzmin, xyzmax);
        glut_viewer_set_region_of_interest(
            float(xyzmin[0]), float(xyzmin[1]), float(xyzmin[2]),
            float(xyzmax[0]), float(xyzmax[1]), float(xyzmax[2])
        );


/*
  // TODO: remove it
        {
            std::cerr << "Computing borders.obj" << std::endl;
            GEO::Mesh M2;
            GEO::MeshIOFlags flags;
            flags.set_attribute(GEO::MESH_FACET_REGION);
            flags.set_element(GEO::MESH_CELLS);
            GEO::mesh_load(filename,M2,flags);
            M2.set_dimension(3);
            M2.compute_cells_boundaries();
            GEO::mesh_save(M2,"borders.obj");
            std::cerr << "Computed borders.obj" << std::endl;            
        }
*/
    }

    /**
     * \brief Loads a mesh from a file icon dropped into the window.
     * \details Specifed as glut_viewer_set_drag_drop_func() callback.
     */
    void dropped_file_cb(char* filename) {
        load_mesh(std::string(filename));
    }

    /**
     * \brief Inverts the normals of a mesh.
     * \details In color mode, this swaps the red and the blue sides.
     */
    void invert_normals() {
        delete_VBOs_if_needed();
        for(GEO::index_t f = 0; f < M.nb_facets(); ++f) {
            GEOGen::MeshMutator<GEO::SinglePrecisionMesh>::flip_facet(M, f);
        }
    }
}

int main(int argc, char** argv) {

    GEO::initialize();
    GEO::Logger::instance()->set_quiet(false);

    GEO::CmdLine::import_arg_group("standard");
    GEO::CmdLine::import_arg_group("algo");
    GEO::CmdLine::declare_arg("full_screen",false,"full screen mode");
   
    std::vector<std::string> filenames;
    if(!GEO::CmdLine::parse(argc, argv, filenames, "<filename>")) {
        return 1;
    }

    if(filenames.size() == 1) {
        load_mesh(filenames[0]);
    } else {
        M.set_dimension(3);
    }

    if(M.nb_cells() != 0) {
        show_volume = GL_TRUE;
        show_mesh = GL_TRUE;
    } else if(M.nb_facets() == 0) {
        show_vertices = GL_TRUE;
    }

    glut_viewer_set_window_title(
        (char*) "[ \\V (O |R |P /A |L |I |N |E ]-[ viewer ]"
    );
    glut_viewer_set_init_func(init);
    glut_viewer_set_display_func(display);
    glut_viewer_set_drag_drop_func(dropped_file_cb);
    glut_viewer_add_toggle('V', &show_volume, "volume");
    glut_viewer_add_toggle('p', &show_vertices, "vertices");
    glut_viewer_add_toggle('S', &show_surface, "surface");
    glut_viewer_add_key_func('L', toggle_lighting, "toggle lighting");
    glut_viewer_add_key_func('n', invert_normals, "invert normals");
   
    if(GEO::CmdLine::get_arg_bool("full_screen")) {
       glut_viewer_enable(GLUT_VIEWER_FULL_SCREEN);
    }
      
    glut_viewer_main_loop(argc, argv);

    // Note: when 'q' is pressed, exit() is called
    // because there is no simple way of exiting from
    // glut's event loop, therefore this line is not
    // reached and memory is not properly freed on exit.
    // TODO: add a function in freeglut to exit event loop.

    return 0;
}

