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

#include <geogram_gfx/mesh/mesh_gfx.h>
#include <geogram_gfx/basic/GLSL.h>
#include <geogram/mesh/mesh.h>
#include <geogram/basic/command_line.h>
#include <geogram/basic/logger.h>

namespace {
    using namespace GEO;

    /**
     * \brief The fragment shader for polygons.
     * \details Used by all polygons and cells programs. 
     *  Does front/back shading in different colors and 
     *  fragment-level mesh display.
     */
    const char* fshader_source =
        "#version 150 compatibility                                         \n"
        "flat in float diffuse;                                             \n"
        "flat in float specular;                                            \n"
        "flat in vec3  edge_mask;                                           \n"
        "in vec2 bary;                                                      \n"
        "uniform bool mesh = true ;                                         \n"
        "uniform vec3 mesh_color = vec3(0.0, 0.0, 0.0) ;                    \n"
        "uniform bool lighting = true ;                                     \n"
        "out vec4 frag_color ;                                              \n"
        "float edge_factor(){                                               \n"
        "    vec3 bary3 = vec3(bary.x, bary.y, 1.0-bary.x-bary.y) ;         \n"
        "    vec3 d = fwidth(bary3);                                        \n"
        "    vec3 a3 = smoothstep(vec3(0.0,0.0,0.0), d*1.3, bary3);         \n"
        "    a3 = vec3(1.0, 1.0, 1.0) - edge_mask + edge_mask*a3;           \n"
        "    return min(min(a3.x, a3.y), a3.z);                             \n"
        "}                                                                  \n"
        "void main() {                                                      \n"
        "    float s = (lighting && gl_FrontFacing) ? 1.0 : -1.0 ;          \n"
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
     * \details Uses vshader_pass_through. Makes the points appear 
     *  as small spheres.
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
    
    /** 
     * \brief Some utility functions for the geometry shaders.
     * \details Provides functions for clipping, projection, and
     *  for generating shaded polygons.
     *  - flat_shaded_triangle(p1,p2,p3,pp1,pp2,pp3,do_clip) where
     *   (p1,p2,p3) are the coordinates in world space, (pp1,pp2,pp3) the
     *   transformed coordinates in clip space and do_clip specifies whether
     *   the triangle should be clipped.
     *  - flat_shaded_quad(p1,p2,p3,p4,pp1,pp2,pp3,pp4,do_clip,edges)
     */
    const char* gshader_utils_source =
        "#version 150 compatibility                                         \n"
        "flat out float diffuse;                                            \n"
        "flat out float specular;                                           \n"
        "flat out vec3 edge_mask;                                           \n"
        "uniform bool lighting=true;                                        \n"
        "uniform bool clipping=true;                                        \n"
        "out float gl_ClipDistance[1];                                      \n"
        "out vec2 bary;                                                     \n"
        "                                                                   \n"
        "vec4 project(vec4 V) {                                             \n"
        "   return gl_ModelViewProjectionMatrix * V ;                       \n"
        "}                                                                  \n"
        "float clip(vec4 V, bool do_clip) {                                 \n"
        "  return do_clip ? dot(gl_ModelViewMatrix*V,gl_ClipPlane[0]) :0.0; \n"
        "}                                                                  \n"
        "float cosangle(vec3 N, vec3 L) {                                   \n"
        "   float s = inversesqrt(dot(N,N)*dot(L,L)) ;                      \n"
        "   return s*dot(N,L) ;                                             \n"
        "}                                                                  \n"
        "void flat_shaded_triangle(                                         \n"
        "     vec4 p1,  vec4 p2,  vec4 p3,                                  \n"
        "     vec4 pp1, vec4 pp2, vec4 pp3,                                 \n"
        "     bool do_clip                                                  \n"
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
        "   edge_mask = vec3(1.0,1.0,1.0);                                  \n"
        "   gl_ClipDistance[0] = clip(p1, do_clip);                         \n"
        "   gl_Position=pp1; bary = vec2(1.0,0.0) ; EmitVertex();           \n"
        "   gl_ClipDistance[0] = clip(p2, do_clip);                         \n"
        "   gl_Position=pp2; bary = vec2(0.0,1.0) ; EmitVertex();           \n"
        "   gl_ClipDistance[0] = clip(p3, do_clip);                         \n"
        "   gl_Position=pp3; bary = vec2(0.0,0.0) ; EmitVertex();           \n"
        "   EndPrimitive();                                                 \n"
        "}                                                                  \n"
        "                                                                   \n"
        "void flat_shaded_quad(                                             \n"
        "     vec4 p1,  vec4 p2,  vec4 p3, vec4 p4,                         \n"
        "     vec4 pp1, vec4 pp2, vec4 pp3, vec4 pp4,                       \n"
        "     bool do_clip                                                  \n"
        "  ) {                                                              \n"
        "   if(lighting) {                                                  \n"
        "      vec3 N = gl_NormalMatrix * (                                 \n" 
        "           cross((p2-p1).xyz,(p4-p1).xyz) +                        \n"
        "           cross((p4-p3).xyz,(p2-p3).xyz)                          \n"
        "           );                                                      \n"
        "      vec3 L = gl_LightSource[0].position.xyz ;                    \n"
        "      diffuse = cosangle(N,L) ;                                    \n"
        "      float NdotH = cosangle(N,gl_LightSource[0].halfVector.xyz) ; \n"
        "      specular = pow(abs(NdotH),gl_FrontMaterial.shininess);       \n"
        "   } else {                                                        \n"
        "       diffuse = -1.0 ; specular = 0.0 ;                           \n"
        "   }                                                               \n"
        "   float cl1 = clip(p1,do_clip);                                   \n"
        "   float cl2 = clip(p2,do_clip);                                   \n"
        "   float cl3 = clip(p3,do_clip);                                   \n"
        "   float cl4 = clip(p4,do_clip);                                   \n"
        "   edge_mask = vec3(0.0, 1.0, 1.0);                                \n"
        "   gl_ClipDistance[0]=cl1; gl_Position=pp1; bary=vec2(1.0,0.0);    \n"
        "   EmitVertex();                                                   \n"
        "   gl_ClipDistance[0]=cl2; gl_Position=pp2; bary=vec2(0.0,1.0);    \n"
        "   EmitVertex();                                                   \n"
        "   gl_ClipDistance[0]=cl4; gl_Position=pp4; bary=vec2(0.0,0.0);    \n"
        "   EmitVertex();                                                   \n"
        "   edge_mask = vec3(0.0, 1.0, 1.0);                                \n"
        "   gl_ClipDistance[0]=cl3; gl_Position=pp3; bary=vec2(1.0,0.0);    \n"
        "   EmitVertex();                                                   \n"
        "   EndPrimitive();                                                 \n"
        "}                                                                  \n"
        ;
    
    /**
     * \brief The pass-through vertex shader.
     * \details Used by points, quads, tets, prisms
     */
    const char* vshader_pass_through_source =
        "#version 150 compatibility                                         \n"
        " void main(void) {                                                 \n"
        "    gl_Position = gl_Vertex;                                       \n"
        " }                                                                 \n";


    /**
     * \brief The tesselation evaluation shader for hexes.
     * \details Gathers the eight vertices into the array passed
     *  to the geometry shader.
     *  Each hex is generated twice, for the second one, the output
     *  boolean discardme is set to one, to indicate to the geometry
     *  shader that the hex should be discarded (sounds stupid, but
     *  I did not find a smarter way...).
     */
    const char* teshader_hex_source =
        "#version 400                                                       \n"
        "layout(isolines, point_mode) in;                                   \n"
        "const int nb_vertices = 8;                                         \n"
        "out Data {                                                         \n"
        "    vec3 p[nb_vertices];                                           \n"
        "    bool discardme;                                                \n"
        "} DataOut;                                                         \n"
        "void main() {                                                      \n"
        "   if(gl_TessCoord.x == 0.0) {                                     \n"
        "      DataOut.discardme = true;                                    \n"
        "      return;                                                      \n"
        "   }                                                               \n"
        "   DataOut.discardme = false;                                      \n"
        "   for(int i=0; i<nb_vertices; ++i) {                              \n"
        "     DataOut.p[i] = gl_in[i].gl_Position.xyz;                      \n"
        "   }                                                               \n"
        "}                                                                  \n";
    
    /**
     * \brief The vertex shader for hexes if tesselation shader cannot be used.
     * \details For hexes, the pass-through vertex shader
     *  cannot be used, since there is no standard OpenGL primitive
     *  with eight vertices (maximum is 6). We use GL_POINTS and 
     *  pass the 8 coordinates as a generic attribute.
     */
    const char* vshader_hex_source =
        "#version 150 compatibility                                         \n"
        "#extension GL_ARB_explicit_attrib_location : enable                \n"
        " const int nb_vertices = 8;                                        \n"
        " layout(location=0) in vec3 p_in[nb_vertices];                     \n"
        " out Data {                                                        \n"
        "     vec3 p[nb_vertices];                                          \n"
        "     bool discardme;                                               \n"
        " } DataOut;                                                        \n"
        " void main(void) {                                                 \n"
        "   DataOut.discardme = false;                                      \n"
        "   for(int i=0; i<nb_vertices; ++i) {                              \n"
        "       DataOut.p[i] = p_in[i];                                     \n"
        "   }                                                               \n"
        " }                                                                 \n";

    /**
     * \brief The geometry shader for hexes.
     * \details Uses vshader_hex and gshader_utils 
     */
    const char* gshader_hex_source =
        "layout(points) in;                                                 \n"
        "layout(triangle_strip, max_vertices = 24) out;                     \n"
        " in Data {                                                         \n"
        "     vec3 p[8];                                                    \n"
        "     bool discardme;                                               \n"
        " } DataIn[] ;                                                      \n"
        "uniform float shrink = 0.0;                                        \n"
        "                                                                   \n"
        "void draw_hex(                                                     \n"
        "   vec4 p0, vec4 p1, vec4 p2, vec4 p3,                             \n"
        "   vec4 p4, vec4 p5, vec4 p6, vec4 p7                              \n"
        ") {                                                                \n"
        "    if(clipping) {                                                 \n"
        "       if(                                                         \n"
        "           clip(p0,true) < 0.0 && clip(p1,true) < 0.0 &&           \n"
        "           clip(p2,true) < 0.0 && clip(p3,true) < 0.0 &&           \n"
        "           clip(p4,true) < 0.0 && clip(p5,true) < 0.0 &&           \n"
        "           clip(p6,true) < 0.0 && clip(p7,true) < 0.0              \n"
        "       ) {                                                         \n"
        "            return ;                                               \n"
        "       }                                                           \n"
        "    }                                                              \n"
        "   if(shrink != 0.0) {                                             \n"
        "      vec4 g = (1.0/8.0)*(p0+p1+p2+p3+p4+p5+p6+p7);                \n"
        "      p0 = shrink*g + (1.0-shrink)*p0;                             \n"
        "      p1 = shrink*g + (1.0-shrink)*p1;                             \n"
        "      p2 = shrink*g + (1.0-shrink)*p2;                             \n"
        "      p3 = shrink*g + (1.0-shrink)*p3;                             \n"
        "      p4 = shrink*g + (1.0-shrink)*p4;                             \n"
        "      p5 = shrink*g + (1.0-shrink)*p5;                             \n"
        "      p6 = shrink*g + (1.0-shrink)*p6;                             \n"
        "      p7 = shrink*g + (1.0-shrink)*p7;                             \n"
        "   }                                                               \n"
        "   vec4 pp0 = project(p0);                                         \n"
        "   vec4 pp1 = project(p1);                                         \n"
        "   vec4 pp2 = project(p2);                                         \n"
        "   vec4 pp3 = project(p3);                                         \n"
        "   vec4 pp4 = project(p4);                                         \n"
        "   vec4 pp5 = project(p5);                                         \n"
        "   vec4 pp6 = project(p6);                                         \n"
        "   vec4 pp7 = project(p7);                                         \n"
        "   flat_shaded_quad(p0,p2,p6,p4,pp0,pp2,pp6,pp4,false);            \n"
        "   flat_shaded_quad(p3,p1,p5,p7,pp3,pp1,pp5,pp7,false);            \n"
        "   flat_shaded_quad(p1,p0,p4,p5,pp1,pp0,pp4,pp5,false);            \n"
        "   flat_shaded_quad(p2,p3,p7,p6,pp2,pp3,pp7,pp6,false);            \n"
        "   flat_shaded_quad(p1,p3,p2,p0,pp1,pp3,pp2,pp0,false);            \n"
        "   flat_shaded_quad(p4,p6,p7,p5,pp4,pp6,pp7,pp5,false);            \n"
        "}                                                                  \n"
        "                                                                   \n"
        "void main() {                                                      \n"
        "   if(DataIn[0].discardme) { return; }                             \n"
        "   draw_hex(                                                       \n"
        "      vec4(DataIn[0].p[0],1.0),                                    \n"
        "      vec4(DataIn[0].p[1],1.0),                                    \n"
        "      vec4(DataIn[0].p[2],1.0),                                    \n"
        "      vec4(DataIn[0].p[3],1.0),                                    \n"
        "      vec4(DataIn[0].p[4],1.0),                                    \n"
        "      vec4(DataIn[0].p[5],1.0),                                    \n"
        "      vec4(DataIn[0].p[6],1.0),                                    \n"
        "      vec4(DataIn[0].p[7],1.0)                                     \n"
        "   );                                                              \n"
        "}                                                                  \n";


    /**
     * \brief The tesselation evaluation shader for pyramids.
     * \details Gathers the eight vertices into the array passed
     *  to the geometry shader.
     *  Each pyramid is generated twice, for the second one, the output
     *  boolean discardme is set to one, to indicate to the geometry
     *  shader that the hex should be discarded (sounds stupid, but
     *  I did not find a smarter way...).
     */
    const char* teshader_pyramid_source =
        "#version 400                                                       \n"
        "layout(isolines, point_mode) in;                                   \n"
        "const int nb_vertices = 5;                                         \n"
        "out Data {                                                         \n"
        "    vec3 p[nb_vertices];                                           \n"
        "    bool discardme;                                                \n"
        "} DataOut;                                                         \n"
        "void main() {                                                      \n"
        "   if(gl_TessCoord.x == 0.0) {                                     \n"
        "      DataOut.discardme = true;                                    \n"
        "      return;                                                      \n"
        "   }                                                               \n"
        "   DataOut.discardme = false;                                      \n"
        "   for(int i=0; i<nb_vertices; ++i) {                              \n"
        "     DataOut.p[i] = gl_in[i].gl_Position.xyz;                      \n"
        "   }                                                               \n"
        "}                                                                  \n";

    
    /**
     * \brief The vertex shader for pyramids if tesselation canot be used.
     * \details For pyramids, the pass-through vertex shader
     *  cannot be used, since there is no standard OpenGL primitive
     *  with five vertices. We use GL_POINTS and 
     *  pass the 5 points as generic attributes.
     */
    const char* vshader_pyramid_source =
        "#version 150 compatibility                                         \n"
        "#extension GL_ARB_explicit_attrib_location : enable                \n"
        " const int nb_vertices = 5;                                        \n"
        " layout(location=0) in vec3 p_in[nb_vertices];                     \n"
        " out Data {                                                        \n"
        "     vec3 p[nb_vertices];                                          \n"
        "     bool discardme;                                               \n"
        " } DataOut ;                                                       \n"
        " void main(void) {                                                 \n"
        "   DataOut.discardme = false;                                      \n"
        "   for(int i=0; i<nb_vertices; ++i) {                              \n"
        "       DataOut.p[i] = p_in[i];                                     \n"
        "   }                                                               \n"
        " }                                                                 \n";
    
    /**
     * \brief The geometry shader for pyramids.
     * \details Uses vshader_pyramid and gshader_utils 
     */
    const char* gshader_pyramid_source = 
        "layout(points) in;                                                 \n"
        "layout(triangle_strip, max_vertices = 28) out;                     \n"
        " in Data {                                                         \n"
        "     vec3 p[5];                                                    \n"
        "     bool discardme;                                               \n"
        " } DataIn[];                                                       \n"
        " uniform float shrink = 0.0;                                       \n"
        "                                                                   \n"
        "void draw_pyramid(                                                 \n"
        "   vec4 p0, vec4 p1, vec4 p2, vec4 p3, vec4 p4                     \n"
        ") {                                                                \n"
        "    if(clipping) {                                                 \n"
        "       if(                                                         \n"
        "           clip(p0,true) < 0.0 &&  clip(p1,true) < 0.0 &&          \n"
        "           clip(p2,true) < 0.0 &&  clip(p3,true) < 0.0 &&          \n"
        "           clip(p4,true) < 0.0                                     \n"
        "       ) {                                                         \n"
        "            return ;                                               \n"
        "       }                                                           \n"
        "    }                                                              \n"
        "   if(shrink != 0.0) {                                             \n"
        "      vec4 g = (1.0/5.0)*(p0+p1+p2+p3+p4);                         \n"
        "      p0 = shrink*g + (1.0-shrink)*p0;                             \n"
        "      p1 = shrink*g + (1.0-shrink)*p1;                             \n"
        "      p2 = shrink*g + (1.0-shrink)*p2;                             \n"
        "      p3 = shrink*g + (1.0-shrink)*p3;                             \n"
        "      p4 = shrink*g + (1.0-shrink)*p4;                             \n"
        "   }                                                               \n"
        "   vec4 pp0 = project(p0);                                         \n"
        "   vec4 pp1 = project(p1);                                         \n"
        "   vec4 pp2 = project(p2);                                         \n"
        "   vec4 pp3 = project(p3);                                         \n"
        "   vec4 pp4 = project(p4);                                         \n"
        "   flat_shaded_quad(p0,p1,p2,p3,pp0,pp1,pp2,pp3,false);            \n"
        "   flat_shaded_triangle(p0,p4,p1,pp0,pp4,pp1,false);               \n"
        "   flat_shaded_triangle(p0,p3,p4,pp0,pp3,pp4,false);               \n"
        "   flat_shaded_triangle(p2,p4,p3,pp2,pp4,pp3,false);               \n"
        "   flat_shaded_triangle(p2,p1,p4,pp2,pp1,pp4,false);               \n"
        "}                                                                  \n"
        "                                                                   \n"
        "void main() {                                                      \n"
        "   if(DataIn[0].discardme) { return; }                             \n"
        "   draw_pyramid(                                                   \n"
        "      vec4(DataIn[0].p[0],1.0),                                    \n"
        "      vec4(DataIn[0].p[1],1.0),                                    \n"
        "      vec4(DataIn[0].p[2],1.0),                                    \n"
        "      vec4(DataIn[0].p[3],1.0),                                    \n"
        "      vec4(DataIn[0].p[4],1.0)                                     \n"
        "   );                                                              \n"
        "}                                                                  \n";

    
    /**
     * \brief The geometry shader for triangles.
     * \details Uses vshader_pass_through and gshader_utils.
     */
    const char* gshader_tri_source =
        "layout(triangles) in;                                              \n"
        "layout(triangle_strip, max_vertices = 3) out;                      \n"
        "                                                                   \n"
        "void draw_triangle(vec4 p1, vec4 p2, vec4 p3) {                    \n"
        "   flat_shaded_triangle(                                           \n"
        "     p1,p2,p3,                                                     \n"
        "     project(p1), project(p2), project(p3),                        \n"
        "     true                                                          \n"
        "   );                                                              \n"
        "}                                                                  \n"
        "                                                                   \n"
        "void main() {                                                      \n"
        "    draw_triangle(                                                 \n"
        "       gl_in[0].gl_Position,                                       \n"
        "       gl_in[1].gl_Position,                                       \n"
        "       gl_in[2].gl_Position                                        \n"
        "    );                                                             \n"
        "}                                                                  \n";


    /**
     * \brief The geometry shader for quads.
     * \details Uses v_shader_pass_through and gshader_utils.
     */
    const char* gshader_quad_source =
        "layout(lines_adjacency) in;                                        \n"
        "layout(triangle_strip, max_vertices = 4) out;                      \n"
        "                                                                   \n"
        "void draw_quad(vec4 p1, vec4 p2, vec4 p3, vec4 p4) {               \n"
        "   flat_shaded_quad(                                               \n"
        "     p1,p2,p3,p4,                                                  \n"
        "     project(p1), project(p2), project(p3), project(p4),           \n"
        "     true                                                          \n"
        "   );                                                              \n"
        "}                                                                  \n"
        "                                                                   \n"
        "void main() {                                                      \n"
        "    draw_quad(                                                     \n"
        "       gl_in[0].gl_Position,                                       \n"
        "       gl_in[1].gl_Position,                                       \n"
        "       gl_in[2].gl_Position,                                       \n"
        "       gl_in[3].gl_Position                                        \n"
        "    );                                                             \n"
        "}                                                                  \n";
    
    /**
     * \brief The geometry shader for tetrahedra.
     * \details Uses v_shader_pass_through and gshader_utils.
     */
    const char* gshader_tet_source =
        "layout(lines_adjacency) in;                                        \n"
        "layout(triangle_strip, max_vertices = 12) out;                     \n"
        "uniform float shrink=0.0;                                          \n"
        "                                                                   \n"
        "void draw_tet(vec4 p0, vec4 p1, vec4 p2, vec4 p3) {                \n"
        "    if(clipping) {                                                 \n"
        "       if(                                                         \n"
        "           clip(p0,true) < 0.0 && clip(p1,true) < 0.0 &&           \n"
        "           clip(p2,true) < 0.0 && clip(p3,true) < 0.0              \n"
        "       ) {                                                         \n"
        "            return ;                                               \n"
        "       }                                                           \n"
        "    }                                                              \n"
        "   if(shrink != 0.0) {                                             \n"
        "      vec4 g = (1.0/4.0)*(p0+p1+p2+p3);                            \n"
        "      p0 = shrink*g + (1.0-shrink)*p0;                             \n"
        "      p1 = shrink*g + (1.0-shrink)*p1;                             \n"
        "      p2 = shrink*g + (1.0-shrink)*p2;                             \n"
        "      p3 = shrink*g + (1.0-shrink)*p3;                             \n"
        "   }                                                               \n"
        "   vec4 pp0 = project(p0);                                         \n"
        "   vec4 pp1 = project(p1);                                         \n"
        "   vec4 pp2 = project(p2);                                         \n"
        "   vec4 pp3 = project(p3);                                         \n"
        "   flat_shaded_triangle(p0,p1,p2,pp0,pp1,pp2,false);               \n"
        "   flat_shaded_triangle(p1,p0,p3,pp1,pp0,pp3,false);               \n"
        "   flat_shaded_triangle(p0,p2,p3,pp0,pp2,pp3,false);               \n"
        "   flat_shaded_triangle(p2,p1,p3,pp2,pp1,pp3,false);               \n"
        "}                                                                  \n"
        "                                                                   \n"
        "void main() {                                                      \n"
        "    draw_tet(                                                      \n"
        "       gl_in[0].gl_Position,                                       \n"
        "       gl_in[1].gl_Position,                                       \n"
        "       gl_in[2].gl_Position,                                       \n"
        "       gl_in[3].gl_Position                                        \n"
        "    );                                                             \n"
        "}                                                                  \n";


    /**
     * \brief The geometry shader for prisms
     * \details Uses vshader_pass_through and gshader_utils.
     */
    const char* gshader_prism_source =
        "layout(triangles_adjacency) in;                                    \n"
        "layout(triangle_strip, max_vertices = 18) out;                     \n"
        "uniform float shrink = 0.0;                                        \n"
        "                                                                   \n"
        "void draw_prism(                                                   \n"
        "   vec4 p0, vec4 p1, vec4 p2,                                      \n"
        "   vec4 p3, vec4 p4, vec4 p5                                       \n"
        ") {                                                                \n"
        "    if(clipping) {                                                 \n"
        "       if(                                                         \n"
        "           clip(p0,true) < 0.0 && clip(p1,true) < 0.0 &&           \n"
        "           clip(p2,true) < 0.0 && clip(p3,true) < 0.0 &&           \n"
        "           clip(p4,true) < 0.0 && clip(p5,true) < 0.0              \n"
        "       ) {                                                         \n"
        "            return ;                                               \n"
        "       }                                                           \n"
        "    }                                                              \n"
        "   if(shrink != 0.0) {                                             \n"
        "      vec4 g = (1.0/6.0)*(p0+p1+p2+p3+p4+p5);                      \n"
        "      p0 = shrink*g + (1.0-shrink)*p0;                             \n"
        "      p1 = shrink*g + (1.0-shrink)*p1;                             \n"
        "      p2 = shrink*g + (1.0-shrink)*p2;                             \n"
        "      p3 = shrink*g + (1.0-shrink)*p3;                             \n"
        "      p4 = shrink*g + (1.0-shrink)*p4;                             \n"
        "      p5 = shrink*g + (1.0-shrink)*p5;                             \n" 
        "   }                                                               \n"
        "   vec4 pp0 = project(p0);                                         \n"
        "   vec4 pp1 = project(p1);                                         \n"
        "   vec4 pp2 = project(p2);                                         \n"
        "   vec4 pp3 = project(p3);                                         \n"
        "   vec4 pp4 = project(p4);                                         \n"
        "   vec4 pp5 = project(p5);                                         \n"
        "   flat_shaded_triangle(p0,p1,p2,pp0,pp1,pp2,false);               \n"
        "   flat_shaded_triangle(p3,p5,p4,pp3,pp5,pp4,false);               \n"
        "   flat_shaded_quad(p0,p3,p4,p1,pp0,pp3,pp4,pp1,false);            \n"
        "   flat_shaded_quad(p0,p2,p5,p3,pp0,pp2,pp5,pp3,false);            \n"
        "   flat_shaded_quad(p1,p4,p5,p2,pp1,pp4,pp5,pp2,false);            \n"
        "}                                                                  \n"
        "                                                                   \n"
        "void main() {                                                      \n"
        "    draw_prism(                                                    \n"
        "       gl_in[0].gl_Position,                                       \n"
        "       gl_in[1].gl_Position,                                       \n"
        "       gl_in[2].gl_Position,                                       \n"
        "       gl_in[3].gl_Position,                                       \n"
        "       gl_in[4].gl_Position,                                       \n"
        "       gl_in[5].gl_Position                                        \n"
        "    );                                                             \n"
        "}                                                                  \n";


    /**
     * \brief Sends a vertex of a mesh to OpenGL
     * \param[in] M a const reference to the mesh
     * \param[in] v the index of the vertex in \p M
     */
    inline void glMeshVertex(const Mesh* M, index_t v) {
        if(M->vertices.single_precision()) {
            glVertex3fv(
                M->vertices.single_precision_point_ptr(v)
            );
        } else {
            glVertex3dv(
                M->vertices.point_ptr(v)
            );
        }
    }

    /**
     * \brief Sends a triangle normal to OpenGL
     * \param[in] pi,pj,pk pointers to the single-precision coordinates
     *  of the free vertices.
     */
    inline void glTriangleNormal(
        const float* pi, const float* pj, const float* pk
    ) {
        float N[3];
        N[0] = -(pi[1]-pj[1])*(pk[2]-pj[2]) + (pi[2]-pj[2])*(pk[1]-pj[1]);
        N[1] = -(pi[2]-pj[2])*(pk[0]-pj[0]) + (pi[0]-pj[0])*(pk[2]-pj[2]);
        N[2] = -(pi[0]-pj[0])*(pk[1]-pj[1]) + (pi[1]-pj[1])*(pk[0]-pj[0]);
        glNormal3fv(N);
    }

    /**
     * \brief Sends a triangle normal to OpenGL
     * \param[in] pi,pj,pk pointers to the double-precision coordinates
     *  of the free vertices.
     */
    inline void glTriangleNormal(
        const double* pi, const double* pj, const double* pk
    ) {
        double N[3];
        N[0] = -(pi[1]-pj[1])*(pk[2]-pj[2]) + (pi[2]-pj[2])*(pk[1]-pj[1]);
        N[1] = -(pi[2]-pj[2])*(pk[0]-pj[0]) + (pi[0]-pj[0])*(pk[2]-pj[2]);
        N[2] = -(pi[0]-pj[0])*(pk[1]-pj[1]) + (pi[1]-pj[1])*(pk[0]-pj[0]);
        glNormal3dv(N);
    }
    
    
    /**
     * \brief Sends a mesh triangle normal to OpenGL
     * \param[in] M a const reference to the mesh
     * \param[in] i,j,k indices of the vertices of the triangle
     */
    inline void glMeshTriangleNormal(
        const Mesh* M, index_t i, index_t j, index_t k
    ) {
        if(M->vertices.single_precision()) {
            const float* pi = M->vertices.single_precision_point_ptr(i);
            const float* pj = M->vertices.single_precision_point_ptr(j);
            const float* pk = M->vertices.single_precision_point_ptr(k);
            glTriangleNormal(pi,pj,pk);
        } else {
            const double* pi = M->vertices.point_ptr(i);
            const double* pj = M->vertices.point_ptr(j);
            const double* pk = M->vertices.point_ptr(k);
            glTriangleNormal(pi,pj,pk);
        }
    }

    /**
     * \brief Sends a mesh surface facet normal to OpenGL
     * \param[in] M a const reference to the mesh
     * \param[in] f index of the facet 
     */
    inline void glMeshFacetNormal(const Mesh* M, index_t f) {
        if(M->vertices.single_precision()) {
            float N[3];
            N[0] = 0.0f;
            N[1] = 0.0f;
            N[2] = 0.0f;
            index_t c1 = M->facets.corners_begin(f);
            const float* p1 = M->vertices.single_precision_point_ptr(
                M->facet_corners.vertex(c1)
            );
            const float* p2 = M->vertices.single_precision_point_ptr(
                M->facet_corners.vertex(c1+1)
            );
            for(index_t c2 = c1+1; c2+1<M->facets.corners_end(f); ++c2) {
                const float* p3 = M->vertices.single_precision_point_ptr(
                    M->facet_corners.vertex(c2+1)                    
                );
                N[0] = N[0]
                    -(p1[1]-p2[1])*(p3[2]-p2[2]) + (p1[2]-p2[2])*(p3[1]-p2[1]);
                N[1] = N[1]
                    -(p1[2]-p2[2])*(p3[0]-p2[0]) + (p1[0]-p2[0])*(p3[2]-p2[2]);
                N[2] = N[2]
                    -(p1[0]-p2[0])*(p3[1]-p2[1]) + (p1[1]-p2[1])*(p3[0]-p2[0]);
                p2 = p3;
            }
            glNormal3fv(N);            
        } else {
            double N[3];
            N[0] = 0.0;
            N[1] = 0.0;
            N[2] = 0.0;
            index_t c1 = M->facets.corners_begin(f);
            const double* p1 = M->vertices.point_ptr(
                M->facet_corners.vertex(c1)
            );
            const double* p2 = M->vertices.point_ptr(
                M->facet_corners.vertex(c1+1)
            );
            for(index_t c2 = c1+1; c2+1<M->facets.corners_end(f); ++c2) {
                const double* p3 = M->vertices.point_ptr(
                    M->facet_corners.vertex(c2+1)                    
                );
                N[0] = N[0]
                    -(p1[1]-p2[1])*(p3[2]-p2[2]) + (p1[2]-p2[2])*(p3[1]-p2[1]);
                N[1] = N[1]
                    -(p1[2]-p2[2])*(p3[0]-p2[0]) + (p1[0]-p2[0])*(p3[2]-p2[2]);
                N[2] = N[2]
                    -(p1[0]-p2[0])*(p3[1]-p2[1]) + (p1[1]-p2[1])*(p3[0]-p2[0]);
                p2 = p3;
            }
            glNormal3dv(N);            
        }
    }

    /**
     * \brief Sends a mesh cell facet normal to OpenGL
     * \param[in] M a const reference to the mesh
     * \param[in] f index of the facet 
     */
    inline void glMeshCellFacetNormal(const Mesh* M, index_t c, index_t lf) {
        glMeshTriangleNormal(
            M,
            M->cells.facet_vertex(c,lf,0),
            M->cells.facet_vertex(c,lf,1),
            M->cells.facet_vertex(c,lf,2)
        );
    }

    /**
     * \brief Tests whether current OpenGL polygon mode is filled.
     * \retval true if OpenGL polygon mode is filled
     * \retval false if OpenGL polygon mode is wireframe
     */
    inline bool glFillsPolygons() {
        GLint polymode[2];
        glGetIntegerv(GL_POLYGON_MODE, polymode);
        return (polymode[0] == GL_FILL);
    }

}

namespace GEO {

    MeshGfx::MeshGfx() : mesh_(nil) {

        vertices_VBO_ = 0;
        facet_indices_VBO_ = 0;
        cell_indices_VBO_ = 0;
        
        show_mesh_ = true;
        shrink_ = 0.0;
        animate_ = false;
        time_ = 0.0;
        lighting_ = true;
        
        triangles_and_quads_ = false;
        GLSL_mode_ = true;
        GLSL_tesselation_ = false;

        for(index_t i=0; i<PRG_NB; ++i) {
            programs_[i] = 0;
            set_color(i, 0.9f, 0.9f, 0.9f);
            set_back_color(i, 0.9f, 0.9f, 0.9f);            
        }

        set_mesh_color(0.0f, 0.0f, 0.0f);
        
        shaders_init_ = false;

        triangles_and_quads_ = false;
        

        for(index_t i=0; i<MESH_NB_CELL_TYPES; ++i) {
            cell_draw_cached_[i] = false;
        }

        for(index_t i=0; i<MESH_NB_CELL_TYPES; ++i) {
            draw_cells_[i] = true;
        }
        
    }

    void MeshGfx::set_mesh(Mesh* M) {
        delete_VBOs();
        clear_cell_draw_cache();
        triangles_and_quads_ = true;
        for(index_t f=0; f<M->facets.nb(); ++f) {
            if(
                M->facets.nb_vertices(f) != 3 &&
                M->facets.nb_vertices(f) != 4
            ) {
                triangles_and_quads_ = false;
                break;
            }
        }
        mesh_ = M;
    }
    
    MeshGfx::~MeshGfx() {
        delete_VBOs();        
        delete_shaders();
    }
    

    void MeshGfx::setup_shaders() {
        if(shaders_init_) {
            return;
        }
        shaders_init_ = true;
        if(CmdLine::get_arg_bool("gfx:GLSL")) {
            Logger::out("GLSL")
                << "Trying to use OpenGL shaders. "
                << "To deactivate, specify use_GLSL=false on the command line"
                << std::endl;
        } else {
            Logger::out("GLSL") << "OpenGL shaders deactivated"
                                     << std::endl;
            GLSL_mode_ = false;
            return;
        }
        const char* shading_language_ver_str = (const char*)glGetString(
            GL_SHADING_LANGUAGE_VERSION
        );
        Logger::out("GLSL") << "version = "
            << shading_language_ver_str << std::endl;
        double shading_language_ver = atof(shading_language_ver_str);
        if(shading_language_ver < 1.5) {
            Logger::out("GLSL")
                << "Deactivated (requires GLSL version >= 1.50)"
                << std::endl;
            GLSL_mode_ = false;
            return;
        } else {
            Logger::out("GLSL") << "Using GLSL shaders" << std::endl;
        }

        GLSL_tesselation_ =
            CmdLine::get_arg_bool("gfx:GLSL_tesselation") &&
            (shading_language_ver >= 4.0);

        if(GLSL_tesselation_) {
            Logger::out("GLSL")
                << "Using Tesselation shaders" << std::endl;
        } else {
            Logger::out("GLSL")
                << "Deactivated Tesselation shaders" << std::endl;
            if(CmdLine::get_arg_bool("gfx:GLSL_tesselation")) {
                Logger::out("GLSL")
                    << "  (unsupported, requires GLSL language ver. >= 4.0)"
                    << std::endl;
            } else {
                Logger::out("GLSL")
                    << "  (specify use_tesselation=true "
                    << "on command line to activate)"
                    << std::endl;
            }
        }
        
        GLuint vshader_pass_through = GLSL::compile_shader(
            GL_VERTEX_SHADER, vshader_pass_through_source
        );
        GLuint gshader_tri = GLSL::compile_shader(
            GL_GEOMETRY_SHADER, gshader_utils_source, gshader_tri_source
        );
        GLuint gshader_tet = GLSL::compile_shader(
            GL_GEOMETRY_SHADER, gshader_utils_source, gshader_tet_source
        );
        GLuint gshader_quad = GLSL::compile_shader(
            GL_GEOMETRY_SHADER, gshader_utils_source, gshader_quad_source
        );
        GLuint gshader_prism = GLSL::compile_shader(
            GL_GEOMETRY_SHADER, gshader_utils_source, gshader_prism_source
        );
        GLuint fshader = GLSL::compile_shader(
            GL_FRAGMENT_SHADER, fshader_source
        );
        GLuint points_fshader = GLSL::compile_shader(
            GL_FRAGMENT_SHADER, points_fshader_source
        );

        GLuint gshader_hex = GLSL::compile_shader(
            GL_GEOMETRY_SHADER, gshader_utils_source, gshader_hex_source
        );
        GLuint gshader_pyramid = GLSL::compile_shader(
            GL_GEOMETRY_SHADER, gshader_utils_source, gshader_pyramid_source
        );
            
        programs_[PRG_POINTS] = GLSL::setup_program(points_fshader, 0);
        programs_[PRG_TRI] = GLSL::setup_program(
            vshader_pass_through, gshader_tri, fshader, 0
        );
        programs_[PRG_QUAD] = GLSL::setup_program(
            vshader_pass_through, gshader_quad, fshader, 0
        );
        programs_[PRG_TET] = GLSL::setup_program(
            vshader_pass_through, gshader_tet, fshader, 0
        );
        programs_[PRG_PRISM] = GLSL::setup_program(
            vshader_pass_through, gshader_prism, fshader, 0
        );
        
        if(GLSL_tesselation_) {
            GLuint teshader_hex = GLSL::compile_shader(
                GL_TESS_EVALUATION_SHADER, teshader_hex_source
            );
            GLuint teshader_pyramid = GLSL::compile_shader(
                GL_TESS_EVALUATION_SHADER, teshader_pyramid_source
            );
            programs_[PRG_HEX] = GLSL::setup_program(
                vshader_pass_through, teshader_hex, gshader_hex, fshader, 0
            );
            programs_[PRG_PYRAMID] = GLSL::setup_program(
                vshader_pass_through, teshader_pyramid,
                gshader_pyramid, fshader, 0
            );                        
        } else {
            GLuint vshader_hex = GLSL::compile_shader(
                GL_VERTEX_SHADER, vshader_hex_source
            );            
            GLuint vshader_pyramid = GLSL::compile_shader(
                GL_VERTEX_SHADER, vshader_pyramid_source
            );
            programs_[PRG_HEX] = GLSL::setup_program(
                vshader_hex, gshader_hex, fshader, 0
            );
            programs_[PRG_PYRAMID] = GLSL::setup_program(
                vshader_pyramid, gshader_pyramid, fshader, 0
            );            
        }
    }

    void MeshGfx::delete_shaders() {
        for(index_t i=0; i<PRG_NB; ++i) {
            if(programs_[i] != 0) {
                glDeleteProgram(programs_[i]);
                programs_[i] = 0;
            }
        }
    }
    
    void MeshGfx::setup_VBOs() {
        if(vertices_VBO_ == 0 && mesh_->vertices.nb() != 0) {
            glGenBuffers(1, &vertices_VBO_);
            glBindBuffer(GL_ARRAY_BUFFER, vertices_VBO_);
            if(mesh_->vertices.single_precision()) {
                glBufferData(
                    GL_ARRAY_BUFFER,
                    GLsizeiptr(
                        mesh_->vertices.nb() *
                        mesh_->vertices.dimension() * sizeof(float)
                    ),
                    mesh_->vertices.single_precision_point_ptr(0),
                    GL_STATIC_DRAW
                );
                glEnableClientState(GL_VERTEX_ARRAY);
                unsigned int stride =
                    (unsigned int) (
                        mesh_->vertices.dimension() * sizeof(float)
                    );
                geo_assert(mesh_->vertices.dimension() >= 3);
                glVertexPointer(3, GL_FLOAT, GLsizei(stride), 0);
            } else {
                glBufferData(
                    GL_ARRAY_BUFFER,
                    GLsizeiptr(
                        mesh_->vertices.nb() *
                        mesh_->vertices.dimension() * sizeof(double)
                    ),
                    mesh_->vertices.point_ptr(0),
                    GL_STATIC_DRAW
                );
                glEnableClientState(GL_VERTEX_ARRAY);
                unsigned int stride =
                    (unsigned int) (
                        mesh_->vertices.dimension() * sizeof(double)
                    );
                geo_assert(mesh_->vertices.dimension() >= 3);
                glVertexPointer(3, GL_DOUBLE, GLsizei(stride), 0);
            }
        }
        if(facet_indices_VBO_ == 0 && mesh_->facets.nb() != 0) {
            glGenBuffers(1, &facet_indices_VBO_);
            glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, facet_indices_VBO_);
            glBufferData(
                GL_ELEMENT_ARRAY_BUFFER,
                GLsizeiptr(mesh_->facet_corners.nb() * sizeof(int)),
                mesh_->facet_corners.vertex_index_ptr(0), GL_STATIC_DRAW
            );
        }
        if(cell_indices_VBO_ == 0 && mesh_->cells.nb() != 0) {
            glGenBuffers(1, &cell_indices_VBO_);
            glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, cell_indices_VBO_);
            glBufferData(
                GL_ELEMENT_ARRAY_BUFFER,
                GLsizeiptr(mesh_->cell_corners.nb() * sizeof(int)),
                mesh_->cell_corners.vertex_index_ptr(0), GL_STATIC_DRAW
            );
        }
    }
    
    void MeshGfx::delete_VBOs() {
        if(vertices_VBO_ != 0) {
            glDeleteBuffers(1, &vertices_VBO_);
            vertices_VBO_ = 0;
        }
        if(facet_indices_VBO_ != 0) {
            glDeleteBuffers(1, &facet_indices_VBO_);
            facet_indices_VBO_ = 0;
        }
        if(cell_indices_VBO_ != 0) {
            glDeleteBuffers(1, &cell_indices_VBO_);
            cell_indices_VBO_ = 0;
        }
    }

    void MeshGfx::begin_draw(MeshElementsFlags what) {
        setup_shaders();
        setup_VBOs();

        glBindBuffer(GL_ARRAY_BUFFER, vertices_VBO_);
            
        switch(what) {
        case MESH_VERTICES: {
            // Nothing else to bind
        } break;
        case MESH_FACETS: {
            glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, facet_indices_VBO_);
        } break;
        case MESH_CELLS: {
            glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, cell_indices_VBO_);
        } break;
        default:
            geo_assert_not_reached;
        }
    }

    void MeshGfx::end_draw() {
        glBindBuffer(GL_ARRAY_BUFFER, 0);
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);        
    }

    void MeshGfx::set_colors(ShaderName name) {
        geo_assert(name < PRG_NB);        
        glMaterialfv(GL_FRONT, GL_DIFFUSE, colors_[name]);        
        glMaterialfv(GL_BACK, GL_DIFFUSE, back_colors_[name]);
        glColor3f(colors_[name][0], colors_[name][1], colors_[name][2]);
    }
    
    void MeshGfx::begin_shader(ShaderName name) {
        geo_assert(name < PRG_NB);
        setup_shaders();
        if(GLSL_mode_) {
            glUseProgram(programs_[name]);
            if(name == PRG_POINTS) {
            } else {
                GLint loc = glGetUniformLocation(programs_[name], "mesh");
                glUniform1i(loc, show_mesh_);
                loc = glGetUniformLocation(programs_[name], "lighting");
                glUniform1i(loc, lighting_);
                loc = glGetUniformLocation(programs_[name], "mesh_color");
                glUniform3f(
                    loc, mesh_color_[0], mesh_color_[1], mesh_color_[2]
                );
                loc = glGetUniformLocation(programs_[name], "shrink");
                glUniform1f(loc,float(shrink_));
            }
        }
        set_colors(name);
    }

    void MeshGfx::end_shader() {
        if(GLSL_mode_) {
            glUseProgram(0);
        }
    }
    
    void MeshGfx::draw_vertices() {
        begin_draw(MESH_VERTICES);
        
        glDisable(GL_LIGHTING);
        glPointSize(15);

        if(GLSL_mode_) {
            glEnable(GL_POINT_SPRITE);
            glTexEnvi(GL_POINT_SPRITE, GL_COORD_REPLACE, GL_TRUE);
            begin_shader(PRG_POINTS);            
        }
        
        glColor3f(0.0f, 1.0f, 0.0f);

        glDrawArrays(GL_POINTS, 0, GLsizei(mesh_->vertices.nb()));
        
        if(GLSL_mode_) {        
            end_shader();
        }

        end_draw();
    }

    void MeshGfx::draw_surface() {
        if(mesh_->facets.nb() == 0) {
            return;
        }

        if(
            show_mesh_ &&
            (!GLSL_mode_ ||
             (!mesh_->facets.are_simplices() && !triangles_and_quads_)
            )
        ) {
            glDisable(GL_LIGHTING);
            glLineWidth(1);
            glColor3f(mesh_color_[0], mesh_color_[1], mesh_color_[2]);
            glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);

            if(animate_) {
                draw_triangles_animation();                
            } else {
                if(mesh_->facets.are_simplices()) {
                    draw_triangles();
                } else {
                    draw_polygons();
                }
            }
            
            glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
        }

        if(animate_) {
            draw_triangles_animation();
        } else {
            if(GLSL_mode_ && mesh_->facets.are_simplices()) {
                draw_triangles();
            } else if(GLSL_mode_ && triangles_and_quads_) {
                draw_triangles_and_quads();
            } else {
                draw_polygons();
            }
        }
    }

    void MeshGfx::draw_surface_borders() {
        if(mesh_->facets.nb() == 0) {
            return;
        }
        glLineWidth(2);
        glColor3f(mesh_color_[0], mesh_color_[1], mesh_color_[2]);
        glBegin(GL_LINES);
        for(index_t f=0; f<mesh_->facets.nb(); ++f) {
            for(
                index_t c1=mesh_->facets.corners_begin(f);
                c1<mesh_->facets.corners_end(f); ++c1
            ) {
                if(mesh_->facet_corners.adjacent_facet(c1) == NO_FACET) {
                    index_t v1 = mesh_->facet_corners.vertex(c1);                    
                    index_t c2 = mesh_->facets.next_corner_around_facet(f,c1);
                    index_t v2 = mesh_->facet_corners.vertex(c2);
                    glMeshVertex(mesh_, v1);
                    glMeshVertex(mesh_, v2);                    
                }
            }
        }
        glEnd();

        glLineWidth(1);        
        Attribute<index_t> region;
        region.bind_if_is_defined(mesh_->facets.attributes(),"region");
        if(region.is_bound()) {
            glBegin(GL_LINES);            
            for(index_t f=0; f<mesh_->facets.nb(); ++f) {
                index_t f_rgn = region[f];
                for(
                    index_t c1=mesh_->facets.corners_begin(f);
                    c1<mesh_->facets.corners_end(f); ++c1
                ) {
                    index_t adj_f = mesh_->facet_corners.adjacent_facet(c1);
                    if(
                        adj_f != NO_FACET && region[adj_f] != f_rgn
                    ) {
                        index_t v1 = mesh_->facet_corners.vertex(c1);                    
                        index_t c2 = mesh_->facets.next_corner_around_facet(f,c1);
                        index_t v2 = mesh_->facet_corners.vertex(c2);
                        glMeshVertex(mesh_, v1);
                        glMeshVertex(mesh_, v2);                    
                    }
                }
            }
            glEnd();            
        }
    }

    
    void MeshGfx::draw_triangles() {
        begin_draw(MESH_FACETS);
        if(glFillsPolygons()) {
            begin_shader(PRG_TRI);
        }
        // Note: the fourth argument (0) corresponds to the bound VBO.        
        glDrawElements(
            GL_TRIANGLES, GLsizei(mesh_->facet_corners.nb()), GL_UNSIGNED_INT, 0
        );
        if(glFillsPolygons()) {
            end_shader();
        }
        end_draw();
    }

    void MeshGfx::draw_triangles_and_quads() {
        begin_draw(MESH_FACETS);

        // TODO: we could probably group primitives of the same type into
        // a single call...
        
        begin_shader(PRG_TRI);
        for(unsigned int f = 0; f < mesh_->facets.nb(); ++f) {
            unsigned int b = mesh_->facets.corners_begin(f);
            unsigned int e = mesh_->facets.corners_end(f);
            // Note: the fourth argument (void*)(b*sizeof(int))
            // is relative to the address of the bound VBO.
            if(e-b == 3) {
                glDrawElements(
                    GL_TRIANGLES, GLsizei(e - b),
                    GL_UNSIGNED_INT, (void*) (b * sizeof(int))
                );
            }
        }
        end_shader();
        
        begin_shader(PRG_QUAD);
        for(unsigned int f = 0; f < mesh_->facets.nb(); f++) {
            unsigned int b = mesh_->facets.corners_begin(f);
            unsigned int e = mesh_->facets.corners_end(f);
            // Note: the fourth argument (void*)(b*sizeof(int))
            // is relative to the address of the bound VBO.
            if(e-b == 4) {
                glDrawElements(
                    GL_LINES_ADJACENCY, GLsizei(e - b),
                    GL_UNSIGNED_INT, (void*) (b * sizeof(int))
                );
            }
        }
        end_shader();

        end_draw();
    }

    void MeshGfx::draw_polygons() {
        bool filled_polygons = glFillsPolygons();
        if(filled_polygons) {
            glEnable(GL_LIGHTING);
            glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);
            set_colors(PRG_TRI);
        }
        begin_draw(MESH_FACETS);
        for(unsigned int f = 0; f < mesh_->facets.nb(); ++f) {
            unsigned int b = mesh_->facets.corners_begin(f);
            unsigned int e = mesh_->facets.corners_end(f);
            // Note: the fourth argument (void*)(b*sizeof(int))
            // is relative to the address of the bound VBO.

            if(filled_polygons) {
                glMeshFacetNormal(mesh_, f);
            }
            
            glDrawElements(
                GL_POLYGON, GLsizei(e - b),
                GL_UNSIGNED_INT, (void*) (b * sizeof(int))
            );
        }
        end_draw();
        if(glFillsPolygons()) {
            glDisable(GL_LIGHTING);
        }
    }

    //   Note: if performance is an issue, we could write a special
    // shader for this one (and do the time interpolation in GLSL...)
    void MeshGfx::draw_triangles_animation() {
        geo_assert(mesh_->facets.are_simplices());

        if(GLSL_mode_) {
            begin_shader(PRG_TRI);
        } else {
            if(glFillsPolygons()) {
                if(lighting_) {
                    glEnable(GL_LIGHTING);
                    glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);
                }
                set_colors(PRG_TRI);
            }
        }

        glBegin(GL_TRIANGLES);

        for(GEO::index_t t=0; t < mesh_->facets.nb(); ++t) {

            if(mesh_->vertices.single_precision()) {
                float p[3][3];                
                for(GEO::index_t i=0; i<3; ++i) {
                    const float* p_t0 = mesh_->vertices.single_precision_point_ptr(mesh_->facets.vertex(t,i));
                    const float* p_t1 = p_t0+3;
                    for(GEO::coord_index_t c=0; c<3; ++c) {
                        p[i][c] = float(time_) * p_t1[c] + float(1.0 - time_) * p_t0[c];
                    }
                }
                if(!GLSL_mode_) {
                    glTriangleNormal(p[0],p[1],p[2]);
                }
                glVertex3fv(p[0]);
                glVertex3fv(p[1]);
                glVertex3fv(p[2]);
            } else {
                double p[3][3];                
                for(GEO::index_t i=0; i<3; ++i) {
                    const double* p_t0 = mesh_->vertices.point_ptr(mesh_->facets.vertex(t,i));
                    const double* p_t1 = p_t0+3;
                    for(GEO::coord_index_t c=0; c<3; ++c) {
                        p[i][c] = time_ * p_t1[c] + (1.0 - time_) * p_t0[c];
                    }
                }
                if(!GLSL_mode_) {
                    glTriangleNormal(p[0],p[1],p[2]);
                }
                glVertex3dv(p[0]);
                glVertex3dv(p[1]);
                glVertex3dv(p[2]);
            }
        }
        
        glEnd();


        if(GLSL_mode_) {
            end_shader();
        } else {
            if(lighting_ && glFillsPolygons()) {
                glDisable(GL_LIGHTING);
                glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_FALSE);
            }
        }
    }
    
    void MeshGfx::draw_volume() {
        if(mesh_->cells.nb() == 0) {
            return;
        }
        if(!GLSL_mode_ && !animate_) {
            glEnable(GL_LIGHTING);
            draw_cells_no_shader();
            glDisable(GL_LIGHTING);
        } else {
            glCullFace(GL_FRONT);
            glEnable(GL_CULL_FACE);
            if(animate_) {
                draw_tets_animation();
            } else {
                if(mesh_->cells.are_simplices()) {
                    draw_tets();
                } else {
                    draw_cells();                
                }
            }
            glDisable(GL_CULL_FACE);
        }
    }

    void MeshGfx::draw_tets() {

        if(!draw_cells_[MESH_TET]) {
            return;
        }

        begin_draw(MESH_CELLS);
        begin_shader(PRG_TET);

        // Note: the fourth argument (0) corresponds to the bound VBO.
        glDrawElements(
            GL_LINES_ADJACENCY,
            GLsizei(mesh_->cells.nb() * 4), GL_UNSIGNED_INT, 0
        );

        end_shader();
        end_draw();
    }

    void MeshGfx::draw_tets_animation() {
        geo_assert(mesh_->cells.are_simplices());
        
        if(GLSL_mode_) {

            begin_shader(PRG_TET);
            
            glBegin(GL_LINES_ADJACENCY);
            for(GEO::index_t t=0; t < mesh_->cells.nb(); ++t) {
                for(GEO::index_t i=0; i<4; ++i) {
                    if(mesh_->vertices.single_precision()) {
                        float p[3];
                        const float* p_t0 = mesh_->vertices.single_precision_point_ptr(
                            mesh_->cells.vertex(t,i)
                        );
                        const float* p_t1 = p_t0+3;
                        for(index_t c=0; c<3; ++c) {
                            p[c] = float(time_) * p_t1[c] + float(1.0 - time_) * p_t0[c];
                        }
                        glVertex3fv(p);
                    } else {
                        double p[3];
                        const double* p_t0 = mesh_->vertices.point_ptr(
                            mesh_->cells.vertex(t,i)
                        );
                        const double* p_t1 = p_t0+3;
                        for(index_t c=0; c<3; ++c) {
                            p[c] = time_ * p_t1[c] + (1.0 - time_) * p_t0[c];
                        }
                        glVertex3dv(p);
                    }
                }
            }
            glEnd();

            end_shader();
        }
    }
    
    
    void MeshGfx::draw_cells_no_shader() {
        static bool cell_type_has_triangles[MESH_NB_CELL_TYPES] = {
            true,  // tets
            false, // hexes
            true,  // prisms
            true,  // pyramids
            true   // connectors
        };

        static bool cell_type_has_quads[MESH_NB_CELL_TYPES] = {
            false, // tets
            true,  // hexes
            true,  // prisms
            true,  // pyramids
            true   // connectors
        };

        static ShaderName cell_type_to_prg[MESH_NB_CELL_TYPES] = {
            PRG_TET,      // tets
            PRG_HEX,      // hexes
            PRG_PRISM,    // prisms
            PRG_PYRAMID,  // pyramids
            PRG_HEX       // connectors (they do not have a prg, using hex)
        };
        
        glCullFace(GL_FRONT);
        glEnable(GL_CULL_FACE);
        
        // Try using loaded vertex array with glDrawElements() 

        set_colors(PRG_TET);
        glBegin(GL_TRIANGLES);
        for(index_t c=0; c<mesh_->cells.nb(); ++c) {
            MeshCellType type = mesh_->cells.type(c);
            if(!draw_cells_[type] || !cell_type_has_triangles[type]) {
                continue;
            }
            if(!mesh_->cells.are_simplices()) {
                set_colors(cell_type_to_prg[type]);
            }
            for(index_t f=0; f<mesh_->cells.nb_facets(c); ++f) {
                if(mesh_->cells.facet_nb_vertices(c,f) == 3) {
                    index_t i = mesh_->cells.facet_vertex(c,f,0);
                    index_t j = mesh_->cells.facet_vertex(c,f,1);
                    index_t k = mesh_->cells.facet_vertex(c,f,2);
                    glMeshTriangleNormal(mesh_,k,j,i);
                    glMeshVertex(mesh_,i);
                    glMeshVertex(mesh_,j);
                    glMeshVertex(mesh_,k);                    
                }
            }
        }
        glEnd();
        
        glBegin(GL_QUADS);
        for(index_t c=0; c<mesh_->cells.nb(); ++c) {
            MeshCellType type = mesh_->cells.type(c);            
            if(!draw_cells_[type] || !cell_type_has_quads[type]) {
                continue;
            }
            set_colors(cell_type_to_prg[type]);            
            for(index_t f=0; f<mesh_->cells.nb_facets(c); ++f) {
                if(mesh_->cells.facet_nb_vertices(c,f) == 4) {
                    index_t i = mesh_->cells.facet_vertex(c,f,0);
                    index_t j = mesh_->cells.facet_vertex(c,f,1);
                    index_t k = mesh_->cells.facet_vertex(c,f,2);
                    index_t l = mesh_->cells.facet_vertex(c,f,3);
                    glMeshTriangleNormal(mesh_,k,j,i);
                    glMeshVertex(mesh_,i);
                    glMeshVertex(mesh_,j);
                    glMeshVertex(mesh_,k);
                    glMeshVertex(mesh_,l);                    
                }
            }
        }
        glEnd();

        glDisable(GL_CULL_FACE);        
    }

    void MeshGfx::draw_cells() {
        begin_draw(MESH_CELLS);

        glCullFace(GL_FRONT);
        glEnable(GL_CULL_FACE);

        // Draw the tets, using GL_LINE_ADJACENCY (sends
        // 4 vertices per primitive)
        begin_shader(PRG_TET);
        draw_mesh_cells_as_opengl_elements(MESH_TET,GL_LINES_ADJACENCY);
        end_shader();

        // Draw the prisms, using GL_TRIANGLES_ADJACENCY (sends
        // 6 vertices per primitive)
        begin_shader(PRG_PRISM);
        draw_mesh_cells_as_opengl_elements(
            MESH_PRISM,GL_TRIANGLES_ADJACENCY
        );
        end_shader();
        
        // Draw the hexes and pyramids using a tesselation shader
        // to lookup the vertices. No standard OpenGL primitive has
        // 8 or 5 vertices, but GL_PATCH (used by tesselation shader)
        // has a configurable number of vertices !

        if(GLSL_tesselation_) {

            // The tesselation shader is just used to lookup
            // hex and pyramid vertices, and group them into
            // a single vertex passed to the geometry shader.
            // This is a way of emulating OpenGL primitives
            // with 8 and 5 vertices (hex and pyramids), since
            // GL_PATCHES has a configurable number of vertices.
            
            // We generate an isoline for each patch, with the
            // minimum tesselation level. This generates two
            // vertices (we discard one of them in the geometry
            // shader).
            static float levels[4] = {1.0, 1.0, 0.0, 0.0};
            glPatchParameterfv(GL_PATCH_DEFAULT_OUTER_LEVEL, levels);

            // Draw the hexes
            glPatchParameteri(GL_PATCH_VERTICES,8);
            begin_shader(PRG_HEX);
            draw_mesh_cells_as_opengl_elements(MESH_HEX, GL_PATCHES);
            end_shader();

            // Draw the pyramids
            glPatchParameteri(GL_PATCH_VERTICES,5);
            begin_shader(PRG_PYRAMID);
            draw_mesh_cells_as_opengl_elements(
                MESH_PYRAMID, GL_PATCHES
            );
            end_shader();
            
        } else {

            // If tesselation shaders are not available,
            // use GL_POINT primitives with 8 vec3's for hexes
            // and 5 vec3's for pyramids. It is probably much
            // much slower, for two reasons:
            // This uses one OpenGL call per point
            // The vertex puller, i.e. the indexing into
            // vertices array performed
            // by glDrawElements(), cannot be used.

            begin_shader(PRG_HEX);
            draw_mesh_cells_as_opengl_points(MESH_HEX);
            end_shader();
            
            begin_shader(PRG_PYRAMID);
            draw_mesh_cells_as_opengl_points(MESH_PYRAMID);
            end_shader();
        }

        glDisable(GL_CULL_FACE);        
        end_draw();
    }

    void MeshGfx::draw_mesh_cells_as_opengl_elements(
        MeshCellType cell_type, GLenum mode
    ) {
        if(!draw_cells_[cell_type]) {
            return;
        }

        // Determine the chunks of contiguous cell indices. Each chunk
        // can be drawn using a single glDrawElements() call.
        if(!cell_draw_cached_[cell_type]) {
            GLsizei vertices_per_cell = GLsizei(
                mesh_->cells.cell_type_to_cell_descriptor(
                    cell_type
                ).nb_vertices
            );
            for(index_t c=0; c<mesh_->cells.nb(); ++c) {
                index_t first_c = c;
                GLsizei nb_v = 0;
                while(mesh_->cells.type(c) == cell_type) {
                    nb_v += vertices_per_cell;
                    ++c;
                }
                if(nb_v != 0) {
                    cell_draw_nb_vertices_[cell_type].push_back(nb_v);
                    cell_draw_start_index_[cell_type].push_back(
                        (void*)(
                            mesh_->cells.corners_begin(first_c) * sizeof(int)
                        )
                    );
                }
            }
            cell_draw_cached_[cell_type] = true;
            Logger::out("GLSL cells")
                << "nb chunks for cell type "
                << cell_type << ":" << cell_draw_nb_vertices_[cell_type].size()
                << std::endl;
        }

        // Issue calls to glDrawElements() using the cached chunks.
        for(
            index_t chunk=0;
            chunk<cell_draw_nb_vertices_[cell_type].size(); ++chunk
        ) {
            GLsizei nb_v  = cell_draw_nb_vertices_[cell_type][chunk];
            void* start = cell_draw_start_index_[cell_type][chunk];
            glDrawElements(mode, nb_v, GL_UNSIGNED_INT, start);
        }
    }

    void MeshGfx::draw_mesh_cells_as_opengl_points(
        MeshCellType cell_type
    ) {
        
        if(!draw_cells_[cell_type]) {
            return;
        }
        
        const index_t nbv =
            mesh_->cells.cell_type_to_cell_descriptor(cell_type).nb_vertices;

        glBegin(GL_POINTS);        
        if(mesh_->vertices.single_precision()) {
            for(index_t c=0; c<mesh_->cells.nb(); ++c) {
                if(mesh_->cells.type(c) == cell_type) {
                    for(index_t i=1; i<nbv; ++i) {
                        index_t v = mesh_->cells.vertex(c,i);
                        glVertexAttrib3fv(
                            i, mesh_->vertices.single_precision_point_ptr(v)
                        );
                    }
                    // Vertex attrib 0 needs to come last, since it is the
                    // one that triggers vertex emit.
                    glVertexAttrib3fv(
                        0,
                        mesh_->vertices.single_precision_point_ptr(
                            mesh_->cells.vertex(c,0)
                        )
                    );
                }
            }
        } else {
            for(index_t c=0; c<mesh_->cells.nb(); ++c) {
                if(mesh_->cells.type(c) == cell_type) {
                    for(index_t i=1; i<nbv; ++i) {
                        index_t v = mesh_->cells.vertex(c,i);
                        glVertexAttrib3dv(i, mesh_->vertices.point_ptr(v));
                    }
                    // Vertex attrib 0 needs to come last, since it is the
                    // one that triggers vertex emit.
                    glVertexAttrib3dv(
                        0,mesh_->vertices.point_ptr(mesh_->cells.vertex(c,0))
                    );
                }
            }
        }
        glEnd();        
    }

    void MeshGfx::clear_cell_draw_cache() {
        for(index_t cell_type=0; cell_type < MESH_NB_CELL_TYPES; ++cell_type) {
            cell_draw_cached_[cell_type] = false;
            cell_draw_nb_vertices_[cell_type].clear();
            cell_draw_start_index_[cell_type].clear();
        }
        
    }
    
}
