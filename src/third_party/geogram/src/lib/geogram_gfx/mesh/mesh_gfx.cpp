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
#include <geogram_gfx/mesh/mesh_gfx_private.h>
#include <geogram/basic/logger.h>
#include <geogram/basic/command_line.h>

namespace GEO {

    MeshGfx::MeshGfx() {
        initialized_ = false;
        // At construction, we create a plain
        // OpenGL implementation. It will be used
        // (at least) to store attributes/drawing modes etc...
        // if some are set by client code before calling
        // any draw() method and before OpenGL is initialized.
        // If supported, it is later replaced with a GLSL-enabled
        // implementation (see setup()), and attributes/drawing modes
        // are copied from this one to the new implementation.
        impl_ = new MeshGfxImplNoShader;
    }

    MeshGfx::~MeshGfx() {
        delete impl_;
        impl_ = nil;
    }

    void MeshGfx::draw_vertices() {
        setup();
        impl_->draw_vertices();
    }

    void MeshGfx::draw_edges() {
        setup();        
        impl_->draw_edges();
    }

    void MeshGfx::draw_surface() {
        setup();        
        impl_->draw_surface();
    }

    void MeshGfx::draw_surface_borders() {
        setup();        
        impl_->draw_surface_borders();
    }

    void MeshGfx::draw_volume() {
        setup();        
        impl_->draw_volume();
    }

    bool MeshGfx::get_show_mesh() const {
        return impl_->get_show_mesh();
    }

    void MeshGfx::set_show_mesh(bool x) {
        impl_->set_show_mesh(x);
    }

    index_t MeshGfx::get_mesh_width() const {
        return impl_->get_mesh_width();
    }

    void MeshGfx::set_mesh_width(index_t x) {
        impl_->set_mesh_width(x);
    }

    index_t MeshGfx::get_mesh_border_width() const {
        return impl_->get_mesh_border_width();
    }

    void MeshGfx::set_mesh_border_width(index_t x) {
        impl_->set_mesh_border_width(x);
    }

    double MeshGfx::get_shrink() const {
        return impl_->get_shrink();
    }

    void MeshGfx::set_shrink(double x) {
        impl_->set_shrink(x);
    }

    bool MeshGfx::get_animate() const {
        return impl_->get_animate();
    }

    void MeshGfx::set_animate(bool x) {
        impl_->set_animate(x);
    }

    double MeshGfx::get_time() const {
        return impl_->get_time();
    }

    void MeshGfx::set_time(double x) {
        impl_->set_time(x);
    }

    bool MeshGfx::get_draw_cells(MeshCellType type) const {
        return impl_->get_draw_cells(type);
    }

    void MeshGfx::set_draw_cells(MeshCellType type, bool x) {
        impl_->set_draw_cells(type,x);
    }

    void MeshGfx::set_points_color(float r, float g, float b) {
        impl_->set_points_color(r,g,b);
    }

    void MeshGfx::get_points_color(float& r, float& g, float& b) const {
        impl_->get_points_color(r,g,b);
    }

    void MeshGfx::set_points_size(float x) {
        impl_->set_points_size(x);
    }

    float MeshGfx::get_points_size() const {
        return impl_->get_points_size();
    }

    void MeshGfx::set_mesh_color(float r, float g, float b) {
        impl_->set_mesh_color(r,g,b);
    }

    void MeshGfx::get_mesh_color(float& r, float& g, float& b) const {
        impl_->get_mesh_color(r,g,b);
    }

    void MeshGfx::set_surface_color(float r, float g, float b) {
        impl_->set_surface_color(r,g,b);
    }

    void MeshGfx::get_surface_color(float& r, float& g, float& b) const {
        impl_->get_surface_color(r,g,b);        
    }

    void MeshGfx::set_backface_surface_color(float r, float g, float b) {
        impl_->set_backface_surface_color(r,g,b);        
    }

    void MeshGfx::set_cells_color(float r, float g, float b) {
        impl_->set_cells_color(r,g,b);        
    }

    void MeshGfx::get_cells_color(float& r, float& g, float& b) const {
        impl_->get_cells_color(r,g,b);                
    }

    void MeshGfx::set_cells_colors_by_type() {
        impl_->set_cells_colors_by_type();
    }

    bool MeshGfx::get_lighting() const {
        return impl_->get_lighting();
    }

    void MeshGfx::set_lighting(bool x) {
        impl_->set_lighting(x);
    }

    void MeshGfx::set_mesh(const Mesh* M) {
        impl_->set_mesh(M);
    }

    const Mesh* MeshGfx::mesh() const {
        return impl_->mesh();
    }

    void MeshGfx::set_picking_mode(MeshElementsFlags what) {
        impl_->set_picking_mode(what);
    }

    MeshElementsFlags MeshGfx::get_picking_mode() const {
        return impl_->get_picking_mode();
    }
    
    void MeshGfx::setup() {
        // First time, we try to "upgrade" implementation
        // to higher GLSL version if supported
        if(!initialized_) {
            double GLSL_version = supported_GLSL_version();
            Logger::out("GLSL")
                << "Selecting among NoShader,GLSL150,GLSL440..." << std::endl;
            if(GLSL_version < 1.5) {
                Logger::out("GLSL")
                    << "Using MeshGfxImplNoShader" << std::endl;
                // We keep the already bound default implementation
                //  (plain OpenGL)
            } else if(GLSL_version < 4.4) {
                Logger::out("GLSL") << "Using MeshGfxImplGLSL150" << std::endl;
                replace_implementation(new MeshGfxImplGLSL150); 
            } else {
                Logger::out("GLSL") << "Using MeshGfxImplGLSL440" << std::endl;
                replace_implementation(new MeshGfxImplGLSL440);
            }
            initialized_ = true;
        }
        
        try {
            impl_->setup();
        } catch(...) {
            Logger::warn("MeshGfx")
                << "Caught an exception, downgrading to plain OpenGL"
                << std::endl;
            // Back to the default implementation (plain OpenGL)
            replace_implementation(new MeshGfxImplNoShader);
            impl_->setup();
        }
    }

    void MeshGfx::replace_implementation(
        MeshGfxImpl* new_impl
    ) {
        new_impl->copy_drawing_attributes(*impl_);
        new_impl->set_mesh(impl_->mesh());
        delete impl_;
        impl_ = new_impl;
    }

    double MeshGfx::supported_GLSL_version() {
        const char* shading_language_ver_str = (const char*)glGetString(
            GL_SHADING_LANGUAGE_VERSION
        );
        const char* vendor = (const char*)glGetString(
            GL_VENDOR
        );
        Logger::out("GLSL") << "vendor = " << vendor << std::endl;
        Logger::out("GLSL") << "version string = "
            << shading_language_ver_str << std::endl;
        double GLSL_version = atof(shading_language_ver_str);
        Logger::out("GLSL") << "version = " << GLSL_version
                            << std::endl;
        if(!CmdLine::get_arg_bool("gfx:GLSL")) {
            Logger::out("GLSL") << "OpenGL shaders deactivated (gfx:GLSL=false)"
                                << std::endl;
            
            GLSL_version = 0.0;
        }
        double forced_version = CmdLine::get_arg_double("gfx:GLSL_version");
        if(forced_version != 0.0) {
            GLSL_version = forced_version;
            Logger::out("GLSL") << "forced to version "
                                << GLSL_version 
                                << " (gfx:GLSL_version)" << std::endl;
        }
        return GLSL_version;
    }
}

