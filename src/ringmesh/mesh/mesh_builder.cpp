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
 *     * Neither the name of ASGA nor the
 *       names of its contributors may be used to endorse or promote products
 *       derived from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL ASGA BE LIABLE FOR ANY
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

/*! \author Francois Bonneau */

#include <ringmesh/mesh/mesh_builder.h>

#include <ringmesh/mesh/geogram_mesh_builder.h>

namespace {
    using namespace RINGMesh ;

    Mesh0DBuilder* create_builder_0d( MeshBase& mesh )
    {
        Mesh0DBuilder* builder = Mesh0DBuilderFactory::create_object(
            mesh.type_name() ) ;
        if( builder ) {
            builder->set_mesh( dynamic_cast< Mesh0D& >( mesh ) ) ;
        }
        return builder ;
    }

    Mesh1DBuilder* create_builder_1d( MeshBase& mesh )
    {
        Mesh1DBuilder* builder = Mesh1DBuilderFactory::create_object(
            mesh.type_name() ) ;
        if( builder ) {
            builder->set_mesh( dynamic_cast< Mesh1D& >( mesh ) ) ;
        }
        return builder ;
    }

    Mesh2DBuilder* create_builder_2d( MeshBase& mesh )
    {
        Mesh2DBuilder* builder = Mesh2DBuilderFactory::create_object(
            mesh.type_name() ) ;
        if( builder ) {
            builder->set_mesh( dynamic_cast< Mesh2D& >( mesh ) ) ;
        }
        return builder ;
    }

    Mesh3DBuilder* create_builder_3d( MeshBase& mesh )
    {
        Mesh3DBuilder* builder = Mesh3DBuilderFactory::create_object(
            mesh.type_name() ) ;
        if( builder ) {
            builder->set_mesh( dynamic_cast< Mesh3D& >( mesh ) ) ;
        }
        return builder ;
    }

    MeshAllDBuilder* create_builder_alld( MeshBase& mesh )
    {
        MeshAllDBuilder* builder = MeshAllDBuilderFactory::create_object(
            mesh.type_name() ) ;
        if( builder ) {
            builder->set_mesh( dynamic_cast< MeshAllD& >( mesh ) ) ;
        }
        return builder ;
    }
}

namespace RINGMesh {

    MeshBaseBuilder* MeshBaseBuilder::create_builder( MeshBase& mesh )
    {
        MeshBaseBuilder* builder = MeshBaseBuilderFactory::create_object(
            mesh.type_name() ) ;
        if( !builder ) {
            builder = create_builder_alld( mesh ) ;
            if( !builder ) {
                builder = create_builder_0d( mesh ) ;
                if( !builder ) {
                    builder = create_builder_1d( mesh ) ;
                    if( !builder ) {
                        builder = create_builder_2d( mesh ) ;
                        if( !builder ) {
                            builder = create_builder_3d( mesh ) ;
                        }
                    }
                }
            }
        }
        if( !builder ) {
            Logger::warn( "MeshBaseBuilder" )
                << "Could not create mesh data structure: " << mesh.type_name()
                << std::endl ;
            Logger::warn( "MeshBaseBuilder" )
                << "Falling back to GeogramMesh0DBuilder data structure"
                << std::endl ;

            builder = new GeogramMeshBaseBuilder ;
        }
        builder->set_mesh( mesh ) ;
        return builder ;
    }

    Mesh0DBuilder* Mesh0DBuilder::create_builder( Mesh0D& mesh )
    {
        Mesh0DBuilder* builder = Mesh0DBuilderFactory::create_object(
            mesh.type_name() ) ;
        if( !builder ) {
            builder = create_builder_alld( mesh ) ;
        }
        if( !builder ) {
            Logger::warn( "Mesh0DBuilder" )
                << "Could not create mesh data structure: " << mesh.type_name()
                << std::endl ;
            Logger::warn( "Mesh0DBuilder" )
                << "Falling back to GeogramMesh0DBuilder data structure"
                << std::endl ;

            builder = new GeogramMesh0DBuilder ;
        }
        builder->set_mesh( mesh ) ;
        return builder ;
    }

	Mesh0D* Mesh0DBuilder::change_mesh_data_structure(Mesh0D* mesh, const MeshType type)
	{
		if (mesh->type_name() == type) {
			return mesh;
		}
		Mesh0D* new_mesh = Mesh0D::create_mesh(type);
		Mesh0DBuilder_var builder = Mesh0DBuilder::create_builder(
			*new_mesh);
		builder->copy(*mesh, true);
		delete mesh;
		return new_mesh;
	}

    Mesh1DBuilder* Mesh1DBuilder::create_builder( Mesh1D& mesh )
    {
        Mesh1DBuilder* builder = Mesh1DBuilderFactory::create_object(
            mesh.type_name() ) ;
        if( !builder ) {
            builder = create_builder_alld( mesh ) ;
        }
        if( !builder ) {
            Logger::warn( "Mesh1DBuilder" )
                << "Could not create mesh data structure: " << mesh.type_name()
                << std::endl ;
            Logger::warn( "Mesh1DBuilder" )
                << "Falling back to GeogramMesh1DBuilder data structure"
                << std::endl ;

            builder = new GeogramMesh1DBuilder ;
        }
        builder->set_mesh( mesh ) ;
        return builder ;
    }

	Mesh1D* Mesh1DBuilder::change_mesh_data_structure(Mesh1D* mesh, const MeshType type)
	{
		if (mesh->type_name() == type) {
			return mesh;
		}
		Mesh1D* new_mesh = Mesh1D::create_mesh(type);
		Mesh1DBuilder_var builder = Mesh1DBuilder::create_builder(
			*new_mesh);
		builder->copy(*mesh, true);
		delete mesh;
		return new_mesh;
	}

    Mesh2DBuilder* Mesh2DBuilder::create_builder( Mesh2D& mesh )
    {
        Mesh2DBuilder* builder = Mesh2DBuilderFactory::create_object(
            mesh.type_name() ) ;
        if( !builder ) {
            builder = create_builder_alld( mesh ) ;
        }
        if( !builder ) {
            Logger::warn( "Mesh2DBuilder" )
                << "Could not create mesh data structure: " << mesh.type_name()
                << std::endl ;
            Logger::warn( "Mesh2DBuilder" )
                << "Falling back to GeogramMesh2DBuilder data structure"
                << std::endl ;

            builder = new GeogramMesh2DBuilder ;
        }
        builder->set_mesh( mesh ) ;
        return builder ;
    }

	Mesh2D* Mesh2DBuilder::change_mesh_data_structure(Mesh2D* mesh, const MeshType type)
	{
		if (mesh->type_name() == type) {
			return mesh;
		}
		Mesh2D* new_mesh = Mesh2D::create_mesh(type);
		Mesh2DBuilder_var builder = Mesh2DBuilder::create_builder(
			*new_mesh);
		builder->copy(*mesh, true);
		delete mesh;
		return new_mesh;
	}

    Mesh3DBuilder* Mesh3DBuilder::create_builder( Mesh3D& mesh )
    {
        Mesh3DBuilder* builder = Mesh3DBuilderFactory::create_object(
            mesh.type_name() ) ;
        if( !builder ) {
            builder = create_builder_alld( mesh ) ;
        }
        if( !builder ) {
            Logger::warn( "Mesh3DBuilder" )
                << "Could not create mesh data structure: " << mesh.type_name()
                << std::endl ;
            Logger::warn( "Mesh3DBuilder" )
                << "Falling back to GeogramMesh3DBuilder data structure"
                << std::endl ;

            builder = new GeogramMesh3DBuilder ;
        }
        builder->set_mesh( mesh ) ;
        return builder ;
    }

	Mesh3D* Mesh3DBuilder::change_mesh_data_structure(Mesh3D* mesh, const MeshType type)
	{
		if (mesh->type_name() == type) {
			return mesh;
		}
		Mesh3D* new_mesh = Mesh3D::create_mesh(type);
		Mesh3DBuilder_var builder = Mesh3DBuilder::create_builder(
			*new_mesh);
		builder->copy(*mesh, true);
		delete mesh;
		return new_mesh;
	}

    MeshAllDBuilder* MeshAllDBuilder::create_builder( MeshAllD& mesh )
    {
        MeshAllDBuilder* builder = MeshAllDBuilderFactory::create_object(
            mesh.type_name() ) ;
        if( !builder ) {
            Logger::warn( "MeshAllDBuilder" )
                << "Could not create mesh data structure: " << mesh.type_name()
                << std::endl ;
            Logger::warn( "MeshAllDBuilder" )
                << "Falling back to GeogramMeshAllDBuilder data structure"
                << std::endl ;

            builder = new GeogramMeshAllDBuilder ;
        }
        builder->set_mesh( mesh ) ;
        return builder ;
    }

	MeshAllD* MeshAllDBuilder::change_mesh_data_structure(MeshAllD* mesh, const MeshType type)
	{
		if (mesh->type_name() == type) {
			return mesh;
		}
		MeshAllD* new_mesh = MeshAllD::create_mesh(type);
		MeshAllDBuilder_var builder = MeshAllDBuilder::create_builder(
			*new_mesh);
		builder->copy(*mesh, true);
		delete mesh;
		return new_mesh;
	}
} // namespace
