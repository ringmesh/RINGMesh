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

/*!
 * @file Initialization of the RINGMesh and geogram library on loading
 * @author Arnaud Botella
 */

#include <ringmesh/basic/common.h>

#include <geogram/basic/command_line.h>
#include <geogram/basic/command_line_args.h>
#include <geogram/basic/common.h>

#include <ringmesh/basic/command_line.h>
#include <ringmesh/basic/logger.h>
#include <ringmesh/basic/plugin_manager.h>

namespace
{
    void configure_geogram()
    {
        RINGMesh::Logger::instance()->unregister_all_clients();
        RINGMesh::Logger::instance()->register_client(
            new RINGMesh::ThreadSafeConsoleLogger );
        GEO::CmdLine::import_arg_group( "sys" );
#ifdef RINGMESH_DEBUG
        GEO::CmdLine::set_arg( "sys:assert", "abort" );
#endif
        GEO::CmdLine::set_arg( "sys:FPE", true );
        GEO::CmdLine::import_arg_group( "algo" );
        GEO::CmdLine::set_arg( "algo:predicates", "exact" );
        GEO::CmdLine::import_arg_group( "log" );
        GEO::CmdLine::set_arg( "sys:use_doubles", true );
    }

    void configure_ringmesh()
    {
        RINGMesh::CmdLine::import_arg_group( "global" );
        RINGMesh::CmdLine::import_arg_group( "validity" );
    }

    RINGMESH_PLUGIN_INITIALIZE(
        RINGMesh_basic,
        // Plugin initialization
        GEO::initialize();
        configure_geogram();
        configure_ringmesh();
    );

} // namespace

namespace RINGMesh
{
    /*!
     * This function configures geogram by setting some geogram options.
     * \pre This function should be call after GEO::initialize().
     */

    void print_header_information()
    {
        Logger::div( "RINGMesh" );
        Logger::out( "", "This project is developed by the RINGMesh",
            " developers team:" );
        Logger::out(
            "", "RINGMesh-dev <georessources-ringmesh-dev@univ-lorraine.fr> " );
        Logger::out( "", "You can have access to the full open-source ",
            "code through its GitHub repository: " );
        Logger::out( "", "https://github.com/ringmesh/RINGMesh" );
        Logger::out( "", "More information on this project and other ",
            "projects of the team: " );
        Logger::out( "", "http://www.ring-team.org" );
    }

    template < index_t DIMENSION >
    vecn< DIMENSION > initialize_vecn_coordinates( double value )
    {
        index_t nb_coords = DIMENSION;
        vecn< DIMENSION > vec;
        for( auto coord : range( nb_coords ) )
        {
            vec[coord] = value;
        }
        return vec;
    }

    template vec2 basic_api initialize_vecn_coordinates( double );
    template vec3 basic_api initialize_vecn_coordinates( double );
} // namespace RINGMesh
