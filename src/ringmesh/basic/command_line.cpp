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
 * @file Defintion of command line functions
 * @author Arnaud Botella
 */

#include <ringmesh/basic/command_line.h>

#include <geogram/basic/command_line.h>
#include <geogram/basic/command_line_args.h>
#include <geogram/basic/file_system.h>

namespace RINGMesh
{
    namespace CmdLine
    {
        void import_arg_group_global()
        {
            GEO::CmdLine::declare_arg( "epsilon", 1e-7,
                "Threshold for numerical precision (ratio of the bbox "
                "diagonal)",
                GEO::CmdLine::ARG_ADVANCED );
            GEO::CmdLine::declare_arg( "validity_save", false,
                "Saves meshes representing geomodel inconsistencies",
                GEO::CmdLine::ARG_ADVANCED );
            GEO::CmdLine::declare_arg( "validity_directory",
                GEO::FileSystem::get_current_working_directory(),
                "Directory to save meshes representing geomodel "
                "inconsistencies",
                GEO::CmdLine::ARG_ADVANCED );
        }

        void import_arg_group_in()
        {
            GEO::CmdLine::declare_arg_group( "in", "Input data" );
            GEO::CmdLine::declare_arg(
                "in:geomodel", "", "Filename of the input geological model" );
            GEO::CmdLine::declare_arg( "in:intersection_check", true,
                "Toggle the surface intersection check at loading",
                GEO::CmdLine::ARG_ADVANCED );
            GEO::CmdLine::declare_arg(
                "in:wells", "", "Filename of the input wells" );
        }

        void import_arg_group_out()
        {
            GEO::CmdLine::declare_arg_group( "out", "Output data" );
            GEO::CmdLine::declare_arg(
                "out:geomodel", "", "Saves the geological model" );
        }

        bool import_arg_group( const std::string& name )
        {
            if( name == "global" )
            {
                import_arg_group_global();
            }
            else if( name == "in" )
            {
                import_arg_group_in();
            }
            else if( name == "out" )
            {
                import_arg_group_out();
            }
            else
            {
                return GEO::CmdLine::import_arg_group( name );
            }
            return true;
        }

    } // namespace CmdLine

} // namespace RINGMesh
