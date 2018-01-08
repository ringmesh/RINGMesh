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
 * @file Defintion of command line functions
 * @author Arnaud Botella
 */

#include <ringmesh/basic/command_line.h>

#include <geogram/basic/command_line.h>
#include <geogram/basic/command_line_args.h>

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
            GEO::CmdLine::declare_arg( "algo:tet", "TetGen",
                "Toggles the tetrahedral mesher (TetGen, MG_Tetra)" );
            GEO::CmdLine::declare_arg( "sys:plugins", "",
                "List of the plugins to load, separated by ;" );
        }

        void import_arg_group_in()
        {
            GEO::CmdLine::declare_arg_group( "in", "Input data" );
            GEO::CmdLine::declare_arg(
                "in:geomodel", "", "Filename of the input geological model" );
            GEO::CmdLine::declare_arg(
                "in:wells", "", "Filename of the input wells" );
        }

        void import_arg_group_out()
        {
            GEO::CmdLine::declare_arg_group( "out", "Output data" );
            GEO::CmdLine::declare_arg(
                "out:geomodel", "", "Saves the geological model" );
        }

        void import_arg_group_validity()
        {
            GEO::CmdLine::declare_arg_group( "validity", "Validity checks" );
            GEO::CmdLine::declare_arg( "validity:save", false,
                "Saves meshes representing geomodel inconsistencies",
                GEO::CmdLine::ARG_ADVANCED );
            GEO::CmdLine::declare_arg( "validity:directory", ".",
                "Directory to save meshes representing geomodel "
                "inconsistencies" );
            GEO::CmdLine::declare_arg( "validity:do_not_check", "0",
                "Toggle off checks at loading:\n"
                "By default all checks are toggled on."
                "'0' to toggle on all checks\n"
                "'A' to toggle off all checks\n"
                "'t' to toggle off all topology checks\n"
                "'g' to toggle off all geometry checks\n"
                "'G' to toggle off all geology checks\n"
                "  Topology checks:\n"
                "    'E' to toggle off checks on finite extension\n"
                "    'c' to toggle off checks on geomodel connectivity\n"
                "  Geometry checks:\n"
                "    's' to toggle off checks on conformity between surfaces "
                "and lines\n"
                "    'r' to toggle off checks on conformity between regions "
                "and surfaces\n"
                "    'm' to toggle off checks on mesh entities\n"
                "    'e' to toggle off checks on non manifold edges\n"
                "    'I' to toggle off checks on polygon intersections\n"
                "  Geology checks:\n"
                "    'f' to toggle off checks on geological entities",
                GEO::CmdLine::ARG_ADVANCED );
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
            else if( name == "validity" )
            {
                import_arg_group_validity();
            }
            else
            {
                return GEO::CmdLine::import_arg_group( name );
            }
            return true;
        }

    } // namespace CmdLine

} // namespace RINGMesh
