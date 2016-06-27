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
 
 /*!
 * @file Defintion of command line functions 
 * @author Arnaud Botella
 */

#include <ringmesh/command_line.h>

#include <geogram/basic/command_line.h>
#include <geogram/basic/command_line_args.h>

namespace RINGMesh {

    namespace CmdLine {

        void import_arg_group_attr()
        {
            GEO::CmdLine::declare_arg_group( "attr", "Attribute handler" ) ;
            GEO::CmdLine::declare_arg( "attr:colormap", "blue/white/red",
                "Colormap with colors separeted with /" ) ;
        }

        void import_arg_group_in()
        {
            GEO::CmdLine::declare_arg_group( "in", "Input data" ) ;
            GEO::CmdLine::declare_arg(
                "in:geomodel", "",
                "Filename of the input geological model" ) ;
            GEO::CmdLine::declare_arg(
                "in:intersection_check", true,
                "Toggle the surface intersection check at loading" ) ;
            GEO::CmdLine::declare_arg(
                "in:wells", "",
                "Filename of the input wells" ) ;
            GEO::CmdLine::declare_arg(
                "in:validity_save", true,
                "Saves meshes representing geomodel inconsistencies" ) ;
        }

        void import_arg_group_out()
        {
            GEO::CmdLine::declare_arg_group( "out", "Output data" ) ;
            GEO::CmdLine::declare_arg(
                "out:geomodel", "",
                "Saves the geological model" ) ;
        }

        void import_arg_group_stats()
        {
            GEO::CmdLine::declare_arg_group( "stats", "Statistics options" ) ;
            GEO::CmdLine::declare_arg(
                "stats:volume", false,
                "Print statistics on the volume" ) ;
            GEO::CmdLine::declare_arg(
                "stats:nb", true,
                "Print statistics on the number of entities" ) ;
        }

        void import_arg_group_remove_gme()
        {
            GEO::CmdLine::declare_arg_group( "remove_gme",
                "Options to remove a GeoModelElement and its dependencies" ) ;
            GEO::CmdLine::declare_arg(
                "remove_gme:type", "REGION",
                "Type of the GeoModelElement to remove" ) ;
            GEO::CmdLine::declare_arg(
                "remove_gme:id", 0,
                "Id of the element to remove (related to remove_gme:type)" ) ;
        }


        void import_arg_group_rotation()
        {
            GEO::CmdLine::declare_arg_group( "rotation",
                "Options to rotate a GeoModel" ) ;
            GEO::CmdLine::declare_arg( "rotation:origin", "0 0 0",
                "Origin of the rotation" ) ;
            GEO::CmdLine::declare_arg( "rotation:axis", "0 0 1",
                "Axis of rotation" ) ;
            GEO::CmdLine::declare_arg( "rotation:angle", 90., "Angle of rotation" ) ;
            GEO::CmdLine::declare_arg( "rotation:unit", "deg",
                "Angle unit (deg for degrees or rad for radians" ) ;
        }

        void import_arg_group_topology_info()
        {
            GEO::CmdLine::declare_arg_group( "topology_info",
                "Options to get topology information a GeoModel" ) ;
            GEO::CmdLine::declare_arg( "topology_info:interface_one_id", 0,
                "TODO" ) ;
        }

        void import_arg_group_fault_offset_info()
        {
            GEO::CmdLine::declare_arg_group( "fault_offset_info",
                "Options to get fault offset informations (throw, heave and dip)" ) ;
            GEO::CmdLine::declare_arg( "fault_offset_info:axis", "x", "TODO" ) ;
            GEO::CmdLine::declare_arg( "fault_offset_info:value", 0., "TODO" ) ;
            GEO::CmdLine::declare_arg( "fault_offset_info:tolerance", 1e-2, "TODO" ) ;
        }

        void import_arg_group_distance()
        {
            GEO::CmdLine::declare_arg_group( "distance",
                "Distance between 2 geomodels" ) ;
            GEO::CmdLine::declare_arg( "distance:geomodel", "",
                "Geomodel to compare to" ) ;
            GEO::CmdLine::declare_arg( "distance:sampling_distance", 10.,
                "Sampling distance (Hausdorff)" ) ;
        }

        bool import_arg_group( const std::string& name )
        {
            if( name == "in" ) {
                import_arg_group_in() ;
            } else if( name == "out" ) {
                import_arg_group_out() ;
            } else if( name == "stats" ) {
                import_arg_group_stats() ;
            } else if( name == "attr" ) {
                import_arg_group_attr() ;
            } else if( name == "remove_gme" ) {
                import_arg_group_remove_gme() ;
            } else if( name == "rotation" ) {
                import_arg_group_rotation() ;
            } else if( name == "topology_info" ) {
                import_arg_group_topology_info() ;
            } else if( name == "fault_offset_info" ) {
                import_arg_group_fault_offset_info() ;
            } else if( name == "distance" ) {
                import_arg_group_distance() ;
            } else {
                return GEO::CmdLine::import_arg_group( name ) ;
            }
            return true ;
        }

    }

}


