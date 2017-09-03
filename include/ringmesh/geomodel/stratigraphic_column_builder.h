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
 * @file ringmesh/geomodel/stratigraphic_column_builder.h
 * @author Marie Sirvent, Pierre Anquez and Francois Bonneau
 */

#pragma once

#include <ringmesh/basic/common.h>

#include <ringmesh/geomodel/geomodel.h>
#include <ringmesh/geomodel/stratigraphic_column.h>

namespace RINGMesh {
    class RINGMESH_API StratigraphicColumnBuilder {
    ringmesh_disable_copy_and_move(StratigraphicColumnBuilder);
    public:
        StratigraphicColumnBuilder(
            StratigraphicColumn& column,
            GeoModel3D& model );
        virtual ~StratigraphicColumnBuilder() = default;
    protected:
        StratigraphicColumn& column_;
        GeoModel3D& model_;
    };

    class RINGMESH_API StratigraphicColumnBuilderFile: public StratigraphicColumnBuilder {
    public:
        StratigraphicColumnBuilderFile(
            StratigraphicColumn& column,
            GeoModel3D& model,
            std::string filename );
        void build_column()
        {
            load_file();
        }
    private:
        virtual void load_file() = 0;

    protected:
        std::string filename_;
    };

    class RINGMESH_API StratigraphicColumnBuilderXML: public StratigraphicColumnBuilderFile {
    public:
        StratigraphicColumnBuilderXML(
            StratigraphicColumn& column,
            GeoModel3D& model,
            const std::string& filename )
            : StratigraphicColumnBuilderFile( column, model, filename )
        {
        }
        void load_file() override;
    };
} // namespace RINGMesh
