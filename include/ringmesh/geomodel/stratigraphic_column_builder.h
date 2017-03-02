/*
 * stratigraphic_column_builder.h
 *
 *  Created on: Mar 1, 2017
 *      Author: sirvent1u
 */

#ifndef INCLUDE_RINGMESH_GEOMODEL_STRATIGRAPHIC_COLUMN_BUILDER_H_
#define INCLUDE_RINGMESH_GEOMODEL_STRATIGRAPHIC_COLUMN_BUILDER_H_

#include <ringmesh/basic/common.h>
#include <ringmesh/geomodel/stratigraphic_column.h>

namespace RINGMesh {
    class StratigraphicColumnBuilder {
    ringmesh_disable_copy(StratigraphicColumnBuilder) ;
    public:
        StratigraphicColumnBuilder( StratigraphicColumn& column ) ;
        virtual ~StratigraphicColumnBuilder()
        {
        }
    protected:
        StratigraphicColumn& column_ ;
    } ;

    class StratigraphicColumnBuilderFile: public StratigraphicColumnBuilder {
    public:
        StratigraphicColumnBuilderFile(
            StratigraphicColumn& column,
            const std::string& filename ) ;
        virtual ~StratigraphicColumnBuilderFile()
        {
        }
        void build_column()
        {
            load_file() ;
        }
    private:
        virtual void load_file() = 0 ;

    protected:
        std::string filename_ ;
    } ;

    class StratigraphicColumnBuilderXML: public StratigraphicColumnBuilderFile {
    public:
        StratigraphicColumnBuilderXML(
            StratigraphicColumn& column,
            const std::string& filename )
            : StratigraphicColumnBuilderFile( column, filename )
        {
        }
        virtual ~StratigraphicColumnBuilderXML()
        {
        }

    private:
        void load_file() ;
        void read_file() ;
        virtual void read_line() ;
    } ;
}

#endif /* INCLUDE_RINGMESH_GEOMODEL_STRATIGRAPHIC_COLUMN_BUILDER_H_ */
