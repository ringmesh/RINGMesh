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

namespace RINGMesh
{
    class StratColBuilder{
    ringmesh_disable_copy(StratColBuilder);
    public:
        StratColBuilder(StratigraphicColumn& column);
        virtual ~StratColBuilder()
        {
        }
    protected:
        StratigraphicColumn& column_;
    };

    class StratColBuilderFile: public StratColBuilder{
    public:
        StratColBuilderFile(StratigraphicColumn& column, const std::string& filename);
        virtual ~StratColBuilderFile()
        {
        }
        void build_column()
        {
            load_file();
        }
    private:
        virtual void load_file() = 0;

    protected:
        std::string filename_;
    };

    class StratColBuilderXML: public StratColBuilderFile{
    public:
        StratColBuilderXML(StratigraphicColumn& column, const std::string& filename)
            : StratColBuilderFile(column, filename)
        {
        }
        virtual ~StratColBuilderXML()
        {
        }

    private:
        void load_file();
        void read_file();
        virtual void read_line();
    };
}


#endif /* INCLUDE_RINGMESH_GEOMODEL_STRATIGRAPHIC_COLUMN_BUILDER_H_ */
