#include <iostream>
#include <geogram/mesh/mesh.h>

void toto_mesh()
{
    using namespace GEO;

    initialize();

    Mesh M;

    vec3 point{};
    M.vertices.create_vertex( point.data() );

    std::cout << MeshCellDescriptors::tet_descriptor.nb_facets << std::endl;

    M.show_stats();
}
