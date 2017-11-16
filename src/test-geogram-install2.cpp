#include <geogram/mesh/mesh.h>

void toto_mesh()
{
    using namespace GEO;

    initialize();

    Mesh M;

    vec3 point{};
    M.vertices.create_vertex( point.data() );

    M.show_stats();
}
