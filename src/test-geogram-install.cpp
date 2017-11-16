#include <geogram/mesh/mesh.h>

int main()
{
    using namespace GEO;

    initialize();

    Mesh M;

    vec3 point{};
    M.vertices.create_vertex( point.data() );

    M.show_stats();


    return 0;
}
