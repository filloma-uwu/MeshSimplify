#include "pqss_proxy_mesh/topology_fill.hpp"

#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>

int main(int argc, char** argv)
{
    try
    {
        if (argc != 3)
            throw std::invalid_argument(
                "usage: pqss-canonicalize-soup input.obj output.obj");
        const auto source = pqss_proxy_mesh::readTriangleSoupObj(argv[1]);
        const auto result = pqss_proxy_mesh::canonicalizeCoplanarTriangleSoup(source);
        std::ofstream output(argv[2]);
        if (!output) throw std::runtime_error("failed to create output OBJ");
        output << std::setprecision(17) << "o coplanar_source_union\n";
        for (const auto point : result.vertices)
            output << "v " << point.x << ' ' << point.y << ' ' << point.z << '\n';
        for (const auto face : result.triangles)
            output << "f " << face[0] + 1 << ' ' << face[1] + 1 << ' '
                   << face[2] + 1 << '\n';
        std::cout << "input_triangles=" << source.triangles.size()
                  << " output_triangles=" << result.triangles.size() << '\n';
        return 0;
    }
    catch (const std::exception& exception)
    {
        std::cerr << exception.what() << '\n';
        return 1;
    }
}
