#include "pqss_proxy_mesh/mesh_pool.hpp"
#include "pqss_proxy_mesh/options.hpp"

#include <iostream>

namespace
{

int check(const bool condition, const char* message)
{
    if (!condition)
    {
        std::cerr << message << '\n';
        return 1;
    }
    return 0;
}

} // namespace

int main()
{
    pqss_proxy_mesh::SimplificationOptions options;
    if (check(!pqss_proxy_mesh::validateOptions(options), "default options must be valid")) return 1;

    options.max_pqss_bvh_depth = 0;
    if (check(pqss_proxy_mesh::validateOptions(options).has_value(), "zero depth must be rejected")) return 1;

    options.max_pqss_bvh_depth = 8;
    options.relative_alpha = 0.0;
    if (check(pqss_proxy_mesh::validateOptions(options).has_value(), "zero alpha must be rejected")) return 1;

    options.relative_alpha = 0.02;
    options.relative_offset = -1.0;
    if (check(pqss_proxy_mesh::validateOptions(options).has_value(), "negative offset must be rejected")) return 1;

    pqss_proxy_mesh::MeshModelPool models;
    if (check(pqss_proxy_mesh::validateModelPool(models).has_value(), "empty pool must be rejected")) return 1;

    pqss_proxy_mesh::MeshModel model;
    model.name = "triangle";
    model.vertices = {{0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}};
    model.triangles = {{{0, 1, 2}}};
    models.push_back(model);
    if (check(!pqss_proxy_mesh::validateModelPool(models), "valid one-model pool must be accepted")) return 1;

    models.push_back(model);
    if (check(!pqss_proxy_mesh::validateModelPool(models), "valid multi-model pool must be accepted")) return 1;

    return 0;
}
