#pragma once

#include "pqss_proxy_mesh/options.hpp"

#include <array>
#include <cstdint>
#include <optional>
#include <string>
#include <vector>

namespace pqss_proxy_mesh
{

struct Position3
{
    double x = 0.0;
    double y = 0.0;
    double z = 0.0;
};

using TriangleIndices = std::array<std::uint32_t, 3>;

struct MeshModel
{
    std::string name;
    std::vector<Position3> vertices;
    std::vector<TriangleIndices> triangles;
};

using MeshModelPool = std::vector<MeshModel>;

struct SimplificationRequest
{
    MeshModelPool source_models;
    SimplificationOptions options;
};

struct SimplificationResult
{
    // proxy_models[i] corresponds to source_models[i].
    MeshModelPool proxy_models;
    ProxyPoolMetrics metrics;
};

[[nodiscard]] std::optional<std::string> validateModelPool(
    const MeshModelPool& models);

} // namespace pqss_proxy_mesh
