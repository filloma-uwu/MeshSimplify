#include "pqss_proxy_mesh/mesh_pool.hpp"

#include <cmath>
#include <limits>

namespace pqss_proxy_mesh
{

std::optional<std::string> validateModelPool(const MeshModelPool& models)
{
    if (models.empty())
    {
        return "source model pool must not be empty";
    }

    for (std::size_t model_index = 0; model_index < models.size(); ++model_index)
    {
        const MeshModel& model = models[model_index];
        const std::string prefix = "model " + std::to_string(model_index) + ": ";
        if (model.vertices.empty())
        {
            return prefix + "vertex list must not be empty";
        }
        if (model.triangles.empty())
        {
            return prefix + "triangle list must not be empty";
        }

        for (const Position3& vertex : model.vertices)
        {
            if (!std::isfinite(vertex.x) || !std::isfinite(vertex.y) || !std::isfinite(vertex.z))
            {
                return prefix + "vertices must have finite coordinates";
            }
        }

        for (const TriangleIndices& triangle : model.triangles)
        {
            for (const std::uint32_t vertex_index : triangle)
            {
                if (vertex_index >= model.vertices.size())
                {
                    return prefix + "triangle index is outside the vertex list";
                }
            }
        }
    }
    return std::nullopt;
}

} // namespace pqss_proxy_mesh
