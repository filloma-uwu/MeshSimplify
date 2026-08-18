#include "pqss_proxy_mesh/hausdorff_simplifier.hpp"

#include <cmath>
#include <cstddef>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>

namespace
{

using pqss_proxy_mesh::Position3;

Position3 subtract(const Position3& first, const Position3& second)
{
    return {first.x - second.x, first.y - second.y, first.z - second.z};
}

Position3 add(const Position3& first, const Position3& second)
{
    return {first.x + second.x, first.y + second.y, first.z + second.z};
}

Position3 scale(const Position3& point, const double factor)
{
    return {point.x * factor, point.y * factor, point.z * factor};
}

double dot(const Position3& first, const Position3& second)
{
    return first.x * second.x + first.y * second.y + first.z * second.z;
}

Position3 closestPointOnTriangle(
    const Position3& point, const Position3& first,
    const Position3& second, const Position3& third)
{
    const Position3 first_edge = subtract(second, first);
    const Position3 second_edge = subtract(third, first);
    const Position3 first_relative = subtract(point, first);
    const double first_projection = dot(first_edge, first_relative);
    const double second_projection = dot(second_edge, first_relative);
    if (first_projection <= 0.0 && second_projection <= 0.0) return first;

    const Position3 second_relative = subtract(point, second);
    const double third_projection = dot(first_edge, second_relative);
    const double fourth_projection = dot(second_edge, second_relative);
    if (third_projection >= 0.0 && fourth_projection <= third_projection) return second;

    const double first_region = first_projection * fourth_projection -
        third_projection * second_projection;
    if (first_region <= 0.0 && first_projection >= 0.0 && third_projection <= 0.0)
        return add(first, scale(first_edge,
            first_projection / (first_projection - third_projection)));

    const Position3 third_relative = subtract(point, third);
    const double fifth_projection = dot(first_edge, third_relative);
    const double sixth_projection = dot(second_edge, third_relative);
    if (sixth_projection >= 0.0 && fifth_projection <= sixth_projection) return third;

    const double second_region = fifth_projection * second_projection -
        first_projection * sixth_projection;
    if (second_region <= 0.0 && second_projection >= 0.0 && sixth_projection <= 0.0)
        return add(first, scale(second_edge,
            second_projection / (second_projection - sixth_projection)));

    const double opposite_region = third_projection * sixth_projection -
        fifth_projection * fourth_projection;
    if (opposite_region <= 0.0 && fourth_projection >= third_projection &&
        fifth_projection >= sixth_projection)
        return add(second, scale(subtract(third, second),
            (fourth_projection - third_projection) /
            ((fourth_projection - third_projection) +
             (fifth_projection - sixth_projection))));

    const double inverse = 1.0 / (opposite_region + second_region + first_region);
    return add(first, add(scale(first_edge, second_region * inverse),
        scale(second_edge, first_region * inverse)));
}

} // namespace

int main(int argc, char** argv)
{
    try
    {
        if (argc != 4)
            throw std::invalid_argument(
                "usage: pqss-inspect-halfedge-vertex phase1_halfedge.bin input.obj vertex_id");
        const auto halfedge = pqss_proxy_mesh::readAnalysisHalfedgeMesh(argv[1]);
        const auto source = pqss_proxy_mesh::readTriangleSoupObj(argv[2]);
        const std::size_t vertex = std::stoull(argv[3]);
        if (vertex >= halfedge.geometry.vertices.size())
            throw std::out_of_range("vertex ID exceeds halfedge vertex count");
        const Position3 point = halfedge.geometry.vertices[vertex];
        std::size_t nearest_face = std::numeric_limits<std::size_t>::max();
        Position3 nearest_point{};
        double nearest_squared = std::numeric_limits<double>::infinity();
        for (std::size_t face = 0; face < source.triangles.size(); ++face)
        {
            const auto triangle = source.triangles[face];
            const Position3 candidate = closestPointOnTriangle(
                point, source.vertices[triangle[0]], source.vertices[triangle[1]],
                source.vertices[triangle[2]]);
            const Position3 delta = subtract(candidate, point);
            const double squared = dot(delta, delta);
            if (squared < nearest_squared)
            {
                nearest_squared = squared;
                nearest_face = face;
                nearest_point = candidate;
            }
        }
        std::cout << std::setprecision(17)
                  << "vertex=" << vertex << '\n'
                  << "position=" << point.x << ',' << point.y << ',' << point.z << '\n'
                  << "nearest_source_face=" << nearest_face << '\n'
                  << "nearest_source_point=" << nearest_point.x << ','
                  << nearest_point.y << ',' << nearest_point.z << '\n'
                  << "nearest_distance=" << std::sqrt(nearest_squared) << '\n';
        const auto triangle = source.triangles[nearest_face];
        for (std::size_t corner = 0; corner < triangle.size(); ++corner)
        {
            const Position3 source_vertex = source.vertices[triangle[corner]];
            std::cout << "nearest_source_triangle_vertex_" << corner << '='
                      << triangle[corner] << ':' << source_vertex.x << ','
                      << source_vertex.y << ',' << source_vertex.z << '\n';
        }
        return 0;
    }
    catch (const std::exception& error)
    {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
