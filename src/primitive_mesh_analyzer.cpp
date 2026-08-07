#include "pqss_proxy_mesh/primitive_mesh_analyzer.hpp"

#include "clipper2/clipper.h"
#include "clipper2/clipper.triangulation.h"
#include "QuickHull.hpp"

#include <algorithm>
#include <array>
#include <bit>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <numeric>
#include <numbers>
#include <optional>
#include <queue>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

namespace pqss_proxy_mesh
{
namespace
{

struct Vec2
{
    double values[2]{};
    Vec2() = default;
    Vec2(const double x, const double y) : values{x, y} {}
    double& x() { return values[0]; }
    double& y() { return values[1]; }
    double x() const { return values[0]; }
    double y() const { return values[1]; }
    double& operator[](const int index) { return values[index]; }
    double operator[](const int index) const { return values[index]; }
    Vec2 operator+(const Vec2& other) const { return {x() + other.x(), y() + other.y()}; }
    Vec2 operator-(const Vec2& other) const { return {x() - other.x(), y() - other.y()}; }
    Vec2 operator*(const double scalar) const { return {x() * scalar, y() * scalar}; }
    Vec2 operator/(const double scalar) const { return {x() / scalar, y() / scalar}; }
    Vec2& operator+=(const Vec2& other) { x() += other.x(); y() += other.y(); return *this; }
    Vec2& operator/=(const double scalar) { x() /= scalar; y() /= scalar; return *this; }
    double dot(const Vec2& other) const { return x() * other.x() + y() * other.y(); }
    double norm() const { return std::sqrt(dot(*this)); }
    Vec2 normalized() const { return *this / norm(); }
    Vec2 cwiseMin(const Vec2& other) const
    { return {std::min(x(), other.x()), std::min(y(), other.y())}; }
    Vec2 cwiseMax(const Vec2& other) const
    { return {std::max(x(), other.x()), std::max(y(), other.y())}; }
    double prod() const { return x() * y(); }
    static Vec2 Zero() { return {}; }
    static Vec2 Constant(const double value) { return {value, value}; }
};

struct Vec3
{
    double values[3]{};
    Vec3() = default;
    Vec3(const double x, const double y, const double z) : values{x, y, z} {}
    double& x() { return values[0]; }
    double& y() { return values[1]; }
    double& z() { return values[2]; }
    double x() const { return values[0]; }
    double y() const { return values[1]; }
    double z() const { return values[2]; }
    double& operator[](const int index) { return values[index]; }
    double operator[](const int index) const { return values[index]; }
    Vec3 operator+(const Vec3& other) const { return {x() + other.x(), y() + other.y(), z() + other.z()}; }
    Vec3 operator-(const Vec3& other) const { return {x() - other.x(), y() - other.y(), z() - other.z()}; }
    Vec3 operator-() const { return {-x(), -y(), -z()}; }
    Vec3 operator*(const double scalar) const { return {x() * scalar, y() * scalar, z() * scalar}; }
    Vec3 operator/(const double scalar) const { return {x() / scalar, y() / scalar, z() / scalar}; }
    Vec3& operator+=(const Vec3& other) { x() += other.x(); y() += other.y(); z() += other.z(); return *this; }
    Vec3& operator*=(const double scalar) { x() *= scalar; y() *= scalar; z() *= scalar; return *this; }
    Vec3& operator/=(const double scalar) { x() /= scalar; y() /= scalar; z() /= scalar; return *this; }
    double dot(const Vec3& other) const { return x() * other.x() + y() * other.y() + z() * other.z(); }
    Vec3 cross(const Vec3& other) const
    {
        return {y() * other.z() - z() * other.y(), z() * other.x() - x() * other.z(),
                x() * other.y() - y() * other.x()};
    }
    double squaredNorm() const { return dot(*this); }
    double norm() const { return std::sqrt(squaredNorm()); }
    Vec3 normalized() const { return *this / norm(); }
    Vec3 cwiseMin(const Vec3& other) const
    { return {std::min(x(), other.x()), std::min(y(), other.y()), std::min(z(), other.z())}; }
    Vec3 cwiseMax(const Vec3& other) const
    { return {std::max(x(), other.x()), std::max(y(), other.y()), std::max(z(), other.z())}; }
    double prod() const { return x() * y() * z(); }
    static Vec3 Zero() { return {}; }
    static Vec3 Constant(const double value) { return {value, value, value}; }
    static Vec3 UnitX() { return {1.0, 0.0, 0.0}; }
    static Vec3 UnitY() { return {0.0, 1.0, 0.0}; }
};

struct Mat3
{
    double values[3][3]{};
    struct Column
    {
        Mat3& matrix;
        int index;
        Column& operator=(const Vec3& vector)
        {
            for (int row = 0; row < 3; ++row) matrix.values[row][index] = vector[row];
            return *this;
        }
        Column& operator*=(const double scalar)
        {
            for (int row = 0; row < 3; ++row) matrix.values[row][index] *= scalar;
            return *this;
        }
        operator Vec3() const
        { return {matrix.values[0][index], matrix.values[1][index], matrix.values[2][index]}; }
        Vec3 operator*(const double scalar) const { return static_cast<Vec3>(*this) * scalar; }
        double dot(const Vec3& vector) const { return static_cast<Vec3>(*this).dot(vector); }
    };
    Column col(const int index) { return {*this, index}; }
    Vec3 col(const int index) const
    { return {values[0][index], values[1][index], values[2][index]}; }
    Vec3 transposeMultiply(const Vec3& vector) const
    { return {col(0).dot(vector), col(1).dot(vector), col(2).dot(vector)}; }
    Vec3 operator*(const Vec3& vector) const
    { return col(0) * vector.x() + col(1) * vector.y() + col(2) * vector.z(); }
    Mat3& operator+=(const Mat3& other)
    {
        for (int row = 0; row < 3; ++row)
            for (int column = 0; column < 3; ++column)
                values[row][column] += other.values[row][column];
        return *this;
    }
    double determinant() const
    {
        return values[0][0] * (values[1][1] * values[2][2] - values[1][2] * values[2][1])
             - values[0][1] * (values[1][0] * values[2][2] - values[1][2] * values[2][0])
             + values[0][2] * (values[1][0] * values[2][1] - values[1][1] * values[2][0]);
    }
    static Mat3 Identity()
    {
        Mat3 result;
        for (int index = 0; index < 3; ++index) result.values[index][index] = 1.0;
        return result;
    }
    static Mat3 Zero() { return {}; }
};

Mat3 outerProduct(const Vec3& vector)
{
    Mat3 result;
    for (int row = 0; row < 3; ++row)
        for (int column = 0; column < 3; ++column)
            result.values[row][column] = vector[row] * vector[column];
    return result;
}

Mat3 symmetricEigenvectors(const Mat3& input)
{
    Mat3 matrix = input;
    Mat3 vectors = Mat3::Identity();
    for (int iteration = 0; iteration < 50; ++iteration)
    {
        int p = 0;
        int q = 1;
        double maximum = std::abs(matrix.values[0][1]);
        for (int row = 0; row < 3; ++row)
            for (int column = row + 1; column < 3; ++column)
                if (std::abs(matrix.values[row][column]) > maximum)
                {
                    p = row;
                    q = column;
                    maximum = std::abs(matrix.values[row][column]);
                }
        if (maximum <= 1.0e-15) break;
        const double angle = 0.5 * std::atan2(2.0 * matrix.values[p][q],
                                             matrix.values[q][q] - matrix.values[p][p]);
        const double cosine = std::cos(angle);
        const double sine = std::sin(angle);
        for (int index = 0; index < 3; ++index)
        {
            const double first = matrix.values[index][p];
            const double second = matrix.values[index][q];
            matrix.values[index][p] = cosine * first - sine * second;
            matrix.values[index][q] = sine * first + cosine * second;
        }
        for (int index = 0; index < 3; ++index)
        {
            const double first = matrix.values[p][index];
            const double second = matrix.values[q][index];
            matrix.values[p][index] = cosine * first - sine * second;
            matrix.values[q][index] = sine * first + cosine * second;
            const double vector_first = vectors.values[index][p];
            const double vector_second = vectors.values[index][q];
            vectors.values[index][p] = cosine * vector_first - sine * vector_second;
            vectors.values[index][q] = sine * vector_first + cosine * vector_second;
        }
    }
    return vectors;
}
using Face = std::array<std::uint32_t, 3>;

struct Mesh
{
    std::vector<Vec3> vertices;
    std::vector<Face> faces;
};

enum class Kind
{
    Rectangle,
    Triangle,
    Polygon,
    Disk,
    Annulus,
    CylindricalBand,
    ConicalBand
};

const char* kindName(Kind kind);

struct Primitive
{
    Kind kind = Kind::Triangle;
    Vec3 center = Vec3::Zero();
    Mat3 axes = Mat3::Identity();
    Vec3 half_size = Vec3::Zero();
    double base_radius = 0.0;
    double top_radius = 0.0;
    double inner_radius = 0.0;
    double height = 0.0;
    std::uint32_t segments = 24;
    std::array<Vec3, 3> triangle{};
    std::vector<Vec3> polygon;
    // A non-empty cyclic range table turns a cylindrical/conical band into a
    // surface clipped in (angle, axial-distance) parameter space.  Entry i is
    // the conservative lower/upper axial boundary at angle 2*pi*i/N.
    std::vector<std::array<double, 2>> band_axial_ranges;
    // Entry i says whether the angular interval [i,i+1] is part of a trimmed
    // band.  Inactive intervals emit no triangles, so a quarter-round or an
    // exposed cylindrical strip stays an open surface instead of becoming a
    // misleading complete cylinder.
    std::vector<std::uint8_t> band_active_segments;
    double volume = 0.0;
};

// A box is only a fitting/certification aid. It is never a primitive: callers
// that accept a box-shaped approximation emit its six rectangle faces directly.
struct BoxFit
{
    Vec3 center = Vec3::Zero();
    Mat3 axes = Mat3::Identity();
    Vec3 half_size = Vec3::Zero();
    double volume = 0.0;
};

struct OutputPrimitive
{
    Primitive primitive;
    std::vector<std::uint32_t> source_faces;
    bool shallow_shell_coalesced = false;
    bool preserves_cavity_opening = false;
    // Nonzero values identify faces emitted together as one complete 3D shell.
    // A planar pass may simplify within one shell, but must not create a partial
    // bridge between different shells or between a shell and an unrelated face.
    std::uint64_t enclosure_group = 0;
};

bool isSurfaceCandidateKind(const Kind kind)
{
    return kind == Kind::Rectangle || kind == Kind::Triangle ||
           kind == Kind::Polygon || kind == Kind::Disk ||
           kind == Kind::Annulus || kind == Kind::CylindricalBand ||
           kind == Kind::ConicalBand;
}

bool isCertifiedRoundSurfaceKind(const Kind kind)
{
    return kind == Kind::Disk || kind == Kind::Annulus ||
           kind == Kind::CylindricalBand || kind == Kind::ConicalBand;
}

void requireSurfaceCandidates(const std::vector<OutputPrimitive>& primitives,
                              const char* stage)
{
    for (const auto& item : primitives)
        if (!isSurfaceCandidateKind(item.primitive.kind))
            throw std::runtime_error(
                std::string(stage) +
                " produced a volumetric primitive; only surface candidates are allowed");
}

struct CertifiedExtrusion
{
    Vec3 origin = Vec3::Zero();
    Mat3 frame = Mat3::Identity();
    std::vector<Vec2> boundary;
    double lower_distance = 0.0;
    double upper_distance = 0.0;
    std::uint64_t enclosure_group = 0;
};

std::vector<CertifiedExtrusion> recognizeCertifiedPrismaticVolumes(
    const std::vector<OutputPrimitive>& primitives,
    double tolerance);

std::vector<CertifiedExtrusion> recognizeEnclosureGroupExtrusions(
    const std::vector<OutputPrimitive>& primitives,
    double tolerance,
    bool include_ungrouped = true);

struct OcclusionClipStats
{
    std::size_t clipped_primitives = 0;
    std::size_t removed_primitives = 0;
    std::size_t input_triangles = 0;
    std::size_t output_triangles = 0;
    double removed_area = 0.0;
};

std::vector<OutputPrimitive> clipPlanarOcclusionByClosedVolumes(
    std::vector<OutputPrimitive> primitives,
    const std::vector<BoxFit>& closed_volumes,
    const std::vector<CertifiedExtrusion>& closed_extrusions,
    double tolerance,
    std::size_t& removed_count,
    std::vector<std::uint32_t>* excluded_faces = nullptr,
    bool preserve_triangle_workload = true,
    OcclusionClipStats* clip_stats = nullptr);

std::vector<OutputPrimitive> clipRegularizedPlanarOcclusion(
    const Mesh& mesh,
    std::vector<OutputPrimitive> primitives,
    double tolerance,
    const std::vector<OutputPrimitive>& coverage_certificates,
    const std::vector<CertifiedExtrusion>* known_owner_extrusions,
    const std::vector<std::uint32_t>* excluded_redundant_faces,
    std::size_t maximum_triangle_workload,
    std::size_t& removed_count,
    std::size_t& recognized_extrusions,
    OcclusionClipStats& clip_stats,
    bool& rolled_back);

void promoteToSemanticPrimitives(std::vector<OutputPrimitive>& primitives);

struct PrimitiveMesh
{
    std::vector<Vec3> vertices;
    std::vector<Face> faces;
};

PrimitiveMesh triangulatePrimitive(const Primitive& primitive);

std::uint64_t edgeKey(std::uint32_t first, std::uint32_t second)
{
    if (first > second) std::swap(first, second);
    return (static_cast<std::uint64_t>(first) << 32U) | second;
}

std::uint32_t parseObjIndex(const std::string& token, const std::size_t vertex_count)
{
    const std::size_t slash = token.find('/');
    const int raw = std::stoi(token.substr(0, slash));
    const long long index = raw > 0 ? raw - 1LL : static_cast<long long>(vertex_count) + raw;
    if (index < 0 || index >= static_cast<long long>(vertex_count))
        throw std::runtime_error("OBJ face index is outside the vertex array");
    return static_cast<std::uint32_t>(index);
}

Mesh readObj(const std::filesystem::path& path)
{
    std::ifstream stream(path);
    if (!stream) throw std::runtime_error("cannot open OBJ: " + path.string());
    Mesh mesh;
    std::string line;
    while (std::getline(stream, line))
    {
        std::istringstream fields(line);
        std::string tag;
        fields >> tag;
        if (tag == "v")
        {
            Vec3 vertex;
            if (!(fields >> vertex.x() >> vertex.y() >> vertex.z()))
                throw std::runtime_error("invalid OBJ vertex in " + path.string());
            mesh.vertices.push_back(vertex);
        }
        else if (tag == "f")
        {
            std::vector<std::uint32_t> polygon;
            std::string token;
            while (fields >> token) polygon.push_back(parseObjIndex(token, mesh.vertices.size()));
            for (std::size_t index = 1; index + 1 < polygon.size(); ++index)
                mesh.faces.push_back({polygon[0], polygon[index], polygon[index + 1]});
        }
    }
    if (mesh.vertices.empty() || mesh.faces.empty())
        throw std::runtime_error("OBJ has no triangle mesh: " + path.string());
    return mesh;
}

std::size_t weldCoincidentVertices(Mesh& mesh)
{
    std::map<std::array<double, 3>, std::uint32_t> canonical;
    std::vector<std::uint32_t> remap(mesh.vertices.size());
    std::size_t welded = 0;
    for (std::uint32_t id = 0; id < mesh.vertices.size(); ++id)
    {
        const Vec3& vertex = mesh.vertices[id];
        const std::array<double, 3> key{vertex.x(), vertex.y(), vertex.z()};
        const auto [iterator, inserted] = canonical.emplace(key, id);
        remap[id] = iterator->second;
        if (!inserted) ++welded;
    }
    for (Face& face : mesh.faces)
        for (auto& vertex : face) vertex = remap[vertex];
    return welded;
}

std::size_t dropDegenerateFaces(Mesh& mesh)
{
    std::vector<Face> kept;
    kept.reserve(mesh.faces.size());
    for (const Face& face : mesh.faces)
    {
        const Vec3 first = mesh.vertices[face[1]] - mesh.vertices[face[0]];
        const Vec3 second = mesh.vertices[face[2]] - mesh.vertices[face[0]];
        const double scale = std::max({first.squaredNorm(), second.squaredNorm(),
                                       (mesh.vertices[face[2]] - mesh.vertices[face[1]]).squaredNorm()});
        // OBJ decimal coordinates frequently encode mathematically collinear
        // triples with a few ulps of residual cross product.  Retaining those
        // numerical slivers creates zero-area semantic polygons later.  The
        // test is scale-relative and still keeps triangles whose altitude is
        // at least 1e-12 of their longest edge.
        if (first.cross(second).norm() > 1.0e-12 * scale)
            kept.push_back(face);
    }
    const std::size_t removed = mesh.faces.size() - kept.size();
    mesh.faces = std::move(kept);
    if (mesh.faces.empty()) throw std::runtime_error("all OBJ triangles are degenerate");
    return removed;
}

std::vector<std::uint32_t> uniqueVertices(const Mesh& mesh,
                                          const std::vector<std::uint32_t>& faces)
{
    std::vector<std::uint32_t> result;
    result.reserve(faces.size() * 3);
    for (const std::uint32_t face_id : faces)
        result.insert(result.end(), mesh.faces[face_id].begin(), mesh.faces[face_id].end());
    std::sort(result.begin(), result.end());
    result.erase(std::unique(result.begin(), result.end()), result.end());
    return result;
}

double faceArea(const Mesh& mesh, const std::uint32_t face_id)
{
    const Face& face = mesh.faces[face_id];
    return 0.5 * (mesh.vertices[face[1]] - mesh.vertices[face[0]])
                     .cross(mesh.vertices[face[2]] - mesh.vertices[face[0]])
                     .norm();
}

double facesArea(const Mesh& mesh, const std::vector<std::uint32_t>& faces)
{
    double area = 0.0;
    for (const auto face_id : faces) area += faceArea(mesh, face_id);
    return area;
}

std::pair<Vec3, Vec3> faceBounds(const Mesh& mesh,
                                 const std::vector<std::uint32_t>& faces)
{
    Vec3 lower = Vec3::Constant(std::numeric_limits<double>::infinity());
    Vec3 upper = Vec3::Constant(-std::numeric_limits<double>::infinity());
    for (const auto vertex : uniqueVertices(mesh, faces))
    {
        lower = lower.cwiseMin(mesh.vertices[vertex]);
        upper = upper.cwiseMax(mesh.vertices[vertex]);
    }
    return {lower, upper};
}

double primitiveSurfaceArea(const Primitive& primitive, const double dimensional_tolerance)
{
    if (primitive.kind == Kind::Polygon)
    {
        const PrimitiveMesh mesh = triangulatePrimitive(primitive);
        double area = 0.0;
        for (const Face& face : mesh.faces)
            area += 0.5 * (mesh.vertices[face[1]] - mesh.vertices[face[0]])
                              .cross(mesh.vertices[face[2]] - mesh.vertices[face[0]]).norm();
        return area;
    }
    if (primitive.kind == Kind::Triangle)
        return 0.5 * (primitive.triangle[1] - primitive.triangle[0])
                         .cross(primitive.triangle[2] - primitive.triangle[0])
                         .norm();
    if (primitive.kind == Kind::Disk)
        return std::numbers::pi * primitive.base_radius * primitive.base_radius;
    if (primitive.kind == Kind::Annulus)
        return std::numbers::pi *
            (primitive.base_radius * primitive.base_radius -
             primitive.inner_radius * primitive.inner_radius);
    if (primitive.kind == Kind::CylindricalBand ||
        primitive.kind == Kind::ConicalBand)
    {
        if (!primitive.band_axial_ranges.empty())
        {
            const PrimitiveMesh mesh = triangulatePrimitive(primitive);
            double area = 0.0;
            for (const Face& face : mesh.faces)
                area += 0.5 *
                    (mesh.vertices[face[1]] - mesh.vertices[face[0]])
                        .cross(mesh.vertices[face[2]] -
                               mesh.vertices[face[0]]).norm();
            return area;
        }
        const double slant = std::hypot(
            primitive.height, primitive.top_radius - primitive.base_radius);
        return std::numbers::pi *
            (primitive.base_radius + primitive.top_radius) * slant;
    }
    std::array<double, 3> extents{
        2.0 * primitive.half_size.x(),
        2.0 * primitive.half_size.y(),
        2.0 * primitive.half_size.z(),
    };
    std::sort(extents.begin(), extents.end());
    if (primitive.kind == Kind::Rectangle || extents[0] <= dimensional_tolerance)
        return extents[1] * extents[2];
    return 2.0 * (extents[0] * extents[1] + extents[1] * extents[2] +
                  extents[2] * extents[0]);
}

double boundsSeparationRatio(const Vec3& first_lower,
                             const Vec3& first_upper,
                             const Vec3& second_lower,
                             const Vec3& second_upper)
{
    Vec3 gap = Vec3::Zero();
    for (int axis = 0; axis < 3; ++axis)
    {
        if (first_upper[axis] < second_lower[axis])
            gap[axis] = second_lower[axis] - first_upper[axis];
        else if (second_upper[axis] < first_lower[axis])
            gap[axis] = first_lower[axis] - second_upper[axis];
    }
    const Vec3 merged_lower = first_lower.cwiseMin(second_lower);
    const Vec3 merged_upper = first_upper.cwiseMax(second_upper);
    return gap.norm() / std::max((merged_upper - merged_lower).norm(), 1.0e-30);
}

Mat3 orthonormalFrame(const Vec3& first_axis)
{
    Vec3 first = first_axis.normalized();
    const Vec3 helper = std::abs(first.x()) < 0.8 ? Vec3::UnitX() : Vec3::UnitY();
    Vec3 second = first.cross(helper).normalized();
    Vec3 third = first.cross(second).normalized();
    Mat3 frame;
    frame.col(0) = first;
    frame.col(1) = second;
    frame.col(2) = third;
    return frame;
}

BoxFit fitBoxOnFrame(const Mesh& mesh,
                     const std::vector<std::uint32_t>& vertex_ids,
                     const Mat3& frame)
{
    Vec3 lower = Vec3::Constant(std::numeric_limits<double>::infinity());
    Vec3 upper = Vec3::Constant(-std::numeric_limits<double>::infinity());
    for (const std::uint32_t id : vertex_ids)
    {
        const Vec3 local = frame.transposeMultiply(mesh.vertices[id]);
        lower = lower.cwiseMin(local);
        upper = upper.cwiseMax(local);
    }
    BoxFit result;
    result.axes = frame;
    result.center = frame * ((lower + upper) * 0.5);
    result.half_size = (upper - lower) * 0.5;
    result.volume = 8.0 * result.half_size.prod();
    return result;
}

std::vector<Mat3> candidateFrames(const Mesh& mesh,
                                  const std::vector<std::uint32_t>& vertex_ids)
{
    Vec3 mean = Vec3::Zero();
    for (const auto id : vertex_ids) mean += mesh.vertices[id];
    mean /= static_cast<double>(vertex_ids.size());
    Mat3 covariance = Mat3::Zero();
    for (const auto id : vertex_ids)
    {
        const Vec3 centered = mesh.vertices[id] - mean;
        covariance += outerProduct(centered);
    }
    Mat3 pca = symmetricEigenvectors(covariance);
    if (pca.determinant() < 0.0) pca.col(0) *= -1.0;
    return {pca, Mat3::Identity()};
}

BoxFit fitBox(const Mesh& mesh, const std::vector<std::uint32_t>& vertex_ids)
{
    BoxFit best;
    best.volume = std::numeric_limits<double>::infinity();
    for (const Mat3& frame : candidateFrames(mesh, vertex_ids))
    {
        BoxFit candidate = fitBoxOnFrame(mesh, vertex_ids, frame);
        const double candidate_area = std::abs(candidate.half_size.x()) +
                                      std::abs(candidate.half_size.y()) +
                                      std::abs(candidate.half_size.z());
        const double best_area = std::abs(best.half_size.x()) +
                                 std::abs(best.half_size.y()) +
                                 std::abs(best.half_size.z());
        if (candidate.volume < best.volume ||
            (candidate.volume == best.volume && candidate_area < best_area))
            best = candidate;
    }
    return best;
}

double boxSurfaceArea(const BoxFit& box)
{
    const Vec3 extent = box.half_size * 2.0;
    return 2.0 * (extent.x() * extent.y() + extent.y() * extent.z() +
                  extent.z() * extent.x());
}

std::optional<Vec2> leastSquaresCircleCenter(
    const std::vector<Vec2>& points)
{
    if (points.size() < 3) return std::nullopt;
    Vec2 mean = Vec2::Zero();
    double mean_squared_radius = 0.0;
    for (const Vec2& point : points)
    {
        mean += point;
        mean_squared_radius += point.dot(point);
    }
    mean /= static_cast<double>(points.size());
    mean_squared_radius /= static_cast<double>(points.size());

    double xx = 0.0;
    double xy = 0.0;
    double yy = 0.0;
    double bx = 0.0;
    double by = 0.0;
    for (const Vec2& point : points)
    {
        const Vec2 centered = point - mean;
        const double radial_equation =
            point.dot(point) - mean_squared_radius;
        xx += 2.0 * centered.x() * centered.x();
        xy += 2.0 * centered.x() * centered.y();
        yy += 2.0 * centered.y() * centered.y();
        bx += centered.x() * radial_equation;
        by += centered.y() * radial_equation;
    }
    const double determinant = xx * yy - xy * xy;
    const double scale = std::max(xx + yy, 1.0e-30);
    if (determinant <= scale * scale * 1.0e-10)
        return std::nullopt;
    return Vec2{(bx * yy - by * xy) / determinant,
                (by * xx - bx * xy) / determinant};
}

Primitive fitConeOnAxis(const Mesh& mesh,
                        const std::vector<std::uint32_t>& vertex_ids,
                        const Vec3& axis,
                        const std::uint32_t segments,
                        const double analytic_tolerance)
{
    const Mat3 frame = orthonormalFrame(axis);
    Vec3 mean = Vec3::Zero();
    for (const auto id : vertex_ids) mean += mesh.vertices[id];
    mean /= static_cast<double>(vertex_ids.size());
    double lower = std::numeric_limits<double>::infinity();
    double upper = -std::numeric_limits<double>::infinity();
    std::vector<Vec3> local_vertices;
    local_vertices.reserve(vertex_ids.size());
    for (const auto id : vertex_ids)
    {
        const Vec3 local = frame.transposeMultiply(mesh.vertices[id] - mean);
        local_vertices.push_back(local);
        lower = std::min(lower, local.x());
        upper = std::max(upper, local.x());
    }
    const double height = upper - lower;
    const double end_tolerance = std::max({height * 1.0e-6,
                                           analytic_tolerance, 1.0e-12});
    std::vector<Vec2> base_points;
    std::vector<Vec2> top_points;
    for (const Vec3& local : local_vertices)
    {
        if (local.x() <= lower + end_tolerance)
            base_points.emplace_back(local.y(), local.z());
        if (local.x() >= upper - end_tolerance)
            top_points.emplace_back(local.y(), local.z());
    }
    const auto base_center = leastSquaresCircleCenter(base_points);
    const auto top_center = leastSquaresCircleCenter(top_points);
    Vec2 radial_center = Vec2::Zero();
    if (base_center && top_center)
        radial_center = (*base_center + *top_center) * 0.5;
    else if (base_center)
        radial_center = *base_center;
    else if (top_center)
        radial_center = *top_center;

    double base_radius = 0.0;
    double top_radius = 0.0;
    for (const Vec3& local : local_vertices)
    {
        const double radius = std::hypot(
            local.y() - radial_center.x(),
            local.z() - radial_center.y());
        if (local.x() <= lower + end_tolerance) base_radius = std::max(base_radius, radius);
        if (local.x() >= upper - end_tolerance) top_radius = std::max(top_radius, radius);
    }
    // Expand both end radii by the largest deficit so the fitted side surface
    // still covers every input vertex, including intermediate axial rings.
    double deficit = 0.0;
    if (height > 0.0)
    {
        for (const Vec3& local : local_vertices)
        {
            const double fraction = std::clamp((local.x() - lower) / height, 0.0, 1.0);
            const double expected = base_radius * (1.0 - fraction) + top_radius * fraction;
            deficit = std::max(deficit, std::hypot(
                local.y() - radial_center.x(),
                local.z() - radial_center.y()) - expected);
        }
    }
    base_radius += std::max(deficit, 0.0);
    top_radius += std::max(deficit, 0.0);
    // The emitted regular polygon is inscribed in its circumcircle. Expand to
    // a circumscribed polygon so every point admitted by the continuous radial
    // fit is still inside the triangulated proxy.
    const double polygon_expansion = segments >= 3
        ? 1.0 / std::cos(std::numbers::pi / static_cast<double>(segments))
        : std::numeric_limits<double>::infinity();
    base_radius *= polygon_expansion * (1.0 + 1.0e-12);
    top_radius *= polygon_expansion * (1.0 + 1.0e-12);
    Primitive result;
    result.kind = std::abs(base_radius - top_radius) <=
            1.0e-8 * std::max({base_radius, top_radius, height})
        ? Kind::CylindricalBand : Kind::ConicalBand;
    result.axes = frame;
    result.height = height;
    result.center = mean + frame.col(1) * radial_center.x() +
        frame.col(2) * radial_center.y() + axis.normalized() * lower;
    result.base_radius = base_radius;
    result.top_radius = top_radius;
    result.segments = segments;
    return result;
}

double conservativeConeProfileSlack(
    const Mesh& mesh, const std::vector<std::uint32_t>& vertex_ids,
    const Primitive& cone)
{
    constexpr std::size_t bin_count = 128;
    std::array<double, bin_count> outer_radii;
    std::array<double, bin_count> outer_fractions{};
    outer_radii.fill(-std::numeric_limits<double>::infinity());
    for (const auto vertex : vertex_ids)
    {
        const Vec3 local = cone.axes.transposeMultiply(
            mesh.vertices[vertex] - cone.center);
        const double fraction = std::clamp(
            local.x() / std::max(cone.height, 1.0e-30), 0.0, 1.0);
        const double radius = std::hypot(local.y(), local.z());
        const std::size_t bin = std::min(
            bin_count - 1,
            static_cast<std::size_t>(fraction * bin_count));
        if (radius > outer_radii[bin])
        {
            outer_radii[bin] = radius;
            outer_fractions[bin] = fraction;
        }
    }
    double maximum_slack = 0.0;
    for (std::size_t bin = 0; bin < bin_count; ++bin)
    {
        if (!std::isfinite(outer_radii[bin])) continue;
        const double expected =
            cone.base_radius * (1.0 - outer_fractions[bin]) +
            cone.top_radius * outer_fractions[bin];
        maximum_slack = std::max(
            maximum_slack, expected - outer_radii[bin]);
    }
    return maximum_slack;
}

double refitConservativeConeProfile(
    const Mesh& mesh, const std::vector<std::uint32_t>& vertex_ids,
    const Vec2& projected_center_offset, Primitive& cone)
{
    cone.center += cone.axes.col(1) * projected_center_offset.x() +
                   cone.axes.col(2) * projected_center_offset.y();
    struct Sample
    {
        double fraction = 0.0;
        double radius = 0.0;
    };
    std::vector<Sample> samples;
    samples.reserve(vertex_ids.size());
    double maximum_radius = 0.0;
    for (const auto vertex : vertex_ids)
    {
        const Vec3 local = cone.axes.transposeMultiply(
            mesh.vertices[vertex] - cone.center);
        const double radius = std::hypot(local.y(), local.z());
        maximum_radius = std::max(maximum_radius, radius);
        samples.push_back({
            std::clamp(local.x() / std::max(cone.height, 1.0e-30),
                       0.0, 1.0),
            radius});
    }
    constexpr std::size_t bin_count = 128;
    std::array<double, bin_count> outer_radii;
    std::array<double, bin_count> outer_fractions{};
    outer_radii.fill(-std::numeric_limits<double>::infinity());
    for (const auto& sample : samples)
    {
        const std::size_t bin = std::min(
            bin_count - 1,
            static_cast<std::size_t>(sample.fraction * bin_count));
        if (sample.radius > outer_radii[bin])
        {
            outer_radii[bin] = sample.radius;
            outer_fractions[bin] = sample.fraction;
        }
    }
    const auto radiiForSlope = [&](const double radius_difference)
    {
        double base_radius = std::max(0.0, -radius_difference);
        for (const auto& sample : samples)
            base_radius = std::max(
                base_radius,
                sample.radius - sample.fraction * radius_difference);
        return std::array<double, 2>{
            base_radius, base_radius + radius_difference};
    };
    const auto objective = [&](const double radius_difference)
    {
        const auto radii = radiiForSlope(radius_difference);
        double maximum_slack = 0.0;
        for (std::size_t bin = 0; bin < bin_count; ++bin)
        {
            if (!std::isfinite(outer_radii[bin])) continue;
            const double expected =
                radii[0] + outer_fractions[bin] * radius_difference;
            maximum_slack = std::max(
                maximum_slack, expected - outer_radii[bin]);
        }
        return maximum_slack +
            1.0e-12 * (radii[0] + radii[1]);
    };
    double lower = -2.0 * maximum_radius;
    double upper = 2.0 * maximum_radius;
    for (int iteration = 0; iteration < 96; ++iteration)
    {
        const double first = lower + (upper - lower) / 3.0;
        const double second = upper - (upper - lower) / 3.0;
        if (objective(first) <= objective(second)) upper = second;
        else lower = first;
    }
    const auto radii = radiiForSlope(0.5 * (lower + upper));
    const double polygon_expansion = cone.segments >= 3
        ? 1.0 / std::cos(
            std::numbers::pi / static_cast<double>(cone.segments))
        : std::numeric_limits<double>::infinity();
    cone.base_radius = radii[0] * polygon_expansion * (1.0 + 1.0e-12);
    cone.top_radius = radii[1] * polygon_expansion * (1.0 + 1.0e-12);
    cone.kind = std::abs(cone.base_radius - cone.top_radius) <=
            1.0e-8 * std::max(
                {cone.base_radius, cone.top_radius, cone.height})
        ? Kind::CylindricalBand : Kind::ConicalBand;
    return conservativeConeProfileSlack(mesh, vertex_ids, cone);
}

std::size_t faceComponentCount(const Mesh& mesh,
                               const std::vector<std::uint32_t>& faces);

std::size_t triangulatedFaceCount(
    const std::vector<OutputPrimitive>& primitives);

void appendBoxRectangles(
    std::vector<OutputPrimitive>& output,
    const BoxFit& box,
    const std::vector<std::uint32_t>& source_faces,
    int covered_face_axis,
    double covered_face_sign,
    std::uint64_t enclosure_group);

std::vector<Vec2> convexHull(std::vector<Vec2> points, double tolerance);

struct RingCircularity
{
    std::size_t point_count = 0;
    double radial_spread_ratio = std::numeric_limits<double>::infinity();
    double angular_gap_ratio = std::numeric_limits<double>::infinity();
};

bool ringLooksCircular(const std::vector<Vec2>& points,
                       const double tolerance,
                       RingCircularity* metrics = nullptr)
{
    if (metrics) metrics->point_count = points.size();
    if (points.empty()) return false;
    std::vector<double> radii;
    std::vector<double> angles;
    radii.reserve(points.size());
    angles.reserve(points.size());
    for (const Vec2& point : points)
    {
        radii.push_back(point.norm());
        angles.push_back(std::atan2(point.y(), point.x()));
    }
    const double mean = std::accumulate(radii.begin(), radii.end(), 0.0) / radii.size();
    if (mean <= 0.0) return false;
    const auto [minimum_radius, maximum_radius] = std::minmax_element(
        radii.begin(), radii.end());
    const double radial_spread_ratio =
        (*maximum_radius - *minimum_radius) / mean;
    std::sort(angles.begin(), angles.end());
    double maximum_gap = angles.front() + 2.0 * std::numbers::pi - angles.back();
    for (std::size_t index = 1; index < angles.size(); ++index)
        maximum_gap = std::max(maximum_gap, angles[index] - angles[index - 1]);
    const double allowed_gap = 2.5 *
        (2.0 * std::numbers::pi / angles.size());
    if (metrics)
    {
        metrics->radial_spread_ratio = radial_spread_ratio;
        metrics->angular_gap_ratio = maximum_gap /
            std::max(allowed_gap, 1.0e-30);
    }
    return points.size() >= 6 && radial_spread_ratio <= tolerance &&
        maximum_gap <= allowed_gap;
}

enum class ConeFitFailure
{
    None,
    Geometry,
    EndRing,
    Normal,
    Radial,
};

struct ConeFitDiagnostics
{
    RingCircularity base_ring;
    RingCircularity top_ring;
    bool base_is_apex = false;
    bool top_is_apex = false;
};

bool coneFitIsCertified(const Mesh& mesh,
                        const std::vector<std::uint32_t>& faces,
                        const Primitive& cone,
                        const double radial_tolerance,
                        const double analytic_tolerance,
                        ConeFitFailure* failure = nullptr,
                        ConeFitDiagnostics* diagnostics = nullptr)
{
    if (failure) *failure = ConeFitFailure::None;
    // Decimal CAD exports often duplicate seam vertices instead of sharing OBJ
    // indices. Analytic certification below is geometric and is stronger than
    // exact index connectivity, so only require exact connectivity when no
    // geometric tolerance was supplied.
    if ((analytic_tolerance <= 0.0 && faceComponentCount(mesh, faces) != 1) ||
        cone.height <= 0.0)
    {
        if (failure) *failure = ConeFitFailure::Geometry;
        return false;
    }
    const auto vertices = uniqueVertices(mesh, faces);
    const double end_tolerance = std::max({cone.height * 1.0e-6,
                                           analytic_tolerance, 1.0e-12});
    std::vector<Vec2> base_ring;
    std::vector<Vec2> top_ring;
    std::vector<std::pair<double, double>> samples;
    for (const auto id : vertices)
    {
        const Vec3 local = cone.axes.transposeMultiply(mesh.vertices[id] - cone.center);
        const double radius = std::hypot(local.y(), local.z());
        samples.emplace_back(local.x(), radius);
        if (local.x() <= end_tolerance &&
            radius >= cone.base_radius * (1.0 - 2.0 * radial_tolerance))
            base_ring.emplace_back(local.y(), local.z());
        if (local.x() >= cone.height - end_tolerance &&
            radius >= cone.top_radius * (1.0 - 2.0 * radial_tolerance))
            top_ring.emplace_back(local.y(), local.z());
    }
    const bool base_is_apex = cone.base_radius <= 1.0e-8 *
        std::max(cone.top_radius, cone.height);
    const bool top_is_apex = cone.top_radius <= 1.0e-8 *
        std::max(cone.base_radius, cone.height);
    RingCircularity base_metrics;
    RingCircularity top_metrics;
    const bool base_is_circular = base_is_apex || ringLooksCircular(
        base_ring, 2.0 * radial_tolerance, &base_metrics);
    const bool top_is_circular = top_is_apex || ringLooksCircular(
        top_ring, 2.0 * radial_tolerance, &top_metrics);
    if (diagnostics)
    {
        diagnostics->base_ring = base_metrics;
        diagnostics->top_ring = top_metrics;
        diagnostics->base_is_apex = base_is_apex;
        diagnostics->top_is_apex = top_is_apex;
    }
    if ((!base_is_apex && !base_is_circular) ||
        (!top_is_apex && !top_is_circular) ||
        (base_is_apex && top_is_apex))
    {
        if (failure) *failure = ConeFitFailure::EndRing;
        return false;
    }

    std::size_t lateral_faces = 0;
    const double expected_axial_normal = std::abs(cone.top_radius - cone.base_radius) /
        std::max(std::hypot(cone.height, cone.top_radius - cone.base_radius), 1.0e-30);
    for (const auto face_id : faces)
    {
        const Face& face = mesh.faces[face_id];
        bool at_base = true;
        bool at_top = true;
        for (const auto vertex : face)
        {
            const double axial = cone.axes.transposeMultiply(
                mesh.vertices[vertex] - cone.center).x();
            at_base &= std::abs(axial) <= end_tolerance;
            at_top &= std::abs(axial - cone.height) <= end_tolerance;
        }
        if (at_base || at_top) continue;
        Vec3 face_normal = (mesh.vertices[face[1]] - mesh.vertices[face[0]])
            .cross(mesh.vertices[face[2]] - mesh.vertices[face[0]]);
        if (face_normal.norm() <= 1.0e-30) continue;
        face_normal = face_normal.normalized();
        const double axial_normal = std::abs(face_normal.dot(cone.axes.col(0)));
        if (std::abs(axial_normal - expected_axial_normal) > 0.15)
        {
            if (failure) *failure = ConeFitFailure::Normal;
            return false;
        }
        ++lateral_faces;
    }
    if (lateral_faces == 0)
    {
        if (failure) *failure = ConeFitFailure::Normal;
        return false;
    }

    // Vertices away from cap planes must lie on the conical side, not merely
    // somewhere inside an enclosing cone.
    for (const auto [axial, radius] : samples)
    {
        if (axial <= end_tolerance || axial >= cone.height - end_tolerance) continue;
        const double fraction = std::clamp(axial / cone.height, 0.0, 1.0);
        const double expected = cone.base_radius * (1.0 - fraction) +
                                cone.top_radius * fraction;
        if (expected > 0.0 && radius < expected * (1.0 - 2.0 * radial_tolerance))
        {
            if (failure) *failure = ConeFitFailure::Radial;
            return false;
        }
    }
    return true;
}

std::optional<Primitive> fitCertifiedRevolvedSurface(
    const Mesh& mesh,
    const std::vector<std::uint32_t>& faces,
    const PrimitiveMeshAnalysisOptions& options,
    const double analytic_tolerance = 0.0)
{
    const auto vertices = uniqueVertices(mesh, faces);
    std::optional<Primitive> best;
    double best_surface_area = std::numeric_limits<double>::infinity();
    if (!options.allow_round_surfaces) return std::nullopt;
    for (const Mat3& frame : candidateFrames(mesh, vertices))
    {
        for (int axis = 0; axis < 3; ++axis)
        {
            Primitive candidate = fitConeOnAxis(
                mesh, vertices, static_cast<Vec3>(frame.col(axis)),
                options.round_surface_segments, analytic_tolerance);
            if (!coneFitIsCertified(
                    mesh, faces, candidate,
                    options.circle_radial_tolerance, analytic_tolerance))
                continue;
            const double surface_area = primitiveSurfaceArea(candidate, 1.0e-12);
            if (!best || surface_area < best_surface_area)
            {
                best = std::move(candidate);
                best_surface_area = surface_area;
            }
        }
    }
    return best;
}

struct ConservativeRevolvedDiagnostics
{
    std::size_t axes_tested = 0;
    std::array<std::size_t, 5> failures{};
    std::size_t relaxed_candidates = 0;
    std::size_t circular_silhouette_candidates = 0;
};

std::optional<Primitive> fitConservativeRevolvedEnvelope(
    const Mesh& mesh,
    const std::vector<std::uint32_t>& faces,
    const PrimitiveMeshAnalysisOptions& options,
    const double analytic_tolerance,
    ConservativeRevolvedDiagnostics* diagnostics = nullptr)
{
    if (!options.allow_round_surfaces || faces.empty()) return std::nullopt;
    const auto vertices = uniqueVertices(mesh, faces);
    std::optional<Primitive> best;
    double best_profile_slack = std::numeric_limits<double>::infinity();
    double best_surface_area = std::numeric_limits<double>::infinity();
    for (const Mat3& frame : candidateFrames(mesh, vertices))
        for (int axis = 0; axis < 3; ++axis)
        {
            if (diagnostics) ++diagnostics->axes_tested;
            Primitive candidate = fitConeOnAxis(
                mesh, vertices, static_cast<Vec3>(frame.col(axis)),
                options.round_surface_segments, analytic_tolerance);
            ConeFitFailure failure = ConeFitFailure::None;
            ConeFitDiagnostics ring_diagnostics;
            const bool exact = coneFitIsCertified(
                mesh, faces, candidate, options.circle_radial_tolerance,
                analytic_tolerance, &failure, &ring_diagnostics);
            if (diagnostics)
                ++diagnostics->failures[static_cast<std::size_t>(failure)];
            std::vector<Vec2> projected;
            projected.reserve(vertices.size());
            for (const auto vertex : vertices)
            {
                const Vec3 local = candidate.axes.transposeMultiply(
                    mesh.vertices[vertex] - candidate.center);
                projected.emplace_back(local.y(), local.z());
            }
            projected = convexHull(std::move(projected), analytic_tolerance);
            const auto silhouette_center =
                leastSquaresCircleCenter(projected);
            std::vector<Vec2> centered_projection = projected;
            if (silhouette_center)
                for (Vec2& point : centered_projection)
                    point = point - *silhouette_center;
            const bool circular_silhouette = silhouette_center &&
                centered_projection.size() >= 6 &&
                ringLooksCircular(centered_projection,
                                  2.0 * options.circle_radial_tolerance);
            // Normal and radial failures occur after both end silhouettes have
            // passed their circular-ring certificates. They indicate interior
            // details or a recessed side, both of which an outer collision
            // proxy may conservatively fill. Geometry and end-ring failures do
            // not establish a circular silhouette and remain hard rejections.
            if (!exact && failure != ConeFitFailure::Normal &&
                failure != ConeFitFailure::Radial &&
                !(failure == ConeFitFailure::EndRing && circular_silhouette))
                continue;
            if (diagnostics && !exact) ++diagnostics->relaxed_candidates;
            if (diagnostics && circular_silhouette)
                ++diagnostics->circular_silhouette_candidates;
            const double profile_slack =
                !exact && circular_silhouette
                ? refitConservativeConeProfile(
                    mesh, vertices, *silhouette_center, candidate)
                : conservativeConeProfileSlack(mesh, vertices, candidate);
            const double surface_area =
                primitiveSurfaceArea(candidate, 1.0e-12);
            if (!best || profile_slack < best_profile_slack - analytic_tolerance ||
                (std::abs(profile_slack - best_profile_slack) <=
                     analytic_tolerance &&
                 surface_area < best_surface_area))
            {
                best = std::move(candidate);
                best_profile_slack = profile_slack;
                best_surface_area = surface_area;
            }
        }
    return best;
}

std::vector<std::vector<std::uint32_t>> coplanarClusters(
    const Mesh& mesh,
    const double tolerance,
    std::vector<std::unordered_set<std::uint32_t>>& cluster_neighbors,
    const std::vector<bool>* included_faces = nullptr)
{
    // CAD OBJ exporters frequently duplicate a vertex for every face even when
    // the coordinates are identical.  Building topology from raw OBJ indices
    // then turns one geometric plate or box shell into thousands of disconnected
    // islands and prevents both coplanar union and protrusion recognition.  Weld
    // coordinates only for topology construction; the source mesh itself stays
    // unchanged, so responsibility and exported geometry retain their original
    // indices and precision.
    struct WeldKey
    {
        std::int64_t x = 0;
        std::int64_t y = 0;
        std::int64_t z = 0;
        bool operator==(const WeldKey&) const = default;
    };
    struct WeldKeyHash
    {
        std::size_t operator()(const WeldKey& key) const noexcept
        {
            std::size_t seed = std::hash<std::int64_t>{}(key.x);
            seed ^= std::hash<std::int64_t>{}(key.y) +
                    0x9e3779b97f4a7c15ULL + (seed << 6U) + (seed >> 2U);
            seed ^= std::hash<std::int64_t>{}(key.z) +
                    0x9e3779b97f4a7c15ULL + (seed << 6U) + (seed >> 2U);
            return seed;
        }
    };
    const double weld_tolerance = std::max(tolerance, 1.0e-12);
    std::unordered_map<WeldKey, std::uint32_t, WeldKeyHash> welded_ids;
    welded_ids.reserve(mesh.vertices.size());
    std::vector<std::uint32_t> welded_vertex(mesh.vertices.size());
    std::uint32_t next_welded_id = 0;
    for (std::size_t vertex = 0; vertex < mesh.vertices.size(); ++vertex)
    {
        const Vec3& point = mesh.vertices[vertex];
        const WeldKey key{
            static_cast<std::int64_t>(std::llround(point.x() / weld_tolerance)),
            static_cast<std::int64_t>(std::llround(point.y() / weld_tolerance)),
            static_cast<std::int64_t>(std::llround(point.z() / weld_tolerance))};
        const auto [found, inserted] = welded_ids.emplace(key, next_welded_id);
        welded_vertex[vertex] = found->second;
        if (inserted) ++next_welded_id;
    }
    std::unordered_map<std::uint64_t, std::vector<std::uint32_t>> edge_faces;
    edge_faces.reserve(mesh.faces.size() * 3);
    std::vector<Vec3> normals(mesh.faces.size());
    for (std::uint32_t id = 0; id < mesh.faces.size(); ++id)
    {
        if (included_faces && !(*included_faces)[id]) continue;
        const Face& face = mesh.faces[id];
        for (int edge = 0; edge < 3; ++edge)
            edge_faces[edgeKey(
                welded_vertex[face[edge]],
                welded_vertex[face[(edge + 1) % 3]])].push_back(id);
        normals[id] = (mesh.vertices[face[1]] - mesh.vertices[face[0]])
                          .cross(mesh.vertices[face[2]] - mesh.vertices[face[0]])
                          .normalized();
    }
    std::vector<std::vector<std::uint32_t>> face_neighbors(mesh.faces.size());
    for (const auto& [key, incident] : edge_faces)
    {
        (void)key;
        for (const auto first : incident)
            for (const auto second : incident)
                if (first != second) face_neighbors[first].push_back(second);
    }

    std::vector<int> owner(mesh.faces.size(), -1);
    if (included_faces)
        for (std::size_t id = 0; id < included_faces->size(); ++id)
            if (!(*included_faces)[id]) owner[id] = -2;
    std::vector<std::vector<std::uint32_t>> clusters;
    for (std::uint32_t seed = 0; seed < mesh.faces.size(); ++seed)
    {
        if (owner[seed] != -1) continue;
        const int cluster_id = static_cast<int>(clusters.size());
        clusters.emplace_back();
        std::queue<std::uint32_t> queue;
        queue.push(seed);
        owner[seed] = cluster_id;
        const Vec3 plane_normal = normals[seed];
        const Vec3 plane_point = mesh.vertices[mesh.faces[seed][0]];
        while (!queue.empty())
        {
            const auto face_id = queue.front();
            queue.pop();
            clusters.back().push_back(face_id);
            for (const auto neighbor : face_neighbors[face_id])
            {
                if (owner[neighbor] >= 0) continue;
                if (std::abs(normals[neighbor].dot(plane_normal)) < 1.0 - 1.0e-8) continue;
                bool in_plane = true;
                for (const auto vertex : mesh.faces[neighbor])
                    in_plane &= std::abs((mesh.vertices[vertex] - plane_point).dot(plane_normal))
                                <= tolerance;
                if (!in_plane) continue;
                owner[neighbor] = cluster_id;
                queue.push(neighbor);
            }
        }
    }

    cluster_neighbors.assign(clusters.size(), {});
    for (const auto& [key, incident] : edge_faces)
    {
        (void)key;
        for (const auto first : incident)
            for (const auto second : incident)
            {
                if (owner[first] < 0 || owner[second] < 0) continue;
                const auto first_cluster = static_cast<std::uint32_t>(owner[first]);
                const auto second_cluster = static_cast<std::uint32_t>(owner[second]);
                if (first_cluster != second_cluster)
                    cluster_neighbors[first_cluster].insert(second_cluster);
            }
    }
    return clusters;
}

struct Bounds2
{
    Vec2 lower = Vec2::Constant(std::numeric_limits<double>::infinity());
    Vec2 upper = Vec2::Constant(-std::numeric_limits<double>::infinity());
};

Bounds2 projectedBounds(const Mesh& mesh,
                        const std::vector<std::uint32_t>& faces,
                        const Vec3& origin,
                        const Mat3& frame)
{
    Bounds2 bounds;
    for (const auto vertex : uniqueVertices(mesh, faces))
    {
        const Vec3 local = frame.transposeMultiply(mesh.vertices[vertex] - origin);
        const Vec2 point(local.y(), local.z());
        bounds.lower = bounds.lower.cwiseMin(point);
        bounds.upper = bounds.upper.cwiseMax(point);
    }
    return bounds;
}

double projectedOverlapRatio(const Bounds2& first, const Bounds2& second)
{
    const Vec2 overlap = first.upper.cwiseMin(second.upper) -
                         first.lower.cwiseMax(second.lower);
    if (overlap.x() <= 0.0 || overlap.y() <= 0.0) return 0.0;
    const double first_area = (first.upper - first.lower).prod();
    const double second_area = (second.upper - second.lower).prod();
    return overlap.prod() / std::max(std::min(first_area, second_area), 1.0e-30);
}

double projectedSeparation(const Bounds2& first, const Bounds2& second)
{
    Vec2 gap = Vec2::Zero();
    for (int axis = 0; axis < 2; ++axis)
    {
        if (first.upper[axis] < second.lower[axis])
            gap[axis] = second.lower[axis] - first.upper[axis];
        else if (second.upper[axis] < first.lower[axis])
            gap[axis] = first.lower[axis] - second.upper[axis];
    }
    return gap.norm();
}

struct RecognizedProtrusion
{
    BoxFit box;
    std::vector<std::uint32_t> faces;
    std::vector<std::uint32_t> clusters;
    int covered_face_axis = -1;
    double covered_face_sign = 0.0;
    // Support protrusions are open five-face shells.  The fitted box can extend
    // behind the finite support surface even when its covered face is omitted,
    // so the four adjacent faces must be clipped to this exact plane.
    bool has_support_plane = false;
    Vec3 support_plane_point = Vec3::Zero();
    Vec3 support_outward = Vec3::Zero();
    // Diagnostics used to validate the scale gate for this generic
    // recognizer; they are geometry-derived and do not affect semantics.
    double footprint_ratio = 0.0;
    double group_depth = 0.0;
};

double boxShellAddedArea(const Mesh& mesh,
                         const std::vector<std::uint32_t>& faces,
                         const BoxFit& box,
                         const int covered_face_axis,
                         const double covered_face_sign,
                         const double inward_band)
{
    double added_area = 0.0;
    constexpr int clipper_precision = 8;
    for (int normal_axis = 0; normal_axis < 3; ++normal_axis)
    {
        const int first_axis = (normal_axis + 1) % 3;
        const int second_axis = (normal_axis + 2) % 3;
        const double face_area = 4.0 * box.half_size[first_axis] *
                                 box.half_size[second_axis];
        for (const double sign : {-1.0, 1.0})
        {
            if (normal_axis == covered_face_axis && sign == covered_face_sign) continue;
            const Vec3 normal = box.axes.col(normal_axis) * sign;
            const Vec3 face_center = box.center + box.axes.col(normal_axis) *
                                                  (sign * box.half_size[normal_axis]);
            Clipper2Lib::PathsD projected;
            for (const auto face_id : faces)
            {
                Clipper2Lib::PathD triangle;
                bool near_shell = true;
                for (const auto vertex_id : mesh.faces[face_id])
                {
                    const Vec3 relative = mesh.vertices[vertex_id] - face_center;
                    const double inward_distance = -relative.dot(normal);
                    if (inward_distance < -inward_band || inward_distance > inward_band)
                    {
                        near_shell = false;
                        break;
                    }
                    triangle.emplace_back(relative.dot(box.axes.col(first_axis)),
                                          relative.dot(box.axes.col(second_axis)));
                }
                if (!near_shell || std::abs(Clipper2Lib::Area(triangle)) <= 1.0e-24)
                    continue;
                if (Clipper2Lib::Area(triangle) < 0.0)
                    std::reverse(triangle.begin(), triangle.end());
                projected.push_back(std::move(triangle));
            }

            double covered_area = 0.0;
            if (!projected.empty())
            {
                const auto united = Clipper2Lib::Union(
                    projected, Clipper2Lib::FillRule::NonZero, clipper_precision);
                const Clipper2Lib::RectD rectangle(
                    -box.half_size[first_axis], -box.half_size[second_axis],
                     box.half_size[first_axis],  box.half_size[second_axis]);
                covered_area = std::abs(Clipper2Lib::Area(
                    Clipper2Lib::RectClip(rectangle, united, clipper_precision)));
            }
            added_area += std::max(face_area - covered_area, 0.0);
        }
    }
    return added_area;
}

std::vector<RecognizedProtrusion> recognizeSupportProtrusions(
    const Mesh& mesh,
    const std::vector<std::vector<std::uint32_t>>& clusters,
    const std::vector<std::unordered_set<std::uint32_t>>& topology_neighbors,
    const PrimitiveMeshAnalysisOptions& options,
    const double model_diagonal,
    const double model_surface_area,
    std::vector<bool>& claimed_clusters,
    double& minimum_candidate_area_excess_ratio)
{
    std::vector<RecognizedProtrusion> result;
    if (!options.recognize_support_protrusions) return result;

    const double dimensional_tolerance = std::max(model_diagonal * 1.0e-9, 1.0e-10);
    std::vector<Vec3> normals(clusters.size());
    std::vector<Vec3> centers(clusters.size());
    std::vector<double> areas(clusters.size());
    std::vector<std::uint32_t> dominant;
    for (std::uint32_t id = 0; id < clusters.size(); ++id)
    {
        const Face& face = mesh.faces[clusters[id].front()];
        normals[id] = (mesh.vertices[face[1]] - mesh.vertices[face[0]])
                          .cross(mesh.vertices[face[2]] - mesh.vertices[face[0]])
                          .normalized();
        const auto [lower, upper] = faceBounds(mesh, clusters[id]);
        centers[id] = (lower + upper) * 0.5;
        areas[id] = facesArea(mesh, clusters[id]);
        if (areas[id] >= 0.03 * model_surface_area) dominant.push_back(id);
    }

    struct SupportCandidate
    {
        std::uint32_t support = 0;
        std::uint32_t partner = 0;
        Vec3 origin = Vec3::Zero();
        Vec3 outward = Vec3::UnitX();
        double outward_depth = 0.0;
    };
    std::vector<SupportCandidate> supports;
    std::unordered_set<std::uint64_t> paired;
    for (const auto first : dominant)
    {
        std::uint32_t best = std::numeric_limits<std::uint32_t>::max();
        double best_separation = std::numeric_limits<double>::infinity();
        const Mat3 frame = orthonormalFrame(normals[first]);
        const Bounds2 first_projection = projectedBounds(
            mesh, clusters[first], centers[first], frame);
        for (const auto second : dominant)
        {
            if (first == second || std::abs(normals[first].dot(normals[second])) < 1.0 - 1.0e-8)
                continue;
            const double separation = std::abs(
                (centers[second] - centers[first]).dot(normals[first]));
            if (separation <= dimensional_tolerance || separation >= best_separation) continue;
            const Bounds2 second_projection = projectedBounds(
                mesh, clusters[second], centers[first], frame);
            if (projectedOverlapRatio(first_projection, second_projection) < 0.5) continue;
            best = second;
            best_separation = separation;
        }
        if (best == std::numeric_limits<std::uint32_t>::max() ||
            !paired.insert(edgeKey(first, best)).second) continue;

        Vec3 first_outward = normals[first];
        if ((centers[best] - centers[first]).dot(first_outward) > 0.0)
            first_outward *= -1.0;
        const Vec3 second_outward = -first_outward;
        const auto outwardDepth = [&](const Vec3& origin, const Vec3& outward)
        {
            double depth = 0.0;
            for (const Vec3& vertex : mesh.vertices)
                depth = std::max(depth, (vertex - origin).dot(outward));
            return depth;
        };
        const double first_depth = outwardDepth(centers[first], first_outward);
        const double second_depth = outwardDepth(centers[best], second_outward);
        if (first_depth <= second_depth)
            supports.push_back({first, best, centers[first], first_outward, first_depth});
        else
            supports.push_back({best, first, centers[best], second_outward, second_depth});
    }
    std::sort(supports.begin(), supports.end(), [](const auto& first, const auto& second)
    { return first.outward_depth < second.outward_depth; });

    const double maximum_gap = options.protrusion_cluster_gap_relative * model_diagonal;
    for (const SupportCandidate& support : supports)
    {
        if (support.outward_depth <= dimensional_tolerance) continue;
        const Mat3 frame = orthonormalFrame(support.outward);
        const Bounds2 support_projection = projectedBounds(
            mesh, clusters[support.support], support.origin, frame);
        // Record how much of the whole-model silhouette the support occupies.
        // Measuring this in the support plane makes the diagnostic invariant
        // to translation and to the choice of world axes.
        Bounds2 model_projection;
        for (const Vec3& vertex : mesh.vertices)
        {
            const Vec3 local = frame.transposeMultiply(vertex - support.origin);
            model_projection.lower = model_projection.lower.cwiseMin(
                {local.x(), local.y()});
            model_projection.upper = model_projection.upper.cwiseMax(
                {local.x(), local.y()});
        }
        const double support_width = std::max(
            support_projection.upper.x() - support_projection.lower.x(), 0.0);
        const double support_height = std::max(
            support_projection.upper.y() - support_projection.lower.y(), 0.0);
        const double model_width = std::max(
            model_projection.upper.x() - model_projection.lower.x(), 0.0);
        const double model_height = std::max(
            model_projection.upper.y() - model_projection.lower.y(), 0.0);
        const double footprint_ratio = support_width * support_height /
            std::max(model_width * model_height, dimensional_tolerance * dimensional_tolerance);
        // The footprint/depth values are retained for diagnostics and for a
        // later caller-selected policy.  Do not reject here: the final
        // directed-distance and coverage certificates are the authoritative
        // acceptance gate, and this recognizer must still discover small
        // protrusions whose support plane happens to be broad in projection.
        std::vector<bool> expanded(clusters.size(), false);
        std::queue<std::uint32_t> expansion_queue;
        for (std::uint32_t id = 0; id < clusters.size(); ++id)
        {
            if (id == support.support || id == support.partner || claimed_clusters[id] ||
                areas[id] >= 0.03 * model_surface_area) continue;
            double minimum = std::numeric_limits<double>::infinity();
            double maximum = -std::numeric_limits<double>::infinity();
            for (const auto vertex : uniqueVertices(mesh, clusters[id]))
            {
                const double distance = (mesh.vertices[vertex] - support.origin).dot(support.outward);
                minimum = std::min(minimum, distance);
                maximum = std::max(maximum, distance);
            }
            // A protruding assembly may have a flange/backing patch exactly on the
            // support plane.  It belongs to the assembly even though it has no
            // positive outward depth; excluding it leaves the flange to be split
            // into hundreds of planar primitives behind the fitted box shell.
            if (maximum < -dimensional_tolerance) continue;
            const Bounds2 projection = projectedBounds(mesh, clusters[id], support.origin, frame);
            if (projectedOverlapRatio(support_projection, projection) <= 0.0) continue;
            expanded[id] = true;
            expansion_queue.push(id);
        }
        while (!expansion_queue.empty())
        {
            const auto current = expansion_queue.front();
            expansion_queue.pop();
            for (const auto neighbor : topology_neighbors[current])
            {
                if (expanded[neighbor] || claimed_clusters[neighbor] ||
                    neighbor == support.support || neighbor == support.partner ||
                    areas[neighbor] >= 0.03 * model_surface_area) continue;
                expanded[neighbor] = true;
                expansion_queue.push(neighbor);
            }
        }
        std::vector<std::uint32_t> selected;
        std::vector<Bounds2> projections;
        for (std::uint32_t id = 0; id < clusters.size(); ++id)
        {
            if (!expanded[id]) continue;
            selected.push_back(id);
            projections.push_back(projectedBounds(mesh, clusters[id], support.origin, frame));
        }
        if (selected.empty()) continue;

        std::vector<std::size_t> parent(selected.size());
        std::iota(parent.begin(), parent.end(), 0);
        const auto find = [&](std::size_t index)
        {
            std::size_t root = index;
            while (parent[root] != root) root = parent[root];
            while (parent[index] != index)
            {
                const std::size_t next = parent[index];
                parent[index] = root;
                index = next;
            }
            return root;
        };
        std::vector<std::size_t> projection_order(selected.size());
        std::iota(projection_order.begin(), projection_order.end(), 0);
        std::sort(projection_order.begin(), projection_order.end(), [&](const auto first,
                                                                        const auto second)
        {
            if (projections[first].lower.x() != projections[second].lower.x())
                return projections[first].lower.x() < projections[second].lower.x();
            return projections[first].lower.y() < projections[second].lower.y();
        });
        for (std::size_t first_position = 0; first_position < projection_order.size();
             ++first_position)
        {
            const std::size_t first = projection_order[first_position];
            for (std::size_t second_position = first_position + 1;
                 second_position < projection_order.size(); ++second_position)
            {
                const std::size_t second = projection_order[second_position];
                if (projections[second].lower.x() >
                    projections[first].upper.x() + maximum_gap) break;
                if (projectedSeparation(projections[first], projections[second]) > maximum_gap)
                    continue;
                const std::size_t first_root = find(first);
                const std::size_t second_root = find(second);
                if (first_root != second_root) parent[second_root] = first_root;
            }
        }
        std::unordered_map<std::size_t, std::vector<std::uint32_t>> groups;
        for (std::size_t index = 0; index < selected.size(); ++index)
            groups[find(index)].push_back(selected[index]);
        for (const auto& [root, group] : groups)
        {
            (void)root;
            if (group.size() < 4) continue;
            std::vector<std::uint32_t> faces;
            for (const auto cluster_id : group)
                faces.insert(faces.end(), clusters[cluster_id].begin(), clusters[cluster_id].end());
            if (faces.size() < 32) continue;
            double group_minimum = std::numeric_limits<double>::infinity();
            double group_maximum = -std::numeric_limits<double>::infinity();
            for (const auto vertex : uniqueVertices(mesh, faces))
            {
                const double distance = (mesh.vertices[vertex] - support.origin).dot(support.outward);
                group_minimum = std::min(group_minimum, distance);
                group_maximum = std::max(group_maximum, distance);
            }
            const double maximum_group_depth = std::max(
                4.0 * support.outward_depth,
                options.protrusion_max_inward_relative * model_diagonal);
            if (group_maximum - group_minimum > maximum_group_depth) continue;
            BoxFit box = fitBox(mesh, uniqueVertices(mesh, faces));
            const Vec3 box_extent = box.half_size * 2.0;
            if (std::max({box_extent.x(), box_extent.y(), box_extent.z()}) >
                options.protrusion_max_extent_relative * model_diagonal)
                continue;
            int covered_axis = 0;
            for (int axis = 1; axis < 3; ++axis)
                if (std::abs(box.axes.col(axis).dot(support.outward)) >
                    std::abs(box.axes.col(covered_axis).dot(support.outward)))
                    covered_axis = axis;
            const double covered_sign = box.axes.col(covered_axis).dot(support.outward) >= 0.0
                ? -1.0 : 1.0;
            // Compare only the exposed Box shell with source triangles close to
            // that same shell. Interior and opposite-side triangles must not
            // cancel a large unsupported face. The final ratio remains global.
            const double shell_added_area = boxShellAddedArea(
                mesh, faces, box, covered_axis, covered_sign,
                std::max(maximum_gap, dimensional_tolerance * 8.0));
            const double area_excess = shell_added_area /
                std::max(model_surface_area, 1.0e-30);
            if (minimum_candidate_area_excess_ratio < 0.0)
                minimum_candidate_area_excess_ratio = area_excess;
            else
                minimum_candidate_area_excess_ratio = std::min(
                    minimum_candidate_area_excess_ratio, area_excess);
            if (area_excess > options.protrusion_max_area_excess_ratio) continue;
            const Face& support_face =
                mesh.faces[clusters[support.support].front()];
            result.push_back(
                {box, std::move(faces), group, covered_axis, covered_sign,
                 true, mesh.vertices[support_face[0]], support.outward,
                 footprint_ratio, group_maximum - group_minimum});
            for (const auto cluster_id : group) claimed_clusters[cluster_id] = true;
        }
    }
    return result;
}

double signedArea(const std::vector<Vec2>& polygon)
{
    double result = 0.0;
    for (std::size_t index = 0; index < polygon.size(); ++index)
    {
        const Vec2& first = polygon[index];
        const Vec2& second = polygon[(index + 1) % polygon.size()];
        result += first.x() * second.y() - first.y() * second.x();
    }
    return 0.5 * result;
}

std::vector<Vec2> convexHull(std::vector<Vec2> points, const double tolerance)
{
    std::sort(points.begin(), points.end(), [](const Vec2& first, const Vec2& second)
    {
        return first.x() < second.x() || (first.x() == second.x() && first.y() < second.y());
    });
    points.erase(std::unique(points.begin(), points.end(), [tolerance](const Vec2& first,
                                                                      const Vec2& second)
    {
        return (first - second).norm() <= tolerance;
    }), points.end());
    if (points.size() <= 3) return points;
    const auto cross = [](const Vec2& origin, const Vec2& first, const Vec2& second)
    {
        const Vec2 a = first - origin;
        const Vec2 b = second - origin;
        return a.x() * b.y() - a.y() * b.x();
    };
    std::vector<Vec2> hull;
    for (const Vec2& point : points)
    {
        while (hull.size() >= 2 &&
               cross(hull[hull.size() - 2], hull.back(), point) <= tolerance)
            hull.pop_back();
        hull.push_back(point);
    }
    const std::size_t lower_size = hull.size();
    for (auto iterator = points.rbegin() + 1; iterator != points.rend(); ++iterator)
    {
        while (hull.size() > lower_size &&
               cross(hull[hull.size() - 2], hull.back(), *iterator) <= tolerance)
            hull.pop_back();
        hull.push_back(*iterator);
    }
    hull.pop_back();
    return hull;
}

std::vector<std::vector<std::uint32_t>> boundaryLoops(
    const Mesh& mesh,
    const std::vector<std::uint32_t>& faces)
{
    std::unordered_map<std::uint64_t, int> counts;
    for (const auto face_id : faces)
    {
        const Face& face = mesh.faces[face_id];
        for (int edge = 0; edge < 3; ++edge)
            ++counts[edgeKey(face[edge], face[(edge + 1) % 3])];
    }
    std::unordered_map<std::uint32_t, std::vector<std::uint32_t>> adjacency;
    for (const auto& [key, count] : counts)
    {
        if (count != 1) continue;
        const auto first = static_cast<std::uint32_t>(key >> 32U);
        const auto second = static_cast<std::uint32_t>(key);
        adjacency[first].push_back(second);
        adjacency[second].push_back(first);
    }
    std::unordered_set<std::uint64_t> visited;
    std::vector<std::vector<std::uint32_t>> loops;
    for (const auto& [start, neighbors] : adjacency)
    {
        for (const auto initial_next : neighbors)
        {
            if (visited.contains(edgeKey(start, initial_next))) continue;
            std::vector<std::uint32_t> loop{start};
            std::uint32_t previous = start;
            std::uint32_t current = initial_next;
            while (current != start && loop.size() <= adjacency.size() + 1)
            {
                loop.push_back(current);
                visited.insert(edgeKey(previous, current));
                const auto& next_candidates = adjacency[current];
                const auto iterator = std::find_if(next_candidates.begin(), next_candidates.end(),
                    [&](const auto candidate) { return candidate != previous; });
                if (iterator == next_candidates.end()) break;
                previous = current;
                current = *iterator;
            }
            if (current == start && loop.size() >= 3)
            {
                visited.insert(edgeKey(previous, current));
                loops.push_back(std::move(loop));
            }
        }
    }
    return loops;
}

std::vector<std::vector<std::uint32_t>> boundaryLoopsApproximate(
    const Mesh& mesh,
    const std::vector<std::uint32_t>& faces,
    const double tolerance)
{
    if (tolerance <= 0.0) return boundaryLoops(mesh, faces);
    struct CellKey
    {
        std::array<std::int64_t, 3> value{};
        bool operator<(const CellKey& other) const { return value < other.value; }
    };
    const auto cellFor = [&](const Vec3& point)
    {
        CellKey key;
        for (int axis = 0; axis < 3; ++axis)
            key.value[axis] = static_cast<std::int64_t>(
                std::floor(point[axis] / tolerance));
        return key;
    };
    std::map<CellKey, std::vector<std::uint32_t>> cells;
    std::vector<Vec3> canonical_points;
    std::vector<std::uint32_t> representatives;
    std::unordered_map<std::uint32_t, std::uint32_t> remap;
    remap.reserve(faces.size() * 3);
    for (const auto face_id : faces)
        for (const auto vertex : mesh.faces[face_id])
        {
            if (remap.contains(vertex)) continue;
            const CellKey center = cellFor(mesh.vertices[vertex]);
            std::uint32_t selected = std::numeric_limits<std::uint32_t>::max();
            for (int dx = -1; dx <= 1 && selected ==
                     std::numeric_limits<std::uint32_t>::max(); ++dx)
                for (int dy = -1; dy <= 1 && selected ==
                         std::numeric_limits<std::uint32_t>::max(); ++dy)
                    for (int dz = -1; dz <= 1 && selected ==
                             std::numeric_limits<std::uint32_t>::max(); ++dz)
                    {
                        CellKey neighbor = center;
                        neighbor.value[0] += dx;
                        neighbor.value[1] += dy;
                        neighbor.value[2] += dz;
                        const auto found = cells.find(neighbor);
                        if (found == cells.end()) continue;
                        for (const auto candidate : found->second)
                            if ((canonical_points[candidate] -
                                 mesh.vertices[vertex]).norm() <= tolerance)
                            {
                                selected = candidate;
                                break;
                            }
                    }
            if (selected == std::numeric_limits<std::uint32_t>::max())
            {
                selected = static_cast<std::uint32_t>(canonical_points.size());
                canonical_points.push_back(mesh.vertices[vertex]);
                representatives.push_back(vertex);
                cells[center].push_back(selected);
            }
            remap.emplace(vertex, selected);
        }

    std::unordered_map<std::uint64_t, int> counts;
    for (const auto face_id : faces)
    {
        const Face& face = mesh.faces[face_id];
        for (int edge = 0; edge < 3; ++edge)
            ++counts[edgeKey(remap.at(face[edge]),
                             remap.at(face[(edge + 1) % 3]))];
    }
    std::unordered_map<std::uint32_t, std::vector<std::uint32_t>> adjacency;
    for (const auto& [key, count] : counts)
    {
        if (count != 1) continue;
        const auto first = static_cast<std::uint32_t>(key >> 32U);
        const auto second = static_cast<std::uint32_t>(key);
        if (first == second) continue;
        adjacency[first].push_back(second);
        adjacency[second].push_back(first);
    }
    std::unordered_set<std::uint64_t> visited;
    std::vector<std::vector<std::uint32_t>> loops;
    for (const auto& [start, neighbors] : adjacency)
        for (const auto initial_next : neighbors)
        {
            if (visited.contains(edgeKey(start, initial_next))) continue;
            std::vector<std::uint32_t> loop{representatives[start]};
            std::uint32_t previous = start;
            std::uint32_t current = initial_next;
            while (current != start && loop.size() <= adjacency.size() + 1)
            {
                loop.push_back(representatives[current]);
                visited.insert(edgeKey(previous, current));
                const auto& next_candidates = adjacency[current];
                const auto iterator = std::find_if(
                    next_candidates.begin(), next_candidates.end(),
                    [&](const auto candidate) { return candidate != previous; });
                if (iterator == next_candidates.end()) break;
                previous = current;
                current = *iterator;
            }
            if (current == start && loop.size() >= 3)
            {
                visited.insert(edgeKey(previous, current));
                loops.push_back(std::move(loop));
            }
        }
    return loops;
}

std::vector<std::vector<std::uint32_t>> faceComponentsFromList(
    const Mesh& mesh, const std::vector<std::uint32_t>& faces);

std::vector<std::vector<std::uint32_t>> faceComponentsFromListApproximate(
    const Mesh& mesh,
    const std::vector<std::uint32_t>& faces,
    const double tolerance)
{
    if (faces.empty()) return {};
    if (tolerance <= 0.0) return faceComponentsFromList(mesh, faces);
    struct CellKey
    {
        std::array<std::int64_t, 3> value{};
        bool operator<(const CellKey& other) const { return value < other.value; }
    };
    const auto cellFor = [&](const Vec3& point)
    {
        CellKey key;
        for (int axis = 0; axis < 3; ++axis)
            key.value[axis] = static_cast<std::int64_t>(
                std::floor(point[axis] / tolerance));
        return key;
    };
    std::map<CellKey, std::vector<std::uint32_t>> cells;
    std::vector<Vec3> canonical_points;
    std::unordered_map<std::uint32_t, std::uint32_t> remap;
    remap.reserve(faces.size() * 3);
    for (const auto face_id : faces)
        for (const auto vertex : mesh.faces[face_id])
        {
            if (remap.contains(vertex)) continue;
            const CellKey center = cellFor(mesh.vertices[vertex]);
            std::uint32_t selected = std::numeric_limits<std::uint32_t>::max();
            for (int dx = -1; dx <= 1 && selected ==
                     std::numeric_limits<std::uint32_t>::max(); ++dx)
                for (int dy = -1; dy <= 1 && selected ==
                         std::numeric_limits<std::uint32_t>::max(); ++dy)
                    for (int dz = -1; dz <= 1 && selected ==
                             std::numeric_limits<std::uint32_t>::max(); ++dz)
                    {
                        CellKey neighbor = center;
                        neighbor.value[0] += dx;
                        neighbor.value[1] += dy;
                        neighbor.value[2] += dz;
                        const auto found = cells.find(neighbor);
                        if (found == cells.end()) continue;
                        for (const auto candidate : found->second)
                            if ((canonical_points[candidate] -
                                 mesh.vertices[vertex]).norm() <= tolerance)
                            {
                                selected = candidate;
                                break;
                            }
                    }
            if (selected == std::numeric_limits<std::uint32_t>::max())
            {
                selected = static_cast<std::uint32_t>(canonical_points.size());
                canonical_points.push_back(mesh.vertices[vertex]);
                cells[center].push_back(selected);
            }
            remap.emplace(vertex, selected);
        }
    std::unordered_map<std::uint64_t, std::vector<std::uint32_t>> edge_faces;
    edge_faces.reserve(faces.size() * 3);
    for (const auto face_id : faces)
    {
        const Face& face = mesh.faces[face_id];
        for (int edge = 0; edge < 3; ++edge)
            edge_faces[edgeKey(remap.at(face[edge]),
                remap.at(face[(edge + 1) % 3]))].push_back(face_id);
    }
    std::unordered_map<std::uint32_t, std::vector<std::uint32_t>> adjacency;
    adjacency.reserve(faces.size());
    for (const auto& [edge, incident] : edge_faces)
    {
        (void)edge;
        for (const auto first : incident)
            for (const auto second : incident)
                if (first != second) adjacency[first].push_back(second);
    }
    std::unordered_set<std::uint32_t> remaining(faces.begin(), faces.end());
    std::vector<std::vector<std::uint32_t>> result;
    while (!remaining.empty())
    {
        const auto seed = *remaining.begin();
        remaining.erase(seed);
        result.emplace_back();
        std::queue<std::uint32_t> queue;
        queue.push(seed);
        while (!queue.empty())
        {
            const auto current = queue.front();
            queue.pop();
            result.back().push_back(current);
            const auto found = adjacency.find(current);
            if (found == adjacency.end()) continue;
            for (const auto neighbor : found->second)
                if (remaining.erase(neighbor) != 0) queue.push(neighbor);
        }
    }
    return result;
}

std::size_t faceComponentCount(const Mesh& mesh,
                               const std::vector<std::uint32_t>& faces)
{
    std::unordered_map<std::uint64_t, std::vector<std::size_t>> edge_faces;
    edge_faces.reserve(faces.size() * 3);
    for (std::size_t local_id = 0; local_id < faces.size(); ++local_id)
    {
        const Face& face = mesh.faces[faces[local_id]];
        for (int edge = 0; edge < 3; ++edge)
            edge_faces[edgeKey(face[edge], face[(edge + 1) % 3])].push_back(local_id);
    }
    std::vector<std::vector<std::size_t>> adjacency(faces.size());
    for (const auto& [key, incident] : edge_faces)
    {
        (void)key;
        for (const auto first : incident)
            for (const auto second : incident)
                if (first != second) adjacency[first].push_back(second);
    }
    std::vector<bool> visited(faces.size(), false);
    std::size_t components = 0;
    for (std::size_t seed = 0; seed < faces.size(); ++seed)
    {
        if (visited[seed]) continue;
        ++components;
        std::queue<std::size_t> queue;
        queue.push(seed);
        visited[seed] = true;
        while (!queue.empty())
        {
            const auto current = queue.front();
            queue.pop();
            for (const auto neighbor : adjacency[current])
            {
                if (visited[neighbor]) continue;
                visited[neighbor] = true;
                queue.push(neighbor);
            }
        }
    }
    return components;
}

std::vector<std::vector<std::uint32_t>> faceComponentsFromList(
    const Mesh& mesh,
    const std::vector<std::uint32_t>& faces)
{
    std::unordered_map<std::uint64_t, std::vector<std::uint32_t>> edge_faces;
    edge_faces.reserve(faces.size() * 3);
    for (const auto face_id : faces)
    {
        const Face& face = mesh.faces[face_id];
        for (int edge = 0; edge < 3; ++edge)
            edge_faces[edgeKey(face[edge], face[(edge + 1) % 3])].push_back(face_id);
    }
    std::unordered_map<std::uint32_t, std::vector<std::uint32_t>> adjacency;
    adjacency.reserve(faces.size());
    for (const auto& [key, incident] : edge_faces)
    {
        (void)key;
        for (const auto first : incident)
            for (const auto second : incident)
                if (first != second) adjacency[first].push_back(second);
    }
    std::unordered_set<std::uint32_t> remaining(faces.begin(), faces.end());
    std::vector<std::vector<std::uint32_t>> result;
    while (!remaining.empty())
    {
        const auto seed = *remaining.begin();
        remaining.erase(seed);
        result.emplace_back();
        std::queue<std::uint32_t> queue;
        queue.push(seed);
        while (!queue.empty())
        {
            const auto current = queue.front();
            queue.pop();
            result.back().push_back(current);
            const auto found = adjacency.find(current);
            if (found == adjacency.end()) continue;
            for (const auto neighbor : found->second)
                if (remaining.erase(neighbor) != 0) queue.push(neighbor);
        }
    }
    return result;
}

std::vector<std::vector<std::uint32_t>> faceComponents(
    const Mesh& mesh,
    const std::vector<bool>& included_faces)
{
    std::unordered_map<std::uint64_t, std::vector<std::uint32_t>> edge_faces;
    edge_faces.reserve(mesh.faces.size() * 3);
    for (std::uint32_t id = 0; id < mesh.faces.size(); ++id)
    {
        if (!included_faces[id]) continue;
        const Face& face = mesh.faces[id];
        for (int edge = 0; edge < 3; ++edge)
            edge_faces[edgeKey(face[edge], face[(edge + 1) % 3])].push_back(id);
    }
    std::vector<std::vector<std::uint32_t>> adjacency(mesh.faces.size());
    for (const auto& [key, incident] : edge_faces)
    {
        (void)key;
        for (const auto first : incident)
            for (const auto second : incident)
                if (first != second) adjacency[first].push_back(second);
    }
    std::vector<bool> visited(mesh.faces.size(), false);
    std::vector<std::vector<std::uint32_t>> result;
    for (std::uint32_t seed = 0; seed < mesh.faces.size(); ++seed)
    {
        if (!included_faces[seed] || visited[seed]) continue;
        result.emplace_back();
        std::queue<std::uint32_t> queue;
        queue.push(seed);
        visited[seed] = true;
        while (!queue.empty())
        {
            const auto current = queue.front();
            queue.pop();
            result.back().push_back(current);
            for (const auto neighbor : adjacency[current])
            {
                if (visited[neighbor]) continue;
                visited[neighbor] = true;
                queue.push(neighbor);
            }
        }
    }
    return result;
}

std::vector<std::vector<std::uint32_t>> faceComponentsApproximate(
    const Mesh& mesh,
    const std::vector<bool>& included_faces,
    const double tolerance)
{
    if (tolerance <= 0.0) return faceComponents(mesh, included_faces);
    struct CellKey
    {
        std::array<std::int64_t, 3> value{};
        bool operator<(const CellKey& other) const { return value < other.value; }
    };
    const auto cellFor = [&](const Vec3& point)
    {
        CellKey key;
        for (int axis = 0; axis < 3; ++axis)
            key.value[axis] = static_cast<std::int64_t>(
                std::floor(point[axis] / tolerance));
        return key;
    };
    std::map<CellKey, std::vector<std::uint32_t>> cells;
    std::vector<Vec3> canonical_points;
    std::vector<std::uint32_t> remap(mesh.vertices.size());
    for (std::uint32_t vertex = 0; vertex < mesh.vertices.size(); ++vertex)
    {
        const CellKey center = cellFor(mesh.vertices[vertex]);
        std::uint32_t selected = std::numeric_limits<std::uint32_t>::max();
        for (int dx = -1; dx <= 1 && selected ==
                 std::numeric_limits<std::uint32_t>::max(); ++dx)
            for (int dy = -1; dy <= 1 && selected ==
                     std::numeric_limits<std::uint32_t>::max(); ++dy)
                for (int dz = -1; dz <= 1 && selected ==
                         std::numeric_limits<std::uint32_t>::max(); ++dz)
                {
                    CellKey neighbor = center;
                    neighbor.value[0] += dx;
                    neighbor.value[1] += dy;
                    neighbor.value[2] += dz;
                    const auto found = cells.find(neighbor);
                    if (found == cells.end()) continue;
                    for (const auto candidate : found->second)
                        if ((canonical_points[candidate] - mesh.vertices[vertex]).norm() <=
                            tolerance)
                        {
                            selected = candidate;
                            break;
                        }
                }
        if (selected == std::numeric_limits<std::uint32_t>::max())
        {
            selected = static_cast<std::uint32_t>(canonical_points.size());
            canonical_points.push_back(mesh.vertices[vertex]);
            cells[center].push_back(selected);
        }
        remap[vertex] = selected;
    }

    std::unordered_map<std::uint64_t, std::vector<std::uint32_t>> edge_faces;
    edge_faces.reserve(mesh.faces.size() * 3);
    for (std::uint32_t id = 0; id < mesh.faces.size(); ++id)
    {
        if (!included_faces[id]) continue;
        const Face& face = mesh.faces[id];
        for (int edge = 0; edge < 3; ++edge)
            edge_faces[edgeKey(remap[face[edge]],
                               remap[face[(edge + 1) % 3]])].push_back(id);
    }
    std::vector<std::vector<std::uint32_t>> adjacency(mesh.faces.size());
    for (const auto& [key, incident] : edge_faces)
    {
        (void)key;
        for (const auto first : incident)
            for (const auto second : incident)
                if (first != second) adjacency[first].push_back(second);
    }
    std::vector<bool> visited(mesh.faces.size(), false);
    std::vector<std::vector<std::uint32_t>> result;
    for (std::uint32_t seed = 0; seed < mesh.faces.size(); ++seed)
    {
        if (!included_faces[seed] || visited[seed]) continue;
        result.emplace_back();
        std::queue<std::uint32_t> queue;
        queue.push(seed);
        visited[seed] = true;
        while (!queue.empty())
        {
            const auto current = queue.front();
            queue.pop();
            result.back().push_back(current);
            for (const auto neighbor : adjacency[current])
                if (!visited[neighbor])
                {
                    visited[neighbor] = true;
                    queue.push(neighbor);
                }
        }
    }
    return result;
}

std::vector<Vec2> simplifyPolygon(std::vector<Vec2> polygon, const double tolerance)
{
    bool changed = true;
    while (changed && polygon.size() > 3)
    {
        changed = false;
        for (std::size_t index = 0; index < polygon.size(); ++index)
        {
            const Vec2 previous = polygon[(index + polygon.size() - 1) % polygon.size()];
            const Vec2 current = polygon[index];
            const Vec2 next = polygon[(index + 1) % polygon.size()];
            const Vec2 first = current - previous;
            const Vec2 second = next - current;
            const double cross = std::abs(first.x() * second.y() - first.y() * second.x());
            if (cross <= tolerance * std::max(first.norm() + second.norm(), 1.0))
            {
                polygon.erase(polygon.begin() + static_cast<std::ptrdiff_t>(index));
                changed = true;
                break;
            }
        }
    }
    return polygon;
}

double cross2(const Vec2& first, const Vec2& second)
{
    return first.x() * second.y() - first.y() * second.x();
}

bool pointOnSegment(const Vec2& point, const Vec2& first, const Vec2& second,
                    const double tolerance)
{
    if (std::abs(cross2(second - first, point - first)) >
        tolerance * std::max((second - first).norm(), 1.0)) return false;
    return point.x() >= std::min(first.x(), second.x()) - tolerance &&
           point.x() <= std::max(first.x(), second.x()) + tolerance &&
           point.y() >= std::min(first.y(), second.y()) - tolerance &&
           point.y() <= std::max(first.y(), second.y()) + tolerance;
}

bool pointInPolygon(const Vec2& point, const std::vector<Vec2>& polygon,
                    const double tolerance)
{
    bool inside = false;
    for (std::size_t index = 0, previous = polygon.size() - 1;
         index < polygon.size(); previous = index++)
    {
        const Vec2& first = polygon[previous];
        const Vec2& second = polygon[index];
        if (pointOnSegment(point, first, second, tolerance)) return true;
        const bool crosses = (first.y() > point.y()) != (second.y() > point.y());
        if (!crosses) continue;
        const double x = first.x() + (point.y() - first.y()) *
                         (second.x() - first.x()) / (second.y() - first.y());
        if (x >= point.x() - tolerance) inside = !inside;
    }
    return inside;
}

bool properSegmentIntersection(const Vec2& a, const Vec2& b,
                               const Vec2& c, const Vec2& d,
                               const double tolerance)
{
    const double first = cross2(b - a, c - a);
    const double second = cross2(b - a, d - a);
    const double third = cross2(d - c, a - c);
    const double fourth = cross2(d - c, b - c);
    return ((first > tolerance && second < -tolerance) ||
            (first < -tolerance && second > tolerance)) &&
           ((third > tolerance && fourth < -tolerance) ||
            (third < -tolerance && fourth > tolerance));
}

bool simplePolygon(const std::vector<Vec2>& polygon, const double tolerance)
{
    for (std::size_t first = 0; first < polygon.size(); ++first)
    {
        const std::size_t first_next = (first + 1) % polygon.size();
        for (std::size_t second = first + 1; second < polygon.size(); ++second)
        {
            const std::size_t second_next = (second + 1) % polygon.size();
            if (first == second || first_next == second || second_next == first) continue;
            if (properSegmentIntersection(polygon[first], polygon[first_next],
                                          polygon[second], polygon[second_next], tolerance))
                return false;
        }
    }
    return true;
}

bool lineIntersection(const Vec2& first_start, const Vec2& first_end,
                      const Vec2& second_start, const Vec2& second_end,
                      Vec2& intersection)
{
    const Vec2 first = first_end - first_start;
    const Vec2 second = second_end - second_start;
    const double denominator = cross2(first, second);
    if (std::abs(denominator) <= 1.0e-14 *
        std::max(first.norm() * second.norm(), 1.0)) return false;
    const double parameter = cross2(second_start - first_start, second) / denominator;
    intersection = first_start + first * parameter;
    return true;
}

std::vector<Vec2> fillShortBoundaryVoids(
    std::vector<Vec2> polygon,
    const double tolerance,
    const double allowed_excess_ratio,
    const double model_surface_area,
    const double maximum_boundary_deviation,
    std::size_t& filled_count,
    double& filled_area)
{
    if (polygon.size() < 4) return polygon;
    const double area_budget = allowed_excess_ratio * model_surface_area;
    while (polygon.size() >= 4)
    {
        Vec2 lower = Vec2::Constant(std::numeric_limits<double>::infinity());
        Vec2 upper = Vec2::Constant(-std::numeric_limits<double>::infinity());
        for (const Vec2& point : polygon)
        {
            lower = lower.cwiseMin(point);
            upper = upper.cwiseMax(point);
        }
        const double maximum_edge = (upper - lower).norm() * 0.03;
        bool accepted = false;
        for (std::size_t edge = 0; edge < polygon.size(); ++edge)
        {
            const std::size_t next = (edge + 1) % polygon.size();
            const std::size_t previous = (edge + polygon.size() - 1) % polygon.size();
            const std::size_t following = (edge + 2) % polygon.size();
            const double edge_length = (polygon[next] - polygon[edge]).norm();
            if (edge_length <= tolerance || edge_length > maximum_edge) continue;
            if (std::min((polygon[edge] - polygon[previous]).norm(),
                         (polygon[following] - polygon[next]).norm()) < 3.0 * edge_length)
                continue;
            Vec2 intersection;
            if (!lineIntersection(polygon[previous], polygon[edge],
                                  polygon[next], polygon[following], intersection)) continue;
            if (std::max((intersection - polygon[edge]).norm(),
                         (intersection - polygon[next]).norm()) >
                std::min(2.0 * edge_length,
                         maximum_boundary_deviation + tolerance))
                continue;

            std::vector<Vec2> candidate;
            candidate.reserve(polygon.size() - 1);
            candidate.push_back(intersection);
            for (std::size_t offset = 2; offset < polygon.size(); ++offset)
                candidate.push_back(polygon[(edge + offset) % polygon.size()]);
            if (!simplePolygon(candidate, tolerance)) continue;
            const double added_area = std::abs(signedArea(candidate)) -
                                      std::abs(signedArea(polygon));
            if (added_area <= tolerance * tolerance ||
                added_area > area_budget) continue;
            const bool contains_source = std::all_of(
                polygon.begin(), polygon.end(), [&](const Vec2& point)
                { return pointInPolygon(point, candidate, tolerance); });
            if (!contains_source) continue;
            polygon = std::move(candidate);
            ++filled_count;
            filled_area += added_area;
            accepted = true;
            break;
        }
        if (!accepted) break;
    }
    return polygon;
}

std::vector<Vec2> absorbTinyPlanarDetails(
    std::vector<Vec2> polygon,
    const double tolerance,
    const double maximum_added_area_ratio,
    const double model_surface_area,
    const std::uint32_t maximum_removed_vertices,
    std::size_t& absorbed_count,
    double& added_area)
{
    if (polygon.size() < 4 || maximum_added_area_ratio <= 0.0 ||
        maximum_removed_vertices == 0) return polygon;
    const std::vector<Vec2> original = polygon;
    const double original_area = std::abs(signedArea(original));
    const double area_budget = model_surface_area * maximum_added_area_ratio;
    while (polygon.size() > 3)
    {
        const double current_signed_area = signedArea(polygon);
        std::vector<Vec2> best;
        std::size_t best_removed = 0;
        double best_added_area = std::numeric_limits<double>::infinity();
        const std::size_t maximum_run = std::min<std::size_t>(
            maximum_removed_vertices, polygon.size() - 3);
        for (std::size_t removed = 1; removed <= maximum_run; ++removed)
        {
            for (std::size_t start = 0; start < polygon.size(); ++start)
            {
                const std::size_t count = polygon.size();
                const std::size_t previous = (start + count - 1) % count;
                const std::size_t next = (start + removed) % count;

                // Replacing one consecutive boundary chain introduces exactly
                // one new edge. Compute its area delta first; almost all CAD
                // boundary positions are convex and fail this constant-time
                // conservative-fill test without constructing a polygon.
                double removed_cross = cross2(polygon[previous], polygon[start]);
                for (std::size_t offset = 0; offset < removed; ++offset)
                {
                    const std::size_t edge = (start + offset) % count;
                    removed_cross += cross2(
                        polygon[edge], polygon[(edge + 1) % count]);
                }
                const double candidate_signed_area = current_signed_area +
                    0.5 * (cross2(polygon[previous], polygon[next]) -
                           removed_cross);
                const double candidate_added_area =
                    std::abs(candidate_signed_area) - original_area;
                if (candidate_signed_area * current_signed_area <= 0.0 ||
                    candidate_added_area <= tolerance * tolerance ||
                    candidate_added_area > area_budget) continue;

                // The input boundary is already certified simple. All retained
                // edges are unchanged, so only the shortcut can introduce a
                // new proper intersection. This is equivalent to rerunning the
                // quadratic all-edge simplePolygon test on every candidate.
                bool intersects = false;
                for (std::size_t edge = 0; edge < count; ++edge)
                {
                    const std::size_t from_previous =
                        (edge + count - previous) % count;
                    if (from_previous <= removed) continue;
                    if (properSegmentIntersection(
                            polygon[previous], polygon[next], polygon[edge],
                            polygon[(edge + 1) % count], tolerance))
                    {
                        intersects = true;
                        break;
                    }
                }
                if (intersects) continue;

                std::vector<Vec2> candidate;
                candidate.reserve(count - removed);
                for (std::size_t offset = removed; offset < count; ++offset)
                    candidate.push_back(polygon[(start + offset) % count]);

                // Retained vertices lie on candidate edges and therefore pass
                // the previous all-source-vertex containment test trivially.
                // Only the removed run needs an explicit point-in-polygon test.
                bool contains_source = true;
                for (std::size_t offset = 0; offset < removed; ++offset)
                    contains_source &= pointInPolygon(
                        polygon[(start + offset) % count], candidate, tolerance);
                if (!contains_source) continue;
                if (removed > best_removed ||
                    (removed == best_removed && candidate_added_area < best_added_area))
                {
                    best = std::move(candidate);
                    best_removed = removed;
                    best_added_area = candidate_added_area;
                }
            }
        }
        if (best.empty()) break;
        polygon = std::move(best);
        ++absorbed_count;
        added_area = std::max(added_area, best_added_area);
    }
    return polygon;
}

bool pointInTriangle(const Vec2& point, const Vec2& a, const Vec2& b, const Vec2& c)
{
    const auto cross = [](const Vec2& u, const Vec2& v)
    { return u.x() * v.y() - u.y() * v.x(); };
    const double first = cross(b - a, point - a);
    const double second = cross(c - b, point - b);
    const double third = cross(a - c, point - c);
    return (first >= -1.0e-12 && second >= -1.0e-12 && third >= -1.0e-12) ||
           (first <= 1.0e-12 && second <= 1.0e-12 && third <= 1.0e-12);
}

std::vector<std::array<std::size_t, 3>> triangulatePolygon(std::vector<Vec2> polygon)
{
    std::vector<std::array<std::size_t, 3>> triangles;
    std::vector<std::size_t> indices(polygon.size());
    std::iota(indices.begin(), indices.end(), 0);
    const bool positive = signedArea(polygon) > 0.0;
    while (indices.size() > 3)
    {
        bool clipped = false;
        for (std::size_t offset = 0; offset < indices.size(); ++offset)
        {
            const std::size_t previous = indices[(offset + indices.size() - 1) % indices.size()];
            const std::size_t current = indices[offset];
            const std::size_t next = indices[(offset + 1) % indices.size()];
            const Vec2 first = polygon[current] - polygon[previous];
            const Vec2 second = polygon[next] - polygon[current];
            const double cross = first.x() * second.y() - first.y() * second.x();
            if ((positive && cross <= 0.0) || (!positive && cross >= 0.0)) continue;
            bool contains = false;
            for (const auto candidate : indices)
            {
                if (candidate == previous || candidate == current || candidate == next) continue;
                contains |= pointInTriangle(polygon[candidate], polygon[previous],
                                            polygon[current], polygon[next]);
            }
            if (contains) continue;
            triangles.push_back({previous, current, next});
            indices.erase(indices.begin() + static_cast<std::ptrdiff_t>(offset));
            clipped = true;
            break;
        }
        if (!clipped) break;
    }
    if (indices.size() == 3) triangles.push_back({indices[0], indices[1], indices[2]});
    return triangles;
}

std::vector<OutputPrimitive> exactTrianglePrimitives(
    const Mesh& mesh, const std::vector<std::uint32_t>& faces)
{
    std::vector<OutputPrimitive> result;
    result.reserve(faces.size());
    for (const auto face_id : faces)
    {
        Primitive triangle;
        triangle.kind = Kind::Triangle;
        for (int corner = 0; corner < 3; ++corner)
            triangle.triangle[corner] = mesh.vertices[mesh.faces[face_id][corner]];
        result.push_back({triangle, {face_id}});
    }
    return result;
}

Primitive polygonPrimitive(const std::vector<Vec2>& boundary,
                           const Vec3& origin,
                           const Mat3& frame)
{
    Primitive polygon;
    polygon.kind = Kind::Polygon;
    polygon.polygon.reserve(boundary.size());
    for (const Vec2& point : boundary)
        polygon.polygon.push_back(origin + frame.col(0) * point.x() +
                                  frame.col(1) * point.y());
    return polygon;
}

Primitive polygonPrimitive(std::initializer_list<Vec3> boundary)
{
    Primitive polygon;
    polygon.kind = Kind::Polygon;
    polygon.polygon.assign(boundary);
    return polygon;
}

std::vector<OutputPrimitive> classifyFinalRegion(
    const Mesh& mesh,
    const std::vector<std::uint32_t>& faces,
    const PrimitiveMeshAnalysisOptions& options,
    const double allowed_excess_ratio,
    const double model_surface_area,
    std::size_t& filled_holes,
    std::size_t& filled_boundary_voids,
    double& filled_boundary_void_area)
{
    const Face& seed = mesh.faces[faces.front()];
    const Vec3 normal = (mesh.vertices[seed[1]] - mesh.vertices[seed[0]])
                            .cross(mesh.vertices[seed[2]] - mesh.vertices[seed[0]])
                            .normalized();
    const Vec3 origin = mesh.vertices[seed[0]];
    Vec3 region_lower = Vec3::Constant(std::numeric_limits<double>::infinity());
    Vec3 region_upper = Vec3::Constant(-std::numeric_limits<double>::infinity());
    for (const auto vertex : uniqueVertices(mesh, faces))
    {
        region_lower = region_lower.cwiseMin(mesh.vertices[vertex]);
        region_upper = region_upper.cwiseMax(mesh.vertices[vertex]);
    }
    const double planar_tolerance = std::max(
        (region_upper - region_lower).norm() * options.coplanar_relative_tolerance,
        1.0e-10);
    bool planar = true;
    for (const auto vertex : uniqueVertices(mesh, faces))
        planar &= std::abs((mesh.vertices[vertex] - origin).dot(normal))
                  <= planar_tolerance;
    if (!planar) return exactTrianglePrimitives(mesh, faces);

    // coplanarClusters already returns geometrically connected regions. Do not
    // rebuild a full-mesh weld table for every cluster here.
    auto loops = boundaryLoopsApproximate(mesh, faces, planar_tolerance);
    if (loops.empty()) return exactTrianglePrimitives(mesh, faces);
    const Mat3 frame = [&]
    {
        Mat3 value = orthonormalFrame(normal);
        Mat3 reordered;
        reordered.col(0) = static_cast<Vec3>(value.col(1));
        reordered.col(1) = static_cast<Vec3>(value.col(2));
        reordered.col(2) = static_cast<Vec3>(value.col(0));
        return reordered;
    }();
    std::vector<Vec2> polygon;
    std::vector<double> loop_areas;
    loop_areas.reserve(loops.size());
    std::size_t outer_loop_index = 0;
    double largest_area = -1.0;
    for (std::size_t loop_index = 0; loop_index < loops.size(); ++loop_index)
    {
        const auto& loop = loops[loop_index];
        std::vector<Vec2> candidate;
        for (const auto vertex : loop)
        {
            const Vec3 local = frame.transposeMultiply(mesh.vertices[vertex] - origin);
            candidate.emplace_back(local.x(), local.y());
        }
        candidate = simplifyPolygon(std::move(candidate), planar_tolerance);
        const double area = std::abs(signedArea(candidate));
        loop_areas.push_back(area);
        if (area > largest_area)
        {
            largest_area = area;
            outer_loop_index = loop_index;
            polygon = std::move(candidate);
        }
    }
    if (polygon.size() < 3) return exactTrianglePrimitives(mesh, faces);
    const auto validSimpleBoundary = [&](const std::vector<Vec2>& candidate)
    {
        return candidate.size() >= 3 &&
               std::abs(signedArea(candidate)) >
                   planar_tolerance * planar_tolerance &&
               simplePolygon(candidate, planar_tolerance) &&
               triangulatePolygon(candidate).size() + 2 == candidate.size();
    };
    if (!validSimpleBoundary(polygon)) return exactTrianglePrimitives(mesh, faces);

    // Ordinary planar inner loops are holes in a surface patch, not 3D cavity
    // boundaries.  Closing them is conservative and does not contribute to the
    // exterior-envelope error.  Volumetric cavity preservation is handled by
    // the separate connected-void pass.
    if (allowed_excess_ratio >= 0.0)
    {
        std::size_t candidate_filled_voids = 0;
        double candidate_filled_area = 0.0;
        auto candidate = fillShortBoundaryVoids(
            polygon, planar_tolerance, allowed_excess_ratio,
            model_surface_area,
            options.maximum_open_error_distance >= 0.0
                ? options.maximum_open_error_distance
                : std::numeric_limits<double>::infinity(),
            candidate_filled_voids, candidate_filled_area);
        if (validSimpleBoundary(candidate))
        {
            polygon = std::move(candidate);
            filled_boundary_voids += candidate_filled_voids;
            filled_boundary_void_area += candidate_filled_area;
        }
    }
    {
        std::size_t candidate_filled_voids = 0;
        double candidate_filled_area = 0.0;
        auto candidate = absorbTinyPlanarDetails(
            polygon, planar_tolerance,
            options.tiny_planar_detail_area_ratio,
            model_surface_area,
            options.tiny_planar_detail_max_vertices,
            candidate_filled_voids, candidate_filled_area);
        if (validSimpleBoundary(candidate))
        {
            polygon = std::move(candidate);
            filled_boundary_voids += candidate_filled_voids;
            filled_boundary_void_area += candidate_filled_area;
        }
    }
    double source_area = 0.0;
    for (const auto face_id : faces)
    {
        const Face& face = mesh.faces[face_id];
        const Vec3 first = frame.transposeMultiply(mesh.vertices[face[1]] - origin);
        const Vec3 second = frame.transposeMultiply(mesh.vertices[face[2]] - origin);
        const Vec3 third = frame.transposeMultiply(mesh.vertices[face[0]] - origin);
        const Vec2 edge_first(first.x() - third.x(), first.y() - third.y());
        const Vec2 edge_second(second.x() - third.x(), second.y() - third.y());
        source_area += 0.5 * std::abs(edge_first.x() * edge_second.y() -
                                      edge_first.y() * edge_second.x());
    }
    const std::size_t hole_count = loops.size() > 1 ? loops.size() - 1 : 0;
    const double effective_source_area = hole_count > 0 ? largest_area : source_area;
    // Phase 1 is the hole-filled reference, not an approximation stage.  An
    // inner boundary loop may be sealed here, but an outer concavity must stay
    // in the reference and compete under the user's directed-distance limit in
    // phase 3.  Fitting the convex hull here turned every U-shaped or notched
    // planar patch with a rectangular hull into a free oversized rectangle.
    if (options.allow_polygon && polygon.size() == 4)
    {
        const Vec2 first = polygon[1] - polygon[0];
        const Vec2 second = polygon[3] - polygon[0];
        const double rectangle_area = std::abs(signedArea(polygon));
        const double rectangle_excess_ratio =
            std::max(rectangle_area - effective_source_area, 0.0) /
            std::max(model_surface_area, 1.0e-30);
        if (std::abs(first.dot(second)) <= 1.0e-8 * first.norm() * second.norm() &&
            rectangle_excess_ratio <= allowed_excess_ratio)
        {
            Primitive rectangle;
            rectangle.kind = Kind::Rectangle;
            const Vec2 first_axis = first.normalized();
            const Vec2 second_axis = second.normalized();
            rectangle.axes.col(0) = frame.col(0) * first_axis.x() + frame.col(1) * first_axis.y();
            rectangle.axes.col(1) = frame.col(0) * second_axis.x() + frame.col(1) * second_axis.y();
            rectangle.axes.col(2) = normal;
            rectangle.half_size = Vec3(first.norm() * 0.5, second.norm() * 0.5, 0.0);
            const Vec2 center = (polygon[0] + polygon[2]) * 0.5;
            rectangle.center = origin + frame.col(0) * center.x() + frame.col(1) * center.y();
            filled_holes += hole_count;
            return {{rectangle, faces}};
        }
    }

    if (options.allow_round_surfaces && polygon.size() >= 8)
    {
        Vec2 center = Vec2::Zero();
        for (const Vec2& point : polygon) center += point;
        center /= static_cast<double>(polygon.size());
        std::vector<double> radii;
        std::vector<Vec2> centered_polygon;
        centered_polygon.reserve(polygon.size());
        double perimeter = 0.0;
        for (std::size_t index = 0; index < polygon.size(); ++index)
        {
            centered_polygon.push_back(polygon[index] - center);
            radii.push_back(centered_polygon.back().norm());
            perimeter += (polygon[(index + 1) % polygon.size()] - polygon[index]).norm();
        }
        const double mean = std::accumulate(radii.begin(), radii.end(), 0.0) / radii.size();
        double maximum_radial_error = 0.0;
        for (const double radius : radii)
            maximum_radial_error = std::max(maximum_radial_error, std::abs(radius - mean));
        const double polygon_area = std::abs(signedArea(polygon));
        const double circularity = perimeter > 0.0
            ? 4.0 * std::numbers::pi * polygon_area / (perimeter * perimeter)
            : 0.0;
        const double maximum_radius = *std::max_element(
            radii.begin(), radii.end());
        const double polygon_expansion = options.round_surface_segments >= 3
            ? 1.0 / std::cos(
                std::numbers::pi /
                static_cast<double>(options.round_surface_segments))
            : std::numeric_limits<double>::infinity();
        const double conservative_radius = maximum_radius *
            polygon_expansion * (1.0 + 1.0e-12);
        const double circle_area = std::numbers::pi *
            conservative_radius * conservative_radius;
        const double circle_excess_ratio =
            std::max(circle_area - effective_source_area, 0.0) /
            std::max(model_surface_area, 1.0e-30);
        if (maximum_radial_error <= options.circle_radial_tolerance * mean &&
            ringLooksCircular(centered_polygon, options.circle_radial_tolerance) &&
            circularity >= 0.9 &&
            circle_excess_ratio <= allowed_excess_ratio)
        {
            Primitive circle;
            circle.kind = Kind::Disk;
            circle.center = origin + frame.col(0) * center.x() + frame.col(1) * center.y();
            circle.axes.col(0) = normal;
            circle.axes.col(1) = frame.col(0);
            circle.axes.col(2) = frame.col(1);
            // The emitted regular polygon is inscribed in its radius.  Expand
            // the fitted circle to the circumscribed radius so the triangle OBJ
            // covers the complete source boundary rather than only its sampled
            // vertices.
            circle.base_radius = conservative_radius;
            circle.top_radius = conservative_radius;
            circle.segments = options.round_surface_segments;
            filled_holes += hole_count;
            return {{circle, faces}};
        }
    }

    filled_holes += hole_count;
    return {{polygonPrimitive(polygon, origin, frame), faces}};
}

struct CircularLoopFit
{
    Vec2 center = Vec2::Zero();
    double minimum_radius = 0.0;
    double maximum_radius = 0.0;
    double mean_radius = 0.0;
    double area = 0.0;
    bool circular = false;
};

CircularLoopFit fitCircularLoop(const std::vector<Vec2>& loop,
                                const double radial_tolerance)
{
    CircularLoopFit fit;
    if (loop.size() < 6) return fit;
    Vec2 lower = Vec2::Constant(std::numeric_limits<double>::infinity());
    Vec2 upper = Vec2::Constant(-std::numeric_limits<double>::infinity());
    for (const Vec2& point : loop)
    {
        lower.x() = std::min(lower.x(), point.x());
        lower.y() = std::min(lower.y(), point.y());
        upper.x() = std::max(upper.x(), point.x());
        upper.y() = std::max(upper.y(), point.y());
    }
    fit.center = (lower + upper) * 0.5;
    fit.minimum_radius = std::numeric_limits<double>::infinity();
    fit.maximum_radius = 0.0;
    double sum = 0.0;
    std::vector<Vec2> centered;
    centered.reserve(loop.size());
    for (const Vec2& point : loop)
    {
        centered.push_back(point - fit.center);
        const double radius = centered.back().norm();
        fit.minimum_radius = std::min(fit.minimum_radius, radius);
        fit.maximum_radius = std::max(fit.maximum_radius, radius);
        sum += radius;
    }
    fit.mean_radius = sum / static_cast<double>(loop.size());
    fit.area = std::abs(signedArea(loop));
    fit.circular = fit.mean_radius > 0.0 &&
        fit.maximum_radius - fit.minimum_radius <=
            2.0 * radial_tolerance * fit.mean_radius &&
        ringLooksCircular(centered, 2.0 * radial_tolerance);
    return fit;
}

bool analyticPlanarSurfaceOwnsFace(
    const Mesh& mesh, const std::uint32_t face_id,
    const Primitive& surface, const double tolerance)
{
    if (face_id >= mesh.faces.size() ||
        (surface.kind != Kind::Disk && surface.kind != Kind::Annulus) ||
        surface.segments < 3)
        return false;
    // The outer regular polygon is emitted with a circumscribed radius, so its
    // inradius is the guaranteed radial coverage.  The inner annulus polygon is
    // inscribed; retaining the continuous inner radius is conservative when
    // deciding which source triangles this candidate may own.
    const double outer_inradius = surface.base_radius * std::cos(
        std::numbers::pi / static_cast<double>(surface.segments));
    const auto covered = [&](const Vec3& point)
    {
        const Vec3 local = surface.axes.transposeMultiply(point - surface.center);
        if (local.x() > tolerance) return false;
        const double radius = std::hypot(local.y(), local.z());
        if (radius > outer_inradius + tolerance) return false;
        return surface.kind != Kind::Annulus ||
               radius >= surface.inner_radius - tolerance;
    };
    const Face& face = mesh.faces[face_id];
    if (!std::all_of(face.begin(), face.end(), [&](const auto vertex)
        { return covered(mesh.vertices[vertex]); }))
        return false;
    const Vec3 centroid = (mesh.vertices[face[0]] + mesh.vertices[face[1]] +
                           mesh.vertices[face[2]]) / 3.0;
    return covered(centroid);
}

std::vector<std::vector<std::uint32_t>> approximatePlanarRegions(
    const Mesh& mesh,
    const std::vector<bool>& included_faces,
    const PrimitiveMeshAnalysisOptions& options,
    const double model_diagonal)
{
    struct FacePlane
    {
        std::uint32_t id = 0;
        Vec3 normal = Vec3::Zero();
        Vec3 centroid = Vec3::Zero();
        double doubled_area = 0.0;
    };
    std::vector<FacePlane> planes;
    planes.reserve(mesh.faces.size());
    for (std::uint32_t id = 0; id < mesh.faces.size(); ++id)
    {
        if (!included_faces[id]) continue;
        const Face& face = mesh.faces[id];
        const Vec3 cross = (mesh.vertices[face[1]] - mesh.vertices[face[0]])
            .cross(mesh.vertices[face[2]] - mesh.vertices[face[0]]);
        const double area = cross.norm();
        if (area <= 1.0e-30) continue;
        planes.push_back({id, cross / area,
            (mesh.vertices[face[0]] + mesh.vertices[face[1]] +
             mesh.vertices[face[2]]) / 3.0, area});
    }
    std::vector<std::size_t> order(planes.size());
    std::iota(order.begin(), order.end(), 0);
    std::sort(order.begin(), order.end(), [&](const auto first, const auto second)
    { return planes[first].doubled_area > planes[second].doubled_area; });

    // Large triangles provide stable plane normals even when small, skinny CAD
    // triangles amplify decimal-coordinate noise. A bounded seed set keeps this
    // pass linear in mesh size with a small fixed multiplier.
    constexpr std::size_t maximum_normal_domains = 256;
    constexpr double minimum_normal_dot = 0.999;
    std::vector<Vec3> seeds;
    for (const auto index : order)
    {
        const Vec3 normal = planes[index].normal;
        if (std::any_of(seeds.begin(), seeds.end(), [&](const Vec3& seed)
            { return std::abs(seed.dot(normal)) >= minimum_normal_dot; })) continue;
        seeds.push_back(normal);
        if (seeds.size() == maximum_normal_domains) break;
    }
    struct Domain
    {
        Vec3 normal_sum = Vec3::Zero();
        Vec3 normal = Vec3::Zero();
        std::vector<std::size_t> faces;
    };
    std::vector<Domain> domains(seeds.size());
    for (std::size_t index = 0; index < seeds.size(); ++index)
        domains[index].normal = seeds[index];
    for (std::size_t index = 0; index < planes.size(); ++index)
    {
        std::size_t best = seeds.size();
        double best_dot = minimum_normal_dot;
        for (std::size_t domain = 0; domain < seeds.size(); ++domain)
        {
            const double dot = std::abs(seeds[domain].dot(planes[index].normal));
            if (dot < best_dot) continue;
            best_dot = dot;
            best = domain;
        }
        if (best == seeds.size()) continue;
        Vec3 weighted = planes[index].normal * planes[index].doubled_area;
        if (weighted.dot(seeds[best]) < 0.0) weighted *= -1.0;
        domains[best].normal_sum += weighted;
        domains[best].faces.push_back(index);
    }

    const double distance_tolerance = std::max(
        model_diagonal * options.analytic_surface_relative_tolerance, 1.0e-10);
    std::vector<std::vector<std::uint32_t>> result;
    for (Domain& domain : domains)
    {
        if (domain.faces.size() < 6 || domain.normal_sum.norm() <= 1.0e-30) continue;
        domain.normal = domain.normal_sum.normalized();
        std::sort(domain.faces.begin(), domain.faces.end(), [&](const auto first,
                                                                const auto second)
        {
            return domain.normal.dot(planes[first].centroid) <
                   domain.normal.dot(planes[second].centroid);
        });
        std::vector<std::uint32_t> plane_group;
        double previous_distance = 0.0;
        const auto flush = [&]
        {
            if (plane_group.size() >= 6)
            {
                auto connected = faceComponentsFromListApproximate(
                    mesh, plane_group, distance_tolerance);
                for (auto& component : connected)
                    if (component.size() >= 6) result.push_back(std::move(component));
            }
            plane_group.clear();
        };
        for (const auto index : domain.faces)
        {
            const double distance = domain.normal.dot(planes[index].centroid);
            if (!plane_group.empty() && distance - previous_distance > distance_tolerance)
                flush();
            plane_group.push_back(planes[index].id);
            previous_distance = distance;
        }
        flush();
    }
    std::sort(result.begin(), result.end(), [](const auto& first, const auto& second)
    {
        if (first.size() != second.size()) return first.size() > second.size();
        return first.front() < second.front();
    });
    return result;
}

std::vector<OutputPrimitive> classifyApproximatePlanarSurface(
    const Mesh& mesh, const std::vector<std::uint32_t>& faces,
    const Vec3& model_center, const double model_diagonal,
    const PrimitiveMeshAnalysisOptions& options)
{
    if (faces.size() < 6) return {};
    Vec3 reference = Vec3::Zero();
    Vec3 normal_sum = Vec3::Zero();
    for (const auto face_id : faces)
    {
        const Face& face = mesh.faces[face_id];
        Vec3 doubled_area = (mesh.vertices[face[1]] - mesh.vertices[face[0]])
            .cross(mesh.vertices[face[2]] - mesh.vertices[face[0]]);
        if (doubled_area.norm() <= 1.0e-30) continue;
        if (reference.norm() == 0.0) reference = doubled_area.normalized();
        if (doubled_area.dot(reference) < 0.0) doubled_area *= -1.0;
        normal_sum += doubled_area;
    }
    if (normal_sum.norm() <= 1.0e-30) return {};
    Vec3 normal = normal_sum.normalized();
    const auto vertices = uniqueVertices(mesh, faces);
    Vec3 center = Vec3::Zero();
    for (const auto vertex : vertices) center += mesh.vertices[vertex];
    center /= static_cast<double>(vertices.size());
    const double tolerance = std::max(
        model_diagonal * options.analytic_surface_relative_tolerance, 1.0e-10);
    double lower = std::numeric_limits<double>::infinity();
    double upper = -std::numeric_limits<double>::infinity();
    for (const auto vertex : vertices)
    {
        const double distance = (mesh.vertices[vertex] - center).dot(normal);
        lower = std::min(lower, distance);
        upper = std::max(upper, distance);
    }
    if (upper - lower > 2.0 * tolerance) return {};

    const auto loops = boundaryLoopsApproximate(mesh, faces, tolerance);
    if (loops.size() != 1 || loops.front().size() < 3) return {};
    if ((center - model_center).dot(normal) < 0.0) normal *= -1.0;
    const Mat3 basis = orthonormalFrame(normal);
    Mat3 frame;
    frame.col(0) = static_cast<Vec3>(basis.col(1));
    frame.col(1) = static_cast<Vec3>(basis.col(2));
    frame.col(2) = normal;
    std::vector<Vec2> boundary;
    boundary.reserve(loops.front().size());
    for (const auto vertex : loops.front())
    {
        const Vec3 local = frame.transposeMultiply(mesh.vertices[vertex] - center);
        boundary.emplace_back(local.x(), local.y());
    }
    boundary = simplifyPolygon(std::move(boundary), tolerance);
    if (boundary.size() < 3 || !simplePolygon(boundary, tolerance) ||
        triangulatePolygon(boundary).size() + 2 != boundary.size())
        return {};
    for (const auto vertex : vertices)
    {
        const Vec3 local = frame.transposeMultiply(mesh.vertices[vertex] - center);
        if (!pointInPolygon({local.x(), local.y()}, boundary, tolerance))
            return {};
    }
    double support = -std::numeric_limits<double>::infinity();
    for (const auto vertex : vertices)
        support = std::max(
            support, (mesh.vertices[vertex] - center).dot(normal));
    Primitive polygon = polygonPrimitive(
        boundary, center + normal * support, frame);
    return {{std::move(polygon), faces}};
}

// Recognize a complete, connected CAD surface before exact coplanar clustering.
// Independent decimal rounding can move vertices of one analytic plane by a few
// micromodel units and amplify face-normal noise on skinny triangles. Requiring
// every triangle to be exactly coplanar therefore fragments an otherwise clear
// disk into thousands of leaves. This certificate tests the complete component
// against one area-weighted plane and emits an outward support surface.
std::vector<OutputPrimitive> classifyApproximateCircularSurface(
    const Mesh& mesh,
    const std::vector<std::uint32_t>& faces,
    const PrimitiveMeshAnalysisOptions& options,
    const Vec3& model_center,
    const double model_diagonal,
    const double model_volume,
    std::size_t& filled_holes,
    double& filled_volume)
{
    if (!options.allow_round_surfaces || faces.size() < 6) return {};

    Vec3 reference = Vec3::Zero();
    Vec3 normal_sum = Vec3::Zero();
    for (const auto face_id : faces)
    {
        const Face& face = mesh.faces[face_id];
        Vec3 doubled_area = (mesh.vertices[face[1]] - mesh.vertices[face[0]])
            .cross(mesh.vertices[face[2]] - mesh.vertices[face[0]]);
        if (doubled_area.norm() <= 1.0e-30) continue;
        if (reference.norm() == 0.0) reference = doubled_area.normalized();
        if (doubled_area.dot(reference) < 0.0) doubled_area *= -1.0;
        normal_sum += doubled_area;
    }
    if (normal_sum.norm() <= 1.0e-30) return {};
    Vec3 normal = normal_sum.normalized();
    const auto vertices = uniqueVertices(mesh, faces);
    Vec3 plane_center = Vec3::Zero();
    for (const auto vertex : vertices) plane_center += mesh.vertices[vertex];
    plane_center /= static_cast<double>(vertices.size());
    double minimum_plane = std::numeric_limits<double>::infinity();
    double maximum_plane = -std::numeric_limits<double>::infinity();
    for (const auto vertex : vertices)
    {
        const double distance = (mesh.vertices[vertex] - plane_center).dot(normal);
        minimum_plane = std::min(minimum_plane, distance);
        maximum_plane = std::max(maximum_plane, distance);
    }
    const double planar_tolerance = std::max(
        model_diagonal * options.analytic_surface_relative_tolerance, 1.0e-10);
    if (maximum_plane - minimum_plane > 2.0 * planar_tolerance)
    {
        return {};
    }

    auto loops3 = boundaryLoopsApproximate(mesh, faces, planar_tolerance);
    if (loops3.empty())
    {
        return {};
    }
    if ((plane_center - model_center).dot(normal) < 0.0) normal *= -1.0;
    const Mat3 basis = orthonormalFrame(normal);
    Mat3 frame;
    frame.col(0) = static_cast<Vec3>(basis.col(1));
    frame.col(1) = static_cast<Vec3>(basis.col(2));
    frame.col(2) = normal;

    struct Loop
    {
        std::vector<Vec2> points;
        CircularLoopFit circle;
    };
    std::vector<Loop> loops;
    loops.reserve(loops3.size());
    for (const auto& loop3 : loops3)
    {
        Loop loop;
        loop.points.reserve(loop3.size());
        for (const auto vertex : loop3)
        {
            const Vec3 local = frame.transposeMultiply(mesh.vertices[vertex] - plane_center);
            loop.points.emplace_back(local.x(), local.y());
        }
        loop.circle = fitCircularLoop(loop.points, options.circle_radial_tolerance);
        loops.push_back(std::move(loop));
    }
    const auto outer_iterator = std::max_element(
        loops.begin(), loops.end(), [](const Loop& first, const Loop& second)
        { return first.circle.area < second.circle.area; });
    if (outer_iterator == loops.end() || !outer_iterator->circle.circular)
    {
        return {};
    }
    const std::size_t outer_index = static_cast<std::size_t>(
        std::distance(loops.begin(), outer_iterator));
    const CircularLoopFit& outer = outer_iterator->circle;
    // Boundary extraction on non-manifold CAD soups may expose only a circular
    // inner loop while omitting the actual outer silhouette.  A disk/annulus is
    // valid only if its fitted outer circle covers every responsibility vertex.
    // Without this certificate a small bore can incorrectly take ownership of
    // an entire rectangular plate.
    const double radial_margin = std::max(
        planar_tolerance,
        2.0 * options.circle_radial_tolerance * outer.maximum_radius);
    for (const auto vertex : vertices)
    {
        const Vec3 local = frame.transposeMultiply(mesh.vertices[vertex] - plane_center);
        const Vec2 radial(local.x() - outer.center.x(),
                          local.y() - outer.center.y());
        if (radial.norm() > outer.maximum_radius + radial_margin) return {};
    }

    double model_lower = std::numeric_limits<double>::infinity();
    double model_upper = -std::numeric_limits<double>::infinity();
    for (const Vec3& vertex : mesh.vertices)
    {
        const double axial = vertex.dot(normal);
        model_lower = std::min(model_lower, axial);
        model_upper = std::max(model_upper, axial);
    }
    const double conservative_depth = model_upper - model_lower;
    const double volume_budget = options.maximum_cavity_added_volume_ratio * model_volume;
    double candidate_accepted_volume = 0.0;
    std::size_t candidate_filled_holes = 0;
    std::vector<std::size_t> retained_holes;
    for (std::size_t index = 0; index < loops.size(); ++index)
    {
        if (index == outer_index) continue;
        // Analytic circular openings use the complete axial model span.  A
        // nearest-hit ray can stop on an unrelated nearby flange and turn a
        // real annulus into a disk, which then fragments the connected conical
        // assembly.  Use the conservative full-span bound for every input.
        const double added_volume = loops[index].circle.area * conservative_depth;
        // The acceptance threshold is per cavity. "Global" means that the
        // denominator is the whole model, not that unrelated holes consume a
        // shared allowance. A large body cavity therefore remains over budget
        // while any number of individually tiny bores may be closed.
        if (added_volume <= volume_budget)
        {
            candidate_accepted_volume += added_volume;
            ++candidate_filled_holes;
        }
        else
            retained_holes.push_back(index);
    }
    if (retained_holes.size() > 1)
    {
        return {};
    }

    Primitive surface;
    surface.kind = retained_holes.empty() ? Kind::Disk : Kind::Annulus;
    surface.axes.col(0) = normal;
    surface.axes.col(1) = static_cast<Vec3>(frame.col(0));
    surface.axes.col(2) = static_cast<Vec3>(frame.col(1));
    // Place the proxy on the outward support plane and use the largest sampled
    // outer radius. Both choices can only add collision volume.
    double support = -std::numeric_limits<double>::infinity();
    for (const auto vertex : vertices)
        support = std::max(support, (mesh.vertices[vertex] - plane_center).dot(normal));
    surface.center = plane_center + normal * support +
        frame.col(0) * outer.center.x() + frame.col(1) * outer.center.y();
    const double polygon_expansion = options.round_surface_segments >= 3
        ? 1.0 / std::cos(
            std::numbers::pi /
            static_cast<double>(options.round_surface_segments))
        : std::numeric_limits<double>::infinity();
    surface.base_radius = outer.maximum_radius * polygon_expansion *
        (1.0 + 1.0e-12);
    surface.top_radius = surface.base_radius;
    surface.segments = options.round_surface_segments;
    if (!retained_holes.empty())
    {
        const CircularLoopFit& inner = loops[retained_holes.front()].circle;
        const double center_offset = (inner.center - outer.center).norm();
        if (!inner.circular || center_offset >
            options.circle_radial_tolerance * outer.mean_radius) return {};
        // Shrinking the open center is the conservative direction.
        surface.inner_radius = std::max(
            inner.minimum_radius - center_offset, 0.0);
    }
    std::vector<std::uint32_t> owned_faces;
    owned_faces.reserve(faces.size());
    for (const auto face_id : faces)
        if (analyticPlanarSurfaceOwnsFace(
                mesh, face_id, surface, planar_tolerance))
            owned_faces.push_back(face_id);
    if (owned_faces.empty()) return {};
    filled_holes += candidate_filled_holes;
    filled_volume += candidate_accepted_volume;
    return {{std::move(surface), std::move(owned_faces)}};
}

void appendRevolvedSurfacePatches(std::vector<OutputPrimitive>& output,
                                  const Mesh& mesh,
                                  const Primitive& fit,
                                  const std::vector<std::uint32_t>& faces,
                                  const double analytic_tolerance)
{
    const double end_tolerance = std::max({fit.height * 1.0e-6,
                                           analytic_tolerance, 1.0e-12});
    std::vector<std::uint32_t> base_faces;
    std::vector<std::uint32_t> top_faces;
    std::vector<std::uint32_t> band_faces;
    for (const auto face_id : faces)
    {
        bool at_base = true;
        bool at_top = true;
        for (const auto vertex : mesh.faces[face_id])
        {
            const double axial = fit.axes.transposeMultiply(
                mesh.vertices[vertex] - fit.center).x();
            at_base &= std::abs(axial) <= end_tolerance;
            at_top &= std::abs(axial - fit.height) <= end_tolerance;
        }
        if (at_base) base_faces.push_back(face_id);
        else if (at_top) top_faces.push_back(face_id);
        else band_faces.push_back(face_id);
    }
    if (!band_faces.empty())
    {
        Primitive band = fit;
        band.kind = std::abs(fit.base_radius - fit.top_radius) <=
                1.0e-8 * std::max({fit.base_radius, fit.top_radius, fit.height})
            ? Kind::CylindricalBand : Kind::ConicalBand;
        output.push_back({std::move(band), std::move(band_faces)});
    }
    const auto append_disk = [&](const Vec3& center, const double radius,
                                 std::vector<std::uint32_t> cap_faces)
    {
        if (cap_faces.empty() || radius <= 0.0) return;
        Primitive disk;
        disk.kind = Kind::Disk;
        disk.center = center;
        disk.axes = fit.axes;
        disk.base_radius = radius;
        disk.top_radius = radius;
        disk.segments = fit.segments;
        output.push_back({std::move(disk), std::move(cap_faces)});
    };
    append_disk(fit.center, fit.base_radius, std::move(base_faces));
    append_disk(fit.center + fit.axes.col(0) * fit.height,
                fit.top_radius, std::move(top_faces));
}

std::uint32_t roundSegmentsForError(const double radius,
                                    const double maximum_error,
                                    const std::uint32_t minimum_segments)
{
    if (radius <= 0.0 || maximum_error <= 0.0)
        return std::max<std::uint32_t>(minimum_segments, 3);
    // A regular polygon with vertex radius R*sec(pi/n) circumscribes a circle
    // of radius R. Its maximum outward radial error is
    // R*(sec(pi/n)-1). Reserve ten percent of the user's total directed error
    // for later arithmetic/canonicalization and solve this bound for n.
    const double budget = maximum_error * 0.9;
    const double cosine = std::clamp(
        radius / (radius + budget), 0.0, 1.0);
    const double angle = std::acos(cosine);
    if (angle <= 1.0e-12) return 4096;
    const auto required = static_cast<std::uint32_t>(
        std::ceil(std::numbers::pi / angle));
    return std::clamp(
        std::max(required, minimum_segments),
        static_cast<std::uint32_t>(3),
        static_cast<std::uint32_t>(4096));
}

void refineRoundSurfaceSegments(std::vector<OutputPrimitive>& surfaces,
                                const double maximum_error)
{
    for (OutputPrimitive& item : surfaces)
    {
        Primitive& primitive = item.primitive;
        if (!isCertifiedRoundSurfaceKind(primitive.kind) ||
            !primitive.band_axial_ranges.empty())
            continue;
        primitive.segments = roundSegmentsForError(
            std::max(primitive.base_radius, primitive.top_radius),
            maximum_error, primitive.segments);
    }
}

struct CylindricalRegionExtractionStats
{
    std::size_t candidate_axes = 0;
    std::size_t smooth_components = 0;
    std::size_t rejected_small = 0;
    std::size_t rejected_circle = 0;
    std::size_t rejected_angular_coverage = 0;
    std::size_t rejected_workload = 0;
    std::size_t accepted_regions = 0;
    std::size_t accepted_source_faces = 0;
};

std::vector<std::vector<std::uint32_t>> smoothFaceComponentsApproximate(
    const Mesh& mesh,
    const std::vector<bool>& included_faces,
    double tolerance,
    double minimum_normal_dot,
    std::vector<std::vector<std::uint32_t>>* adjacency_output);

// Extract complete cylindrical side patches before planar clustering.  The old
// pipeline tested one whole connected model at a time; a cap, rib, or attached
// bracket then made the end-ring certificate fail and every curved strip fell
// through as an unrelated planar polygon.  This routine grows only smooth
// lateral regions around edge-supported axis hypotheses.  It emits a surface,
// never a closed solid and never a box fallback.
std::vector<OutputPrimitive> extractCylindricalSurfaceRegions(
    const Mesh& mesh,
    const std::vector<bool>& included_faces,
    const PrimitiveMeshAnalysisOptions& options,
    const double geometric_tolerance,
    CylindricalRegionExtractionStats& stats)
{
    struct AxisBucket
    {
        Vec3 sum = Vec3::Zero();
        std::size_t count = 0;
    };
    std::map<std::array<int, 3>, AxisBucket> histogram;
    const auto canonicalAxis = [](Vec3 direction)
    {
        direction = direction.normalized();
        for (int axis = 0; axis < 3; ++axis)
            if (std::abs(direction[axis]) > 1.0e-12)
            {
                if (direction[axis] < 0.0) direction *= -1.0;
                break;
            }
        return direction;
    };
    for (std::uint32_t face_id = 0; face_id < mesh.faces.size(); ++face_id)
    {
        if (!included_faces[face_id]) continue;
        const Face& face = mesh.faces[face_id];
        for (int edge = 0; edge < 3; ++edge)
        {
            Vec3 direction = mesh.vertices[face[(edge + 1) % 3]] -
                             mesh.vertices[face[edge]];
            if (direction.norm() <= geometric_tolerance) continue;
            direction = canonicalAxis(direction);
            constexpr double quantization = 32.0;
            std::array<int, 3> key{};
            for (int axis = 0; axis < 3; ++axis)
                key[axis] = static_cast<int>(
                    std::llround(direction[axis] * quantization));
            AxisBucket& bucket = histogram[key];
            bucket.sum += direction;
            ++bucket.count;
        }
    }
    std::vector<std::pair<std::size_t, Vec3>> ranked_axes;
    ranked_axes.reserve(histogram.size() + 6);
    for (const auto& [key, bucket] : histogram)
    {
        (void)key;
        if (bucket.count < 4 || bucket.sum.norm() <= 1.0e-30) continue;
        ranked_axes.emplace_back(bucket.count,
                                 canonicalAxis(bucket.sum));
    }
    std::sort(ranked_axes.begin(), ranked_axes.end(),
        [](const auto& first, const auto& second)
        { return first.first > second.first; });
    std::vector<Vec3> axes{
        Vec3{1.0, 0.0, 0.0}, Vec3{0.0, 1.0, 0.0}, Vec3{0.0, 0.0, 1.0}};
    const std::size_t maximum_histogram_axes = 24;
    for (std::size_t index = 0;
         index < std::min(ranked_axes.size(), maximum_histogram_axes); ++index)
    {
        const Vec3 axis = ranked_axes[index].second;
        if (std::none_of(axes.begin(), axes.end(), [&](const Vec3& existing)
            { return std::abs(existing.dot(axis)) >= 0.9995; }))
            axes.push_back(axis);
    }
    stats.candidate_axes = axes.size();

    struct Candidate
    {
        Primitive surface;
        std::vector<std::uint32_t> faces;
        std::size_t output_triangles = 0;
    };
    std::vector<Candidate> candidates;
    const double maximum_axial_normal = 0.20;
    const double minimum_smooth_normal_dot =
        std::cos(50.0 * std::numbers::pi / 180.0);
    // Geometric welding and edge adjacency are independent of the candidate
    // axis. Build them once; rebuilding the same CAD topology for every one of
    // roughly two dozen axes dominated runtime on the real 4/5/16 models.
    std::vector<Vec3> face_normals(mesh.faces.size(), Vec3::Zero());
    for (std::uint32_t face_id = 0; face_id < mesh.faces.size(); ++face_id)
    {
        if (!included_faces[face_id]) continue;
        const Face& face = mesh.faces[face_id];
        Vec3 normal = (mesh.vertices[face[1]] - mesh.vertices[face[0]])
            .cross(mesh.vertices[face[2]] - mesh.vertices[face[0]]);
        if (normal.norm() > 1.0e-30)
            face_normals[face_id] = normal.normalized();
    }
    std::vector<std::vector<std::uint32_t>> smooth_adjacency;
    const auto ignored_smooth_components = smoothFaceComponentsApproximate(
        mesh, included_faces, geometric_tolerance,
        minimum_smooth_normal_dot, &smooth_adjacency);
    (void)ignored_smooth_components;
    for (const Vec3& axis : axes)
    {
        std::vector<bool> lateral_faces(mesh.faces.size(), false);
        for (std::uint32_t face_id = 0; face_id < mesh.faces.size(); ++face_id)
        {
            if (!included_faces[face_id]) continue;
            if (face_normals[face_id].norm() <= 1.0e-30)
                continue;
            if (std::abs(face_normals[face_id].dot(axis)) >
                maximum_axial_normal) continue;
            lateral_faces[face_id] = true;
        }
        // Remove the interiors of large planar patches before circle fitting.
        // A tessellated cylinder has a small but persistent normal rotation
        // between neighboring strips; a flat wall does not.  Include the
        // coplanar triangle mate of every curved seed so split quads remain
        // complete, then rebuild components in the same geometric adjacency.
        std::vector<bool> curved_faces(mesh.faces.size(), false);
        constexpr double planar_neighbor_dot = 0.99995;
        for (std::uint32_t face_id = 0; face_id < mesh.faces.size(); ++face_id)
        {
            if (!lateral_faces[face_id]) continue;
            const bool curved_seed = std::any_of(
                smooth_adjacency[face_id].begin(),
                smooth_adjacency[face_id].end(),
                [&](const auto neighbor)
                {
                    return lateral_faces[neighbor] &&
                        std::abs(face_normals[face_id].dot(
                        face_normals[neighbor])) < planar_neighbor_dot;
                });
            if (!curved_seed) continue;
            curved_faces[face_id] = true;
            for (const auto neighbor : smooth_adjacency[face_id])
                if (lateral_faces[neighbor]) curved_faces[neighbor] = true;
        }
        std::vector<std::vector<std::uint32_t>> components;
        std::vector<bool> visited(mesh.faces.size(), false);
        for (std::uint32_t seed = 0; seed < mesh.faces.size(); ++seed)
        {
            if (!curved_faces[seed] || visited[seed]) continue;
            components.emplace_back();
            std::queue<std::uint32_t> queue;
            queue.push(seed);
            visited[seed] = true;
            while (!queue.empty())
            {
                const auto current = queue.front();
                queue.pop();
                components.back().push_back(current);
                for (const auto neighbor : smooth_adjacency[current])
                    if (curved_faces[neighbor] && !visited[neighbor])
                    {
                        visited[neighbor] = true;
                        queue.push(neighbor);
                    }
            }
        }
        stats.smooth_components += components.size();
        for (auto& component : components)
        {
            if (component.size() < 8)
            {
                ++stats.rejected_small;
                continue;
            }
            const auto vertices = uniqueVertices(mesh, component);
            if (vertices.size() < 8)
            {
                ++stats.rejected_small;
                continue;
            }
            const Mat3 frame = orthonormalFrame(axis);
            Vec3 origin = Vec3::Zero();
            for (const auto vertex : vertices) origin += mesh.vertices[vertex];
            origin /= static_cast<double>(vertices.size());
            std::vector<Vec2> projected;
            projected.reserve(vertices.size());
            double lower = std::numeric_limits<double>::infinity();
            double upper = -std::numeric_limits<double>::infinity();
            for (const auto vertex : vertices)
            {
                const Vec3 local = frame.transposeMultiply(
                    mesh.vertices[vertex] - origin);
                lower = std::min(lower, local.x());
                upper = std::max(upper, local.x());
                projected.emplace_back(local.y(), local.z());
            }
            const auto center = leastSquaresCircleCenter(projected);
            if (!center)
            {
                ++stats.rejected_circle;
                continue;
            }
            double minimum_radius = std::numeric_limits<double>::infinity();
            double maximum_radius = 0.0;
            double mean_radius = 0.0;
            std::vector<double> angles;
            angles.reserve(projected.size());
            for (Vec2 point : projected)
            {
                point = point - *center;
                const double radius = point.norm();
                minimum_radius = std::min(minimum_radius, radius);
                maximum_radius = std::max(maximum_radius, radius);
                mean_radius += radius;
                double angle = std::atan2(point.y(), point.x());
                if (angle < 0.0) angle += 2.0 * std::numbers::pi;
                angles.push_back(angle);
            }
            mean_radius /= static_cast<double>(projected.size());
            const double allowed_radial_spread = std::max(
                2.0 * options.circle_radial_tolerance,
                4.0 * geometric_tolerance /
                    std::max(mean_radius, geometric_tolerance));
            if (mean_radius <= geometric_tolerance ||
                (maximum_radius - minimum_radius) / mean_radius >
                    allowed_radial_spread)
            {
                ++stats.rejected_circle;
                continue;
            }
            const Vec3 fitted_axis_origin = origin +
                frame.col(1) * center->x() + frame.col(2) * center->y();
            bool radial_normals = true;
            const double minimum_radial_normal_dot =
                std::cos(20.0 * std::numbers::pi / 180.0);
            for (const auto face_id : component)
            {
                const Face& face = mesh.faces[face_id];
                const Vec3 centroid = (mesh.vertices[face[0]] +
                    mesh.vertices[face[1]] + mesh.vertices[face[2]]) / 3.0;
                Vec3 normal = (mesh.vertices[face[1]] - mesh.vertices[face[0]])
                    .cross(mesh.vertices[face[2]] - mesh.vertices[face[0]]);
                if (normal.norm() <= 1.0e-30) continue;
                normal = normal.normalized();
                Vec3 radial = centroid - fitted_axis_origin;
                radial = radial - axis * radial.dot(axis);
                if (radial.norm() <= geometric_tolerance ||
                    std::abs(normal.dot(radial.normalized())) <
                        minimum_radial_normal_dot)
                {
                    radial_normals = false;
                    break;
                }
            }
            if (!radial_normals)
            {
                ++stats.rejected_circle;
                continue;
            }
            std::sort(angles.begin(), angles.end());
            angles.erase(std::unique(angles.begin(), angles.end(),
                [](const double first, const double second)
                { return std::abs(first - second) <= 1.0e-5; }), angles.end());
            double maximum_gap = angles.empty()
                ? 2.0 * std::numbers::pi
                : angles.front() + 2.0 * std::numbers::pi - angles.back();
            for (std::size_t index = 1; index < angles.size(); ++index)
                maximum_gap = std::max(maximum_gap,
                                       angles[index] - angles[index - 1]);
            const double covered_angle = 2.0 * std::numbers::pi - maximum_gap;
            if (angles.size() < 3 || covered_angle <
                    20.0 * std::numbers::pi / 180.0)
            {
                ++stats.rejected_angular_coverage;
                continue;
            }
            const double height = upper - lower;
            if (height <= geometric_tolerance)
            {
                ++stats.rejected_small;
                continue;
            }
            Primitive surface;
            surface.kind = Kind::CylindricalBand;
            surface.axes = frame;
            surface.center = origin + frame.col(0) * lower +
                frame.col(1) * center->x() + frame.col(2) * center->y();
            surface.height = height;
            const double polygon_expansion = options.round_surface_segments >= 3
                ? 1.0 / std::cos(std::numbers::pi /
                    static_cast<double>(options.round_surface_segments))
                : std::numeric_limits<double>::infinity();
            surface.base_radius = surface.top_radius = maximum_radius *
                polygon_expansion * (1.0 + 1.0e-12);
            surface.segments = options.round_surface_segments;
            surface.band_axial_ranges.assign(
                surface.segments, std::array<double, 2>{0.0, height});
            surface.band_active_segments.assign(surface.segments, 0);
            const double angular_step = 2.0 * std::numbers::pi /
                static_cast<double>(surface.segments);
            bool spans_too_wide = false;
            for (const auto face_id : component)
            {
                std::array<double, 3> face_angles{};
                for (int corner = 0; corner < 3; ++corner)
                {
                    const Vec3 local = frame.transposeMultiply(
                        mesh.vertices[mesh.faces[face_id][corner]] -
                        fitted_axis_origin);
                    face_angles[corner] = std::atan2(local.z(), local.y());
                }
                for (int corner = 1; corner < 3; ++corner)
                {
                    while (face_angles[corner] - face_angles[0] >
                           std::numbers::pi)
                        face_angles[corner] -= 2.0 * std::numbers::pi;
                    while (face_angles[corner] - face_angles[0] <
                          -std::numbers::pi)
                        face_angles[corner] += 2.0 * std::numbers::pi;
                }
                const auto [minimum_angle, maximum_angle] =
                    std::minmax_element(face_angles.begin(), face_angles.end());
                if (*maximum_angle - *minimum_angle > angular_step * 1.05)
                {
                    spans_too_wide = true;
                    break;
                }
                const int first_segment = static_cast<int>(std::floor(
                    (*minimum_angle - 1.0e-9) / angular_step));
                const int last_segment = static_cast<int>(std::floor(
                    (*maximum_angle + 1.0e-9) / angular_step));
                for (int segment = first_segment; segment <= last_segment;
                     ++segment)
                {
                    const int cyclic = (segment %
                        static_cast<int>(surface.segments) +
                        static_cast<int>(surface.segments)) %
                        static_cast<int>(surface.segments);
                    surface.band_active_segments[cyclic] = 1;
                }
            }
            if (spans_too_wide)
            {
                ++stats.rejected_angular_coverage;
                continue;
            }
            const std::size_t output_triangles =
                triangulatePrimitive(surface).faces.size();
            if (output_triangles > component.size())
            {
                ++stats.rejected_workload;
                continue;
            }
            candidates.push_back(
                {std::move(surface), std::move(component), output_triangles});
        }
    }
    std::sort(candidates.begin(), candidates.end(),
        [](const Candidate& first, const Candidate& second)
        {
            const std::size_t first_gain =
                first.faces.size() - first.output_triangles;
            const std::size_t second_gain =
                second.faces.size() - second.output_triangles;
            if (first_gain != second_gain) return first_gain > second_gain;
            return first.faces.size() > second.faces.size();
        });
    std::vector<bool> claimed(mesh.faces.size(), false);
    std::vector<OutputPrimitive> result;
    for (Candidate& candidate : candidates)
    {
        if (std::any_of(candidate.faces.begin(), candidate.faces.end(),
            [&](const auto face) { return claimed[face]; }))
            continue;
        for (const auto face : candidate.faces) claimed[face] = true;
        stats.accepted_source_faces += candidate.faces.size();
        ++stats.accepted_regions;
        result.push_back(
            {std::move(candidate.surface), std::move(candidate.faces)});
    }
    return result;
}

PrimitiveMesh triangulatePrimitive(const Primitive& primitive)
{
    PrimitiveMesh result;
    if (primitive.kind == Kind::Polygon)
    {
        if (primitive.polygon.size() < 3)
            throw std::runtime_error("polygon primitive has fewer than three vertices");
        Vec3 polygon_lower = Vec3::Constant(
            std::numeric_limits<double>::infinity());
        Vec3 polygon_upper = Vec3::Constant(
            -std::numeric_limits<double>::infinity());
        for (const Vec3& vertex : primitive.polygon)
        {
            polygon_lower = polygon_lower.cwiseMin(vertex);
            polygon_upper = polygon_upper.cwiseMax(vertex);
        }
        const double duplicate_tolerance = std::max(
            (polygon_upper - polygon_lower).norm() * 1.0e-12,
            1.0e-12);
        result.vertices.reserve(primitive.polygon.size());
        for (const Vec3& vertex : primitive.polygon)
            if (std::none_of(
                    result.vertices.begin(), result.vertices.end(),
                    [&](const Vec3& existing)
                    { return (existing - vertex).norm() <= duplicate_tolerance; }))
                result.vertices.push_back(vertex);
        if (result.vertices.size() < 3)
        {
            result.vertices.clear();
            return result;
        }

        Vec3 normal = Vec3::Zero();
        for (std::size_t second = 1;
             second + 1 < result.vertices.size() &&
             normal.norm() <= duplicate_tolerance * duplicate_tolerance;
             ++second)
            for (std::size_t third = second + 1;
                 third < result.vertices.size(); ++third)
            {
                normal = (result.vertices[second] - result.vertices[0]).cross(
                    result.vertices[third] - result.vertices[0]);
                if (normal.norm() >
                    duplicate_tolerance * duplicate_tolerance)
                    break;
            }
        if (normal.norm() <= duplicate_tolerance * duplicate_tolerance)
        {
            result.vertices.clear();
            return result;
        }
        Mat3 basis = orthonormalFrame(normal.normalized());
        Mat3 frame;
        frame.col(0) = static_cast<Vec3>(basis.col(1));
        frame.col(1) = static_cast<Vec3>(basis.col(2));
        frame.col(2) = normal.normalized();
        std::vector<Vec2> boundary;
        boundary.reserve(result.vertices.size());
        for (const Vec3& vertex : result.vertices)
        {
            const Vec3 local = frame.transposeMultiply(vertex - result.vertices.front());
            boundary.emplace_back(local.x(), local.y());
        }
        const auto indices = triangulatePolygon(std::move(boundary));
        if (indices.size() + 2 != result.vertices.size())
        {
            std::ostringstream message;
            message << "simple polygon triangulation did not produce n-2 triangles; vertices=";
            for (const Vec3& vertex : primitive.polygon)
                message << " (" << vertex.x() << ',' << vertex.y() << ',' << vertex.z() << ')';
            throw std::runtime_error(message.str());
        }
        for (const auto& triangle : indices)
            result.faces.push_back({static_cast<std::uint32_t>(triangle[0]),
                                    static_cast<std::uint32_t>(triangle[1]),
                                    static_cast<std::uint32_t>(triangle[2])});
        return result;
    }
    if (primitive.kind == Kind::Triangle)
    {
        result.vertices.assign(primitive.triangle.begin(), primitive.triangle.end());
        result.faces.push_back({0, 1, 2});
        return result;
    }
    if (primitive.kind == Kind::Rectangle)
    {
        for (const double x : {-1.0, 1.0})
            for (const double y : {-1.0, 1.0})
                result.vertices.push_back(primitive.center + primitive.axes.col(0) *
                    (x * primitive.half_size.x()) + primitive.axes.col(1) *
                    (y * primitive.half_size.y()));
        result.faces = {{0, 2, 3}, {0, 3, 1}};
        return result;
    }

    const std::uint32_t segments = std::max<std::uint32_t>(primitive.segments, 3);
    if (primitive.kind == Kind::Disk)
    {
        result.vertices.push_back(primitive.center);
        for (std::uint32_t index = 0; index < segments; ++index)
        {
            const double angle = 2.0 * std::numbers::pi * index / segments;
            result.vertices.push_back(primitive.center + primitive.axes.col(1) *
                (primitive.base_radius * std::cos(angle)) + primitive.axes.col(2) *
                (primitive.base_radius * std::sin(angle)));
            result.faces.push_back({0, index + 1, (index + 1) % segments + 1});
        }
        return result;
    }
    if (primitive.kind == Kind::Annulus)
    {
        for (std::uint32_t index = 0; index < segments; ++index)
        {
            const double angle = 2.0 * std::numbers::pi * index / segments;
            const Vec3 radial = primitive.axes.col(1) * std::cos(angle) +
                                primitive.axes.col(2) * std::sin(angle);
            result.vertices.push_back(primitive.center + radial * primitive.inner_radius);
            result.vertices.push_back(primitive.center + radial * primitive.base_radius);
        }
        for (std::uint32_t index = 0; index < segments; ++index)
        {
            const std::uint32_t next = (index + 1) % segments;
            result.faces.push_back({2 * index, 2 * index + 1, 2 * next + 1});
            result.faces.push_back({2 * index, 2 * next + 1, 2 * next});
        }
        return result;
    }
    const bool band = primitive.kind == Kind::CylindricalBand ||
                      primitive.kind == Kind::ConicalBand;
    if (band && !primitive.band_axial_ranges.empty())
    {
        const std::uint32_t range_count = static_cast<std::uint32_t>(
            primitive.band_axial_ranges.size());
        if (range_count < 3)
            throw std::runtime_error(
                "trimmed analytic band has fewer than three angular ranges");
        result.vertices.reserve(2 * range_count);
        result.faces.reserve(2 * range_count);
        for (std::uint32_t index = 0; index < range_count; ++index)
        {
            const double angle = 2.0 * std::numbers::pi * index /
                                 static_cast<double>(range_count);
            const Vec3 radial_direction =
                primitive.axes.col(1) * std::cos(angle) +
                primitive.axes.col(2) * std::sin(angle);
            for (const double axial : primitive.band_axial_ranges[index])
            {
                const double fraction = primitive.height > 0.0
                    ? std::clamp(axial / primitive.height, 0.0, 1.0)
                    : 0.0;
                const double radius =
                    primitive.base_radius * (1.0 - fraction) +
                    primitive.top_radius * fraction;
                result.vertices.push_back(
                    primitive.center + primitive.axes.col(0) * axial +
                    radial_direction * radius);
            }
        }
        for (std::uint32_t index = 0; index < range_count; ++index)
        {
            const std::uint32_t next = (index + 1) % range_count;
            if (!primitive.band_active_segments.empty() &&
                (index >= primitive.band_active_segments.size() ||
                 primitive.band_active_segments[index] == 0))
                continue;
            result.faces.push_back(
                {2 * index, 2 * next, 2 * next + 1});
            result.faces.push_back(
                {2 * index, 2 * next + 1, 2 * index + 1});
        }
        return result;
    }
    if (!band)
        throw std::logic_error("unsupported non-surface primitive kind");
    const Vec3 top_center = primitive.center +
        primitive.axes.col(0) * primitive.height;
    for (std::uint32_t index = 0; index < segments; ++index)
    {
        const double angle = 2.0 * std::numbers::pi * index / segments;
        const Vec3 radial = primitive.axes.col(1) * std::cos(angle) +
                            primitive.axes.col(2) * std::sin(angle);
        result.vertices.push_back(primitive.center + radial * primitive.base_radius);
        result.vertices.push_back(top_center + radial * primitive.top_radius);
    }
    for (std::uint32_t index = 0; index < segments; ++index)
    {
        const std::uint32_t next = (index + 1) % segments;
        const std::uint32_t base = 2 * index;
        const std::uint32_t top = base + 1;
        const std::uint32_t next_base = 2 * next;
        const std::uint32_t next_top = next_base + 1;
        result.faces.push_back({base, next_base, next_top});
        result.faces.push_back({base, next_top, top});
    }
    return result;
}

struct Bounds
{
    Vec3 lower = Vec3::Constant(std::numeric_limits<double>::infinity());
    Vec3 upper = Vec3::Constant(-std::numeric_limits<double>::infinity());
};

Bounds primitiveBounds(const Primitive& primitive)
{
    Bounds bounds;
    for (const Vec3& vertex : triangulatePrimitive(primitive).vertices)
    {
        bounds.lower = bounds.lower.cwiseMin(vertex);
        bounds.upper = bounds.upper.cwiseMax(vertex);
    }
    return bounds;
}

bool boundsContain(const Bounds& outer, const Bounds& inner, const double tolerance)
{
    for (int axis = 0; axis < 3; ++axis)
        if (inner.lower[axis] < outer.lower[axis] - tolerance ||
            inner.upper[axis] > outer.upper[axis] + tolerance) return false;
    return true;
}

bool containsPoint(const Primitive& primitive, const Vec3& point, const double tolerance)
{
    if (primitive.kind == Kind::Disk || primitive.kind == Kind::Annulus)
    {
        const Vec3 local = primitive.axes.transposeMultiply(point - primitive.center);
        if (std::abs(local.x()) > tolerance) return false;
        const double radius = std::hypot(local.y(), local.z());
        return radius <= primitive.base_radius + tolerance &&
               (primitive.kind != Kind::Annulus ||
                radius >= primitive.inner_radius - tolerance);
    }
    if (primitive.kind == Kind::Polygon)
    {
        const PrimitiveMesh mesh = triangulatePrimitive(primitive);
        for (const Face& face : mesh.faces)
        {
            Primitive triangle;
            triangle.kind = Kind::Triangle;
            for (int corner = 0; corner < 3; ++corner)
                triangle.triangle[corner] = mesh.vertices[face[corner]];
            if (containsPoint(triangle, point, tolerance)) return true;
        }
        return false;
    }
    if (primitive.kind == Kind::Rectangle)
    {
        const Vec3 local = primitive.axes.transposeMultiply(point - primitive.center);
        return std::abs(local.x()) <= primitive.half_size.x() + tolerance &&
               std::abs(local.y()) <= primitive.half_size.y() + tolerance &&
               std::abs(local.z()) <= primitive.half_size.z() + tolerance;
    }
    const Vec3 first = primitive.triangle[1] - primitive.triangle[0];
    const Vec3 second = primitive.triangle[2] - primitive.triangle[0];
    const Vec3 normal = first.cross(second);
    const double normal_length = normal.norm();
    if (normal_length <= tolerance) return false;
    if (std::abs((point - primitive.triangle[0]).dot(normal)) > tolerance * normal_length)
        return false;
    const Vec3 relative = point - primitive.triangle[0];
    const double first_first = first.dot(first);
    const double first_second = first.dot(second);
    const double second_second = second.dot(second);
    const double relative_first = relative.dot(first);
    const double relative_second = relative.dot(second);
    const double denominator = first_first * second_second - first_second * first_second;
    if (std::abs(denominator) <= 1.0e-30) return false;
    const double u = (second_second * relative_first - first_second * relative_second) / denominator;
    const double v = (first_first * relative_second - first_second * relative_first) / denominator;
    return u >= -tolerance && v >= -tolerance && u + v <= 1.0 + tolerance;
}

bool planarNormal(const Primitive& primitive, double tolerance, Vec3& normal);
Vec3 planarPoint(const Primitive& primitive);

bool trimmedAnalyticBandContainsTriangle(
    const Primitive& primitive, const std::array<Vec3, 3>& triangle,
    const double tolerance)
{
    if ((primitive.kind != Kind::CylindricalBand &&
         primitive.kind != Kind::ConicalBand) ||
        primitive.height <= 0.0)
        return false;
    const bool trimmed = !primitive.band_axial_ranges.empty();
    const std::size_t segments = trimmed
        ? primitive.band_axial_ranges.size()
        : static_cast<std::size_t>(primitive.segments);
    if (segments < 3) return false;
    if (!primitive.band_active_segments.empty() &&
        primitive.band_active_segments.size() != segments)
        return false;
    const double step = 2.0 * std::numbers::pi /
                        static_cast<double>(segments);
    std::array<double, 3> angles{};
    std::array<Vec3, 3> local{};
    for (std::size_t corner = 0; corner < triangle.size(); ++corner)
    {
        local[corner] = primitive.axes.transposeMultiply(
            triangle[corner] - primitive.center);
        angles[corner] = std::atan2(local[corner].z(), local[corner].y());
        if (angles[corner] < 0.0)
            angles[corner] += 2.0 * std::numbers::pi;
    }
    // A certified source triangle must lie in one angular polygon sector.
    // Unwrap around its first vertex so the 0/2*pi seam is not special.
    for (std::size_t corner = 1; corner < angles.size(); ++corner)
    {
        while (angles[corner] - angles[0] > std::numbers::pi)
            angles[corner] -= 2.0 * std::numbers::pi;
        while (angles[corner] - angles[0] < -std::numbers::pi)
            angles[corner] += 2.0 * std::numbers::pi;
    }
    const auto [minimum_angle, maximum_angle] = std::minmax_element(
        angles.begin(), angles.end());
    if (*maximum_angle - *minimum_angle > step * 1.05)
        return false;

    for (std::size_t corner = 0; corner < triangle.size(); ++corner)
    {
        double cyclic_angle = std::fmod(
            angles[corner], 2.0 * std::numbers::pi);
        if (cyclic_angle < 0.0) cyclic_angle += 2.0 * std::numbers::pi;
        const double sample = cyclic_angle / step;
        const std::size_t first = static_cast<std::size_t>(
            std::floor(sample)) % segments;
        const std::size_t second = (first + 1) % segments;
        if (!primitive.band_active_segments.empty() &&
            primitive.band_active_segments[first] == 0)
            return false;
        const double fraction = sample - std::floor(sample);
        const double lower = trimmed
            ? primitive.band_axial_ranges[first][0] * (1.0 - fraction) +
                primitive.band_axial_ranges[second][0] * fraction
            : 0.0;
        const double upper = trimmed
            ? primitive.band_axial_ranges[first][1] * (1.0 - fraction) +
                primitive.band_axial_ranges[second][1] * fraction
            : primitive.height;
        if (local[corner].x() < lower - tolerance ||
            local[corner].x() > upper + tolerance)
            return false;
        const double axial_fraction = std::clamp(
            local[corner].x() / primitive.height, 0.0, 1.0);
        const double vertex_radius =
            primitive.base_radius * (1.0 - axial_fraction) +
            primitive.top_radius * axial_fraction;
        const double sector_midpoint =
            (static_cast<double>(first) + 0.5) * step;
        double midpoint_offset = cyclic_angle - sector_midpoint;
        if (midpoint_offset > std::numbers::pi)
            midpoint_offset -= 2.0 * std::numbers::pi;
        if (midpoint_offset < -std::numbers::pi)
            midpoint_offset += 2.0 * std::numbers::pi;
        const double polygon_radius = vertex_radius *
            std::cos(0.5 * step) /
            std::max(std::cos(midpoint_offset), 1.0e-12);
        if (std::hypot(local[corner].y(), local[corner].z()) >
            polygon_radius + tolerance)
            return false;
    }
    return true;
}

bool outwardPlanarSupportContainsTriangle(
    const Primitive& primitive, const std::array<Vec3, 3>& triangle,
    const double tolerance)
{
    Vec3 normal;
    if (!planarNormal(primitive, tolerance, normal)) return false;
    const Vec3 point = planarPoint(primitive);
    const PrimitiveMesh proxy = triangulatePrimitive(primitive);
    if (proxy.faces.empty()) return false;
    for (const Vec3& vertex : triangle)
    {
        const double signed_distance = (vertex - point).dot(normal);
        // Analytic CAD planes are fitted on their outward support plane.  The
        // source may lie microscopically behind that plane after decimal
        // quantization, but it may never protrude through the proxy.
        if (signed_distance > tolerance) return false;
        const Vec3 projected = vertex - normal * signed_distance;
        bool covered = false;
        for (const Face& face : proxy.faces)
        {
            Primitive proxy_triangle;
            proxy_triangle.kind = Kind::Triangle;
            for (int corner = 0; corner < 3; ++corner)
                proxy_triangle.triangle[corner] = proxy.vertices[face[corner]];
            if (containsPoint(proxy_triangle, projected, tolerance))
            {
                covered = true;
                break;
            }
        }
        if (!covered) return false;
    }
    return true;
}

bool containsPointCached(const Primitive& primitive,
                         const PrimitiveMesh& triangulated,
                         const Vec3& point,
                         const double tolerance)
{
    if (primitive.kind != Kind::Polygon)
        return containsPoint(primitive, point, tolerance);
    for (const Face& face : triangulated.faces)
    {
        Primitive triangle;
        triangle.kind = Kind::Triangle;
        for (int corner = 0; corner < 3; ++corner)
            triangle.triangle[corner] = triangulated.vertices[face[corner]];
        if (containsPoint(triangle, point, tolerance)) return true;
    }
    if (primitive.kind == Kind::Disk || primitive.kind == Kind::Annulus ||
        primitive.kind == Kind::CylindricalBand ||
        primitive.kind == Kind::ConicalBand)
    {
        const PrimitiveMesh mesh = triangulatePrimitive(primitive);
        for (const Face& face : mesh.faces)
        {
            Primitive triangle;
            triangle.kind = Kind::Triangle;
            for (int corner = 0; corner < 3; ++corner)
                triangle.triangle[corner] = mesh.vertices[face[corner]];
            if (containsPoint(triangle, point, tolerance)) return true;
        }
        return false;
    }
    return false;
}

bool planarNormal(const Primitive& primitive, const double tolerance, Vec3& normal)
{
    if (primitive.kind == Kind::Disk || primitive.kind == Kind::Annulus)
    {
        normal = primitive.axes.col(0);
        return true;
    }
    if (primitive.kind == Kind::Polygon)
    {
        if (primitive.polygon.size() < 3) return false;
        normal = (primitive.polygon[1] - primitive.polygon[0])
                     .cross(primitive.polygon[2] - primitive.polygon[0]);
        if (normal.norm() <= tolerance) return false;
        normal = normal.normalized();
        return true;
    }
    if (primitive.kind == Kind::Rectangle)
    {
        normal = primitive.axes.col(2);
        return true;
    }
    if (primitive.kind == Kind::Triangle)
    {
        normal = (primitive.triangle[1] - primitive.triangle[0])
                     .cross(primitive.triangle[2] - primitive.triangle[0]);
        if (normal.norm() <= tolerance) return false;
        normal = normal.normalized();
        return true;
    }
    return false;
}

Vec3 planarPoint(const Primitive& primitive)
{
    if (primitive.kind == Kind::Polygon) return primitive.polygon.front();
    return primitive.kind == Kind::Triangle ? primitive.triangle[0] : primitive.center;
}

std::vector<OutputPrimitive> removePlanarSurfacesOccludedByDisks(
    std::vector<OutputPrimitive> primitives, const Vec3& model_center,
    const double tolerance, std::size_t& removed_count)
{
    removed_count = 0;
    std::vector<std::uint8_t> removed(primitives.size(), 0);
    std::vector<PrimitiveMesh> meshes(primitives.size());
    std::vector<std::size_t> disks;
    for (std::size_t index = 0; index < primitives.size(); ++index)
    {
        meshes[index] = triangulatePrimitive(primitives[index].primitive);
        if (primitives[index].primitive.kind == Kind::Disk)
            disks.push_back(index);
    }
    const auto pointInMesh = [&](const PrimitiveMesh& mesh, const Vec3& point)
    {
        for (const Face& face : mesh.faces)
        {
            Primitive triangle;
            triangle.kind = Kind::Triangle;
            for (int corner = 0; corner < 3; ++corner)
                triangle.triangle[corner] = mesh.vertices[face[corner]];
            if (containsPoint(triangle, point, tolerance)) return true;
        }
        return false;
    };

    for (std::size_t candidate = 0; candidate < primitives.size(); ++candidate)
    {
        if (primitives[candidate].preserves_cavity_opening) continue;
        Vec3 candidate_normal;
        if (!planarNormal(
                primitives[candidate].primitive, tolerance, candidate_normal))
            continue;
        const Vec3 candidate_point = planarPoint(
            primitives[candidate].primitive);
        if ((candidate_point - model_center).dot(candidate_normal) < 0.0)
            candidate_normal *= -1.0;
        const double candidate_outward =
            (candidate_point - model_center).dot(candidate_normal);

        for (const std::size_t coverer : disks)
        {
            if (coverer == candidate || removed[coverer]) continue;
            Vec3 coverer_normal = primitives[coverer].primitive.axes.col(0);
            const Vec3 coverer_point = primitives[coverer].primitive.center;
            if ((coverer_point - model_center).dot(coverer_normal) < 0.0)
                coverer_normal *= -1.0;
            if (coverer_normal.dot(candidate_normal) < 1.0 - 1.0e-8)
                continue;
            const double coverer_outward =
                (coverer_point - model_center).dot(coverer_normal);
            if (coverer_outward <= candidate_outward + tolerance) continue;

            bool covered = true;
            for (const Vec3& vertex : meshes[candidate].vertices)
            {
                const double distance =
                    (vertex - coverer_point).dot(coverer_normal);
                const Vec3 projected = vertex - coverer_normal * distance;
                if (!pointInMesh(meshes[coverer], projected))
                {
                    covered = false;
                    break;
                }
            }
            if (!covered) continue;

            auto& responsibility = primitives[coverer].source_faces;
            responsibility.insert(
                responsibility.end(),
                primitives[candidate].source_faces.begin(),
                primitives[candidate].source_faces.end());
            removed[candidate] = 1;
            ++removed_count;
            break;
        }
    }

    std::vector<OutputPrimitive> result;
    result.reserve(primitives.size() - removed_count);
    for (std::size_t index = 0; index < primitives.size(); ++index)
        if (!removed[index])
        {
            auto& responsibility = primitives[index].source_faces;
            std::sort(responsibility.begin(), responsibility.end());
            responsibility.erase(
                std::unique(responsibility.begin(), responsibility.end()),
                responsibility.end());
            result.push_back(std::move(primitives[index]));
        }
    return result;
}

std::vector<Vec3> planarCorners(const Primitive& primitive)
{
    if (primitive.kind == Kind::Polygon) return primitive.polygon;
    if (primitive.kind == Kind::Triangle)
        return {primitive.triangle.begin(), primitive.triangle.end()};
    if (primitive.kind == Kind::Rectangle)
    {
        std::vector<Vec3> corners;
        for (const double first : {-1.0, 1.0})
            for (const double second : {-1.0, 1.0})
                corners.push_back(primitive.center + primitive.axes.col(0) *
                    (first * primitive.half_size.x()) + primitive.axes.col(1) *
                    (second * primitive.half_size.y()));
        return corners;
    }
    return {};
}

double planarPrimitiveArea(const Primitive& primitive)
{
    if (primitive.kind == Kind::Polygon)
        return primitiveSurfaceArea(primitive, 0.0);
    if (primitive.kind == Kind::Triangle)
        return 0.5 * (primitive.triangle[1] - primitive.triangle[0])
                         .cross(primitive.triangle[2] - primitive.triangle[0]).norm();
    if (primitive.kind == Kind::Rectangle)
        return 4.0 * primitive.half_size.x() * primitive.half_size.y();
    return 0.0;
}

double maximumFilledSurfaceDistance(
    const Mesh& filled_surface_mesh,
    const std::vector<OutputPrimitive>& output,
    double sample_spacing,
    double maximum_distance);

struct FilledSurfaceDistanceCertificate
{
    bool exceeded = false;
    bool certified = false;
    double observed_maximum = 0.0;
    std::size_t samples = 0;
};

FilledSurfaceDistanceCertificate certifyFilledSurfaceDistance(
    const Mesh& filled_surface_mesh,
    const std::vector<OutputPrimitive>& output,
    double maximum_distance);

Vec3 closestPointOnTriangle(const Vec3& point,
                            const Vec3& first,
                            const Vec3& second,
                            const Vec3& third);

std::vector<OutputPrimitive> mergeAdjacentSurfacePrimitives(
    const Mesh& responsibility_mesh,
    const Mesh& distance_reference,
    std::vector<OutputPrimitive> primitives,
    const double tolerance,
    const double maximum_open_error_distance,
    const double error_sample_spacing,
    std::size_t& merged_count,
    const std::filesystem::path& profile_path)
{
    const auto started = std::chrono::steady_clock::now();
    struct Item
    {
        OutputPrimitive output;
        std::vector<std::uint32_t> source_vertices;
        std::unordered_set<std::size_t> neighbors;
        std::uint64_t version = 0;
        bool active = true;
    };
    struct Fit
    {
        Primitive surface;
        std::size_t triangles = 0;
        double error = 0.0;
        std::vector<std::size_t> absorbed_neighbors;
    };
    struct Candidate
    {
        std::size_t first = 0;
        std::size_t second = 0;
        std::uint64_t first_version = 0;
        std::uint64_t second_version = 0;
    };
    struct Profile
    {
        std::size_t input_primitives = 0;
        std::size_t initial_adjacencies = 0;
        std::size_t responsibility_adjacencies = 0;
        std::size_t boundary_adjacencies = 0;
        std::size_t segment_contact_adjacencies = 0;
        std::size_t candidate_evaluations = 0;
        std::size_t cache_hits = 0;
        std::size_t fit_attempts = 0;
        std::size_t containment_rejections = 0;
        std::size_t connectivity_rejections = 0;
        std::size_t error_rejections = 0;
        std::size_t workload_rejections = 0;
        std::size_t protected_cutout_rejections = 0;
        std::size_t unsupported_rejections = 0;
        std::size_t accepted_merges = 0;
        std::size_t absorbed_neighbors = 0;
        std::size_t stale_neighbors_removed = 0;
        std::size_t accepted_triangle_increases = 0;
        std::size_t remaining_acceptable_candidates = 0;
        std::size_t certified_distance_accepts = 0;
        std::size_t certified_distance_rejections = 0;
        std::size_t certificate_samples = 0;
        std::size_t fine_distance_evaluations = 0;
        double certificate_seconds = 0.0;
        double fine_distance_seconds = 0.0;
        double distance_seconds = 0.0;
    } profile;
    profile.input_primitives = primitives.size();
    std::vector<Item> items;
    items.reserve(primitives.size());
    for (auto& primitive : primitives)
    {
        Item item;
        item.output = std::move(primitive);
        for (const auto face_id : item.output.source_faces)
        {
            if (face_id >= responsibility_mesh.faces.size()) continue;
            const Face& face = responsibility_mesh.faces[face_id];
            item.source_vertices.insert(
                item.source_vertices.end(), face.begin(), face.end());
        }
        std::sort(item.source_vertices.begin(), item.source_vertices.end());
        item.source_vertices.erase(
            std::unique(item.source_vertices.begin(), item.source_vertices.end()),
            item.source_vertices.end());
        items.push_back(std::move(item));
    }

    const auto connect = [&](const std::size_t first, const std::size_t second)
    {
        if (first == second) return false;
        const bool inserted = items[first].neighbors.insert(second).second;
        items[second].neighbors.insert(first);
        return inserted;
    };
    std::vector<std::vector<std::size_t>> face_owners(
        responsibility_mesh.faces.size());
    for (std::size_t primitive = 0; primitive < items.size(); ++primitive)
        for (const auto face : items[primitive].output.source_faces)
            if (face < face_owners.size()) face_owners[face].push_back(primitive);
    std::unordered_map<std::uint64_t, std::vector<std::uint32_t>> edge_faces;
    edge_faces.reserve(responsibility_mesh.faces.size() * 3);
    for (std::uint32_t face = 0; face < responsibility_mesh.faces.size(); ++face)
    {
        const Face& triangle = responsibility_mesh.faces[face];
        for (int edge = 0; edge < 3; ++edge)
            edge_faces[edgeKey(triangle[edge], triangle[(edge + 1) % 3])]
                .push_back(face);
        const auto& owners = face_owners[face];
        for (std::size_t first = 0; first < owners.size(); ++first)
            for (std::size_t second = first + 1; second < owners.size(); ++second)
                connect(owners[first], owners[second]);
    }
    for (const auto& [edge, faces] : edge_faces)
    {
        (void)edge;
        for (std::size_t first_face = 0; first_face < faces.size(); ++first_face)
            for (std::size_t second_face = first_face + 1;
                 second_face < faces.size(); ++second_face)
                for (const auto first : face_owners[faces[first_face]])
                    for (const auto second : face_owners[faces[second_face]])
                        connect(first, second);
    }
    for (const auto& item : items)
        profile.responsibility_adjacencies += item.neighbors.size();
    profile.responsibility_adjacencies /= 2;

    // Adjacency belongs to the current stage-2 surfaces, not merely to their
    // source-triangle responsibility records. Hole caps, restored cavity
    // surfaces, and canonicalized union pieces can legitimately have no source
    // faces. If adjacency is inferred only from source-face ownership, those
    // primitives can never participate in the fixed-point merge. Connect
    // primitives whose actual boundary vertices touch instead. Quantized
    // neighboring cells make this tolerant to the tiny coordinate differences
    // introduced by independent surface fitting, while the source-edge graph
    // above remains a conservative fallback for analytically equivalent edges.
    struct QuantizedPoint
    {
        std::int64_t x = 0;
        std::int64_t y = 0;
        std::int64_t z = 0;
        bool operator==(const QuantizedPoint&) const = default;
    };
    struct QuantizedPointHash
    {
        std::size_t operator()(const QuantizedPoint& point) const noexcept
        {
            std::size_t seed = std::hash<std::int64_t>{}(point.x);
            seed ^= std::hash<std::int64_t>{}(point.y) +
                    0x9e3779b97f4a7c15ULL + (seed << 6U) + (seed >> 2U);
            seed ^= std::hash<std::int64_t>{}(point.z) +
                    0x9e3779b97f4a7c15ULL + (seed << 6U) + (seed >> 2U);
            return seed;
        }
    };
    const double boundary_cell = std::max(tolerance * 8.0, 1.0e-12);
    const auto quantizePoint = [&](const Vec3& point)
    {
        return QuantizedPoint{
            static_cast<std::int64_t>(std::llround(point.x() / boundary_cell)),
            static_cast<std::int64_t>(std::llround(point.y() / boundary_cell)),
            static_cast<std::int64_t>(std::llround(point.z() / boundary_cell))};
    };
    std::unordered_map<QuantizedPoint, std::vector<std::size_t>,
                       QuantizedPointHash> boundary_owners;
    struct BoundarySegment
    {
        Vec3 first;
        Vec3 second;
        Bounds bounds;
        std::size_t primitive = 0;
    };
    std::vector<BoundarySegment> boundary_segments;
    for (std::size_t primitive = 0; primitive < items.size(); ++primitive)
    {
        const PrimitiveMesh surface =
            triangulatePrimitive(items[primitive].output.primitive);
        std::unordered_map<std::uint64_t, std::size_t> edge_counts;
        edge_counts.reserve(surface.faces.size() * 3);
        for (const Face& face : surface.faces)
            for (int edge = 0; edge < 3; ++edge)
                ++edge_counts[edgeKey(face[edge], face[(edge + 1) % 3])];
        std::vector<std::uint32_t> boundary_vertices;
        for (const Face& face : surface.faces)
            for (int edge = 0; edge < 3; ++edge)
                if (edge_counts[edgeKey(face[edge], face[(edge + 1) % 3])] == 1)
                {
                    const auto first = face[edge];
                    const auto second = face[(edge + 1) % 3];
                    boundary_vertices.push_back(first);
                    boundary_vertices.push_back(second);
                    BoundarySegment segment;
                    segment.first = surface.vertices[first];
                    segment.second = surface.vertices[second];
                    segment.bounds.lower =
                        segment.first.cwiseMin(segment.second);
                    segment.bounds.upper =
                        segment.first.cwiseMax(segment.second);
                    segment.primitive = primitive;
                    boundary_segments.push_back(std::move(segment));
                }
        std::sort(boundary_vertices.begin(), boundary_vertices.end());
        boundary_vertices.erase(
            std::unique(boundary_vertices.begin(), boundary_vertices.end()),
            boundary_vertices.end());
        std::unordered_set<QuantizedPoint, QuantizedPointHash> item_cells;
        item_cells.reserve(boundary_vertices.size());
        for (const auto vertex : boundary_vertices)
        {
            const QuantizedPoint cell = quantizePoint(surface.vertices[vertex]);
            for (std::int64_t dx = -1; dx <= 1; ++dx)
                for (std::int64_t dy = -1; dy <= 1; ++dy)
                    for (std::int64_t dz = -1; dz <= 1; ++dz)
                    {
                        const QuantizedPoint neighbor{
                            cell.x + dx, cell.y + dy, cell.z + dz};
                        const auto found = boundary_owners.find(neighbor);
                        if (found == boundary_owners.end()) continue;
                        for (const auto owner : found->second)
                            connect(primitive, owner);
                    }
            item_cells.insert(cell);
        }
        for (const QuantizedPoint& cell : item_cells)
            boundary_owners[cell].push_back(primitive);
    }

    // Boundary vertices alone do not define surface adjacency. A T junction,
    // a partially overlapping edge, or an endpoint landing in the interior of
    // another edge has no coincident vertex pair, although the two current
    // primitives are geometrically adjacent and must receive a merge attempt.
    // Find all touching boundary-segment pairs with a compact static AABB tree.
    // This keeps the adjacency construction close to O(E log E) for the large
    // real-scene models instead of materializing every segment pair.
    struct SegmentBvhNode
    {
        Bounds bounds;
        std::size_t begin = 0;
        std::size_t end = 0;
        std::size_t left = std::numeric_limits<std::size_t>::max();
        std::size_t right = std::numeric_limits<std::size_t>::max();
    };
    std::vector<std::size_t> segment_order(boundary_segments.size());
    std::iota(segment_order.begin(), segment_order.end(), 0);
    std::vector<SegmentBvhNode> segment_bvh;
    segment_bvh.reserve(boundary_segments.size() * 2);
    const auto mergeBounds = [](Bounds& target, const Bounds& source)
    {
        target.lower = target.lower.cwiseMin(source.lower);
        target.upper = target.upper.cwiseMax(source.upper);
    };
    std::function<std::size_t(std::size_t, std::size_t)> buildSegmentBvh =
        [&](const std::size_t begin, const std::size_t end)
    {
        SegmentBvhNode node;
        node.begin = begin;
        node.end = end;
        Bounds centroid_bounds;
        for (std::size_t offset = begin; offset < end; ++offset)
        {
            const Bounds& bounds =
                boundary_segments[segment_order[offset]].bounds;
            mergeBounds(node.bounds, bounds);
            const Vec3 centroid = (bounds.lower + bounds.upper) * 0.5;
            centroid_bounds.lower = centroid_bounds.lower.cwiseMin(centroid);
            centroid_bounds.upper = centroid_bounds.upper.cwiseMax(centroid);
        }
        const std::size_t node_index = segment_bvh.size();
        segment_bvh.push_back(node);
        if (end - begin <= 8) return node_index;
        const Vec3 extent = centroid_bounds.upper - centroid_bounds.lower;
        int axis = 0;
        if (extent.y() > extent.x()) axis = 1;
        if (extent.z() > extent[axis]) axis = 2;
        const std::size_t middle = begin + (end - begin) / 2;
        std::nth_element(
            segment_order.begin() + begin,
            segment_order.begin() + middle,
            segment_order.begin() + end,
            [&](const std::size_t first, const std::size_t second)
            {
                const Bounds& first_bounds = boundary_segments[first].bounds;
                const Bounds& second_bounds = boundary_segments[second].bounds;
                return first_bounds.lower[axis] + first_bounds.upper[axis] <
                       second_bounds.lower[axis] + second_bounds.upper[axis];
            });
        const std::size_t left = buildSegmentBvh(begin, middle);
        const std::size_t right = buildSegmentBvh(middle, end);
        segment_bvh[node_index].left = left;
        segment_bvh[node_index].right = right;
        return node_index;
    };
    const auto boundsOverlap = [](const Bounds& first, const Bounds& second,
                                  const double expansion)
    {
        for (int axis = 0; axis < 3; ++axis)
            if (first.upper[axis] + expansion < second.lower[axis] ||
                second.upper[axis] + expansion < first.lower[axis])
                return false;
        return true;
    };
    const auto segmentDistanceSquared = [](const Vec3& first_begin,
                                           const Vec3& first_end,
                                           const Vec3& second_begin,
                                           const Vec3& second_end)
    {
        const Vec3 first_direction = first_end - first_begin;
        const Vec3 second_direction = second_end - second_begin;
        const Vec3 offset = first_begin - second_begin;
        const double first_length = first_direction.dot(first_direction);
        const double second_length = second_direction.dot(second_direction);
        const double second_projection = second_direction.dot(offset);
        constexpr double epsilon = 1.0e-30;
        double first_parameter = 0.0;
        double second_parameter = 0.0;
        if (first_length <= epsilon && second_length <= epsilon)
            return offset.dot(offset);
        if (first_length <= epsilon)
            second_parameter = std::clamp(
                second_projection / second_length, 0.0, 1.0);
        else
        {
            const double first_projection = first_direction.dot(offset);
            if (second_length <= epsilon)
                first_parameter = std::clamp(
                    -first_projection / first_length, 0.0, 1.0);
            else
            {
                const double coupling =
                    first_direction.dot(second_direction);
                const double denominator =
                    first_length * second_length - coupling * coupling;
                if (std::abs(denominator) > epsilon)
                    first_parameter = std::clamp(
                        (coupling * second_projection -
                         first_projection * second_length) / denominator,
                        0.0, 1.0);
                second_parameter =
                    (coupling * first_parameter + second_projection) /
                    second_length;
                if (second_parameter < 0.0)
                {
                    second_parameter = 0.0;
                    first_parameter = std::clamp(
                        -first_projection / first_length, 0.0, 1.0);
                }
                else if (second_parameter > 1.0)
                {
                    second_parameter = 1.0;
                    first_parameter = std::clamp(
                        (coupling - first_projection) / first_length,
                        0.0, 1.0);
                }
            }
        }
        const Vec3 first_point =
            first_begin + first_direction * first_parameter;
        const Vec3 second_point =
            second_begin + second_direction * second_parameter;
        const Vec3 difference = first_point - second_point;
        return difference.dot(difference);
    };
    if (!boundary_segments.empty())
    {
        buildSegmentBvh(0, boundary_segments.size());
        const double contact_tolerance = boundary_cell * 2.0;
        const double contact_tolerance_squared =
            contact_tolerance * contact_tolerance;
        std::function<void(std::size_t, std::size_t)> querySegmentBvh =
            [&](const std::size_t segment_index, const std::size_t node_index)
        {
            const BoundarySegment& segment = boundary_segments[segment_index];
            const SegmentBvhNode& node = segment_bvh[node_index];
            if (!boundsOverlap(segment.bounds, node.bounds, contact_tolerance))
                return;
            if (node.left == std::numeric_limits<std::size_t>::max())
            {
                for (std::size_t offset = node.begin; offset < node.end; ++offset)
                {
                    const std::size_t other_index = segment_order[offset];
                    if (other_index <= segment_index) continue;
                    const BoundarySegment& other = boundary_segments[other_index];
                    if (segment.primitive == other.primitive ||
                        items[segment.primitive].neighbors.contains(
                            other.primitive))
                        continue;
                    if (!boundsOverlap(
                            segment.bounds, other.bounds, contact_tolerance))
                        continue;
                    if (segmentDistanceSquared(
                            segment.first, segment.second,
                            other.first, other.second) >
                        contact_tolerance_squared)
                        continue;
                    if (connect(segment.primitive, other.primitive))
                        ++profile.segment_contact_adjacencies;
                }
                return;
            }
            querySegmentBvh(segment_index, node.left);
            querySegmentBvh(segment_index, node.right);
        };
        for (std::size_t segment = 0; segment < boundary_segments.size(); ++segment)
            querySegmentBvh(segment, 0);
    }
    for (const auto& item : items)
        profile.initial_adjacencies += item.neighbors.size();
    profile.initial_adjacencies /= 2;
    profile.boundary_adjacencies =
        profile.initial_adjacencies - profile.responsibility_adjacencies;

    Bounds model_bounds;
    for (const Vec3& vertex : responsibility_mesh.vertices)
    {
        model_bounds.lower = model_bounds.lower.cwiseMin(vertex);
        model_bounds.upper = model_bounds.upper.cwiseMax(vertex);
    }
    const Vec3 model_center = (model_bounds.lower + model_bounds.upper) * 0.5;
    const auto representativeNormal = [&](const OutputPrimitive& item,
                                          const PrimitiveMesh& surface)
        -> std::optional<Vec3>
    {
        Vec3 normal;
        if (planarNormal(item.primitive, tolerance, normal))
        {
            if ((planarPoint(item.primitive) - model_center).dot(normal) < 0.0)
                normal *= -1.0;
            return normal;
        }
        Vec3 sum = Vec3::Zero();
        Vec3 point = Vec3::Zero();
        std::size_t point_count = 0;
        for (const Face& face : surface.faces)
        {
            Vec3 doubled_area =
                (surface.vertices[face[1]] - surface.vertices[face[0]]).cross(
                surface.vertices[face[2]] - surface.vertices[face[0]]);
            sum += doubled_area;
            for (const auto vertex : face)
            {
                point += surface.vertices[vertex];
                ++point_count;
            }
        }
        if (sum.norm() <= tolerance || point_count == 0) return std::nullopt;
        normal = sum.normalized();
        point /= static_cast<double>(point_count);
        if ((point - model_center).dot(normal) < 0.0) normal *= -1.0;
        return normal;
    };
    const auto responsibilityFaces = [&](const std::size_t first,
                                         const std::size_t second)
    {
        std::vector<std::uint32_t> faces = items[first].output.source_faces;
        faces.insert(faces.end(), items[second].output.source_faces.begin(),
                     items[second].output.source_faces.end());
        std::sort(faces.begin(), faces.end());
        faces.erase(std::unique(faces.begin(), faces.end()), faces.end());
        return faces;
    };
    const double connection_tolerance = boundary_cell * 2.0;
    const double connection_tolerance_squared =
        connection_tolerance * connection_tolerance;
    const auto surfacesTouch = [&](const Primitive& first_primitive,
                                   const PrimitiveMesh& first_surface,
                                   const Primitive& second_primitive,
                                   const PrimitiveMesh& second_surface)
    {
        for (const Vec3& vertex : first_surface.vertices)
            if (containsPointCached(
                    second_primitive, second_surface,
                    vertex, connection_tolerance))
                return true;
        for (const Vec3& vertex : second_surface.vertices)
            if (containsPointCached(
                    first_primitive, first_surface,
                    vertex, connection_tolerance))
                return true;
        for (const Face& first_face : first_surface.faces)
            for (int first_edge = 0; first_edge < 3; ++first_edge)
            {
                const Vec3& first_begin =
                    first_surface.vertices[first_face[first_edge]];
                const Vec3& first_end = first_surface.vertices[
                    first_face[(first_edge + 1) % 3]];
                for (const Face& second_face : second_surface.faces)
                    for (int second_edge = 0; second_edge < 3; ++second_edge)
                    {
                        const Vec3& second_begin = second_surface.vertices[
                            second_face[second_edge]];
                        const Vec3& second_end = second_surface.vertices[
                            second_face[(second_edge + 1) % 3]];
                        if (segmentDistanceSquared(
                                first_begin, first_end,
                                second_begin, second_end) <=
                            connection_tolerance_squared)
                            return true;
                    }
            }
        const auto segmentTouchesTriangle = [&](const Vec3& begin,
                                                const Vec3& end,
                                                const Vec3& first,
                                                const Vec3& second,
                                                const Vec3& third)
        {
            const Vec3 area_normal =
                (second - first).cross(third - first);
            const double area_length = area_normal.norm();
            if (area_length <= tolerance) return false;
            const Vec3 normal = area_normal / area_length;
            const double begin_distance = (begin - first).dot(normal);
            const double end_distance = (end - first).dot(normal);
            if ((begin_distance > connection_tolerance &&
                 end_distance > connection_tolerance) ||
                (begin_distance < -connection_tolerance &&
                 end_distance < -connection_tolerance))
                return false;
            const double denominator = begin_distance - end_distance;
            if (std::abs(denominator) <= tolerance) return false;
            const double parameter = begin_distance / denominator;
            if (parameter < -1.0e-10 || parameter > 1.0 + 1.0e-10)
                return false;
            const Vec3 intersection = begin + (end - begin) *
                std::clamp(parameter, 0.0, 1.0);
            return (closestPointOnTriangle(
                intersection, first, second, third) - intersection).norm() <=
                connection_tolerance;
        };
        // Vertex and edge-to-edge tests miss the common case where an edge of
        // one finite patch crosses the interior of a triangle in the other.
        // Test both edge/triangle directions so adjacency reconstruction uses
        // actual surface contact instead of only coincident boundaries.
        for (const Face& first_face : first_surface.faces)
            for (int first_edge = 0; first_edge < 3; ++first_edge)
            {
                const Vec3& begin =
                    first_surface.vertices[first_face[first_edge]];
                const Vec3& end = first_surface.vertices[
                    first_face[(first_edge + 1) % 3]];
                for (const Face& second_face : second_surface.faces)
                    if (segmentTouchesTriangle(
                            begin, end,
                            second_surface.vertices[second_face[0]],
                            second_surface.vertices[second_face[1]],
                            second_surface.vertices[second_face[2]]))
                        return true;
            }
        for (const Face& second_face : second_surface.faces)
            for (int second_edge = 0; second_edge < 3; ++second_edge)
            {
                const Vec3& begin =
                    second_surface.vertices[second_face[second_edge]];
                const Vec3& end = second_surface.vertices[
                    second_face[(second_edge + 1) % 3]];
                for (const Face& first_face : first_surface.faces)
                    if (segmentTouchesTriangle(
                            begin, end,
                            first_surface.vertices[first_face[0]],
                            first_surface.vertices[first_face[1]],
                            first_surface.vertices[first_face[2]]))
                        return true;
            }
        return false;
    };
    const auto fitPair = [&](const std::size_t first,
                             const std::size_t second) -> std::optional<Fit>
    {
        ++profile.fit_attempts;
        // A certified analytic round surface is already a tighter and cheaper
        // representation than an arbitrary planar support replacement. Keep it
        // atomic here; otherwise a correctly recognized disk or cylinder can be
        // flattened and later rebuilt as a visibly looser collection of boxes.
        if (isCertifiedRoundSurfaceKind(
                items[first].output.primitive.kind) ||
            isCertifiedRoundSurfaceKind(
                items[second].output.primitive.kind))
        {
            ++profile.unsupported_rejections;
            return std::nullopt;
        }
        // A directed candidate-to-reference distance cannot detect geometry
        // deleted behind another source surface. Preserve explicitly restored
        // openings until a complete enclosure candidate replaces the whole
        // responsibility set.
        if (items[first].output.preserves_cavity_opening ||
            items[second].output.preserves_cavity_opening)
        {
            ++profile.protected_cutout_rejections;
            return std::nullopt;
        }
        std::array<PrimitiveMesh, 2> current_surfaces{
            triangulatePrimitive(items[first].output.primitive),
            triangulatePrimitive(items[second].output.primitive)};
        const auto first_normal = representativeNormal(
            items[first].output, current_surfaces[0]);
        const auto second_normal = representativeNormal(
            items[second].output, current_surfaces[1]);
        // This pass merges planar layers. Flattening differently oriented
        // surfaces onto an arbitrary support plane creates large crossing
        // panels. Non-planar unions must be recognized later as a certified
        // analytic surface or remain unchanged.
        if (!first_normal || !second_normal ||
            first_normal->dot(*second_normal) < 1.0 - 1.0e-8)
        {
            ++profile.unsupported_rejections;
            return std::nullopt;
        }
        struct ExternalNeighborSurface
        {
            std::size_t id = 0;
            PrimitiveMesh surface;
        };
        std::vector<ExternalNeighborSurface> external_neighbors;
        std::unordered_set<std::size_t> external_neighbor_ids =
            items[first].neighbors;
        external_neighbor_ids.insert(
            items[second].neighbors.begin(), items[second].neighbors.end());
        external_neighbor_ids.erase(first);
        external_neighbor_ids.erase(second);
        for (const auto neighbor : external_neighbor_ids)
        {
            if (!items[neighbor].active) continue;
            PrimitiveMesh neighbor_surface =
                triangulatePrimitive(items[neighbor].output.primitive);
            if (!surfacesTouch(
                    items[first].output.primitive, current_surfaces[0],
                    items[neighbor].output.primitive, neighbor_surface) &&
                !surfacesTouch(
                    items[second].output.primitive, current_surfaces[1],
                    items[neighbor].output.primitive, neighbor_surface))
                continue;
            external_neighbors.push_back({neighbor, std::move(neighbor_surface)});
        }
        std::vector<Vec3> points;
        // The current primitives are the geometry being replaced. Always add
        // their vertices so a fitted surface cannot silently drop a stage-2 cap
        // or restored cavity face. Every current primitive already carries a
        // strict source-triangle certificate. Proving that the replacement
        // contains these two current surfaces is therefore transitive and avoids
        // copying and rescanning their ever-growing source-face sets for every
        // rejected neighbor candidate.
        for (const PrimitiveMesh& surface : current_surfaces)
            points.insert(points.end(), surface.vertices.begin(),
                          surface.vertices.end());
        for (const auto index : {first, second})
            for (const auto vertex : items[index].source_vertices)
                if (vertex < responsibility_mesh.vertices.size())
                    points.push_back(responsibility_mesh.vertices[vertex]);
        if (points.size() < 3)
        {
            ++profile.unsupported_rejections;
            return std::nullopt;
        }
        std::vector<Vec3> directions;
        const bool first_is_larger =
            primitiveSurfaceArea(items[first].output.primitive, tolerance) >=
            primitiveSurfaceArea(items[second].output.primitive, tolerance);
        if (first_normal && second_normal)
        {
            directions.push_back(first_is_larger ? *first_normal : *second_normal);
            if ((*first_normal + *second_normal).norm() > tolerance)
                directions.push_back((*first_normal + *second_normal).normalized());
            directions.push_back(first_is_larger ? *second_normal : *first_normal);
        }
        else
        {
            if (first_normal) directions.push_back(*first_normal);
            if (second_normal) directions.push_back(*second_normal);
        }
        // All face directions below are parallel for a certified planar pair;
        // retaining them preserves deterministic direction ordering.
        for (const PrimitiveMesh& surface : current_surfaces)
        {
            struct FaceDirection
            {
                double area = 0.0;
                Vec3 normal = Vec3::Zero();
            };
            std::vector<FaceDirection> face_directions;
            face_directions.reserve(surface.faces.size());
            for (const Face& face : surface.faces)
            {
                const Vec3& a = surface.vertices[face[0]];
                const Vec3& b = surface.vertices[face[1]];
                const Vec3& c = surface.vertices[face[2]];
                Vec3 area_normal = (b - a).cross(c - a);
                const double area = area_normal.norm();
                if (area <= tolerance) continue;
                Vec3 normal = area_normal / area;
                const Vec3 centroid = (a + b + c) / 3.0;
                if ((centroid - model_center).dot(normal) < 0.0) normal *= -1.0;
                face_directions.push_back({area, normal});
            }
            std::sort(face_directions.begin(), face_directions.end(),
                [](const FaceDirection& left, const FaceDirection& right)
                { return left.area > right.area; });
            const std::size_t limit = std::min<std::size_t>(
                face_directions.size(), 8);
            for (std::size_t index = 0; index < limit; ++index)
                directions.push_back(face_directions[index].normal);
        }
        std::vector<Vec3> unique_directions;
        for (const Vec3& direction : directions)
        {
            if (std::any_of(unique_directions.begin(), unique_directions.end(),
                [&](const Vec3& existing)
                { return existing.dot(direction) >= 1.0 - 1.0e-10; }))
                continue;
            unique_directions.push_back(direction);
        }
        struct GeometricCandidate
        {
            Primitive surface;
            std::size_t triangles = 0;
            double estimated_error = 0.0;
            std::vector<std::size_t> absorbed_neighbors;
        };
        std::vector<GeometricCandidate> geometric_candidates;
        for (Vec3 normal : unique_directions)
        {
            if (normal.norm() <= tolerance) continue;
            normal = normal.normalized();
            double support = -std::numeric_limits<double>::infinity();
            for (const Vec3& point : points)
                support = std::max(support, (point - model_center).dot(normal));
            const Vec3 origin = model_center + normal * support;
            double maximum_removed_surface_distance = 0.0;
            for (const Vec3& point : points)
                maximum_removed_surface_distance = std::max(
                    maximum_removed_surface_distance,
                    support - (point - model_center).dot(normal));
            // A local replacement may not move either current patch farther
            // behind its support plane than the permitted open-surface error.
            // The one-sided candidate-to-reference audit below cannot detect a
            // dropped face when the candidate happens to coincide with some
            // other source surface; this symmetric local bound prevents that
            // false certificate and keeps the proxy shell near the geometry it
            // replaces.
            if (maximum_removed_surface_distance >
                maximum_open_error_distance + tolerance)
            {
                ++profile.error_rejections;
                continue;
            }
            std::vector<Vec3> candidate_points = points;
            for (const auto& neighbor : external_neighbors)
            {
                std::vector<Vec3> intersections;
                for (const Face& face : neighbor.surface.faces)
                {
                    std::array<Vec3, 3> triangle{
                        neighbor.surface.vertices[face[0]],
                        neighbor.surface.vertices[face[1]],
                        neighbor.surface.vertices[face[2]]};
                    std::array<double, 3> signed_distances{};
                    bool coplanar = true;
                    for (int corner = 0; corner < 3; ++corner)
                    {
                        signed_distances[corner] =
                            (triangle[corner] - origin).dot(normal);
                        coplanar &= std::abs(signed_distances[corner]) <=
                                    connection_tolerance;
                    }
                    if (coplanar)
                    {
                        intersections.insert(
                            intersections.end(), triangle.begin(), triangle.end());
                        continue;
                    }
                    for (int edge = 0; edge < 3; ++edge)
                    {
                        const int next = (edge + 1) % 3;
                        const double first_distance = signed_distances[edge];
                        const double second_distance = signed_distances[next];
                        if (std::abs(first_distance) <= connection_tolerance)
                            intersections.push_back(triangle[edge]);
                        if ((first_distance < -connection_tolerance &&
                             second_distance > connection_tolerance) ||
                            (first_distance > connection_tolerance &&
                             second_distance < -connection_tolerance))
                        {
                            const double fraction = first_distance /
                                (first_distance - second_distance);
                            intersections.push_back(
                                triangle[edge] +
                                (triangle[next] - triangle[edge]) * fraction);
                        }
                    }
                }
                candidate_points.insert(
                    candidate_points.end(),
                    intersections.begin(), intersections.end());
            }
            const Mat3 basis = orthonormalFrame(normal);
            Mat3 frame;
            frame.col(0) = basis.col(1);
            frame.col(1) = basis.col(2);
            frame.col(2) = normal;
            std::vector<Vec2> projected;
            projected.reserve(candidate_points.size());
            for (const Vec3& point : candidate_points)
            {
                const Vec3 local = frame.transposeMultiply(point - origin);
                projected.emplace_back(local.x(), local.y());
            }
            auto hull = simplifyPolygon(
                convexHull(std::move(projected), tolerance), tolerance);
            if (hull.size() < 3 || !simplePolygon(hull, tolerance)) continue;
            std::vector<std::vector<Vec2>> enclosing_boundaries;
            enclosing_boundaries.push_back(std::move(hull));
            for (auto& boundary : enclosing_boundaries)
            {
                Primitive surface = polygonPrimitive(boundary, origin, frame);
                const PrimitiveMesh candidate_surface =
                    triangulatePrimitive(surface);
                bool contains = true;
                Vec3 candidate_normal;
                if (!planarNormal(surface, tolerance, candidate_normal))
                    contains = false;
                const Vec3 candidate_point = planarPoint(surface);
                // A replacement surface must retain a two-dimensional
                // projection of both current surfaces. Merely placing every
                // vertex behind a support plane is insufficient: an
                // orthogonal box face would collapse to a line, allowing a
                // six-face enclosure to lose one side at zero directed error.
                // This is a coverage check, not an error heuristic.
                for (const PrimitiveMesh& current : current_surfaces)
                {
                    if (!contains) break;
                    for (const Face& face : current.faces)
                    {
                        const Vec3 area_normal =
                            (current.vertices[face[1]] -
                             current.vertices[face[0]]).cross(
                                current.vertices[face[2]] -
                                current.vertices[face[0]]);
                        const double area = area_normal.norm();
                        if (area <= tolerance) continue;
                        if (std::abs(area_normal.dot(candidate_normal)) <=
                            area * 1.0e-10)
                        {
                            contains = false;
                            break;
                        }
                    }
                }
                for (const auto index : {first, second})
                {
                    if (!contains) break;
                    for (const auto vertex_id : items[index].source_vertices)
                    {
                        if (vertex_id >= responsibility_mesh.vertices.size())
                            continue;
                        const Vec3& vertex =
                            responsibility_mesh.vertices[vertex_id];
                        const double signed_distance =
                            (vertex - candidate_point).dot(candidate_normal);
                        if (signed_distance > tolerance ||
                            !containsPointCached(
                                surface, candidate_surface,
                                vertex - candidate_normal * signed_distance,
                                tolerance))
                        {
                            contains = false;
                            break;
                        }
                    }
                }
                if (!contains)
                {
                    ++profile.containment_rejections;
                    continue;
                }
                // Replacing an adjacent pair must not leave an exposed seam.
                // A former neighbor is valid after the replacement when the
                // new finite surface still touches it.  The only exception is
                // a same-facing patch that lies completely behind, and whose
                // entire projection is covered by, the new outward support
                // surface.  Such a patch has become internal and is absorbed
                // immediately, together with its responsibility triangles.
                // Requiring an intersection with every historical neighbor was
                // too strong, while blindly inheriting every old graph edge
                // preserved adjacencies that no longer existed geometrically.
                std::vector<std::size_t> absorbed_neighbors;
                std::unordered_set<std::size_t> visited_neighbors;
                std::queue<std::size_t> neighbor_queue;
                for (const auto& neighbor : external_neighbors)
                    neighbor_queue.push(neighbor.id);
                double absorbed_error = maximum_removed_surface_distance;
                bool preserves_connections = true;
                while (!neighbor_queue.empty())
                {
                    const std::size_t neighbor_id = neighbor_queue.front();
                    neighbor_queue.pop();
                    if (!visited_neighbors.insert(neighbor_id).second ||
                        !items[neighbor_id].active || neighbor_id == first ||
                        neighbor_id == second)
                        continue;
                    const PrimitiveMesh neighbor_surface = triangulatePrimitive(
                        items[neighbor_id].output.primitive);
                    if (surfacesTouch(
                            surface, candidate_surface,
                            items[neighbor_id].output.primitive,
                            neighbor_surface))
                        continue;

                    const auto neighbor_normal = representativeNormal(
                        items[neighbor_id].output, neighbor_surface);
                    bool occluded = neighbor_normal &&
                        candidate_normal.dot(*neighbor_normal) >=
                            1.0 - 1.0e-8;
                    double neighbor_depth = 0.0;
                    if (occluded)
                        for (const Vec3& vertex : neighbor_surface.vertices)
                        {
                            const double signed_distance =
                                (vertex - candidate_point).dot(candidate_normal);
                            neighbor_depth = std::max(
                                neighbor_depth, -signed_distance);
                            if (signed_distance > connection_tolerance ||
                                neighbor_depth >
                                    maximum_open_error_distance + tolerance ||
                                !containsPointCached(
                                    surface, candidate_surface,
                                    vertex - candidate_normal * signed_distance,
                                    connection_tolerance))
                            {
                                occluded = false;
                                break;
                            }
                        }
                    if (!occluded)
                    {
                        preserves_connections = false;
                        break;
                    }
                    absorbed_error = std::max(absorbed_error, neighbor_depth);
                    absorbed_neighbors.push_back(neighbor_id);
                    for (const auto next : items[neighbor_id].neighbors)
                        neighbor_queue.push(next);
                }
                if (!preserves_connections)
                {
                    ++profile.connectivity_rejections;
                    continue;
                }
                std::size_t replaced_triangles =
                    current_surfaces[0].faces.size() +
                    current_surfaces[1].faces.size();
                for (const auto absorbed : absorbed_neighbors)
                    replaced_triangles += triangulatePrimitive(
                        items[absorbed].output.primitive).faces.size();
                if (boundary.size() - 2 > replaced_triangles)
                {
                    ++profile.workload_rejections;
                    continue;
                }
                geometric_candidates.push_back({
                    std::move(surface), boundary.size() - 2,
                    absorbed_error, std::move(absorbed_neighbors)});
            }
        }
        std::sort(geometric_candidates.begin(), geometric_candidates.end(),
            [](const GeometricCandidate& first,
               const GeometricCandidate& second)
            {
                if (first.triangles != second.triangles)
                    return first.triangles < second.triangles;
                return first.estimated_error < second.estimated_error;
            });
        for (auto& candidate : geometric_candidates)
        {
            OutputPrimitive merged;
            merged.primitive = candidate.surface;
            const auto distance_started = std::chrono::steady_clock::now();
            const FilledSurfaceDistanceCertificate certificate =
                certifyFilledSurfaceDistance(
                distance_reference, {merged},
                maximum_open_error_distance + tolerance);
            const double certificate_seconds = std::chrono::duration<double>(
                std::chrono::steady_clock::now() - distance_started).count();
            profile.certificate_seconds += certificate_seconds;
            profile.distance_seconds += certificate_seconds;
            profile.certificate_samples += certificate.samples;
            if (certificate.exceeded)
            {
                ++profile.certified_distance_rejections;
                ++profile.error_rejections;
                continue;
            }
            if (certificate.certified)
            {
                ++profile.certified_distance_accepts;
                return Fit{std::move(candidate.surface), candidate.triangles,
                           candidate.estimated_error,
                           std::move(candidate.absorbed_neighbors)};
            }
            ++profile.fine_distance_evaluations;
            const auto fine_distance_started = std::chrono::steady_clock::now();
            const double error = maximumFilledSurfaceDistance(
                distance_reference, {merged}, error_sample_spacing,
                maximum_open_error_distance + tolerance);
            const double fine_seconds = std::chrono::duration<double>(
                std::chrono::steady_clock::now() - fine_distance_started).count();
            profile.fine_distance_seconds += fine_seconds;
            profile.distance_seconds += fine_seconds;
            if (error > maximum_open_error_distance + tolerance)
            {
                ++profile.error_rejections;
                continue;
            }
            return Fit{std::move(candidate.surface), candidate.triangles,
                       candidate.estimated_error,
                       std::move(candidate.absorbed_neighbors)};
        }
        if (geometric_candidates.empty())
            ++profile.unsupported_rejections;
        return std::nullopt;
    };

    struct CacheEntry
    {
        std::uint64_t first_version = 0;
        std::uint64_t second_version = 0;
        std::optional<Fit> fit;
    };
    std::unordered_map<std::uint64_t, CacheEntry> cache;
    const auto pairKey = [](std::size_t first, std::size_t second)
    {
        if (first > second) std::swap(first, second);
        return (static_cast<std::uint64_t>(first) << 32U) |
               static_cast<std::uint64_t>(second);
    };
    const auto evaluate = [&](std::size_t first, std::size_t second)
        -> std::optional<Fit>
    {
        ++profile.candidate_evaluations;
        if (first > second) std::swap(first, second);
        const auto key = pairKey(first, second);
        const auto found = cache.find(key);
        if (found != cache.end() &&
            found->second.first_version == items[first].version &&
            found->second.second_version == items[second].version)
        {
            ++profile.cache_hits;
            return found->second.fit;
        }
        auto fit = fitPair(first, second);
        cache[key] = {items[first].version, items[second].version, fit};
        return fit;
    };
    // Process the current adjacency graph as a streaming greedy work list.
    // Do not pre-fit and rank every edge: most of those expensive candidates
    // become stale as soon as an earlier accepted merge changes an endpoint.
    // The user contract is local and immediate: evaluate the current adjacent
    // pair, merge it when the measured error passes, otherwise cache that
    // failure until either endpoint version changes.
    std::queue<Candidate> queue;
    const auto enqueue = [&](std::size_t first, std::size_t second)
    {
        if (first > second) std::swap(first, second);
        if (!items[first].active || !items[second].active ||
            !items[first].neighbors.contains(second)) return;
        queue.push({first, second, items[first].version, items[second].version});
    };
    for (std::size_t first = 0; first < items.size(); ++first)
        for (const auto second : items[first].neighbors)
            if (first < second) enqueue(first, second);

    merged_count = 0;
    while (!queue.empty())
    {
        const Candidate candidate = queue.front();
        queue.pop();
        if (!items[candidate.first].active || !items[candidate.second].active ||
            items[candidate.first].version != candidate.first_version ||
            items[candidate.second].version != candidate.second_version ||
            !items[candidate.first].neighbors.contains(candidate.second))
            continue;
        const auto fit = evaluate(candidate.first, candidate.second);
        if (!fit) continue;
        std::size_t replaced_triangles =
            triangulatePrimitive(items[candidate.first].output.primitive).faces.size() +
            triangulatePrimitive(items[candidate.second].output.primitive).faces.size();
        for (const auto absorbed : fit->absorbed_neighbors)
            if (absorbed != candidate.first && absorbed != candidate.second &&
                items[absorbed].active)
                replaced_triangles += triangulatePrimitive(
                    items[absorbed].output.primitive).faces.size();
        if (fit->triangles > replaced_triangles)
        {
            ++profile.workload_rejections;
            continue;
        }
        std::vector<std::size_t> consumed{
            candidate.first, candidate.second};
        consumed.insert(consumed.end(), fit->absorbed_neighbors.begin(),
                        fit->absorbed_neighbors.end());
        std::sort(consumed.begin(), consumed.end());
        consumed.erase(std::unique(consumed.begin(), consumed.end()),
                       consumed.end());
        OutputPrimitive replacement;
        replacement.primitive = fit->surface;
        std::uint64_t common_enclosure =
            items[candidate.first].output.enclosure_group;
        bool common_enclosure_valid = true;
        std::vector<std::uint32_t> merged_source_vertices;
        for (const auto index : consumed)
        {
            replacement.source_faces.insert(
                replacement.source_faces.end(),
                items[index].output.source_faces.begin(),
                items[index].output.source_faces.end());
            merged_source_vertices.insert(
                merged_source_vertices.end(),
                items[index].source_vertices.begin(),
                items[index].source_vertices.end());
            replacement.shallow_shell_coalesced |=
                items[index].output.shallow_shell_coalesced;
            replacement.preserves_cavity_opening |=
                items[index].output.preserves_cavity_opening;
            common_enclosure_valid &=
                items[index].output.enclosure_group == common_enclosure;
        }
        std::sort(replacement.source_faces.begin(),
                  replacement.source_faces.end());
        replacement.source_faces.erase(
            std::unique(replacement.source_faces.begin(),
                        replacement.source_faces.end()),
            replacement.source_faces.end());
        replacement.enclosure_group =
            common_enclosure_valid ? common_enclosure : 0;
        std::sort(merged_source_vertices.begin(), merged_source_vertices.end());
        merged_source_vertices.erase(
            std::unique(merged_source_vertices.begin(),
                        merged_source_vertices.end()),
            merged_source_vertices.end());
        std::vector<std::size_t> affected;
        for (const auto index : consumed)
            affected.insert(affected.end(), items[index].neighbors.begin(),
                            items[index].neighbors.end());
        std::sort(affected.begin(), affected.end());
        affected.erase(std::unique(affected.begin(), affected.end()), affected.end());
        items[candidate.first].output = std::move(replacement);
        items[candidate.first].source_vertices =
            std::move(merged_source_vertices);
        ++items[candidate.first].version;
        for (const auto index : consumed)
        {
            if (index == candidate.first) continue;
            items[index].active = false;
            std::vector<std::uint32_t>().swap(items[index].source_vertices);
            ++items[index].version;
        }
        items[candidate.first].neighbors.clear();
        const PrimitiveMesh replacement_surface = triangulatePrimitive(
            items[candidate.first].output.primitive);
        for (const auto neighbor : affected)
        {
            for (const auto index : consumed)
                items[neighbor].neighbors.erase(index);
            if (std::binary_search(consumed.begin(), consumed.end(), neighbor) ||
                !items[neighbor].active)
                continue;
            const PrimitiveMesh neighbor_surface = triangulatePrimitive(
                items[neighbor].output.primitive);
            if (surfacesTouch(
                    items[candidate.first].output.primitive,
                    replacement_surface,
                    items[neighbor].output.primitive, neighbor_surface))
            {
                items[neighbor].neighbors.insert(candidate.first);
                items[candidate.first].neighbors.insert(neighbor);
            }
            else
                ++profile.stale_neighbors_removed;
        }
        for (const auto index : consumed)
            if (index != candidate.first) items[index].neighbors.clear();
        merged_count += consumed.size() - 1;
        ++profile.accepted_merges;
        profile.absorbed_neighbors += consumed.size() - 2;
        if (fit->triangles > replaced_triangles)
            ++profile.accepted_triangle_increases;
        for (const auto neighbor : items[candidate.first].neighbors)
            enqueue(candidate.first, neighbor);
    }

    // A fixed-point merge is complete only when every still-adjacent pair has
    // been evaluated at its current endpoint versions and none satisfies the
    // error limit. This audit guards against stale-queue and missed-neighbor
    // bugs: silently returning while an acceptable adjacent merge remains would
    // violate the stage-3 contract.
    for (std::size_t first = 0; first < items.size(); ++first)
    {
        if (!items[first].active) continue;
        for (const auto second : items[first].neighbors)
        {
            if (first >= second || !items[second].active) continue;
            if (evaluate(first, second))
                ++profile.remaining_acceptable_candidates;
        }
    }
    if (profile.remaining_acceptable_candidates != 0)
        throw std::runtime_error(
            "adjacent surface merge stopped before reaching a fixed point");

    std::vector<OutputPrimitive> result;
    result.reserve(items.size() - merged_count);
    for (auto& item : items)
        if (item.active) result.push_back(std::move(item.output));

    std::ofstream profile_stream(profile_path);
    profile_stream << std::setprecision(17)
        << "{\"complete\":true"
        << ",\"strategy\":\"streaming_greedy_adjacent_fixed_point\""
        << ",\"input_primitives\":" << profile.input_primitives
        << ",\"initial_adjacencies\":" << profile.initial_adjacencies
        << ",\"responsibility_adjacencies\":"
        << profile.responsibility_adjacencies
        << ",\"boundary_adjacencies\":" << profile.boundary_adjacencies
        << ",\"segment_contact_adjacencies\":"
        << profile.segment_contact_adjacencies
        << ",\"candidate_evaluations\":" << profile.candidate_evaluations
        << ",\"cache_hits\":" << profile.cache_hits
        << ",\"fit_attempts\":" << profile.fit_attempts
        << ",\"containment_rejections\":" << profile.containment_rejections
        << ",\"connectivity_rejections\":"
        << profile.connectivity_rejections
        << ",\"error_rejections\":" << profile.error_rejections
        << ",\"workload_rejections\":" << profile.workload_rejections
        << ",\"protected_cutout_rejections\":"
        << profile.protected_cutout_rejections
        << ",\"unsupported_rejections\":" << profile.unsupported_rejections
        << ",\"accepted_merges\":" << profile.accepted_merges
        << ",\"absorbed_neighbors\":" << profile.absorbed_neighbors
        << ",\"stale_neighbors_removed\":"
        << profile.stale_neighbors_removed
        << ",\"accepted_triangle_increases\":"
        << profile.accepted_triangle_increases
        << ",\"remaining_acceptable_candidates\":"
        << profile.remaining_acceptable_candidates
        << ",\"certified_distance_accepts\":"
        << profile.certified_distance_accepts
        << ",\"certified_distance_rejections\":"
        << profile.certified_distance_rejections
        << ",\"certificate_samples\":" << profile.certificate_samples
        << ",\"fine_distance_evaluations\":"
        << profile.fine_distance_evaluations
        << ",\"certificate_seconds\":"
        << profile.certificate_seconds
        << ",\"fine_distance_seconds\":"
        << profile.fine_distance_seconds
        << ",\"distance_seconds\":" << profile.distance_seconds
        << ",\"output_primitives\":" << result.size()
        << ",\"total_seconds\":" << std::chrono::duration<double>(
            std::chrono::steady_clock::now() - started).count() << "}\n";
    return result;
}

std::vector<OutputPrimitive> mergeLocalCoplanarPrimitives(
    const Mesh& source_mesh,
    std::vector<OutputPrimitive> primitives,
    const double tolerance,
    const double maximum_open_error_distance,
    const double error_sample_spacing,
    std::size_t& merged_count,
    const std::filesystem::path& profile_path)
{
    const auto merge_started = std::chrono::steady_clock::now();
    struct MergeProfile
    {
        std::size_t input_primitives = 0;
        std::size_t plane_groups = 0;
        std::size_t normal_groups = 0;
        std::size_t maximum_plane_group_size = 0;
        std::size_t sweep_pairs_visited = 0;
        std::size_t adjacency_tests = 0;
        std::size_t adjacency_passes = 0;
        std::size_t evaluate_calls = 0;
        std::size_t evaluation_cache_hits = 0;
        std::size_t evaluation_cache_misses = 0;
        std::size_t hull_attempts = 0;
        std::size_t surface_hulls = 0;
        std::size_t hausdorff_calls = 0;
        std::size_t rejected_cavity_opening = 0;
        std::size_t rejected_enclosure_group = 0;
        std::size_t rejected_nonplanar = 0;
        std::size_t rejected_plane_mismatch = 0;
        std::size_t rejected_adjacency = 0;
        std::size_t rejected_nonquadrilateral_hull = 0;
        std::size_t rejected_nonrectangle_hull = 0;
        std::size_t rejected_hausdorff = 0;
        std::size_t recompute_calls = 0;
        std::size_t recompute_partner_scans = 0;
        std::size_t accepted_merges = 0;
        double plane_grouping_seconds = 0.0;
        double hausdorff_seconds = 0.0;
    } profile;
    profile.input_primitives = primitives.size();
    auto last_profile_flush = merge_started - std::chrono::seconds(2);
    const auto writeProfile = [&](const bool force = false,
                                  const bool complete = false)
    {
        const auto now = std::chrono::steady_clock::now();
        if (!force && now - last_profile_flush < std::chrono::seconds(1)) return;
        last_profile_flush = now;
        std::ofstream stream(profile_path);
        stream << std::setprecision(17)
               << "{\"complete\":" << (complete ? "true" : "false")
               << ",\"input_primitives\":" << profile.input_primitives
               << ",\"plane_groups\":" << profile.plane_groups
               << ",\"normal_groups\":" << profile.normal_groups
               << ",\"maximum_plane_group_size\":"
               << profile.maximum_plane_group_size
               << ",\"plane_grouping_seconds\":"
               << profile.plane_grouping_seconds
               << ",\"sweep_pairs_visited\":"
               << profile.sweep_pairs_visited
               << ",\"adjacency_tests\":" << profile.adjacency_tests
               << ",\"adjacency_passes\":" << profile.adjacency_passes
               << ",\"evaluate_calls\":" << profile.evaluate_calls
               << ",\"evaluation_cache_hits\":"
               << profile.evaluation_cache_hits
               << ",\"evaluation_cache_misses\":"
               << profile.evaluation_cache_misses
               << ",\"hull_attempts\":" << profile.hull_attempts
               << ",\"surface_hulls\":" << profile.surface_hulls
               << ",\"hausdorff_calls\":" << profile.hausdorff_calls
               << ",\"rejected_cavity_opening\":"
               << profile.rejected_cavity_opening
               << ",\"rejected_enclosure_group\":"
               << profile.rejected_enclosure_group
               << ",\"rejected_nonplanar\":" << profile.rejected_nonplanar
               << ",\"rejected_plane_mismatch\":"
               << profile.rejected_plane_mismatch
               << ",\"rejected_adjacency\":" << profile.rejected_adjacency
               << ",\"rejected_nonquadrilateral_hull\":"
               << profile.rejected_nonquadrilateral_hull
               << ",\"rejected_nonrectangle_hull\":"
               << profile.rejected_nonrectangle_hull
               << ",\"rejected_hausdorff\":" << profile.rejected_hausdorff
               << ",\"hausdorff_seconds\":" << profile.hausdorff_seconds
               << ",\"recompute_calls\":" << profile.recompute_calls
               << ",\"recompute_partner_scans\":"
               << profile.recompute_partner_scans
               << ",\"accepted_merges\":" << profile.accepted_merges
               << ",\"total_seconds\":"
               << std::chrono::duration<double>(now - merge_started).count()
               << "}\n";
    };
    writeProfile();
    merged_count = 0;
    struct Item
    {
        OutputPrimitive output;
        std::uint64_t version = 0;
        bool active = true;
    };
    struct MergeFit
    {
        double excess = 0.0;
        Primitive surface;
    };
    struct MergeCandidate
    {
        double excess = 0.0;
        std::size_t first = 0;
        std::size_t second = 0;
        std::uint64_t first_version = 0;
        std::uint64_t second_version = 0;
        std::size_t owner = 0;
        std::uint64_t owner_generation = 0;
    };
    struct MergeCandidateGreater
    {
        bool operator()(const MergeCandidate& left,
                        const MergeCandidate& right) const
        {
            if (left.excess != right.excess) return left.excess > right.excess;
            if (left.first != right.first) return left.first > right.first;
            return left.second > right.second;
        }
    };

    std::vector<Item> items;
    items.reserve(primitives.size());
    for (auto& primitive : primitives)
        items.push_back({std::move(primitive)});

    const auto primitiveBounds = [&](const Primitive& primitive)
    {
        Bounds bounds;
        for (const Vec3& corner : planarCorners(primitive))
        {
            bounds.lower = bounds.lower.cwiseMin(corner);
            bounds.upper = bounds.upper.cwiseMax(corner);
        }
        return bounds;
    };
    std::vector<Bounds> item_bounds(items.size());
    for (std::size_t index = 0; index < items.size(); ++index)
        item_bounds[index] = primitiveBounds(items[index].output.primitive);
    const auto potentiallyAdjacent = [&](const std::size_t first,
                                         const std::size_t second)
    {
        ++profile.adjacency_tests;
        const Bounds& first_bounds = item_bounds[first];
        const Bounds& second_bounds = item_bounds[second];
        Vec3 separation = Vec3::Zero();
        for (int axis = 0; axis < 3; ++axis)
        {
            if (first_bounds.upper[axis] < second_bounds.lower[axis])
                separation[axis] =
                    second_bounds.lower[axis] - first_bounds.upper[axis];
            else if (second_bounds.upper[axis] < first_bounds.lower[axis])
                separation[axis] =
                    first_bounds.lower[axis] - second_bounds.upper[axis];
        }
        const double error_controlled_distance = std::max(
            tolerance * 8.0, 2.0 * maximum_open_error_distance + tolerance);
        const bool adjacent = separation.norm() <= error_controlled_distance;
        if (adjacent) ++profile.adjacency_passes;
        return adjacent;
    };

    const auto evaluateUncached = [&](const std::size_t first,
                                       const std::size_t second) -> std::optional<MergeFit>
    {
        if (first == second || !items[first].active || !items[second].active)
            return std::nullopt;
        // A retained cavity opening is already represented by an exact
        // constrained planar triangulation.  Replacing any of those triangles
        // with a pairwise bounding rectangle can bridge the opening and can
        // extend far beyond the source silhouette.
        if (items[first].output.preserves_cavity_opening ||
            items[second].output.preserves_cavity_opening)
        {
            ++profile.rejected_cavity_opening;
            return std::nullopt;
        }
        const std::uint64_t first_enclosure =
            items[first].output.enclosure_group;
        const std::uint64_t second_enclosure =
            items[second].output.enclosure_group;
        if (first_enclosure != second_enclosure &&
            (first_enclosure != 0 || second_enclosure != 0))
        {
            ++profile.rejected_enclosure_group;
            return std::nullopt;
        }
        const Primitive& first_primitive = items[first].output.primitive;
        const Primitive& second_primitive = items[second].output.primitive;
        Vec3 normal;
        Vec3 second_normal;
        if (!planarNormal(first_primitive, tolerance, normal) ||
            !planarNormal(second_primitive, tolerance, second_normal))
        {
            ++profile.rejected_nonplanar;
            return std::nullopt;
        }
        if (std::abs(normal.dot(second_normal)) < 1.0 - 1.0e-8)
        {
            ++profile.rejected_plane_mismatch;
            return std::nullopt;
        }
        const Vec3 origin = planarPoint(first_primitive);
        if (std::abs((planarPoint(second_primitive) - origin).dot(normal)) > tolerance)
        {
            ++profile.rejected_plane_mismatch;
            return std::nullopt;
        }

        if (!potentiallyAdjacent(first, second))
        {
            ++profile.rejected_adjacency;
            return std::nullopt;
        }
        const Mat3 basis = orthonormalFrame(normal);
        Mat3 frame;
        frame.col(0) = basis.col(1);
        frame.col(1) = basis.col(2);
        frame.col(2) = normal;
        std::vector<Vec2> points;
        for (const auto index : {first, second})
            for (const Vec3& corner : planarCorners(items[index].output.primitive))
            {
                const Vec3 local = frame.transposeMultiply(corner - origin);
                points.emplace_back(local.x(), local.y());
            }
        ++profile.hull_attempts;
        auto hull = simplifyPolygon(convexHull(std::move(points), tolerance), tolerance);
        if (hull.size() < 3 || !simplePolygon(hull, tolerance))
        {
            ++profile.rejected_nonquadrilateral_hull;
            return std::nullopt;
        }
        ++profile.surface_hulls;
        Primitive surface = polygonPrimitive(hull, origin, frame);
        OutputPrimitive merged;
        merged.primitive = surface;
        merged.source_faces = items[first].output.source_faces;
        merged.source_faces.insert(
            merged.source_faces.end(),
            items[second].output.source_faces.begin(),
            items[second].output.source_faces.end());
        const auto hausdorff_started = std::chrono::steady_clock::now();
        ++profile.hausdorff_calls;
        const double maximum_distance = maximumFilledSurfaceDistance(
            source_mesh, {merged}, error_sample_spacing,
            maximum_open_error_distance + tolerance);
        profile.hausdorff_seconds += std::chrono::duration<double>(
            std::chrono::steady_clock::now() - hausdorff_started).count();
        if (maximum_distance > maximum_open_error_distance + tolerance)
        {
            ++profile.rejected_hausdorff;
            return std::nullopt;
        }
        return MergeFit{maximum_distance, std::move(surface)};
    };
    struct CachedEvaluation
    {
        std::uint64_t first_version = 0;
        std::uint64_t second_version = 0;
        std::optional<MergeFit> fit;
    };
    std::unordered_map<std::uint64_t, CachedEvaluation> evaluation_cache;
    evaluation_cache.reserve(items.size() * 2);
    const auto evaluate = [&](std::size_t first,
                              std::size_t second) -> std::optional<MergeFit>
    {
        ++profile.evaluate_calls;
        writeProfile();
        if (first > second) std::swap(first, second);
        const std::uint64_t key =
            (static_cast<std::uint64_t>(first) << 32U) |
            static_cast<std::uint64_t>(second);
        const auto found = evaluation_cache.find(key);
        if (found != evaluation_cache.end() &&
            found->second.first_version == items[first].version &&
            found->second.second_version == items[second].version)
        {
            ++profile.evaluation_cache_hits;
            return found->second.fit;
        }
        ++profile.evaluation_cache_misses;
        auto fit = evaluateUncached(first, second);
        evaluation_cache[key] = {
            items[first].version, items[second].version, fit};
        return fit;
    };

    // Coplanar planes are independent merge domains. Quantized analytic plane
    // keys replace the former all-pairs disjoint-set construction; splitting a
    // numerically borderline plane is safe, while quadratic grouping is not.
    using PlaneKey = std::array<std::int64_t, 4>;
    const auto grouping_started = std::chrono::steady_clock::now();
    std::map<PlaneKey, std::vector<std::size_t>> plane_groups;
    std::set<std::array<std::int64_t, 3>> normal_groups;
    constexpr double normal_quantization = 1.0e8;
    const double distance_quantization = 1.0 /
        std::max(tolerance * 2.0, 1.0e-12);
    for (std::size_t index = 0; index < items.size(); ++index)
    {
        Vec3 normal;
        if (!planarNormal(items[index].output.primitive, tolerance, normal))
            continue;
        int dominant = 0;
        for (int axis = 1; axis < 3; ++axis)
            if (std::abs(normal[axis]) > std::abs(normal[dominant]))
                dominant = axis;
        if (normal[dominant] < 0.0) normal *= -1.0;
        const double distance = normal.dot(
            planarPoint(items[index].output.primitive));
        const PlaneKey key{{
            static_cast<std::int64_t>(std::llround(
                normal.x() * normal_quantization)),
            static_cast<std::int64_t>(std::llround(
                normal.y() * normal_quantization)),
            static_cast<std::int64_t>(std::llround(
                normal.z() * normal_quantization)),
            static_cast<std::int64_t>(std::llround(
                distance * distance_quantization))}};
        plane_groups[key].push_back(index);
        normal_groups.insert({key[0], key[1], key[2]});
    }
    profile.plane_groups = plane_groups.size();
    profile.normal_groups = normal_groups.size();
    for (const auto& [key, group] : plane_groups)
    {
        (void)key;
        profile.maximum_plane_group_size = std::max(
            profile.maximum_plane_group_size, group.size());
    }
    profile.plane_grouping_seconds = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - grouping_started).count();
    writeProfile(true);

    for (auto& [plane_key, group] : plane_groups)
    {
        (void)plane_key;
        if (group.size() < 2) continue;
        std::priority_queue<MergeCandidate, std::vector<MergeCandidate>,
                            MergeCandidateGreater> candidates;
        const std::size_t no_partner = std::numeric_limits<std::size_t>::max();
        std::vector<std::size_t> best_partner(items.size(), no_partner);
        std::vector<double> best_excess(
            items.size(), std::numeric_limits<double>::infinity());
        std::vector<std::uint64_t> best_generation(items.size(), 0);
        const auto canonicalPair = [](std::size_t first, std::size_t second)
        {
            if (first > second) std::swap(first, second);
            return std::pair{first, second};
        };
        const auto isBetter = [&](const std::size_t owner,
                                  const std::size_t partner,
                                  const double excess)
        {
            if (excess != best_excess[owner])
                return excess < best_excess[owner];
            if (best_partner[owner] == no_partner) return true;
            return canonicalPair(owner, partner) <
                   canonicalPair(owner, best_partner[owner]);
        };
        const auto consider = [&](const std::size_t owner,
                                  const std::size_t partner,
                                  const MergeFit& fit)
        {
            if (!isBetter(owner, partner, fit.excess)) return false;
            best_partner[owner] = partner;
            best_excess[owner] = fit.excess;
            return true;
        };
        const auto pushCurrent = [&](const std::size_t owner)
        {
            const std::size_t partner = best_partner[owner];
            if (!items[owner].active || partner == no_partner ||
                !items[partner].active)
                return;
            const auto [first, second] = canonicalPair(owner, partner);
            candidates.push({best_excess[owner], first, second,
                             items[first].version, items[second].version,
                             owner, best_generation[owner]});
        };
        int sweep_axis = 0;
        Vec3 center_lower = Vec3::Constant(
            std::numeric_limits<double>::infinity());
        Vec3 center_upper = Vec3::Constant(
            -std::numeric_limits<double>::infinity());
        for (const std::size_t id : group)
        {
            const Vec3 center =
                (item_bounds[id].lower + item_bounds[id].upper) * 0.5;
            center_lower = center_lower.cwiseMin(center);
            center_upper = center_upper.cwiseMax(center);
        }
        for (int axis = 1; axis < 3; ++axis)
            if (center_upper[axis] - center_lower[axis] >
                center_upper[sweep_axis] - center_lower[sweep_axis])
                sweep_axis = axis;
        std::sort(group.begin(), group.end(), [&](const auto first,
                                                  const auto second)
        {
            return item_bounds[first].lower[sweep_axis] <
                   item_bounds[second].lower[sweep_axis];
        });
        for (std::size_t first = 0; first < group.size(); ++first)
        {
            const std::size_t first_id = group[first];
            const double sweep_gap = std::max(
                tolerance * 8.0,
                2.0 * maximum_open_error_distance + tolerance);
            for (std::size_t second = first + 1;
                 second < group.size(); ++second)
            {
                ++profile.sweep_pairs_visited;
                const std::size_t second_id = group[second];
                if (item_bounds[second_id].lower[sweep_axis] >
                    item_bounds[first_id].upper[sweep_axis] + sweep_gap)
                    break;
                if (!potentiallyAdjacent(first_id, second_id)) continue;
                const auto fit = evaluate(first_id, second_id);
                if (!fit) continue;
                (void)consider(first_id, second_id, *fit);
                (void)consider(second_id, first_id, *fit);
            }
        }
        for (const std::size_t owner : group) pushCurrent(owner);

        const auto recompute = [&](const std::size_t owner)
        {
            ++profile.recompute_calls;
            best_partner[owner] = no_partner;
            best_excess[owner] = std::numeric_limits<double>::infinity();
            ++best_generation[owner];
            if (!items[owner].active) return;
            for (const std::size_t partner : group)
            {
                ++profile.recompute_partner_scans;
                if (partner == owner || !items[partner].active) continue;
                if (!potentiallyAdjacent(owner, partner)) continue;
                const auto [first, second] = canonicalPair(owner, partner);
                const auto fit = evaluate(first, second);
                if (fit) (void)consider(owner, partner, *fit);
            }
            pushCurrent(owner);
        };

        while (!candidates.empty())
        {
            const MergeCandidate candidate = candidates.top();
            candidates.pop();
            if (best_generation[candidate.owner] != candidate.owner_generation ||
                best_partner[candidate.owner] == no_partner ||
                canonicalPair(candidate.owner, best_partner[candidate.owner]) !=
                    std::pair{candidate.first, candidate.second} ||
                !items[candidate.first].active || !items[candidate.second].active ||
                items[candidate.first].version != candidate.first_version ||
                items[candidate.second].version != candidate.second_version)
                continue;
            const auto fit = evaluate(candidate.first, candidate.second);
            if (!fit) continue;
            OutputPrimitive replacement;
            replacement.primitive = fit->surface;
            replacement.source_faces = items[candidate.first].output.source_faces;
            replacement.source_faces.insert(
                replacement.source_faces.end(),
                items[candidate.second].output.source_faces.begin(),
                items[candidate.second].output.source_faces.end());
            replacement.shallow_shell_coalesced =
                items[candidate.first].output.shallow_shell_coalesced ||
                items[candidate.second].output.shallow_shell_coalesced;
            replacement.preserves_cavity_opening =
                items[candidate.first].output.preserves_cavity_opening ||
                items[candidate.second].output.preserves_cavity_opening;
            replacement.enclosure_group =
                items[candidate.first].output.enclosure_group;
            items[candidate.first].output = std::move(replacement);
            item_bounds[candidate.first] = primitiveBounds(
                items[candidate.first].output.primitive);
            ++items[candidate.first].version;
            items[candidate.second].active = false;
            ++items[candidate.second].version;
            ++merged_count;
            ++profile.accepted_merges;
            ++best_generation[candidate.second];
            best_partner[candidate.second] = no_partner;
            std::size_t active_count = 0;
            for (const std::size_t other : group)
            {
                if (!items[other].active) continue;
                ++active_count;
                if (other == candidate.first ||
                    best_partner[other] == candidate.first ||
                    best_partner[other] == candidate.second)
                {
                    recompute(other);
                    continue;
                }
                if (!potentiallyAdjacent(other, candidate.first)) continue;
                const auto [first, second] = canonicalPair(other, candidate.first);
                const auto new_fit = evaluate(first, second);
                if (new_fit && consider(other, candidate.first, *new_fit))
                {
                    ++best_generation[other];
                    pushCurrent(other);
                }
            }
            if (candidates.size() > 4 * active_count + 64)
            {
                decltype(candidates) compact;
                candidates.swap(compact);
                for (const std::size_t owner : group)
                    if (items[owner].active) pushCurrent(owner);
            }
        }
    }
    std::vector<OutputPrimitive> result;
    result.reserve(items.size() - merged_count);
    for (auto& item : items)
        if (item.active) result.push_back(std::move(item.output));
    writeProfile(true, true);
    return result;
}

struct CoplanarCanonicalizationStats
{
    std::size_t groups = 0;
    std::size_t removed_primitives = 0;
    double removed_overlap_area = 0.0;
};

std::vector<OutputPrimitive> canonicalizeCoplanarPrimitiveUnion(
    const Mesh& source_mesh,
    std::vector<OutputPrimitive> primitives,
    const double tolerance,
    CoplanarCanonicalizationStats& stats)
{
    stats = {};
    std::vector<bool> consumed(primitives.size(), false);
    std::vector<OutputPrimitive> result;
    result.reserve(primitives.size());
    constexpr int clipper_precision = 8;

    struct PlaneKey
    {
        std::array<std::int64_t, 4> value{};
        bool operator<(const PlaneKey& other) const { return value < other.value; }
    };
    const double angular_quantum = 1.0e-8;
    const double distance_quantum = std::max(tolerance, 1.0e-12);
    std::map<PlaneKey, std::vector<std::size_t>> plane_buckets;
    std::vector<Vec3> plane_normals(primitives.size());
    for (std::size_t index = 0; index < primitives.size(); ++index)
    {
        Vec3 normal;
        if (!planarNormal(primitives[index].primitive, tolerance, normal)) continue;
        int dominant_axis = 0;
        for (int axis = 1; axis < 3; ++axis)
            if (std::abs(normal[axis]) > std::abs(normal[dominant_axis])) dominant_axis = axis;
        if (normal[dominant_axis] < 0.0) normal = -normal;
        plane_normals[index] = normal;
        const double distance = normal.dot(planarPoint(primitives[index].primitive));
        PlaneKey key;
        for (int axis = 0; axis < 3; ++axis)
            key.value[axis] = static_cast<std::int64_t>(
                std::llround(normal[axis] / angular_quantum));
        key.value[3] = static_cast<std::int64_t>(
            std::llround(distance / distance_quantum));
        plane_buckets[key].push_back(index);
    }
    constexpr std::size_t no_component = std::numeric_limits<std::size_t>::max();
    std::vector<std::size_t> planar_component(primitives.size(), no_component);
    std::vector<std::vector<std::size_t>> planar_components;

    while (true)
    {
        const auto seed_iterator = std::find(consumed.begin(), consumed.end(), false);
        if (seed_iterator == consumed.end()) break;
        const std::size_t seed = static_cast<std::size_t>(
            std::distance(consumed.begin(), seed_iterator));
        if (plane_normals[seed].norm() == 0.0)
        {
            consumed[seed] = true;
            result.push_back(std::move(primitives[seed]));
            continue;
        }
        const Vec3 normal = plane_normals[seed];
        const Vec3 origin = planarPoint(primitives[seed].primitive);
        const Mat3 basis = orthonormalFrame(normal);
        Mat3 frame;
        frame.col(0) = static_cast<Vec3>(basis.col(1));
        frame.col(1) = static_cast<Vec3>(basis.col(2));
        frame.col(2) = normal;

        std::vector<std::size_t> group;
        Clipper2Lib::PathsD subjects;
        std::size_t input_triangle_count = 0;
        double input_area = 0.0;
        std::vector<std::uint32_t> source_faces;
        bool shallow_shell_coalesced = false;
        bool preserves_cavity_opening = false;
        std::uint64_t enclosure_group = 0;
        bool enclosure_initialized = false;
        bool enclosure_consistent = true;
        if (planar_component[seed] == no_component)
        {
            const double distance = normal.dot(origin);
            PlaneKey center_key;
            for (int axis = 0; axis < 3; ++axis)
                center_key.value[axis] = static_cast<std::int64_t>(
                    std::llround(normal[axis] / angular_quantum));
            center_key.value[3] = static_cast<std::int64_t>(
                std::llround(distance / distance_quantum));
            std::vector<std::size_t> candidates;
            for (int first_delta = -1; first_delta <= 1; ++first_delta)
                for (int second_delta = -1; second_delta <= 1; ++second_delta)
                    for (int third_delta = -1; third_delta <= 1; ++third_delta)
                        for (int distance_delta = -1; distance_delta <= 1; ++distance_delta)
                        {
                            PlaneKey key = center_key;
                            key.value[0] += first_delta;
                            key.value[1] += second_delta;
                            key.value[2] += third_delta;
                            key.value[3] += distance_delta;
                            const auto found = plane_buckets.find(key);
                            if (found != plane_buckets.end())
                                candidates.insert(candidates.end(), found->second.begin(),
                                                  found->second.end());
                        }
            std::sort(candidates.begin(), candidates.end());
            candidates.erase(std::unique(candidates.begin(), candidates.end()), candidates.end());
            struct ProjectedCandidate
            {
                std::size_t primitive_index = 0;
                Bounds2 bounds;
            };
            std::vector<ProjectedCandidate> projected_candidates;
            projected_candidates.reserve(candidates.size());
            for (const std::size_t index : candidates)
            {
                if (consumed[index] || planar_component[index] != no_component) continue;
                const Vec3 candidate_normal = plane_normals[index];
                if (std::abs(candidate_normal.dot(normal)) < 1.0 - 1.0e-8 ||
                    std::abs((planarPoint(primitives[index].primitive) - origin).dot(normal)) >
                        tolerance) continue;

                const PrimitiveMesh mesh = triangulatePrimitive(primitives[index].primitive);
                Bounds2 bounds;
                for (const Vec3& vertex : mesh.vertices)
                {
                    const Vec3 local = frame.transposeMultiply(vertex - origin);
                    bounds.lower = bounds.lower.cwiseMin({local.x(), local.y()});
                    bounds.upper = bounds.upper.cwiseMax({local.x(), local.y()});
                }
                projected_candidates.push_back({index, bounds});
            }

            // A quantized plane bucket can contain thousands of mutually distant
            // CAD patches. Sending all of them through one polygon union makes the
            // Boolean engine combine unrelated path arrangements and can require
            // quadratic memory. Build and cache planar connected components first.
            // The sweep only compares rectangles whose x intervals overlap; the y
            // test then joins patches that overlap or touch within the tolerance.
            std::vector<std::size_t> parent(projected_candidates.size());
            std::iota(parent.begin(), parent.end(), 0);
            const auto findRoot = [&](std::size_t item)
            {
                std::size_t root = item;
                while (parent[root] != root) root = parent[root];
                while (parent[item] != item)
                {
                    const std::size_t next = parent[item];
                    parent[item] = root;
                    item = next;
                }
                return root;
            };
            const auto unite = [&](const std::size_t first, const std::size_t second)
            {
                const std::size_t first_root = findRoot(first);
                const std::size_t second_root = findRoot(second);
                if (first_root != second_root) parent[second_root] = first_root;
            };
            std::vector<std::size_t> order(projected_candidates.size());
            std::iota(order.begin(), order.end(), 0);
            std::sort(order.begin(), order.end(), [&](const std::size_t first,
                                                      const std::size_t second)
            {
                return projected_candidates[first].bounds.lower.x() <
                       projected_candidates[second].bounds.lower.x();
            });
            std::vector<std::size_t> active;
            for (const std::size_t current : order)
            {
                const Bounds2& current_bounds = projected_candidates[current].bounds;
                active.erase(std::remove_if(active.begin(), active.end(),
                    [&](const std::size_t other)
                    {
                        return projected_candidates[other].bounds.upper.x() + tolerance <
                               current_bounds.lower.x();
                    }), active.end());
                for (const std::size_t other : active)
                {
                    const Bounds2& other_bounds = projected_candidates[other].bounds;
                    if (other_bounds.upper.y() + tolerance < current_bounds.lower.y() ||
                        current_bounds.upper.y() + tolerance < other_bounds.lower.y()) continue;
                    unite(current, other);
                }
                active.push_back(current);
            }

            std::map<std::size_t, std::size_t> component_for_root;
            for (std::size_t local = 0; local < projected_candidates.size(); ++local)
            {
                const std::size_t root = findRoot(local);
                auto [iterator, inserted] = component_for_root.emplace(
                    root, planar_components.size());
                if (inserted) planar_components.emplace_back();
                const std::size_t component = iterator->second;
                const std::size_t index = projected_candidates[local].primitive_index;
                planar_component[index] = component;
                planar_components[component].push_back(index);
            }
        }

        if (planar_component[seed] == no_component)
        {
            consumed[seed] = true;
            result.push_back(std::move(primitives[seed]));
            continue;
        }
        group = planar_components[planar_component[seed]];
        for (const std::size_t index : group)
        {
            if (consumed[index]) continue;

            const PrimitiveMesh mesh = triangulatePrimitive(primitives[index].primitive);
            input_triangle_count += mesh.faces.size();
            source_faces.insert(source_faces.end(), primitives[index].source_faces.begin(),
                                primitives[index].source_faces.end());
            shallow_shell_coalesced |= primitives[index].shallow_shell_coalesced;
            preserves_cavity_opening |= primitives[index].preserves_cavity_opening;
            if (!enclosure_initialized)
            {
                enclosure_group = primitives[index].enclosure_group;
                enclosure_initialized = true;
            }
            else if (enclosure_group != primitives[index].enclosure_group)
                enclosure_consistent = false;
            for (const Face& face : mesh.faces)
            {
                Clipper2Lib::PathD path;
                for (const auto vertex : face)
                {
                    const Vec3 local = frame.transposeMultiply(mesh.vertices[vertex] - origin);
                    path.emplace_back(local.x(), local.y());
                }
                const double area = Clipper2Lib::Area(path);
                if (std::abs(area) <= tolerance * tolerance) continue;
                if (area < 0.0) std::reverse(path.begin(), path.end());
                input_area += std::abs(area);
                subjects.push_back(std::move(path));
            }
        }
        if (group.size() < 2 || subjects.empty())
        {
            for (const std::size_t index : group)
            {
                consumed[index] = true;
                result.push_back(std::move(primitives[index]));
            }
            if (group.empty())
            {
                consumed[seed] = true;
                result.push_back(std::move(primitives[seed]));
            }
            continue;
        }

        const Clipper2Lib::PathsD covered = Clipper2Lib::Union(
            subjects, Clipper2Lib::FillRule::NonZero, clipper_precision);
        const double union_area = std::abs(Clipper2Lib::Area(covered));
        const double overlap_area = std::max(input_area - union_area, 0.0);
        std::vector<OutputPrimitive> replacements;
        bool invalid_simple_union = false;
        // A connected Boolean-union component is one semantic polygon. Its
        // triangulation belongs to the later collision-mesh stage.
        if (!covered.empty())
        {
            const bool has_holes = std::any_of(covered.begin(), covered.end(),
                [](const auto& path) { return !Clipper2Lib::IsPositive(path); });
            if (has_holes)
            {
                // A clean support-contact partition is deliberately split into
                // local polygons so PQSS receives compact leaf triangles.  Its
                // union has holes by construction; globally retriangulating the
                // union would recreate long triangles spanning between distant
                // contacts.  With no measured overlap there is nothing for the
                // canonicalizer to remove, so preserve the local partition.
                if (preserves_cavity_opening &&
                    overlap_area <= tolerance * tolerance)
                {
                    for (const std::size_t index : group)
                    {
                        consumed[index] = true;
                        result.push_back(std::move(primitives[index]));
                    }
                    continue;
                }
                // A polygon-with-hole is represented as a set of non-overlapping
                // triangular Polygon primitives. This preserves the hole and
                // removes overlap without introducing a polygon-with-hole type.
                Clipper2Lib::PathsD triangles;
                if (Clipper2Lib::Triangulate(covered, clipper_precision, triangles, false) !=
                    Clipper2Lib::TriangulateResult::success)
                {
                    for (const std::size_t index : group)
                    {
                        consumed[index] = true;
                        result.push_back(std::move(primitives[index]));
                    }
                    continue;
                }
                for (const auto& triangle : triangles)
                {
                    std::vector<Vec2> polygon;
                    for (const auto& point : triangle)
                        polygon.emplace_back(point.x, point.y);
                    if (polygon.size() != 3) continue;
                    Primitive replacement = polygonPrimitive(polygon, origin, frame);
                    const Vec3 doubled_area =
                        (replacement.polygon[1] - replacement.polygon[0]).cross(
                            replacement.polygon[2] - replacement.polygon[0]);
                    // Clipper triangulation can emit a zero-area numerical
                    // sliver at a coincident boundary event.  It has no covered
                    // support and must not become a semantic primitive.
                    if (doubled_area.norm() <= tolerance * tolerance) continue;
                    replacements.push_back({std::move(replacement), {}});
                }
            }
            else for (const auto& path : covered)
            {
                std::vector<Vec2> polygon;
                polygon.reserve(path.size());
                for (const auto& point : path) polygon.emplace_back(point.x, point.y);
                polygon = simplifyPolygon(std::move(polygon), tolerance);
                if (polygon.size() < 3 ||
                    triangulatePolygon(polygon).size() + 2 != polygon.size())
                {
                    replacements.clear();
                    invalid_simple_union = true;
                    break;
                }
                Primitive replacement = polygonPrimitive(polygon, origin, frame);
                Vec3 doubled_area = Vec3::Zero();
                const Vec3 area_origin = replacement.polygon.front();
                for (std::size_t index = 1; index + 1 < replacement.polygon.size(); ++index)
                    doubled_area +=
                        (replacement.polygon[index] - area_origin).cross(
                            replacement.polygon[index + 1] - area_origin);
                if (doubled_area.norm() <= tolerance * tolerance) continue;
                replacements.push_back({std::move(replacement), {}});
            }
        }
        if (invalid_simple_union)
        {
            // A self-touching Boolean boundary can be a valid planar union even
            // though it is not representable as one simple Polygon primitive.
            // Keeping the overlapping inputs in that case defeats the purpose
            // of canonicalization. Use Clipper's constrained triangulation as
            // the non-overlapping surface partition, exactly as for a union
            // with holes.
            Clipper2Lib::PathsD triangles;
            if (Clipper2Lib::Triangulate(
                    covered, clipper_precision, triangles, false) ==
                Clipper2Lib::TriangulateResult::success)
                for (const auto& triangle : triangles)
                {
                    std::vector<Vec2> polygon;
                    for (const auto& point : triangle)
                        polygon.emplace_back(point.x, point.y);
                    if (polygon.size() != 3) continue;
                    Primitive replacement =
                        polygonPrimitive(polygon, origin, frame);
                    const Vec3 doubled_area =
                        (replacement.polygon[1] -
                         replacement.polygon[0]).cross(
                            replacement.polygon[2] -
                            replacement.polygon[0]);
                    if (doubled_area.norm() <= tolerance * tolerance)
                        continue;
                    replacements.push_back(
                        {std::move(replacement), {}});
                }
        }
        for (auto& replacement : replacements)
        {
            replacement.shallow_shell_coalesced = shallow_shell_coalesced;
            replacement.preserves_cavity_opening = preserves_cavity_opening;
            replacement.enclosure_group = enclosure_consistent
                ? enclosure_group : 0;
        }
        if (replacements.empty() ||
            (overlap_area <= tolerance * tolerance &&
             replacements.size() >= group.size()))
        {
            for (const std::size_t index : group)
            {
                consumed[index] = true;
                result.push_back(std::move(primitives[index]));
            }
            continue;
        }

        std::sort(source_faces.begin(), source_faces.end());
        source_faces.erase(std::unique(source_faces.begin(), source_faces.end()),
                           source_faces.end());
        const std::size_t replacement_begin = result.size();
        result.insert(result.end(), std::make_move_iterator(replacements.begin()),
                      std::make_move_iterator(replacements.end()));

        struct ReplacementProjection
        {
            Bounds2 bounds;
        };
        std::vector<ReplacementProjection> replacement_projections(
            result.size() - replacement_begin);
        for (std::size_t index = replacement_begin; index < result.size(); ++index)
        {
            const PrimitiveMesh replacement_mesh =
                triangulatePrimitive(result[index].primitive);
            Bounds2& bounds = replacement_projections[
                index - replacement_begin].bounds;
            for (const Vec3& vertex : replacement_mesh.vertices)
            {
                const Vec3 local = frame.transposeMultiply(vertex - origin);
                bounds.lower = bounds.lower.cwiseMin({local.x(), local.y()});
                bounds.upper = bounds.upper.cwiseMax({local.x(), local.y()});
            }
        }

        for (const auto face_id : source_faces)
        {
            const Face& face = source_mesh.faces[face_id];
            Bounds2 source_bounds;
            for (const auto vertex : face)
            {
                const Vec3 local = frame.transposeMultiply(
                    source_mesh.vertices[vertex] - origin);
                source_bounds.lower = source_bounds.lower.cwiseMin(
                    {local.x(), local.y()});
                source_bounds.upper = source_bounds.upper.cwiseMax(
                    {local.x(), local.y()});
            }

            // A source triangle can straddle several polygons created by the
            // planar Boolean union.  Assigning it only to the polygon that
            // contains its centroid loses the other pieces; worse, when no
            // polygon contains the (slightly offset) centroid the old code
            // silently assigned it to the first, unrelated replacement.  Give
            // every locally overlapping replacement responsibility instead.
            // The final union audit then proves coverage jointly and no full
            // source triangle has to be appended over an already-covered area.
            for (std::size_t index = replacement_begin; index < result.size(); ++index)
            {
                const Bounds2& candidate = replacement_projections[
                    index - replacement_begin].bounds;
                if (candidate.upper.x() + tolerance < source_bounds.lower.x() ||
                    source_bounds.upper.x() + tolerance < candidate.lower.x() ||
                    candidate.upper.y() + tolerance < source_bounds.lower.y() ||
                    source_bounds.upper.y() + tolerance < candidate.lower.y())
                    continue;
                result[index].source_faces.push_back(face_id);
            }
        }
        for (const auto index : group) consumed[index] = true;
        ++stats.groups;
        // One union component replaces all overlapping input patches even when
        // its non-convex boundary needs more output triangles than the inputs.
        stats.removed_primitives += group.size() > covered.size()
            ? group.size() - covered.size() : 0;
        stats.removed_overlap_area += overlap_area;
    }
    for (std::size_t index = 0; index < primitives.size(); ++index)
        if (!consumed[index]) result.push_back(std::move(primitives[index]));
    return result;
}

bool coveredByPlanarUnion(const std::vector<OutputPrimitive>& primitives,
                          const std::size_t excluded,
                          const Vec3& point,
                          const Vec3& normal,
                          const double tolerance)
{
    for (std::size_t index = 0; index < primitives.size(); ++index)
    {
        if (index == excluded) continue;
        Vec3 candidate_normal;
        if (!planarNormal(primitives[index].primitive, tolerance, candidate_normal)) continue;
        if (std::abs(candidate_normal.dot(normal)) < 1.0 - 1.0e-8) continue;
        if (std::abs((point - planarPoint(primitives[index].primitive)).dot(normal)) > tolerance)
            continue;
        if (containsPoint(primitives[index].primitive, point, tolerance)) return true;
    }
    return false;
}

std::vector<std::uint32_t> removeSealedVoidWalls(
    std::vector<OutputPrimitive>& primitives,
    const double model_diagonal,
    const double tolerance,
    std::size_t& removed_count)
{
    std::vector<bool> removed(primitives.size(), false);
    std::vector<std::uint32_t> excluded_faces;
    std::vector<Vec3> plane_normals(primitives.size(), Vec3::Zero());
    std::vector<Vec3> plane_points(primitives.size(), Vec3::Zero());
    std::vector<PrimitiveMesh> triangulated;
    triangulated.reserve(primitives.size());
    for (std::size_t index = 0; index < primitives.size(); ++index)
    {
        (void)planarNormal(primitives[index].primitive, tolerance,
                           plane_normals[index]);
        plane_points[index] = planarPoint(primitives[index].primitive);
        triangulated.push_back(triangulatePrimitive(primitives[index].primitive));
    }
    const auto covered = [&](const std::size_t excluded,
                             const Vec3& point,
                             const Vec3& normal)
    {
        for (std::size_t index = 0; index < primitives.size(); ++index)
        {
            if (index == excluded || plane_normals[index].norm() == 0.0 ||
                std::abs(plane_normals[index].dot(normal)) < 1.0 - 1.0e-8 ||
                std::abs((point - plane_points[index]).dot(normal)) > tolerance * 8.0)
                continue;
            if (containsPointCached(primitives[index].primitive,
                                    triangulated[index], point,
                                    tolerance * 8.0))
                return true;
        }
        return false;
    };
    for (std::size_t index = 0; index < primitives.size(); ++index)
    {
        const Primitive& wall = primitives[index].primitive;
        if (wall.kind != Kind::Rectangle) continue;
        const std::array<Vec3, 2> axes{wall.axes.col(0), wall.axes.col(1)};
        const std::array<double, 2> half_sizes{wall.half_size.x(), wall.half_size.y()};
        for (int span_axis = 0; span_axis < 2 && !removed[index]; ++span_axis)
        {
            const int line_axis = 1 - span_axis;
            if (half_sizes[span_axis] <= tolerance || half_sizes[line_axis] <= tolerance)
                continue;
            // A void wall bridges a locally thin depth interval.  Without this
            // guard a large cap can be interpreted sideways and delete itself
            // merely because other exterior faces exist at its distant edges.
            if (2.0 * half_sizes[span_axis] > 0.25 * model_diagonal) continue;
            const Vec3 cap_normal = axes[span_axis].normalized();
            const Vec3 line_direction = axes[line_axis].normalized();
            Vec3 side_direction = cap_normal.cross(line_direction);
            if (side_direction.norm() <= tolerance) continue;
            side_direction = side_direction.normalized();
            const double offset = std::max(
                tolerance * 16.0,
                std::min(half_sizes[line_axis] * 0.05, model_diagonal * 1.0e-4));
            bool closed = true;
            for (const double end_sign : {-1.0, 1.0})
            {
                const Vec3 end_center = wall.center + cap_normal *
                    (end_sign * half_sizes[span_axis]);
                for (const double along : {-0.75, 0.0, 0.75})
                {
                    const Vec3 line_point = end_center + line_direction *
                        (along * half_sizes[line_axis]);
                    for (const double side : {-1.0, 1.0})
                    {
                        const Vec3 sample = line_point + side_direction * (side * offset);
                        if (!covered(index, sample, cap_normal))
                        {
                            closed = false;
                            break;
                        }
                    }
                    if (!closed) break;
                }
                if (!closed) break;
            }
            if (!closed) continue;
            removed[index] = true;
            excluded_faces.insert(excluded_faces.end(),
                                  primitives[index].source_faces.begin(),
                                  primitives[index].source_faces.end());
        }
    }

    std::vector<OutputPrimitive> kept;
    kept.reserve(primitives.size());
    for (std::size_t index = 0; index < primitives.size(); ++index)
        if (!removed[index]) kept.push_back(std::move(primitives[index]));
    removed_count = primitives.size() - kept.size();
    primitives = std::move(kept);
    std::sort(excluded_faces.begin(), excluded_faces.end());
    excluded_faces.erase(std::unique(excluded_faces.begin(), excluded_faces.end()),
                         excluded_faces.end());
    return excluded_faces;
}

std::vector<std::uint32_t> removeBlindCavitySurfaces(
    std::vector<OutputPrimitive>& primitives,
    const Vec3& model_center,
    const double model_diagonal,
    const double tolerance,
    std::size_t& removed_count)
{
    std::vector<bool> removed(primitives.size(), false);
    std::vector<std::uint32_t> excluded_faces;
    std::vector<PrimitiveMesh> triangulated;
    std::vector<double> areas;
    std::vector<std::size_t> cap_indices;
    triangulated.reserve(primitives.size());
    areas.reserve(primitives.size());
    for (std::size_t index = 0; index < primitives.size(); ++index)
    {
        triangulated.push_back(triangulatePrimitive(primitives[index].primitive));
        double area = 0.0;
        for (const Face& face : triangulated.back().faces)
            area += 0.5 *
                (triangulated.back().vertices[face[1]] -
                 triangulated.back().vertices[face[0]]).cross(
                    triangulated.back().vertices[face[2]] -
                    triangulated.back().vertices[face[0]]).norm();
        areas.push_back(area);
        if (primitives[index].primitive.kind == Kind::Rectangle)
            cap_indices.push_back(index);
    }
    for (std::size_t candidate = 0; candidate < primitives.size(); ++candidate)
    {
        const PrimitiveMesh& candidate_mesh = triangulated[candidate];
        if (candidate_mesh.vertices.empty()) continue;
        for (const std::size_t cap_index : cap_indices)
        {
            if (candidate == cap_index) continue;
            const Primitive& cap = primitives[cap_index].primitive;
            const double cap_area = 4.0 * cap.half_size.x() * cap.half_size.y();
            if (cap_area <= 4.0 * areas[candidate]) continue;
            Vec3 inward_normal = cap.axes.col(2);
            if ((model_center - cap.center).dot(inward_normal) < 0.0)
                inward_normal *= -1.0;
            double minimum_depth = std::numeric_limits<double>::infinity();
            double maximum_depth = -std::numeric_limits<double>::infinity();
            bool strictly_inside = true;
            const double margin = std::max(
                tolerance * 16.0,
                std::min({cap.half_size.x(), cap.half_size.y(), model_diagonal * 1.0e-3}));
            for (const Vec3& vertex : candidate_mesh.vertices)
            {
                const Vec3 local = cap.axes.transposeMultiply(vertex - cap.center);
                const double inward_depth = (vertex - cap.center).dot(inward_normal);
                minimum_depth = std::min(minimum_depth, inward_depth);
                maximum_depth = std::max(maximum_depth, inward_depth);
                strictly_inside &= std::abs(local.x()) <= cap.half_size.x() - margin &&
                                   std::abs(local.y()) <= cap.half_size.y() - margin;
            }
            if (!strictly_inside) continue;
            // A cap can only hide a cavity on its inward side. Projected
            // containment on the outward side describes a protrusion (for
            // example a foot below a base plate), which must remain part of the
            // conservative responsibility set until an enclosing 3D proxy
            // replaces it.
            if (minimum_depth < -tolerance ||
                maximum_depth > 0.03 * model_diagonal ||
                maximum_depth <= tolerance) continue;
            // A fragmented opposite skin of a thin shell is not a blind
            // cavity.  Assess its complete coplanar family before deleting an
            // individual fragment; otherwise a filled support plane can erase
            // the other side of an ordinary plate one patch at a time.
            Vec3 candidate_normal;
            if (planarNormal(primitives[candidate].primitive, tolerance,
                             candidate_normal) &&
                std::abs(candidate_normal.dot(cap.axes.col(2))) >= 1.0 - 1.0e-8)
            {
                const double candidate_coordinate =
                    planarPoint(primitives[candidate].primitive).dot(cap.axes.col(2));
                double coplanar_supported_area = 0.0;
                for (std::size_t peer = 0; peer < primitives.size(); ++peer)
                {
                    Vec3 peer_normal;
                    if (!planarNormal(primitives[peer].primitive, tolerance,
                                      peer_normal) ||
                        std::abs(peer_normal.dot(cap.axes.col(2))) < 1.0 - 1.0e-8 ||
                        std::abs(planarPoint(primitives[peer].primitive).dot(
                                     cap.axes.col(2)) - candidate_coordinate) >
                            tolerance * 16.0)
                        continue;
                    bool within_cap_projection = true;
                    for (const Vec3& vertex : triangulated[peer].vertices)
                    {
                        const Vec3 local = cap.axes.transposeMultiply(vertex - cap.center);
                        if (std::abs(local.x()) > cap.half_size.x() + margin ||
                            std::abs(local.y()) > cap.half_size.y() + margin)
                        {
                            within_cap_projection = false;
                            break;
                        }
                    }
                    if (within_cap_projection)
                        coplanar_supported_area += areas[peer];
                }
                if (coplanar_supported_area >= 0.25 * cap_area) continue;
            }
            removed[candidate] = true;
            excluded_faces.insert(excluded_faces.end(),
                                  primitives[candidate].source_faces.begin(),
                                  primitives[candidate].source_faces.end());
            break;
        }
    }
    std::vector<OutputPrimitive> kept;
    kept.reserve(primitives.size());
    for (std::size_t index = 0; index < primitives.size(); ++index)
        if (!removed[index]) kept.push_back(std::move(primitives[index]));
    removed_count = primitives.size() - kept.size();
    primitives = std::move(kept);
    std::sort(excluded_faces.begin(), excluded_faces.end());
    excluded_faces.erase(std::unique(excluded_faces.begin(), excluded_faces.end()),
                         excluded_faces.end());
    return excluded_faces;
}

bool axisAlignedRectangle(const Primitive& primitive, const double tolerance,
                          int& normal_axis, Bounds& bounds)
{
    if (primitive.kind != Kind::Rectangle) return false;
    bounds = primitiveBounds(primitive);
    const Vec3 extent = bounds.upper - bounds.lower;
    normal_axis = 0;
    if (extent.y() < extent[normal_axis]) normal_axis = 1;
    if (extent.z() < extent[normal_axis]) normal_axis = 2;
    if (extent[normal_axis] > tolerance) return false;
    return true;
}

void recognizeClosedAxisAlignedBoxes(std::vector<OutputPrimitive>& primitives,
                                     const double tolerance,
                                     std::size_t& recognized_count,
                                     std::vector<RecognizedProtrusion>* recognized_boxes = nullptr)
{
    std::vector<bool> consumed(primitives.size(), false);
    for (std::size_t first = 0; first < primitives.size(); ++first)
    {
        if (consumed[first]) continue;
        int first_axis = 0;
        Bounds first_bounds;
        if (!axisAlignedRectangle(primitives[first].primitive, tolerance,
                                  first_axis, first_bounds)) continue;
        for (std::size_t opposite = first + 1; opposite < primitives.size(); ++opposite)
        {
            if (consumed[opposite]) continue;
            int opposite_axis = 0;
            Bounds opposite_bounds;
            if (!axisAlignedRectangle(primitives[opposite].primitive, tolerance,
                                      opposite_axis, opposite_bounds) ||
                opposite_axis != first_axis) continue;
            bool same_projection = true;
            for (int axis = 0; axis < 3; ++axis)
            {
                if (axis == first_axis) continue;
                same_projection &= std::abs(first_bounds.lower[axis] -
                                            opposite_bounds.lower[axis]) <= tolerance &&
                                   std::abs(first_bounds.upper[axis] -
                                            opposite_bounds.upper[axis]) <= tolerance;
            }
            if (!same_projection || std::abs(first_bounds.lower[first_axis] -
                                              opposite_bounds.lower[first_axis]) <= tolerance)
                continue;
            Bounds box_bounds;
            box_bounds.lower = first_bounds.lower.cwiseMin(opposite_bounds.lower);
            box_bounds.upper = first_bounds.upper.cwiseMax(opposite_bounds.upper);
            std::array<std::size_t, 6> faces;
            faces.fill(std::numeric_limits<std::size_t>::max());
            for (std::size_t index = 0; index < primitives.size(); ++index)
            {
                if (consumed[index]) continue;
                int axis = 0;
                Bounds bounds;
                if (!axisAlignedRectangle(primitives[index].primitive, tolerance,
                                          axis, bounds)) continue;
                bool matches = true;
                for (int tangent = 0; tangent < 3; ++tangent)
                {
                    if (tangent == axis) continue;
                    matches &= std::abs(bounds.lower[tangent] - box_bounds.lower[tangent]) <= tolerance &&
                               std::abs(bounds.upper[tangent] - box_bounds.upper[tangent]) <= tolerance;
                }
                if (!matches) continue;
                const double coordinate = bounds.lower[axis];
                if (std::abs(coordinate - box_bounds.lower[axis]) <= tolerance)
                    faces[2 * axis] = index;
                else if (std::abs(coordinate - box_bounds.upper[axis]) <= tolerance)
                    faces[2 * axis + 1] = index;
            }
            if (std::any_of(faces.begin(), faces.end(), [](const auto index)
                { return index == std::numeric_limits<std::size_t>::max(); })) continue;
            BoxFit box;
            box.center = (box_bounds.lower + box_bounds.upper) * 0.5;
            box.axes = Mat3::Identity();
            box.half_size = (box_bounds.upper - box_bounds.lower) * 0.5;
            box.volume = 8.0 * box.half_size.prod();
            std::vector<std::uint32_t> source_faces;
            for (const auto index : faces)
            {
                consumed[index] = true;
                source_faces.insert(source_faces.end(), primitives[index].source_faces.begin(),
                                    primitives[index].source_faces.end());
            }
            if (recognized_boxes != nullptr)
                recognized_boxes->push_back({box, std::move(source_faces), {}, -1, 0.0});
            ++recognized_count;
            break;
        }
    }
}

void appendBoxRectangles(std::vector<OutputPrimitive>& output,
                         const BoxFit& box,
                         const std::vector<std::uint32_t>& source_faces,
                         const int covered_face_axis = -1,
                         const double covered_face_sign = 0.0,
                         const std::uint64_t enclosure_group = 0)
{
    for (int normal_axis = 0; normal_axis < 3; ++normal_axis)
    {
        const int first_axis = (normal_axis + 1) % 3;
        const int second_axis = (normal_axis + 2) % 3;
        for (const double sign : {-1.0, 1.0})
        {
            if (normal_axis == covered_face_axis && sign == covered_face_sign) continue;
            Primitive rectangle;
            rectangle.kind = Kind::Rectangle;
            rectangle.center = box.center + box.axes.col(normal_axis) *
                (sign * box.half_size[normal_axis]);
            rectangle.axes.col(0) = box.axes.col(first_axis);
            rectangle.axes.col(1) = box.axes.col(second_axis);
            rectangle.axes.col(2) = box.axes.col(normal_axis) * sign;
            rectangle.half_size = {
                box.half_size[first_axis], box.half_size[second_axis], 0.0};
            OutputPrimitive item;
            item.primitive = std::move(rectangle);
            item.source_faces = source_faces;
            item.enclosure_group = enclosure_group;
            output.push_back(std::move(item));
        }
    }
}

struct SupportClippedBoxStats
{
    std::size_t clipped_faces = 0;
    std::size_t discarded_faces = 0;
    double maximum_removed_penetration = 0.0;
};

void appendSupportClippedBoxSurfaces(
    std::vector<OutputPrimitive>& output,
    const RecognizedProtrusion& candidate,
    const double tolerance,
    SupportClippedBoxStats& stats)
{
    std::vector<OutputPrimitive> box_surfaces;
    appendBoxRectangles(
        box_surfaces, candidate.box, candidate.faces,
        candidate.covered_face_axis, candidate.covered_face_sign, 0);
    if (!candidate.has_support_plane)
    {
        output.insert(output.end(),
                      std::make_move_iterator(box_surfaces.begin()),
                      std::make_move_iterator(box_surfaces.end()));
        return;
    }

    const Vec3 outward = candidate.support_outward.normalized();
    const auto signedDistance = [&](const Vec3& point)
    {
        return (point - candidate.support_plane_point).dot(outward);
    };
    const auto clipToOutsideHalfspace = [&](const std::vector<Vec3>& boundary)
    {
        std::vector<Vec3> clipped;
        if (boundary.empty()) return clipped;
        for (std::size_t index = 0; index < boundary.size(); ++index)
        {
            const Vec3& first = boundary[index];
            const Vec3& second = boundary[(index + 1) % boundary.size()];
            const double first_distance = signedDistance(first);
            const double second_distance = signedDistance(second);
            const bool first_inside = first_distance >= -tolerance;
            const bool second_inside = second_distance >= -tolerance;
            stats.maximum_removed_penetration = std::max(
                stats.maximum_removed_penetration,
                std::max(-first_distance, 0.0));
            if (first_inside) clipped.push_back(first);
            if (first_inside == second_inside) continue;
            const double denominator = first_distance - second_distance;
            if (std::abs(denominator) <= 1.0e-30) continue;
            const double parameter = std::clamp(
                first_distance / denominator, 0.0, 1.0);
            clipped.push_back(first + (second - first) * parameter);
        }
        std::vector<Vec3> cleaned;
        cleaned.reserve(clipped.size());
        for (const Vec3& point : clipped)
            if (cleaned.empty() ||
                (point - cleaned.back()).norm() > tolerance)
                cleaned.push_back(point);
        if (cleaned.size() > 1 &&
            (cleaned.front() - cleaned.back()).norm() <= tolerance)
            cleaned.pop_back();
        bool changed = true;
        while (changed && cleaned.size() >= 3)
        {
            changed = false;
            for (std::size_t index = 0; index < cleaned.size(); ++index)
            {
                const Vec3 first = cleaned[index] -
                    cleaned[(index + cleaned.size() - 1) % cleaned.size()];
                const Vec3 second = cleaned[(index + 1) % cleaned.size()] -
                    cleaned[index];
                const double scale = std::max(first.norm() * second.norm(),
                                              1.0e-30);
                if (first.cross(second).norm() > tolerance * scale)
                    continue;
                cleaned.erase(cleaned.begin() + index);
                changed = true;
                break;
            }
        }
        return cleaned;
    };

    for (OutputPrimitive& surface : box_surfaces)
    {
        const PrimitiveMesh rectangle = triangulatePrimitive(surface.primitive);
        const std::vector<Vec3> boundary{
            rectangle.vertices[0], rectangle.vertices[2],
            rectangle.vertices[3], rectangle.vertices[1]};
        const bool penetrates = std::any_of(
            boundary.begin(), boundary.end(), [&](const Vec3& point)
            { return signedDistance(point) < -tolerance; });
        if (!penetrates)
        {
            output.push_back(std::move(surface));
            continue;
        }
        std::vector<Vec3> clipped = clipToOutsideHalfspace(boundary);
        if (clipped.size() < 3)
        {
            ++stats.discarded_faces;
            continue;
        }
        surface.primitive = Primitive{};
        surface.primitive.kind = Kind::Polygon;
        surface.primitive.polygon = std::move(clipped);
        output.push_back(std::move(surface));
        ++stats.clipped_faces;
    }
}

std::vector<Vec2> boxPlaneCrossSection(const BoxFit& box,
                                       const Vec3& origin,
                                       const Mat3& frame,
                                       double tolerance,
                                       bool& occupies_negative_side,
                                       bool& occupies_positive_side);

struct SupportContactClipStats
{
    std::size_t clipped_primitives = 0;
    std::size_t removed_primitives = 0;
    std::size_t output_fragments = 0;
    double removed_area = 0.0;
};

struct ParallelOcclusionStats
{
    std::size_t candidate_count = 0;
    std::size_t accepted_candidates = 0;
    std::size_t clipped_primitives = 0;
    std::size_t removed_primitives = 0;
    std::size_t output_fragments = 0;
    double removed_area = 0.0;
    double maximum_covered_ratio = 0.0;
    std::vector<std::pair<std::size_t, double>> candidate_ratios;
    std::vector<std::size_t> accepted_ids;
};

// A face removed as an internal parallel layer is also an occlusion witness
// for the portions of orthogonal neighbour faces that continue through that
// layer.  Treating the witness as a 2-D coverer only removes the layer itself;
// the neighbour then leaves a thin, fully internal skirt in the proxy.  Clip
// that skirt against the witness plane whenever the two polygons actually
// meet.  The decision is purely geometric (plane intersection + footprint
// containment) and therefore applies equally to any combination of panels.
std::vector<OutputPrimitive> clipAdjacentFacesAtOcclusionPlanes(
    std::vector<OutputPrimitive> primitives,
    const std::vector<OutputPrimitive>& occlusion_certificates,
    const Vec3& model_center,
    const double tolerance,
    ParallelOcclusionStats& stats)
{
    const auto polygonArea = [](const std::vector<Vec3>& polygon)
    {
        if (polygon.size() < 3) return 0.0;
        Vec3 vector_area = Vec3::Zero();
        for (std::size_t index = 0; index < polygon.size(); ++index)
            vector_area += polygon[index].cross(
                polygon[(index + 1) % polygon.size()]);
        return 0.5 * vector_area.norm();
    };
    const auto convexProjected = [](const std::vector<Vec2>& polygon)
    {
        int sign = 0;
        for (std::size_t index = 0; index < polygon.size(); ++index)
        {
            const Vec2 first = polygon[(index + polygon.size() - 1) % polygon.size()];
            const Vec2 second = polygon[index];
            const Vec2 third = polygon[(index + 1) % polygon.size()];
            const Vec2 first_edge = second - first;
            const Vec2 second_edge = third - second;
            const double cross = first_edge.x() * second_edge.y() -
                first_edge.y() * second_edge.x();
            if (std::abs(cross) <= 1.0e-12) continue;
            const int current_sign = cross > 0.0 ? 1 : -1;
            if (sign == 0) sign = current_sign;
            else if (sign != current_sign) return false;
        }
        return sign != 0;
    };
    auto clipToHalfspace = [&](const std::vector<Vec3>& polygon,
                               const Vec3& normal, const double offset,
                               const double keep_sign)
    {
        std::vector<Vec3> result;
        if (polygon.empty()) return result;
        const auto signed_value = [&](const Vec3& point)
        { return keep_sign * (normal.dot(point) - offset); };
        for (std::size_t index = 0; index < polygon.size(); ++index)
        {
            const Vec3& first = polygon[index];
            const Vec3& second = polygon[(index + 1) % polygon.size()];
            const double first_value = signed_value(first);
            const double second_value = signed_value(second);
            const bool first_inside = first_value >= -tolerance;
            const bool second_inside = second_value >= -tolerance;
            if (first_inside) result.push_back(first);
            if (first_inside != second_inside)
            {
                const double denominator = first_value - second_value;
                if (std::abs(denominator) > 1.0e-30)
                {
                    const double weight = first_value / denominator;
                    result.push_back(first + (second - first) * weight);
                }
            }
        }
        std::vector<Vec3> cleaned;
        for (const Vec3& point : result)
            if (cleaned.empty() || (point - cleaned.back()).norm() > tolerance)
                cleaned.push_back(point);
        if (cleaned.size() > 1 &&
            (cleaned.front() - cleaned.back()).norm() <= tolerance)
            cleaned.pop_back();
        bool removed_collinear = true;
        while (removed_collinear && cleaned.size() >= 3)
        {
            removed_collinear = false;
            for (std::size_t index = 0; index < cleaned.size(); ++index)
            {
                const Vec3 previous = cleaned[
                    (index + cleaned.size() - 1) % cleaned.size()];
                const Vec3 current = cleaned[index];
                const Vec3 next = cleaned[(index + 1) % cleaned.size()];
                const Vec3 first = current - previous;
                const Vec3 second = next - current;
                if (first.norm() <= tolerance || second.norm() <= tolerance ||
                    first.cross(second).norm() <= tolerance *
                        std::max(first.norm() * second.norm(), 1.0))
                {
                    cleaned.erase(cleaned.begin() + index);
                    removed_collinear = true;
                    break;
                }
            }
        }
        return cleaned;
    };

    auto planeIntersectionPoints = [&](const std::vector<Vec3>& polygon,
                                       const Vec3& normal, const double offset)
    {
        std::vector<Vec3> points;
        for (std::size_t index = 0; index < polygon.size(); ++index)
        {
            const Vec3& first = polygon[index];
            const Vec3& second = polygon[(index + 1) % polygon.size()];
            const double first_value = normal.dot(first) - offset;
            const double second_value = normal.dot(second) - offset;
            if (std::abs(first_value) <= tolerance)
                points.push_back(first);
            if ((first_value < -tolerance && second_value > tolerance) ||
                (first_value > tolerance && second_value < -tolerance))
            {
                const double weight = first_value /
                    (first_value - second_value);
                points.push_back(first + (second - first) * weight);
            }
        }
        std::vector<Vec3> unique;
        for (const Vec3& point : points)
            if (std::none_of(unique.begin(), unique.end(),
                [&](const Vec3& other)
                { return (point - other).norm() <= tolerance * 8.0; }))
                unique.push_back(point);
        return unique;
    };

    std::vector<OutputPrimitive> result;
    result.reserve(primitives.size());
    for (OutputPrimitive& item : primitives)
    {
        if (item.primitive.kind != Kind::Polygon ||
            item.primitive.polygon.size() < 3)
        {
            result.push_back(std::move(item));
            continue;
        }
        std::vector<Vec3> current = item.primitive.polygon;
        bool changed = false;
        for (const OutputPrimitive& certificate : occlusion_certificates)
        {
            if (certificate.primitive.kind != Kind::Polygon ||
                certificate.primitive.polygon.size() < 3)
                continue;
            Vec3 certificate_normal;
            if (!planarNormal(certificate.primitive, tolerance,
                              certificate_normal))
                continue;
            const Vec3 certificate_point =
                certificate.primitive.polygon.front();
            const double offset = certificate_normal.dot(certificate_point);
            std::vector<Vec3> intersection = planeIntersectionPoints(
                current, certificate_normal, offset);
            if (intersection.size() < 2)
            {
                // A few CAD polygons are exported with a seam that skips the
                // geometric boundary edge. Recover the plane cut from all
                // opposite-side vertex pairs instead of trusting that seam's
                // winding order.
                for (std::size_t first = 0; first < current.size(); ++first)
                    for (std::size_t second = first + 1;
                         second < current.size(); ++second)
                    {
                        const double a = certificate_normal.dot(current[first]) - offset;
                        const double b = certificate_normal.dot(current[second]) - offset;
                        if ((a < -tolerance && b > tolerance) ||
                            (a > tolerance && b < -tolerance))
                        {
                            const double weight = a / (a - b);
                            intersection.push_back(current[first] +
                                (current[second] - current[first]) * weight);
                        }
                    }
                if (intersection.size() > 2)
                {
                    std::vector<Vec3> unique;
                    for (const Vec3& point : intersection)
                        if (std::none_of(unique.begin(), unique.end(),
                            [&](const Vec3& other)
                            { return (point - other).norm() <= tolerance * 8.0; }))
                            unique.push_back(point);
                    intersection = std::move(unique);
                }
            }
            if (intersection.size() < 2)
            {
                // Fall back to the actual triangulated surface.  This handles
                // polygon seams whose winding order does not describe the
                // geometric boundary, while keeping every retained point on
                // the requested half-space.
                const PrimitiveMesh triangle_mesh = triangulatePrimitive(
                    item.primitive);
                std::vector<Vec3> kept_points;
                for (const Face& face : triangle_mesh.faces)
                {
                    std::vector<Vec3> triangle{
                        triangle_mesh.vertices[face[0]],
                        triangle_mesh.vertices[face[1]],
                        triangle_mesh.vertices[face[2]]};
                    const double first_side = certificate_normal.dot(triangle[0]) - offset;
                    const double second_side = certificate_normal.dot(triangle[1]) - offset;
                    const double third_side = certificate_normal.dot(triangle[2]) - offset;
                    if (!((first_side < -tolerance &&
                           (second_side > tolerance || third_side > tolerance)) ||
                          (second_side < -tolerance &&
                           (first_side > tolerance || third_side > tolerance)) ||
                          (third_side < -tolerance &&
                           (first_side > tolerance || second_side > tolerance))))
                        continue;
                    const auto clipped_triangle = clipToHalfspace(
                        triangle, certificate_normal, offset,
                        certificate_normal.dot(model_center - certificate_point) >
                            0.0 ? 1.0 : -1.0);
                    kept_points.insert(kept_points.end(),
                        clipped_triangle.begin(), clipped_triangle.end());
                }
                if (kept_points.size() >= 3)
                {
                    Vec3 local_normal;
                    if (!planarNormal(item.primitive, tolerance, local_normal))
                        continue;
                    const Mat3 basis = orthonormalFrame(local_normal);
                    std::vector<Vec2> projected_points;
                    projected_points.reserve(kept_points.size());
                    for (const Vec3& point : kept_points)
                    {
                        const Vec3 local = basis.transposeMultiply(
                            point - item.primitive.polygon.front());
                        projected_points.emplace_back(local.x(), local.y());
                    }
                    const auto hull = convexHull(projected_points, tolerance);
                    if (hull.size() >= 3)
                    {
                        current.clear();
                        current.reserve(hull.size());
                        for (const Vec2& point : hull)
                            current.push_back(item.primitive.polygon.front() +
                                basis * Vec3(point.x(), point.y(), 0.0));
                        changed = true;
                        ++stats.clipped_primitives;
                        continue;
                    }
                }
                continue;
            }
            const Bounds certificate_bounds = primitiveBounds(
                certificate.primitive);
            const bool footprint_match = std::any_of(
                intersection.begin(), intersection.end(), [&](const Vec3& point)
                {
                    return point.x() >= certificate_bounds.lower.x() - tolerance * 16.0 &&
                           point.x() <= certificate_bounds.upper.x() + tolerance * 16.0 &&
                           point.y() >= certificate_bounds.lower.y() - tolerance * 16.0 &&
                           point.y() <= certificate_bounds.upper.y() + tolerance * 16.0 &&
                           point.z() >= certificate_bounds.lower.z() - tolerance * 16.0 &&
                           point.z() <= certificate_bounds.upper.z() + tolerance * 16.0;
                });
            if (!footprint_match && !certificate.preserves_cavity_opening)
                continue;

            const double orientation = certificate_normal.dot(
                model_center - certificate_point);
            if (std::abs(orientation) <= tolerance) continue;
            const double keep_sign = orientation > 0.0 ? 1.0 : -1.0;
            std::vector<Vec3> clipped = clipToHalfspace(
                current, certificate_normal, offset, keep_sign);
            if (clipped.size() < 3) { current.clear(); changed = true; break; }
            const double before_area = polygonArea(current);
            const double after_area = polygonArea(clipped);
            double minimum_side = std::numeric_limits<double>::infinity();
            double maximum_side = -std::numeric_limits<double>::infinity();
            for (const Vec3& point : current)
            {
                const double side = keep_sign * (certificate_normal.dot(point) - offset);
                minimum_side = std::min(minimum_side, side);
                maximum_side = std::max(maximum_side, side);
            }
            const bool crosses_plane = minimum_side < -tolerance &&
                                       maximum_side > tolerance;
            if (after_area + tolerance * tolerance < before_area || crosses_plane)
            {
                Vec3 clipped_normal = Vec3::Zero();
                for (std::size_t index = 1; index + 1 < clipped.size(); ++index)
                {
                    clipped_normal = (clipped[index] - clipped.front()).cross(
                        clipped[index + 1] - clipped.front());
                    if (clipped_normal.norm() > tolerance * tolerance)
                    {
                        clipped_normal /= clipped_normal.norm();
                        break;
                    }
                }
                if (clipped_normal.norm() <= 0.0) continue;
                const Mat3 clipped_basis = orthonormalFrame(clipped_normal);
                std::vector<Vec2> clipped_projected;
                clipped_projected.reserve(clipped.size());
                bool coplanar = true;
                for (const Vec3& point : clipped)
                {
                    if (std::abs(clipped_normal.dot(point - clipped.front())) >
                        tolerance * 16.0)
                    {
                        coplanar = false;
                        break;
                    }
                    const Vec3 local = clipped_basis.transposeMultiply(
                        point - clipped.front());
                    clipped_projected.emplace_back(local.x(), local.y());
                }
                if (coplanar && signedArea(clipped_projected) < 0.0)
                    std::reverse(clipped_projected.begin(),
                                clipped_projected.end());
                const auto clipped_triangles = triangulatePolygon(clipped_projected);
                if (!coplanar || (clipped_triangles.size() + 2 !=
                    clipped_projected.size() &&
                    !convexProjected(clipped_projected)))
                {
                    // Some recognised CAD polygons contain a collinear seam
                    // that is harmless before clipping but makes their vertex
                    // order non-simple after a plane cut.  Reconstruct only
                    // this clipped fragment from its own boundary hull; this
                    // never broadens the kept half-space.
                    const auto hull = convexHull(clipped_projected, tolerance);
                    if (hull.size() < 3) continue;
                    const Vec3 clip_origin = current.front();
                    clipped.clear();
                    clipped.reserve(hull.size());
                    for (const Vec2& point : hull)
                        clipped.push_back(clip_origin + clipped_basis *
                            Vec3(point.x(), point.y(), 0.0));
                }
                current = clipped;
                changed = true;
                ++stats.clipped_primitives;
                stats.removed_area += std::max(before_area - after_area, 0.0);
            }
        }
        if (current.size() < 3 || polygonArea(current) <= tolerance * tolerance)
        {
            ++stats.removed_primitives;
            continue;
        }
        Vec3 candidate_normal = Vec3::Zero();
        for (std::size_t index = 1; index + 1 < current.size(); ++index)
        {
            candidate_normal = (current[index] - current.front()).cross(
                current[index + 1] - current.front());
            if (candidate_normal.norm() > tolerance * tolerance)
            {
                candidate_normal /= candidate_normal.norm();
                break;
            }
        }
        std::vector<Vec2> projected;
        if (candidate_normal.norm() <= 0.0)
        {
            result.push_back(std::move(item));
            continue;
        }
        const Mat3 candidate_basis = orthonormalFrame(candidate_normal);
        projected.reserve(current.size());
        for (const Vec3& point : current)
        {
            if (std::abs(candidate_normal.dot(point - current.front())) >
                tolerance * 16.0)
            {
                projected.clear();
                break;
            }
            const Vec3 local = candidate_basis.transposeMultiply(
                point - current.front());
            projected.emplace_back(local.x(), local.y());
        }
        if (!projected.empty() && signedArea(projected) < 0.0)
            std::reverse(projected.begin(), projected.end());
        if (projected.size() != current.size() ||
            (triangulatePolygon(projected).size() + 2 != projected.size() &&
             !convexProjected(projected)))
        {
            // The clipping planes can meet at a polygon vertex.  If numerical
            // noise makes the resulting ordering non-simple, retain the
            // original face rather than exporting an invalid/self-crossing
            // polygon.
            result.push_back(std::move(item));
            continue;
        }
        if (!changed)
        {
            result.push_back(std::move(item));
            continue;
        }
        OutputPrimitive clipped = item;
        clipped.primitive.polygon = std::move(current);
        clipped.preserves_cavity_opening = true;
        result.push_back(std::move(clipped));
        ++stats.output_fragments;
    }
    return result;
}

std::vector<OutputPrimitive> clipAdjacentFaceTrianglesAtOcclusionPlanes(
    std::vector<OutputPrimitive> primitives,
    const std::vector<OutputPrimitive>& occlusion_certificates,
    const Vec3& model_center,
    const double tolerance,
    ParallelOcclusionStats& stats)
{
    constexpr int clipper_precision = 8;
    struct Plane
    {
        const Primitive* surface = nullptr;
        Vec3 normal = Vec3::Zero();
        double offset = 0.0;
        double keep_sign = 1.0;
        Bounds bounds;
        bool seam_tolerant = false;
    };
    std::vector<Plane> planes;
    for (const auto& certificate : occlusion_certificates)
    {
        if (certificate.primitive.kind != Kind::Polygon ||
            certificate.primitive.polygon.size() < 3)
            continue;
        Vec3 normal;
        if (!planarNormal(certificate.primitive, tolerance, normal)) continue;
        const Vec3 point = certificate.primitive.polygon.front();
        const double orientation = normal.dot(model_center - point);
        if (std::abs(orientation) <= tolerance) continue;
        planes.push_back({
            &certificate.primitive, normal, normal.dot(point),
            orientation > 0.0 ? 1.0 : -1.0,
            primitiveBounds(certificate.primitive),
            certificate.preserves_cavity_opening});
    }
    if (planes.empty()) return primitives;

    const auto clippedHalfspace = [&](const std::vector<Vec3>& polygon,
                                      const Plane& plane)
    {
        std::vector<Vec3> clipped;
        for (std::size_t index = 0; index < polygon.size(); ++index)
        {
            const Vec3& first = polygon[index];
            const Vec3& second = polygon[(index + 1) % polygon.size()];
            const double first_value = plane.keep_sign *
                (plane.normal.dot(first) - plane.offset);
            const double second_value = plane.keep_sign *
                (plane.normal.dot(second) - plane.offset);
            const bool first_inside = first_value >= -tolerance;
            const bool second_inside = second_value >= -tolerance;
            if (first_inside) clipped.push_back(first);
            if (first_inside != second_inside)
            {
                const double denominator = first_value - second_value;
                if (std::abs(denominator) > 1.0e-30)
                    clipped.push_back(first + (second - first) *
                        (first_value / denominator));
            }
        }
        std::vector<Vec3> unique;
        for (const Vec3& point : clipped)
            if (unique.empty() || (point - unique.back()).norm() > tolerance)
                unique.push_back(point);
        if (unique.size() > 1 &&
            (unique.front() - unique.back()).norm() <= tolerance)
            unique.pop_back();
        return unique;
    };
    const auto candidateMeetsFootprint = [&](const Primitive& candidate,
                                             const Plane& plane)
    {
        if (plane.seam_tolerant) return true;
        std::vector<Vec3> points;
        for (std::size_t index = 0; index < candidate.polygon.size(); ++index)
        {
            const Vec3& first = candidate.polygon[index];
            const Vec3& second = candidate.polygon[
                (index + 1) % candidate.polygon.size()];
            const double a = plane.normal.dot(first) - plane.offset;
            const double b = plane.normal.dot(second) - plane.offset;
            if (std::abs(a) <= tolerance) points.push_back(first);
            if ((a < -tolerance && b > tolerance) ||
                (a > tolerance && b < -tolerance))
                points.push_back(first + (second - first) * (a / (a - b)));
        }
        return std::any_of(points.begin(), points.end(), [&](const Vec3& point)
        {
            return point.x() >= plane.bounds.lower.x() - tolerance * 16.0 &&
                   point.x() <= plane.bounds.upper.x() + tolerance * 16.0 &&
                   point.y() >= plane.bounds.lower.y() - tolerance * 16.0 &&
                   point.y() <= plane.bounds.upper.y() + tolerance * 16.0 &&
                   point.z() >= plane.bounds.lower.z() - tolerance * 16.0 &&
                   point.z() <= plane.bounds.upper.z() + tolerance * 16.0;
        });
    };

    std::vector<OutputPrimitive> result;
    result.reserve(primitives.size());
    for (OutputPrimitive& item : primitives)
    {
        if (item.primitive.kind != Kind::Polygon ||
            item.primitive.polygon.size() < 3)
        {
            result.push_back(std::move(item));
            continue;
        }
        std::vector<const Plane*> active_planes;
        for (const Plane& plane : planes)
        {
            double minimum = std::numeric_limits<double>::infinity();
            double maximum = -std::numeric_limits<double>::infinity();
            for (const Vec3& point : item.primitive.polygon)
            {
                const double side = plane.keep_sign *
                    (plane.normal.dot(point) - plane.offset);
                minimum = std::min(minimum, side);
                maximum = std::max(maximum, side);
            }
            if (minimum >= -tolerance || maximum <= tolerance) continue;
            if (!candidateMeetsFootprint(item.primitive, plane)) continue;
            active_planes.push_back(&plane);
        }
        if (active_planes.empty())
        {
            result.push_back(std::move(item));
            continue;
        }

        const PrimitiveMesh mesh = triangulatePrimitive(item.primitive);
        std::vector<std::vector<Vec3>> pieces;
        pieces.reserve(mesh.faces.size());
        for (const Face& face : mesh.faces)
            pieces.push_back({mesh.vertices[face[0]], mesh.vertices[face[1]],
                              mesh.vertices[face[2]]});
        const std::size_t input_triangles = pieces.size();
        for (const Plane* plane : active_planes)
        {
            std::vector<std::vector<Vec3>> next;
            next.reserve(pieces.size() * 2);
            for (const auto& piece : pieces)
            {
                double minimum = std::numeric_limits<double>::infinity();
                double maximum = -std::numeric_limits<double>::infinity();
                for (const Vec3& point : piece)
                {
                    const double side = plane->keep_sign *
                        (plane->normal.dot(point) - plane->offset);
                    minimum = std::min(minimum, side);
                    maximum = std::max(maximum, side);
                }
                if (minimum >= -tolerance)
                {
                    next.push_back(piece);
                    continue;
                }
                if (maximum <= tolerance) continue;
                auto clipped = clippedHalfspace(piece, *plane);
                if (clipped.size() >= 3) next.push_back(std::move(clipped));
            }
            pieces = std::move(next);
            if (pieces.empty()) break;
        }
        if (pieces.empty())
        {
            ++stats.removed_primitives;
            continue;
        }

        Vec3 normal;
        if (!planarNormal(item.primitive, tolerance, normal))
        {
            result.push_back(std::move(item));
            continue;
        }
        const Vec3 origin = item.primitive.polygon.front();
        const Mat3 basis = orthonormalFrame(normal);
        Mat3 frame;
        frame.col(0) = basis.col(1);
        frame.col(1) = basis.col(2);
        frame.col(2) = normal;
        Clipper2Lib::PathsD paths;
        for (const auto& piece : pieces)
        {
            Clipper2Lib::PathD path;
            for (const Vec3& point : piece)
            {
                const Vec3 local = frame.transposeMultiply(point - origin);
                path.emplace_back(local.x(), local.y());
            }
            if (std::abs(Clipper2Lib::Area(path)) <= tolerance * tolerance)
                continue;
            if (Clipper2Lib::Area(path) < 0.0) std::reverse(path.begin(), path.end());
            paths.push_back(std::move(path));
        }
        auto united = Clipper2Lib::Union(
            paths, Clipper2Lib::FillRule::NonZero, clipper_precision);
        if (united.empty())
        {
            ++stats.removed_primitives;
            continue;
        }
        if (std::any_of(united.begin(), united.end(),
            [](const auto& path) { return !Clipper2Lib::IsPositive(path); }))
        {
            Clipper2Lib::PathsD triangles;
            if (Clipper2Lib::Triangulate(
                    united, clipper_precision, triangles, false) !=
                Clipper2Lib::TriangulateResult::success)
            {
                result.push_back(std::move(item));
                continue;
            }
            united = std::move(triangles);
        }
        std::size_t emitted = 0;
        for (const auto& path : united)
        {
            std::vector<Vec2> boundary;
            for (const auto& point : path)
                boundary.emplace_back(point.x, point.y);
            boundary = simplifyPolygon(std::move(boundary), tolerance);
            if (boundary.size() < 3 ||
                triangulatePolygon(boundary).size() + 2 != boundary.size())
                continue;
            OutputPrimitive fragment = item;
            fragment.primitive = polygonPrimitive(boundary, origin, frame);
            fragment.preserves_cavity_opening = true;
            result.push_back(std::move(fragment));
            ++emitted;
        }
        if (emitted == 0)
        {
            result.push_back(std::move(item));
            continue;
        }
        ++stats.clipped_primitives;
        stats.output_fragments += emitted;
        (void)input_triangles;
    }
    return result;
}

// Build geometric face components while refusing sharp creases.  This is used
// for analytic-surface region growing: a cylindrical side must be allowed to
// cross the small angle between adjacent tessellation strips, but it must not
// absorb an attached cap, rib, or box wall merely because the OBJ happens to be
// one connected mesh.  OBJ indices are not authoritative here because CAD
// exporters commonly duplicate vertices at every patch boundary.
std::vector<std::vector<std::uint32_t>> smoothFaceComponentsApproximate(
    const Mesh& mesh,
    const std::vector<bool>& included_faces,
    const double tolerance,
    const double minimum_normal_dot,
    std::vector<std::vector<std::uint32_t>>* adjacency_output)
{
    struct CellKey
    {
        std::array<std::int64_t, 3> value{};
        bool operator<(const CellKey& other) const { return value < other.value; }
    };
    const double weld_tolerance = std::max(tolerance, 1.0e-12);
    const auto cellFor = [&](const Vec3& point)
    {
        CellKey key;
        for (int axis = 0; axis < 3; ++axis)
            key.value[axis] = static_cast<std::int64_t>(
                std::floor(point[axis] / weld_tolerance));
        return key;
    };
    std::map<CellKey, std::vector<std::uint32_t>> cells;
    std::vector<Vec3> canonical_points;
    std::vector<std::uint32_t> remap(mesh.vertices.size());
    for (std::uint32_t vertex = 0; vertex < mesh.vertices.size(); ++vertex)
    {
        const CellKey center = cellFor(mesh.vertices[vertex]);
        std::uint32_t selected = std::numeric_limits<std::uint32_t>::max();
        for (int dx = -1; dx <= 1 && selected ==
                 std::numeric_limits<std::uint32_t>::max(); ++dx)
            for (int dy = -1; dy <= 1 && selected ==
                     std::numeric_limits<std::uint32_t>::max(); ++dy)
                for (int dz = -1; dz <= 1 && selected ==
                         std::numeric_limits<std::uint32_t>::max(); ++dz)
                {
                    CellKey neighbor = center;
                    neighbor.value[0] += dx;
                    neighbor.value[1] += dy;
                    neighbor.value[2] += dz;
                    const auto found = cells.find(neighbor);
                    if (found == cells.end()) continue;
                    for (const auto candidate : found->second)
                        if ((canonical_points[candidate] -
                             mesh.vertices[vertex]).norm() <= weld_tolerance)
                        {
                            selected = candidate;
                            break;
                        }
                }
        if (selected == std::numeric_limits<std::uint32_t>::max())
        {
            selected = static_cast<std::uint32_t>(canonical_points.size());
            canonical_points.push_back(mesh.vertices[vertex]);
            cells[center].push_back(selected);
        }
        remap[vertex] = selected;
    }

    std::vector<Vec3> normals(mesh.faces.size(), Vec3::Zero());
    std::unordered_map<std::uint64_t, std::vector<std::uint32_t>> edge_faces;
    edge_faces.reserve(mesh.faces.size() * 3);
    for (std::uint32_t face_id = 0; face_id < mesh.faces.size(); ++face_id)
    {
        if (!included_faces[face_id]) continue;
        const Face& face = mesh.faces[face_id];
        Vec3 normal = (mesh.vertices[face[1]] - mesh.vertices[face[0]])
            .cross(mesh.vertices[face[2]] - mesh.vertices[face[0]]);
        if (normal.norm() > 1.0e-30) normals[face_id] = normal.normalized();
        for (int edge = 0; edge < 3; ++edge)
            edge_faces[edgeKey(remap[face[edge]],
                               remap[face[(edge + 1) % 3]])].push_back(face_id);
    }
    std::vector<std::vector<std::uint32_t>> adjacency(mesh.faces.size());
    for (const auto& [key, incident] : edge_faces)
    {
        (void)key;
        for (std::size_t first = 0; first < incident.size(); ++first)
            for (std::size_t second = first + 1; second < incident.size(); ++second)
            {
                const auto a = incident[first];
                const auto b = incident[second];
                if (std::abs(normals[a].dot(normals[b])) < minimum_normal_dot)
                    continue;
                adjacency[a].push_back(b);
                adjacency[b].push_back(a);
            }
    }
    std::vector<bool> visited(mesh.faces.size(), false);
    std::vector<std::vector<std::uint32_t>> result;
    for (std::uint32_t seed = 0; seed < mesh.faces.size(); ++seed)
    {
        if (!included_faces[seed] || visited[seed]) continue;
        result.emplace_back();
        std::queue<std::uint32_t> queue;
        queue.push(seed);
        visited[seed] = true;
        while (!queue.empty())
        {
            const auto current = queue.front();
            queue.pop();
            result.back().push_back(current);
            for (const auto neighbor : adjacency[current])
                if (!visited[neighbor])
                {
                    visited[neighbor] = true;
                    queue.push(neighbor);
                }
        }
    }
    if (adjacency_output != nullptr)
        *adjacency_output = std::move(adjacency);
    return result;
}

struct ConvexHullSurface
{
    std::vector<OutputPrimitive> shell;
    std::vector<Vec3> vertices;
    std::vector<Face> faces;
    std::vector<Vec3> plane_normals;
    std::vector<double> plane_offsets;
    Bounds bounds;
};

std::optional<ConvexHullSurface> fitConvexHullSurface(
    std::vector<Vec3> points,
    const std::vector<std::uint32_t>& responsibility,
    const std::uint64_t enclosure_group,
    const double tolerance)
{
    std::sort(points.begin(), points.end(),
        [](const Vec3& first, const Vec3& second)
        {
            if (first.x() != second.x()) return first.x() < second.x();
            if (first.y() != second.y()) return first.y() < second.y();
            return first.z() < second.z();
        });
    points.erase(std::unique(points.begin(), points.end(),
        [&](const Vec3& first, const Vec3& second)
        { return (first - second).norm() <= tolerance; }), points.end());
    if (points.size() < 4) return std::nullopt;

    Mesh point_mesh;
    point_mesh.vertices = points;
    std::vector<std::uint32_t> point_ids(points.size());
    std::iota(point_ids.begin(), point_ids.end(), 0);
    const BoxFit dimensionality = fitBox(point_mesh, point_ids);
    if (std::min({dimensionality.half_size.x(), dimensionality.half_size.y(),
                  dimensionality.half_size.z()}) <= tolerance)
        return std::nullopt;

    Bounds input_bounds;
    std::vector<quickhull::Vector3<double>> cloud;
    cloud.reserve(points.size());
    for (const Vec3& point : points)
    {
        input_bounds.lower = input_bounds.lower.cwiseMin(point);
        input_bounds.upper = input_bounds.upper.cwiseMax(point);
        cloud.emplace_back(point.x(), point.y(), point.z());
    }
    const double scale = std::max(
        (input_bounds.upper - input_bounds.lower).norm(), tolerance);
    quickhull::QuickHull<double> builder;
    const auto hull = builder.getConvexHull(
        cloud, true, false, std::max(tolerance / scale, 1.0e-12));
    if (hull.getVertexBuffer().size() < 4 ||
        hull.getIndexBuffer().size() < 12 ||
        hull.getIndexBuffer().size() % 3 != 0)
        return std::nullopt;

    ConvexHullSurface result;
    result.vertices.reserve(hull.getVertexBuffer().size());
    Vec3 center = Vec3::Zero();
    for (const auto& vertex : hull.getVertexBuffer())
    {
        const Vec3 converted(vertex.x, vertex.y, vertex.z);
        result.vertices.push_back(converted);
        result.bounds.lower = result.bounds.lower.cwiseMin(converted);
        result.bounds.upper = result.bounds.upper.cwiseMax(converted);
        center += converted;
    }
    center /= static_cast<double>(result.vertices.size());
    const auto& indices = hull.getIndexBuffer();
    result.faces.reserve(indices.size() / 3);
    result.shell.reserve(indices.size() / 3);
    for (std::size_t offset = 0; offset < indices.size(); offset += 3)
    {
        Face face{
            static_cast<std::uint32_t>(indices[offset]),
            static_cast<std::uint32_t>(indices[offset + 1]),
            static_cast<std::uint32_t>(indices[offset + 2])};
        Vec3 normal = (result.vertices[face[1]] - result.vertices[face[0]])
                          .cross(result.vertices[face[2]] -
                                 result.vertices[face[0]]);
        if (normal.norm() <= tolerance * tolerance) continue;
        if (normal.dot(result.vertices[face[0]] - center) < 0.0)
            std::swap(face[1], face[2]);
        normal = (result.vertices[face[1]] - result.vertices[face[0]])
                     .cross(result.vertices[face[2]] - result.vertices[face[0]])
                     .normalized();
        result.faces.push_back(face);
        result.plane_normals.push_back(normal);
        result.plane_offsets.push_back(
            normal.dot(result.vertices[face[0]]));
        Primitive triangle;
        triangle.kind = Kind::Polygon;
        triangle.polygon = {result.vertices[face[0]], result.vertices[face[1]],
                            result.vertices[face[2]]};
        OutputPrimitive item;
        item.primitive = std::move(triangle);
        item.source_faces = responsibility;
        item.enclosure_group = enclosure_group;
        result.shell.push_back(std::move(item));
    }
    if (result.faces.size() < 4) return std::nullopt;
    return result;
}

bool convexBoundsOverlap(const Bounds& first, const Bounds& second,
                         const double tolerance)
{
    for (int axis = 0; axis < 3; ++axis)
        if (first.upper[axis] < second.lower[axis] - tolerance ||
            second.upper[axis] < first.lower[axis] - tolerance)
            return false;
    return true;
}

bool boundsEntirelyOutsideHull(const ConvexHullSurface& hull,
                               const Bounds& bounds,
                               const double tolerance)
{
    for (std::size_t face = 0; face < hull.faces.size(); ++face)
    {
        const Vec3& normal = hull.plane_normals[face];
        double minimum_projection = 0.0;
        for (int axis = 0; axis < 3; ++axis)
            minimum_projection += normal[axis] >= 0.0
                ? normal[axis] * bounds.lower[axis]
                : normal[axis] * bounds.upper[axis];
        if (minimum_projection - hull.plane_offsets[face] > tolerance)
            return true;
    }
    return false;
}

bool pointStrictlyInsideHull(const ConvexHullSurface& hull,
                             const Vec3& point,
                             const double tolerance)
{
    for (std::size_t face = 0; face < hull.faces.size(); ++face)
        if (hull.plane_normals[face].dot(point) -
                hull.plane_offsets[face] > -tolerance)
            return false;
    return true;
}

bool segmentIntersectsHullInterior(const ConvexHullSurface& hull,
                                   const Vec3& first,
                                   const Vec3& second,
                                   const double tolerance)
{
    double lower = 0.0;
    double upper = 1.0;
    for (std::size_t face = 0; face < hull.faces.size(); ++face)
    {
        const Vec3& normal = hull.plane_normals[face];
        const double first_distance =
            normal.dot(first) - hull.plane_offsets[face];
        const double delta = normal.dot(second - first);
        const double limit = -tolerance;
        if (std::abs(delta) <= 1.0e-30)
        {
            if (first_distance > limit) return false;
            continue;
        }
        const double crossing = (limit - first_distance) / delta;
        if (delta > 0.0) upper = std::min(upper, crossing);
        else lower = std::max(lower, crossing);
        if (lower > upper) return false;
    }
    return upper >= 0.0 && lower <= 1.0 &&
           std::max(lower, 0.0) <= std::min(upper, 1.0);
}

bool segmentIntersectsTriangleInterior(
    const Vec3& first,
    const Vec3& second,
    const std::array<Vec3, 3>& triangle,
    const double tolerance)
{
    const Vec3 direction = second - first;
    const Vec3 edge1 = triangle[1] - triangle[0];
    const Vec3 edge2 = triangle[2] - triangle[0];
    const Vec3 cross = direction.cross(edge2);
    const double determinant = edge1.dot(cross);
    const double scale = std::max(
        {direction.norm(), edge1.norm(), edge2.norm(), 1.0});
    if (std::abs(determinant) <= tolerance * scale * scale) return false;
    const double inverse = 1.0 / determinant;
    const Vec3 relative = first - triangle[0];
    const double first_barycentric = relative.dot(cross) * inverse;
    const Vec3 second_cross = relative.cross(edge1);
    const double second_barycentric =
        direction.dot(second_cross) * inverse;
    const double segment_parameter = edge2.dot(second_cross) * inverse;
    const double barycentric_tolerance = tolerance / scale;
    return first_barycentric > barycentric_tolerance &&
           second_barycentric > barycentric_tolerance &&
           first_barycentric + second_barycentric <
               1.0 - barycentric_tolerance &&
           segment_parameter > barycentric_tolerance &&
           segment_parameter < 1.0 - barycentric_tolerance;
}

struct PlanarProjectionBounds
{
    Vec3 first_axis = Vec3::Zero();
    Vec3 second_axis = Vec3::Zero();
    double first_lower = std::numeric_limits<double>::infinity();
    double first_upper = -std::numeric_limits<double>::infinity();
    double second_lower = std::numeric_limits<double>::infinity();
    double second_upper = -std::numeric_limits<double>::infinity();
};

bool meshIntrudesConvexHull(const ConvexHullSurface& hull,
                            const PrimitiveMesh& mesh,
                            const Bounds& mesh_bounds,
                            const double tolerance,
                            bool* halfspace_rejected = nullptr,
                            const Vec3* mesh_plane_normal = nullptr,
                            const double mesh_plane_offset = 0.0,
                            bool* mesh_plane_rejected = nullptr,
                            const PlanarProjectionBounds* projection = nullptr,
                            bool* projection_rejected = nullptr)
{
    if (halfspace_rejected) *halfspace_rejected = false;
    if (mesh_plane_rejected) *mesh_plane_rejected = false;
    if (projection_rejected) *projection_rejected = false;
    if (!convexBoundsOverlap(hull.bounds, mesh_bounds, tolerance))
        return false;
    // A hull face plane that places every mesh vertex on its exterior side is
    // a separating plane for every triangle in the mesh. This proves that the
    // surface cannot enter the strict hull interior and avoids the much more
    // expensive edge/triangle intersection tests below.
    for (std::size_t plane = 0; plane < hull.faces.size(); ++plane)
    {
        bool all_outside = true;
        for (const Vec3& vertex : mesh.vertices)
            if (hull.plane_normals[plane].dot(vertex) -
                    hull.plane_offsets[plane] < -tolerance)
            {
                all_outside = false;
                break;
            }
        if (all_outside)
        {
            if (halfspace_rejected) *halfspace_rejected = true;
            return false;
        }
    }
    if (projection != nullptr)
    {
        double first_lower = std::numeric_limits<double>::infinity();
        double first_upper = -std::numeric_limits<double>::infinity();
        double second_lower = std::numeric_limits<double>::infinity();
        double second_upper = -std::numeric_limits<double>::infinity();
        for (const Vec3& vertex : hull.vertices)
        {
            const double first = projection->first_axis.dot(vertex);
            const double second = projection->second_axis.dot(vertex);
            first_lower = std::min(first_lower, first);
            first_upper = std::max(first_upper, first);
            second_lower = std::min(second_lower, second);
            second_upper = std::max(second_upper, second);
        }
        if (first_upper < projection->first_lower - tolerance ||
            projection->first_upper < first_lower - tolerance ||
            second_upper < projection->second_lower - tolerance ||
            projection->second_upper < second_lower - tolerance)
        {
            if (projection_rejected) *projection_rejected = true;
            return false;
        }
    }
    if (mesh_plane_normal != nullptr)
    {
        double minimum = std::numeric_limits<double>::infinity();
        double maximum = -std::numeric_limits<double>::infinity();
        for (const Vec3& vertex : hull.vertices)
        {
            const double distance =
                mesh_plane_normal->dot(vertex) - mesh_plane_offset;
            minimum = std::min(minimum, distance);
            maximum = std::max(maximum, distance);
        }
        if (minimum > tolerance || maximum < -tolerance)
        {
            if (mesh_plane_rejected) *mesh_plane_rejected = true;
            return false;
        }
    }
    for (const Vec3& vertex : mesh.vertices)
        if (pointStrictlyInsideHull(hull, vertex, tolerance)) return true;
    for (const Face& face : mesh.faces)
        for (int edge = 0; edge < 3; ++edge)
            if (segmentIntersectsHullInterior(
                    hull, mesh.vertices[face[edge]],
                    mesh.vertices[face[(edge + 1) % 3]], tolerance))
                return true;
    for (const Face& hull_face : hull.faces)
        for (int edge = 0; edge < 3; ++edge)
        {
            const Vec3& first = hull.vertices[hull_face[edge]];
            const Vec3& second =
                hull.vertices[hull_face[(edge + 1) % 3]];
            for (const Face& face : mesh.faces)
            {
                const std::array<Vec3, 3> triangle{
                    mesh.vertices[face[0]], mesh.vertices[face[1]],
                    mesh.vertices[face[2]]};
                if (segmentIntersectsTriangleInterior(
                        first, second, triangle, tolerance))
                    return true;
            }
        }
    return false;
}

bool primitiveIntrudesConvexHull(const ConvexHullSurface& hull,
                                 const OutputPrimitive& primitive,
                                 const double tolerance)
{
    const PrimitiveMesh mesh = triangulatePrimitive(primitive.primitive);
    return meshIntrudesConvexHull(
        hull, mesh, primitiveBounds(primitive.primitive), tolerance);
}

std::vector<Vec3> clipPolygonInsideConvexHull(
    std::vector<Vec3> polygon,
    const ConvexHullSurface& hull,
    const double tolerance)
{
    for (const Face& face : hull.faces)
    {
        const Vec3& plane_point = hull.vertices[face[0]];
        Vec3 normal = (hull.vertices[face[1]] - plane_point).cross(
            hull.vertices[face[2]] - plane_point);
        const double length = normal.norm();
        if (length <= 1.0e-30) continue;
        normal /= length;
        std::vector<Vec3> clipped;
        if (polygon.empty()) break;
        clipped.reserve(polygon.size() + 2);
        for (std::size_t index = 0; index < polygon.size(); ++index)
        {
            const Vec3& first = polygon[index];
            const Vec3& second = polygon[(index + 1) % polygon.size()];
            const double first_distance = normal.dot(first - plane_point);
            const double second_distance = normal.dot(second - plane_point);
            const bool first_inside = first_distance <= tolerance;
            const bool second_inside = second_distance <= tolerance;
            if (first_inside) clipped.push_back(first);
            if (first_inside == second_inside) continue;
            const double denominator = first_distance - second_distance;
            if (std::abs(denominator) <= 1.0e-30) continue;
            clipped.push_back(first + (second - first) *
                std::clamp(first_distance / denominator, 0.0, 1.0));
        }
        polygon = std::move(clipped);
    }
    return polygon;
}

std::vector<OutputPrimitive> convexHullUnionOuterSurface(
    const std::vector<ConvexHullSurface>& hulls,
    const double tolerance)
{
    constexpr int clipper_precision = 8;
    std::vector<OutputPrimitive> result;
    for (std::size_t hull_index = 0; hull_index < hulls.size(); ++hull_index)
    {
        const ConvexHullSurface& hull = hulls[hull_index];
        for (std::size_t face_index = 0;
             face_index < hull.faces.size(); ++face_index)
        {
            const Face& face = hull.faces[face_index];
            const Vec3& origin = hull.vertices[face[0]];
            Vec3 normal = (hull.vertices[face[1]] - origin).cross(
                hull.vertices[face[2]] - origin);
            const double normal_length = normal.norm();
            if (normal_length <= tolerance * tolerance) continue;
            normal /= normal_length;
            const Mat3 basis = orthonormalFrame(normal);
            Mat3 frame;
            frame.col(0) = basis.col(1);
            frame.col(1) = basis.col(2);
            frame.col(2) = normal;
            Clipper2Lib::PathD original_path;
            for (const auto vertex : face)
            {
                const Vec3 local = frame.transposeMultiply(
                    hull.vertices[vertex] - origin);
                original_path.emplace_back(local.x(), local.y());
            }
            if (Clipper2Lib::Area(original_path) < 0.0)
                std::reverse(original_path.begin(), original_path.end());
            Clipper2Lib::PathsD remaining{original_path};

            for (std::size_t other_index = 0;
                 other_index < hulls.size() && !remaining.empty();
                 ++other_index)
            {
                if (other_index == hull_index ||
                    !convexBoundsOverlap(
                        hull.bounds, hulls[other_index].bounds,
                        tolerance * 8.0))
                    continue;
                const ConvexHullSurface& other = hulls[other_index];
                bool same_facing_coplanar_boundary = false;
                for (const Face& other_face : other.faces)
                {
                    const Vec3& other_point = other.vertices[other_face[0]];
                    Vec3 other_normal =
                        (other.vertices[other_face[1]] - other_point).cross(
                            other.vertices[other_face[2]] - other_point);
                    const double length = other_normal.norm();
                    if (length <= 1.0e-30) continue;
                    other_normal /= length;
                    if (normal.dot(other_normal) < 1.0 - 1.0e-8 ||
                        std::abs(normal.dot(other_point - origin)) >
                            tolerance * 8.0)
                        continue;
                    same_facing_coplanar_boundary = true;
                    break;
                }
                // Coincident external faces are a duplicate, not an internal
                // interface. Keep the lower-index owner and subtract it from
                // later hulls; opposite-facing contact faces are removed from
                // both sides by the ordinary inside clipping below.
                if (same_facing_coplanar_boundary &&
                    hull_index < other_index)
                    continue;

                std::vector<Vec3> intersection = clipPolygonInsideConvexHull(
                    {hull.vertices[face[0]], hull.vertices[face[1]],
                     hull.vertices[face[2]]},
                    other, tolerance * 8.0);
                if (intersection.size() < 3) continue;
                Clipper2Lib::PathD clip_path;
                for (const Vec3& point : intersection)
                {
                    const Vec3 local =
                        frame.transposeMultiply(point - origin);
                    clip_path.emplace_back(local.x(), local.y());
                }
                if (std::abs(Clipper2Lib::Area(clip_path)) <=
                    tolerance * tolerance)
                    continue;
                if (Clipper2Lib::Area(clip_path) < 0.0)
                    std::reverse(clip_path.begin(), clip_path.end());
                remaining = Clipper2Lib::Difference(
                    remaining, Clipper2Lib::PathsD{clip_path},
                    Clipper2Lib::FillRule::NonZero, clipper_precision);
            }
            if (remaining.empty()) continue;
            Clipper2Lib::PathsD triangles;
            if (Clipper2Lib::Triangulate(
                    remaining, clipper_precision, triangles, false) !=
                Clipper2Lib::TriangulateResult::success)
                continue;
            for (const auto& triangle_path : triangles)
            {
                if (triangle_path.size() != 3) continue;
                std::vector<Vec2> polygon;
                polygon.reserve(3);
                for (const auto& point : triangle_path)
                    polygon.emplace_back(point.x, point.y);
                Primitive surface = polygonPrimitive(polygon, origin, frame);
                if (primitiveSurfaceArea(surface, tolerance) <=
                    tolerance * tolerance)
                    continue;
                OutputPrimitive item;
                item.primitive = std::move(surface);
                item.source_faces = hull.shell[face_index].source_faces;
                result.push_back(std::move(item));
            }
        }
    }
    return result;
}

struct ConvexPartitionDiagnostics
{
    std::size_t levels_tested = 0;
    std::size_t degenerate_levels = 0;
    std::size_t workload_rejections = 0;
    std::size_t overlap_rejections = 0;
    std::size_t error_rejections = 0;
    std::size_t accepted_parts = 0;
    std::size_t accepted_triangles = 0;
    double accepted_error = 0.0;
    std::vector<std::size_t> tested_triangles;
    std::vector<double> tested_errors;
};

struct ConvexPartitionResult
{
    std::vector<ConvexHullSurface> hulls;
    std::vector<OutputPrimitive> outer_surface;
};

struct ConvexPartitionFragment
{
    std::vector<Vec3> polygon;
    std::uint32_t source_face = 0;
};

struct ConvexPartitionCell
{
    std::vector<ConvexPartitionFragment> fragments;
};

std::vector<Vec3> clipPolygonToAxisHalfspace(
    const std::vector<Vec3>& polygon, const int axis,
    const double coordinate, const bool keep_lower)
{
    std::vector<Vec3> clipped;
    clipped.reserve(polygon.size() + 2);
    const auto inside = [&](const Vec3& point)
    {
        return keep_lower ? point[axis] <= coordinate
                          : point[axis] >= coordinate;
    };
    for (std::size_t index = 0; index < polygon.size(); ++index)
    {
        const Vec3& first = polygon[index];
        const Vec3& second = polygon[(index + 1) % polygon.size()];
        const bool first_inside = inside(first);
        const bool second_inside = inside(second);
        if (first_inside) clipped.push_back(first);
        if (first_inside == second_inside) continue;
        const double denominator = second[axis] - first[axis];
        if (std::abs(denominator) <= 1.0e-30) continue;
        const double parameter =
            std::clamp((coordinate - first[axis]) / denominator, 0.0, 1.0);
        Vec3 intersection = first + (second - first) * parameter;
        intersection[axis] = coordinate;
        clipped.push_back(intersection);
    }
    return clipped;
}

std::optional<std::array<ConvexPartitionCell, 2>> splitConvexPartitionCell(
    const ConvexPartitionCell& cell, const int axis, const double coordinate,
    const double tolerance)
{
    std::array<ConvexPartitionCell, 2> children;
    children[0].fragments.reserve(cell.fragments.size());
    children[1].fragments.reserve(cell.fragments.size());
    for (const auto& fragment : cell.fragments)
    {
        double lower = std::numeric_limits<double>::infinity();
        double upper = -std::numeric_limits<double>::infinity();
        for (const Vec3& point : fragment.polygon)
        {
            lower = std::min(lower, point[axis]);
            upper = std::max(upper, point[axis]);
        }
        if (upper <= coordinate + tolerance)
        {
            children[0].fragments.push_back(fragment);
            continue;
        }
        if (lower >= coordinate - tolerance)
        {
            children[1].fragments.push_back(fragment);
            continue;
        }
        auto lower_polygon = clipPolygonToAxisHalfspace(
            fragment.polygon, axis, coordinate, true);
        auto upper_polygon = clipPolygonToAxisHalfspace(
            fragment.polygon, axis, coordinate, false);
        if (lower_polygon.size() >= 3)
            children[0].fragments.push_back(
                {std::move(lower_polygon), fragment.source_face});
        if (upper_polygon.size() >= 3)
            children[1].fragments.push_back(
                {std::move(upper_polygon), fragment.source_face});
    }
    if (children[0].fragments.empty() || children[1].fragments.empty())
        return std::nullopt;
    return children;
}

std::optional<ConvexPartitionResult>
fitNonOverlappingConvexPartition(
    const Mesh& responsibility_mesh,
    const Mesh& distance_reference,
    const std::vector<std::uint32_t>& component,
    const std::size_t input_workload,
    const std::uint64_t first_enclosure_group,
    const double tolerance,
    const double maximum_error_distance,
    const double error_sample_spacing,
    ConvexPartitionDiagnostics& diagnostics)
{
    diagnostics = {};
    if (component.size() < 8) return std::nullopt;
    ConvexPartitionCell initial;
    initial.fragments.reserve(component.size());
    for (const auto face_id : component)
    {
        const Face& face = responsibility_mesh.faces[face_id];
        initial.fragments.push_back({
            {responsibility_mesh.vertices[face[0]],
             responsibility_mesh.vertices[face[1]],
             responsibility_mesh.vertices[face[2]]},
            face_id});
    }
    std::vector<ConvexPartitionCell> partitions;
    partitions.push_back(std::move(initial));

    const auto fitCells = [&](const std::vector<ConvexPartitionCell>& cells)
        -> std::optional<ConvexPartitionResult>
    {
        ConvexPartitionResult result;
        result.hulls.reserve(cells.size());
        for (std::size_t index = 0; index < cells.size(); ++index)
        {
            std::vector<Vec3> points;
            std::vector<std::uint32_t> responsibility;
            for (const auto& fragment : cells[index].fragments)
            {
                points.insert(points.end(), fragment.polygon.begin(),
                              fragment.polygon.end());
                responsibility.push_back(fragment.source_face);
            }
            std::sort(responsibility.begin(), responsibility.end());
            responsibility.erase(
                std::unique(responsibility.begin(), responsibility.end()),
                responsibility.end());
            auto hull = fitConvexHullSurface(
                std::move(points), responsibility,
                first_enclosure_group + index, tolerance);
            if (!hull) return std::nullopt;
            result.hulls.push_back(std::move(*hull));
        }
        result.outer_surface =
            convexHullUnionOuterSurface(result.hulls, tolerance);
        if (result.outer_surface.empty()) return std::nullopt;
        return result;
    };

    // The first cut determines whether a coarse decomposition becomes a useful
    // proxy or merely exposes a distant artificial section. Search a bounded,
    // geometry-only family of axis/quantile cuts instead of committing to one
    // arbitrary longest-axis median. Later levels refine the best-error cut.
    std::optional<ConvexPartitionResult> best_accepted;
    std::vector<ConvexPartitionCell> best_seed_partition;
    std::size_t best_accepted_workload =
        std::numeric_limits<std::size_t>::max();
    double best_accepted_error = std::numeric_limits<double>::infinity();
    double best_seed_error = std::numeric_limits<double>::infinity();
    Bounds initial_bounds;
    std::array<std::vector<double>, 3> centroid_coordinates;
    for (const auto& fragment : partitions.front().fragments)
    {
        Vec3 centroid = Vec3::Zero();
        for (const Vec3& point : fragment.polygon)
        {
            initial_bounds.lower = initial_bounds.lower.cwiseMin(point);
            initial_bounds.upper = initial_bounds.upper.cwiseMax(point);
            centroid += point;
        }
        centroid /= static_cast<double>(fragment.polygon.size());
        for (int axis = 0; axis < 3; ++axis)
            centroid_coordinates[axis].push_back(centroid[axis]);
    }
    constexpr std::array<double, 5> quantiles{
        0.2, 0.35, 0.5, 0.65, 0.8};
    for (int axis = 0; axis < 3; ++axis)
    {
        if (initial_bounds.upper[axis] - initial_bounds.lower[axis] <=
            tolerance * 8.0)
            continue;
        auto coordinates = centroid_coordinates[axis];
        std::sort(coordinates.begin(), coordinates.end());
        for (const double quantile : quantiles)
        {
            const std::size_t coordinate_index = std::min(
                coordinates.size() - 1,
                static_cast<std::size_t>(
                    quantile * static_cast<double>(coordinates.size() - 1)));
            double coordinate = coordinates[coordinate_index];
            if (coordinate <= initial_bounds.lower[axis] + tolerance * 8.0 ||
                coordinate >= initial_bounds.upper[axis] - tolerance * 8.0)
                continue;
            auto children = splitConvexPartitionCell(
                partitions.front(), axis, coordinate, tolerance);
            if (!children) continue;
            std::vector<ConvexPartitionCell> cells;
            cells.push_back(std::move((*children)[0]));
            cells.push_back(std::move((*children)[1]));
            auto candidate = fitCells(cells);
            if (!candidate) continue;
            const std::size_t workload =
                triangulatedFaceCount(candidate->outer_surface);
            if (workload >= input_workload) continue;
            const double error = maximumFilledSurfaceDistance(
                distance_reference, candidate->outer_surface,
                error_sample_spacing,
                std::numeric_limits<double>::infinity());
            diagnostics.tested_triangles.push_back(workload);
            diagnostics.tested_errors.push_back(error);
            if (error < best_seed_error)
            {
                best_seed_error = error;
                best_seed_partition = cells;
            }
            if (error <= maximum_error_distance + tolerance &&
                (workload < best_accepted_workload ||
                 (workload == best_accepted_workload &&
                  error < best_accepted_error)))
            {
                best_accepted_workload = workload;
                best_accepted_error = error;
                best_accepted = std::move(candidate);
            }
        }
    }
    if (best_accepted)
    {
        diagnostics.levels_tested = 1;
        diagnostics.accepted_parts = best_accepted->hulls.size();
        diagnostics.accepted_triangles = best_accepted_workload;
        diagnostics.accepted_error = best_accepted_error;
        return best_accepted;
    }
    if (best_seed_partition.empty()) return std::nullopt;
    partitions = std::move(best_seed_partition);

    constexpr std::size_t maximum_levels = 6;
    for (std::size_t level = 0; level < maximum_levels; ++level)
    {
        std::vector<ConvexPartitionCell> split_partitions;
        split_partitions.reserve(partitions.size() * 2);
        bool split_any = false;
        for (auto& cell : partitions)
        {
            if (cell.fragments.size() < 8)
            {
                split_partitions.push_back(std::move(cell));
                continue;
            }
            Bounds point_bounds;
            std::vector<Vec3> centroids;
            centroids.reserve(cell.fragments.size());
            for (const auto& fragment : cell.fragments)
            {
                Vec3 centroid = Vec3::Zero();
                for (const Vec3& point : fragment.polygon)
                {
                    point_bounds.lower = point_bounds.lower.cwiseMin(point);
                    point_bounds.upper = point_bounds.upper.cwiseMax(point);
                    centroid += point;
                }
                centroid /= static_cast<double>(fragment.polygon.size());
                centroids.push_back(centroid);
            }
            const Vec3 extent = point_bounds.upper - point_bounds.lower;
            int axis = 0;
            if (extent.y() > extent.x()) axis = 1;
            if (extent.z() > extent[axis]) axis = 2;
            if (extent[axis] <= tolerance * 8.0)
            {
                split_partitions.push_back(std::move(cell));
                continue;
            }
            const std::size_t middle = centroids.size() / 2;
            std::nth_element(
                centroids.begin(), centroids.begin() + middle,
                centroids.end(), [&](const Vec3& first, const Vec3& second)
                { return first[axis] < second[axis]; });
            double split = centroids[middle][axis];
            if (split <= point_bounds.lower[axis] + tolerance * 8.0 ||
                split >= point_bounds.upper[axis] - tolerance * 8.0)
                split = 0.5 * (point_bounds.lower[axis] +
                               point_bounds.upper[axis]);

            ConvexPartitionCell lower_cell;
            ConvexPartitionCell upper_cell;
            lower_cell.fragments.reserve(cell.fragments.size());
            upper_cell.fragments.reserve(cell.fragments.size());
            for (auto& fragment : cell.fragments)
            {
                double lower = std::numeric_limits<double>::infinity();
                double upper = -std::numeric_limits<double>::infinity();
                for (const Vec3& point : fragment.polygon)
                {
                    lower = std::min(lower, point[axis]);
                    upper = std::max(upper, point[axis]);
                }
                if (upper <= split + tolerance)
                {
                    lower_cell.fragments.push_back(std::move(fragment));
                    continue;
                }
                if (lower >= split - tolerance)
                {
                    upper_cell.fragments.push_back(std::move(fragment));
                    continue;
                }
                auto lower_polygon = clipPolygonToAxisHalfspace(
                    fragment.polygon, axis, split, true);
                auto upper_polygon = clipPolygonToAxisHalfspace(
                    fragment.polygon, axis, split, false);
                if (lower_polygon.size() >= 3)
                    lower_cell.fragments.push_back(
                        {std::move(lower_polygon), fragment.source_face});
                if (upper_polygon.size() >= 3)
                    upper_cell.fragments.push_back(
                        {std::move(upper_polygon), fragment.source_face});
            }
            if (lower_cell.fragments.empty() || upper_cell.fragments.empty())
            {
                split_partitions.push_back(std::move(cell));
                continue;
            }
            split_partitions.push_back(std::move(lower_cell));
            split_partitions.push_back(std::move(upper_cell));
            split_any = true;
        }
        partitions = std::move(split_partitions);
        if (!split_any) break;
        ++diagnostics.levels_tested;

        std::vector<ConvexHullSurface> hulls;
        hulls.reserve(partitions.size());
        bool degenerate = false;
        for (std::size_t index = 0; index < partitions.size(); ++index)
        {
            std::vector<Vec3> points;
            std::vector<std::uint32_t> responsibility;
            for (const auto& fragment : partitions[index].fragments)
            {
                points.insert(points.end(), fragment.polygon.begin(),
                              fragment.polygon.end());
                responsibility.push_back(fragment.source_face);
            }
            std::sort(responsibility.begin(), responsibility.end());
            responsibility.erase(
                std::unique(responsibility.begin(), responsibility.end()),
                responsibility.end());
            auto hull = fitConvexHullSurface(
                std::move(points), responsibility,
                first_enclosure_group + index, tolerance);
            if (!hull)
            {
                degenerate = true;
                break;
            }
            hulls.push_back(std::move(*hull));
        }
        if (degenerate)
        {
            ++diagnostics.degenerate_levels;
            continue;
        }
        std::vector<OutputPrimitive> candidate =
            convexHullUnionOuterSurface(hulls, tolerance);
        const std::size_t candidate_workload =
            triangulatedFaceCount(candidate);
        if (candidate.empty() || candidate_workload >= input_workload)
        {
            ++diagnostics.workload_rejections;
            continue;
        }
        const FilledSurfaceDistanceCertificate certificate =
            certifyFilledSurfaceDistance(
                distance_reference, candidate,
                maximum_error_distance + tolerance);
        bool accepted = certificate.certified && !certificate.exceeded;
        double observed_error = certificate.observed_maximum;
        if (!accepted && !certificate.exceeded)
        {
            observed_error = maximumFilledSurfaceDistance(
                distance_reference, candidate, error_sample_spacing,
                maximum_error_distance + tolerance);
            accepted = observed_error <= maximum_error_distance + tolerance;
        }
        diagnostics.tested_triangles.push_back(candidate_workload);
        diagnostics.tested_errors.push_back(observed_error);
        if (!accepted)
        {
            ++diagnostics.error_rejections;
            continue;
        }
        diagnostics.accepted_parts = hulls.size();
        diagnostics.accepted_triangles = candidate_workload;
        diagnostics.accepted_error = observed_error;
        return ConvexPartitionResult{
            std::move(hulls), std::move(candidate)};
    }
    return std::nullopt;
}

struct EnvelopeMergeGroup
{
    struct IntrusionSurface
    {
        PrimitiveMesh mesh;
        Bounds bounds;
        Vec3 plane_normal = Vec3::Zero();
        double plane_offset = 0.0;
        PlanarProjectionBounds projection;
        bool planar = false;
    };
    std::vector<OutputPrimitive> shell;
    std::vector<IntrusionSurface> intrusion_surfaces;
    std::vector<std::uint32_t> responsibility;
    std::vector<std::uint32_t> source_vertices;
    Bounds bounds;
    bool contains_round_surface = false;
    std::size_t triangle_workload = 0;
    std::unordered_set<std::size_t> neighbors;
    std::uint64_t version = 0;
    bool active = true;
};

struct EnvelopeFitResult
{
    std::vector<OutputPrimitive> shell;
    std::vector<std::uint32_t> responsibility;
    std::vector<std::uint32_t> source_vertices;
    Bounds bounds;
    bool contains_round_surface = false;
    std::size_t triangle_workload = 0;
    std::vector<std::size_t> consumed_groups;
};

std::vector<OutputPrimitive> mergeAdjacentEnvelopeGroups(
    const Mesh& responsibility_mesh,
    const Mesh& distance_reference,
    std::vector<OutputPrimitive> primitives,
    const PrimitiveMeshAnalysisOptions& options,
    const double tolerance,
    const double maximum_error_distance,
    const double error_sample_spacing,
    std::size_t& merged_group_count,
    std::vector<OutputPrimitive>& convex_enclosure_certificates,
    const std::filesystem::path& profile_path)
{
    const auto started = std::chrono::steady_clock::now();
    using Group = EnvelopeMergeGroup;
    struct Candidate
    {
        std::size_t first = 0;
        std::size_t second = 0;
        std::uint64_t first_version = 0;
        std::uint64_t second_version = 0;
        std::size_t priority = 0;
    };
    struct CandidateCompare
    {
        bool operator()(const Candidate& left,
                        const Candidate& right) const noexcept
        {
            if (left.priority != right.priority)
                return left.priority < right.priority;
            if (left.first != right.first) return left.first > right.first;
            return left.second > right.second;
        }
    };
    using FitResult = EnvelopeFitResult;
    struct Profile
    {
        std::size_t input_primitives = 0;
        std::size_t initial_adjacencies = 0;
        std::size_t candidate_evaluations = 0;
        std::size_t cache_hits = 0;
        std::size_t degenerate_rejections = 0;
        std::size_t error_rejections = 0;
        std::size_t certified_round_candidates = 0;
        std::size_t accepted_round_groups = 0;
        std::size_t conservative_round_component_attempts = 0;
        std::size_t accepted_conservative_round_components = 0;
        std::size_t convex_component_attempts = 0;
        std::size_t accepted_convex_components = 0;
        std::size_t accepted_convex_partition_components = 0;
        std::size_t invalid_convex_components = 0;
        std::size_t convex_workload_rejections = 0;
        std::size_t convex_error_rejections = 0;
        std::size_t convex_intersection_rejections = 0;
        std::size_t accepted_groups = 0;
        std::size_t coarse_distance_evaluations = 0;
        std::size_t coarse_distance_rejections = 0;
        std::size_t convex_face_cache_hits = 0;
        std::size_t convex_hull_fits = 0;
        std::size_t convex_growth_iterations = 0;
        std::size_t convex_debt_growth_steps = 0;
        std::size_t convex_debt_rejections = 0;
        std::size_t convex_error_cache_hits = 0;
        std::size_t convex_transition_cache_hits = 0;
        std::size_t intrusion_tests = 0;
        std::size_t intrusion_index_candidates = 0;
        std::size_t intrusion_index_rebuilds = 0;
        std::size_t intrusion_index_halfspace_rejections = 0;
        std::size_t intrusion_halfspace_rejections = 0;
        std::size_t intrusion_surface_plane_rejections = 0;
        std::size_t intrusion_planar_projection_rejections = 0;
        std::size_t certificate_samples = 0;
        std::size_t fine_distance_evaluations = 0;
        double coarse_distance_seconds = 0.0;
        double convex_hull_seconds = 0.0;
        double intrusion_seconds = 0.0;
        double distance_seconds = 0.0;
    } profile;
    profile.input_primitives = primitives.size();
    auto last_profile_flush = started - std::chrono::seconds(2);
    const auto flushProfile = [&](const bool force = false)
    {
        const auto now = std::chrono::steady_clock::now();
        if (!force && now - last_profile_flush < std::chrono::seconds(1))
            return;
        last_profile_flush = now;
        std::ofstream stream(profile_path);
        stream << std::setprecision(17)
               << "{\"complete\":false"
               << ",\"input_primitives\":" << profile.input_primitives
               << ",\"initial_adjacencies\":"
               << profile.initial_adjacencies
               << ",\"candidate_evaluations\":"
               << profile.candidate_evaluations
               << ",\"accepted_groups\":" << profile.accepted_groups
               << ",\"error_rejections\":" << profile.error_rejections
               << ",\"coarse_distance_evaluations\":"
               << profile.coarse_distance_evaluations
               << ",\"coarse_distance_rejections\":"
               << profile.coarse_distance_rejections
               << ",\"convex_face_cache_hits\":"
               << profile.convex_face_cache_hits
               << ",\"convex_hull_fits\":" << profile.convex_hull_fits
               << ",\"convex_growth_iterations\":"
               << profile.convex_growth_iterations
               << ",\"convex_debt_growth_steps\":"
               << profile.convex_debt_growth_steps
               << ",\"convex_debt_rejections\":"
               << profile.convex_debt_rejections
               << ",\"convex_error_cache_hits\":"
               << profile.convex_error_cache_hits
               << ",\"convex_transition_cache_hits\":"
               << profile.convex_transition_cache_hits
               << ",\"intrusion_tests\":" << profile.intrusion_tests
               << ",\"intrusion_index_candidates\":"
               << profile.intrusion_index_candidates
               << ",\"intrusion_index_rebuilds\":"
               << profile.intrusion_index_rebuilds
               << ",\"intrusion_index_halfspace_rejections\":"
               << profile.intrusion_index_halfspace_rejections
               << ",\"intrusion_halfspace_rejections\":"
               << profile.intrusion_halfspace_rejections
               << ",\"intrusion_surface_plane_rejections\":"
               << profile.intrusion_surface_plane_rejections
               << ",\"intrusion_planar_projection_rejections\":"
               << profile.intrusion_planar_projection_rejections
               << ",\"certificate_samples\":"
               << profile.certificate_samples
               << ",\"coarse_distance_seconds\":"
               << profile.coarse_distance_seconds
               << ",\"convex_hull_seconds\":"
               << profile.convex_hull_seconds
               << ",\"intrusion_seconds\":" << profile.intrusion_seconds
               << ",\"distance_seconds\":" << profile.distance_seconds
               << ",\"elapsed_seconds\":"
               << std::chrono::duration<double>(now - started).count()
               << "}\n";
    };
    merged_group_count = 0;
    if (primitives.size() < 2) return primitives;
    // The user-provided directed distance is the sole approximation budget.
    // Candidate certification and the final full-surface audit both use that
    // same limit; an undocumented local safety factor must not suppress an
    // otherwise admissible lower-work candidate.
    const double candidate_error_limit = maximum_error_distance;

    // Test complete geometric components before pairwise growth. A conservative
    // revolved envelope is admissible only when both end silhouettes are
    // certified circular, it reduces triangle work, and its directed error is
    // within the same user limit used by every later merge.
    std::uint64_t next_component_enclosure_group = 1;
    for (const auto& primitive : primitives)
        next_component_enclosure_group = std::max(
            next_component_enclosure_group, primitive.enclosure_group + 1);
    std::vector<bool> all_source_faces(responsibility_mesh.faces.size(), true);
    Bounds responsibility_bounds;
    for (const Vec3& vertex : responsibility_mesh.vertices)
    {
        responsibility_bounds.lower =
            responsibility_bounds.lower.cwiseMin(vertex);
        responsibility_bounds.upper =
            responsibility_bounds.upper.cwiseMax(vertex);
    }
    const double component_analytic_tolerance = std::max(
        tolerance * 8.0,
        (responsibility_bounds.upper - responsibility_bounds.lower).norm() *
            options.analytic_surface_relative_tolerance);
    auto source_components = faceComponentsApproximate(
        responsibility_mesh, all_source_faces,
        component_analytic_tolerance);
    if (source_components.size() > 1)
    {
        std::vector<std::uint32_t> complete_model_faces(
            responsibility_mesh.faces.size());
        std::iota(complete_model_faces.begin(), complete_model_faces.end(), 0);
        source_components.insert(source_components.begin(),
                                 std::move(complete_model_faces));
    }
    std::vector<bool> consumed_primitives(primitives.size(), false);
    std::vector<OutputPrimitive> component_replacements;
    std::ofstream convex_component_diagnostics(
        profile_path.parent_path() / "convex_component_profile.jsonl");
    for (const auto& component : source_components)
    {
        ++profile.conservative_round_component_attempts;
        ConservativeRevolvedDiagnostics fit_diagnostics;
        const auto revolved = fitConservativeRevolvedEnvelope(
            responsibility_mesh, component, options,
            component_analytic_tolerance,
            &fit_diagnostics);
        {
            std::ofstream diagnostic(profile_path.parent_path() /
                                     "conservative_round_fit_profile.jsonl",
                                     std::ios::app);
            diagnostic << "{\"component_faces\":" << component.size()
                       << ",\"axes_tested\":" << fit_diagnostics.axes_tested
                       << ",\"failure_geometry\":"
                       << fit_diagnostics.failures[0]
                       << ",\"failure_end_ring\":"
                       << fit_diagnostics.failures[2]
                       << ",\"failure_normal\":"
                       << fit_diagnostics.failures[3]
                       << ",\"failure_radial\":"
                       << fit_diagnostics.failures[4]
                       << ",\"relaxed_candidates\":"
                       << fit_diagnostics.relaxed_candidates
                       << ",\"circular_silhouette_candidates\":"
                       << fit_diagnostics.circular_silhouette_candidates
                       << "}\n";
        }
        if (!revolved) continue;
        std::unordered_set<std::uint32_t> component_faces(
            component.begin(), component.end());
        std::vector<std::size_t> owners;
        std::size_t input_workload = 0;
        bool partial_owner = false;
        for (std::size_t index = 0; index < primitives.size(); ++index)
        {
            if (consumed_primitives[index]) continue;
            if (primitives[index].source_faces.empty()) continue;
            const bool touches = std::any_of(
                primitives[index].source_faces.begin(),
                primitives[index].source_faces.end(),
                [&](const auto face) { return component_faces.contains(face); });
            if (!touches) continue;
            const bool fully_owned = std::all_of(
                primitives[index].source_faces.begin(),
                primitives[index].source_faces.end(),
                [&](const auto face) { return component_faces.contains(face); });
            if (!fully_owned)
            {
                partial_owner = true;
                break;
            }
            owners.push_back(index);
            input_workload +=
                triangulatePrimitive(primitives[index].primitive).faces.size();
        }
        if (partial_owner || owners.empty()) continue;

        std::vector<OutputPrimitive> candidate;
        Primitive side = *revolved;
        const std::uint32_t adaptive_segments = roundSegmentsForError(
            std::max(side.base_radius, side.top_radius),
            candidate_error_limit, side.segments);
        side.segments = adaptive_segments;
        side.kind = std::abs(side.base_radius - side.top_radius) <=
                1.0e-8 * std::max({side.base_radius, side.top_radius,
                                    side.height})
            ? Kind::CylindricalBand : Kind::ConicalBand;
        candidate.push_back({std::move(side), component});
        const auto append_cap = [&](const Vec3& center, const double radius)
        {
            Primitive disk;
            disk.kind = Kind::Disk;
            disk.center = center;
            disk.axes = revolved->axes;
            disk.base_radius = radius;
            disk.top_radius = radius;
            disk.segments = adaptive_segments;
            candidate.push_back({std::move(disk), component});
        };
        append_cap(revolved->center, revolved->base_radius);
        append_cap(revolved->center +
                       revolved->axes.col(0) * revolved->height,
                   revolved->top_radius);
        const std::size_t candidate_workload =
            triangulatedFaceCount(candidate);
        if (candidate.size() < 2 || candidate_workload >= input_workload)
        {
            std::ofstream diagnostic(profile_path.parent_path() /
                                     "conservative_round_candidate_profile.jsonl",
                                     std::ios::app);
            diagnostic << std::setprecision(17)
                       << "{\"component_faces\":" << component.size()
                       << ",\"height\":" << revolved->height
                       << ",\"base_radius\":" << revolved->base_radius
                       << ",\"top_radius\":" << revolved->top_radius
                       << ",\"candidate_primitives\":" << candidate.size()
                       << ",\"candidate_triangles\":" << candidate_workload
                       << ",\"input_triangles\":" << input_workload
                       << ",\"rejection\":\"workload_or_open_shell\"}\n";
            continue;
        }
        for (auto& item : candidate)
        {
            item.source_faces = component;
            item.enclosure_group = next_component_enclosure_group;
        }
        const FilledSurfaceDistanceCertificate certificate =
            certifyFilledSurfaceDistance(
                distance_reference, candidate,
                candidate_error_limit + tolerance);
        bool accepted = certificate.certified && !certificate.exceeded;
        double measured_distance = certificate.observed_maximum;
        if (!accepted && !certificate.exceeded)
        {
            measured_distance = maximumFilledSurfaceDistance(
                distance_reference, candidate, error_sample_spacing,
                candidate_error_limit + tolerance);
            accepted = measured_distance <= candidate_error_limit + tolerance;
        }
        {
            std::ofstream diagnostic(profile_path.parent_path() /
                                     "conservative_round_candidate_profile.jsonl",
                                     std::ios::app);
            diagnostic << std::setprecision(17)
                       << "{\"component_faces\":" << component.size()
                       << ",\"height\":" << revolved->height
                       << ",\"base_radius\":" << revolved->base_radius
                       << ",\"top_radius\":" << revolved->top_radius
                       << ",\"candidate_primitives\":" << candidate.size()
                       << ",\"candidate_triangles\":" << candidate_workload
                       << ",\"input_triangles\":" << input_workload
                       << ",\"certificate_exceeded\":"
                       << (certificate.exceeded ? "true" : "false")
                       << ",\"certificate_certified\":"
                       << (certificate.certified ? "true" : "false")
                       << ",\"observed_distance\":" << measured_distance
                       << ",\"limit\":" << candidate_error_limit
                       << ",\"accepted\":" << (accepted ? "true" : "false")
                       << "}\n";
        }
        if (!accepted) continue;
        for (const auto owner : owners) consumed_primitives[owner] = true;
        component_replacements.insert(
            component_replacements.end(),
            std::make_move_iterator(candidate.begin()),
            std::make_move_iterator(candidate.end()));
        ++next_component_enclosure_group;
        ++profile.accepted_conservative_round_components;
        merged_group_count += owners.size() - 1;
    }

    // A complete connected responsibility component may also compete as its
    // conservative convex surface.  This is deliberately a component-level
    // candidate only: the former pairwise convex growth greedily produced many
    // intersecting local hulls on noisy CAD surfaces.  Analytic round surfaces
    // remain preferred and are never flattened into a polyhedral hull.
    for (const auto& component : source_components)
    {
        std::unordered_set<std::uint32_t> component_faces(
            component.begin(), component.end());
        std::vector<std::size_t> owners;
        std::size_t input_workload = 0;
        bool partial_owner = false;
        for (std::size_t index = 0; index < primitives.size(); ++index)
        {
            if (consumed_primitives[index] ||
                primitives[index].source_faces.empty())
                continue;
            const bool touches = std::any_of(
                primitives[index].source_faces.begin(),
                primitives[index].source_faces.end(),
                [&](const auto face)
                { return component_faces.contains(face); });
            if (!touches) continue;
            const bool fully_owned = std::all_of(
                primitives[index].source_faces.begin(),
                primitives[index].source_faces.end(),
                [&](const auto face)
                { return component_faces.contains(face); });
            if (!fully_owned)
            {
                partial_owner = true;
                break;
            }
            owners.push_back(index);
            input_workload += triangulatePrimitive(
                primitives[index].primitive).faces.size();
        }
        // A component accepted by the preceding analytic-envelope pass has
        // already been marked consumed, so it cannot reach this candidate.
        // Do not veto a genuinely complex component merely because one small
        // disk or band is among its many current surface patches.
        if (partial_owner || owners.empty()) continue;
        ++profile.convex_component_attempts;

        // Convexification must grow from adjacent current primitives. A
        // complete-component hull is a shortcut that bypasses the merge fixed
        // point and is therefore not an admissible candidate.
        continue;

        // The bundled QuickHull implementation becomes quadratic on large,
        // nearly coplanar CAD point clouds. Do not spend minutes constructing
        // a global candidate when the already-recognized surface still has an
        // industrial-scale workload: those candidates overwhelmingly fail the
        // same directed-error audit and previously produced no output change.
        // A coarse component whose current proxy is already small (for example
        // a rail assembly) remains eligible regardless of source face count.
        constexpr std::size_t maximum_dense_hull_component_faces = 20000;
        constexpr std::size_t maximum_dense_hull_input_workload = 7000;
        if (component.size() > maximum_dense_hull_component_faces &&
            input_workload > maximum_dense_hull_input_workload)
        {
            ++profile.convex_workload_rejections;
            convex_component_diagnostics
                << "{\"component_faces\":" << component.size()
                << ",\"input_triangles\":" << input_workload
                << ",\"rejection\":\"dense_hull_cost\"}\n";
            continue;
        }

        std::vector<Vec3> points;
        for (const auto vertex : uniqueVertices(responsibility_mesh, component))
            points.push_back(responsibility_mesh.vertices[vertex]);
        auto hull = fitConvexHullSurface(
            std::move(points), component, next_component_enclosure_group,
            tolerance);
        if (!hull)
        {
            ++profile.invalid_convex_components;
            convex_component_diagnostics
                << "{\"component_faces\":" << component.size()
                << ",\"input_triangles\":" << input_workload
                << ",\"rejection\":\"degenerate_hull\"}\n";
            continue;
        }
        if (hull->faces.size() >= input_workload)
        {
            ++profile.convex_workload_rejections;
            convex_component_diagnostics
                << "{\"component_faces\":" << component.size()
                << ",\"input_triangles\":" << input_workload
                << ",\"hull_triangles\":" << hull->faces.size()
                << ",\"rejection\":\"workload\"}\n";
            continue;
        }

        const FilledSurfaceDistanceCertificate certificate =
            certifyFilledSurfaceDistance(
                distance_reference, hull->shell,
                candidate_error_limit + tolerance);
        bool accepted = certificate.certified && !certificate.exceeded;
        double observed_error = certificate.observed_maximum;
        if (!accepted && !certificate.exceeded)
        {
            observed_error = maximumFilledSurfaceDistance(
                distance_reference, hull->shell, error_sample_spacing,
                candidate_error_limit + tolerance);
            accepted = observed_error <= candidate_error_limit + tolerance;
        }

        std::vector<ConvexHullSurface> accepted_hulls;
        std::vector<OutputPrimitive> accepted_outer_surface;
        bool used_partition = false;
        if (accepted)
        {
            accepted_outer_surface = hull->shell;
            accepted_hulls.push_back(std::move(*hull));
        }
        else
        {
            ++profile.convex_error_rejections;
            convex_component_diagnostics << std::setprecision(17)
                << "{\"component_faces\":" << component.size()
                << ",\"input_triangles\":" << input_workload
                << ",\"whole_hull_triangles\":" << hull->faces.size()
                << ",\"whole_hull_observed_distance\":"
                << observed_error
                << ",\"limit\":" << candidate_error_limit
                << ",\"rejection\":\"maximum_error\"}\n";
            continue;
        }

        std::unordered_set<std::size_t> owner_set(
            owners.begin(), owners.end());
        std::vector<std::size_t> enclosed_generated_primitives;
        bool intersects_unowned = false;
        for (std::size_t index = 0; index < primitives.size(); ++index)
        {
            if (consumed_primitives[index] || owner_set.contains(index))
                continue;
            bool intrudes = false;
            for (const auto& accepted_hull : accepted_hulls)
                if (primitiveIntrudesConvexHull(
                        accepted_hull, primitives[index], tolerance * 8.0))
                {
                    intrudes = true;
                    break;
                }
            if (!intrudes) continue;
            if (primitives[index].source_faces.empty())
            {
                enclosed_generated_primitives.push_back(index);
                continue;
            }
            intersects_unowned = true;
            break;
        }
        if (intersects_unowned)
        {
            ++profile.convex_intersection_rejections;
            convex_component_diagnostics
                << "{\"component_faces\":" << component.size()
                << ",\"input_triangles\":" << input_workload
                << ",\"hull_parts\":" << accepted_hulls.size()
                << ",\"rejection\":\"unowned_intersection\"}\n";
            continue;
        }
        for (const auto owner : owners)
            consumed_primitives[owner] = true;
        for (const auto enclosed : enclosed_generated_primitives)
            consumed_primitives[enclosed] = true;
        const std::size_t accepted_triangle_count =
            triangulatedFaceCount(accepted_outer_surface);
        if (used_partition)
        {
            for (const auto& accepted_hull : accepted_hulls)
                convex_enclosure_certificates.insert(
                    convex_enclosure_certificates.end(),
                    accepted_hull.shell.begin(), accepted_hull.shell.end());
            // The clipped union is generally non-convex, so its triangles must
            // not share a convex enclosure id. The complete per-part hulls
            // above remain the conservative offline coverage certificates.
            for (auto& surface : accepted_outer_surface)
                surface.enclosure_group = 0;
        }
        component_replacements.insert(
            component_replacements.end(),
            std::make_move_iterator(accepted_outer_surface.begin()),
            std::make_move_iterator(accepted_outer_surface.end()));
        next_component_enclosure_group += accepted_hulls.size();
        ++profile.accepted_convex_components;
        if (used_partition)
            ++profile.accepted_convex_partition_components;
        merged_group_count += owners.size() - 1;
        if (!used_partition)
            convex_component_diagnostics << std::setprecision(17)
                << "{\"component_faces\":" << component.size()
                << ",\"input_triangles\":" << input_workload
                << ",\"hull_triangles\":" << accepted_triangle_count
                << ",\"observed_distance\":" << observed_error
                << ",\"limit\":" << candidate_error_limit
                << ",\"accepted\":true}\n";
    }
    if (!component_replacements.empty())
    {
        std::vector<OutputPrimitive> replaced;
        replaced.reserve(primitives.size() + component_replacements.size());
        for (std::size_t index = 0; index < primitives.size(); ++index)
            if (!consumed_primitives[index])
                replaced.push_back(std::move(primitives[index]));
        replaced.insert(replaced.end(),
                        std::make_move_iterator(component_replacements.begin()),
                        std::make_move_iterator(component_replacements.end()));
        primitives = std::move(replaced);
    }

    std::vector<Group> groups;
    groups.reserve(primitives.size());
    const auto makeIntrusionSurface = [&](const OutputPrimitive& output)
    {
        Group::IntrusionSurface result;
        result.mesh = triangulatePrimitive(output.primitive);
        result.bounds = primitiveBounds(output.primitive);
        result.planar = planarNormal(
            output.primitive, tolerance * 8.0, result.plane_normal);
        if (result.planar && !result.mesh.vertices.empty())
        {
            result.plane_offset = result.plane_normal.dot(
                result.mesh.vertices.front());
            const Mat3 frame = orthonormalFrame(result.plane_normal);
            result.projection.first_axis = frame.col(1);
            result.projection.second_axis = frame.col(2);
            for (const Vec3& vertex : result.mesh.vertices)
            {
                const double first =
                    result.projection.first_axis.dot(vertex);
                const double second =
                    result.projection.second_axis.dot(vertex);
                result.projection.first_lower = std::min(
                    result.projection.first_lower, first);
                result.projection.first_upper = std::max(
                    result.projection.first_upper, first);
                result.projection.second_lower = std::min(
                    result.projection.second_lower, second);
                result.projection.second_upper = std::max(
                    result.projection.second_upper, second);
            }
        }
        return result;
    };
    for (auto& primitive : primitives)
    {
        Group group;
        group.responsibility = primitive.source_faces;
        for (const auto face : group.responsibility)
            if (face < responsibility_mesh.faces.size())
                group.source_vertices.insert(
                    group.source_vertices.end(),
                    responsibility_mesh.faces[face].begin(),
                    responsibility_mesh.faces[face].end());
        std::sort(group.source_vertices.begin(), group.source_vertices.end());
        group.source_vertices.erase(
            std::unique(group.source_vertices.begin(),
                        group.source_vertices.end()),
            group.source_vertices.end());
        group.bounds = primitiveBounds(primitive.primitive);
        group.contains_round_surface =
            isCertifiedRoundSurfaceKind(primitive.primitive.kind);
        group.triangle_workload =
            triangulatePrimitive(primitive.primitive).faces.size();
        group.shell.push_back(std::move(primitive));
        group.intrusion_surfaces.push_back(
            makeIntrusionSurface(group.shell.front()));
        groups.push_back(std::move(group));
    }
    const auto connect = [&](const std::size_t first, const std::size_t second)
    {
        if (first == second) return;
        groups[first].neighbors.insert(second);
        groups[second].neighbors.insert(first);
    };

    // Responsibility adjacency is stable under a group merge and remains the
    // authoritative topology. It lets a closed envelope grow through the same
    // connected surface without requiring its newly fitted faces to retain all
    // historical pairwise intersections.
    std::vector<std::vector<std::size_t>> face_owners(
        responsibility_mesh.faces.size());
    for (std::size_t group = 0; group < groups.size(); ++group)
        for (const auto face : groups[group].shell.front().source_faces)
            if (face < face_owners.size()) face_owners[face].push_back(group);
    std::unordered_map<std::uint64_t, std::vector<std::uint32_t>> edge_faces;
    edge_faces.reserve(responsibility_mesh.faces.size() * 3);
    for (std::uint32_t face = 0; face < responsibility_mesh.faces.size(); ++face)
    {
        const Face& triangle = responsibility_mesh.faces[face];
        for (int edge = 0; edge < 3; ++edge)
            edge_faces[edgeKey(triangle[edge], triangle[(edge + 1) % 3])]
                .push_back(face);
        const auto& owners = face_owners[face];
        for (std::size_t first = 0; first < owners.size(); ++first)
            for (std::size_t second = first + 1; second < owners.size(); ++second)
                connect(owners[first], owners[second]);
    }
    for (const auto& [edge, faces] : edge_faces)
    {
        (void)edge;
        for (std::size_t first_face = 0; first_face < faces.size(); ++first_face)
            for (std::size_t second_face = first_face + 1;
                 second_face < faces.size(); ++second_face)
                for (const auto first : face_owners[faces[first_face]])
                    for (const auto second : face_owners[faces[second_face]])
                        connect(first, second);
    }

    // Restored caps and fitted surfaces may carry no source face. Add geometric
    // vertex adjacency so they still participate in the same fixed point.
    struct PointKey
    {
        std::int64_t x = 0;
        std::int64_t y = 0;
        std::int64_t z = 0;
        bool operator==(const PointKey&) const = default;
    };
    struct PointKeyHash
    {
        std::size_t operator()(const PointKey& point) const noexcept
        {
            std::size_t seed = std::hash<std::int64_t>{}(point.x);
            seed ^= std::hash<std::int64_t>{}(point.y) +
                    0x9e3779b97f4a7c15ULL + (seed << 6U) + (seed >> 2U);
            seed ^= std::hash<std::int64_t>{}(point.z) +
                    0x9e3779b97f4a7c15ULL + (seed << 6U) + (seed >> 2U);
            return seed;
        }
    };
    const double cell_size = std::max(tolerance * 8.0, 1.0e-12);
    const auto pointKey = [&](const Vec3& point)
    {
        return PointKey{
            static_cast<std::int64_t>(std::llround(point.x() / cell_size)),
            static_cast<std::int64_t>(std::llround(point.y() / cell_size)),
            static_cast<std::int64_t>(std::llround(point.z() / cell_size))};
    };
    std::unordered_map<PointKey, std::vector<std::size_t>, PointKeyHash> owners;
    for (std::size_t group = 0; group < groups.size(); ++group)
    {
        std::unordered_set<PointKey, PointKeyHash> keys;
        for (const auto& primitive : groups[group].shell)
        {
            const PrimitiveMesh surface = triangulatePrimitive(primitive.primitive);
            for (const Vec3& vertex : surface.vertices) keys.insert(pointKey(vertex));
        }
        for (const PointKey& key : keys)
        {
            for (std::int64_t dx = -1; dx <= 1; ++dx)
                for (std::int64_t dy = -1; dy <= 1; ++dy)
                    for (std::int64_t dz = -1; dz <= 1; ++dz)
                    {
                        const auto found = owners.find(
                            {key.x + dx, key.y + dy, key.z + dz});
                        if (found == owners.end()) continue;
                        for (const auto other : found->second) connect(group, other);
                    }
            owners[key].push_back(group);
        }
    }
    for (const Group& group : groups)
        profile.initial_adjacencies += group.neighbors.size();
    profile.initial_adjacencies /= 2;

    struct IndexedSurface
    {
        std::size_t group = 0;
        std::size_t surface = 0;
        std::uint64_t version = 0;
        Bounds bounds;
    };
    struct SurfaceBoundsNode
    {
        Bounds bounds;
        std::size_t begin = 0;
        std::size_t end = 0;
        std::size_t first_child = std::numeric_limits<std::size_t>::max();
        std::size_t second_child = std::numeric_limits<std::size_t>::max();
    };
    std::vector<IndexedSurface> surface_bounds_entries;
    std::vector<SurfaceBoundsNode> surface_bounds_nodes;
    std::vector<IndexedSurface> surface_bounds_delta;
    const auto rebuildSurfaceBoundsIndex = [&]
    {
        surface_bounds_entries.clear();
        for (std::size_t group = 0; group < groups.size(); ++group)
        {
            if (!groups[group].active) continue;
            for (std::size_t surface = 0;
                 surface < groups[group].intrusion_surfaces.size(); ++surface)
                surface_bounds_entries.push_back({
                    group, surface, groups[group].version,
                    groups[group].intrusion_surfaces[surface].bounds});
        }
        surface_bounds_nodes.clear();
        surface_bounds_nodes.reserve(surface_bounds_entries.size() * 2);
        const std::function<std::size_t(std::size_t, std::size_t)> build =
            [&](const std::size_t begin, const std::size_t end)
        {
            const std::size_t node_index = surface_bounds_nodes.size();
            surface_bounds_nodes.emplace_back();
            Bounds bounds;
            Bounds center_bounds;
            for (std::size_t index = begin; index < end; ++index)
            {
                const Bounds& item = surface_bounds_entries[index].bounds;
                bounds.lower = bounds.lower.cwiseMin(item.lower);
                bounds.upper = bounds.upper.cwiseMax(item.upper);
                const Vec3 center = (item.lower + item.upper) * 0.5;
                center_bounds.lower = center_bounds.lower.cwiseMin(center);
                center_bounds.upper = center_bounds.upper.cwiseMax(center);
            }
            surface_bounds_nodes[node_index].bounds = bounds;
            surface_bounds_nodes[node_index].begin = begin;
            surface_bounds_nodes[node_index].end = end;
            if (end - begin <= 8) return node_index;
            const Vec3 extent = center_bounds.upper - center_bounds.lower;
            int axis = 0;
            if (extent.y() > extent.x()) axis = 1;
            if (extent.z() > extent[axis]) axis = 2;
            const std::size_t middle = begin + (end - begin) / 2;
            std::nth_element(
                surface_bounds_entries.begin() + begin,
                surface_bounds_entries.begin() + middle,
                surface_bounds_entries.begin() + end,
                [axis](const IndexedSurface& left,
                       const IndexedSurface& right)
                {
                    const double left_center =
                        left.bounds.lower[axis] + left.bounds.upper[axis];
                    const double right_center =
                        right.bounds.lower[axis] + right.bounds.upper[axis];
                    if (left_center != right_center)
                        return left_center < right_center;
                    if (left.group != right.group)
                        return left.group < right.group;
                    return left.surface < right.surface;
                });
            const std::size_t first_child = build(begin, middle);
            const std::size_t second_child = build(middle, end);
            surface_bounds_nodes[node_index].first_child = first_child;
            surface_bounds_nodes[node_index].second_child = second_child;
            return node_index;
        };
        if (!surface_bounds_entries.empty())
            build(0, surface_bounds_entries.size());
        surface_bounds_delta.clear();
        ++profile.intrusion_index_rebuilds;
    };
    const auto indexGroupSurfaces = [&](const std::size_t group)
    {
        for (std::size_t surface = 0;
             surface < groups[group].intrusion_surfaces.size(); ++surface)
            surface_bounds_delta.push_back({
                group, surface, groups[group].version,
                groups[group].intrusion_surfaces[surface].bounds});
    };
    const auto eraseGroupSurfaces = [&](const std::size_t)
    {
    };
    const auto maybeRebuildSurfaceBoundsIndex = [&]
    {
        const std::size_t threshold = std::max<std::size_t>(
            1024, surface_bounds_entries.size() / 8);
        if (surface_bounds_delta.size() >= threshold)
            rebuildSurfaceBoundsIndex();
    };
    const auto querySurfaceBoundsIndex =
        [&](const ConvexHullSurface& query_hull, const auto& visit)
    {
        const Bounds& query = query_hull.bounds;
        if (!surface_bounds_nodes.empty())
        {
            std::vector<std::size_t> pending{0};
            while (!pending.empty())
            {
                const std::size_t node_index = pending.back();
                pending.pop_back();
                const SurfaceBoundsNode& node =
                    surface_bounds_nodes[node_index];
                if (!convexBoundsOverlap(
                        query, node.bounds, tolerance * 8.0))
                    continue;
                if (boundsEntirelyOutsideHull(
                        query_hull, node.bounds, tolerance * 8.0))
                {
                    ++profile.intrusion_index_halfspace_rejections;
                    continue;
                }
                if (node.first_child ==
                    std::numeric_limits<std::size_t>::max())
                {
                    for (std::size_t index = node.begin;
                         index < node.end; ++index)
                        if (convexBoundsOverlap(
                                query, surface_bounds_entries[index].bounds,
                                tolerance * 8.0))
                            visit(surface_bounds_entries[index]);
                    continue;
                }
                pending.push_back(node.first_child);
                pending.push_back(node.second_child);
            }
        }
        for (const IndexedSurface& indexed : surface_bounds_delta)
            if (convexBoundsOverlap(
                    query, indexed.bounds, tolerance * 8.0))
                visit(indexed);
    };
    rebuildSurfaceBoundsIndex();

    struct FaceVertexKey
    {
        std::uint64_t x = 0;
        std::uint64_t y = 0;
        std::uint64_t z = 0;
        bool operator==(const FaceVertexKey&) const = default;
    };
    struct ConvexFaceKey
    {
        std::array<FaceVertexKey, 3> vertices{};
        bool operator==(const ConvexFaceKey&) const = default;
    };
    struct ConvexFaceKeyHash
    {
        std::size_t operator()(const ConvexFaceKey& key) const noexcept
        {
            std::size_t seed = 0;
            for (const auto& vertex : key.vertices)
                for (const auto coordinate : {vertex.x, vertex.y, vertex.z})
                    seed ^= std::hash<std::uint64_t>{}(coordinate) +
                        0x9e3779b97f4a7c15ULL +
                        (seed << 6U) + (seed >> 2U);
            return seed;
        }
    };
    std::unordered_map<ConvexFaceKey, bool, ConvexFaceKeyHash>
        convex_face_error_cache;
    convex_face_error_cache.reserve(16384);

    struct ClosureStateKey
    {
        std::vector<std::uint64_t> members;
        bool operator==(const ClosureStateKey&) const = default;
    };
    struct ClosureStateKeyHash
    {
        std::size_t operator()(const ClosureStateKey& key) const noexcept
        {
            std::size_t seed = 0;
            for (const auto member : key.members)
                seed ^= std::hash<std::uint64_t>{}(member) +
                    0x9e3779b97f4a7c15ULL +
                    (seed << 6U) + (seed >> 2U);
            return seed;
        }
    };
    std::unordered_set<ClosureStateKey, ClosureStateKeyHash>
        rejected_closure_states;
    std::unordered_map<ClosureStateKey, std::vector<std::uint64_t>,
                       ClosureStateKeyHash> closure_transitions;
    std::size_t rejected_closure_members = 0;
    std::size_t transition_closure_members = 0;
    constexpr std::size_t maximum_cached_closure_members = 2000000;
    const auto closureStateKey = [&](const std::vector<std::size_t>& consumed)
    {
        ClosureStateKey key;
        key.members.reserve(consumed.size());
        for (const auto group : consumed)
            key.members.push_back(
                (groups[group].version << 32U) |
                static_cast<std::uint64_t>(group));
        std::sort(key.members.begin(), key.members.end());
        return key;
    };
    const auto cacheRejectedClosure = [&](ClosureStateKey key)
    {
        if (rejected_closure_members + key.members.size() >
            maximum_cached_closure_members)
        {
            rejected_closure_states.clear();
            rejected_closure_members = 0;
        }
        rejected_closure_members += key.members.size();
        rejected_closure_states.insert(std::move(key));
    };

    const std::function<std::optional<FitResult>(std::size_t, std::size_t)>
        fitPair = [&](std::size_t first, std::size_t second)
        -> std::optional<FitResult>
    {
        ++profile.candidate_evaluations;
        flushProfile();
        if (first > second) std::swap(first, second);
        const auto candidatePassesError =
            [&](const std::vector<OutputPrimitive>& candidate)
        {
            const auto distance_started = std::chrono::steady_clock::now();
            const FilledSurfaceDistanceCertificate certificate =
                certifyFilledSurfaceDistance(
                    distance_reference, candidate,
                    candidate_error_limit + tolerance);
            profile.distance_seconds += std::chrono::duration<double>(
                std::chrono::steady_clock::now() - distance_started).count();
            profile.certificate_samples += certificate.samples;
            if (certificate.exceeded) return false;
            if (certificate.certified) return true;
            ++profile.fine_distance_evaluations;
            const auto fine_started = std::chrono::steady_clock::now();
            const bool accepted = maximumFilledSurfaceDistance(
                distance_reference, candidate, error_sample_spacing,
                candidate_error_limit + tolerance) <=
                candidate_error_limit + tolerance;
            profile.distance_seconds += std::chrono::duration<double>(
                std::chrono::steady_clock::now() - fine_started).count();
            return accepted;
        };
        const auto convexFaceKey = [&](const OutputPrimitive& face)
            -> std::optional<ConvexFaceKey>
        {
            const PrimitiveMesh surface =
                triangulatePrimitive(face.primitive);
            if (surface.faces.size() != 1 || surface.vertices.size() != 3)
                return std::nullopt;
            ConvexFaceKey key;
            for (std::size_t index = 0; index < 3; ++index)
            {
                const Vec3& vertex = surface.vertices[index];
                const auto bits = [](const double value)
                {
                    return std::bit_cast<std::uint64_t>(
                        value == 0.0 ? 0.0 : value);
                };
                key.vertices[index] = {
                    bits(vertex.x()), bits(vertex.y()), bits(vertex.z())};
            }
            std::sort(key.vertices.begin(), key.vertices.end(),
                [](const FaceVertexKey& left, const FaceVertexKey& right)
                {
                    if (left.x != right.x) return left.x < right.x;
                    if (left.y != right.y) return left.y < right.y;
                    return left.z < right.z;
                });
            return key;
        };
        const auto convexCandidatePassesError =
            [&](const std::vector<OutputPrimitive>& candidate)
        {
            for (const auto& face : candidate)
            {
                const auto key = convexFaceKey(face);
                if (!key)
                {
                    if (!candidatePassesError({face})) return false;
                    continue;
                }
                const auto cached = convex_face_error_cache.find(*key);
                if (cached != convex_face_error_cache.end())
                {
                    ++profile.convex_face_cache_hits;
                    if (!cached->second) return false;
                    continue;
                }
                // A triangle's vertices and centroid are exact points on the
                // proposed proxy. Any one of them exceeding the directed
                // distance limit proves that the whole candidate fails. This
                // cheap rejection pass never accepts a face; passing faces
                // still go through the continuous certificate below.
                ++profile.coarse_distance_evaluations;
                const auto coarse_started = std::chrono::steady_clock::now();
                const bool coarse_accepted = maximumFilledSurfaceDistance(
                    distance_reference, {face},
                    std::numeric_limits<double>::infinity(),
                    candidate_error_limit + tolerance) <=
                    candidate_error_limit + tolerance;
                const double coarse_seconds = std::chrono::duration<double>(
                    std::chrono::steady_clock::now() - coarse_started).count();
                profile.coarse_distance_seconds += coarse_seconds;
                profile.distance_seconds += coarse_seconds;
                if (!coarse_accepted)
                {
                    ++profile.coarse_distance_rejections;
                    if (convex_face_error_cache.size() >= 500000)
                        convex_face_error_cache.clear();
                    convex_face_error_cache.emplace(*key, false);
                    return false;
                }
                const bool accepted = candidatePassesError({face});
                if (convex_face_error_cache.size() >= 500000)
                    convex_face_error_cache.clear();
                convex_face_error_cache.emplace(*key, accepted);
                if (!accepted) return false;
            }
            return true;
        };

        std::vector<std::size_t> consumed{first, second};
        std::vector<std::uint32_t> responsibility;
        std::vector<std::uint32_t> source_vertices;
        std::size_t input_triangle_workload = 0;
        for (const auto group : consumed)
        {
            responsibility.insert(
                responsibility.end(), groups[group].responsibility.begin(),
                groups[group].responsibility.end());
            source_vertices.insert(
                source_vertices.end(),
                groups[group].source_vertices.begin(),
                groups[group].source_vertices.end());
            input_triangle_workload += groups[group].triangle_workload;
        }
        std::sort(responsibility.begin(), responsibility.end());
        responsibility.erase(
            std::unique(responsibility.begin(), responsibility.end()),
            responsibility.end());
        std::sort(source_vertices.begin(), source_vertices.end());
        source_vertices.erase(
            std::unique(source_vertices.begin(), source_vertices.end()),
            source_vertices.end());

        // Non-coplanar stage-3 candidates must be certified analytic
        // surfaces. An arbitrary 3-D convex hull is a solid fallback, not a
        // recognized surface, and it creates the large overlapping shells
        // that are particularly harmful to PQSS traversal.
        if (const auto revolved = fitCertifiedRevolvedSurface(
                responsibility_mesh, responsibility, options,
                tolerance * 8.0))
        {
            ++profile.certified_round_candidates;
            std::vector<OutputPrimitive> round_candidate;
            appendRevolvedSurfacePatches(
                round_candidate, responsibility_mesh, *revolved,
                responsibility, tolerance * 8.0);
            refineRoundSurfaceSegments(
                round_candidate, candidate_error_limit);
            const std::size_t round_workload =
                triangulatedFaceCount(round_candidate);
            if (!round_candidate.empty() &&
                round_workload < input_triangle_workload &&
                candidatePassesError(round_candidate))
            {
                ++profile.accepted_round_groups;
                Bounds bounds;
                for (const auto& item : round_candidate)
                {
                    const Bounds item_bounds = primitiveBounds(item.primitive);
                    bounds.lower = bounds.lower.cwiseMin(item_bounds.lower);
                    bounds.upper = bounds.upper.cwiseMax(item_bounds.upper);
                }
                FitResult result;
                result.shell = std::move(round_candidate);
                result.responsibility = std::move(responsibility);
                result.source_vertices = std::move(source_vertices);
                result.bounds = bounds;
                result.contains_round_surface = true;
                result.triangle_workload = round_workload;
                result.consumed_groups = consumed;
                return result;
            }
        }

        if (groups[first].contains_round_surface ||
            groups[second].contains_round_surface)
        {
            ++profile.degenerate_rejections;
            return std::nullopt;
        }

        // Grow a convex candidate only from an adjacent pair. If the fitted
        // volume intrudes another active surface, absorb that owner and refit
        // before acceptance; emitting both shells would leave internal or
        // intersecting faces for PQSS to traverse.
        std::unordered_set<std::size_t> consumed_set(consumed.begin(),
                                                      consumed.end());
        std::optional<ConvexHullSurface> convex_candidate;
        for (;;)
        {
            ++profile.convex_growth_iterations;
            ClosureStateKey closure_key = closureStateKey(consumed);
            if (rejected_closure_states.contains(closure_key))
            {
                ++profile.convex_error_cache_hits;
                ++profile.error_rejections;
                return std::nullopt;
            }
            const auto cached_transition =
                closure_transitions.find(closure_key);
            if (cached_transition != closure_transitions.end())
            {
                bool valid = !cached_transition->second.empty();
                std::vector<std::size_t> additions;
                for (const auto token : cached_transition->second)
                {
                    const std::size_t group =
                        static_cast<std::uint32_t>(token);
                    const std::uint64_t version = token >> 32U;
                    if (group >= groups.size() || !groups[group].active ||
                        groups[group].version != version ||
                        consumed_set.contains(group))
                    {
                        valid = false;
                        break;
                    }
                    additions.push_back(group);
                }
                if (valid)
                {
                    ++profile.convex_transition_cache_hits;
                    for (const auto group : additions)
                    {
                        consumed_set.insert(group);
                        consumed.push_back(group);
                    }
                    continue;
                }
                transition_closure_members -=
                    cached_transition->first.members.size() +
                    cached_transition->second.size();
                closure_transitions.erase(cached_transition);
            }
            responsibility.clear();
            source_vertices.clear();
            input_triangle_workload = 0;
            std::vector<Vec3> points;
            for (const auto group : consumed)
            {
                responsibility.insert(
                    responsibility.end(), groups[group].responsibility.begin(),
                    groups[group].responsibility.end());
                source_vertices.insert(
                    source_vertices.end(),
                    groups[group].source_vertices.begin(),
                    groups[group].source_vertices.end());
                input_triangle_workload += groups[group].triangle_workload;
                for (const auto& surface : groups[group].intrusion_surfaces)
                    points.insert(points.end(), surface.mesh.vertices.begin(),
                                  surface.mesh.vertices.end());
            }
            std::sort(responsibility.begin(), responsibility.end());
            responsibility.erase(
                std::unique(responsibility.begin(), responsibility.end()),
                responsibility.end());
            std::sort(source_vertices.begin(), source_vertices.end());
            source_vertices.erase(
                std::unique(source_vertices.begin(), source_vertices.end()),
                source_vertices.end());
            for (const auto vertex : source_vertices)
                points.push_back(responsibility_mesh.vertices[vertex]);
            ++profile.convex_hull_fits;
            const auto hull_started = std::chrono::steady_clock::now();
            convex_candidate = fitConvexHullSurface(
                std::move(points), responsibility,
                next_component_enclosure_group, tolerance);
            profile.convex_hull_seconds += std::chrono::duration<double>(
                std::chrono::steady_clock::now() - hull_started).count();
            if (!convex_candidate)
            {
                cacheRejectedClosure(std::move(closure_key));
                ++profile.degenerate_rejections;
                return std::nullopt;
            }
            if (!convexCandidatePassesError(convex_candidate->shell))
            {
                cacheRejectedClosure(std::move(closure_key));
                ++profile.error_rejections;
                return std::nullopt;
            }

            bool expanded = false;
            std::vector<std::size_t> iteration_additions;
            const auto intrusion_scan_started =
                std::chrono::steady_clock::now();
            std::unordered_set<std::size_t> intruded_groups;
            querySurfaceBoundsIndex(*convex_candidate,
                [&](const IndexedSurface& indexed)
            {
                const std::size_t group = indexed.group;
                if (!groups[group].active ||
                    groups[group].version != indexed.version ||
                    consumed_set.contains(group) ||
                    intruded_groups.contains(group) ||
                    indexed.surface >= groups[group].intrusion_surfaces.size())
                    return;
                ++profile.intrusion_index_candidates;
                const auto& surface =
                    groups[group].intrusion_surfaces[indexed.surface];
                if (!convexBoundsOverlap(
                        convex_candidate->bounds, groups[group].bounds,
                        tolerance * 8.0))
                    return;
                ++profile.intrusion_tests;
                bool halfspace_rejected = false;
                bool surface_plane_rejected = false;
                bool planar_projection_rejected = false;
                const bool intrudes = meshIntrudesConvexHull(
                    *convex_candidate, surface.mesh, surface.bounds,
                    tolerance * 8.0, &halfspace_rejected,
                    surface.planar ? &surface.plane_normal : nullptr,
                    surface.plane_offset, &surface_plane_rejected,
                    surface.planar ? &surface.projection : nullptr,
                    &planar_projection_rejected);
                if (halfspace_rejected)
                    ++profile.intrusion_halfspace_rejections;
                if (surface_plane_rejected)
                    ++profile.intrusion_surface_plane_rejections;
                if (planar_projection_rejected)
                    ++profile.intrusion_planar_projection_rejections;
                if (!intrudes) return;
                intruded_groups.insert(group);
            });
            std::vector<std::size_t> ordered_intruded_groups(
                intruded_groups.begin(), intruded_groups.end());
            std::sort(ordered_intruded_groups.begin(),
                      ordered_intruded_groups.end(),
                [&](const std::size_t left, const std::size_t right)
                {
                    if (groups[left].bounds.lower.x() !=
                        groups[right].bounds.lower.x())
                        return groups[left].bounds.lower.x() <
                               groups[right].bounds.lower.x();
                    return left < right;
                });
            for (const auto group : ordered_intruded_groups)
            {
                consumed_set.insert(group);
                consumed.push_back(group);
                iteration_additions.push_back(group);
                expanded = true;
            }
            profile.intrusion_seconds += std::chrono::duration<double>(
                std::chrono::steady_clock::now() -
                intrusion_scan_started).count();
            if (!expanded &&
                triangulatedFaceCount(convex_candidate->shell) >
                    input_triangle_workload)
            {
                // A pair of non-coplanar surface patches often needs closure
                // faces before it can grow into a useful envelope. Do not
                // commit that temporary triangle increase. Continue the same
                // adjacent-growth sequence until the candidate repays its
                // triangle debt, or reject it if no adjacent continuation is
                // available. Every intermediate hull has already passed the
                // same directed-distance audit above.
                std::vector<std::size_t> frontier;
                for (const auto group : consumed)
                    for (const auto neighbor : groups[group].neighbors)
                        if (groups[neighbor].active &&
                            !consumed_set.contains(neighbor))
                            frontier.push_back(neighbor);
                std::sort(frontier.begin(), frontier.end());
                frontier.erase(std::unique(frontier.begin(), frontier.end()),
                               frontier.end());
                std::sort(frontier.begin(), frontier.end(),
                    [&](const std::size_t left, const std::size_t right)
                    {
                        if (groups[left].triangle_workload !=
                            groups[right].triangle_workload)
                            return groups[left].triangle_workload >
                                   groups[right].triangle_workload;
                        return left < right;
                    });
                if (frontier.empty())
                {
                    ++profile.convex_debt_rejections;
                    ++profile.convex_workload_rejections;
                    return std::nullopt;
                }
                consumed_set.insert(frontier.front());
                consumed.push_back(frontier.front());
                iteration_additions.push_back(frontier.front());
                ++profile.convex_debt_growth_steps;
                expanded = true;
            }
            if (expanded)
            {
                if (transition_closure_members + closure_key.members.size() +
                        iteration_additions.size() >
                    maximum_cached_closure_members)
                {
                    closure_transitions.clear();
                    transition_closure_members = 0;
                }
                std::vector<std::uint64_t> addition_tokens;
                addition_tokens.reserve(iteration_additions.size());
                for (const auto group : iteration_additions)
                    addition_tokens.push_back(
                        (groups[group].version << 32U) |
                        static_cast<std::uint64_t>(group));
                transition_closure_members +=
                    closure_key.members.size() + addition_tokens.size();
                closure_transitions.insert_or_assign(
                    std::move(closure_key), std::move(addition_tokens));
            }
            if (!expanded) break;
        }
        const std::size_t convex_workload =
            triangulatedFaceCount(convex_candidate->shell);
        FitResult convex_result;
        convex_result.shell = std::move(convex_candidate->shell);
        convex_result.responsibility = std::move(responsibility);
        convex_result.source_vertices = std::move(source_vertices);
        convex_result.bounds = convex_candidate->bounds;
        convex_result.contains_round_surface = false;
        convex_result.triangle_workload = convex_workload;
        convex_result.consumed_groups = std::move(consumed);
        return convex_result;

        ++profile.degenerate_rejections;
        return std::nullopt;
    };

    std::priority_queue<
        Candidate, std::vector<Candidate>, CandidateCompare> queue;
    struct PairVersionKey
    {
        std::size_t first = 0;
        std::size_t second = 0;
        std::uint64_t first_version = 0;
        std::uint64_t second_version = 0;
        bool operator==(const PairVersionKey&) const = default;
    };
    struct PairVersionKeyHash
    {
        std::size_t operator()(const PairVersionKey& key) const noexcept
        {
            std::size_t seed = std::hash<std::size_t>{}(key.first);
            for (const auto value : {
                     static_cast<std::uint64_t>(key.second),
                     key.first_version, key.second_version})
                seed ^= std::hash<std::uint64_t>{}(value) +
                    0x9e3779b97f4a7c15ULL +
                    (seed << 6U) + (seed >> 2U);
            return seed;
        }
    };
    const auto pairVersionKey = [&](std::size_t first, std::size_t second)
    {
        if (first > second) std::swap(first, second);
        return PairVersionKey{first, second, groups[first].version,
                              groups[second].version};
    };
    std::unordered_set<PairVersionKey, PairVersionKeyHash> scheduled_pairs;
    std::unordered_set<PairVersionKey, PairVersionKeyHash> attempted_pairs;
    const auto enqueueBestNeighbor = [&](const std::size_t group)
    {
        if (group >= groups.size() || !groups[group].active) return;
        std::optional<Candidate> best;
        for (const auto neighbor : groups[group].neighbors)
        {
            if (!groups[neighbor].active) continue;
            const PairVersionKey key = pairVersionKey(group, neighbor);
            if (scheduled_pairs.contains(key) || attempted_pairs.contains(key))
                continue;
            const std::size_t priority =
                groups[key.first].triangle_workload +
                groups[key.second].triangle_workload;
            Candidate candidate{key.first, key.second,
                                key.first_version, key.second_version,
                                priority};
            if (!best || CandidateCompare{}(*best, candidate))
                best = candidate;
        }
        if (!best) return;
        const PairVersionKey key{best->first, best->second,
                                 best->first_version,
                                 best->second_version};
        scheduled_pairs.insert(key);
        queue.push(*best);
    };
    for (std::size_t first = 0; first < groups.size(); ++first)
        enqueueBestNeighbor(first);

    while (!queue.empty())
    {
        const Candidate pair = queue.top();
        queue.pop();
        const PairVersionKey key{pair.first, pair.second,
                                 pair.first_version, pair.second_version};
        scheduled_pairs.erase(key);
        if (!groups[pair.first].active || !groups[pair.second].active ||
            groups[pair.first].version != pair.first_version ||
            groups[pair.second].version != pair.second_version ||
            !groups[pair.first].neighbors.contains(pair.second))
        {
            enqueueBestNeighbor(pair.first);
            enqueueBestNeighbor(pair.second);
            continue;
        }
        attempted_pairs.insert(key);
        auto candidate = fitPair(pair.first, pair.second);
        if (!candidate)
        {
            enqueueBestNeighbor(pair.first);
            enqueueBestNeighbor(pair.second);
            continue;
        }
        closure_transitions.clear();
        transition_closure_members = 0;
        const bool accepted_convex_enclosure = std::any_of(
            candidate->shell.begin(), candidate->shell.end(),
            [&](const OutputPrimitive& item)
            { return item.enclosure_group == next_component_enclosure_group; });

        std::vector<std::size_t> affected;
        for (const auto consumed : candidate->consumed_groups)
            affected.insert(affected.end(), groups[consumed].neighbors.begin(),
                            groups[consumed].neighbors.end());
        std::sort(affected.begin(), affected.end());
        affected.erase(std::unique(affected.begin(), affected.end()), affected.end());
        for (const auto consumed : candidate->consumed_groups)
            eraseGroupSurfaces(consumed);
        groups[pair.first].shell = std::move(candidate->shell);
        groups[pair.first].intrusion_surfaces.clear();
        groups[pair.first].intrusion_surfaces.reserve(
            groups[pair.first].shell.size());
        for (const auto& surface : groups[pair.first].shell)
            groups[pair.first].intrusion_surfaces.push_back(
                makeIntrusionSurface(surface));
        groups[pair.first].responsibility = std::move(candidate->responsibility);
        groups[pair.first].source_vertices =
            std::move(candidate->source_vertices);
        groups[pair.first].bounds = candidate->bounds;
        groups[pair.first].contains_round_surface =
            candidate->contains_round_surface;
        groups[pair.first].triangle_workload = candidate->triangle_workload;
        groups[pair.first].neighbors.clear();
        ++groups[pair.first].version;
        indexGroupSurfaces(pair.first);
        for (const auto consumed : candidate->consumed_groups)
            if (consumed != pair.first)
            {
                groups[consumed].active = false;
                ++groups[consumed].version;
            }
        maybeRebuildSurfaceBoundsIndex();
        for (const auto neighbor : affected)
        {
            for (const auto consumed : candidate->consumed_groups)
                groups[neighbor].neighbors.erase(consumed);
            if (neighbor == pair.first || !groups[neighbor].active) continue;
            connect(pair.first, neighbor);
        }
        for (const auto consumed : candidate->consumed_groups)
            if (consumed != pair.first) groups[consumed].neighbors.clear();
        merged_group_count += candidate->consumed_groups.size() - 1;
        ++profile.accepted_groups;
        if (accepted_convex_enclosure)
            ++next_component_enclosure_group;
        enqueueBestNeighbor(pair.first);
        for (const auto neighbor : groups[pair.first].neighbors)
            enqueueBestNeighbor(neighbor);
    }

    std::vector<OutputPrimitive> result;
    for (auto& group : groups)
        if (group.active)
            result.insert(result.end(),
                          std::make_move_iterator(group.shell.begin()),
                          std::make_move_iterator(group.shell.end()));
    std::ofstream stream(profile_path);
    stream << std::setprecision(17)
           << "{\"complete\":true"
           << ",\"strategy\":\"certified_component_envelope_fixed_point\""
           << ",\"input_primitives\":" << profile.input_primitives
           << ",\"initial_adjacencies\":" << profile.initial_adjacencies
           << ",\"candidate_evaluations\":" << profile.candidate_evaluations
           << ",\"cache_hits\":" << profile.cache_hits
           << ",\"degenerate_rejections\":" << profile.degenerate_rejections
           << ",\"error_rejections\":" << profile.error_rejections
           << ",\"certified_round_candidates\":"
           << profile.certified_round_candidates
           << ",\"accepted_round_groups\":"
           << profile.accepted_round_groups
           << ",\"conservative_round_component_attempts\":"
           << profile.conservative_round_component_attempts
           << ",\"accepted_conservative_round_components\":"
           << profile.accepted_conservative_round_components
           << ",\"convex_component_attempts\":"
           << profile.convex_component_attempts
           << ",\"accepted_convex_components\":"
           << profile.accepted_convex_components
           << ",\"accepted_convex_partition_components\":"
           << profile.accepted_convex_partition_components
           << ",\"invalid_convex_components\":"
           << profile.invalid_convex_components
           << ",\"convex_workload_rejections\":"
           << profile.convex_workload_rejections
           << ",\"convex_error_rejections\":"
           << profile.convex_error_rejections
           << ",\"convex_intersection_rejections\":"
           << profile.convex_intersection_rejections
           << ",\"pairwise_enabled\":"
           << "true"
           << ",\"accepted_groups\":" << profile.accepted_groups
           << ",\"coarse_distance_evaluations\":"
           << profile.coarse_distance_evaluations
           << ",\"coarse_distance_rejections\":"
           << profile.coarse_distance_rejections
           << ",\"convex_face_cache_hits\":"
           << profile.convex_face_cache_hits
           << ",\"convex_hull_fits\":" << profile.convex_hull_fits
           << ",\"convex_growth_iterations\":"
           << profile.convex_growth_iterations
           << ",\"convex_debt_growth_steps\":"
           << profile.convex_debt_growth_steps
           << ",\"convex_debt_rejections\":"
           << profile.convex_debt_rejections
           << ",\"convex_error_cache_hits\":"
           << profile.convex_error_cache_hits
           << ",\"convex_transition_cache_hits\":"
           << profile.convex_transition_cache_hits
           << ",\"intrusion_tests\":" << profile.intrusion_tests
           << ",\"intrusion_index_candidates\":"
           << profile.intrusion_index_candidates
           << ",\"intrusion_index_rebuilds\":"
           << profile.intrusion_index_rebuilds
           << ",\"intrusion_index_halfspace_rejections\":"
           << profile.intrusion_index_halfspace_rejections
           << ",\"intrusion_halfspace_rejections\":"
           << profile.intrusion_halfspace_rejections
           << ",\"intrusion_surface_plane_rejections\":"
           << profile.intrusion_surface_plane_rejections
           << ",\"intrusion_planar_projection_rejections\":"
           << profile.intrusion_planar_projection_rejections
           << ",\"certificate_samples\":" << profile.certificate_samples
           << ",\"fine_distance_evaluations\":"
           << profile.fine_distance_evaluations
           << ",\"coarse_distance_seconds\":"
           << profile.coarse_distance_seconds
           << ",\"convex_hull_seconds\":" << profile.convex_hull_seconds
           << ",\"intrusion_seconds\":" << profile.intrusion_seconds
           << ",\"distance_seconds\":" << profile.distance_seconds
           << ",\"output_primitives\":" << result.size()
           << ",\"output_triangles\":" << triangulatedFaceCount(result)
           << ",\"total_seconds\":" << std::chrono::duration<double>(
                std::chrono::steady_clock::now() - started).count() << "}\n";
    return result;
}

std::size_t protectLargeOpposingCavityWalls(
    std::vector<OutputPrimitive>& primitives,
    const double cavity_volume_limit,
    const double model_diagonal,
    const double tolerance)
{
    if (!std::isfinite(cavity_volume_limit)) return 0;
    struct Wall
    {
        std::size_t primitive = 0;
        Vec3 normal = Vec3::Zero();
        double distance = 0.0;
        Bounds2 projection;
    };
    std::vector<Wall> walls;
    for (std::size_t index = 0; index < primitives.size(); ++index)
    {
        OutputPrimitive& item = primitives[index];
        if (item.primitive.kind != Kind::Polygon ||
            item.primitive.polygon.size() < 3)
            continue;
        const double area = primitiveSurfaceArea(item.primitive, tolerance);
        if (area * model_diagonal <= cavity_volume_limit + tolerance)
            continue;
        Vec3 normal;
        if (!planarNormal(item.primitive, tolerance, normal)) continue;
        for (int axis = 0; axis < 3; ++axis)
        {
            if (std::abs(normal[axis]) <= 1.0e-12) continue;
            if (normal[axis] < 0.0) normal = -normal;
            break;
        }
        const Mat3 frame = orthonormalFrame(normal);
        Bounds2 projection;
        for (const Vec3& vertex : item.primitive.polygon)
        {
            projection.lower = projection.lower.cwiseMin({
                vertex.dot(frame.col(1)), vertex.dot(frame.col(2))});
            projection.upper = projection.upper.cwiseMax({
                vertex.dot(frame.col(1)), vertex.dot(frame.col(2))});
        }
        walls.push_back({index, normal,
                         normal.dot(item.primitive.polygon.front()),
                         projection});
    }

    std::vector<bool> protected_wall(primitives.size(), false);
    for (std::size_t first = 0; first < walls.size(); ++first)
        for (std::size_t second = first + 1; second < walls.size(); ++second)
        {
            const Wall& a = walls[first];
            const Wall& b = walls[second];
            if (a.normal.dot(b.normal) < 1.0 - 1.0e-8) continue;
            const double separation = std::abs(a.distance - b.distance);
            if (separation <= tolerance) continue;
            const Vec2 overlap =
                a.projection.upper.cwiseMin(b.projection.upper) -
                a.projection.lower.cwiseMax(b.projection.lower);
            if (overlap.x() <= tolerance || overlap.y() <= tolerance)
                continue;
            if (overlap.prod() * separation <=
                cavity_volume_limit + tolerance)
                continue;
            protected_wall[a.primitive] = true;
            protected_wall[b.primitive] = true;
        }

    std::size_t count = 0;
    for (std::size_t index = 0; index < primitives.size(); ++index)
        if (protected_wall[index])
        {
            primitives[index].preserves_cavity_opening = true;
            ++count;
        }
    return count;
}

std::vector<OutputPrimitive> clipParallelInternalSurfaceOcclusion(
    std::vector<OutputPrimitive> primitives,
    const Vec3& model_center,
    const double maximum_depth,
    const double tolerance,
    std::vector<OutputPrimitive>& coverage_certificates,
    ParallelOcclusionStats& stats)
{
    constexpr int clipper_precision = 8;
    // Occlusion decisions in one pass must see the same geometry.  If a
    // coverer is moved out of `primitives` before a later candidate is tested,
    // a face hidden by several coverers can incorrectly survive merely because
    // the first coverers were already consumed.  Keep only the immutable
    // primitive geometry as the decision snapshot; responsibility metadata is
    // still taken from the live candidate below.
    std::vector<Primitive> source_geometry;
    source_geometry.reserve(primitives.size());
    for (const auto& item : primitives)
        source_geometry.push_back(item.primitive);
    std::vector<OutputPrimitive> result;
    result.reserve(primitives.size());
    for (std::size_t candidate_id = 0; candidate_id < primitives.size();
         ++candidate_id)
    {
        OutputPrimitive& candidate = primitives[candidate_id];
        if (candidate.preserves_cavity_opening ||
            candidate.primitive.kind != Kind::Polygon ||
            candidate.primitive.polygon.size() < 3)
        {
            result.push_back(std::move(candidate));
            continue;
        }
        Vec3 normal;
        if (!planarNormal(candidate.primitive, tolerance, normal))
        {
            result.push_back(std::move(candidate));
            continue;
        }
        const Vec3 origin = candidate.primitive.polygon.front();
        // Orient the candidate normal toward the exterior of the complete
        // model.  A parallel coverer may remove a candidate only when it is
        // farther along this outward direction; an inner skin must never be
        // allowed to erase the actual outer shell.
        if (normal.dot(origin - model_center) < 0.0)
            normal = -normal;
        const double candidate_distance = normal.dot(origin);
        const Mat3 basis = orthonormalFrame(normal);
        Mat3 frame;
        frame.col(0) = basis.col(1);
        frame.col(1) = basis.col(2);
        frame.col(2) = normal;
        Clipper2Lib::PathD candidate_path;
        for (const Vec3& vertex : candidate.primitive.polygon)
        {
            const Vec3 local = frame.transposeMultiply(vertex - origin);
            candidate_path.emplace_back(local.x(), local.y());
        }
        if (Clipper2Lib::Area(candidate_path) < 0.0)
            std::reverse(candidate_path.begin(), candidate_path.end());
        const double candidate_area = std::abs(
            Clipper2Lib::Area(candidate_path));
        if (candidate_area <= tolerance * tolerance)
        {
            result.push_back(std::move(candidate));
            continue;
        }

        ++stats.candidate_count;
        Clipper2Lib::PathsD outward_coverers;
        for (std::size_t coverer_id = 0; coverer_id < primitives.size();
             ++coverer_id)
        {
            if (coverer_id == candidate_id ||
                source_geometry[coverer_id].kind != Kind::Polygon ||
                source_geometry[coverer_id].polygon.size() < 3)
                continue;
            Vec3 coverer_normal;
            if (!planarNormal(
                    source_geometry[coverer_id], tolerance,
                    coverer_normal))
                continue;
            if (std::abs(coverer_normal.dot(normal)) < 1.0 - 1.0e-8)
                continue;
            const double coverer_distance = normal.dot(
                source_geometry[coverer_id].polygon.front());
            const double signed_depth = coverer_distance - candidate_distance;
            if (signed_depth <= tolerance || signed_depth > maximum_depth)
                continue;
            Clipper2Lib::PathD path;
            for (const Vec3& vertex : source_geometry[coverer_id].polygon)
            {
                const Vec3 local = frame.transposeMultiply(vertex - origin);
                path.emplace_back(local.x(), local.y());
            }
            if (std::abs(Clipper2Lib::Area(path)) <= tolerance * tolerance)
                continue;
            if (Clipper2Lib::Area(path) < 0.0)
                std::reverse(path.begin(), path.end());
            outward_coverers.push_back(std::move(path));
        }
        if (outward_coverers.empty())
        {
            result.push_back(std::move(candidate));
            continue;
        }
        const auto cover_union = Clipper2Lib::Union(
            outward_coverers, Clipper2Lib::FillRule::NonZero,
            clipper_precision);
        const double covered_area = std::abs(Clipper2Lib::Area(
            Clipper2Lib::Intersect(
                Clipper2Lib::PathsD{candidate_path}, cover_union,
                Clipper2Lib::FillRule::NonZero, clipper_precision)));
        const double covered_ratio = covered_area / candidate_area;
        stats.candidate_ratios.emplace_back(candidate_id, covered_ratio);
        stats.maximum_covered_ratio = std::max(
            stats.maximum_covered_ratio, covered_ratio);
        // A high same-facing outward coverage ratio is the geometric signature
        // of an internal layer.  A lower threshold permits partial inner caps
        // such as a wall hidden by several outer strips, while the direction
        // and bounded depth reject opposing exterior faces.
        if (covered_ratio < 0.80)
        {
            result.push_back(std::move(candidate));
            continue;
        }
        ++stats.accepted_candidates;
        stats.accepted_ids.push_back(candidate_id);
        const auto remaining = Clipper2Lib::Difference(
            Clipper2Lib::PathsD{candidate_path}, cover_union,
            Clipper2Lib::FillRule::NonZero, clipper_precision);
        const double removed_area = std::max(
            candidate_area - std::abs(Clipper2Lib::Area(remaining)), 0.0);
        coverage_certificates.push_back(candidate);
        ++stats.clipped_primitives;
        stats.removed_area += removed_area;
        if (remaining.empty())
        {
            ++stats.removed_primitives;
            continue;
        }
        Clipper2Lib::PathsD output_paths = remaining;
        if (std::any_of(remaining.begin(), remaining.end(),
            [](const auto& path) { return !Clipper2Lib::IsPositive(path); }))
        {
            Clipper2Lib::PathsD triangles;
            if (Clipper2Lib::Triangulate(
                    remaining, clipper_precision, triangles, false) !=
                Clipper2Lib::TriangulateResult::success)
            {
                result.push_back(std::move(candidate));
                --stats.clipped_primitives;
                stats.removed_area -= removed_area;
                coverage_certificates.pop_back();
                continue;
            }
            output_paths = std::move(triangles);
        }
        for (const auto& path : output_paths)
        {
            if (std::abs(Clipper2Lib::Area(path)) <= tolerance * tolerance)
                continue;
            std::vector<Vec2> boundary;
            for (const auto& point : path)
                boundary.emplace_back(point.x, point.y);
            boundary = simplifyPolygon(std::move(boundary), tolerance);
            if (boundary.size() < 3 ||
                triangulatePolygon(boundary).size() + 2 != boundary.size())
                continue;
            OutputPrimitive fragment = candidate;
            fragment.primitive = polygonPrimitive(boundary, origin, frame);
            fragment.preserves_cavity_opening = true;
            result.push_back(std::move(fragment));
            ++stats.output_fragments;
        }
    }
    return result;
}

std::vector<OutputPrimitive> clipSupportContactOcclusion(
    std::vector<OutputPrimitive> primitives,
    const std::vector<RecognizedProtrusion>& protrusions,
    const double tolerance,
    std::vector<OutputPrimitive>& coverage_certificates,
    SupportContactClipStats& stats)
{
    constexpr int clipper_precision = 8;
    std::vector<OutputPrimitive> result;
    result.reserve(primitives.size());
    for (OutputPrimitive& candidate : primitives)
    {
        if (candidate.primitive.kind != Kind::Polygon ||
            candidate.primitive.polygon.size() < 3)
        {
            result.push_back(std::move(candidate));
            continue;
        }
        Vec3 normal;
        if (!planarNormal(candidate.primitive, tolerance, normal))
        {
            result.push_back(std::move(candidate));
            continue;
        }
        const Vec3 origin = candidate.primitive.polygon.front();
        const Mat3 basis = orthonormalFrame(normal);
        Mat3 frame;
        frame.col(0) = basis.col(1);
        frame.col(1) = basis.col(2);
        frame.col(2) = normal;
        Clipper2Lib::PathD candidate_path;
        for (const Vec3& vertex : candidate.primitive.polygon)
        {
            const Vec3 local = frame.transposeMultiply(vertex - origin);
            candidate_path.emplace_back(local.x(), local.y());
        }
        if (Clipper2Lib::Area(candidate_path) < 0.0)
            std::reverse(candidate_path.begin(), candidate_path.end());

        Clipper2Lib::PathsD contact_paths;
        for (const RecognizedProtrusion& protrusion : protrusions)
        {
            if (!protrusion.has_support_plane ||
                std::abs(normal.dot(protrusion.support_outward)) <
                    1.0 - 1.0e-8 ||
                std::abs((origin - protrusion.support_plane_point).dot(
                    protrusion.support_outward)) > tolerance * 16.0)
                continue;
            bool occupies_negative_side = false;
            bool occupies_positive_side = false;
            std::vector<Vec2> section = boxPlaneCrossSection(
                protrusion.box, origin, frame, tolerance * 8.0,
                occupies_negative_side, occupies_positive_side);
            if (section.size() < 3) continue;
            Clipper2Lib::PathD path;
            path.reserve(section.size());
            for (const Vec2& point : section)
                path.emplace_back(point.x(), point.y());
            if (Clipper2Lib::Area(path) < 0.0)
                std::reverse(path.begin(), path.end());
            contact_paths.push_back(std::move(path));
        }
        if (contact_paths.empty())
        {
            result.push_back(std::move(candidate));
            continue;
        }
        const auto contacts = Clipper2Lib::Union(
            contact_paths, Clipper2Lib::FillRule::NonZero,
            clipper_precision);
        const auto remaining = Clipper2Lib::Difference(
            Clipper2Lib::PathsD{candidate_path}, contacts,
            Clipper2Lib::FillRule::NonZero, clipper_precision);
        const double removed_area = std::max(
            std::abs(Clipper2Lib::Area(candidate_path)) -
                std::abs(Clipper2Lib::Area(remaining)),
            0.0);
        if (removed_area <= tolerance * tolerance)
        {
            result.push_back(std::move(candidate));
            continue;
        }

        // The uncut owner remains an audit-only certificate.  It is never
        // exported or inserted into the PQSS BVH; the query geometry contains
        // only the visible support fragments around the contact footprint.
        coverage_certificates.push_back(candidate);
        stats.removed_area += removed_area;
        ++stats.clipped_primitives;
        if (remaining.empty())
        {
            ++stats.removed_primitives;
            continue;
        }

        Clipper2Lib::PathsD output_paths = remaining;
        if (std::any_of(remaining.begin(), remaining.end(),
            [](const auto& path) { return !Clipper2Lib::IsPositive(path); }))
        {
            // Do not triangulate one large polygon around every contact hole.
            // Constrained triangulation tends to connect a hole corner to a
            // distant outer corner, producing long triangles with large RSS
            // bounds.  Slice at every contact-vertex ordinate instead.  Within
            // one open strip a contact crosses the strip boundary, so the
            // remaining pieces are local simple polygons rather than a polygon
            // with holes.
            double lower_x = std::numeric_limits<double>::infinity();
            double upper_x = -std::numeric_limits<double>::infinity();
            double lower_y = std::numeric_limits<double>::infinity();
            double upper_y = -std::numeric_limits<double>::infinity();
            for (const auto& point : candidate_path)
            {
                lower_x = std::min(lower_x, point.x);
                upper_x = std::max(upper_x, point.x);
                lower_y = std::min(lower_y, point.y);
                upper_y = std::max(upper_y, point.y);
            }
            std::vector<double> cuts{lower_y, upper_y};
            for (const auto& path : contacts)
                for (const auto& point : path)
                    if (point.y > lower_y + tolerance &&
                        point.y < upper_y - tolerance)
                        cuts.push_back(point.y);
            std::sort(cuts.begin(), cuts.end());
            cuts.erase(std::unique(cuts.begin(), cuts.end(),
                [&](const double first, const double second)
                { return std::abs(first - second) <= tolerance; }), cuts.end());
            output_paths.clear();
            const double margin = std::max(
                upper_x - lower_x, tolerance * 32.0);
            for (std::size_t strip = 0; strip + 1 < cuts.size(); ++strip)
            {
                if (cuts[strip + 1] - cuts[strip] <= tolerance) continue;
                Clipper2Lib::PathD slab{
                    {lower_x - margin, cuts[strip]},
                    {upper_x + margin, cuts[strip]},
                    {upper_x + margin, cuts[strip + 1]},
                    {lower_x - margin, cuts[strip + 1]},
                };
                const auto subject_piece = Clipper2Lib::Intersect(
                    Clipper2Lib::PathsD{candidate_path},
                    Clipper2Lib::PathsD{slab},
                    Clipper2Lib::FillRule::NonZero, clipper_precision);
                if (subject_piece.empty()) continue;
                auto strip_remaining = Clipper2Lib::Difference(
                    subject_piece, contacts, Clipper2Lib::FillRule::NonZero,
                    clipper_precision);
                if (std::any_of(strip_remaining.begin(), strip_remaining.end(),
                    [](const auto& path) { return !Clipper2Lib::IsPositive(path); }))
                {
                    Clipper2Lib::PathsD triangles;
                    if (Clipper2Lib::Triangulate(
                            strip_remaining, clipper_precision, triangles,
                            false) != Clipper2Lib::TriangulateResult::success)
                        throw std::runtime_error(
                            "failed to triangulate a sliced support contact subtraction");
                    strip_remaining = std::move(triangles);
                }
                output_paths.insert(
                    output_paths.end(),
                    std::make_move_iterator(strip_remaining.begin()),
                    std::make_move_iterator(strip_remaining.end()));
            }
        }
        for (const auto& path : output_paths)
        {
            if (std::abs(Clipper2Lib::Area(path)) <=
                tolerance * tolerance) continue;
            std::vector<Vec2> boundary;
            boundary.reserve(path.size());
            for (const auto& point : path)
                boundary.emplace_back(point.x, point.y);
            boundary = simplifyPolygon(std::move(boundary), tolerance);
            if (boundary.size() < 3 ||
                triangulatePolygon(boundary).size() + 2 != boundary.size())
                throw std::runtime_error(
                    "support contact subtraction produced an invalid polygon");
            OutputPrimitive fragment = candidate;
            fragment.primitive = polygonPrimitive(boundary, origin, frame);
            fragment.preserves_cavity_opening = true;
            result.push_back(std::move(fragment));
            ++stats.output_fragments;
        }
    }
    return result;
}

std::vector<BoxFit> selectActiveBoxCertificates(
    const std::vector<OutputPrimitive>& active_primitives,
    const std::vector<BoxFit>& historical_boxes,
    const double tolerance)
{
    std::unordered_map<std::uint64_t, std::vector<Vec3>> group_vertices;
    for (const auto& item : active_primitives)
    {
        if (item.enclosure_group == 0) continue;
        const PrimitiveMesh mesh = triangulatePrimitive(item.primitive);
        auto& vertices = group_vertices[item.enclosure_group];
        vertices.insert(vertices.end(), mesh.vertices.begin(), mesh.vertices.end());
    }

    std::vector<BoxFit> active;
    active.reserve(group_vertices.size());
    const double match_tolerance = tolerance * 32.0;
    for (const auto& [group, vertices] : group_vertices)
    {
        (void)group;
        if (vertices.size() < 8) continue;
        for (const BoxFit& box : historical_boxes)
        {
            const bool vertices_inside = std::all_of(
                vertices.begin(), vertices.end(), [&](const Vec3& vertex)
                {
                    const Vec3 local = box.axes.transposeMultiply(
                        vertex - box.center);
                    return std::abs(local.x()) <=
                               box.half_size.x() + match_tolerance &&
                           std::abs(local.y()) <=
                               box.half_size.y() + match_tolerance &&
                           std::abs(local.z()) <=
                               box.half_size.z() + match_tolerance;
                });
            if (!vertices_inside) continue;

            bool all_corners_present = true;
            for (int corner = 0; corner < 8 && all_corners_present; ++corner)
            {
                Vec3 position = box.center;
                for (int axis = 0; axis < 3; ++axis)
                    position += box.axes.col(axis) *
                        ((corner & (1 << axis)) ? box.half_size[axis]
                                                : -box.half_size[axis]);
                all_corners_present = std::any_of(
                    vertices.begin(), vertices.end(), [&](const Vec3& vertex)
                    { return (vertex - position).norm() <= match_tolerance; });
            }
            if (!all_corners_present) continue;
            active.push_back(box);
            break;
        }
    }
    return active;
}

std::vector<std::uint32_t> removeContainedPrimitives(
    std::vector<OutputPrimitive>& primitives,
    const double tolerance,
    std::size_t& removed_count)
{
    std::vector<Bounds> bounds;
    std::vector<PrimitiveMesh> meshes;
    bounds.reserve(primitives.size());
    meshes.reserve(primitives.size());
    for (const auto& item : primitives)
    {
        bounds.push_back(primitiveBounds(item.primitive));
        meshes.push_back(triangulatePrimitive(item.primitive));
    }
    std::vector<bool> removed(primitives.size(), false);
    std::vector<std::uint32_t> excluded;
    for (std::size_t candidate = 0; candidate < primitives.size(); ++candidate)
    {
        if (primitives[candidate].preserves_cavity_opening) continue;
        const Vec3 candidate_extent = bounds[candidate].upper - bounds[candidate].lower;
        for (std::size_t coverer = 0; coverer < primitives.size(); ++coverer)
        {
            if (candidate == coverer || removed[coverer] ||
                !boundsContain(bounds[coverer], bounds[candidate], tolerance)) continue;
            const Vec3 coverer_extent = bounds[coverer].upper - bounds[coverer].lower;
            const double candidate_scale = candidate_extent.norm();
            const double coverer_scale = coverer_extent.norm();
            if (coverer_scale + tolerance < candidate_scale ||
                (std::abs(coverer_scale - candidate_scale) <= tolerance && coverer > candidate))
                continue;
            const bool contained = std::all_of(
                meshes[candidate].vertices.begin(), meshes[candidate].vertices.end(),
                [&](const Vec3& vertex)
                {
                    return containsPointCached(
                        primitives[coverer].primitive, meshes[coverer],
                        vertex, tolerance);
                });
            if (!contained) continue;
            removed[candidate] = true;
            excluded.insert(excluded.end(), primitives[candidate].source_faces.begin(),
                            primitives[candidate].source_faces.end());
            break;
        }
    }
    std::vector<OutputPrimitive> kept;
    kept.reserve(primitives.size());
    for (std::size_t index = 0; index < primitives.size(); ++index)
        if (!removed[index]) kept.push_back(std::move(primitives[index]));
    removed_count = primitives.size() - kept.size();
    primitives = std::move(kept);
    std::sort(excluded.begin(), excluded.end());
    excluded.erase(std::unique(excluded.begin(), excluded.end()), excluded.end());
    return excluded;
}

std::uint64_t enclosureGroupForFaces(
    const std::vector<std::uint32_t>& faces)
{
    std::vector<std::uint32_t> sorted = faces;
    std::sort(sorted.begin(), sorted.end());
    sorted.erase(std::unique(sorted.begin(), sorted.end()), sorted.end());
    std::uint64_t hash = 1469598103934665603ull;
    for (const auto face : sorted)
    {
        hash ^= static_cast<std::uint64_t>(face) + 1;
        hash *= 1099511628211ull;
    }
    return hash == 0 ? 1 : hash;
}

std::vector<OutputPrimitive> removePrimitivesInsideEnclosureGroups(
    const Mesh& source_mesh,
    std::vector<OutputPrimitive> primitives,
    const double tolerance,
    std::size_t& removed_count)
{
    struct Plane
    {
        Vec3 point;
        Vec3 normal;
    };
    struct Enclosure
    {
        std::uint64_t group = 0;
        std::vector<Plane> planes;
        Bounds bounds;
    };

    std::unordered_map<std::uint64_t, std::vector<std::size_t>> groups;
    for (std::size_t index = 0; index < primitives.size(); ++index)
        if (primitives[index].enclosure_group != 0)
            groups[primitives[index].enclosure_group].push_back(index);

    std::vector<Enclosure> enclosures;
    for (const auto& [group, members] : groups)
    {
        std::vector<PrimitiveMesh> meshes;
        std::vector<Vec3> vertices;
        for (const auto member : members)
        {
            meshes.push_back(triangulatePrimitive(primitives[member].primitive));
            vertices.insert(vertices.end(), meshes.back().vertices.begin(),
                            meshes.back().vertices.end());
        }
        if (vertices.size() < 4) continue;
        Vec3 center = Vec3::Zero();
        for (const Vec3& vertex : vertices) center += vertex;
        center /= static_cast<double>(vertices.size());

        Enclosure enclosure;
        enclosure.group = group;
        for (const Vec3& vertex : vertices)
        {
            enclosure.bounds.lower = enclosure.bounds.lower.cwiseMin(vertex);
            enclosure.bounds.upper = enclosure.bounds.upper.cwiseMax(vertex);
        }
        for (const auto& mesh : meshes)
            for (const Face& face : mesh.faces)
            {
                const Vec3& first = mesh.vertices[face[0]];
                Vec3 normal = (mesh.vertices[face[1]] - first).cross(
                    mesh.vertices[face[2]] - first);
                const double length = normal.norm();
                if (length <= 1.0e-30) continue;
                normal /= length;
                if (normal.dot(center - first) > 0.0) normal = -normal;
                enclosure.planes.push_back({first, normal});
            }
        if (enclosure.planes.size() < 4) continue;
        const bool convex = std::all_of(vertices.begin(), vertices.end(),
            [&](const Vec3& vertex)
            {
                return std::all_of(enclosure.planes.begin(), enclosure.planes.end(),
                    [&](const Plane& plane)
                    {
                        return plane.normal.dot(vertex - plane.point) <= tolerance;
                    });
            });
        if (convex) enclosures.push_back(std::move(enclosure));
    }

    std::vector<bool> removed(primitives.size(), false);
    std::vector<std::uint64_t> covering_group(primitives.size(), 0);
    std::unordered_set<std::uint64_t> removed_groups;
    std::unordered_map<std::uint64_t, std::vector<std::uint32_t>>
        transferred_faces;

    // Remove a complete enclosing candidate when all source triangles for
    // which it is responsible are already contained by geometrically
    // dominating enclosure groups. Doing this before per-face removal avoids
    // leaving four side walls behind after the two coplanar caps are merged
    // into larger surfaces.
    std::unordered_map<std::uint64_t, const Enclosure*> enclosure_by_group;
    std::unordered_map<std::uint64_t, double> enclosure_order_volume;
    for (const Enclosure& enclosure : enclosures)
    {
        enclosure_by_group[enclosure.group] = &enclosure;
        enclosure_order_volume[enclosure.group] =
            (enclosure.bounds.upper - enclosure.bounds.lower)
                .cwiseMax(Vec3::Zero()).prod();
    }
    std::vector<std::uint64_t> group_order;
    group_order.reserve(groups.size());
    for (const auto& [group, members] : groups)
        if (enclosure_by_group.contains(group) && !members.empty())
            group_order.push_back(group);
    std::sort(group_order.begin(), group_order.end(),
        [&](const std::uint64_t first, const std::uint64_t second)
        {
            if (enclosure_order_volume[first] != enclosure_order_volume[second])
                return enclosure_order_volume[first] <
                       enclosure_order_volume[second];
            return first > second;
        });
    const auto enclosureContainsFace = [&](const Enclosure& enclosure,
                                           const std::uint32_t face_id)
    {
        const Face& face = source_mesh.faces[face_id];
        Bounds triangle_bounds;
        for (const auto vertex_id : face)
        {
            const Vec3& vertex = source_mesh.vertices[vertex_id];
            triangle_bounds.lower = triangle_bounds.lower.cwiseMin(vertex);
            triangle_bounds.upper = triangle_bounds.upper.cwiseMax(vertex);
        }
        if (!boundsContain(enclosure.bounds, triangle_bounds, tolerance))
            return false;
        return std::all_of(face.begin(), face.end(), [&](const auto vertex_id)
        {
            const Vec3& vertex = source_mesh.vertices[vertex_id];
            return std::all_of(
                enclosure.planes.begin(), enclosure.planes.end(),
                [&](const Plane& plane)
                {
                    return plane.normal.dot(vertex - plane.point) <= tolerance;
                });
        });
    };
    for (const std::uint64_t candidate_group : group_order)
    {
        if (removed_groups.contains(candidate_group)) continue;
        std::vector<std::uint32_t> responsibility;
        for (const auto member : groups[candidate_group])
            responsibility.insert(
                responsibility.end(), primitives[member].source_faces.begin(),
                primitives[member].source_faces.end());
        if (const auto transferred = transferred_faces.find(candidate_group);
            transferred != transferred_faces.end())
            responsibility.insert(
                responsibility.end(), transferred->second.begin(),
                transferred->second.end());
        std::sort(responsibility.begin(), responsibility.end());
        responsibility.erase(
            std::unique(responsibility.begin(), responsibility.end()),
            responsibility.end());
        if (responsibility.empty()) continue;

        std::unordered_map<std::uint64_t, std::vector<std::uint32_t>>
            candidate_transfers;
        bool fully_covered = true;
        for (const auto face_id : responsibility)
        {
            std::uint64_t coverer_group = 0;
            for (const Enclosure& coverer : enclosures)
            {
                if (coverer.group == candidate_group ||
                    removed_groups.contains(coverer.group))
                    continue;
                const double candidate_volume =
                    enclosure_order_volume[candidate_group];
                const double coverer_volume =
                    enclosure_order_volume[coverer.group];
                if (coverer_volume < candidate_volume ||
                    (coverer_volume == candidate_volume &&
                     coverer.group > candidate_group))
                    continue;
                if (enclosureContainsFace(coverer, face_id))
                {
                    coverer_group = coverer.group;
                    break;
                }
            }
            if (coverer_group == 0)
            {
                fully_covered = false;
                break;
            }
            candidate_transfers[coverer_group].push_back(face_id);
        }
        if (!fully_covered) continue;
        removed_groups.insert(candidate_group);
        for (const auto member : groups[candidate_group]) removed[member] = true;
        for (auto& [coverer, faces] : candidate_transfers)
        {
            auto& target = transferred_faces[coverer];
            target.insert(target.end(), faces.begin(), faces.end());
        }
    }

    std::vector<Bounds> bounds;
    std::vector<PrimitiveMesh> meshes;
    bounds.reserve(primitives.size());
    meshes.reserve(primitives.size());
    for (const auto& item : primitives)
    {
        bounds.push_back(primitiveBounds(item.primitive));
        meshes.push_back(triangulatePrimitive(item.primitive));
    }
    for (std::size_t candidate = 0; candidate < primitives.size(); ++candidate)
    {
        if (removed[candidate] ||
            primitives[candidate].preserves_cavity_opening)
            continue;
        for (const Enclosure& enclosure : enclosures)
        {
            if (removed_groups.contains(enclosure.group)) continue;
            if (primitives[candidate].enclosure_group == enclosure.group ||
                !boundsContain(enclosure.bounds, bounds[candidate], tolerance))
                continue;
            const bool inside = std::all_of(
                meshes[candidate].vertices.begin(), meshes[candidate].vertices.end(),
                [&](const Vec3& vertex)
                {
                    return std::all_of(enclosure.planes.begin(), enclosure.planes.end(),
                        [&](const Plane& plane)
                        {
                            return plane.normal.dot(vertex - plane.point) <= tolerance;
                        });
                });
            if (inside)
            {
                removed[candidate] = true;
                covering_group[candidate] = enclosure.group;
                break;
            }
        }
    }

    // Geometric containment proves that the enclosing group covers the removed
    // primitive. Preserve that responsibility explicitly so the final export
    // audit and the region viewer do not lose the removed source faces.
    for (std::size_t index = 0; index < primitives.size(); ++index)
        if (removed[index] && covering_group[index] != 0)
        {
            auto& faces = transferred_faces[covering_group[index]];
            faces.insert(faces.end(), primitives[index].source_faces.begin(),
                         primitives[index].source_faces.end());
        }
    for (auto& [group, faces] : transferred_faces)
    {
        std::sort(faces.begin(), faces.end());
        faces.erase(std::unique(faces.begin(), faces.end()), faces.end());
        for (std::size_t index = 0; index < primitives.size(); ++index)
            if (!removed[index] &&
                primitives[index].enclosure_group == group)
            {
                auto& item = primitives[index];
                item.source_faces.insert(item.source_faces.end(),
                                         faces.begin(), faces.end());
                std::sort(item.source_faces.begin(), item.source_faces.end());
                item.source_faces.erase(
                    std::unique(item.source_faces.begin(), item.source_faces.end()),
                    item.source_faces.end());
            }
    }

    std::vector<OutputPrimitive> kept;
    kept.reserve(primitives.size());
    for (std::size_t index = 0; index < primitives.size(); ++index)
        if (!removed[index]) kept.push_back(std::move(primitives[index]));
    removed_count = primitives.size() - kept.size();
    return kept;
}

struct FinalCoverageAudit
{
    std::size_t assigned_source_faces = 0;
    std::size_t enclosure_source_faces = 0;
    std::size_t planar_source_faces = 0;
    std::size_t unassigned_source_faces = 0;
    std::size_t failed_source_faces = 0;
    std::size_t failed_with_enclosure_owner = 0;
    std::size_t failed_with_planar_owner = 0;
    std::size_t failed_with_certificate_planar_owner = 0;
    std::size_t failed_with_certificate_band_owner = 0;
    std::vector<std::uint32_t> failed_face_ids;
};

FinalCoverageAudit auditFinalConservativeCoverage(
    const Mesh& source_mesh,
    const std::vector<OutputPrimitive>& primitives,
    const double tolerance,
    const std::vector<OutputPrimitive>* enclosure_certificate_primitives = nullptr,
    const std::vector<CertifiedExtrusion>* extrusion_certificates = nullptr,
    const std::vector<std::uint32_t>* excluded_redundant_faces = nullptr)
{
    struct Plane
    {
        Vec3 point;
        Vec3 normal;
    };
    struct Enclosure
    {
        Bounds bounds;
        std::vector<Plane> planes;
    };

    // Coplanar canonicalization preserves the exact surface union, but a merged
    // face may combine surfaces from different enclosing candidates and thus
    // cannot retain one enclosure_group id. Build convex-volume certificates
    // from the selected pre-canonicalization candidates while checking planar
    // coverage against the actual final primitives below.
    const auto& certificate_primitives = enclosure_certificate_primitives
        ? *enclosure_certificate_primitives : primitives;
    std::unordered_map<std::uint64_t, std::vector<std::size_t>> groups;
    for (std::size_t index = 0; index < certificate_primitives.size(); ++index)
        if (certificate_primitives[index].enclosure_group != 0)
            groups[certificate_primitives[index].enclosure_group].push_back(index);

    std::vector<Enclosure> enclosures;
    enclosures.reserve(groups.size());
    for (const auto& [group, members] : groups)
    {
        (void)group;
        std::vector<PrimitiveMesh> meshes;
        std::vector<Vec3> vertices;
        for (const auto member : members)
        {
            meshes.push_back(triangulatePrimitive(
                certificate_primitives[member].primitive));
            vertices.insert(vertices.end(), meshes.back().vertices.begin(),
                            meshes.back().vertices.end());
        }
        if (vertices.size() < 4) continue;
        Vec3 center = Vec3::Zero();
        for (const Vec3& vertex : vertices) center += vertex;
        center /= static_cast<double>(vertices.size());

        Enclosure enclosure;
        for (const Vec3& vertex : vertices)
        {
            enclosure.bounds.lower = enclosure.bounds.lower.cwiseMin(vertex);
            enclosure.bounds.upper = enclosure.bounds.upper.cwiseMax(vertex);
        }
        for (const PrimitiveMesh& mesh : meshes)
            for (const Face& face : mesh.faces)
            {
                const Vec3& first = mesh.vertices[face[0]];
                Vec3 normal = (mesh.vertices[face[1]] - first).cross(
                    mesh.vertices[face[2]] - first);
                const double length = normal.norm();
                if (length <= 1.0e-30) continue;
                normal /= length;
                if (normal.dot(center - first) > 0.0) normal = -normal;
                enclosure.planes.push_back({first, normal});
            }
        if (enclosure.planes.size() < 4) continue;
        const bool convex = std::all_of(vertices.begin(), vertices.end(),
            [&](const Vec3& vertex)
            {
                return std::all_of(enclosure.planes.begin(), enclosure.planes.end(),
                    [&](const Plane& plane)
                    {
                        return plane.normal.dot(vertex - plane.point) <= tolerance;
                    });
            });
        if (convex) enclosures.push_back(std::move(enclosure));
    }

    if (extrusion_certificates != nullptr)
        for (const auto& extrusion : *extrusion_certificates)
        {
            const auto hull = convexHull(extrusion.boundary, tolerance);
            if (hull.size() < 3) continue;
            const double boundary_area = std::abs(signedArea(extrusion.boundary));
            const double hull_area = std::abs(signedArea(hull));
            if (hull_area - boundary_area > tolerance * tolerance * 16.0)
                continue;

            Enclosure enclosure;
            std::vector<Vec3> vertices;
            vertices.reserve(2 * hull.size());
            for (const double depth : {extrusion.lower_distance,
                                       extrusion.upper_distance})
                for (const Vec2& point : hull)
                {
                    const Vec3 vertex = extrusion.origin +
                        extrusion.frame.col(0) * point.x() +
                        extrusion.frame.col(1) * point.y() +
                        extrusion.frame.col(2) * depth;
                    vertices.push_back(vertex);
                    enclosure.bounds.lower = enclosure.bounds.lower.cwiseMin(vertex);
                    enclosure.bounds.upper = enclosure.bounds.upper.cwiseMax(vertex);
                }
            Vec3 center = Vec3::Zero();
            for (const Vec3& vertex : vertices) center += vertex;
            center /= static_cast<double>(vertices.size());
            const auto addPlane = [&](const Vec3& point, Vec3 normal)
            {
                const double length = normal.norm();
                if (length <= 1.0e-30) return;
                normal /= length;
                if (normal.dot(center - point) > 0.0) normal = -normal;
                enclosure.planes.push_back({point, normal});
            };
            addPlane(vertices.front(), -extrusion.frame.col(2));
            addPlane(vertices[hull.size()], extrusion.frame.col(2));
            for (std::size_t edge = 0; edge < hull.size(); ++edge)
            {
                const Vec3 first = vertices[edge];
                const Vec3 second = vertices[(edge + 1) % hull.size()];
                addPlane(first, (second - first).cross(extrusion.frame.col(2)));
            }
            if (enclosure.planes.size() >= 4)
                enclosures.push_back(std::move(enclosure));
        }

    std::vector<std::vector<std::size_t>> owners(source_mesh.faces.size());
    for (std::size_t primitive = 0; primitive < primitives.size(); ++primitive)
        for (const auto face : primitives[primitive].source_faces)
            if (face < owners.size()) owners[face].push_back(primitive);
    std::vector<std::uint8_t> responsibility_certified(
        source_mesh.faces.size(), 0);
    std::vector<std::vector<std::size_t>> certificate_owners(
        source_mesh.faces.size());
    for (std::size_t certificate_index = 0;
         certificate_index < certificate_primitives.size();
         ++certificate_index)
        for (const auto face :
             certificate_primitives[certificate_index].source_faces)
            if (face < responsibility_certified.size())
            {
                responsibility_certified[face] = 1;
                certificate_owners[face].push_back(certificate_index);
            }
    std::vector<std::uint8_t> explicitly_excluded(
        source_mesh.faces.size(), 0);
    if (excluded_redundant_faces != nullptr)
        for (const auto face : *excluded_redundant_faces)
            if (face < explicitly_excluded.size()) explicitly_excluded[face] = 1;

    constexpr int clipper_precision = 8;
    FinalCoverageAudit audit;
    std::vector<std::size_t> all_primitives(primitives.size());
    std::iota(all_primitives.begin(), all_primitives.end(), 0);
    for (std::size_t face_id = 0; face_id < source_mesh.faces.size(); ++face_id)
    {
        if (explicitly_excluded[face_id])
        {
            ++audit.assigned_source_faces;
            continue;
        }
        const Face& face = source_mesh.faces[face_id];
        const std::array<Vec3, 3> triangle{{
            source_mesh.vertices[face[0]], source_mesh.vertices[face[1]],
            source_mesh.vertices[face[2]]}};
        if (owners[face_id].empty() && !responsibility_certified[face_id])
        {
            ++audit.unassigned_source_faces;
            audit.failed_face_ids.push_back(static_cast<std::uint32_t>(face_id));
            continue;
        }
        ++audit.assigned_source_faces;

        Bounds triangle_bounds;
        for (const Vec3& vertex : triangle)
        {
            triangle_bounds.lower = triangle_bounds.lower.cwiseMin(vertex);
            triangle_bounds.upper = triangle_bounds.upper.cwiseMax(vertex);
        }
        bool covered = false;
        for (const Enclosure& enclosure : enclosures)
        {
            if (!boundsContain(enclosure.bounds, triangle_bounds, tolerance)) continue;
            covered = std::all_of(triangle.begin(), triangle.end(),
                [&](const Vec3& vertex)
                {
                    return std::all_of(enclosure.planes.begin(), enclosure.planes.end(),
                        [&](const Plane& plane)
                        {
                            return plane.normal.dot(vertex - plane.point) <= tolerance;
                        });
                });
            if (covered) break;
        }
        if (covered)
        {
            ++audit.enclosure_source_faces;
            continue;
        }

        // A trimmed analytic band is an outward, circumscribed radial
        // enclosure over an explicitly certified angular/axial domain.  Its
        // source triangle need not be coplanar with the proxy facets, so audit
        // it in the band's parameter space before the planar-union fallback.
        covered = std::any_of(
            owners[face_id].begin(), owners[face_id].end(),
            [&](const std::size_t owner)
            {
                return trimmedAnalyticBandContainsTriangle(
                    primitives[owner].primitive, triangle, tolerance);
            });
        if (!covered)
            covered = std::any_of(
                certificate_owners[face_id].begin(),
                certificate_owners[face_id].end(),
                [&](const std::size_t owner)
                {
                    return trimmedAnalyticBandContainsTriangle(
                        certificate_primitives[owner].primitive,
                        triangle, tolerance);
                });
        if (covered)
        {
            ++audit.enclosure_source_faces;
            continue;
        }

        // Approximate planar recognition places disks, annuli, and polygons on
        // an outward support plane.  They need not be exactly coplanar with a
        // quantized source triangle, so certify the one-sided projection before
        // the exact coplanar-union path below.
        covered = std::any_of(
            owners[face_id].begin(), owners[face_id].end(),
            [&](const std::size_t owner)
            {
                return outwardPlanarSupportContainsTriangle(
                    primitives[owner].primitive, triangle, tolerance);
            });
        if (!covered)
            covered = std::any_of(
                certificate_owners[face_id].begin(),
                certificate_owners[face_id].end(),
                [&](const std::size_t owner)
                {
                    return outwardPlanarSupportContainsTriangle(
                        certificate_primitives[owner].primitive,
                        triangle, tolerance);
                });
        if (covered)
        {
            ++audit.planar_source_faces;
            continue;
        }

        // Exact local repair triangles are a stronger certificate than the
        // projected polygon-union fallback below. At large CAD coordinates,
        // Clipper can report a microscopic remainder even when the retained
        // triangle has the same three vertices as the source face.
        for (const auto owner : owners[face_id])
        {
            const OutputPrimitive& item = primitives[owner];
            if (item.source_faces.size() != 1 ||
                item.source_faces.front() != face_id)
                continue;
            const PrimitiveMesh proxy = triangulatePrimitive(item.primitive);
            for (const Face& proxy_face : proxy.faces)
            {
                const bool same_vertices = std::all_of(
                    triangle.begin(), triangle.end(), [&](const Vec3& vertex)
                    {
                        return std::any_of(
                            proxy_face.begin(), proxy_face.end(),
                            [&](const auto proxy_vertex)
                            {
                                return (proxy.vertices[proxy_vertex] - vertex)
                                    .norm() <= tolerance;
                            });
                    });
                if (same_vertices)
                {
                    covered = true;
                    break;
                }
            }
            if (covered) break;
        }
        if (covered)
        {
            ++audit.planar_source_faces;
            continue;
        }

        Vec3 normal = (triangle[1] - triangle[0]).cross(triangle[2] - triangle[0]);
        const double doubled_area = normal.norm();
        if (doubled_area <= tolerance * tolerance)
        {
            ++audit.failed_source_faces;
            continue;
        }
        normal /= doubled_area;
        const Mat3 basis = orthonormalFrame(normal);
        Mat3 frame;
        frame.col(0) = static_cast<Vec3>(basis.col(1));
        frame.col(1) = static_cast<Vec3>(basis.col(2));
        frame.col(2) = normal;
        Clipper2Lib::PathD source_path;
        for (const Vec3& vertex : triangle)
        {
            const Vec3 local = frame.transposeMultiply(vertex - triangle[0]);
            source_path.emplace_back(local.x(), local.y());
        }
        if (Clipper2Lib::Area(source_path) < 0.0)
            std::reverse(source_path.begin(), source_path.end());

        const auto coveredByPlanarPrimitives = [&](
            const std::vector<std::size_t>& candidates)
        {
            Clipper2Lib::PathsD cover_paths;
            for (const std::size_t index : candidates)
            {
                Vec3 candidate_normal;
                if (!planarNormal(primitives[index].primitive, tolerance,
                                  candidate_normal) ||
                    std::abs(candidate_normal.dot(normal)) < 1.0 - 1.0e-8)
                    continue;
                if (std::abs((planarPoint(primitives[index].primitive) -
                              triangle[0]).dot(normal)) > tolerance)
                    continue;
                const PrimitiveMesh mesh = triangulatePrimitive(
                    primitives[index].primitive);
                for (const Face& proxy_face : mesh.faces)
                {
                    Clipper2Lib::PathD path;
                    for (const auto vertex : proxy_face)
                    {
                        const Vec3 local = frame.transposeMultiply(
                            mesh.vertices[vertex] - triangle[0]);
                        path.emplace_back(local.x(), local.y());
                    }
                    if (std::abs(Clipper2Lib::Area(path)) <=
                        tolerance * tolerance) continue;
                    if (Clipper2Lib::Area(path) < 0.0)
                        std::reverse(path.begin(), path.end());
                    cover_paths.push_back(std::move(path));
                }
            }
            // A clipped source triangle can be covered jointly by the retained
            // coplanar proxy and by the convex volume that swallowed its buried
            // portion. Intersect every certified volume with this source plane
            // and include those exact sections in the same 2D union audit.
            for (const Enclosure& enclosure : enclosures)
            {
                std::vector<Vec2> section;
                section.reserve(source_path.size());
                for (const auto& point : source_path)
                    section.emplace_back(point.x, point.y);
                for (const Plane& plane : enclosure.planes)
                {
                    if (section.empty()) break;
                    std::vector<Vec2> clipped;
                    clipped.reserve(section.size() + 1);
                    for (std::size_t edge = 0; edge < section.size(); ++edge)
                    {
                        const Vec2 first = section[edge];
                        const Vec2 second = section[(edge + 1) % section.size()];
                        const auto distance = [&](const Vec2& point)
                        {
                            const Vec3 position = triangle[0] +
                                frame.col(0) * point.x() +
                                frame.col(1) * point.y();
                            return plane.normal.dot(position - plane.point);
                        };
                        const double first_distance = distance(first);
                        const double second_distance = distance(second);
                        const bool first_inside = first_distance <= tolerance;
                        const bool second_inside = second_distance <= tolerance;
                        if (first_inside) clipped.push_back(first);
                        if (first_inside == second_inside) continue;
                        const double denominator = first_distance - second_distance;
                        if (std::abs(denominator) <= 1.0e-30) continue;
                        clipped.push_back(first + (second - first) *
                            std::clamp(first_distance / denominator, 0.0, 1.0));
                    }
                    section = simplifyPolygon(std::move(clipped), tolerance);
                }
                if (section.size() < 3 ||
                    std::abs(signedArea(section)) <= tolerance * tolerance)
                    continue;
                Clipper2Lib::PathD path;
                for (const Vec2& point : section)
                    path.emplace_back(point.x(), point.y());
                if (Clipper2Lib::Area(path) < 0.0)
                    std::reverse(path.begin(), path.end());
                cover_paths.push_back(std::move(path));
            }
            if (cover_paths.empty()) return false;
            const Clipper2Lib::PathsD covered_paths = Clipper2Lib::Union(
                cover_paths, Clipper2Lib::FillRule::NonZero, clipper_precision);
            const Clipper2Lib::PathsD remaining = Clipper2Lib::Difference(
                Clipper2Lib::PathsD{source_path}, covered_paths,
                Clipper2Lib::FillRule::NonZero, clipper_precision);
            double remaining_area = 0.0;
            for (const auto& path : remaining)
                remaining_area += std::abs(Clipper2Lib::Area(path));
            return remaining_area <= std::max(
                tolerance * tolerance * 64.0, 0.5 * doubled_area * 1.0e-10);
        };

        covered = coveredByPlanarPrimitives(owners[face_id]);
        if (!covered && owners[face_id].size() != primitives.size())
            covered = coveredByPlanarPrimitives(all_primitives);
        if (covered) ++audit.planar_source_faces;
        else
        {
            ++audit.failed_source_faces;
            audit.failed_face_ids.push_back(static_cast<std::uint32_t>(face_id));
            audit.failed_with_enclosure_owner += std::any_of(
                owners[face_id].begin(), owners[face_id].end(),
                [&](const auto owner)
                { return primitives[owner].enclosure_group != 0; });
            audit.failed_with_planar_owner += std::any_of(
                owners[face_id].begin(), owners[face_id].end(),
                [&](const auto owner)
                {
                    Vec3 ignored;
                    return planarNormal(primitives[owner].primitive,
                                        tolerance, ignored);
                });
            audit.failed_with_certificate_planar_owner += std::any_of(
                certificate_owners[face_id].begin(),
                certificate_owners[face_id].end(),
                [&](const auto owner)
                {
                    Vec3 ignored;
                    return planarNormal(
                        certificate_primitives[owner].primitive,
                        tolerance, ignored);
                });
            audit.failed_with_certificate_band_owner += std::any_of(
                certificate_owners[face_id].begin(),
                certificate_owners[face_id].end(),
                [&](const auto owner)
                {
                    const Kind kind =
                        certificate_primitives[owner].primitive.kind;
                    return kind == Kind::CylindricalBand ||
                           kind == Kind::ConicalBand;
                });
        }
    }
    return audit;
}

std::vector<OutputPrimitive> buildNonOverlappingPlanarCoverageRepair(
    const Mesh& source_mesh,
    const std::vector<OutputPrimitive>& retained,
    const std::vector<std::uint32_t>& failed_face_ids,
    const double tolerance,
    std::size_t& merged_plane_groups)
{
    merged_plane_groups = 0;
    constexpr int clipper_precision = 8;
    struct PlaneKey
    {
        std::array<std::int64_t, 4> values{};
        bool operator<(const PlaneKey& other) const
        {
            return values < other.values;
        }
    };
    struct PlaneGroup
    {
        Vec3 normal = Vec3::Zero();
        Vec3 origin = Vec3::Zero();
        Mat3 frame = Mat3::Identity();
        std::vector<std::uint32_t> faces;
        Clipper2Lib::PathsD source_paths;
    };
    const double angular_quantum = 1.0e-8;
    const double distance_quantum = std::max(tolerance, 1.0e-12);
    const auto canonicalNormal = [](Vec3 normal)
    {
        int dominant = 0;
        for (int axis = 1; axis < 3; ++axis)
            if (std::abs(normal[axis]) > std::abs(normal[dominant]))
                dominant = axis;
        if (normal[dominant] < 0.0) normal *= -1.0;
        return normal;
    };
    const auto planeKey = [&](const Vec3& normal, const Vec3& point)
    {
        PlaneKey key;
        for (int axis = 0; axis < 3; ++axis)
            key.values[axis] = static_cast<std::int64_t>(
                std::llround(normal[axis] / angular_quantum));
        key.values[3] = static_cast<std::int64_t>(
            std::llround(normal.dot(point) / distance_quantum));
        return key;
    };

    std::map<PlaneKey, PlaneGroup> groups;
    for (const auto face_id : failed_face_ids)
    {
        if (face_id >= source_mesh.faces.size()) continue;
        const Face& face = source_mesh.faces[face_id];
        const std::array<Vec3, 3> triangle{
            source_mesh.vertices[face[0]],
            source_mesh.vertices[face[1]],
            source_mesh.vertices[face[2]]};
        Vec3 normal = (triangle[1] - triangle[0]).cross(
            triangle[2] - triangle[0]);
        if (normal.norm() <= tolerance * tolerance) continue;
        normal = canonicalNormal(normal.normalized());
        const PlaneKey key = planeKey(normal, triangle[0]);
        auto [iterator, inserted] = groups.try_emplace(key);
        PlaneGroup& group = iterator->second;
        if (inserted)
        {
            group.normal = normal;
            group.origin = triangle[0];
            const Mat3 basis = orthonormalFrame(normal);
            group.frame.col(0) = basis.col(1);
            group.frame.col(1) = basis.col(2);
            group.frame.col(2) = normal;
        }
        Clipper2Lib::PathD path;
        for (const Vec3& vertex : triangle)
        {
            const Vec3 local =
                group.frame.transposeMultiply(vertex - group.origin);
            path.emplace_back(local.x(), local.y());
        }
        if (Clipper2Lib::Area(path) < 0.0)
            std::reverse(path.begin(), path.end());
        if (std::abs(Clipper2Lib::Area(path)) <= tolerance * tolerance)
            continue;
        group.faces.push_back(face_id);
        group.source_paths.push_back(std::move(path));
    }

    struct PlanarRetained
    {
        Vec3 normal = Vec3::Zero();
        Vec3 point = Vec3::Zero();
        PrimitiveMesh mesh;
    };
    std::vector<PlanarRetained> planar_retained;
    planar_retained.reserve(retained.size());
    for (const OutputPrimitive& item : retained)
    {
        Vec3 normal;
        if (!planarNormal(item.primitive, tolerance, normal)) continue;
        normal = canonicalNormal(normal);
        PrimitiveMesh mesh = triangulatePrimitive(item.primitive);
        if (mesh.faces.empty()) continue;
        planar_retained.push_back(
            {normal, planarPoint(item.primitive), std::move(mesh)});
    }

    std::vector<OutputPrimitive> result;
    const auto appendExactFallback = [&](const PlaneGroup& group)
    {
        for (const auto face_id : group.faces)
        {
            const Face& face = source_mesh.faces[face_id];
            Primitive triangle;
            triangle.kind = Kind::Triangle;
            for (int corner = 0; corner < 3; ++corner)
                triangle.triangle[corner] =
                    source_mesh.vertices[face[corner]];
            result.push_back({std::move(triangle), {face_id}});
        }
    };

    for (const auto& [key, group] : groups)
    {
        (void)key;
        if (group.source_paths.empty()) continue;
        const Clipper2Lib::PathsD source_union = Clipper2Lib::Union(
            group.source_paths, Clipper2Lib::FillRule::NonZero,
            clipper_precision);
        Clipper2Lib::PathsD cover_paths;
        for (const PlanarRetained& candidate : planar_retained)
        {
            if (std::abs(candidate.normal.dot(group.normal)) <
                    1.0 - 1.0e-8 ||
                std::abs((candidate.point - group.origin).dot(group.normal)) >
                    tolerance)
                continue;
            for (const Face& face : candidate.mesh.faces)
            {
                Clipper2Lib::PathD path;
                for (const auto vertex : face)
                {
                    const Vec3 local = group.frame.transposeMultiply(
                        candidate.mesh.vertices[vertex] - group.origin);
                    path.emplace_back(local.x(), local.y());
                }
                const double area = Clipper2Lib::Area(path);
                if (std::abs(area) <= tolerance * tolerance) continue;
                if (area < 0.0) std::reverse(path.begin(), path.end());
                cover_paths.push_back(std::move(path));
            }
        }
        const Clipper2Lib::PathsD cover_union = cover_paths.empty()
            ? Clipper2Lib::PathsD{}
            : Clipper2Lib::Union(
                cover_paths, Clipper2Lib::FillRule::NonZero,
                clipper_precision);
        Clipper2Lib::PathsD remaining = source_union;
        if (!cover_union.empty())
            remaining = Clipper2Lib::Difference(
                source_union, cover_union,
                Clipper2Lib::FillRule::NonZero, clipper_precision);
        if (remaining.empty()) continue;

        const auto appendPerFaceResidual = [&]()
        {
            for (std::size_t face_index = 0;
                 face_index < group.faces.size(); ++face_index)
            {
                Clipper2Lib::PathsD residual{
                    group.source_paths[face_index]};
                if (!cover_union.empty())
                    residual = Clipper2Lib::Difference(
                        residual, cover_union,
                        Clipper2Lib::FillRule::NonZero,
                        clipper_precision);
                if (residual.empty()) continue;
                Clipper2Lib::PathsD paths = residual;
                if (std::any_of(residual.begin(), residual.end(),
                    [](const auto& path)
                    { return !Clipper2Lib::IsPositive(path); }))
                {
                    Clipper2Lib::PathsD triangles;
                    if (Clipper2Lib::Triangulate(
                            residual, clipper_precision, triangles, false) !=
                        Clipper2Lib::TriangulateResult::success)
                    {
                        PlaneGroup one_face = group;
                        one_face.faces = {group.faces[face_index]};
                        appendExactFallback(one_face);
                        continue;
                    }
                    paths = std::move(triangles);
                }
                bool emitted = false;
                for (const auto& path : paths)
                {
                    if (!Clipper2Lib::IsPositive(path) ||
                        std::abs(Clipper2Lib::Area(path)) <=
                            tolerance * tolerance)
                        continue;
                    std::vector<Vec2> boundary;
                    for (const auto& point : path)
                        boundary.emplace_back(point.x, point.y);
                    boundary = simplifyPolygon(
                        std::move(boundary), tolerance);
                    if (boundary.size() < 3 ||
                        !simplePolygon(boundary, tolerance) ||
                        triangulatePolygon(boundary).size() + 2 !=
                            boundary.size())
                        continue;
                    OutputPrimitive output;
                    output.primitive = polygonPrimitive(
                        boundary, group.origin, group.frame);
                    output.source_faces = {group.faces[face_index]};
                    result.push_back(std::move(output));
                    emitted = true;
                }
                if (!emitted)
                {
                    PlaneGroup one_face = group;
                    one_face.faces = {group.faces[face_index]};
                    appendExactFallback(one_face);
                }
            }
        };

        Clipper2Lib::PathsD output_paths = remaining;
        if (std::any_of(remaining.begin(), remaining.end(),
            [](const auto& path) { return !Clipper2Lib::IsPositive(path); }))
        {
            Clipper2Lib::PathsD triangles;
            if (Clipper2Lib::Triangulate(
                    remaining, clipper_precision, triangles, false) !=
                Clipper2Lib::TriangulateResult::success)
            {
                appendPerFaceResidual();
                continue;
            }
            output_paths = std::move(triangles);
        }

        struct RepairPiece
        {
            OutputPrimitive output;
            Bounds2 bounds;
        };
        std::vector<RepairPiece> pieces;
        bool conversion_failed = false;
        for (const auto& path : output_paths)
        {
            if (!Clipper2Lib::IsPositive(path) ||
                std::abs(Clipper2Lib::Area(path)) <=
                    tolerance * tolerance)
                continue;
            std::vector<Vec2> boundary;
            Bounds2 bounds;
            boundary.reserve(path.size());
            for (const auto& point : path)
            {
                const Vec2 value(point.x, point.y);
                boundary.push_back(value);
                bounds.lower = bounds.lower.cwiseMin(value);
                bounds.upper = bounds.upper.cwiseMax(value);
            }
            boundary = simplifyPolygon(std::move(boundary), tolerance);
            if (boundary.size() < 3 || !simplePolygon(boundary, tolerance) ||
                triangulatePolygon(boundary).size() + 2 != boundary.size())
            {
                conversion_failed = true;
                break;
            }
            OutputPrimitive output;
            output.primitive = polygonPrimitive(
                boundary, group.origin, group.frame);
            pieces.push_back({std::move(output), bounds});
        }
        if (conversion_failed || pieces.empty())
        {
            appendPerFaceResidual();
            continue;
        }

        for (const auto face_id : group.faces)
        {
            Bounds2 face_bounds;
            for (const auto vertex : source_mesh.faces[face_id])
            {
                const Vec3 local = group.frame.transposeMultiply(
                    source_mesh.vertices[vertex] - group.origin);
                face_bounds.lower = face_bounds.lower.cwiseMin(
                    {local.x(), local.y()});
                face_bounds.upper = face_bounds.upper.cwiseMax(
                    {local.x(), local.y()});
            }
            for (RepairPiece& piece : pieces)
                if (!(piece.bounds.upper.x() + tolerance <
                          face_bounds.lower.x() ||
                      face_bounds.upper.x() + tolerance <
                          piece.bounds.lower.x() ||
                      piece.bounds.upper.y() + tolerance <
                          face_bounds.lower.y() ||
                      face_bounds.upper.y() + tolerance <
                          piece.bounds.lower.y()))
                    piece.output.source_faces.push_back(face_id);
        }
        for (RepairPiece& piece : pieces)
        {
            if (piece.output.source_faces.empty())
                piece.output.source_faces = group.faces;
            result.push_back(std::move(piece.output));
        }
        if (pieces.size() < group.faces.size()) ++merged_plane_groups;
    }
    return result;
}

std::vector<OutputPrimitive> removeRedundantEnclosureGroupsByUnionCoverage(
    const Mesh& source_mesh,
    std::vector<OutputPrimitive> primitives,
    const double tolerance,
    std::size_t& removed_count)
{
    removed_count = 0;
    for (;;)
    {
        std::unordered_map<std::uint64_t, std::vector<std::size_t>> groups;
        for (std::size_t index = 0; index < primitives.size(); ++index)
            if (primitives[index].enclosure_group != 0)
                groups[primitives[index].enclosure_group].push_back(index);
        if (groups.size() < 2) break;

        std::vector<std::uint64_t> order;
        order.reserve(groups.size());
        for (const auto& [group, members] : groups)
        {
            (void)members;
            order.push_back(group);
        }
        std::sort(order.begin(), order.end(), [&](const auto first,
                                                  const auto second)
        {
            const auto groupVolume = [&](const std::uint64_t group)
            {
                Bounds bounds;
                for (const auto member : groups[group])
                {
                    const Bounds item = primitiveBounds(
                        primitives[member].primitive);
                    bounds.lower = bounds.lower.cwiseMin(item.lower);
                    bounds.upper = bounds.upper.cwiseMax(item.upper);
                }
                return (bounds.upper - bounds.lower)
                    .cwiseMax(Vec3::Zero()).prod();
            };
            const double first_volume = groupVolume(first);
            const double second_volume = groupVolume(second);
            if (first_volume != second_volume)
                return first_volume < second_volume;
            return first < second;
        });

        bool changed = false;
        for (const std::uint64_t group : order)
        {
            const auto members_iterator = groups.find(group);
            if (members_iterator == groups.end()) continue;
            std::vector<std::uint8_t> member_mask(primitives.size(), 0);
            std::vector<std::uint32_t> responsibility;
            for (const auto member : members_iterator->second)
            {
                member_mask[member] = 1;
                responsibility.insert(
                    responsibility.end(),
                    primitives[member].source_faces.begin(),
                    primitives[member].source_faces.end());
            }
            std::sort(responsibility.begin(), responsibility.end());
            responsibility.erase(
                std::unique(responsibility.begin(), responsibility.end()),
                responsibility.end());
            if (responsibility.empty()) continue;

            std::vector<OutputPrimitive> tentative;
            tentative.reserve(primitives.size() - members_iterator->second.size());
            for (std::size_t index = 0; index < primitives.size(); ++index)
                if (!member_mask[index]) tentative.push_back(primitives[index]);
            if (tentative.empty()) continue;
            // Responsibility is bookkeeping only here. The audit first checks
            // every retained convex enclosure, then the exact union of all
            // retained coplanar proxy faces, independently of this owner.
            tentative.front().source_faces.insert(
                tentative.front().source_faces.end(), responsibility.begin(),
                responsibility.end());
            std::sort(tentative.front().source_faces.begin(),
                      tentative.front().source_faces.end());
            tentative.front().source_faces.erase(
                std::unique(tentative.front().source_faces.begin(),
                            tentative.front().source_faces.end()),
                tentative.front().source_faces.end());

            const FinalCoverageAudit audit = auditFinalConservativeCoverage(
                source_mesh, tentative, tolerance);
            if (audit.unassigned_source_faces != 0 ||
                audit.failed_source_faces != 0)
                continue;
            removed_count += members_iterator->second.size();
            primitives = std::move(tentative);
            changed = true;
            break;
        }
        if (!changed) break;
    }
    return primitives;
}

void writeCoverageFailureDiagnostics(
    const std::filesystem::path& directory,
    const Mesh& source_mesh,
    const std::vector<OutputPrimitive>& primitives,
    const FinalCoverageAudit& audit)
{
    std::ofstream obj(directory / "coverage_failed_faces.obj");
    obj << std::setprecision(17);
    for (const Vec3& vertex : source_mesh.vertices)
        obj << "v " << vertex.x() << ' ' << vertex.y() << ' '
            << vertex.z() << '\n';
    for (const auto face_id : audit.failed_face_ids)
    {
        const Face& face = source_mesh.faces[face_id];
        obj << "f " << face[0] + 1 << ' ' << face[1] + 1 << ' '
            << face[2] + 1 << '\n';
    }

    std::ofstream report(directory / "coverage_failure.txt");
    report << "unassigned_source_faces=" << audit.unassigned_source_faces << '\n'
           << "failed_source_faces=" << audit.failed_source_faces << '\n'
           << "failed_with_enclosure_owner="
           << audit.failed_with_enclosure_owner << '\n'
           << "failed_with_planar_owner=" << audit.failed_with_planar_owner << '\n'
           << "failed_face_ids=";
    for (std::size_t index = 0; index < audit.failed_face_ids.size(); ++index)
    {
        if (index != 0) report << ',';
        report << audit.failed_face_ids[index];
    }
    report << '\n';

    std::vector<std::vector<std::size_t>> owners(source_mesh.faces.size());
    for (std::size_t primitive = 0; primitive < primitives.size(); ++primitive)
        for (const auto face : primitives[primitive].source_faces)
            if (face < owners.size()) owners[face].push_back(primitive);
    for (const auto face_id : audit.failed_face_ids)
    {
        report << "face " << face_id << " owners=";
        for (std::size_t owner_index = 0; owner_index < owners[face_id].size();
             ++owner_index)
        {
            if (owner_index != 0) report << ',';
            const auto owner = owners[face_id][owner_index];
            report << owner << ':' << kindName(primitives[owner].primitive.kind)
                   << ":group=" << primitives[owner].enclosure_group;
        }
        report << '\n';
    }
}

bool pointInTriangleInclusive(const Vec2& point, const std::array<Vec2, 3>& triangle)
{
    const double first = cross2(triangle[1] - triangle[0], point - triangle[0]);
    const double second = cross2(triangle[2] - triangle[1], point - triangle[1]);
    const double third = cross2(triangle[0] - triangle[2], point - triangle[2]);
    return (first >= -1.0e-10 && second >= -1.0e-10 && third >= -1.0e-10) ||
           (first <= 1.0e-10 && second <= 1.0e-10 && third <= 1.0e-10);
}

bool projectedTriangleIntersectsCell(const std::array<Vec2, 3>& triangle,
                                     const int x, const int y)
{
    const double lower_x = static_cast<double>(x);
    const double lower_y = static_cast<double>(y);
    const double upper_x = lower_x + 1.0;
    const double upper_y = lower_y + 1.0;
    for (const Vec2& point : triangle)
        if (point.x() >= lower_x && point.x() <= upper_x &&
            point.y() >= lower_y && point.y() <= upper_y) return true;
    const std::array<Vec2, 4> corners{{{lower_x, lower_y}, {upper_x, lower_y},
                                      {upper_x, upper_y}, {lower_x, upper_y}}};
    for (const Vec2& corner : corners)
        if (pointInTriangleInclusive(corner, triangle)) return true;
    for (int edge = 0; edge < 3; ++edge)
        for (int side = 0; side < 4; ++side)
            if (properSegmentIntersection(
                    triangle[edge], triangle[(edge + 1) % 3],
                    corners[side], corners[(side + 1) % 4], 1.0e-12)) return true;
    return false;
}

std::vector<bool> rasterizeProjection(const Mesh& mesh,
                                      const int first_axis,
                                      const int second_axis,
                                      const Vec3& lower,
                                      const Vec3& step,
                                      const int first_size,
                                      const int second_size)
{
    std::vector<bool> mask(static_cast<std::size_t>(first_size) * second_size, false);
    for (const Face& face : mesh.faces)
    {
        std::array<Vec2, 3> triangle;
        for (int corner = 0; corner < 3; ++corner)
        {
            const Vec3& vertex = mesh.vertices[face[corner]];
            triangle[corner] = {
                (vertex[first_axis] - lower[first_axis]) / step[first_axis],
                (vertex[second_axis] - lower[second_axis]) / step[second_axis]};
        }
        double minimum_x = triangle[0].x();
        double maximum_x = triangle[0].x();
        double minimum_y = triangle[0].y();
        double maximum_y = triangle[0].y();
        for (int corner = 1; corner < 3; ++corner)
        {
            minimum_x = std::min(minimum_x, triangle[corner].x());
            maximum_x = std::max(maximum_x, triangle[corner].x());
            minimum_y = std::min(minimum_y, triangle[corner].y());
            maximum_y = std::max(maximum_y, triangle[corner].y());
        }
        const int begin_x = std::clamp(static_cast<int>(std::floor(minimum_x)), 0,
                                       first_size - 1);
        const int end_x = std::clamp(static_cast<int>(std::floor(maximum_x)), 0,
                                     first_size - 1);
        const int begin_y = std::clamp(static_cast<int>(std::floor(minimum_y)), 0,
                                       second_size - 1);
        const int end_y = std::clamp(static_cast<int>(std::floor(maximum_y)), 0,
                                     second_size - 1);
        for (int y = begin_y; y <= end_y; ++y)
            for (int x = begin_x; x <= end_x; ++x)
                if (projectedTriangleIntersectsCell(triangle, x, y))
                    mask[static_cast<std::size_t>(y) * first_size + x] = true;
    }
    return mask;
}

void fillProjectionHoles(std::vector<bool>& mask, const int width, const int height)
{
    std::vector<bool> exterior(mask.size(), false);
    std::queue<std::pair<int, int>> queue;
    const auto push = [&](const int x, const int y)
    {
        const std::size_t index = static_cast<std::size_t>(y) * width + x;
        if (mask[index] || exterior[index]) return;
        exterior[index] = true;
        queue.emplace(x, y);
    };
    for (int x = 0; x < width; ++x) { push(x, 0); push(x, height - 1); }
    for (int y = 0; y < height; ++y) { push(0, y); push(width - 1, y); }
    constexpr std::array<std::pair<int, int>, 4> offsets{{{1, 0}, {-1, 0}, {0, 1}, {0, -1}}};
    while (!queue.empty())
    {
        const auto [x, y] = queue.front();
        queue.pop();
        for (const auto [dx, dy] : offsets)
        {
            const int next_x = x + dx;
            const int next_y = y + dy;
            if (next_x >= 0 && next_x < width && next_y >= 0 && next_y < height)
                push(next_x, next_y);
        }
    }
    for (std::size_t index = 0; index < mask.size(); ++index)
        if (!mask[index] && !exterior[index]) mask[index] = true;
}

double boxVisualHullAddedVolumeRatio(
    const Mesh& mesh,
    const std::vector<std::uint32_t>& faces,
    const BoxFit& box,
    const double model_volume,
    const int resolution)
{
    if (resolution < 4 || box.volume <= 0.0) return 0.0;
    const Vec3 extent = box.half_size * 2.0;
    if (extent.x() <= 0.0 || extent.y() <= 0.0 || extent.z() <= 0.0)
        return 0.0;

    const auto rasterize = [&](const int first_axis, const int second_axis)
    {
        std::vector<bool> mask(
            static_cast<std::size_t>(resolution) * resolution, false);
        for (const auto face_id : faces)
        {
            std::array<Vec2, 3> triangle;
            for (int corner = 0; corner < 3; ++corner)
            {
                const Vec3 local = box.axes.transposeMultiply(
                    mesh.vertices[mesh.faces[face_id][corner]] - box.center);
                triangle[corner] = {
                    (local[first_axis] + box.half_size[first_axis]) /
                        extent[first_axis] * resolution,
                    (local[second_axis] + box.half_size[second_axis]) /
                        extent[second_axis] * resolution,
                };
            }
            double minimum_x = triangle[0].x();
            double maximum_x = triangle[0].x();
            double minimum_y = triangle[0].y();
            double maximum_y = triangle[0].y();
            for (int corner = 1; corner < 3; ++corner)
            {
                minimum_x = std::min(minimum_x, triangle[corner].x());
                maximum_x = std::max(maximum_x, triangle[corner].x());
                minimum_y = std::min(minimum_y, triangle[corner].y());
                maximum_y = std::max(maximum_y, triangle[corner].y());
            }
            const int begin_x = std::clamp(
                static_cast<int>(std::floor(minimum_x)), 0, resolution - 1);
            const int end_x = std::clamp(
                static_cast<int>(std::floor(maximum_x)), 0, resolution - 1);
            const int begin_y = std::clamp(
                static_cast<int>(std::floor(minimum_y)), 0, resolution - 1);
            const int end_y = std::clamp(
                static_cast<int>(std::floor(maximum_y)), 0, resolution - 1);
            for (int y = begin_y; y <= end_y; ++y)
                for (int x = begin_x; x <= end_x; ++x)
                    if (projectedTriangleIntersectsCell(triangle, x, y))
                        mask[static_cast<std::size_t>(y) * resolution + x] = true;
        }
        fillProjectionHoles(mask, resolution, resolution);
        return mask;
    };

    const auto xy = rasterize(0, 1);
    const auto xz = rasterize(0, 2);
    const auto yz = rasterize(1, 2);
    std::uint64_t visual_hull_cells = 0;
    for (int z = 0; z < resolution; ++z)
        for (int y = 0; y < resolution; ++y)
            for (int x = 0; x < resolution; ++x)
                visual_hull_cells +=
                    xy[static_cast<std::size_t>(y) * resolution + x] &&
                    xz[static_cast<std::size_t>(z) * resolution + x] &&
                    yz[static_cast<std::size_t>(z) * resolution + y];
    const double occupied_fraction = static_cast<double>(visual_hull_cells) /
        static_cast<double>(resolution * resolution * resolution);
    return box.volume * std::max(1.0 - occupied_fraction, 0.0) /
        std::max(model_volume, 1.0e-30);
}

double boxFloodFillAddedVolumeRatio(
    const Mesh& mesh,
    const std::vector<std::uint32_t>& faces,
    const BoxFit& box,
    const double model_volume,
    const int resolution)
{
    if (resolution < 8 || box.volume <= 0.0) return 0.0;
    const Vec3 extent = box.half_size * 2.0;
    if (extent.x() <= 0.0 || extent.y() <= 0.0 || extent.z() <= 0.0)
        return 0.0;
    const std::size_t cell_count = static_cast<std::size_t>(resolution) *
        resolution * resolution;
    std::vector<std::uint8_t> surface(cell_count, 0);
    const auto indexOf = [&](const int x, const int y, const int z)
    {
        return (static_cast<std::size_t>(z) * resolution + y) * resolution + x;
    };
    const auto mark = [&](const int x, const int y, const int z)
    {
        if (x >= 0 && x < resolution && y >= 0 && y < resolution &&
            z >= 0 && z < resolution)
            surface[indexOf(x, y, z)] = 1;
    };
    const auto gridPoint = [&](const Vec3& world)
    {
        const Vec3 local = box.axes.transposeMultiply(world - box.center);
        return Vec3(
            (local.x() + box.half_size.x()) / extent.x() * resolution,
            (local.y() + box.half_size.y()) / extent.y() * resolution,
            (local.z() + box.half_size.z()) / extent.z() * resolution);
    };

    for (const auto face_id : faces)
    {
        std::array<Vec3, 3> triangle;
        for (int corner = 0; corner < 3; ++corner)
            triangle[corner] = gridPoint(mesh.vertices[mesh.faces[face_id][corner]]);
        const Vec3 normal = (triangle[1] - triangle[0]).cross(
            triangle[2] - triangle[0]);
        int normal_axis = 0;
        for (int axis = 1; axis < 3; ++axis)
            if (std::abs(normal[axis]) > std::abs(normal[normal_axis]))
                normal_axis = axis;
        if (std::abs(normal[normal_axis]) <= 1.0e-12) continue;
        const int first_axis = (normal_axis + 1) % 3;
        const int second_axis = (normal_axis + 2) % 3;
        std::array<Vec2, 3> projected;
        for (int corner = 0; corner < 3; ++corner)
            projected[corner] = {
                triangle[corner][first_axis], triangle[corner][second_axis]};
        double minimum_first = projected[0].x();
        double maximum_first = projected[0].x();
        double minimum_second = projected[0].y();
        double maximum_second = projected[0].y();
        for (int corner = 1; corner < 3; ++corner)
        {
            minimum_first = std::min(minimum_first, projected[corner].x());
            maximum_first = std::max(maximum_first, projected[corner].x());
            minimum_second = std::min(minimum_second, projected[corner].y());
            maximum_second = std::max(maximum_second, projected[corner].y());
        }
        const int begin_first = std::clamp(
            static_cast<int>(std::floor(minimum_first)), 0, resolution - 1);
        const int end_first = std::clamp(
            static_cast<int>(std::floor(maximum_first)), 0, resolution - 1);
        const int begin_second = std::clamp(
            static_cast<int>(std::floor(minimum_second)), 0, resolution - 1);
        const int end_second = std::clamp(
            static_cast<int>(std::floor(maximum_second)), 0, resolution - 1);
        for (int second = begin_second; second <= end_second; ++second)
            for (int first = begin_first; first <= end_first; ++first)
            {
                if (!projectedTriangleIntersectsCell(
                        projected, first, second)) continue;
                Vec3 sample = Vec3::Zero();
                sample[first_axis] = static_cast<double>(first) + 0.5;
                sample[second_axis] = static_cast<double>(second) + 0.5;
                sample[normal_axis] = triangle[0][normal_axis] -
                    (normal[first_axis] *
                         (sample[first_axis] - triangle[0][first_axis]) +
                     normal[second_axis] *
                         (sample[second_axis] - triangle[0][second_axis])) /
                        normal[normal_axis];
                const int normal_cell = std::clamp(
                    static_cast<int>(std::floor(sample[normal_axis])),
                    0, resolution - 1);
                std::array<int, 3> cell{};
                cell[first_axis] = first;
                cell[second_axis] = second;
                for (int offset = -1; offset <= 1; ++offset)
                {
                    cell[normal_axis] = normal_cell + offset;
                    mark(cell[0], cell[1], cell[2]);
                }
            }
    }

    std::vector<std::uint8_t> exterior(cell_count, 0);
    std::queue<std::array<int, 3>> queue;
    const auto pushExterior = [&](const int x, const int y, const int z)
    {
        const std::size_t index = indexOf(x, y, z);
        if (surface[index] || exterior[index]) return;
        exterior[index] = 1;
        queue.push({x, y, z});
    };
    for (int z = 0; z < resolution; ++z)
        for (int y = 0; y < resolution; ++y)
        {
            pushExterior(0, y, z);
            pushExterior(resolution - 1, y, z);
        }
    for (int z = 0; z < resolution; ++z)
        for (int x = 0; x < resolution; ++x)
        {
            pushExterior(x, 0, z);
            pushExterior(x, resolution - 1, z);
        }
    for (int y = 0; y < resolution; ++y)
        for (int x = 0; x < resolution; ++x)
        {
            pushExterior(x, y, 0);
            pushExterior(x, y, resolution - 1);
        }
    constexpr std::array<std::array<int, 3>, 6> offsets{{
        {{1, 0, 0}}, {{-1, 0, 0}}, {{0, 1, 0}},
        {{0, -1, 0}}, {{0, 0, 1}}, {{0, 0, -1}},
    }};
    while (!queue.empty())
    {
        const auto cell = queue.front();
        queue.pop();
        for (const auto& offset : offsets)
        {
            const int x = cell[0] + offset[0];
            const int y = cell[1] + offset[1];
            const int z = cell[2] + offset[2];
            if (x >= 0 && x < resolution && y >= 0 && y < resolution &&
                z >= 0 && z < resolution)
                pushExterior(x, y, z);
        }
    }
    const std::size_t exterior_cells = static_cast<std::size_t>(std::count(
        exterior.begin(), exterior.end(), std::uint8_t{1}));
    const double exterior_fraction = static_cast<double>(exterior_cells) /
        static_cast<double>(cell_count);
    return box.volume * exterior_fraction / std::max(model_volume, 1.0e-30);
}

struct FinalOpenErrorAudit
{
    std::size_t distance_sample_count = 0;
    double mean_distance = 0.0;
    double maximum_distance = 0.0;
    double distance_area_integral = 0.0;
    Vec3 maximum_proxy_point = Vec3::Zero();
    Vec3 maximum_source_point = Vec3::Zero();
};
PrimitiveMesh triangulateOutputPrimitives(
    const std::vector<OutputPrimitive>& output)
{
    PrimitiveMesh result;
    for (const OutputPrimitive& item : output)
    {
        PrimitiveMesh part = triangulatePrimitive(item.primitive);
        const std::uint32_t base = static_cast<std::uint32_t>(
            result.vertices.size());
        result.vertices.insert(result.vertices.end(), part.vertices.begin(),
                               part.vertices.end());
        for (Face face : part.faces)
        {
            for (auto& vertex : face) vertex += base;
            result.faces.push_back(face);
        }
    }
    return result;
}

std::size_t triangulatedFaceCount(
    const std::vector<OutputPrimitive>& output)
{
    std::size_t count = 0;
    for (const OutputPrimitive& item : output)
        count += triangulatePrimitive(item.primitive).faces.size();
    return count;
}

Vec3 closestPointOnTriangle(const Vec3& point,
                            const Vec3& first,
                            const Vec3& second,
                            const Vec3& third)
{
    const Vec3 first_edge = second - first;
    const Vec3 second_edge = third - first;
    const Vec3 first_relative = point - first;
    const double first_projection = first_edge.dot(first_relative);
    const double second_projection = second_edge.dot(first_relative);
    if (first_projection <= 0.0 && second_projection <= 0.0) return first;

    const Vec3 second_relative = point - second;
    const double third_projection = first_edge.dot(second_relative);
    const double fourth_projection = second_edge.dot(second_relative);
    if (third_projection >= 0.0 && fourth_projection <= third_projection)
        return second;

    const double first_edge_region = first_projection * fourth_projection -
        third_projection * second_projection;
    if (first_edge_region <= 0.0 && first_projection >= 0.0 &&
        third_projection <= 0.0)
        return first + first_edge *
            (first_projection / (first_projection - third_projection));

    const Vec3 third_relative = point - third;
    const double fifth_projection = first_edge.dot(third_relative);
    const double sixth_projection = second_edge.dot(third_relative);
    if (sixth_projection >= 0.0 && fifth_projection <= sixth_projection)
        return third;

    const double second_edge_region = fifth_projection * second_projection -
        first_projection * sixth_projection;
    if (second_edge_region <= 0.0 && second_projection >= 0.0 &&
        sixth_projection <= 0.0)
        return first + second_edge *
            (second_projection / (second_projection - sixth_projection));

    const double opposite_edge_region = third_projection * sixth_projection -
        fifth_projection * fourth_projection;
    if (opposite_edge_region <= 0.0 &&
        fourth_projection - third_projection >= 0.0 &&
        fifth_projection - sixth_projection >= 0.0)
        return second + (third - second) *
            ((fourth_projection - third_projection) /
             ((fourth_projection - third_projection) +
              (fifth_projection - sixth_projection)));

    const double denominator = 1.0 /
        (opposite_edge_region + second_edge_region + first_edge_region);
    const double second_weight = second_edge_region * denominator;
    const double third_weight = first_edge_region * denominator;
    return first + first_edge * second_weight + second_edge * third_weight;
}

class SourceTriangleBvh
{
public:
    explicit SourceTriangleBvh(const Mesh& mesh) : mesh_(mesh)
    {
        face_ids_.resize(mesh.faces.size());
        std::iota(face_ids_.begin(), face_ids_.end(), 0);
        if (!face_ids_.empty()) build(0, face_ids_.size());
    }

    struct ClosestPointResult
    {
        Vec3 point = Vec3::Zero();
        std::size_t face_id = std::numeric_limits<std::size_t>::max();
    };

    [[nodiscard]] ClosestPointResult closestPoint(const Vec3& point) const
    {
        ClosestPointResult result;
        result.point = point;
        if (nodes_.empty()) return result;
        double best_distance_squared = std::numeric_limits<double>::infinity();
        query(0, point, best_distance_squared, result);
        return result;
    }

    [[nodiscard]] bool rayIntersectsAny(const Vec3& origin,
                                        const Vec3& direction,
                                        const std::size_t excluded_face,
                                        const double minimum_distance) const
    {
        return !nodes_.empty() && rayQuery(
            0, origin, direction, excluded_face, minimum_distance);
    }

private:
    struct Node
    {
        Bounds bounds;
        std::size_t begin = 0;
        std::size_t end = 0;
        std::size_t left = 0;
        std::size_t right = 0;
        bool leaf = false;
    };

    [[nodiscard]] Bounds faceBounds(const std::size_t face_id) const
    {
        Bounds bounds;
        for (const auto vertex : mesh_.faces[face_id])
        {
            bounds.lower = bounds.lower.cwiseMin(mesh_.vertices[vertex]);
            bounds.upper = bounds.upper.cwiseMax(mesh_.vertices[vertex]);
        }
        return bounds;
    }

    [[nodiscard]] Vec3 faceCentroid(const std::size_t face_id) const
    {
        const Face& face = mesh_.faces[face_id];
        return (mesh_.vertices[face[0]] + mesh_.vertices[face[1]] +
                mesh_.vertices[face[2]]) / 3.0;
    }

    std::size_t build(const std::size_t begin, const std::size_t end)
    {
        Node node;
        node.begin = begin;
        node.end = end;
        Bounds centroid_bounds;
        for (std::size_t index = begin; index < end; ++index)
        {
            const Bounds bounds = faceBounds(face_ids_[index]);
            node.bounds.lower = node.bounds.lower.cwiseMin(bounds.lower);
            node.bounds.upper = node.bounds.upper.cwiseMax(bounds.upper);
            const Vec3 centroid = faceCentroid(face_ids_[index]);
            centroid_bounds.lower = centroid_bounds.lower.cwiseMin(centroid);
            centroid_bounds.upper = centroid_bounds.upper.cwiseMax(centroid);
        }
        const std::size_t node_index = nodes_.size();
        nodes_.push_back(node);
        if (end - begin <= 8)
        {
            nodes_[node_index].leaf = true;
            return node_index;
        }
        const Vec3 centroid_extent = centroid_bounds.upper - centroid_bounds.lower;
        int axis = 0;
        if (centroid_extent.y() > centroid_extent.x()) axis = 1;
        if (centroid_extent.z() > centroid_extent[axis]) axis = 2;
        const std::size_t middle = begin + (end - begin) / 2;
        std::nth_element(face_ids_.begin() + begin, face_ids_.begin() + middle,
                         face_ids_.begin() + end,
            [&](const std::size_t first, const std::size_t second)
            {
                return faceCentroid(first)[axis] < faceCentroid(second)[axis];
            });
        nodes_[node_index].left = build(begin, middle);
        nodes_[node_index].right = build(middle, end);
        return node_index;
    }

    [[nodiscard]] static double distanceSquaredToBounds(
        const Vec3& point, const Bounds& bounds)
    {
        double result = 0.0;
        for (int axis = 0; axis < 3; ++axis)
        {
            const double delta = point[axis] < bounds.lower[axis]
                ? bounds.lower[axis] - point[axis]
                : (point[axis] > bounds.upper[axis]
                    ? point[axis] - bounds.upper[axis] : 0.0);
            result += delta * delta;
        }
        return result;
    }

    void query(const std::size_t node_index,
               const Vec3& point,
               double& best_distance_squared,
               ClosestPointResult& best) const
    {
        const Node& node = nodes_[node_index];
        if (distanceSquaredToBounds(point, node.bounds) >= best_distance_squared)
            return;
        if (node.leaf)
        {
            for (std::size_t index = node.begin; index < node.end; ++index)
            {
                const Face& face = mesh_.faces[face_ids_[index]];
                const Vec3 candidate = closestPointOnTriangle(
                    point, mesh_.vertices[face[0]], mesh_.vertices[face[1]],
                    mesh_.vertices[face[2]]);
                const double distance_squared = (candidate - point).squaredNorm();
                if (distance_squared < best_distance_squared)
                {
                    best_distance_squared = distance_squared;
                    best.point = candidate;
                    best.face_id = face_ids_[index];
                }
            }
            return;
        }
        const double left_distance = distanceSquaredToBounds(
            point, nodes_[node.left].bounds);
        const double right_distance = distanceSquaredToBounds(
            point, nodes_[node.right].bounds);
        if (left_distance <= right_distance)
        {
            query(node.left, point, best_distance_squared, best);
            query(node.right, point, best_distance_squared, best);
        }
        else
        {
            query(node.right, point, best_distance_squared, best);
            query(node.left, point, best_distance_squared, best);
        }
    }

    [[nodiscard]] static bool rayIntersectsBounds(
        const Vec3& origin, const Vec3& direction, const Bounds& bounds,
        const double minimum_distance)
    {
        double entry = minimum_distance;
        double exit = std::numeric_limits<double>::infinity();
        for (int axis = 0; axis < 3; ++axis)
        {
            if (std::abs(direction[axis]) <= 1.0e-30)
            {
                if (origin[axis] < bounds.lower[axis] ||
                    origin[axis] > bounds.upper[axis]) return false;
                continue;
            }
            double first = (bounds.lower[axis] - origin[axis]) / direction[axis];
            double second = (bounds.upper[axis] - origin[axis]) / direction[axis];
            if (first > second) std::swap(first, second);
            entry = std::max(entry, first);
            exit = std::min(exit, second);
            if (entry > exit) return false;
        }
        return exit >= minimum_distance;
    }

    [[nodiscard]] bool rayIntersectsTriangle(
        const Vec3& origin, const Vec3& direction, const std::size_t face_id,
        const double minimum_distance) const
    {
        const Face& face = mesh_.faces[face_id];
        const Vec3& first = mesh_.vertices[face[0]];
        const Vec3 first_edge = mesh_.vertices[face[1]] - first;
        const Vec3 second_edge = mesh_.vertices[face[2]] - first;
        const Vec3 cross = direction.cross(second_edge);
        const double determinant = first_edge.dot(cross);
        if (std::abs(determinant) <= 1.0e-15) return false;
        const double inverse = 1.0 / determinant;
        const Vec3 relative = origin - first;
        const double u = relative.dot(cross) * inverse;
        if (u < -1.0e-10 || u > 1.0 + 1.0e-10) return false;
        const Vec3 second_cross = relative.cross(first_edge);
        const double v = direction.dot(second_cross) * inverse;
        if (v < -1.0e-10 || u + v > 1.0 + 1.0e-10) return false;
        const double distance = second_edge.dot(second_cross) * inverse;
        return distance > minimum_distance;
    }

    [[nodiscard]] bool rayQuery(
        const std::size_t node_index, const Vec3& origin,
        const Vec3& direction, const std::size_t excluded_face,
        const double minimum_distance) const
    {
        const Node& node = nodes_[node_index];
        if (!rayIntersectsBounds(
                origin, direction, node.bounds, minimum_distance)) return false;
        if (node.leaf)
        {
            for (std::size_t index = node.begin; index < node.end; ++index)
                if (face_ids_[index] != excluded_face &&
                    rayIntersectsTriangle(origin, direction, face_ids_[index],
                                          minimum_distance)) return true;
            return false;
        }
        return rayQuery(node.left, origin, direction, excluded_face,
                        minimum_distance) ||
               rayQuery(node.right, origin, direction, excluded_face,
                        minimum_distance);
    }

    const Mesh& mesh_;
    std::vector<std::size_t> face_ids_;
    std::vector<Node> nodes_;
};

struct SourceOpenSurfaceReference
{
    explicit SourceOpenSurfaceReference(const Mesh& source_mesh,
                                        const double model_diagonal)
        : mesh(source_mesh), bvh(source_mesh), open_faces(source_mesh.faces.size(), 0)
    {
        const double ray_epsilon = std::max(model_diagonal * 1.0e-9, 1.0e-10);
        for (std::size_t face_id = 0; face_id < mesh.faces.size(); ++face_id)
        {
            const Face& face = mesh.faces[face_id];
            const Vec3& first = mesh.vertices[face[0]];
            const Vec3& second = mesh.vertices[face[1]];
            const Vec3& third = mesh.vertices[face[2]];
            Vec3 normal = (second - first).cross(third - first);
            const double length = normal.norm();
            if (length <= ray_epsilon * ray_epsilon) continue;
            normal /= length;
            const Vec3 centroid = (first + second + third) / 3.0;
            const bool positive_blocked = bvh.rayIntersectsAny(
                centroid + normal * ray_epsilon, normal, face_id, ray_epsilon);
            const bool negative_blocked = bvh.rayIntersectsAny(
                centroid - normal * ray_epsilon, -normal, face_id, ray_epsilon);
            open_faces[face_id] = !(positive_blocked && negative_blocked);
        }
    }

    const Mesh& mesh;
    SourceTriangleBvh bvh;
    std::vector<std::uint8_t> open_faces;
};

FilledSurfaceDistanceCertificate certifyFilledSurfaceDistance(
    const Mesh& filled_surface_mesh,
    const std::vector<OutputPrimitive>& output,
    const double maximum_distance)
{
    struct DistancePointKey
    {
        std::uint64_t x = 0;
        std::uint64_t y = 0;
        std::uint64_t z = 0;
        bool operator==(const DistancePointKey&) const = default;
    };
    struct DistancePointKeyHash
    {
        std::size_t operator()(const DistancePointKey& key) const noexcept
        {
            std::size_t seed = std::hash<std::uint64_t>{}(key.x);
            seed ^= std::hash<std::uint64_t>{}(key.y) +
                0x9e3779b97f4a7c15ULL + (seed << 6U) + (seed >> 2U);
            seed ^= std::hash<std::uint64_t>{}(key.z) +
                0x9e3779b97f4a7c15ULL + (seed << 6U) + (seed >> 2U);
            return seed;
        }
    };
    struct CachedReference
    {
        const Vec3* vertices = nullptr;
        const Face* faces = nullptr;
        std::size_t vertex_count = 0;
        std::size_t face_count = 0;
        std::unique_ptr<SourceTriangleBvh> reference;
        std::unordered_map<DistancePointKey, double, DistancePointKeyHash>
            point_distances;
    };
    static thread_local CachedReference cache;
    if (cache.vertices != filled_surface_mesh.vertices.data() ||
        cache.faces != filled_surface_mesh.faces.data() ||
        cache.vertex_count != filled_surface_mesh.vertices.size() ||
        cache.face_count != filled_surface_mesh.faces.size())
    {
        cache.vertices = filled_surface_mesh.vertices.data();
        cache.faces = filled_surface_mesh.faces.data();
        cache.vertex_count = filled_surface_mesh.vertices.size();
        cache.face_count = filled_surface_mesh.faces.size();
        cache.reference = std::make_unique<SourceTriangleBvh>(
            filled_surface_mesh);
        cache.point_distances.clear();
        cache.point_distances.reserve(65536);
    }

    FilledSurfaceDistanceCertificate result;
    if (!cache.reference)
    {
        result.certified = output.empty();
        result.exceeded = !output.empty();
        return result;
    }
    const PrimitiveMesh proxy = triangulateOutputPrimitives(output);
    constexpr int maximum_depth = 7;
    constexpr std::size_t maximum_samples = 4096;
    bool uncertain = false;
    const auto distanceAt = [&](const Vec3& point)
    {
        ++result.samples;
        const auto bits = [](const double value)
        {
            return std::bit_cast<std::uint64_t>(
                value == 0.0 ? 0.0 : value);
        };
        const DistancePointKey key{
            bits(point.x()), bits(point.y()), bits(point.z())};
        double distance = 0.0;
        const auto cached_distance = cache.point_distances.find(key);
        if (cached_distance != cache.point_distances.end())
            distance = cached_distance->second;
        else
        {
            const Vec3 closest = cache.reference->closestPoint(point).point;
            distance = (closest - point).norm();
            if (cache.point_distances.size() >= 1000000)
                cache.point_distances.clear();
            cache.point_distances.emplace(key, distance);
        }
        result.observed_maximum = std::max(
            result.observed_maximum, distance);
        if (distance > maximum_distance) result.exceeded = true;
        return distance;
    };
    std::function<bool(const Vec3&, const Vec3&, const Vec3&,
                       double, double, double, int)> certifyTriangle;
    certifyTriangle = [&](const Vec3& first, const Vec3& second,
                          const Vec3& third, const double first_distance,
                          const double second_distance,
                          const double third_distance, const int depth)
    {
        if (result.exceeded) return false;
        const double maximum_vertex_distance = std::max({
            first_distance, second_distance, third_distance});
        const double maximum_edge = std::max({
            (second - first).norm(), (third - second).norm(),
            (first - third).norm()});
        // Distance to a closed triangle mesh is 1-Lipschitz. Every point in
        // this subtriangle lies no farther than maximum_edge from at least one
        // sampled vertex, so this is a conservative continuous upper bound.
        if (maximum_vertex_distance + maximum_edge <= maximum_distance)
            return true;
        if (depth >= maximum_depth || result.samples + 3 > maximum_samples)
        {
            uncertain = true;
            return false;
        }
        const Vec3 first_second = (first + second) * 0.5;
        const Vec3 second_third = (second + third) * 0.5;
        const Vec3 third_first = (third + first) * 0.5;
        const double first_second_distance = distanceAt(first_second);
        if (result.exceeded) return false;
        const double second_third_distance = distanceAt(second_third);
        if (result.exceeded) return false;
        const double third_first_distance = distanceAt(third_first);
        if (result.exceeded) return false;
        return certifyTriangle(
                   first, first_second, third_first,
                   first_distance, first_second_distance,
                   third_first_distance, depth + 1) &&
               certifyTriangle(
                   first_second, second, second_third,
                   first_second_distance, second_distance,
                   second_third_distance, depth + 1) &&
               certifyTriangle(
                   third_first, second_third, third,
                   third_first_distance, second_third_distance,
                   third_distance, depth + 1) &&
               certifyTriangle(
                   first_second, second_third, third_first,
                   first_second_distance, second_third_distance,
                   third_first_distance, depth + 1);
    };

    for (const Face& face : proxy.faces)
    {
        const Vec3& first = proxy.vertices[face[0]];
        const Vec3& second = proxy.vertices[face[1]];
        const Vec3& third = proxy.vertices[face[2]];
        const double first_distance = distanceAt(first);
        if (result.exceeded) break;
        const double second_distance = distanceAt(second);
        if (result.exceeded) break;
        const double third_distance = distanceAt(third);
        if (result.exceeded) break;
        if (!certifyTriangle(
                first, second, third, first_distance, second_distance,
                third_distance, 0)) break;
    }
    result.certified = !result.exceeded && !uncertain;
    return result;
}

FinalOpenErrorAudit measureFilledSurfaceDistance(
    const SourceTriangleBvh& filled_surface,
    const PrimitiveMesh& proxy,
    const double sample_spacing,
    const double early_exit_distance = std::numeric_limits<double>::infinity())
{
    FinalOpenErrorAudit audit;
    double weighted_distance = 0.0;
    double total_weight = 0.0;
    bool exceeded = false;
    for (const Face& face : proxy.faces)
    {
        const Vec3& first = proxy.vertices[face[0]];
        const Vec3& second = proxy.vertices[face[1]];
        const Vec3& third = proxy.vertices[face[2]];
        const double area = 0.5 * (second - first).cross(third - first).norm();
        if (area <= 1.0e-30) continue;
        const double maximum_edge = std::max({
            (second - first).norm(), (third - second).norm(),
            (first - third).norm()});
        const int subdivisions = std::clamp(
            static_cast<int>(std::ceil(maximum_edge /
                std::max(sample_spacing, 1.0e-30))), 1, 128);
        const double sample_weight = area /
            static_cast<double>(subdivisions * subdivisions);
        const auto sample = [&](const double second_weight,
                                const double third_weight,
                                const double area_weight)
        {
            const Vec3 proxy_point = first + (second - first) * second_weight +
                (third - first) * third_weight;
            const auto closest = filled_surface.closestPoint(proxy_point);
            const double distance = (closest.point - proxy_point).norm();
            weighted_distance += distance * area_weight;
            total_weight += area_weight;
            ++audit.distance_sample_count;
            if (distance > audit.maximum_distance)
            {
                audit.maximum_distance = distance;
                audit.maximum_proxy_point = proxy_point;
                audit.maximum_source_point = closest.point;
            }
            exceeded = distance > early_exit_distance;
        };

        sample(0.0, 0.0, 0.0);
        if (exceeded) break;
        sample(1.0, 0.0, 0.0);
        if (exceeded) break;
        sample(0.0, 1.0, 0.0);
        if (exceeded) break;
        // Merges that bridge a recess or cavity commonly agree with the
        // reference along all three boundary vertices while reaching their
        // maximum error in the triangle interior. Test the centroid before the
        // edge lattice so rejected candidates can stop after four queries.
        sample(1.0 / 3.0, 1.0 / 3.0, 0.0);
        if (exceeded) break;
        for (int step = 1; step < subdivisions && !exceeded; ++step)
        {
            const double weight = static_cast<double>(step) / subdivisions;
            sample(weight, 0.0, 0.0);
            if (exceeded) break;
            sample(0.0, weight, 0.0);
            if (exceeded) break;
            sample(1.0 - weight, weight, 0.0);
        }
        for (int row = 0; row < subdivisions && !exceeded; ++row)
            for (int column = 0; column < subdivisions - row && !exceeded;
                 ++column)
            {
                sample((static_cast<double>(row) + 1.0 / 3.0) / subdivisions,
                       (static_cast<double>(column) + 1.0 / 3.0) / subdivisions,
                       sample_weight);
                if (exceeded) break;
                if (row + column + 1 < subdivisions)
                    sample((static_cast<double>(row) + 2.0 / 3.0) / subdivisions,
                           (static_cast<double>(column) + 2.0 / 3.0) / subdivisions,
                           sample_weight);
            }
        if (exceeded) break;
    }
    if (total_weight > 0.0)
        audit.mean_distance = weighted_distance / total_weight;
    audit.distance_area_integral = weighted_distance;
    return audit;
}

FinalOpenErrorAudit measureOpenSurfaceDistance(
    const SourceOpenSurfaceReference& source,
    const PrimitiveMesh& proxy,
    const double sample_spacing,
    const std::vector<CertifiedExtrusion>* ignored_planar_holes = nullptr)
{
    FinalOpenErrorAudit audit;
    double weighted_distance = 0.0;
    double total_weight = 0.0;
    for (const Face& face : proxy.faces)
    {
        const Vec3& first = proxy.vertices[face[0]];
        const Vec3& second = proxy.vertices[face[1]];
        const Vec3& third = proxy.vertices[face[2]];
        const double area = 0.5 * (second - first).cross(third - first).norm();
        if (area <= 1.0e-30) continue;
        const double maximum_edge = std::max({
            (second - first).norm(), (third - second).norm(),
            (first - third).norm()});
        const int subdivisions = std::clamp(
            static_cast<int>(std::ceil(maximum_edge /
                std::max(sample_spacing, 1.0e-30))), 1, 128);
        const double sample_weight = area /
            static_cast<double>(subdivisions * subdivisions);
        const auto sample = [&](const double second_weight,
                                 const double third_weight,
                                 const double area_weight)
        {
            const Vec3 proxy_point = first + (second - first) * second_weight +
                (third - first) * third_weight;
            if (ignored_planar_holes != nullptr &&
                std::any_of(
                    ignored_planar_holes->begin(),
                    ignored_planar_holes->end(),
                    [&](const CertifiedExtrusion& hole)
                    {
                        const Vec3 local = hole.frame.transposeMultiply(
                            proxy_point - hole.origin);
                        return std::abs(local.z()) <= hole.upper_distance &&
                               pointInPolygon(
                                   Vec2(local.x(), local.y()), hole.boundary,
                                   hole.lower_distance);
                    }))
                return;
            const auto closest = source.bvh.closestPoint(proxy_point);
            if (closest.face_id >= source.open_faces.size() ||
                !source.open_faces[closest.face_id]) return;
            const double distance = (closest.point - proxy_point).norm();
            // Boundary samples participate in the maximum but carry no area
            // weight. This mirrors surface-Hausdorff tools that sample mesh
            // vertices, edges, and faces while keeping the mean surface-area
            // weighted.
            weighted_distance += distance * area_weight;
            total_weight += area_weight;
            ++audit.distance_sample_count;
            if (distance > audit.maximum_distance)
            {
                audit.maximum_distance = distance;
                audit.maximum_proxy_point = proxy_point;
                audit.maximum_source_point = closest.point;
            }
        };

        // A face-interior grid alone can miss the actual sampled maximum at a
        // proxy corner or along a polygon boundary. Include all three vertices
        // and deterministic, uniformly spaced edge samples before sampling the
        // face interiors. Adjacent proxy faces may repeat a boundary point;
        // that is harmless for the maximum and avoids topology assumptions.
        sample(0.0, 0.0, 0.0);
        sample(1.0, 0.0, 0.0);
        sample(0.0, 1.0, 0.0);
        for (int step = 1; step < subdivisions; ++step)
        {
            const double weight = static_cast<double>(step) / subdivisions;
            sample(weight, 0.0, 0.0);
            sample(0.0, weight, 0.0);
            sample(1.0 - weight, weight, 0.0);
        }
        for (int row = 0; row < subdivisions; ++row)
            for (int column = 0; column < subdivisions - row; ++column)
            {
                sample((static_cast<double>(row) + 1.0 / 3.0) / subdivisions,
                       (static_cast<double>(column) + 1.0 / 3.0) / subdivisions,
                       sample_weight);
                if (row + column + 1 < subdivisions)
                    sample((static_cast<double>(row) + 2.0 / 3.0) / subdivisions,
                           (static_cast<double>(column) + 2.0 / 3.0) / subdivisions,
                           sample_weight);
            }
    }
    if (total_weight > 0.0)
        audit.mean_distance = weighted_distance / total_weight;
    audit.distance_area_integral = weighted_distance;
    return audit;
}

std::vector<CertifiedExtrusion> certifiedPlanarHoleExclusions(
    const Mesh& source_mesh,
    const std::vector<OutputPrimitive>& output,
    const double tolerance)
{
    std::vector<CertifiedExtrusion> exclusions;
    for (const OutputPrimitive& item : output)
    {
        if (item.source_faces.size() <= 1) continue;
        Vec3 normal;
        if (!planarNormal(item.primitive, tolerance, normal)) continue;
        const Vec3 origin = planarPoint(item.primitive);
        bool coplanar = true;
        for (const auto face_id : item.source_faces)
        {
            if (face_id >= source_mesh.faces.size())
            {
                coplanar = false;
                break;
            }
            for (const auto vertex_id : source_mesh.faces[face_id])
                if (std::abs((source_mesh.vertices[vertex_id] - origin).dot(normal)) >
                    tolerance)
                {
                    coplanar = false;
                    break;
                }
            if (!coplanar) break;
        }
        if (!coplanar) continue;

        const auto loops = boundaryLoops(source_mesh, item.source_faces);
        if (loops.size() <= 1) continue;
        const Mat3 basis = orthonormalFrame(normal);
        Mat3 frame;
        frame.col(0) = static_cast<Vec3>(basis.col(1));
        frame.col(1) = static_cast<Vec3>(basis.col(2));
        frame.col(2) = normal;
        std::vector<std::vector<Vec2>> polygons;
        polygons.reserve(loops.size());
        std::size_t outer = 0;
        double outer_area = 0.0;
        for (const auto& loop : loops)
        {
            std::vector<Vec2> polygon;
            polygon.reserve(loop.size());
            for (const auto vertex_id : loop)
            {
                const Vec3 local = frame.transposeMultiply(
                    source_mesh.vertices[vertex_id] - origin);
                polygon.emplace_back(local.x(), local.y());
            }
            polygon = simplifyPolygon(std::move(polygon), tolerance);
            const double area = std::abs(signedArea(polygon));
            if (area > outer_area)
            {
                outer_area = area;
                outer = polygons.size();
            }
            polygons.push_back(std::move(polygon));
        }
        if (polygons[outer].size() < 3) continue;
        for (std::size_t index = 0; index < polygons.size(); ++index)
        {
            if (index == outer || polygons[index].size() < 3 ||
                std::abs(signedArea(polygons[index])) <= tolerance * tolerance)
                continue;
            Vec2 centroid = Vec2::Zero();
            for (const Vec2& point : polygons[index]) centroid += point;
            centroid /= static_cast<double>(polygons[index].size());
            if (!pointInPolygon(centroid, polygons[outer], tolerance)) continue;
            CertifiedExtrusion exclusion;
            exclusion.origin = origin;
            exclusion.frame = frame;
            exclusion.boundary = std::move(polygons[index]);
            // Reuse the scalar fields as the certified plane and 2D tolerances;
            // these records never participate in volume/occlusion processing.
            exclusion.lower_distance = tolerance;
            exclusion.upper_distance = tolerance;
            exclusions.push_back(std::move(exclusion));
        }
    }
    return exclusions;
}

FinalOpenErrorAudit measureChargeableOpenSurfaceDistance(
    const SourceOpenSurfaceReference& source,
    const Mesh& source_mesh,
    const std::vector<OutputPrimitive>& output,
    const double sample_spacing,
    const double tolerance,
    const std::vector<OutputPrimitive>* planar_hole_certificates = nullptr)
{
    const auto exclusions = certifiedPlanarHoleExclusions(
        source_mesh,
        planar_hole_certificates ? *planar_hole_certificates : output,
        tolerance);
    return measureOpenSurfaceDistance(
        source, triangulateOutputPrimitives(output), sample_spacing,
        exclusions.empty() ? nullptr : &exclusions);
}

double maximumFilledSurfaceDistance(
    const Mesh& filled_surface_mesh,
    const std::vector<OutputPrimitive>& output,
    const double sample_spacing,
    const double maximum_distance)
{
    struct CachedReference
    {
        const Vec3* vertices = nullptr;
        const Face* faces = nullptr;
        std::size_t vertex_count = 0;
        std::size_t face_count = 0;
        std::unique_ptr<SourceTriangleBvh> reference;
    };
    static thread_local CachedReference cache;
    if (cache.vertices != filled_surface_mesh.vertices.data() ||
        cache.faces != filled_surface_mesh.faces.data() ||
        cache.vertex_count != filled_surface_mesh.vertices.size() ||
        cache.face_count != filled_surface_mesh.faces.size())
    {
        cache.vertices = filled_surface_mesh.vertices.data();
        cache.faces = filled_surface_mesh.faces.data();
        cache.vertex_count = filled_surface_mesh.vertices.size();
        cache.face_count = filled_surface_mesh.faces.size();
        cache.reference = std::make_unique<SourceTriangleBvh>(
            filled_surface_mesh);
    }
    return measureFilledSurfaceDistance(
        *cache.reference, triangulateOutputPrimitives(output),
        sample_spacing, maximum_distance).maximum_distance;
}

std::vector<OutputPrimitive> mergeSpatialSurfaceGroups(
    const Mesh& source_mesh,
    const Mesh& filled_surface_mesh,
    std::vector<OutputPrimitive> primitives,
    const double tolerance,
    const double maximum_error_distance,
    const double error_sample_spacing,
    std::size_t& merged_group_count,
    std::vector<BoxFit>& accepted_boxes,
    const std::filesystem::path& profile_path)
{
    merged_group_count = 0;
    if (primitives.size() < 2) return primitives;
    const auto started = std::chrono::steady_clock::now();
    struct Record
    {
        std::size_t primitive = 0;
        Bounds bounds;
        Vec3 center = Vec3::Zero();
        std::size_t triangles = 0;
        bool polygon_surface = false;
    };
    std::vector<Record> records;
    records.reserve(primitives.size());
    for (std::size_t index = 0; index < primitives.size(); ++index)
    {
        const PrimitiveMesh surface = triangulatePrimitive(
            primitives[index].primitive);
        Bounds bounds;
        for (const Vec3& vertex : surface.vertices)
        {
            bounds.lower = bounds.lower.cwiseMin(vertex);
            bounds.upper = bounds.upper.cwiseMax(vertex);
        }
        const Kind kind = primitives[index].primitive.kind;
        records.push_back({index, bounds, (bounds.lower + bounds.upper) * 0.5,
                           surface.faces.size(),
                           kind == Kind::Polygon || kind == Kind::Rectangle ||
                           kind == Kind::Triangle});
    }
    std::vector<std::size_t> order(records.size());
    std::iota(order.begin(), order.end(), 0);
    std::uint64_t next_enclosure_group = 1;
    for (const OutputPrimitive& primitive : primitives)
        next_enclosure_group = std::max(
            next_enclosure_group, primitive.enclosure_group + 1);

    struct Profile
    {
        std::size_t nodes = 0;
        std::size_t box_candidates = 0;
        std::size_t analytic_rejections = 0;
        std::size_t workload_rejections = 0;
        std::size_t error_rejections = 0;
        std::size_t accepted_boxes = 0;
        std::size_t certificate_samples = 0;
        double distance_seconds = 0.0;
    } profile;
    const auto writeProfile = [&](const bool complete)
    {
        std::ofstream stream(profile_path);
        stream << std::setprecision(17)
               << "{\"complete\":" << (complete ? "true" : "false")
               << ",\"surface_candidates_only\":true"
               << ",\"input_primitives\":" << primitives.size()
               << ",\"hierarchy_nodes\":" << profile.nodes
               << ",\"box_surface_set_candidates\":"
               << profile.box_candidates
               << ",\"analytic_rejections\":"
               << profile.analytic_rejections
               << ",\"workload_rejections\":"
               << profile.workload_rejections
               << ",\"error_rejections\":" << profile.error_rejections
               << ",\"accepted_box_surface_sets\":"
               << profile.accepted_boxes
               << ",\"certificate_samples\":"
               << profile.certificate_samples
               << ",\"distance_seconds\":" << profile.distance_seconds
               << ",\"total_seconds\":" << std::chrono::duration<double>(
                    std::chrono::steady_clock::now() - started).count()
               << "}\n";
    };
    writeProfile(false);

    const auto solve = [&](auto&& self, const std::size_t begin,
                           const std::size_t end)
        -> std::vector<OutputPrimitive>
    {
        ++profile.nodes;
        const std::size_t count = end - begin;
        if (count == 1)
            return {primitives[records[order[begin]].primitive]};

        std::size_t original_triangles = 0;
        bool polygon_surfaces_only = true;
        Bounds center_bounds;
        for (std::size_t offset = begin; offset < end; ++offset)
        {
            const Record& record = records[order[offset]];
            center_bounds.lower = center_bounds.lower.cwiseMin(record.center);
            center_bounds.upper = center_bounds.upper.cwiseMax(record.center);
            original_triangles += record.triangles;
            polygon_surfaces_only &= record.polygon_surface;
        }

        if (!polygon_surfaces_only)
            ++profile.analytic_rejections;
        else if (original_triangles <= 12)
            ++profile.workload_rejections;
        else
        {
            Mesh fitting_mesh;
            std::vector<std::uint32_t> fitting_vertices;
            std::vector<std::uint32_t> responsibility;
            for (std::size_t offset = begin; offset < end; ++offset)
            {
                const OutputPrimitive& item =
                    primitives[records[order[offset]].primitive];
                const PrimitiveMesh surface =
                    triangulatePrimitive(item.primitive);
                for (const Vec3& vertex : surface.vertices)
                {
                    fitting_vertices.push_back(static_cast<std::uint32_t>(
                        fitting_mesh.vertices.size()));
                    fitting_mesh.vertices.push_back(vertex);
                }
                responsibility.insert(responsibility.end(),
                                      item.source_faces.begin(),
                                      item.source_faces.end());
            }
            std::sort(responsibility.begin(), responsibility.end());
            responsibility.erase(
                std::unique(responsibility.begin(), responsibility.end()),
                responsibility.end());
            BoxFit box = fitBox(fitting_mesh, fitting_vertices);
            Vec3 local_lower = -box.half_size;
            Vec3 local_upper = box.half_size;
            for (const std::uint32_t face_id : responsibility)
            {
                if (face_id >= source_mesh.faces.size()) continue;
                for (const std::uint32_t vertex_id : source_mesh.faces[face_id])
                {
                    if (vertex_id >= source_mesh.vertices.size()) continue;
                    const Vec3 local = box.axes.transposeMultiply(
                        source_mesh.vertices[vertex_id] - box.center);
                    local_lower = local_lower.cwiseMin(local);
                    local_upper = local_upper.cwiseMax(local);
                }
            }
            const Vec3 local_center = (local_lower + local_upper) * 0.5;
            box.center += box.axes * local_center;
            box.half_size = (local_upper - local_lower) * 0.5;
            box.volume = 8.0 * box.half_size.prod();
            const Vec3 extent = box.half_size * 2.0;
            if (extent.x() > tolerance && extent.y() > tolerance &&
                extent.z() > tolerance)
            {
                std::vector<OutputPrimitive> candidate;
                appendBoxRectangles(candidate, box, responsibility, -1, 0.0,
                                    next_enclosure_group);
                ++profile.box_candidates;
                const auto distance_started = std::chrono::steady_clock::now();
                const FilledSurfaceDistanceCertificate certificate =
                    certifyFilledSurfaceDistance(
                        filled_surface_mesh, candidate,
                        maximum_error_distance + tolerance);
                profile.certificate_samples += certificate.samples;
                bool accepted = certificate.certified && !certificate.exceeded;
                if (!accepted && !certificate.exceeded)
                    accepted = maximumFilledSurfaceDistance(
                        filled_surface_mesh, candidate, error_sample_spacing,
                        maximum_error_distance + tolerance) <=
                        maximum_error_distance + tolerance;
                profile.distance_seconds += std::chrono::duration<double>(
                    std::chrono::steady_clock::now() - distance_started).count();
                if (accepted)
                {
                    ++next_enclosure_group;
                    ++profile.accepted_boxes;
                    ++merged_group_count;
                    accepted_boxes.push_back(box);
                    return candidate;
                }
                ++profile.error_rejections;
            }
        }

        const Vec3 center_extent = center_bounds.upper - center_bounds.lower;
        int axis = 0;
        if (center_extent.y() > center_extent.x()) axis = 1;
        if (center_extent.z() > center_extent[axis]) axis = 2;
        const std::size_t middle = begin + count / 2;
        std::nth_element(order.begin() + begin, order.begin() + middle,
                         order.begin() + end,
            [&](const std::size_t first, const std::size_t second)
            {
                return records[first].center[axis] < records[second].center[axis];
            });
        auto left = self(self, begin, middle);
        auto right = self(self, middle, end);
        left.insert(left.end(), std::make_move_iterator(right.begin()),
                    std::make_move_iterator(right.end()));
        return left;
    };

    auto result = solve(solve, 0, order.size());
    writeProfile(true);
    return result;
}

std::vector<OutputPrimitive> mergeSupportProtrusionSurfaceSets(
    const Mesh& source_mesh,
    const Mesh& filled_surface_mesh,
    std::vector<OutputPrimitive> primitives,
    const PrimitiveMeshAnalysisOptions& options,
    const double model_diagonal,
    const double model_surface_area,
    const double tolerance,
    const double maximum_error_distance,
    const double error_sample_spacing,
    std::size_t& accepted_count,
    std::vector<OutputPrimitive>& occluded_support_certificates,
    const std::filesystem::path& profile_path)
{
    const auto started = std::chrono::steady_clock::now();
    accepted_count = 0;
    std::vector<bool> represented_faces(source_mesh.faces.size(), false);
    for (const OutputPrimitive& primitive : primitives)
        for (const auto face : primitive.source_faces)
            if (face < represented_faces.size()) represented_faces[face] = true;

    std::vector<std::unordered_set<std::uint32_t>> adjacency;
    const auto clusters = coplanarClusters(
        source_mesh, model_diagonal * options.coplanar_relative_tolerance,
        adjacency, &represented_faces);
    std::vector<bool> claimed_clusters(clusters.size(), false);
    PrimitiveMeshAnalysisOptions candidate_options = options;
    // Candidate generation may be permissive; acceptance below is controlled
    // solely by the user's directed-distance limit and triangle workload.
    candidate_options.recognize_support_protrusions = true;
    candidate_options.protrusion_max_area_excess_ratio =
        std::numeric_limits<double>::infinity();
    double minimum_area_excess = -1.0;
    const auto candidates = recognizeSupportProtrusions(
        source_mesh, clusters, adjacency, candidate_options,
        model_diagonal, model_surface_area, claimed_clusters,
        minimum_area_excess);

    struct Profile
    {
        std::size_t generated = 0;
        std::size_t overlap_rejections = 0;
        std::size_t ownership_rejections = 0;
        std::size_t workload_rejections = 0;
        std::size_t error_rejections = 0;
        std::size_t accepted = 0;
        std::size_t consumed_primitives = 0;
        std::size_t certificate_samples = 0;
        std::size_t support_clipped_candidates = 0;
        std::size_t support_clipped_faces = 0;
        std::size_t support_discarded_faces = 0;
        double maximum_removed_penetration = 0.0;
        std::size_t contact_clipped_primitives = 0;
        std::size_t contact_removed_primitives = 0;
        std::size_t contact_output_fragments = 0;
        double contact_removed_area = 0.0;
        double distance_seconds = 0.0;
    } profile;
    profile.generated = candidates.size();
    std::vector<bool> consumed(primitives.size(), false);
    std::vector<bool> claimed_faces(source_mesh.faces.size(), false);
    std::vector<OutputPrimitive> replacements;
    std::vector<RecognizedProtrusion> accepted_protrusions;
    for (const RecognizedProtrusion& candidate : candidates)
    {
        if (std::any_of(candidate.faces.begin(), candidate.faces.end(),
            [&](const auto face)
            { return face < claimed_faces.size() && claimed_faces[face]; }))
        {
            ++profile.overlap_rejections;
            continue;
        }
        std::unordered_set<std::uint32_t> face_set(
            candidate.faces.begin(), candidate.faces.end());
        std::vector<std::size_t> owners;
        std::size_t original_triangles = 0;
        for (std::size_t index = 0; index < primitives.size(); ++index)
        {
            if (consumed[index] || primitives[index].source_faces.empty()) continue;
            const bool has_candidate_face = std::any_of(
                primitives[index].source_faces.begin(),
                primitives[index].source_faces.end(),
                [&](const auto face) { return face_set.contains(face); });
            if (!has_candidate_face) continue;
            const bool fully_owned = std::all_of(
                primitives[index].source_faces.begin(),
                primitives[index].source_faces.end(),
                [&](const auto face) { return face_set.contains(face); });
            if (!fully_owned) continue;
            owners.push_back(index);
            original_triangles += triangulatePrimitive(
                primitives[index].primitive).faces.size();
        }
        if (owners.empty())
        {
            ++profile.ownership_rejections;
            continue;
        }
        std::vector<OutputPrimitive> surface_set;
        SupportClippedBoxStats clip_stats;
        appendSupportClippedBoxSurfaces(
            surface_set, candidate, tolerance, clip_stats);
        if (clip_stats.clipped_faces != 0 ||
            clip_stats.discarded_faces != 0)
            ++profile.support_clipped_candidates;
        profile.support_clipped_faces += clip_stats.clipped_faces;
        profile.support_discarded_faces += clip_stats.discarded_faces;
        profile.maximum_removed_penetration = std::max(
            profile.maximum_removed_penetration,
            clip_stats.maximum_removed_penetration);
        if (surface_set.empty())
        {
            ++profile.ownership_rejections;
            continue;
        }
        const std::size_t candidate_triangles =
            triangulatedFaceCount(surface_set);
        if (candidate_triangles >= original_triangles)
        {
            ++profile.workload_rejections;
            continue;
        }
        const auto distance_started = std::chrono::steady_clock::now();
        const FilledSurfaceDistanceCertificate certificate =
            certifyFilledSurfaceDistance(
                filled_surface_mesh, surface_set,
                maximum_error_distance + tolerance);
        profile.certificate_samples += certificate.samples;
        bool accepted = certificate.certified && !certificate.exceeded;
        if (!accepted && !certificate.exceeded)
            accepted = maximumFilledSurfaceDistance(
                filled_surface_mesh, surface_set, error_sample_spacing,
                maximum_error_distance + tolerance) <=
                maximum_error_distance + tolerance;
        profile.distance_seconds += std::chrono::duration<double>(
            std::chrono::steady_clock::now() - distance_started).count();
        if (!accepted)
        {
            ++profile.error_rejections;
            continue;
        }
        for (const auto owner : owners) consumed[owner] = true;
        for (const auto face : candidate.faces)
            if (face < claimed_faces.size()) claimed_faces[face] = true;
        profile.consumed_primitives += owners.size();
        ++profile.accepted;
        ++accepted_count;
        accepted_protrusions.push_back(candidate);
        replacements.insert(
            replacements.end(),
            std::make_move_iterator(surface_set.begin()),
            std::make_move_iterator(surface_set.end()));
    }

    std::vector<OutputPrimitive> result;
    result.reserve(primitives.size() + replacements.size());
    for (std::size_t index = 0; index < primitives.size(); ++index)
        if (!consumed[index]) result.push_back(std::move(primitives[index]));
    result.insert(result.end(),
                  std::make_move_iterator(replacements.begin()),
                  std::make_move_iterator(replacements.end()));
    SupportContactClipStats contact_stats;
    result = clipSupportContactOcclusion(
        std::move(result), accepted_protrusions, tolerance,
        occluded_support_certificates, contact_stats);
    profile.contact_clipped_primitives = contact_stats.clipped_primitives;
    profile.contact_removed_primitives = contact_stats.removed_primitives;
    profile.contact_output_fragments = contact_stats.output_fragments;
    profile.contact_removed_area = contact_stats.removed_area;
    std::ofstream stream(profile_path);
    stream << std::setprecision(17)
           << "{\"generated_candidates\":" << profile.generated
           << ",\"overlap_rejections\":" << profile.overlap_rejections
           << ",\"ownership_rejections\":" << profile.ownership_rejections
           << ",\"workload_rejections\":" << profile.workload_rejections
           << ",\"error_rejections\":" << profile.error_rejections
           << ",\"accepted\":" << profile.accepted
           << ",\"consumed_primitives\":" << profile.consumed_primitives
           << ",\"certificate_samples\":" << profile.certificate_samples
           << ",\"support_clipped_candidates\":"
           << profile.support_clipped_candidates
           << ",\"support_clipped_faces\":" << profile.support_clipped_faces
           << ",\"support_discarded_faces\":"
           << profile.support_discarded_faces
           << ",\"maximum_removed_penetration\":"
           << profile.maximum_removed_penetration
           << ",\"contact_clipped_primitives\":"
           << profile.contact_clipped_primitives
           << ",\"contact_removed_primitives\":"
           << profile.contact_removed_primitives
           << ",\"contact_output_fragments\":"
           << profile.contact_output_fragments
           << ",\"contact_removed_area\":"
           << profile.contact_removed_area
           << ",\"distance_seconds\":" << profile.distance_seconds
           << ",\"accepted_geometry\":[";
    for (std::size_t index = 0; index < accepted_protrusions.size(); ++index)
    {
        if (index != 0) stream << ',';
        const auto& item = accepted_protrusions[index];
        stream << "{\"footprint_ratio\":" << item.footprint_ratio
               << ",\"group_depth\":" << item.group_depth
               << ",\"box_extent\":["
               << item.box.half_size.x() * 2.0 << ','
               << item.box.half_size.y() * 2.0 << ','
               << item.box.half_size.z() * 2.0 << "]}";
    }
    stream << "]"
           << ",\"output_primitives\":" << result.size()
           << ",\"total_seconds\":" << std::chrono::duration<double>(
                std::chrono::steady_clock::now() - started).count()
           << "}\n";
    return result;
}

std::vector<OutputPrimitive> fillCertifiedIntercomponentGaps(
    const Mesh& mesh,
    std::vector<OutputPrimitive> primitives,
    const double model_volume,
    const double maximum_gap_volume_ratio,
    const double tolerance,
    std::size_t& filled_gap_count,
    std::vector<BoxFit>& certified_volumes,
    const std::filesystem::path& profile_path)
{
    filled_gap_count = 0;
    if (mesh.faces.empty() || maximum_gap_volume_ratio <= 0.0)
        return primitives;
    std::vector<std::uint32_t> all_vertices(mesh.vertices.size());
    std::iota(all_vertices.begin(), all_vertices.end(), 0);
    const BoxFit model_box = fitBox(mesh, all_vertices);
    const Mat3 frame = model_box.axes;

    struct ComponentBox
    {
        Bounds local;
        BoxFit world;
        std::vector<std::uint32_t> faces;
    };
    std::vector<ComponentBox> components;
    std::unordered_map<std::uint64_t, std::vector<std::uint32_t>> edge_faces;
    edge_faces.reserve(mesh.faces.size() * 3);
    for (std::uint32_t face_id = 0; face_id < mesh.faces.size(); ++face_id)
    {
        const Face& face = mesh.faces[face_id];
        for (int edge = 0; edge < 3; ++edge)
            edge_faces[edgeKey(face[edge], face[(edge + 1) % 3])]
                .push_back(face_id);
    }
    std::vector<std::vector<std::uint32_t>> manifold_neighbors(mesh.faces.size());
    for (const auto& [edge, incident] : edge_faces)
    {
        (void)edge;
        if (incident.size() != 2) continue;
        manifold_neighbors[incident[0]].push_back(incident[1]);
        manifold_neighbors[incident[1]].push_back(incident[0]);
    }
    std::vector<std::uint8_t> visited(mesh.faces.size(), 0);
    std::vector<std::vector<std::uint32_t>> manifold_components;
    for (std::uint32_t seed = 0; seed < mesh.faces.size(); ++seed)
    {
        if (visited[seed]) continue;
        manifold_components.emplace_back();
        std::queue<std::uint32_t> queue;
        queue.push(seed);
        visited[seed] = 1;
        while (!queue.empty())
        {
            const auto face = queue.front();
            queue.pop();
            manifold_components.back().push_back(face);
            for (const auto neighbor : manifold_neighbors[face])
                if (!visited[neighbor])
                {
                    visited[neighbor] = 1;
                    queue.push(neighbor);
                }
        }
    }
    for (auto faces : manifold_components)
    {
        if (faces.size() < 12) continue;
        std::unordered_map<std::uint64_t, std::uint32_t> edge_counts;
        for (const auto face_id : faces)
        {
            const Face& face = mesh.faces[face_id];
            for (int edge = 0; edge < 3; ++edge)
                ++edge_counts[edgeKey(face[edge], face[(edge + 1) % 3])];
        }
        if (std::any_of(edge_counts.begin(), edge_counts.end(),
            [](const auto& edge) { return edge.second != 2; }))
            continue;
        const auto vertex_ids = uniqueVertices(mesh, faces);
        Bounds local;
        for (const auto vertex_id : vertex_ids)
        {
            const Vec3 point = frame.transposeMultiply(mesh.vertices[vertex_id]);
            local.lower = local.lower.cwiseMin(point);
            local.upper = local.upper.cwiseMax(point);
        }
        const Vec3 extent = local.upper - local.lower;
        const double volume = extent.cwiseMax(Vec3::Zero()).prod();
        if (volume <= std::max(model_volume * 1.0e-4,
                               tolerance * tolerance * tolerance))
            continue;
        const double minimum_extent = std::min({extent.x(), extent.y(), extent.z()});
        if (minimum_extent <= tolerance) continue;
        double maximum_inward = 0.0;
        for (const auto vertex_id : vertex_ids)
        {
            const Vec3 point = frame.transposeMultiply(mesh.vertices[vertex_id]);
            maximum_inward = std::max(maximum_inward, std::min({
                point.x() - local.lower.x(), local.upper.x() - point.x(),
                point.y() - local.lower.y(), local.upper.y() - point.y(),
                point.z() - local.lower.z(), local.upper.z() - point.z()}));
        }
        if (maximum_inward > 0.1 * minimum_extent + tolerance) continue;
        BoxFit box;
        box.axes = frame;
        box.center = frame * ((local.lower + local.upper) * 0.5);
        box.half_size = extent * 0.5;
        box.volume = volume;
        components.push_back({local, box, std::move(faces)});
        certified_volumes.push_back(box);
    }
    if (components.size() < 2) return primitives;

    std::vector<BoxFit> bridges;
    std::set<std::array<std::int64_t, 6>> bridge_keys;
    const double quantization = 1.0 / std::max(tolerance * 16.0, 1.0e-9);
    for (int separation_axis = 0; separation_axis < 3; ++separation_axis)
    {
        std::vector<std::size_t> order(components.size());
        std::iota(order.begin(), order.end(), 0);
        std::sort(order.begin(), order.end(), [&](const auto first,
                                                  const auto second)
        {
            return components[first].local.lower[separation_axis] <
                   components[second].local.lower[separation_axis];
        });
        const int first_tangent = (separation_axis + 1) % 3;
        const int second_tangent = (separation_axis + 2) % 3;
        for (std::size_t position = 0; position < order.size(); ++position)
            for (std::size_t next = position + 1;
                 next < std::min(order.size(), position + 9); ++next)
            {
                const Bounds& first = components[order[position]].local;
                const Bounds& second = components[order[next]].local;
                const double gap = second.lower[separation_axis] -
                                   first.upper[separation_axis];
                if (gap <= tolerance) continue;
                const double first_extent_a =
                    first.upper[first_tangent] - first.lower[first_tangent];
                const double second_extent_a =
                    second.upper[first_tangent] - second.lower[first_tangent];
                const double first_extent_b =
                    first.upper[second_tangent] - first.lower[second_tangent];
                const double second_extent_b =
                    second.upper[second_tangent] - second.lower[second_tangent];
                if (std::max(first_extent_a, second_extent_a) >
                        1.25 * std::min(first_extent_a, second_extent_a) ||
                    std::max(first_extent_b, second_extent_b) >
                        1.25 * std::min(first_extent_b, second_extent_b))
                    continue;
                const double lower_a = std::max(
                    first.lower[first_tangent], second.lower[first_tangent]);
                const double upper_a = std::min(
                    first.upper[first_tangent], second.upper[first_tangent]);
                const double lower_b = std::max(
                    first.lower[second_tangent], second.lower[second_tangent]);
                const double upper_b = std::min(
                    first.upper[second_tangent], second.upper[second_tangent]);
                const double overlap_a = upper_a - lower_a;
                const double overlap_b = upper_b - lower_b;
                if (overlap_a < 0.8 * std::min(first_extent_a, second_extent_a) ||
                    overlap_b < 0.8 * std::min(first_extent_b, second_extent_b))
                    continue;
                const double gap_volume = gap * overlap_a * overlap_b;
                if (gap_volume / std::max(model_volume, 1.0e-30) >
                    maximum_gap_volume_ratio)
                    continue;
                Vec3 lower = Vec3::Zero();
                Vec3 upper = Vec3::Zero();
                lower[separation_axis] = first.upper[separation_axis];
                upper[separation_axis] = second.lower[separation_axis];
                lower[first_tangent] = lower_a;
                upper[first_tangent] = upper_a;
                lower[second_tangent] = lower_b;
                upper[second_tangent] = upper_b;
                const std::array<std::int64_t, 6> key{{
                    static_cast<std::int64_t>(std::llround(lower.x() * quantization)),
                    static_cast<std::int64_t>(std::llround(lower.y() * quantization)),
                    static_cast<std::int64_t>(std::llround(lower.z() * quantization)),
                    static_cast<std::int64_t>(std::llround(upper.x() * quantization)),
                    static_cast<std::int64_t>(std::llround(upper.y() * quantization)),
                    static_cast<std::int64_t>(std::llround(upper.z() * quantization))}};
                if (!bridge_keys.insert(key).second) continue;
                BoxFit bridge;
                bridge.axes = frame;
                bridge.center = frame * ((lower + upper) * 0.5);
                bridge.half_size = (upper - lower) * 0.5;
                bridge.volume = gap_volume;
                bridges.push_back(bridge);
            }
    }
    if (bridges.empty()) return primitives;

    // Several independently matched component pairs can describe the same
    // filled space at different cross-section widths. Keep the maximal bridge
    // whenever another candidate is completely contained in it. This is a
    // geometric canonicalization, not a model-specific count rule, and avoids
    // exporting coincident internal bridge faces.
    const std::size_t raw_bridge_count = bridges.size();
    std::vector<Bounds> bridge_bounds(bridges.size());
    for (std::size_t index = 0; index < bridges.size(); ++index)
    {
        const Vec3 center = frame.transposeMultiply(bridges[index].center);
        bridge_bounds[index] = {
            center - bridges[index].half_size,
            center + bridges[index].half_size};
    }
    std::vector<std::uint8_t> redundant_bridge(bridges.size(), 0);
    for (std::size_t candidate = 0; candidate < bridges.size(); ++candidate)
        for (std::size_t container = 0; container < bridges.size(); ++container)
        {
            if (candidate == container) continue;
            bool contained = true;
            bool strictly_smaller = false;
            for (int axis = 0; axis < 3; ++axis)
            {
                contained &= bridge_bounds[candidate].lower[axis] >=
                                 bridge_bounds[container].lower[axis] - tolerance &&
                             bridge_bounds[candidate].upper[axis] <=
                                 bridge_bounds[container].upper[axis] + tolerance;
                strictly_smaller |= bridge_bounds[candidate].lower[axis] >
                                        bridge_bounds[container].lower[axis] + tolerance ||
                                    bridge_bounds[candidate].upper[axis] <
                                        bridge_bounds[container].upper[axis] - tolerance;
            }
            if (contained && strictly_smaller)
            {
                redundant_bridge[candidate] = 1;
                break;
            }
        }
    std::vector<BoxFit> canonical_bridges;
    canonical_bridges.reserve(bridges.size());
    for (std::size_t index = 0; index < bridges.size(); ++index)
        if (!redundant_bridge[index])
            canonical_bridges.push_back(bridges[index]);
    bridges = std::move(canonical_bridges);

    {
        std::ofstream profile(profile_path);
        profile << std::setprecision(17)
                << "{\"certified_components\":" << components.size()
                << ",\"raw_bridge_candidates\":" << raw_bridge_count
                << ",\"canonical_bridges\":" << bridges.size()
                << ",\"components\":[";
        for (std::size_t index = 0; index < components.size(); ++index)
        {
            if (index) profile << ',';
            const Bounds& bounds = components[index].local;
            profile << "{\"lower\":[" << bounds.lower.x() << ','
                    << bounds.lower.y() << ',' << bounds.lower.z()
                    << "],\"upper\":[" << bounds.upper.x() << ','
                    << bounds.upper.y() << ',' << bounds.upper.z() << "]}";
        }
        profile << "],\"bridges\":[";
        for (std::size_t index = 0; index < bridges.size(); ++index)
        {
            if (index) profile << ',';
            const BoxFit& bridge = bridges[index];
            const Vec3 center = frame.transposeMultiply(bridge.center);
            const Vec3 lower = center - bridge.half_size;
            const Vec3 upper = center + bridge.half_size;
            profile << "{\"lower\":[" << lower.x() << ',' << lower.y() << ','
                    << lower.z() << "],\"upper\":[" << upper.x() << ','
                    << upper.y() << ',' << upper.z()
                    << "],\"volume_ratio\":"
                    << bridge.volume / std::max(model_volume, 1.0e-30) << '}';
        }
        profile << "]}\n";
    }

    std::uint64_t next_group = 1;
    for (const OutputPrimitive& primitive : primitives)
        next_group = std::max(next_group, primitive.enclosure_group + 1);
    const std::size_t first_bridge_primitive = primitives.size();
    for (const BoxFit& bridge : bridges)
    {
        appendBoxRectangles(primitives, bridge, {}, -1, 0.0, next_group++);
        certified_volumes.push_back(bridge);
    }
    std::size_t removed = 0;
    std::vector<std::uint32_t> absorbed_faces;
    primitives = clipPlanarOcclusionByClosedVolumes(
        std::move(primitives), certified_volumes, {}, tolerance, removed,
        &absorbed_faces, false);
    if (!absorbed_faces.empty())
    {
        for (OutputPrimitive& primitive : primitives)
            if (primitive.enclosure_group >= next_group - bridges.size())
            {
                primitive.source_faces.insert(
                    primitive.source_faces.end(), absorbed_faces.begin(),
                    absorbed_faces.end());
                break;
            }
    }
    (void)first_bridge_primitive;
    filled_gap_count = bridges.size();
    return primitives;
}

void writeOpenErrorVisualization(
    const std::filesystem::path& path,
    const FinalOpenErrorAudit& audit)
{
    std::ofstream stream(path);
    stream << std::setprecision(17)
           << "{\"boundary_points\":[],\"maximum_pair\":{\"proxy\":["
           << audit.maximum_proxy_point.x() << ','
           << audit.maximum_proxy_point.y() << ','
           << audit.maximum_proxy_point.z() << "],\"source\":["
           << audit.maximum_source_point.x() << ','
           << audit.maximum_source_point.y() << ','
           << audit.maximum_source_point.z() << "],\"distance\":"
           << audit.maximum_distance << "}}\n";
}

std::vector<RecognizedProtrusion> recognizeLowErrorBoxEnvelopes(
    const Mesh& mesh,
    const PrimitiveMeshAnalysisOptions& options,
    const Vec3& model_lower,
    const Vec3& model_upper,
    const double model_diagonal,
    const double model_volume,
    const double model_surface_area,
    std::vector<bool>& responsibility_faces)
{
    struct CandidateBounds
    {
        Bounds bounds;
        std::size_t estimated_faces = 0;
        std::vector<std::uint32_t> component_faces;
        std::size_t component_count = 1;
    };
    const double tolerance = std::max(model_diagonal * 1.0e-9, 1.0e-10);
    const std::vector<bool> original_responsibility = responsibility_faces;
    const std::size_t original_face_count = static_cast<std::size_t>(std::count(
        original_responsibility.begin(), original_responsibility.end(), true));
    std::vector<std::unordered_set<std::uint32_t>> neighbors;
    const auto clusters = coplanarClusters(
        mesh, model_diagonal * options.coplanar_relative_tolerance,
        neighbors, &responsibility_faces);

    struct PlaneCluster
    {
        std::uint32_t id = 0;
        int normal_axis = 0;
        double coordinate = 0.0;
        Bounds bounds;
        double area = 0.0;
    };
    std::vector<PlaneCluster> planes;
    for (std::uint32_t id = 0; id < clusters.size(); ++id)
    {
        if (clusters[id].empty()) continue;
        const Face& face = mesh.faces[clusters[id].front()];
        const Vec3 normal = (mesh.vertices[face[1]] - mesh.vertices[face[0]])
                                .cross(mesh.vertices[face[2]] - mesh.vertices[face[0]])
                                .normalized();
        int axis = 0;
        for (int candidate = 1; candidate < 3; ++candidate)
            if (std::abs(normal[candidate]) > std::abs(normal[axis])) axis = candidate;
        const double area = facesArea(mesh, clusters[id]);
        if (std::abs(normal[axis]) < 1.0 - 1.0e-8 ||
            area < 0.002 * model_surface_area) continue;
        const auto [lower, upper] = faceBounds(mesh, clusters[id]);
        planes.push_back({id, axis,
                          0.5 * (lower[axis] + upper[axis]),
                          {lower, upper}, area});
    }

    const auto facesInside = [&](const Bounds& bounds,
                                 const std::vector<bool>& included_faces)
    {
        std::vector<std::uint32_t> faces;
        for (std::uint32_t face_id = 0; face_id < mesh.faces.size(); ++face_id)
        {
            if (!included_faces[face_id]) continue;
            bool inside = true;
            for (const auto vertex_id : mesh.faces[face_id])
                for (int axis = 0; axis < 3; ++axis)
                    inside &= mesh.vertices[vertex_id][axis] >= bounds.lower[axis] - tolerance &&
                              mesh.vertices[vertex_id][axis] <= bounds.upper[axis] + tolerance;
            if (inside) faces.push_back(face_id);
        }
        return faces;
    };
    const auto appendCandidate = [&](std::vector<CandidateBounds>& candidates,
                                     const Bounds& bounds,
                                     std::vector<std::uint32_t> component_faces = {},
                                     const std::size_t component_count = 1)
    {
        const Vec3 extent = bounds.upper - bounds.lower;
        if (extent.x() <= tolerance || extent.y() <= tolerance || extent.z() <= tolerance)
            return;
        const std::size_t face_count = component_faces.empty()
            ? facesInside(bounds, responsibility_faces).size()
            : component_faces.size();
        if (face_count < 12) return;
        candidates.push_back(
            {bounds, face_count, std::move(component_faces), component_count});
    };

    std::vector<CandidateBounds> candidates;
    struct SourceComponent
    {
        std::vector<std::uint32_t> faces;
        Bounds bounds;
        double bounds_volume = 0.0;
    };
    const std::vector<bool> all_source_faces(mesh.faces.size(), true);
    std::vector<SourceComponent> source_components;
    for (auto faces : faceComponents(mesh, all_source_faces))
    {
        const auto [component_lower, component_upper] = faceBounds(mesh, faces);
        const Vec3 extent = (component_upper - component_lower).cwiseMax(Vec3::Zero());
        source_components.push_back({std::move(faces),
                                     {component_lower, component_upper},
                                     extent.prod()});
    }
    for (std::size_t index = 0; index < source_components.size(); ++index)
        appendCandidate(candidates, source_components[index].bounds,
                        source_components[index].faces);

    const auto mergedBounds = [](const Bounds& first, const Bounds& second)
    {
        return Bounds{first.lower.cwiseMin(second.lower),
                      first.upper.cwiseMax(second.upper)};
    };
    const auto boundsVolume = [](const Bounds& bounds)
    {
        return (bounds.upper - bounds.lower).cwiseMax(Vec3::Zero()).prod();
    };
    const auto intersectionVolume = [&](const Bounds& first, const Bounds& second)
    {
        const Vec3 overlap = (first.upper.cwiseMin(second.upper) -
                              first.lower.cwiseMax(second.lower)).cwiseMax(Vec3::Zero());
        return overlap.prod();
    };
    const auto boundsBox = [](const Bounds& bounds)
    {
        BoxFit box;
        box.axes = Mat3::Identity();
        box.center = (bounds.lower + bounds.upper) * 0.5;
        box.half_size = (bounds.upper - bounds.lower) * 0.5;
        box.volume = 8.0 * box.half_size.prod();
        return box;
    };

    // Build a model-independent hierarchy over disconnected source components.
    // A merge is legal only when the final shared AABB passes the globally
    // normalized added-volume budget. Workload reduction affects ordering,
    // never the error certificate itself.
    struct ComponentGroup
    {
        std::vector<std::uint32_t> faces;
        Bounds bounds;
        double bounds_volume = 0.0;
        std::size_t component_count = 1;
        bool active = true;
    };
    std::vector<ComponentGroup> groups;
    groups.reserve(source_components.size() * 2);
    for (const auto& component : source_components)
        groups.push_back({component.faces, component.bounds,
                          component.bounds_volume});
    std::set<std::pair<std::size_t, std::size_t>> rejected_group_pairs;
    while (true)
    {
        std::size_t best_first = groups.size();
        std::size_t best_second = groups.size();
        double best_score = std::numeric_limits<double>::infinity();
        for (std::size_t first = 0; first < groups.size(); ++first)
        {
            if (!groups[first].active) continue;
            for (std::size_t second = first + 1; second < groups.size(); ++second)
            {
                if (!groups[second].active ||
                    rejected_group_pairs.contains({first, second})) continue;
                const Bounds merged = mergedBounds(groups[first].bounds,
                                                   groups[second].bounds);
                const double merged_volume = boundsVolume(merged);
                const double union_bounds_volume = groups[first].bounds_volume +
                    groups[second].bounds_volume -
                    intersectionVolume(groups[first].bounds, groups[second].bounds);
                const double lower_bound_added_volume_ratio = std::max(
                    merged_volume - union_bounds_volume, 0.0) /
                    std::max(model_volume, 1.0e-30);
                if (lower_bound_added_volume_ratio >
                    options.maximum_group_box_added_volume_ratio) continue;
                const std::size_t eliminated_shells =
                    groups[first].component_count + groups[second].component_count - 1;
                const double score = lower_bound_added_volume_ratio /
                    static_cast<double>(eliminated_shells);
                if (score < best_score)
                {
                    best_score = score;
                    best_first = first;
                    best_second = second;
                }
            }
        }
        if (best_first == groups.size()) break;

        const Bounds merged = mergedBounds(groups[best_first].bounds,
                                           groups[best_second].bounds);
        std::vector<std::uint32_t> merged_faces = groups[best_first].faces;
        merged_faces.insert(merged_faces.end(), groups[best_second].faces.begin(),
                            groups[best_second].faces.end());
        const BoxFit box = boundsBox(merged);
        const int resolution = std::clamp(
            static_cast<int>(options.projection_envelope_resolution / 8), 16, 32);
        const double added_volume_ratio = boxVisualHullAddedVolumeRatio(
            mesh, merged_faces, box, model_volume, resolution);
        if (added_volume_ratio > options.maximum_group_box_added_volume_ratio)
        {
            rejected_group_pairs.emplace(best_first, best_second);
            continue;
        }

        const std::size_t component_count = groups[best_first].component_count +
                                            groups[best_second].component_count;
        groups[best_first].active = false;
        groups[best_second].active = false;
        groups.push_back({merged_faces, merged, box.volume,
                          component_count, true});
        appendCandidate(candidates, merged, std::move(merged_faces), component_count);
    }

    const Vec3 model_extent = model_upper - model_lower;
    for (const auto& plane : planes)
    {
        for (const int sign : {-1, 1})
        {
            Bounds bounds = plane.bounds;
            if (sign < 0)
            {
                bounds.lower[plane.normal_axis] = model_lower[plane.normal_axis];
                bounds.upper[plane.normal_axis] = plane.coordinate;
            }
            else
            {
                bounds.lower[plane.normal_axis] = plane.coordinate;
                bounds.upper[plane.normal_axis] = model_upper[plane.normal_axis];
            }
            const double thickness = bounds.upper[plane.normal_axis] -
                                     bounds.lower[plane.normal_axis];
            if (thickness <= tolerance ||
                thickness > 0.10 * model_extent[plane.normal_axis]) continue;
            appendCandidate(candidates, bounds);
        }
    }

    std::sort(candidates.begin(), candidates.end(), [](const auto& first, const auto& second)
    {
        return first.estimated_faces > second.estimated_faces;
    });
    std::vector<RecognizedProtrusion> accepted;
    for (const auto& candidate : candidates)
    {
        const auto certification_faces = candidate.component_faces.empty()
            ? facesInside(candidate.bounds, original_responsibility)
            : candidate.component_faces;
        auto faces = candidate.component_faces.empty()
            ? facesInside(candidate.bounds, responsibility_faces)
            : candidate.component_faces;
        if (!candidate.component_faces.empty())
            std::erase_if(faces, [&](const auto face_id)
            { return !responsibility_faces[face_id]; });
        if (faces.size() < 12) continue;
        BoxFit box;
        box.axes = Mat3::Identity();
        box.center = (candidate.bounds.lower + candidate.bounds.upper) * 0.5;
        box.half_size = (candidate.bounds.upper - candidate.bounds.lower) * 0.5;
        box.volume = 8.0 * box.half_size.prod();
        const bool grouped_components = candidate.component_count > 1;
        if (!grouped_components && options.uniform_structure_policy &&
            options.forbid_main_body_box_approximation &&
            certification_faces.size() == original_face_count) continue;
        if (!grouped_components && box.volume / model_volume >
            4.0 * options.maximum_cavity_added_volume_ratio) continue;
        const int resolution = std::clamp(
            static_cast<int>(options.projection_envelope_resolution / 8), 16, 32);
        const double added_volume_ratio = boxVisualHullAddedVolumeRatio(
            mesh, certification_faces, box, model_volume, resolution);
        const double added_area_ratio = std::max(
            boxSurfaceArea(box) - facesArea(mesh, certification_faces), 0.0) /
            std::max(model_surface_area, 1.0e-30);
        const double volume_limit = grouped_components
            ? options.maximum_group_box_added_volume_ratio
            : options.maximum_cavity_added_volume_ratio;
        if (added_volume_ratio > volume_limit ||
            (!grouped_components &&
             added_area_ratio > options.maximum_box_envelope_added_area_ratio)) continue;
        for (const auto face_id : faces) responsibility_faces[face_id] = false;
        accepted.push_back({box, std::move(faces), {}, -1, 0.0});
    }
    return accepted;
}

std::vector<std::vector<std::size_t>> maskComponents(const std::vector<bool>& mask,
                                                     int width,
                                                     int height);

std::size_t projectionHoleCount(const std::vector<bool>& mask,
                                const int width,
                                const int height)
{
    std::vector<bool> filled = mask;
    fillProjectionHoles(filled, width, height);
    std::vector<bool> holes(mask.size(), false);
    for (std::size_t index = 0; index < mask.size(); ++index)
        holes[index] = filled[index] && !mask[index];
    return maskComponents(holes, width, height).size();
}

std::vector<bool> squareMorphology(const std::vector<bool>& input,
                                   const int width, const int height,
                                   const int radius, const bool dilate)
{
    if (radius <= 0) return input;
    std::vector<int> integral(static_cast<std::size_t>(width + 1) * (height + 1), 0);
    for (int y = 0; y < height; ++y)
        for (int x = 0; x < width; ++x)
            integral[static_cast<std::size_t>(y + 1) * (width + 1) + x + 1] =
                static_cast<int>(input[static_cast<std::size_t>(y) * width + x]) +
                integral[static_cast<std::size_t>(y) * (width + 1) + x + 1] +
                integral[static_cast<std::size_t>(y + 1) * (width + 1) + x] -
                integral[static_cast<std::size_t>(y) * (width + 1) + x];
    std::vector<bool> result(input.size(), false);
    for (int y = 0; y < height; ++y)
        for (int x = 0; x < width; ++x)
        {
            const int x0 = std::max(0, x - radius);
            const int y0 = std::max(0, y - radius);
            const int x1 = std::min(width - 1, x + radius);
            const int y1 = std::min(height - 1, y + radius);
            const int sum = integral[static_cast<std::size_t>(y1 + 1) * (width + 1) + x1 + 1] -
                            integral[static_cast<std::size_t>(y0) * (width + 1) + x1 + 1] -
                            integral[static_cast<std::size_t>(y1 + 1) * (width + 1) + x0] +
                            integral[static_cast<std::size_t>(y0) * (width + 1) + x0];
            const int area = (x1 - x0 + 1) * (y1 - y0 + 1);
            result[static_cast<std::size_t>(y) * width + x] = dilate ? sum > 0 : sum == area;
        }
    return result;
}

void closeShallowProjectionConcavities(std::vector<bool>& mask,
                                       const int width, const int height,
                                       const double ratio)
{
    const std::vector<bool> original = mask;
    const int radius = std::max(0, static_cast<int>(
        std::lround(ratio * std::min(width, height))));
    mask = squareMorphology(squareMorphology(mask, width, height, radius, true),
                            width, height, radius, false);
    for (std::size_t index = 0; index < mask.size(); ++index)
        mask[index] = mask[index] || original[index];
}

std::vector<std::vector<std::size_t>> maskComponents(const std::vector<bool>& mask,
                                                     const int width,
                                                     const int height)
{
    std::vector<bool> visited(mask.size(), false);
    std::vector<std::vector<std::size_t>> result;
    constexpr std::array<std::pair<int, int>, 4> offsets{{{1, 0}, {-1, 0}, {0, 1}, {0, -1}}};
    for (int y = 0; y < height; ++y)
        for (int x = 0; x < width; ++x)
        {
            const std::size_t seed = static_cast<std::size_t>(y) * width + x;
            if (!mask[seed] || visited[seed]) continue;
            std::vector<std::size_t> component;
            std::queue<std::pair<int, int>> queue;
            queue.emplace(x, y);
            visited[seed] = true;
            while (!queue.empty())
            {
                const auto [current_x, current_y] = queue.front();
                queue.pop();
                component.push_back(static_cast<std::size_t>(current_y) * width + current_x);
                for (const auto [dx, dy] : offsets)
                {
                    const int next_x = current_x + dx;
                    const int next_y = current_y + dy;
                    if (next_x < 0 || next_x >= width || next_y < 0 || next_y >= height)
                        continue;
                    const std::size_t next = static_cast<std::size_t>(next_y) * width + next_x;
                    if (!mask[next] || visited[next]) continue;
                    visited[next] = true;
                    queue.emplace(next_x, next_y);
                }
            }
            result.push_back(std::move(component));
        }
    return result;
}

void acceptProjectionCavities(std::vector<bool>& mask,
                              const int width,
                              const int height,
                              const double maximum_added_volume_ratio,
                              const double closing_ratio)
{
    const auto accept = [&](const std::vector<bool>& proposed,
                            const double minimum_added_volume_ratio = 0.0,
                            const bool require_interior_core = false)
    {
        std::vector<bool> added(mask.size(), false);
        for (std::size_t index = 0; index < mask.size(); ++index)
            added[index] = proposed[index] && !mask[index];
        for (const auto& component : maskComponents(added, width, height))
        {
            // Extruding one projection cell through the AABB contributes exactly
            // 1/(width*height) of its volume, independent of projection axis.
            const double volume_ratio = static_cast<double>(component.size()) /
                static_cast<double>(width * height);
            if (volume_ratio < minimum_added_volume_ratio ||
                volume_ratio > maximum_added_volume_ratio) continue;
            if (require_interior_core)
            {
                const std::unordered_set<std::size_t> members(
                    component.begin(), component.end());
                bool has_core = false;
                for (const auto index : component)
                {
                    const int x = static_cast<int>(index % width);
                    const int y = static_cast<int>(index / width);
                    if (x == 0 || x + 1 == width || y == 0 || y + 1 == height) continue;
                    bool full = true;
                    for (int dy = -1; dy <= 1; ++dy)
                        for (int dx = -1; dx <= 1; ++dx)
                            full &= members.contains(
                                static_cast<std::size_t>(y + dy) * width + x + dx);
                    if (full) { has_core = true; break; }
                }
                if (!has_core) continue;
            }
            for (const auto index : component) mask[index] = true;
        }
    };

    std::vector<bool> holes_filled = mask;
    fillProjectionHoles(holes_filled, width, height);
    accept(holes_filled);

    std::vector<bool> shallow_closed = mask;
    closeShallowProjectionConcavities(shallow_closed, width, height, closing_ratio);
    accept(shallow_closed);

    std::vector<Vec2> occupied;
    occupied.reserve(mask.size());
    for (int y = 0; y < height; ++y)
        for (int x = 0; x < width; ++x)
            if (mask[static_cast<std::size_t>(y) * width + x])
                occupied.emplace_back(x + 0.5, y + 0.5);
    if (occupied.size() >= 3)
    {
        const std::vector<Vec2> hull = convexHull(std::move(occupied), 1.0e-12);
        std::vector<bool> hull_fill(mask.size(), false);
        for (int y = 0; y < height; ++y)
            for (int x = 0; x < width; ++x)
                hull_fill[static_cast<std::size_t>(y) * width + x] =
                    pointInPolygon({x + 0.5, y + 0.5}, hull, 1.0e-12);
        // Each connected missing region is evaluated separately. Small shallow
        // recesses pass; a large main cavity remains because its prism volume
        // exceeds the same model-independent threshold.
        // Ignore one-cell raster cracks along diagonal source edges. Meaningful
        // shallow recesses remain well above this resolution-scaled threshold.
        accept(hull_fill, 0.00025, true);
    }

    // A convex hull cannot remove shallow exterior notches when another part of
    // the model keeps those notches on the hull boundary.  Such notches create
    // many long, one-cell-wide side rectangles after extrusion.  Propose the
    // projection bounding rectangle as one more fill source, but continue to
    // decide every connected missing region independently by added volume.  A
    // small rail recess can therefore be absorbed while a large body cavity is
    // rejected by exactly the same test.
    int minimum_x = width;
    int minimum_y = height;
    int maximum_x = -1;
    int maximum_y = -1;
    for (int y = 0; y < height; ++y)
        for (int x = 0; x < width; ++x)
            if (mask[static_cast<std::size_t>(y) * width + x])
            {
                minimum_x = std::min(minimum_x, x);
                minimum_y = std::min(minimum_y, y);
                maximum_x = std::max(maximum_x, x);
                maximum_y = std::max(maximum_y, y);
            }
    if (maximum_x >= minimum_x && maximum_y >= minimum_y)
    {
        std::vector<bool> rectangle_fill(mask.size(), false);
        for (int y = minimum_y; y <= maximum_y; ++y)
            for (int x = minimum_x; x <= maximum_x; ++x)
                rectangle_fill[static_cast<std::size_t>(y) * width + x] = true;
        accept(rectangle_fill, 0.00025, true);
    }
}

std::vector<OutputPrimitive> buildProjectionEnvelope(
    const Mesh& mesh,
    const PrimitiveMeshAnalysisOptions& options)
{
    const int resolution = static_cast<int>(options.projection_envelope_resolution);
    if (resolution < 16) return {};
    Vec3 lower = Vec3::Constant(std::numeric_limits<double>::infinity());
    Vec3 upper = Vec3::Constant(-std::numeric_limits<double>::infinity());
    for (const Vec3& vertex : mesh.vertices)
    {
        lower = lower.cwiseMin(vertex);
        upper = upper.cwiseMax(vertex);
    }
    Vec3 step;
    for (int axis = 0; axis < 3; ++axis)
    {
        const double extent = std::max(upper[axis] - lower[axis], 1.0e-9);
        step[axis] = extent / resolution;
    }
    auto xy = rasterizeProjection(mesh, 0, 1, lower, step, resolution, resolution);
    auto xz = rasterizeProjection(mesh, 0, 2, lower, step, resolution, resolution);
    auto yz = rasterizeProjection(mesh, 1, 2, lower, step, resolution, resolution);
    for (auto* mask : {&xy, &xz, &yz})
    {
        fillProjectionHoles(*mask, resolution, resolution);
        closeShallowProjectionConcavities(*mask, resolution, resolution,
                                          options.shallow_concavity_closing_ratio);
    }
    const auto occupied = [&](const int x, const int y, const int z)
    {
        return xy[static_cast<std::size_t>(y) * resolution + x] &&
               xz[static_cast<std::size_t>(z) * resolution + x] &&
               yz[static_cast<std::size_t>(z) * resolution + y];
    };
    std::vector<OutputPrimitive> output;
    const std::array<int, 3> tangent_first{{1, 0, 0}};
    const std::array<int, 3> tangent_second{{2, 2, 1}};
    for (int normal_axis = 0; normal_axis < 3; ++normal_axis)
        for (const int sign : {-1, 1})
            for (int slice = 0; slice < resolution; ++slice)
            {
                const int first_axis = tangent_first[normal_axis];
                const int second_axis = tangent_second[normal_axis];
                std::vector<bool> faces(static_cast<std::size_t>(resolution) * resolution, false);
                for (int second = 0; second < resolution; ++second)
                    for (int first = 0; first < resolution; ++first)
                    {
                        std::array<int, 3> cell{};
                        cell[normal_axis] = slice;
                        cell[first_axis] = first;
                        cell[second_axis] = second;
                        if (!occupied(cell[0], cell[1], cell[2])) continue;
                        cell[normal_axis] += sign;
                        const bool neighbor = cell[normal_axis] >= 0 &&
                            cell[normal_axis] < resolution &&
                            occupied(cell[0], cell[1], cell[2]);
                        if (!neighbor)
                            faces[static_cast<std::size_t>(second) * resolution + first] = true;
                    }
                for (int second = 0; second < resolution; ++second)
                    for (int first = 0; first < resolution; ++first)
                    {
                        const std::size_t seed = static_cast<std::size_t>(second) * resolution + first;
                        if (!faces[seed]) continue;
                        int width = 1;
                        while (first + width < resolution &&
                               faces[static_cast<std::size_t>(second) * resolution + first + width])
                            ++width;
                        int height = 1;
                        for (;;)
                        {
                            if (second + height >= resolution) break;
                            bool full = true;
                            for (int offset = 0; offset < width; ++offset)
                                full &= faces[static_cast<std::size_t>(second + height) * resolution +
                                              first + offset];
                            if (!full) break;
                            ++height;
                        }
                        for (int row = 0; row < height; ++row)
                            for (int column = 0; column < width; ++column)
                                faces[static_cast<std::size_t>(second + row) * resolution +
                                      first + column] = false;
                        Primitive rectangle;
                        rectangle.kind = Kind::Rectangle;
                        rectangle.axes = Mat3::Identity();
                        rectangle.axes.col(0) = Vec3::Zero();
                        rectangle.axes.col(1) = Vec3::Zero();
                        rectangle.axes.col(2) = Vec3::Zero();
                        rectangle.axes.values[first_axis][0] = 1.0;
                        rectangle.axes.values[second_axis][1] = 1.0;
                        rectangle.axes.values[normal_axis][2] = static_cast<double>(sign);
                        rectangle.center = lower;
                        rectangle.center[first_axis] += (first + width * 0.5) * step[first_axis];
                        rectangle.center[second_axis] += (second + height * 0.5) * step[second_axis];
                        rectangle.center[normal_axis] +=
                            (slice + (sign > 0 ? 1.0 : 0.0)) * step[normal_axis];
                        rectangle.half_size = {width * step[first_axis] * 0.5,
                                               height * step[second_axis] * 0.5, 0.0};
                        output.push_back({rectangle, {}});
                    }
            }
    return output;
}

Primitive makeAxisAlignedRectangle(const Vec3& center,
                                   const int first_axis,
                                   const int second_axis,
                                   const int normal_axis,
                                   const double first_length,
                                   const double second_length,
                                   const int normal_sign)
{
    Primitive rectangle;
    rectangle.kind = Kind::Rectangle;
    rectangle.axes = Mat3::Zero();
    rectangle.axes.values[first_axis][0] = 1.0;
    rectangle.axes.values[second_axis][1] = 1.0;
    rectangle.axes.values[normal_axis][2] = static_cast<double>(normal_sign);
    rectangle.center = center;
    rectangle.half_size = {first_length * 0.5, second_length * 0.5, 0.0};
    return rectangle;
}

std::vector<OutputPrimitive> buildAxisEnvelopeCandidate(
    const Mesh& mesh,
    const PrimitiveMeshAnalysisOptions& options,
    const Vec3& lower,
    const Vec3& upper,
    const int main_axis,
    double& added_volume_ratio,
    std::size_t& remaining_cavities)
{
    const int resolution = static_cast<int>(options.projection_envelope_resolution);
    const Vec3 extent = upper - lower;
    const int first_axis = (main_axis + 1) % 3;
    const int second_axis = (main_axis + 2) % 3;
    Vec3 step = Vec3::Zero();
    step[first_axis] = std::max(extent[first_axis], 1.0e-9) / resolution;
    step[second_axis] = std::max(extent[second_axis], 1.0e-9) / resolution;
    step[main_axis] = std::max(extent[main_axis], 1.0e-9);
    auto mask = rasterizeProjection(
        mesh, first_axis, second_axis, lower, step, resolution, resolution);
    acceptProjectionCavities(mask, resolution, resolution,
                             options.maximum_cavity_added_volume_ratio,
                             options.shallow_concavity_closing_ratio);
    remaining_cavities = projectionHoleCount(mask, resolution, resolution);

    Vec3 grid_step;
    for (int axis = 0; axis < 3; ++axis)
        grid_step[axis] = std::max(extent[axis], 1.0e-9) / resolution;
    const auto xy = rasterizeProjection(mesh, 0, 1, lower, grid_step,
                                        resolution, resolution);
    const auto xz = rasterizeProjection(mesh, 0, 2, lower, grid_step,
                                        resolution, resolution);
    const auto yz = rasterizeProjection(mesh, 1, 2, lower, grid_step,
                                        resolution, resolution);
    std::uint64_t added_cells = 0;
    for (int second = 0; second < resolution; ++second)
        for (int first = 0; first < resolution; ++first)
        {
            if (!mask[static_cast<std::size_t>(second) * resolution + first]) continue;
            for (int main = 0; main < resolution; ++main)
            {
                std::array<int, 3> cell{};
                cell[main_axis] = main;
                cell[first_axis] = first;
                cell[second_axis] = second;
                const bool visual_hull =
                    xy[static_cast<std::size_t>(cell[1]) * resolution + cell[0]] &&
                    xz[static_cast<std::size_t>(cell[2]) * resolution + cell[0]] &&
                    yz[static_cast<std::size_t>(cell[2]) * resolution + cell[1]];
                if (!visual_hull) ++added_cells;
            }
        }
    added_volume_ratio = static_cast<double>(added_cells) /
        (static_cast<double>(resolution) * resolution * resolution);

    std::vector<OutputPrimitive> output;
    // End caps: greedily tile the occupied cross-section with large rectangles.
    for (const int sign : {-1, 1})
    {
        std::vector<bool> remaining = mask;
        for (int second = 0; second < resolution; ++second)
            for (int first = 0; first < resolution; ++first)
            {
                const std::size_t seed = static_cast<std::size_t>(second) * resolution + first;
                if (!remaining[seed]) continue;
                int width = 1;
                while (first + width < resolution &&
                       remaining[static_cast<std::size_t>(second) * resolution + first + width])
                    ++width;
                int height = 1;
                while (second + height < resolution)
                {
                    bool full = true;
                    for (int offset = 0; offset < width; ++offset)
                        full &= remaining[static_cast<std::size_t>(second + height) * resolution +
                                          first + offset];
                    if (!full) break;
                    ++height;
                }
                for (int row = 0; row < height; ++row)
                    for (int column = 0; column < width; ++column)
                        remaining[static_cast<std::size_t>(second + row) * resolution +
                                  first + column] = false;
                Vec3 center = lower;
                center[first_axis] += (first + width * 0.5) * step[first_axis];
                center[second_axis] += (second + height * 0.5) * step[second_axis];
                center[main_axis] = sign < 0 ? lower[main_axis] : upper[main_axis];
                output.push_back({makeAxisAlignedRectangle(
                    center, first_axis, second_axis, main_axis,
                    width * step[first_axis], height * step[second_axis], sign), {}});
            }
    }

    const auto occupied = [&](const int first, const int second)
    {
        return first >= 0 && first < resolution && second >= 0 && second < resolution &&
               mask[static_cast<std::size_t>(second) * resolution + first];
    };
    // Side walls normal to the first cross-section axis. Merge consecutive
    // boundary cells along the second axis into one full-length rectangle.
    for (const int sign : {-1, 1})
        for (int boundary = 0; boundary <= resolution; ++boundary)
        {
            int second = 0;
            while (second < resolution)
            {
                const int cell_first = sign < 0 ? boundary : boundary - 1;
                const int neighbor_first = sign < 0 ? boundary - 1 : boundary;
                const bool exposed = occupied(cell_first, second) &&
                                     !occupied(neighbor_first, second);
                if (!exposed) { ++second; continue; }
                const int begin = second++;
                while (second < resolution && occupied(cell_first, second) &&
                       !occupied(neighbor_first, second)) ++second;
                Vec3 center = lower;
                center[main_axis] = (lower[main_axis] + upper[main_axis]) * 0.5;
                center[first_axis] += boundary * step[first_axis];
                center[second_axis] += (begin + (second - begin) * 0.5) * step[second_axis];
                output.push_back({makeAxisAlignedRectangle(
                    center, main_axis, second_axis, first_axis,
                    extent[main_axis], (second - begin) * step[second_axis], sign), {}});
            }
        }
    // Side walls normal to the second cross-section axis.
    for (const int sign : {-1, 1})
        for (int boundary = 0; boundary <= resolution; ++boundary)
        {
            int first = 0;
            while (first < resolution)
            {
                const int cell_second = sign < 0 ? boundary : boundary - 1;
                const int neighbor_second = sign < 0 ? boundary - 1 : boundary;
                const bool exposed = occupied(first, cell_second) &&
                                     !occupied(first, neighbor_second);
                if (!exposed) { ++first; continue; }
                const int begin = first++;
                while (first < resolution && occupied(first, cell_second) &&
                       !occupied(first, neighbor_second)) ++first;
                Vec3 center = lower;
                center[main_axis] = (lower[main_axis] + upper[main_axis]) * 0.5;
                center[first_axis] += (begin + (first - begin) * 0.5) * step[first_axis];
                center[second_axis] += boundary * step[second_axis];
                output.push_back({makeAxisAlignedRectangle(
                    center, main_axis, first_axis, second_axis,
                    extent[main_axis], (first - begin) * step[first_axis], sign), {}});
            }
        }
    return output;
}

struct StructuralCleanup
{
    struct CavityGroup
    {
        // Oriented sweep coordinates: column 0 separates the opposing walls;
        // columns 1 and 2 span their common projected responsibility.
        Mat3 frame = Mat3::Identity();
        Vec3 lower = Vec3::Zero();
        Vec3 upper = Vec3::Zero();
    };
    // Every face recognized as internal under the fully closed candidate.  The
    // main exterior analysis excludes this complete set so restored cavities
    // cannot perturb exterior region merging.
    std::vector<std::uint32_t> candidate_excluded_faces;
    // Over-budget connected voids are emitted later through an isolated cavity
    // surface pass.
    std::vector<std::uint32_t> restored_cavity_faces;
    std::vector<std::vector<std::uint32_t>> restored_cavity_face_groups;
    std::vector<CavityGroup> restored_cavities;
    std::vector<std::uint32_t> excluded_faces;
    std::vector<std::uint32_t> sealed_void_wall_faces;
    std::vector<std::uint32_t> blind_cavity_faces;
    std::vector<std::uint32_t> contained_faces;
    std::size_t sealed_void_wall_primitives = 0;
    std::size_t sealed_void_wall_triangles = 0;
    std::size_t blind_cavity_primitives = 0;
    std::size_t contained_primitives = 0;
};

StructuralCleanup identifyStructuralRedundantFaces(
    const Mesh& mesh,
    const PrimitiveMeshAnalysisOptions& options,
    const double model_diagonal,
    const double model_surface_area,
    const double model_volume)
{
    StructuralCleanup result;
    PrimitiveMeshAnalysisOptions closure_options = options;
    // This pass needs the closed candidate shell in order to identify its
    // interior walls.  Whether a cavity is actually allowed to close is decided
    // below per connected cavity against the caller's model-scale volume budget.
    closure_options.maximum_cavity_added_volume_ratio =
        std::numeric_limits<double>::infinity();
    closure_options.maximum_open_error_distance =
        options.maximum_open_error_distance >= 0.0
            ? options.maximum_open_error_distance
            : 100.0;
    std::vector<std::unordered_set<std::uint32_t>> ignored_neighbors;
    const auto clusters = coplanarClusters(
        mesh, model_diagonal * options.coplanar_relative_tolerance, ignored_neighbors);

    // This fixed, strength-independent surface partition is used only to identify
    // geometry that becomes redundant after an exterior hole/cavity is closed.
    // Performing the test on the final strength cut made the responsibility set
    // change with strength and could therefore make primitive counts decrease.
    std::vector<OutputPrimitive> structural_primitives;
    structural_primitives.reserve(clusters.size());
    std::size_t ignored_holes = 0;
    std::size_t ignored_boundary_voids = 0;
    double ignored_boundary_void_area = 0.0;
    for (const auto& cluster : clusters)
    {
        auto classified = classifyFinalRegion(
            mesh, cluster, closure_options, std::numeric_limits<double>::infinity(),
            model_surface_area,
            ignored_holes, ignored_boundary_voids, ignored_boundary_void_area);
        structural_primitives.insert(
            structural_primitives.end(), std::make_move_iterator(classified.begin()),
            std::make_move_iterator(classified.end()));
    }

    const double tolerance = std::max(model_diagonal * 1.0e-9, 1.0e-10);
    Vec3 lower = Vec3::Constant(std::numeric_limits<double>::infinity());
    Vec3 upper = Vec3::Constant(-std::numeric_limits<double>::infinity());
    for (const Vec3& vertex : mesh.vertices)
    {
        lower = lower.cwiseMin(vertex);
        upper = upper.cwiseMax(vertex);
    }
    const Vec3 model_center = (lower + upper) * 0.5;
    auto append = [&](std::vector<std::uint32_t> faces)
    {
        result.excluded_faces.insert(result.excluded_faces.end(), faces.begin(), faces.end());
    };

    std::size_t sealed_void_wall_primitives = 0;
    auto sealed = removeSealedVoidWalls(
        structural_primitives, model_diagonal, tolerance,
        sealed_void_wall_primitives);
    result.sealed_void_wall_primitives = sealed_void_wall_primitives;
    result.sealed_void_wall_triangles = sealed.size();
    result.sealed_void_wall_faces = sealed;
    append(std::move(sealed));
    auto blind = removeBlindCavitySurfaces(
        structural_primitives, model_center, model_diagonal, tolerance,
        result.blind_cavity_primitives);
    result.blind_cavity_faces = blind;
    append(std::move(blind));

    // Boundary-void closure modifies the cap silhouette before ordinary region
    // analysis. Any bevel wall that is then strictly inside the resulting
    // certified closed volume belongs to that filled void and must lose its
    // source responsibility here, not in a late visual cleanup pass. Run this
    // after the broad sealed/blind cleanup so the certification graph remains
    // small even for triangle soups with thousands of cavity-wall patches.
    std::vector<RecognizedProtrusion> structural_boxes;
    std::size_t ignored_structural_box_count = 0;
    recognizeClosedAxisAlignedBoxes(
        structural_primitives, tolerance, ignored_structural_box_count,
        &structural_boxes);
    std::vector<BoxFit> structural_box_volumes;
    for (const auto& box : structural_boxes)
        structural_box_volumes.push_back(box.box);
    promoteToSemanticPrimitives(structural_primitives);
    const auto structural_extrusions = recognizeCertifiedPrismaticVolumes(
        structural_primitives, tolerance);
    std::size_t boundary_void_wall_primitives = 0;
    std::vector<std::uint32_t> boundary_void_wall_faces;
    structural_primitives = clipPlanarOcclusionByClosedVolumes(
        std::move(structural_primitives), structural_box_volumes,
        structural_extrusions, tolerance, boundary_void_wall_primitives,
        &boundary_void_wall_faces);
    result.sealed_void_wall_primitives += boundary_void_wall_primitives;
    result.sealed_void_wall_triangles += boundary_void_wall_faces.size();
    result.sealed_void_wall_faces.insert(
        result.sealed_void_wall_faces.end(),
        boundary_void_wall_faces.begin(), boundary_void_wall_faces.end());
    append(std::move(boundary_void_wall_faces));

    auto contained = removeContainedPrimitives(
        structural_primitives, tolerance, result.contained_primitives);
    result.contained_faces = contained;
    append(std::move(contained));
    std::sort(result.excluded_faces.begin(), result.excluded_faces.end());
    result.excluded_faces.erase(
        std::unique(result.excluded_faces.begin(), result.excluded_faces.end()),
        result.excluded_faces.end());
    result.candidate_excluded_faces = result.excluded_faces;

    // Build geometric void groups. CAD triangle soups often duplicate vertices
    // across the two opposing walls of one cavity, so index connectivity alone
    // splits a single void. Pair thin, parallel components only when their
    // tangent projections nearly coincide and cover a macroscopic fraction of
    // the model cross-section. Small through-holes remain surface holes and are
    // intentionally filled by planar classification.
    const double cavity_volume_limit =
        options.maximum_cavity_added_volume_ratio * model_volume;
    struct Component
    {
        std::vector<std::uint32_t> faces;
        Vec3 lower;
        Vec3 upper;
        Vec3 extent;
        Vec3 normal = Vec3::Zero();
        Mat3 frame = Mat3::Identity();
        Bounds2 tangent_bounds;
        double normal_lower = 0.0;
        double normal_upper = 0.0;
        bool planar = false;
    };
    std::vector<Component> components;
    for (auto& faces : faceComponentsFromList(mesh, result.excluded_faces))
    {
        const auto [component_lower, component_upper] = faceBounds(mesh, faces);
        Component component;
        component.faces = std::move(faces);
        component.lower = component_lower;
        component.upper = component_upper;
        component.extent = (component_upper - component_lower).cwiseMax(Vec3::Zero());
        const auto component_vertices = uniqueVertices(mesh, component.faces);
        if (component_vertices.size() >= 3)
        {
            const BoxFit oriented_bounds = fitBox(mesh, component_vertices);
            int thin_axis = 0;
            for (int axis = 1; axis < 3; ++axis)
                if (oriented_bounds.half_size[axis] <
                    oriented_bounds.half_size[thin_axis]) thin_axis = axis;
            component.normal = oriented_bounds.axes.col(thin_axis);
            // Canonicalize an unoriented plane normal lexicographically. Using
            // the numerically largest component makes a 45-degree rotation
            // unstable: opposite normals can choose different tied axes and
            // consequently different tangent frames.
            for (int axis = 0; axis < 3; ++axis)
            {
                if (std::abs(component.normal[axis]) <= 1.0e-12) continue;
                if (component.normal[axis] < 0.0)
                    component.normal = -component.normal;
                break;
            }
            const Mat3 basis = orthonormalFrame(component.normal);
            component.frame.col(0) = component.normal;
            component.frame.col(1) = static_cast<Vec3>(basis.col(1));
            component.frame.col(2) = static_cast<Vec3>(basis.col(2));
            component.normal_lower = std::numeric_limits<double>::infinity();
            component.normal_upper = -std::numeric_limits<double>::infinity();
            for (const auto vertex_id : component_vertices)
            {
                const Vec3& vertex = mesh.vertices[vertex_id];
                const double normal_coordinate = vertex.dot(component.normal);
                const Vec2 tangent_coordinate(
                    vertex.dot(component.frame.col(1)),
                    vertex.dot(component.frame.col(2)));
                component.normal_lower = std::min(
                    component.normal_lower, normal_coordinate);
                component.normal_upper = std::max(
                    component.normal_upper, normal_coordinate);
                component.tangent_bounds.lower =
                    component.tangent_bounds.lower.cwiseMin(tangent_coordinate);
                component.tangent_bounds.upper =
                    component.tangent_bounds.upper.cwiseMax(tangent_coordinate);
            }
            const Vec2 tangent_extent =
                component.tangent_bounds.upper - component.tangent_bounds.lower;
            const double tangent_scale = std::max(
                {tangent_extent.x(), tangent_extent.y(), tolerance});
            component.planar = component.normal_upper - component.normal_lower <=
                0.1 * tangent_scale;
        }
        components.push_back(std::move(component));
    }
    std::vector<std::size_t> parent(components.size());
    std::iota(parent.begin(), parent.end(), 0);
    const auto find = [&](std::size_t index)
    {
        std::size_t root = index;
        while (parent[root] != root) root = parent[root];
        while (parent[index] != index)
        {
            const std::size_t next = parent[index];
            parent[index] = root;
            index = next;
        }
        return root;
    };
    for (std::size_t first = 0; first < components.size(); ++first)
        for (std::size_t second = first + 1; second < components.size(); ++second)
        {
            const Component& a = components[first];
            const Component& b = components[second];
            if (!a.planar || !b.planar ||
                std::abs(a.normal.dot(b.normal)) < 1.0 - 1.0e-8) continue;
            // Canonicalized parallel normals produce the same deterministic
            // tangent frame, so their 2D bounds are directly comparable.
            const Vec2 a_tangent_extent =
                (a.tangent_bounds.upper - a.tangent_bounds.lower).cwiseMax(Vec2::Zero());
            const Vec2 b_tangent_extent =
                (b.tangent_bounds.upper - b.tangent_bounds.lower).cwiseMax(Vec2::Zero());
            const double tangent_scale = std::max(
                {a_tangent_extent.x(), a_tangent_extent.y(),
                 b_tangent_extent.x(), b_tangent_extent.y(), tolerance});
            const double a_thickness = a.normal_upper - a.normal_lower;
            const double b_thickness = b.normal_upper - b.normal_lower;
            if (a_thickness > 0.1 * tangent_scale ||
                b_thickness > 0.1 * tangent_scale) continue;
            const Vec2 overlap =
                a.tangent_bounds.upper.cwiseMin(b.tangent_bounds.upper) -
                a.tangent_bounds.lower.cwiseMax(b.tangent_bounds.lower);
            if (overlap.x() < 0.9 * std::min(a_tangent_extent.x(), b_tangent_extent.x()) ||
                overlap.y() < 0.9 * std::min(a_tangent_extent.y(), b_tangent_extent.y()))
                continue;
            const double projection_area = overlap.prod();
            Bounds2 model_projection;
            for (const Vec3& vertex : mesh.vertices)
            {
                const Vec2 projected(vertex.dot(a.frame.col(1)),
                                     vertex.dot(a.frame.col(2)));
                model_projection.lower = model_projection.lower.cwiseMin(projected);
                model_projection.upper = model_projection.upper.cwiseMax(projected);
            }
            const double model_projection_area =
                (model_projection.upper - model_projection.lower).prod();
            if (projection_area < 0.25 * model_projection_area) continue;
            const double gap = std::max({b.normal_lower - a.normal_upper,
                                         a.normal_lower - b.normal_upper, 0.0});
            if (gap <= 2.0 * std::max(a_thickness, b_thickness)) continue;
            parent[find(second)] = find(first);
        }
    std::unordered_map<std::size_t, std::vector<std::size_t>> groups;
    for (std::size_t index = 0; index < components.size(); ++index)
        groups[find(index)].push_back(index);
    std::vector<std::uint32_t> certified_excluded;
    for (const auto& [root, members] : groups)
    {
        (void)root;
        Vec3 group_lower = Vec3::Constant(std::numeric_limits<double>::infinity());
        Vec3 group_upper = Vec3::Constant(-std::numeric_limits<double>::infinity());
        for (const auto member : members)
        {
            group_lower = group_lower.cwiseMin(components[member].lower);
            group_upper = group_upper.cwiseMax(components[member].upper);
        }
        double swept_volume = (group_upper - group_lower).prod();
        if (components[members.front()].planar)
        {
            swept_volume = 0.0;
            const Component& reference = components[members.front()];
            Bounds2 tangent_overlap;
            tangent_overlap.lower = Vec2::Constant(-std::numeric_limits<double>::infinity());
            tangent_overlap.upper = Vec2::Constant(std::numeric_limits<double>::infinity());
            double normal_lower = std::numeric_limits<double>::infinity();
            double normal_upper = -std::numeric_limits<double>::infinity();
            bool common_plane_family = true;
            for (const auto member : members)
            {
                const Component& component = components[member];
                common_plane_family &= component.planar &&
                    std::abs(component.normal.dot(reference.normal)) >= 1.0 - 1.0e-8;
                tangent_overlap.lower = tangent_overlap.lower.cwiseMax(
                    component.tangent_bounds.lower);
                tangent_overlap.upper = tangent_overlap.upper.cwiseMin(
                    component.tangent_bounds.upper);
                normal_lower = std::min(normal_lower, component.normal_lower);
                normal_upper = std::max(normal_upper, component.normal_upper);
            }
            const Vec2 tangent_extent = tangent_overlap.upper - tangent_overlap.lower;
            if (common_plane_family && members.size() > 1 &&
                tangent_extent.x() > tolerance && tangent_extent.y() > tolerance)
            {
                swept_volume = (normal_upper - normal_lower) *
                    tangent_extent.prod();
                if (swept_volume > cavity_volume_limit)
                {
                    StructuralCleanup::CavityGroup cavity;
                    cavity.frame = reference.frame;
                    cavity.lower = Vec3(normal_lower, tangent_overlap.lower.x(),
                                        tangent_overlap.lower.y());
                    cavity.upper = Vec3(normal_upper, tangent_overlap.upper.x(),
                                        tangent_overlap.upper.y());
                    result.restored_cavities.push_back(cavity);
                }
            }
        }
        if (swept_volume > cavity_volume_limit)
        {
            std::vector<std::uint32_t> restored_group;
            for (const auto member : members)
            {
                result.restored_cavity_faces.insert(
                    result.restored_cavity_faces.end(),
                    components[member].faces.begin(), components[member].faces.end());
                restored_group.insert(
                    restored_group.end(), components[member].faces.begin(),
                    components[member].faces.end());
            }
            std::sort(restored_group.begin(), restored_group.end());
            restored_group.erase(
                std::unique(restored_group.begin(), restored_group.end()),
                restored_group.end());
            result.restored_cavity_face_groups.push_back(
                std::move(restored_group));
        }
        else for (const auto member : members)
            certified_excluded.insert(certified_excluded.end(),
                components[member].faces.begin(), components[member].faces.end());
    }
    result.excluded_faces = std::move(certified_excluded);
    if (result.excluded_faces.empty())
    {
        result.sealed_void_wall_primitives = 0;
        result.sealed_void_wall_triangles = 0;
        result.blind_cavity_primitives = 0;
        result.contained_primitives = 0;
    }
    return result;
}

std::vector<OutputPrimitive> reopenOverBudgetPlanarHoles(
    const Mesh& mesh,
    std::vector<OutputPrimitive> primitives,
    const double maximum_cavity_added_volume_ratio,
    const double model_volume,
    const double model_diagonal,
    const double tolerance)
{
    if (maximum_cavity_added_volume_ratio < 0.0 ||
        std::isinf(maximum_cavity_added_volume_ratio)) return primitives;

    constexpr int clipper_precision = 8;
    const double volume_limit = maximum_cavity_added_volume_ratio * model_volume;
    std::vector<OutputPrimitive> result;
    result.reserve(primitives.size());

    for (auto& item : primitives)
    {
        Vec3 normal;
        if (!planarNormal(item.primitive, tolerance, normal))
        {
            result.push_back(std::move(item));
            continue;
        }
        // No subset of this proxy can exceed the per-hole budget. This cheap
        // upper bound avoids planar Boolean work on thousands of small patches.
        if (primitiveSurfaceArea(item.primitive, tolerance) * model_diagonal <=
            volume_limit)
        {
            result.push_back(std::move(item));
            continue;
        }

        const Vec3 origin = planarPoint(item.primitive);
        const Mat3 basis = orthonormalFrame(normal);
        Mat3 frame;
        frame.col(0) = static_cast<Vec3>(basis.col(1));
        frame.col(1) = static_cast<Vec3>(basis.col(2));
        frame.col(2) = normal;
        const auto project = [&](const Vec3& point)
        {
            const Vec3 local = frame.transposeMultiply(point - origin);
            return Clipper2Lib::PointD(local.x(), local.y());
        };

        Clipper2Lib::PathsD subjects;
        const PrimitiveMesh primitive_mesh = triangulatePrimitive(item.primitive);
        for (const Face& face : primitive_mesh.faces)
        {
            Clipper2Lib::PathD triangle;
            for (const auto vertex : face)
                triangle.push_back(project(primitive_mesh.vertices[vertex]));
            if (std::abs(Clipper2Lib::Area(triangle)) <= tolerance * tolerance) continue;
            if (Clipper2Lib::Area(triangle) < 0.0)
                std::reverse(triangle.begin(), triangle.end());
            subjects.push_back(std::move(triangle));
        }
        Clipper2Lib::PathsD source_paths;
        for (const auto face_id : item.source_faces)
        {
            const Face& face = mesh.faces[face_id];
            Clipper2Lib::PathD triangle;
            bool coplanar = true;
            for (const auto vertex_id : face)
            {
                const Vec3& vertex = mesh.vertices[vertex_id];
                coplanar &= std::abs((vertex - origin).dot(normal)) <= tolerance * 8.0;
                triangle.push_back(project(vertex));
            }
            if (!coplanar ||
                std::abs(Clipper2Lib::Area(triangle)) <= tolerance * tolerance) continue;
            if (Clipper2Lib::Area(triangle) < 0.0)
                std::reverse(triangle.begin(), triangle.end());
            source_paths.push_back(std::move(triangle));
        }
        if (subjects.empty() || source_paths.empty())
        {
            result.push_back(std::move(item));
            continue;
        }
        const auto subject_union = Clipper2Lib::Union(
            subjects, Clipper2Lib::FillRule::NonZero, clipper_precision);
        const auto source_union = Clipper2Lib::Union(
            source_paths, Clipper2Lib::FillRule::NonZero, clipper_precision);
        const auto added_regions = Clipper2Lib::Difference(
            subject_union, source_union, Clipper2Lib::FillRule::NonZero,
            clipper_precision);

        double model_lower = std::numeric_limits<double>::infinity();
        double model_upper = -std::numeric_limits<double>::infinity();
        for (const Vec3& vertex : mesh.vertices)
        {
            const double coordinate = vertex.dot(normal);
            model_lower = std::min(model_lower, coordinate);
            model_upper = std::max(model_upper, coordinate);
        }
        const double conservative_depth = model_upper - model_lower;
        const double boundary_tolerance_squared =
            tolerance * tolerance * 64.0 * 64.0;
        const auto squaredPointSegmentDistance = [](const Clipper2Lib::PointD& point,
                                                    const Clipper2Lib::PointD& first,
                                                    const Clipper2Lib::PointD& second)
        {
            const double dx = second.x - first.x;
            const double dy = second.y - first.y;
            const double length_squared = dx * dx + dy * dy;
            const double parameter = length_squared > 0.0
                ? std::clamp(((point.x - first.x) * dx + (point.y - first.y) * dy) /
                             length_squared, 0.0, 1.0)
                : 0.0;
            const double offset_x = point.x - (first.x + parameter * dx);
            const double offset_y = point.y - (first.y + parameter * dy);
            return offset_x * offset_x + offset_y * offset_y;
        };
        Clipper2Lib::PathsD clips;
        for (const auto& region : added_regions)
        {
            if (!Clipper2Lib::IsPositive(region) ||
                std::abs(Clipper2Lib::Area(region)) * conservative_depth <= volume_limit)
                continue;
            bool touches_outer_boundary = false;
            for (const auto& point : region)
            {
                for (const auto& subject : subject_union)
                {
                    if (!Clipper2Lib::IsPositive(subject)) continue;
                    for (std::size_t edge = 0; edge < subject.size(); ++edge)
                    {
                        if (squaredPointSegmentDistance(
                                point, subject[edge], subject[(edge + 1) % subject.size()]) <=
                            boundary_tolerance_squared)
                        {
                            touches_outer_boundary = true;
                            break;
                        }
                    }
                    if (touches_outer_boundary) break;
                }
                if (touches_outer_boundary) break;
            }
            // Added regions touching the proxy's exterior boundary are ordinary
            // silhouette/concavity approximations, not enclosed planar holes.
            if (!touches_outer_boundary) clips.push_back(region);
        }
        if (clips.empty())
        {
            result.push_back(std::move(item));
            continue;
        }
        const auto clip_union = Clipper2Lib::Union(
            clips, Clipper2Lib::FillRule::NonZero, clipper_precision);
        auto remaining = Clipper2Lib::Difference(
            subject_union, clip_union, Clipper2Lib::FillRule::NonZero,
            clipper_precision);
        remaining.insert(remaining.end(), source_union.begin(), source_union.end());
        remaining = Clipper2Lib::Union(
            remaining, Clipper2Lib::FillRule::NonZero, clipper_precision);

        std::vector<OutputPrimitive> replacements;
        bool conversion_failed = false;
        const bool has_holes = std::any_of(remaining.begin(), remaining.end(),
            [](const auto& path) { return !Clipper2Lib::IsPositive(path); });
        if (has_holes)
        {
            Clipper2Lib::PathsD triangles;
            if (Clipper2Lib::Triangulate(
                    remaining, clipper_precision, triangles, false) !=
                Clipper2Lib::TriangulateResult::success)
            {
                conversion_failed = true;
            }
            else for (const auto& triangle : triangles)
            {
                std::vector<Vec2> boundary;
                for (const auto& point : triangle)
                    boundary.emplace_back(point.x, point.y);
                if (boundary.size() != 3 ||
                    std::abs(signedArea(boundary)) <= tolerance * tolerance)
                    continue;
                OutputPrimitive replacement;
                replacement.primitive = polygonPrimitive(boundary, origin, frame);
                replacement.preserves_cavity_opening = true;
                replacements.push_back(std::move(replacement));
            }
        }
        else for (const auto& path : remaining)
        {
            std::vector<Vec2> boundary;
            boundary.reserve(path.size());
            for (const auto& point : path) boundary.emplace_back(point.x, point.y);
            boundary = simplifyPolygon(std::move(boundary), tolerance);
            if (boundary.size() < 3 || !simplePolygon(boundary, tolerance) ||
                triangulatePolygon(boundary).size() + 2 != boundary.size())
            {
                conversion_failed = true;
                break;
            }
            OutputPrimitive replacement;
            replacement.primitive = polygonPrimitive(boundary, origin, frame);
            replacement.preserves_cavity_opening = true;
            replacements.push_back(std::move(replacement));
        }
        if (conversion_failed)
        {
            result.push_back(std::move(item));
            continue;
        }

        // The replacement is (proxy minus over-budget holes) union the complete
        // coplanar source projection. This is a planar coverage certificate and
        // discards the source triangulation's internal edges. Non-planar cavity
        // walls remain the responsibility of the isolated cavity pass.
        if (!replacements.empty())
            replacements.front().source_faces = item.source_faces;
        result.insert(result.end(), std::make_move_iterator(replacements.begin()),
                      std::make_move_iterator(replacements.end()));
    }
    return result;
}

std::vector<OutputPrimitive> reopenRestoredCavityVolumes(
    const Mesh& mesh,
    std::vector<OutputPrimitive> primitives,
    const StructuralCleanup& cleanup,
    const double maximum_cavity_added_volume_ratio,
    const double model_volume,
    const double tolerance)
{
    if (cleanup.restored_cavities.empty()) return primitives;
    constexpr int clipper_precision = 8;
    constexpr std::array<std::array<int, 2>, 12> box_edges{{
        {{0, 1}}, {{0, 2}}, {{0, 4}}, {{1, 3}}, {{1, 5}}, {{2, 3}},
        {{2, 6}}, {{3, 7}}, {{4, 5}}, {{4, 6}}, {{5, 7}}, {{6, 7}}
    }};
    std::vector<OutputPrimitive> result;
    result.reserve(primitives.size());
    const double cavity_volume_limit =
        maximum_cavity_added_volume_ratio * model_volume;

    for (auto& item : primitives)
    {
        Vec3 normal;
        if (!planarNormal(item.primitive, tolerance, normal))
        {
            result.push_back(std::move(item));
            continue;
        }
        const Vec3 origin = planarPoint(item.primitive);

        // The two end walls of an over-budget cavity are part of the opening,
        // even when they survived the structural pass and therefore were not
        // emitted through restored_cavity_output. Mark them here so later
        // parallel/volume occlusion cannot mistake an exposed cavity wall for
        // a redundant inner skin hidden by the corresponding exterior wall.
        bool bounds_restored_cavity = false;
        for (const auto& cavity : cleanup.restored_cavities)
        {
            const Vec3 separation = cavity.frame.col(0);
            if (std::abs(normal.dot(separation)) < 1.0 - 1.0e-8)
                continue;
            const double coordinate = origin.dot(separation);
            if (std::abs(coordinate - cavity.lower.x()) > tolerance * 16.0 &&
                std::abs(coordinate - cavity.upper.x()) > tolerance * 16.0)
                continue;

            Bounds2 projected;
            const PrimitiveMesh planar_mesh =
                triangulatePrimitive(item.primitive);
            for (const Vec3& vertex : planar_mesh.vertices)
            {
                projected.lower = projected.lower.cwiseMin({
                    vertex.dot(cavity.frame.col(1)),
                    vertex.dot(cavity.frame.col(2))});
                projected.upper = projected.upper.cwiseMax({
                    vertex.dot(cavity.frame.col(1)),
                    vertex.dot(cavity.frame.col(2))});
            }
            const Vec2 overlap =
                projected.upper.cwiseMin(
                    Vec2(cavity.upper.y(), cavity.upper.z())) -
                projected.lower.cwiseMax(
                    Vec2(cavity.lower.y(), cavity.lower.z()));
            if (overlap.x() > tolerance && overlap.y() > tolerance)
            {
                bounds_restored_cavity = true;
                break;
            }
        }
        if (bounds_restored_cavity)
        {
            item.preserves_cavity_opening = true;
            result.push_back(std::move(item));
            continue;
        }

        const Mat3 basis = orthonormalFrame(normal);
        Mat3 frame;
        frame.col(0) = static_cast<Vec3>(basis.col(1));
        frame.col(1) = static_cast<Vec3>(basis.col(2));
        frame.col(2) = normal;
        const auto project = [&](const Vec3& point)
        {
            const Vec3 local = frame.transposeMultiply(point - origin);
            return Clipper2Lib::PointD(local.x(), local.y());
        };

        struct SectionCandidate
        {
            Clipper2Lib::PathD path;
            double depth = 0.0;
        };
        std::vector<SectionCandidate> sections;
        for (const auto& cavity : cleanup.restored_cavities)
        {
            const Vec3 separation = cavity.frame.col(0);
            if (std::abs(normal.dot(separation)) >= 1.0 - 1.0e-8)
            {
                const double coordinate = origin.dot(separation);
                if (std::abs(coordinate - cavity.lower.x()) <= tolerance * 16.0 ||
                    std::abs(coordinate - cavity.upper.x()) <= tolerance * 16.0)
                    continue;
            }

            std::array<Vec3, 8> corners;
            for (int index = 0; index < 8; ++index)
            {
                const Vec3 local(
                    (index & 1) ? cavity.upper.x() : cavity.lower.x(),
                    (index & 2) ? cavity.upper.y() : cavity.lower.y(),
                    (index & 4) ? cavity.upper.z() : cavity.lower.z());
                corners[index] = cavity.frame * local;
            }
            std::vector<Vec3> section_points;
            const auto appendUnique = [&](const Vec3& point)
            {
                if (std::none_of(section_points.begin(), section_points.end(),
                    [&](const Vec3& existing)
                    { return (existing - point).norm() <= tolerance * 8.0; }))
                    section_points.push_back(point);
            };
            std::array<double, 8> distances{};
            for (int index = 0; index < 8; ++index)
            {
                distances[index] = (corners[index] - origin).dot(normal);
                if (std::abs(distances[index]) <= tolerance * 8.0)
                    appendUnique(corners[index]);
            }
            for (const auto& edge : box_edges)
            {
                const double first_distance = distances[edge[0]];
                const double second_distance = distances[edge[1]];
                if ((first_distance < -tolerance && second_distance < -tolerance) ||
                    (first_distance > tolerance && second_distance > tolerance) ||
                    std::abs(first_distance - second_distance) <= tolerance)
                    continue;
                const double parameter = std::clamp(
                    first_distance / (first_distance - second_distance), 0.0, 1.0);
                appendUnique(corners[edge[0]] +
                    (corners[edge[1]] - corners[edge[0]]) * parameter);
            }
            if (section_points.size() < 3) continue;
            std::vector<Vec2> section;
            section.reserve(section_points.size());
            for (const Vec3& point : section_points)
            {
                const auto projected = project(point);
                section.emplace_back(projected.x, projected.y);
            }
            section = simplifyPolygon(
                convexHull(std::move(section), tolerance), tolerance);
            if (section.size() < 3 ||
                std::abs(signedArea(section)) <= tolerance * tolerance)
                continue;
            Clipper2Lib::PathD clip;
            for (const Vec2& point : section) clip.emplace_back(point.x(), point.y());
            if (Clipper2Lib::Area(clip) < 0.0)
                std::reverse(clip.begin(), clip.end());
            const Vec3 cavity_extent = cavity.upper - cavity.lower;
            const double depth =
                std::abs(normal.dot(cavity.frame.col(0))) * cavity_extent.x() +
                std::abs(normal.dot(cavity.frame.col(1))) * cavity_extent.y() +
                std::abs(normal.dot(cavity.frame.col(2))) * cavity_extent.z();
            sections.push_back({std::move(clip), depth});
        }
        if (sections.empty())
        {
            result.push_back(std::move(item));
            continue;
        }

        Clipper2Lib::PathsD subjects;
        const PrimitiveMesh item_mesh = triangulatePrimitive(item.primitive);
        for (const Face& face : item_mesh.faces)
        {
            Clipper2Lib::PathD triangle;
            for (const auto vertex : face)
                triangle.push_back(project(item_mesh.vertices[vertex]));
            if (std::abs(Clipper2Lib::Area(triangle)) <= tolerance * tolerance) continue;
            if (Clipper2Lib::Area(triangle) < 0.0)
                std::reverse(triangle.begin(), triangle.end());
            subjects.push_back(std::move(triangle));
        }
        const auto subject_union = Clipper2Lib::Union(
            subjects, Clipper2Lib::FillRule::NonZero, clipper_precision);
        const double subject_area = std::abs(Clipper2Lib::Area(subject_union));
        Clipper2Lib::PathsD source_paths;
        for (std::uint32_t face_id = 0; face_id < mesh.faces.size(); ++face_id)
        {
            const Face& face = mesh.faces[face_id];
            Clipper2Lib::PathD triangle;
            bool coplanar = true;
            for (const auto vertex_id : face)
            {
                const Vec3& vertex = mesh.vertices[vertex_id];
                coplanar &= std::abs((vertex - origin).dot(normal)) <= tolerance * 8.0;
                triangle.push_back(project(vertex));
            }
            if (!coplanar ||
                std::abs(Clipper2Lib::Area(triangle)) <= tolerance * tolerance) continue;
            if (Clipper2Lib::Area(triangle) < 0.0)
                std::reverse(triangle.begin(), triangle.end());
            source_paths.push_back(std::move(triangle));
        }
        Clipper2Lib::PathsD source_union;
        if (!source_paths.empty())
            source_union = Clipper2Lib::Union(
                source_paths, Clipper2Lib::FillRule::NonZero, clipper_precision);

        // Estimate the local shell thickness from the nearest parallel source
        // surface with substantial projected overlap.  A restored cavity can
        // be hundreds of times deeper than the plate containing a hole; using
        // the cavity depth would incorrectly reopen inexpensive through-holes.
        double local_shell_depth = std::numeric_limits<double>::infinity();
        std::map<double, Clipper2Lib::PathsD> parallel_planes;
        for (const Face& face : mesh.faces)
        {
            const Vec3& first = mesh.vertices[face[0]];
            const Vec3& second = mesh.vertices[face[1]];
            const Vec3& third = mesh.vertices[face[2]];
            Vec3 face_normal = (second - first).cross(third - first);
            const double face_normal_length = face_normal.norm();
            if (face_normal_length <= tolerance * tolerance) continue;
            face_normal = face_normal / face_normal_length;
            if (std::abs(face_normal.dot(normal)) < 1.0 - 1.0e-8) continue;
            const double depth = std::abs((first - origin).dot(normal));
            if (depth <= tolerance * 16.0) continue;
            Clipper2Lib::PathD triangle;
            triangle.reserve(3);
            for (const auto vertex_id : face)
                triangle.push_back(project(mesh.vertices[vertex_id]));
            if (std::abs(Clipper2Lib::Area(triangle)) <= tolerance * tolerance)
                continue;
            if (Clipper2Lib::Area(triangle) < 0.0)
                std::reverse(triangle.begin(), triangle.end());
            auto plane = parallel_planes.lower_bound(depth - tolerance * 16.0);
            if (plane == parallel_planes.end() ||
                std::abs(plane->first - depth) > tolerance * 16.0)
                plane = parallel_planes.emplace(depth, Clipper2Lib::PathsD{}).first;
            plane->second.push_back(std::move(triangle));
        }
        for (const auto& [depth, paths] : parallel_planes)
        {
            const auto plane_union = Clipper2Lib::Union(
                paths, Clipper2Lib::FillRule::NonZero, clipper_precision);
            const auto overlap = Clipper2Lib::Intersect(
                subject_union, plane_union, Clipper2Lib::FillRule::NonZero,
                clipper_precision);
            if (std::abs(Clipper2Lib::Area(overlap)) >= 0.25 * subject_area)
            {
                local_shell_depth = depth;
                break;
            }
        }

        Clipper2Lib::PathsD clips;
        for (const auto& section : sections)
        {
            const auto candidate = Clipper2Lib::Intersect(
                subject_union, Clipper2Lib::PathsD{section.path},
                Clipper2Lib::FillRule::NonZero, clipper_precision);
            const auto unsupported = source_union.empty()
                ? candidate
                : Clipper2Lib::Difference(
                    candidate, source_union, Clipper2Lib::FillRule::NonZero,
                    clipper_precision);
            for (const auto& path : unsupported)
            {
                if (!Clipper2Lib::IsPositive(path) || path.size() < 3) continue;
                const double fill_depth = std::min(section.depth, local_shell_depth);
                if (std::abs(Clipper2Lib::Area(path)) * fill_depth <=
                    cavity_volume_limit) continue;
                clips.push_back(path);
            }
        }
        if (clips.empty())
        {
            result.push_back(std::move(item));
            continue;
        }
        const auto clip_union = Clipper2Lib::Union(
            clips, Clipper2Lib::FillRule::NonZero, clipper_precision);
        const auto removed = Clipper2Lib::Intersect(
            subject_union, clip_union, Clipper2Lib::FillRule::NonZero,
                clipper_precision);
        if (std::abs(Clipper2Lib::Area(removed)) <= tolerance * tolerance)
        {
            result.push_back(std::move(item));
            continue;
        }
        const auto remaining = Clipper2Lib::Difference(
            subject_union, clip_union, Clipper2Lib::FillRule::NonZero,
            clipper_precision);

        std::vector<OutputPrimitive> replacements;
        bool conversion_failed = false;
        const bool has_holes = std::any_of(remaining.begin(), remaining.end(),
            [](const auto& path) { return !Clipper2Lib::IsPositive(path); });
        if (has_holes)
        {
            Clipper2Lib::PathsD triangles;
            if (Clipper2Lib::Triangulate(
                    remaining, clipper_precision, triangles, false) !=
                Clipper2Lib::TriangulateResult::success)
                conversion_failed = true;
            else for (const auto& triangle : triangles)
            {
                std::vector<Vec2> boundary;
                for (const auto& point : triangle)
                    boundary.emplace_back(point.x, point.y);
                if (boundary.size() != 3 ||
                    std::abs(signedArea(boundary)) <= tolerance * tolerance) continue;
                OutputPrimitive replacement;
                replacement.primitive = polygonPrimitive(boundary, origin, frame);
                replacement.preserves_cavity_opening = true;
                replacements.push_back(std::move(replacement));
            }
        }
        else for (const auto& path : remaining)
        {
            std::vector<Vec2> boundary;
            for (const auto& point : path) boundary.emplace_back(point.x, point.y);
            boundary = simplifyPolygon(std::move(boundary), tolerance);
            if (boundary.size() < 3 || !simplePolygon(boundary, tolerance) ||
                triangulatePolygon(boundary).size() + 2 != boundary.size())
            {
                conversion_failed = true;
                break;
            }
            OutputPrimitive replacement;
            replacement.primitive = polygonPrimitive(boundary, origin, frame);
            replacement.preserves_cavity_opening = true;
            replacements.push_back(std::move(replacement));
        }
        if (conversion_failed)
        {
            result.push_back(std::move(item));
            continue;
        }
        if (!replacements.empty()) replacements.front().source_faces = item.source_faces;
        result.insert(result.end(), std::make_move_iterator(replacements.begin()),
                      std::make_move_iterator(replacements.end()));
    }
    return result;
}

const char* kindName(const Kind kind)
{
    switch (kind)
    {
    case Kind::Rectangle: return "rectangle";
    case Kind::Triangle: return "triangle";
    case Kind::Polygon: return "polygon";
    case Kind::Disk: return "disk";
    case Kind::Annulus: return "annulus";
    case Kind::CylindricalBand: return "cylindricalband";
    case Kind::ConicalBand: return "conicalband";
    }
    return "unknown";
}

void promoteToSemanticPrimitives(std::vector<OutputPrimitive>& primitives)
{
    for (OutputPrimitive& item : primitives)
    {
        Primitive& primitive = item.primitive;
        if (primitive.kind == Kind::Triangle)
        {
            primitive.polygon.assign(primitive.triangle.begin(), primitive.triangle.end());
            primitive.kind = Kind::Polygon;
            continue;
        }
        if (primitive.kind != Kind::Rectangle) continue;

        const Vec3 first = primitive.axes.col(0) * primitive.half_size.x();
        const Vec3 second = primitive.axes.col(1) * primitive.half_size.y();
        primitive.polygon = {
            primitive.center - first - second,
            primitive.center + first - second,
            primitive.center + first + second,
            primitive.center - first + second,
        };
        primitive.kind = Kind::Polygon;
    }
}

std::vector<OutputPrimitive> coalesceShallowParallelShells(
    std::vector<OutputPrimitive> primitives,
    const Vec3& model_center,
    const double model_diagonal,
    const double model_volume,
    const PrimitiveMeshAnalysisOptions& options,
    const double tolerance,
    const bool broad_terrace_phase)
{
    if (primitives.size() < 2) return primitives;
    std::vector<std::size_t> parent(primitives.size());
    std::iota(parent.begin(), parent.end(), 0);
    const auto find = [&](std::size_t index)
    {
        std::size_t root = index;
        while (parent[root] != root) root = parent[root];
        while (parent[index] != index)
        {
            const std::size_t next = parent[index];
            parent[index] = root;
            index = next;
        }
        return root;
    };
    struct ShellEdge
    {
        std::size_t first = 0;
        std::size_t second = 0;
        double swept_volume_ratio = 0.0;
        double exterior_priority = 0.0;
    };
    std::vector<ShellEdge> shell_edges;
    const auto projectedBoundsOf = [](const Primitive& polygon,
                                      const Vec3& origin,
                                      const Mat3& frame)
    {
        Bounds2 bounds;
        for (const Vec3& vertex : polygon.polygon)
        {
            const Vec3 local = frame.transposeMultiply(vertex - origin);
            bounds.lower = bounds.lower.cwiseMin({local.y(), local.z()});
            bounds.upper = bounds.upper.cwiseMax({local.y(), local.z()});
        }
        return bounds;
    };
    const auto boundsOverlap = [tolerance](const Bounds2& first,
                                            const Bounds2& second)
    {
        return first.upper.x() >= second.lower.x() - tolerance &&
               first.lower.x() <= second.upper.x() + tolerance &&
               first.upper.y() >= second.lower.y() - tolerance &&
               first.lower.y() <= second.upper.y() + tolerance;
    };
    for (std::size_t first = 0; first < primitives.size(); ++first)
    {
        if (primitives[first].enclosure_group != 0 ||
            primitives[first].preserves_cavity_opening ||
            primitives[first].primitive.kind != Kind::Polygon) continue;
        Vec3 first_normal;
        if (!planarNormal(primitives[first].primitive, tolerance, first_normal)) continue;
        int dominant_axis = 0;
        for (int axis = 1; axis < 3; ++axis)
            if (std::abs(first_normal[axis]) > std::abs(first_normal[dominant_axis]))
                dominant_axis = axis;
        if (first_normal[dominant_axis] < 0.0) first_normal *= -1.0;
        const Vec3 origin = primitives[first].primitive.polygon.front();
        const Mat3 frame = orthonormalFrame(first_normal);
        const Bounds2 first_bounds = projectedBoundsOf(
            primitives[first].primitive, origin, frame);
        const double first_distance = first_normal.dot(origin);
        for (std::size_t second = first + 1; second < primitives.size(); ++second)
        {
            if (primitives[second].enclosure_group != 0 ||
                primitives[second].preserves_cavity_opening ||
                primitives[second].primitive.kind != Kind::Polygon) continue;
            Vec3 second_normal;
            if (!planarNormal(primitives[second].primitive, tolerance, second_normal) ||
                std::abs(second_normal.dot(first_normal)) < 1.0 - 1.0e-8) continue;
            const double second_distance = first_normal.dot(
                primitives[second].primitive.polygon.front());
            const double depth = std::abs(second_distance - first_distance);
            if (depth <= tolerance) continue;
            const double shallow_depth =
                options.shallow_parallel_merge_depth_relative * model_diagonal;
            if ((!broad_terrace_phase && depth > shallow_depth) ||
                (broad_terrace_phase && depth <= shallow_depth))
                continue;
            const Bounds2 second_bounds = projectedBoundsOf(
                primitives[second].primitive, origin, frame);
            const double overlap_first = std::max(
                0.0, std::min(first_bounds.upper.x(), second_bounds.upper.x()) -
                     std::max(first_bounds.lower.x(), second_bounds.lower.x()));
            const double overlap_second = std::max(
                0.0, std::min(first_bounds.upper.y(), second_bounds.upper.y()) -
                     std::max(first_bounds.lower.y(), second_bounds.lower.y()));
            if (overlap_first <= tolerance || overlap_second <= tolerance) continue;

            std::size_t bridge_count = 0;
            const double lower_distance = std::min(first_distance, second_distance);
            const double upper_distance = std::max(first_distance, second_distance);
            for (std::size_t bridge = 0; bridge < primitives.size(); ++bridge)
            {
                if (bridge == first || bridge == second ||
                    primitives[bridge].primitive.kind != Kind::Polygon) continue;
                Vec3 bridge_normal;
                if (!planarNormal(primitives[bridge].primitive, tolerance, bridge_normal) ||
                    std::abs(bridge_normal.dot(first_normal)) > 1.0e-6) continue;
                double minimum = std::numeric_limits<double>::infinity();
                double maximum = -std::numeric_limits<double>::infinity();
                for (const Vec3& vertex : primitives[bridge].primitive.polygon)
                {
                    const double distance = first_normal.dot(vertex);
                    minimum = std::min(minimum, distance);
                    maximum = std::max(maximum, distance);
                }
                if (minimum > lower_distance + tolerance ||
                    maximum < upper_distance - tolerance) continue;
                const Bounds2 bridge_bounds = projectedBoundsOf(
                    primitives[bridge].primitive, origin, frame);
                if (boundsOverlap(bridge_bounds, first_bounds) &&
                    boundsOverlap(bridge_bounds, second_bounds))
                    ++bridge_count;
            }
            if (bridge_count < 1) continue;
            const double swept_volume_ratio =
                std::max(primitiveSurfaceArea(primitives[first].primitive, tolerance),
                         primitiveSurfaceArea(primitives[second].primitive, tolerance)) *
                depth / std::max(model_volume, 1.0e-30);
            if (swept_volume_ratio > options.maximum_cavity_added_volume_ratio) continue;
            Bounds2 pair_bounds;
            pair_bounds.lower = first_bounds.lower.cwiseMin(second_bounds.lower);
            pair_bounds.upper = first_bounds.upper.cwiseMax(second_bounds.upper);
            double support_minimum = std::min(first_distance, second_distance);
            double support_maximum = std::max(first_distance, second_distance);
            for (const OutputPrimitive& support : primitives)
            {
                if (support.primitive.kind != Kind::Polygon) continue;
                Vec3 support_normal;
                if (!planarNormal(support.primitive, tolerance, support_normal) ||
                    std::abs(support_normal.dot(first_normal)) < 1.0 - 1.0e-8) continue;
                const Bounds2 support_bounds = projectedBoundsOf(
                    support.primitive, origin, frame);
                const double overlap_x =
                    std::min(pair_bounds.upper.x(), support_bounds.upper.x()) -
                    std::max(pair_bounds.lower.x(), support_bounds.lower.x());
                const double overlap_y =
                    std::min(pair_bounds.upper.y(), support_bounds.upper.y()) -
                    std::max(pair_bounds.lower.y(), support_bounds.lower.y());
                if (overlap_x <= tolerance || overlap_y <= tolerance) continue;
                const double distance = first_normal.dot(
                    support.primitive.polygon.front());
                support_minimum = std::min(support_minimum, distance);
                support_maximum = std::max(support_maximum, distance);
            }
            const double support_span = support_maximum - support_minimum;
            const double local_center = support_span > 2.0 * depth + tolerance
                ? 0.5 * (support_minimum + support_maximum)
                : first_normal.dot(model_center);
            if ((first_distance - local_center) *
                    (second_distance - local_center) <= 0.0 &&
                depth > options.shallow_parallel_merge_depth_relative * model_diagonal)
                continue;
            const double exterior_priority = std::max(
                std::abs(first_distance - local_center),
                std::abs(second_distance - local_center)) /
                std::max(support_span, tolerance);
            shell_edges.push_back(
                {first, second, swept_volume_ratio, exterior_priority});
        }
    }

    std::sort(shell_edges.begin(), shell_edges.end(), [](const auto& first, const auto& second)
    {
        if (first.exterior_priority != second.exterior_priority)
            return first.exterior_priority > second.exterior_priority;
        if (first.swept_volume_ratio != second.swept_volume_ratio)
            return first.swept_volume_ratio < second.swept_volume_ratio;
        if (first.first != second.first) return first.first < second.first;
        return first.second < second.second;
    });
    std::vector<std::vector<std::size_t>> members(primitives.size());
    for (std::size_t index = 0; index < primitives.size(); ++index)
        members[index].push_back(index);
    const auto groupSweepRatioToTarget = [&primitives, tolerance, model_volume](
        const std::vector<std::size_t>& group,
        const double requested_target)
    {
        Vec3 normal;
        if (group.empty() ||
            !planarNormal(primitives[group.front()].primitive, tolerance, normal))
            return std::numeric_limits<double>::infinity();
        int dominant_axis = 0;
        for (int axis = 1; axis < 3; ++axis)
            if (std::abs(normal[axis]) > std::abs(normal[dominant_axis]))
                dominant_axis = axis;
        if (normal[dominant_axis] < 0.0) normal *= -1.0;
        const Vec3 origin = primitives[group.front()].primitive.polygon.front();
        const Mat3 frame = orthonormalFrame(normal);
        struct Layer
        {
            double distance = 0.0;
            Clipper2Lib::PathD path;
        };
        std::vector<Layer> layers;
        std::vector<double> levels;
        for (const auto index : group)
        {
            const double distance = normal.dot(
                primitives[index].primitive.polygon.front());
            Clipper2Lib::PathD path;
            for (const Vec3& vertex : primitives[index].primitive.polygon)
            {
                const Vec3 local = frame.transposeMultiply(vertex - origin);
                path.emplace_back(local.y(), local.z());
            }
            if (Clipper2Lib::Area(path) < 0.0) std::reverse(path.begin(), path.end());
            layers.push_back({distance, std::move(path)});
            levels.push_back(distance);
        }
        std::sort(levels.begin(), levels.end());
        levels.erase(std::unique(levels.begin(), levels.end(), [tolerance](
            const double first, const double second)
            { return std::abs(first - second) <= tolerance; }), levels.end());
        if (levels.size() < 2) return 0.0;
        const bool target_is_minimum =
            std::abs(requested_target - levels.front()) <=
            std::abs(requested_target - levels.back());
        constexpr int clipper_precision = 8;
        double volume = 0.0;
        for (std::size_t level = 0; level + 1 < levels.size(); ++level)
        {
            const double lower = levels[level];
            const double upper = levels[level + 1];
            Clipper2Lib::PathsD active;
            for (const Layer& layer : layers)
            {
                const bool crosses = target_is_minimum
                    ? layer.distance >= upper - tolerance
                    : layer.distance <= lower + tolerance;
                if (crosses) active.push_back(layer.path);
            }
            if (active.empty()) continue;
            const double area = std::abs(Clipper2Lib::Area(Clipper2Lib::Union(
                active, Clipper2Lib::FillRule::NonZero, clipper_precision)));
            volume += area * (upper - lower);
        }
        return volume / std::max(model_volume, 1.0e-30);
    };
    const auto groupEnvelopeRatio = [&](const std::vector<std::size_t>& group)
    {
        Vec3 normal;
        if (group.empty() ||
            !planarNormal(primitives[group.front()].primitive, tolerance, normal))
            return std::numeric_limits<double>::infinity();
        int dominant_axis = 0;
        for (int axis = 1; axis < 3; ++axis)
            if (std::abs(normal[axis]) > std::abs(normal[dominant_axis]))
                dominant_axis = axis;
        if (normal[dominant_axis] < 0.0) normal *= -1.0;
        const Vec3 origin = primitives[group.front()].primitive.polygon.front();
        const Mat3 frame = orthonormalFrame(normal);
        double minimum = std::numeric_limits<double>::infinity();
        double maximum = -std::numeric_limits<double>::infinity();
        struct ProjectedLayer
        {
            double distance = 0.0;
            Clipper2Lib::PathD path;
        };
        std::vector<ProjectedLayer> projected;
        for (const auto index : group)
        {
            const double distance = normal.dot(
                primitives[index].primitive.polygon.front());
            minimum = std::min(minimum, distance);
            maximum = std::max(maximum, distance);
            Clipper2Lib::PathD path;
            for (const Vec3& vertex : primitives[index].primitive.polygon)
            {
                const Vec3 local = frame.transposeMultiply(vertex - origin);
                path.emplace_back(local.y(), local.z());
            }
            if (Clipper2Lib::Area(path) < 0.0) std::reverse(path.begin(), path.end());
            projected.push_back({distance, std::move(path)});
        }
        if (maximum - minimum >
            options.shallow_parallel_merge_depth_relative * model_diagonal)
        {
            Clipper2Lib::PathsD all_paths;
            Clipper2Lib::PathsD minimum_paths;
            Clipper2Lib::PathsD maximum_paths;
            for (const auto& layer : projected)
            {
                all_paths.push_back(layer.path);
                if (std::abs(layer.distance - minimum) <= tolerance)
                    minimum_paths.push_back(layer.path);
                if (std::abs(layer.distance - maximum) <= tolerance)
                    maximum_paths.push_back(layer.path);
            }
            constexpr int clipper_precision = 8;
            const auto all_union = Clipper2Lib::Union(
                all_paths, Clipper2Lib::FillRule::NonZero, clipper_precision);
            const auto minimum_union = Clipper2Lib::Union(
                minimum_paths, Clipper2Lib::FillRule::NonZero, clipper_precision);
            const auto maximum_union = Clipper2Lib::Union(
                maximum_paths, Clipper2Lib::FillRule::NonZero, clipper_precision);
            const double all_area = std::abs(Clipper2Lib::Area(all_union));
            const double minimum_area = std::abs(Clipper2Lib::Area(minimum_union));
            const double maximum_area = std::abs(Clipper2Lib::Area(maximum_union));
            const double intersection_area = std::abs(Clipper2Lib::Area(
                Clipper2Lib::Intersect(minimum_union, maximum_union,
                    Clipper2Lib::FillRule::NonZero, clipper_precision)));
            if (minimum_area / std::max(all_area, 1.0e-30) >= 0.85 &&
                maximum_area / std::max(all_area, 1.0e-30) >= 0.85 &&
                intersection_area /
                    std::max(std::max(minimum_area, maximum_area), 1.0e-30) >= 0.85)
                return std::numeric_limits<double>::infinity();
        }
        return std::min(groupSweepRatioToTarget(group, minimum),
                        groupSweepRatioToTarget(group, maximum));
    };
    for (const ShellEdge& edge : shell_edges)
    {
        const std::size_t first_root = find(edge.first);
        const std::size_t second_root = find(edge.second);
        if (first_root == second_root) continue;
        std::vector<std::size_t> proposed = members[first_root];
        proposed.insert(proposed.end(), members[second_root].begin(),
                        members[second_root].end());
        const double proposed_ratio = groupEnvelopeRatio(proposed);
        if (proposed_ratio > options.maximum_cavity_added_volume_ratio)
            continue;
        parent[second_root] = first_root;
        members[first_root] = std::move(proposed);
        members[second_root].clear();
    }

    std::unordered_map<std::size_t, std::vector<std::size_t>> groups;
    for (std::size_t index = 0; index < primitives.size(); ++index)
        groups[find(index)].push_back(index);
    for (const auto& [root, group] : groups)
    {
        (void)root;
        if (group.size() < 2) continue;
        Vec3 normal;
        if (!planarNormal(primitives[group.front()].primitive, tolerance, normal)) continue;
        int dominant_axis = 0;
        for (int axis = 1; axis < 3; ++axis)
            if (std::abs(normal[axis]) > std::abs(normal[dominant_axis])) dominant_axis = axis;
        if (normal[dominant_axis] < 0.0) normal *= -1.0;
        const Vec3 origin = primitives[group.front()].primitive.polygon.front();
        const Mat3 frame = orthonormalFrame(normal);
        Bounds2 group_bounds;
        double group_minimum = std::numeric_limits<double>::infinity();
        double group_maximum = -std::numeric_limits<double>::infinity();
        for (const auto index : group)
        {
            const double distance = normal.dot(
                primitives[index].primitive.polygon.front());
            group_minimum = std::min(group_minimum, distance);
            group_maximum = std::max(group_maximum, distance);
            const Bounds2 bounds = projectedBoundsOf(
                primitives[index].primitive, origin, frame);
            group_bounds.lower = group_bounds.lower.cwiseMin(bounds.lower);
            group_bounds.upper = group_bounds.upper.cwiseMax(bounds.upper);
        }
        double support_minimum = group_minimum;
        double support_maximum = group_maximum;
        for (const OutputPrimitive& candidate : primitives)
        {
            if (candidate.primitive.kind != Kind::Polygon) continue;
            Vec3 candidate_normal;
            if (!planarNormal(candidate.primitive, tolerance, candidate_normal) ||
                std::abs(candidate_normal.dot(normal)) < 1.0 - 1.0e-8) continue;
            const Bounds2 bounds = projectedBoundsOf(candidate.primitive, origin, frame);
            const double overlap_x = std::min(group_bounds.upper.x(), bounds.upper.x()) -
                                     std::max(group_bounds.lower.x(), bounds.lower.x());
            const double overlap_y = std::min(group_bounds.upper.y(), bounds.upper.y()) -
                                     std::max(group_bounds.lower.y(), bounds.lower.y());
            if (overlap_x <= tolerance || overlap_y <= tolerance) continue;
            const double distance = normal.dot(candidate.primitive.polygon.front());
            support_minimum = std::min(support_minimum, distance);
            support_maximum = std::max(support_maximum, distance);
        }
        const double group_depth = group_maximum - group_minimum;
        constexpr int clipper_precision = 8;
        Clipper2Lib::PathsD all_paths;
        Clipper2Lib::PathsD minimum_paths;
        Clipper2Lib::PathsD maximum_paths;
        for (const auto index : group)
        {
            const double distance = normal.dot(
                primitives[index].primitive.polygon.front());
            Clipper2Lib::PathD path;
            for (const Vec3& vertex : primitives[index].primitive.polygon)
            {
                const Vec3 local = frame.transposeMultiply(vertex - origin);
                path.emplace_back(local.y(), local.z());
            }
            if (Clipper2Lib::Area(path) < 0.0) std::reverse(path.begin(), path.end());
            all_paths.push_back(path);
            if (std::abs(distance - group_minimum) <= tolerance)
                minimum_paths.push_back(path);
            if (std::abs(distance - group_maximum) <= tolerance)
                maximum_paths.push_back(std::move(path));
        }
        const Clipper2Lib::PathsD all_union = Clipper2Lib::Union(
            all_paths, Clipper2Lib::FillRule::NonZero, clipper_precision);
        const Clipper2Lib::PathsD minimum_union = Clipper2Lib::Union(
            minimum_paths, Clipper2Lib::FillRule::NonZero, clipper_precision);
        const Clipper2Lib::PathsD maximum_union = Clipper2Lib::Union(
            maximum_paths, Clipper2Lib::FillRule::NonZero, clipper_precision);
        const double all_area = std::abs(Clipper2Lib::Area(all_union));
        const double minimum_area = std::abs(Clipper2Lib::Area(minimum_union));
        const double maximum_area = std::abs(Clipper2Lib::Area(maximum_union));
        const double extreme_intersection_area = std::abs(Clipper2Lib::Area(
            Clipper2Lib::Intersect(minimum_union, maximum_union,
                Clipper2Lib::FillRule::NonZero, clipper_precision)));
        const double extreme_similarity = extreme_intersection_area /
            std::max(std::max(minimum_area, maximum_area), 1.0e-30);
        const bool opposing_full_faces =
            minimum_area / std::max(all_area, 1.0e-30) >= 0.85 &&
            maximum_area / std::max(all_area, 1.0e-30) >= 0.85 &&
            extreme_similarity >= 0.85;
        if (opposing_full_faces && group_depth >
            options.shallow_parallel_merge_depth_relative * model_diagonal)
            continue;
        const bool has_local_depth_context =
            support_maximum - support_minimum > 2.0 * group_depth + tolerance;
        const double center_distance = has_local_depth_context
            ? 0.5 * (support_minimum + support_maximum)
            : normal.dot(model_center);
        double target_distance = normal.dot(
            primitives[group.front()].primitive.polygon.front());
        for (const auto index : group)
        {
            const double distance = normal.dot(
                primitives[index].primitive.polygon.front());
            const double target_exposure = std::abs(target_distance - center_distance);
            const double candidate_exposure = std::abs(distance - center_distance);
            if (candidate_exposure > target_exposure + tolerance ||
                (std::abs(candidate_exposure - target_exposure) <= tolerance &&
                 distance > target_distance))
                target_distance = distance;
        }
        const double group_swept_volume_ratio =
            groupSweepRatioToTarget(group, target_distance);
        if (group_swept_volume_ratio > options.maximum_cavity_added_volume_ratio)
            continue;
        for (const auto index : group)
        {
            for (Vec3& vertex : primitives[index].primitive.polygon)
                vertex += normal * (target_distance - normal.dot(vertex));
            primitives[index].shallow_shell_coalesced = true;
        }
    }
    return primitives;
}

void simplifyConservativePolygonDetails(
    std::vector<OutputPrimitive>& primitives,
    const double model_surface_area,
    const double maximum_added_area_ratio,
    const double tolerance)
{
    constexpr int clipper_precision = 8;
    const auto canonicalizeLoop = [](std::vector<Vec2>& polygon)
    {
        if (polygon.empty()) return;
        Clipper2Lib::PathD path;
        path.reserve(polygon.size());
        for (const Vec2& point : polygon) path.emplace_back(point.x(), point.y());
        if (Clipper2Lib::Area(path) < 0.0) std::reverse(polygon.begin(), polygon.end());
        const auto first = std::min_element(
            polygon.begin(), polygon.end(), [](const Vec2& left, const Vec2& right)
            {
                if (left.x() != right.x()) return left.x() < right.x();
                return left.y() < right.y();
            });
        std::rotate(polygon.begin(), first, polygon.end());
    };
    for (OutputPrimitive& item : primitives)
    {
        Primitive& primitive = item.primitive;
        if (item.preserves_cavity_opening || !item.shallow_shell_coalesced ||
            primitive.kind != Kind::Polygon ||
            primitive.polygon.size() <= 3) continue;
        Vec3 normal;
        if (!planarNormal(primitive, tolerance, normal)) continue;
        const Vec3 origin = primitive.polygon.front();
        const Mat3 frame = orthonormalFrame(normal);
        std::vector<Vec2> polygon;
        polygon.reserve(primitive.polygon.size());
        for (const Vec3& vertex : primitive.polygon)
        {
            const Vec3 local = frame.transposeMultiply(vertex - origin);
            polygon.emplace_back(local.y(), local.z());
        }
        canonicalizeLoop(polygon);
        double spent_area = 0.0;
        for (;;)
        {
            Clipper2Lib::PathD original_path;
            for (const Vec2& point : polygon)
                original_path.emplace_back(point.x(), point.y());
            if (Clipper2Lib::Area(original_path) < 0.0)
                std::reverse(original_path.begin(), original_path.end());
            const double original_area = std::abs(Clipper2Lib::Area(original_path));
            std::vector<Vec2> best;
            double best_cost = std::numeric_limits<double>::infinity();
            std::size_t best_removed = 0;
            for (std::size_t first = 0; first < polygon.size(); ++first)
                for (std::size_t second = first + 2; second < polygon.size(); ++second)
                {
                    const std::size_t removed = second - first - 1;
                    if (polygon.size() - removed < 3) continue;
                    std::vector<Vec2> candidate;
                    candidate.reserve(polygon.size() - removed);
                    candidate.insert(candidate.end(), polygon.begin(),
                                     polygon.begin() + static_cast<std::ptrdiff_t>(first + 1));
                    candidate.insert(candidate.end(),
                                     polygon.begin() + static_cast<std::ptrdiff_t>(second),
                                     polygon.end());
                    if (!simplePolygon(candidate, tolerance) ||
                        triangulatePolygon(candidate).size() + 2 != candidate.size()) continue;
                    Clipper2Lib::PathD candidate_path;
                    for (const Vec2& point : candidate)
                        candidate_path.emplace_back(point.x(), point.y());
                    if (Clipper2Lib::Area(candidate_path) < 0.0)
                        std::reverse(candidate_path.begin(), candidate_path.end());
                    const double missing_area = std::abs(Clipper2Lib::Area(
                        Clipper2Lib::Difference(
                            Clipper2Lib::PathsD{original_path},
                            Clipper2Lib::PathsD{candidate_path},
                            Clipper2Lib::FillRule::NonZero, clipper_precision)));
                    if (missing_area > tolerance * tolerance) continue;
                    const double added_area = std::max(
                        std::abs(Clipper2Lib::Area(candidate_path)) - original_area, 0.0);
                    if ((spent_area + added_area) /
                            std::max(model_surface_area, 1.0e-30) >
                        maximum_added_area_ratio) continue;
                    const double cost = added_area / static_cast<double>(removed);
                    if (cost < best_cost ||
                        (cost == best_cost && removed > best_removed))
                    {
                        best = std::move(candidate);
                        best_cost = cost;
                        best_removed = removed;
                    }
                }
            // A short bevel between two meaningful boundary lines should be
            // simplified by extending those lines to their intersection, not
            // by inventing an axis-aligned bounding-box corner. This candidate
            // is completely geometric and is still accepted only after exact
            // conservative containment and global added-area checks.
            const auto cross2 = [](const Vec2& first, const Vec2& second)
            { return first.x() * second.y() - first.y() * second.x(); };
            for (std::size_t first = 0; first + 3 < polygon.size(); ++first)
            {
                const Vec2 a = polygon[first];
                const Vec2 b = polygon[first + 1];
                const Vec2 c = polygon[first + 2];
                const Vec2 d = polygon[first + 3];
                const Vec2 first_direction = b - a;
                const Vec2 second_direction = d - c;
                const double bevel_length = (c - b).norm();
                if (bevel_length > 0.1 * std::min(
                        first_direction.norm(), second_direction.norm()))
                    continue;
                const double denominator = cross2(first_direction, second_direction);
                if (std::abs(denominator) <= tolerance *
                    std::max(first_direction.norm() * second_direction.norm(), 1.0))
                    continue;
                const double parameter = cross2(c - a, second_direction) / denominator;
                const Vec2 intersection = a + first_direction * parameter;
                if ((intersection - b).norm() + (intersection - c).norm() >
                    2.0 * (first_direction.norm() + second_direction.norm()))
                    continue;
                std::vector<Vec2> candidate;
                candidate.reserve(polygon.size() - 1);
                candidate.insert(candidate.end(), polygon.begin(),
                                 polygon.begin() + static_cast<std::ptrdiff_t>(first + 1));
                candidate.push_back(intersection);
                candidate.insert(candidate.end(),
                                 polygon.begin() + static_cast<std::ptrdiff_t>(first + 3),
                                 polygon.end());
                candidate = simplifyPolygon(std::move(candidate), tolerance);
                if (candidate.size() >= polygon.size() || candidate.size() < 3 ||
                    !simplePolygon(candidate, tolerance) ||
                    triangulatePolygon(candidate).size() + 2 != candidate.size())
                    continue;
                Clipper2Lib::PathD candidate_path;
                for (const Vec2& point : candidate)
                    candidate_path.emplace_back(point.x(), point.y());
                if (Clipper2Lib::Area(candidate_path) < 0.0)
                    std::reverse(candidate_path.begin(), candidate_path.end());
                const double missing_area = std::abs(Clipper2Lib::Area(
                    Clipper2Lib::Difference(
                        Clipper2Lib::PathsD{original_path},
                        Clipper2Lib::PathsD{candidate_path},
                        Clipper2Lib::FillRule::NonZero, clipper_precision)));
                if (missing_area > tolerance * tolerance) continue;
                const double added_area = std::max(
                    std::abs(Clipper2Lib::Area(candidate_path)) - original_area, 0.0);
                if ((spent_area + added_area) /
                        std::max(model_surface_area, 1.0e-30) >
                    maximum_added_area_ratio) continue;
                const std::size_t removed = polygon.size() - candidate.size();
                const double cost = added_area / static_cast<double>(removed);
                if (cost < best_cost ||
                    (cost == best_cost && removed > best_removed))
                {
                    best = std::move(candidate);
                    best_cost = cost;
                    best_removed = removed;
                }
            }
            if (best.empty()) break;
            spent_area += best_cost * static_cast<double>(best_removed);
            polygon = std::move(best);
            canonicalizeLoop(polygon);
        }
        canonicalizeLoop(polygon);
        primitive.polygon.clear();
        primitive.polygon.reserve(polygon.size());
        for (const Vec2& point : polygon)
            primitive.polygon.push_back(
                origin + frame.col(1) * point.x() + frame.col(2) * point.y());
    }
}

void synchronizeOpposingShellSilhouettes(
    std::vector<OutputPrimitive>& primitives,
    const double model_surface_area,
    const double maximum_added_area_ratio,
    const double tolerance)
{
    struct Candidate
    {
        std::size_t first = 0;
        std::size_t second = 0;
        double similarity = 0.0;
        double added_area = 0.0;
        Vec3 normal = Vec3::Zero();
        Vec3 origin = Vec3::Zero();
        Mat3 frame = Mat3::Identity();
        double second_distance = 0.0;
        std::vector<Vec2> boundary;
    };
    constexpr int clipper_precision = 8;
    std::vector<Candidate> candidates;
    for (std::size_t first = 0; first < primitives.size(); ++first)
    {
        // Only a shell moved by the preceding parallel-layer phases can have
        // acquired a silhouette that needs synchronization.  Restricting the
        // pair search to that event set changes the all-polygon cubic scan into
        // a narrow post-process without changing geometric acceptance.
        if (primitives[first].preserves_cavity_opening ||
            !primitives[first].shallow_shell_coalesced ||
            primitives[first].primitive.kind != Kind::Polygon) continue;
        Vec3 normal;
        if (!planarNormal(primitives[first].primitive, tolerance, normal)) continue;
        int dominant_axis = 0;
        for (int axis = 1; axis < 3; ++axis)
            if (std::abs(normal[axis]) > std::abs(normal[dominant_axis]))
                dominant_axis = axis;
        if (normal[dominant_axis] < 0.0) normal *= -1.0;
        const Vec3 origin = primitives[first].primitive.polygon.front();
        const Mat3 basis = orthonormalFrame(normal);
        Mat3 frame;
        frame.col(0) = basis.col(1);
        frame.col(1) = basis.col(2);
        frame.col(2) = normal;
        const auto project = [&](const Primitive& primitive)
        {
            Clipper2Lib::PathD path;
            for (const Vec3& vertex : primitive.polygon)
            {
                const Vec3 local = frame.transposeMultiply(vertex - origin);
                path.emplace_back(local.x(), local.y());
            }
            if (Clipper2Lib::Area(path) < 0.0) std::reverse(path.begin(), path.end());
            return path;
        };
        const auto boundsOf = [&](const Primitive& primitive)
        {
            Bounds2 bounds;
            for (const Vec3& vertex : primitive.polygon)
            {
                const Vec3 local = frame.transposeMultiply(vertex - origin);
                bounds.lower = bounds.lower.cwiseMin({local.x(), local.y()});
                bounds.upper = bounds.upper.cwiseMax({local.x(), local.y()});
            }
            return bounds;
        };
        const auto overlaps = [tolerance](const Bounds2& first_bounds,
                                           const Bounds2& second_bounds)
        {
            return first_bounds.upper.x() >= second_bounds.lower.x() - tolerance &&
                   first_bounds.lower.x() <= second_bounds.upper.x() + tolerance &&
                   first_bounds.upper.y() >= second_bounds.lower.y() - tolerance &&
                   first_bounds.lower.y() <= second_bounds.upper.y() + tolerance;
        };
        const auto first_path = project(primitives[first].primitive);
        const double first_area = std::abs(Clipper2Lib::Area(first_path));
        const Bounds2 first_bounds = boundsOf(primitives[first].primitive);
        for (std::size_t second = first + 1; second < primitives.size(); ++second)
        {
            if (primitives[second].preserves_cavity_opening ||
                !primitives[second].shallow_shell_coalesced ||
                primitives[second].primitive.kind != Kind::Polygon) continue;
            Vec3 second_normal;
            if (!planarNormal(primitives[second].primitive, tolerance, second_normal) ||
                std::abs(second_normal.dot(normal)) < 1.0 - 1.0e-8) continue;
            const double second_distance = normal.dot(
                primitives[second].primitive.polygon.front() - origin);
            if (std::abs(second_distance) <= tolerance) continue;
            const auto second_path = project(primitives[second].primitive);
            const auto united = Clipper2Lib::Union(
                Clipper2Lib::PathsD{first_path, second_path},
                Clipper2Lib::FillRule::NonZero, clipper_precision);
            if (united.size() != 1 || !Clipper2Lib::IsPositive(united.front())) continue;
            const double first_second_union_area = std::abs(Clipper2Lib::Area(united));
            const double intersection_area = std::abs(Clipper2Lib::Area(
                Clipper2Lib::Intersect(
                    Clipper2Lib::PathsD{first_path},
                    Clipper2Lib::PathsD{second_path},
                    Clipper2Lib::FillRule::NonZero, clipper_precision)));
            const double similarity = intersection_area /
                std::max(first_second_union_area, 1.0e-30);
            if (similarity < 0.95) continue;

            const Bounds2 second_bounds = boundsOf(primitives[second].primitive);
            std::size_t bridge_count = 0;
            const double lower_distance = std::min(0.0, second_distance);
            const double upper_distance = std::max(0.0, second_distance);
            for (std::size_t bridge = 0; bridge < primitives.size(); ++bridge)
            {
                if (bridge == first || bridge == second ||
                    primitives[bridge].primitive.kind != Kind::Polygon) continue;
                Vec3 bridge_normal;
                if (!planarNormal(primitives[bridge].primitive, tolerance, bridge_normal) ||
                    std::abs(bridge_normal.dot(normal)) > 1.0e-6) continue;
                double minimum = std::numeric_limits<double>::infinity();
                double maximum = -std::numeric_limits<double>::infinity();
                for (const Vec3& vertex : primitives[bridge].primitive.polygon)
                {
                    const double distance = normal.dot(vertex - origin);
                    minimum = std::min(minimum, distance);
                    maximum = std::max(maximum, distance);
                }
                if (minimum > lower_distance + tolerance ||
                    maximum < upper_distance - tolerance) continue;
                const Bounds2 bridge_bounds = boundsOf(primitives[bridge].primitive);
                if (overlaps(bridge_bounds, first_bounds) &&
                    overlaps(bridge_bounds, second_bounds))
                    ++bridge_count;
            }
            if (bridge_count < 2) continue;
            const double second_area = std::abs(Clipper2Lib::Area(second_path));
            const double added_area = std::max(
                2.0 * first_second_union_area - first_area - second_area, 0.0);
            if (added_area / std::max(model_surface_area, 1.0e-30) >
                maximum_added_area_ratio) continue;
            std::vector<Vec2> boundary;
            for (const auto& point : united.front())
                boundary.emplace_back(point.x, point.y);
            boundary = simplifyPolygon(std::move(boundary), tolerance);
            if (boundary.size() < 3 || !simplePolygon(boundary, tolerance) ||
                triangulatePolygon(boundary).size() + 2 != boundary.size()) continue;
            candidates.push_back({first, second, similarity, added_area, normal,
                                  origin, frame, second_distance, std::move(boundary)});
        }
    }
    std::sort(candidates.begin(), candidates.end(), [](const auto& first, const auto& second)
    {
        if (first.similarity != second.similarity)
            return first.similarity > second.similarity;
        if (first.added_area != second.added_area)
            return first.added_area < second.added_area;
        if (first.first != second.first) return first.first < second.first;
        return first.second < second.second;
    });
    std::vector<bool> used(primitives.size(), false);
    for (const Candidate& candidate : candidates)
    {
        if (used[candidate.first] || used[candidate.second]) continue;
        primitives[candidate.first].primitive = polygonPrimitive(
            candidate.boundary, candidate.origin, candidate.frame);
        primitives[candidate.second].primitive = polygonPrimitive(
            candidate.boundary,
            candidate.origin + candidate.normal * candidate.second_distance,
            candidate.frame);
        used[candidate.first] = true;
        used[candidate.second] = true;
    }
}

// Once two opposing outer patches have been assigned the same conservative
// silhouette, their old bridge patches may still describe the wider silhouette
// that existed before synchronization.  Keeping those patches is conservative,
// but it needlessly preserves a prism of empty space.  Rebuild a certified
// polygonal extrusion from the synchronized caps when all source responsibility
// carried by the replaced bridge patches lies inside that extrusion.
//
// Closure is not the acceptance criterion here.  The useful operation is the
// reduction of the implied swept volume at a small triangle-count cost; a closed
// side shell is merely the natural representation of that tighter extrusion.
std::vector<OutputPrimitive> rebuildCertifiedExtrudedShells(
    const Mesh& source_mesh,
    std::vector<OutputPrimitive> primitives,
    const double tolerance,
    std::vector<CertifiedExtrusion>* certified_extrusions = nullptr)
{
    constexpr int clipper_precision = 8;
    const auto projectionContained = [&](const Face& face,
                                         const Vec3& origin,
                                         const Mat3& frame,
                                         const std::vector<Vec2>& boundary,
                                         const Clipper2Lib::PathD& boundary_path,
                                         const double lower_distance,
                                         const double upper_distance)
    {
        std::array<Vec2, 3> projected;
        Clipper2Lib::PathD triangle;
        for (int corner = 0; corner < 3; ++corner)
        {
            const Vec3 local = frame.transposeMultiply(
                source_mesh.vertices[face[corner]] - origin);
            if (local.z() < lower_distance - tolerance ||
                local.z() > upper_distance + tolerance)
                return false;
            projected[corner] = {local.x(), local.y()};
            if (!pointInPolygon(projected[corner], boundary, tolerance)) return false;
            triangle.emplace_back(local.x(), local.y());
        }
        const double projected_area = std::abs(Clipper2Lib::Area(triangle));
        if (projected_area > tolerance * tolerance)
        {
            if (Clipper2Lib::Area(triangle) < 0.0) std::reverse(triangle.begin(), triangle.end());
            const double outside_area = std::abs(Clipper2Lib::Area(
                Clipper2Lib::Difference(
                    Clipper2Lib::PathsD{triangle},
                    Clipper2Lib::PathsD{boundary_path},
                    Clipper2Lib::FillRule::NonZero, clipper_precision)));
            return outside_area <= tolerance * tolerance;
        }
        // A wall triangle can collapse to a segment in cap projection.  End
        // points alone are insufficient for a concave cap, so reject any edge
        // that crosses the silhouette and also check its midpoint.
        for (int edge = 0; edge < 3; ++edge)
        {
            const Vec2 first = projected[edge];
            const Vec2 second = projected[(edge + 1) % 3];
            if (!pointInPolygon((first + second) * 0.5, boundary, tolerance))
                return false;
            for (std::size_t boundary_edge = 0;
                 boundary_edge < boundary.size(); ++boundary_edge)
                if (properSegmentIntersection(
                        first, second,
                        boundary[boundary_edge],
                        boundary[(boundary_edge + 1) % boundary.size()],
                        tolerance))
                    return false;
        }
        return true;
    };

    for (;;)
    {
        bool changed = false;
        for (std::size_t first = 0; first < primitives.size() && !changed; ++first)
        {
            if (primitives[first].preserves_cavity_opening ||
                !primitives[first].shallow_shell_coalesced ||
                primitives[first].primitive.kind != Kind::Polygon) continue;
            Vec3 normal;
            if (!planarNormal(primitives[first].primitive, tolerance, normal)) continue;
            int dominant_axis = 0;
            for (int axis = 1; axis < 3; ++axis)
                if (std::abs(normal[axis]) > std::abs(normal[dominant_axis]))
                    dominant_axis = axis;
            if (normal[dominant_axis] < 0.0) normal = -normal;
            const Vec3 origin = primitives[first].primitive.polygon.front();
            const Mat3 basis = orthonormalFrame(normal);
            Mat3 frame;
            frame.col(0) = basis.col(1);
            frame.col(1) = basis.col(2);
            frame.col(2) = normal;
            std::vector<Vec2> first_boundary;
            Clipper2Lib::PathD first_path;
            for (const Vec3& vertex : primitives[first].primitive.polygon)
            {
                const Vec3 local = frame.transposeMultiply(vertex - origin);
                first_boundary.emplace_back(local.x(), local.y());
                first_path.emplace_back(local.x(), local.y());
            }
            if (Clipper2Lib::Area(first_path) < 0.0)
            {
                std::reverse(first_path.begin(), first_path.end());
                std::reverse(first_boundary.begin(), first_boundary.end());
            }
            first_boundary = simplifyPolygon(std::move(first_boundary), tolerance);
            if (first_boundary.size() < 3 ||
                !simplePolygon(first_boundary, tolerance)) continue;
            first_path.clear();
            for (const Vec2& point : first_boundary)
                first_path.emplace_back(point.x(), point.y());

            for (std::size_t second = first + 1;
                 second < primitives.size() && !changed; ++second)
            {
                if (primitives[second].preserves_cavity_opening ||
                    !primitives[second].shallow_shell_coalesced ||
                    primitives[second].primitive.kind != Kind::Polygon) continue;
                Vec3 second_normal;
                if (!planarNormal(primitives[second].primitive, tolerance, second_normal) ||
                    std::abs(second_normal.dot(normal)) < 1.0 - 1.0e-8)
                    continue;
                const double cap_distance = normal.dot(
                    primitives[second].primitive.polygon.front() - origin);
                if (std::abs(cap_distance) <= tolerance) continue;
                Clipper2Lib::PathD second_path;
                for (const Vec3& vertex : primitives[second].primitive.polygon)
                {
                    const Vec3 local = frame.transposeMultiply(vertex - origin);
                    second_path.emplace_back(local.x(), local.y());
                }
                if (Clipper2Lib::Area(second_path) < 0.0)
                    std::reverse(second_path.begin(), second_path.end());
                const double union_area = std::abs(Clipper2Lib::Area(
                    Clipper2Lib::Union(
                        Clipper2Lib::PathsD{first_path, second_path},
                        Clipper2Lib::FillRule::NonZero, clipper_precision)));
                const double intersection_area = std::abs(Clipper2Lib::Area(
                    Clipper2Lib::Intersect(
                        Clipper2Lib::PathsD{first_path},
                        Clipper2Lib::PathsD{second_path},
                        Clipper2Lib::FillRule::NonZero, clipper_precision)));
                if (union_area <= tolerance * tolerance ||
                    union_area - intersection_area > tolerance * tolerance * 16.0)
                    continue;

                const double lower_distance = std::min(0.0, cap_distance);
                const double upper_distance = std::max(0.0, cap_distance);
                std::vector<std::size_t> bridges;
                std::vector<Vec2> old_projection_points = first_boundary;
                std::size_t old_triangles = 0;
                std::vector<std::uint32_t> bridge_source_faces;
                bool has_outside_projection = false;
                for (std::size_t bridge = 0; bridge < primitives.size(); ++bridge)
                {
                    if (bridge == first || bridge == second ||
                        primitives[bridge].primitive.kind != Kind::Polygon)
                        continue;
                    Vec3 bridge_normal;
                    if (!planarNormal(primitives[bridge].primitive, tolerance, bridge_normal) ||
                        std::abs(bridge_normal.dot(normal)) > 1.0e-6)
                        continue;
                    double minimum = std::numeric_limits<double>::infinity();
                    double maximum = -std::numeric_limits<double>::infinity();
                    bool projection_near_cap = true;
                    std::vector<Vec2> bridge_projection_points;
                    bool bridge_has_outside_projection = false;
                    for (const Vec3& vertex : primitives[bridge].primitive.polygon)
                    {
                        const Vec3 local = frame.transposeMultiply(vertex - origin);
                        minimum = std::min(minimum, local.z());
                        maximum = std::max(maximum, local.z());
                        const Vec2 projected(local.x(), local.y());
                        bridge_projection_points.push_back(projected);
                        if (!pointInPolygon(projected, first_boundary, tolerance))
                            bridge_has_outside_projection = true;
                    }
                    if (minimum > lower_distance + tolerance ||
                        maximum < upper_distance - tolerance)
                        continue;
                    if (primitives[bridge].source_faces.empty()) continue;
                    for (const auto face_id : primitives[bridge].source_faces)
                        if (!projectionContained(
                                source_mesh.faces[face_id], origin, frame,
                                first_boundary, first_path,
                                lower_distance, upper_distance))
                        {
                            projection_near_cap = false;
                            break;
                        }
                    if (!projection_near_cap) continue;
                    bridges.push_back(bridge);
                    old_projection_points.insert(
                        old_projection_points.end(), bridge_projection_points.begin(),
                        bridge_projection_points.end());
                    has_outside_projection |= bridge_has_outside_projection;
                    old_triangles += primitives[bridge].primitive.polygon.size() - 2;
                    bridge_source_faces.insert(
                        bridge_source_faces.end(),
                        primitives[bridge].source_faces.begin(),
                        primitives[bridge].source_faces.end());
                }
                if (bridges.size() < 3 || !has_outside_projection) continue;
                const std::size_t new_triangles = 2 * first_boundary.size();
                if (new_triangles > old_triangles + 2) continue;
                const auto old_hull = simplifyPolygon(
                    convexHull(std::move(old_projection_points), tolerance), tolerance);
                const double removed_cross_section = std::max(
                    std::abs(signedArea(old_hull)) -
                    std::abs(signedArea(first_boundary)), 0.0);
                if (removed_cross_section <= tolerance * tolerance) continue;

                std::sort(bridges.begin(), bridges.end());
                std::sort(bridge_source_faces.begin(), bridge_source_faces.end());
                bridge_source_faces.erase(
                    std::unique(bridge_source_faces.begin(), bridge_source_faces.end()),
                    bridge_source_faces.end());
                std::vector<OutputPrimitive> rebuilt;
                rebuilt.reserve(first_boundary.size());
                std::vector<std::vector<std::uint32_t>> assigned(first_boundary.size());
                for (const auto face_id : bridge_source_faces)
                {
                    const Face& face = source_mesh.faces[face_id];
                    Vec2 centroid;
                    for (const auto vertex_id : face)
                    {
                        const Vec3 local = frame.transposeMultiply(
                            source_mesh.vertices[vertex_id] - origin);
                        centroid += Vec2(local.x(), local.y());
                    }
                    centroid /= 3.0;
                    std::size_t nearest = 0;
                    double nearest_distance = std::numeric_limits<double>::infinity();
                    for (std::size_t edge = 0; edge < first_boundary.size(); ++edge)
                    {
                        const Vec2 a = first_boundary[edge];
                        const Vec2 b = first_boundary[(edge + 1) % first_boundary.size()];
                        const Vec2 direction = b - a;
                        const double denominator = std::max(direction.dot(direction), 1.0e-30);
                        const double parameter = std::clamp(
                            (centroid - a).dot(direction) / denominator, 0.0, 1.0);
                        const double distance = (centroid - (a + direction * parameter)).norm();
                        if (distance < nearest_distance)
                        {
                            nearest_distance = distance;
                            nearest = edge;
                        }
                    }
                    assigned[nearest].push_back(face_id);
                }
                for (std::size_t edge = 0; edge < first_boundary.size(); ++edge)
                {
                    const Vec2 a = first_boundary[edge];
                    const Vec2 b = first_boundary[(edge + 1) % first_boundary.size()];
                    const Vec3 first_a = origin + frame.col(0) * a.x() + frame.col(1) * a.y();
                    const Vec3 first_b = origin + frame.col(0) * b.x() + frame.col(1) * b.y();
                    const Vec3 second_b = first_b + normal * cap_distance;
                    const Vec3 second_a = first_a + normal * cap_distance;
                    rebuilt.push_back({polygonPrimitive(
                        {first_a, first_b, second_b, second_a}),
                        std::move(assigned[edge]), true});
                }
                if (certified_extrusions != nullptr)
                    certified_extrusions->push_back(
                        {origin, frame, first_boundary,
                         lower_distance, upper_distance});
                std::vector<OutputPrimitive> result;
                result.reserve(primitives.size() - bridges.size() + rebuilt.size());
                for (std::size_t index = 0; index < primitives.size(); ++index)
                    if (!std::binary_search(bridges.begin(), bridges.end(), index))
                        result.push_back(std::move(primitives[index]));
                result.insert(result.end(),
                    std::make_move_iterator(rebuilt.begin()),
                    std::make_move_iterator(rebuilt.end()));
                primitives = std::move(result);
                changed = true;
            }
        }
        if (!changed) break;
    }
    return primitives;
}

std::vector<OutputPrimitive> clipParallelOuterOcclusion(
    std::vector<OutputPrimitive> primitives,
    const Vec3& model_center,
    const double maximum_depth,
    const double tolerance)
{
    if (maximum_depth <= tolerance) return primitives;
    constexpr int clipper_precision = 8;
    std::vector<OutputPrimitive> result;
    result.reserve(primitives.size());
    for (std::size_t candidate_id = 0; candidate_id < primitives.size(); ++candidate_id)
    {
        OutputPrimitive& candidate = primitives[candidate_id];
        if (!candidate.preserves_cavity_opening ||
            candidate.primitive.kind != Kind::Polygon ||
            candidate.primitive.polygon.size() < 3)
        {
            result.push_back(candidate);
            continue;
        }
        Vec3 normal;
        if (!planarNormal(candidate.primitive, tolerance, normal))
        {
            result.push_back(candidate);
            continue;
        }
        int dominant_axis = 0;
        for (int axis = 1; axis < 3; ++axis)
            if (std::abs(normal[axis]) > std::abs(normal[dominant_axis]))
                dominant_axis = axis;
        if (normal[dominant_axis] < 0.0) normal = -normal;
        const Vec3 origin = candidate.primitive.polygon.front();
        const double candidate_distance = normal.dot(origin);
        const double center_distance = normal.dot(model_center);
        const double outward_sign = candidate_distance < center_distance ? -1.0 : 1.0;
        const Mat3 basis = orthonormalFrame(normal);
        Mat3 frame;
        frame.col(0) = static_cast<Vec3>(basis.col(1));
        frame.col(1) = static_cast<Vec3>(basis.col(2));
        frame.col(2) = normal;
        const auto project = [&](const Primitive& polygon)
        {
            Clipper2Lib::PathD path;
            path.reserve(polygon.polygon.size());
            for (const Vec3& vertex : polygon.polygon)
            {
                const Vec3 local = frame.transposeMultiply(vertex - origin);
                path.emplace_back(local.x(), local.y());
            }
            if (Clipper2Lib::Area(path) < 0.0) std::reverse(path.begin(), path.end());
            return path;
        };
        const Clipper2Lib::PathD candidate_path = project(candidate.primitive);
        Vec2 circular_center = Vec2::Zero();
        for (const auto& point : candidate_path)
            circular_center += Vec2(point.x, point.y);
        circular_center /= static_cast<double>(candidate_path.size());
        std::vector<Vec2> centered_candidate;
        centered_candidate.reserve(candidate_path.size());
        for (const auto& point : candidate_path)
            centered_candidate.push_back(
                Vec2(point.x, point.y) - circular_center);
        if (candidate_path.size() < 8 ||
            !ringLooksCircular(centered_candidate, 0.05))
        {
            result.push_back(candidate);
            continue;
        }

        Clipper2Lib::PathsD outward_coverers;
        for (std::size_t coverer_id = 0; coverer_id < primitives.size(); ++coverer_id)
        {
            if (coverer_id == candidate_id) continue;
            const Primitive& coverer = primitives[coverer_id].primitive;
            if (coverer.kind != Kind::Polygon || coverer.polygon.size() < 3) continue;
            Vec3 coverer_normal;
            if (!planarNormal(coverer, tolerance, coverer_normal) ||
                std::abs(coverer_normal.dot(normal)) < 1.0 - 1.0e-8) continue;
            const double coverer_distance = normal.dot(coverer.polygon.front());
            const double outward_depth =
                outward_sign * (coverer_distance - candidate_distance);
            if (outward_depth <= tolerance || outward_depth > maximum_depth) continue;
            outward_coverers.push_back(project(coverer));
        }
        if (outward_coverers.empty())
        {
            result.push_back(candidate);
            continue;
        }

        const Clipper2Lib::PathsD remaining = Clipper2Lib::Difference(
            Clipper2Lib::PathsD{candidate_path},
            Clipper2Lib::Union(outward_coverers, Clipper2Lib::FillRule::NonZero,
                               clipper_precision),
            Clipper2Lib::FillRule::NonZero, clipper_precision);
        if (std::any_of(remaining.begin(), remaining.end(),
            [](const auto& path) { return !Clipper2Lib::IsPositive(path); }))
        {
            // A polygon-with-hole cannot be represented by one semantic simple
            // polygon. Retain it until a certified hole-splitting path exists.
            result.push_back(candidate);
            continue;
        }
        for (const auto& path : remaining)
        {
            if (std::abs(Clipper2Lib::Area(path)) <= tolerance * tolerance) continue;
            std::vector<Vec2> boundary;
            boundary.reserve(path.size());
            for (const auto& point : path) boundary.emplace_back(point.x, point.y);
            boundary = simplifyPolygon(std::move(boundary), tolerance);
            if (boundary.size() < 3 ||
                triangulatePolygon(boundary).size() + 2 != boundary.size()) continue;
            result.push_back({polygonPrimitive(boundary, origin, frame),
                              candidate.source_faces});
        }
    }
    return result;
}

std::vector<Vec2> boxPlaneCrossSection(const BoxFit& box,
                                       const Vec3& origin,
                                       const Mat3& frame,
                                       const double tolerance,
                                       bool& occupies_negative_side,
                                       bool& occupies_positive_side)
{
    std::array<Vec3, 8> corners;
    std::array<double, 8> distances{};
    double minimum_distance = std::numeric_limits<double>::infinity();
    double maximum_distance = -std::numeric_limits<double>::infinity();
    for (int index = 0; index < 8; ++index)
    {
        Vec3 corner = box.center;
        for (int axis = 0; axis < 3; ++axis)
            corner += box.axes.col(axis) *
                ((index & (1 << axis)) ? box.half_size[axis] : -box.half_size[axis]);
        corners[index] = corner;
        distances[index] = frame.col(2).dot(corner - origin);
        minimum_distance = std::min(minimum_distance, distances[index]);
        maximum_distance = std::max(maximum_distance, distances[index]);
    }
    occupies_negative_side = minimum_distance < -tolerance;
    occupies_positive_side = maximum_distance > tolerance;
    if (minimum_distance > tolerance || maximum_distance < -tolerance ||
        (!occupies_negative_side && !occupies_positive_side)) return {};

    std::vector<Vec2> points;
    for (int axis = 0; axis < 3; ++axis)
        for (int index = 0; index < 8; ++index)
        {
            if (index & (1 << axis)) continue;
            const int other = index | (1 << axis);
            const double first_distance = distances[index];
            const double second_distance = distances[other];
            if ((first_distance < -tolerance && second_distance < -tolerance) ||
                (first_distance > tolerance && second_distance > tolerance)) continue;
            const double denominator = first_distance - second_distance;
            if (std::abs(denominator) <= 1.0e-30) continue;
            const double parameter = std::clamp(
                first_distance / denominator, 0.0, 1.0);
            const Vec3 point = corners[index] +
                (corners[other] - corners[index]) * parameter;
            const Vec3 local = frame.transposeMultiply(point - origin);
            const Vec2 projected(local.x(), local.y());
            if (std::none_of(points.begin(), points.end(), [&](const Vec2& existing)
                { return (existing - projected).norm() <= tolerance; }))
                points.push_back(projected);
        }
    if (points.size() < 3) return {};
    return simplifyPolygon(convexHull(std::move(points), tolerance), tolerance);
}

struct ClassifiedPlaneSections
{
    Clipper2Lib::PathsD negative;
    Clipper2Lib::PathsD positive;
    Clipper2Lib::PathsD crossing;
};

ClassifiedPlaneSections extrusionPlaneCrossSections(
    const CertifiedExtrusion& extrusion,
    const Vec3& plane_origin,
    const Mat3& plane_frame,
    const double tolerance)
{
    ClassifiedPlaneSections sections;
    const auto triangles = triangulatePolygon(extrusion.boundary);
    constexpr std::array<std::pair<int, int>, 9> edges{{
        {0, 1}, {1, 2}, {2, 0}, {3, 4}, {4, 5}, {5, 3},
        {0, 3}, {1, 4}, {2, 5},
    }};
    for (const auto& triangle : triangles)
    {
        std::array<Vec3, 6> vertices;
        std::array<double, 6> distances{};
        double minimum_distance = std::numeric_limits<double>::infinity();
        double maximum_distance = -std::numeric_limits<double>::infinity();
        for (int corner = 0; corner < 3; ++corner)
            for (int end = 0; end < 2; ++end)
            {
                const Vec2 point = extrusion.boundary[triangle[corner]];
                const double depth = end == 0 ? extrusion.lower_distance
                                              : extrusion.upper_distance;
                const int index = corner + 3 * end;
                vertices[index] = extrusion.origin +
                    extrusion.frame.col(0) * point.x() +
                    extrusion.frame.col(1) * point.y() +
                    extrusion.frame.col(2) * depth;
                distances[index] = plane_frame.col(2).dot(
                    vertices[index] - plane_origin);
                minimum_distance = std::min(minimum_distance, distances[index]);
                maximum_distance = std::max(maximum_distance, distances[index]);
            }
        const bool occupies_negative_side = minimum_distance < -tolerance;
        const bool occupies_positive_side = maximum_distance > tolerance;
        if (minimum_distance > tolerance || maximum_distance < -tolerance ||
            (!occupies_negative_side && !occupies_positive_side))
            continue;
        std::vector<Vec2> points;
        for (const auto [first, second] : edges)
        {
            const double first_distance = distances[first];
            const double second_distance = distances[second];
            if ((first_distance < -tolerance && second_distance < -tolerance) ||
                (first_distance > tolerance && second_distance > tolerance)) continue;
            const double denominator = first_distance - second_distance;
            if (std::abs(denominator) <= 1.0e-30) continue;
            const double parameter = std::clamp(
                first_distance / denominator, 0.0, 1.0);
            const Vec3 point = vertices[first] +
                (vertices[second] - vertices[first]) * parameter;
            const Vec3 local = plane_frame.transposeMultiply(point - plane_origin);
            const Vec2 projected(local.x(), local.y());
            if (std::none_of(points.begin(), points.end(), [&](const Vec2& existing)
                { return (existing - projected).norm() <= tolerance; }))
                points.push_back(projected);
        }
        if (points.size() < 3) continue;
        const auto polygon = simplifyPolygon(
            convexHull(std::move(points), tolerance), tolerance);
        if (polygon.size() < 3) continue;
        Clipper2Lib::PathD path;
        for (const Vec2& point : polygon) path.emplace_back(point.x(), point.y());
        if (Clipper2Lib::Area(path) < 0.0) std::reverse(path.begin(), path.end());
        if (occupies_negative_side && occupies_positive_side)
            sections.crossing.push_back(std::move(path));
        else if (occupies_negative_side)
            sections.negative.push_back(std::move(path));
        else
            sections.positive.push_back(std::move(path));
    }
    return sections;
}

std::vector<CertifiedExtrusion> recognizeCertifiedPrismaticVolumes(
    const std::vector<OutputPrimitive>& primitives,
    const double tolerance)
{
    constexpr int clipper_precision = 8;
    std::vector<CertifiedExtrusion> result;
    for (std::size_t first = 0; first < primitives.size(); ++first)
    {
        if (primitives[first].primitive.kind != Kind::Polygon) continue;
        Vec3 normal;
        if (!planarNormal(primitives[first].primitive, tolerance, normal)) continue;
        const Vec3 origin = primitives[first].primitive.polygon.front();
        const Mat3 basis = orthonormalFrame(normal);
        Mat3 frame;
        frame.col(0) = basis.col(1);
        frame.col(1) = basis.col(2);
        frame.col(2) = normal;
        const auto projectCap = [&](const Primitive& primitive)
        {
            Clipper2Lib::PathD path;
            for (const Vec3& vertex : primitive.polygon)
            {
                const Vec3 local = frame.transposeMultiply(vertex - origin);
                path.emplace_back(local.x(), local.y());
            }
            if (Clipper2Lib::Area(path) < 0.0) std::reverse(path.begin(), path.end());
            return path;
        };
        const auto first_path = projectCap(primitives[first].primitive);
        for (std::size_t second = first + 1; second < primitives.size(); ++second)
        {
            if (primitives[second].primitive.kind != Kind::Polygon) continue;
            Vec3 second_normal;
            if (!planarNormal(primitives[second].primitive, tolerance, second_normal) ||
                std::abs(second_normal.dot(normal)) < 1.0 - 1.0e-8) continue;
            const double cap_distance = normal.dot(
                primitives[second].primitive.polygon.front() - origin);
            if (std::abs(cap_distance) <= tolerance) continue;
            const auto second_path = projectCap(primitives[second].primitive);
            const auto intersections = Clipper2Lib::Intersect(
                Clipper2Lib::PathsD{first_path},
                Clipper2Lib::PathsD{second_path},
                Clipper2Lib::FillRule::NonZero, clipper_precision);
            for (const auto& intersection : intersections)
            {
                if (!Clipper2Lib::IsPositive(intersection) ||
                    std::abs(Clipper2Lib::Area(intersection)) <=
                        tolerance * tolerance) continue;
                std::vector<Vec2> boundary;
                for (const auto& point : intersection)
                    boundary.emplace_back(point.x, point.y);
                boundary = simplifyPolygon(std::move(boundary), tolerance);
                if (boundary.size() < 3 || !simplePolygon(boundary, tolerance)) continue;

                bool closed = true;
                for (std::size_t edge = 0; edge < boundary.size() && closed; ++edge)
                {
                    const Vec2 a = boundary[edge];
                    const Vec2 b = boundary[(edge + 1) % boundary.size()];
                    const Vec3 first_a = origin + frame.col(0) * a.x() +
                                         frame.col(1) * a.y();
                    const Vec3 first_b = origin + frame.col(0) * b.x() +
                                         frame.col(1) * b.y();
                    const Vec3 second_a = first_a + normal * cap_distance;
                    const Vec3 second_b = first_b + normal * cap_distance;
                    const Vec3 edge_vector = first_b - first_a;
                    const double edge_length = edge_vector.norm();
                    const double depth = std::abs(cap_distance);
                    if (edge_length <= tolerance || depth <= tolerance)
                    {
                        closed = false;
                        break;
                    }
                    const Vec3 edge_axis = edge_vector / edge_length;
                    const Vec3 depth_axis = (second_a - first_a) / depth;
                    const Vec3 side_normal = edge_axis.cross(depth_axis).normalized();
                    Clipper2Lib::PathsD side_coverers;
                    for (std::size_t side = 0; side < primitives.size(); ++side)
                    {
                        if (side == first || side == second ||
                            primitives[side].primitive.kind != Kind::Polygon) continue;
                        Vec3 side_candidate_normal;
                        if (!planarNormal(primitives[side].primitive, tolerance,
                                          side_candidate_normal) ||
                            std::abs(side_candidate_normal.dot(side_normal)) <
                                1.0 - 1.0e-8) continue;
                        Clipper2Lib::PathD path;
                        bool coplanar = true;
                        for (const Vec3& vertex : primitives[side].primitive.polygon)
                        {
                            const Vec3 relative = vertex - first_a;
                            if (std::abs(relative.dot(side_normal)) > tolerance * 16.0)
                            {
                                coplanar = false;
                                break;
                            }
                            path.emplace_back(relative.dot(edge_axis),
                                              relative.dot(depth_axis));
                        }
                        if (!coplanar || std::abs(Clipper2Lib::Area(path)) <=
                            tolerance * tolerance) continue;
                        if (Clipper2Lib::Area(path) < 0.0)
                            std::reverse(path.begin(), path.end());
                        side_coverers.push_back(std::move(path));
                    }
                    if (side_coverers.empty())
                    {
                        closed = false;
                        break;
                    }
                    Clipper2Lib::PathD target{
                        {0.0, 0.0}, {edge_length, 0.0},
                        {edge_length, depth}, {0.0, depth},
                    };
                    const auto uncovered = Clipper2Lib::Difference(
                        Clipper2Lib::PathsD{target},
                        Clipper2Lib::Union(
                            side_coverers, Clipper2Lib::FillRule::NonZero,
                            clipper_precision),
                        Clipper2Lib::FillRule::NonZero, clipper_precision);
                    if (std::abs(Clipper2Lib::Area(uncovered)) >
                        tolerance * tolerance * 16.0) closed = false;
                }
                if (closed)
                    result.push_back({origin, frame, boundary,
                                      std::min(0.0, cap_distance),
                                      std::max(0.0, cap_distance)});
            }
        }
    }
    return result;
}

std::vector<CertifiedExtrusion> recognizeEnclosureGroupExtrusions(
    const std::vector<OutputPrimitive>& primitives,
    const double tolerance,
    const bool include_ungrouped)
{
    std::unordered_map<std::uint64_t, std::vector<OutputPrimitive>> groups;
    std::vector<OutputPrimitive> ungrouped;
    for (const auto& item : primitives)
        if (item.enclosure_group != 0)
            groups[item.enclosure_group].push_back(item);
        else if (include_ungrouped)
            ungrouped.push_back(item);
    std::vector<CertifiedExtrusion> result;
    for (const auto& [group, members] : groups)
    {
        (void)group;
        auto recognized = recognizeCertifiedPrismaticVolumes(members, tolerance);
        for (auto& extrusion : recognized)
            extrusion.enclosure_group = group;
        result.insert(result.end(),
                      std::make_move_iterator(recognized.begin()),
                      std::make_move_iterator(recognized.end()));
    }
    if (include_ungrouped)
    {
        auto recognized_ungrouped = recognizeCertifiedPrismaticVolumes(
            ungrouped, tolerance);
        result.insert(result.end(),
                      std::make_move_iterator(recognized_ungrouped.begin()),
                      std::make_move_iterator(recognized_ungrouped.end()));
    }
    // A rectangular box is a valid extrusion along all three axes. The
    // recognizer therefore returns three certificates for the same solid;
    // submitting all of their nearly coincident plane sections to Clipper can
    // create a huge Boolean arrangement. Canonicalize by the quantized 3D
    // vertex set, which is orientation-independent and preserves distinct
    // overlapping solids.
    const double quantum = std::max(tolerance * 8.0, 1.0e-9);
    using VertexKey = std::array<std::int64_t, 3>;
    std::set<std::vector<VertexKey>> signatures;
    std::vector<CertifiedExtrusion> unique;
    unique.reserve(result.size());
    for (auto& extrusion : result)
    {
        std::vector<VertexKey> signature;
        signature.reserve(2 * extrusion.boundary.size());
        for (const double depth : {extrusion.lower_distance,
                                   extrusion.upper_distance})
            for (const Vec2& point : extrusion.boundary)
            {
                const Vec3 vertex = extrusion.origin +
                    extrusion.frame.col(0) * point.x() +
                    extrusion.frame.col(1) * point.y() +
                    extrusion.frame.col(2) * depth;
                signature.push_back({
                    static_cast<std::int64_t>(std::llround(vertex.x() / quantum)),
                    static_cast<std::int64_t>(std::llround(vertex.y() / quantum)),
                    static_cast<std::int64_t>(std::llround(vertex.z() / quantum)),
                });
            }
        std::sort(signature.begin(), signature.end());
        signature.erase(std::unique(signature.begin(), signature.end()),
                        signature.end());
        if (signatures.insert(std::move(signature)).second)
            unique.push_back(std::move(extrusion));
    }
    return unique;
}

std::vector<OutputPrimitive> clipPlanarOcclusionByClosedVolumes(
    std::vector<OutputPrimitive> primitives,
    const std::vector<BoxFit>& closed_volumes,
    const std::vector<CertifiedExtrusion>& closed_extrusions,
    const double tolerance,
    std::size_t& removed_count,
    std::vector<std::uint32_t>* excluded_faces,
    const bool preserve_triangle_workload,
    OcclusionClipStats* clip_stats)
{
    removed_count = 0;
    if (clip_stats) *clip_stats = {};
    if (closed_volumes.empty() && closed_extrusions.empty()) return primitives;
    if (clip_stats) clip_stats->input_triangles = triangulatedFaceCount(primitives);
    constexpr int clipper_precision = 8;
    std::vector<OutputPrimitive> result;
    result.reserve(primitives.size());
    for (auto& candidate : primitives)
    {
        if (candidate.preserves_cavity_opening ||
            candidate.primitive.kind != Kind::Polygon ||
            candidate.primitive.polygon.size() < 3)
        {
            result.push_back(std::move(candidate));
            continue;
        }
        Vec3 normal;
        if (!planarNormal(candidate.primitive, tolerance, normal))
        {
            result.push_back(std::move(candidate));
            continue;
        }
        const Vec3 origin = candidate.primitive.polygon.front();
        const Mat3 basis = orthonormalFrame(normal);
        Mat3 frame;
        frame.col(0) = basis.col(1);
        frame.col(1) = basis.col(2);
        frame.col(2) = normal;
        Clipper2Lib::PathD candidate_path;
        for (const Vec3& vertex : candidate.primitive.polygon)
        {
            const Vec3 local = frame.transposeMultiply(vertex - origin);
            candidate_path.emplace_back(local.x(), local.y());
        }
        if (Clipper2Lib::Area(candidate_path) < 0.0)
            std::reverse(candidate_path.begin(), candidate_path.end());

        Clipper2Lib::PathsD crossing_coverers;
        Clipper2Lib::PathsD owner_negative;
        Clipper2Lib::PathsD owner_positive;
        Clipper2Lib::PathsD other_negative;
        Clipper2Lib::PathsD other_positive;
        for (const BoxFit& volume : closed_volumes)
        {
            bool occupies_negative_side = false;
            bool occupies_positive_side = false;
            const auto cross_section = boxPlaneCrossSection(
                volume, origin, frame, tolerance * 8.0,
                occupies_negative_side, occupies_positive_side);
            if (cross_section.size() < 3) continue;
            Clipper2Lib::PathD path;
            for (const Vec2& point : cross_section)
                path.emplace_back(point.x(), point.y());
            if (Clipper2Lib::Area(path) < 0.0) std::reverse(path.begin(), path.end());
            if (occupies_negative_side && occupies_positive_side)
                crossing_coverers.push_back(std::move(path));
        }
        for (const CertifiedExtrusion& extrusion : closed_extrusions)
        {
            auto sections = extrusionPlaneCrossSections(
                extrusion, origin, frame, tolerance * 8.0);
            const bool is_owner = candidate.enclosure_group != 0 &&
                candidate.enclosure_group == extrusion.enclosure_group;
            const auto append = [](Clipper2Lib::PathsD& destination,
                                   Clipper2Lib::PathsD& source)
            {
                destination.insert(destination.end(),
                    std::make_move_iterator(source.begin()),
                    std::make_move_iterator(source.end()));
            };
            if (is_owner)
            {
                // A one-sided owner section identifies a true boundary side.
                // A crossing owner section proves that this candidate lies
                // strictly inside the already certified solid and is therefore
                // redundant even when both carry the same enclosure id.
                append(crossing_coverers, sections.crossing);
                append(owner_negative, sections.negative);
                append(owner_positive, sections.positive);
            }
            else
            {
                append(crossing_coverers, sections.crossing);
                append(other_negative, sections.negative);
                append(other_positive, sections.positive);
            }
        }
        using PointKey = std::pair<std::int64_t, std::int64_t>;
        constexpr double clipper_scale = 100000000.0;
        const auto unionUniqueSections = [&](Clipper2Lib::PathsD paths)
        {
            std::set<std::vector<PointKey>> section_signatures;
            Clipper2Lib::PathsD unique;
            unique.reserve(paths.size());
            for (auto& path : paths)
            {
                std::vector<PointKey> signature;
                signature.reserve(path.size());
                for (const auto& point : path)
                    signature.emplace_back(
                        static_cast<std::int64_t>(std::llround(
                            point.x * clipper_scale)),
                        static_cast<std::int64_t>(std::llround(
                            point.y * clipper_scale)));
                std::sort(signature.begin(), signature.end());
                signature.erase(std::unique(signature.begin(), signature.end()),
                                signature.end());
                if (signature.size() >= 3 &&
                    section_signatures.insert(std::move(signature)).second)
                    unique.push_back(std::move(path));
            }
            if (unique.empty()) return Clipper2Lib::PathsD{};
            return Clipper2Lib::Union(
                unique, Clipper2Lib::FillRule::NonZero, clipper_precision);
        };

        Clipper2Lib::PathsD internally_covered =
            unionUniqueSections(std::move(crossing_coverers));
        if (candidate.enclosure_group != 0)
        {
            const auto owner_negative_union =
                unionUniqueSections(std::move(owner_negative));
            const auto owner_positive_union =
                unionUniqueSections(std::move(owner_positive));
            const auto other_negative_union =
                unionUniqueSections(std::move(other_negative));
            const auto other_positive_union =
                unionUniqueSections(std::move(other_positive));
            const auto appendInterface = [&](const Clipper2Lib::PathsD& owner,
                                             const Clipper2Lib::PathsD& other)
            {
                if (owner.empty() || other.empty()) return;
                auto overlap = Clipper2Lib::Intersect(
                    owner, other, Clipper2Lib::FillRule::NonZero,
                    clipper_precision);
                internally_covered.insert(internally_covered.end(),
                    std::make_move_iterator(overlap.begin()),
                    std::make_move_iterator(overlap.end()));
            };
            appendInterface(owner_negative_union, other_positive_union);
            appendInterface(owner_positive_union, other_negative_union);
        }
        if (internally_covered.empty())
        {
            result.push_back(std::move(candidate));
            continue;
        }
        const auto covered = unionUniqueSections(std::move(internally_covered));
        const auto remaining = Clipper2Lib::Difference(
            Clipper2Lib::PathsD{candidate_path}, covered,
            Clipper2Lib::FillRule::NonZero, clipper_precision);
        const double removed_area = std::abs(Clipper2Lib::Area(candidate_path)) -
                                    std::abs(Clipper2Lib::Area(remaining));
        if (removed_area <= tolerance * tolerance)
        {
            result.push_back(std::move(candidate));
            continue;
        }
        if (remaining.empty())
        {
            if (excluded_faces != nullptr)
                excluded_faces->insert(excluded_faces->end(),
                    candidate.source_faces.begin(), candidate.source_faces.end());
            ++removed_count;
            if (clip_stats)
            {
                ++clip_stats->removed_primitives;
                ++clip_stats->clipped_primitives;
                clip_stats->removed_area += std::abs(
                    Clipper2Lib::Area(candidate_path));
            }
            continue;
        }
        // A hole or a split generally makes PQSS do more triangle work. Keep
        // the original conservative face unless clipping is no more expensive.
        const bool has_holes = std::any_of(
            remaining.begin(), remaining.end(),
            [](const auto& path) { return !Clipper2Lib::IsPositive(path); });
        if (has_holes && preserve_triangle_workload)
        {
            result.push_back(std::move(candidate));
            continue;
        }
        Clipper2Lib::PathsD output_paths = remaining;
        if (has_holes)
        {
            Clipper2Lib::PathsD triangles;
            if (Clipper2Lib::Triangulate(
                    remaining, clipper_precision, triangles, false) !=
                Clipper2Lib::TriangulateResult::success)
            {
                result.push_back(std::move(candidate));
                continue;
            }
            output_paths = std::move(triangles);
        }
        std::vector<std::vector<Vec2>> boundaries;
        std::size_t new_triangle_count = 0;
        bool valid = true;
        for (const auto& path : output_paths)
        {
            if (std::abs(Clipper2Lib::Area(path)) <= tolerance * tolerance) continue;
            std::vector<Vec2> boundary;
            for (const auto& point : path) boundary.emplace_back(point.x, point.y);
            boundary = simplifyPolygon(std::move(boundary), tolerance);
            const auto triangles = triangulatePolygon(boundary);
            if (boundary.size() < 3 || triangles.size() + 2 != boundary.size())
            {
                valid = false;
                break;
            }
            new_triangle_count += triangles.size();
            boundaries.push_back(std::move(boundary));
        }
        const std::size_t old_triangle_count =
            triangulatePolygon([&]
            {
                std::vector<Vec2> boundary;
                for (const auto& point : candidate_path)
                    boundary.emplace_back(point.x, point.y);
                return boundary;
            }()).size();
        if (!valid || boundaries.empty() ||
            (preserve_triangle_workload &&
             new_triangle_count > old_triangle_count))
        {
            result.push_back(std::move(candidate));
            continue;
        }
        for (const auto& boundary : boundaries)
        {
            OutputPrimitive clipped;
            clipped.primitive = polygonPrimitive(boundary, origin, frame);
            clipped.source_faces = candidate.source_faces;
            clipped.shallow_shell_coalesced = candidate.shallow_shell_coalesced;
            clipped.preserves_cavity_opening = candidate.preserves_cavity_opening;
            clipped.enclosure_group = candidate.enclosure_group;
            result.push_back(std::move(clipped));
        }
        if (clip_stats)
        {
            ++clip_stats->clipped_primitives;
            clip_stats->removed_area += std::max(removed_area, 0.0);
        }
    }
    if (clip_stats) clip_stats->output_triangles = triangulatedFaceCount(result);
    return result;
}

std::vector<OutputPrimitive> clipRegularizedPlanarOcclusion(
    const Mesh& mesh,
    std::vector<OutputPrimitive> primitives,
    const double tolerance,
    const std::vector<OutputPrimitive>& coverage_certificates,
    const std::vector<CertifiedExtrusion>* known_owner_extrusions,
    const std::vector<std::uint32_t>* excluded_redundant_faces,
    const std::size_t maximum_triangle_workload,
    std::size_t& removed_count,
    std::size_t& recognized_extrusions,
    OcclusionClipStats& clip_stats,
    bool& rolled_back)
{
    removed_count = 0;
    recognized_extrusions = 0;
    clip_stats = {};
    rolled_back = false;
    if (triangulatedFaceCount(primitives) > maximum_triangle_workload)
        return primitives;

    promoteToSemanticPrimitives(primitives);
    auto extrusions = recognizeCertifiedPrismaticVolumes(
        primitives, tolerance);
    if (known_owner_extrusions != nullptr)
        extrusions.insert(extrusions.end(), known_owner_extrusions->begin(),
                          known_owner_extrusions->end());
    recognized_extrusions = extrusions.size();
    if (extrusions.empty()) return primitives;

    std::vector<OutputPrimitive> before = primitives;
    primitives = clipPlanarOcclusionByClosedVolumes(
        std::move(primitives), {}, extrusions, tolerance, removed_count,
        nullptr, true, &clip_stats);
    if (clip_stats.clipped_primitives == 0) return primitives;

    const FinalCoverageAudit coverage = auditFinalConservativeCoverage(
        mesh, primitives, tolerance, &coverage_certificates, nullptr,
        excluded_redundant_faces);
    if (coverage.unassigned_source_faces == 0 &&
        coverage.failed_source_faces == 0)
        return primitives;

    removed_count = 0;
    clip_stats = {};
    rolled_back = true;
    return before;
}

std::vector<OutputPrimitive> fillLowErrorPlanarEnvelopes(
    std::vector<OutputPrimitive> primitives,
    const double model_surface_area,
    const double maximum_added_area_ratio,
    const double tolerance)
{
    constexpr int clipper_precision = 8;
    std::vector<bool> consumed(primitives.size(), false);
    std::vector<OutputPrimitive> result;
    result.reserve(primitives.size());
    for (std::size_t seed = 0; seed < primitives.size(); ++seed)
    {
        if (consumed[seed] || primitives[seed].primitive.kind != Kind::Polygon)
            continue;
        Vec3 normal;
        if (!planarNormal(primitives[seed].primitive, tolerance, normal)) continue;
        int dominant_axis = 0;
        for (int axis = 1; axis < 3; ++axis)
            if (std::abs(normal[axis]) > std::abs(normal[dominant_axis]))
                dominant_axis = axis;
        if (normal[dominant_axis] < 0.0) normal = -normal;
        const Vec3 origin = primitives[seed].primitive.polygon.front();
        const Mat3 basis = orthonormalFrame(normal);
        Mat3 frame;
        frame.col(0) = static_cast<Vec3>(basis.col(1));
        frame.col(1) = static_cast<Vec3>(basis.col(2));
        frame.col(2) = normal;

        std::vector<std::size_t> group;
        Clipper2Lib::PathsD subjects;
        std::vector<std::uint32_t> source_faces;
        for (std::size_t index = seed; index < primitives.size(); ++index)
        {
            if (consumed[index] || primitives[index].primitive.kind != Kind::Polygon)
                continue;
            Vec3 candidate_normal;
            if (!planarNormal(primitives[index].primitive, tolerance, candidate_normal) ||
                std::abs(candidate_normal.dot(normal)) < 1.0 - 1.0e-8 ||
                std::abs((primitives[index].primitive.polygon.front() - origin).dot(normal)) >
                    tolerance) continue;
            Clipper2Lib::PathD path;
            for (const Vec3& vertex : primitives[index].primitive.polygon)
            {
                const Vec3 local = frame.transposeMultiply(vertex - origin);
                path.emplace_back(local.x(), local.y());
            }
            if (std::abs(Clipper2Lib::Area(path)) <= tolerance * tolerance) continue;
            if (Clipper2Lib::Area(path) < 0.0) std::reverse(path.begin(), path.end());
            group.push_back(index);
            subjects.push_back(std::move(path));
            source_faces.insert(source_faces.end(),
                primitives[index].source_faces.begin(),
                primitives[index].source_faces.end());
        }
        if (subjects.empty()) continue;
        const Clipper2Lib::PathsD covered = Clipper2Lib::Union(
            subjects, Clipper2Lib::FillRule::NonZero, clipper_precision);
        const Clipper2Lib::RectD bounds = Clipper2Lib::GetBounds(covered);
        const double rectangle_area = bounds.Width() * bounds.Height();
        const double covered_area = std::abs(Clipper2Lib::Area(covered));
        const double added_ratio = std::max(rectangle_area - covered_area, 0.0) /
            std::max(model_surface_area, 1.0e-30);
        if (rectangle_area <= tolerance * tolerance ||
            added_ratio > maximum_added_area_ratio)
            continue;

        std::vector<Vec2> rectangle{
            {bounds.left, bounds.top},
            {bounds.left, bounds.bottom},
            {bounds.right, bounds.bottom},
            {bounds.right, bounds.top},
        };
        std::sort(source_faces.begin(), source_faces.end());
        source_faces.erase(std::unique(source_faces.begin(), source_faces.end()),
                           source_faces.end());
        result.push_back({polygonPrimitive(rectangle, origin, frame),
                          std::move(source_faces)});
        for (const auto index : group) consumed[index] = true;
    }
    for (std::size_t index = 0; index < primitives.size(); ++index)
        if (!consumed[index]) result.push_back(std::move(primitives[index]));
    return result;
}

std::vector<OutputPrimitive> regularizeShallowNearRectangles(
    std::vector<OutputPrimitive> primitives,
    const double model_surface_area,
    const double maximum_added_area_ratio,
    const double maximum_open_error_distance,
    const double tolerance)
{
    constexpr double maximum_local_added_area_ratio = 0.03;
    for (auto& item : primitives)
    {
        if (item.primitive.kind != Kind::Polygon ||
            item.primitive.polygon.size() <= 4 ||
            item.preserves_cavity_opening)
            continue;
        Vec3 normal;
        if (!planarNormal(item.primitive, tolerance, normal)) continue;
        const Vec3 origin = item.primitive.polygon.front();
        const Mat3 basis = orthonormalFrame(normal);
        Mat3 frame;
        frame.col(0) = basis.col(1);
        frame.col(1) = basis.col(2);
        frame.col(2) = normal;
        std::vector<Vec2> polygon;
        polygon.reserve(item.primitive.polygon.size());
        for (const Vec3& vertex : item.primitive.polygon)
        {
            const Vec3 local = frame.transposeMultiply(vertex - origin);
            polygon.emplace_back(local.x(), local.y());
        }
        std::size_t filled_boundary_voids = 0;
        double filled_boundary_void_area = 0.0;
        auto locally_closed = fillShortBoundaryVoids(
            polygon, tolerance, maximum_added_area_ratio, model_surface_area,
            maximum_open_error_distance, filled_boundary_voids,
            filled_boundary_void_area);
        if (filled_boundary_voids != 0 &&
            locally_closed.size() >= 3 &&
            simplePolygon(locally_closed, tolerance))
        {
            polygon = std::move(locally_closed);
            item.primitive = polygonPrimitive(polygon, origin, frame);
        }
        const double polygon_area = std::abs(signedArea(polygon));
        const auto hull = convexHull(polygon, tolerance);
        if (polygon_area <= tolerance * tolerance || hull.size() < 3) continue;

        std::vector<Vec2> best_rectangle;
        double best_area = std::numeric_limits<double>::infinity();
        double best_maximum_deviation = std::numeric_limits<double>::infinity();
        for (std::size_t edge = 0; edge < hull.size(); ++edge)
        {
            const Vec2 direction = hull[(edge + 1) % hull.size()] - hull[edge];
            if (direction.norm() <= tolerance) continue;
            const Vec2 first_axis = direction.normalized();
            const Vec2 second_axis(-first_axis.y(), first_axis.x());
            double minimum_first = std::numeric_limits<double>::infinity();
            double maximum_first = -std::numeric_limits<double>::infinity();
            double minimum_second = std::numeric_limits<double>::infinity();
            double maximum_second = -std::numeric_limits<double>::infinity();
            for (const Vec2& point : polygon)
            {
                const double first = point.dot(first_axis);
                const double second = point.dot(second_axis);
                minimum_first = std::min(minimum_first, first);
                maximum_first = std::max(maximum_first, first);
                minimum_second = std::min(minimum_second, second);
                maximum_second = std::max(maximum_second, second);
            }
            const double width = maximum_first - minimum_first;
            const double height = maximum_second - minimum_second;
            const double area = width * height;
            if (width <= tolerance || height <= tolerance || area >= best_area)
                continue;
            double maximum_deviation = 0.0;
            for (const Vec2& point : polygon)
            {
                const double first = point.dot(first_axis);
                const double second = point.dot(second_axis);
                maximum_deviation = std::max(maximum_deviation, std::min({
                    first - minimum_first, maximum_first - first,
                    second - minimum_second, maximum_second - second}));
            }
            best_area = area;
            best_maximum_deviation = maximum_deviation;
            best_rectangle = {
                first_axis * minimum_first + second_axis * minimum_second,
                first_axis * maximum_first + second_axis * minimum_second,
                first_axis * maximum_first + second_axis * maximum_second,
                first_axis * minimum_first + second_axis * maximum_second,
            };
        }
        if (best_rectangle.empty()) continue;
        const double added_area = std::max(best_area - polygon_area, 0.0);
        if (added_area <= tolerance * tolerance ||
            added_area / std::max(model_surface_area, 1.0e-30) >
                maximum_added_area_ratio ||
            added_area / polygon_area > maximum_local_added_area_ratio ||
            best_maximum_deviation >
                maximum_open_error_distance + tolerance)
            continue;

        OutputPrimitive replacement;
        replacement.primitive = polygonPrimitive(best_rectangle, origin, frame);
        replacement.source_faces = std::move(item.source_faces);
        replacement.shallow_shell_coalesced = item.shallow_shell_coalesced;
        replacement.preserves_cavity_opening = item.preserves_cavity_opening;
        replacement.enclosure_group = item.enclosure_group;
        item = std::move(replacement);
    }
    return primitives;
}

void writeSourceObj(const std::filesystem::path& path, const Mesh& mesh)
{
    std::ofstream stream(path);
    stream << std::setprecision(17);
    for (const Vec3& vertex : mesh.vertices)
        stream << "v " << vertex.x() << ' ' << vertex.y() << ' ' << vertex.z() << '\n';
    for (const Face& face : mesh.faces)
        stream << "f " << face[0] + 1 << ' ' << face[1] + 1 << ' ' << face[2] + 1 << '\n';
}

void writeRegionsObj(const std::filesystem::path& path,
                     const Mesh& mesh,
                     const std::vector<OutputPrimitive>& primitives,
                     const std::vector<std::uint32_t>& excluded_faces)
{
    std::ofstream stream(path);
    stream << std::setprecision(17);
    for (const Vec3& vertex : mesh.vertices)
        stream << "v " << vertex.x() << ' ' << vertex.y() << ' ' << vertex.z() << '\n';
    for (std::size_t id = 0; id < primitives.size(); ++id)
    {
        stream << "g region_" << std::setw(5) << std::setfill('0') << id << '_'
               << kindName(primitives[id].primitive.kind) << std::setfill(' ') << '\n';
        for (const auto face_id : primitives[id].source_faces)
        {
            const Face& face = mesh.faces[face_id];
            stream << "f " << face[0] + 1 << ' ' << face[1] + 1 << ' ' << face[2] + 1 << '\n';
        }
    }
    if (!excluded_faces.empty())
    {
        stream << "g excluded_redundant_surface\n";
        for (const auto face_id : excluded_faces)
        {
            const Face& face = mesh.faces[face_id];
            stream << "f " << face[0] + 1 << ' ' << face[1] + 1 << ' ' << face[2] + 1 << '\n';
        }
    }
}

void writeSemanticPrimitiveObj(const std::filesystem::path& path,
                               const std::vector<OutputPrimitive>& primitives)
{
    std::ofstream stream(path);
    stream << std::setprecision(17);
    std::uint32_t offset = 0;
    for (std::size_t id = 0; id < primitives.size(); ++id)
    {
        const Primitive& primitive = primitives[id].primitive;
        stream << "g primitive_" << std::setw(5) << std::setfill('0') << id << '_'
               << kindName(primitive.kind) << std::setfill(' ') << '\n';
        if (primitive.kind == Kind::Polygon)
        {
            for (const Vec3& vertex : primitive.polygon)
                stream << "v " << vertex.x() << ' ' << vertex.y() << ' ' << vertex.z() << '\n';
            stream << 'f';
            for (std::uint32_t index = 0; index < primitive.polygon.size(); ++index)
                stream << ' ' << offset + index + 1;
            stream << '\n';
            offset += static_cast<std::uint32_t>(primitive.polygon.size());
            continue;
        }
        const PrimitiveMesh mesh = triangulatePrimitive(primitive);
        for (const Vec3& vertex : mesh.vertices)
            stream << "v " << vertex.x() << ' ' << vertex.y() << ' ' << vertex.z() << '\n';
        for (const Face& face : mesh.faces)
            stream << "f " << offset + face[0] + 1 << ' ' << offset + face[1] + 1 << ' '
                   << offset + face[2] + 1 << '\n';
        offset += static_cast<std::uint32_t>(mesh.vertices.size());
    }
}

std::size_t writeTriangulatedObj(const std::filesystem::path& path,
                                 const std::vector<OutputPrimitive>& primitives)
{
    std::ofstream stream(path);
    stream << std::setprecision(17);
    std::uint32_t offset = 0;
    std::size_t triangle_count = 0;
    Bounds bounds;
    for (const auto& item : primitives)
        for (const Vec3& vertex : triangulatePrimitive(item.primitive).vertices)
        {
            bounds.lower = bounds.lower.cwiseMin(vertex);
            bounds.upper = bounds.upper.cwiseMax(vertex);
        }
    const double diagonal = (bounds.upper - bounds.lower).norm();
    const double area_epsilon_squared = std::max(
        std::pow(diagonal, 4.0) * 1.0e-28, 1.0e-48);
    const double quantum = std::max(diagonal * 1.0e-10, 1.0e-12);
    using Coordinate = std::array<std::int64_t, 3>;
    using TriangleKey = std::array<Coordinate, 3>;
    std::set<TriangleKey> emitted_triangles;
    for (std::size_t id = 0; id < primitives.size(); ++id)
    {
        const PrimitiveMesh mesh = triangulatePrimitive(primitives[id].primitive);
        stream << "g primitive_" << std::setw(5) << std::setfill('0') << id << '_'
               << kindName(primitives[id].primitive.kind) << std::setfill(' ') << '\n';
        for (const Vec3& vertex : mesh.vertices)
            stream << "v " << vertex.x() << ' ' << vertex.y() << ' ' << vertex.z() << '\n';
        for (const Face& face : mesh.faces)
        {
            const Vec3& first = mesh.vertices[face[0]];
            const Vec3& second = mesh.vertices[face[1]];
            const Vec3& third = mesh.vertices[face[2]];
            const Vec3 cross = (second - first).cross(third - first);
            if (cross.squaredNorm() <= area_epsilon_squared)
                continue;
            TriangleKey key{{
                Coordinate{{std::llround(first.x() / quantum),
                            std::llround(first.y() / quantum),
                            std::llround(first.z() / quantum)}},
                Coordinate{{std::llround(second.x() / quantum),
                            std::llround(second.y() / quantum),
                            std::llround(second.z() / quantum)}},
                Coordinate{{std::llround(third.x() / quantum),
                            std::llround(third.y() / quantum),
                            std::llround(third.z() / quantum)}},
            }};
            std::sort(key.begin(), key.end());
            if (!emitted_triangles.insert(key).second) continue;
            stream << "f " << offset + face[0] + 1 << ' '
                   << offset + face[1] + 1 << ' '
                   << offset + face[2] + 1 << '\n';
            ++triangle_count;
        }
        offset += static_cast<std::uint32_t>(mesh.vertices.size());
    }
    return triangle_count;
}

void writeMetadata(const std::filesystem::path& directory,
                   const std::filesystem::path& source_path,
                   const std::vector<OutputPrimitive>& primitives,
                   const PrimitiveMeshAnalysisStats& stats,
                   const PrimitiveMeshAnalysisOptions& options)
{
    std::ofstream model(directory / "model.json");
    model << std::setprecision(17)
          << "{\"stats\":{\"model\":\"" << source_path.filename().string()
          << "\",\"source_triangles\":" << stats.source_triangles
          << ",\"discarded_degenerate_triangles\":" << stats.discarded_degenerate_triangles
          << ",\"primitive_count\":" << stats.primitive_count
          << ",\"primitive_types\":{\"polygon\":" << stats.polygon_count
          << ",\"disk\":" << stats.disk_count
          << ",\"annulus\":" << stats.annulus_count
          << ",\"cylindricalband\":" << stats.cylindrical_band_count
          << ",\"conicalband\":" << stats.conical_band_count
          << "},\"triangulated_proxy_triangles\":" << stats.proxy_triangles
          << ",\"proxy_triangles\":" << stats.proxy_triangles
          << ",\"filled_planar_holes\":" << stats.filled_planar_holes
          << ",\"filled_cavity_volume_ratio\":" << stats.filled_cavity_volume_ratio
          << ",\"filled_boundary_voids\":" << stats.filled_boundary_voids
          << ",\"filled_intercomponent_gaps\":"
          << stats.filled_intercomponent_gaps
          << ",\"filled_boundary_void_area\":" << stats.filled_boundary_void_area
          << ",\"removed_contained_primitives\":" << stats.removed_contained_primitives
          << ",\"removed_sealed_void_wall_primitives\":"
          << stats.removed_sealed_void_wall_primitives
          << ",\"excluded_sealed_void_wall_triangles\":"
          << stats.excluded_sealed_void_wall_triangles
          << ",\"removed_blind_cavity_primitives\":"
          << stats.removed_blind_cavity_primitives
          << ",\"recognized_closed_box_shells\":" << stats.recognized_closed_box_shells
          << ",\"recognized_protrusion_box_shells\":"
          << stats.recognized_protrusion_box_shells
          << ",\"merged_local_planar_primitives\":"
          << stats.merged_local_planar_primitives
          << ",\"merged_spatial_primitive_groups\":"
          << stats.merged_spatial_primitive_groups
          << ",\"canonicalized_coplanar_groups\":"
          << stats.canonicalized_coplanar_groups
          << ",\"removed_coplanar_redundant_primitives\":"
          << stats.removed_coplanar_redundant_primitives
          << ",\"removed_coplanar_overlap_area\":"
          << stats.removed_coplanar_overlap_area
          << ",\"coverage_audit\":{\"assigned_source_faces\":"
          << stats.coverage_assigned_source_faces
          << ",\"enclosure_source_faces\":"
          << stats.coverage_enclosure_source_faces
          << ",\"planar_source_faces\":"
          << stats.coverage_planar_source_faces
          << ",\"unassigned_source_faces\":"
          << stats.coverage_unassigned_source_faces
          << ",\"failed_source_faces\":"
          << stats.coverage_failed_source_faces
          << ",\"seconds\":" << stats.coverage_audit_seconds << '}'
          << ",\"containment_validation\":{\"passed\":"
          << (stats.containment_validation_passed ? "true" : "false")
          << ",\"method\":\"triangle_coverage_certificate\""
          << ",\"source_triangles\":" << stats.source_triangles
          << ",\"assigned_source_triangles\":"
          << stats.coverage_assigned_source_faces
          << ",\"unassigned_source_triangles\":"
          << stats.coverage_unassigned_source_faces
          << ",\"failed_source_triangles\":"
          << stats.coverage_failed_source_faces << '}'
           << ",\"simplification_error\":{\"distance_method\":"
              "\"sampled_directed_phase3_to_phase1_surface_distance\""
           << ",\"sampling_method\":"
              "\"deterministic_area_surface_with_vertices_and_edges\""
           << ",\"reference\":\"phase1_hole_filled.obj\""
           << ",\"maximum_is_sample_estimate\":true"
           << ",\"maximum_distance_limit\":"
          << stats.maximum_open_error_distance_limit
          << ",\"limit_was_default\":"
          << (options.maximum_open_error_distance < 0.0 ? "true" : "false")
           << ",\"distance_sample_count\":"
           << stats.open_error_distance_sample_count
           << ",\"mean_distance\":" << stats.open_mean_distance
          << ",\"maximum_distance\":" << stats.open_max_distance
          << ",\"mean_distance_ratio\":"
          << stats.open_mean_distance_ratio
          << ",\"maximum_distance_ratio\":"
          << stats.open_max_distance_ratio
          << ",\"maximum_pair\":{\"proxy\":["
          << stats.open_max_proxy_point[0] << ','
          << stats.open_max_proxy_point[1] << ','
          << stats.open_max_proxy_point[2] << "],\"source\":["
          << stats.open_max_source_point[0] << ','
          << stats.open_max_source_point[1] << ','
          << stats.open_max_source_point[2] << "]},\"seconds\":"
          << stats.open_error_audit_seconds << '}'
          << ",\"minimum_protrusion_candidate_area_excess_ratio\":"
          << stats.minimum_protrusion_candidate_area_excess_ratio
          << ",\"selected_envelope_added_volume_ratio\":"
          << stats.selected_envelope_added_volume_ratio
          << ",\"selected_envelope_axis\":" << stats.selected_envelope_axis
          << ",\"envelope_candidate_added_volume_ratios\":["
          << stats.envelope_candidate_added_volume_ratios[0] << ','
          << stats.envelope_candidate_added_volume_ratios[1] << ','
          << stats.envelope_candidate_added_volume_ratios[2] << ']'
          << ",\"envelope_candidate_primitive_counts\":["
          << stats.envelope_candidate_primitive_counts[0] << ','
          << stats.envelope_candidate_primitive_counts[1] << ','
          << stats.envelope_candidate_primitive_counts[2] << ']'
          << ",\"envelope_candidate_remaining_cavities\":["
          << stats.envelope_candidate_remaining_cavities[0] << ','
          << stats.envelope_candidate_remaining_cavities[1] << ','
          << stats.envelope_candidate_remaining_cavities[2] << ']'
          << ",\"uniform_structure_policy\":"
          << (options.uniform_structure_policy ? "true" : "false")
          << ",\"timings_seconds\":{\"total\":"
          << stats.analysis_seconds << "}},\"source\":\"source.obj\""
          << ",\"phase1_hole_filled\":\"phase1_hole_filled.obj\""
          << ",\"phase2_recognized_surfaces\":\"phase2_recognized_surfaces.obj\""
          << ",\"phase3_simplified_surfaces\":\"primitives.obj\""
          << ",\"phase4_triangulated\":\"proxy.obj\""
          << ",\"regions\":\"regions.obj\""
          << ",\"primitive_analysis\":\"primitives.obj\""
          << ",\"triangulated_proxy\":\"proxy.obj\""
          << ",\"proxy\":\"proxy.obj\""
          << ",\"open_error_visualization\":\"open_error.json\""
          << ",\"proxy_components\":[";
    for (std::size_t id = 0; id < primitives.size(); ++id)
    {
        if (id) model << ',';
        const Primitive& primitive = primitives[id].primitive;
        const PrimitiveMesh triangulated = triangulatePrimitive(primitive);
        model << "{\"id\":" << id << ",\"type\":\""
              << kindName(primitive.kind) << "\",\"vertex_count\":"
              << (primitive.kind == Kind::Polygon ? primitive.polygon.size() : 0)
              << ",\"triangulated_face_count\":" << triangulated.faces.size()
              << ",\"source_face_count\":" << primitives[id].source_faces.size()
              << ",\"enclosure_group\":" << primitives[id].enclosure_group
              << "}";
    }
    model << "]}\n";

    std::ofstream manifest(directory / "viewer_manifest.json");
    std::string model_id = source_path.stem().string();
    manifest << "{\"algorithm\":\"CppProxyMeshGeneration\"," 
             << "\"complete\":true,\"model_count\":1,\"models\":[{\"id\":"
             << model_id << ",\"metadata\":\"model.json\"}]}\n";
}

} // namespace

double analysisStrengthToAllowedExcessRatio(const double strength)
{
    if (strength < 0.0 || strength > 1.0)
        throw std::invalid_argument("analysis_strength must be in [0, 1]");

    // IEEE-754 division makes the exponent +infinity at strength == 0.
    return std::pow(10.0, 3.0 / strength - 8.0);
}

PrimitiveMeshAnalysisStats analyzePrimitiveMeshObj(
    const std::filesystem::path& input_obj,
    const std::filesystem::path& output_directory,
    const PrimitiveMeshAnalysisOptions& options)
{
    const auto started = std::chrono::steady_clock::now();
    std::filesystem::create_directories(output_directory);
    const auto markStage = [&](const char* name)
    {
        std::ofstream progress(output_directory / "analysis_stages.txt", std::ios::app);
        progress << name << ' ' << std::setprecision(10)
                 << std::chrono::duration<double>(
                        std::chrono::steady_clock::now() - started).count()
                 << '\n';
    };
    {
        std::ofstream progress(output_directory / "analysis_stages.txt");
        progress << "start 0\n";
    }
    {
        std::ofstream profile(output_directory / "stage_error_profile.jsonl");
    }
    PrimitiveMeshAnalysisOptions effective_options = options;

    if (!effective_options.allow_polygon)
        throw std::invalid_argument("polygon must be enabled for exact fallback");
    if (effective_options.round_surface_segments < 3)
        throw std::invalid_argument("round_surface_segments must be at least 3");
    if (effective_options.analysis_strength < 0.0 || effective_options.analysis_strength > 1.0)
        throw std::invalid_argument("analysis_strength must be in [0, 1]");

    Mesh mesh = readObj(input_obj);
    weldCoincidentVertices(mesh);
    PrimitiveMeshAnalysisStats stats;
    stats.source_triangles = mesh.faces.size();
    stats.discarded_degenerate_triangles = dropDegenerateFaces(mesh);
    stats.source_triangles = mesh.faces.size();
    markStage("input_preprocessed");

    Vec3 lower = Vec3::Constant(std::numeric_limits<double>::infinity());
    Vec3 upper = Vec3::Constant(-std::numeric_limits<double>::infinity());
    for (const Vec3& vertex : mesh.vertices)
    {
        lower = lower.cwiseMin(vertex);
        upper = upper.cwiseMax(vertex);
    }
    const Vec3 extent = upper - lower;
    const double diagonal = extent.norm();
    const double model_volume = std::max(extent.prod(), diagonal * diagonal * diagonal * 1.0e-12);
    const double model_surface_area = facesArea(mesh, [&]
    {
        std::vector<std::uint32_t> ids(mesh.faces.size());
        std::iota(ids.begin(), ids.end(), 0);
        return ids;
    }());
    const double threshold = effective_options.maximum_added_volume_ratio >= 0.0
        ? effective_options.maximum_added_volume_ratio
        : analysisStrengthToAllowedExcessRatio(
            effective_options.uniform_structure_policy ? 0.24 : effective_options.analysis_strength);
    const double maximum_open_error_distance =
        effective_options.maximum_open_error_distance >= 0.0
            ? effective_options.maximum_open_error_distance
            : 100.0;
    stats.maximum_open_error_distance_limit = maximum_open_error_distance;

    // There is deliberately no enclosing-volume candidate phase.  Every
    // candidate created below is one surface patch.
    markStage("surface_candidate_setup");
    StructuralCleanup structural_cleanup = identifyStructuralRedundantFaces(
        mesh, effective_options, diagonal, model_surface_area, model_volume);
    stats.removed_sealed_void_wall_primitives =
        structural_cleanup.sealed_void_wall_primitives;
    stats.excluded_sealed_void_wall_triangles =
        structural_cleanup.sealed_void_wall_triangles;
    stats.removed_blind_cavity_primitives =
        structural_cleanup.blind_cavity_primitives;
    stats.removed_contained_primitives = structural_cleanup.contained_primitives;
    markStage("structural_cleanup");

    std::vector<bool> responsibility_faces(mesh.faces.size(), true);
    // Keep the exterior partition identical to the fully closed candidate.
    // Over-budget cavity faces are analyzed later in an isolated pass; letting
    // them re-enter this graph changes otherwise-correct outer side merges.
    for (const auto face : structural_cleanup.candidate_excluded_faces)
        responsibility_faces[face] = false;

    std::vector<OutputPrimitive> output;
    std::vector<BoxFit> certified_closed_volumes;
    std::vector<CertifiedExtrusion> certified_closed_extrusions;
    double analytic_filled_cavity_volume = 0.0;
    if (effective_options.allow_round_surfaces)
    {
        // Certify a complete revolved CAD component before hole filling or
        // planar partitioning changes its responsibility graph.  This is a
        // strict geometric test: protrusions, internal walls, or an off-axis
        // recess make the certificate fail and the component continues through
        // the ordinary surface pipeline.  A true disk-like or frustum-like
        // body is therefore represented by its analytic side and end surfaces
        // instead of thousands of nearly coplanar wedges.
        std::vector<bool> all_faces(mesh.faces.size(), true);
        const double analytic_tolerance = diagonal *
            effective_options.analytic_surface_relative_tolerance;
        for (const auto& component : faceComponentsApproximate(
                 mesh, all_faces, analytic_tolerance))
        {
            const auto fit = fitCertifiedRevolvedSurface(
                mesh, component, effective_options, analytic_tolerance);
            if (!fit) continue;
            const std::size_t previous_size = output.size();
            appendRevolvedSurfacePatches(
                output, mesh, *fit, component, analytic_tolerance);
            if (output.size() == previous_size) continue;
            for (const auto face : component)
                responsibility_faces[face] = false;
        }
        CylindricalRegionExtractionStats cylindrical_stats;
        auto cylindrical_regions = extractCylindricalSurfaceRegions(
            mesh, responsibility_faces, effective_options,
            std::max(analytic_tolerance, 1.0e-10), cylindrical_stats);
        for (const OutputPrimitive& region : cylindrical_regions)
            for (const auto face : region.source_faces)
                responsibility_faces[face] = false;
        output.insert(output.end(),
            std::make_move_iterator(cylindrical_regions.begin()),
            std::make_move_iterator(cylindrical_regions.end()));
        {
            std::ofstream profile(output_directory /
                                  "cylindrical_region_profile.json");
            profile << "{\"candidate_axes\":"
                    << cylindrical_stats.candidate_axes
                    << ",\"smooth_components\":"
                    << cylindrical_stats.smooth_components
                    << ",\"rejected_small\":"
                    << cylindrical_stats.rejected_small
                    << ",\"rejected_circle\":"
                    << cylindrical_stats.rejected_circle
                    << ",\"rejected_angular_coverage\":"
                    << cylindrical_stats.rejected_angular_coverage
                    << ",\"rejected_workload\":"
                    << cylindrical_stats.rejected_workload
                    << ",\"accepted_regions\":"
                    << cylindrical_stats.accepted_regions
                    << ",\"accepted_source_faces\":"
                    << cylindrical_stats.accepted_source_faces << "}\n";
        }
        // Surface recognition precedes exact coplanar clustering. It sees the
        // complete connected CAD patch, including all of its inner boundary
        // loops, so small bores cannot disappear into thousands of fallback
        // triangles before the global fill-volume decision is made.
        for (const auto& component : approximatePlanarRegions(
                 mesh, responsibility_faces, effective_options, diagonal))
        {
            auto surface = classifyApproximateCircularSurface(
                mesh, component, effective_options, (lower + upper) * 0.5,
                diagonal, model_volume, stats.filled_planar_holes,
                analytic_filled_cavity_volume);
            if (surface.empty())
                surface = classifyApproximatePlanarSurface(
                    mesh, component, (lower + upper) * 0.5, diagonal,
                    effective_options);
            if (surface.empty()) continue;
            std::vector<std::uint32_t> recognized_faces;
            for (const OutputPrimitive& item : surface)
                recognized_faces.insert(recognized_faces.end(),
                                        item.source_faces.begin(),
                                        item.source_faces.end());
            output.insert(output.end(), std::make_move_iterator(surface.begin()),
                          std::make_move_iterator(surface.end()));
            for (const auto face : recognized_faces)
                responsibility_faces[face] = false;
        }
    }
    markStage("analytic_surface_recognition");
    stats.filled_cavity_volume_ratio = analytic_filled_cavity_volume /
        std::max(model_volume, 1.0e-30);
    if (effective_options.allow_round_surfaces)
    {
        // Recognize only analytic surface patches.  A shared axis/radius fit may
        // identify a lateral band and its planar end faces, but no closed
        // closed round object is ever inserted into the candidate set.
        const double analytic_tolerance = diagonal *
            effective_options.analytic_surface_relative_tolerance;
        for (const auto& component : faceComponentsApproximate(
                 mesh, responsibility_faces, analytic_tolerance))
        {
            const auto fit = fitCertifiedRevolvedSurface(
                mesh, component, effective_options, analytic_tolerance);
            if (!fit) continue;
            appendRevolvedSurfacePatches(
                output, mesh, *fit, component, analytic_tolerance);
            for (const auto face : component) responsibility_faces[face] = false;
        }
    }
    requireSurfaceCandidates(output, "analytic surface recognition");
    std::vector<std::unordered_set<std::uint32_t>> adjacency;
    auto clusters = coplanarClusters(
        mesh, diagonal * effective_options.coplanar_relative_tolerance,
        adjacency, &responsibility_faces);
    for (const auto& cluster : clusters)
    {
        auto classified = classifyFinalRegion(
            mesh, cluster, effective_options, threshold, model_surface_area,
            stats.filled_planar_holes, stats.filled_boundary_voids,
            stats.filled_boundary_void_area);
        output.insert(output.end(),
            std::make_move_iterator(classified.begin()),
            std::make_move_iterator(classified.end()));
    }
    refineRoundSurfaceSegments(output, maximum_open_error_distance);
    decltype(adjacency){}.swap(adjacency);
    decltype(clusters){}.swap(clusters);
    output = fillCertifiedIntercomponentGaps(
        mesh, std::move(output), model_volume,
        effective_options.maximum_cavity_added_volume_ratio,
        std::max(diagonal * 1.0e-9, 1.0e-10),
        stats.filled_intercomponent_gaps, certified_closed_volumes,
        output_directory / "intercomponent_gap_profile.json");
    markStage("region_classification");

    // Phase 1 reference: make every accepted hole/cavity decision explicit.
    // Large cavities are restored here, while accepted fills keep their cap and
    // omit the now-occluded inner wall. Later simplification error is measured
    // against this surface, so filling a hole itself is intentionally free.
    std::vector<OutputPrimitive> restored_cavity_output;
    if (!structural_cleanup.restored_cavity_faces.empty())
    {
        std::vector<bool> cavity_mask(mesh.faces.size(), false);
        for (const auto face : structural_cleanup.restored_cavity_faces)
            cavity_mask[face] = true;
        std::vector<std::unordered_set<std::uint32_t>> cavity_adjacency;
        const auto cavity_clusters = coplanarClusters(
            mesh, diagonal * effective_options.coplanar_relative_tolerance,
            cavity_adjacency, &cavity_mask);
        for (const auto& cluster : cavity_clusters)
        {
            auto classified = classifyFinalRegion(
                mesh, cluster, effective_options, threshold,
                model_surface_area,
                stats.filled_planar_holes, stats.filled_boundary_voids,
                stats.filled_boundary_void_area);
            for (OutputPrimitive& item : classified)
                item.preserves_cavity_opening = true;
            restored_cavity_output.insert(restored_cavity_output.end(),
                std::make_move_iterator(classified.begin()),
                std::make_move_iterator(classified.end()));
        }
        CoplanarCanonicalizationStats cavity_coplanar_stats;
        restored_cavity_output = canonicalizeCoplanarPrimitiveUnion(
            mesh, std::move(restored_cavity_output),
            std::max(diagonal * 1.0e-9, 1.0e-10), cavity_coplanar_stats);
        promoteToSemanticPrimitives(restored_cavity_output);
    }
    const auto withRestoredCavities = [&](std::vector<OutputPrimitive> stage)
    {
        stage.insert(stage.end(), restored_cavity_output.begin(),
                     restored_cavity_output.end());
        stage = reopenRestoredCavityVolumes(
            mesh, std::move(stage), structural_cleanup,
            effective_options.maximum_cavity_added_volume_ratio, model_volume,
            std::max(diagonal * 1.0e-9, 1.0e-10));
        CoplanarCanonicalizationStats stage_coplanar_stats;
        stage = canonicalizeCoplanarPrimitiveUnion(
            mesh, std::move(stage), std::max(diagonal * 1.0e-9, 1.0e-10),
            stage_coplanar_stats);
        return stage;
    };
    std::vector<OutputPrimitive> phase1_hole_filled =
        withRestoredCavities(output);
    PrimitiveMesh phase1_triangulated =
        triangulateOutputPrimitives(phase1_hole_filled);
    (void)writeTriangulatedObj(
        output_directory / "phase1_hole_filled.obj", phase1_hole_filled);
    Mesh filled_surface_mesh;
    filled_surface_mesh.vertices = std::move(phase1_triangulated.vertices);
    filled_surface_mesh.faces = std::move(phase1_triangulated.faces);
    std::vector<OutputPrimitive>().swap(phase1_hole_filled);
    markStage("phase1_hole_filled_reference");

    CoplanarCanonicalizationStats coplanar_stats;
    output = canonicalizeCoplanarPrimitiveUnion(
        mesh, std::move(output), std::max(diagonal * 1.0e-9, 1.0e-10),
        coplanar_stats);
    stats.canonicalized_coplanar_groups = coplanar_stats.groups;
    stats.removed_coplanar_redundant_primitives = coplanar_stats.removed_primitives;
    stats.removed_coplanar_overlap_area = coplanar_stats.removed_overlap_area;
    std::vector<OutputPrimitive> phase2_recognized_surfaces =
        withRestoredCavities(output);
    promoteToSemanticPrimitives(phase2_recognized_surfaces);
    requireSurfaceCandidates(
        phase2_recognized_surfaces, "phase 2 surface recognition");
    writeSemanticPrimitiveObj(
        output_directory / "phase2_recognized_surfaces.obj",
        phase2_recognized_surfaces);
    // Stage 3 must simplify the complete stage-2 surface, including every
    // cavity opening that phase 1 deliberately retained.  Keeping those
    // restored surfaces on a side channel and appending them after Hausdorff
    // merging made them bypass the user's error limit: even an effectively
    // unlimited limit could not replace the whole model by one box shell.
    output = std::move(phase2_recognized_surfaces);
    std::size_t disk_occluded_surfaces = 0;
    output = removePlanarSurfacesOccludedByDisks(
        std::move(output), (lower + upper) * 0.5,
        std::max(diagonal * 1.0e-9, 1.0e-10), disk_occluded_surfaces);
    stats.removed_contained_primitives += disk_occluded_surfaces;
    std::vector<std::uint32_t> excluded_faces =
        std::move(structural_cleanup.excluded_faces);
    structural_cleanup = StructuralCleanup{};
    std::vector<OutputPrimitive>().swap(restored_cavity_output);
    std::vector<bool>().swap(responsibility_faces);
    markStage("exact_surface_union");

    const double merge_tolerance =
        std::max(diagonal * 1.0e-9, 1.0e-10);

    // A non-coplanar protruding component cannot in general be represented by
    // one planar patch.  Detect a complete component from a dominant support
    // plane and emit only its five exposed box faces.  This avoids arbitrary
    // spatial subtrees cutting the main body into unrelated local boxes.
    std::size_t support_surface_set_merges = 0;
    std::vector<OutputPrimitive> occluded_support_certificates;
    const double support_surface_error_limit = std::min(
        maximum_open_error_distance, diagonal * 0.02);
    output = mergeSupportProtrusionSurfaceSets(
        mesh, filled_surface_mesh, std::move(output), effective_options,
        diagonal, model_surface_area, merge_tolerance,
        support_surface_error_limit,
        std::max(diagonal / 192.0, 1.0e-30),
        support_surface_set_merges,
        occluded_support_certificates,
        output_directory / "support_protrusion_merge_profile.json");
    stats.merged_spatial_primitive_groups = support_surface_set_merges;
    {
        std::ofstream profile(
            output_directory / "spatial_group_fixed_point_profile.json");
        profile << "{\"passes\":1,\"surface_only\":true"
                << ",\"requested_limit\":" << std::setprecision(17)
                << maximum_open_error_distance
                << ",\"accepted_groups\":" << support_surface_set_merges
                << ",\"removed_enclosed_primitives\":0"
                << ",\"output_primitives\":" << output.size() << "}\n";
    }
    output = mergeAdjacentSurfacePrimitives(
        mesh, filled_surface_mesh, std::move(output),
        std::max(diagonal * 1.0e-9, 1.0e-10),
        maximum_open_error_distance,
        std::max(diagonal / 192.0, 1.0e-30),
        stats.merged_local_planar_primitives,
        output_directory / "surface_merge_profile.json");
    requireSurfaceCandidates(output, "phase 3 surface merging");
    std::size_t adjacent_envelope_group_merges = 0;
    std::vector<OutputPrimitive> convex_enclosure_certificates;
    output = mergeAdjacentEnvelopeGroups(
        mesh, filled_surface_mesh, std::move(output), effective_options,
        merge_tolerance,
        maximum_open_error_distance,
        std::max(diagonal / 192.0, 1.0e-30),
        adjacent_envelope_group_merges,
        convex_enclosure_certificates,
        output_directory / "adjacent_envelope_group_profile.json");
    stats.merged_spatial_primitive_groups += adjacent_envelope_group_merges;
    // Boolean clipping may leave a zero-area loop that carries duplicate or
    // collinear coordinates. It has no collision surface and cannot cover a
    // nondegenerate source triangle, so remove it before the occlusion passes;
    // the final responsibility audit remains authoritative.
    output.erase(std::remove_if(output.begin(), output.end(),
        [](const OutputPrimitive& item)
        { return triangulatePrimitive(item.primitive).faces.empty(); }),
        output.end());
    requireSurfaceCandidates(output, "phase 3 closed envelope merging");

    // Preserve the exact responsibility owners produced by the fixed-point
    // merge before deleting surfaces hidden inside enclosure groups. A removed
    // surface can prove coverage of its own source faces; transferring those
    // faces to an arbitrary surviving enclosure face does not imply that the
    // survivor's finite polygon covers the same projection.
    const std::size_t protected_cavity_walls =
        protectLargeOpposingCavityWalls(
            output,
            effective_options.maximum_cavity_added_volume_ratio * model_volume,
            diagonal, merge_tolerance);
    {
        std::ofstream profile(
            output_directory / "cavity_wall_protection_profile.json");
        profile << "{\"protected_primitives\":"
                << protected_cavity_walls << "}\n";
    }
    std::vector<OutputPrimitive> coverage_certificate_primitives = output;
    coverage_certificate_primitives.insert(
        coverage_certificate_primitives.end(),
        std::make_move_iterator(occluded_support_certificates.begin()),
        std::make_move_iterator(occluded_support_certificates.end()));
    coverage_certificate_primitives.insert(
        coverage_certificate_primitives.end(),
        std::make_move_iterator(convex_enclosure_certificates.begin()),
        std::make_move_iterator(convex_enclosure_certificates.end()));

    // Remove inner parallel layers from the same current surface set before
    // any of their outer coverers are discarded.  Coverage is computed from
    // the union of all same-facing outward patches, so a bottom face hidden by
    // three side-panel bottoms is handled as one Boolean occlusion event, just
    // like the support face hidden by several protrusion contacts.
    ParallelOcclusionStats parallel_occlusion_stats;
    const std::vector<OutputPrimitive> pre_parallel_occlusion_output = output;
    const std::size_t parallel_certificate_begin =
        coverage_certificate_primitives.size();
    output = clipParallelInternalSurfaceOcclusion(
        std::move(output), (lower + upper) * 0.5,
        diagonal * 0.03, merge_tolerance,
        coverage_certificate_primitives, parallel_occlusion_stats);
    // The newly appended certificates are the actual faces removed in this
    // pass.  Use only those witnesses for the orthogonal propagation step;
    // feeding every active primitive here would turn ordinary adjacent panels
    // into clipping planes and could erase legitimate exterior surfaces.
    std::vector<OutputPrimitive> parallel_certificates;
    if (parallel_certificate_begin < coverage_certificate_primitives.size())
        parallel_certificates.insert(
            parallel_certificates.end(),
            coverage_certificate_primitives.begin() + parallel_certificate_begin,
            coverage_certificate_primitives.end());
    // Keep the surviving parallel coverers paired with the removed witnesses.
    // A removed lower skin may be separated from its surviving support plane
    // by a thin wall; orthogonal neighbours must be clipped at that surviving
    // plane, not at the already removed skin itself.
    std::vector<OutputPrimitive> orthogonal_occlusion_planes =
        parallel_certificates;
    const auto active_end = coverage_certificate_primitives.begin() +
        static_cast<std::ptrdiff_t>(parallel_certificate_begin);
    for (auto active = coverage_certificate_primitives.begin();
         active != active_end; ++active)
    {
        Vec3 active_normal;
        if (!planarNormal(active->primitive, merge_tolerance, active_normal))
            continue;
        const double active_distance = active_normal.dot(
            planarPoint(active->primitive));
        bool paired = false;
        for (const auto& certificate : parallel_certificates)
        {
            Vec3 certificate_normal;
            if (!planarNormal(certificate.primitive, merge_tolerance,
                              certificate_normal) ||
                std::abs(active_normal.dot(certificate_normal)) < 1.0 - 1.0e-8)
                continue;
            if (active_normal.dot(certificate_normal) < 0.0)
                certificate_normal = -certificate_normal;
            const double certificate_distance = certificate_normal.dot(
                planarPoint(certificate.primitive));
            if (std::abs(active_distance - certificate_distance) > merge_tolerance &&
                std::abs(active_distance - certificate_distance) <= diagonal * 0.03)
            {
                paired = true;
                break;
            }
        }
        if (paired)
        {
            OutputPrimitive plane = *active;
            // Internal marker for the orthogonal pass: this plane is a
            // surviving coverer paired with a removed witness, so its finite
            // footprint need not contain the intersection endpoint exactly
            // (the endpoint may lie on a shared CAD seam).
            plane.preserves_cavity_opening = true;
            orthogonal_occlusion_planes.push_back(std::move(plane));
        }
    }
    const std::size_t orthogonal_input_count = output.size();
    std::vector<OutputPrimitive> orthogonal_input = output;
    std::vector<OutputPrimitive> orthogonal_output =
        clipAdjacentFaceTrianglesAtOcclusionPlanes(
            std::move(output), orthogonal_occlusion_planes,
            (lower + upper) * 0.5, merge_tolerance,
            parallel_occlusion_stats);
    // Orthogonal propagation is deliberately transactional.  A malformed CAD
    // seam can make a plane witness appear to cover unrelated patches; if one
    // pass would erase more than a fifth of the active surface set, retain the
    // pre-pass geometry and let the conservative coverage audit decide the
    // remaining cleanup.  This keeps the generic rule from becoming a hidden
    // model-specific aggressive simplifier.
    if (orthogonal_output.size() * 5 < orthogonal_input_count * 4)
        output = std::move(orthogonal_input);
    else
        output = std::move(orthogonal_output);
    // Occlusion is an optimization, never the authority for conservative
    // coverage. Roll the entire pass back if it loses even one source face;
    // this avoids reaching the expensive per-face exact repair path.
    const FinalCoverageAudit post_occlusion_coverage =
        auditFinalConservativeCoverage(
            mesh, output, merge_tolerance, &coverage_certificate_primitives,
            nullptr, &excluded_faces);
    if (post_occlusion_coverage.unassigned_source_faces != 0 ||
        post_occlusion_coverage.failed_source_faces != 0)
        output = pre_parallel_occlusion_output;
    {
        std::ofstream parallel_profile(
            output_directory / "parallel_occlusion_profile.json");
        parallel_profile << "{\"candidate_count\":"
            << parallel_occlusion_stats.candidate_count
            << ",\"accepted_candidates\":"
            << parallel_occlusion_stats.accepted_candidates
            << ",\"clipped_primitives\":"
            << parallel_occlusion_stats.clipped_primitives
            << ",\"removed_primitives\":"
            << parallel_occlusion_stats.removed_primitives
            << ",\"output_fragments\":"
            << parallel_occlusion_stats.output_fragments
            << ",\"removed_area\":"
            << parallel_occlusion_stats.removed_area
            << ",\"maximum_covered_ratio\":"
            << parallel_occlusion_stats.maximum_covered_ratio
            << ",\"accepted_ids\":[";
        for (std::size_t index = 0; index <
             parallel_occlusion_stats.accepted_ids.size(); ++index)
        {
            if (index != 0) parallel_profile << ',';
            parallel_profile << parallel_occlusion_stats.accepted_ids[index];
        }
        parallel_profile << "],\"candidate_ratios\":[";
        for (std::size_t index = 0; index <
             parallel_occlusion_stats.candidate_ratios.size(); ++index)
        {
            if (index != 0) parallel_profile << ',';
            parallel_profile << "["
                << parallel_occlusion_stats.candidate_ratios[index].first
                << "," << parallel_occlusion_stats.candidate_ratios[index].second
                << "]";
        }
        parallel_profile << "]}\n";
    }

    // Run the generic prismatic-occlusion pass before enclosure bookkeeping.
    // This is intentionally model-independent: if the current surface set
    // contains a certified open/closed extrusion, planar caps hidden inside it
    // are removed using the same coverage certificate used by the final audit.
    // The workload ceiling prevents this cleanup from trading a small amount of
    // redundant surface for an unbounded triangle increase.
    std::size_t regularized_occluded = 0;
    std::size_t regularized_extrusions = 0;
    OcclusionClipStats regularized_clip_stats;
    bool regularized_rolled_back = false;
    const std::size_t regularized_workload = triangulatedFaceCount(output);
    output = clipRegularizedPlanarOcclusion(
        mesh, std::move(output), merge_tolerance,
        coverage_certificate_primitives, nullptr, &excluded_faces,
        regularized_workload, regularized_occluded,
        regularized_extrusions, regularized_clip_stats,
        regularized_rolled_back);
    stats.removed_contained_primitives += regularized_occluded;
    {
        std::ofstream regularized_profile(
            output_directory / "regularized_occlusion_profile.json");
        regularized_profile << "{\"recognized_extrusions\":"
            << regularized_extrusions
            << ",\"removed_primitives\":" << regularized_occluded
            << ",\"rolled_back\":"
            << (regularized_rolled_back ? "true" : "false")
            << ",\"clipped_primitives\":"
            << regularized_clip_stats.clipped_primitives
            << ",\"input_triangles\":"
            << regularized_clip_stats.input_triangles
            << ",\"output_triangles\":"
            << regularized_clip_stats.output_triangles << "}\n";
    }

    std::size_t locally_enclosed_primitives = 0;
    output = removePrimitivesInsideEnclosureGroups(
        mesh, std::move(output), merge_tolerance,
        locally_enclosed_primitives);
    stats.removed_contained_primitives += locally_enclosed_primitives;

    CoplanarCanonicalizationStats post_merge_coplanar_stats;
    output = canonicalizeCoplanarPrimitiveUnion(
        mesh, std::move(output), std::max(diagonal * 1.0e-9, 1.0e-10),
        post_merge_coplanar_stats);
    stats.canonicalized_coplanar_groups += post_merge_coplanar_stats.groups;
    stats.removed_coplanar_redundant_primitives +=
        post_merge_coplanar_stats.removed_primitives;
    stats.removed_coplanar_overlap_area +=
        post_merge_coplanar_stats.removed_overlap_area;
    markStage("hausdorff_surface_merging");

    // Closed-box recognition is bookkeeping only. Its six rectangles remain the
    // actual primitives; there is no box primitive in the type system.
    std::vector<RecognizedProtrusion> discovered_closed_boxes;
    recognizeClosedAxisAlignedBoxes(
        output, std::max(diagonal * 1.0e-9, 1.0e-10),
        stats.recognized_closed_box_shells, &discovered_closed_boxes);
    for (const auto& box : discovered_closed_boxes)
        certified_closed_volumes.push_back(box.box);
    std::size_t post_merge_contained = 0;
    (void)removeContainedPrimitives(
        output, std::max(diagonal * 1.0e-9, 1.0e-10), post_merge_contained);
    stats.removed_contained_primitives += post_merge_contained;
    // Containment cleanup can expose coplanar shell fragments that were not in
    // the same earlier canonicalization group. The exported semantic layer is
    // required to be overlap-free, so canonicalize once at the final boundary.
    CoplanarCanonicalizationStats final_coplanar_stats;
    output = canonicalizeCoplanarPrimitiveUnion(
        mesh, std::move(output), std::max(diagonal * 1.0e-9, 1.0e-10),
        final_coplanar_stats);
    stats.canonicalized_coplanar_groups += final_coplanar_stats.groups;
    stats.removed_coplanar_redundant_primitives +=
        final_coplanar_stats.removed_primitives;
    stats.removed_coplanar_overlap_area +=
        final_coplanar_stats.removed_overlap_area;
    markStage("containment_and_internal_canonicalization");
    stats.excluded_redundant_triangles = excluded_faces.size();
    promoteToSemanticPrimitives(output);
    // From this boundary onward only topology-preserving cleanup is allowed.
    // Approximate terrace, silhouette, and extrusion rewrites previously used
    // area/volume budgets and could silently override the Hausdorff decision.
    CoplanarCanonicalizationStats export_coplanar_stats;
    output = canonicalizeCoplanarPrimitiveUnion(
        mesh, std::move(output), std::max(diagonal * 1.0e-9, 1.0e-10),
        export_coplanar_stats);
    stats.canonicalized_coplanar_groups += export_coplanar_stats.groups;
    stats.removed_coplanar_redundant_primitives +=
        export_coplanar_stats.removed_primitives;
    stats.removed_coplanar_overlap_area +=
        export_coplanar_stats.removed_overlap_area;
    markStage("export_canonicalization");

    // Restored cavity surfaces already participated in stage 3. Do not append
    // or reopen them here: doing so would undo a valid error-bounded merge and
    // reintroduce geometry that the user explicitly allowed us to remove.
    CoplanarCanonicalizationStats final_surface_coplanar_stats;
    output = canonicalizeCoplanarPrimitiveUnion(
        mesh, std::move(output), std::max(diagonal * 1.0e-9, 1.0e-10),
        final_surface_coplanar_stats);
    markStage("final_surface_union");
    // Rebuild the occlusion certificates from the active fixed-point result.
    // Historical accepted boxes include candidates that were later swallowed
    // by a larger enclosure. Feeding both generations to Clipper creates a
    // dense arrangement of overlapping sections and can exhaust the 2 GB job
    // limit on large CAD meshes.
    certified_closed_extrusions = recognizeEnclosureGroupExtrusions(
        coverage_certificate_primitives,
        std::max(diagonal * 1.0e-9, 1.0e-10), false);
    const std::vector<BoxFit> active_closed_volumes =
        selectActiveBoxCertificates(
            coverage_certificate_primitives, certified_closed_volumes,
            std::max(diagonal * 1.0e-9, 1.0e-10));
    {
        std::size_t boundary_vertices = 0;
        for (const auto& extrusion : certified_closed_extrusions)
            boundary_vertices += extrusion.boundary.size();
        std::ofstream profile(output_directory /
                              "final_occlusion_certificate_profile.json");
        profile << "{\"historical_box_count\":"
                << certified_closed_volumes.size()
                << ",\"active_box_count\":"
                << active_closed_volumes.size()
                << ",\"active_extrusion_count\":"
                << certified_closed_extrusions.size()
                << ",\"active_boundary_vertices\":"
                << boundary_vertices << "}\n";
    }
    markStage("final_volume_recognition");
    std::size_t volume_occluded_primitives = 0;
    const std::vector<OutputPrimitive> pre_volume_occlusion_output = output;
    output = clipPlanarOcclusionByClosedVolumes(
        std::move(output), active_closed_volumes, certified_closed_extrusions,
        std::max(diagonal * 1.0e-9, 1.0e-10), volume_occluded_primitives,
        nullptr, false);
    const FinalCoverageAudit post_volume_coverage =
        auditFinalConservativeCoverage(
            mesh, output, merge_tolerance, &coverage_certificate_primitives,
            &certified_closed_extrusions, &excluded_faces);
    if (post_volume_coverage.unassigned_source_faces != 0 ||
        post_volume_coverage.failed_source_faces != 0)
    {
        output = pre_volume_occlusion_output;
        volume_occluded_primitives = 0;
    }
    stats.removed_contained_primitives += volume_occluded_primitives;
    markStage("final_volume_occlusion");
    const std::vector<OutputPrimitive> pre_outer_occlusion_output = output;
    output = clipParallelOuterOcclusion(
        std::move(output), (lower + upper) * 0.5, diagonal * 0.03,
        std::max(diagonal * 1.0e-9, 1.0e-10));
    const FinalCoverageAudit post_outer_coverage =
        auditFinalConservativeCoverage(
            mesh, output, merge_tolerance, &coverage_certificate_primitives,
            &certified_closed_extrusions, &excluded_faces);
    if (post_outer_coverage.unassigned_source_faces != 0 ||
        post_outer_coverage.failed_source_faces != 0)
        output = pre_outer_occlusion_output;
    std::vector<OutputPrimitive> simplification_error_primitives = output;
    markStage("final_surface_canonicalization");

    const double final_tolerance = std::max(
        diagonal * 1.0e-9, 1.0e-10);
    FinalCoverageAudit coverage = auditFinalConservativeCoverage(
        mesh, output, final_tolerance, &coverage_certificate_primitives,
        &certified_closed_extrusions, &excluded_faces);
    markStage("coverage_audit_pre_repair");
    const std::size_t pre_repair_primitives = output.size();
    const std::size_t pre_repair_triangles = triangulatedFaceCount(output);
    const std::size_t pre_repair_unassigned = coverage.unassigned_source_faces;
    const std::size_t pre_repair_failed = coverage.failed_source_faces;
    const std::size_t pre_repair_failed_with_planar_owner =
        coverage.failed_with_planar_owner;
    const std::size_t pre_repair_failed_with_certificate_planar_owner =
        coverage.failed_with_certificate_planar_owner;
    const std::size_t pre_repair_failed_with_certificate_band_owner =
        coverage.failed_with_certificate_band_owner;
    std::size_t repair_merged_groups = 0;
    std::size_t repair_output_primitives = 0;
    std::size_t repair_output_triangles = 0;
    std::unordered_set<std::uint32_t> repaired_face_set;
    std::unordered_set<std::uint32_t> exact_repaired_face_set;
    std::size_t repair_pass = 0;
    while (!coverage.failed_face_ids.empty())
    {
        std::vector<std::uint32_t> newly_failed_faces;
        newly_failed_faces.reserve(coverage.failed_face_ids.size());
        for (const auto face_id : coverage.failed_face_ids)
        {
            if (face_id >= mesh.faces.size()) continue;
            repaired_face_set.insert(face_id);
            newly_failed_faces.push_back(face_id);
        }
        std::size_t merged_groups = 0;
        std::vector<OutputPrimitive> repair_output;
        if (repair_pass++ == 0)
        {
            repair_output = buildNonOverlappingPlanarCoverageRepair(
                mesh, output, newly_failed_faces, final_tolerance,
                merged_groups);
            repair_merged_groups += merged_groups;
        }
        else
            for (const auto face_id : newly_failed_faces)
            {
                if (!exact_repaired_face_set.insert(face_id).second) continue;
                const Face& face = mesh.faces[face_id];
                Primitive triangle;
                triangle.kind = Kind::Triangle;
                for (int corner = 0; corner < 3; ++corner)
                    triangle.triangle[corner] =
                        mesh.vertices[face[corner]];
                repair_output.push_back({std::move(triangle), {face_id}});
            }
        if (repair_output.empty()) break;
        CoplanarCanonicalizationStats repair_coplanar_stats;
        repair_output = canonicalizeCoplanarPrimitiveUnion(
            mesh, std::move(repair_output), final_tolerance,
            repair_coplanar_stats);
        repair_merged_groups += repair_coplanar_stats.groups;
        // Repair the union of failed coplanar source faces, subtracting the
        // already retained exact-coplanar surface before emitting the residual.
        // This is still an exact surface fallback, never a fitted box shell,
        // but it avoids stacking every complete source triangle over a mostly
        // covered proxy patch.

        std::uint64_t group_offset = 0;
        for (const auto& item : output)
            group_offset = std::max(group_offset, item.enclosure_group);
        for (auto& item : repair_output)
            if (item.enclosure_group != 0)
                item.enclosure_group += group_offset;

        repair_output_primitives += repair_output.size();
        repair_output_triangles += triangulatedFaceCount(repair_output);
        for (const auto& item : repair_output)
        {
            coverage_certificate_primitives.push_back(item);
            // Approximate repair enclosures are part of stage-3 error. Exact
            // source triangles have zero error to the original mesh but may be
            // absent from the intentionally hole-filled phase-1 reference.
            if (item.enclosure_group != 0)
                simplification_error_primitives.push_back(item);
        }
        output.insert(output.end(),
            std::make_move_iterator(repair_output.begin()),
            std::make_move_iterator(repair_output.end()));
        // A union replacement can land in the neighboring quantized plane
        // bucket of a component already visited in this pass. Iterate the
        // canonicalizer to a fixed point so those numerically equivalent
        // groups cannot survive as a second overlapping layer.
        for (std::size_t pass = 0; pass < 8; ++pass)
        {
            CoplanarCanonicalizationStats post_repair_coplanar_stats;
            output = canonicalizeCoplanarPrimitiveUnion(
                mesh, std::move(output), final_tolerance,
                post_repair_coplanar_stats);
            repair_merged_groups += post_repair_coplanar_stats.groups;
            if (post_repair_coplanar_stats.groups == 0) break;
        }
        promoteToSemanticPrimitives(output);
        requireSurfaceCandidates(output, "coverage repair");
        coverage = auditFinalConservativeCoverage(
            mesh, output, final_tolerance, &coverage_certificate_primitives,
            &certified_closed_extrusions, &excluded_faces);
    }
    std::vector<std::uint32_t> repair_face_ids(
        repaired_face_set.begin(), repaired_face_set.end());
    std::sort(repair_face_ids.begin(), repair_face_ids.end());
    {
        std::ofstream profile(output_directory /
                              "coverage_audit_pre_repair.json");
        profile << "{\"output_primitives\":" << pre_repair_primitives
                << ",\"output_triangles\":" << pre_repair_triangles
                << ",\"unassigned_source_faces\":"
                << pre_repair_unassigned
                << ",\"failed_source_faces\":" << pre_repair_failed
                << ",\"failed_with_planar_owner\":"
                << pre_repair_failed_with_planar_owner
                << ",\"failed_with_certificate_planar_owner\":"
                << pre_repair_failed_with_certificate_planar_owner
                << ",\"failed_with_certificate_band_owner\":"
                << pre_repair_failed_with_certificate_band_owner
                << ",\"repair_face_count\":" << repair_face_ids.size()
                << ",\"repair_merged_groups\":" << repair_merged_groups
                << ",\"repair_output_primitives\":"
                << repair_output_primitives
                << ",\"repair_output_triangles\":"
                << repair_output_triangles
                << ",\"repair_face_ids\":[";
        for (std::size_t index = 0; index < repair_face_ids.size(); ++index)
        {
            if (index != 0) profile << ',';
            profile << repair_face_ids[index];
        }
        profile << "]}\n";
    }
    const PrimitiveMesh simplification_error_proxy =
        triangulateOutputPrimitives(simplification_error_primitives);
    std::vector<OutputPrimitive>().swap(simplification_error_primitives);
    markStage("coverage_audit_final");
    stats.coverage_assigned_source_faces = coverage.assigned_source_faces;
    stats.coverage_enclosure_source_faces = coverage.enclosure_source_faces;
    stats.coverage_planar_source_faces = coverage.planar_source_faces;
    stats.coverage_unassigned_source_faces = coverage.unassigned_source_faces;
    stats.coverage_failed_source_faces = coverage.failed_source_faces;
    stats.containment_validation_passed =
        coverage.unassigned_source_faces == 0 &&
        coverage.failed_source_faces == 0;
    if (!stats.containment_validation_passed)
    {
        writeCoverageFailureDiagnostics(output_directory, mesh, output, coverage);
        throw std::runtime_error(
            "staged surface pipeline failed conservative coverage audit");
    }

    const auto error_audit_started = std::chrono::steady_clock::now();
    const SourceTriangleBvh filled_surface_reference(filled_surface_mesh);
    const FinalOpenErrorAudit open_error = measureFilledSurfaceDistance(
        filled_surface_reference, simplification_error_proxy,
        std::max(diagonal / 192.0, 1.0e-30));
    markStage("open_error_audit");
    stats.open_error_audit_seconds = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - error_audit_started).count();
    stats.open_error_distance_sample_count = open_error.distance_sample_count;
    stats.open_mean_distance = open_error.mean_distance;
    stats.open_max_distance = open_error.maximum_distance;
    stats.open_mean_distance_ratio = open_error.mean_distance /
        std::max(diagonal, 1.0e-30);
    stats.open_max_distance_ratio = open_error.maximum_distance /
        std::max(diagonal, 1.0e-30);
    stats.open_max_proxy_point = {{open_error.maximum_proxy_point.x(),
                                   open_error.maximum_proxy_point.y(),
                                   open_error.maximum_proxy_point.z()}};
    stats.open_max_source_point = {{open_error.maximum_source_point.x(),
                                    open_error.maximum_source_point.y(),
                                    open_error.maximum_source_point.z()}};
    if (open_error.maximum_distance >
        maximum_open_error_distance + final_tolerance)
    {
        writeOpenErrorVisualization(
            output_directory / "open_error.json", open_error);
        throw std::runtime_error(
            "staged surface pipeline exceeded maximum open error distance");
    }
    std::filesystem::create_directories(output_directory);
    requireSurfaceCandidates(output, "final export");
    writeSourceObj(output_directory / "source.obj", mesh);
    writeRegionsObj(
        output_directory / "regions.obj", mesh, output, excluded_faces);
    stats.primitive_count = output.size();
    for (const auto& item : output)
    {
        switch (item.primitive.kind)
        {
        case Kind::Polygon: ++stats.polygon_count; break;
        case Kind::Disk: ++stats.disk_count; break;
        case Kind::Annulus: ++stats.annulus_count; break;
        case Kind::CylindricalBand: ++stats.cylindrical_band_count; break;
        case Kind::ConicalBand: ++stats.conical_band_count; break;
        case Kind::Rectangle:
        case Kind::Triangle: throw std::logic_error("internal planar primitive was not promoted");
        }
    }
    writeSemanticPrimitiveObj(output_directory / "primitives.obj", output);
    stats.proxy_triangles = writeTriangulatedObj(output_directory / "proxy.obj", output);
    writeOpenErrorVisualization(output_directory / "open_error.json", open_error);
    stats.analysis_seconds = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - started).count();
    writeMetadata(output_directory, input_obj, output, stats, effective_options);
    return stats;
}

} // namespace pqss_proxy_mesh
