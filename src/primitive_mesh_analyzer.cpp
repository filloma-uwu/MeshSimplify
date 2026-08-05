#include "pqss_proxy_mesh/primitive_mesh_analyzer.hpp"

#include "clipper2/clipper.h"
#include "clipper2/clipper.triangulation.h"

#include <algorithm>
#include <array>
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
    Frustum,
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

struct VolumeFit
{
    BoxFit box;
    Primitive frustum;
    bool use_frustum = false;
    double volume = 0.0;
    double surface_area = 0.0;
};

VolumeFit boxVolumeFit(const BoxFit& box)
{
    VolumeFit fit;
    fit.box = box;
    fit.volume = box.volume;
    const Vec3 extent = box.half_size * 2.0;
    fit.surface_area = 2.0 * (extent.x() * extent.y() + extent.y() * extent.z() +
                              extent.z() * extent.x());
    return fit;
}

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

struct Node
{
    std::vector<std::uint32_t> faces;
    std::unordered_set<std::uint32_t> neighbors;
    VolumeFit fit;
    Vec3 lower = Vec3::Constant(std::numeric_limits<double>::infinity());
    Vec3 upper = Vec3::Constant(-std::numeric_limits<double>::infinity());
    double source_area = 0.0;
    double activation_threshold = 0.0;
    std::uint32_t parent = std::numeric_limits<std::uint32_t>::max();
    std::uint64_t version = 0;
    bool active = true;
    bool recognized_protrusion = false;
    bool recognized_closed_box = false;
    int covered_box_face_axis = -1;
    double covered_box_face_sign = 0.0;
};

struct Candidate
{
    double added_ratio = 0.0;
    double volume = 0.0;
    std::uint32_t first = 0;
    std::uint32_t second = 0;
    std::uint64_t first_version = 0;
    std::uint64_t second_version = 0;
    VolumeFit fit;
};

struct CandidateGreater
{
    bool operator()(const Candidate& a, const Candidate& b) const
    {
        if (a.added_ratio != b.added_ratio) return a.added_ratio > b.added_ratio;
        return a.volume > b.volume;
    }
};

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
    if (primitive.kind == Kind::Frustum)
    {
        if (primitive.height <= dimensional_tolerance)
            return std::numbers::pi * primitive.base_radius * primitive.base_radius;
        const double slant = std::hypot(
            primitive.height, primitive.top_radius - primitive.base_radius);
        return std::numbers::pi * (primitive.base_radius + primitive.top_radius) * slant +
               std::numbers::pi * (primitive.base_radius * primitive.base_radius +
                                   primitive.top_radius * primitive.top_radius);
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
    // Expand both end radii by the largest deficit so the fitted frustum still
    // contains every input vertex, including intermediate axial rings.
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
    result.kind = Kind::Frustum;
    result.axes = frame;
    result.height = height;
    result.center = mean + frame.col(1) * radial_center.x() +
        frame.col(2) * radial_center.y() + axis.normalized() * lower;
    result.base_radius = base_radius;
    result.top_radius = top_radius;
    result.segments = segments;
    result.volume = std::numbers::pi * height / 3.0 *
        (base_radius * base_radius + base_radius * top_radius + top_radius * top_radius);
    return result;
}

std::size_t faceComponentCount(const Mesh& mesh,
                               const std::vector<std::uint32_t>& faces);

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

double boxSurfaceArea(const BoxFit& box)
{
    const Vec3 extent = box.half_size * 2.0;
    return 2.0 * (extent.x() * extent.y() + extent.y() * extent.z() +
                  extent.z() * extent.x());
}

VolumeFit fitBestVolume(const Mesh& mesh,
                        const std::vector<std::uint32_t>& faces,
                        const PrimitiveMeshAnalysisOptions& options,
                        const double analytic_tolerance = 0.0,
                        const bool allow_conservative_enclosing_frustum = false)
{
    const auto vertices = uniqueVertices(mesh, faces);
    VolumeFit best;
    best.box = fitBox(mesh, vertices);
    best.volume = best.box.volume;
    best.surface_area = boxSurfaceArea(best.box);
    if (options.allow_frustum)
    {
        for (const Mat3& frame : candidateFrames(mesh, vertices))
        {
            for (int axis = 0; axis < 3; ++axis)
            {
                Primitive candidate = fitConeOnAxis(mesh, vertices, static_cast<Vec3>(frame.col(axis)),
                                                      options.frustum_segments,
                                                      analytic_tolerance);
                if (candidate.volume < best.volume &&
                    (allow_conservative_enclosing_frustum ||
                     coneFitIsCertified(mesh, faces, candidate,
                                        options.circle_radial_tolerance,
                                        analytic_tolerance)))
                {
                    best.frustum = candidate;
                    best.use_frustum = true;
                    best.volume = candidate.volume;
                    best.surface_area = primitiveSurfaceArea(candidate, 1.0e-12);
                }
            }
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
    std::unordered_map<std::uint64_t, std::vector<std::uint32_t>> edge_faces;
    edge_faces.reserve(mesh.faces.size() * 3);
    std::vector<Vec3> normals(mesh.faces.size());
    for (std::uint32_t id = 0; id < mesh.faces.size(); ++id)
    {
        if (included_faces && !(*included_faces)[id]) continue;
        const Face& face = mesh.faces[id];
        for (int edge = 0; edge < 3; ++edge)
            edge_faces[edgeKey(face[edge], face[(edge + 1) % 3])].push_back(id);
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
            result.push_back(
                {box, std::move(faces), group, covered_axis, covered_sign});
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

    // Spatially merged planar islands are separate objects, not holes in one patch.
    if (faceComponentCount(mesh, faces) > 1)
        return exactTrianglePrimitives(mesh, faces);

    auto loops = boundaryLoops(mesh, faces);
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
    const std::vector<Vec2> hull = simplifyPolygon(convexHull(polygon, planar_tolerance),
                                                    planar_tolerance);
    const std::size_t hole_count = loops.size() > 1 ? loops.size() - 1 : 0;
    const double effective_source_area = hole_count > 0 ? largest_area : source_area;
    if (options.allow_polygon && hull.size() == 4)
    {
        const Vec2 first = hull[1] - hull[0];
        const Vec2 second = hull[3] - hull[0];
        const double rectangle_area = std::abs(signedArea(hull));
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
            const Vec2 center = (hull[0] + hull[2]) * 0.5;
            rectangle.center = origin + frame.col(0) * center.x() + frame.col(1) * center.y();
            filled_holes += hole_count;
            return {{rectangle, faces}};
        }
    }

    if (options.allow_frustum && polygon.size() >= 8)
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
        const double circle_area = std::numbers::pi * mean * mean;
        const double circle_excess_ratio =
            std::max(circle_area - effective_source_area, 0.0) /
            std::max(model_surface_area, 1.0e-30);
        if (maximum_radial_error <= options.circle_radial_tolerance * mean &&
            ringLooksCircular(centered_polygon, options.circle_radial_tolerance) &&
            circularity >= 0.9 &&
            circle_excess_ratio <= allowed_excess_ratio)
        {
            Primitive circle;
            circle.kind = Kind::Frustum;
            circle.center = origin + frame.col(0) * center.x() + frame.col(1) * center.y();
            circle.axes.col(0) = normal;
            circle.axes.col(1) = frame.col(0);
            circle.axes.col(2) = frame.col(1);
            circle.base_radius = mean;
            circle.top_radius = mean;
            circle.segments = options.frustum_segments;
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
                auto connected = faceComponentsFromList(mesh, plane_group);
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
    if (!options.allow_frustum || faces.size() < 6 ||
        faceComponentCount(mesh, faces) != 1) return {};

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

    auto loops3 = boundaryLoops(mesh, faces);
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
    surface.base_radius = outer.maximum_radius;
    surface.top_radius = surface.base_radius;
    surface.segments = options.frustum_segments;
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
    filled_holes += candidate_filled_holes;
    filled_volume += candidate_accepted_volume;
    return {{std::move(surface), faces}};
}

void appendConicalSurfaceAssembly(std::vector<OutputPrimitive>& output,
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

void appendEnclosingConicalSurfaceAssembly(
    std::vector<OutputPrimitive>& output,
    const Primitive& fit,
    const std::vector<std::uint32_t>& faces,
    const std::uint64_t enclosure_group)
{
    if (fit.height <= 0.0 || fit.segments < 3) return;
    Primitive band = fit;
    band.kind = std::abs(fit.base_radius - fit.top_radius) <=
            1.0e-8 * std::max({fit.base_radius, fit.top_radius, fit.height})
        ? Kind::CylindricalBand : Kind::ConicalBand;
    OutputPrimitive band_item{std::move(band), faces};
    band_item.enclosure_group = enclosure_group;
    output.push_back(std::move(band_item));

    const auto append_disk = [&](const Vec3& center, const double radius)
    {
        if (radius <= 0.0) return;
        Primitive disk;
        disk.kind = Kind::Disk;
        disk.center = center;
        disk.axes = fit.axes;
        disk.base_radius = radius;
        disk.top_radius = radius;
        disk.segments = fit.segments;
        OutputPrimitive disk_item{std::move(disk), faces};
        disk_item.enclosure_group = enclosure_group;
        output.push_back(std::move(disk_item));
    };
    append_disk(fit.center, fit.base_radius);
    append_disk(fit.center + fit.axes.col(0) * fit.height, fit.top_radius);
}

PrimitiveMesh triangulatePrimitive(const Primitive& primitive)
{
    PrimitiveMesh result;
    if (primitive.kind == Kind::Polygon)
    {
        if (primitive.polygon.size() < 3)
            throw std::runtime_error("polygon primitive has fewer than three vertices");
        result.vertices = primitive.polygon;
        const Vec3 normal = (primitive.polygon[1] - primitive.polygon[0])
                                .cross(primitive.polygon[2] - primitive.polygon[0]);
        if (normal.norm() <= 1.0e-30)
        {
            std::ostringstream message;
            message << "polygon primitive has a degenerate plane; vertices=";
            for (const Vec3& vertex : primitive.polygon)
                message << " (" << vertex.x() << ',' << vertex.y() << ',' << vertex.z() << ')';
            throw std::runtime_error(message.str());
        }
        Mat3 basis = orthonormalFrame(normal.normalized());
        Mat3 frame;
        frame.col(0) = static_cast<Vec3>(basis.col(1));
        frame.col(1) = static_cast<Vec3>(basis.col(2));
        frame.col(2) = normal.normalized();
        std::vector<Vec2> boundary;
        boundary.reserve(primitive.polygon.size());
        for (const Vec3& vertex : primitive.polygon)
        {
            const Vec3 local = frame.transposeMultiply(vertex - primitive.polygon.front());
            boundary.emplace_back(local.x(), local.y());
        }
        const auto indices = triangulatePolygon(std::move(boundary));
        if (indices.size() + 2 != primitive.polygon.size())
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
    if (primitive.kind == Kind::Disk ||
        (primitive.kind == Kind::Frustum && primitive.height == 0.0))
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
    const Vec3 top_center = primitive.center + primitive.axes.col(0) * primitive.height;
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
            result.faces.push_back(
                {2 * index, 2 * next, 2 * next + 1});
            result.faces.push_back(
                {2 * index, 2 * next + 1, 2 * index + 1});
        }
        return result;
    }
    if (!band)
    {
        result.vertices.push_back(primitive.center);
        result.vertices.push_back(top_center);
    }
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
        const std::uint32_t vertex_offset = band ? 0 : 2;
        const std::uint32_t base = vertex_offset + 2 * index;
        const std::uint32_t top = base + 1;
        const std::uint32_t next_base = vertex_offset + 2 * next;
        const std::uint32_t next_top = next_base + 1;
        result.faces.push_back({base, next_base, next_top});
        result.faces.push_back({base, next_top, top});
        if (!band)
        {
            result.faces.push_back({0, next_base, base});
            if (primitive.top_radius > 0.0) result.faces.push_back({1, top, next_top});
        }
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
    if (primitive.kind == Kind::Frustum)
    {
        const Vec3 local = primitive.axes.transposeMultiply(point - primitive.center);
        if (primitive.height <= tolerance)
            return std::abs(local.x()) <= tolerance &&
                   std::hypot(local.y(), local.z()) <= primitive.base_radius + tolerance;
        if (local.x() < -tolerance || local.x() > primitive.height + tolerance) return false;
        const double parameter = std::clamp(local.x() / primitive.height, 0.0, 1.0);
        const double radius = primitive.base_radius +
                              parameter * (primitive.top_radius - primitive.base_radius);
        return std::hypot(local.y(), local.z()) <= radius + tolerance;
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

bool trimmedAnalyticBandContainsTriangle(
    const Primitive& primitive, const std::array<Vec3, 3>& triangle,
    const double tolerance)
{
    if ((primitive.kind != Kind::CylindricalBand &&
         primitive.kind != Kind::ConicalBand) ||
        primitive.band_axial_ranges.size() < 3 || primitive.height <= 0.0)
        return false;
    const std::size_t segments = primitive.band_axial_ranges.size();
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
        const double fraction = sample - std::floor(sample);
        const double lower =
            primitive.band_axial_ranges[first][0] * (1.0 - fraction) +
            primitive.band_axial_ranges[second][0] * fraction;
        const double upper =
            primitive.band_axial_ranges[first][1] * (1.0 - fraction) +
            primitive.band_axial_ranges[second][1] * fraction;
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
    if (primitive.kind == Kind::Frustum && primitive.height <= tolerance)
    {
        normal = primitive.axes.col(0);
        return true;
    }
    return false;
}

Vec3 planarPoint(const Primitive& primitive)
{
    if (primitive.kind == Kind::Polygon) return primitive.polygon.front();
    return primitive.kind == Kind::Triangle ? primitive.triangle[0] : primitive.center;
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

double maximumChargeableOpenSurfaceDistance(
    const Mesh& source_mesh,
    const std::vector<OutputPrimitive>& output,
    double sample_spacing,
    double tolerance);

std::vector<OutputPrimitive> mergeLocalCoplanarPrimitives(
    const Mesh& source_mesh,
    std::vector<OutputPrimitive> primitives,
    const double tolerance,
    const double maximum_open_error_distance,
    const double error_sample_spacing,
    std::size_t& merged_count)
{
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
        Primitive rectangle;
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

    const auto evaluate = [&](const std::size_t first,
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
            return std::nullopt;
        const std::uint64_t first_enclosure =
            items[first].output.enclosure_group;
        const std::uint64_t second_enclosure =
            items[second].output.enclosure_group;
        if (first_enclosure != second_enclosure &&
            (first_enclosure != 0 || second_enclosure != 0))
            return std::nullopt;
        const Primitive& first_primitive = items[first].output.primitive;
        const Primitive& second_primitive = items[second].output.primitive;
        Vec3 normal;
        Vec3 second_normal;
        if (first_primitive.kind == Kind::Frustum ||
            second_primitive.kind == Kind::Frustum ||
            !planarNormal(first_primitive, tolerance, normal) ||
            !planarNormal(second_primitive, tolerance, second_normal) ||
            std::abs(normal.dot(second_normal)) < 1.0 - 1.0e-8)
            return std::nullopt;
        const Vec3 origin = planarPoint(first_primitive);
        if (std::abs((planarPoint(second_primitive) - origin).dot(normal)) > tolerance)
            return std::nullopt;
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
        auto hull = simplifyPolygon(convexHull(std::move(points), tolerance), tolerance);
        // This pass deliberately performs only 2-to-1 rectangle merges.
        // General n-gon replacement belongs in the region classifier.
        if (hull.size() != 4) return std::nullopt;
        const Vec2 first_edge = hull[1] - hull[0];
        const Vec2 second_edge = hull[3] - hull[0];
        if (first_edge.norm() <= tolerance || second_edge.norm() <= tolerance ||
            std::abs(first_edge.dot(second_edge)) >
                1.0e-8 * first_edge.norm() * second_edge.norm())
            return std::nullopt;
        Primitive rectangle;
        rectangle.kind = Kind::Rectangle;
        rectangle.axes.col(0) = frame.col(0) * (first_edge.x() / first_edge.norm()) +
                                frame.col(1) * (first_edge.y() / first_edge.norm());
        rectangle.axes.col(1) = frame.col(0) * (second_edge.x() / second_edge.norm()) +
                                frame.col(1) * (second_edge.y() / second_edge.norm());
        rectangle.axes.col(2) = normal;
        rectangle.half_size = {0.5 * first_edge.norm(), 0.5 * second_edge.norm(), 0.0};
        const Vec2 center = (hull[0] + hull[2]) * 0.5;
        rectangle.center = origin + frame.col(0) * center.x() +
                           frame.col(1) * center.y();
        OutputPrimitive merged;
        merged.primitive = rectangle;
        merged.source_faces = items[first].output.source_faces;
        merged.source_faces.insert(
            merged.source_faces.end(),
            items[second].output.source_faces.begin(),
            items[second].output.source_faces.end());
        const double maximum_distance = maximumChargeableOpenSurfaceDistance(
            source_mesh, {merged}, error_sample_spacing, tolerance);
        if (maximum_distance > maximum_open_error_distance + tolerance)
            return std::nullopt;
        return MergeFit{maximum_distance, std::move(rectangle)};
    };

    // Coplanar planes are independent merge domains.  Build those domains with
    // a linear-memory disjoint set, then allocate a candidate heap for only one
    // plane at a time.  This preserves the exact pairwise acceptance/order
    // within every plane while avoiding a global O(n^2) resident candidate set.
    std::vector<std::size_t> parent(items.size());
    std::iota(parent.begin(), parent.end(), 0);
    const auto rootOf = [&](std::size_t index)
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
    std::vector<Vec3> plane_normals(items.size(), Vec3::Zero());
    std::vector<Vec3> plane_points(items.size(), Vec3::Zero());
    for (std::size_t index = 0; index < items.size(); ++index)
    {
        (void)planarNormal(items[index].output.primitive, tolerance,
                           plane_normals[index]);
        plane_points[index] = planarPoint(items[index].output.primitive);
    }
    for (std::size_t first = 0; first < items.size(); ++first)
    {
        if (plane_normals[first].norm() == 0.0) continue;
        for (std::size_t second = first + 1; second < items.size(); ++second)
        {
            if (plane_normals[second].norm() == 0.0 ||
                std::abs(plane_normals[first].dot(plane_normals[second])) <
                    1.0 - 1.0e-8 ||
                std::abs((plane_points[second] - plane_points[first]).dot(
                    plane_normals[first])) > tolerance)
                continue;
            const std::size_t first_root = rootOf(first);
            const std::size_t second_root = rootOf(second);
            if (first_root != second_root) parent[second_root] = first_root;
        }
    }
    std::map<std::size_t, std::vector<std::size_t>> plane_groups;
    for (std::size_t index = 0; index < items.size(); ++index)
        plane_groups[rootOf(index)].push_back(index);

    for (const auto& [root, group] : plane_groups)
    {
        (void)root;
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
        for (std::size_t first = 0; first < group.size(); ++first)
            for (std::size_t second = first + 1; second < group.size(); ++second)
            {
                const std::size_t first_id = group[first];
                const std::size_t second_id = group[second];
                const auto fit = evaluate(first_id, second_id);
                if (!fit) continue;
                (void)consider(first_id, second_id, *fit);
                (void)consider(second_id, first_id, *fit);
            }
        for (const std::size_t owner : group) pushCurrent(owner);

        const auto recompute = [&](const std::size_t owner)
        {
            best_partner[owner] = no_partner;
            best_excess[owner] = std::numeric_limits<double>::infinity();
            ++best_generation[owner];
            if (!items[owner].active) return;
            for (const std::size_t partner : group)
            {
                if (partner == owner || !items[partner].active) continue;
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
            replacement.primitive = fit->rectangle;
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
            ++items[candidate.first].version;
            items[candidate.second].active = false;
            ++items[candidate.second].version;
            ++merged_count;
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
        // A connected Boolean-union component is one semantic polygon. Its
        // triangulation belongs to the later collision-mesh stage.
        if (!covered.empty())
        {
            const bool has_holes = std::any_of(covered.begin(), covered.end(),
                [](const auto& path) { return !Clipper2Lib::IsPositive(path); });
            if (has_holes)
            {
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

        for (const auto face_id : source_faces)
        {
            const Face& face = source_mesh.faces[face_id];
            const Vec3 centroid3 = (source_mesh.vertices[face[0]] +
                                    source_mesh.vertices[face[1]] +
                                    source_mesh.vertices[face[2]]) / 3.0;
            std::size_t selected = replacement_begin;
            for (std::size_t index = replacement_begin; index < result.size(); ++index)
                if (containsPoint(result[index].primitive, centroid3, tolerance * 8.0))
                {
                    selected = index;
                    break;
                }
            result[selected].source_faces.push_back(face_id);
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
        if (removed[candidate]) continue;
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
    for (const OutputPrimitive& certificate : certificate_primitives)
        for (const auto face : certificate.source_faces)
            if (face < responsibility_certified.size())
                responsibility_certified[face] = 1;
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
        if (covered)
        {
            ++audit.enclosure_source_faces;
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

        const auto coveredByPlanarPrimitives = [&](const std::vector<std::size_t>& candidates)
        {
            Clipper2Lib::PathsD cover_paths;
            for (const std::size_t index : candidates)
            {
                Vec3 candidate_normal;
                if (!planarNormal(primitives[index].primitive, tolerance,
                                  candidate_normal) ||
                    std::abs(candidate_normal.dot(normal)) < 1.0 - 1.0e-8 ||
                    std::abs((planarPoint(primitives[index].primitive) -
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
        }
    }
    return audit;
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

double maximumChargeableOpenSurfaceDistance(
    const Mesh& source_mesh,
    const std::vector<OutputPrimitive>& output,
    const double sample_spacing,
    const double tolerance)
{
    Vec3 lower = Vec3::Constant(std::numeric_limits<double>::infinity());
    Vec3 upper = Vec3::Constant(-std::numeric_limits<double>::infinity());
    for (const Vec3& vertex : source_mesh.vertices)
    {
        lower = lower.cwiseMin(vertex);
        upper = upper.cwiseMax(vertex);
    }
    const SourceOpenSurfaceReference source(
        source_mesh, (upper - lower).norm());
    return measureChargeableOpenSurfaceDistance(
        source, source_mesh, output, sample_spacing, tolerance)
        .maximum_distance;
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

struct UnifiedProxyChoiceNode
{
    std::vector<OutputPrimitive> direct_output;
    std::shared_ptr<const UnifiedProxyChoiceNode> left;
    std::shared_ptr<const UnifiedProxyChoiceNode> right;
    std::shared_ptr<const UnifiedProxyChoiceNode> refinement;
};

struct UnifiedProxyChoice
{
    std::shared_ptr<const UnifiedProxyChoiceNode> node;
    // Combined candidates remain lazy until they survive frontier pruning.
    // Materializing every temporary Cartesian-product entry causes the Windows
    // heap to retain gigabytes of discarded choice nodes on large CAD meshes.
    std::shared_ptr<const UnifiedProxyChoiceNode> pending_left;
    std::shared_ptr<const UnifiedProxyChoiceNode> pending_right;
    double workload = std::numeric_limits<double>::infinity();
    double open_mean_distance_ratio = 0.0;
    double open_maximum_distance = 0.0;
    double open_maximum_distance_ratio = 0.0;
    Vec3 open_maximum_proxy_point = Vec3::Zero();
    Vec3 open_maximum_source_point = Vec3::Zero();
    std::size_t final_enclosure_overlap_pairs = 0;
    double final_enclosure_overlap_score = 0.0;
    bool final_coverage_passed = false;
    std::size_t final_unassigned_source_faces = 0;
    std::size_t final_failed_source_faces = 0;
    std::size_t analytic_responsibility_faces = 0;
    double analytic_responsibility_area = 0.0;
    std::size_t enclosure_responsibility_faces = 0;
    double enclosure_responsibility_area = 0.0;
    bool preserve_pre_removal_enclosures = true;
    // An enclosing candidate may temporarily exceed the global error limit
    // while still owning a certified, finer replacement.  Keep that state in
    // the hierarchy so the final audit can expand only the branch containing
    // the measured maximum-error point.
    bool can_refine = false;
    std::vector<std::uint32_t> coverage_patch_faces;
};

void flattenUnifiedChoice(const std::shared_ptr<const UnifiedProxyChoiceNode>& node,
                          std::vector<OutputPrimitive>& output)
{
    if (!node) return;
    if (node->left || node->right)
    {
        flattenUnifiedChoice(node->left, output);
        flattenUnifiedChoice(node->right, output);
        return;
    }
    output.insert(output.end(), node->direct_output.begin(), node->direct_output.end());
}

bool unifiedChoiceSurfaceContainsPoint(
    const std::shared_ptr<const UnifiedProxyChoiceNode>& node,
    const Vec3& point, const double tolerance)
{
    if (!node) return false;
    const double tolerance_squared = tolerance * tolerance;
    return std::any_of(
        node->direct_output.begin(), node->direct_output.end(),
        [&](const OutputPrimitive& item)
        {
            const PrimitiveMesh mesh = triangulatePrimitive(item.primitive);
            return std::any_of(mesh.faces.begin(), mesh.faces.end(),
                [&](const Face& face)
                {
                    return (closestPointOnTriangle(
                        point, mesh.vertices[face[0]],
                        mesh.vertices[face[1]],
                        mesh.vertices[face[2]]) - point).squaredNorm() <=
                        tolerance_squared;
                });
        });
}

double unifiedChoiceSurfaceDistanceSquared(
    const std::shared_ptr<const UnifiedProxyChoiceNode>& node,
    const Vec3& point)
{
    if (!node) return std::numeric_limits<double>::infinity();
    double distance_squared = std::numeric_limits<double>::infinity();
    for (const OutputPrimitive& item : node->direct_output)
    {
        const PrimitiveMesh mesh = triangulatePrimitive(item.primitive);
        for (const Face& face : mesh.faces)
            distance_squared = std::min(
                distance_squared,
                (closestPointOnTriangle(
                    point, mesh.vertices[face[0]], mesh.vertices[face[1]],
                    mesh.vertices[face[2]]) - point).squaredNorm());
    }
    return distance_squared;
}

const UnifiedProxyChoiceNode* findNearestRefinableUnifiedChoice(
    const std::shared_ptr<const UnifiedProxyChoiceNode>& node,
    const Vec3& point, double& best_distance_squared)
{
    if (!node) return nullptr;
    const UnifiedProxyChoiceNode* best = nullptr;
    if (node->refinement)
    {
        const double distance_squared =
            unifiedChoiceSurfaceDistanceSquared(node, point);
        if (distance_squared < best_distance_squared)
        {
            best_distance_squared = distance_squared;
            best = node.get();
        }
    }
    for (const auto& child : {node->left, node->right})
    {
        const UnifiedProxyChoiceNode* candidate =
            findNearestRefinableUnifiedChoice(
                child, point, best_distance_squared);
        if (candidate != nullptr) best = candidate;
    }
    return best;
}

std::shared_ptr<const UnifiedProxyChoiceNode> refineUnifiedChoiceNode(
    const std::shared_ptr<const UnifiedProxyChoiceNode>& node,
    const UnifiedProxyChoiceNode* target, bool& changed)
{
    if (!node || changed) return node;
    if (node.get() == target && node->refinement)
    {
        changed = true;
        return node->refinement;
    }
    if (!node->left && !node->right) return node;
    const auto left = refineUnifiedChoiceNode(node->left, target, changed);
    const auto right = refineUnifiedChoiceNode(node->right, target, changed);
    if (left == node->left && right == node->right) return node;
    auto refined = std::make_shared<UnifiedProxyChoiceNode>(*node);
    refined->left = left;
    refined->right = right;
    return refined;
}

std::shared_ptr<const UnifiedProxyChoiceNode> refineUnifiedChoiceAtPoint(
    const std::shared_ptr<const UnifiedProxyChoiceNode>& node,
    const Vec3& point, const double tolerance, bool& changed)
{
    if (!node || changed) return node;
    if (node->refinement &&
        unifiedChoiceSurfaceContainsPoint(node, point, tolerance))
    {
        changed = true;
        return node->refinement;
    }
    if (!node->left && !node->right) return node;
    const auto left = refineUnifiedChoiceAtPoint(
        node->left, point, tolerance, changed);
    const auto right = refineUnifiedChoiceAtPoint(
        node->right, point, tolerance, changed);
    if (left == node->left && right == node->right) return node;
    auto refined = std::make_shared<UnifiedProxyChoiceNode>(*node);
    refined->left = left;
    refined->right = right;
    return refined;
}

std::size_t triangulatedFaceCount(const std::vector<OutputPrimitive>& output)
{
    std::size_t count = 0;
    for (const auto& item : output)
        count += triangulatePrimitive(item.primitive).faces.size();
    return count;
}

std::vector<OutputPrimitive> regularizeShallowNearRectangles(
    std::vector<OutputPrimitive> primitives,
    double model_surface_area,
    double maximum_added_area_ratio,
    double maximum_open_error_distance,
    double tolerance);

std::vector<OutputPrimitive> optimizeUnifiedProxyCandidates(
    const Mesh& mesh,
    const PrimitiveMeshAnalysisOptions& options,
    const double model_diagonal,
    const double model_volume,
    const double model_surface_area,
    PrimitiveMeshAnalysisStats& stats,
    const std::filesystem::path& performance_path,
    const std::vector<std::uint32_t>* excluded_redundant_faces = nullptr,
    std::vector<OutputPrimitive>* selected_coverage_certificates = nullptr)
{
    const double maximum_open_error_distance =
        options.maximum_open_error_distance >= 0.0
            ? options.maximum_open_error_distance
            : model_diagonal * 0.08;
    stats.maximum_open_error_distance_limit = maximum_open_error_distance;
    const auto diagnostic_started = std::chrono::steady_clock::now();
    const auto diagnostic = [&](const char* stage, const std::size_t first = 0,
                                const std::size_t second = 0)
    {
        const double seconds = std::chrono::duration<double>(
            std::chrono::steady_clock::now() - diagnostic_started).count();
        std::ofstream output(performance_path);
        output << std::setprecision(10)
               << "stage=" << stage << '\n'
               << "elapsed_seconds=" << seconds << '\n'
               << "stage_first=" << first << '\n'
               << "stage_second=" << second << '\n'
               << "approximate_planar_regions="
               << stats.unified_approximate_planar_regions << '\n'
               << "exact_coplanar_regions="
               << stats.unified_exact_coplanar_regions << '\n'
               << "hierarchy_regions=" << stats.unified_hierarchy_regions << '\n'
               << "connected_components="
               << stats.unified_connected_components << '\n'
               << "lateral_sweep_axes="
               << stats.unified_lateral_sweep_axes << '\n'
               << "lateral_sweep_patches="
               << stats.unified_lateral_sweep_patches << '\n'
               << "lateral_sweep_certified="
               << stats.unified_lateral_sweep_certified << '\n'
               << "lateral_sweep_candidates="
               << stats.unified_lateral_sweep_candidates << '\n'
               << "lateral_sweep_max_faces="
               << stats.unified_lateral_sweep_max_faces << '\n'
               << "lateral_significant_candidates="
               << stats.unified_lateral_significant_candidates << '\n'
               << "lateral_sweep_max_axial_aspect="
               << stats.unified_lateral_sweep_max_axial_aspect << '\n'
               << "lateral_sweep_max_area_ratio="
               << stats.unified_lateral_sweep_max_area_ratio << '\n'
               << "lateral_reject_geometry="
               << stats.unified_lateral_reject_geometry << '\n'
               << "lateral_reject_end_ring="
               << stats.unified_lateral_reject_end_ring << '\n'
               << "lateral_reject_normal="
               << stats.unified_lateral_reject_normal << '\n'
               << "lateral_reject_radial="
               << stats.unified_lateral_reject_radial << '\n'
               << "lateral_best_end_ring_min_points="
               << stats.unified_lateral_best_end_ring_min_points << '\n'
               << "lateral_best_end_ring_spread_ratio="
               << stats.unified_lateral_best_end_ring_spread_ratio << '\n'
               << "lateral_best_end_ring_gap_ratio="
               << stats.unified_lateral_best_end_ring_gap_ratio << '\n'
               << "lateral_best_end_ring_score="
               << stats.unified_lateral_best_end_ring_score << '\n'
               << "classify_calls=" << stats.unified_classify_calls << '\n'
               << "classify_faces=" << stats.unified_classify_faces << '\n'
               << "classify_seconds=" << stats.unified_classify_seconds << '\n'
               << "frontier_prune_calls="
               << stats.unified_frontier_prune_calls << '\n'
               << "frontier_prune_choices="
               << stats.unified_frontier_prune_choices << '\n'
               << "frontier_prune_max_choices="
               << stats.unified_frontier_prune_max_choices << '\n'
               << "frontier_prune_comparisons="
               << stats.unified_frontier_prune_comparisons << '\n'
               << "frontier_prune_seconds="
               << stats.unified_frontier_prune_seconds << '\n'
               << "planar_region_seconds="
               << stats.unified_planar_region_seconds << '\n'
               << "connected_component_seconds="
               << stats.unified_connected_component_seconds << '\n'
               << "analytic_recognition_seconds="
               << stats.unified_analytic_recognition_seconds << '\n'
               << "hierarchy_seconds=" << stats.unified_hierarchy_seconds << '\n';
    };
    const auto planar_started = std::chrono::steady_clock::now();
    std::vector<bool> included(mesh.faces.size(), true);
    if (excluded_redundant_faces != nullptr)
        for (const auto face : *excluded_redundant_faces)
            if (face < included.size()) included[face] = false;
    const auto approximate_regions = approximatePlanarRegions(
        mesh, included, options, model_diagonal);
    std::vector<std::unordered_set<std::uint32_t>> adjacency;
    const auto exact_regions = coplanarClusters(
        mesh, model_diagonal * std::max(
            options.coplanar_relative_tolerance,
            options.analytic_surface_relative_tolerance),
        adjacency, &included);
    stats.unified_approximate_planar_regions = approximate_regions.size();
    stats.unified_exact_coplanar_regions = exact_regions.size();

    // CAD exporters frequently duplicate vertices and perturb the normals of
    // skinny triangles from one analytic plane.  Using only index-connected,
    // pairwise-coplanar clusters turns such a plane into tens of thousands of
    // hierarchy leaves.  The area-weighted approximate domains recover the
    // complete plane; exact clusters then partition every still-unassigned
    // source face so coverage remains exhaustive.  classifyFinalRegion still
    // performs the final whole-domain plane certificate and falls back to the
    // original triangles if an approximate domain is not planar enough.
    std::vector<std::vector<std::uint32_t>> regions;
    regions.reserve(approximate_regions.size() + exact_regions.size());
    std::vector<bool> assigned(mesh.faces.size(), false);
    for (const auto& region : approximate_regions)
    {
        std::vector<std::uint32_t> unassigned;
        unassigned.reserve(region.size());
        for (const auto face : region)
            if (!assigned[face]) unassigned.push_back(face);
        if (unassigned.empty()) continue;
        for (const auto face : unassigned) assigned[face] = true;
        regions.push_back(std::move(unassigned));
    }
    for (const auto& exact : exact_regions)
    {
        std::vector<std::uint32_t> unassigned;
        unassigned.reserve(exact.size());
        for (const auto face : exact)
            if (!assigned[face]) unassigned.push_back(face);
        if (unassigned.empty()) continue;
        for (auto& component : faceComponentsFromList(mesh, unassigned))
        {
            for (const auto face : component) assigned[face] = true;
            regions.push_back(std::move(component));
        }
    }
    for (std::uint32_t face = 0; face < mesh.faces.size(); ++face)
        if (included[face] && !assigned[face]) regions.push_back({face});
    if (regions.empty()) return {};
    stats.unified_hierarchy_regions = regions.size();
    stats.unified_planar_region_seconds = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - planar_started).count();
    diagnostic("coplanar_regions", regions.size(), mesh.faces.size());

    const double tolerance = std::max(model_diagonal * 1.0e-9, 1.0e-10);
    const auto timedClassifyFinalRegion = [&](
        const std::vector<std::uint32_t>& faces,
        const PrimitiveMeshAnalysisOptions& final_options,
        const double allowed_excess_ratio = 0.0)
    {
        const auto started = std::chrono::steady_clock::now();
        ++stats.unified_classify_calls;
        stats.unified_classify_faces += faces.size();
        if (faces.size() >= 256)
            diagnostic("classify_begin", stats.unified_classify_calls,
                       faces.size());
        std::size_t ignored_holes = 0;
        std::size_t ignored_voids = 0;
        double ignored_void_area = 0.0;
        auto result = classifyFinalRegion(
            mesh, faces, final_options, allowed_excess_ratio, model_surface_area,
            ignored_holes, ignored_voids, ignored_void_area);
        stats.unified_classify_seconds += std::chrono::duration<double>(
            std::chrono::steady_clock::now() - started).count();
        if (faces.size() >= 256)
            diagnostic("classify_end", stats.unified_classify_calls,
                       faces.size());
        return result;
    };
    const SourceOpenSurfaceReference open_surface_reference(
        mesh, model_diagonal);
    const double local_open_error_sample_spacing =
        std::max(model_diagonal / 16.0, 1.0e-30);
    const double final_open_error_sample_spacing =
        std::max(model_diagonal / 192.0, 1.0e-30);
    const auto makeChoice = [&](std::vector<OutputPrimitive> output)
    {
        UnifiedProxyChoice choice;
        choice.workload = static_cast<double>(triangulatedFaceCount(output));
        const FinalOpenErrorAudit open_error = measureChargeableOpenSurfaceDistance(
            open_surface_reference, mesh, output,
            local_open_error_sample_spacing, tolerance);
        choice.open_mean_distance_ratio =
            open_error.distance_area_integral /
            std::max(model_surface_area * model_diagonal, 1.0e-30);
        choice.open_maximum_distance = open_error.maximum_distance;
        choice.open_maximum_distance_ratio =
            open_error.maximum_distance / std::max(model_diagonal, 1.0e-30);
        choice.open_maximum_proxy_point = open_error.maximum_proxy_point;
        choice.open_maximum_source_point = open_error.maximum_source_point;
        if (choice.open_maximum_distance >
            maximum_open_error_distance + tolerance)
            choice.workload = std::numeric_limits<double>::infinity();
        if (std::any_of(output.begin(), output.end(),
                        [](const OutputPrimitive& item)
                        {
                            return item.enclosure_group != 0;
                        }))
        {
            std::unordered_set<std::uint32_t> enclosure_faces;
            for (const auto& item : output)
                if (item.enclosure_group != 0)
                    enclosure_faces.insert(item.source_faces.begin(),
                                           item.source_faces.end());
            choice.enclosure_responsibility_faces = enclosure_faces.size();
            std::vector<std::uint32_t> unique_enclosure_faces(
                enclosure_faces.begin(), enclosure_faces.end());
            choice.enclosure_responsibility_area = facesArea(
                mesh, unique_enclosure_faces);
        }
        if (std::any_of(output.begin(), output.end(),
                        [](const OutputPrimitive& item)
                        {
                            return item.primitive.kind == Kind::CylindricalBand ||
                                   item.primitive.kind == Kind::ConicalBand;
                        }))
        {
            std::unordered_set<std::uint32_t> analytic_faces;
            for (const OutputPrimitive& item : output)
            {
                if (item.primitive.kind != Kind::CylindricalBand &&
                    item.primitive.kind != Kind::ConicalBand)
                    continue;
                analytic_faces.insert(item.source_faces.begin(),
                                      item.source_faces.end());
            }
            choice.analytic_responsibility_faces = analytic_faces.size();
            std::vector<std::uint32_t> unique_analytic_faces(
                analytic_faces.begin(), analytic_faces.end());
            choice.analytic_responsibility_area = facesArea(
                mesh, unique_analytic_faces);
        }
        auto node = std::make_shared<UnifiedProxyChoiceNode>();
        node->direct_output = std::move(output);
        choice.node = std::move(node);
        return choice;
    };
    const auto makeExactFallback = [&]()
    {
        std::vector<OutputPrimitive> exact;
        exact.reserve(mesh.faces.size());
        for (std::uint32_t face_id = 0; face_id < mesh.faces.size(); ++face_id)
        {
            if (!included[face_id]) continue;
            const Face& face = mesh.faces[face_id];
            Primitive triangle;
            triangle.kind = Kind::Triangle;
            for (int corner = 0; corner < 3; ++corner)
                triangle.triangle[corner] = mesh.vertices[face[corner]];
            exact.push_back({std::move(triangle), {face_id}});
        }
        promoteToSemanticPrimitives(exact);
        if (selected_coverage_certificates != nullptr)
            *selected_coverage_certificates = exact;
        stats.unified_exact_fallback_selected = true;
        stats.unified_selected_workload = static_cast<double>(exact.size());
        return exact;
    };
    const auto combineChoices = [](const UnifiedProxyChoice& left,
                                   const UnifiedProxyChoice& right)
    {
        UnifiedProxyChoice choice;
        choice.pending_left = left.node;
        choice.pending_right = right.node;
        choice.workload = left.workload + right.workload;
        choice.open_mean_distance_ratio = left.open_mean_distance_ratio +
                                          right.open_mean_distance_ratio;
        choice.open_maximum_distance = std::max(
            left.open_maximum_distance, right.open_maximum_distance);
        choice.open_maximum_distance_ratio = std::max(
            left.open_maximum_distance_ratio,
            right.open_maximum_distance_ratio);
        choice.analytic_responsibility_faces =
            left.analytic_responsibility_faces +
            right.analytic_responsibility_faces;
        choice.analytic_responsibility_area =
            left.analytic_responsibility_area +
            right.analytic_responsibility_area;
        choice.enclosure_responsibility_faces =
            left.enclosure_responsibility_faces +
            right.enclosure_responsibility_faces;
        choice.enclosure_responsibility_area =
            left.enclosure_responsibility_area +
            right.enclosure_responsibility_area;
        choice.can_refine = left.can_refine || right.can_refine;
        return choice;
    };
    const auto materializeChoice = [](UnifiedProxyChoice& choice)
    {
        if (choice.node) return;
        if (!choice.pending_left || !choice.pending_right) return;
        auto node = std::make_shared<UnifiedProxyChoiceNode>();
        node->left = std::move(choice.pending_left);
        node->right = std::move(choice.pending_right);
        choice.node = std::move(node);
    };
    bool preserve_significant_analytic_choices = false;
    bool final_frontier_pruning = false;
    const auto pruneFrontier = [&](std::vector<UnifiedProxyChoice> choices)
    {
        const auto prune_started = std::chrono::steady_clock::now();
        ++stats.unified_frontier_prune_calls;
        stats.unified_frontier_prune_choices += choices.size();
        stats.unified_frontier_prune_max_choices = std::max(
            stats.unified_frontier_prune_max_choices, choices.size());
        if (choices.size() >= 1024)
            diagnostic("frontier_prune_begin", choices.size(),
                       stats.unified_frontier_prune_calls);
        // Error is a Pareto dimension, so exact and approximate alternatives
        // must both survive. Bound that tradeoff curve aggressively: recursive
        // combination is quadratic in frontier size, while ten uniformly
        // sampled representatives retain the zero-error and near-limit ends.
        const std::size_t maximum_frontier_size = 10;
        std::erase_if(choices, [&](const UnifiedProxyChoice& choice)
        {
            return (!choice.node &&
                    (!choice.pending_left || !choice.pending_right)) ||
                   !std::isfinite(choice.workload) ||
                   (!choice.can_refine &&
                    choice.open_maximum_distance >
                        maximum_open_error_distance + tolerance);
        });
        std::vector<bool> dominated(choices.size(), false);
        stats.unified_frontier_prune_comparisons +=
            static_cast<std::uint64_t>(choices.size()) * choices.size();
        for (std::size_t candidate = 0; candidate < choices.size(); ++candidate)
            for (std::size_t other = 0; other < choices.size(); ++other)
            {
                if (candidate == other) continue;
                const bool no_worse =
                    choices[other].workload <= choices[candidate].workload &&
                    choices[other].open_maximum_distance <=
                        choices[candidate].open_maximum_distance + tolerance &&
                    choices[other].enclosure_responsibility_area >=
                        choices[candidate].enclosure_responsibility_area &&
                    (!preserve_significant_analytic_choices ||
                     choices[other].analytic_responsibility_area >=
                          choices[candidate].analytic_responsibility_area);
                const bool strictly_better =
                    choices[other].workload < choices[candidate].workload ||
                    choices[other].open_maximum_distance + tolerance <
                        choices[candidate].open_maximum_distance ||
                    choices[other].enclosure_responsibility_area >
                        choices[candidate].enclosure_responsibility_area ||
                    (preserve_significant_analytic_choices &&
                      choices[other].analytic_responsibility_area >
                          choices[candidate].analytic_responsibility_area);
                if (no_worse && strictly_better)
                {
                    dominated[candidate] = true;
                    break;
                }
            }
        std::vector<UnifiedProxyChoice> frontier;
        frontier.reserve(choices.size());
        for (std::size_t index = 0; index < choices.size(); ++index)
            if (!dominated[index]) frontier.push_back(std::move(choices[index]));
        std::sort(frontier.begin(), frontier.end(), [](const auto& first,
                                                       const auto& second)
        {
            if (first.open_maximum_distance != second.open_maximum_distance)
                return first.open_maximum_distance <
                       second.open_maximum_distance;
            if (first.enclosure_responsibility_area !=
                second.enclosure_responsibility_area)
                return first.enclosure_responsibility_area >
                       second.enclosure_responsibility_area;
            return first.workload < second.workload;
        });
        const auto finishPrune = [&]
        {
            stats.unified_frontier_prune_seconds += std::chrono::duration<double>(
                std::chrono::steady_clock::now() - prune_started).count();
            if (choices.size() >= 1024)
                diagnostic("frontier_prune_end", choices.size(), frontier.size());
        };
        if (frontier.size() <= maximum_frontier_size)
        {
            for (auto& choice : frontier) materializeChoice(choice);
            finishPrune();
            return frontier;
        }
        std::vector<UnifiedProxyChoice> sampled;
        sampled.reserve(maximum_frontier_size);
        std::set<std::size_t> selected_indices;
        const auto selectBest = [&](const auto& better)
        {
            const auto iterator = std::min_element(
                frontier.begin(), frontier.end(), better);
            selected_indices.insert(static_cast<std::size_t>(
                std::distance(frontier.begin(), iterator)));
        };
        selected_indices.insert(0);
        selected_indices.insert(frontier.size() - 1);
        selectBest([](const auto& first, const auto& second)
        {
            return first.workload < second.workload;
        });
        selectBest([&](const auto& first, const auto& second)
        {
            const bool first_feasible =
                first.open_maximum_distance <=
                    maximum_open_error_distance + tolerance;
            const bool second_feasible =
                second.open_maximum_distance <=
                    maximum_open_error_distance + tolerance;
            if (first_feasible != second_feasible) return first_feasible;
            if (first.workload != second.workload)
                return first.workload < second.workload;
            return first.open_maximum_distance < second.open_maximum_distance;
        });
        selectBest([](const auto& first, const auto& second)
        {
            if (first.analytic_responsibility_area !=
                second.analytic_responsibility_area)
                return first.analytic_responsibility_area >
                       second.analytic_responsibility_area;
            return first.workload < second.workload;
        });
        selectBest([](const auto& first, const auto& second)
        {
            const bool first_analytic =
                first.analytic_responsibility_area > 0.0;
            const bool second_analytic =
                second.analytic_responsibility_area > 0.0;
            if (first_analytic != second_analytic) return first_analytic;
            if (!first_analytic) return first.workload < second.workload;
            if (first.workload != second.workload)
                return first.workload < second.workload;
            return first.open_maximum_distance <
                   second.open_maximum_distance;
        });
        selectBest([](const auto& first, const auto& second)
        {
            if (first.enclosure_responsibility_area !=
                second.enclosure_responsibility_area)
                return first.enclosure_responsibility_area >
                       second.enclosure_responsibility_area;
            return first.workload < second.workload;
        });
        // Final union cleanup can increase the sampled error relative to this
        // local estimate. Reserve representatives densely below the high-error
        // end, where the minimum-work feasible candidate normally lies.
        if (final_frontier_pruning)
            for (const std::size_t sixteenth : {15u, 14u, 13u, 12u, 8u, 4u})
            {
                if (selected_indices.size() >= maximum_frontier_size) break;
                selected_indices.insert(
                    sixteenth * (frontier.size() - 1) / 16);
            }
        // Keep every mandatory extreme above, then repeatedly choose the point
        // farthest from the retained indices. The previous increasing 1/31,
        // 2/31, ... scan stopped as soon as eight entries existed and therefore
        // sampled almost exclusively from the low-error end of large fronts.
        while (selected_indices.size() < maximum_frontier_size)
        {
            std::size_t best_index = 0;
            std::size_t best_separation = 0;
            for (std::size_t index = 0; index < frontier.size(); ++index)
            {
                if (selected_indices.contains(index)) continue;
                std::size_t separation = frontier.size();
                for (const auto selected_index : selected_indices)
                    separation = std::min(
                        separation,
                        index > selected_index ? index - selected_index
                                               : selected_index - index);
                if (separation > best_separation)
                {
                    best_separation = separation;
                    best_index = index;
                }
            }
            if (best_separation == 0) break;
            selected_indices.insert(best_index);
        }
        for (const std::size_t index : selected_indices)
        {
            if (sampled.size() == maximum_frontier_size) break;
            sampled.push_back(std::move(frontier[index]));
        }
        for (auto& choice : sampled) materializeChoice(choice);
        finishPrune();
        return sampled;
    };

    PrimitiveMeshAnalysisOptions planar_options = options;
    planar_options.coplanar_relative_tolerance = std::max(
        options.coplanar_relative_tolerance,
        options.analytic_surface_relative_tolerance);

    const auto connected_components_started = std::chrono::steady_clock::now();
    const auto connected_components = faceComponentsApproximate(
        mesh, included,
        model_diagonal * options.analytic_surface_relative_tolerance);
    stats.unified_connected_components = connected_components.size();
    stats.unified_connected_component_seconds = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - connected_components_started).count();
    diagnostic("connected_components", connected_components.size());
    std::vector<std::uint32_t> face_component(
        mesh.faces.size(), std::numeric_limits<std::uint32_t>::max());
    for (std::uint32_t component = 0;
         component < connected_components.size(); ++component)
        for (const auto face : connected_components[component])
            face_component[face] = component;
    std::vector<std::vector<std::uint32_t>> component_regions(
        connected_components.size());
    for (std::uint32_t region = 0; region < regions.size(); ++region)
        component_regions[face_component[regions[region].front()]].push_back(region);
    std::vector<Vec3> region_centers(regions.size(), Vec3::Zero());
    std::vector<Bounds> region_bounds(regions.size());
    std::vector<std::size_t> region_weights(regions.size(), 0);
    std::vector<std::uint32_t> face_region(
        mesh.faces.size(), std::numeric_limits<std::uint32_t>::max());
    for (std::size_t region = 0; region < regions.size(); ++region)
    {
        const auto [lower, upper] = faceBounds(mesh, regions[region]);
        region_centers[region] = (lower + upper) * 0.5;
        region_bounds[region] = {lower, upper};
        region_weights[region] = regions[region].size();
        for (const auto face : regions[region])
            face_region[face] = static_cast<std::uint32_t>(region);
    }

    struct AnalyticSweepCandidate
    {
        Primitive fit;
        std::vector<std::uint32_t> regions;
        std::vector<std::uint32_t> faces;
        double axial_aspect = 0.0;
        double responsibility_area = 0.0;
    };
    std::vector<AnalyticSweepCandidate> sweep_candidates;
    const auto analyticSweepIsSignificant = [&](
        const AnalyticSweepCandidate& sweep)
    {
        const double area_ratio = sweep.responsibility_area /
            std::max(model_surface_area, 1.0e-30);
        // A long shaft can be important with comparatively little surface
        // area, while a short drum is important because it owns a substantial
        // part of the model.  Treat these as one continuous criterion instead
        // of requiring every useful sweep to be either slender or dominant.
        // Keeping such a sweep in the model-level frontier prevents a coarse
        // enclosing hierarchy from replacing its curved body with intersecting
        // box shells before final workload and error are evaluated.
        const double aspect_credit = std::max(1.0, sweep.axial_aspect);
        return area_ratio * aspect_credit >= 0.125;
    };
    const auto analytic_started = std::chrono::steady_clock::now();
    if (options.allow_frustum)
    {
        Vec3 model_lower = Vec3::Constant(
            std::numeric_limits<double>::infinity());
        Vec3 model_upper = Vec3::Constant(
            -std::numeric_limits<double>::infinity());
        for (const Vec3& vertex : mesh.vertices)
        {
            model_lower = model_lower.cwiseMin(vertex);
            model_upper = model_upper.cwiseMax(vertex);
        }
        const Vec3 model_center = (model_lower + model_upper) * 0.5;
        std::vector<OutputPrimitive> circular_caps;
        std::size_t ignored_holes = 0;
        double ignored_volume = 0.0;
        diagnostic("approximate_planar_regions", approximate_regions.size());
        for (const auto& planar_faces : approximate_regions)
        {
            auto cap = classifyApproximateCircularSurface(
                mesh, planar_faces, options, model_center, model_diagonal,
                model_volume, ignored_holes, ignored_volume);
            if (cap.size() == 1 &&
                (cap.front().primitive.kind == Kind::Disk ||
                 cap.front().primitive.kind == Kind::Annulus))
                circular_caps.push_back(std::move(cap.front()));
        }
        diagnostic("circular_caps", circular_caps.size());

        const double axial_tolerance = std::max(
            model_diagonal * options.analytic_surface_relative_tolerance,
            tolerance);
        const auto appendCertifiedSweepCandidate = [&]
        (const std::vector<std::uint32_t>& certified_faces, const Vec3& axis)
        {
            if (certified_faces.size() < 12) return;
            const auto certified_vertices = uniqueVertices(
                mesh, certified_faces);
            if (certified_vertices.size() < 12) return;
            Primitive fit = fitConeOnAxis(
                mesh, certified_vertices, axis, options.frustum_segments,
                axial_tolerance);
            if (fit.height <= 4.0 * axial_tolerance)
            {
                ++stats.unified_lateral_reject_geometry;
                return;
            }
            ConeFitFailure failure = ConeFitFailure::None;
            ConeFitDiagnostics fit_diagnostics;
            if (!coneFitIsCertified(
                    mesh, certified_faces, fit,
                    options.circle_radial_tolerance, axial_tolerance,
                    &failure, &fit_diagnostics))
            {
                switch (failure)
                {
                case ConeFitFailure::Geometry:
                    ++stats.unified_lateral_reject_geometry;
                    break;
                case ConeFitFailure::EndRing:
                    ++stats.unified_lateral_reject_end_ring;
                    {
                        const std::size_t minimum_points = std::min(
                            fit_diagnostics.base_ring.point_count,
                            fit_diagnostics.top_ring.point_count);
                        const double spread = std::max(
                            fit_diagnostics.base_ring.radial_spread_ratio,
                            fit_diagnostics.top_ring.radial_spread_ratio);
                        const double gap = std::max(
                            fit_diagnostics.base_ring.angular_gap_ratio,
                            fit_diagnostics.top_ring.angular_gap_ratio);
                        const double point_score = minimum_points == 0
                            ? std::numeric_limits<double>::infinity()
                            : std::max(6.0 /
                                static_cast<double>(minimum_points), 1.0);
                        const double spread_score = spread /
                            std::max(2.0 * options.circle_radial_tolerance,
                                     1.0e-30);
                        const double score = std::max(
                            {point_score, spread_score, gap});
                        if (std::isfinite(score) &&
                            (stats.unified_lateral_best_end_ring_score < 0.0 ||
                             score < stats.unified_lateral_best_end_ring_score))
                        {
                            stats.unified_lateral_best_end_ring_score = score;
                            stats.unified_lateral_best_end_ring_min_points =
                                minimum_points;
                            stats.unified_lateral_best_end_ring_spread_ratio =
                                spread;
                            stats.unified_lateral_best_end_ring_gap_ratio = gap;
                        }
                    }
                    break;
                case ConeFitFailure::Normal:
                    ++stats.unified_lateral_reject_normal;
                    break;
                case ConeFitFailure::Radial:
                    ++stats.unified_lateral_reject_radial;
                    break;
                case ConeFitFailure::None:
                    break;
                }
                return;
            }
            ++stats.unified_lateral_sweep_certified;

            std::vector<std::uint32_t> claimed_regions;
            std::vector<std::uint32_t> claimed_faces;
            for (std::uint32_t region = 0; region < regions.size(); ++region)
            {
                bool contained = true;
                for (const auto vertex : uniqueVertices(mesh, regions[region]))
                {
                    const Vec3 local = fit.axes.transposeMultiply(
                        mesh.vertices[vertex] - fit.center);
                    if (local.x() < -axial_tolerance ||
                        local.x() > fit.height + axial_tolerance)
                    {
                        contained = false;
                        break;
                    }
                    const double fraction = std::clamp(
                        local.x() / fit.height, 0.0, 1.0);
                    const double radius =
                        fit.base_radius * (1.0 - fraction) +
                        fit.top_radius * fraction;
                    if (std::hypot(local.y(), local.z()) >
                        radius * (1.0 + 2.0 *
                            options.circle_radial_tolerance) + axial_tolerance)
                    {
                        contained = false;
                        break;
                    }
                }
                if (!contained) continue;
                claimed_regions.push_back(region);
                claimed_faces.insert(claimed_faces.end(),
                                     regions[region].begin(),
                                     regions[region].end());
            }
            if (claimed_regions.empty() || claimed_faces.size() < 12) return;
            std::sort(claimed_regions.begin(), claimed_regions.end());
            const bool duplicate = std::any_of(
                sweep_candidates.begin(), sweep_candidates.end(),
                [&](const AnalyticSweepCandidate& candidate)
                { return candidate.regions == claimed_regions; });
            if (duplicate) return;
            const double maximum_radius = std::max(
                fit.base_radius, fit.top_radius);
            const double axial_aspect = fit.height /
                std::max(maximum_radius, axial_tolerance);
            const double responsibility_area = facesArea(mesh, claimed_faces);
            sweep_candidates.push_back(
                {std::move(fit), std::move(claimed_regions),
                  std::move(claimed_faces),
                  axial_aspect, responsibility_area});
            ++stats.unified_lateral_sweep_candidates;
        };
        const auto appendTrimmedCylinderCandidate = [&]
        (const std::vector<std::uint32_t>& support_faces, const Vec3& axis)
        {
            if (support_faces.size() < 12) return;
            std::vector<std::uint8_t> supported(mesh.faces.size(), 0);
            for (const auto face : support_faces) supported[face] = 1;
            std::vector<std::uint32_t> claimed_regions;
            std::vector<std::uint32_t> claimed_faces;
            std::vector<std::uint8_t> considered_region(regions.size(), 0);
            for (const auto face : support_faces)
            {
                const auto region = face_region[face];
                if (region >= regions.size() || considered_region[region])
                    continue;
                considered_region[region] = 1;
                if (!std::all_of(
                        regions[region].begin(), regions[region].end(),
                        [&](const auto region_face)
                        { return supported[region_face] != 0; }))
                    continue;
                claimed_regions.push_back(region);
                claimed_faces.insert(claimed_faces.end(),
                                     regions[region].begin(),
                                     regions[region].end());
            }
            if (claimed_faces.size() < 12) return;
            const auto vertices = uniqueVertices(mesh, claimed_faces);
            if (vertices.size() < 12) return;

            const Mat3 frame = orthonormalFrame(axis);
            Vec3 mean = Vec3::Zero();
            for (const auto vertex : vertices) mean += mesh.vertices[vertex];
            mean /= static_cast<double>(vertices.size());
            std::vector<Vec2> projected;
            projected.reserve(vertices.size());
            double lower = std::numeric_limits<double>::infinity();
            double upper = -std::numeric_limits<double>::infinity();
            for (const auto vertex : vertices)
            {
                const Vec3 local = frame.transposeMultiply(
                    mesh.vertices[vertex] - mean);
                projected.emplace_back(local.y(), local.z());
                lower = std::min(lower, local.x());
                upper = std::max(upper, local.x());
            }
            const double height = upper - lower;
            if (height <= 4.0 * axial_tolerance) return;
            const auto circle_center = leastSquaresCircleCenter(projected);
            if (!circle_center) return;
            std::vector<double> radii;
            std::vector<double> angles;
            radii.reserve(projected.size());
            angles.reserve(projected.size());
            for (const Vec2& point : projected)
            {
                const Vec2 radial = point - *circle_center;
                radii.push_back(radial.norm());
                double angle = std::atan2(radial.y(), radial.x());
                if (angle < 0.0) angle += 2.0 * std::numbers::pi;
                angles.push_back(angle);
            }
            const double mean_radius = std::accumulate(
                radii.begin(), radii.end(), 0.0) / radii.size();
            const auto [minimum_radius, maximum_radius] =
                std::minmax_element(radii.begin(), radii.end());
            if (mean_radius <= axial_tolerance ||
                (*maximum_radius - *minimum_radius) / mean_radius >
                    2.0 * options.circle_radial_tolerance)
                return;
            std::sort(angles.begin(), angles.end());
            double maximum_gap = angles.front() + 2.0 * std::numbers::pi -
                                 angles.back();
            for (std::size_t index = 1; index < angles.size(); ++index)
                maximum_gap = std::max(
                    maximum_gap, angles[index] - angles[index - 1]);
            constexpr double minimum_angular_coverage = 0.65;
            if (maximum_gap >
                (1.0 - minimum_angular_coverage) * 2.0 * std::numbers::pi)
                return;

            for (const auto face_id : claimed_faces)
            {
                const Face& face = mesh.faces[face_id];
                Vec3 normal =
                    (mesh.vertices[face[1]] - mesh.vertices[face[0]])
                        .cross(mesh.vertices[face[2]] -
                               mesh.vertices[face[0]]);
                const double normal_length = normal.norm();
                if (normal_length <= 1.0e-30) return;
                normal /= normal_length;
                Vec3 centroid = Vec3::Zero();
                for (const auto vertex : face) centroid += mesh.vertices[vertex];
                centroid /= 3.0;
                const Vec3 centroid_local = frame.transposeMultiply(
                    centroid - mean);
                Vec2 radial2{
                    centroid_local.y() - circle_center->x(),
                    centroid_local.z() - circle_center->y()};
                if (radial2.norm() <= axial_tolerance) return;
                radial2 /= radial2.norm();
                const Vec3 radial_world = frame.col(1) * radial2.x() +
                                          frame.col(2) * radial2.y();
                if (std::abs(normal.dot(axis)) > 0.20 ||
                    std::abs(normal.dot(radial_world)) < 0.85)
                    return;
            }

            Primitive fit;
            fit.kind = Kind::CylindricalBand;
            fit.axes = frame;
            fit.height = height;
            fit.center = mean + frame.col(1) * circle_center->x() +
                frame.col(2) * circle_center->y() +
                axis.normalized() * lower;
            fit.segments = std::max<std::uint32_t>(
                options.frustum_segments, 12);
            const double polygon_expansion = 1.0 / std::cos(
                std::numbers::pi / static_cast<double>(fit.segments));
            fit.base_radius = *maximum_radius * polygon_expansion *
                              (1.0 + 1.0e-12);
            fit.top_radius = fit.base_radius;
            fit.band_axial_ranges.assign(
                fit.segments,
                {std::numeric_limits<double>::infinity(),
                 -std::numeric_limits<double>::infinity()});
            const double angle_step = 2.0 * std::numbers::pi /
                                      static_cast<double>(fit.segments);
            for (const auto face_id : claimed_faces)
            {
                const Face& face = mesh.faces[face_id];
                std::array<double, 3> face_angles{};
                double face_lower = std::numeric_limits<double>::infinity();
                double face_upper = -std::numeric_limits<double>::infinity();
                for (std::size_t corner = 0; corner < face.size(); ++corner)
                {
                    const Vec3 local = fit.axes.transposeMultiply(
                        mesh.vertices[face[corner]] - fit.center);
                    face_lower = std::min(face_lower, local.x());
                    face_upper = std::max(face_upper, local.x());
                    face_angles[corner] = std::atan2(local.z(), local.y());
                    if (face_angles[corner] < 0.0)
                        face_angles[corner] += 2.0 * std::numbers::pi;
                }
                for (std::size_t corner = 1; corner < face_angles.size(); ++corner)
                {
                    while (face_angles[corner] - face_angles[0] >
                           std::numbers::pi)
                        face_angles[corner] -= 2.0 * std::numbers::pi;
                    while (face_angles[corner] - face_angles[0] <
                           -std::numbers::pi)
                        face_angles[corner] += 2.0 * std::numbers::pi;
                }
                const auto [face_minimum_angle, face_maximum_angle] =
                    std::minmax_element(
                        face_angles.begin(), face_angles.end());
                if (*face_maximum_angle - *face_minimum_angle >
                    1.25 * angle_step)
                    return;
                const int first_sample = static_cast<int>(
                    std::floor(*face_minimum_angle / angle_step));
                const int last_sample = static_cast<int>(
                    std::ceil(*face_maximum_angle / angle_step));
                for (int sample = first_sample; sample <= last_sample; ++sample)
                {
                    const std::size_t cyclic = static_cast<std::size_t>(
                        (sample % static_cast<int>(fit.segments) +
                         static_cast<int>(fit.segments)) %
                        static_cast<int>(fit.segments));
                    fit.band_axial_ranges[cyclic][0] = std::min(
                        fit.band_axial_ranges[cyclic][0], face_lower);
                    fit.band_axial_ranges[cyclic][1] = std::max(
                        fit.band_axial_ranges[cyclic][1], face_upper);
                }
            }
            std::size_t valid_ranges = 0;
            for (const auto& range : fit.band_axial_ranges)
                valid_ranges += std::isfinite(range[0]) &&
                                std::isfinite(range[1]);
            if (valid_ranges < static_cast<std::size_t>(
                    std::ceil(minimum_angular_coverage * fit.segments)))
                return;
            for (std::size_t index = 0;
                 index < fit.band_axial_ranges.size(); ++index)
            {
                auto& range = fit.band_axial_ranges[index];
                if (!std::isfinite(range[0]) || !std::isfinite(range[1]))
                {
                    std::size_t previous = index;
                    std::size_t next = index;
                    do previous = (previous + fit.segments - 1) % fit.segments;
                    while (!std::isfinite(
                        fit.band_axial_ranges[previous][0]));
                    do next = (next + 1) % fit.segments;
                    while (!std::isfinite(fit.band_axial_ranges[next][0]));
                    range[0] = std::min(
                        fit.band_axial_ranges[previous][0],
                        fit.band_axial_ranges[next][0]);
                    range[1] = std::max(
                        fit.band_axial_ranges[previous][1],
                        fit.band_axial_ranges[next][1]);
                }
                range[0] = std::max(0.0, range[0] - axial_tolerance);
                range[1] = std::min(height, range[1] + axial_tolerance);
                if (range[1] <= range[0] + axial_tolerance) return;
            }

            std::sort(claimed_regions.begin(), claimed_regions.end());
            const bool duplicate = std::any_of(
                sweep_candidates.begin(), sweep_candidates.end(),
                [&](const AnalyticSweepCandidate& candidate)
                {
                    return candidate.regions == claimed_regions &&
                        !candidate.fit.band_axial_ranges.empty();
                });
            if (duplicate) return;
            ++stats.unified_lateral_sweep_certified;
            ++stats.unified_lateral_sweep_candidates;
            const double responsibility_area = facesArea(
                mesh, claimed_faces);
            sweep_candidates.push_back(
                {std::move(fit), std::move(claimed_regions),
                 std::move(claimed_faces), height / mean_radius,
                 responsibility_area});
        };
        for (std::size_t first = 0; first < circular_caps.size(); ++first)
            for (std::size_t second = first + 1;
                 second < circular_caps.size(); ++second)
            {
                const Primitive& first_cap = circular_caps[first].primitive;
                const Primitive& second_cap = circular_caps[second].primitive;
                Vec3 delta = second_cap.center - first_cap.center;
                const double height = delta.norm();
                if (height <= 4.0 * axial_tolerance) continue;
                const Vec3 axis = delta / height;
                if (std::abs(first_cap.axes.col(0).dot(axis)) < 0.995 ||
                    std::abs(second_cap.axes.col(0).dot(axis)) < 0.995)
                    continue;
                const double maximum_radius = std::max(
                    first_cap.base_radius, second_cap.base_radius);
                const double minimum_radius = std::min(
                    first_cap.base_radius, second_cap.base_radius);
                if (minimum_radius <= axial_tolerance ||
                    maximum_radius > 4.0 * minimum_radius)
                    continue;

                std::vector<std::uint32_t> claimed_regions;
                std::vector<std::uint32_t> claimed_faces;
                for (std::uint32_t region = 0; region < regions.size(); ++region)
                {
                    bool contained = true;
                    for (const auto vertex : uniqueVertices(mesh, regions[region]))
                    {
                        const Vec3 offset =
                            mesh.vertices[vertex] - first_cap.center;
                        const double axial = offset.dot(axis);
                        if (axial < -axial_tolerance ||
                            axial > height + axial_tolerance)
                        {
                            contained = false;
                            break;
                        }
                        const double parameter =
                            std::clamp(axial / height, 0.0, 1.0);
                        const double radius =
                            first_cap.base_radius * (1.0 - parameter) +
                            second_cap.base_radius * parameter;
                        const double radial =
                            (offset - axis * axial).norm();
                        if (radial > radius *
                                (1.0 + 2.0 * options.circle_radial_tolerance) +
                            axial_tolerance)
                        {
                            contained = false;
                            break;
                        }
                    }
                    if (!contained) continue;
                    claimed_regions.push_back(region);
                    claimed_faces.insert(claimed_faces.end(),
                                         regions[region].begin(),
                                         regions[region].end());
                }
                if (claimed_regions.size() < 3 || claimed_faces.size() < 12)
                    continue;
                Primitive fit = fitConeOnAxis(
                    mesh, uniqueVertices(mesh, claimed_faces), axis,
                    options.frustum_segments, axial_tolerance);
                if (fit.height <= 4.0 * axial_tolerance) continue;
                std::sort(claimed_regions.begin(), claimed_regions.end());
                const bool duplicate = std::any_of(
                    sweep_candidates.begin(), sweep_candidates.end(),
                    [&](const AnalyticSweepCandidate& candidate)
                    { return candidate.regions == claimed_regions; });
                if (!duplicate)
                {
                    const double responsibility_area = facesArea(
                        mesh, claimed_faces);
                    sweep_candidates.push_back(
                        {std::move(fit), std::move(claimed_regions),
                          std::move(claimed_faces),
                          height / std::max(maximum_radius, axial_tolerance),
                          responsibility_area});
                }
            }

        // A closed or open cylindrical/conical side already carries a stronger
        // analytic certificate than a pair of inferred planar caps. Detect it
        // directly so missing, perforated, or separately indexed caps do not
        // force the side back to dozens of planar leaves. Candidate axes come
        // from the component frame and crosses of representative face normals;
        // the latter recovers a local sweep axis even inside a larger assembly.
        for (const auto& component_faces : connected_components)
        {
            if (component_faces.size() < 12) continue;
            const auto component_vertices = uniqueVertices(
                mesh, component_faces);
            if (component_vertices.size() < 12) continue;

            std::vector<Vec3> candidate_axes;
            const auto addAxis = [&](Vec3 axis)
            {
                const double norm = axis.norm();
                if (norm <= 1.0e-12) return;
                axis /= norm;
                for (int coordinate = 0; coordinate < 3; ++coordinate)
                    if (std::abs(axis[coordinate]) > 1.0e-12)
                    {
                        if (axis[coordinate] < 0.0) axis *= -1.0;
                        break;
                    }
                if (std::any_of(
                        candidate_axes.begin(), candidate_axes.end(),
                        [&](const Vec3& existing)
                        { return std::abs(existing.dot(axis)) >= 0.997; }))
                    return;
                if (candidate_axes.size() < 48)
                    candidate_axes.push_back(axis);
            };
            for (const Mat3& frame : candidateFrames(mesh, component_vertices))
                for (int axis = 0; axis < 3; ++axis)
                    addAxis(static_cast<Vec3>(frame.col(axis)));

            struct NormalSample
            {
                double doubled_area = 0.0;
                Vec3 normal = Vec3::Zero();
            };
            std::vector<NormalSample> normal_samples;
            normal_samples.reserve(component_faces.size());
            for (const auto face_id : component_faces)
            {
                const Face& face = mesh.faces[face_id];
                const Vec3 cross =
                    (mesh.vertices[face[1]] - mesh.vertices[face[0]])
                        .cross(mesh.vertices[face[2]] - mesh.vertices[face[0]]);
                const double area = cross.norm();
                if (area > 1.0e-30)
                    normal_samples.push_back({area, cross / area});
            }
            constexpr std::size_t maximum_normal_samples = 32;
            if (normal_samples.size() > maximum_normal_samples)
            {
                std::partial_sort(
                    normal_samples.begin(),
                    normal_samples.begin() + maximum_normal_samples,
                    normal_samples.end(),
                    [](const NormalSample& first, const NormalSample& second)
                    { return first.doubled_area > second.doubled_area; });
                normal_samples.resize(maximum_normal_samples);
            }
            for (std::size_t first = 0;
                 first < normal_samples.size() && candidate_axes.size() < 48;
                 ++first)
                for (std::size_t second = first + 1;
                     second < normal_samples.size() && candidate_axes.size() < 48;
                     ++second)
                    addAxis(normal_samples[first].normal.cross(
                        normal_samples[second].normal));
            for (const NormalSample& sample : normal_samples)
                addAxis(sample.normal);

            constexpr double axial_normal_bin_width = 0.05;
            constexpr std::size_t axial_normal_bin_count = 17;
            for (const Vec3& axis : candidate_axes)
            {
                ++stats.unified_lateral_sweep_axes;
                std::array<std::vector<std::uint32_t>,
                           axial_normal_bin_count> lateral_bins;
                for (const auto face_id : component_faces)
                {
                    const Face& face = mesh.faces[face_id];
                    Vec3 normal =
                        (mesh.vertices[face[1]] - mesh.vertices[face[0]])
                            .cross(mesh.vertices[face[2]] - mesh.vertices[face[0]]);
                    const double norm = normal.norm();
                    if (norm <= 1.0e-30) continue;
                    normal /= norm;
                    const double axial_normal = std::abs(normal.dot(axis));
                    const std::size_t bin = static_cast<std::size_t>(
                        std::llround(axial_normal / axial_normal_bin_width));
                    if (bin < axial_normal_bin_count)
                        lateral_bins[bin].push_back(face_id);
                }
                for (const auto& bin_faces : lateral_bins)
                {
                    if (bin_faces.size() < 12) continue;
                    auto lateral_components = faceComponentsFromList(
                        mesh, bin_faces);
                    std::sort(
                        lateral_components.begin(), lateral_components.end(),
                        [](const auto& first, const auto& second)
                        { return first.size() > second.size(); });

                    // A flange, seam, or attached detail can disconnect one
                    // analytic side into several angular patches. Grow each
                    // sizeable seed by the fitted radial equation, then run the
                    // normal full-ring certificate on the union. A local arc
                    // is never accepted on its own merely because it fits a
                    // circle.
                    constexpr std::size_t maximum_aggregate_seeds = 16;
                    std::size_t aggregate_seeds = 0;
                    for (const auto& seed_faces : lateral_components)
                    {
                        if (seed_faces.size() < 6 ||
                            aggregate_seeds++ >= maximum_aggregate_seeds)
                            break;
                        const auto seed_vertices = uniqueVertices(
                            mesh, seed_faces);
                        if (seed_vertices.size() < 6) continue;
                        const Primitive seed_fit = fitConeOnAxis(
                            mesh, seed_vertices, axis,
                            options.frustum_segments, axial_tolerance);
                        if (seed_fit.height <= 4.0 * axial_tolerance) continue;

                        std::vector<std::uint32_t> radial_support;
                        for (const auto face_id : bin_faces)
                        {
                            bool on_side = true;
                            for (const auto vertex : mesh.faces[face_id])
                            {
                                const Vec3 local =
                                    seed_fit.axes.transposeMultiply(
                                        mesh.vertices[vertex] - seed_fit.center);
                                if (local.x() < -axial_tolerance ||
                                    local.x() > seed_fit.height + axial_tolerance)
                                {
                                    on_side = false;
                                    break;
                                }
                                const double fraction = std::clamp(
                                    local.x() / seed_fit.height, 0.0, 1.0);
                                const double radius =
                                    seed_fit.base_radius * (1.0 - fraction) +
                                    seed_fit.top_radius * fraction;
                                const double radial = std::hypot(
                                    local.y(), local.z());
                                const double radial_tolerance = std::max(
                                    2.0 * options.circle_radial_tolerance *
                                        std::max(radius, axial_tolerance),
                                    2.0 * axial_tolerance);
                                if (std::abs(radial - radius) > radial_tolerance)
                                {
                                    on_side = false;
                                    break;
                                }
                            }
                            if (on_side) radial_support.push_back(face_id);
                        }
                        if (radial_support.size() > seed_faces.size())
                        {
                            ++stats.unified_lateral_sweep_patches;
                            appendTrimmedCylinderCandidate(
                                radial_support, axis);
                            appendCertifiedSweepCandidate(
                                radial_support, axis);
                        }
                    }

                    for (const auto& lateral_component : lateral_components)
                    {
                        if (lateral_component.size() >= 12)
                            ++stats.unified_lateral_sweep_patches;
                        appendTrimmedCylinderCandidate(
                            lateral_component, axis);
                        appendCertifiedSweepCandidate(
                            lateral_component, axis);
                    }
                }
            }
        }
        std::sort(sweep_candidates.begin(), sweep_candidates.end(),
                  [](const auto& first, const auto& second)
                  {
                      if (first.responsibility_area != second.responsibility_area)
                          return first.responsibility_area >
                                 second.responsibility_area;
                      return first.faces.size() > second.faces.size();
                  });
        if (sweep_candidates.size() > 96) sweep_candidates.resize(96);
        for (const AnalyticSweepCandidate& sweep : sweep_candidates)
        {
            stats.unified_lateral_sweep_max_faces = std::max(
                stats.unified_lateral_sweep_max_faces, sweep.faces.size());
            stats.unified_lateral_sweep_max_axial_aspect = std::max(
                stats.unified_lateral_sweep_max_axial_aspect,
                sweep.axial_aspect);
            stats.unified_lateral_sweep_max_area_ratio = std::max(
                stats.unified_lateral_sweep_max_area_ratio,
                sweep.responsibility_area /
                    std::max(model_surface_area, 1.0e-30));
            if (analyticSweepIsSignificant(sweep))
                ++stats.unified_lateral_significant_candidates;
        }
        diagnostic("analytic_sweeps", sweep_candidates.size());
        preserve_significant_analytic_choices = std::any_of(
            sweep_candidates.begin(), sweep_candidates.end(),
            [&](const AnalyticSweepCandidate& sweep)
            {
                return analyticSweepIsSignificant(sweep);
            });
    }
    stats.unified_analytic_recognition_seconds = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - analytic_started).count();

    const auto splitByOptimizedSah = [&](
        std::vector<std::uint32_t> ids,
        const auto& item_bounds,
        const auto& item_centers,
        const auto& item_weights)
    {
        constexpr int bin_count = 16;
        struct Bin
        {
            Bounds bounds;
            std::size_t weight = 0;
        };
        const auto extend = [](Bounds& target, const Bounds& source)
        {
            target.lower = target.lower.cwiseMin(source.lower);
            target.upper = target.upper.cwiseMax(source.upper);
        };
        const auto area = [](const Bounds& bounds)
        {
            const Vec3 extent = (bounds.upper - bounds.lower).cwiseMax(Vec3::Zero());
            return 2.0 * (extent.x() * extent.y() +
                          extent.y() * extent.z() +
                          extent.z() * extent.x());
        };

        Vec3 center_lower = Vec3::Constant(std::numeric_limits<double>::infinity());
        Vec3 center_upper = Vec3::Constant(-std::numeric_limits<double>::infinity());
        for (const auto id : ids)
        {
            center_lower = center_lower.cwiseMin(item_centers[id]);
            center_upper = center_upper.cwiseMax(item_centers[id]);
        }

        int best_axis = -1;
        int best_split = -1;
        double best_cost = std::numeric_limits<double>::infinity();
        for (int axis = 0; axis < 3; ++axis)
        {
            const double span = center_upper[axis] - center_lower[axis];
            if (span <= tolerance) continue;
            std::array<Bin, bin_count> bins;
            for (const auto id : ids)
            {
                const int bin = std::clamp(static_cast<int>(
                    (item_centers[id][axis] - center_lower[axis]) /
                    span * bin_count), 0, bin_count - 1);
                extend(bins[bin].bounds, item_bounds[id]);
                bins[bin].weight += item_weights[id];
            }
            std::array<Bounds, bin_count> left_bounds;
            std::array<Bounds, bin_count> right_bounds;
            std::array<std::size_t, bin_count> left_weight{};
            std::array<std::size_t, bin_count> right_weight{};
            for (int bin = 0; bin < bin_count; ++bin)
            {
                if (bin > 0)
                {
                    left_bounds[bin] = left_bounds[bin - 1];
                    left_weight[bin] = left_weight[bin - 1];
                }
                if (bins[bin].weight > 0)
                {
                    extend(left_bounds[bin], bins[bin].bounds);
                    left_weight[bin] += bins[bin].weight;
                }
            }
            for (int bin = bin_count - 1; bin >= 0; --bin)
            {
                if (bin + 1 < bin_count)
                {
                    right_bounds[bin] = right_bounds[bin + 1];
                    right_weight[bin] = right_weight[bin + 1];
                }
                if (bins[bin].weight > 0)
                {
                    extend(right_bounds[bin], bins[bin].bounds);
                    right_weight[bin] += bins[bin].weight;
                }
            }
            for (int split = 0; split + 1 < bin_count; ++split)
            {
                if (left_weight[split] == 0 || right_weight[split + 1] == 0)
                    continue;
                const double cost = area(left_bounds[split]) * left_weight[split] +
                    area(right_bounds[split + 1]) * right_weight[split + 1];
                if (cost < best_cost)
                {
                    best_cost = cost;
                    best_axis = axis;
                    best_split = split;
                }
            }
        }

        std::vector<std::uint32_t> left;
        std::vector<std::uint32_t> right;
        if (best_axis >= 0)
        {
            const double span = center_upper[best_axis] - center_lower[best_axis];
            for (const auto id : ids)
            {
                const int bin = std::clamp(static_cast<int>(
                    (item_centers[id][best_axis] - center_lower[best_axis]) /
                    span * bin_count), 0, bin_count - 1);
                (bin <= best_split ? left : right).push_back(id);
            }
        }
        if (left.empty() || right.empty())
        {
            int axis = 0;
            for (int candidate = 1; candidate < 3; ++candidate)
                if (center_upper[candidate] - center_lower[candidate] >
                    center_upper[axis] - center_lower[axis]) axis = candidate;
            std::sort(ids.begin(), ids.end(), [&](const auto first, const auto second)
            {
                return item_centers[first][axis] < item_centers[second][axis];
            });
            const std::size_t middle = ids.size() / 2;
            left.assign(ids.begin(), ids.begin() + middle);
            right.assign(ids.begin() + middle, ids.end());
        }
        return std::pair{std::move(left), std::move(right)};
    };

    const auto splitByAxisMedian = [&](std::vector<std::uint32_t> ids,
                                       const auto& item_centers,
                                       const int axis)
    {
        std::stable_sort(ids.begin(), ids.end(), [&](const auto first,
                                                     const auto second)
        {
            if (item_centers[first][axis] != item_centers[second][axis])
                return item_centers[first][axis] < item_centers[second][axis];
            for (int offset = 1; offset < 3; ++offset)
            {
                const int other_axis = (axis + offset) % 3;
                if (item_centers[first][other_axis] !=
                    item_centers[second][other_axis])
                    return item_centers[first][other_axis] <
                           item_centers[second][other_axis];
            }
            return first < second;
        });
        const std::size_t middle = ids.size() / 2;
        std::vector<std::uint32_t> left(ids.begin(), ids.begin() + middle);
        std::vector<std::uint32_t> right(ids.begin() + middle, ids.end());
        return std::pair{std::move(left), std::move(right)};
    };

    std::vector<std::uint32_t> face_membership_generation(mesh.faces.size(), 0);
    std::uint32_t next_face_membership_generation = 0;
    const auto containsSignificantAnalyticSweep = [&](
        const std::vector<std::uint32_t>& faces)
    {
        if (sweep_candidates.empty()) return false;
        if (++next_face_membership_generation == 0)
        {
            std::fill(face_membership_generation.begin(),
                      face_membership_generation.end(), 0);
            ++next_face_membership_generation;
        }
        for (const auto face : faces)
            face_membership_generation[face] =
                next_face_membership_generation;
        for (const AnalyticSweepCandidate& sweep : sweep_candidates)
        {
            if (sweep.faces.size() < 12 ||
                !analyticSweepIsSignificant(sweep))
                continue;
            if (std::all_of(
                    sweep.faces.begin(), sweep.faces.end(),
                    [&](const auto face)
                    {
                        return face_membership_generation[face] ==
                               next_face_membership_generation;
                    }))
                return true;
        }
        return false;
    };

    const auto considerEnclosingCandidates = [&](
        const std::vector<std::uint32_t>& faces,
        std::vector<UnifiedProxyChoice> frontier)
    {
        if (faces.size() < 12) return pruneFrontier(std::move(frontier));
        const auto vertices = uniqueVertices(mesh, faces);
        const BoxFit box = fitBox(mesh, vertices);
        const Vec3 box_extent = box.half_size * 2.0;
        if (box_extent.x() <= tolerance || box_extent.y() <= tolerance ||
            box_extent.z() <= tolerance) return pruneFrontier(std::move(frontier));
        const int resolution = std::clamp(
            static_cast<int>(options.projection_envelope_resolution / 6), 24, 40);
        (void)resolution;
        std::shared_ptr<const UnifiedProxyChoiceNode> refinement_node;
        if (!frontier.empty())
        {
            const auto refinement = std::min_element(
                frontier.begin(), frontier.end(),
                [](const UnifiedProxyChoice& first,
                   const UnifiedProxyChoice& second)
                {
                    if (first.workload != second.workload)
                        return first.workload < second.workload;
                    return first.open_maximum_distance <
                           second.open_maximum_distance;
                });
            materializeChoice(*refinement);
            refinement_node = refinement->node;
        }
        const auto attachRefinement = [&](UnifiedProxyChoice& choice)
        {
            if (!choice.node || !refinement_node) return;
            auto node = std::make_shared<UnifiedProxyChoiceNode>(*choice.node);
            node->refinement = refinement_node;
            choice.node = std::move(node);
            choice.can_refine = true;
            if (!std::isfinite(choice.workload))
                choice.workload = static_cast<double>(
                    triangulatedFaceCount(choice.node->direct_output));
        };
        // Candidate acceptance is normalized only by the complete model. A
        // local empty-fraction gate makes the same geometric error legal or
        // illegal depending on how the hierarchy happened to partition it,
        // and causes open triangle soups to collapse back to one primitive per
        // source face. The global error term below already makes inaccurate
        // boxes compete against their saved PQSS leaf work.
        if (!containsSignificantAnalyticSweep(faces))
        {
            std::vector<OutputPrimitive> box_output;
            appendBoxRectangles(
                box_output, box, faces, -1, 0.0,
                enclosureGroupForFaces(faces));
            UnifiedProxyChoice box_choice = makeChoice(std::move(box_output));
            const bool box_passed_local_budget =
                std::isfinite(box_choice.workload);
            attachRefinement(box_choice);
            if (faces.size() == mesh.faces.size())
            {
                stats.unified_whole_model_box_workload =
                    box_passed_local_budget
                        ? box_choice.workload : -1.0;
                stats.unified_whole_model_box_passed_local_budget =
                    box_passed_local_budget;
            }
            if (std::isfinite(box_choice.workload))
                frontier.push_back(std::move(box_choice));
        }

        // Do not fit a generic enclosing frustum to an arbitrary hierarchy
        // node.  Curved output must come from the certified sweep recognizer
        // above; otherwise unrelated planar attachments can be absorbed into
        // spurious cylinders/cones that interleave with the box hierarchy.
        return pruneFrontier(std::move(frontier));
    };

    const auto combineFrontiers = [&](const std::vector<UnifiedProxyChoice>& left,
                                      const std::vector<UnifiedProxyChoice>& right)
    {
        std::vector<UnifiedProxyChoice> combined;
        combined.reserve(left.size() * right.size());
        for (const auto& first : left)
            for (const auto& second : right)
            {
                UnifiedProxyChoice choice = combineChoices(first, second);
                if (std::isfinite(choice.workload) &&
                    (choice.can_refine ||
                     choice.open_maximum_distance <=
                         maximum_open_error_distance + tolerance))
                    combined.push_back(std::move(choice));
            }
        return pruneFrontier(std::move(combined));
    };

    PrimitiveMeshAnalysisOptions strict_planar_options = planar_options;
    strict_planar_options.tiny_planar_detail_area_ratio = 0.0;
    strict_planar_options.maximum_local_planar_fill_area_ratio = 0.0;
    std::ofstream region_error_profile(
        performance_path.parent_path() / "region_leaf_error_profile.csv");
    region_error_profile <<
        "region,source_faces,standard_workload,standard_maximum_distance,"
        "strict_workload,strict_maximum_distance,surviving_choices\n";
    std::size_t written_region_diagnostics = 0;
    constexpr std::size_t maximum_region_diagnostics = 16;
    const auto writeRegionSourceObj = [&](const std::filesystem::path& path,
                                          const std::vector<std::uint32_t>& faces)
    {
        std::ofstream stream(path);
        stream << std::setprecision(17);
        for (const Vec3& vertex : mesh.vertices)
            stream << "v " << vertex.x() << ' ' << vertex.y() << ' '
                   << vertex.z() << '\n';
        for (const auto face_id : faces)
        {
            const Face& face = mesh.faces[face_id];
            stream << "f " << face[0] + 1 << ' ' << face[1] + 1 << ' '
                   << face[2] + 1 << '\n';
        }
    };
    const auto writeRegionCandidateObj = [&](const std::filesystem::path& path,
                                             const std::vector<OutputPrimitive>& output)
    {
        const PrimitiveMesh proxy = triangulateOutputPrimitives(output);
        std::ofstream stream(path);
        stream << std::setprecision(17);
        for (const Vec3& vertex : proxy.vertices)
            stream << "v " << vertex.x() << ' ' << vertex.y() << ' '
                   << vertex.z() << '\n';
        for (const Face& face : proxy.faces)
            stream << "f " << face[0] + 1 << ' ' << face[1] + 1 << ' '
                   << face[2] + 1 << '\n';
    };
    const auto writeRegionError = [&](const std::filesystem::path& path,
                                      const FinalOpenErrorAudit& error)
    {
        std::ofstream stream(path);
        stream << std::setprecision(17)
               << "{\"maximum_distance\":" << error.maximum_distance
               << ",\"mean_distance\":" << error.mean_distance
               << ",\"maximum_pair\":{\"proxy\":["
               << error.maximum_proxy_point.x() << ','
               << error.maximum_proxy_point.y() << ','
               << error.maximum_proxy_point.z() << "],\"source\":["
               << error.maximum_source_point.x() << ','
               << error.maximum_source_point.y() << ','
               << error.maximum_source_point.z() << "]}}\n";
    };
    std::vector<std::vector<UnifiedProxyChoice>> region_leaf_choices;
    region_leaf_choices.reserve(regions.size());
    for (std::size_t region_index = 0;
         region_index < regions.size(); ++region_index)
    {
        const auto& region = regions[region_index];
        auto standard_output = timedClassifyFinalRegion(region, planar_options);
        auto strict_output = timedClassifyFinalRegion(
            region, strict_planar_options, -1.0);
        std::vector<UnifiedProxyChoice> choices;
        choices.push_back(makeChoice(standard_output));
        choices.push_back(makeChoice(strict_output));
        const double standard_workload = choices[0].workload;
        const double standard_maximum_distance =
            choices[0].open_maximum_distance;
        const double strict_workload = choices[1].workload;
        const double strict_maximum_distance =
            choices[1].open_maximum_distance;
        auto surviving = pruneFrontier(std::move(choices));
        const bool failed = surviving.empty() ||
            standard_maximum_distance > maximum_open_error_distance ||
            strict_maximum_distance > maximum_open_error_distance;
        if (failed)
            region_error_profile << region_index << ',' << region.size() << ','
                << standard_workload << ',' << standard_maximum_distance << ','
                << strict_workload << ',' << strict_maximum_distance << ','
                << surviving.size() << '\n';
        if (failed && options.write_unified_region_diagnostics &&
            written_region_diagnostics < maximum_region_diagnostics)
        {
            std::ostringstream name;
            name << "region_" << std::setw(4) << std::setfill('0')
                 << region_index;
            const auto directory = performance_path.parent_path() /
                "region_diagnostics" / name.str();
            std::filesystem::create_directories(directory);
            writeRegionSourceObj(directory / "source.obj", region);
            writeRegionCandidateObj(directory / "standard.obj", standard_output);
            writeRegionCandidateObj(directory / "strict.obj", strict_output);
            writeRegionError(
                directory / "standard_error.json",
                measureChargeableOpenSurfaceDistance(
                    open_surface_reference, mesh, standard_output,
                    local_open_error_sample_spacing, tolerance));
            writeRegionError(
                directory / "strict_error.json",
                measureChargeableOpenSurfaceDistance(
                    open_surface_reference, mesh, strict_output,
                    local_open_error_sample_spacing, tolerance));
            ++written_region_diagnostics;
        }
        region_leaf_choices.push_back(std::move(surviving));
    }

    std::function<std::vector<UnifiedProxyChoice>(
        std::vector<std::uint32_t>, int)>
        solve_regions;
    solve_regions = [&](std::vector<std::uint32_t> region_ids,
                        const int split_policy)
        -> std::vector<UnifiedProxyChoice>
    {
        if (region_ids.size() == 1)
            return region_leaf_choices[region_ids.front()];
        auto [left_ids, right_ids] = split_policy < 0
            ? splitByOptimizedSah(
                  region_ids, region_bounds, region_centers, region_weights)
            : splitByAxisMedian(region_ids, region_centers, split_policy);
        auto left = solve_regions(std::move(left_ids), split_policy);
        auto right = solve_regions(std::move(right_ids), split_policy);
        auto baseline = combineFrontiers(left, right);
        std::vector<std::uint32_t> faces;
        for (const auto region : region_ids)
            faces.insert(faces.end(), regions[region].begin(), regions[region].end());
        return considerEnclosingCandidates(
            faces, std::move(baseline));
    };

    const auto hierarchy_started = std::chrono::steady_clock::now();
    std::vector<std::vector<UnifiedProxyChoice>> component_choices;
    std::vector<Vec3> component_centers;
    std::vector<Bounds> component_bounds;
    std::vector<std::size_t> component_weights;
    component_choices.reserve(connected_components.size());
    component_centers.reserve(connected_components.size());
    component_bounds.reserve(connected_components.size());
    component_weights.reserve(connected_components.size());
    for (std::size_t component = 0;
         component < connected_components.size(); ++component)
    {
        diagnostic("component_begin", component,
                   component_regions[component].size());
        std::vector<UnifiedProxyChoice> baseline;
        for (int split_policy = -1; split_policy < 3; ++split_policy)
        {
            auto alternative = solve_regions(
                component_regions[component], split_policy);
            baseline.insert(baseline.end(),
                std::make_move_iterator(alternative.begin()),
                std::make_move_iterator(alternative.end()));
        }
        baseline = pruneFrontier(std::move(baseline));
        component_choices.push_back(considerEnclosingCandidates(
            connected_components[component], std::move(baseline)));
        const auto [lower, upper] = faceBounds(
            mesh, connected_components[component]);
        component_centers.push_back((lower + upper) * 0.5);
        component_bounds.push_back({lower, upper});
        component_weights.push_back(connected_components[component].size());
        if ((component + 1) % 8 == 0 || component + 1 == connected_components.size())
            diagnostic("component_hierarchy", component + 1,
                       connected_components.size());
    }

    std::function<std::vector<UnifiedProxyChoice>(
        std::vector<std::uint32_t>, int)>
        solve_components;
    solve_components = [&](std::vector<std::uint32_t> component_ids,
                           const int split_policy)
        -> std::vector<UnifiedProxyChoice>
    {
        if (component_ids.size() == 1)
            return component_choices[component_ids.front()];
        auto [left_ids, right_ids] = split_policy < 0
            ? splitByOptimizedSah(
                  component_ids, component_bounds, component_centers,
                  component_weights)
            : splitByAxisMedian(component_ids, component_centers, split_policy);
        auto left = solve_components(std::move(left_ids), split_policy);
        auto right = solve_components(std::move(right_ids), split_policy);
        auto baseline = combineFrontiers(left, right);

        std::vector<std::uint32_t> faces;
        for (const auto component : component_ids)
            faces.insert(faces.end(), connected_components[component].begin(),
                         connected_components[component].end());
        return considerEnclosingCandidates(
            faces, std::move(baseline));
    };

    std::vector<std::uint32_t> component_ids(connected_components.size());
    std::iota(component_ids.begin(), component_ids.end(), 0);
    std::vector<UnifiedProxyChoice> frontier;
    for (int split_policy = -1; split_policy < 3; ++split_policy)
    {
        auto alternative = solve_components(component_ids, split_policy);
        frontier.insert(frontier.end(),
            std::make_move_iterator(alternative.begin()),
            std::make_move_iterator(alternative.end()));
    }
    frontier = pruneFrontier(std::move(frontier));
    diagnostic("component_frontier", frontier.size());

    // Circular recognizers only propose atomic leaves. The original hierarchy
    // remains in the frontier, while a second hierarchy is built from accepted
    // non-overlapping sweep leaves plus every still-unclaimed planar region.
    // Consequently a sweep can replace one shaft without forcing the unrelated
    // remainder of the model into one loose envelope.
    std::vector<std::vector<UnifiedProxyChoice>> atom_choices;
    std::vector<std::vector<std::uint32_t>> atom_faces;
    std::vector<Bounds> atom_bounds;
    std::vector<Vec3> atom_centers;
    std::vector<std::size_t> atom_weights;
    std::vector<bool> claimed_regions(regions.size(), false);
    for (const AnalyticSweepCandidate& sweep : sweep_candidates)
    {
        // Certified local arcs remain useful recognition diagnostics, but only
        // a model-significant sweep may become a model-level atom.  Greedily
        // atomizing every small certified arc lets incidental round details
        // accumulate more responsibility than the actual cylindrical body and
        // produces interleaved analytic/enclosing shells.
        if (!analyticSweepIsSignificant(sweep)) continue;
        if (std::any_of(sweep.regions.begin(), sweep.regions.end(),
                        [&](const auto region)
                        { return claimed_regions[region]; }))
            continue;
        std::vector<OutputPrimitive> sweep_output;
        if (!sweep.fit.band_axial_ranges.empty())
            sweep_output.push_back({sweep.fit, sweep.faces});
        else
            appendEnclosingConicalSurfaceAssembly(
                sweep_output, sweep.fit, sweep.faces,
                enclosureGroupForFaces(sweep.faces));
        UnifiedProxyChoice sweep_choice = makeChoice(std::move(sweep_output));
        if (!std::isfinite(sweep_choice.workload))
            continue;
        for (const auto region : sweep.regions) claimed_regions[region] = true;
        atom_choices.push_back({std::move(sweep_choice)});
        atom_faces.push_back(sweep.faces);
        const auto [lower, upper] = faceBounds(mesh, sweep.faces);
        atom_bounds.push_back({lower, upper});
        atom_centers.push_back((lower + upper) * 0.5);
        atom_weights.push_back(sweep.faces.size());
    }
    const std::size_t sweep_atom_count = atom_choices.size();
    for (std::uint32_t region = 0; region < regions.size(); ++region)
    {
        if (claimed_regions[region]) continue;
        atom_choices.push_back(region_leaf_choices[region]);
        atom_faces.push_back(regions[region]);
        atom_bounds.push_back(region_bounds[region]);
        atom_centers.push_back(region_centers[region]);
        atom_weights.push_back(region_weights[region]);
    }
    if (sweep_atom_count > 0)
    {
        std::function<std::vector<UnifiedProxyChoice>(
            std::vector<std::uint32_t>)> solve_atoms;
        solve_atoms = [&](std::vector<std::uint32_t> ids)
            -> std::vector<UnifiedProxyChoice>
        {
            if (ids.size() == 1) return atom_choices[ids.front()];
            auto [left_ids, right_ids] = splitByOptimizedSah(
                ids, atom_bounds, atom_centers, atom_weights);
            auto left = solve_atoms(std::move(left_ids));
            auto right = solve_atoms(std::move(right_ids));
            auto baseline = combineFrontiers(left, right);
            std::vector<std::uint32_t> faces;
            for (const auto id : ids)
                faces.insert(faces.end(), atom_faces[id].begin(),
                             atom_faces[id].end());
            return considerEnclosingCandidates(
                faces, std::move(baseline));
        };
        std::vector<std::uint32_t> atom_ids(atom_choices.size());
        std::iota(atom_ids.begin(), atom_ids.end(), 0);
        auto atomized = solve_atoms(std::move(atom_ids));
        frontier.insert(frontier.end(),
                        std::make_move_iterator(atomized.begin()),
                        std::make_move_iterator(atomized.end()));

        // Also optimize the non-analytic remainder independently before
        // attaching the certified sweeps.  A central drum otherwise becomes
        // an ancestor of many unrelated planar regions in the SAH tree; every
        // such ancestor is intentionally forbidden from becoming a box, so a
        // useful "drum + coarse attachments" solution degenerates toward the
        // exact mesh.  This alternative preserves the sweep as an atomic leaf
        // without preventing aggressive hierarchy choices for the remainder.
        std::vector<UnifiedProxyChoice> separated;
        if (sweep_atom_count < atom_choices.size())
        {
            std::vector<std::uint32_t> remainder_ids(
                atom_choices.size() - sweep_atom_count);
            std::iota(remainder_ids.begin(), remainder_ids.end(),
                      static_cast<std::uint32_t>(sweep_atom_count));
            separated = solve_atoms(std::move(remainder_ids));
        }
        else
            separated = atom_choices.front();
        const std::size_t first_sweep_to_attach =
            sweep_atom_count < atom_choices.size() ? 0 : 1;
        for (std::size_t sweep = first_sweep_to_attach;
             sweep < sweep_atom_count; ++sweep)
            separated = combineFrontiers(
                separated, atom_choices[sweep]);
        frontier.insert(frontier.end(),
                        std::make_move_iterator(separated.begin()),
                        std::make_move_iterator(separated.end()));
    }
    final_frontier_pruning = true;
    frontier = pruneFrontier(std::move(frontier));
    final_frontier_pruning = false;
    stats.unified_hierarchy_seconds = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - hierarchy_started).count();
    std::stable_sort(frontier.begin(), frontier.end(),
        [](const UnifiedProxyChoice& first, const UnifiedProxyChoice& second)
        {
            if (first.workload != second.workload)
                return first.workload < second.workload;
            return first.open_maximum_distance < second.open_maximum_distance;
        });
    diagnostic("final_frontier", frontier.size());
    if (frontier.empty()) return makeExactFallback();
    std::ofstream final_candidate_profile(
        performance_path.parent_path() / "final_candidate_profile.csv");
    final_candidate_profile <<
        "candidate,stage,primitives,triangles,aux_count,raw_workload,"
        "analytic_faces,analytic_area_ratio\n";
    struct FinalCandidateEvaluation
    {
        double workload = 0.0;
        std::size_t enclosure_overlap_pairs = 0;
        double enclosure_overlap_score = 0.0;
        std::vector<OutputPrimitive> output;
        std::vector<OutputPrimitive> coverage_certificates;
    };
    const auto evaluateFinalUnionBoundary = [&](
        std::vector<OutputPrimitive> output,
        const std::size_t candidate_index,
        const double raw_workload,
        const bool preserve_pre_removal_enclosures,
        const bool write_profile = true)
    {
        const auto profileStage = [&](const char* stage,
                                      const std::size_t auxiliary_count = 0)
        {
            if (!write_profile) return;
            final_candidate_profile << candidate_index << ',' << stage << ','
                << output.size() << ',' << triangulatedFaceCount(output) << ','
                << auxiliary_count << ',' << raw_workload << ','
                << frontier[candidate_index % frontier.size()]
                       .analytic_responsibility_faces << ','
                << frontier[candidate_index % frontier.size()]
                       .analytic_responsibility_area /
                       std::max(model_surface_area, 1.0e-30) << '\n';
            final_candidate_profile.flush();
            if (options.write_unified_region_diagnostics &&
                raw_workload <= 1000.0)
            {
                const std::string_view stage_name(stage);
                if (stage_name == "flattened" ||
                    stage_name == "contained_removed" ||
                    stage_name == "occlusion_clipped" ||
                    stage_name == "coplanar_canonicalized" ||
                    stage_name == "near_rectangles_regularized")
                {
                    std::ostringstream name;
                    name << "candidate_" << std::setw(2)
                         << std::setfill('0') << candidate_index;
                    const auto directory = performance_path.parent_path() /
                        "final_candidate_diagnostics" / name.str() / "stages";
                    std::filesystem::create_directories(directory);
                    writeRegionCandidateObj(
                        directory / (std::string(stage) + ".obj"), output);
                }
            }
        };
        profileStage("flattened");
        const auto finalization_started = std::chrono::steady_clock::now();
        std::size_t redundant_enclosure_primitives = 0;
        if (preserve_pre_removal_enclosures)
            output = removeRedundantEnclosureGroupsByUnionCoverage(
                mesh, std::move(output), tolerance,
                redundant_enclosure_primitives);
        profileStage("redundant_enclosures_removed",
                     redundant_enclosure_primitives);
        std::vector<OutputPrimitive> coverage_certificates = output;
        promoteToSemanticPrimitives(coverage_certificates);
        std::size_t contained_primitives = 0;
        output = removePrimitivesInsideEnclosureGroups(
            mesh, std::move(output), tolerance, contained_primitives);
        profileStage("contained_removed", contained_primitives);
        promoteToSemanticPrimitives(output);
        const auto extrusions = recognizeEnclosureGroupExtrusions(
            preserve_pre_removal_enclosures ? coverage_certificates : output,
            tolerance, false);
        std::vector<Bounds> extrusion_bounds;
        extrusion_bounds.reserve(extrusions.size());
        for (const auto& extrusion : extrusions)
        {
            Bounds bounds;
            for (const double depth : {extrusion.lower_distance,
                                       extrusion.upper_distance})
                for (const Vec2& point : extrusion.boundary)
                {
                    const Vec3 vertex = extrusion.origin +
                        extrusion.frame.col(0) * point.x() +
                        extrusion.frame.col(1) * point.y() +
                        extrusion.frame.col(2) * depth;
                    bounds.lower = bounds.lower.cwiseMin(vertex);
                    bounds.upper = bounds.upper.cwiseMax(vertex);
                }
            extrusion_bounds.push_back(bounds);
        }
        std::size_t enclosure_overlap_pairs = 0;
        double enclosure_overlap_score = 0.0;
        for (std::size_t first = 0; first < extrusion_bounds.size(); ++first)
            for (std::size_t second = first + 1;
                 second < extrusion_bounds.size(); ++second)
            {
                const Vec3 overlap =
                    (extrusion_bounds[first].upper.cwiseMin(
                         extrusion_bounds[second].upper) -
                     extrusion_bounds[first].lower.cwiseMax(
                         extrusion_bounds[second].lower))
                        .cwiseMax(Vec3::Zero());
                const double overlap_volume = overlap.prod();
                if (overlap_volume <= tolerance * tolerance * tolerance)
                    continue;
                const double first_volume =
                    (extrusion_bounds[first].upper -
                     extrusion_bounds[first].lower)
                        .cwiseMax(Vec3::Zero()).prod();
                const double second_volume =
                    (extrusion_bounds[second].upper -
                     extrusion_bounds[second].lower)
                        .cwiseMax(Vec3::Zero()).prod();
                ++enclosure_overlap_pairs;
                enclosure_overlap_score += overlap_volume /
                    std::max(std::min(first_volume, second_volume), 1.0e-30);
            }
        std::size_t extrusion_boundary_vertices = 0;
        for (const auto& extrusion : extrusions)
            extrusion_boundary_vertices += extrusion.boundary.size();
        profileStage("extrusions_recognized", extrusions.size());
        if (write_profile)
        {
            final_candidate_profile << candidate_index
                << ",extrusion_boundary_vertices," << output.size() << ','
                << triangulatedFaceCount(output) << ','
                << extrusion_boundary_vertices << ',' << raw_workload << ','
                << frontier[candidate_index % frontier.size()]
                       .analytic_responsibility_faces << ','
                << frontier[candidate_index % frontier.size()]
                       .analytic_responsibility_area /
                       std::max(model_surface_area, 1.0e-30) << '\n';
            final_candidate_profile.flush();
        }
        std::vector<OutputPrimitive> before_occlusion = output;
        std::size_t occlusion_removed = 0;
        OcclusionClipStats occlusion_stats;
        constexpr std::size_t maximum_occlusion_input_workload = 2000;
        const bool skip_occlusion =
            triangulatedFaceCount(output) > maximum_occlusion_input_workload;
        if (!skip_occlusion)
            output = clipPlanarOcclusionByClosedVolumes(
                std::move(output), {}, extrusions, tolerance,
                occlusion_removed, nullptr, false, &occlusion_stats);
        else
            profileStage("occlusion_skipped_workload");
        const FinalCoverageAudit clipped_coverage =
            auditFinalConservativeCoverage(
                mesh, output, tolerance, &coverage_certificates,
                extrusions.empty() ? nullptr : &extrusions,
                excluded_redundant_faces);
        bool occlusion_rolled_back = false;
        if (clipped_coverage.unassigned_source_faces != 0 ||
            clipped_coverage.failed_source_faces != 0)
        {
            output = std::move(before_occlusion);
            occlusion_stats = {};
            occlusion_removed = 0;
            occlusion_rolled_back = true;
            profileStage("occlusion_rollback");
        }
        else
            profileStage("occlusion_clipped");
        CoplanarCanonicalizationStats canonicalization;
        output = canonicalizeCoplanarPrimitiveUnion(
            mesh, std::move(output), tolerance, canonicalization);
        profileStage("coplanar_canonicalized");
        output = regularizeShallowNearRectangles(
            std::move(output), model_surface_area,
            options.tiny_planar_detail_area_ratio,
            maximum_open_error_distance, tolerance);
        CoplanarCanonicalizationStats regularized_canonicalization;
        output = canonicalizeCoplanarPrimitiveUnion(
            mesh, std::move(output), tolerance,
            regularized_canonicalization);
        std::size_t post_regularization_removed = 0;
        std::size_t post_regularization_extrusions = 0;
        OcclusionClipStats post_regularization_clip_stats;
        bool post_regularization_rolled_back = false;
        output = clipRegularizedPlanarOcclusion(
            mesh, std::move(output), tolerance, coverage_certificates,
            &extrusions, excluded_redundant_faces, 1000,
            post_regularization_removed,
            post_regularization_extrusions,
            post_regularization_clip_stats,
            post_regularization_rolled_back);
        CoplanarCanonicalizationStats clipped_canonicalization;
        if (post_regularization_clip_stats.clipped_primitives != 0)
        {
            output = canonicalizeCoplanarPrimitiveUnion(
                mesh, std::move(output), tolerance,
                clipped_canonicalization);
        }
        profileStage(post_regularization_rolled_back
            ? "post_regularization_occlusion_rollback"
            : "post_regularization_occlusion",
            post_regularization_extrusions);
        profileStage("near_rectangles_regularized");
        if (!write_profile)
        {
            stats.removed_contained_primitives =
                redundant_enclosure_primitives + contained_primitives;
            stats.unified_enclosure_extrusions = extrusions.size();
            stats.unified_occlusion_clipped_primitives =
                occlusion_stats.clipped_primitives +
                post_regularization_clip_stats.clipped_primitives;
            stats.unified_occlusion_removed_primitives =
                occlusion_stats.removed_primitives +
                post_regularization_clip_stats.removed_primitives;
            stats.unified_occlusion_input_triangles =
                occlusion_stats.input_triangles +
                post_regularization_clip_stats.input_triangles;
            stats.unified_occlusion_output_triangles =
                occlusion_stats.output_triangles +
                post_regularization_clip_stats.output_triangles;
            stats.unified_occlusion_removed_area =
                occlusion_stats.removed_area +
                post_regularization_clip_stats.removed_area;
            stats.unified_occlusion_rolled_back =
                occlusion_rolled_back || post_regularization_rolled_back;
            stats.unified_occlusion_seconds = std::chrono::duration<double>(
                std::chrono::steady_clock::now() - finalization_started).count();
            stats.canonicalized_coplanar_groups =
                canonicalization.groups +
                regularized_canonicalization.groups +
                clipped_canonicalization.groups;
            stats.removed_coplanar_redundant_primitives =
                canonicalization.removed_primitives +
                regularized_canonicalization.removed_primitives +
                clipped_canonicalization.removed_primitives;
            stats.removed_coplanar_overlap_area =
                canonicalization.removed_overlap_area +
                regularized_canonicalization.removed_overlap_area +
                clipped_canonicalization.removed_overlap_area;
        }
        return FinalCandidateEvaluation{
            static_cast<double>(triangulatedFaceCount(output)),
            enclosure_overlap_pairs, enclosure_overlap_score,
            std::move(output), std::move(coverage_certificates)};
    };
    const auto auditFinalEvaluation = [&](const FinalCandidateEvaluation& evaluation)
    {
        FinalOpenErrorAudit open_error = measureChargeableOpenSurfaceDistance(
            open_surface_reference, mesh, evaluation.output,
            final_open_error_sample_spacing, tolerance,
            &evaluation.coverage_certificates);
        FinalCoverageAudit coverage = auditFinalConservativeCoverage(
            mesh, evaluation.output, tolerance,
            &evaluation.coverage_certificates,
            nullptr, excluded_redundant_faces);
        return std::pair{std::move(open_error), std::move(coverage)};
    };
    const auto patchCoverageFailures =
        [&](FinalCandidateEvaluation evaluation,
            const FinalCoverageAudit& coverage)
    {
        std::vector<OutputPrimitive> exact_patches;
        exact_patches.reserve(coverage.failed_face_ids.size());
        for (const auto face_id : coverage.failed_face_ids)
        {
            if (face_id >= mesh.faces.size()) continue;
            const Face& face = mesh.faces[face_id];
            Primitive triangle;
            triangle.kind = Kind::Triangle;
            for (int corner = 0; corner < 3; ++corner)
                triangle.triangle[corner] = mesh.vertices[face[corner]];
            exact_patches.push_back({std::move(triangle), {face_id}});
        }
        promoteToSemanticPrimitives(exact_patches);
        evaluation.output.insert(
            evaluation.output.end(), exact_patches.begin(), exact_patches.end());
        evaluation.coverage_certificates.insert(
            evaluation.coverage_certificates.end(),
            std::make_move_iterator(exact_patches.begin()),
            std::make_move_iterator(exact_patches.end()));
        evaluation.workload = static_cast<double>(
            triangulatedFaceCount(evaluation.output));
        return evaluation;
    };
    double best_feasible_workload =
        std::numeric_limits<double>::infinity();
    for (std::size_t candidate_index = 0;
         candidate_index < frontier.size(); ++candidate_index)
    {
        UnifiedProxyChoice& choice = frontier[candidate_index];
        std::vector<OutputPrimitive> output;
        flattenUnifiedChoice(choice.node, output);
        const double raw_workload = choice.workload;
        constexpr double industrial_fallback_workload = 10000.0;
        constexpr double acceptable_incumbent_workload = 2000.0;
        if (best_feasible_workload < acceptable_incumbent_workload &&
            raw_workload > industrial_fallback_workload &&
            (!preserve_significant_analytic_choices ||
             choice.analytic_responsibility_faces == 0))
        {
            choice.final_coverage_passed = false;
            ++stats.unified_skipped_industrial_fallback_candidates;
            final_candidate_profile << candidate_index
                << ",industrial_fallback_skipped," << output.size() << ','
                << raw_workload << ",0," << raw_workload << ','
                << choice.analytic_responsibility_faces << ','
                << choice.analytic_responsibility_area /
                       std::max(model_surface_area, 1.0e-30) << '\n';
            final_candidate_profile.flush();
            continue;
        }
        FinalCandidateEvaluation evaluation = evaluateFinalUnionBoundary(
            output, candidate_index, raw_workload, true);
        // Candidates are ordered by their estimated triangle workload.  Once
        // a fully audited feasible result exists, a candidate whose actual
        // finalized workload is already larger cannot win the workload-first
        // objective.  Avoid the expensive sampled-distance and coverage audit;
        // geometry finalization still ran, so this compares export workloads
        // rather than pre-union estimates.
        if (evaluation.workload > best_feasible_workload &&
            (!preserve_significant_analytic_choices ||
             choice.analytic_responsibility_faces == 0))
        {
            choice.workload = evaluation.workload;
            choice.final_coverage_passed = false;
            ++stats.unified_skipped_final_audits;
            final_candidate_profile << candidate_index
                << ",audit_skipped_incumbent," << evaluation.output.size()
                << ',' << evaluation.workload << ",0," << raw_workload << ','
                << choice.analytic_responsibility_faces << ','
                << choice.analytic_responsibility_area /
                       std::max(model_surface_area, 1.0e-30) << '\n';
            final_candidate_profile.flush();
            continue;
        }
        auto [open_error, coverage] = auditFinalEvaluation(evaluation);
        const auto isFeasible = [&](const FinalOpenErrorAudit& error,
                                    const FinalCoverageAudit& audit)
        {
            return audit.unassigned_source_faces == 0 &&
                   audit.failed_source_faces == 0 &&
                   error.maximum_distance <=
                       maximum_open_error_distance + tolerance;
        };
        choice.preserve_pre_removal_enclosures = true;
        if (!isFeasible(open_error, coverage))
        {
            FinalCandidateEvaluation conservative_evaluation =
                evaluateFinalUnionBoundary(
                    std::move(output), frontier.size() + candidate_index,
                    raw_workload, false);
            auto [conservative_error, conservative_coverage] =
                auditFinalEvaluation(conservative_evaluation);
            const bool conservative_feasible = isFeasible(
                conservative_error, conservative_coverage);
            const bool current_coverage_passed =
                coverage.unassigned_source_faces == 0 &&
                coverage.failed_source_faces == 0;
            const bool conservative_coverage_passed =
                conservative_coverage.unassigned_source_faces == 0 &&
                conservative_coverage.failed_source_faces == 0;
            const bool prefer_conservative =
                conservative_feasible ||
                (!isFeasible(open_error, coverage) &&
                 ((conservative_coverage_passed != current_coverage_passed &&
                   conservative_coverage_passed) ||
                  (conservative_coverage_passed == current_coverage_passed &&
                   (conservative_error.maximum_distance <
                        open_error.maximum_distance ||
                    (conservative_error.maximum_distance ==
                         open_error.maximum_distance &&
                     conservative_evaluation.workload <
                         evaluation.workload)))));
            if (prefer_conservative)
            {
                evaluation = std::move(conservative_evaluation);
                open_error = std::move(conservative_error);
                coverage = std::move(conservative_coverage);
                choice.preserve_pre_removal_enclosures = false;
            }
        }
        constexpr std::size_t maximum_final_refinements = 96;
        for (std::size_t refinement_iteration = 0;
             open_error.maximum_distance >
                 maximum_open_error_distance + tolerance &&
             refinement_iteration < maximum_final_refinements;
             ++refinement_iteration)
        {
            bool changed = false;
            const auto refined_node = refineUnifiedChoiceAtPoint(
                choice.node, open_error.maximum_proxy_point,
                std::max(tolerance * 64.0, 1.0e-12), changed);
            auto attributed_refined_node = refined_node;
            if (!changed)
            {
                double nearest_distance_squared =
                    std::numeric_limits<double>::infinity();
                const UnifiedProxyChoiceNode* nearest =
                    findNearestRefinableUnifiedChoice(
                        choice.node, open_error.maximum_proxy_point,
                        nearest_distance_squared);
                if (nearest != nullptr)
                    attributed_refined_node = refineUnifiedChoiceNode(
                        choice.node, nearest, changed);
            }
            if (!changed) break;
            std::vector<OutputPrimitive> refined_output;
            flattenUnifiedChoice(attributed_refined_node, refined_output);
            FinalCandidateEvaluation refined_evaluation =
                evaluateFinalUnionBoundary(
                    std::move(refined_output), candidate_index,
                    evaluation.workload,
                    choice.preserve_pre_removal_enclosures);
            auto [refined_error, refined_coverage] =
                auditFinalEvaluation(refined_evaluation);
            // Replacing an enclosure by its certified child is geometrically
            // finer, but planar union and occlusion cleanup can move the
            // sampled global maximum to another branch and make this scalar
            // temporarily non-monotone.  Keep refining infeasible candidates;
            // feasibility is still decided only by the complete final audit.
            if (choice.analytic_responsibility_faces == 0 &&
                refined_error.maximum_distance >
                    open_error.maximum_distance + tolerance)
                break;
            choice.node = attributed_refined_node;
            evaluation = std::move(refined_evaluation);
            open_error = std::move(refined_error);
            coverage = std::move(refined_coverage);
            ++stats.unified_adaptive_refinements;
        }
        std::size_t patched_coverage_faces = 0;
        while (open_error.maximum_distance <=
                   maximum_open_error_distance + tolerance &&
               (coverage.unassigned_source_faces != 0 ||
                coverage.failed_source_faces != 0) &&
               !coverage.failed_face_ids.empty() &&
               patched_coverage_faces + coverage.failed_face_ids.size() <= 128 &&
               evaluation.workload <= 2000.0)
        {
            patched_coverage_faces += coverage.failed_face_ids.size();
            choice.coverage_patch_faces.insert(
                choice.coverage_patch_faces.end(),
                coverage.failed_face_ids.begin(),
                coverage.failed_face_ids.end());
            evaluation = patchCoverageFailures(
                std::move(evaluation), coverage);
            std::tie(open_error, coverage) =
                auditFinalEvaluation(evaluation);
        }
        choice.workload = evaluation.workload;
        choice.open_mean_distance_ratio =
            open_error.distance_area_integral /
            std::max(model_surface_area * model_diagonal, 1.0e-30);
        choice.open_maximum_distance = open_error.maximum_distance;
        choice.open_maximum_distance_ratio =
            open_error.maximum_distance / std::max(model_diagonal, 1.0e-30);
        choice.open_maximum_proxy_point = open_error.maximum_proxy_point;
        choice.open_maximum_source_point = open_error.maximum_source_point;
        choice.final_enclosure_overlap_pairs =
            evaluation.enclosure_overlap_pairs;
        choice.final_enclosure_overlap_score =
            evaluation.enclosure_overlap_score;
        choice.final_coverage_passed =
            coverage.unassigned_source_faces == 0 &&
            coverage.failed_source_faces == 0;
        choice.final_unassigned_source_faces =
            coverage.unassigned_source_faces;
        choice.final_failed_source_faces = coverage.failed_source_faces;
        if (choice.final_coverage_passed &&
            choice.open_maximum_distance <=
                maximum_open_error_distance + tolerance)
            best_feasible_workload = std::min(
                best_feasible_workload, choice.workload);
        if (options.write_unified_region_diagnostics &&
            coverage.failed_source_faces > 0 &&
            coverage.failed_source_faces <= 64)
        {
            std::ostringstream name;
            name << "candidate_" << std::setw(2) << std::setfill('0')
                 << candidate_index;
            const auto directory = performance_path.parent_path() /
                "final_candidate_diagnostics" / name.str();
            std::filesystem::create_directories(directory);
            writeRegionSourceObj(
                directory / "failed_source_faces.obj",
                coverage.failed_face_ids);
            writeRegionCandidateObj(
                directory / "final_candidate.obj", evaluation.output);
            writeRegionError(directory / "error.json", open_error);
        }
    }
    std::ofstream final_error_profile(
        performance_path.parent_path() / "final_candidate_error_profile.csv");
    final_error_profile <<
        "candidate,triangles,maximum_distance,maximum_distance_ratio,"
        "global_mean_distance_ratio,coverage_passed,unassigned_source_faces,"
        "failed_source_faces,maximum_proxy_x,maximum_proxy_y,maximum_proxy_z,"
        "maximum_source_x,maximum_source_y,maximum_source_z,"
        "enclosure_overlap_pairs,enclosure_overlap_score\n";
    for (std::size_t index = 0; index < frontier.size(); ++index)
    {
        const UnifiedProxyChoice& choice = frontier[index];
        final_error_profile << index << ',' << choice.workload << ','
            << choice.open_maximum_distance << ','
            << choice.open_maximum_distance_ratio << ','
               << choice.open_mean_distance_ratio << ','
               << (choice.final_coverage_passed ? 1 : 0) << ','
               << choice.final_unassigned_source_faces << ','
               << choice.final_failed_source_faces << ','
               << choice.open_maximum_proxy_point.x() << ','
               << choice.open_maximum_proxy_point.y() << ','
               << choice.open_maximum_proxy_point.z() << ','
               << choice.open_maximum_source_point.x() << ','
               << choice.open_maximum_source_point.y() << ','
               << choice.open_maximum_source_point.z() << ','
               << choice.final_enclosure_overlap_pairs << ','
               << choice.final_enclosure_overlap_score << '\n';
    }
    std::erase_if(frontier, [&](const UnifiedProxyChoice& choice)
    {
        return !choice.final_coverage_passed ||
               choice.open_maximum_distance >
                   maximum_open_error_distance + tolerance;
    });
    if (frontier.empty())
        return makeExactFallback();
    stats.unified_final_frontier_size = frontier.size();
    stats.unified_final_frontier_workloads.clear();
    for (const UnifiedProxyChoice& choice : frontier)
        stats.unified_final_frontier_workloads.push_back(choice.workload);
    stats.unified_final_analytic_choices = static_cast<std::size_t>(std::count_if(
        frontier.begin(), frontier.end(), [](const UnifiedProxyChoice& choice)
        { return choice.analytic_responsibility_faces > 0; }));
    const auto staticWorkload = [](const UnifiedProxyChoice& choice)
    {
        // One overlapping pair of independently enclosed shells can survive a
        // binary BVTT node test and expand into four child pairs.  Use that
        // minimum expansion as a pose-independent PQSS workload surrogate;
        // the continuous overlap score remains telemetry for validation.
        return choice.workload +
            4.0 * static_cast<double>(
                choice.final_enclosure_overlap_pairs);
    };
    const auto betterChoice = [&](const UnifiedProxyChoice& first,
                                  const UnifiedProxyChoice& second)
    {
        // The user-provided maximum distance is a hard feasibility condition.
        // Among feasible candidates, minimize the static PQSS workload rather
        // than triangle count alone: interpenetrating enclosure shells keep
        // both RSS subtrees alive and multiply BVTT work.
        const double first_static_workload = staticWorkload(first);
        const double second_static_workload = staticWorkload(second);
        if (first_static_workload != second_static_workload)
            return first_static_workload < second_static_workload;
        if (first.workload != second.workload)
            return first.workload < second.workload;
        if (first.enclosure_responsibility_area !=
            second.enclosure_responsibility_area)
            return first.enclosure_responsibility_area >
                   second.enclosure_responsibility_area;
        return first.analytic_responsibility_area > second.analytic_responsibility_area;
    };
    const auto selected = std::min_element(
        frontier.begin(), frontier.end(), betterChoice);
    const auto best_analytic = std::min_element(
        frontier.begin(), frontier.end(), [&](const auto& first, const auto& second)
        {
            const bool first_analytic = first.analytic_responsibility_faces > 0;
            const bool second_analytic = second.analytic_responsibility_faces > 0;
            if (first_analytic != second_analytic) return first_analytic;
            if (!first_analytic) return false;
            return betterChoice(first, second);
        });
    stats.unified_selected_workload = selected->workload;
    stats.unified_selected_static_workload = staticWorkload(*selected);
    stats.unified_selected_enclosure_overlap_pairs =
        selected->final_enclosure_overlap_pairs;
    stats.unified_selected_enclosure_overlap_score =
        selected->final_enclosure_overlap_score;
    stats.unified_selected_analytic_faces =
        selected->analytic_responsibility_faces;
    stats.unified_selected_analytic_area_ratio =
        selected->analytic_responsibility_area /
        std::max(model_surface_area, 1.0e-30);
    if (best_analytic != frontier.end() &&
        best_analytic->analytic_responsibility_faces > 0)
    {
        stats.unified_best_analytic_workload = best_analytic->workload;
        stats.unified_best_analytic_faces =
            best_analytic->analytic_responsibility_faces;
        stats.unified_best_analytic_area_ratio =
            best_analytic->analytic_responsibility_area /
            std::max(model_surface_area, 1.0e-30);
    }
    // Export exactly the post-union geometry that competes in the final
    // workload/error audit. The raw choice is retained separately as its
    // conservative coverage certificate.
    std::vector<OutputPrimitive> output;
    flattenUnifiedChoice(selected->node, output);
    FinalCandidateEvaluation selected_evaluation = evaluateFinalUnionBoundary(
        std::move(output), 0, selected->workload,
        selected->preserve_pre_removal_enclosures, false);
    if (!selected->coverage_patch_faces.empty())
    {
        FinalCoverageAudit selected_patch;
        selected_patch.failed_face_ids = selected->coverage_patch_faces;
        selected_evaluation = patchCoverageFailures(
            std::move(selected_evaluation), selected_patch);
    }
    std::size_t final_patch_faces = 0;
    for (;;)
    {
        const FinalCoverageAudit replay_coverage =
            auditFinalConservativeCoverage(
                mesh, selected_evaluation.output, tolerance,
                &selected_evaluation.coverage_certificates,
                nullptr, excluded_redundant_faces);
        if (replay_coverage.unassigned_source_faces == 0 &&
            replay_coverage.failed_source_faces == 0)
            break;
        if (replay_coverage.failed_face_ids.empty() ||
            final_patch_faces + replay_coverage.failed_face_ids.size() > 128)
            break;
        final_patch_faces += replay_coverage.failed_face_ids.size();
        selected_evaluation = patchCoverageFailures(
            std::move(selected_evaluation), replay_coverage);
    }
    if (selected_coverage_certificates != nullptr)
        *selected_coverage_certificates =
            std::move(selected_evaluation.coverage_certificates);
    stats.unified_selected_workload = selected_evaluation.workload;
    return std::move(selected_evaluation.output);
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
            : model_diagonal * 0.08;
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
    case Kind::Frustum: return "frustum";
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
    for (std::size_t id = 0; id < primitives.size(); ++id)
    {
        const PrimitiveMesh mesh = triangulatePrimitive(primitives[id].primitive);
        stream << "g primitive_" << std::setw(5) << std::setfill('0') << id << '_'
               << kindName(primitives[id].primitive.kind) << std::setfill(' ') << '\n';
        for (const Vec3& vertex : mesh.vertices)
            stream << "v " << vertex.x() << ' ' << vertex.y() << ' ' << vertex.z() << '\n';
        for (const Face& face : mesh.faces)
            stream << "f " << offset + face[0] + 1 << ' ' << offset + face[1] + 1 << ' '
                   << offset + face[2] + 1 << '\n';
        offset += static_cast<std::uint32_t>(mesh.vertices.size());
        triangle_count += mesh.faces.size();
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
          << ",\"frustum\":" << stats.frustum_count
          << ",\"disk\":" << stats.disk_count
          << ",\"annulus\":" << stats.annulus_count
          << ",\"cylindricalband\":" << stats.cylindrical_band_count
          << ",\"conicalband\":" << stats.conical_band_count
          << "},\"triangulated_proxy_triangles\":" << stats.proxy_triangles
          << ",\"proxy_triangles\":" << stats.proxy_triangles
          << ",\"filled_planar_holes\":" << stats.filled_planar_holes
          << ",\"filled_cavity_volume_ratio\":" << stats.filled_cavity_volume_ratio
          << ",\"filled_boundary_voids\":" << stats.filled_boundary_voids
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
           << ",\"open_region_error\":{\"distance_method\":"
              "\"sampled_closest_source_open_surface_directed_distance\""
           << ",\"sampling_method\":"
              "\"deterministic_area_surface_with_vertices_and_edges\""
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
          << ",\"optimizer_profile\":{"
          << "\"approximate_planar_regions\":"
          << stats.unified_approximate_planar_regions
          << ",\"exact_coplanar_regions\":"
          << stats.unified_exact_coplanar_regions
          << ",\"hierarchy_regions\":" << stats.unified_hierarchy_regions
          << ",\"connected_components\":"
          << stats.unified_connected_components
          << ",\"lateral_sweep_axes\":"
          << stats.unified_lateral_sweep_axes
          << ",\"lateral_sweep_patches\":"
          << stats.unified_lateral_sweep_patches
          << ",\"lateral_sweep_certified\":"
          << stats.unified_lateral_sweep_certified
          << ",\"lateral_sweep_candidates\":"
          << stats.unified_lateral_sweep_candidates
          << ",\"lateral_sweep_max_faces\":"
          << stats.unified_lateral_sweep_max_faces
          << ",\"lateral_significant_candidates\":"
          << stats.unified_lateral_significant_candidates
          << ",\"lateral_sweep_max_axial_aspect\":"
          << stats.unified_lateral_sweep_max_axial_aspect
          << ",\"lateral_sweep_max_area_ratio\":"
          << stats.unified_lateral_sweep_max_area_ratio
          << ",\"lateral_reject_geometry\":"
          << stats.unified_lateral_reject_geometry
          << ",\"lateral_reject_end_ring\":"
          << stats.unified_lateral_reject_end_ring
          << ",\"lateral_reject_normal\":"
          << stats.unified_lateral_reject_normal
          << ",\"lateral_reject_radial\":"
          << stats.unified_lateral_reject_radial
          << ",\"lateral_best_end_ring_min_points\":"
          << stats.unified_lateral_best_end_ring_min_points
          << ",\"lateral_best_end_ring_spread_ratio\":"
          << stats.unified_lateral_best_end_ring_spread_ratio
          << ",\"lateral_best_end_ring_gap_ratio\":"
          << stats.unified_lateral_best_end_ring_gap_ratio
          << ",\"lateral_best_end_ring_score\":"
          << stats.unified_lateral_best_end_ring_score
          << ",\"final_frontier_size\":"
          << stats.unified_final_frontier_size
          << ",\"exact_fallback_selected\":"
          << (stats.unified_exact_fallback_selected ? "true" : "false")
          << ",\"final_analytic_choices\":"
          << stats.unified_final_analytic_choices
          << ",\"selected_workload\":"
          << stats.unified_selected_workload
          << ",\"selected_static_workload\":"
          << stats.unified_selected_static_workload
          << ",\"selected_enclosure_overlap_pairs\":"
          << stats.unified_selected_enclosure_overlap_pairs
          << ",\"selected_enclosure_overlap_score\":"
          << stats.unified_selected_enclosure_overlap_score
           << ",\"selected_analytic_faces\":"
          << stats.unified_selected_analytic_faces
          << ",\"selected_analytic_area_ratio\":"
          << stats.unified_selected_analytic_area_ratio
           << ",\"best_analytic_workload\":"
           << stats.unified_best_analytic_workload
           << ",\"best_analytic_faces\":"
          << stats.unified_best_analytic_faces
          << ",\"best_analytic_area_ratio\":"
          << stats.unified_best_analytic_area_ratio
           << ",\"whole_model_box_workload\":"
           << stats.unified_whole_model_box_workload
           << ",\"whole_model_box_passed_local_budget\":"
          << (stats.unified_whole_model_box_passed_local_budget
                  ? "true" : "false")
          << ",\"adaptive_refinements\":"
          << stats.unified_adaptive_refinements
          << ",\"skipped_final_audits\":"
          << stats.unified_skipped_final_audits
          << ",\"skipped_industrial_fallback_candidates\":"
          << stats.unified_skipped_industrial_fallback_candidates
          << ",\"final_frontier_choices\":[";
    for (std::size_t index = 0;
         index < stats.unified_final_frontier_workloads.size(); ++index)
    {
        if (index != 0) model << ',';
        model << "{\"workload\":"
              << stats.unified_final_frontier_workloads[index] << '}';
    }
    model << ']'
          << ",\"enclosure_occlusion\":{\"certified_extrusions\":"
          << stats.unified_enclosure_extrusions
          << ",\"clipped_primitives\":"
          << stats.unified_occlusion_clipped_primitives
          << ",\"removed_primitives\":"
          << stats.unified_occlusion_removed_primitives
          << ",\"input_triangles\":"
          << stats.unified_occlusion_input_triangles
          << ",\"output_triangles\":"
          << stats.unified_occlusion_output_triangles
          << ",\"removed_area\":"
          << stats.unified_occlusion_removed_area
          << ",\"rolled_back\":"
          << (stats.unified_occlusion_rolled_back ? "true" : "false")
          << ",\"seconds\":" << stats.unified_occlusion_seconds << '}'
           << ",\"classify_calls\":" << stats.unified_classify_calls
           << ",\"classify_faces\":" << stats.unified_classify_faces
           << ",\"frontier_prune_calls\":"
          << stats.unified_frontier_prune_calls
          << ",\"frontier_prune_choices\":"
          << stats.unified_frontier_prune_choices
          << ",\"frontier_prune_max_choices\":"
          << stats.unified_frontier_prune_max_choices
          << ",\"frontier_prune_comparisons\":"
          << stats.unified_frontier_prune_comparisons
           << ",\"timings_seconds\":{"
           << "\"planar_regions\":" << stats.unified_planar_region_seconds
           << ",\"connected_components\":"
          << stats.unified_connected_component_seconds
           << ",\"analytic_recognition\":"
           << stats.unified_analytic_recognition_seconds
           << ",\"classify\":" << stats.unified_classify_seconds
           << ",\"frontier_prune\":"
          << stats.unified_frontier_prune_seconds
          << ",\"hierarchy\":" << stats.unified_hierarchy_seconds
          << "}}"
          << ",\"uniform_structure_policy\":"
          << (options.uniform_structure_policy ? "true" : "false")
          << ",\"timings_seconds\":{\"total\":"
          << stats.analysis_seconds << "}},\"source\":\"source.obj\""
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
    PrimitiveMeshAnalysisOptions effective_options = options;

    if (!effective_options.allow_polygon)
        throw std::invalid_argument("polygon must be enabled for exact fallback");
    if (effective_options.frustum_segments < 3)
        throw std::invalid_argument("frustum_segments must be at least 3");
    if (effective_options.analysis_strength < 0.0 || effective_options.analysis_strength > 1.0)
        throw std::invalid_argument("analysis_strength must be in [0, 1]");

    // Production generation is staged: repair openings, recognize complete
    // surfaces, simplify the cleaned surface graph, then triangulate.  The
    // former unified frontier mixed geometry reconstruction with a PQSS
    // workload objective and is retained only as implementation history.
    if (effective_options.use_staged_surface_pipeline)
        effective_options.use_unified_candidate_optimizer = false;

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
            : diagonal * 0.08;
    stats.maximum_open_error_distance_limit = maximum_open_error_distance;

    std::vector<OutputPrimitive> envelope_output;
    double envelope_added_ratio = std::numeric_limits<double>::infinity();
    int envelope_axis = -1;
    std::size_t envelope_remaining_cavities = std::numeric_limits<std::size_t>::max();
    std::array<std::vector<OutputPrimitive>, 3> envelope_candidates;
    if (effective_options.uniform_structure_policy &&
        effective_options.enable_volume_evaluated_envelope)
    {
        for (int axis = 0; axis < 3; ++axis)
        {
            double added_ratio = 0.0;
            std::size_t remaining_cavities = 0;
            auto candidate = buildAxisEnvelopeCandidate(
                mesh, effective_options, lower, upper, axis, added_ratio,
                remaining_cavities);
            stats.envelope_candidate_added_volume_ratios[axis] = added_ratio;
            stats.envelope_candidate_primitive_counts[axis] = candidate.size();
            stats.envelope_candidate_remaining_cavities[axis] = remaining_cavities;
            envelope_candidates[axis] = candidate;
            if (candidate.empty() ||
                added_ratio > effective_options.maximum_envelope_added_volume_ratio)
                continue;
            if (envelope_output.empty() ||
                remaining_cavities < envelope_remaining_cavities ||
                (remaining_cavities == envelope_remaining_cavities &&
                 candidate.size() < envelope_output.size()) ||
                (remaining_cavities == envelope_remaining_cavities &&
                 candidate.size() == envelope_output.size() &&
                 added_ratio < envelope_added_ratio))
            {
                envelope_output = candidate;
                envelope_added_ratio = added_ratio;
                envelope_axis = axis;
                envelope_remaining_cavities = remaining_cavities;
            }
        }

        // Once an enclosing approximation has already spent a volume-error
        // budget, evaluate a coarser envelope by its incremental cost.  Comparing
        // every candidate only with the source rejected useful second-stage
        // merges: for example, replacing dozens of raster stair faces by six box
        // faces can add less volume than the cavity fill already accepted.  The
        // rule is geometry-only and applies to every OBJ; it also requires a
        // substantial primitive-count reduction and cannot hide more cavities.
        if (!envelope_output.empty())
        {
            for (int axis = 0; axis < 3; ++axis)
            {
                const auto& candidate = envelope_candidates[axis];
                const double candidate_added =
                    stats.envelope_candidate_added_volume_ratios[axis];
                const std::size_t candidate_cavities =
                    stats.envelope_candidate_remaining_cavities[axis];
                if (candidate.empty() || candidate.size() >= envelope_output.size() ||
                    candidate_cavities > envelope_remaining_cavities ||
                    envelope_output.size() <
                        effective_options.minimum_envelope_primitive_reduction_ratio *
                            candidate.size())
                    continue;
                const double incremental_added =
                    std::max(candidate_added - envelope_added_ratio, 0.0);
                if (incremental_added > envelope_added_ratio) continue;
                envelope_output = candidate;
                envelope_added_ratio = candidate_added;
                envelope_axis = axis;
                envelope_remaining_cavities = candidate_cavities;
            }
        }
    }
    markStage("envelope_candidates");

    if (effective_options.use_unified_candidate_optimizer)
    {
        const StructuralCleanup structural_cleanup =
            identifyStructuralRedundantFaces(
                mesh, effective_options, diagonal, model_surface_area,
                model_volume);
        if (effective_options.write_unified_region_diagnostics)
        {
            std::ofstream cleanup_faces(
                output_directory / "structural_cleanup_faces.csv");
            cleanup_faces << "category,face\n";
            const auto write_faces = [&](const char* category,
                                         const auto& faces)
            {
                for (const auto face : faces)
                    cleanup_faces << category << ',' << face << '\n';
            };
            write_faces("sealed_void_wall",
                        structural_cleanup.sealed_void_wall_faces);
            write_faces("blind_cavity",
                        structural_cleanup.blind_cavity_faces);
            write_faces("contained", structural_cleanup.contained_faces);
            write_faces("restored_cavity",
                        structural_cleanup.restored_cavity_faces);
        }
        std::unordered_set<std::uint32_t> restored(
            structural_cleanup.restored_cavity_faces.begin(),
            structural_cleanup.restored_cavity_faces.end());
        const std::unordered_set<std::uint32_t> sealed_faces(
            structural_cleanup.sealed_void_wall_faces.begin(),
            structural_cleanup.sealed_void_wall_faces.end());
        const std::unordered_set<std::uint32_t> contained_faces(
            structural_cleanup.contained_faces.begin(),
            structural_cleanup.contained_faces.end());
        for (const auto& group :
             structural_cleanup.restored_cavity_face_groups)
        {
            const bool contains_support_surface = std::any_of(
                group.begin(), group.end(), [&](const auto face)
                {
                    return contained_faces.contains(face);
                });
            if (!contains_support_surface) continue;
            // A mixed group can have a large AABB solely because long exterior
            // support surfaces meet a tiny boundary bevel. The support surfaces
            // still need restoration, but a wall independently certified between
            // the filled cap planes remains redundant; restoring it would reopen
            // the already accepted hole boundary.
            for (const auto face : group)
                if (sealed_faces.contains(face)) restored.erase(face);
        }
        std::vector<std::uint32_t> redundant_faces;
        redundant_faces.reserve(
            structural_cleanup.sealed_void_wall_faces.size() +
            structural_cleanup.blind_cavity_faces.size());
        redundant_faces.insert(
            redundant_faces.end(),
            structural_cleanup.sealed_void_wall_faces.begin(),
            structural_cleanup.sealed_void_wall_faces.end());
        redundant_faces.insert(
            redundant_faces.end(),
            structural_cleanup.blind_cavity_faces.begin(),
            structural_cleanup.blind_cavity_faces.end());
        std::erase_if(redundant_faces, [&](const auto face)
        {
            return restored.contains(face);
        });
        std::sort(redundant_faces.begin(), redundant_faces.end());
        redundant_faces.erase(
            std::unique(redundant_faces.begin(), redundant_faces.end()),
            redundant_faces.end());
        stats.excluded_redundant_triangles = redundant_faces.size();
        stats.removed_sealed_void_wall_primitives =
            structural_cleanup.sealed_void_wall_primitives;
        stats.excluded_sealed_void_wall_triangles =
            structural_cleanup.sealed_void_wall_triangles;
        stats.removed_blind_cavity_primitives =
            structural_cleanup.blind_cavity_primitives;
        stats.removed_contained_primitives = 0;
        markStage("structural_cleanup");
        std::vector<OutputPrimitive> enclosure_certificates;
        std::vector<OutputPrimitive> unified = optimizeUnifiedProxyCandidates(
            mesh, effective_options, diagonal, model_volume, model_surface_area,
            stats, output_directory / "unified_performance.txt",
            &redundant_faces, &enclosure_certificates);
        const double unified_tolerance = std::max(
            diagonal * 1.0e-9, 1.0e-10);
        std::vector<CertifiedExtrusion> enclosure_extrusions;
        if (stats.unified_exact_fallback_selected)
        {
            // Every fallback primitive is one unchanged source triangle with a
            // direct responsibility certificate.  Global canonicalization is
            // not needed for correctness here, and its all-pairs planar union
            // is precisely the unsafe path when an industrial mesh has tens of
            // thousands of fallback triangles.
            promoteToSemanticPrimitives(unified);
        }
        else
        {
            // The optimizer already returned the post-union boundary that won
            // its final audit. Keep the raw owners only as certificates; a
            // caller-side geometry pass would invalidate the selected error
            // and workload.
            promoteToSemanticPrimitives(enclosure_certificates);
            promoteToSemanticPrimitives(unified);
            enclosure_extrusions = recognizeEnclosureGroupExtrusions(
                enclosure_certificates, unified_tolerance, false);
            stats.unified_enclosure_extrusions = enclosure_extrusions.size();
        }
        const auto coverage_started = std::chrono::steady_clock::now();
        const FinalCoverageAudit coverage = auditFinalConservativeCoverage(
            mesh, unified, std::max(diagonal * 1.0e-9, 1.0e-10),
            enclosure_certificates.empty() ? nullptr : &enclosure_certificates,
            nullptr,
            &redundant_faces);
        stats.coverage_assigned_source_faces = coverage.assigned_source_faces;
        stats.coverage_enclosure_source_faces = coverage.enclosure_source_faces;
        stats.coverage_planar_source_faces = coverage.planar_source_faces;
        stats.coverage_unassigned_source_faces = coverage.unassigned_source_faces;
        stats.coverage_failed_source_faces = coverage.failed_source_faces;
        stats.containment_validation_passed =
            coverage.unassigned_source_faces == 0 &&
            coverage.failed_source_faces == 0;
        stats.coverage_audit_seconds = std::chrono::duration<double>(
            std::chrono::steady_clock::now() - coverage_started).count();
        if (coverage.unassigned_source_faces != 0 ||
            coverage.failed_source_faces != 0)
        {
            writeCoverageFailureDiagnostics(
                output_directory, mesh, unified, coverage);
            std::ostringstream message;
            message << "final conservative coverage audit failed: unassigned="
                    << coverage.unassigned_source_faces << ", uncovered="
                    << coverage.failed_source_faces
                    << ", failed_with_enclosure_owner="
                    << coverage.failed_with_enclosure_owner
                    << ", failed_with_planar_owner="
                    << coverage.failed_with_planar_owner;
            throw std::runtime_error(message.str());
        }
        markStage("unified_candidate_optimization");

        std::filesystem::create_directories(output_directory);
        writeSourceObj(output_directory / "source.obj", mesh);
        writeRegionsObj(
            output_directory / "regions.obj", mesh, unified,
            redundant_faces);
        stats.primitive_count = unified.size();
        for (const auto& item : unified)
        {
            switch (item.primitive.kind)
            {
            case Kind::Polygon: ++stats.polygon_count; break;
            case Kind::Frustum: ++stats.frustum_count; break;
            case Kind::Disk: ++stats.disk_count; break;
            case Kind::Annulus: ++stats.annulus_count; break;
            case Kind::CylindricalBand: ++stats.cylindrical_band_count; break;
            case Kind::ConicalBand: ++stats.conical_band_count; break;
            case Kind::Rectangle:
            case Kind::Triangle:
                throw std::logic_error("unified optimizer emitted an internal primitive");
            }
        }
        writeSemanticPrimitiveObj(output_directory / "primitives.obj", unified);
        stats.proxy_triangles = writeTriangulatedObj(
            output_directory / "proxy.obj", unified);
        const auto error_started = std::chrono::steady_clock::now();
        const SourceOpenSurfaceReference final_open_surface_reference(
            mesh, diagonal);
        {
            std::ofstream openness(output_directory / "region_openness.csv");
            openness << "primitive,source_faces,open_source_faces\n";
            for (std::size_t primitive = 0; primitive < unified.size();
                 ++primitive)
            {
                std::unordered_set<std::uint32_t> faces(
                    unified[primitive].source_faces.begin(),
                    unified[primitive].source_faces.end());
                const std::size_t open_faces = static_cast<std::size_t>(
                    std::count_if(faces.begin(), faces.end(),
                        [&](const auto face)
                        {
                            return face < final_open_surface_reference.open_faces.size() &&
                                final_open_surface_reference.open_faces[face] != 0;
                        }));
                openness << primitive << ',' << faces.size() << ','
                         << open_faces << '\n';
            }
        }
        const FinalOpenErrorAudit error = measureChargeableOpenSurfaceDistance(
            final_open_surface_reference, mesh, unified,
            std::max(diagonal / 192.0, 1.0e-30), unified_tolerance,
            enclosure_certificates.empty() ? nullptr : &enclosure_certificates);
        stats.open_error_distance_sample_count = error.distance_sample_count;
        stats.open_mean_distance = error.mean_distance;
        stats.open_max_distance = error.maximum_distance;
        stats.open_mean_distance_ratio = error.mean_distance /
            std::max(diagonal, 1.0e-30);
        stats.open_max_distance_ratio = error.maximum_distance /
            std::max(diagonal, 1.0e-30);
        stats.open_max_proxy_point = {{error.maximum_proxy_point.x(),
                                       error.maximum_proxy_point.y(),
                                       error.maximum_proxy_point.z()}};
        stats.open_max_source_point = {{error.maximum_source_point.x(),
                                        error.maximum_source_point.y(),
                                        error.maximum_source_point.z()}};
        stats.open_error_audit_seconds = std::chrono::duration<double>(
            std::chrono::steady_clock::now() - error_started).count();
        if (error.maximum_distance >
            stats.maximum_open_error_distance_limit + unified_tolerance)
        {
            std::ostringstream message;
            message << "final proxy exceeds maximum open error distance: maximum="
                    << error.maximum_distance << ", limit="
                    << stats.maximum_open_error_distance_limit;
            throw std::runtime_error(message.str());
        }
        writeOpenErrorVisualization(
            output_directory / "open_error.json", error);
        stats.analysis_seconds = std::chrono::duration<double>(
            std::chrono::steady_clock::now() - started).count();
        writeMetadata(output_directory, input_obj, unified, stats, effective_options);
        return stats;
    }

    const StructuralCleanup structural_cleanup = identifyStructuralRedundantFaces(
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
    std::uint64_t next_enclosure_group = 1;
    double analytic_filled_cavity_volume = 0.0;
    if (effective_options.allow_frustum)
    {
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
            if (surface.empty()) continue;
            output.insert(output.end(), std::make_move_iterator(surface.begin()),
                          std::make_move_iterator(surface.end()));
            for (const auto face : component) responsibility_faces[face] = false;
        }
    }
    markStage("analytic_surface_recognition");
    stats.filled_cavity_volume_ratio = analytic_filled_cavity_volume /
        std::max(model_volume, 1.0e-30);
    if (effective_options.allow_frustum)
    {
        // A round volume must be recognized from a complete connected surface.
        // Greedy pairwise merging cannot certify a circle while it only sees one
        // polygonal side strip, and accepting those partial fits caused unrelated
        // strips to bridge a large empty region.
        for (const auto& component : faceComponents(mesh, responsibility_faces))
        {
            const double analytic_tolerance = diagonal *
                effective_options.analytic_surface_relative_tolerance;
            VolumeFit fit = fitBestVolume(
                mesh, component, effective_options, analytic_tolerance);
            if (!fit.use_frustum) continue;
            appendConicalSurfaceAssembly(
                output, mesh, fit.frustum, component, analytic_tolerance);
            for (const auto face : component) responsibility_faces[face] = false;
        }
    }
    std::vector<std::unordered_set<std::uint32_t>> adjacency;
        auto clusters = coplanarClusters(mesh, diagonal * effective_options.coplanar_relative_tolerance,
                                         adjacency, &responsibility_faces);
        std::vector<bool> claimed_clusters(clusters.size(), false);
        // Component boxes and support-protrusion boxes belonged to the old
        // volume-first hierarchy.  Surface recognition must finish before any
        // error-bounded simplification, so planar clusters stay as surface
        // regions here.
        std::vector<RecognizedProtrusion> protrusions;
        std::vector<Vec3> cluster_centers(clusters.size(), Vec3::Zero());
        std::vector<Vec3> cluster_lowers(clusters.size());
        std::vector<Vec3> cluster_uppers(clusters.size());
        for (std::size_t id = 0; id < clusters.size(); ++id)
        {
            Vec3 cluster_lower = Vec3::Constant(std::numeric_limits<double>::infinity());
            Vec3 cluster_upper = Vec3::Constant(-std::numeric_limits<double>::infinity());
            for (const auto vertex : uniqueVertices(mesh, clusters[id]))
            {
                cluster_lower = cluster_lower.cwiseMin(mesh.vertices[vertex]);
                cluster_upper = cluster_upper.cwiseMax(mesh.vertices[vertex]);
            }
            cluster_lowers[id] = cluster_lower;
            cluster_uppers[id] = cluster_upper;
            cluster_centers[id] = (cluster_lower + cluster_upper) * 0.5;
        }
        std::vector<std::uint32_t> spatial_order(clusters.size());
        std::iota(spatial_order.begin(), spatial_order.end(), 0);
        for (int axis = 0; axis < 3; ++axis)
        {
            std::sort(spatial_order.begin(), spatial_order.end(), [&](const auto first,
                                                                      const auto second)
            {
                if (cluster_centers[first][axis] != cluster_centers[second][axis])
                    return cluster_centers[first][axis] < cluster_centers[second][axis];
                const int second_axis = (axis + 1) % 3;
                return cluster_centers[first][second_axis] < cluster_centers[second][second_axis];
            });
            for (std::size_t index = 1; index < spatial_order.size(); ++index)
            {
                const auto first = spatial_order[index - 1];
                const auto second = spatial_order[index];
                adjacency[first].insert(second);
                adjacency[second].insert(first);
            }
        }
        std::vector<Node> nodes;
        nodes.reserve(2 * clusters.size() + protrusions.size());
        std::vector<std::uint32_t> replacement(clusters.size());
        std::iota(replacement.begin(), replacement.end(), 0);
        for (std::size_t index = 0; index < protrusions.size(); ++index)
            for (const auto cluster_id : protrusions[index].clusters)
                replacement[cluster_id] = static_cast<std::uint32_t>(clusters.size() + index);
        for (std::size_t id = 0; id < clusters.size(); ++id)
        {
            Node node;
            node.faces = clusters[id];
            node.fit = fitBestVolume(mesh, clusters[id], effective_options);
            node.lower = cluster_lowers[id];
            node.upper = cluster_uppers[id];
            node.source_area = facesArea(mesh, clusters[id]);
            node.active = !claimed_clusters[id];
            nodes.push_back(std::move(node));
        }
        for (auto& protrusion : protrusions)
        {
            Node node;
            node.faces = std::move(protrusion.faces);
            node.fit = boxVolumeFit(protrusion.box);
            const auto [node_lower, node_upper] = faceBounds(mesh, node.faces);
            node.lower = node_lower;
            node.upper = node_upper;
            node.source_area = facesArea(mesh, node.faces);
            node.recognized_protrusion = true;
            node.covered_box_face_axis = protrusion.covered_face_axis;
            node.covered_box_face_sign = protrusion.covered_face_sign;
            nodes.push_back(std::move(node));
        }
        for (std::uint32_t id = 0; id < clusters.size(); ++id)
        {
            if (!nodes[id].active) continue;
            for (const auto neighbor : adjacency[id])
            {
                const auto mapped = replacement[neighbor];
                if (mapped != id) nodes[id].neighbors.insert(mapped);
            }
        }
        for (std::size_t index = 0; index < protrusions.size(); ++index)
        {
            const auto node_id = static_cast<std::uint32_t>(clusters.size() + index);
            for (const auto cluster_id : protrusions[index].clusters)
            {
                for (const auto neighbor : adjacency[cluster_id])
                {
                    const auto mapped = replacement[neighbor];
                    if (mapped != node_id) nodes[node_id].neighbors.insert(mapped);
                }
            }
        }
        std::priority_queue<Candidate, std::vector<Candidate>, CandidateGreater> heap;
        const auto push = [&](const std::uint32_t first, const std::uint32_t second)
        {
            // The later surface-primitive pass owns all approximate merging.
            // Do not fit arbitrary source-region pairs to a closed box/frustum.
            return;
            if (first == second || !nodes[first].active || !nodes[second].active) return;
            // A recognized independent assembly is already a semantic leaf.
            // Letting the generic region merger absorb it destroys the component
            // boundary that recognition was meant to preserve.
            if (nodes[first].recognized_protrusion || nodes[second].recognized_protrusion)
                return;
            std::vector<std::uint32_t> merged = nodes[first].faces;
            merged.insert(merged.end(), nodes[second].faces.begin(), nodes[second].faces.end());
            VolumeFit fit = fitBestVolume(mesh, merged, effective_options);
            if (effective_options.uniform_structure_policy &&
                effective_options.forbid_main_body_box_approximation &&
                !fit.use_frustum) return;
            const double added = fit.volume - nodes[first].fit.volume - nodes[second].fit.volume;
            const double source_area = nodes[first].source_area + nodes[second].source_area;
            const double area_excess_ratio = std::max(
                fit.surface_area - source_area, 0.0) /
                std::max(model_surface_area, 1.0e-30);
            const double volume_excess_ratio = std::max(added, 0.0) / model_volume;
            const double separation_ratio = boundsSeparationRatio(
                nodes[first].lower, nodes[first].upper,
                nodes[second].lower, nodes[second].upper);
            const double object_boundary_ratio = 0.0;
            const double score = std::max({area_excess_ratio, volume_excess_ratio,
                                           separation_ratio, object_boundary_ratio});
            heap.push({score, fit.volume, first, second,
                       nodes[first].version, nodes[second].version, fit});
        };
        for (std::uint32_t first = 0; first < nodes.size(); ++first)
            for (const auto second : nodes[first].neighbors)
                if (first < second) push(first, second);

        while (!heap.empty())
        {
            Candidate candidate = heap.top();
            heap.pop();
            if (candidate.added_ratio > threshold) break;
            if (candidate.first >= nodes.size() || candidate.second >= nodes.size()) continue;
            Node& first = nodes[candidate.first];
            Node& second = nodes[candidate.second];
            if (!first.active || !second.active || first.version != candidate.first_version ||
                second.version != candidate.second_version ||
                !first.neighbors.contains(candidate.second)) continue;
            Node merged;
            merged.faces = first.faces;
            merged.faces.insert(merged.faces.end(), second.faces.begin(), second.faces.end());
            merged.fit = candidate.fit;
            merged.lower = first.lower.cwiseMin(second.lower);
            merged.upper = first.upper.cwiseMax(second.upper);
            merged.source_area = first.source_area + second.source_area;
            merged.neighbors = first.neighbors;
            merged.neighbors.insert(second.neighbors.begin(), second.neighbors.end());
            merged.neighbors.erase(candidate.first);
            merged.neighbors.erase(candidate.second);
            first.active = false;
            second.active = false;
            ++first.version;
            ++second.version;
            const std::uint32_t merged_id = static_cast<std::uint32_t>(nodes.size());
            nodes.push_back(std::move(merged));
            for (const auto neighbor : nodes.back().neighbors)
            {
                nodes[neighbor].neighbors.erase(candidate.first);
                nodes[neighbor].neighbors.erase(candidate.second);
                nodes[neighbor].neighbors.insert(merged_id);
            }
            for (const auto neighbor : nodes.back().neighbors) push(merged_id, neighbor);
        }
    for (const Node& node : nodes)
    {
        if (!node.active) continue;
        if (node.recognized_protrusion)
        {
            appendBoxRectangles(output, node.fit.box, node.faces,
                                node.covered_box_face_axis,
                                node.covered_box_face_sign,
                                next_enclosure_group++);
            ++stats.recognized_protrusion_box_shells;
            continue;
        }
        auto classified = classifyFinalRegion(mesh, node.faces, effective_options, threshold,
                                              model_surface_area,
                                              stats.filled_planar_holes,
                                              stats.filled_boundary_voids,
                                              stats.filled_boundary_void_area);
        output.insert(output.end(), std::make_move_iterator(classified.begin()),
                      std::make_move_iterator(classified.end()));
    }
    markStage("region_classification");

    output = mergeLocalCoplanarPrimitives(
        mesh, std::move(output),
        std::max(diagonal * 1.0e-9, 1.0e-10),
        maximum_open_error_distance,
        std::max(diagonal / 192.0, 1.0e-30),
        stats.merged_local_planar_primitives);
    CoplanarCanonicalizationStats coplanar_stats;
    output = canonicalizeCoplanarPrimitiveUnion(
        mesh, std::move(output), std::max(diagonal * 1.0e-9, 1.0e-10),
        coplanar_stats);
    stats.canonicalized_coplanar_groups = coplanar_stats.groups;
    stats.removed_coplanar_redundant_primitives = coplanar_stats.removed_primitives;
    stats.removed_coplanar_overlap_area = coplanar_stats.removed_overlap_area;
    markStage("initial_coplanar_processing");

    const std::vector<std::uint32_t>& excluded_faces = structural_cleanup.excluded_faces;
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

    if (!structural_cleanup.restored_cavity_faces.empty())
    {
        // Restore only the interior surfaces of over-budget cavities.  They do
        // not participate in exterior box/protrusion merging, shallow-shell
        // coalescing, or silhouette synchronization.  The exterior pass above
        // therefore keeps a stable partition independent of whether a candidate
        // cavity is ultimately accepted or restored.
        std::vector<bool> cavity_mask(mesh.faces.size(), false);
        for (const auto face : structural_cleanup.restored_cavity_faces)
            cavity_mask[face] = true;
        std::vector<std::unordered_set<std::uint32_t>> cavity_adjacency;
        const auto cavity_clusters = coplanarClusters(
            mesh, diagonal * effective_options.coplanar_relative_tolerance,
            cavity_adjacency, &cavity_mask);
        std::vector<OutputPrimitive> cavity_output;
        for (const auto& cluster : cavity_clusters)
        {
            auto classified = classifyFinalRegion(
                mesh, cluster, effective_options, threshold,
                model_surface_area,
                stats.filled_planar_holes, stats.filled_boundary_voids,
                stats.filled_boundary_void_area);
            for (OutputPrimitive& item : classified)
                item.preserves_cavity_opening = true;
            cavity_output.insert(cavity_output.end(),
                std::make_move_iterator(classified.begin()),
                std::make_move_iterator(classified.end()));
        }
        CoplanarCanonicalizationStats cavity_coplanar_stats;
        cavity_output = canonicalizeCoplanarPrimitiveUnion(
            mesh, std::move(cavity_output),
            std::max(diagonal * 1.0e-9, 1.0e-10), cavity_coplanar_stats);
        promoteToSemanticPrimitives(cavity_output);
        output.insert(output.end(),
            std::make_move_iterator(cavity_output.begin()),
            std::make_move_iterator(cavity_output.end()));
    }
    output = reopenRestoredCavityVolumes(
        mesh, std::move(output), structural_cleanup,
        effective_options.maximum_cavity_added_volume_ratio, model_volume,
        std::max(diagonal * 1.0e-9, 1.0e-10));
    CoplanarCanonicalizationStats final_surface_coplanar_stats;
    output = canonicalizeCoplanarPrimitiveUnion(
        mesh, std::move(output), std::max(diagonal * 1.0e-9, 1.0e-10),
        final_surface_coplanar_stats);
    auto final_closed_extrusions = recognizeCertifiedPrismaticVolumes(
        output, std::max(diagonal * 1.0e-9, 1.0e-10));
    certified_closed_extrusions.insert(
        certified_closed_extrusions.end(),
        std::make_move_iterator(final_closed_extrusions.begin()),
        std::make_move_iterator(final_closed_extrusions.end()));
    std::size_t volume_occluded_primitives = 0;
    output = clipPlanarOcclusionByClosedVolumes(
        std::move(output), certified_closed_volumes, certified_closed_extrusions,
        std::max(diagonal * 1.0e-9, 1.0e-10), volume_occluded_primitives,
        nullptr, false);
    stats.removed_contained_primitives += volume_occluded_primitives;
    output = clipParallelOuterOcclusion(
        std::move(output), (lower + upper) * 0.5, diagonal * 0.03,
        std::max(diagonal * 1.0e-9, 1.0e-10));
    markStage("final_surface_canonicalization");

    const double final_tolerance = std::max(
        diagonal * 1.0e-9, 1.0e-10);
    FinalCoverageAudit coverage = auditFinalConservativeCoverage(
        mesh, output, final_tolerance, &output, nullptr, &excluded_faces);
    if (!coverage.failed_face_ids.empty())
    {
        for (const auto face_id : coverage.failed_face_ids)
        {
            if (face_id >= mesh.faces.size()) continue;
            Primitive triangle;
            triangle.kind = Kind::Triangle;
            for (int corner = 0; corner < 3; ++corner)
                triangle.triangle[corner] =
                    mesh.vertices[mesh.faces[face_id][corner]];
            output.push_back({std::move(triangle), {face_id}});
        }
        promoteToSemanticPrimitives(output);
        coverage = auditFinalConservativeCoverage(
            mesh, output, final_tolerance, &output, nullptr, &excluded_faces);
    }
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

    const SourceOpenSurfaceReference open_surface_reference(mesh, diagonal);
    const FinalOpenErrorAudit open_error = measureChargeableOpenSurfaceDistance(
        open_surface_reference, mesh, output,
        std::max(diagonal / 192.0, 1.0e-30), final_tolerance, &output);
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
    writeSourceObj(output_directory / "source.obj", mesh);
    writeRegionsObj(
        output_directory / "regions.obj", mesh, output, excluded_faces);
    stats.primitive_count = output.size();
    for (const auto& item : output)
    {
        switch (item.primitive.kind)
        {
        case Kind::Polygon: ++stats.polygon_count; break;
        case Kind::Frustum: ++stats.frustum_count; break;
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

std::vector<PrimitiveMeshAnalysisStats> analyzePrimitiveMeshObjStrengthSweep(
    const std::filesystem::path& input_obj,
    const std::filesystem::path& output_root,
    const PrimitiveMeshAnalysisOptions& options)
{
    const auto sweep_started = std::chrono::steady_clock::now();
    if (!options.allow_polygon)
        throw std::invalid_argument("polygon must be enabled for exact fallback");
    if (options.frustum_segments < 3)
        throw std::invalid_argument("frustum_segments must be at least 3");

    Mesh mesh = readObj(input_obj);
    weldCoincidentVertices(mesh);
    PrimitiveMeshAnalysisStats base_stats;
    base_stats.source_triangles = mesh.faces.size();
    base_stats.discarded_degenerate_triangles = dropDegenerateFaces(mesh);
    base_stats.source_triangles = mesh.faces.size();

    Vec3 lower = Vec3::Constant(std::numeric_limits<double>::infinity());
    Vec3 upper = Vec3::Constant(-std::numeric_limits<double>::infinity());
    for (const Vec3& vertex : mesh.vertices)
    {
        lower = lower.cwiseMin(vertex);
        upper = upper.cwiseMax(vertex);
    }
    const Vec3 extent = upper - lower;
    const double diagonal = extent.norm();
    const double model_volume = std::max(
        extent.prod(), diagonal * diagonal * diagonal * 1.0e-12);
    std::vector<std::uint32_t> all_faces(mesh.faces.size());
    std::iota(all_faces.begin(), all_faces.end(), 0);
    const double model_surface_area = facesArea(mesh, all_faces);
    const StructuralCleanup structural_cleanup = identifyStructuralRedundantFaces(
        mesh, options, diagonal, model_surface_area, model_volume);
    base_stats.removed_sealed_void_wall_primitives =
        structural_cleanup.sealed_void_wall_primitives;
    base_stats.excluded_sealed_void_wall_triangles =
        structural_cleanup.sealed_void_wall_triangles;
    base_stats.removed_blind_cavity_primitives =
        structural_cleanup.blind_cavity_primitives;
    base_stats.removed_contained_primitives = structural_cleanup.contained_primitives;
    base_stats.excluded_redundant_triangles = structural_cleanup.excluded_faces.size();
    std::vector<bool> responsibility_faces(mesh.faces.size(), true);
    for (const auto face : structural_cleanup.excluded_faces)
        responsibility_faces[face] = false;

    std::vector<std::unordered_set<std::uint32_t>> adjacency;
    auto clusters = coplanarClusters(
        mesh, diagonal * options.coplanar_relative_tolerance,
        adjacency, &responsibility_faces);
    std::vector<bool> claimed_clusters(clusters.size(), false);

    // Closed boxes are geometric identities, so recognize them once on the
    // finest fixed partition and keep them atomic at every strength.  Doing
    // this after each strength cut could make the primitive count decrease.
    std::vector<OutputPrimitive> finest_partition;
    std::size_t ignored_holes = 0;
    std::size_t ignored_boundary_voids = 0;
    double ignored_boundary_void_area = 0.0;
    const double finest_threshold = analysisStrengthToAllowedExcessRatio(1.0);
    for (const auto& cluster : clusters)
    {
        auto classified = classifyFinalRegion(
            mesh, cluster, options, finest_threshold, model_surface_area,
            ignored_holes,
            ignored_boundary_voids, ignored_boundary_void_area);
        finest_partition.insert(finest_partition.end(),
            std::make_move_iterator(classified.begin()),
            std::make_move_iterator(classified.end()));
    }
    std::vector<RecognizedProtrusion> closed_boxes;
    std::size_t ignored_closed_box_count = 0;
    recognizeClosedAxisAlignedBoxes(
        finest_partition, std::max(diagonal * 1.0e-9, 1.0e-10),
        ignored_closed_box_count, &closed_boxes);

    std::vector<RecognizedProtrusion> atomic_closed_boxes;
    for (auto& box : closed_boxes)
    {
        std::unordered_set<std::uint32_t> face_set(
            box.faces.begin(), box.faces.end());
        std::vector<std::uint32_t> box_clusters;
        std::vector<std::uint32_t> box_faces;
        bool overlaps_claimed = false;
        for (std::uint32_t cluster_id = 0; cluster_id < clusters.size(); ++cluster_id)
        {
            const bool intersects = std::any_of(
                clusters[cluster_id].begin(), clusters[cluster_id].end(),
                [&](const auto face) { return face_set.contains(face); });
            if (!intersects) continue;
            const bool contained = std::all_of(
                clusters[cluster_id].begin(), clusters[cluster_id].end(),
                [&](const auto face) { return face_set.contains(face); });
            if (!contained || claimed_clusters[cluster_id])
            {
                overlaps_claimed = true;
                break;
            }
            box_clusters.push_back(cluster_id);
            box_faces.insert(box_faces.end(), clusters[cluster_id].begin(),
                             clusters[cluster_id].end());
        }
        if (overlaps_claimed || box_faces.size() != face_set.size()) continue;
        for (const auto cluster_id : box_clusters) claimed_clusters[cluster_id] = true;
        atomic_closed_boxes.push_back(
            {box.box, std::move(box_faces), std::move(box_clusters)});
    }
    auto protrusions = recognizeSupportProtrusions(
        mesh, clusters, adjacency, options, diagonal, model_surface_area,
        claimed_clusters, base_stats.minimum_protrusion_candidate_area_excess_ratio);
    const std::size_t support_protrusion_count = protrusions.size();
    protrusions.insert(protrusions.end(),
        std::make_move_iterator(atomic_closed_boxes.begin()),
        std::make_move_iterator(atomic_closed_boxes.end()));

    std::vector<Vec3> cluster_centers(clusters.size(), Vec3::Zero());
    std::vector<Vec3> cluster_lowers(clusters.size());
    std::vector<Vec3> cluster_uppers(clusters.size());
    for (std::size_t id = 0; id < clusters.size(); ++id)
    {
        Vec3 cluster_lower = Vec3::Constant(std::numeric_limits<double>::infinity());
        Vec3 cluster_upper = Vec3::Constant(-std::numeric_limits<double>::infinity());
        for (const auto vertex : uniqueVertices(mesh, clusters[id]))
        {
            cluster_lower = cluster_lower.cwiseMin(mesh.vertices[vertex]);
            cluster_upper = cluster_upper.cwiseMax(mesh.vertices[vertex]);
        }
        cluster_lowers[id] = cluster_lower;
        cluster_uppers[id] = cluster_upper;
        cluster_centers[id] = (cluster_lower + cluster_upper) * 0.5;
    }
    std::vector<std::uint32_t> spatial_order(clusters.size());
    std::iota(spatial_order.begin(), spatial_order.end(), 0);
    for (int axis = 0; axis < 3; ++axis)
    {
        std::sort(spatial_order.begin(), spatial_order.end(), [&](const auto first,
                                                                  const auto second)
        {
            if (cluster_centers[first][axis] != cluster_centers[second][axis])
                return cluster_centers[first][axis] < cluster_centers[second][axis];
            const int second_axis = (axis + 1) % 3;
            return cluster_centers[first][second_axis] < cluster_centers[second][second_axis];
        });
        for (std::size_t index = 1; index < spatial_order.size(); ++index)
        {
            const auto first = spatial_order[index - 1];
            const auto second = spatial_order[index];
            adjacency[first].insert(second);
            adjacency[second].insert(first);
        }
    }

    std::vector<Node> nodes;
    nodes.reserve(2 * clusters.size() + protrusions.size());
    std::vector<std::uint32_t> replacement(clusters.size());
    std::iota(replacement.begin(), replacement.end(), 0);
    for (std::size_t index = 0; index < protrusions.size(); ++index)
        for (const auto cluster_id : protrusions[index].clusters)
            replacement[cluster_id] = static_cast<std::uint32_t>(clusters.size() + index);
    for (std::size_t id = 0; id < clusters.size(); ++id)
    {
        Node node;
        node.faces = clusters[id];
        node.fit = fitBestVolume(mesh, clusters[id], options);
        node.lower = cluster_lowers[id];
        node.upper = cluster_uppers[id];
        node.source_area = facesArea(mesh, clusters[id]);
        node.active = !claimed_clusters[id];
        nodes.push_back(std::move(node));
    }
    for (std::size_t index = 0; index < protrusions.size(); ++index)
    {
        auto& protrusion = protrusions[index];
        Node node;
        node.faces = std::move(protrusion.faces);
        node.fit = boxVolumeFit(protrusion.box);
        const auto [node_lower, node_upper] = faceBounds(mesh, node.faces);
        node.lower = node_lower;
        node.upper = node_upper;
        node.source_area = facesArea(mesh, node.faces);
        node.recognized_protrusion = index < support_protrusion_count;
        node.recognized_closed_box = index >= support_protrusion_count;
        node.covered_box_face_axis = protrusion.covered_face_axis;
        node.covered_box_face_sign = protrusion.covered_face_sign;
        nodes.push_back(std::move(node));
    }
    for (std::uint32_t id = 0; id < clusters.size(); ++id)
    {
        if (!nodes[id].active) continue;
        for (const auto neighbor : adjacency[id])
        {
            const auto mapped = replacement[neighbor];
            if (mapped != id) nodes[id].neighbors.insert(mapped);
        }
    }
    for (std::size_t index = 0; index < protrusions.size(); ++index)
    {
        const auto node_id = static_cast<std::uint32_t>(clusters.size() + index);
        for (const auto cluster_id : protrusions[index].clusters)
            for (const auto neighbor : adjacency[cluster_id])
            {
                const auto mapped = replacement[neighbor];
                if (mapped != node_id) nodes[node_id].neighbors.insert(mapped);
            }
    }

    std::size_t active_count = static_cast<std::size_t>(std::count_if(
        nodes.begin(), nodes.end(), [](const Node& node) { return node.active; }));
    while (active_count > 1)
    {
        std::vector<Candidate> candidates;
        for (std::uint32_t first = 0; first < nodes.size(); ++first)
        {
            if (!nodes[first].active) continue;
            for (const auto second : nodes[first].neighbors)
            {
                if (first >= second || second >= nodes.size() || !nodes[second].active) continue;
                std::vector<std::uint32_t> merged_faces = nodes[first].faces;
                merged_faces.insert(merged_faces.end(), nodes[second].faces.begin(),
                                    nodes[second].faces.end());
                VolumeFit fit = fitBestVolume(mesh, merged_faces, options);
                const double added = fit.volume - nodes[first].fit.volume -
                                     nodes[second].fit.volume;
                const double source_area = nodes[first].source_area + nodes[second].source_area;
                const double area_excess_ratio = std::max(
                    fit.surface_area - source_area, 0.0) /
                    std::max(model_surface_area, 1.0e-30);
                const double volume_excess_ratio = std::max(added, 0.0) / model_volume;
                const double separation_ratio = boundsSeparationRatio(
                    nodes[first].lower, nodes[first].upper,
                    nodes[second].lower, nodes[second].upper);
                const double object_boundary_ratio =
                    (nodes[first].recognized_protrusion || nodes[second].recognized_protrusion)
                        ? options.protrusion_max_area_excess_ratio : 0.0;
                const double score = std::max({area_excess_ratio, volume_excess_ratio,
                                               separation_ratio, object_boundary_ratio});
                candidates.push_back({score, fit.volume, first, second, 0, 0, fit});
            }
        }
        if (candidates.empty()) break;
        std::sort(candidates.begin(), candidates.end(), [](const Candidate& first,
                                                            const Candidate& second)
        {
            if (first.added_ratio != second.added_ratio)
                return first.added_ratio < second.added_ratio;
            if (first.volume != second.volume) return first.volume < second.volume;
            if (first.first != second.first) return first.first < second.first;
            return first.second < second.second;
        });

        std::vector<bool> matched(nodes.size(), false);
        std::size_t merged_count = 0;
        for (const Candidate& candidate : candidates)
        {
            if (candidate.first >= nodes.size() || candidate.second >= nodes.size() ||
                matched[candidate.first] || matched[candidate.second] ||
                !nodes[candidate.first].active || !nodes[candidate.second].active) continue;
            Node& first = nodes[candidate.first];
            Node& second = nodes[candidate.second];
            Node merged;
            merged.faces = first.faces;
            merged.faces.insert(merged.faces.end(), second.faces.begin(), second.faces.end());
            merged.fit = candidate.fit;
            merged.lower = first.lower.cwiseMin(second.lower);
            merged.upper = first.upper.cwiseMax(second.upper);
            merged.source_area = first.source_area + second.source_area;
            merged.activation_threshold = std::max(
                {candidate.added_ratio, first.activation_threshold,
                 second.activation_threshold});
            merged.neighbors = first.neighbors;
            merged.neighbors.insert(second.neighbors.begin(), second.neighbors.end());
            merged.neighbors.erase(candidate.first);
            merged.neighbors.erase(candidate.second);
            first.active = false;
            second.active = false;
            const std::uint32_t merged_id = static_cast<std::uint32_t>(nodes.size());
            first.parent = merged_id;
            second.parent = merged_id;
            matched[candidate.first] = true;
            matched[candidate.second] = true;
            nodes.push_back(std::move(merged));
            for (const auto neighbor : nodes.back().neighbors)
            {
                nodes[neighbor].neighbors.erase(candidate.first);
                nodes[neighbor].neighbors.erase(candidate.second);
                nodes[neighbor].neighbors.insert(merged_id);
            }
            ++merged_count;
        }
        if (merged_count == 0) break;
        active_count -= merged_count;
    }
    const double hierarchy_seconds = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - sweep_started).count();

    std::filesystem::create_directories(output_root);
    const auto common_source = output_root / "source.obj";
    writeSourceObj(common_source, mesh);
    std::vector<PrimitiveMeshAnalysisStats> results;
    results.reserve(101);
    for (int strength_index = 0; strength_index <= 100; ++strength_index)
    {
        const auto emit_started = std::chrono::steady_clock::now();
        PrimitiveMeshAnalysisOptions strength_options = options;
        strength_options.analysis_strength = static_cast<double>(strength_index) / 100.0;
        const double threshold = analysisStrengthToAllowedExcessRatio(
            strength_options.analysis_strength);
        PrimitiveMeshAnalysisStats stats = base_stats;
        std::vector<OutputPrimitive> output;
        for (std::uint32_t node_id = 0; node_id < nodes.size(); ++node_id)
        {
            const Node& node = nodes[node_id];
            if (node.activation_threshold > threshold) continue;
            if (node.parent != std::numeric_limits<std::uint32_t>::max() &&
                nodes[node.parent].activation_threshold <= threshold) continue;
            if (node.recognized_closed_box)
            {
                appendBoxRectangles(output, node.fit.box, node.faces,
                                    node.covered_box_face_axis,
                                    node.covered_box_face_sign);
                ++stats.recognized_closed_box_shells;
                continue;
            }
            if (node.recognized_protrusion)
            {
                appendBoxRectangles(output, node.fit.box, node.faces,
                                    node.covered_box_face_axis,
                                    node.covered_box_face_sign);
                ++stats.recognized_protrusion_box_shells;
                continue;
            }
            auto classified = classifyFinalRegion(
                mesh, node.faces, strength_options, threshold, model_surface_area,
                stats.filled_planar_holes, stats.filled_boundary_voids,
                stats.filled_boundary_void_area);
            output.insert(output.end(), std::make_move_iterator(classified.begin()),
                          std::make_move_iterator(classified.end()));
        }
        output = mergeLocalCoplanarPrimitives(
            mesh, std::move(output),
            std::max(diagonal * 1.0e-9, 1.0e-10),
            strength_options.maximum_open_error_distance >= 0.0
                ? strength_options.maximum_open_error_distance
                : diagonal * 0.08,
            std::max(diagonal / 192.0, 1.0e-30),
            stats.merged_local_planar_primitives);
        CoplanarCanonicalizationStats coplanar_stats;
        output = canonicalizeCoplanarPrimitiveUnion(
            mesh, std::move(output), std::max(diagonal * 1.0e-9, 1.0e-10),
            coplanar_stats);
        stats.canonicalized_coplanar_groups = coplanar_stats.groups;
        stats.removed_coplanar_redundant_primitives = coplanar_stats.removed_primitives;
        stats.removed_coplanar_overlap_area = coplanar_stats.removed_overlap_area;
        promoteToSemanticPrimitives(output);
        output = coalesceShallowParallelShells(
            std::move(output), (lower + upper) * 0.5, diagonal, model_volume,
            strength_options, std::max(diagonal * 1.0e-9, 1.0e-10), false);
        // Keep the same conservative export contract as the production path:
        // exact coplanar union is allowed, bounding-rectangle silhouette
        // replacement and projected-only parallel deletion are not.
        CoplanarCanonicalizationStats export_coplanar_stats;
        output = canonicalizeCoplanarPrimitiveUnion(
            mesh, std::move(output), std::max(diagonal * 1.0e-9, 1.0e-10),
            export_coplanar_stats);
        stats.canonicalized_coplanar_groups += export_coplanar_stats.groups;
        stats.removed_coplanar_redundant_primitives +=
            export_coplanar_stats.removed_primitives;
        stats.removed_coplanar_overlap_area +=
            export_coplanar_stats.removed_overlap_area;
        output = coalesceShallowParallelShells(
            std::move(output), (lower + upper) * 0.5, diagonal, model_volume,
            strength_options, std::max(diagonal * 1.0e-9, 1.0e-10), true);
        CoplanarCanonicalizationStats terrace_coplanar_stats;
        output = canonicalizeCoplanarPrimitiveUnion(
            mesh, std::move(output), std::max(diagonal * 1.0e-9, 1.0e-10),
            terrace_coplanar_stats);
        stats.canonicalized_coplanar_groups += terrace_coplanar_stats.groups;
        stats.removed_coplanar_redundant_primitives +=
            terrace_coplanar_stats.removed_primitives;
        stats.removed_coplanar_overlap_area +=
            terrace_coplanar_stats.removed_overlap_area;
        simplifyConservativePolygonDetails(
            output, model_surface_area,
            strength_options.shallow_shell_silhouette_added_area_ratio,
            std::max(diagonal * 1.0e-9, 1.0e-10));
        synchronizeOpposingShellSilhouettes(
            output, model_surface_area,
            strength_options.shallow_shell_silhouette_added_area_ratio,
            std::max(diagonal * 1.0e-9, 1.0e-10));
        output = rebuildCertifiedExtrudedShells(
            mesh, std::move(output), std::max(diagonal * 1.0e-9, 1.0e-10));
        CoplanarCanonicalizationStats extrusion_coplanar_stats;
        output = canonicalizeCoplanarPrimitiveUnion(
            mesh, std::move(output), std::max(diagonal * 1.0e-9, 1.0e-10),
            extrusion_coplanar_stats);
        stats.canonicalized_coplanar_groups += extrusion_coplanar_stats.groups;
        stats.removed_coplanar_redundant_primitives +=
            extrusion_coplanar_stats.removed_primitives;
        stats.removed_coplanar_overlap_area +=
            extrusion_coplanar_stats.removed_overlap_area;
        stats.primitive_count = output.size();
        for (const auto& item : output)
        {
            switch (item.primitive.kind)
            {
            case Kind::Polygon: ++stats.polygon_count; break;
            case Kind::Frustum: ++stats.frustum_count; break;
            case Kind::Disk: ++stats.disk_count; break;
            case Kind::Annulus: ++stats.annulus_count; break;
            case Kind::CylindricalBand: ++stats.cylindrical_band_count; break;
            case Kind::ConicalBand: ++stats.conical_band_count; break;
            case Kind::Rectangle:
            case Kind::Triangle: throw std::logic_error("internal planar primitive was not promoted");
            }
        }

        std::ostringstream directory_name;
        directory_name << 's' << std::setw(3) << std::setfill('0') << strength_index;
        const auto directory = output_root / directory_name.str();
        std::filesystem::create_directories(directory);
        std::error_code link_error;
        const auto source_link = directory / "source.obj";
        if (!std::filesystem::exists(source_link))
        {
            std::filesystem::create_hard_link(common_source, source_link, link_error);
            if (link_error)
                std::filesystem::copy_file(
                    common_source, source_link,
                    std::filesystem::copy_options::overwrite_existing);
        }
        writeRegionsObj(directory / "regions.obj", mesh, output,
                        structural_cleanup.excluded_faces);
        writeSemanticPrimitiveObj(directory / "primitives.obj", output);
        stats.proxy_triangles = writeTriangulatedObj(directory / "proxy.obj", output);
        stats.analysis_seconds = hierarchy_seconds + std::chrono::duration<double>(
            std::chrono::steady_clock::now() - emit_started).count();
        writeMetadata(directory, input_obj, output, stats, strength_options);
        results.push_back(stats);
    }

    std::ofstream manifest(output_root / "viewer_manifest.json");
    const std::string model_id = input_obj.stem().string();
    manifest << "{\"algorithm\":\"CppHierarchicalPrimitiveMeshAnalysis\","
             << "\"complete\":true,\"model_count\":1,\"models\":[{\"id\":"
             << model_id << ",\"metadata\":\"s040/model.json\"}],"
             << "\"strength_variants\":[";
    for (int index = 0; index <= 100; ++index)
    {
        if (index) manifest << ',';
        manifest << "{\"model_id\":" << model_id << ",\"strength\":"
                 << std::setprecision(17) << static_cast<double>(index) / 100.0
                 << ",\"metadata\":\"s" << std::setw(3) << std::setfill('0') << index
                 << "/model.json\"}" << std::setfill(' ');
    }
    manifest << "]}\n";
    return results;
}

} // namespace pqss_proxy_mesh
