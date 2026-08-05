#include "pqss_proxy_mesh/primitive_bvh_pool.hpp"

#include <algorithm>
#include <array>
#include <chrono>
#include <charconv>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <limits>
#include <numbers>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <unordered_map>
#include <utility>
#include <vector>

#include "pqss/dist.hpp"
#include "pqss/model.hpp"
#include "pqss/utils/mat_vec.hpp"

namespace pqss_proxy_mesh
{
namespace
{

using Clock = std::chrono::steady_clock;
constexpr std::size_t kInvalidIndex = std::numeric_limits<std::size_t>::max();
constexpr int kSahBuckets = 16;
constexpr std::uint32_t kPrimitiveBvhCacheVersion = 2;
constexpr std::array<char, 8> kPrimitiveBvhCacheMagic = {'P', 'Q', 'P', 'B', 'V', 'H', '0', '1'};

struct Bounds
{
    pqss::Vec3 minimum = {
        std::numeric_limits<pqss::Real>::max(),
        std::numeric_limits<pqss::Real>::max(),
        std::numeric_limits<pqss::Real>::max(),
    };
    pqss::Vec3 maximum = {
        std::numeric_limits<pqss::Real>::lowest(),
        std::numeric_limits<pqss::Real>::lowest(),
        std::numeric_limits<pqss::Real>::lowest(),
    };

    void expand(const pqss::Vec3& point)
    {
        for (int axis = 0; axis < 3; ++axis)
        {
            minimum[axis] = std::min(minimum[axis], point[axis]);
            maximum[axis] = std::max(maximum[axis], point[axis]);
        }
    }

    void expand(const Bounds& other)
    {
        expand(other.minimum);
        expand(other.maximum);
    }

    [[nodiscard]] pqss::Real surfaceArea() const
    {
        const pqss::Vec3 extent = {
            maximum[0] - minimum[0],
            maximum[1] - minimum[1],
            maximum[2] - minimum[2],
        };
        if (extent[0] < pqss::k_zero || extent[1] < pqss::k_zero || extent[2] < pqss::k_zero)
            return pqss::k_zero;
        return static_cast<pqss::Real>(2) *
               (extent[0] * extent[1] + extent[1] * extent[2] + extent[2] * extent[0]);
    }
};

struct Shape
{
    pqss::Mat3 R{};
    pqss::Vec3 T{};
    std::array<pqss::Real, 2> lengths{};
    pqss::Real radius = pqss::k_zero;
    pqss::Real size = pqss::k_zero;
    AnalyticPrimitiveType type = AnalyticPrimitiveType::Rss;
};

struct Primitive
{
    Shape shape;
    std::vector<std::size_t> triangle_ids;
    Bounds bounds;
    pqss::Vec3 centroid{};
};

struct Node
{
    Shape shape;
    std::size_t left = kInvalidIndex;
    std::size_t right = kInvalidIndex;
    std::size_t primitive = kInvalidIndex;

    [[nodiscard]] bool leaf() const
    {
        return primitive != kInvalidIndex;
    }
};

struct ModelTree
{
    std::vector<pqss::build::Tri> source_triangles;
    std::vector<Primitive> primitives;
    std::vector<std::size_t> primitive_order;
    std::vector<Node> nodes;
    PrimitiveBvhModelStats stats;
};

std::vector<std::string_view> Split(const std::string& line, const char delimiter)
{
    std::vector<std::string_view> result;
    std::size_t begin = 0;
    while (begin <= line.size())
    {
        const std::size_t end = line.find(delimiter, begin);
        result.emplace_back(line.data() + begin,
                            (end == std::string::npos ? line.size() : end) - begin);
        if (end == std::string::npos) break;
        begin = end + 1;
    }
    return result;
}

std::size_t ParseIndex(const std::string_view value, const std::string_view field)
{
    std::size_t result = 0;
    const auto parsed = std::from_chars(value.data(), value.data() + value.size(), result);
    if (parsed.ec != std::errc{} || parsed.ptr != value.data() + value.size())
        throw std::runtime_error("invalid " + std::string(field) + ": " + std::string(value));
    return result;
}

pqss::Real ParseReal(const std::string_view value, const std::string_view field)
{
    std::string text(value);
    std::size_t parsed = 0;
    const double result = std::stod(text, &parsed);
    if (parsed != text.size() || !std::isfinite(result))
        throw std::runtime_error("invalid " + std::string(field) + ": " + text);
    return static_cast<pqss::Real>(result);
}

AnalyticPrimitiveType ParseType(const std::string_view value)
{
    if (value == "sphere") return AnalyticPrimitiveType::Sphere;
    if (value == "capsule") return AnalyticPrimitiveType::Capsule;
    if (value == "rss") return AnalyticPrimitiveType::Rss;
    throw std::runtime_error("unsupported analytic primitive type: " + std::string(value));
}

pqss::Real ShapeSize(const Shape& shape)
{
    return shape.lengths[0] * shape.lengths[1] +
           std::numbers::pi_v<pqss::Real> * shape.radius *
               (shape.lengths[0] + shape.lengths[1] + static_cast<pqss::Real>(2) * shape.radius);
}

pqss::Mat3 ComposeRotation(const pqss::Mat3& parent, const pqss::Mat3& relative)
{
    return pqss::MxM(parent, relative);
}

pqss::Vec3 ComposeTranslation(const pqss::Mat3& parent_R,
                              const pqss::Vec3& parent_T,
                              const pqss::Vec3& relative_T)
{
    return pqss::MxVpV(parent_R, relative_T, parent_T);
}

pqss::Vec3 RectangleCorner(const Shape& shape, const int x, const int y)
{
    pqss::Vec3 result = shape.T;
    for (int row = 0; row < 3; ++row)
    {
        result[row] += shape.R[row][0] * shape.lengths[0] * static_cast<pqss::Real>(x);
        result[row] += shape.R[row][1] * shape.lengths[1] * static_cast<pqss::Real>(y);
    }
    return result;
}

Bounds ShapeBounds(const Shape& shape)
{
    Bounds bounds;
    for (int x = 0; x <= 1; ++x)
    {
        for (int y = 0; y <= 1; ++y)
        {
            pqss::Vec3 minimum = RectangleCorner(shape, x, y);
            pqss::Vec3 maximum = minimum;
            for (int axis = 0; axis < 3; ++axis)
            {
                minimum[axis] -= shape.radius;
                maximum[axis] += shape.radius;
            }
            bounds.expand(minimum);
            bounds.expand(maximum);
        }
    }
    return bounds;
}

pqss::Real PointRectangleDistanceSq(const pqss::Vec3& point, const Shape& shape)
{
    const pqss::Vec3 local = pqss::MTxVmV(shape.R, point, shape.T);
    const pqss::Real x = std::clamp(local[0], pqss::k_zero, shape.lengths[0]);
    const pqss::Real y = std::clamp(local[1], pqss::k_zero, shape.lengths[1]);
    const pqss::Vec3 delta = {local[0] - x, local[1] - y, local[2]};
    return pqss::VdotV(delta, delta);
}

pqss::Mat3 InitialFrame(const std::vector<pqss::build::SubTri>& triangles)
{
    pqss::Vec3 first_moment{};
    pqss::Mat3 second_moment{};
    for (const auto& triangle : triangles)
    {
        for (int vertex = 0; vertex < 3; ++vertex)
        {
            const pqss::Vec3& point = triangle[static_cast<std::size_t>(vertex)];
            for (int row = 0; row < 3; ++row)
            {
                first_moment[row] += point[row];
                for (int column = 0; column < 3; ++column)
                    second_moment[row][column] += point[row] * point[column];
            }
        }
    }
    const pqss::Real count = static_cast<pqss::Real>(3 * triangles.size());
    pqss::Mat3 covariance{};
    for (int row = 0; row < 3; ++row)
        for (int column = 0; column < 3; ++column)
            covariance[row][column] =
                second_moment[row][column] - first_moment[row] * first_moment[column] / count;

    auto [eigenvectors, eigenvalues] = pqss::Meigen(covariance);
    int minimum = 0;
    int middle = 1;
    int maximum = 2;
    if (eigenvalues[maximum] < eigenvalues[middle]) std::swap(maximum, middle);
    if (eigenvalues[maximum] < eigenvalues[minimum]) std::swap(maximum, minimum);
    if (eigenvalues[middle] < eigenvalues[minimum]) std::swap(middle, minimum);
    pqss::Mat3 frame{};
    for (int row = 0; row < 3; ++row)
    {
        frame[row][0] = eigenvectors[row][maximum];
        frame[row][1] = eigenvectors[row][middle];
    }
    frame[0][2] = frame[1][0] * frame[2][1] - frame[2][0] * frame[1][1];
    frame[1][2] = frame[2][0] * frame[0][1] - frame[0][0] * frame[2][1];
    frame[2][2] = frame[0][0] * frame[1][1] - frame[1][0] * frame[0][1];
    return frame;
}

Shape FitRss(const std::vector<pqss::build::SubTri>& triangles,
             const pqss::BuildStrategy strategy)
{
    if (triangles.empty()) throw std::logic_error("cannot fit RSS to an empty triangle set");
    pqss::build::BV fitted;
    fitted.FitToTris(InitialFrame(triangles), triangles, 0, triangles.size(), strategy);
    Shape shape;
    shape.R = fitted.R();
    shape.T = fitted.Tr();
    shape.lengths = fitted.L();
    shape.radius = fitted.Radius();
    shape.size = fitted.Size();
    shape.type = AnalyticPrimitiveType::Rss;

    const pqss::Real scale = std::max({pqss::Real{1}, shape.lengths[0], shape.lengths[1], shape.radius,
                                       std::abs(shape.T[0]), std::abs(shape.T[1]), std::abs(shape.T[2])});
    const pqss::Real numerical_tolerance =
        scale * scale * std::numeric_limits<pqss::Real>::epsilon() * static_cast<pqss::Real>(256);
    const pqss::Real radius_sq = shape.radius * shape.radius;
    pqss::Real maximum_distance_sq = pqss::k_zero;
    pqss::Vec3 maximum_point{};
    for (const auto& triangle : triangles)
    {
        for (int vertex = 0; vertex < 3; ++vertex)
        {
            const pqss::Vec3& point = triangle[static_cast<std::size_t>(vertex)];
            const pqss::Real distance_sq = PointRectangleDistanceSq(point, shape);
            if (distance_sq > maximum_distance_sq)
            {
                maximum_distance_sq = distance_sq;
                maximum_point = point;
            }
        }
    }
    if (maximum_distance_sq > radius_sq + numerical_tolerance)
    {
        const pqss::Vec3 maximum_local = pqss::MTxVmV(shape.R, maximum_point, shape.T);
        const std::array<pqss::Real, 2> point_lengths{};
        const pqss::Real pqss_point_distance_sq = pqss::detail::RectDistSq(
            pqss::MTxM(shape.R, pqss::Meident()), maximum_local,
            shape.lengths, point_lengths);

        pqss::build::BV fast_fitted;
        fast_fitted.FitToTris(
            InitialFrame(triangles), triangles, 0, triangles.size(), pqss::BuildStrategy::Fast);
        Shape fast_shape;
        fast_shape.R = fast_fitted.R();
        fast_shape.T = fast_fitted.Tr();
        fast_shape.lengths = fast_fitted.L();
        fast_shape.radius = fast_fitted.Radius();
        pqss::Real fast_maximum_distance_sq = pqss::k_zero;
        for (const auto& triangle : triangles)
            for (int vertex = 0; vertex < 3; ++vertex)
                fast_maximum_distance_sq = std::max(
                    fast_maximum_distance_sq,
                    PointRectangleDistanceSq(triangle[static_cast<std::size_t>(vertex)], fast_shape));

        if (const char* diagnostic_path = std::getenv("PQSS_RSS_DIAGNOSTIC"))
        {
            std::ofstream diagnostic(diagnostic_path, std::ios::trunc);
            diagnostic.precision(17);
            diagnostic << "PQSS_RSS_CONTAINMENT_DIAGNOSTIC_V1\n";
            diagnostic << "triangles\t" << triangles.size() << '\n';
            diagnostic << "optimized_T\t" << shape.T[0] << '\t' << shape.T[1] << '\t' << shape.T[2] << '\n';
            for (int row = 0; row < 3; ++row)
                diagnostic << "optimized_R" << row << '\t' << shape.R[row][0] << '\t'
                           << shape.R[row][1] << '\t' << shape.R[row][2] << '\n';
            diagnostic << "optimized_LR\t" << shape.lengths[0] << '\t' << shape.lengths[1]
                       << '\t' << shape.radius << '\n';
            diagnostic << "maximum_point\t" << maximum_point[0] << '\t' << maximum_point[1]
                       << '\t' << maximum_point[2] << '\n';
            diagnostic << "maximum_local\t" << maximum_local[0] << '\t' << maximum_local[1]
                       << '\t' << maximum_local[2] << '\n';
            diagnostic << "direct_distance_sq\t" << maximum_distance_sq << '\n';
            diagnostic << "pqss_rect_dist_sq\t" << pqss_point_distance_sq << '\n';
            diagnostic << "radius_sq\t" << radius_sq << '\n';
            diagnostic << "fast_radius_sq\t" << fast_shape.radius * fast_shape.radius << '\n';
            diagnostic << "fast_maximum_distance_sq\t" << fast_maximum_distance_sq << '\n';
            for (const auto& triangle : triangles)
                diagnostic << "triangle\t" << triangle.p1[0] << '\t' << triangle.p1[1] << '\t'
                           << triangle.p1[2] << '\t' << triangle.p2[0] << '\t' << triangle.p2[1]
                           << '\t' << triangle.p2[2] << '\t' << triangle.p3[0] << '\t'
                           << triangle.p3[1] << '\t' << triangle.p3[2] << '\n';
        }
        std::ostringstream message;
        message << "PQSS Optimized RSS failed internal-node containment certification: triangles="
                << triangles.size() << " radius_sq=" << radius_sq
                << " maximum_distance_sq=" << maximum_distance_sq
                << " pqss_point_distance_sq=" << pqss_point_distance_sq
                << " comparison_tolerance=" << numerical_tolerance
                << " excess=" << maximum_distance_sq - radius_sq
                << " fast_radius_sq=" << fast_shape.radius * fast_shape.radius
                << " fast_maximum_distance_sq=" << fast_maximum_distance_sq;
        throw std::runtime_error(message.str());
    }
    return shape;
}

std::vector<pqss::build::Tri> ReadObj(const std::filesystem::path& path)
{
    std::ifstream stream(path);
    if (!stream) throw std::runtime_error("cannot open source OBJ: " + path.string());
    std::vector<pqss::Vec3> vertices;
    std::vector<pqss::build::Tri> triangles;
    std::string line;
    std::size_t line_number = 0;
    while (std::getline(stream, line))
    {
        ++line_number;
        std::istringstream fields(line);
        std::string kind;
        fields >> kind;
        if (kind.empty() || kind[0] == '#') continue;
        if (kind == "v")
        {
            pqss::Vec3 vertex{};
            if (!(fields >> vertex[0] >> vertex[1] >> vertex[2]))
                throw std::runtime_error(path.string() + ": invalid vertex at line " + std::to_string(line_number));
            vertices.push_back(vertex);
        }
        else if (kind == "f")
        {
            std::vector<std::size_t> polygon;
            std::string token;
            while (fields >> token)
            {
                const std::size_t slash = token.find('/');
                const std::string index_text = token.substr(0, slash);
                const long long raw = std::stoll(index_text);
                const long long index = raw > 0 ? raw - 1 : static_cast<long long>(vertices.size()) + raw;
                if (index < 0 || static_cast<std::size_t>(index) >= vertices.size())
                    throw std::runtime_error(path.string() + ": face index out of range at line " +
                                             std::to_string(line_number));
                polygon.push_back(static_cast<std::size_t>(index));
            }
            if (polygon.size() < 3)
                throw std::runtime_error(path.string() + ": incomplete face at line " +
                                         std::to_string(line_number));
            for (std::size_t index = 1; index + 1 < polygon.size(); ++index)
                triangles.emplace_back(vertices[polygon[0]], vertices[polygon[index]], vertices[polygon[index + 1]]);
        }
    }
    if (triangles.empty()) throw std::runtime_error("source OBJ has no triangles: " + path.string());
    return triangles;
}

pqss::Real RectangleDistanceSq(const Shape& first,
                               const Shape& second,
                               const pqss::Mat3& object_R,
                               const pqss::Vec3& object_T)
{
    const pqss::Mat3 relative_R = pqss::MTxM(first.R, pqss::MxM(object_R, second.R));
    const pqss::Vec3 relative_T =
        pqss::MTxVmV(first.R, pqss::MxVpV(object_R, second.T, object_T), first.T);
    return pqss::detail::RectDistSq(relative_R, relative_T, first.lengths, second.lengths);
}

bool Separated(const Shape& first,
               const Shape& second,
               const pqss::Mat3& object_R,
               const pqss::Vec3& object_T,
               const pqss::Real tolerance,
               pqss::Real& rectangle_distance_sq)
{
    rectangle_distance_sq = RectangleDistanceSq(first, second, object_R, object_T);
    const pqss::Real threshold = tolerance + first.radius + second.radius;
    return rectangle_distance_sq > threshold * threshold;
}

struct PrimitiveBvhCacheHeader
{
    char magic[8]{};
    std::uint32_t version = 0;
    std::uint32_t real_size = 0;
    std::uint32_t build_strategy = 0;
    std::uint32_t reserved = 0;
    std::uint64_t fingerprint = 0;
    std::uint64_t model_id = 0;
    std::uint64_t source_triangles = 0;
    std::uint64_t primitives = 0;
    std::uint64_t nodes = 0;
    std::uint64_t max_depth = 0;
    std::uint64_t primitive_types[3]{};
};

template <typename Value>
bool ReadExact(std::istream& stream, Value& value)
{
    return static_cast<bool>(stream.read(reinterpret_cast<char*>(&value),
                                         static_cast<std::streamsize>(sizeof(Value))));
}

template <typename Value>
bool WriteExact(std::ostream& stream, const Value& value)
{
    return static_cast<bool>(stream.write(reinterpret_cast<const char*>(&value),
                                          static_cast<std::streamsize>(sizeof(Value))));
}

bool ReadShape(std::istream& stream, Shape& shape)
{
    for (int row = 0; row < 3; ++row)
        for (int column = 0; column < 3; ++column)
            if (!ReadExact(stream, shape.R[row][column])) return false;
    for (int axis = 0; axis < 3; ++axis)
        if (!ReadExact(stream, shape.T[axis])) return false;
    std::uint32_t type = 0;
    if (!ReadExact(stream, shape.lengths[0]) || !ReadExact(stream, shape.lengths[1]) ||
        !ReadExact(stream, shape.radius) || !ReadExact(stream, shape.size) || !ReadExact(stream, type) ||
        type >= static_cast<std::uint32_t>(AnalyticPrimitiveType::Count))
        return false;
    shape.type = static_cast<AnalyticPrimitiveType>(type);
    return true;
}

bool WriteShape(std::ostream& stream, const Shape& shape)
{
    for (int row = 0; row < 3; ++row)
        for (int column = 0; column < 3; ++column)
            if (!WriteExact(stream, shape.R[row][column])) return false;
    for (int axis = 0; axis < 3; ++axis)
        if (!WriteExact(stream, shape.T[axis])) return false;
    const auto type = static_cast<std::uint32_t>(shape.type);
    return WriteExact(stream, shape.lengths[0]) && WriteExact(stream, shape.lengths[1]) &&
           WriteExact(stream, shape.radius) && WriteExact(stream, shape.size) && WriteExact(stream, type);
}

void HashBytes(std::uint64_t& hash, const char* data, const std::size_t size)
{
    constexpr std::uint64_t prime = 1099511628211ull;
    for (std::size_t index = 0; index < size; ++index)
    {
        hash ^= static_cast<std::uint64_t>(static_cast<unsigned char>(data[index]));
        hash *= prime;
    }
}

void HashFile(std::uint64_t& hash, const std::filesystem::path& path)
{
    std::ifstream stream(path, std::ios::binary);
    if (!stream) throw std::runtime_error("cannot fingerprint cache input: " + path.string());
    std::array<char, 64 * 1024> buffer{};
    while (stream)
    {
        stream.read(buffer.data(), static_cast<std::streamsize>(buffer.size()));
        HashBytes(hash, buffer.data(), static_cast<std::size_t>(stream.gcount()));
    }
}

std::filesystem::path SourcePathFromModelFile(const std::filesystem::path& model_path)
{
    std::ifstream stream(model_path);
    if (!stream) throw std::runtime_error("cannot open primitive model: " + model_path.string());
    std::string line;
    while (std::getline(stream, line))
    {
        const auto fields = Split(line, '\t');
        if (fields.size() == 2 && fields[0] == "source")
            return (model_path.parent_path() / std::filesystem::path(fields[1])).lexically_normal();
    }
    throw std::runtime_error("primitive model has no source row: " + model_path.string());
}

std::uint64_t CacheFingerprint(const std::filesystem::path& model_path,
                               const pqss::BuildStrategy strategy)
{
    std::uint64_t hash = 14695981039346656037ull;
    HashFile(hash, model_path);
    HashFile(hash, SourcePathFromModelFile(model_path));
    HashBytes(hash, reinterpret_cast<const char*>(&kPrimitiveBvhCacheVersion),
              sizeof(kPrimitiveBvhCacheVersion));
    const auto strategy_value = static_cast<std::uint32_t>(strategy);
    HashBytes(hash, reinterpret_cast<const char*>(&strategy_value), sizeof(strategy_value));
    const std::uint32_t real_size = sizeof(pqss::Real);
    HashBytes(hash, reinterpret_cast<const char*>(&real_size), sizeof(real_size));
    return hash;
}

bool ValidateCachedTree(const ModelTree& tree)
{
    if (tree.nodes.empty() || tree.stats.primitives == 0 ||
        tree.nodes.size() != 2 * tree.stats.primitives - 1)
        return false;
    for (const Node& node : tree.nodes)
    {
        if (!std::isfinite(node.shape.radius) || node.shape.radius < pqss::k_zero ||
            !std::isfinite(node.shape.lengths[0]) || node.shape.lengths[0] < pqss::k_zero ||
            !std::isfinite(node.shape.lengths[1]) || node.shape.lengths[1] < pqss::k_zero)
            return false;
        if (node.leaf())
        {
            if (node.left != kInvalidIndex || node.right != kInvalidIndex ||
                node.primitive >= tree.stats.primitives)
                return false;
        }
        else if (node.left >= tree.nodes.size() || node.right >= tree.nodes.size() || node.left == node.right)
            return false;
    }
    return true;
}

bool LoadTreeCache(const std::filesystem::path& path,
                   const std::size_t model_id,
                   const pqss::BuildStrategy strategy,
                   const std::uint64_t fingerprint,
                   ModelTree& tree)
{
    std::ifstream stream(path, std::ios::binary);
    if (!stream) return false;
    PrimitiveBvhCacheHeader header;
    if (!ReadExact(stream, header) ||
        !std::equal(kPrimitiveBvhCacheMagic.begin(), kPrimitiveBvhCacheMagic.end(), std::begin(header.magic)) ||
        header.version != kPrimitiveBvhCacheVersion || header.real_size != sizeof(pqss::Real) ||
        header.build_strategy != static_cast<std::uint32_t>(strategy) || header.fingerprint != fingerprint ||
        header.model_id != model_id || header.primitives == 0 || header.nodes != 2 * header.primitives - 1)
        return false;
    tree.stats.model_id = model_id;
    tree.stats.source_triangles = static_cast<std::size_t>(header.source_triangles);
    tree.stats.primitives = static_cast<std::size_t>(header.primitives);
    tree.stats.nodes = static_cast<std::size_t>(header.nodes);
    tree.stats.max_depth = static_cast<std::size_t>(header.max_depth);
    for (std::size_t index = 0; index < 3; ++index)
        tree.stats.primitive_types[index] = static_cast<std::size_t>(header.primitive_types[index]);
    tree.stats.analyzed = true;
    tree.stats.internal_containment_certified = true;
    try
    {
        tree.nodes.resize(tree.stats.nodes);
    }
    catch (const std::bad_alloc&)
    {
        return false;
    }
    for (Node& node : tree.nodes)
    {
        std::uint64_t left = 0;
        std::uint64_t right = 0;
        std::uint64_t primitive = 0;
        if (!ReadShape(stream, node.shape) || !ReadExact(stream, left) ||
            !ReadExact(stream, right) || !ReadExact(stream, primitive))
            return false;
        node.left = left == std::numeric_limits<std::uint64_t>::max()
            ? kInvalidIndex : static_cast<std::size_t>(left);
        node.right = right == std::numeric_limits<std::uint64_t>::max()
            ? kInvalidIndex : static_cast<std::size_t>(right);
        node.primitive = primitive == std::numeric_limits<std::uint64_t>::max()
            ? kInvalidIndex : static_cast<std::size_t>(primitive);
    }
    return ValidateCachedTree(tree);
}

bool SaveTreeCache(const std::filesystem::path& path,
                   const std::size_t model_id,
                   const pqss::BuildStrategy strategy,
                   const std::uint64_t fingerprint,
                   const ModelTree& tree)
{
    std::error_code error;
    std::filesystem::create_directories(path.parent_path(), error);
    if (error) return false;
    const std::filesystem::path temporary = path.string() + ".tmp";
    std::ofstream stream(temporary, std::ios::binary | std::ios::trunc);
    if (!stream) return false;
    PrimitiveBvhCacheHeader header;
    std::copy(kPrimitiveBvhCacheMagic.begin(), kPrimitiveBvhCacheMagic.end(), std::begin(header.magic));
    header.version = kPrimitiveBvhCacheVersion;
    header.real_size = sizeof(pqss::Real);
    header.build_strategy = static_cast<std::uint32_t>(strategy);
    header.fingerprint = fingerprint;
    header.model_id = model_id;
    header.source_triangles = tree.stats.source_triangles;
    header.primitives = tree.stats.primitives;
    header.nodes = tree.nodes.size();
    header.max_depth = tree.stats.max_depth;
    for (std::size_t index = 0; index < 3; ++index)
        header.primitive_types[index] = tree.stats.primitive_types[index];
    if (!WriteExact(stream, header)) return false;
    for (const Node& node : tree.nodes)
    {
        const std::uint64_t left = node.left;
        const std::uint64_t right = node.right;
        const std::uint64_t primitive = node.primitive;
        if (!WriteShape(stream, node.shape) || !WriteExact(stream, left) ||
            !WriteExact(stream, right) || !WriteExact(stream, primitive))
            return false;
    }
    stream.close();
    if (!stream) return false;
    std::filesystem::remove(path, error);
    error.clear();
    std::filesystem::rename(temporary, path, error);
    return !error;
}

} // namespace

struct PrimitiveBvhPool::Impl
{
    std::unordered_map<std::size_t, ModelTree> models;
    PrimitiveBvhCacheStats cache_stats;

    static void PreparePrimitive(Primitive& primitive)
    {
        primitive.shape.size = ShapeSize(primitive.shape);
        primitive.bounds = ShapeBounds(primitive.shape);
        primitive.centroid = primitive.shape.T;
        for (int row = 0; row < 3; ++row)
        {
            primitive.centroid[row] +=
                static_cast<pqss::Real>(0.5) * primitive.shape.R[row][0] * primitive.shape.lengths[0];
            primitive.centroid[row] +=
                static_cast<pqss::Real>(0.5) * primitive.shape.R[row][1] * primitive.shape.lengths[1];
        }
    }

    std::size_t ChooseSplit(ModelTree& tree, const std::size_t begin, const std::size_t end) const
    {
        struct Bucket
        {
            Bounds bounds;
            std::size_t count = 0;
        };
        pqss::Real best_cost = std::numeric_limits<pqss::Real>::max();
        int best_axis = -1;
        int best_bucket = -1;
        Bounds centroid_bounds;
        for (std::size_t offset = begin; offset < end; ++offset)
            centroid_bounds.expand(tree.primitives[tree.primitive_order[offset]].centroid);

        for (int axis = 0; axis < 3; ++axis)
        {
            const pqss::Real extent = centroid_bounds.maximum[axis] - centroid_bounds.minimum[axis];
            if (extent <= pqss::k_zero) continue;
            std::array<Bucket, kSahBuckets> buckets{};
            for (std::size_t offset = begin; offset < end; ++offset)
            {
                const Primitive& primitive = tree.primitives[tree.primitive_order[offset]];
                int bucket = static_cast<int>(kSahBuckets *
                    (primitive.centroid[axis] - centroid_bounds.minimum[axis]) / extent);
                bucket = std::clamp(bucket, 0, kSahBuckets - 1);
                ++buckets[static_cast<std::size_t>(bucket)].count;
                buckets[static_cast<std::size_t>(bucket)].bounds.expand(primitive.bounds);
            }
            std::array<Bounds, kSahBuckets> prefix_bounds{};
            std::array<Bounds, kSahBuckets> suffix_bounds{};
            std::array<std::size_t, kSahBuckets> prefix_count{};
            std::array<std::size_t, kSahBuckets> suffix_count{};
            for (int bucket = 0; bucket < kSahBuckets; ++bucket)
            {
                if (bucket > 0)
                {
                    prefix_bounds[static_cast<std::size_t>(bucket)] =
                        prefix_bounds[static_cast<std::size_t>(bucket - 1)];
                    prefix_count[static_cast<std::size_t>(bucket)] =
                        prefix_count[static_cast<std::size_t>(bucket - 1)];
                }
                if (buckets[static_cast<std::size_t>(bucket)].count > 0)
                    prefix_bounds[static_cast<std::size_t>(bucket)].expand(
                        buckets[static_cast<std::size_t>(bucket)].bounds);
                prefix_count[static_cast<std::size_t>(bucket)] += buckets[static_cast<std::size_t>(bucket)].count;
            }
            for (int bucket = kSahBuckets - 1; bucket >= 0; --bucket)
            {
                if (bucket + 1 < kSahBuckets)
                {
                    suffix_bounds[static_cast<std::size_t>(bucket)] =
                        suffix_bounds[static_cast<std::size_t>(bucket + 1)];
                    suffix_count[static_cast<std::size_t>(bucket)] =
                        suffix_count[static_cast<std::size_t>(bucket + 1)];
                }
                if (buckets[static_cast<std::size_t>(bucket)].count > 0)
                    suffix_bounds[static_cast<std::size_t>(bucket)].expand(
                        buckets[static_cast<std::size_t>(bucket)].bounds);
                suffix_count[static_cast<std::size_t>(bucket)] += buckets[static_cast<std::size_t>(bucket)].count;
            }
            for (int bucket = 0; bucket + 1 < kSahBuckets; ++bucket)
            {
                const std::size_t left_count = prefix_count[static_cast<std::size_t>(bucket)];
                const std::size_t right_count = suffix_count[static_cast<std::size_t>(bucket + 1)];
                if (left_count == 0 || right_count == 0) continue;
                const pqss::Real cost =
                    prefix_bounds[static_cast<std::size_t>(bucket)].surfaceArea() *
                        static_cast<pqss::Real>(left_count) +
                    suffix_bounds[static_cast<std::size_t>(bucket + 1)].surfaceArea() *
                        static_cast<pqss::Real>(right_count);
                if (cost < best_cost)
                {
                    best_cost = cost;
                    best_axis = axis;
                    best_bucket = bucket;
                }
            }
        }

        const std::size_t middle = begin + (end - begin) / 2;
        if (best_axis < 0)
            return middle;
        const pqss::Real extent = centroid_bounds.maximum[best_axis] - centroid_bounds.minimum[best_axis];
        auto split = std::stable_partition(
            tree.primitive_order.begin() + static_cast<std::ptrdiff_t>(begin),
            tree.primitive_order.begin() + static_cast<std::ptrdiff_t>(end),
            [&](const std::size_t primitive_index)
            {
                int bucket = static_cast<int>(kSahBuckets *
                    (tree.primitives[primitive_index].centroid[best_axis] - centroid_bounds.minimum[best_axis]) /
                    extent);
                bucket = std::clamp(bucket, 0, kSahBuckets - 1);
                return bucket <= best_bucket;
            });
        const std::size_t split_index =
            static_cast<std::size_t>(std::distance(tree.primitive_order.begin(), split));
        if (split_index != begin && split_index != end)
            return split_index;

        std::nth_element(
            tree.primitive_order.begin() + static_cast<std::ptrdiff_t>(begin),
            tree.primitive_order.begin() + static_cast<std::ptrdiff_t>(middle),
            tree.primitive_order.begin() + static_cast<std::ptrdiff_t>(end),
            [&](const std::size_t first, const std::size_t second)
            {
                return tree.primitives[first].centroid[best_axis] < tree.primitives[second].centroid[best_axis];
            });
        return middle;
    }

    Shape FitRange(const ModelTree& tree,
                   const std::size_t begin,
                   const std::size_t end,
                   const pqss::BuildStrategy strategy) const
    {
        std::size_t count = 0;
        for (std::size_t offset = begin; offset < end; ++offset)
            count += tree.primitives[tree.primitive_order[offset]].triangle_ids.size();
        std::vector<pqss::build::SubTri> triangles;
        triangles.reserve(count);
        for (std::size_t offset = begin; offset < end; ++offset)
        {
            const Primitive& primitive = tree.primitives[tree.primitive_order[offset]];
            for (const std::size_t triangle_id : primitive.triangle_ids)
            {
                const auto& triangle = tree.source_triangles.at(triangle_id);
                triangles.emplace_back(triangle.p1, triangle.p2, triangle.p3, triangle_id);
            }
        }
        return FitRss(triangles, strategy);
    }

    std::size_t BuildRange(ModelTree& tree,
                           const std::size_t begin,
                           const std::size_t end,
                           const std::size_t depth,
                           const pqss::BuildStrategy strategy) const
    {
        const std::size_t node_index = tree.nodes.size();
        tree.nodes.emplace_back();
        tree.stats.max_depth = std::max(tree.stats.max_depth, depth);
        if (end - begin == 1)
        {
            const std::size_t primitive_index = tree.primitive_order[begin];
            tree.nodes[node_index].shape = tree.primitives[primitive_index].shape;
            tree.nodes[node_index].primitive = primitive_index;
            return node_index;
        }
        const std::size_t split = ChooseSplit(tree, begin, end);
        tree.nodes[node_index].shape = FitRange(tree, begin, end, strategy);
        tree.nodes[node_index].left = BuildRange(tree, begin, split, depth + 1, strategy);
        tree.nodes[node_index].right = BuildRange(tree, split, end, depth + 1, strategy);
        return node_index;
    }

    ModelTree LoadAnalyzedModel(const std::size_t model_id,
                                const std::filesystem::path& path,
                                const pqss::BuildStrategy strategy) const
    {
        std::ifstream stream(path);
        if (!stream) throw std::runtime_error("cannot open primitive model: " + path.string());
        std::string line;
        if (!std::getline(stream, line) || line != "PQSS_PRIMITIVE_BVH_MODEL_V1")
            throw std::runtime_error("invalid primitive model header: " + path.string());
        ModelTree tree;
        std::size_t declared_triangles = 0;
        std::size_t declared_primitives = 0;
        std::filesystem::path source_path;
        while (std::getline(stream, line))
        {
            if (line.empty()) continue;
            const auto fields = Split(line, '\t');
            if (fields[0] == "source")
            {
                if (fields.size() != 2) throw std::runtime_error("invalid source row in " + path.string());
                source_path = path.parent_path() / std::filesystem::path(fields[1]);
            }
            else if (fields[0] == "source_triangles")
            {
                declared_triangles = ParseIndex(fields.at(1), "source triangle count");
            }
            else if (fields[0] == "primitive_count")
            {
                declared_primitives = ParseIndex(fields.at(1), "primitive count");
                tree.primitives.reserve(declared_primitives);
            }
            else if (fields[0] == "primitive")
            {
                if (fields.size() != 19) throw std::runtime_error("invalid primitive row in " + path.string());
                const std::size_t primitive_id = ParseIndex(fields[1], "primitive id");
                if (primitive_id != tree.primitives.size())
                    throw std::runtime_error("primitive IDs must be contiguous in " + path.string());
                Primitive primitive;
                primitive.shape.type = ParseType(fields[2]);
                for (int axis = 0; axis < 3; ++axis)
                    primitive.shape.T[axis] = ParseReal(fields[static_cast<std::size_t>(3 + axis)], "origin");
                std::size_t field = 6;
                for (int row = 0; row < 3; ++row)
                    for (int column = 0; column < 3; ++column)
                        primitive.shape.R[row][column] = ParseReal(fields[field++], "axis");
                primitive.shape.lengths[0] = ParseReal(fields[15], "length 0");
                primitive.shape.lengths[1] = ParseReal(fields[16], "length 1");
                primitive.shape.radius = ParseReal(fields[17], "radius");
                const std::string triangle_list(fields[18]);
                for (const std::string_view triangle : Split(triangle_list, ','))
                    primitive.triangle_ids.push_back(ParseIndex(triangle, "triangle id"));
                PreparePrimitive(primitive);
                tree.primitives.push_back(std::move(primitive));
            }
        }
        tree.source_triangles = ReadObj(source_path.lexically_normal());
        if (tree.source_triangles.size() != declared_triangles)
            throw std::runtime_error("source triangle count mismatch for model " + std::to_string(model_id));
        if (tree.primitives.size() != declared_primitives || tree.primitives.empty())
            throw std::runtime_error("primitive count mismatch for model " + std::to_string(model_id));
        std::vector<std::uint8_t> assigned(tree.source_triangles.size(), 0);
        for (const Primitive& primitive : tree.primitives)
        {
            for (const std::size_t triangle_id : primitive.triangle_ids)
            {
                if (triangle_id >= assigned.size() || assigned[triangle_id] != 0)
                    throw std::runtime_error("triangle responsibility is not a strict partition for model " +
                                             std::to_string(model_id));
                assigned[triangle_id] = 1;
            }
        }
        if (std::find(assigned.begin(), assigned.end(), 0) != assigned.end())
            throw std::runtime_error("triangle responsibility is incomplete for model " + std::to_string(model_id));

        tree.primitive_order.resize(tree.primitives.size());
        std::iota(tree.primitive_order.begin(), tree.primitive_order.end(), 0);
        tree.nodes.reserve(2 * tree.primitives.size() - 1);
        tree.stats.model_id = model_id;
        tree.stats.source_triangles = tree.source_triangles.size();
        tree.stats.primitives = tree.primitives.size();
        tree.stats.analyzed = true;
        tree.stats.internal_containment_certified = true;
        for (const Primitive& primitive : tree.primitives)
            ++tree.stats.primitive_types[static_cast<std::size_t>(primitive.shape.type)];
        BuildRange(tree, 0, tree.primitive_order.size(), 0, strategy);
        tree.stats.nodes = tree.nodes.size();
        tree.source_triangles.clear();
        tree.source_triangles.shrink_to_fit();
        tree.primitives.clear();
        tree.primitives.shrink_to_fit();
        tree.primitive_order.clear();
        tree.primitive_order.shrink_to_fit();
        return tree;
    }

    void ImportNode(ModelTree& tree,
                    const pqss::build::Model& source_model,
                    const std::size_t index,
                    const pqss::Mat3& parent_R,
                    const pqss::Vec3& parent_T,
                    const std::size_t depth) const
    {
        const pqss::build::BV& source = source_model.Bvs().at(index);
        Node& node = tree.nodes.at(index);
        node.shape.R = ComposeRotation(parent_R, source.R());
        node.shape.T = ComposeTranslation(parent_R, parent_T, source.Tr());
        node.shape.lengths = source.L();
        node.shape.radius = source.Radius();
        node.shape.size = source.Size();
        node.shape.type = AnalyticPrimitiveType::Rss;
        tree.stats.max_depth = std::max(tree.stats.max_depth, depth);
        if (source.IsLeaf())
        {
            node.primitive = tree.stats.primitives++;
            ++tree.stats.primitive_types[static_cast<std::size_t>(AnalyticPrimitiveType::Rss)];
            return;
        }
        node.left = source.FirstChild();
        node.right = source.FirstChild() + 1;
        ImportNode(tree, source_model, node.left, node.shape.R, node.shape.T, depth + 1);
        ImportNode(tree, source_model, node.right, node.shape.R, node.shape.T, depth + 1);
    }

    struct QueryState
    {
        PrimitiveBvhQueryResult result;
        pqss::Mat3 object_R{};
        pqss::Vec3 object_T{};
        pqss::Real tolerance = pqss::k_zero;
    };

    bool TestPair(const Node& first,
                  const Node& second,
                  QueryState& state,
                  const bool root,
                  pqss::Real& rectangle_distance_sq) const
    {
        if (root) ++state.result.root_bv_tests;
        else ++state.result.internal_bv_tests;
        return !Separated(first.shape, second.shape, state.object_R, state.object_T,
                          state.tolerance, rectangle_distance_sq);
    }

    void Traverse(const ModelTree& first_tree,
                  const std::size_t first_index,
                  const ModelTree& second_tree,
                  const std::size_t second_index,
                  QueryState& state) const
    {
        if (state.result.closer_than_tolerance) return;
        const Node& first = first_tree.nodes[first_index];
        const Node& second = second_tree.nodes[second_index];
        if (first.leaf() && second.leaf())
        {
            ++state.result.leaf_pair_tests;
            ++state.result.type_pair_tests[static_cast<std::size_t>(first.shape.type)]
                                                [static_cast<std::size_t>(second.shape.type)];
            state.result.closer_than_tolerance = true;
            return;
        }
        struct ChildPair
        {
            std::size_t first = 0;
            std::size_t second = 0;
            pqss::Real distance_sq = pqss::k_zero;
            bool passed = false;
        };
        std::array<ChildPair, 2> children{};
        if (second.leaf() || (!first.leaf() && first.shape.size > second.shape.size))
        {
            children[0].first = first.left;
            children[0].second = second_index;
            children[1].first = first.right;
            children[1].second = second_index;
        }
        else
        {
            children[0].first = first_index;
            children[0].second = second.left;
            children[1].first = first_index;
            children[1].second = second.right;
        }
        for (ChildPair& child : children)
        {
            child.passed = TestPair(first_tree.nodes[child.first], second_tree.nodes[child.second],
                                    state, false, child.distance_sq);
        }
        if (children[1].distance_sq < children[0].distance_sq)
            std::swap(children[0], children[1]);
        for (const ChildPair& child : children)
        {
            if (!child.passed) continue;
            Traverse(first_tree, child.first, second_tree, child.second, state);
            if (state.result.closer_than_tolerance) return;
        }
    }
};

const char* analyticPrimitiveTypeName(const AnalyticPrimitiveType type)
{
    switch (type)
    {
    case AnalyticPrimitiveType::Sphere: return "sphere";
    case AnalyticPrimitiveType::Capsule: return "capsule";
    case AnalyticPrimitiveType::Rss: return "rss";
    case AnalyticPrimitiveType::Count: break;
    }
    return "unknown";
}

PrimitiveBvhPool::PrimitiveBvhPool() : m_impl(std::make_unique<Impl>()) {}
PrimitiveBvhPool::~PrimitiveBvhPool() = default;
PrimitiveBvhPool::PrimitiveBvhPool(PrimitiveBvhPool&&) noexcept = default;
PrimitiveBvhPool& PrimitiveBvhPool::operator=(PrimitiveBvhPool&&) noexcept = default;

void PrimitiveBvhPool::importReferencePool(const pqss::ModelPool& reference_pool,
                                           const std::vector<std::size_t>& object_ids)
{
    const auto& source_models = reference_pool.Models();
    if (source_models.size() != object_ids.size())
        throw std::invalid_argument("reference object ID count does not match PQSS pool size");
    for (std::size_t index = 0; index < source_models.size(); ++index)
    {
        const auto& source = source_models[index];
        if (source.Bvs().empty()) throw std::runtime_error("reference PQSS model has no BVH");
        ModelTree tree;
        tree.nodes.resize(source.Bvs().size());
        tree.stats.model_id = object_ids[index];
        tree.stats.source_triangles = source.OriginalTris().size();
        tree.stats.analyzed = false;
        tree.stats.internal_containment_certified = true;
        m_impl->ImportNode(tree, source, 0, pqss::Meident(), pqss::Veident(), 0);
        tree.stats.nodes = tree.nodes.size();
        m_impl->models[object_ids[index]] = std::move(tree);
    }
}

void PrimitiveBvhPool::loadAnalyzedPool(const std::filesystem::path& manifest,
                                        const pqss::BuildStrategy internal_fit_strategy,
                                        const std::filesystem::path& cache_directory)
{
    m_impl->cache_stats = {};
    std::ifstream stream(manifest);
    if (!stream) throw std::runtime_error("cannot open primitive pool manifest: " + manifest.string());
    std::string line;
    if (!std::getline(stream, line) || line != "PQSS_PRIMITIVE_BVH_POOL_V1")
        throw std::runtime_error("invalid primitive pool manifest header");
    if (!std::getline(stream, line) || line != "model_id\tmodel_file")
        throw std::runtime_error("invalid primitive pool manifest columns");
    while (std::getline(stream, line))
    {
        if (line.empty()) continue;
        const auto fields = Split(line, '\t');
        if (fields.size() != 2) throw std::runtime_error("invalid primitive pool manifest row: " + line);
        const std::size_t model_id = ParseIndex(fields[0], "model id");
        const std::filesystem::path model_path =
            manifest.parent_path() / std::filesystem::path(std::string(fields[1]));
        if (cache_directory.empty())
        {
            m_impl->models[model_id] =
                m_impl->LoadAnalyzedModel(model_id, model_path, internal_fit_strategy);
            continue;
        }

        ++m_impl->cache_stats.eligible_models;
        const std::uint64_t fingerprint = CacheFingerprint(model_path, internal_fit_strategy);
        const std::filesystem::path cache_path = cache_directory / (std::to_string(model_id) + ".bin");
        ModelTree tree;
        const auto load_started = Clock::now();
        const bool loaded = LoadTreeCache(cache_path, model_id, internal_fit_strategy, fingerprint, tree);
        m_impl->cache_stats.cache_load_ms +=
            std::chrono::duration<double, std::milli>(Clock::now() - load_started).count();
        if (loaded)
        {
            ++m_impl->cache_stats.cache_hits;
            m_impl->models[model_id] = std::move(tree);
            continue;
        }

        ++m_impl->cache_stats.cache_misses;
        const auto build_started = Clock::now();
        try
        {
            tree = m_impl->LoadAnalyzedModel(model_id, model_path, internal_fit_strategy);
        }
        catch (const std::exception& error)
        {
            throw std::runtime_error("primitive BVH model " + std::to_string(model_id) +
                                     " build failed: " + error.what());
        }
        m_impl->cache_stats.cache_build_ms +=
            std::chrono::duration<double, std::milli>(Clock::now() - build_started).count();
        const auto save_started = Clock::now();
        if (!SaveTreeCache(cache_path, model_id, internal_fit_strategy, fingerprint, tree))
            throw std::runtime_error("cannot save primitive BVH cache: " + cache_path.string());
        m_impl->cache_stats.cache_save_ms +=
            std::chrono::duration<double, std::milli>(Clock::now() - save_started).count();
        m_impl->models[model_id] = std::move(tree);
    }
}

PrimitiveBvhQueryResult PrimitiveBvhPool::query(const std::size_t model_id_1,
                                                const pqss::Mat3& R1,
                                                const pqss::Vec3& T1,
                                                const std::size_t model_id_2,
                                                const pqss::Mat3& R2,
                                                const pqss::Vec3& T2,
                                                const pqss::Real tolerance) const
{
    if (tolerance < pqss::k_zero) throw std::invalid_argument("collision tolerance must be non-negative");
    const auto first_found = m_impl->models.find(model_id_1);
    const auto second_found = m_impl->models.find(model_id_2);
    if (first_found == m_impl->models.end() || second_found == m_impl->models.end())
        throw std::out_of_range("query model ID is not loaded in PrimitiveBvhPool");
    const ModelTree& first = first_found->second;
    const ModelTree& second = second_found->second;
    Impl::QueryState state;
    state.object_R = pqss::MTxM(R1, R2);
    state.object_T = pqss::MTxVmV(R1, T2, T1);
    state.tolerance = tolerance;
    const auto started = Clock::now();
    pqss::Real root_distance_sq = pqss::k_zero;
    if (m_impl->TestPair(first.nodes[0], second.nodes[0], state, true, root_distance_sq))
        m_impl->Traverse(first, 0, second, 0, state);
    state.result.query_time_ms =
        std::chrono::duration<double, std::milli>(Clock::now() - started).count();
    return state.result;
}

bool PrimitiveBvhPool::hasModel(const std::size_t model_id) const
{
    return m_impl->models.contains(model_id);
}

const PrimitiveBvhModelStats& PrimitiveBvhPool::modelStats(const std::size_t model_id) const
{
    return m_impl->models.at(model_id).stats;
}

std::vector<std::size_t> PrimitiveBvhPool::modelIds() const
{
    std::vector<std::size_t> result;
    result.reserve(m_impl->models.size());
    for (const auto& [model_id, _] : m_impl->models)
        result.push_back(model_id);
    std::sort(result.begin(), result.end());
    return result;
}

const PrimitiveBvhCacheStats& PrimitiveBvhPool::cacheStats() const
{
    return m_impl->cache_stats;
}

void PrimitiveBvhPool::installIntoPqssPool(pqss::ModelPool& pool,
                                           const std::vector<std::size_t>& object_ids) const
{
    if (pool.PoolSize() != object_ids.size())
        throw std::invalid_argument("PQSS pool size does not match object ID mapping");
    for (std::size_t pool_index = 0; pool_index < object_ids.size(); ++pool_index)
    {
        const auto found = m_impl->models.find(object_ids[pool_index]);
        if (found == m_impl->models.end() || !found->second.stats.analyzed)
            continue;
        const ModelTree& tree = found->second;
        std::vector<pqss::RssBvhNode> nodes;
        nodes.reserve(tree.nodes.size());
        for (const Node& source : tree.nodes)
        {
            pqss::RssBvhNode node;
            node.R = source.shape.R;
            node.T = source.shape.T;
            node.lengths = source.shape.lengths;
            node.radius = source.shape.radius;
            node.size = source.shape.size;
            node.left = source.left;
            node.right = source.right;
            node.leaf = source.leaf();
            nodes.push_back(node);
        }
        const pqss::Return status = pool.ReplaceModelBvhWithRssLeaves(pool_index, nodes);
        if (status != pqss::Return::Success)
            throw std::runtime_error("cannot install primitive BVH for model " +
                                     std::to_string(object_ids[pool_index]) +
                                     "; Return=" + std::to_string(static_cast<int>(status)));
    }
}

} // namespace pqss_proxy_mesh
