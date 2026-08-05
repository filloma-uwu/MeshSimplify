#include "pqss_proxy_mesh/independent_primitive_pool.hpp"

#include <algorithm>
#include <array>
#include <chrono>
#include <charconv>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <limits>
#include <stdexcept>
#include <string_view>
#include <unordered_set>
#include <utility>
#include <vector>

#include "pqss/dist.hpp"
#include "pqss/utils/mat_vec.hpp"

namespace pqss_proxy_mesh
{
namespace
{

using Clock = std::chrono::steady_clock;

struct RssNode
{
    pqss::Mat3 R{};
    pqss::Vec3 T{};
    std::array<pqss::Real, 2> lengths{};
    pqss::Real radius = pqss::k_zero;
    pqss::Real size = pqss::k_zero;
    std::size_t first_child = 0;
    std::size_t triangle = 0;
    std::uint32_t global_id = 0;
    bool leaf = false;
    std::vector<std::uint32_t> containers;
    pqss::Vec3 bounds_min{};
    pqss::Vec3 bounds_max{};
};

struct PhysicalTree
{
    std::size_t physical_model = 0;
    std::size_t logical_model = 0;
    std::vector<RssNode> nodes;
};

struct PairKey
{
    std::uint32_t first = 0;
    std::uint32_t second = 0;

    bool operator==(const PairKey&) const = default;
};

struct PairHash
{
    std::size_t operator()(const PairKey pair) const noexcept
    {
        const std::uint64_t packed = (static_cast<std::uint64_t>(pair.first) << 32) | pair.second;
        return std::hash<std::uint64_t>{}(packed);
    }
};

std::vector<std::string_view> SplitTabs(const std::string& line)
{
    std::vector<std::string_view> result;
    std::size_t begin = 0;
    while (begin <= line.size())
    {
        const std::size_t end = line.find('\t', begin);
        result.emplace_back(line.data() + begin,
                            (end == std::string::npos ? line.size() : end) - begin);
        if (end == std::string::npos)
            break;
        begin = end + 1;
    }
    return result;
}

std::size_t ParseIndex(const std::string_view value, const std::string_view field)
{
    std::size_t result = 0;
    const char* begin = value.data();
    const char* end = value.data() + value.size();
    const auto parsed = std::from_chars(begin, end, result);
    if (parsed.ec != std::errc{} || parsed.ptr != end)
        throw std::runtime_error("invalid " + std::string(field) + ": " + std::string(value));
    return result;
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

pqss::Vec3 RectangleCorner(const RssNode& node, const int x, const int y)
{
    pqss::Vec3 result = node.T;
    for (int row = 0; row < 3; ++row)
    {
        result[row] += node.R[row][0] * node.lengths[0] * static_cast<pqss::Real>(x);
        result[row] += node.R[row][1] * node.lengths[1] * static_cast<pqss::Real>(y);
    }
    return result;
}

void ComputeBounds(RssNode& node)
{
    const pqss::Real limit = std::numeric_limits<pqss::Real>::max();
    node.bounds_min = {limit, limit, limit};
    node.bounds_max = {-limit, -limit, -limit};
    for (int x = 0; x <= 1; ++x)
    {
        for (int y = 0; y <= 1; ++y)
        {
            const pqss::Vec3 corner = RectangleCorner(node, x, y);
            for (int axis = 0; axis < 3; ++axis)
            {
                node.bounds_min[axis] = std::min(node.bounds_min[axis], corner[axis] - node.radius);
                node.bounds_max[axis] = std::max(node.bounds_max[axis], corner[axis] + node.radius);
            }
        }
    }
}

pqss::Real PointRectangleDistanceSq(const pqss::Vec3& point, const RssNode& rectangle)
{
    const pqss::Vec3 local = pqss::MTxVmV(rectangle.R, point, rectangle.T);
    const pqss::Real x = std::clamp(local[0], pqss::k_zero, rectangle.lengths[0]);
    const pqss::Real y = std::clamp(local[1], pqss::k_zero, rectangle.lengths[1]);
    const pqss::Vec3 delta = {local[0] - x, local[1] - y, local[2]};
    return pqss::VdotV(delta, delta);
}

bool RssContains(const RssNode& outer, const RssNode& inner)
{
    const pqss::Real scale = std::max({pqss::Real{1}, outer.lengths[0], outer.lengths[1],
                                       outer.radius, inner.lengths[0], inner.lengths[1], inner.radius});
    const pqss::Real margin = scale * static_cast<pqss::Real>(1.0e-10);
    const pqss::Real available = outer.radius - inner.radius;
    if (available <= margin)
        return false;
    for (int axis = 0; axis < 3; ++axis)
    {
        if (outer.bounds_min[axis] > inner.bounds_min[axis] + margin ||
            outer.bounds_max[axis] < inner.bounds_max[axis] - margin)
            return false;
    }
    const pqss::Real certified = available - margin;
    const pqss::Real limit_sq = certified * certified;
    for (int x = 0; x <= 1; ++x)
        for (int y = 0; y <= 1; ++y)
            if (PointRectangleDistanceSq(RectangleCorner(inner, x, y), outer) > limit_sq)
                return false;
    return true;
}

pqss::Real RssRectangleDistanceSq(const RssNode& first,
                                  const RssNode& second,
                                  const pqss::Mat3& object_R,
                                  const pqss::Vec3& object_T)
{
    const pqss::Mat3 relative_R = pqss::MTxM(first.R, pqss::MxM(object_R, second.R));
    const pqss::Vec3 relative_T =
        pqss::MTxVmV(first.R, pqss::MxVpV(object_R, second.T, object_T), first.T);
    return pqss::detail::RectDistSq(relative_R, relative_T, first.lengths, second.lengths);
}

bool RssSeparated(const RssNode& first,
                  const RssNode& second,
                  const pqss::Mat3& object_R,
                  const pqss::Vec3& object_T,
                  const pqss::Real tolerance,
                  pqss::Real* rectangle_distance_sq = nullptr)
{
    const pqss::Real distance_sq = RssRectangleDistanceSq(first, second, object_R, object_T);
    if (rectangle_distance_sq != nullptr)
        *rectangle_distance_sq = distance_sq;
    const pqss::Real threshold = tolerance + first.radius + second.radius;
    return distance_sq > threshold * threshold;
}

} // namespace

struct IndependentPrimitivePool::Impl
{
    pqss::ModelPool pool;
    std::vector<std::vector<std::size_t>> logical_trees;
    std::vector<PhysicalTree> trees;
    std::vector<std::size_t> physical_to_tree;
    std::size_t next_physical_model = 0;
    ContainmentStats stats;
    bool base_built = false;
    bool finalized = false;

    Impl() : logical_trees(36) {}

    static void Require(const pqss::Return result, const std::string& operation)
    {
        if (result != pqss::Return::Success)
            throw std::runtime_error(operation + " failed with Return=" +
                                     std::to_string(static_cast<int>(result)));
    }

    void LoadPhysical(const std::size_t logical_model,
                      const std::filesystem::path& path,
                      const pqss::BuildStrategy strategy)
    {
        if (logical_model >= logical_trees.size())
            throw std::out_of_range("logical model index is out of range");
        const std::size_t physical_model = next_physical_model++;
        Require(pool.LoadModelFromFile(path.string(), strategy), "LoadModelFromFile(" + path.string() + ")");
        PhysicalTree tree;
        tree.physical_model = physical_model;
        tree.logical_model = logical_model;
        logical_trees[logical_model].push_back(trees.size());
        trees.push_back(std::move(tree));
        finalized = false;
    }

    void BuildTreeNodes(PhysicalTree& tree,
                        const pqss::build::Model& model,
                        const std::size_t node_index,
                        const pqss::Mat3& parent_R,
                        const pqss::Vec3& parent_T,
                        const std::vector<std::uint32_t>& ancestors,
                        std::uint32_t& next_global)
    {
        const pqss::build::BV& source = model.Bvs().at(node_index);
        RssNode& node = tree.nodes.at(node_index);
        node.R = ComposeRotation(parent_R, source.R());
        node.T = ComposeTranslation(parent_R, parent_T, source.Tr());
        node.lengths = source.L();
        node.radius = source.Radius();
        node.size = source.Size();
        node.leaf = source.IsLeaf();
        node.first_child = source.IsLeaf() ? 0 : source.FirstChild();
        node.triangle = source.IsLeaf() ? source.TriIndex() : 0;
        if (next_global == std::numeric_limits<std::uint32_t>::max())
            throw std::overflow_error("too many RSS nodes for 32-bit query keys");
        node.global_id = next_global++;
        node.containers = ancestors;
        node.containers.push_back(node.global_id);
        stats.tree_ancestor_links += ancestors.size();
        ComputeBounds(node);

        if (!node.leaf)
        {
            BuildTreeNodes(tree, model, node.first_child, node.R, node.T, node.containers, next_global);
            BuildTreeNodes(tree, model, node.first_child + 1, node.R, node.T, node.containers, next_global);
        }
    }

    void Finalize()
    {
        if (!base_built)
            throw std::logic_error("base pool must be built before containment finalization");
        stats = {};
        const auto& models = pool.Models();
        if (models.size() != trees.size())
            throw std::runtime_error("physical model mapping does not match PQSS pool");
        physical_to_tree.assign(models.size(), std::numeric_limits<std::size_t>::max());
        std::uint32_t next_global = 0;
        const pqss::Mat3 identity_R = pqss::Meident();
        const pqss::Vec3 identity_T = pqss::Veident();
        for (std::size_t tree_index = 0; tree_index < trees.size(); ++tree_index)
        {
            PhysicalTree& tree = trees[tree_index];
            const auto& model = models.at(tree.physical_model);
            tree.nodes.assign(model.Bvs().size(), RssNode{});
            BuildTreeNodes(tree, model, 0, identity_R, identity_T, {}, next_global);
            physical_to_tree.at(tree.physical_model) = tree_index;
        }

        for (const auto& group : logical_trees)
        {
            for (const std::size_t target_tree_index : group)
            {
                PhysicalTree& target_tree = trees[target_tree_index];
                for (RssNode& target : target_tree.nodes)
                {
                    for (const std::size_t container_tree_index : group)
                    {
                        if (container_tree_index == target_tree_index)
                            continue;
                        for (const RssNode& candidate : trees[container_tree_index].nodes)
                        {
                            if (RssContains(candidate, target))
                            {
                                target.containers.push_back(candidate.global_id);
                                ++stats.cross_tree_containment_links;
                            }
                        }
                    }
                    std::sort(target.containers.begin(), target.containers.end());
                    target.containers.erase(
                        std::unique(target.containers.begin(), target.containers.end()),
                        target.containers.end());
                }
            }
        }
        stats.physical_models = trees.size();
        stats.bvh_nodes = next_global;
        finalized = true;
    }

    struct SeparationAxis
    {
        pqss::Vec3 direction{};
    };

    struct QueryState
    {
        IndependentQueryResult result;
        pqss::Mat3 object_R{};
        pqss::Vec3 object_T{};
        pqss::Real tolerance = pqss::k_zero;
        DistanceFilterMode filter_mode = DistanceFilterMode::Hybrid;
        std::unordered_set<PairKey, PairHash> pruned;
        std::vector<SeparationAxis> separating_axes;
    };

    static bool UsesContainment(const DistanceFilterMode mode)
    {
        return mode == DistanceFilterMode::Containment || mode == DistanceFilterMode::Hybrid;
    }

    static bool UsesSeparatingAxes(const DistanceFilterMode mode)
    {
        return mode == DistanceFilterMode::SeparatingAxis || mode == DistanceFilterMode::Hybrid;
    }

    static pqss::Vec3 Column(const pqss::Mat3& matrix, const int column)
    {
        return {matrix[0][column], matrix[1][column], matrix[2][column]};
    }

    static pqss::Vec3 Cross(const pqss::Vec3& first, const pqss::Vec3& second)
    {
        return {
            first[1] * second[2] - first[2] * second[1],
            first[2] * second[0] - first[0] * second[2],
            first[0] * second[1] - first[1] * second[0],
        };
    }

    static std::pair<pqss::Real, pqss::Real> ProjectionInterval(
        const RssNode& node,
        const pqss::Vec3& direction,
        const pqss::Mat3* object_R = nullptr,
        const pqss::Vec3* object_T = nullptr)
    {
        const pqss::Mat3 frame = object_R == nullptr ? node.R : pqss::MxM(*object_R, node.R);
        const pqss::Vec3 origin = object_R == nullptr
            ? node.T
            : pqss::MxVpV(*object_R, node.T, *object_T);
        const pqss::Real first_extent = pqss::VdotV(direction, Column(frame, 0)) * node.lengths[0];
        const pqss::Real second_extent = pqss::VdotV(direction, Column(frame, 1)) * node.lengths[1];
        const pqss::Real base = pqss::VdotV(direction, origin);
        const pqss::Real minimum = base + std::min(pqss::k_zero, first_extent) +
                                   std::min(pqss::k_zero, second_extent) - node.radius;
        const pqss::Real maximum = base + std::max(pqss::k_zero, first_extent) +
                                   std::max(pqss::k_zero, second_extent) + node.radius;
        return {minimum, maximum};
    }

    static pqss::Real ProjectionGap(const RssNode& first,
                                    const RssNode& second,
                                    const pqss::Vec3& direction,
                                    const QueryState& state)
    {
        const auto first_interval = ProjectionInterval(first, direction);
        const auto second_interval = ProjectionInterval(
            second, direction, &state.object_R, &state.object_T);
        return std::max({pqss::k_zero,
                         second_interval.first - first_interval.second,
                         first_interval.first - second_interval.second});
    }

    bool ContainmentFiltered(const RssNode& first,
                             const RssNode& second,
                             QueryState& state) const
    {
        if (!UsesContainment(state.filter_mode))
            return false;
        for (const std::uint32_t first_container : first.containers)
        {
            for (const std::uint32_t second_container : second.containers)
            {
                ++state.result.containment_filter_checks;
                if (state.pruned.contains({first_container, second_container}))
                {
                    ++state.result.containment_filter_skips;
                    return true;
                }
            }
        }
        return false;
    }

    bool AxisFiltered(const RssNode& first,
                      const RssNode& second,
                      QueryState& state) const
    {
        if (!UsesSeparatingAxes(state.filter_mode))
            return false;
        for (const SeparationAxis& axis : state.separating_axes)
        {
            ++state.result.axis_filter_checks;
            if (ProjectionGap(first, second, axis.direction, state) > state.tolerance)
            {
                ++state.result.axis_filter_skips;
                return true;
            }
        }
        return false;
    }

    bool DistanceFiltered(const RssNode& first,
                          const RssNode& second,
                          QueryState& state) const
    {
        // Projection tests are cheaper and usually have much smaller state
        // than the Cartesian product of containment ancestors.
        const bool filtered = AxisFiltered(first, second, state) ||
                              ContainmentFiltered(first, second, state);
        if (filtered)
            ++state.result.distance_filter_skips;
        return filtered;
    }

    bool TryRememberAxis(const RssNode& first,
                         const RssNode& second,
                         pqss::Vec3 direction,
                         QueryState& state) const
    {
        constexpr std::size_t max_cached_axes = 4;
        if (state.separating_axes.size() >= max_cached_axes)
            return false;
        const pqss::Real length_sq = pqss::VdotV(direction, direction);
        if (length_sq <= static_cast<pqss::Real>(1.0e-20))
            return false;
        const pqss::Real inverse_length = pqss::k_one / std::sqrt(length_sq);
        for (pqss::Real& value : direction)
            value *= inverse_length;
        for (int axis = 0; axis < 3; ++axis)
        {
            if (std::abs(direction[axis]) <= static_cast<pqss::Real>(1.0e-12))
                continue;
            if (direction[axis] < pqss::k_zero)
                for (pqss::Real& value : direction) value = -value;
            break;
        }

        ++state.result.axis_generation_checks;
        if (ProjectionGap(first, second, direction, state) <= state.tolerance)
            return false;
        for (const SeparationAxis& existing : state.separating_axes)
            if (std::abs(pqss::VdotV(existing.direction, direction)) >
                static_cast<pqss::Real>(0.9995))
                return false;
        state.separating_axes.push_back({direction});
        return true;
    }

    void RememberSeparationAxes(const RssNode& first,
                                const RssNode& second,
                                QueryState& state) const
    {
        if (!UsesSeparatingAxes(state.filter_mode))
            return;
        constexpr std::size_t max_cached_axes = 4;
        if (state.separating_axes.size() >= max_cached_axes)
            return;

        const pqss::Mat3 second_frame = pqss::MxM(state.object_R, second.R);
        pqss::Vec3 first_center = first.T;
        pqss::Vec3 second_center = pqss::MxVpV(state.object_R, second.T, state.object_T);
        for (int row = 0; row < 3; ++row)
        {
            first_center[row] += static_cast<pqss::Real>(0.5) *
                (first.R[row][0] * first.lengths[0] + first.R[row][1] * first.lengths[1]);
            second_center[row] += static_cast<pqss::Real>(0.5) *
                (second_frame[row][0] * second.lengths[0] + second_frame[row][1] * second.lengths[1]);
        }
        pqss::Vec3 center_delta{};
        for (int axis = 0; axis < 3; ++axis)
            center_delta[axis] = second_center[axis] - first_center[axis];
        if (TryRememberAxis(first, second, center_delta, state))
            return;

        std::array<pqss::Vec3, 3> first_axes = {
            Column(first.R, 0), Column(first.R, 1), Column(first.R, 2)};
        std::array<pqss::Vec3, 3> second_axes = {
            Column(second_frame, 0), Column(second_frame, 1), Column(second_frame, 2)};
        for (const pqss::Vec3& axis : first_axes)
            if (TryRememberAxis(first, second, axis, state)) return;
        for (const pqss::Vec3& axis : second_axes)
            if (TryRememberAxis(first, second, axis, state)) return;
        for (const pqss::Vec3& first_axis : first_axes)
            for (const pqss::Vec3& second_axis : second_axes)
                if (TryRememberAxis(first, second, Cross(first_axis, second_axis), state)) return;
    }

    bool TestAndRecord(const RssNode& first,
                       const RssNode& second,
                       QueryState& state,
                       const bool primitive_root,
                       pqss::Real* rectangle_distance_sq = nullptr) const
    {
        if (DistanceFiltered(first, second, state))
            return false;
        if (primitive_root)
            ++state.result.primitive_pair_tests;
        else
            ++state.result.internal_bv_tests;
        if (RssSeparated(first, second, state.object_R, state.object_T,
                         state.tolerance, rectangle_distance_sq))
        {
            if (UsesContainment(state.filter_mode))
                state.pruned.insert({first.global_id, second.global_id});
            RememberSeparationAxes(first, second, state);
            return false;
        }
        return true;
    }

    void TraversePassed(const PhysicalTree& first_tree,
                        const std::size_t first_index,
                        const PhysicalTree& second_tree,
                        const std::size_t second_index,
                        QueryState& state) const
    {
        if (state.result.closer_than_tolerance)
            return;
        const RssNode& first = first_tree.nodes[first_index];
        const RssNode& second = second_tree.nodes[second_index];
        if (first.leaf && second.leaf)
        {
            ++state.result.triangle_tests;
            const auto& first_tri = pool.Models()[first_tree.physical_model].OriginalTris().at(first.triangle);
            const auto& second_tri = pool.Models()[second_tree.physical_model].OriginalTris().at(second.triangle);
            const std::array<pqss::Vec3, 3> first_points = {first_tri[0], first_tri[1], first_tri[2]};
            const std::array<pqss::Vec3, 3> second_points = {
                pqss::MxVpV(state.object_R, second_tri[0], state.object_T),
                pqss::MxVpV(state.object_R, second_tri[1], state.object_T),
                pqss::MxVpV(state.object_R, second_tri[2], state.object_T),
            };
            pqss::Vec3 p{};
            pqss::Vec3 q{};
            const pqss::Real distance_sq = pqss::detail::TriDistSq(p, q, first_points, second_points);
            state.result.closer_than_tolerance =
                distance_sq <= state.tolerance * state.tolerance;
            return;
        }

        struct ChildPair
        {
            std::size_t first = 0;
            std::size_t second = 0;
            pqss::Real rectangle_distance_sq = pqss::k_zero;
        };
        std::array<ChildPair, 2> children{};
        if (second.leaf || (!first.leaf && first.size > second.size))
        {
            children[0] = {first.first_child, second_index, pqss::k_zero};
            children[1] = {first.first_child + 1, second_index, pqss::k_zero};
        }
        else
        {
            children[0] = {first_index, second.first_child, pqss::k_zero};
            children[1] = {first_index, second.first_child + 1, pqss::k_zero};
        }

        std::array<bool, 2> passed{};
        for (std::size_t index = 0; index < children.size(); ++index)
        {
            const RssNode& child_first = first_tree.nodes[children[index].first];
            const RssNode& child_second = second_tree.nodes[children[index].second];
            passed[index] = TestAndRecord(child_first, child_second, state, false,
                                          &children[index].rectangle_distance_sq);
        }
        if (children[1].rectangle_distance_sq < children[0].rectangle_distance_sq)
        {
            std::swap(children[0], children[1]);
            std::swap(passed[0], passed[1]);
        }
        for (std::size_t index = 0; index < children.size(); ++index)
        {
            if (!passed[index])
                continue;
            TraversePassed(first_tree, children[index].first,
                           second_tree, children[index].second, state);
            if (state.result.closer_than_tolerance)
                return;
        }
    }

    IndependentQueryResult Query(const std::size_t logical_model_1,
                                 const pqss::Mat3& R1,
                                 const pqss::Vec3& T1,
                                 const std::size_t logical_model_2,
                                 const pqss::Mat3& R2,
                                 const pqss::Vec3& T2,
                                 const pqss::Real tolerance,
                                 const DistanceFilterMode filter_mode) const
    {
        if (!finalized)
            throw std::logic_error("containment metadata has not been finalized");
        if (logical_model_1 >= logical_trees.size() || logical_model_2 >= logical_trees.size())
            throw std::out_of_range("logical model index is out of range");
        if (logical_trees[logical_model_1].empty() || logical_trees[logical_model_2].empty())
            throw std::runtime_error("logical model has no physical primitives");

        QueryState state;
        state.object_R = pqss::MTxM(R1, R2);
        state.object_T = pqss::MTxVmV(R1, T2, T1);
        state.tolerance = tolerance;
        state.filter_mode = filter_mode;
        const auto started = Clock::now();

        struct RootPair
        {
            std::size_t first_tree = 0;
            std::size_t second_tree = 0;
            pqss::Real size = pqss::k_zero;
        };
        std::vector<RootPair> roots;
        for (const std::size_t first_tree : logical_trees[logical_model_1])
            for (const std::size_t second_tree : logical_trees[logical_model_2])
                roots.push_back({first_tree, second_tree,
                                 trees[first_tree].nodes[0].size + trees[second_tree].nodes[0].size});
        std::stable_sort(roots.begin(), roots.end(),
                         [](const RootPair& first, const RootPair& second)
                         {
                             return first.size > second.size;
                         });

        for (const RootPair& root : roots)
        {
            const PhysicalTree& first_tree = trees[root.first_tree];
            const PhysicalTree& second_tree = trees[root.second_tree];
            if (!TestAndRecord(first_tree.nodes[0], second_tree.nodes[0], state, true))
                continue;
            TraversePassed(first_tree, 0, second_tree, 0, state);
            if (state.result.closer_than_tolerance)
                break;
        }
        state.result.query_time_ms =
            std::chrono::duration<double, std::milli>(Clock::now() - started).count();
        state.result.cached_separating_axes = state.separating_axes.size();
        return state.result;
    }
};

const char* distanceFilterModeName(const DistanceFilterMode mode)
{
    switch (mode)
    {
    case DistanceFilterMode::None: return "none";
    case DistanceFilterMode::Containment: return "containment";
    case DistanceFilterMode::SeparatingAxis: return "separating-axis";
    case DistanceFilterMode::Hybrid: return "hybrid";
    }
    return "unknown";
}

IndependentPrimitivePool::IndependentPrimitivePool() : m_impl(std::make_unique<Impl>()) {}
IndependentPrimitivePool::~IndependentPrimitivePool() = default;
IndependentPrimitivePool::IndependentPrimitivePool(IndependentPrimitivePool&&) noexcept = default;
IndependentPrimitivePool& IndependentPrimitivePool::operator=(IndependentPrimitivePool&&) noexcept = default;

void IndependentPrimitivePool::setBuildCacheTag(const std::string& cache_tag)
{
    Impl::Require(m_impl->pool.SetBuildCacheTag(cache_tag), "SetBuildCacheTag");
}

void IndependentPrimitivePool::loadBaseManifest(const std::filesystem::path& manifest,
                                                const pqss::BuildStrategy strategy)
{
    if (m_impl->base_built)
        throw std::logic_error("base pool is already built");
    std::ifstream stream(manifest);
    if (!stream)
        throw std::runtime_error("cannot open primitive manifest: " + manifest.string());
    const std::filesystem::path root = manifest.parent_path();
    std::string line;
    if (!std::getline(stream, line) || line != "logical_model_id\tpart_index\trelative_path\tselected")
        throw std::runtime_error("invalid primitive manifest header");
    std::array<std::size_t, 23> expected_part{};
    while (std::getline(stream, line))
    {
        if (line.empty())
            continue;
        const auto fields = SplitTabs(line);
        if (fields.size() != 4)
            throw std::runtime_error("invalid primitive manifest row: " + line);
        const std::size_t object_id = ParseIndex(fields[0], "logical_model_id");
        const std::size_t part_index = ParseIndex(fields[1], "part_index");
        if (object_id < 1 || object_id > 23)
            throw std::runtime_error("base logical_model_id must be in 1-23");
        if (part_index != expected_part[object_id - 1]++)
            throw std::runtime_error("part indices must be contiguous for logical model " +
                                     std::to_string(object_id));
        m_impl->LoadPhysical(object_id - 1, root / std::filesystem::path(fields[2]), strategy);
    }
    for (std::size_t logical = 0; logical < expected_part.size(); ++logical)
        if (expected_part[logical] == 0)
            throw std::runtime_error("manifest is missing logical model " + std::to_string(logical + 1));
}

void IndependentPrimitivePool::endBasePool()
{
    Impl::Require(m_impl->pool.EndPool(), "EndPool");
    m_impl->base_built = true;
}

void IndependentPrimitivePool::appendLogicalModel(const std::size_t logical_pool_index,
                                                  const std::filesystem::path& model,
                                                  const pqss::BuildStrategy strategy)
{
    if (!m_impl->base_built)
        throw std::logic_error("base pool must be built before appending models");
    m_impl->LoadPhysical(logical_pool_index, model, strategy);
}

void IndependentPrimitivePool::finalizeContainment()
{
    m_impl->Finalize();
}

IndependentQueryResult IndependentPrimitivePool::query(const std::size_t logical_model_1,
                                                       const pqss::Mat3& R1,
                                                       const pqss::Vec3& T1,
                                                       const std::size_t logical_model_2,
                                                       const pqss::Mat3& R2,
                                                       const pqss::Vec3& T2,
                                                       const pqss::Real tolerance,
                                                       const DistanceFilterMode filter_mode) const
{
    return m_impl->Query(logical_model_1, R1, T1, logical_model_2, R2, T2,
                         tolerance, filter_mode);
}

const pqss::ModelPool& IndependentPrimitivePool::physicalPool() const
{
    return m_impl->pool;
}

const ContainmentStats& IndependentPrimitivePool::containmentStats() const
{
    return m_impl->stats;
}

std::size_t IndependentPrimitivePool::logicalModelCount() const
{
    return m_impl->logical_trees.size();
}

std::size_t IndependentPrimitivePool::primitiveCount(const std::size_t logical_model) const
{
    if (logical_model >= m_impl->logical_trees.size())
        throw std::out_of_range("logical model index is out of range");
    return m_impl->logical_trees[logical_model].size();
}

std::size_t IndependentPrimitivePool::physicalModelId(const std::size_t logical_model,
                                                      const std::size_t primitive_index) const
{
    if (logical_model >= m_impl->logical_trees.size())
        throw std::out_of_range("logical model index is out of range");
    const auto& group = m_impl->logical_trees[logical_model];
    if (primitive_index >= group.size())
        throw std::out_of_range("primitive index is out of range");
    return m_impl->trees[group[primitive_index]].physical_model;
}

} // namespace pqss_proxy_mesh
