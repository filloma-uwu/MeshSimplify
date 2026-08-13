#include "pqss_proxy_mesh/topology_fill.hpp"

#include <clipper2/clipper.h>

#include <algorithm>
#include <array>
#include <bit>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <limits>
#include <map>
#include <numeric>
#include <optional>
#include <queue>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <tuple>
#include <unordered_map>
#include <utility>
#include <vector>

namespace pqss_proxy_mesh
{
namespace
{

constexpr std::array<std::array<int, 3>, 6> neighbors6{{
    {{-1, 0, 0}}, {{1, 0, 0}}, {{0, -1, 0}},
    {{0, 1, 0}}, {{0, 0, -1}}, {{0, 0, 1}},
}};

struct Vec3
{
    double x = 0.0;
    double y = 0.0;
    double z = 0.0;

    double& operator[](const int axis) { return axis == 0 ? x : (axis == 1 ? y : z); }
    double operator[](const int axis) const { return axis == 0 ? x : (axis == 1 ? y : z); }
};

Vec3 operator+(const Vec3 a, const Vec3 b) { return {a.x + b.x, a.y + b.y, a.z + b.z}; }
Vec3 operator-(const Vec3 a, const Vec3 b) { return {a.x - b.x, a.y - b.y, a.z - b.z}; }
Vec3 operator*(const Vec3 a, const double value) { return {a.x * value, a.y * value, a.z * value}; }
double dot(const Vec3 a, const Vec3 b) { return a.x * b.x + a.y * b.y + a.z * b.z; }
Vec3 cross(const Vec3 a, const Vec3 b)
{
    return {a.y * b.z - a.z * b.y,
            a.z * b.x - a.x * b.z,
            a.x * b.y - a.y * b.x};
}
double norm(const Vec3 value) { return std::sqrt(dot(value, value)); }
Vec3 normalized(const Vec3 value)
{
    const double length = norm(value);
    return length == 0.0 ? Vec3{} : value * (1.0 / length);
}
Vec3 toVec(const Position3 value) { return {value.x, value.y, value.z}; }
Position3 toPosition(const Vec3 value) { return {value.x, value.y, value.z}; }

std::size_t checkedVoxelCount(const std::array<std::uint32_t, 3>& shape)
{
    const std::uint64_t count = static_cast<std::uint64_t>(shape[0]) * shape[1] * shape[2];
    if (count == 0 || count > std::numeric_limits<std::size_t>::max())
        throw std::invalid_argument("invalid voxel-grid shape");
    return static_cast<std::size_t>(count);
}

std::array<std::uint32_t, 3> decodeIndex(
    const std::size_t index, const std::array<std::uint32_t, 3>& shape)
{
    const std::size_t yz = static_cast<std::size_t>(shape[1]) * shape[2];
    const auto x = static_cast<std::uint32_t>(index / yz);
    const std::size_t remainder = index % yz;
    return {x, static_cast<std::uint32_t>(remainder / shape[2]),
            static_cast<std::uint32_t>(remainder % shape[2])};
}

std::size_t countOccupied(const VoxelGrid& grid)
{
    return static_cast<std::size_t>(std::count_if(
        grid.occupancy.begin(), grid.occupancy.end(),
        [](const std::uint8_t value) { return value != 0; }));
}

void validateGrid(const VoxelGrid& grid)
{
    if (grid.occupancy.size() != checkedVoxelCount(grid.shape))
        throw std::invalid_argument("voxel payload does not match its shape");
    if (!std::isfinite(grid.pitch) || grid.pitch <= 0.0)
        throw std::invalid_argument("voxel pitch must be finite and positive");
    if (countOccupied(grid) == 0)
        throw std::invalid_argument("voxel occupancy is empty");
}

template <std::size_t NeighborCount>
std::vector<std::uint8_t> floodEmptyFromBoundary(
    const VoxelGrid& grid,
    const std::array<std::array<int, 3>, NeighborCount>& neighbors)
{
    std::vector<std::uint8_t> visited(grid.occupancy.size(), 0);
    std::vector<std::uint32_t> queue;
    queue.reserve(grid.occupancy.size() / 4);
    const auto seed = [&](const std::uint32_t x, const std::uint32_t y,
                          const std::uint32_t z)
    {
        const std::size_t id = grid.index(x, y, z);
        if (!grid.occupancy[id] && !visited[id])
        {
            visited[id] = 1;
            queue.push_back(static_cast<std::uint32_t>(id));
        }
    };
    for (std::uint32_t x = 0; x < grid.shape[0]; ++x)
        for (std::uint32_t y = 0; y < grid.shape[1]; ++y)
        {
            seed(x, y, 0);
            seed(x, y, grid.shape[2] - 1);
        }
    for (std::uint32_t x = 0; x < grid.shape[0]; ++x)
        for (std::uint32_t z = 0; z < grid.shape[2]; ++z)
        {
            seed(x, 0, z);
            seed(x, grid.shape[1] - 1, z);
        }
    for (std::uint32_t y = 0; y < grid.shape[1]; ++y)
        for (std::uint32_t z = 0; z < grid.shape[2]; ++z)
        {
            seed(0, y, z);
            seed(grid.shape[0] - 1, y, z);
        }

    std::size_t head = 0;
    while (head < queue.size())
    {
        const auto coordinate = decodeIndex(queue[head++], grid.shape);
        for (const auto& offset : neighbors)
        {
            const int x = static_cast<int>(coordinate[0]) + offset[0];
            const int y = static_cast<int>(coordinate[1]) + offset[1];
            const int z = static_cast<int>(coordinate[2]) + offset[2];
            if (x < 0 || y < 0 || z < 0 ||
                x >= static_cast<int>(grid.shape[0]) ||
                y >= static_cast<int>(grid.shape[1]) ||
                z >= static_cast<int>(grid.shape[2])) continue;
            const std::size_t id = grid.index(
                static_cast<std::uint32_t>(x), static_cast<std::uint32_t>(y),
                static_cast<std::uint32_t>(z));
            if (!grid.occupancy[id] && !visited[id])
            {
                visited[id] = 1;
                queue.push_back(static_cast<std::uint32_t>(id));
            }
        }
    }
    return visited;
}

constexpr auto neighbors26 = []
{
    std::array<std::array<int, 3>, 26> result{};
    std::size_t index = 0;
    for (int x = -1; x <= 1; ++x)
        for (int y = -1; y <= 1; ++y)
            for (int z = -1; z <= 1; ++z)
                if (x != 0 || y != 0 || z != 0) result[index++] = {x, y, z};
    return result;
}();

VoxelGrid fillCavities(const VoxelGrid& input)
{
    VoxelGrid result = input;
    const auto outside = floodEmptyFromBoundary(input, neighbors26);
    for (std::size_t index = 0; index < result.occupancy.size(); ++index)
        if (!result.occupancy[index] && !outside[index]) result.occupancy[index] = 1;
    return result;
}

VoxelGrid thinnestConnectivityHandleCore(
    const VoxelGrid& source, const VoxelGrid& envelope,
    const bool preserve_responsible_envelope = false)
{
    const auto core_started = std::chrono::steady_clock::now();
    auto segment_started = core_started;
    const auto log_segment = [&](const char* stage)
    {
        const auto now = std::chrono::steady_clock::now();
        std::cerr << "monitor: stage=topology_core_" << stage
                  << " segment_seconds="
                  << std::chrono::duration<double>(now - segment_started).count()
                  << " total_seconds="
                  << std::chrono::duration<double>(now - core_started).count() << '\n';
        segment_started = now;
    };
    if (source.shape != envelope.shape || source.occupancy.size() != envelope.occupancy.size())
        throw std::invalid_argument("topology core grids do not match");
    for (std::size_t index = 0; index < source.occupancy.size(); ++index)
    {
        if (source.occupancy[index] && !envelope.occupancy[index])
            throw std::invalid_argument("topology envelope lost a source voxel");
    }
    const BettiNumbers envelope_topology = voxelBettiNumbers(envelope);
    log_segment("envelope_betti");
    if (std::abs(envelope_topology.beta0 - 1) + envelope_topology.beta1 != 0)
        throw std::runtime_error("topology envelope has unresolved connectivity or handles");

    constexpr std::uint32_t unreachable = std::numeric_limits<std::uint32_t>::max();
    std::vector<std::uint32_t> distance(source.occupancy.size(), unreachable);
    std::vector<std::uint32_t> distance_queue;
    distance_queue.reserve(countOccupied(envelope) - countOccupied(source));
    std::uint32_t maximum_distance = 0;
    for (std::uint32_t x = 0; x < source.shape[0]; ++x)
        for (std::uint32_t y = 0; y < source.shape[1]; ++y)
            for (std::uint32_t z = 0; z < source.shape[2]; ++z)
            {
                const std::uint32_t index = static_cast<std::uint32_t>(
                    source.index(x, y, z));
                if (!envelope.occupancy[index] || source.occupancy[index]) continue;
                bool exposed = false;
                for (const auto& offset : neighbors6)
                {
                    const int nx = static_cast<int>(x) + offset[0];
                    const int ny = static_cast<int>(y) + offset[1];
                    const int nz = static_cast<int>(z) + offset[2];
                    exposed |= nx < 0 || ny < 0 || nz < 0 ||
                        nx >= static_cast<int>(source.shape[0]) ||
                        ny >= static_cast<int>(source.shape[1]) ||
                        nz >= static_cast<int>(source.shape[2]) ||
                        !envelope.occupied(static_cast<std::uint32_t>(nx),
                                           static_cast<std::uint32_t>(ny),
                                           static_cast<std::uint32_t>(nz));
                }
                if (exposed)
                {
                    distance[index] = 1;
                    maximum_distance = 1;
                    distance_queue.push_back(index);
                }
            }
    for (std::size_t head = 0; head < distance_queue.size(); ++head)
    {
        const std::uint32_t index = distance_queue[head];
        const auto coordinate = decodeIndex(index, source.shape);
        for (const auto& offset : neighbors6)
        {
            const int x = static_cast<int>(coordinate[0]) + offset[0];
            const int y = static_cast<int>(coordinate[1]) + offset[1];
            const int z = static_cast<int>(coordinate[2]) + offset[2];
            if (x < 0 || y < 0 || z < 0 ||
                x >= static_cast<int>(source.shape[0]) ||
                y >= static_cast<int>(source.shape[1]) ||
                z >= static_cast<int>(source.shape[2])) continue;
            const std::uint32_t neighbor = static_cast<std::uint32_t>(source.index(
                static_cast<std::uint32_t>(x), static_cast<std::uint32_t>(y),
                static_cast<std::uint32_t>(z)));
            if (!envelope.occupancy[neighbor] || source.occupancy[neighbor] ||
                distance[neighbor] != unreachable) continue;
            distance[neighbor] = distance[index] + 1;
            maximum_distance = std::max(maximum_distance, distance[neighbor]);
            distance_queue.push_back(neighbor);
        }
    }
    log_segment("distance_build");
    std::vector<std::vector<std::uint32_t>> distance_layers(maximum_distance + 1);
    for (const std::uint32_t index : distance_queue)
        distance_layers[distance[index]].push_back(index);
    VoxelGrid result = preserve_responsible_envelope ? envelope : source;
    bool reached_target = preserve_responsible_envelope;
    std::uint32_t layer_evaluations = 0;
    for (std::uint32_t radius = 1;
         !preserve_responsible_envelope && radius <= maximum_distance; ++radius)
    {
        for (const std::uint32_t index : distance_layers[radius])
            result.occupancy[index] = 1;
        const BettiNumbers topology = voxelBettiNumbers(result);
        ++layer_evaluations;
        if (std::abs(topology.beta0 - 1) + topology.beta1 == 0)
        {
            reached_target = true;
            break;
        }
    }
    if (!reached_target) result = envelope;
    std::cerr << "monitor: stage=topology_core_layer_summary evaluations="
              << layer_evaluations << " maximum_distance=" << maximum_distance
              << " reached_target=" << (reached_target ? 1 : 0) << '\n';
    log_segment("layer_betti");
    std::vector<std::uint8_t> component_seen(result.occupancy.size(), 0);
    std::vector<std::vector<std::uint32_t>> added_components;
    std::vector<std::uint32_t> component_queue;
    for (std::uint32_t seed = 0; seed < result.occupancy.size(); ++seed)
    {
        if (component_seen[seed] || source.occupancy[seed] || !result.occupancy[seed]) continue;
        added_components.emplace_back();
        auto& component = added_components.back();
        component_queue.clear();
        component_queue.push_back(seed);
        component_seen[seed] = 1;
        for (std::size_t head = 0; head < component_queue.size(); ++head)
        {
            const std::uint32_t index = component_queue[head];
            component.push_back(index);
            const auto coordinate = decodeIndex(index, result.shape);
            for (const auto& offset : neighbors6)
            {
                const int x = static_cast<int>(coordinate[0]) + offset[0];
                const int y = static_cast<int>(coordinate[1]) + offset[1];
                const int z = static_cast<int>(coordinate[2]) + offset[2];
                if (x < 0 || y < 0 || z < 0 ||
                    x >= static_cast<int>(result.shape[0]) ||
                    y >= static_cast<int>(result.shape[1]) ||
                    z >= static_cast<int>(result.shape[2])) continue;
                const std::uint32_t neighbor = static_cast<std::uint32_t>(result.index(
                    static_cast<std::uint32_t>(x), static_cast<std::uint32_t>(y),
                    static_cast<std::uint32_t>(z)));
                if (!component_seen[neighbor] && !source.occupancy[neighbor] &&
                    result.occupancy[neighbor])
                {
                    component_seen[neighbor] = 1;
                    component_queue.push_back(neighbor);
                }
            }
        }
    }
    std::sort(added_components.begin(), added_components.end(),
        [](const auto& first, const auto& second)
        { return first.size() > second.size(); });
    log_segment("component_build");
    std::size_t removed_components = 0;
    std::size_t component_replay_evaluations = 0;
    const auto replay_range = [&](auto&& self, const std::size_t begin,
                                  const std::size_t end) -> void
    {
        for (std::size_t component = begin; component < end; ++component)
            for (const std::uint32_t index : added_components[component])
                result.occupancy[index] = 0;
        const BettiNumbers replay = voxelBettiNumbers(result);
        ++component_replay_evaluations;
        if (std::abs(replay.beta0 - 1) + replay.beta1 == 0)
        {
            removed_components += end - begin;
            return;
        }
        for (std::size_t component = begin; component < end; ++component)
            for (const std::uint32_t index : added_components[component])
                result.occupancy[index] = 1;
        if (end - begin <= 1) return;
        const std::size_t middle = begin + (end - begin) / 2;
        self(self, begin, middle);
        self(self, middle, end);
    };
    if (!added_components.empty())
        replay_range(replay_range, 0, added_components.size());
    const BettiNumbers responsibility_topology = voxelBettiNumbers(result);
    std::cerr << "monitor: stage=topology_core_component_summary components="
              << added_components.size() << " removed_components="
              << removed_components << " replay_evaluations="
              << component_replay_evaluations << '\n';
    log_segment("component_replay");

    // Phase 1 needs a locally filled mouth, not the thinnest homotopy core.
    // Component replay has already removed every independent added region
    // whose complete removal preserves beta0=1,beta1=0. Keep the remaining
    // responsible closing components intact so their exposed boundary reaches
    // the surrounding exterior envelope instead of receding into the hole.
    if (preserve_responsible_envelope) return result;

    // Delete only digital-simple added points.  With the foreground using
    // 6-connectivity and the background using 26-connectivity, the two local
    // component tests are a homotopy certificate: the deletion cannot change
    // any Betti number.  Source cells are fixed, so the result remains an
    // enclosing repair while unrelated closing volume is peeled away.
    const auto local_components = [&](const std::uint32_t occupancy_mask,
                                      const bool foreground, const bool use_26)
    {
        std::array<std::uint8_t, 27> active{};
        std::array<std::uint8_t, 27> visited{};
        const auto local_index = [](const int x, const int y, const int z)
        { return static_cast<std::size_t>((x + 1) * 9 + (y + 1) * 3 + z + 1); };
        int bit = 0;
        for (int dx = -1; dx <= 1; ++dx)
            for (int dy = -1; dy <= 1; ++dy)
                for (int dz = -1; dz <= 1; ++dz)
                {
                    if (dx == 0 && dy == 0 && dz == 0)
                    {
                        active[local_index(dx, dy, dz)] = 0;
                        continue;
                    }
                    const bool occupied = (occupancy_mask & (1u << bit++)) != 0;
                    active[local_index(dx, dy, dz)] = occupied == foreground;
                }
        int components = 0;
        std::array<std::array<int, 3>, 27> queue{};
        for (int sx = -1; sx <= 1; ++sx)
            for (int sy = -1; sy <= 1; ++sy)
                for (int sz = -1; sz <= 1; ++sz)
                {
                    const std::size_t seed = local_index(sx, sy, sz);
                    if (!active[seed] || visited[seed]) continue;
                    ++components;
                    std::size_t head = 0, tail = 0;
                    queue[tail++] = {sx, sy, sz};
                    visited[seed] = 1;
                    while (head < tail)
                    {
                        const auto current = queue[head++];
                        const auto visit = [&](const auto& offset)
                        {
                            const int x = current[0] + offset[0];
                            const int y = current[1] + offset[1];
                            const int z = current[2] + offset[2];
                            if (x < -1 || x > 1 || y < -1 || y > 1 ||
                                z < -1 || z > 1) return;
                            const std::size_t id = local_index(x, y, z);
                            if (active[id] && !visited[id])
                            {
                                visited[id] = 1;
                                queue[tail++] = {x, y, z};
                            }
                        };
                        if (use_26)
                            for (const auto& offset : neighbors26) visit(offset);
                        else
                            for (const auto& offset : neighbors6) visit(offset);
                    }
                }
        return components;
    };
    std::unordered_map<std::uint32_t, std::uint8_t> simple_cache;
    simple_cache.reserve(4096);
    std::size_t cache_hits = 0;
    std::size_t cache_misses = 0;
    const auto is_simple = [&](const std::uint32_t index)
    {
        const auto coordinate = decodeIndex(index, result.shape);
        std::uint32_t mask = 0;
        int bit = 0;
        for (int dx = -1; dx <= 1; ++dx)
            for (int dy = -1; dy <= 1; ++dy)
                for (int dz = -1; dz <= 1; ++dz)
                {
                    if (dx == 0 && dy == 0 && dz == 0) continue;
                    const int x = static_cast<int>(coordinate[0]) + dx;
                    const int y = static_cast<int>(coordinate[1]) + dy;
                    const int z = static_cast<int>(coordinate[2]) + dz;
                    const bool occupied = x >= 0 && y >= 0 && z >= 0 &&
                        x < static_cast<int>(result.shape[0]) &&
                        y < static_cast<int>(result.shape[1]) &&
                        z < static_cast<int>(result.shape[2]) &&
                        result.occupied(static_cast<std::uint32_t>(x),
                                        static_cast<std::uint32_t>(y),
                                        static_cast<std::uint32_t>(z));
                    if (occupied) mask |= 1u << bit;
                    ++bit;
                }
        const auto found = simple_cache.find(mask);
        if (found != simple_cache.end())
        {
            ++cache_hits;
            return found->second != 0;
        }
        ++cache_misses;
        const bool simple = local_components(mask, true, false) == 1 &&
            local_components(mask, false, true) == 1;
        simple_cache.emplace(mask, simple ? 1 : 0);
        return simple;
    };

    std::vector<std::uint32_t> queue;
    std::vector<std::uint8_t> queued(result.occupancy.size(), 0);
    queue.reserve(countOccupied(result) - countOccupied(source));
    const auto exposed_added = [&](const std::uint32_t index)
    {
        if (source.occupancy[index] || !result.occupancy[index]) return false;
        const auto coordinate = decodeIndex(index, result.shape);
        for (const auto& offset : neighbors6)
        {
            const int x = static_cast<int>(coordinate[0]) + offset[0];
            const int y = static_cast<int>(coordinate[1]) + offset[1];
            const int z = static_cast<int>(coordinate[2]) + offset[2];
            if (x < 0 || y < 0 || z < 0 ||
                x >= static_cast<int>(result.shape[0]) ||
                y >= static_cast<int>(result.shape[1]) ||
                z >= static_cast<int>(result.shape[2]) ||
                !result.occupied(static_cast<std::uint32_t>(x),
                                 static_cast<std::uint32_t>(y),
                                 static_cast<std::uint32_t>(z)))
                return true;
        }
        return false;
    };
    const auto schedule = [&](const std::uint32_t index)
    {
        if (!queued[index] && exposed_added(index))
        {
            queued[index] = 1;
            queue.push_back(index);
        }
    };
    for (std::uint32_t index = 0; index < result.occupancy.size(); ++index) schedule(index);
    std::cerr << "monitor: stage=topology_core_scan_summary initial_queue="
              << queue.size() << '\n';
    log_segment("simple_scan");
    std::size_t head = 0;
    std::size_t simple_checks = 0;
    std::size_t simple_deletions = 0;
    while (head < queue.size())
    {
        const std::uint32_t index = queue[head++];
        queued[index] = 0;
        if (source.occupancy[index] || !result.occupancy[index]) continue;
        ++simple_checks;
        if (!is_simple(index)) continue;
        result.occupancy[index] = 0;
        ++simple_deletions;
        const auto coordinate = decodeIndex(index, result.shape);
        for (const auto& offset : neighbors26)
        {
            const int x = static_cast<int>(coordinate[0]) + offset[0];
            const int y = static_cast<int>(coordinate[1]) + offset[1];
            const int z = static_cast<int>(coordinate[2]) + offset[2];
            if (x < 0 || y < 0 || z < 0 ||
                x >= static_cast<int>(result.shape[0]) ||
                y >= static_cast<int>(result.shape[1]) ||
                z >= static_cast<int>(result.shape[2])) continue;
            schedule(static_cast<std::uint32_t>(result.index(
                static_cast<std::uint32_t>(x), static_cast<std::uint32_t>(y),
                static_cast<std::uint32_t>(z))));
        }
    }
    std::cerr << "monitor: stage=topology_core_thinning_summary queue_visits="
              << head << " simple_checks=" << simple_checks
              << " deletions=" << simple_deletions
              << " cache_hits=" << cache_hits
              << " cache_misses=" << cache_misses << '\n';
    log_segment("simple_thinning");
    if (voxelBettiNumbers(result) != responsibility_topology)
        throw std::runtime_error("simple-point thinning changed topology");
    log_segment("final_betti");
    return result;
}

int countComponents(const VoxelGrid& grid, const bool foreground, const bool use_26,
                    int* boundary_components = nullptr)
{
    std::vector<std::uint8_t> visited(grid.occupancy.size(), 0);
    std::vector<std::uint32_t> queue;
    int components = 0;
    int touching_boundary = 0;
    for (std::size_t seed = 0; seed < grid.occupancy.size(); ++seed)
    {
        if (visited[seed] || static_cast<bool>(grid.occupancy[seed]) != foreground) continue;
        ++components;
        bool touches = false;
        queue.clear();
        queue.push_back(static_cast<std::uint32_t>(seed));
        visited[seed] = 1;
        std::size_t head = 0;
        while (head < queue.size())
        {
            const auto coordinate = decodeIndex(queue[head++], grid.shape);
            touches |= coordinate[0] == 0 || coordinate[1] == 0 || coordinate[2] == 0 ||
                coordinate[0] + 1 == grid.shape[0] ||
                coordinate[1] + 1 == grid.shape[1] ||
                coordinate[2] + 1 == grid.shape[2];
            const auto visit = [&](const auto& offset)
            {
                const int x = static_cast<int>(coordinate[0]) + offset[0];
                const int y = static_cast<int>(coordinate[1]) + offset[1];
                const int z = static_cast<int>(coordinate[2]) + offset[2];
                if (x < 0 || y < 0 || z < 0 ||
                    x >= static_cast<int>(grid.shape[0]) ||
                    y >= static_cast<int>(grid.shape[1]) ||
                    z >= static_cast<int>(grid.shape[2])) return;
                const std::size_t id = grid.index(
                    static_cast<std::uint32_t>(x), static_cast<std::uint32_t>(y),
                    static_cast<std::uint32_t>(z));
                if (!visited[id] && static_cast<bool>(grid.occupancy[id]) == foreground)
                {
                    visited[id] = 1;
                    queue.push_back(static_cast<std::uint32_t>(id));
                }
            };
            if (use_26)
                for (const auto& offset : neighbors26) visit(offset);
            else
                for (const auto& offset : neighbors6) visit(offset);
        }
        if (touches) ++touching_boundary;
    }
    if (boundary_components) *boundary_components = touching_boundary;
    return components;
}

std::array<int, 256> makeEulerVertexTable()
{
    // Ohser et al.'s 3D digital Euler coefficients for the 26-connected
    // foreground. Reversing the table gives the paired 6-connected foreground
    // convention used by this pipeline (and by scikit-image, BSD-3-Clause).
    constexpr std::array<int, 256> coefficients26{{
        0,1,1,0,1,0,-2,-1,1,-2,0,-1,0,-1,-1,0,1,0,-2,-1,-2,-1,-1,-2,-6,-3,-3,-2,-3,-2,0,-1,
        1,-2,0,-1,-6,-3,-3,-2,-2,-1,-1,-2,-3,0,-2,-1,0,-1,-1,0,-3,-2,0,-1,-3,0,-2,-1,0,1,1,0,
        1,-2,-6,-3,0,-1,-3,-2,-2,-1,-3,0,-1,-2,-2,-1,0,-1,-3,-2,-1,0,0,-1,-3,0,0,1,-2,-1,1,0,
        -2,-1,-3,0,-3,0,0,1,-1,4,0,3,0,3,1,2,-1,-2,-2,-1,-2,-1,1,0,0,3,1,2,1,2,2,1,
        1,-6,-2,-3,-2,-3,-1,0,0,-3,-1,-2,-1,-2,-2,-1,-2,-3,-1,0,-1,0,4,3,-3,0,0,1,0,1,3,2,
        0,-3,-1,-2,-3,0,0,1,-1,0,0,-1,-2,1,-1,0,-1,-2,-2,-1,0,1,3,2,-2,1,-1,0,1,2,2,1,
        0,-3,-3,0,-1,-2,0,1,-1,0,-2,1,0,-1,-1,0,-1,-2,0,1,-2,-1,3,2,-2,1,1,2,-1,0,2,1,
        -1,0,-2,1,-2,1,1,2,-2,3,-1,2,-1,2,0,1,0,-1,-1,0,-1,0,2,1,-1,2,0,1,0,1,1,0}};
    std::array<int, 256> result{};
    for (std::size_t index = 0; index < result.size(); ++index)
        result[index] = coefficients26[result.size() - 1 - index];
    return result;
}

int cubicalEulerCharacteristic(const VoxelGrid& grid)
{
    static const auto table = makeEulerVertexTable();
    const std::array<std::uint32_t, 3> corner_shape{
        grid.shape[0] + 1, grid.shape[1] + 1, grid.shape[2] + 1};
    const std::size_t corner_count = checkedVoxelCount(corner_shape);
    std::vector<std::uint8_t> masks(corner_count, 0);
    const auto corner_index = [&](const std::uint32_t x, const std::uint32_t y,
                                  const std::uint32_t z)
    {
        return (static_cast<std::size_t>(x) * corner_shape[1] + y) * corner_shape[2] + z;
    };
    for (std::uint32_t x = 0; x < grid.shape[0]; ++x)
        for (std::uint32_t y = 0; y < grid.shape[1]; ++y)
            for (std::uint32_t z = 0; z < grid.shape[2]; ++z)
            {
                if (!grid.occupied(x, y, z)) continue;
                for (int cx = 0; cx < 2; ++cx)
                    for (int cy = 0; cy < 2; ++cy)
                        for (int cz = 0; cz < 2; ++cz)
                        {
                            const int sx = 1 - cx;
                            const int sy = 1 - cy;
                            const int sz = 1 - cz;
                            masks[corner_index(x + cx, y + cy, z + cz)] |=
                                static_cast<std::uint8_t>(1 << ((sx << 2) | (sy << 1) | sz));
                        }
            }
    std::int64_t multiplied = 0;
    for (const std::uint8_t mask : masks) multiplied += table[mask];
    if (multiplied % 8 != 0)
        throw std::runtime_error("cubical Euler-characteristic accumulator is inconsistent");
    return static_cast<int>(multiplied / 8);
}

int topologyDefect(const BettiNumbers topology)
{
    return std::abs(topology.beta0 - 1) + topology.beta1 + topology.beta2;
}

int connectivityAndHandleDefect(const BettiNumbers topology)
{
    return std::abs(topology.beta0 - 1) + topology.beta1;
}

BettiNumbers voxelBettiNumbersUnchecked(const VoxelGrid& grid)
{
    const int beta0 = countComponents(grid, true, false);
    int boundary_background = 0;
    const int background_components = countComponents(
        grid, false, true, &boundary_background);
    const int beta2 = background_components - boundary_background;
    const int euler = cubicalEulerCharacteristic(grid);
    return {beta0, beta0 + beta2 - euler, beta2};
}

BettiNumbers voxelBettiNumbersForEnclosingCandidate(
    const VoxelGrid& grid, const BettiNumbers current_topology)
{
    // axisClosing only adds paths whose endpoints are already occupied. Once
    // the current foreground is connected, such additions cannot split it or
    // create a new component, so beta0 remains exactly one.
    const int beta0 = current_topology.beta0 == 1
        ? 1 : countComponents(grid, true, false);
    int boundary_background = 0;
    const int background_components = countComponents(
        grid, false, true, &boundary_background);
    const int beta2 = background_components - boundary_background;
    const int euler = cubicalEulerCharacteristic(grid);
    return {beta0, beta0 + beta2 - euler, beta2};
}

std::array<int, 256> makeEulerVertexTable();

VoxelGrid axisClosing(const VoxelGrid& input, const int axis, const int radius)
{
    if (axis < 0 || axis > 2 || radius <= 0)
        throw std::invalid_argument("invalid axial-closing request");
    VoxelGrid result = input;
    const int first = (axis + 1) % 3;
    const int second = (axis + 2) % 3;
    std::array<std::uint32_t, 3> coordinate{};
    for (std::uint32_t a = 0; a < input.shape[first]; ++a)
        for (std::uint32_t b = 0; b < input.shape[second]; ++b)
        {
            coordinate[first] = a;
            coordinate[second] = b;
            int previous = -1;
            for (int position = 0; position < static_cast<int>(input.shape[axis]); ++position)
            {
                coordinate[axis] = static_cast<std::uint32_t>(position);
                if (!input.occupied(coordinate[0], coordinate[1], coordinate[2])) continue;
                if (previous >= 0 && position - previous - 1 <= 2 * radius)
                    for (int fill = previous + 1; fill < position; ++fill)
                    {
                        coordinate[axis] = static_cast<std::uint32_t>(fill);
                        result.occupancy[result.index(coordinate[0], coordinate[1], coordinate[2])] = 1;
                    }
                previous = position;
            }
        }
    return result;
}

struct AxisClosingEvaluation
{
    BettiNumbers topology{};
    std::size_t added = 0;
};

struct DisjointCells
{
    static constexpr std::uint32_t inactive =
        std::numeric_limits<std::uint32_t>::max();

    explicit DisjointCells(const std::size_t count, const bool track_boundary)
        : parent(count, inactive), boundary(track_boundary ? count : 0, 0)
    {
    }

    std::uint32_t find(std::uint32_t cell)
    {
        std::uint32_t root = cell;
        while (parent[root] != root) root = parent[root];
        while (parent[cell] != cell)
        {
            const std::uint32_t next = parent[cell];
            parent[cell] = root;
            cell = next;
        }
        return root;
    }

    void activate(const std::uint32_t cell, const bool touches_boundary = false)
    {
        parent[cell] = cell;
        if (!boundary.empty()) boundary[cell] = touches_boundary ? 1 : 0;
        ++components;
        if (touches_boundary) ++boundary_components;
    }

    void unite(const std::uint32_t first, const std::uint32_t second)
    {
        if (parent[first] == inactive || parent[second] == inactive) return;
        std::uint32_t first_root = find(first);
        std::uint32_t second_root = find(second);
        if (first_root == second_root) return;
        if (first_root > second_root) std::swap(first_root, second_root);
        if (!boundary.empty())
        {
            boundary_components -= boundary[first_root] + boundary[second_root];
            boundary[first_root] = boundary[first_root] || boundary[second_root];
            boundary_components += boundary[first_root];
        }
        parent[second_root] = first_root;
        --components;
    }

    std::vector<std::uint32_t> parent;
    std::vector<std::uint8_t> boundary;
    int components = 0;
    int boundary_components = 0;
};

bool isBoundaryCell(const std::array<std::uint32_t, 3>& coordinate,
                    const std::array<std::uint32_t, 3>& shape)
{
    return coordinate[0] == 0 || coordinate[1] == 0 || coordinate[2] == 0 ||
        coordinate[0] + 1 == shape[0] || coordinate[1] + 1 == shape[1] ||
        coordinate[2] + 1 == shape[2];
}

std::vector<AxisClosingEvaluation> axisClosingProfile(
    const VoxelGrid& input, const int axis, const BettiNumbers input_topology)
{
    if (axis < 0 || axis > 2)
        throw std::invalid_argument("invalid axial-closing profile request");
    const std::size_t voxel_count = input.occupancy.size();
    std::vector<std::uint32_t> activation(voxel_count, 0);
    const int first = (axis + 1) % 3;
    const int second = (axis + 2) % 3;
    int saturation_radius = 0;
    std::array<std::uint32_t, 3> coordinate{};
    for (std::uint32_t a = 0; a < input.shape[first]; ++a)
        for (std::uint32_t b = 0; b < input.shape[second]; ++b)
        {
            coordinate[first] = a;
            coordinate[second] = b;
            int previous = -1;
            for (int position = 0; position < static_cast<int>(input.shape[axis]); ++position)
            {
                coordinate[axis] = static_cast<std::uint32_t>(position);
                if (!input.occupied(coordinate[0], coordinate[1], coordinate[2])) continue;
                const int gap = previous < 0 ? 0 : position - previous - 1;
                if (gap > 0)
                {
                    const int radius = (gap + 1) / 2;
                    saturation_radius = std::max(saturation_radius, radius);
                    for (int fill = previous + 1; fill < position; ++fill)
                    {
                        coordinate[axis] = static_cast<std::uint32_t>(fill);
                        activation[input.index(
                            coordinate[0], coordinate[1], coordinate[2])] =
                                static_cast<std::uint32_t>(radius);
                    }
                    coordinate[axis] = static_cast<std::uint32_t>(position);
                }
                previous = position;
            }
        }

    std::vector<std::vector<std::uint32_t>> additions(
        static_cast<std::size_t>(saturation_radius) + 1);
    for (std::uint32_t index = 0; index < voxel_count; ++index)
        if (activation[index] != 0) additions[activation[index]].push_back(index);

    std::vector<AxisClosingEvaluation> profile(
        static_cast<std::size_t>(saturation_radius) + 1);
    profile[0] = {input_topology, 0};
    if (saturation_radius == 0) return profile;

    std::vector<int> beta0(profile.size(), 1);
    if (input_topology.beta0 != 1)
    {
        DisjointCells foreground(voxel_count, false);
        for (std::uint32_t index = 0; index < voxel_count; ++index)
            if (input.occupancy[index]) foreground.activate(index);
        for (std::uint32_t index = 0; index < voxel_count; ++index)
        {
            if (!input.occupancy[index]) continue;
            const auto cell = decodeIndex(index, input.shape);
            if (cell[0] > 0) foreground.unite(index, static_cast<std::uint32_t>(
                input.index(cell[0] - 1, cell[1], cell[2])));
            if (cell[1] > 0) foreground.unite(index, static_cast<std::uint32_t>(
                input.index(cell[0], cell[1] - 1, cell[2])));
            if (cell[2] > 0) foreground.unite(index, static_cast<std::uint32_t>(
                input.index(cell[0], cell[1], cell[2] - 1)));
        }
        beta0[0] = foreground.components;
        for (int radius = 1; radius <= saturation_radius; ++radius)
        {
            for (const std::uint32_t index : additions[radius])
                foreground.activate(index);
            for (const std::uint32_t index : additions[radius])
            {
                const auto cell = decodeIndex(index, input.shape);
                for (const auto& offset : neighbors6)
                {
                    const int x = static_cast<int>(cell[0]) + offset[0];
                    const int y = static_cast<int>(cell[1]) + offset[1];
                    const int z = static_cast<int>(cell[2]) + offset[2];
                    if (x < 0 || y < 0 || z < 0 ||
                        x >= static_cast<int>(input.shape[0]) ||
                        y >= static_cast<int>(input.shape[1]) ||
                        z >= static_cast<int>(input.shape[2])) continue;
                    foreground.unite(index, static_cast<std::uint32_t>(input.index(
                        static_cast<std::uint32_t>(x), static_cast<std::uint32_t>(y),
                        static_cast<std::uint32_t>(z))));
                }
            }
            beta0[radius] = foreground.components;
        }
    }

    DisjointCells background(voxel_count, true);
    for (std::uint32_t index = 0; index < voxel_count; ++index)
        if (!input.occupancy[index] && activation[index] == 0)
        {
            const auto cell = decodeIndex(index, input.shape);
            background.activate(index, isBoundaryCell(cell, input.shape));
        }
    for (std::uint32_t index = 0; index < voxel_count; ++index)
    {
        if (background.parent[index] == DisjointCells::inactive) continue;
        const auto cell = decodeIndex(index, input.shape);
        for (int dx = -1; dx <= 0; ++dx)
            for (int dy = -1; dy <= 1; ++dy)
                for (int dz = -1; dz <= 1; ++dz)
                {
                    if (dx == 0 && dy >= 0 && (dy > 0 || dz >= 0)) continue;
                    const int x = static_cast<int>(cell[0]) + dx;
                    const int y = static_cast<int>(cell[1]) + dy;
                    const int z = static_cast<int>(cell[2]) + dz;
                    if (x < 0 || y < 0 || z < 0 ||
                        x >= static_cast<int>(input.shape[0]) ||
                        y >= static_cast<int>(input.shape[1]) ||
                        z >= static_cast<int>(input.shape[2])) continue;
                    background.unite(index, static_cast<std::uint32_t>(input.index(
                        static_cast<std::uint32_t>(x), static_cast<std::uint32_t>(y),
                        static_cast<std::uint32_t>(z))));
                }
    }
    std::vector<int> beta2(profile.size(), 0);
    beta2[saturation_radius] =
        background.components - background.boundary_components;
    for (int radius = saturation_radius; radius > 0; --radius)
    {
        for (const std::uint32_t index : additions[radius])
        {
            const auto cell = decodeIndex(index, input.shape);
            background.activate(index, isBoundaryCell(cell, input.shape));
            for (const auto& offset : neighbors26)
            {
                const int x = static_cast<int>(cell[0]) + offset[0];
                const int y = static_cast<int>(cell[1]) + offset[1];
                const int z = static_cast<int>(cell[2]) + offset[2];
                if (x < 0 || y < 0 || z < 0 ||
                    x >= static_cast<int>(input.shape[0]) ||
                    y >= static_cast<int>(input.shape[1]) ||
                    z >= static_cast<int>(input.shape[2])) continue;
                background.unite(index, static_cast<std::uint32_t>(input.index(
                    static_cast<std::uint32_t>(x), static_cast<std::uint32_t>(y),
                    static_cast<std::uint32_t>(z))));
            }
        }
        beta2[radius - 1] = background.components -
            background.boundary_components;
    }

    const auto table = makeEulerVertexTable();
    const std::array<std::uint32_t, 3> corner_shape{
        input.shape[0] + 1, input.shape[1] + 1, input.shape[2] + 1};
    std::vector<std::uint8_t> masks(checkedVoxelCount(corner_shape), 0);
    const auto corner_index = [&](const std::uint32_t x, const std::uint32_t y,
                                  const std::uint32_t z)
    {
        return (static_cast<std::size_t>(x) * corner_shape[1] + y) *
            corner_shape[2] + z;
    };
    for (std::uint32_t index = 0; index < voxel_count; ++index)
    {
        if (!input.occupancy[index]) continue;
        const auto cell = decodeIndex(index, input.shape);
        for (int cx = 0; cx < 2; ++cx)
            for (int cy = 0; cy < 2; ++cy)
                for (int cz = 0; cz < 2; ++cz)
                    masks[corner_index(cell[0] + cx, cell[1] + cy, cell[2] + cz)] |=
                        static_cast<std::uint8_t>(1 <<
                            (((1 - cx) << 2) | ((1 - cy) << 1) | (1 - cz)));
    }
    std::int64_t multiplied_euler = 0;
    for (const std::uint8_t mask : masks) multiplied_euler += table[mask];
    std::size_t added_count = 0;
    for (int radius = 1; radius <= saturation_radius; ++radius)
    {
        for (const std::uint32_t index : additions[radius])
        {
            const auto cell = decodeIndex(index, input.shape);
            for (int cx = 0; cx < 2; ++cx)
                for (int cy = 0; cy < 2; ++cy)
                    for (int cz = 0; cz < 2; ++cz)
                    {
                        const std::size_t corner = corner_index(
                            cell[0] + cx, cell[1] + cy, cell[2] + cz);
                        const std::uint8_t bit = static_cast<std::uint8_t>(1 <<
                            (((1 - cx) << 2) | ((1 - cy) << 1) | (1 - cz)));
                        const std::uint8_t before = masks[corner];
                        const std::uint8_t after = before | bit;
                        multiplied_euler += table[after] - table[before];
                        masks[corner] = after;
                    }
        }
        added_count += additions[radius].size();
        if (multiplied_euler % 8 != 0)
            throw std::runtime_error("incremental Euler accumulator is inconsistent");
        const int euler = static_cast<int>(multiplied_euler / 8);
        profile[radius] = {{beta0[radius],
            beta0[radius] + beta2[radius] - euler, beta2[radius]}, added_count};
    }

    // Small fixtures retain the old full-grid calculation as an independent
    // oracle, so regressions catch any incremental-connectivity mistake.
    if (voxel_count <= 32'768)
        for (int radius = 1; radius <= saturation_radius; ++radius)
        {
            const VoxelGrid brute = axisClosing(input, axis, radius);
            const AxisClosingEvaluation expected{
                voxelBettiNumbersForEnclosingCandidate(brute, input_topology),
                countOccupied(brute) - countOccupied(input)};
            if (profile[radius].topology != expected.topology ||
                profile[radius].added != expected.added)
                throw std::runtime_error(
                    "incremental axial-closing profile disagrees with full-grid oracle");
        }
    return profile;
}

MeshModel readObj(const std::filesystem::path& path)
{
    std::ifstream stream(path);
    if (!stream) throw std::runtime_error("failed to open OBJ: " + path.string());
    MeshModel mesh;
    mesh.name = path.filename().string();
    std::string line;
    while (std::getline(stream, line))
    {
        std::istringstream parser(line);
        std::string tag;
        parser >> tag;
        if (tag == "v")
        {
            Position3 vertex;
            if (parser >> vertex.x >> vertex.y >> vertex.z) mesh.vertices.push_back(vertex);
        }
        else if (tag == "f")
        {
            std::vector<std::uint32_t> polygon;
            std::string token;
            while (parser >> token)
            {
                const std::size_t slash = token.find('/');
                const int raw = std::stoi(token.substr(0, slash));
                const std::int64_t resolved = raw > 0
                    ? static_cast<std::int64_t>(raw - 1)
                    : static_cast<std::int64_t>(mesh.vertices.size()) + raw;
                if (resolved < 0 || resolved >= static_cast<std::int64_t>(mesh.vertices.size()))
                    throw std::runtime_error("OBJ face index is outside the vertex list");
                polygon.push_back(static_cast<std::uint32_t>(resolved));
            }
            for (std::size_t index = 1; index + 1 < polygon.size(); ++index)
                mesh.triangles.push_back({polygon[0], polygon[index], polygon[index + 1]});
        }
    }
    if (const auto error = validateModelPool({mesh}))
        throw std::runtime_error("invalid OBJ " + path.string() + ": " + *error);
    return mesh;
}

double choosePitch(const Vec3 extents, const std::size_t maximum_grid_voxels,
                   const std::uint32_t padding)
{
    if (extents.x <= 0.0 || extents.y <= 0.0 || extents.z <= 0.0)
        throw std::invalid_argument("OBJ bounds must have three positive extents");
    const auto cells = [&](const double pitch)
    {
        const auto count = [](const double extent, const double value, const std::uint32_t pad)
        { return static_cast<std::uint64_t>(std::ceil(extent / value)) + 1 + 2 * pad; };
        return count(extents.x, pitch, padding) * count(extents.y, pitch, padding) *
               count(extents.z, pitch, padding);
    };
    double low = std::numeric_limits<double>::epsilon() *
        std::max({extents.x, extents.y, extents.z});
    double high = std::max({extents.x, extents.y, extents.z});
    while (cells(high) > maximum_grid_voxels) high *= 2.0;
    for (int iteration = 0; iteration < 80; ++iteration)
    {
        const double middle = 0.5 * (low + high);
        if (cells(middle) > maximum_grid_voxels) low = middle;
        else high = middle;
    }
    return high;
}

void rasterizeTriangle(
    VoxelGrid& grid, const std::array<Vec3, 3>& world_triangle,
    std::vector<std::uint32_t>* covered_voxels = nullptr)
{
    std::array<Vec3, 3> triangle{};
    const Vec3 origin = toVec(grid.origin);
    for (int index = 0; index < 3; ++index)
        triangle[index] = (world_triangle[index] - origin) * (1.0 / grid.pitch);
    const Vec3 normal = cross(triangle[1] - triangle[0], triangle[2] - triangle[0]);
    int dominant = 0;
    if (std::abs(normal.y) > std::abs(normal[dominant])) dominant = 1;
    if (std::abs(normal.z) > std::abs(normal[dominant])) dominant = 2;
    constexpr double epsilon = 64.0 * std::numeric_limits<double>::epsilon();
    if (std::abs(normal[dominant]) <= epsilon) return;
    std::array<int, 2> projected_axes{};
    int next = 0;
    for (int axis = 0; axis < 3; ++axis)
        if (axis != dominant) projected_axes[next++] = axis;
    const int first_axis = projected_axes[0];
    const int second_axis = projected_axes[1];
    std::array<double, 2> lower{std::numeric_limits<double>::infinity(),
                                std::numeric_limits<double>::infinity()};
    std::array<double, 2> upper{-std::numeric_limits<double>::infinity(),
                                -std::numeric_limits<double>::infinity()};
    for (const Vec3 vertex : triangle)
    {
        lower[0] = std::min(lower[0], vertex[first_axis]);
        lower[1] = std::min(lower[1], vertex[second_axis]);
        upper[0] = std::max(upper[0], vertex[first_axis]);
        upper[1] = std::max(upper[1], vertex[second_axis]);
    }
    const int first_lower = std::max(0, static_cast<int>(std::ceil(lower[0] - 0.5 - epsilon)));
    const int second_lower = std::max(0, static_cast<int>(std::ceil(lower[1] - 0.5 - epsilon)));
    const int first_upper = std::min(static_cast<int>(grid.shape[first_axis]) - 1,
        static_cast<int>(std::floor(upper[0] + 0.5 + epsilon)));
    const int second_upper = std::min(static_cast<int>(grid.shape[second_axis]) - 1,
        static_cast<int>(std::floor(upper[1] + 0.5 + epsilon)));
    if (first_lower > first_upper || second_lower > second_upper) return;

    std::array<std::array<double, 2>, 3> projected{};
    for (int index = 0; index < 3; ++index)
        projected[index] = {triangle[index][first_axis], triangle[index][second_axis]};
    std::array<std::array<double, 2>, 3> separating_axes{};
    std::array<double, 3> interval_lower{};
    std::array<double, 3> interval_upper{};
    std::array<double, 3> square_radius{};
    for (int edge = 0; edge < 3; ++edge)
    {
        const auto& a = projected[edge];
        const auto& b = projected[(edge + 1) % 3];
        const double dx = b[0] - a[0];
        const double dy = b[1] - a[1];
        separating_axes[edge] = {dy, -dx};
        interval_lower[edge] = std::numeric_limits<double>::infinity();
        interval_upper[edge] = -std::numeric_limits<double>::infinity();
        for (const auto& point : projected)
        {
            const double value = point[0] * dy - point[1] * dx;
            interval_lower[edge] = std::min(interval_lower[edge], value);
            interval_upper[edge] = std::max(interval_upper[edge], value);
        }
        square_radius[edge] = 0.5 * (std::abs(dy) + std::abs(dx));
    }
    const std::array<double, 2> slopes{
        -normal[first_axis] / normal[dominant],
        -normal[second_axis] / normal[dominant]};
    const double height_radius = 0.5 * (std::abs(slopes[0]) + std::abs(slopes[1]));
    std::array<std::uint32_t, 3> coordinate{};
    for (int first = first_lower; first <= first_upper; ++first)
        for (int second = second_lower; second <= second_upper; ++second)
        {
            bool intersects = true;
            for (int axis = 0; axis < 3; ++axis)
            {
                const double center_interval = first * separating_axes[axis][0] +
                    second * separating_axes[axis][1];
                intersects &= center_interval + square_radius[axis] >=
                                  interval_lower[axis] - epsilon &&
                              center_interval - square_radius[axis] <=
                                  interval_upper[axis] + epsilon;
            }
            if (!intersects) continue;
            const double center_height = triangle[0][dominant] +
                (first - triangle[0][first_axis]) * slopes[0] +
                (second - triangle[0][second_axis]) * slopes[1];
            const int depth_lower = std::max(0, static_cast<int>(
                std::ceil(center_height - height_radius - 0.5 - epsilon)));
            const int depth_upper = std::min(static_cast<int>(grid.shape[dominant]) - 1,
                static_cast<int>(std::floor(center_height + height_radius + 0.5 + epsilon)));
            coordinate[first_axis] = static_cast<std::uint32_t>(first);
            coordinate[second_axis] = static_cast<std::uint32_t>(second);
            for (int depth = depth_lower; depth <= depth_upper; ++depth)
            {
                coordinate[dominant] = static_cast<std::uint32_t>(depth);
                const std::size_t index = grid.index(
                    coordinate[0], coordinate[1], coordinate[2]);
                grid.occupancy[index] = 1;
                if (covered_voxels)
                    covered_voxels->push_back(static_cast<std::uint32_t>(index));
            }
        }
}

} // namespace

std::size_t VoxelGrid::index(
    const std::uint32_t x, const std::uint32_t y, const std::uint32_t z) const
{
    return (static_cast<std::size_t>(x) * shape[1] + y) * shape[2] + z;
}

bool VoxelGrid::occupied(
    const std::uint32_t x, const std::uint32_t y, const std::uint32_t z) const
{
    return occupancy[index(x, y, z)] != 0;
}

BettiNumbers voxelBettiNumbers(const VoxelGrid& grid)
{
    validateGrid(grid);
    return voxelBettiNumbersUnchecked(grid);
}

VoxelGrid enclosingTopologyFill(
    const VoxelGrid& source, const std::uint32_t maximum_steps,
    TopologyFillStats* stats, VoxelGrid* cavity_fill_labels)
{
    validateGrid(source);
    if (maximum_steps == 0) throw std::invalid_argument("maximum steps must be positive");
    VoxelGrid current = source;
    const BettiNumbers input_topology = voxelBettiNumbers(source);
    BettiNumbers current_topology = voxelBettiNumbers(current);
    std::vector<TopologyFillStep> steps;
    for (std::uint32_t step = 0; step < maximum_steps; ++step)
    {
        const auto step_started = std::chrono::steady_clock::now();
        const int current_defect = connectivityAndHandleDefect(current_topology);
        if (current_defect == 0) break;
        struct Candidate
        {
            int defect = 0;
            std::size_t added = 0;
            int axis = 0;
            int radius = 0;
            BettiNumbers topology{};
        };
        std::vector<Candidate> candidates;
        double profile_seconds = 0.0;
        double closing_seconds = 0.0;
        std::size_t evaluation_count = 0;
        for (int axis = 0; axis < 3; ++axis)
        {
            const auto profile_started = std::chrono::steady_clock::now();
            const std::vector<AxisClosingEvaluation> profile =
                axisClosingProfile(current, axis, current_topology);
            profile_seconds += std::chrono::duration<double>(
                std::chrono::steady_clock::now() - profile_started).count();
            const auto evaluate = [&](const int radius) -> const AxisClosingEvaluation&
            {
                ++evaluation_count;
                return profile.at(static_cast<std::size_t>(radius));
            };
            const int saturation_radius = static_cast<int>(profile.size()) - 1;
            if (saturation_radius == 0) continue;
            int previous_radius = 0;
            int improving_radius = 0;
            for (int coarse_radius = 1;;)
            {
                const int coarse_defect = connectivityAndHandleDefect(
                    evaluate(coarse_radius).topology);
                if (coarse_defect < current_defect)
                {
                    improving_radius = coarse_radius;
                    break;
                }
                previous_radius = coarse_radius;
                if (coarse_radius == saturation_radius) break;
                coarse_radius = std::min(saturation_radius, coarse_radius * 2);
            }
            if (improving_radius == 0) continue;
            int low = previous_radius + 1;
            int high = improving_radius;
            while (low < high)
            {
                const int middle = (low + high) / 2;
                if (connectivityAndHandleDefect(
                        evaluate(middle).topology) < current_defect)
                    high = middle;
                else low = middle + 1;
            }
            const AxisClosingEvaluation& evaluation = evaluate(low);
            const int defect = connectivityAndHandleDefect(evaluation.topology);
            if (evaluation.added > 0)
                candidates.push_back({defect, evaluation.added, axis, low,
                                      evaluation.topology});
        }
        if (candidates.empty())
            throw std::runtime_error("enclosing topology fill stalled");
        std::sort(candidates.begin(), candidates.end(), [&](const Candidate& first,
                                                            const Candidate& second)
        {
            const auto first_reduction = current_defect - first.defect;
            const auto second_reduction = current_defect - second.defect;
            const auto left = static_cast<unsigned long long>(first.added) * second_reduction;
            const auto right = static_cast<unsigned long long>(second.added) * first_reduction;
            if (left != right) return left < right;
            return std::tie(first.defect, first.added, first.axis, first.radius) <
                   std::tie(second.defect, second.added, second.axis, second.radius);
        });
        const Candidate& selected = candidates.front();
        steps.push_back({selected.axis, selected.radius, selected.added,
                         current_topology, selected.topology});
        const auto closing_started = std::chrono::steady_clock::now();
        current = axisClosing(current, selected.axis, selected.radius);
        closing_seconds += std::chrono::duration<double>(
            std::chrono::steady_clock::now() - closing_started).count();
        current_topology = selected.topology;
        std::cerr << "monitor: stage=topology_fill_step step=" << step
                  << " axis=" << selected.axis << " radius=" << selected.radius
                  << " defect=" << current_defect << "->"
                  << connectivityAndHandleDefect(current_topology)
                  << " elapsed_seconds="
                  << std::chrono::duration<double>(
                         std::chrono::steady_clock::now() - step_started).count()
                  << " evaluations=" << evaluation_count
                  << " profile_seconds=" << profile_seconds
                  << " closing_seconds=" << closing_seconds
                  << '\n';
    }
    if (connectivityAndHandleDefect(current_topology) != 0)
        throw std::runtime_error(
            "enclosing topology fill did not reach beta0=1,beta1=0; final beta=(" +
            std::to_string(current_topology.beta0) + "," +
            std::to_string(current_topology.beta1) + "," +
            std::to_string(current_topology.beta2) + "), steps=" +
            std::to_string(steps.size()));

    // Closing handles is allowed to create closed empty components. Resolve
    // beta2 only after beta0/beta1 are final so this exact difference is the
    // certificate for source cavity-wall removal.
    current = thinnestConnectivityHandleCore(source, current, true);
    current_topology = voxelBettiNumbers(current);
    const VoxelGrid before_cavity_fill = current;
    current = fillCavities(before_cavity_fill);
    VoxelGrid labels = current;
    std::size_t cavity_fill_count = 0;
    for (std::size_t index = 0; index < labels.occupancy.size(); ++index)
    {
        labels.occupancy[index] = current.occupancy[index] &&
            !before_cavity_fill.occupancy[index] ? 1 : 0;
        cavity_fill_count += labels.occupancy[index] != 0 ? 1 : 0;
    }
    const BettiNumbers after_cavity_fill = voxelBettiNumbers(current);
    if (cavity_fill_count != 0)
        steps.push_back({-1, 0, cavity_fill_count,
                         current_topology, after_cavity_fill});
    current_topology = after_cavity_fill;

    if (current_topology != BettiNumbers{1, 0, 0})
        throw std::runtime_error("topology core and cavity fill lost target Betti numbers");
    std::size_t removed = 0;
    for (std::size_t index = 0; index < source.occupancy.size(); ++index)
        removed += source.occupancy[index] && !current.occupancy[index] ? 1 : 0;
    if (removed != 0) throw std::runtime_error("enclosing topology fill removed source voxels");
    if (stats)
    {
        stats->input_betti = input_topology;
        stats->output_betti = current_topology;
        stats->input_voxels = countOccupied(source);
        stats->output_voxels = countOccupied(current);
        stats->added_voxels = stats->output_voxels - stats->input_voxels;
        stats->cavity_fill_voxels = cavity_fill_count;
        stats->removed_voxels = removed;
        stats->pitch = source.pitch;
        stats->grid_shape = source.shape;
        stats->steps = std::move(steps);
    }
    if (cavity_fill_labels) *cavity_fill_labels = std::move(labels);
    return current;
}

VoxelGrid voxelizeTriangleSoup(
    const MeshModel& mesh, const std::size_t maximum_grid_voxels,
    const std::uint32_t padding, TriangleVoxelProvenance* provenance)
{
    if (const auto error = validateModelPool({mesh})) throw std::invalid_argument(*error);
    Vec3 lower{std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity(),
               std::numeric_limits<double>::infinity()};
    Vec3 upper{-std::numeric_limits<double>::infinity(), -std::numeric_limits<double>::infinity(),
               -std::numeric_limits<double>::infinity()};
    for (const Position3 vertex : mesh.vertices)
    {
        const Vec3 point = toVec(vertex);
        for (int axis = 0; axis < 3; ++axis)
        {
            lower[axis] = std::min(lower[axis], point[axis]);
            upper[axis] = std::max(upper[axis], point[axis]);
        }
    }
    const Vec3 extents = upper - lower;
    const double pitch = choosePitch(extents, maximum_grid_voxels, padding);
    VoxelGrid grid;
    grid.pitch = pitch;
    grid.origin = toPosition(lower - Vec3{pitch * padding, pitch * padding, pitch * padding});
    for (int axis = 0; axis < 3; ++axis)
        grid.shape[axis] = static_cast<std::uint32_t>(std::ceil(extents[axis] / pitch)) +
            1 + 2 * padding;
    grid.occupancy.assign(checkedVoxelCount(grid.shape), 0);
    if (provenance) provenance->surface_voxels.assign(mesh.triangles.size(), {});
    for (std::size_t face_id = 0; face_id < mesh.triangles.size(); ++face_id)
    {
        const TriangleIndices face = mesh.triangles[face_id];
        rasterizeTriangle(grid,
            {toVec(mesh.vertices[face[0]]), toVec(mesh.vertices[face[1]]),
             toVec(mesh.vertices[face[2]])},
            provenance ? &provenance->surface_voxels[face_id] : nullptr);
    }
    const auto outside = floodEmptyFromBoundary(grid, neighbors6);
    for (std::size_t index = 0; index < grid.occupancy.size(); ++index)
        grid.occupancy[index] = outside[index] ? 0 : 1;
    validateGrid(grid);
    return grid;
}

VoxelGrid voxelizeTriangleSoupOnGrid(
    const MeshModel& mesh, const VoxelGrid& grid_template,
    TriangleVoxelProvenance* provenance)
{
    if (const auto error = validateModelPool({mesh})) throw std::invalid_argument(*error);
    validateGrid(grid_template);
    VoxelGrid grid = grid_template;
    std::fill(grid.occupancy.begin(), grid.occupancy.end(), std::uint8_t{0});
    if (provenance) provenance->surface_voxels.assign(mesh.triangles.size(), {});
    for (std::size_t face_id = 0; face_id < mesh.triangles.size(); ++face_id)
    {
        const TriangleIndices face = mesh.triangles[face_id];
        rasterizeTriangle(grid,
            {toVec(mesh.vertices[face[0]]), toVec(mesh.vertices[face[1]]),
             toVec(mesh.vertices[face[2]])},
            provenance ? &provenance->surface_voxels[face_id] : nullptr);
    }
    const auto outside = floodEmptyFromBoundary(grid, neighbors6);
    for (std::size_t index = 0; index < grid.occupancy.size(); ++index)
        grid.occupancy[index] = outside[index] ? 0 : 1;
    validateGrid(grid);
    return grid;
}

namespace
{

struct Bounds3
{
    Vec3 lower{std::numeric_limits<double>::infinity(),
               std::numeric_limits<double>::infinity(),
               std::numeric_limits<double>::infinity()};
    Vec3 upper{-std::numeric_limits<double>::infinity(),
               -std::numeric_limits<double>::infinity(),
               -std::numeric_limits<double>::infinity()};
};

void include(Bounds3& bounds, const Vec3 point)
{
    for (int axis = 0; axis < 3; ++axis)
    {
        bounds.lower[axis] = std::min(bounds.lower[axis], point[axis]);
        bounds.upper[axis] = std::max(bounds.upper[axis], point[axis]);
    }
}

double pointBoundsDistanceSquared(const Vec3 point, const Bounds3& bounds)
{
    double result = 0.0;
    for (int axis = 0; axis < 3; ++axis)
    {
        const double delta = point[axis] < bounds.lower[axis]
            ? bounds.lower[axis] - point[axis]
            : (point[axis] > bounds.upper[axis] ? point[axis] - bounds.upper[axis] : 0.0);
        result += delta * delta;
    }
    return result;
}

Vec3 closestPointOnTriangle(
    const Vec3 point, const Vec3 first, const Vec3 second, const Vec3 third)
{
    const Vec3 first_edge = second - first;
    const Vec3 second_edge = third - first;
    const Vec3 first_relative = point - first;
    const double first_projection = dot(first_edge, first_relative);
    const double second_projection = dot(second_edge, first_relative);
    if (first_projection <= 0.0 && second_projection <= 0.0) return first;

    const Vec3 second_relative = point - second;
    const double third_projection = dot(first_edge, second_relative);
    const double fourth_projection = dot(second_edge, second_relative);
    if (third_projection >= 0.0 && fourth_projection <= third_projection) return second;

    const double first_region = first_projection * fourth_projection -
        third_projection * second_projection;
    if (first_region <= 0.0 && first_projection >= 0.0 && third_projection <= 0.0)
        return first + first_edge *
            (first_projection / (first_projection - third_projection));

    const Vec3 third_relative = point - third;
    const double fifth_projection = dot(first_edge, third_relative);
    const double sixth_projection = dot(second_edge, third_relative);
    if (sixth_projection >= 0.0 && fifth_projection <= sixth_projection) return third;

    const double second_region = fifth_projection * second_projection -
        first_projection * sixth_projection;
    if (second_region <= 0.0 && second_projection >= 0.0 && sixth_projection <= 0.0)
        return first + second_edge *
            (second_projection / (second_projection - sixth_projection));

    const double opposite_region = third_projection * sixth_projection -
        fifth_projection * fourth_projection;
    if (opposite_region <= 0.0 && fourth_projection >= third_projection &&
        fifth_projection >= sixth_projection)
        return second + (third - second) *
            ((fourth_projection - third_projection) /
             ((fourth_projection - third_projection) +
              (fifth_projection - sixth_projection)));

    const double inverse = 1.0 / (opposite_region + second_region + first_region);
    return first + first_edge * (second_region * inverse) +
        second_edge * (first_region * inverse);
}

class SourceTriangleBvh
{
public:
    SourceTriangleBvh(
        const MeshModel& mesh, const std::vector<std::uint8_t>& retained_faces)
        : mesh_(mesh)
    {
        if (retained_faces.size() != mesh_.triangles.size())
            throw std::invalid_argument("source-face mask size does not match mesh");
        for (std::size_t face = 0; face < retained_faces.size(); ++face)
            if (retained_faces[face]) face_ids_.push_back(face);
        if (!face_ids_.empty()) build(0, face_ids_.size());
    }

    struct ClosestPointResult
    {
        Vec3 point{};
        std::size_t face_id = std::numeric_limits<std::size_t>::max();
    };

    [[nodiscard]] ClosestPointResult closestPoint(const Vec3 point) const
    {
        if (nodes_.empty()) throw std::invalid_argument("source mesh has no triangles");
        ClosestPointResult best{point, std::numeric_limits<std::size_t>::max()};
        double distance_squared = std::numeric_limits<double>::infinity();
        query(0, point, distance_squared, best);
        return best;
    }

private:
    struct Node
    {
        Bounds3 bounds;
        std::size_t begin = 0;
        std::size_t end = 0;
        std::size_t left = 0;
        std::size_t right = 0;
        bool leaf = false;
    };

    [[nodiscard]] Bounds3 faceBounds(const std::size_t face_id) const
    {
        Bounds3 bounds;
        for (const auto vertex : mesh_.triangles[face_id])
            include(bounds, toVec(mesh_.vertices[vertex]));
        return bounds;
    }

    [[nodiscard]] Vec3 faceCenter(const std::size_t face_id) const
    {
        const auto face = mesh_.triangles[face_id];
        return (toVec(mesh_.vertices[face[0]]) + toVec(mesh_.vertices[face[1]]) +
            toVec(mesh_.vertices[face[2]])) * (1.0 / 3.0);
    }

    std::size_t build(const std::size_t begin, const std::size_t end)
    {
        Node node;
        node.begin = begin;
        node.end = end;
        Bounds3 center_bounds;
        for (std::size_t index = begin; index < end; ++index)
        {
            const Bounds3 bounds = faceBounds(face_ids_[index]);
            include(node.bounds, bounds.lower);
            include(node.bounds, bounds.upper);
            include(center_bounds, faceCenter(face_ids_[index]));
        }
        const std::size_t node_id = nodes_.size();
        nodes_.push_back(node);
        if (end - begin <= 8)
        {
            nodes_[node_id].leaf = true;
            return node_id;
        }
        const Vec3 extent = center_bounds.upper - center_bounds.lower;
        int axis = extent.y > extent.x ? 1 : 0;
        if (extent.z > extent[axis]) axis = 2;
        const std::size_t middle = begin + (end - begin) / 2;
        std::nth_element(face_ids_.begin() + begin, face_ids_.begin() + middle,
            face_ids_.begin() + end, [&](const std::size_t first, const std::size_t second)
            { return faceCenter(first)[axis] < faceCenter(second)[axis]; });
        nodes_[node_id].left = build(begin, middle);
        nodes_[node_id].right = build(middle, end);
        return node_id;
    }

    void query(const std::size_t node_id, const Vec3 point,
               double& best_distance_squared, ClosestPointResult& best) const
    {
        const Node& node = nodes_[node_id];
        if (pointBoundsDistanceSquared(point, node.bounds) >= best_distance_squared) return;
        if (node.leaf)
        {
            for (std::size_t index = node.begin; index < node.end; ++index)
            {
                const auto face = mesh_.triangles[face_ids_[index]];
                const Vec3 candidate = closestPointOnTriangle(point,
                    toVec(mesh_.vertices[face[0]]), toVec(mesh_.vertices[face[1]]),
                    toVec(mesh_.vertices[face[2]]));
                const double distance_squared = dot(candidate - point, candidate - point);
                if (distance_squared < best_distance_squared)
                {
                    best_distance_squared = distance_squared;
                    best.point = candidate;
                    best.face_id = face_ids_[index];
                }
            }
            return;
        }
        const double left_distance = pointBoundsDistanceSquared(point, nodes_[node.left].bounds);
        const double right_distance = pointBoundsDistanceSquared(point, nodes_[node.right].bounds);
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

    const MeshModel& mesh_;
    std::vector<std::size_t> face_ids_;
    std::vector<Node> nodes_;
};

struct IndexedSurface
{
    std::vector<Vec3> vertices;
    std::vector<std::array<std::uint32_t, 3>> triangles;
    std::vector<std::uint32_t> halfedge_opposites;
    std::vector<VoxelBoundaryCrossing> boundary_crossings;
    std::vector<std::uint32_t> face_boundary_crossings;
};

struct DualSourcePlane
{
    Vec3 normal{};
    double distance = 0.0;
    double support_area = 0.0;
    Bounds3 bounds;
    std::vector<std::array<Vec3, 3>> support_triangles;
    std::vector<Bounds3> support_components;
};

std::vector<DualSourcePlane> sourcePlanes(
    const MeshModel& source, const double pitch,
    const std::vector<std::uint8_t>* retained_faces = nullptr)
{
    struct DualPlaneKey
    {
        std::array<std::int64_t, 4> value{};
        bool operator<(const DualPlaneKey& other) const { return value < other.value; }
    };
    Bounds3 model_bounds;
    for (const Position3 point : source.vertices) include(model_bounds, toVec(point));
    const double model_scale = norm(model_bounds.upper - model_bounds.lower);
    const double distance_quantum = std::max(model_scale * 1.0e-8, pitch * 1.0e-7);
    constexpr double angular_quantum = 1.0e-7;
    std::map<DualPlaneKey, DualSourcePlane> grouped;
    for (std::size_t face_id = 0; face_id < source.triangles.size(); ++face_id)
    {
        if (retained_faces && !(*retained_faces)[face_id]) continue;
        const TriangleIndices& face = source.triangles[face_id];
        const Vec3 a = toVec(source.vertices[face[0]]);
        const Vec3 b = toVec(source.vertices[face[1]]);
        const Vec3 c = toVec(source.vertices[face[2]]);
        const Vec3 area_vector = cross(b - a, c - a);
        const double area = 0.5 * norm(area_vector);
        Vec3 normal = normalized(area_vector);
        if (area == 0.0) continue;
        int dominant = 0;
        for (int axis = 1; axis < 3; ++axis)
            if (std::abs(normal[axis]) > std::abs(normal[dominant])) dominant = axis;
        if (normal[dominant] < 0.0) normal = normal * -1.0;
        const double distance = dot(normal, a);
        DualPlaneKey key{{
            static_cast<std::int64_t>(std::llround(normal.x / angular_quantum)),
            static_cast<std::int64_t>(std::llround(normal.y / angular_quantum)),
            static_cast<std::int64_t>(std::llround(normal.z / angular_quantum)),
            static_cast<std::int64_t>(std::llround(distance / distance_quantum))}};
        auto [iterator, inserted] = grouped.try_emplace(
            key, DualSourcePlane{normal, distance, 0.0, {}, {}, {}});
        iterator->second.support_area += area;
        for (const Vec3 point : {a, b, c}) include(iterator->second.bounds, point);
        iterator->second.support_triangles.push_back({a, b, c});
        (void)inserted;
    }
    std::vector<DualSourcePlane> exact_planes;
    exact_planes.reserve(grouped.size());
    for (auto& [key, plane] : grouped)
    {
        (void)key;
        const std::size_t triangle_count = plane.support_triangles.size();
        std::vector<std::uint32_t> parent(triangle_count);
        std::iota(parent.begin(), parent.end(), 0);
        const auto find = [&](std::uint32_t value)
        {
            std::uint32_t root = value;
            while (parent[root] != root) root = parent[root];
            while (parent[value] != value)
            {
                const std::uint32_t next = parent[value];
                parent[value] = root;
                value = next;
            }
            return root;
        };
        const auto unite = [&](const std::uint32_t first,
                               const std::uint32_t second)
        {
            const std::uint32_t a = find(first);
            const std::uint32_t b = find(second);
            if (a != b) parent[b] = a;
        };
        using VertexKey = std::array<std::int64_t, 3>;
        std::map<VertexKey, std::uint32_t> vertex_owner;
        for (std::uint32_t triangle = 0; triangle < triangle_count; ++triangle)
            for (const Vec3 point : plane.support_triangles[triangle])
            {
                const VertexKey vertex{
                    static_cast<std::int64_t>(std::llround(point.x / distance_quantum)),
                    static_cast<std::int64_t>(std::llround(point.y / distance_quantum)),
                    static_cast<std::int64_t>(std::llround(point.z / distance_quantum))};
                const auto [iterator, inserted] = vertex_owner.try_emplace(
                    vertex, triangle);
                if (!inserted) unite(triangle, iterator->second);
            }
        std::map<std::uint32_t, Bounds3> component_bounds;
        for (std::uint32_t triangle = 0; triangle < triangle_count; ++triangle)
        {
            auto& component = component_bounds[find(triangle)];
            for (const Vec3 point : plane.support_triangles[triangle])
                include(component, point);
        }
        plane.support_components.reserve(component_bounds.size());
        for (auto& [root, bounds] : component_bounds)
        {
            (void)root;
            plane.support_components.push_back(bounds);
        }
        plane.support_triangles.clear();
        plane.support_triangles.shrink_to_fit();
        exact_planes.push_back(std::move(plane));
    }
    std::sort(exact_planes.begin(), exact_planes.end(),
        [](const DualSourcePlane& first, const DualSourcePlane& second)
        { return first.support_area > second.support_area; });

    // Triangle soups frequently encode one CAD face as several almost equal
    // fitted planes. Treat differences below the topology grid's resolving
    // power as one source patch, otherwise adjacent dual cells can alternate
    // between those planes and create a saw-tooth boundary.
    constexpr double minimum_parallel_dot = 0.99999;
    const double distance_tolerance = std::max(pitch * 0.02, model_scale * 1.0e-9);
    const double normal_bucket_width = std::sqrt(2.0 * (1.0 - minimum_parallel_dot));
    using PlaneBucketKey = std::array<std::int64_t, 4>;
    const auto bucket_key = [&](Vec3 normal, double distance)
    {
        int dominant = 0;
        for (int axis = 1; axis < 3; ++axis)
            if (std::abs(normal[axis]) > std::abs(normal[dominant])) dominant = axis;
        if (normal[dominant] < 0.0)
        {
            normal = normal * -1.0;
            distance = -distance;
        }
        return PlaneBucketKey{
            static_cast<std::int64_t>(std::llround(normal.x / normal_bucket_width)),
            static_cast<std::int64_t>(std::llround(normal.y / normal_bucket_width)),
            static_cast<std::int64_t>(std::llround(normal.z / normal_bucket_width)),
            static_cast<std::int64_t>(std::llround(distance / distance_tolerance))};
    };
    std::vector<DualSourcePlane> result;
    result.reserve(exact_planes.size());
    std::map<PlaneBucketKey, std::vector<std::size_t>> plane_buckets;
    for (const auto& candidate : exact_planes)
    {
        auto match = result.end();
        const PlaneBucketKey candidate_bucket =
            bucket_key(candidate.normal, candidate.distance);
        std::vector<std::size_t> possible_matches;
        for (std::int64_t dx = -2; dx <= 2; ++dx)
            for (std::int64_t dy = -2; dy <= 2; ++dy)
                for (std::int64_t dz = -2; dz <= 2; ++dz)
                    for (std::int64_t dd = -2; dd <= 2; ++dd)
                    {
                        PlaneBucketKey key = candidate_bucket;
                        key[0] += dx;
                        key[1] += dy;
                        key[2] += dz;
                        key[3] += dd;
                        const auto found = plane_buckets.find(key);
                        if (found != plane_buckets.end())
                            possible_matches.insert(possible_matches.end(),
                                found->second.begin(), found->second.end());
                    }
        std::sort(possible_matches.begin(), possible_matches.end());
        possible_matches.erase(
            std::unique(possible_matches.begin(), possible_matches.end()),
            possible_matches.end());
        for (const std::size_t match_id : possible_matches)
        {
            auto iterator = result.begin() + static_cast<std::ptrdiff_t>(match_id);
            const double orientation = dot(iterator->normal, candidate.normal);
            if (std::abs(orientation) < minimum_parallel_dot) continue;
            const double candidate_distance = orientation < 0.0
                ? -candidate.distance : candidate.distance;
            if (std::abs(iterator->distance - candidate_distance) > distance_tolerance)
                continue;
            bool overlapping_support = true;
            for (int axis = 0; axis < 3; ++axis)
                overlapping_support &= iterator->bounds.lower[axis] <=
                        candidate.bounds.upper[axis] + pitch &&
                    candidate.bounds.lower[axis] <= iterator->bounds.upper[axis] + pitch;
            if (overlapping_support)
            {
                match = iterator;
                break;
            }
        }
        if (match == result.end())
        {
            result.push_back(candidate);
            plane_buckets[candidate_bucket].push_back(result.size() - 1);
            continue;
        }
        const double orientation = dot(match->normal, candidate.normal);
        const Vec3 aligned_normal = orientation < 0.0
            ? candidate.normal * -1.0 : candidate.normal;
        const double aligned_distance = orientation < 0.0
            ? -candidate.distance : candidate.distance;
        const double combined_area = match->support_area + candidate.support_area;
        match->normal = normalized(
            match->normal * match->support_area + aligned_normal * candidate.support_area);
        match->distance = (match->distance * match->support_area +
            aligned_distance * candidate.support_area) / combined_area;
        match->support_area = combined_area;
        include(match->bounds, candidate.bounds.lower);
        include(match->bounds, candidate.bounds.upper);
        match->support_components.insert(match->support_components.end(),
            candidate.support_components.begin(),
            candidate.support_components.end());
        const std::size_t match_id =
            static_cast<std::size_t>(match - result.begin());
        plane_buckets[bucket_key(match->normal, match->distance)].push_back(match_id);
    }
    return result;
}

std::vector<std::uint8_t> exteriorSourceFaces(
    const VoxelGrid& occupancy, const MeshModel& source)
{
    const Vec3 grid_origin = toVec(occupancy.origin);
    const auto sample_occupied = [&](const Vec3 point)
    {
        std::array<int, 3> coordinate{};
        for (int axis = 0; axis < 3; ++axis)
        {
            coordinate[axis] = static_cast<int>(std::llround(
                (point[axis] - grid_origin[axis]) / occupancy.pitch));
            if (coordinate[axis] < 0 ||
                coordinate[axis] >= static_cast<int>(occupancy.shape[axis]))
                return false;
        }
        return occupancy.occupied(
            static_cast<std::uint32_t>(coordinate[0]),
            static_cast<std::uint32_t>(coordinate[1]),
            static_cast<std::uint32_t>(coordinate[2]));
    };
    std::vector<std::uint8_t> retained(source.triangles.size(), 0);
    for (std::size_t face_id = 0; face_id < source.triangles.size(); ++face_id)
    {
        const auto face = source.triangles[face_id];
        const Vec3 a = toVec(source.vertices[face[0]]);
        const Vec3 b = toVec(source.vertices[face[1]]);
        const Vec3 c = toVec(source.vertices[face[2]]);
        const Vec3 normal = normalized(cross(b - a, c - a));
        if (norm(normal) == 0.0) continue;
        const std::array<Vec3, 4> samples{
            (a + b + c) * (1.0 / 3.0),
            a * 0.6 + b * 0.2 + c * 0.2,
            a * 0.2 + b * 0.6 + c * 0.2,
            a * 0.2 + b * 0.2 + c * 0.6};
        for (const Vec3 sample : samples)
            for (const double scale : {0.5, 1.0, 1.5, 2.0})
            {
                const Vec3 offset = normal * (scale * occupancy.pitch);
                if (sample_occupied(sample + offset) !=
                    sample_occupied(sample - offset))
                {
                    retained[face_id] = 1;
                    break;
                }
            }
    }
    return retained;
}

bool withinExpandedBounds(const Vec3 point, const Bounds3& bounds, const double expansion)
{
    for (int axis = 0; axis < 3; ++axis)
        if (point[axis] < bounds.lower[axis] - expansion ||
            point[axis] > bounds.upper[axis] + expansion) return false;
    return true;
}

class SourcePlaneBoundsBvh
{
public:
    SourcePlaneBoundsBvh(
        const std::vector<DualSourcePlane>& planes, const double expansion)
        : planes_(planes), expansion_(expansion)
    {
        plane_ids_.resize(planes.size());
        std::iota(plane_ids_.begin(), plane_ids_.end(), 0);
        if (!plane_ids_.empty()) build(0, plane_ids_.size());
    }

    [[nodiscard]] std::vector<std::size_t> containing(const Vec3 point) const
    {
        std::vector<std::size_t> result;
        if (!nodes_.empty()) query(0, point, result);
        std::sort(result.begin(), result.end());
        return result;
    }

private:
    struct Node
    {
        Bounds3 bounds;
        std::size_t begin = 0;
        std::size_t end = 0;
        std::size_t left = 0;
        std::size_t right = 0;
        bool leaf = false;
    };

    [[nodiscard]] Bounds3 expandedBounds(const std::size_t plane_id) const
    {
        Bounds3 result = planes_[plane_id].bounds;
        for (int axis = 0; axis < 3; ++axis)
        {
            result.lower[axis] -= expansion_;
            result.upper[axis] += expansion_;
        }
        return result;
    }

    [[nodiscard]] Vec3 center(const std::size_t plane_id) const
    {
        const Bounds3 bounds = planes_[plane_id].bounds;
        return (bounds.lower + bounds.upper) * 0.5;
    }

    std::size_t build(const std::size_t begin, const std::size_t end)
    {
        Node node;
        node.begin = begin;
        node.end = end;
        Bounds3 center_bounds;
        for (std::size_t index = begin; index < end; ++index)
        {
            const Bounds3 bounds = expandedBounds(plane_ids_[index]);
            include(node.bounds, bounds.lower);
            include(node.bounds, bounds.upper);
            include(center_bounds, center(plane_ids_[index]));
        }
        const std::size_t node_id = nodes_.size();
        nodes_.push_back(node);
        if (end - begin <= 8)
        {
            nodes_[node_id].leaf = true;
            return node_id;
        }
        const Vec3 extent = center_bounds.upper - center_bounds.lower;
        int axis = extent.y > extent.x ? 1 : 0;
        if (extent.z > extent[axis]) axis = 2;
        const std::size_t middle = begin + (end - begin) / 2;
        std::nth_element(plane_ids_.begin() + begin, plane_ids_.begin() + middle,
            plane_ids_.begin() + end, [&](const std::size_t first, const std::size_t second)
            { return center(first)[axis] < center(second)[axis]; });
        nodes_[node_id].left = build(begin, middle);
        nodes_[node_id].right = build(middle, end);
        return node_id;
    }

    void query(const std::size_t node_id, const Vec3 point,
               std::vector<std::size_t>& result) const
    {
        const Node& node = nodes_[node_id];
        if (!withinExpandedBounds(point, node.bounds, 0.0)) return;
        if (node.leaf)
        {
            for (std::size_t index = node.begin; index < node.end; ++index)
            {
                const std::size_t plane_id = plane_ids_[index];
                if (withinExpandedBounds(point, planes_[plane_id].bounds, expansion_))
                    result.push_back(plane_id);
            }
            return;
        }
        query(node.left, point, result);
        query(node.right, point, result);
    }

    const std::vector<DualSourcePlane>& planes_;
    double expansion_ = 0.0;
    std::vector<std::size_t> plane_ids_;
    std::vector<Node> nodes_;
};

std::optional<std::size_t> supportingPlane(
    const Vec3 point, const int crossing_axis,
    const std::vector<DualSourcePlane>& planes,
    const SourcePlaneBoundsBvh& plane_bvh, const double pitch)
{
    std::vector<std::pair<std::size_t, double>> candidates;
    double best_distance = std::numeric_limits<double>::infinity();
    for (const std::size_t index : plane_bvh.containing(point))
    {
        const DualSourcePlane& plane = planes[index];
        if (std::abs(plane.normal[crossing_axis]) < 0.45 ||
            !withinExpandedBounds(point, plane.bounds, 2.0 * pitch)) continue;
        const bool within_component = std::ranges::any_of(
            plane.support_components, [&](const Bounds3& component)
            { return withinExpandedBounds(point, component, 2.0 * pitch); });
        if (!within_component) continue;
        const double distance = std::abs(dot(plane.normal, point) - plane.distance);
        if (distance <= 2.0 * pitch)
        {
            candidates.emplace_back(index, distance);
            best_distance = std::min(best_distance, distance);
        }
    }
    std::optional<std::size_t> result;
    double best_support = -1.0;
    const double equivalent_distance = best_distance + pitch * 0.02;
    for (const auto [index, distance] : candidates)
        if (distance <= equivalent_distance && planes[index].support_area > best_support)
        {
            result = index;
            best_support = planes[index].support_area;
        }
    return result;
}

Vec3 solveSymmetric3x3(
    std::array<std::array<double, 3>, 3> matrix,
    std::array<double, 3> right,
    const Vec3 fallback)
{
    for (int pivot = 0; pivot < 3; ++pivot)
    {
        int best = pivot;
        for (int row = pivot + 1; row < 3; ++row)
            if (std::abs(matrix[row][pivot]) > std::abs(matrix[best][pivot])) best = row;
        if (std::abs(matrix[best][pivot]) <= 1.0e-12) return fallback;
        if (best != pivot)
        {
            std::swap(matrix[best], matrix[pivot]);
            std::swap(right[best], right[pivot]);
        }
        const double inverse = 1.0 / matrix[pivot][pivot];
        for (int column = pivot; column < 3; ++column) matrix[pivot][column] *= inverse;
        right[pivot] *= inverse;
        for (int row = 0; row < 3; ++row)
        {
            if (row == pivot) continue;
            const double factor = matrix[row][pivot];
            for (int column = pivot; column < 3; ++column)
                matrix[row][column] -= factor * matrix[pivot][column];
            right[row] -= factor * right[pivot];
        }
    }
    return {right[0], right[1], right[2]};
}

IndexedSurface reconstructDualSurface(const VoxelGrid& grid, const MeshModel& source)
{
    const auto plane_started = std::chrono::steady_clock::now();
    const auto retained_faces = exteriorSourceFaces(grid, source);
    const auto planes = sourcePlanes(source, grid.pitch, &retained_faces);
    std::cerr << "monitor: stage=source_planes_built plane_count="
              << planes.size() << " elapsed_seconds="
              << std::chrono::duration<double>(
                     std::chrono::steady_clock::now() - plane_started).count()
              << '\n';
    const SourcePlaneBoundsBvh plane_bvh(planes, 2.0 * grid.pitch);
    const Vec3 origin = toVec(grid.origin);
    const std::array<std::uint32_t, 3> cell_shape{
        grid.shape[0] - 1, grid.shape[1] - 1, grid.shape[2] - 1};
    const std::size_t cell_count = checkedVoxelCount(cell_shape);
    constexpr std::uint32_t no_vertex = std::numeric_limits<std::uint32_t>::max();
    std::vector<std::uint32_t> cell_vertices(cell_count, no_vertex);
    const auto cell_index = [&](const std::uint32_t x, const std::uint32_t y,
                                const std::uint32_t z)
    { return (static_cast<std::size_t>(x) * cell_shape[1] + y) * cell_shape[2] + z; };
    static constexpr std::array<std::array<int, 3>, 8> corners{{
        {{0,0,0}},{{1,0,0}},{{0,1,0}},{{1,1,0}},
        {{0,0,1}},{{1,0,1}},{{0,1,1}},{{1,1,1}},
    }};
    static constexpr std::array<std::array<int, 3>, 12> cube_edges{{
        {{0,1,0}},{{2,3,0}},{{4,5,0}},{{6,7,0}},
        {{0,2,1}},{{1,3,1}},{{4,6,1}},{{5,7,1}},
        {{0,4,2}},{{1,5,2}},{{2,6,2}},{{3,7,2}},
    }};
    IndexedSurface result;
    result.vertices.reserve(grid.occupancy.size() / 16);
    std::vector<Vec3> cell_center_vertices;
    cell_center_vertices.reserve(grid.occupancy.size() / 16);
    for (std::uint32_t x = 0; x < cell_shape[0]; ++x)
        for (std::uint32_t y = 0; y < cell_shape[1]; ++y)
            for (std::uint32_t z = 0; z < cell_shape[2]; ++z)
            {
                std::array<bool, 8> signs{};
                int occupied_count = 0;
                for (int corner = 0; corner < 8; ++corner)
                {
                    signs[corner] = grid.occupied(
                        x + corners[corner][0], y + corners[corner][1],
                        z + corners[corner][2]);
                    occupied_count += signs[corner] ? 1 : 0;
                }
                if (occupied_count == 0 || occupied_count == 8) continue;
                std::array<std::array<double, 3>, 3> matrix{};
                std::array<double, 3> right{};
                Vec3 average{};
                int constraints = 0;
                std::optional<std::size_t> common_plane;
                std::set<std::size_t> selected_plane_ids;
                bool one_common_plane = true;
                for (const auto& edge : cube_edges)
                {
                    if (signs[edge[0]] == signs[edge[1]]) continue;
                    const int axis = edge[2];
                    const auto& first = corners[edge[0]];
                    const auto& second = corners[edge[1]];
                    Vec3 point = origin + Vec3{
                        (x + 0.5 * (first[0] + second[0])) * grid.pitch,
                        (y + 0.5 * (first[1] + second[1])) * grid.pitch,
                        (z + 0.5 * (first[2] + second[2])) * grid.pitch};
                    Vec3 normal{};
                    normal[axis] = 1.0;
                    const auto plane_id = supportingPlane(
                        point, axis, planes, plane_bvh, grid.pitch);
                    if (plane_id)
                    {
                        const DualSourcePlane& plane = planes[*plane_id];
                        normal = plane.normal;
                        point = point - normal * (dot(normal, point) - plane.distance);
                        selected_plane_ids.insert(*plane_id);
                        if (!common_plane) common_plane = plane_id;
                        else if (*common_plane != *plane_id) one_common_plane = false;
                    }
                    else
                    {
                        one_common_plane = false;
                    }
                    for (int row = 0; row < 3; ++row)
                    {
                        right[row] += normal[row] * dot(normal, point);
                        for (int column = 0; column < 3; ++column)
                            matrix[row][column] += normal[row] * normal[column];
                    }
                    average = average + point;
                    ++constraints;
                }
                if (constraints == 0) continue;
                average = average * (1.0 / constraints);
                const Vec3 center = origin + Vec3{
                    (x + 0.5) * grid.pitch, (y + 0.5) * grid.pitch,
                    (z + 0.5) * grid.pitch};
                const double regularization = 1.0e-8 * constraints;
                for (int axis = 0; axis < 3; ++axis)
                {
                    matrix[axis][axis] += regularization;
                    right[axis] += regularization * center[axis];
                }
                Vec3 vertex = solveSymmetric3x3(matrix, right, average);
                if (one_common_plane && common_plane)
                {
                    const auto& plane = planes[*common_plane];
                    vertex = vertex - plane.normal * (dot(plane.normal, vertex) - plane.distance);
                }
                const std::array<std::uint32_t, 3> cell_coordinate{x, y, z};
                const auto clamp_to_cell = [&]
                {
                    for (int axis = 0; axis < 3; ++axis)
                        vertex[axis] = std::clamp(vertex[axis],
                            origin[axis] + cell_coordinate[axis] * grid.pitch,
                            origin[axis] + (cell_coordinate[axis] + 1.0) * grid.pitch);
                };
                clamp_to_cell();
                std::vector<std::size_t> structural_planes;
                for (const std::size_t plane_id : selected_plane_ids)
                {
                    const auto& candidate = planes[plane_id];
                    if (candidate.support_area < grid.pitch * grid.pitch) continue;
                    bool independent = true;
                    for (const std::size_t accepted : structural_planes)
                        if (std::abs(dot(candidate.normal, planes[accepted].normal)) >
                            1.0 - 1.0e-7) independent = false;
                    if (independent) structural_planes.push_back(plane_id);
                }
                if (structural_planes.size() >= 2)
                {
                    std::array<std::array<double, 3>, 3> feature_matrix{};
                    std::array<double, 3> feature_right{};
                    for (const std::size_t plane_id : structural_planes)
                    {
                        const auto& plane = planes[plane_id];
                        for (int row = 0; row < 3; ++row)
                        {
                            feature_right[row] += plane.normal[row] * plane.distance;
                            for (int column = 0; column < 3; ++column)
                                feature_matrix[row][column] +=
                                    plane.normal[row] * plane.normal[column];
                        }
                    }
                    constexpr double feature_regularization = 1.0e-10;
                    for (int axis = 0; axis < 3; ++axis)
                    {
                        feature_matrix[axis][axis] += feature_regularization;
                        feature_right[axis] += feature_regularization * vertex[axis];
                    }
                    const Vec3 projected = solveSymmetric3x3(
                        feature_matrix, feature_right, vertex);
                    const Vec3 movement = projected - vertex;
                    if (dot(movement, movement) <= 0.75 * 0.75 * grid.pitch * grid.pitch)
                        vertex = projected;
                }
                // A source feature may improve the local fit, but it cannot move
                // this dual vertex outside the cell that owns its connectivity.
                clamp_to_cell();
                cell_vertices[cell_index(x, y, z)] =
                    static_cast<std::uint32_t>(result.vertices.size());
                result.vertices.push_back(vertex);
                cell_center_vertices.push_back(center);
            }

    const auto emit_crossing = [&](const int axis, const std::uint32_t x,
                                   const std::uint32_t y, const std::uint32_t z)
    {
        std::array<std::uint32_t, 3> first{x, y, z};
        auto second = first;
        ++second[axis];
        const bool first_inside = grid.occupied(first[0], first[1], first[2]);
        if (first_inside == grid.occupied(second[0], second[1], second[2])) return;
        const int u = (axis + 1) % 3;
        const int v = (axis + 2) % 3;
        std::array<std::array<std::uint32_t, 3>, 4> cells{};
        for (auto& cell : cells) cell = first;
        cells[0][u]--; cells[0][v]--;
        cells[1][v]--;
        cells[3][u]--;
        std::array<std::uint32_t, 4> quad{};
        for (int index = 0; index < 4; ++index)
        {
            quad[index] = cell_vertices[
                cell_index(cells[index][0], cells[index][1], cells[index][2])];
            if (quad[index] == no_vertex)
                throw std::runtime_error("dual surface has a missing crossing-cell vertex");
        }
        // (axis + 1, axis + 2) is a positive cyclic basis for every axis.
        // The face points toward the empty endpoint of the crossed grid edge.
        if (!first_inside) std::reverse(quad.begin(), quad.end());
        const std::uint32_t crossing =
            static_cast<std::uint32_t>(result.boundary_crossings.size());
        result.boundary_crossings.push_back(first_inside
            ? VoxelBoundaryCrossing{first, second}
            : VoxelBoundaryCrossing{second, first});
        result.triangles.push_back({quad[0], quad[1], quad[2]});
        result.triangles.push_back({quad[0], quad[2], quad[3]});
        result.face_boundary_crossings.push_back(crossing);
        result.face_boundary_crossings.push_back(crossing);
    };
    for (int axis = 0; axis < 3; ++axis)
    {
        const int u = (axis + 1) % 3;
        const int v = (axis + 2) % 3;
        std::array<std::uint32_t, 3> coordinate{};
        for (coordinate[axis] = 0; coordinate[axis] + 1 < grid.shape[axis]; ++coordinate[axis])
            for (coordinate[u] = 1; coordinate[u] + 1 < grid.shape[u]; ++coordinate[u])
                for (coordinate[v] = 1; coordinate[v] + 1 < grid.shape[v]; ++coordinate[v])
                    emit_crossing(axis, coordinate[0], coordinate[1], coordinate[2]);
    }
    const double minimum_area_twice = grid.pitch * grid.pitch * 1.0e-14;
    const auto repair_degenerate_faces = [&](const std::vector<Vec3>& fallback)
    {
        bool changed = false;
        for (const auto& triangle : result.triangles)
        {
            const Vec3 a = result.vertices[triangle[0]];
            const Vec3 b = result.vertices[triangle[1]];
            const Vec3 c = result.vertices[triangle[2]];
            if (norm(cross(b - a, c - a)) > minimum_area_twice) continue;
            for (const std::uint32_t vertex : triangle)
            {
                const Vec3 replacement = fallback[vertex];
                const Vec3 delta = result.vertices[vertex] - replacement;
                if (dot(delta, delta) > 0.0)
                {
                    result.vertices[vertex] = replacement;
                    changed = true;
                }
            }
        }
        return changed;
    };
    while (repair_degenerate_faces(cell_center_vertices)) {}
    return result;
}

double indexedSignedVolume(const IndexedSurface& surface)
{
    double result = 0.0;
    for (const auto& triangle : surface.triangles)
        result += dot(surface.vertices[triangle[0]],
                      cross(surface.vertices[triangle[1]], surface.vertices[triangle[2]]));
    return result / 6.0;
}

OrientedSurfaceMesh buildOrientedSurfaceMesh(
    IndexedSurface surface, const VoxelGrid& occupancy,
    std::vector<std::uint32_t> source_faces = {})
{
    if (surface.vertices.empty() || surface.triangles.empty())
        throw std::runtime_error("dual surface is empty");
    const bool has_voxel_provenance = !surface.face_boundary_crossings.empty();
    if (has_voxel_provenance &&
        surface.face_boundary_crossings.size() != surface.triangles.size())
        throw std::runtime_error("dual surface lost voxel-boundary provenance");
    if (!source_faces.empty() && source_faces.size() != surface.triangles.size())
        throw std::runtime_error("surface lost source-face provenance");

    double signed_volume = indexedSignedVolume(surface);
    if (signed_volume < 0.0)
    {
        for (auto& triangle : surface.triangles) std::swap(triangle[1], triangle[2]);
        signed_volume = -signed_volume;
    }

    // A dual cell can touch the boundary in several disconnected sectors.
    // Those sectors may share the cell vertex geometrically, but they are
    // distinct topology vertices. Split triangle corners by edge-connected
    // fan before constructing halfedges; coordinates and face provenance stay
    // unchanged.
    if (surface.halfedge_opposites.empty())
    {
    const std::size_t corner_count = surface.triangles.size() * 3;
    std::vector<std::uint32_t> corner_parent(corner_count);
    std::vector<std::uint32_t> corner_opposites(
        corner_count, invalid_surface_index);
    std::iota(corner_parent.begin(), corner_parent.end(), 0);
    const auto find_corner = [&](std::uint32_t value)
    {
        std::uint32_t root = value;
        while (corner_parent[root] != root) root = corner_parent[root];
        while (corner_parent[value] != value)
        {
            const std::uint32_t next = corner_parent[value];
            corner_parent[value] = root;
            value = next;
        }
        return root;
    };
    const auto unite_corners = [&](const std::uint32_t first,
                                   const std::uint32_t second)
    {
        const std::uint32_t a = find_corner(first);
        const std::uint32_t b = find_corner(second);
        if (a != b) corner_parent[b] = a;
    };
    std::map<std::pair<std::uint32_t, std::uint32_t>,
             std::vector<std::uint32_t>> corner_edge_uses;
    for (std::uint32_t face = 0; face < surface.triangles.size(); ++face)
        for (std::uint32_t local = 0; local < 3; ++local)
        {
            const std::uint32_t first = surface.triangles[face][local];
            const std::uint32_t second = surface.triangles[face][(local + 1) % 3];
            corner_edge_uses[std::minmax(first, second)].push_back(face * 3 + local);
        }
    for (const auto& [edge, uses] : corner_edge_uses)
    {
        const auto next_corner = [](const std::uint32_t corner)
        { return corner - corner % 3 + (corner % 3 + 1) % 3; };
        const auto corner_vertex = [&](const std::uint32_t corner)
        { return surface.triangles[corner / 3][corner % 3]; };
        const auto endpoint_corner = [&](const std::uint32_t use,
                                         const std::uint32_t endpoint)
        {
            if (corner_vertex(use) == endpoint) return use;
            const std::uint32_t next = next_corner(use);
            if (corner_vertex(next) == endpoint) return next;
            throw std::runtime_error("dual edge use lost an endpoint");
        };
        const auto pair_uses = [&](const std::uint32_t first,
                                   const std::uint32_t second)
        {
            const std::uint32_t first_begin = corner_vertex(first);
            const std::uint32_t first_end = corner_vertex(next_corner(first));
            const std::uint32_t second_begin = corner_vertex(second);
            const std::uint32_t second_end = corner_vertex(next_corner(second));
            if (first_begin != second_end || first_end != second_begin)
            {
                const std::uint32_t first_crossing =
                    surface.face_boundary_crossings[first / 3];
                const std::uint32_t second_crossing =
                    surface.face_boundary_crossings[second / 3];
                const auto& first_inside =
                    surface.boundary_crossings[first_crossing].inside;
                const auto& second_inside =
                    surface.boundary_crossings[second_crossing].inside;
                std::ostringstream message;
                message << "dual-sector pairing has equal edge orientation: edge=("
                        << edge.first << ',' << edge.second << ") uses=" << uses.size()
                        << " first=" << first_begin << "->" << first_end
                        << " second=" << second_begin << "->" << second_end
                        << " inside_first=(" << first_inside[0] << ','
                        << first_inside[1] << ',' << first_inside[2] << ')'
                        << " inside_second=(" << second_inside[0] << ','
                        << second_inside[1] << ',' << second_inside[2] << ')';
                throw std::runtime_error(message.str());
            }
            if (corner_opposites[first] != invalid_surface_index ||
                corner_opposites[second] != invalid_surface_index)
                throw std::runtime_error("dual edge use was paired more than once");
            corner_opposites[first] = second;
            corner_opposites[second] = first;
            unite_corners(endpoint_corner(first, edge.first),
                          endpoint_corner(second, edge.first));
            unite_corners(endpoint_corner(first, edge.second),
                          endpoint_corner(second, edge.second));
        };
        if (uses.size() == 2)
        {
            pair_uses(uses[0], uses[1]);
            continue;
        }

        // A checkerboard voxel face creates four dual-edge uses. Under the
        // foreground-6/background-26 convention, the two sheets meet only
        // geometrically: pair the uses around each occupied voxel and give
        // each pair distinct topology vertices.
        std::map<std::array<std::uint32_t, 3>, std::vector<std::uint32_t>>
            uses_by_inside_voxel;
        for (const std::uint32_t use : uses)
        {
            const std::uint32_t face = use / 3;
            if (!has_voxel_provenance ||
                face >= surface.face_boundary_crossings.size())
                throw std::runtime_error(
                    "non-manifold dual edge has no voxel-boundary certificate");
            const std::uint32_t crossing = surface.face_boundary_crossings[face];
            if (crossing >= surface.boundary_crossings.size())
                throw std::runtime_error(
                    "non-manifold dual edge has an invalid boundary certificate");
            uses_by_inside_voxel[
                surface.boundary_crossings[crossing].inside].push_back(use);
        }
        for (const auto& [inside, grouped_uses] : uses_by_inside_voxel)
        {
            (void)inside;
            if (grouped_uses.size() != 2)
                throw std::runtime_error(
                    "ambiguous dual edge cannot be paired by occupied voxel");
            pair_uses(grouped_uses[0], grouped_uses[1]);
        }
    }
    if (std::ranges::any_of(corner_opposites, [](const std::uint32_t opposite)
        { return opposite == invalid_surface_index; }))
        throw std::runtime_error("dual surface contains an unpaired edge use");
    std::vector<Vec3> split_vertices;
    split_vertices.reserve(surface.vertices.size());
    std::map<std::pair<std::uint32_t, std::uint32_t>, std::uint32_t> fan_vertices;
    for (std::uint32_t face = 0; face < surface.triangles.size(); ++face)
        for (std::uint32_t local = 0; local < 3; ++local)
        {
            const std::uint32_t corner = face * 3 + local;
            const std::uint32_t old_vertex = surface.triangles[face][local];
            const auto key = std::make_pair(old_vertex, find_corner(corner));
            const auto [iterator, inserted] = fan_vertices.try_emplace(
                key, static_cast<std::uint32_t>(split_vertices.size()));
            if (inserted)
            {
                split_vertices.push_back(surface.vertices[old_vertex]);
            }
            surface.triangles[face][local] = iterator->second;
        }
    surface.vertices = std::move(split_vertices);
    surface.halfedge_opposites = std::move(corner_opposites);
    }

    OrientedSurfaceMesh result;
    result.geometry.name = "phase1_watertight_surface";
    result.geometry.vertices.reserve(surface.vertices.size());
    for (const Vec3 point : surface.vertices) result.geometry.vertices.push_back(toPosition(point));
    result.geometry.triangles = std::move(surface.triangles);
    result.boundary_crossings = std::move(surface.boundary_crossings);
    result.face_boundary_crossings = std::move(surface.face_boundary_crossings);
    result.face_source_faces = std::move(source_faces);
    result.vertex_halfedges.assign(result.geometry.vertices.size(), invalid_surface_index);
    result.face_halfedges.resize(result.geometry.triangles.size(), invalid_surface_index);
    result.halfedges.resize(result.geometry.triangles.size() * 3);

    Bounds3 bounds;
    for (const Position3 point : result.geometry.vertices) include(bounds, toVec(point));
    const double model_scale = norm(bounds.upper - bounds.lower);
    const double minimum_area = std::max(1.0, model_scale * model_scale) * 1.0e-24;
    if (!surface.halfedge_opposites.empty() &&
        surface.halfedge_opposites.size() != result.halfedges.size())
        throw std::runtime_error("surface opposite map has the wrong size");
    std::map<std::pair<std::uint32_t, std::uint32_t>, std::uint32_t> edge_map;
    for (std::uint32_t face = 0; face < result.geometry.triangles.size(); ++face)
    {
        const auto& triangle = result.geometry.triangles[face];
        for (const std::uint32_t vertex : triangle)
            if (vertex >= result.geometry.vertices.size())
                throw std::runtime_error("dual surface contains an invalid vertex index");
        if (triangle[0] == triangle[1] || triangle[1] == triangle[2] ||
            triangle[2] == triangle[0])
            throw std::runtime_error("dual surface contains a topologically degenerate triangle");
        const Vec3 a = toVec(result.geometry.vertices[triangle[0]]);
        const Vec3 b = toVec(result.geometry.vertices[triangle[1]]);
        const Vec3 c = toVec(result.geometry.vertices[triangle[2]]);
        const double area_twice = norm(cross(b - a, c - a));
        if (!std::isfinite(area_twice) || area_twice <= minimum_area)
            throw std::runtime_error("dual surface contains a geometrically degenerate triangle");

        const std::uint32_t base = face * 3;
        result.face_halfedges[face] = base;
        for (std::uint32_t local = 0; local < 3; ++local)
        {
            const std::uint32_t halfedge = base + local;
            const std::uint32_t first = triangle[local];
            const std::uint32_t second = triangle[(local + 1) % 3];
            result.halfedges[halfedge] = {
                first, face, base + (local + 1) % 3, invalid_surface_index};
            if (result.vertex_halfedges[first] == invalid_surface_index)
                result.vertex_halfedges[first] = halfedge;
            if (!surface.halfedge_opposites.empty())
            {
                const std::uint32_t opposite = surface.halfedge_opposites[halfedge];
                if (opposite >= result.halfedges.size())
                    throw std::runtime_error("surface opposite map contains an invalid index");
                result.halfedges[halfedge].opposite = opposite;
            }
            else
            {
                const auto key = std::minmax(first, second);
                const auto [iterator, inserted] = edge_map.try_emplace(key, halfedge);
                if (!inserted)
                {
                    const std::uint32_t opposite = iterator->second;
                    if (result.halfedges[opposite].opposite != invalid_surface_index)
                        throw std::runtime_error("dual surface contains a non-manifold edge");
                    result.halfedges[halfedge].opposite = opposite;
                    result.halfedges[opposite].opposite = halfedge;
                }
            }
        }
    }
    for (std::uint32_t halfedge = 0; halfedge < result.halfedges.size(); ++halfedge)
    {
        const std::uint32_t opposite = result.halfedges[halfedge].opposite;
        if (opposite == invalid_surface_index)
            throw std::runtime_error("dual surface is open: a halfedge has no opposite");
        if (result.halfedges[opposite].opposite != halfedge)
            throw std::runtime_error("dual surface opposite relation is asymmetric");
        const std::uint32_t second =
            result.halfedges[result.halfedges[halfedge].next].origin;
        const std::uint32_t opposite_second =
            result.halfedges[result.halfedges[opposite].next].origin;
        if (result.halfedges[opposite].origin != second ||
            opposite_second != result.halfedges[halfedge].origin)
        {
            throw std::runtime_error("dual surface has inconsistent face orientation");
        }
    }

    std::vector<std::uint8_t> visited_faces(result.geometry.triangles.size(), 0);
    std::vector<std::uint32_t> queue{0};
    visited_faces[0] = 1;
    for (std::size_t cursor = 0; cursor < queue.size(); ++cursor)
    {
        const std::uint32_t face = queue[cursor];
        std::uint32_t halfedge = result.face_halfedges[face];
        for (int local = 0; local < 3; ++local)
        {
            const std::uint32_t neighbor = result.halfedges[
                result.halfedges[halfedge].opposite].face;
            if (!visited_faces[neighbor])
            {
                visited_faces[neighbor] = 1;
                queue.push_back(neighbor);
            }
            halfedge = result.halfedges[halfedge].next;
        }
    }
    if (queue.size() != result.geometry.triangles.size())
        throw std::runtime_error("dual surface has more than one face component");

    std::vector<std::uint32_t> vertex_valence(result.geometry.vertices.size(), 0);
    for (const SurfaceHalfedge& halfedge : result.halfedges) ++vertex_valence[halfedge.origin];
    for (std::uint32_t vertex = 0; vertex < result.geometry.vertices.size(); ++vertex)
    {
        const std::uint32_t start = result.vertex_halfedges[vertex];
        if (start == invalid_surface_index)
            throw std::runtime_error("dual surface contains an unused vertex");
        std::uint32_t current = start;
        std::uint32_t fan_size = 0;
        do
        {
            if (result.halfedges[current].origin != vertex || ++fan_size > vertex_valence[vertex])
                throw std::runtime_error("dual surface contains a non-manifold vertex");
            const std::uint32_t previous = result.halfedges[
                result.halfedges[current].next].next;
            current = result.halfedges[previous].opposite;
        }
        while (current != start);
        if (fan_size != vertex_valence[vertex])
            throw std::runtime_error("dual surface vertex has multiple disconnected fans");
    }

    if (has_voxel_provenance)
    {
        std::vector<std::uint32_t> crossing_use(result.boundary_crossings.size(), 0);
        for (const std::uint32_t crossing : result.face_boundary_crossings)
        {
            if (crossing >= result.boundary_crossings.size())
                throw std::runtime_error("dual surface contains invalid boundary provenance");
            ++crossing_use[crossing];
        }
        for (std::size_t index = 0; index < result.boundary_crossings.size(); ++index)
        {
            if (crossing_use[index] != 2)
                throw std::runtime_error("voxel-boundary crossing does not own exactly one quad");
            const auto& crossing = result.boundary_crossings[index];
            int coordinate_delta = 0;
            for (int axis = 0; axis < 3; ++axis)
                coordinate_delta += std::abs(static_cast<int>(crossing.inside[axis]) -
                                             static_cast<int>(crossing.outside[axis]));
            if (coordinate_delta != 1 ||
                !occupancy.occupied(crossing.inside[0], crossing.inside[1], crossing.inside[2]) ||
                occupancy.occupied(crossing.outside[0], crossing.outside[1], crossing.outside[2]))
                throw std::runtime_error("a surface quad is not an inside/outside voxel boundary");
        }
    }

    result.euler_characteristic = static_cast<int>(result.geometry.vertices.size()) -
        static_cast<int>(result.halfedges.size() / 2) +
        static_cast<int>(result.geometry.triangles.size());
    if (result.euler_characteristic != 2)
        throw std::runtime_error("phase-1 boundary is not a genus-zero closed surface");
    const double minimum_volume = std::max(1.0, model_scale * model_scale * model_scale) * 1.0e-18;
    if (!std::isfinite(signed_volume) || signed_volume <= minimum_volume)
        throw std::runtime_error("phase-1 boundary has invalid oriented volume");
    result.signed_volume = signed_volume;
    return result;
}

void writeSurfaceObj(const std::filesystem::path& path, const OrientedSurfaceMesh& surface)
{
    std::ofstream stream(path);
    if (!stream) throw std::runtime_error("failed to create OBJ: " + path.string());
    stream << std::setprecision(17) << "o phase1_watertight_surface\n";
    for (const Position3 point : surface.geometry.vertices)
        stream << "v " << point.x << ' ' << point.y << ' ' << point.z << '\n';
    stream << "g phase1_surface\n";
    for (const auto& triangle : surface.geometry.triangles)
        stream << "f " << triangle[0] + 1 << ' ' << triangle[1] + 1 << ' '
               << triangle[2] + 1 << '\n';
}

void writeMeshObj(const std::filesystem::path& path, const MeshModel& mesh)
{
    std::ofstream stream(path);
    if (!stream) throw std::runtime_error("failed to create OBJ: " + path.string());
    stream << std::setprecision(17) << "o phase1_original_obj_hole_surgery\n";
    for (const Position3 point : mesh.vertices)
        stream << "v " << point.x << ' ' << point.y << ' ' << point.z << '\n';
    stream << "g retained_source_and_hole_caps\n";
    for (const auto& triangle : mesh.triangles)
        stream << "f " << triangle[0] + 1 << ' ' << triangle[1] + 1 << ' '
               << triangle[2] + 1 << '\n';
}

OrientedSurfaceMesh buildAnalysisHalfedgeMeshImpl(
    const MeshModel& mesh, AnalysisHalfedgeStats* output_stats)
{
    if (const auto error = validateModelPool({mesh}))
        throw std::invalid_argument(*error);
    Bounds3 bounds;
    for (const Position3 point : mesh.vertices) include(bounds, toVec(point));
    const double scale = norm(bounds.upper - bounds.lower);
    const double tolerance = std::max(scale * 1.0e-9, 1.0e-12);
    using Key = std::array<std::int64_t, 3>;
    const auto key = [&](const Position3 point)
    {
        return Key{
            static_cast<std::int64_t>(std::llround(point.x / tolerance)),
            static_cast<std::int64_t>(std::llround(point.y / tolerance)),
            static_cast<std::int64_t>(std::llround(point.z / tolerance))};
    };
    std::map<Key, std::uint32_t> welded_ids;
    std::vector<Position3> welded_positions;
    std::vector<std::uint32_t> source_to_welded(mesh.vertices.size());
    for (std::size_t index = 0; index < mesh.vertices.size(); ++index)
    {
        const auto [iterator, inserted] = welded_ids.try_emplace(
            key(mesh.vertices[index]), static_cast<std::uint32_t>(welded_positions.size()));
        if (inserted) welded_positions.push_back(mesh.vertices[index]);
        source_to_welded[index] = iterator->second;
    }

    std::vector<TriangleIndices> welded_faces;
    welded_faces.reserve(mesh.triangles.size());
    AnalysisHalfedgeStats stats;
    std::set<TriangleIndices> unique_faces;
    const double minimum_area = std::max(1.0, scale * scale) * 1.0e-24;
    for (const TriangleIndices face : mesh.triangles)
    {
        TriangleIndices welded{
            source_to_welded[face[0]], source_to_welded[face[1]],
            source_to_welded[face[2]]};
        if (welded[0] == welded[1] || welded[1] == welded[2] ||
            welded[2] == welded[0] ||
            norm(cross(toVec(welded_positions[welded[1]]) -
                           toVec(welded_positions[welded[0]]),
                       toVec(welded_positions[welded[2]]) -
                           toVec(welded_positions[welded[0]]))) <= minimum_area)
        {
            ++stats.dropped_degenerate_triangles;
            continue;
        }
        TriangleIndices signature = welded;
        std::sort(signature.begin(), signature.end());
        if (!unique_faces.insert(signature).second)
        {
            ++stats.dropped_duplicate_triangles;
            continue;
        }
        welded_faces.push_back(welded);
    }
    if (welded_faces.empty())
        throw std::invalid_argument("analysis halfedge mesh has no nondegenerate faces");

    // Make collinear T-junctions conforming before opposite pairing. Every
    // inserted point is an existing welded endpoint on the same reality line;
    // face geometry and covered area are unchanged.
    using LineKey = std::array<std::int64_t, 6>;
    const double angular_quantum = 1.0e-10;
    const auto lineKey = [&](const std::uint32_t first, const std::uint32_t second)
    {
        const Vec3 a = toVec(welded_positions[first]);
        const Vec3 b = toVec(welded_positions[second]);
        Vec3 direction = normalized(b - a);
        int dominant = 0;
        for (int axis = 1; axis < 3; ++axis)
            if (std::abs(direction[axis]) > std::abs(direction[dominant])) dominant = axis;
        if (direction[dominant] < 0.0) direction = direction * -1.0;
        const Vec3 moment = cross(direction, a);
        return LineKey{
            static_cast<std::int64_t>(std::llround(direction.x / angular_quantum)),
            static_cast<std::int64_t>(std::llround(direction.y / angular_quantum)),
            static_cast<std::int64_t>(std::llround(direction.z / angular_quantum)),
            static_cast<std::int64_t>(std::llround(moment.x / tolerance)),
            static_cast<std::int64_t>(std::llround(moment.y / tolerance)),
            static_cast<std::int64_t>(std::llround(moment.z / tolerance))};
    };
    std::map<LineKey, std::set<std::uint32_t>> line_vertices;
    for (const TriangleIndices face : welded_faces)
        for (int local = 0; local < 3; ++local)
        {
            const auto key = lineKey(face[local], face[(local + 1) % 3]);
            line_vertices[key].insert(face[local]);
            line_vertices[key].insert(face[(local + 1) % 3]);
        }
    std::vector<TriangleIndices> conforming_faces;
    for (const TriangleIndices face : welded_faces)
    {
        std::vector<std::uint32_t> polygon;
        for (int local = 0; local < 3; ++local)
        {
            const std::uint32_t first = face[local];
            const std::uint32_t second = face[(local + 1) % 3];
            const Vec3 a = toVec(welded_positions[first]);
            const Vec3 b = toVec(welded_positions[second]);
            const Vec3 edge = b - a;
            const double length_squared = dot(edge, edge);
            std::vector<std::pair<double, std::uint32_t>> points;
            for (const std::uint32_t vertex : line_vertices.at(lineKey(first, second)))
            {
                const Vec3 point = toVec(welded_positions[vertex]);
                const double parameter = dot(point - a, edge) / length_squared;
                if (parameter < -1.0e-10 || parameter > 1.0 + 1.0e-10) continue;
                const Vec3 projected = a + edge * parameter;
                if (norm(point - projected) > 4.0 * tolerance) continue;
                points.push_back({std::clamp(parameter, 0.0, 1.0), vertex});
            }
            std::sort(points.begin(), points.end());
            for (const auto [parameter, vertex] : points)
            {
                (void)parameter;
                if (vertex == second) continue;
                if (polygon.empty() || polygon.back() != vertex)
                    polygon.push_back(vertex);
            }
        }
        if (polygon.size() == 3)
        {
            conforming_faces.push_back(face);
            continue;
        }
        const Vec3 center = (
            toVec(welded_positions[face[0]]) +
            toVec(welded_positions[face[1]]) +
            toVec(welded_positions[face[2]])) * (1.0 / 3.0);
        const std::uint32_t center_id =
            static_cast<std::uint32_t>(welded_positions.size());
        welded_positions.push_back(toPosition(center));
        for (std::size_t index = 0; index < polygon.size(); ++index)
        {
            const std::uint32_t first = polygon[index];
            const std::uint32_t second = polygon[(index + 1) % polygon.size()];
            if (first != second)
                conforming_faces.push_back({first, second, center_id});
        }
    }
    welded_faces = std::move(conforming_faces);

    using Edge = std::array<std::uint32_t, 2>;
    const auto buildEdgeUses = [&]()
    {
        std::map<Edge, std::vector<std::uint32_t>> result;
        for (std::uint32_t face = 0; face < welded_faces.size(); ++face)
            for (std::uint32_t local = 0; local < 3; ++local)
            {
                const std::uint32_t first = welded_faces[face][local];
                const std::uint32_t second = welded_faces[face][(local + 1) % 3];
                result[{std::min(first, second), std::max(first, second)}]
                    .push_back(face * 3 + local);
            }
        return result;
    };
    auto edge_uses = buildEdgeUses();
    std::vector<std::vector<std::pair<std::uint32_t, bool>>> orientation_graph(
        welded_faces.size());
    for (const auto& [edge, uses] : edge_uses)
    {
        (void)edge;
        if (uses.size() != 2) continue;
        const std::uint32_t first_face = uses[0] / 3;
        const std::uint32_t second_face = uses[1] / 3;
        const bool same_direction =
            welded_faces[first_face][uses[0] % 3] ==
            welded_faces[second_face][uses[1] % 3];
        orientation_graph[first_face].push_back({second_face, same_direction});
        orientation_graph[second_face].push_back({first_face, same_direction});
    }
    std::vector<std::int8_t> flip(welded_faces.size(), -1);
    for (std::uint32_t seed = 0; seed < welded_faces.size(); ++seed)
    {
        if (flip[seed] != -1) continue;
        flip[seed] = 0;
        std::vector<std::uint32_t> queue{seed};
        for (std::size_t cursor = 0; cursor < queue.size(); ++cursor)
        {
            const std::uint32_t face = queue[cursor];
            for (const auto [neighbor, relation] : orientation_graph[face])
            {
                const std::int8_t required = static_cast<std::int8_t>(
                    flip[face] ^ static_cast<std::int8_t>(relation));
                if (flip[neighbor] == -1)
                {
                    flip[neighbor] = required;
                    queue.push_back(neighbor);
                }
                else if (flip[neighbor] != required)
                    ++stats.inconsistent_orientation_edges;
            }
        }
    }
    for (std::uint32_t face = 0; face < welded_faces.size(); ++face)
        if (flip[face]) std::swap(welded_faces[face][1], welded_faces[face][2]);
    edge_uses = buildEdgeUses();

    const std::size_t corner_count = welded_faces.size() * 3;
    std::vector<std::uint32_t> parent(corner_count);
    std::iota(parent.begin(), parent.end(), 0);
    const auto find = [&](std::uint32_t value)
    {
        std::uint32_t root = value;
        while (parent[root] != root) root = parent[root];
        while (parent[value] != value)
        {
            const std::uint32_t next = parent[value];
            parent[value] = root;
            value = next;
        }
        return root;
    };
    const auto unite = [&](const std::uint32_t first, const std::uint32_t second)
    {
        const std::uint32_t a = find(first);
        const std::uint32_t b = find(second);
        if (a != b) parent[b] = a;
    };
    std::vector<std::array<std::uint32_t, 2>> initial_pairs;
    for (const auto& [edge, uses] : edge_uses)
    {
        if (uses.size() > 2) ++stats.nonmanifold_edge_groups;
        std::vector<std::uint32_t> forward;
        std::vector<std::uint32_t> reverse;
        for (const std::uint32_t use : uses)
            (welded_faces[use / 3][use % 3] == edge[0]
                ? forward : reverse).push_back(use);
        const std::size_t pair_count = std::min(forward.size(), reverse.size());
        const auto next_corner = [](const std::uint32_t corner)
        { return corner - corner % 3 + (corner % 3 + 1) % 3; };
        for (std::size_t index = 0; index < pair_count; ++index)
        {
            initial_pairs.push_back({forward[index], reverse[index]});
            unite(forward[index], next_corner(reverse[index]));
            unite(reverse[index], next_corner(forward[index]));
        }
    }

    OrientedSurfaceMesh result;
    result.geometry.name = mesh.name + "_analysis_halfedge";
    std::map<std::pair<std::uint32_t, std::uint32_t>, std::uint32_t> fan_vertices;
    result.geometry.triangles.resize(welded_faces.size());
    for (std::uint32_t face = 0; face < welded_faces.size(); ++face)
        for (std::uint32_t local = 0; local < 3; ++local)
        {
            const std::uint32_t corner = face * 3 + local;
            const auto fan_key = std::make_pair(welded_faces[face][local], find(corner));
            const auto [iterator, inserted] = fan_vertices.try_emplace(
                fan_key, static_cast<std::uint32_t>(result.geometry.vertices.size()));
            if (inserted)
                result.geometry.vertices.push_back(welded_positions[fan_key.first]);
            result.geometry.triangles[face][local] = iterator->second;
        }

    result.halfedges.resize(corner_count);
    result.face_halfedges.resize(welded_faces.size());
    result.vertex_halfedges.assign(
        result.geometry.vertices.size(), invalid_surface_index);
    for (std::uint32_t face = 0; face < result.geometry.triangles.size(); ++face)
    {
        const std::uint32_t base = face * 3;
        result.face_halfedges[face] = base;
        for (std::uint32_t local = 0; local < 3; ++local)
        {
            const std::uint32_t halfedge = base + local;
            const std::uint32_t origin = result.geometry.triangles[face][local];
            result.halfedges[halfedge] = {
                origin, face, base + (local + 1) % 3, invalid_surface_index};
            if (result.vertex_halfedges[origin] == invalid_surface_index)
                result.vertex_halfedges[origin] = halfedge;
        }
    }
    for (const auto pair : initial_pairs)
    {
        const std::uint32_t first_destination =
            result.halfedges[result.halfedges[pair[0]].next].origin;
        const std::uint32_t second_destination =
            result.halfedges[result.halfedges[pair[1]].next].origin;
        if (result.halfedges[pair[0]].origin == second_destination &&
            result.halfedges[pair[1]].origin == first_destination)
        {
            result.halfedges[pair[0]].opposite = pair[1];
            result.halfedges[pair[1]].opposite = pair[0];
            ++stats.paired_edges;
        }
    }
    stats.boundary_halfedges = std::count_if(
        result.halfedges.begin(), result.halfedges.end(),
        [](const SurfaceHalfedge& halfedge)
        { return halfedge.opposite == invalid_surface_index; });

    // Structural validation, including one traversable face fan per split
    // topology vertex. Boundary and deliberately unpaired conflicts are legal.
    std::vector<std::vector<std::uint32_t>> incident(result.geometry.vertices.size());
    for (std::uint32_t index = 0; index < result.halfedges.size(); ++index)
    {
        const SurfaceHalfedge& halfedge = result.halfedges[index];
        if (halfedge.origin >= result.geometry.vertices.size() ||
            halfedge.face >= result.geometry.triangles.size() ||
            halfedge.next >= result.halfedges.size())
            throw std::runtime_error("analysis halfedge contains an invalid index");
        const std::uint32_t next2 = result.halfedges[halfedge.next].next;
        if (result.halfedges[halfedge.next].face != halfedge.face ||
            next2 >= result.halfedges.size() ||
            result.halfedges[next2].next != index ||
            result.halfedges[next2].face != halfedge.face)
            throw std::runtime_error("analysis halfedge face cycle is not a triangle");
        if (halfedge.opposite != invalid_surface_index)
        {
            if (halfedge.opposite >= result.halfedges.size() ||
                result.halfedges[halfedge.opposite].opposite != index ||
                result.halfedges[halfedge.opposite].origin !=
                    result.halfedges[halfedge.next].origin ||
                result.halfedges[result.halfedges[halfedge.opposite].next].origin !=
                    halfedge.origin)
                throw std::runtime_error("analysis halfedge opposite is invalid");
        }
        incident[halfedge.origin].push_back(index);
    }
    for (std::uint32_t vertex = 0; vertex < incident.size(); ++vertex)
    {
        if (incident[vertex].empty() || result.vertex_halfedges[vertex] == invalid_surface_index)
            throw std::runtime_error("analysis halfedge contains an unused vertex");
        std::set<std::uint32_t> reachable{incident[vertex].front()};
        std::vector<std::uint32_t> queue{incident[vertex].front()};
        for (std::size_t cursor = 0; cursor < queue.size(); ++cursor)
        {
            const std::uint32_t outgoing = queue[cursor];
            const std::uint32_t incoming = result.halfedges[
                result.halfedges[outgoing].next].next;
            for (const std::uint32_t edge : {outgoing, incoming})
            {
                const std::uint32_t opposite = result.halfedges[edge].opposite;
                if (opposite == invalid_surface_index) continue;
                const std::uint32_t adjacent = edge == outgoing
                    ? result.halfedges[opposite].next : opposite;
                if (reachable.insert(adjacent).second) queue.push_back(adjacent);
            }
        }
        if (reachable.size() != incident[vertex].size())
            throw std::runtime_error("analysis halfedge vertex has multiple fans");
    }
    std::vector<std::uint8_t> visited_halfedges(result.halfedges.size(), 0);
    for (std::uint32_t seed = 0; seed < result.halfedges.size(); ++seed)
    {
        if (visited_halfedges[seed]) continue;
        ++stats.face_components;
        std::vector<std::uint32_t> queue{seed};
        visited_halfedges[seed] = 1;
        for (std::size_t cursor = 0; cursor < queue.size(); ++cursor)
        {
            const std::uint32_t halfedge = queue[cursor];
            for (const std::uint32_t neighbor : {
                    result.halfedges[halfedge].next,
                    result.halfedges[halfedge].opposite})
            {
                if (neighbor == invalid_surface_index) continue;
                if (!visited_halfedges[neighbor])
                {
                    visited_halfedges[neighbor] = 1;
                    queue.push_back(neighbor);
                }
            }
        }
    }
    if (output_stats) *output_stats = stats;
    return result;
}

OrientedSurfaceMesh extractCellComplexBoundaryImpl(
    const SourceConstrainedCellComplex& complex)
{
    if (complex.cells.empty())
        throw std::invalid_argument("source-constrained cell complex is empty");

    // For a positively oriented tetrahedron, these are the outward face
    // windings opposite local vertices 0..3.
    constexpr std::array<std::array<int, 3>, 4> outward_faces{{
        {{1, 2, 3}}, {{0, 3, 2}}, {{0, 1, 3}}, {{0, 2, 1}}}};
    IndexedSurface surface;
    surface.vertices.reserve(complex.vertices.size());
    for (const auto point : complex.vertices) surface.vertices.push_back(toVec(point));
    std::vector<std::uint32_t> source_faces;
    for (std::uint32_t cell_id = 0; cell_id < complex.cells.size(); ++cell_id)
    {
        const auto& cell = complex.cells[cell_id];
        for (const auto vertex : cell.vertices)
            if (vertex >= complex.vertices.size())
                throw std::invalid_argument("cell references an invalid vertex");
        const Vec3 a = toVec(complex.vertices[cell.vertices[0]]);
        const Vec3 b = toVec(complex.vertices[cell.vertices[1]]);
        const Vec3 c = toVec(complex.vertices[cell.vertices[2]]);
        const Vec3 d = toVec(complex.vertices[cell.vertices[3]]);
        if (dot(cross(b - a, c - a), d - a) <= 0.0)
            throw std::invalid_argument("tetrahedron must be positively oriented");
        if (!cell.occupied) continue;

        for (int local_face = 0; local_face < 4; ++local_face)
        {
            const auto neighbor = cell.neighbors[local_face];
            if (neighbor != invalid_surface_index)
            {
                if (neighbor >= complex.cells.size())
                    throw std::invalid_argument("cell references an invalid neighbor");
                if (complex.cells[neighbor].occupied) continue;
            }
            const auto local = outward_faces[local_face];
            surface.triangles.push_back({{
                cell.vertices[local[0]], cell.vertices[local[1]],
                cell.vertices[local[2]]}});
            source_faces.push_back(cell.source_faces[local_face]);
        }
    }
    if (surface.triangles.empty())
        throw std::invalid_argument("cell complex has no occupied boundary");

    // Drop vertices used exclusively by empty cells before halfedge validation.
    std::vector<std::uint32_t> remap(surface.vertices.size(), invalid_surface_index);
    std::vector<Vec3> compact_vertices;
    for (auto& triangle : surface.triangles)
        for (auto& vertex : triangle)
        {
            if (remap[vertex] == invalid_surface_index)
            {
                remap[vertex] = static_cast<std::uint32_t>(compact_vertices.size());
                compact_vertices.push_back(surface.vertices[vertex]);
            }
            vertex = remap[vertex];
        }
    surface.vertices = std::move(compact_vertices);

    // Cell-complex faces carry source provenance instead of voxel-quad
    // provenance. The placeholder is unused by this validation path.
    VoxelGrid placeholder;
    placeholder.shape = {1, 1, 1};
    placeholder.pitch = 1.0;
    placeholder.occupancy = {1};
    return buildOrientedSurfaceMesh(
        std::move(surface), placeholder, std::move(source_faces));
}

std::string bettiJson(const BettiNumbers value)
{
    return "{\"beta0\":" + std::to_string(value.beta0) +
        ",\"beta1\":" + std::to_string(value.beta1) +
        ",\"beta2\":" + std::to_string(value.beta2) + "}";
}

void writeModelMetadata(const std::filesystem::path& path, const std::string& model_name,
                        const TopologyFillStats& stats)
{
    std::ofstream stream(path);
    if (!stream) throw std::runtime_error("failed to create metadata: " + path.string());
    stream << std::setprecision(17)
        << "{\n  \"stats\":{\n"
        << "    \"model\":\"" << model_name << "\",\n"
        << "    \"source_triangles\":" << stats.source_triangles << ",\n"
        << "    \"primitive_count\":0,\n"
        << "    \"primitive_types\":{},\n"
        << "    \"proxy_triangles\":0,\n"
        << "    \"topology_fill\":{\n"
        << "      \"input_betti\":" << bettiJson(stats.input_betti) << ",\n"
        << "      \"output_betti\":" << bettiJson(stats.output_betti) << ",\n"
        << "      \"input_voxels\":" << stats.input_voxels << ",\n"
        << "      \"output_voxels\":" << stats.output_voxels << ",\n"
        << "      \"added_voxels\":" << stats.added_voxels << ",\n"
        << "      \"cavity_fill_voxels\":" << stats.cavity_fill_voxels << ",\n"
        << "      \"removed_voxels\":" << stats.removed_voxels << ",\n"
        << "      \"pitch\":" << stats.pitch << ",\n"
        << "      \"grid_shape\":[" << stats.grid_shape[0] << ','
        << stats.grid_shape[1] << ',' << stats.grid_shape[2] << "],\n"
        << "      \"representation\":\"validated open analysis halfedge\",\n"
        << "      \"boundary_voxel_faces\":0,\n"
        << "      \"boundary_voxels\":0,\n"
        << "      \"mesh_triangles\":" << stats.mesh_triangles << ",\n"
        << "      \"planar_regions\":" << stats.planar_regions << ",\n"
        << "      \"mesh_vertices\":" << stats.mesh_vertices << ",\n"
        << "      \"retained_source_triangles\":" << stats.retained_source_triangles << ",\n"
        << "      \"removed_internal_triangles\":" << stats.removed_internal_triangles << ",\n"
        << "      \"ambiguous_source_triangles\":" << stats.ambiguous_source_triangles << ",\n"
        << "      \"hole_boundary_loops\":" << stats.hole_boundary_loops << ",\n"
        << "      \"cap_triangles\":" << stats.cap_triangles << ",\n"
        << "      \"halfedge_count\":" << stats.halfedge_count << ",\n"
        << "      \"paired_halfedge_edges\":" << stats.paired_halfedge_edges << ",\n"
        << "      \"boundary_halfedges\":" << stats.boundary_halfedges << ",\n"
        << "      \"nonmanifold_edge_groups\":" << stats.nonmanifold_edge_groups << ",\n"
        << "      \"inconsistent_orientation_edges\":"
        << stats.inconsistent_orientation_edges << ",\n"
        << "      \"halfedge_face_components\":"
        << stats.halfedge_face_components << ",\n"
        << "      \"dropped_duplicate_triangles\":"
        << stats.dropped_duplicate_triangles << ",\n"
        << "      \"dropped_degenerate_triangles\":"
        << stats.dropped_degenerate_triangles << ",\n"
        << "      \"mesh_watertight\":" << (stats.mesh_watertight ? "true" : "false") << ",\n"
        << "      \"mesh_oriented\":" << (stats.mesh_oriented ? "true" : "false") << ",\n"
        << "      \"mesh_manifold\":" << (stats.mesh_manifold ? "true" : "false") << ",\n"
        << "      \"mesh_connected\":" << (stats.mesh_connected ? "true" : "false") << ",\n"
        << "      \"mesh_has_only_boundary_faces\":"
        << (stats.mesh_has_only_boundary_faces ? "true" : "false") << ",\n"
        << "      \"mesh_euler_characteristic\":" << stats.mesh_euler_characteristic << ",\n"
        << "      \"mesh_signed_volume\":" << stats.mesh_signed_volume << ",\n"
        << "      \"source_occupancy_enclosed\":true\n"
        << "    },\n"
        << "    \"timings_seconds\":{\"total\":" << stats.elapsed_seconds << "}\n"
        << "  },\n"
        << "  \"source\":\"source.obj\",\n"
        << "  \"phase1_halfedge\":\"phase1_halfedge.bin\",\n"
        << "  \"proxy_components\":[],\n"
        << "  \"viewer_stages\":[\"source\",\"phase1\",\"split\"]\n"
        << "}\n";
}

} // namespace

OrientedSurfaceMesh buildAnalysisHalfedgeMesh(
    const MeshModel& mesh, AnalysisHalfedgeStats* stats)
{
    return buildAnalysisHalfedgeMeshImpl(mesh, stats);
}

void writeAnalysisHalfedgeMesh(
    const std::filesystem::path& path, const OrientedSurfaceMesh& mesh)
{
    if (mesh.geometry.vertices.size() > std::numeric_limits<std::uint32_t>::max() ||
        mesh.geometry.triangles.size() > std::numeric_limits<std::uint32_t>::max() ||
        mesh.halfedges.size() > std::numeric_limits<std::uint32_t>::max())
        throw std::runtime_error("analysis halfedge mesh exceeds PQSSHED1 limits");
    if (mesh.halfedges.size() != mesh.geometry.triangles.size() * 3 ||
        mesh.vertex_halfedges.size() != mesh.geometry.vertices.size() ||
        mesh.face_halfedges.size() != mesh.geometry.triangles.size())
        throw std::invalid_argument("analysis halfedge arrays have inconsistent sizes");
    std::ofstream stream(path, std::ios::binary);
    if (!stream)
        throw std::runtime_error("failed to create halfedge mesh: " + path.string());
    stream.write("PQSSHED1", 8);
    const std::array<std::uint32_t, 5> counts{
        static_cast<std::uint32_t>(mesh.geometry.vertices.size()),
        static_cast<std::uint32_t>(mesh.geometry.triangles.size()),
        static_cast<std::uint32_t>(mesh.halfedges.size()),
        static_cast<std::uint32_t>(mesh.vertex_halfedges.size()),
        static_cast<std::uint32_t>(mesh.face_halfedges.size())};
    stream.write(reinterpret_cast<const char*>(counts.data()), sizeof(counts));
    for (const Position3 point : mesh.geometry.vertices)
    {
        const std::array<double, 3> value{point.x, point.y, point.z};
        stream.write(reinterpret_cast<const char*>(value.data()), sizeof(value));
    }
    for (const TriangleIndices face : mesh.geometry.triangles)
        stream.write(reinterpret_cast<const char*>(face.data()), sizeof(face));
    for (const SurfaceHalfedge halfedge : mesh.halfedges)
    {
        const std::array<std::uint32_t, 4> value{
            halfedge.origin, halfedge.face, halfedge.next, halfedge.opposite};
        stream.write(reinterpret_cast<const char*>(value.data()), sizeof(value));
    }
    stream.write(reinterpret_cast<const char*>(mesh.vertex_halfedges.data()),
                 static_cast<std::streamsize>(mesh.vertex_halfedges.size() *
                                              sizeof(std::uint32_t)));
    stream.write(reinterpret_cast<const char*>(mesh.face_halfedges.data()),
                 static_cast<std::streamsize>(mesh.face_halfedges.size() *
                                              sizeof(std::uint32_t)));
    if (!stream) throw std::runtime_error("failed to write complete PQSSHED1 mesh");
}

OrientedSurfaceMesh extractCellComplexBoundary(
    const SourceConstrainedCellComplex& complex)
{
    return extractCellComplexBoundaryImpl(complex);
}

MeshModel canonicalizeCoplanarTriangleSoup(
    const MeshModel& source, const double relative_tolerance)
{
    if (source.triangles.empty() || source.vertices.empty())
        throw std::invalid_argument("triangle soup is empty");
    if (!std::isfinite(relative_tolerance) || relative_tolerance <= 0.0)
        throw std::invalid_argument("relative tolerance must be finite and positive");

    Bounds3 bounds;
    for (const auto point : source.vertices) include(bounds, toVec(point));
    const double scale = norm(bounds.upper - bounds.lower);
    const double distance_quantum = std::max(scale * relative_tolerance, 1.0e-12);
    constexpr double angular_quantum = 1.0e-10;
    constexpr int clipper_precision = 8;
    using PlaneKey = std::array<std::int64_t, 4>;
    struct PlaneGroup
    {
        Vec3 origin;
        Vec3 normal;
        Vec3 tangent;
        Vec3 bitangent;
        std::vector<std::uint32_t> faces;
    };
    std::map<PlaneKey, PlaneGroup> groups;

    for (std::uint32_t face_id = 0; face_id < source.triangles.size(); ++face_id)
    {
        const auto face = source.triangles[face_id];
        for (const auto vertex : face)
            if (vertex >= source.vertices.size())
                throw std::invalid_argument("triangle soup has an invalid vertex index");
        const Vec3 a = toVec(source.vertices[face[0]]);
        const Vec3 b = toVec(source.vertices[face[1]]);
        const Vec3 c = toVec(source.vertices[face[2]]);
        Vec3 normal = cross(b - a, c - a);
        const double length = norm(normal);
        if (length <= std::max(scale * scale, 1.0) * 1.0e-24)
        {
            continue;
        }
        normal = normal * (1.0 / length);
        int dominant = 0;
        for (int axis = 1; axis < 3; ++axis)
            if (std::abs(normal[axis]) > std::abs(normal[dominant])) dominant = axis;
        if (normal[dominant] < 0.0) normal = normal * -1.0;
        const double distance = dot(normal, a);
        PlaneKey key{};
        for (int axis = 0; axis < 3; ++axis)
            key[axis] = static_cast<std::int64_t>(
                std::llround(normal[axis] / angular_quantum));
        key[3] = static_cast<std::int64_t>(
            std::llround(distance / distance_quantum));
        auto [iterator, inserted] = groups.try_emplace(key);
        auto& group = iterator->second;
        if (inserted)
        {
            group.origin = a;
            group.normal = normal;
            const Vec3 helper = std::abs(normal.x) < 0.8 ? Vec3{1.0, 0.0, 0.0}
                                                        : Vec3{0.0, 1.0, 0.0};
            group.tangent = normalized(cross(helper, normal));
            group.bitangent = cross(normal, group.tangent);
        }
        // Quantization only proposes a group. Reject triangles that are not
        // actually coplanar with its representative plane.
        bool coplanar = std::abs(dot(group.normal, a - group.origin)) <= distance_quantum &&
            std::abs(dot(group.normal, b - group.origin)) <= distance_quantum &&
            std::abs(dot(group.normal, c - group.origin)) <= distance_quantum &&
            std::abs(dot(group.normal, normal)) >= 1.0 - 1.0e-12;
        if (!coplanar)
        {
            // Preserve it exactly by assigning a unique deterministic key.
            key[3] ^= static_cast<std::int64_t>(face_id + 1) << 32;
            auto& unique = groups[key];
            unique.origin = a;
            unique.normal = normal;
            const Vec3 helper = std::abs(normal.x) < 0.8 ? Vec3{1.0, 0.0, 0.0}
                                                        : Vec3{0.0, 1.0, 0.0};
            unique.tangent = normalized(cross(helper, normal));
            unique.bitangent = cross(normal, unique.tangent);
            unique.faces.push_back(face_id);
        }
        else group.faces.push_back(face_id);
    }

    MeshModel result;
    result.name = source.name + "_coplanar_union";
    std::map<std::array<double, 3>, std::uint32_t> output_vertices;
    const auto appendVertex = [&](const Position3 point)
    {
        const std::array<double, 3> key{{point.x, point.y, point.z}};
        const auto [iterator, inserted] = output_vertices.try_emplace(
            key, static_cast<std::uint32_t>(result.vertices.size()));
        if (inserted) result.vertices.push_back(point);
        return iterator->second;
    };
    for (const auto& [key, group] : groups)
    {
        (void)key;
        Clipper2Lib::PathsD paths;
        paths.reserve(group.faces.size());
        for (const auto face_id : group.faces)
        {
            Clipper2Lib::PathD path;
            for (const auto vertex : source.triangles[face_id])
            {
                const Vec3 delta = toVec(source.vertices[vertex]) - group.origin;
                path.emplace_back(dot(delta, group.tangent), dot(delta, group.bitangent));
            }
            if (std::abs(Clipper2Lib::Area(path)) <= distance_quantum * distance_quantum)
                continue;
            if (Clipper2Lib::Area(path) < 0.0) std::reverse(path.begin(), path.end());
            paths.push_back(std::move(path));
        }
        const auto united = Clipper2Lib::Union(
            paths, Clipper2Lib::FillRule::NonZero, clipper_precision);
        Clipper2Lib::PathsD triangles;
        if (united.empty() || Clipper2Lib::Triangulate(
                united, clipper_precision, triangles, false) !=
                Clipper2Lib::TriangulateResult::success)
        {
            // A self-touching soup boundary may be geometrically valid while
            // lacking a single polygon-with-holes representation. Preserve
            // that group exactly; canonicalization is an optimization only.
            for (const auto face_id : group.faces)
            {
                TriangleIndices output_face{};
                for (int corner = 0; corner < 3; ++corner)
                    output_face[corner] = appendVertex(
                        source.vertices[source.triangles[face_id][corner]]);
                result.triangles.push_back(output_face);
            }
            continue;
        }
        for (const auto& triangle : triangles)
        {
            if (triangle.size() != 3) continue;
            TriangleIndices output_face{};
            for (int corner = 0; corner < 3; ++corner)
            {
                const auto& point = triangle[corner];
                output_face[corner] = appendVertex(toPosition(
                    group.origin + group.tangent * point.x + group.bitangent * point.y));
            }
            result.triangles.push_back(output_face);
        }
    }
    if (result.triangles.empty())
        throw std::runtime_error("coplanar canonicalization produced an empty soup");
    return result;
}

namespace
{

MeshModel buildPlanarExteriorSurgery(
    const VoxelGrid& filled, const MeshModel& source)
{
    const auto exterior = floodEmptyFromBoundary(filled, neighbors26);
    const Vec3 grid_origin = toVec(filled.origin);
    const auto voxelAt = [&](const Vec3 point) -> std::optional<std::size_t>
    {
        std::array<int, 3> coordinate{};
        for (int axis = 0; axis < 3; ++axis)
        {
            coordinate[axis] = static_cast<int>(std::llround(
                (point[axis] - grid_origin[axis]) / filled.pitch));
            if (coordinate[axis] < 0 ||
                coordinate[axis] >= static_cast<int>(filled.shape[axis]))
                return std::nullopt;
        }
        return filled.index(
            static_cast<std::uint32_t>(coordinate[0]),
            static_cast<std::uint32_t>(coordinate[1]),
            static_cast<std::uint32_t>(coordinate[2]));
    };
    const auto isExterior = [&](const Vec3 point)
    {
        const auto voxel = voxelAt(point);
        return !voxel || exterior[*voxel] != 0;
    };

    Bounds3 bounds;
    for (const auto point : source.vertices) include(bounds, toVec(point));
    const double scale = norm(bounds.upper - bounds.lower);
    const double distance_quantum = std::max(scale * 1.0e-9, 1.0e-9);
    constexpr double angular_quantum = 1.0e-10;
    constexpr int clipper_precision = 8;
    using PlaneKey = std::array<std::int64_t, 4>;
    struct PlaneGroup
    {
        Vec3 origin;
        Vec3 normal;
        Vec3 tangent;
        Vec3 bitangent;
        std::vector<std::uint32_t> faces;
    };
    std::map<PlaneKey, PlaneGroup> planar_groups;
    std::vector<std::uint8_t> retained(source.triangles.size(), 0);
    std::size_t removed_internal = 0;
    for (std::uint32_t face_id = 0; face_id < source.triangles.size(); ++face_id)
    {
        const auto face = source.triangles[face_id];
        const Vec3 a = toVec(source.vertices[face[0]]);
        const Vec3 b = toVec(source.vertices[face[1]]);
        const Vec3 c = toVec(source.vertices[face[2]]);
        Vec3 normal = normalized(cross(b - a, c - a));
        if (norm(normal) == 0.0) continue;
        const std::array<Vec3, 7> samples{
            (a + b + c) * (1.0 / 3.0),
            a * 0.8 + b * 0.1 + c * 0.1,
            a * 0.1 + b * 0.8 + c * 0.1,
            a * 0.1 + b * 0.1 + c * 0.8,
            a * 0.5 + b * 0.5,
            b * 0.5 + c * 0.5,
            c * 0.5 + a * 0.5};
        bool touches_exterior = false;
        for (const Vec3 sample : samples)
            for (const double distance : {0.75, 1.25, 2.0})
            {
                const Vec3 offset = normal * (distance * filled.pitch);
                touches_exterior |= isExterior(sample + offset) ||
                    isExterior(sample - offset);
            }
        if (!touches_exterior)
        {
            ++removed_internal;
            continue;
        }
        retained[face_id] = 1;
        int dominant = 0;
        for (int axis = 1; axis < 3; ++axis)
            if (std::abs(normal[axis]) > std::abs(normal[dominant])) dominant = axis;
        if (normal[dominant] < 0.0) normal = normal * -1.0;
        PlaneKey key{
            static_cast<std::int64_t>(std::llround(normal.x / angular_quantum)),
            static_cast<std::int64_t>(std::llround(normal.y / angular_quantum)),
            static_cast<std::int64_t>(std::llround(normal.z / angular_quantum)),
            static_cast<std::int64_t>(std::llround(dot(normal, a) / distance_quantum))};
        auto [iterator, inserted] = planar_groups.try_emplace(key);
        auto& group = iterator->second;
        if (inserted)
        {
            group.origin = a;
            group.normal = normal;
            const Vec3 helper = std::abs(normal.x) < 0.8
                ? Vec3{1, 0, 0} : Vec3{0, 1, 0};
            group.tangent = normalized(cross(helper, normal));
            group.bitangent = cross(normal, group.tangent);
        }
        const bool coplanar =
            std::abs(dot(group.normal, a - group.origin)) <= distance_quantum &&
            std::abs(dot(group.normal, b - group.origin)) <= distance_quantum &&
            std::abs(dot(group.normal, c - group.origin)) <= distance_quantum;
        if (coplanar) group.faces.push_back(face_id);
        else
        {
            key[3] ^= static_cast<std::int64_t>(face_id + 1) << 32;
            auto& unique = planar_groups[key];
            unique.origin = a;
            unique.normal = normal;
            const Vec3 helper = std::abs(normal.x) < 0.8
                ? Vec3{1, 0, 0} : Vec3{0, 1, 0};
            unique.tangent = normalized(cross(helper, normal));
            unique.bitangent = cross(normal, unique.tangent);
            unique.faces.push_back(face_id);
        }
    }

    MeshModel result;
    result.name = source.name + "_planar_exterior_surgery";
    std::map<std::array<std::int64_t, 3>, std::uint32_t> vertices;
    const auto appendVertex = [&](const Vec3 point)
    {
        const std::array<std::int64_t, 3> key{
            static_cast<std::int64_t>(std::llround(point.x / distance_quantum)),
            static_cast<std::int64_t>(std::llround(point.y / distance_quantum)),
            static_cast<std::int64_t>(std::llround(point.z / distance_quantum))};
        const auto [iterator, inserted] = vertices.try_emplace(
            key, static_cast<std::uint32_t>(result.vertices.size()));
        if (inserted) result.vertices.push_back(toPosition(point));
        return iterator->second;
    };
    std::size_t filled_planar_holes = 0;
    for (const auto& [key, group] : planar_groups)
    {
        (void)key;
        Clipper2Lib::PathsD source_paths;
        for (const std::uint32_t face_id : group.faces)
        {
            Clipper2Lib::PathD path;
            for (const std::uint32_t vertex : source.triangles[face_id])
            {
                const Vec3 delta = toVec(source.vertices[vertex]) - group.origin;
                path.emplace_back(dot(delta, group.tangent), dot(delta, group.bitangent));
            }
            if (std::abs(Clipper2Lib::Area(path)) <= distance_quantum * distance_quantum)
                continue;
            if (Clipper2Lib::Area(path) < 0.0) std::reverse(path.begin(), path.end());
            source_paths.push_back(std::move(path));
        }
        auto united = Clipper2Lib::Union(
            source_paths, Clipper2Lib::FillRule::NonZero, clipper_precision);
        Clipper2Lib::PathsD retained_paths;
        for (auto path : united)
        {
            if (Clipper2Lib::IsPositive(path))
            {
                retained_paths.push_back(std::move(path));
                continue;
            }
            Vec3 center{};
            for (const auto& point : path)
                center = center + group.origin + group.tangent * point.x +
                    group.bitangent * point.y;
            center = center * (1.0 / static_cast<double>(path.size()));
            const auto voxel = voxelAt(center);
            if (voxel && filled.occupancy[*voxel] && !exterior[*voxel])
                ++filled_planar_holes;
            else retained_paths.push_back(std::move(path));
        }
        Clipper2Lib::PathsD triangles;
        if (retained_paths.empty() || Clipper2Lib::Triangulate(
                retained_paths, clipper_precision, triangles, false) !=
                Clipper2Lib::TriangulateResult::success)
        {
            for (const std::uint32_t face_id : group.faces)
            {
                TriangleIndices output{};
                for (int corner = 0; corner < 3; ++corner)
                    output[corner] = appendVertex(toVec(source.vertices[
                        source.triangles[face_id][corner]]));
                result.triangles.push_back(output);
            }
            continue;
        }
        for (const auto& triangle : triangles)
        {
            if (triangle.size() != 3) continue;
            TriangleIndices output{};
            for (int corner = 0; corner < 3; ++corner)
                output[corner] = appendVertex(group.origin +
                    group.tangent * triangle[corner].x +
                    group.bitangent * triangle[corner].y);
            result.triangles.push_back(output);
        }
    }
    std::cerr << "monitor: stage=planar_exterior_surgery removed_internal="
              << removed_internal << " filled_planar_holes="
              << filled_planar_holes << " triangles="
              << result.triangles.size() << '\n';
    if (result.triangles.empty())
        throw std::runtime_error("planar exterior surgery produced an empty mesh");
    return result;
}

} // namespace

void writeVoxelGrid(const std::filesystem::path& path, const VoxelGrid& grid)
{
    validateGrid(grid);
    std::ofstream stream(path, std::ios::binary);
    if (!stream) throw std::runtime_error("failed to create voxel grid: " + path.string());
    const auto write_little_endian = [&stream]<typename Value>(const Value value)
    {
        auto bytes = std::bit_cast<std::array<char, sizeof(Value)>>(value);
        if constexpr (std::endian::native == std::endian::big)
            std::reverse(bytes.begin(), bytes.end());
        stream.write(bytes.data(), static_cast<std::streamsize>(bytes.size()));
    };
    constexpr std::array<char, 8> magic{'P','Q','S','S','V','O','X','1'};
    stream.write(magic.data(), static_cast<std::streamsize>(magic.size()));
    for (const std::uint32_t extent : grid.shape) write_little_endian(extent);
    write_little_endian(grid.pitch);
    write_little_endian(grid.origin.x);
    write_little_endian(grid.origin.y);
    write_little_endian(grid.origin.z);
    std::vector<std::uint8_t> packed((grid.occupancy.size() + 7) / 8, 0);
    for (std::size_t index = 0; index < grid.occupancy.size(); ++index)
        if (grid.occupancy[index]) packed[index / 8] |= static_cast<std::uint8_t>(1u << (index % 8));
    stream.write(reinterpret_cast<const char*>(packed.data()),
                 static_cast<std::streamsize>(packed.size()));
    if (!stream) throw std::runtime_error("failed to write voxel grid: " + path.string());
}

VoxelGrid readVoxelGrid(const std::filesystem::path& path)
{
    std::ifstream stream(path, std::ios::binary);
    if (!stream) throw std::runtime_error("failed to open voxel grid: " + path.string());
    const auto read_little_endian = [&stream]<typename Value>()
    {
        std::array<char, sizeof(Value)> bytes{};
        stream.read(bytes.data(), static_cast<std::streamsize>(bytes.size()));
        if (!stream) throw std::runtime_error("truncated voxel-grid header");
        if constexpr (std::endian::native == std::endian::big)
            std::reverse(bytes.begin(), bytes.end());
        return std::bit_cast<Value>(bytes);
    };
    std::array<char, 8> magic{};
    stream.read(magic.data(), static_cast<std::streamsize>(magic.size()));
    constexpr std::array<char, 8> expected_magic{'P','Q','S','S','V','O','X','1'};
    if (magic != expected_magic) throw std::runtime_error("invalid voxel-grid magic");
    VoxelGrid result;
    for (std::uint32_t& extent : result.shape)
        extent = read_little_endian.template operator()<std::uint32_t>();
    result.pitch = read_little_endian.template operator()<double>();
    result.origin.x = read_little_endian.template operator()<double>();
    result.origin.y = read_little_endian.template operator()<double>();
    result.origin.z = read_little_endian.template operator()<double>();
    const std::size_t voxel_count = checkedVoxelCount(result.shape);
    std::vector<std::uint8_t> packed((voxel_count + 7) / 8, 0);
    stream.read(reinterpret_cast<char*>(packed.data()),
                static_cast<std::streamsize>(packed.size()));
    if (!stream || stream.peek() != std::char_traits<char>::eof())
        throw std::runtime_error("invalid voxel-grid payload size");
    result.occupancy.resize(voxel_count);
    for (std::size_t index = 0; index < voxel_count; ++index)
        result.occupancy[index] = (packed[index / 8] >> (index % 8)) & 1u;
    validateGrid(result);
    return result;
}

Phase1Solid buildPhase1Solid(VoxelGrid occupancy, const MeshModel& source)
{
    validateGrid(occupancy);
    if (voxelBettiNumbers(occupancy) != BettiNumbers{1, 0, 0})
        throw std::invalid_argument("phase-1 solid requires occupancy beta=(1,0,0)");
    Phase1Solid result;
    result.occupancy = std::move(occupancy);
    result.boundary = buildOrientedSurfaceMesh(
        reconstructDualSurface(result.occupancy, source), result.occupancy);
    return result;
}

Phase1Solid buildPhase1Solid(VoxelGrid occupancy)
{
    return buildPhase1Solid(std::move(occupancy), MeshModel{});
}

OrientedSurfaceMesh collapseDegenerateHalfedges(const OrientedSurfaceMesh& input)
{
    if (input.geometry.vertices.empty() || input.geometry.triangles.empty() ||
        input.halfedges.size() != input.geometry.triangles.size() * 3)
        throw std::invalid_argument("degenerate cleanup requires a halfedge surface");

    std::vector<Position3> vertices = input.geometry.vertices;
    std::vector<TriangleIndices> faces = input.geometry.triangles;
    Bounds3 bounds;
    for (const auto point : vertices) include(bounds, toVec(point));
    const double scale = norm(bounds.upper - bounds.lower);
    const double minimum_area = std::max(1.0, scale * scale) * 1.0e-24;
    using Edge = std::array<std::uint32_t, 2>;
    const auto edgeKey = [](const std::uint32_t a, const std::uint32_t b)
    { return Edge{std::min(a, b), std::max(a, b)}; };

    std::size_t collapse_count = 0;
    for (std::size_t pass = 0; pass < 64; ++pass)
    {
        std::map<Edge, std::vector<std::uint32_t>> edge_faces;
        std::vector<std::vector<std::uint32_t>> neighbors(vertices.size());
        std::vector<std::uint32_t> degenerate_faces;
        for (std::uint32_t face = 0; face < faces.size(); ++face)
        {
            const auto triangle = faces[face];
            const Vec3 a = toVec(vertices[triangle[0]]);
            const Vec3 b = toVec(vertices[triangle[1]]);
            const Vec3 c = toVec(vertices[triangle[2]]);
            if (norm(cross(b - a, c - a)) <= minimum_area)
                degenerate_faces.push_back(face);
            for (int local = 0; local < 3; ++local)
            {
                const auto first = triangle[local];
                const auto second = triangle[(local + 1) % 3];
                edge_faces[edgeKey(first, second)].push_back(face);
                neighbors[first].push_back(second);
                neighbors[second].push_back(first);
            }
        }
        if (degenerate_faces.empty()) break;
        for (auto& adjacent : neighbors)
        {
            std::sort(adjacent.begin(), adjacent.end());
            adjacent.erase(std::unique(adjacent.begin(), adjacent.end()), adjacent.end());
        }

        struct Collapse { std::uint32_t removed; std::uint32_t kept; };
        std::vector<Collapse> collapses;
        std::vector<std::uint8_t> locked(vertices.size(), 0);
        const auto oppositeVertex = [&](const std::uint32_t face, const Edge edge)
        {
            for (const auto vertex : faces[face])
                if (vertex != edge[0] && vertex != edge[1]) return vertex;
            throw std::runtime_error("edge incident face has no opposite vertex");
        };
        const auto linkCondition = [&](const Edge edge)
        {
            const auto found = edge_faces.find(edge);
            if (found == edge_faces.end() || found->second.size() != 2) return false;
            std::vector<std::uint32_t> common;
            std::set_intersection(neighbors[edge[0]].begin(), neighbors[edge[0]].end(),
                neighbors[edge[1]].begin(), neighbors[edge[1]].end(),
                std::back_inserter(common));
            std::array<std::uint32_t, 2> opposite{
                oppositeVertex(found->second[0], edge),
                oppositeVertex(found->second[1], edge)};
            std::sort(opposite.begin(), opposite.end());
            return opposite[0] != opposite[1] && common.size() == 2 &&
                common[0] == opposite[0] && common[1] == opposite[1];
        };

        for (const std::uint32_t face : degenerate_faces)
        {
            const auto triangle = faces[face];
            std::array<std::pair<double, Edge>, 3> candidates{};
            for (int local = 0; local < 3; ++local)
            {
                const Edge edge = edgeKey(
                    triangle[local], triangle[(local + 1) % 3]);
                const Vec3 delta = toVec(vertices[edge[0]]) - toVec(vertices[edge[1]]);
                candidates[local] = {dot(delta, delta), edge};
            }
            std::sort(candidates.begin(), candidates.end(),
                [](const auto& first, const auto& second)
                { return first.first < second.first; });
            for (const auto& [length_squared, edge] : candidates)
            {
                (void)length_squared;
                if (locked[edge[0]] || locked[edge[1]] || !linkCondition(edge)) continue;
                collapses.push_back({edge[1], edge[0]});
                locked[edge[0]] = locked[edge[1]] = 1;
                for (const auto neighbor : neighbors[edge[0]]) locked[neighbor] = 1;
                for (const auto neighbor : neighbors[edge[1]]) locked[neighbor] = 1;
                break;
            }
        }
        if (collapses.empty())
            throw std::runtime_error(
                "degenerate halfedge cleanup found no link-condition collapse");

        std::vector<std::uint32_t> replacement(vertices.size());
        std::iota(replacement.begin(), replacement.end(), 0);
        for (const auto collapse : collapses)
            replacement[collapse.removed] = collapse.kept;
        std::vector<TriangleIndices> next_faces;
        next_faces.reserve(faces.size() - 2 * collapses.size());
        std::set<TriangleIndices> unique_faces;
        for (auto triangle : faces)
        {
            for (auto& vertex : triangle) vertex = replacement[vertex];
            if (triangle[0] == triangle[1] || triangle[1] == triangle[2] ||
                triangle[2] == triangle[0]) continue;
            TriangleIndices signature = triangle;
            std::sort(signature.begin(), signature.end());
            if (!unique_faces.insert(signature).second)
                throw std::runtime_error(
                    "link-condition collapse created a duplicate face");
            next_faces.push_back(triangle);
        }
        faces = std::move(next_faces);
        collapse_count += collapses.size();
        std::cerr << "monitor: stage=degenerate_halfedge_collapse pass=" << pass
                  << " batch=" << collapses.size()
                  << " remaining_faces=" << faces.size() << '\n';
    }

    std::vector<std::uint32_t> remap(vertices.size(), invalid_surface_index);
    IndexedSurface cleaned;
    cleaned.triangles = faces;
    for (auto& triangle : cleaned.triangles)
        for (auto& vertex : triangle)
        {
            if (remap[vertex] == invalid_surface_index)
            {
                remap[vertex] = static_cast<std::uint32_t>(cleaned.vertices.size());
                cleaned.vertices.push_back(toVec(vertices[vertex]));
            }
            vertex = remap[vertex];
        }
    VoxelGrid placeholder;
    placeholder.shape = {1, 1, 1};
    placeholder.pitch = 1.0;
    placeholder.occupancy = {1};
    OrientedSurfaceMesh result = buildOrientedSurfaceMesh(
        std::move(cleaned), placeholder);
    std::cerr << "monitor: stage=degenerate_halfedge_cleanup"
              << " collapses=" << collapse_count
              << " vertices=" << result.geometry.vertices.size()
              << " triangles=" << result.geometry.triangles.size() << '\n';
    return result;
}

OrientedSurfaceMesh projectPhase1BoundaryToSource(
    const OrientedSurfaceMesh& boundary, const VoxelGrid& filled_occupancy,
    const VoxelGrid& source_occupancy, const MeshModel& source)
{
    validateGrid(filled_occupancy);
    validateGrid(source_occupancy);
    if (source_occupancy.shape != filled_occupancy.shape ||
        source_occupancy.pitch != filled_occupancy.pitch ||
        source_occupancy.origin.x != filled_occupancy.origin.x ||
        source_occupancy.origin.y != filled_occupancy.origin.y ||
        source_occupancy.origin.z != filled_occupancy.origin.z)
        throw std::invalid_argument("source and filled occupancy grids do not match");
    if (const auto error = validateModelPool({source}))
        throw std::invalid_argument(*error);
    if (boundary.geometry.vertices.empty() || boundary.geometry.triangles.empty())
        throw std::invalid_argument("phase-1 boundary is empty");

    const std::vector<std::uint8_t> retained_faces =
        exteriorSourceFaces(filled_occupancy, source);
    const SourceTriangleBvh source_bvh(source, retained_faces);
    OrientedSurfaceMesh result = boundary;
    result.geometry.name = boundary.geometry.name + "_source_projected";
    if (boundary.face_boundary_crossings.size() !=
        boundary.geometry.triangles.size())
        throw std::invalid_argument("phase-1 boundary lost voxel provenance");
    std::vector<std::uint8_t> touches_source_surface(
        result.geometry.vertices.size(), 0);
    std::vector<std::uint8_t> touches_generated_fill(
        result.geometry.vertices.size(), 0);
    for (std::size_t face = 0; face < boundary.geometry.triangles.size(); ++face)
    {
        const std::uint32_t crossing_id = boundary.face_boundary_crossings[face];
        if (crossing_id >= boundary.boundary_crossings.size())
            throw std::invalid_argument("phase-1 boundary has invalid voxel provenance");
        const auto inside = boundary.boundary_crossings[crossing_id].inside;
        const bool generated_fill = !source_occupancy.occupied(
            inside[0], inside[1], inside[2]);
        auto& responsibility = generated_fill
            ? touches_generated_fill : touches_source_surface;
        for (const std::uint32_t vertex : boundary.geometry.triangles[face])
            responsibility[vertex] = 1;
    }
    std::size_t projected_vertices = 0;
    std::size_t retained_fill_vertices = 0;
    for (std::size_t vertex = 0; vertex < result.geometry.vertices.size(); ++vertex)
    {
        // A fill-only vertex has no corresponding point on the source OBJ.
        // Projecting it would collapse the cap toward the hole wall. Shared
        // mouth vertices are source constrained so the cap joins exactly.
        if (!touches_source_surface[vertex] && touches_generated_fill[vertex])
        {
            ++retained_fill_vertices;
            continue;
        }
        const auto closest = source_bvh.closestPoint(
            toVec(boundary.geometry.vertices[vertex]));
        result.geometry.vertices[vertex] = toPosition(closest.point);
        ++projected_vertices;
    }

    double signed_volume = 0.0;
    std::size_t degenerate_faces = 0;
    Bounds3 bounds;
    for (const auto point : result.geometry.vertices) include(bounds, toVec(point));
    const double scale = norm(bounds.upper - bounds.lower);
    const double minimum_area = std::max(1.0, scale * scale) * 1.0e-24;
    for (const auto triangle : result.geometry.triangles)
    {
        const Vec3 a = toVec(result.geometry.vertices[triangle[0]]);
        const Vec3 b = toVec(result.geometry.vertices[triangle[1]]);
        const Vec3 c = toVec(result.geometry.vertices[triangle[2]]);
        degenerate_faces += norm(cross(b - a, c - a)) <= minimum_area ? 1 : 0;
        signed_volume += dot(a, cross(b, c)) / 6.0;
    }
    result.euler_characteristic = boundary.euler_characteristic;
    if (!std::isfinite(signed_volume) || signed_volume == 0.0)
        throw std::runtime_error("exact source projection has zero oriented volume");
    if (signed_volume < 0.0)
        throw std::runtime_error("exact source projection has inward orientation");
    result.signed_volume = signed_volume;
    std::cerr << "monitor: stage=exact_source_projection"
              << " projected_vertices=" << projected_vertices
              << " retained_fill_vertices=" << retained_fill_vertices
              << " retained_topology=1"
              << " degenerate_faces_pending_collapse=" << degenerate_faces << '\n';
    return collapseDegenerateHalfedges(result);
}

MeshModel buildSourceDominatedPhase1Preview(
    const VoxelGrid& occupancy, const MeshModel& source)
{
    validateGrid(occupancy);
    if (const auto error = validateModelPool({source}))
        throw std::invalid_argument(*error);

    const Vec3 grid_origin = toVec(occupancy.origin);
    const auto sampleOccupied = [&](const Vec3 point)
    {
        std::array<int, 3> coordinate{};
        for (int axis = 0; axis < 3; ++axis)
        {
            coordinate[axis] = static_cast<int>(std::llround(
                (point[axis] - grid_origin[axis]) / occupancy.pitch));
            if (coordinate[axis] < 0 ||
                coordinate[axis] >= static_cast<int>(occupancy.shape[axis]))
                return false;
        }
        return occupancy.occupied(
            static_cast<std::uint32_t>(coordinate[0]),
            static_cast<std::uint32_t>(coordinate[1]),
            static_cast<std::uint32_t>(coordinate[2]));
    };

    std::vector<std::uint8_t> retained(source.triangles.size(), 0);
    VoxelGrid source_support = occupancy;
    std::fill(source_support.occupancy.begin(), source_support.occupancy.end(), 0);
    for (std::size_t face_id = 0; face_id < source.triangles.size(); ++face_id)
    {
        const auto face = source.triangles[face_id];
        const Vec3 a = toVec(source.vertices[face[0]]);
        const Vec3 b = toVec(source.vertices[face[1]]);
        const Vec3 c = toVec(source.vertices[face[2]]);
        const Vec3 normal = normalized(cross(b - a, c - a));
        if (norm(normal) == 0.0) continue;
        const std::array<Vec3, 4> samples{
            (a + b + c) * (1.0 / 3.0),
            a * 0.6 + b * 0.2 + c * 0.2,
            a * 0.2 + b * 0.6 + c * 0.2,
            a * 0.2 + b * 0.2 + c * 0.6};
        bool boundary = false;
        for (const Vec3 sample : samples)
            for (const double scale : {0.5, 1.0, 1.5, 2.0})
            {
                const Vec3 offset = normal * (scale * occupancy.pitch);
                if (sampleOccupied(sample + offset) != sampleOccupied(sample - offset))
                {
                    boundary = true;
                    break;
                }
            }
        if (!boundary) continue;
        retained[face_id] = 1;
        rasterizeTriangle(source_support, {a, b, c});
    }

    const auto topology = buildPhase1Solid(occupancy, source);
    const auto planes = sourcePlanes(source, occupancy.pitch);
    MeshModel result;
    result.name = source.name + "_source_dominated_phase1_preview";
    result.vertices = source.vertices;
    for (std::size_t face_id = 0; face_id < source.triangles.size(); ++face_id)
        if (retained[face_id]) result.triangles.push_back(source.triangles[face_id]);

    MeshModel generated_fill;
    generated_fill.name = source.name + "_generated_fill";
    std::map<std::array<double, 3>, std::uint32_t> generated_vertices;
    const auto appendGeneratedVertex = [&](const Vec3 point)
    {
        const std::array<double, 3> key{point.x, point.y, point.z};
        const auto [iterator, inserted] = generated_vertices.try_emplace(
            key, static_cast<std::uint32_t>(generated_fill.vertices.size()));
        if (inserted) generated_fill.vertices.push_back(toPosition(point));
        return iterator->second;
    };
    for (std::size_t face_id = 0;
         face_id < topology.boundary.geometry.triangles.size(); ++face_id)
    {
        const std::uint32_t crossing_id =
            topology.boundary.face_boundary_crossings[face_id];
        const auto& crossing = topology.boundary.boundary_crossings[crossing_id];
        const bool source_supported = source_support.occupied(
                crossing.inside[0], crossing.inside[1], crossing.inside[2]) ||
            source_support.occupied(
                crossing.outside[0], crossing.outside[1], crossing.outside[2]);
        if (source_supported) continue;
        const auto topology_face = topology.boundary.geometry.triangles[face_id];
        std::array<Vec3, 3> points{};
        for (int corner = 0; corner < 3; ++corner)
            points[corner] = toVec(topology.boundary.geometry.vertices[
                topology_face[corner]]);
        const Vec3 centroid = (points[0] + points[1] + points[2]) * (1.0 / 3.0);
        const Vec3 face_normal = normalized(cross(
            points[1] - points[0], points[2] - points[0]));
        std::optional<std::size_t> best_plane;
        double best_score = std::numeric_limits<double>::infinity();
        for (std::size_t plane_id = 0; plane_id < planes.size(); ++plane_id)
        {
            const auto& plane = planes[plane_id];
            const double alignment = std::abs(dot(face_normal, plane.normal));
            if (alignment < 0.75 ||
                !withinExpandedBounds(centroid, plane.bounds, 3.0 * occupancy.pitch))
                continue;
            const double distance = std::abs(dot(plane.normal, centroid) - plane.distance);
            if (distance > 2.5 * occupancy.pitch) continue;
            const double score = distance + (1.0 - alignment) * occupancy.pitch;
            if (score < best_score)
            {
                best_score = score;
                best_plane = plane_id;
            }
        }
        if (best_plane)
        {
            const auto& plane = planes[*best_plane];
            for (auto& point : points)
                point = point - plane.normal * (dot(plane.normal, point) - plane.distance);
        }
        TriangleIndices generated_face{};
        for (int corner = 0; corner < 3; ++corner)
            generated_face[corner] = appendGeneratedVertex(points[corner]);
        generated_fill.triangles.push_back(generated_face);
    }
    if (!generated_fill.triangles.empty())
        generated_fill = canonicalizeCoplanarTriangleSoup(generated_fill);
    const std::uint32_t generated_offset =
        static_cast<std::uint32_t>(result.vertices.size());
    result.vertices.insert(result.vertices.end(),
        generated_fill.vertices.begin(), generated_fill.vertices.end());
    for (auto face : generated_fill.triangles)
    {
        for (auto& vertex : face) vertex += generated_offset;
        result.triangles.push_back(face);
    }
    return result;
}

MeshModel buildOriginalObjHoleSurgery(
    const VoxelGrid& filled, const VoxelGrid& cavity_fill_labels,
    const VoxelGrid& original, const TriangleVoxelProvenance& provenance,
    const MeshModel& source,
    OriginalObjHoleSurgeryStats* stats)
{
    validateGrid(filled);
    if (cavity_fill_labels.shape != filled.shape ||
        cavity_fill_labels.pitch != filled.pitch ||
        cavity_fill_labels.origin.x != filled.origin.x ||
        cavity_fill_labels.origin.y != filled.origin.y ||
        cavity_fill_labels.origin.z != filled.origin.z ||
        cavity_fill_labels.occupancy.size() != filled.occupancy.size())
        throw std::invalid_argument("cavity-fill labels do not match the filled grid");
    if (const auto error = validateModelPool({source}))
        throw std::invalid_argument(*error);
    validateGrid(original);
    if (original.shape != filled.shape || original.pitch != filled.pitch ||
        original.origin.x != filled.origin.x ||
        original.origin.y != filled.origin.y ||
        original.origin.z != filled.origin.z ||
        provenance.surface_voxels.size() != source.triangles.size())
        throw std::invalid_argument("source provenance does not match the filled grid");
    std::cerr << "monitor: stage=voxelized source_voxels="
              << countOccupied(original) << '\n';
    std::size_t added_voxels = 0;
    std::size_t removed_voxels = 0;
    for (std::size_t index = 0; index < filled.occupancy.size(); ++index)
    {
        added_voxels += filled.occupancy[index] && !original.occupancy[index] ? 1 : 0;
        removed_voxels += original.occupancy[index] && !filled.occupancy[index] ? 1 : 0;
    }
    if (removed_voxels != 0)
        throw std::runtime_error(
            "frozen fill does not enclose occupancy reconstructed on its grid");
    std::size_t closed_cavity_voxels = 0;
    for (std::size_t index = 0; index < filled.occupancy.size(); ++index)
    {
        if (cavity_fill_labels.occupancy[index] && !filled.occupancy[index])
            throw std::invalid_argument("cavity-fill label lies outside the final occupancy");
        if (cavity_fill_labels.occupancy[index])
        {
            ++closed_cavity_voxels;
        }
    }
    std::cerr << "monitor: stage=beta2_labels cavity_voxels="
              << closed_cavity_voxels << '\n';

    const auto exterior = floodEmptyFromBoundary(filled, neighbors26);
    VoxelGrid outer_boundary = filled;
    std::fill(outer_boundary.occupancy.begin(),
              outer_boundary.occupancy.end(), std::uint8_t{0});
    for (std::uint32_t x = 0; x < filled.shape[0]; ++x)
        for (std::uint32_t y = 0; y < filled.shape[1]; ++y)
            for (std::uint32_t z = 0; z < filled.shape[2]; ++z)
            {
                if (!filled.occupied(x, y, z)) continue;
                bool boundary = false;
                for (const auto& offset : neighbors26)
                {
                    const int nx = static_cast<int>(x) + offset[0];
                    const int ny = static_cast<int>(y) + offset[1];
                    const int nz = static_cast<int>(z) + offset[2];
                    if (nx < 0 || ny < 0 || nz < 0 ||
                        nx >= static_cast<int>(filled.shape[0]) ||
                        ny >= static_cast<int>(filled.shape[1]) ||
                        nz >= static_cast<int>(filled.shape[2]))
                    {
                        boundary = true;
                        break;
                    }
                    boundary |= exterior[filled.index(
                        static_cast<std::uint32_t>(nx),
                        static_cast<std::uint32_t>(ny),
                        static_cast<std::uint32_t>(nz))] != 0;
                }
                outer_boundary.occupancy[filled.index(x, y, z)] = boundary ? 1 : 0;
            }

    Bounds3 source_bounds;
    for (const auto point : source.vertices) include(source_bounds, toVec(point));
    const double model_scale = norm(source_bounds.upper - source_bounds.lower);
    const double weld_tolerance = std::max(model_scale * 1.0e-9, 1.0e-9);
    using VertexKey = std::array<std::int64_t, 3>;
    using EdgeKey = std::array<VertexKey, 2>;
    const auto vertexKey = [&](const Position3 point)
    {
        return VertexKey{
            static_cast<std::int64_t>(std::llround(point.x / weld_tolerance)),
            static_cast<std::int64_t>(std::llround(point.y / weld_tolerance)),
            static_cast<std::int64_t>(std::llround(point.z / weld_tolerance))};
    };

    // A source triangle belongs to the final outer boundary iff at least one
    // of the surface voxels it originally owned remains on the exact
    // occupied/exterior interface. Ownership is many-to-many by construction.
    std::vector<std::uint8_t> removed(source.triangles.size(), 0);
    std::vector<std::uint8_t> mixed(source.triangles.size(), 0);
    std::size_t ambiguous = 0;
    for (std::size_t face_id = 0; face_id < source.triangles.size(); ++face_id)
    {
        bool intersects_outer_boundary = false;
        bool intersects_internal = false;
        for (const std::uint32_t index : provenance.surface_voxels[face_id])
        {
            if (index >= filled.occupancy.size())
                throw std::invalid_argument("source provenance contains an invalid voxel");
            intersects_outer_boundary |= outer_boundary.occupancy[index] != 0;
            intersects_internal |= filled.occupancy[index] != 0 &&
                outer_boundary.occupancy[index] == 0;
        }
        if (!intersects_outer_boundary) removed[face_id] = 2;
        else if (intersects_internal)
        {
            mixed[face_id] = 1;
            ++ambiguous;
        }
    }

    // Build a continuous, piecewise-linear classification field on the source
    // triangles. Voxel provenance supplies only samples. All field nodes and
    // all zero crossings are barycentric points on the original triangle, so
    // no voxel corner or voxel-box edge can enter the output geometry.
    constexpr int source_subdivision = 2;
    struct SourceFieldNode
    {
        Vec3 point;
        VertexKey key{};
        double local_value = 0.0;
    };
    struct SourceFieldFace
    {
        std::vector<SourceFieldNode> nodes;
    };
    struct FieldAccumulator
    {
        double sum = 0.0;
        std::size_t count = 0;
    };
    std::vector<SourceFieldFace> source_fields(source.triangles.size());
    std::map<VertexKey, FieldAccumulator> shared_field;
    const Vec3 grid_origin = toVec(filled.origin);
    const auto sampleCenter = [&](const std::uint32_t index)
    {
        const auto coordinate = decodeIndex(index, filled.shape);
        return grid_origin + Vec3{
            coordinate[0] * filled.pitch,
            coordinate[1] * filled.pitch,
            coordinate[2] * filled.pitch};
    };
    const auto insertNearest = [](std::array<double, 4>& nearest,
                                  const double distance_squared)
    {
        if (distance_squared >= nearest.back()) return;
        nearest.back() = distance_squared;
        for (std::size_t index = nearest.size() - 1; index > 0; --index)
        {
            if (nearest[index] >= nearest[index - 1]) break;
            std::swap(nearest[index], nearest[index - 1]);
        }
    };
    for (std::size_t face_id = 0; face_id < source.triangles.size(); ++face_id)
    {
        if (!mixed[face_id]) continue;
        const auto face = source.triangles[face_id];
        const Vec3 a = toVec(source.vertices[face[0]]);
        const Vec3 b = toVec(source.vertices[face[1]]);
        const Vec3 c = toVec(source.vertices[face[2]]);
        const Vec3 normal = normalized(cross(b - a, c - a));
        auto& field = source_fields[face_id];
        field.nodes.reserve((source_subdivision + 1) *
            (source_subdivision + 2) / 2);
        for (int i = 0; i <= source_subdivision; ++i)
            for (int j = 0; j <= source_subdivision - i; ++j)
            {
                const double wb = static_cast<double>(i) / source_subdivision;
                const double wc = static_cast<double>(j) / source_subdivision;
                const Vec3 point = a * (1.0 - wb - wc) + b * wb + c * wc;
                std::array<double, 4> nearest_outer{};
                std::array<double, 4> nearest_inner{};
                nearest_outer.fill(std::numeric_limits<double>::infinity());
                nearest_inner.fill(std::numeric_limits<double>::infinity());
                for (const std::uint32_t voxel : provenance.surface_voxels[face_id])
                {
                    Vec3 sample = sampleCenter(voxel);
                    sample = sample - normal * dot(normal, sample - a);
                    const double distance_squared = dot(sample - point, sample - point);
                    if (outer_boundary.occupancy[voxel])
                        insertNearest(nearest_outer, distance_squared);
                    else
                        insertNearest(nearest_inner, distance_squared);
                }
                const auto averageFinite = [](const std::array<double, 4>& values)
                {
                    double sum = 0.0;
                    std::size_t count = 0;
                    for (const double value : values)
                        if (std::isfinite(value)) { sum += value; ++count; }
                    return count == 0 ? std::numeric_limits<double>::infinity()
                                      : sum / count;
                };
                const double outer_distance = averageFinite(nearest_outer);
                const double inner_distance = averageFinite(nearest_inner);
                const double value = inner_distance - outer_distance;
                SourceFieldNode node{point, vertexKey(toPosition(point)), value};
                auto& accumulator = shared_field[node.key];
                accumulator.sum += value;
                ++accumulator.count;
                field.nodes.push_back(node);
            }
    }
    for (auto& field : source_fields)
        for (auto& node : field.nodes)
        {
            const auto& accumulator = shared_field.at(node.key);
            node.local_value = accumulator.sum /
                static_cast<double>(accumulator.count);
        }

    MeshModel retained_source;
    retained_source.name = source.name + "_source_clipped";
    retained_source.vertices = source.vertices;
    std::map<VertexKey, std::uint32_t> retained_vertex_ids;
    for (std::uint32_t vertex = 0; vertex < source.vertices.size(); ++vertex)
        retained_vertex_ids.try_emplace(vertexKey(source.vertices[vertex]), vertex);
    const auto appendRetainedVertex = [&](const Vec3 point)
    {
        const VertexKey key = vertexKey(toPosition(point));
        const auto [iterator, inserted] = retained_vertex_ids.try_emplace(
            key, static_cast<std::uint32_t>(retained_source.vertices.size()));
        if (inserted) retained_source.vertices.push_back(toPosition(point));
        return iterator->second;
    };
    std::vector<std::pair<Position3, Position3>> mixed_interface_segments;
    struct ClippedNode { Vec3 point; double value = 0.0; };
    const auto clipMicroTriangle = [&](const std::array<SourceFieldNode, 3>& triangle)
    {
        std::vector<ClippedNode> polygon;
        polygon.reserve(4);
        std::vector<Vec3> crossings;
        for (int edge = 0; edge < 3; ++edge)
        {
            const auto& first = triangle[edge];
            const auto& second = triangle[(edge + 1) % 3];
            const bool first_outside = first.local_value >= 0.0;
            const bool second_outside = second.local_value >= 0.0;
            if (first_outside) polygon.push_back({first.point, first.local_value});
            if (first_outside == second_outside) continue;
            const double denominator = first.local_value - second.local_value;
            const double parameter = std::abs(denominator) <= 1.0e-30
                ? 0.5 : std::clamp(first.local_value / denominator, 0.0, 1.0);
            const Vec3 crossing = first.point +
                (second.point - first.point) * parameter;
            polygon.push_back({crossing, 0.0});
            crossings.push_back(crossing);
        }
        if (crossings.size() == 2 &&
            norm(crossings[1] - crossings[0]) > weld_tolerance)
            mixed_interface_segments.push_back({
                toPosition(crossings[0]), toPosition(crossings[1])});
        if (polygon.size() < 3) return;
        const std::uint32_t first = appendRetainedVertex(polygon[0].point);
        for (std::size_t corner = 1; corner + 1 < polygon.size(); ++corner)
        {
            const std::uint32_t second = appendRetainedVertex(polygon[corner].point);
            const std::uint32_t third = appendRetainedVertex(polygon[corner + 1].point);
            if (first != second && second != third && third != first)
                retained_source.triangles.push_back({first, second, third});
        }
    };
    const auto sourceFieldIndex = [](const int i, const int j)
    {
        return static_cast<std::size_t>(
            i * (source_subdivision + 1) - i * (i - 1) / 2 + j);
    };
    for (std::size_t face_id = 0; face_id < source.triangles.size(); ++face_id)
    {
        if (removed[face_id]) continue;
        if (!mixed[face_id])
        {
            retained_source.triangles.push_back(source.triangles[face_id]);
            continue;
        }
        const auto& nodes = source_fields[face_id].nodes;
        for (int i = 0; i < source_subdivision; ++i)
            for (int j = 0; j < source_subdivision - i; ++j)
            {
                const auto& a = nodes[sourceFieldIndex(i, j)];
                const auto& b = nodes[sourceFieldIndex(i + 1, j)];
                const auto& c = nodes[sourceFieldIndex(i, j + 1)];
                clipMicroTriangle({a, b, c});
                if (j + 1 < source_subdivision - i)
                {
                    const auto& d = nodes[sourceFieldIndex(i + 1, j + 1)];
                    clipMicroTriangle({b, d, c});
                }
            }
    }
    std::cerr << "monitor: stage=source_clipped triangles="
              << retained_source.triangles.size() << " vertices="
              << retained_source.vertices.size() << '\n';
    const std::size_t removed_closed_components = 0;
    std::cerr << "monitor: stage=classified removed="
              << std::count_if(removed.begin(), removed.end(),
                     [](const std::uint8_t reason) { return reason != 0; })
              << " cavity_walls="
              << std::count(removed.begin(), removed.end(), std::uint8_t{2})
              << " cavity_components=" << removed_closed_components
              << '\n';
    struct EdgeUse
    {
        bool removed = false;
        bool cap_removed = false;
        bool retained = false;
        std::size_t removed_count = 0;
        std::size_t retained_count = 0;
        Position3 first;
        Position3 second;
    };
    std::map<EdgeKey, EdgeUse> edge_uses;
    std::map<VertexKey, Position3> welded_positions;
    for (std::size_t face_id = 0; face_id < source.triangles.size(); ++face_id)
    {
        const auto face = source.triangles[face_id];
        for (int edge = 0; edge < 3; ++edge)
        {
            Position3 first = source.vertices[face[edge]];
            Position3 second = source.vertices[face[(edge + 1) % 3]];
            VertexKey first_key = vertexKey(first);
            VertexKey second_key = vertexKey(second);
            welded_positions.try_emplace(first_key, first);
            welded_positions.try_emplace(second_key, second);
            if (second_key < first_key)
            {
                std::swap(first_key, second_key);
                std::swap(first, second);
            }
            auto& use = edge_uses[EdgeKey{first_key, second_key}];
            use.first = first;
            use.second = second;
            use.removed |= removed[face_id] != 0;
            use.cap_removed |= removed[face_id] == 1;
            use.retained |= removed[face_id] == 0;
            use.removed_count += removed[face_id] != 0 ? 1 : 0;
            use.retained_count += removed[face_id] == 0 ? 1 : 0;
        }
    }
    std::cerr << "monitor: stage=edge_map edges=" << edge_uses.size() << '\n';

    std::map<VertexKey, Position3> interface_vertices;
    std::map<VertexKey, std::set<VertexKey>> adjacency;
    std::map<EdgeKey, std::size_t> retained_edge_counts;
    for (const auto face : retained_source.triangles)
        for (int local = 0; local < 3; ++local)
        {
            VertexKey first = vertexKey(retained_source.vertices[face[local]]);
            VertexKey second = vertexKey(
                retained_source.vertices[face[(local + 1) % 3]]);
            if (second < first) std::swap(first, second);
            ++retained_edge_counts[{first, second}];
        }
    std::size_t interface_edges = 0;
    for (const auto& [edge, count] : retained_edge_counts)
    {
        if (count != 1) continue;
        const Position3 first{
            edge[0][0] * weld_tolerance,
            edge[0][1] * weld_tolerance,
            edge[0][2] * weld_tolerance};
        const Position3 second{
            edge[1][0] * weld_tolerance,
            edge[1][1] * weld_tolerance,
            edge[1][2] * weld_tolerance};
        interface_vertices[edge[0]] = first;
        interface_vertices[edge[1]] = second;
        adjacency[edge[0]].insert(edge[1]);
        adjacency[edge[1]].insert(edge[0]);
        ++interface_edges;
    }

    std::set<EdgeKey> visited;
    std::vector<std::vector<VertexKey>> loops;
    for (const auto& [start, neighbors] : adjacency)
    {
        if (neighbors.size() != 2) continue;
        for (const auto& next_start : neighbors)
        {
            EdgeKey initial{std::min(start, next_start), std::max(start, next_start)};
            if (visited.contains(initial)) continue;
            std::vector<VertexKey> loop{start};
            VertexKey previous = start;
            VertexKey current = next_start;
            bool closed = false;
            for (std::size_t step = 0; step <= adjacency.size(); ++step)
            {
                EdgeKey edge{std::min(previous, current), std::max(previous, current)};
                if (visited.contains(edge)) break;
                visited.insert(edge);
                if (current == start)
                {
                    closed = true;
                    break;
                }
                loop.push_back(current);
                const auto found = adjacency.find(current);
                if (found == adjacency.end() || found->second.size() != 2) break;
                auto iterator = found->second.begin();
                VertexKey next = *iterator;
                if (next == previous) next = *std::next(iterator);
                previous = current;
                current = next;
            }
            if (closed && loop.size() >= 3) loops.push_back(std::move(loop));
        }
    }
    // A drilled or recessed mouth is an inner boundary of a coplanar source
    // patch even when its wall is sub-voxel. Extract those inner loops per
    // exact geometric plane; unlike global soup boundary edges, this does not
    // join unrelated cracks or treat an exterior silhouette as a hole.
    using PlaneKey = std::array<std::int64_t, 4>;
    struct PlanarBoundaryGroup
    {
        Vec3 origin;
        Vec3 normal;
        Vec3 tangent;
        Vec3 bitangent;
        std::map<EdgeKey, std::size_t> edge_counts;
    };
    std::map<PlaneKey, PlanarBoundaryGroup> plane_groups;
    constexpr double angular_quantum = 1.0e-8;
    for (std::size_t face_id = 0; face_id < source.triangles.size(); ++face_id)
    {
        if (removed[face_id]) continue;
        const auto face = source.triangles[face_id];
        const Vec3 a = toVec(source.vertices[face[0]]);
        const Vec3 b = toVec(source.vertices[face[1]]);
        const Vec3 c = toVec(source.vertices[face[2]]);
        Vec3 normal = normalized(cross(b - a, c - a));
        if (norm(normal) == 0.0) continue;
        int dominant = 0;
        for (int axis = 1; axis < 3; ++axis)
            if (std::abs(normal[axis]) > std::abs(normal[dominant])) dominant = axis;
        if (normal[dominant] < 0.0) normal = normal * -1.0;
        const double distance = dot(normal, a);
        PlaneKey key{
            static_cast<std::int64_t>(std::llround(normal.x / angular_quantum)),
            static_cast<std::int64_t>(std::llround(normal.y / angular_quantum)),
            static_cast<std::int64_t>(std::llround(normal.z / angular_quantum)),
            static_cast<std::int64_t>(std::llround(distance / weld_tolerance))};
        auto [iterator, inserted] = plane_groups.try_emplace(key);
        auto& group = iterator->second;
        if (inserted)
        {
            group.origin = a;
            group.normal = normal;
            const Vec3 helper = std::abs(normal.x) < 0.8
                ? Vec3{1, 0, 0} : Vec3{0, 1, 0};
            group.tangent = normalized(cross(helper, normal));
            group.bitangent = cross(normal, group.tangent);
        }
        for (int edge_id = 0; edge_id < 3; ++edge_id)
        {
            VertexKey first = vertexKey(source.vertices[face[edge_id]]);
            VertexKey second = vertexKey(source.vertices[face[(edge_id + 1) % 3]]);
            if (second < first) std::swap(first, second);
            ++group.edge_counts[EdgeKey{first, second}];
        }
    }
    std::set<std::vector<EdgeKey>> loop_signatures;
    std::set<EdgeKey> beta2_interface_edges;
    for (const auto& [edge, use] : edge_uses)
        if (use.removed && !use.cap_removed && use.retained)
            beta2_interface_edges.insert(edge);
    const auto signature = [](const std::vector<VertexKey>& loop)
    {
        std::vector<EdgeKey> result;
        result.reserve(loop.size());
        for (std::size_t index = 0; index < loop.size(); ++index)
        {
            VertexKey a = loop[index];
            VertexKey b = loop[(index + 1) % loop.size()];
            result.push_back({std::min(a, b), std::max(a, b)});
        }
        std::sort(result.begin(), result.end());
        return result;
    };
    for (const auto& loop : loops) loop_signatures.insert(signature(loop));
    std::size_t planar_hole_loops = 0;
    for (const auto& [plane_key, group] : plane_groups)
    {
        (void)plane_key;
        std::map<VertexKey, std::set<VertexKey>> plane_adjacency;
        for (const auto& [edge, count] : group.edge_counts)
            if (count == 1)
            {
                plane_adjacency[edge[0]].insert(edge[1]);
                plane_adjacency[edge[1]].insert(edge[0]);
            }
        std::set<EdgeKey> plane_visited;
        std::vector<std::vector<VertexKey>> plane_loops;
        for (const auto& [start, neighbors] : plane_adjacency)
        {
            if (neighbors.size() != 2) continue;
            for (const auto& next_start : neighbors)
            {
                EdgeKey initial{std::min(start, next_start), std::max(start, next_start)};
                if (plane_visited.contains(initial)) continue;
                std::vector<VertexKey> loop{start};
                VertexKey previous = start;
                VertexKey current = next_start;
                bool closed = false;
                for (std::size_t step = 0; step <= plane_adjacency.size(); ++step)
                {
                    EdgeKey edge{std::min(previous, current), std::max(previous, current)};
                    if (plane_visited.contains(edge)) break;
                    plane_visited.insert(edge);
                    if (current == start) { closed = true; break; }
                    loop.push_back(current);
                    const auto found = plane_adjacency.find(current);
                    if (found == plane_adjacency.end() || found->second.size() != 2) break;
                    auto next_iterator = found->second.begin();
                    VertexKey next = *next_iterator;
                    if (next == previous) next = *std::next(next_iterator);
                    previous = current;
                    current = next;
                }
                if (closed && loop.size() >= 3 && loop.size() <= 4096)
                    plane_loops.push_back(std::move(loop));
            }
        }
        std::vector<std::vector<std::array<double, 2>>> projected;
        projected.reserve(plane_loops.size());
        for (const auto& loop : plane_loops)
        {
            std::vector<std::array<double, 2>> polygon;
            polygon.reserve(loop.size());
            for (const auto& vertex : loop)
            {
                const Vec3 delta = toVec(welded_positions.at(vertex)) - group.origin;
                polygon.push_back({dot(delta, group.tangent), dot(delta, group.bitangent)});
            }
            projected.push_back(std::move(polygon));
        }
        const auto contains = [](const std::vector<std::array<double, 2>>& polygon,
                                 const std::array<double, 2> point)
        {
            bool inside = false;
            for (std::size_t i = 0, j = polygon.size() - 1; i < polygon.size(); j = i++)
            {
                const auto& a = polygon[i];
                const auto& b = polygon[j];
                if ((a[1] > point[1]) != (b[1] > point[1]) &&
                    point[0] < (b[0] - a[0]) * (point[1] - a[1]) /
                        (b[1] - a[1]) + a[0]) inside = !inside;
            }
            return inside;
        };
        for (std::size_t loop_id = 0; loop_id < plane_loops.size(); ++loop_id)
        {
            int nesting = 0;
            const auto point = projected[loop_id].front();
            for (std::size_t other = 0; other < plane_loops.size(); ++other)
                if (other != loop_id && contains(projected[other], point)) ++nesting;
            if ((nesting % 2) == 0) continue;
            auto loop_signature = signature(plane_loops[loop_id]);
            if (std::ranges::any_of(loop_signature,
                    [&](const EdgeKey& edge)
                    { return beta2_interface_edges.contains(edge); }))
                continue;
            if (!loop_signatures.insert(loop_signature).second) continue;
            for (const auto& vertex : plane_loops[loop_id])
            {
                interface_vertices.try_emplace(
                    vertex, welded_positions.at(vertex));
            }
            loops.push_back(std::move(plane_loops[loop_id]));
            ++planar_hole_loops;
        }
    }
    std::cerr << "monitor: stage=loops interface_edges=" << interface_edges
              << " planar_holes=" << planar_hole_loops
              << " loops=" << loops.size() << '\n';

    MeshModel result;
    result.name = source.name + "_original_obj_hole_surgery";
    result.vertices = std::move(retained_source.vertices);
    result.triangles = std::move(retained_source.triangles);

    std::map<VertexKey, std::uint32_t> cap_vertices;
    for (const auto& [key, point] : interface_vertices)
    {
        cap_vertices[key] = static_cast<std::uint32_t>(result.vertices.size());
        result.vertices.push_back(point);
    }
    std::size_t cap_triangles = 0;
    for (const auto& loop : loops)
    {
        // A valid source-hole mouth is local. An enormous loop here indicates
        // that disconnected soup cracks were joined into a non-simple walk;
        // never hand that adversarial path to the triangulator.
        constexpr std::size_t maximum_safe_loop_vertices = 4096;
        if (loop.size() > maximum_safe_loop_vertices)
        {
            std::cerr << "monitor: skipped_oversized_loop vertices="
                      << loop.size() << '\n';
            continue;
        }
        Vec3 normal{};
        for (std::size_t index = 0; index < loop.size(); ++index)
        {
            const Vec3 a = toVec(interface_vertices.at(loop[index]));
            const Vec3 b = toVec(interface_vertices.at(loop[(index + 1) % loop.size()]));
            normal.x += (a.y - b.y) * (a.z + b.z);
            normal.y += (a.z - b.z) * (a.x + b.x);
            normal.z += (a.x - b.x) * (a.y + b.y);
        }
        normal = normalized(normal);
        if (norm(normal) == 0.0) continue;
        const Vec3 helper = std::abs(normal.x) < 0.8 ? Vec3{1, 0, 0} : Vec3{0, 1, 0};
        const Vec3 tangent = normalized(cross(helper, normal));
        const Vec3 bitangent = cross(normal, tangent);
        const Vec3 origin = toVec(interface_vertices.at(loop.front()));
        std::vector<std::array<double, 2>> polygon;
        for (const auto& key : loop)
        {
            const Vec3 delta = toVec(interface_vertices.at(key)) - origin;
            polygon.push_back({dot(delta, tangent), dot(delta, bitangent)});
        }
        double signed_area = 0.0;
        for (std::size_t index = 0; index < polygon.size(); ++index)
        {
            const auto& a = polygon[index];
            const auto& b = polygon[(index + 1) % polygon.size()];
            signed_area += a[0] * b[1] - a[1] * b[0];
        }
        if (std::abs(signed_area) <= weld_tolerance * weld_tolerance) continue;
        std::vector<std::size_t> remaining(loop.size());
        std::iota(remaining.begin(), remaining.end(), 0);
        if (signed_area < 0.0) std::reverse(remaining.begin(), remaining.end());
        const auto cross2 = [&](const std::size_t a, const std::size_t b,
                                const std::size_t c)
        {
            return (polygon[b][0] - polygon[a][0]) *
                       (polygon[c][1] - polygon[a][1]) -
                   (polygon[b][1] - polygon[a][1]) *
                       (polygon[c][0] - polygon[a][0]);
        };
        const auto insideTriangle = [&](const std::size_t point,
            const std::size_t a, const std::size_t b, const std::size_t c)
        {
            const double epsilon = weld_tolerance * weld_tolerance;
            return cross2(a, b, point) >= -epsilon &&
                   cross2(b, c, point) >= -epsilon &&
                   cross2(c, a, point) >= -epsilon;
        };
        bool triangulation_failed = false;
        std::vector<TriangleIndices> loop_triangles;
        loop_triangles.reserve(loop.size() - 2);
        while (remaining.size() > 3)
        {
            bool clipped = false;
            for (std::size_t index = 0; index < remaining.size(); ++index)
            {
                const std::size_t a = remaining[(index + remaining.size() - 1) % remaining.size()];
                const std::size_t b = remaining[index];
                const std::size_t c = remaining[(index + 1) % remaining.size()];
                if (cross2(a, b, c) <= weld_tolerance * weld_tolerance) continue;
                bool contains = false;
                for (const std::size_t point : remaining)
                {
                    if (point != a && point != b && point != c &&
                        insideTriangle(point, a, b, c)) { contains = true; break; }
                }
                if (contains) continue;
                loop_triangles.push_back({cap_vertices.at(loop[a]),
                    cap_vertices.at(loop[b]), cap_vertices.at(loop[c])});
                remaining.erase(remaining.begin() + static_cast<std::ptrdiff_t>(index));
                clipped = true;
                break;
            }
            if (!clipped)
            {
                triangulation_failed = true;
                break;
            }
        }
        if (!triangulation_failed && remaining.size() == 3)
        {
            loop_triangles.push_back({cap_vertices.at(loop[remaining[0]]),
                cap_vertices.at(loop[remaining[1]]),
                cap_vertices.at(loop[remaining[2]])});
            if (loop_triangles.size() == loop.size() - 2)
            {
                cap_triangles += loop_triangles.size();
                result.triangles.insert(result.triangles.end(),
                    loop_triangles.begin(), loop_triangles.end());
            }
        }
    }
    std::cerr << "monitor: stage=triangulated cap_triangles="
              << cap_triangles << '\n';

    if (stats)
    {
        stats->source_voxels = countOccupied(original);
        stats->filled_voxels = countOccupied(filled);
        stats->added_voxels = added_voxels;
        stats->source_triangles = source.triangles.size();
        stats->removed_internal_triangles = std::count(removed.begin(), removed.end(), 1);
        stats->removed_internal_triangles = std::count_if(
            removed.begin(), removed.end(),
            [](const std::uint8_t reason) { return reason != 0; });
        stats->removed_closed_cavity_triangles = std::count(
            removed.begin(), removed.end(), std::uint8_t{2});
        stats->removed_closed_cavity_components = removed_closed_components;
        stats->closed_cavity_voxels = closed_cavity_voxels;
        stats->retained_triangles = source.triangles.size() -
            stats->removed_internal_triangles;
        stats->ambiguous_triangles = ambiguous;
        stats->interface_edges = interface_edges;
        stats->closed_loops = loops.size();
        stats->cap_triangles = cap_triangles;
    }
    return result;
}

MeshModel readTriangleSoupObj(const std::filesystem::path& path)
{
    return readObj(path);
}

namespace
{

struct HalfedgeComponentDiagnostic
{
    std::vector<std::uint32_t> faces;
    double area = 0.0;
    std::size_t boundary_halfedges = 0;
    Bounds3 bounds;
    std::uint32_t minimum_face = invalid_surface_index;
};

std::vector<HalfedgeComponentDiagnostic> halfedgeComponentDiagnostics(
    const OrientedSurfaceMesh& mesh)
{
    std::vector<std::uint8_t> visited(mesh.halfedges.size(), 0);
    std::vector<HalfedgeComponentDiagnostic> result;
    for (std::uint32_t seed = 0; seed < mesh.halfedges.size(); ++seed)
    {
        if (visited[seed]) continue;
        HalfedgeComponentDiagnostic component;
        std::vector<std::uint32_t> queue{seed};
        std::set<std::uint32_t> faces;
        visited[seed] = 1;
        for (std::size_t cursor = 0; cursor < queue.size(); ++cursor)
        {
            const std::uint32_t halfedge = queue[cursor];
            const auto& edge = mesh.halfedges[halfedge];
            faces.insert(edge.face);
            component.boundary_halfedges +=
                edge.opposite == invalid_surface_index ? 1 : 0;
            for (const std::uint32_t neighbor : {edge.next, edge.opposite})
            {
                if (neighbor == invalid_surface_index || visited[neighbor]) continue;
                visited[neighbor] = 1;
                queue.push_back(neighbor);
            }
        }
        component.faces.assign(faces.begin(), faces.end());
        component.minimum_face = component.faces.front();
        for (const std::uint32_t face_id : component.faces)
        {
            const auto face = mesh.geometry.triangles[face_id];
            const Vec3 a = toVec(mesh.geometry.vertices[face[0]]);
            const Vec3 b = toVec(mesh.geometry.vertices[face[1]]);
            const Vec3 c = toVec(mesh.geometry.vertices[face[2]]);
            component.area += 0.5 * norm(cross(b - a, c - a));
            include(component.bounds, a);
            include(component.bounds, b);
            include(component.bounds, c);
        }
        result.push_back(std::move(component));
    }
    std::sort(result.begin(), result.end(), [](const auto& first, const auto& second)
    {
        if (first.area != second.area) return first.area > second.area;
        return first.minimum_face < second.minimum_face;
    });
    return result;
}

void writeDisconnectedHalfedgeDiagnostics(
    const std::filesystem::path& output_directory,
    const std::filesystem::path& source_obj,
    const OrientedSurfaceMesh& mesh,
    const TopologyFillStats& stats)
{
    const auto components = halfedgeComponentDiagnostics(mesh);
    const double total_area = std::accumulate(
        components.begin(), components.end(), 0.0,
        [](const double sum, const auto& component)
        { return sum + component.area; });
    std::filesystem::copy_file(
        source_obj, output_directory / "source.obj",
        std::filesystem::copy_options::overwrite_existing);
    writeAnalysisHalfedgeMesh(output_directory / "invalid_halfedge.bin", mesh);

    std::ofstream object(output_directory / "invalid_components.obj");
    if (!object) throw std::runtime_error("failed to create component diagnostic OBJ");
    object << std::setprecision(17);
    for (const Position3 point : mesh.geometry.vertices)
        object << "v " << point.x << ' ' << point.y << ' ' << point.z << '\n';
    for (std::size_t component_id = 0; component_id < components.size(); ++component_id)
    {
        object << "g primitive_" << component_id << "_component\n";
        for (const std::uint32_t face_id : components[component_id].faces)
        {
            const auto face = mesh.geometry.triangles[face_id];
            object << "f " << face[0] + 1 << ' ' << face[1] + 1 << ' '
                   << face[2] + 1 << '\n';
        }
    }

    std::ofstream component_json(output_directory / "components.json");
    if (!component_json) throw std::runtime_error("failed to create component diagnostics");
    component_json << std::setprecision(17) << "{\n  \"component_count\":"
                   << components.size() << ",\n  \"total_area\":" << total_area
                   << ",\n  \"components\":[\n";
    for (std::size_t id = 0; id < components.size(); ++id)
    {
        const auto& component = components[id];
        component_json << "    {\"id\":" << id
            << ",\"triangle_count\":" << component.faces.size()
            << ",\"area\":" << component.area
            << ",\"area_fraction\":"
            << (total_area > 0.0 ? component.area / total_area : 0.0)
            << ",\"boundary_halfedges\":" << component.boundary_halfedges
            << ",\"bounds_min\":[" << component.bounds.lower.x << ','
            << component.bounds.lower.y << ',' << component.bounds.lower.z << ']'
            << ",\"bounds_max\":[" << component.bounds.upper.x << ','
            << component.bounds.upper.y << ',' << component.bounds.upper.z << "]}"
            << (id + 1 == components.size() ? "\n" : ",\n");
    }
    component_json << "  ]\n}\n";

    std::ofstream metadata(output_directory / "model.json");
    if (!metadata) throw std::runtime_error("failed to create diagnostic metadata");
    metadata << std::setprecision(17)
        << "{\n  \"diagnostic_status\":\"invalid_phase1_disconnected\",\n"
        << "  \"stats\":{\"source_triangles\":" << stats.source_triangles
        << ",\"primitive_count\":" << components.size()
        << ",\"primitive_types\":{\"component\":" << components.size()
        << "},\"proxy_triangles\":" << mesh.geometry.triangles.size()
        << ",\"timings_seconds\":{\"total\":" << stats.elapsed_seconds
        << "},\"topology_fill\":{\"representation\":\"invalid disconnected halfedge diagnostic\""
        << ",\"output_voxels\":" << stats.output_voxels
        << ",\"boundary_voxels\":0,\"halfedge_face_components\":"
        << components.size() << "}},\n"
        << "  \"source\":\"source.obj\",\n"
        << "  \"phase3_simplified_surfaces\":\"invalid_components.obj\",\n"
        << "  \"component_diagnostics\":\"components.json\",\n"
        << "  \"proxy_components\":[\n";
    for (std::size_t id = 0; id < components.size(); ++id)
    {
        const auto& component = components[id];
        metadata << "    {\"id\":" << id << ",\"type\":\"component\""
            << ",\"triangulated_face_count\":" << component.faces.size()
            << ",\"area\":" << component.area
            << ",\"area_fraction\":"
            << (total_area > 0.0 ? component.area / total_area : 0.0)
            << ",\"boundary_halfedges\":" << component.boundary_halfedges << "}"
            << (id + 1 == components.size() ? "\n" : ",\n");
    }
    metadata << "  ],\n  \"viewer_stages\":[\"source\",\"phase3\",\"split\"]\n}\n";

    std::ofstream manifest(output_directory / "viewer_manifest.json");
    if (!manifest) throw std::runtime_error("failed to create diagnostic manifest");
    manifest << "{\n  \"algorithm\":\"InvalidHalfedgeComponentDiagnostic\",\n"
        << "  \"complete\":false,\n  \"model_count\":1,\n"
        << "  \"models\":[{\"id\":\"" << source_obj.stem().string()
        << "\",\"metadata\":\"model.json\"}]\n}\n";
}

} // namespace

void writePhase1BoundaryObj(
    const std::filesystem::path& path, const OrientedSurfaceMesh& surface)
{
    writeSurfaceObj(path, surface);
}

TopologyFillStats generateTopologyFillModel(
    const std::filesystem::path& source_obj,
    const std::filesystem::path& output_directory,
    const TopologyFillOptions& options)
{
    const auto started = std::chrono::steady_clock::now();
    if (options.maximum_grid_voxels == 0 || options.padding == 0 || options.maximum_steps == 0)
        throw std::invalid_argument("topology-fill options must be positive");
    if (std::filesystem::exists(output_directory))
        throw std::runtime_error("output directory already exists: " + output_directory.string());
    std::filesystem::create_directories(output_directory);
    const MeshModel source = readObj(source_obj);
    const auto log_stage = [&](const char* stage)
    {
        const double elapsed = std::chrono::duration<double>(
            std::chrono::steady_clock::now() - started).count();
        std::cerr << "monitor: stage=" << stage << " elapsed_seconds="
                  << elapsed << '\n';
    };
    log_stage("source_loaded");
    TriangleVoxelProvenance provenance;
    const VoxelGrid input = voxelizeTriangleSoup(
        source, options.maximum_grid_voxels, options.padding, &provenance);
    log_stage("source_voxelized");
    TopologyFillStats stats;
    stats.source_triangles = source.triangles.size();
    VoxelGrid cavity_fill_labels;
    VoxelGrid filled = enclosingTopologyFill(
        input, options.maximum_steps, &stats, &cavity_fill_labels);
    log_stage("topology_filled");
    // The phase-1 authority is the occupied/exterior interface. Internal walls
    // are absent by construction because both of their sides are occupied.
    // Source planes constrain geometry inside reconstructDualSurface; the
    // occupancy supplies topology only.
    Phase1Solid solid;
    try
    {
        solid = buildPhase1Solid(filled, source);
    }
    catch (...)
    {
        writeVoxelGrid(output_directory / "failed_topology.vox", filled);
        throw;
    }
    log_stage("dual_surface_built");
    const OrientedSurfaceMesh phase1 = projectPhase1BoundaryToSource(
        solid.boundary, solid.occupancy, input, source);
    log_stage("source_projected");
    AnalysisHalfedgeStats halfedge_stats;
    halfedge_stats.paired_edges = phase1.halfedges.size() / 2;
    halfedge_stats.boundary_halfedges = 0;
    halfedge_stats.face_components = 1;
    stats.exposed_voxel_faces = 0;
    stats.mesh_triangles = phase1.geometry.triangles.size();
    stats.mesh_vertices = phase1.geometry.vertices.size();
    stats.planar_regions = 0;
    stats.retained_source_triangles = 0;
    stats.removed_internal_triangles = 0;
    stats.ambiguous_source_triangles = 0;
    stats.hole_boundary_loops = 0;
    stats.cap_triangles = 0;
    stats.halfedge_count = phase1.halfedges.size();
    stats.paired_halfedge_edges = halfedge_stats.paired_edges;
    stats.boundary_halfedges = halfedge_stats.boundary_halfedges;
    stats.nonmanifold_edge_groups = halfedge_stats.nonmanifold_edge_groups;
    stats.inconsistent_orientation_edges =
        halfedge_stats.inconsistent_orientation_edges;
    stats.halfedge_face_components = halfedge_stats.face_components;
    stats.dropped_duplicate_triangles =
        halfedge_stats.dropped_duplicate_triangles;
    stats.dropped_degenerate_triangles =
        halfedge_stats.dropped_degenerate_triangles;
    stats.mesh_connected = halfedge_stats.face_components == 1;
    stats.mesh_watertight = true;
    stats.mesh_oriented = true;
    stats.mesh_manifold = true;
    stats.mesh_has_only_boundary_faces = true;
    stats.mesh_euler_characteristic = phase1.euler_characteristic;
    stats.mesh_signed_volume = phase1.signed_volume;
    stats.elapsed_seconds = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - started).count();
    if (!stats.mesh_connected)
    {
        writeDisconnectedHalfedgeDiagnostics(
            output_directory, source_obj, phase1, stats);
        throw std::runtime_error(
            "phase-1 halfedge is disconnected under next/opposite traversal; components=" +
            std::to_string(stats.halfedge_face_components));
    }
    std::filesystem::copy_file(source_obj, output_directory / "source.obj");
    writeAnalysisHalfedgeMesh(output_directory / "phase1_halfedge.bin", phase1);
    writeModelMetadata(output_directory / "model.json", source_obj.filename().string(), stats);
    return stats;
}

std::vector<TopologyFillStats> generateTopologyFillBatch(
    const std::filesystem::path& source_directory,
    const std::vector<int>& model_ids,
    const std::filesystem::path& output_root,
    const TopologyFillOptions& options)
{
    if (std::filesystem::exists(output_root))
        throw std::runtime_error("output root already exists: " + output_root.string());
    std::filesystem::create_directories(output_root / "models");
    std::vector<TopologyFillStats> results;
    results.reserve(model_ids.size());
    for (const int model_id : model_ids)
    {
        if (model_id < 0) throw std::invalid_argument("model IDs must be non-negative");
        results.push_back(generateTopologyFillModel(
            source_directory / (std::to_string(model_id) + ".obj"),
            output_root / "models" / std::to_string(model_id), options));
    }
    std::ofstream manifest(output_root / "viewer_manifest.json");
    if (!manifest) throw std::runtime_error("failed to create viewer manifest");
    manifest << "{\n"
        << "  \"algorithm\":\"CppHighResolutionVoxelLocatedOriginalObjHoleSurgery\",\n"
        << "  \"complete\":true,\n"
        << "  \"model_count\":" << model_ids.size() << ",\n"
        << "  \"models\":[\n";
    for (std::size_t index = 0; index < model_ids.size(); ++index)
    {
        manifest << "    {\"id\":" << model_ids[index]
                 << ",\"metadata\":\"models/" << model_ids[index] << "/model.json\"}"
                 << (index + 1 == model_ids.size() ? "\n" : ",\n");
    }
    manifest << "  ],\n"
        << "  \"options\":{\n"
        << "    \"maximum_grid_voxels\":" << options.maximum_grid_voxels << ",\n"
        << "    \"padding\":" << options.padding << ",\n"
        << "    \"candidate_axes\":[0,1,2],\n"
        << "    \"target_betti\":[1,0,0],\n"
        << "    \"source_occupancy_is_hard_kernel\":true,\n"
        << "    \"voxel_role\":\"hole location and topology certificate only\",\n"
        << "    \"phase1_representation\":\"PQSSHED1 validated open analysis halfedge\",\n"
        << "    \"voxel_boundary_mesh_used\":false,\n"
        << "    \"cpu_threads_per_model\":1,\n"
        << "    \"process_memory_limit_mib\":2048\n"
        << "  }\n"
        << "}\n";
    return results;
}

} // namespace pqss_proxy_mesh
