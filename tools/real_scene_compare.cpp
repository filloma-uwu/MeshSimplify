#include <algorithm>
#include <array>
#include <barrier>
#include <cctype>
#include <charconv>
#include <chrono>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <latch>
#include <limits>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <thread>
#include <vector>

#include "pqss/model_pool.hpp"
#include "pqss_proxy_mesh/proxy_bvh_build.hpp"

namespace
{

using Clock = std::chrono::steady_clock;

constexpr std::array<int, 13> kAppendedModelIds = {24, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37};
constexpr std::size_t kModelCount = 36;

struct QueryPair
{
    std::size_t pool_id_1{};
    std::size_t pool_id_2{};
    pqss::Mat3 R1{};
    pqss::Vec3 T1{};
    pqss::Mat3 R2{};
    pqss::Vec3 T2{};
    pqss::Real tolerance{};
};

struct AdditionTiming
{
    int model_id = 0;
    double build_ms = 0.0;
};

struct PoolBuildResult
{
    std::unique_ptr<pqss::ModelPool> pool;
    std::string cache_tag;
    double base_load_ms = 0.0;
    double base_end_pool_ms = 0.0;
    double base_total_ms = 0.0;
    double total_with_appends_ms = 0.0;
    std::vector<AdditionTiming> additions;
    std::size_t memory_bytes = 0;
    pqss::Real area_threshold = 0.0;
    std::vector<pqss::Real> model_area_thresholds;
    std::vector<pqss::ModelPool::ModelStats> model_stats;
    std::vector<pqss::ModelPool::BvhSignature> bvh_signatures;
    pqss::ModelPool::BuildCacheStats cache_stats;
};

struct ConfusionMatrix
{
    std::size_t true_positive = 0;
    std::size_t true_negative = 0;
    std::size_t false_positive = 0;
    std::size_t false_negative = 0;
};

struct QueryMetrics
{
    struct Pair
    {
        std::size_t queries = 0;
        std::size_t bv_tests = 0;
        std::size_t tri_tests = 0;
        std::size_t positives = 0;
        ConfusionMatrix confusion;
    };

    double parallel_wall_ms = 0.0;
    double summed_call_ms = 0.0;
    std::size_t total_bv_tests = 0;
    std::size_t total_tri_tests = 0;
    std::size_t positives = 0;
    std::array<Pair, kModelCount * kModelCount> pairs{};
};

struct ModelAggregate
{
    std::size_t input_triangles = 0;
    std::size_t built_triangles = 0;
    std::size_t subdivisions = 0;
    std::size_t bvs = 0;
    std::size_t leaves = 0;
    std::size_t internal_nodes = 0;
    std::size_t max_depth = 0;
    pqss::Real root_size_sum = 0.0;
    pqss::Real total_size = 0.0;
    pqss::Real leaf_size_sum = 0.0;
    pqss::Real internal_size_sum = 0.0;
};

std::filesystem::path ModelPath(const std::filesystem::path& directory, const int id)
{
    return directory / (std::to_string(id) + ".obj");
}

std::size_t PoolIndex(const std::size_t object_id)
{
    if (object_id >= 1 && object_id <= 23)
        return object_id - 1;
    if (object_id == 24)
        return 23;
    if (object_id >= 26 && object_id <= 28)
        return 24 + (object_id - 26);
    if (object_id >= 29 && object_id <= 32)
        return 27 + (object_id - 29);
    if (object_id >= 33 && object_id <= 37)
        return 31 + (object_id - 33);
    return std::numeric_limits<std::size_t>::max();
}

std::size_t ObjectIdFromName(const std::string_view name)
{
    const std::size_t dot = name.find('.');
    if (dot == std::string_view::npos)
        return std::numeric_limits<std::size_t>::max();
    std::size_t id = 0;
    const auto [end, error] = std::from_chars(name.data(), name.data() + dot, id);
    if (error != std::errc{} || end != name.data() + dot)
        return std::numeric_limits<std::size_t>::max();
    return id;
}

void SkipSpaces(const char*& cursor, const char* end)
{
    while (cursor < end && std::isspace(static_cast<unsigned char>(*cursor)))
        ++cursor;
}

std::string_view ReadToken(const char*& cursor, const char* end)
{
    SkipSpaces(cursor, end);
    const char* begin = cursor;
    while (cursor < end && !std::isspace(static_cast<unsigned char>(*cursor)))
        ++cursor;
    return {begin, cursor};
}

template <typename Value>
bool ParseNumber(const char*& cursor, const char* end, Value& value)
{
    SkipSpaces(cursor, end);
    const auto [next, error] = std::from_chars(cursor, end, value);
    if (error != std::errc{})
        return false;
    cursor = next;
    return true;
}

bool ParseRotation(const char*& cursor, const char* end, pqss::Mat3& rotation)
{
    for (int row = 0; row < 3; ++row)
        for (int column = 0; column < 3; ++column)
            if (!ParseNumber(cursor, end, rotation[row][column]))
                return false;
    return true;
}

bool ParseTranslation(const char*& cursor, const char* end, pqss::Vec3& translation)
{
    for (int index = 0; index < 3; ++index)
        if (!ParseNumber(cursor, end, translation[index]))
            return false;
    return true;
}

std::vector<QueryPair> LoadQueryPairs(const std::filesystem::path& filename)
{
    std::ifstream stream(filename);
    if (!stream)
        throw std::runtime_error("cannot open query file: " + filename.string());

    std::vector<QueryPair> pairs;
    pairs.reserve(900'000);
    std::string line;
    std::size_t line_number = 0;
    while (std::getline(stream, line))
    {
        ++line_number;
        if (line.empty())
            continue;
        const char* cursor = line.data();
        const char* end = cursor + line.size();
        const std::size_t id1 = ObjectIdFromName(ReadToken(cursor, end));
        QueryPair pair;
        if (!ParseRotation(cursor, end, pair.R1) || !ParseTranslation(cursor, end, pair.T1))
            throw std::runtime_error("invalid first pose at query line " + std::to_string(line_number));
        const std::size_t id2 = ObjectIdFromName(ReadToken(cursor, end));
        if (!ParseRotation(cursor, end, pair.R2) || !ParseTranslation(cursor, end, pair.T2) ||
            !ParseNumber(cursor, end, pair.tolerance))
            throw std::runtime_error("invalid second pose at query line " + std::to_string(line_number));
        pair.pool_id_1 = PoolIndex(id1);
        pair.pool_id_2 = PoolIndex(id2);
        if (pair.pool_id_1 == std::numeric_limits<std::size_t>::max() ||
            pair.pool_id_2 == std::numeric_limits<std::size_t>::max())
            throw std::runtime_error("unknown model id at query line " + std::to_string(line_number));
        pairs.push_back(pair);
    }
    return pairs;
}

double ElapsedMs(const Clock::time_point begin, const Clock::time_point end)
{
    return std::chrono::duration<double, std::milli>(end - begin).count();
}

void RequireSuccess(const pqss::Return result, const std::string& operation, const pqss::ModelPool& pool)
{
    static_cast<void>(pool);
    if (result != pqss::Return::Success)
        throw std::runtime_error(operation + " failed with Return=" +
                                 std::to_string(static_cast<int>(result)));
}

PoolBuildResult BuildPool(const std::filesystem::path& base_directory,
                          const std::filesystem::path& workpiece_directory,
                          const std::string& cache_tag,
                          const bool use_proxy_bvh,
                          const pqss::Real area_threshold_override = pqss::k_zero)
{
    pqss_proxy_mesh::setProxyBvhBuildEnabled(use_proxy_bvh);
    PoolBuildResult result;
    result.pool = std::make_unique<pqss::ModelPool>();
    result.cache_tag = cache_tag;
    RequireSuccess(result.pool->SetBuildCacheTag(cache_tag), "SetBuildCacheTag", *result.pool);
    RequireSuccess(result.pool->SetAreaThresholdOverride(area_threshold_override),
                   "SetAreaThresholdOverride", *result.pool);

    const auto base_begin = Clock::now();
    const auto load_begin = Clock::now();
    for (int id = 1; id <= 23; ++id)
    {
        const auto path = ModelPath(base_directory, id);
        RequireSuccess(result.pool->LoadModelFromFile(path.string(), pqss::BuildStrategy::Optimized),
                       "LoadModelFromFile(" + path.string() + ")", *result.pool);
    }
    const auto load_end = Clock::now();
    RequireSuccess(result.pool->EndPool(), "EndPool", *result.pool);
    const auto base_end = Clock::now();
    result.base_load_ms = ElapsedMs(load_begin, load_end);
    result.base_end_pool_ms = ElapsedMs(load_end, base_end);
    result.base_total_ms = ElapsedMs(base_begin, base_end);

    for (const int id : kAppendedModelIds)
    {
        const auto path = ModelPath(workpiece_directory, id);
        const auto begin = Clock::now();
        RequireSuccess(result.pool->LoadModelFromFile(path.string(), pqss::BuildStrategy::Optimized),
                       "LoadModelFromFile(" + path.string() + ")", *result.pool);
        result.additions.push_back({id, ElapsedMs(begin, Clock::now())});
    }
    result.total_with_appends_ms = ElapsedMs(base_begin, Clock::now());
    result.memory_bytes = result.pool->Info(false);
    result.area_threshold = result.pool->MaxArea();
    result.model_area_thresholds = result.pool->ModelAreaThresholds();
    result.model_stats = result.pool->GetModelStats();
    result.bvh_signatures = result.pool->BvhSignatures();
    result.cache_stats = result.pool->GetBuildCacheStats();
    pqss_proxy_mesh::setProxyBvhBuildEnabled(false);
    return result;
}

struct WorkerResult
{
    double call_ms = 0.0;
    std::size_t bv_tests = 0;
    std::size_t tri_tests = 0;
    std::size_t positives = 0;
    ConfusionMatrix confusion;
    std::array<QueryMetrics::Pair, kModelCount * kModelCount> pairs{};
};

QueryMetrics RunQueries(const pqss::ModelPool& pool,
                        const std::vector<QueryPair>& pairs,
                        const std::size_t thread_count,
                        std::vector<std::uint8_t>* decisions,
                        const std::vector<std::uint8_t>* baseline,
                        ConfusionMatrix* confusion)
{
    const std::size_t active_threads = std::min(thread_count, pairs.size());
    std::vector<WorkerResult> worker_results(active_threads);
    std::vector<std::thread> workers;
    workers.reserve(active_threads);
    std::latch workers_ready(active_threads);
    std::latch workers_finished(active_threads);
    std::barrier phase_start(static_cast<std::ptrdiff_t>(active_threads + 1));

    for (std::size_t thread_index = 0; thread_index < active_threads; ++thread_index)
    {
        const std::size_t begin = pairs.size() * thread_index / active_threads;
        const std::size_t end = pairs.size() * (thread_index + 1) / active_threads;
        workers.emplace_back([&, thread_index, begin, end]()
        {
            WorkerResult& worker = worker_results[thread_index];
            workers_ready.count_down();
            phase_start.arrive_and_wait();
            const auto call_begin = Clock::now();
            for (std::size_t index = begin; index < end; ++index)
            {
                const QueryPair& pair = pairs[index];
                pqss::QueryResult query_result;
                const pqss::Return status = pool.Query(query_result,
                                                       pair.pool_id_1, pair.R1, pair.T1,
                                                       pair.pool_id_2, pair.R2, pair.T2,
                                                       pair.tolerance);
                if (status != pqss::Return::Success)
                    throw std::runtime_error("PQSS query failed");
                const bool positive = query_result.CloserThanTolerance();
                const std::size_t pair_index = pair.pool_id_1 * kModelCount + pair.pool_id_2;
                QueryMetrics::Pair& pair_metrics = worker.pairs[pair_index];
                ++pair_metrics.queries;
                pair_metrics.bv_tests += query_result.NumBvTests();
                pair_metrics.tri_tests += query_result.NumTriTests();
                pair_metrics.positives += positive ? 1u : 0u;
                if (decisions != nullptr)
                    (*decisions)[index] = positive ? 1u : 0u;
                if (baseline != nullptr)
                {
                    const bool expected = (*baseline)[index] != 0;
                    if (positive && expected) ++worker.confusion.true_positive;
                    else if (!positive && !expected) ++worker.confusion.true_negative;
                    else if (positive) ++worker.confusion.false_positive;
                    else ++worker.confusion.false_negative;
                    if (positive && expected) ++pair_metrics.confusion.true_positive;
                    else if (!positive && !expected) ++pair_metrics.confusion.true_negative;
                    else if (positive) ++pair_metrics.confusion.false_positive;
                    else ++pair_metrics.confusion.false_negative;
                }
                worker.positives += positive ? 1u : 0u;
                worker.bv_tests += query_result.NumBvTests();
                worker.tri_tests += query_result.NumTriTests();
            }
            worker.call_ms = ElapsedMs(call_begin, Clock::now());
            workers_finished.count_down();
        });
    }

    workers_ready.wait();
    const auto wall_begin = Clock::now();
    phase_start.arrive_and_wait();
    workers_finished.wait();
    const auto wall_end = Clock::now();
    for (std::thread& worker : workers)
        worker.join();

    QueryMetrics metrics;
    metrics.parallel_wall_ms = ElapsedMs(wall_begin, wall_end);
    for (const WorkerResult& worker : worker_results)
    {
        metrics.summed_call_ms += worker.call_ms;
        metrics.total_bv_tests += worker.bv_tests;
        metrics.total_tri_tests += worker.tri_tests;
        metrics.positives += worker.positives;
        for (std::size_t index = 0; index < metrics.pairs.size(); ++index)
        {
            metrics.pairs[index].queries += worker.pairs[index].queries;
            metrics.pairs[index].bv_tests += worker.pairs[index].bv_tests;
            metrics.pairs[index].tri_tests += worker.pairs[index].tri_tests;
            metrics.pairs[index].positives += worker.pairs[index].positives;
            metrics.pairs[index].confusion.true_positive += worker.pairs[index].confusion.true_positive;
            metrics.pairs[index].confusion.true_negative += worker.pairs[index].confusion.true_negative;
            metrics.pairs[index].confusion.false_positive += worker.pairs[index].confusion.false_positive;
            metrics.pairs[index].confusion.false_negative += worker.pairs[index].confusion.false_negative;
        }
        if (confusion != nullptr)
        {
            confusion->true_positive += worker.confusion.true_positive;
            confusion->true_negative += worker.confusion.true_negative;
            confusion->false_positive += worker.confusion.false_positive;
            confusion->false_negative += worker.confusion.false_negative;
        }
    }
    return metrics;
}

std::vector<int> PoolObjectIds()
{
    std::vector<int> result;
    for (int id = 1; id <= 23; ++id)
        result.push_back(id);
    result.insert(result.end(), kAppendedModelIds.begin(), kAppendedModelIds.end());
    return result;
}

ModelAggregate AggregateBaseModels(const PoolBuildResult& pool)
{
    ModelAggregate result;
    for (std::size_t index = 0; index < 23; ++index)
    {
        const auto& model = pool.model_stats.at(index);
        const auto& bvh = pool.bvh_signatures.at(index);
        result.input_triangles += model.original_tri_count;
        result.built_triangles += model.built_tri_count;
        result.subdivisions += model.subdivision_count;
        result.bvs += model.bv_count;
        result.leaves += bvh.leaf_count;
        result.internal_nodes += bvh.internal_count;
        result.max_depth = std::max(result.max_depth, bvh.max_depth);
        result.root_size_sum += bvh.root_size;
        result.total_size += bvh.total_size;
        result.leaf_size_sum += bvh.leaf_size_sum;
        result.internal_size_sum += bvh.internal_size_sum;
    }
    return result;
}

void WriteModelCsv(const std::filesystem::path& path,
                   const PoolBuildResult& original,
                   const PoolBuildResult& proxy)
{
    std::ofstream stream(path);
    if (!stream)
        throw std::runtime_error("cannot write " + path.string());
    stream << "model_id,category,original_triangles,proxy_input_triangles,original_built_triangles,"
              "proxy_built_triangles,original_subdivisions,proxy_subdivisions,original_bv_count,proxy_bv_count,"
              "original_leaf_count,proxy_leaf_count,original_internal_count,proxy_internal_count,"
              "original_max_depth,proxy_max_depth,original_root_size,proxy_root_size,original_total_size,"
              "proxy_total_size,original_leaf_size_sum,proxy_leaf_size_sum,original_internal_size_sum,"
              "proxy_internal_size_sum\n";
    const std::vector<int> object_ids = PoolObjectIds();
    for (std::size_t index = 0; index < object_ids.size(); ++index)
    {
        const auto& om = original.model_stats.at(index);
        const auto& pm = proxy.model_stats.at(index);
        const auto& ob = original.bvh_signatures.at(index);
        const auto& pb = proxy.bvh_signatures.at(index);
        stream << object_ids[index] << ',' << (object_ids[index] <= 23 ? "base" : "workpiece") << ','
               << om.original_tri_count << ',' << pm.original_tri_count << ','
               << om.built_tri_count << ',' << pm.built_tri_count << ','
               << om.subdivision_count << ',' << pm.subdivision_count << ','
               << om.bv_count << ',' << pm.bv_count << ','
               << ob.leaf_count << ',' << pb.leaf_count << ','
               << ob.internal_count << ',' << pb.internal_count << ','
               << ob.max_depth << ',' << pb.max_depth << ','
               << ob.root_size << ',' << pb.root_size << ','
               << ob.total_size << ',' << pb.total_size << ','
               << ob.leaf_size_sum << ',' << pb.leaf_size_sum << ','
               << ob.internal_size_sum << ',' << pb.internal_size_sum << '\n';
    }
}

void WritePairCsv(const std::filesystem::path& path,
                  const QueryMetrics& original,
                  const QueryMetrics& proxy)
{
    std::ofstream stream(path);
    if (!stream)
        throw std::runtime_error("cannot write " + path.string());
    stream << "model_id_1,model_id_2,queries,original_bv_tests,proxy_bv_tests,original_tri_tests,"
              "proxy_tri_tests,original_positives,proxy_positives,tp,tn,fp,fn\n";
    const std::vector<int> object_ids = PoolObjectIds();
    for (std::size_t first = 0; first < kModelCount; ++first)
    {
        for (std::size_t second = 0; second < kModelCount; ++second)
        {
            const std::size_t index = first * kModelCount + second;
            if (original.pairs[index].queries == 0 && proxy.pairs[index].queries == 0)
                continue;
            const auto& baseline = original.pairs[index];
            const auto& candidate = proxy.pairs[index];
            stream << object_ids[first] << ',' << object_ids[second] << ',' << baseline.queries << ','
                   << baseline.bv_tests << ',' << candidate.bv_tests << ','
                   << baseline.tri_tests << ',' << candidate.tri_tests << ','
                   << baseline.positives << ',' << candidate.positives << ','
                   << candidate.confusion.true_positive << ',' << candidate.confusion.true_negative << ','
                   << candidate.confusion.false_positive << ',' << candidate.confusion.false_negative << '\n';
        }
    }
}

void WriteAdditionTable(std::ostream& stream, const PoolBuildResult& original, const PoolBuildResult& proxy)
{
    stream << "| Model | Original append ms | Proxy-pool append ms |\n"
              "|---:|---:|---:|\n";
    for (std::size_t index = 0; index < original.additions.size(); ++index)
        stream << '|' << original.additions[index].model_id << '|'
               << original.additions[index].build_ms << '|'
               << proxy.additions[index].build_ms << "|\n";
}

void WriteReport(const std::filesystem::path& path,
                 const std::size_t query_count,
                 const std::size_t thread_count,
                 const PoolBuildResult& original,
                 const PoolBuildResult& proxy,
                 const QueryMetrics& original_query,
                 const QueryMetrics& proxy_query,
                 const ConfusionMatrix& confusion,
                 const std::size_t requested_max_depth,
                 const bool use_proxy_bvh)
{
    std::ofstream stream(path);
    if (!stream)
        throw std::runtime_error("cannot write " + path.string());
    const std::size_t negatives = confusion.false_positive + confusion.true_negative;
    const double false_positive_rate = negatives == 0 ? 0.0 :
        static_cast<double>(confusion.false_positive) / static_cast<double>(negatives);
    const ModelAggregate original_models = AggregateBaseModels(original);
    const ModelAggregate proxy_models = AggregateBaseModels(proxy);
    stream << std::fixed << std::setprecision(6);
    stream << "# PQSS Real Scene Conservative Proxy Comparison\n\n"
           << "- Build strategy: `Optimized`\n"
           << "- BVH builder: `" << (use_proxy_bvh ? "project six-axis" : "PQSS standard")
           << "`\n"
           << "- Threads: " << thread_count << "\n"
           << "- Queries: " << query_count << "\n"
           << "- Base models: proxy `1-23` versus original `1-23`\n"
           << "- Appended workpieces: original geometry in both pools\n\n"
           << "- Requested proxy max depth: " << requested_max_depth << "\n"
           << "- Actual proxy max depth: " << proxy_models.max_depth << "\n"
           << "- Depth constraint satisfied: `"
           << (proxy_models.max_depth <= requested_max_depth ? "true" : "false") << "`\n\n"
           << "## Build\n\n"
           << "| Metric | Original pool | Proxy pool |\n|---|---:|---:|\n"
           << "| Base OBJ load ms |" << original.base_load_ms << '|' << proxy.base_load_ms << "|\n"
           << "| Base EndPool ms |" << original.base_end_pool_ms << '|' << proxy.base_end_pool_ms << "|\n"
           << "| Base total ms |" << original.base_total_ms << '|' << proxy.base_total_ms << "|\n"
           << "| Total including appends ms |" << original.total_with_appends_ms << '|'
           << proxy.total_with_appends_ms << "|\n"
           << "| Pool memory bytes |" << original.memory_bytes << '|' << proxy.memory_bytes << "|\n"
           << "| Area threshold |" << original.area_threshold << '|' << proxy.area_threshold << "|\n"
           << "| Cache hits |" << original.cache_stats.cache_hits << '|' << proxy.cache_stats.cache_hits << "|\n"
           << "| Cache misses |" << original.cache_stats.cache_misses << '|' << proxy.cache_stats.cache_misses << "|\n\n"
           << "- Original cache tag: `" << original.cache_tag << "`\n"
           << "- Proxy cache tag: `" << proxy.cache_tag << "`\n\n"
           << "### Cache Phases\n\n"
           << "| Metric | Original pool | Proxy pool |\n|---|---:|---:|\n"
           << "| Eligible models |" << original.cache_stats.eligible_models << '|'
           << proxy.cache_stats.eligible_models << "|\n"
           << "| Area threshold ms |" << original.cache_stats.area_threshold_ms << '|'
           << proxy.cache_stats.area_threshold_ms << "|\n"
           << "| Cache load ms |" << original.cache_stats.cache_load_ms << '|'
           << proxy.cache_stats.cache_load_ms << "|\n"
           << "| Cache build ms |" << original.cache_stats.cache_build_ms << '|'
           << proxy.cache_stats.cache_build_ms << "|\n"
           << "| Cache save ms |" << original.cache_stats.cache_save_ms << '|'
           << proxy.cache_stats.cache_save_ms << "|\n\n"
           << "### Appended Workpieces\n\n";
    WriteAdditionTable(stream, original, proxy);
    stream << "\n## Base Model And BVH Data\n\n"
           << "| Metric | Original `1-23` | Proxy `1-23` |\n|---|---:|---:|\n"
           << "| Input triangles |" << original_models.input_triangles << '|' << proxy_models.input_triangles << "|\n"
           << "| Built triangles |" << original_models.built_triangles << '|' << proxy_models.built_triangles << "|\n"
           << "| Subdivisions |" << original_models.subdivisions << '|' << proxy_models.subdivisions << "|\n"
           << "| BV nodes |" << original_models.bvs << '|' << proxy_models.bvs << "|\n"
           << "| Leaf nodes |" << original_models.leaves << '|' << proxy_models.leaves << "|\n"
           << "| Internal nodes |" << original_models.internal_nodes << '|' << proxy_models.internal_nodes << "|\n"
           << "| Maximum depth |" << original_models.max_depth << '|' << proxy_models.max_depth << "|\n"
           << "| Root size sum |" << original_models.root_size_sum << '|' << proxy_models.root_size_sum << "|\n"
           << "| Total size |" << original_models.total_size << '|' << proxy_models.total_size << "|\n"
           << "| Leaf size sum |" << original_models.leaf_size_sum << '|' << proxy_models.leaf_size_sum << "|\n"
           << "| Internal size sum |" << original_models.internal_size_sum << '|'
           << proxy_models.internal_size_sum << "|\n";
    stream << "\n## Queries\n\n"
           << "| Metric | Original pool | Proxy pool |\n|---|---:|---:|\n"
           << "| Parallel wall ms |" << original_query.parallel_wall_ms << '|' << proxy_query.parallel_wall_ms << "|\n"
           << "| Summed call ms |" << original_query.summed_call_ms << '|' << proxy_query.summed_call_ms << "|\n"
           << "| BV tests |" << original_query.total_bv_tests << '|' << proxy_query.total_bv_tests << "|\n"
           << "| Triangle tests |" << original_query.total_tri_tests << '|' << proxy_query.total_tri_tests << "|\n"
           << "| Positive results |" << original_query.positives << '|' << proxy_query.positives << "|\n\n"
           << "## Correctness\n\n"
           << "| TP | TN | FP | FN | False-positive rate |\n|---:|---:|---:|---:|---:|\n|"
           << confusion.true_positive << '|' << confusion.true_negative << '|'
           << confusion.false_positive << '|' << confusion.false_negative << '|'
           << false_positive_rate << "|\n\n"
           << "The false-positive denominator is `FP + TN`. A valid conservative result requires `FN = 0`.\n\n"
           << "Per-model geometry and BVH data are in `model_bvh_stats.csv`; query-pair data are in "
              "`query_pair_stats.csv`.\n";
}

std::string UniqueCacheTag(const std::string_view prefix)
{
    const auto ticks = Clock::now().time_since_epoch().count();
    return std::string(prefix) + '_' + std::to_string(ticks);
}

} // namespace

int main(int argc, char** argv)
{
    try
    {
        if (argc < 6 || argc > 11)
        {
            std::cerr << "Usage: pqss-real-scene-compare <original-base-dir> <proxy-base-dir> "
                         "<original-workpiece-dir> <collision-pairs.txt> <report-dir> [threads] "
                         "[original-cache-tag] [proxy-cache-tag] [requested-max-depth] "
                         "[standard|six-axis]\n";
            return 2;
        }
        const std::filesystem::path original_base(argv[1]);
        const std::filesystem::path proxy_base(argv[2]);
        const std::filesystem::path workpieces(argv[3]);
        const std::filesystem::path query_file(argv[4]);
        const std::filesystem::path report_directory(argv[5]);
        std::size_t thread_count = 1;
        if (argc >= 7)
        {
            const std::string_view text = argv[6];
            const auto [end, error] = std::from_chars(text.data(), text.data() + text.size(), thread_count);
            if (error != std::errc{} || end != text.data() + text.size() || thread_count == 0)
                throw std::runtime_error("threads must be a positive integer");
        }
        const std::string original_cache_tag = argc >= 8 ? argv[7] : UniqueCacheTag("original");
        const std::string proxy_cache_tag = argc >= 9 ? argv[8] : UniqueCacheTag("proxy");
        std::size_t requested_max_depth = 8;
        if (argc >= 10)
        {
            const std::string_view text = argv[9];
            const auto [end, error] = std::from_chars(text.data(), text.data() + text.size(), requested_max_depth);
            if (error != std::errc{} || end != text.data() + text.size() || requested_max_depth == 0)
                throw std::runtime_error("requested max depth must be a positive integer");
        }
        bool use_proxy_bvh = true;
        if (argc >= 11)
        {
            const std::string_view mode(argv[10]);
            if (mode == "standard") use_proxy_bvh = false;
            else if (mode == "six-axis") use_proxy_bvh = true;
            else throw std::runtime_error("BVH mode must be standard or six-axis");
        }
        std::cout << "Loading real-scene query poses...\n";
        const std::vector<QueryPair> pairs = LoadQueryPairs(query_file);
        std::cout << "queries=" << pairs.size() << "\nBuilding original Optimized pool...\n" << std::flush;
        PoolBuildResult original = BuildPool(original_base, workpieces, original_cache_tag, false);
        std::cout << "original_base_build_ms=" << original.base_total_ms
                  << "\nBuilding proxy Optimized pool...\n" << std::flush;
        PoolBuildResult proxy = BuildPool(
            proxy_base, workpieces, proxy_cache_tag, use_proxy_bvh,
            original.area_threshold);
        std::cout << "proxy_base_build_ms=" << proxy.base_total_ms << "\nRunning original queries...\n" << std::flush;

        std::vector<std::uint8_t> original_decisions(pairs.size());
        const QueryMetrics original_query =
            RunQueries(*original.pool, pairs, thread_count, &original_decisions, nullptr, nullptr);
        std::cout << "Running proxy queries...\n" << std::flush;
        ConfusionMatrix confusion;
        const QueryMetrics proxy_query =
            RunQueries(*proxy.pool, pairs, thread_count, nullptr, &original_decisions, &confusion);

        std::filesystem::create_directories(report_directory);
        WriteModelCsv(report_directory / "model_bvh_stats.csv", original, proxy);
        WritePairCsv(report_directory / "query_pair_stats.csv", original_query, proxy_query);
        WriteReport(report_directory / "comparison_report.md", pairs.size(), thread_count,
                    original, proxy, original_query, proxy_query, confusion,
                    requested_max_depth, use_proxy_bvh);

        std::cout << "original_query_wall_ms=" << original_query.parallel_wall_ms
                  << "\nproxy_query_wall_ms=" << proxy_query.parallel_wall_ms
                  << "\nfalse_positive=" << confusion.false_positive
                  << "\nfalse_negative=" << confusion.false_negative
                  << "\nreport=" << (report_directory / "comparison_report.md").string() << '\n';
        const ModelAggregate proxy_models = AggregateBaseModels(proxy);
        if (confusion.false_negative != 0)
            return 3;
        return proxy_models.max_depth <= requested_max_depth ? 0 : 4;
    }
    catch (const std::exception& error)
    {
        std::cerr << "error: " << error.what() << '\n';
        return 1;
    }
}
