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
#include <mutex>
#include <stdexcept>
#include <string>
#include <string_view>
#include <thread>
#include <vector>

#include "pqss/model_pool.hpp"
#include "pqss_proxy_mesh/independent_primitive_pool.hpp"
#include "pqss_proxy_mesh/proxy_bvh_build.hpp"

namespace
{

using Clock = std::chrono::steady_clock;
constexpr std::size_t kLogicalModelCount = 36;
constexpr std::array<int, 13> kWorkpieces = {24, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37};

struct QueryPair
{
    std::size_t first = 0;
    std::size_t second = 0;
    pqss::Mat3 R1{};
    pqss::Vec3 T1{};
    pqss::Mat3 R2{};
    pqss::Vec3 T2{};
    pqss::Real tolerance = pqss::k_zero;
};

struct FilterMetrics
{
    std::size_t skips = 0;
    std::size_t containment_checks = 0;
    std::size_t containment_skips = 0;
    std::size_t axis_checks = 0;
    std::size_t axis_skips = 0;
    std::size_t axis_generation_checks = 0;
    std::size_t cached_axes = 0;
};

void AddFilter(FilterMetrics& target, const FilterMetrics& source)
{
    target.skips += source.skips;
    target.containment_checks += source.containment_checks;
    target.containment_skips += source.containment_skips;
    target.axis_checks += source.axis_checks;
    target.axis_skips += source.axis_skips;
    target.axis_generation_checks += source.axis_generation_checks;
    target.cached_axes += source.cached_axes;
}

struct Sample
{
    bool positive = false;
    std::size_t root_tests = 0;
    std::size_t internal_tests = 0;
    std::size_t triangle_tests = 0;
    FilterMetrics filter;
    double call_ms = 0.0;
};

struct Confusion
{
    std::size_t tp = 0;
    std::size_t tn = 0;
    std::size_t fp = 0;
    std::size_t fn = 0;
};

struct PairMetrics
{
    std::size_t queries = 0;
    std::size_t positives = 0;
    std::size_t root_tests = 0;
    std::size_t internal_tests = 0;
    std::size_t triangle_tests = 0;
    FilterMetrics filter;
    Confusion confusion;
};

struct Metrics
{
    double wall_ms = 0.0;
    double summed_call_ms = 0.0;
    std::size_t positives = 0;
    std::size_t root_tests = 0;
    std::size_t internal_tests = 0;
    std::size_t triangle_tests = 0;
    FilterMetrics filter;
    Confusion confusion;
    std::array<PairMetrics, kLogicalModelCount * kLogicalModelCount> pairs{};
};

struct StandardBuild
{
    std::unique_ptr<pqss::ModelPool> pool;
    double base_ms = 0.0;
    double total_ms = 0.0;
    std::size_t memory_bytes = 0;
};

struct IndependentBuild
{
    std::unique_ptr<pqss_proxy_mesh::IndependentPrimitivePool> pool;
    double base_ms = 0.0;
    double containment_ms = 0.0;
    double total_ms = 0.0;
    std::size_t memory_bytes = 0;
};

std::size_t PoolIndex(const std::size_t object_id)
{
    if (object_id >= 1 && object_id <= 23) return object_id - 1;
    if (object_id == 24) return 23;
    if (object_id >= 26 && object_id <= 28) return 24 + object_id - 26;
    if (object_id >= 29 && object_id <= 32) return 27 + object_id - 29;
    if (object_id >= 33 && object_id <= 37) return 31 + object_id - 33;
    return std::numeric_limits<std::size_t>::max();
}

std::vector<int> ObjectIds()
{
    std::vector<int> result;
    for (int id = 1; id <= 23; ++id) result.push_back(id);
    result.insert(result.end(), kWorkpieces.begin(), kWorkpieces.end());
    return result;
}

std::filesystem::path ModelPath(const std::filesystem::path& directory, const int id)
{
    return directory / (std::to_string(id) + ".obj");
}

std::size_t ObjectIdFromToken(const std::string_view token)
{
    const std::size_t dot = token.find('.');
    if (dot == std::string_view::npos) return std::numeric_limits<std::size_t>::max();
    std::size_t result = 0;
    const auto parsed = std::from_chars(token.data(), token.data() + dot, result);
    return parsed.ec == std::errc{} && parsed.ptr == token.data() + dot
        ? result : std::numeric_limits<std::size_t>::max();
}

void SkipSpaces(const char*& cursor, const char* end)
{
    while (cursor < end && std::isspace(static_cast<unsigned char>(*cursor))) ++cursor;
}

std::string_view ReadToken(const char*& cursor, const char* end)
{
    SkipSpaces(cursor, end);
    const char* begin = cursor;
    while (cursor < end && !std::isspace(static_cast<unsigned char>(*cursor))) ++cursor;
    return {begin, cursor};
}

template <typename T>
bool ParseNumber(const char*& cursor, const char* end, T& value)
{
    SkipSpaces(cursor, end);
    const auto parsed = std::from_chars(cursor, end, value);
    if (parsed.ec != std::errc{}) return false;
    cursor = parsed.ptr;
    return true;
}

bool ParsePose(const char*& cursor, const char* end, pqss::Mat3& R, pqss::Vec3& T)
{
    for (int row = 0; row < 3; ++row)
        for (int column = 0; column < 3; ++column)
            if (!ParseNumber(cursor, end, R[row][column])) return false;
    for (int axis = 0; axis < 3; ++axis)
        if (!ParseNumber(cursor, end, T[axis])) return false;
    return true;
}

std::vector<QueryPair> LoadQueries(const std::filesystem::path& path)
{
    std::ifstream stream(path);
    if (!stream) throw std::runtime_error("cannot open query file: " + path.string());
    std::vector<QueryPair> result;
    result.reserve(750'000);
    std::string line;
    std::size_t line_number = 0;
    while (std::getline(stream, line))
    {
        ++line_number;
        if (line.empty()) continue;
        const char* cursor = line.data();
        const char* end = cursor + line.size();
        const std::size_t first_id = ObjectIdFromToken(ReadToken(cursor, end));
        QueryPair query;
        if (!ParsePose(cursor, end, query.R1, query.T1))
            throw std::runtime_error("invalid first pose on line " + std::to_string(line_number));
        const std::size_t second_id = ObjectIdFromToken(ReadToken(cursor, end));
        if (!ParsePose(cursor, end, query.R2, query.T2) ||
            !ParseNumber(cursor, end, query.tolerance))
            throw std::runtime_error("invalid second pose on line " + std::to_string(line_number));
        query.first = PoolIndex(first_id);
        query.second = PoolIndex(second_id);
        if (query.first >= kLogicalModelCount || query.second >= kLogicalModelCount)
            throw std::runtime_error("unknown model id on line " + std::to_string(line_number));
        result.push_back(query);
    }
    return result;
}

void Require(const pqss::Return result, const std::string& operation)
{
    if (result != pqss::Return::Success)
        throw std::runtime_error(operation + " failed with Return=" +
                                 std::to_string(static_cast<int>(result)));
}

double Milliseconds(const Clock::time_point begin, const Clock::time_point end)
{
    return std::chrono::duration<double, std::milli>(end - begin).count();
}

StandardBuild BuildStandard(const std::filesystem::path& base,
                            const std::filesystem::path& workpieces,
                            const std::string& cache_tag)
{
    pqss_proxy_mesh::setProxyBvhBuildEnabled(false);
    StandardBuild result;
    result.pool = std::make_unique<pqss::ModelPool>();
    Require(result.pool->SetBuildCacheTag(cache_tag), "SetBuildCacheTag");
    const auto started = Clock::now();
    for (int id = 1; id <= 23; ++id)
        Require(result.pool->LoadModelFromFile(ModelPath(base, id).string(), pqss::BuildStrategy::Optimized),
                "LoadModelFromFile");
    Require(result.pool->EndPool(), "EndPool");
    result.base_ms = Milliseconds(started, Clock::now());
    for (const int id : kWorkpieces)
        Require(result.pool->LoadModelFromFile(ModelPath(workpieces, id).string(), pqss::BuildStrategy::Optimized),
                "append workpiece");
    result.total_ms = Milliseconds(started, Clock::now());
    result.memory_bytes = result.pool->Info(false);
    return result;
}

IndependentBuild BuildIndependent(const std::filesystem::path& manifest,
                                  const std::filesystem::path& workpieces,
                                  const std::string& cache_tag)
{
    pqss_proxy_mesh::setProxyBvhBuildEnabled(true);
    IndependentBuild result;
    result.pool = std::make_unique<pqss_proxy_mesh::IndependentPrimitivePool>();
    result.pool->setBuildCacheTag(cache_tag);
    const auto started = Clock::now();
    result.pool->loadBaseManifest(manifest);
    result.pool->endBasePool();
    result.base_ms = Milliseconds(started, Clock::now());
    for (const int id : kWorkpieces)
        result.pool->appendLogicalModel(PoolIndex(id), ModelPath(workpieces, id));
    const auto containment_started = Clock::now();
    result.pool->finalizeContainment();
    result.containment_ms = Milliseconds(containment_started, Clock::now());
    result.total_ms = Milliseconds(started, Clock::now());
    result.memory_bytes = result.pool->physicalPool().Info(false);
    pqss_proxy_mesh::setProxyBvhBuildEnabled(false);
    return result;
}

void AddConfusion(Confusion& confusion, const bool actual, const bool expected)
{
    if (actual && expected) ++confusion.tp;
    else if (!actual && !expected) ++confusion.tn;
    else if (actual) ++confusion.fp;
    else ++confusion.fn;
}

template <typename QueryFunction>
Metrics RunQueries(const std::vector<QueryPair>& queries,
                   const std::size_t requested_threads,
                   QueryFunction query_function,
                   std::vector<std::uint8_t>* decisions,
                   const std::vector<std::uint8_t>* truth)
{
    const std::size_t thread_count = std::min(requested_threads, queries.size());
    std::vector<Metrics> workers(thread_count);
    std::vector<std::thread> threads;
    std::latch ready(thread_count);
    std::latch finished(thread_count);
    std::barrier start(static_cast<std::ptrdiff_t>(thread_count + 1));
    std::mutex error_mutex;
    std::exception_ptr error;
    for (std::size_t thread_index = 0; thread_index < thread_count; ++thread_index)
    {
        const std::size_t begin = queries.size() * thread_index / thread_count;
        const std::size_t end = queries.size() * (thread_index + 1) / thread_count;
        threads.emplace_back([&, thread_index, begin, end]
        {
            ready.count_down();
            start.arrive_and_wait();
            const auto call_started = Clock::now();
            try
            {
                Metrics& worker = workers[thread_index];
                for (std::size_t index = begin; index < end; ++index)
                {
                    const QueryPair& query = queries[index];
                    const Sample sample = query_function(query);
                    if (decisions != nullptr) (*decisions)[index] = sample.positive ? 1 : 0;
                    worker.positives += sample.positive;
                    worker.root_tests += sample.root_tests;
                    worker.internal_tests += sample.internal_tests;
                    worker.triangle_tests += sample.triangle_tests;
                    AddFilter(worker.filter, sample.filter);
                    worker.summed_call_ms += sample.call_ms;
                    PairMetrics& pair = worker.pairs[query.first * kLogicalModelCount + query.second];
                    ++pair.queries;
                    pair.positives += sample.positive;
                    pair.root_tests += sample.root_tests;
                    pair.internal_tests += sample.internal_tests;
                    pair.triangle_tests += sample.triangle_tests;
                    AddFilter(pair.filter, sample.filter);
                    if (truth != nullptr)
                    {
                        const bool expected = (*truth)[index] != 0;
                        AddConfusion(worker.confusion, sample.positive, expected);
                        AddConfusion(pair.confusion, sample.positive, expected);
                    }
                }
                worker.summed_call_ms = Milliseconds(call_started, Clock::now());
            }
            catch (...)
            {
                std::lock_guard lock(error_mutex);
                if (error == nullptr) error = std::current_exception();
            }
            finished.count_down();
        });
    }
    ready.wait();
    const auto wall_started = Clock::now();
    start.arrive_and_wait();
    finished.wait();
    const double wall_ms = Milliseconds(wall_started, Clock::now());
    for (auto& thread : threads) thread.join();
    if (error != nullptr) std::rethrow_exception(error);

    Metrics result;
    result.wall_ms = wall_ms;
    for (const Metrics& worker : workers)
    {
        result.summed_call_ms += worker.summed_call_ms;
        result.positives += worker.positives;
        result.root_tests += worker.root_tests;
        result.internal_tests += worker.internal_tests;
        result.triangle_tests += worker.triangle_tests;
        AddFilter(result.filter, worker.filter);
        result.confusion.tp += worker.confusion.tp;
        result.confusion.tn += worker.confusion.tn;
        result.confusion.fp += worker.confusion.fp;
        result.confusion.fn += worker.confusion.fn;
        for (std::size_t index = 0; index < result.pairs.size(); ++index)
        {
            auto& out = result.pairs[index];
            const auto& in = worker.pairs[index];
            out.queries += in.queries;
            out.positives += in.positives;
            out.root_tests += in.root_tests;
            out.internal_tests += in.internal_tests;
            out.triangle_tests += in.triangle_tests;
            AddFilter(out.filter, in.filter);
            out.confusion.tp += in.confusion.tp;
            out.confusion.tn += in.confusion.tn;
            out.confusion.fp += in.confusion.fp;
            out.confusion.fn += in.confusion.fn;
        }
    }
    return result;
}

Sample StandardSample(const pqss::ModelPool& pool, const QueryPair& query)
{
    pqss::QueryResult result;
    Require(pool.Query(result, query.first, query.R1, query.T1,
                       query.second, query.R2, query.T2, query.tolerance), "Query");
    return {result.CloserThanTolerance(), 1, result.NumBvTests(), result.NumTriTests(), {},
            result.QueryTimeMs()};
}

Sample IndependentSample(const pqss_proxy_mesh::IndependentPrimitivePool& pool,
                         const QueryPair& query,
                         const pqss_proxy_mesh::DistanceFilterMode filter_mode)
{
    const auto result = pool.query(query.first, query.R1, query.T1,
                                   query.second, query.R2, query.T2, query.tolerance,
                                   filter_mode);
    return {result.closer_than_tolerance, result.primitive_pair_tests,
            result.internal_bv_tests, result.triangle_tests,
            {
                result.distance_filter_skips,
                result.containment_filter_checks,
                result.containment_filter_skips,
                result.axis_filter_checks,
                result.axis_filter_skips,
                result.axis_generation_checks,
                result.cached_separating_axes,
            },
            result.query_time_ms};
}

double FalsePositiveRate(const Confusion& confusion)
{
    const std::size_t negatives = confusion.fp + confusion.tn;
    return negatives == 0 ? 0.0 : static_cast<double>(confusion.fp) / negatives;
}

pqss_proxy_mesh::DistanceFilterMode ParseDistanceFilter(const std::string_view value)
{
    using pqss_proxy_mesh::DistanceFilterMode;
    if (value == "none") return DistanceFilterMode::None;
    if (value == "containment") return DistanceFilterMode::Containment;
    if (value == "separating-axis" || value == "axis") return DistanceFilterMode::SeparatingAxis;
    if (value == "hybrid") return DistanceFilterMode::Hybrid;
    throw std::runtime_error("distance filter must be none, containment, separating-axis, or hybrid");
}

void WritePairCsv(const std::filesystem::path& path,
                  const Metrics& official,
                  const Metrics& candidate)
{
    std::ofstream stream(path);
    if (!stream) throw std::runtime_error("cannot write " + path.string());
    stream << "model_id_1,model_id_2,queries,official_internal_bv,candidate_internal_bv,"
              "official_root_bv,candidate_primitive_root_bv,official_total_rss,candidate_total_rss,"
              "official_triangles,candidate_triangles,distance_filter_skips,containment_filter_checks,"
              "containment_filter_skips,axis_filter_checks,axis_filter_skips,official_fp,official_fn,"
              "candidate_fp,candidate_fn\n";
    const auto ids = ObjectIds();
    for (std::size_t first = 0; first < kLogicalModelCount; ++first)
    {
        for (std::size_t second = 0; second < kLogicalModelCount; ++second)
        {
            const std::size_t index = first * kLogicalModelCount + second;
            const auto& o = official.pairs[index];
            const auto& c = candidate.pairs[index];
            if (o.queries == 0 && c.queries == 0) continue;
            stream << ids[first] << ',' << ids[second] << ',' << o.queries << ','
                   << o.internal_tests << ',' << c.internal_tests << ','
                   << o.root_tests << ',' << c.root_tests << ','
                   << o.root_tests + o.internal_tests << ',' << c.root_tests + c.internal_tests << ','
                   << o.triangle_tests << ',' << c.triangle_tests << ',' << c.filter.skips << ','
                   << c.filter.containment_checks << ',' << c.filter.containment_skips << ','
                   << c.filter.axis_checks << ',' << c.filter.axis_skips << ','
                   << o.confusion.fp << ',' << o.confusion.fn << ','
                   << c.confusion.fp << ',' << c.confusion.fn << '\n';
        }
    }
}

void WriteModelCsv(const std::filesystem::path& path,
                   const StandardBuild& official,
                   const IndependentBuild& candidate)
{
    std::ofstream stream(path);
    if (!stream) throw std::runtime_error("cannot write " + path.string());
    stream << "model_id,official_primitives,candidate_primitives,official_triangles,candidate_triangles,"
              "official_bvs,candidate_bvs,official_max_depth,candidate_max_depth\n";
    const auto ids = ObjectIds();
    const auto official_models = official.pool->GetModelStats();
    const auto official_bvhs = official.pool->BvhSignatures();
    const auto candidate_models = candidate.pool->physicalPool().GetModelStats();
    const auto candidate_bvhs = candidate.pool->physicalPool().BvhSignatures();
    for (std::size_t logical = 0; logical < kLogicalModelCount; ++logical)
    {
        std::size_t triangles = 0;
        std::size_t bvs = 0;
        std::size_t depth = 0;
        for (std::size_t primitive = 0; primitive < candidate.pool->primitiveCount(logical); ++primitive)
        {
            const std::size_t physical = candidate.pool->physicalModelId(logical, primitive);
            triangles += candidate_models.at(physical).built_tri_count;
            bvs += candidate_models.at(physical).bv_count;
            depth = std::max(depth, candidate_bvhs.at(physical).max_depth);
        }
        stream << ids[logical] << ",1," << candidate.pool->primitiveCount(logical) << ','
               << official_models.at(logical).built_tri_count << ',' << triangles << ','
               << official_models.at(logical).bv_count << ',' << bvs << ','
               << official_bvhs.at(logical).max_depth << ',' << depth << '\n';
    }
}

void WriteReport(const std::filesystem::path& path,
                 const std::size_t query_count,
                 const std::size_t threads,
                 const StandardBuild& original,
                 const StandardBuild& official,
                 const IndependentBuild& candidate,
                 const Metrics& original_metrics,
                 const Metrics& official_metrics,
                 const Metrics& candidate_metrics,
                 const pqss_proxy_mesh::DistanceFilterMode filter_mode)
{
    std::ofstream stream(path);
    if (!stream) throw std::runtime_error("cannot write " + path.string());
    const auto& containment = candidate.pool->containmentStats();
    stream << std::fixed << std::setprecision(6)
           << "# Independent Primitive RSS BVH Comparison\n\n"
           << "- Queries: " << query_count << "\n"
           << "- Threads: " << threads << "\n"
           << "- Truth: original geometry, original one-model/one-BVH query\n"
           << "- Official: official simplified geometry, original one-model/one-BVH query\n"
           << "- Candidate: independent primitive BLASes, cross-tree RSS containment pruning\n"
           << "- Candidate builder: project-local six-axis Optimized builder\n\n"
           << "- Distance filter: `" << pqss_proxy_mesh::distanceFilterModeName(filter_mode) << "`\n\n"
           << "## Build\n\n"
           << "| Metric | Original | Official | Candidate |\n|---|---:|---:|---:|\n"
           << "| Base build ms |" << original.base_ms << '|' << official.base_ms << '|' << candidate.base_ms << "|\n"
           << "| Total build ms |" << original.total_ms << '|' << official.total_ms << '|' << candidate.total_ms << "|\n"
           << "| Containment analysis ms |0|0|" << candidate.containment_ms << "|\n"
           << "| PQSS pool memory bytes |" << original.memory_bytes << '|' << official.memory_bytes << '|'
           << candidate.memory_bytes << "|\n"
           << "| Physical models |" << kLogicalModelCount << '|' << kLogicalModelCount << '|'
           << containment.physical_models << "|\n"
           << "| RSS BVH nodes |-| -|" << containment.bvh_nodes << "|\n"
           << "| Tree ancestor links |0|0|" << containment.tree_ancestor_links << "|\n"
           << "| Cross-tree containment links |0|0|" << containment.cross_tree_containment_links << "|\n\n"
           << "## Query Work\n\n"
           << "The original PQSS counter excludes each model-pair root RSS test. Both views are reported.\n\n"
           << "| Metric | Original truth | Official simplified | Candidate |\n|---|---:|---:|---:|\n"
           << "| PQSS internal BV tests |" << original_metrics.internal_tests << '|'
           << official_metrics.internal_tests << '|' << candidate_metrics.internal_tests << "|\n"
           << "| Root/scheduler RSS tests |" << original_metrics.root_tests << '|'
           << official_metrics.root_tests << '|' << candidate_metrics.root_tests << "|\n"
           << "| All RSS distance tests |" << original_metrics.root_tests + original_metrics.internal_tests << '|'
           << official_metrics.root_tests + official_metrics.internal_tests << '|'
           << candidate_metrics.root_tests + candidate_metrics.internal_tests << "|\n"
           << "| Triangle tests |" << original_metrics.triangle_tests << '|'
           << official_metrics.triangle_tests << '|' << candidate_metrics.triangle_tests << "|\n"
           << "| Distance-filter skips |0|0|" << candidate_metrics.filter.skips << "|\n"
           << "| Containment checks |0|0|" << candidate_metrics.filter.containment_checks << "|\n"
           << "| Containment skips |0|0|" << candidate_metrics.filter.containment_skips << "|\n"
           << "| Separating-axis checks |0|0|" << candidate_metrics.filter.axis_checks << "|\n"
           << "| Separating-axis skips |0|0|" << candidate_metrics.filter.axis_skips << "|\n"
           << "| Separating-axis generation checks |0|0|"
           << candidate_metrics.filter.axis_generation_checks << "|\n"
           << "| Positive results |" << original_metrics.positives << '|'
           << official_metrics.positives << '|' << candidate_metrics.positives << "|\n"
           << "| Parallel wall ms |" << original_metrics.wall_ms << '|'
           << official_metrics.wall_ms << '|' << candidate_metrics.wall_ms << "|\n\n"
           << "## Correctness Against Original Geometry\n\n"
           << "| Candidate | TP | TN | FP | FN | False-positive rate |\n|---|---:|---:|---:|---:|---:|\n"
           << "| Official simplified |" << official_metrics.confusion.tp << '|'
           << official_metrics.confusion.tn << '|' << official_metrics.confusion.fp << '|'
           << official_metrics.confusion.fn << '|' << FalsePositiveRate(official_metrics.confusion) << "|\n"
           << "| Independent primitives |" << candidate_metrics.confusion.tp << '|'
           << candidate_metrics.confusion.tn << '|' << candidate_metrics.confusion.fp << '|'
           << candidate_metrics.confusion.fn << '|' << FalsePositiveRate(candidate_metrics.confusion) << "|\n\n"
           << "A conservative candidate must have `FN = 0`. Detailed data are in "
              "`model_bvh_stats.csv` and `query_pair_stats.csv`.\n";
}

} // namespace

int main(int argc, char** argv)
{
    try
    {
        if (argc < 7 || argc > 12)
        {
            std::cerr << "Usage: pqss-independent-real-scene-compare <original-dir> <official-dir> "
                         "<primitive-manifest.tsv> <workpiece-dir> <queries.txt> <report-dir> [threads] "
                         "[original-cache-tag] [official-cache-tag] [candidate-cache-tag] [distance-filter]\n";
            return 2;
        }
        std::size_t threads = 1;
        if (argc >= 8)
        {
            const std::string_view value = argv[7];
            const auto parsed = std::from_chars(value.data(), value.data() + value.size(), threads);
            if (parsed.ec != std::errc{} || parsed.ptr != value.data() + value.size() || threads == 0)
                throw std::runtime_error("threads must be a positive integer");
        }
        const std::string original_cache = argc >= 9 ? argv[8] : "original_91194221491100";
        const std::string official_cache = argc >= 10 ? argv[9] : "official_simplified_v1";
        const std::string candidate_cache =
            argc >= 11 ? argv[10] : "independent_primitives_strength8_proxy_bvh_v3";
        const pqss_proxy_mesh::DistanceFilterMode filter_mode =
            argc >= 12 ? ParseDistanceFilter(argv[11]) : pqss_proxy_mesh::DistanceFilterMode::SeparatingAxis;
        const auto queries = LoadQueries(argv[5]);
        std::cout << "queries=" << queries.size() << "\nBuilding original truth pool...\n" << std::flush;
        const auto original = BuildStandard(argv[1], argv[4], original_cache);
        std::cout << "Building official simplified pool...\n" << std::flush;
        const auto official = BuildStandard(argv[2], argv[4], official_cache);
        std::cout << "Building independent primitive pool...\n" << std::flush;
        const auto candidate = BuildIndependent(argv[3], argv[4], candidate_cache);

        std::vector<std::uint8_t> truth(queries.size());
        std::cout << "Running original truth queries...\n" << std::flush;
        const Metrics original_metrics = RunQueries(
            queries, threads, [&](const QueryPair& query) { return StandardSample(*original.pool, query); },
            &truth, nullptr);
        std::cout << "Running official simplified queries...\n" << std::flush;
        const Metrics official_metrics = RunQueries(
            queries, threads, [&](const QueryPair& query) { return StandardSample(*official.pool, query); },
            nullptr, &truth);
        std::cout << "Running independent primitive queries...\n" << std::flush;
        const Metrics candidate_metrics = RunQueries(
            queries, threads, [&](const QueryPair& query)
            {
                return IndependentSample(*candidate.pool, query, filter_mode);
            },
            nullptr, &truth);

        const std::filesystem::path report_dir = argv[6];
        std::filesystem::create_directories(report_dir);
        WriteModelCsv(report_dir / "model_bvh_stats.csv", official, candidate);
        WritePairCsv(report_dir / "query_pair_stats.csv", official_metrics, candidate_metrics);
        WriteReport(report_dir / "comparison_report.md", queries.size(), threads,
                    original, official, candidate,
                    original_metrics, official_metrics, candidate_metrics, filter_mode);
        std::cout << "official_internal_bv=" << official_metrics.internal_tests
                  << "\ncandidate_internal_bv=" << candidate_metrics.internal_tests
                  << "\ncandidate_scheduler_bv=" << candidate_metrics.root_tests
                  << "\ncandidate_triangle_tests=" << candidate_metrics.triangle_tests
                  << "\ndistance_filter=" << pqss_proxy_mesh::distanceFilterModeName(filter_mode)
                  << "\ncandidate_filter_skips=" << candidate_metrics.filter.skips
                  << "\ncandidate_containment_skips=" << candidate_metrics.filter.containment_skips
                  << "\ncandidate_axis_skips=" << candidate_metrics.filter.axis_skips
                  << "\ncandidate_fp=" << candidate_metrics.confusion.fp
                  << "\ncandidate_fn=" << candidate_metrics.confusion.fn
                  << "\nreport=" << (report_dir / "comparison_report.md").string() << '\n';
        return candidate_metrics.confusion.fn == 0 ? 0 : 3;
    }
    catch (const std::exception& error)
    {
        std::cerr << "error: " << error.what() << '\n';
        return 1;
    }
}
