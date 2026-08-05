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
#include "pqss_proxy_mesh/primitive_bvh_pool.hpp"
#include "pqss_proxy_mesh/proxy_bvh_build.hpp"

namespace
{

using Clock = std::chrono::steady_clock;
constexpr std::size_t kModelCount = 36;
constexpr std::array<int, 13> kWorkpieces = {24, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37};

struct Query
{
    std::size_t first_pool = 0;
    std::size_t second_pool = 0;
    std::size_t first_object = 0;
    std::size_t second_object = 0;
    pqss::Mat3 R1{};
    pqss::Vec3 T1{};
    pqss::Mat3 R2{};
    pqss::Vec3 T2{};
    pqss::Real tolerance = pqss::k_zero;
};

struct Confusion
{
    std::size_t tp = 0;
    std::size_t tn = 0;
    std::size_t fp = 0;
    std::size_t fn = 0;
};

struct Sample
{
    bool positive = false;
    std::size_t root_bv_tests = 0;
    std::size_t internal_bv_tests = 0;
    std::size_t leaf_tests = 0;
    std::array<std::array<std::size_t, 3>, 3> type_tests{};
};

struct PairMetrics
{
    std::size_t queries = 0;
    std::size_t positives = 0;
    std::size_t root_bv_tests = 0;
    std::size_t internal_bv_tests = 0;
    std::size_t leaf_tests = 0;
    Confusion confusion;
};

struct Metrics
{
    double wall_ms = 0.0;
    std::size_t positives = 0;
    std::size_t root_bv_tests = 0;
    std::size_t internal_bv_tests = 0;
    std::size_t leaf_tests = 0;
    std::array<std::array<std::size_t, 3>, 3> type_tests{};
    Confusion confusion;
    std::array<PairMetrics, kModelCount * kModelCount> pairs{};
};

struct BuildResult
{
    std::unique_ptr<pqss::ModelPool> reference;
    std::unique_ptr<pqss::ModelPool> official;
    std::unique_ptr<pqss::ModelPool> candidate_pqss;
    std::unique_ptr<pqss_proxy_mesh::PrimitiveBvhPool> candidate;
    double reference_ms = 0.0;
    double official_ms = 0.0;
    double candidate_import_ms = 0.0;
    double analyzed_bvh_ms = 0.0;
    double pqss_install_ms = 0.0;
    std::size_t reference_memory_bytes = 0;
    std::size_t official_memory_bytes = 0;
    pqss_proxy_mesh::PrimitiveBvhCacheStats cache_stats;
};

double Milliseconds(const Clock::time_point begin, const Clock::time_point end)
{
    return std::chrono::duration<double, std::milli>(end - begin).count();
}

std::size_t PoolIndex(const std::size_t object_id)
{
    if (object_id >= 1 && object_id <= 23) return object_id - 1;
    if (object_id == 24) return 23;
    if (object_id >= 26 && object_id <= 28) return 24 + object_id - 26;
    if (object_id >= 29 && object_id <= 32) return 27 + object_id - 29;
    if (object_id >= 33 && object_id <= 37) return 31 + object_id - 33;
    return std::numeric_limits<std::size_t>::max();
}

std::vector<std::size_t> ObjectIds()
{
    std::vector<std::size_t> result;
    for (std::size_t id = 1; id <= 23; ++id) result.push_back(id);
    for (const int id : kWorkpieces) result.push_back(static_cast<std::size_t>(id));
    return result;
}

std::filesystem::path ModelPath(const std::filesystem::path& directory, const int id)
{
    return directory / (std::to_string(id) + ".obj");
}

void Require(const pqss::Return result, const std::string& operation)
{
    if (result != pqss::Return::Success)
        throw std::runtime_error(operation + " failed with Return=" +
                                 std::to_string(static_cast<int>(result)));
}

std::size_t ParseObjectId(const std::string_view token)
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

template <typename Value>
bool ParseNumber(const char*& cursor, const char* end, Value& value)
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

std::vector<Query> LoadQueries(const std::filesystem::path& path)
{
    std::ifstream stream(path);
    if (!stream) throw std::runtime_error("cannot open query file: " + path.string());
    std::vector<Query> queries;
    queries.reserve(750'000);
    std::string line;
    std::size_t line_number = 0;
    while (std::getline(stream, line))
    {
        ++line_number;
        if (line.empty()) continue;
        const char* cursor = line.data();
        const char* end = cursor + line.size();
        Query query;
        query.first_object = ParseObjectId(ReadToken(cursor, end));
        if (!ParsePose(cursor, end, query.R1, query.T1))
            throw std::runtime_error("invalid first pose on line " + std::to_string(line_number));
        query.second_object = ParseObjectId(ReadToken(cursor, end));
        if (!ParsePose(cursor, end, query.R2, query.T2) ||
            !ParseNumber(cursor, end, query.tolerance))
            throw std::runtime_error("invalid second pose on line " + std::to_string(line_number));
        query.first_pool = PoolIndex(query.first_object);
        query.second_pool = PoolIndex(query.second_object);
        if (query.first_pool >= kModelCount || query.second_pool >= kModelCount)
            throw std::runtime_error("unknown model ID on line " + std::to_string(line_number));
        queries.push_back(query);
    }
    return queries;
}

std::unique_ptr<pqss::ModelPool> BuildStandardPool(const std::filesystem::path& base,
                                                    const std::filesystem::path& workpieces,
                                                    const std::string& cache_tag)
{
    auto pool = std::make_unique<pqss::ModelPool>();
    Require(pool->SetBuildCacheTag(cache_tag), "SetBuildCacheTag");
    for (int id = 1; id <= 23; ++id)
        Require(pool->LoadModelFromFile(ModelPath(base, id).string(), pqss::BuildStrategy::Optimized),
                "LoadModelFromFile(" + std::to_string(id) + ")");
    Require(pool->EndPool(), "EndPool");
    for (const int id : kWorkpieces)
        Require(pool->LoadModelFromFile(ModelPath(workpieces, id).string(), pqss::BuildStrategy::Optimized),
                "append workpiece " + std::to_string(id));
    return pool;
}

BuildResult Build(const std::filesystem::path& original_base,
                  const std::filesystem::path& official_base,
                  const std::filesystem::path& workpieces,
                  const std::filesystem::path& primitive_manifest,
                  const std::string& reference_cache_tag,
                  const std::string& official_cache_tag)
{
    pqss_proxy_mesh::setProxyBvhBuildEnabled(false);
    BuildResult result;
    const auto reference_started = Clock::now();
    result.reference = BuildStandardPool(original_base, workpieces, reference_cache_tag);
    result.reference_ms = Milliseconds(reference_started, Clock::now());
    result.reference_memory_bytes = result.reference->Info(false);

    const auto official_started = Clock::now();
    result.official = BuildStandardPool(official_base, workpieces, official_cache_tag);
    result.official_ms = Milliseconds(official_started, Clock::now());
    result.official_memory_bytes = result.official->Info(false);

    result.candidate_pqss = std::make_unique<pqss::ModelPool>(*result.reference);
    result.candidate = std::make_unique<pqss_proxy_mesh::PrimitiveBvhPool>();
    const auto import_started = Clock::now();
    result.candidate->importReferencePool(*result.reference, ObjectIds());
    result.candidate_import_ms = Milliseconds(import_started, Clock::now());
    const auto analyzed_started = Clock::now();
    result.candidate->loadAnalyzedPool(primitive_manifest, pqss::BuildStrategy::Optimized,
                                       primitive_manifest.parent_path() / "bvh_cache");
    result.analyzed_bvh_ms = Milliseconds(analyzed_started, Clock::now());
    result.cache_stats = result.candidate->cacheStats();
    const auto install_started = Clock::now();
    result.candidate->installIntoPqssPool(*result.candidate_pqss, ObjectIds());
    result.pqss_install_ms = Milliseconds(install_started, Clock::now());
    return result;
}

void AddConfusion(Confusion& result, const bool actual, const bool expected)
{
    if (actual && expected) ++result.tp;
    else if (!actual && !expected) ++result.tn;
    else if (actual) ++result.fp;
    else ++result.fn;
}

template <typename Function>
Metrics RunQueries(const std::vector<Query>& queries,
                   const std::size_t requested_threads,
                   Function function,
                   std::vector<std::uint8_t>* decisions,
                   const std::vector<std::uint8_t>* truth)
{
    const std::size_t thread_count = std::min(requested_threads, queries.size());
    if (thread_count == 0) return {};
    std::vector<Metrics> workers(thread_count);
    std::vector<std::thread> threads;
    threads.reserve(thread_count);
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
            try
            {
                Metrics& worker = workers[thread_index];
                for (std::size_t index = begin; index < end; ++index)
                {
                    const Query& query = queries[index];
                    const Sample sample = function(query);
                    if (decisions != nullptr) (*decisions)[index] = sample.positive ? 1 : 0;
                    worker.positives += sample.positive;
                    worker.root_bv_tests += sample.root_bv_tests;
                    worker.internal_bv_tests += sample.internal_bv_tests;
                    worker.leaf_tests += sample.leaf_tests;
                    for (std::size_t first = 0; first < 3; ++first)
                        for (std::size_t second = 0; second < 3; ++second)
                            worker.type_tests[first][second] += sample.type_tests[first][second];
                    PairMetrics& pair = worker.pairs[query.first_pool * kModelCount + query.second_pool];
                    ++pair.queries;
                    pair.positives += sample.positive;
                    pair.root_bv_tests += sample.root_bv_tests;
                    pair.internal_bv_tests += sample.internal_bv_tests;
                    pair.leaf_tests += sample.leaf_tests;
                    if (truth != nullptr)
                    {
                        const bool expected = (*truth)[index] != 0;
                        AddConfusion(worker.confusion, sample.positive, expected);
                        AddConfusion(pair.confusion, sample.positive, expected);
                    }
                }
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
    const auto started = Clock::now();
    start.arrive_and_wait();
    finished.wait();
    const double wall_ms = Milliseconds(started, Clock::now());
    for (auto& thread : threads) thread.join();
    if (error != nullptr) std::rethrow_exception(error);

    Metrics result;
    result.wall_ms = wall_ms;
    for (const Metrics& worker : workers)
    {
        result.positives += worker.positives;
        result.root_bv_tests += worker.root_bv_tests;
        result.internal_bv_tests += worker.internal_bv_tests;
        result.leaf_tests += worker.leaf_tests;
        result.confusion.tp += worker.confusion.tp;
        result.confusion.tn += worker.confusion.tn;
        result.confusion.fp += worker.confusion.fp;
        result.confusion.fn += worker.confusion.fn;
        for (std::size_t first = 0; first < 3; ++first)
            for (std::size_t second = 0; second < 3; ++second)
                result.type_tests[first][second] += worker.type_tests[first][second];
        for (std::size_t index = 0; index < result.pairs.size(); ++index)
        {
            PairMetrics& out = result.pairs[index];
            const PairMetrics& in = worker.pairs[index];
            out.queries += in.queries;
            out.positives += in.positives;
            out.root_bv_tests += in.root_bv_tests;
            out.internal_bv_tests += in.internal_bv_tests;
            out.leaf_tests += in.leaf_tests;
            out.confusion.tp += in.confusion.tp;
            out.confusion.tn += in.confusion.tn;
            out.confusion.fp += in.confusion.fp;
            out.confusion.fn += in.confusion.fn;
        }
    }
    return result;
}

Sample ReferenceSample(const pqss::ModelPool& pool, const Query& query)
{
    pqss::QueryResult result;
    Require(pool.Query(result, query.first_pool, query.R1, query.T1,
                       query.second_pool, query.R2, query.T2, query.tolerance), "Query");
    return {result.CloserThanTolerance(), 1, result.NumBvTests(), result.NumTriTests(), {}};
}

Sample CandidateSample(const pqss::ModelPool& pool, const Query& query)
{
    pqss::QueryResult result;
    Require(pool.Query(result, query.first_pool, query.R1, query.T1,
                       query.second_pool, query.R2, query.T2, query.tolerance), "candidate Query");
    return {result.CloserThanTolerance(), 1, result.NumBvTests(), result.NumTriTests(), {}};
}

double FalsePositiveRate(const Confusion& confusion)
{
    const std::size_t negatives = confusion.tn + confusion.fp;
    return negatives == 0 ? 0.0 : static_cast<double>(confusion.fp) / static_cast<double>(negatives);
}

void WriteModels(const std::filesystem::path& path,
                 const BuildResult& build)
{
    std::ofstream stream(path);
    if (!stream) throw std::runtime_error("cannot write " + path.string());
    stream << "model_id,analyzed,source_triangles,primitive_leaves,primitive_nodes,primitive_max_depth,"
              "official_input_triangles,official_built_triangles,official_nodes,official_max_depth,"
              "spheres,capsules,rss,containment_certified\n";
    const auto official_models = build.official->GetModelStats();
    const auto official_bvhs = build.official->BvhSignatures();
    const auto ids = ObjectIds();
    for (std::size_t index = 0; index < ids.size(); ++index)
    {
        const std::size_t id = ids[index];
        const auto& stats = build.candidate->modelStats(id);
        stream << id << ',' << (stats.analyzed ? 1 : 0) << ',' << stats.source_triangles << ','
               << stats.primitives << ',' << stats.nodes << ',' << stats.max_depth << ','
               << official_models[index].original_tri_count << ',' << official_models[index].built_tri_count << ','
               << official_models[index].bv_count << ',' << official_bvhs[index].max_depth << ','
               << stats.primitive_types[0] << ',' << stats.primitive_types[1] << ','
               << stats.primitive_types[2] << ',' << (stats.internal_containment_certified ? 1 : 0) << '\n';
    }
}

void WritePairs(const std::filesystem::path& path,
                const Metrics& reference,
                const Metrics& official,
                const Metrics& candidate)
{
    std::ofstream stream(path);
    if (!stream) throw std::runtime_error("cannot write " + path.string());
    stream << "model_id_1,model_id_2,queries,reference_bv_tests,official_bv_tests,"
              "candidate_root_bv_tests,candidate_internal_bv_tests,reference_triangle_tests,"
              "official_triangle_tests,candidate_leaf_tests,official_tp,official_tn,official_fp,official_fn,"
              "candidate_tp,candidate_tn,candidate_fp,candidate_fn\n";
    const auto ids = ObjectIds();
    for (std::size_t first = 0; first < kModelCount; ++first)
    {
        for (std::size_t second = 0; second < kModelCount; ++second)
        {
            const std::size_t index = first * kModelCount + second;
            const auto& r = reference.pairs[index];
            const auto& o = official.pairs[index];
            const auto& c = candidate.pairs[index];
            if (r.queries == 0 && c.queries == 0) continue;
            stream << ids[first] << ',' << ids[second] << ',' << r.queries << ','
                   << r.internal_bv_tests << ',' << o.internal_bv_tests << ',' << c.root_bv_tests << ','
                   << c.internal_bv_tests << ',' << r.leaf_tests << ',' << o.leaf_tests << ','
                   << c.leaf_tests << ','
                   << o.confusion.tp << ',' << o.confusion.tn << ','
                   << o.confusion.fp << ',' << o.confusion.fn << ','
                   << c.confusion.tp << ',' << c.confusion.tn << ','
                   << c.confusion.fp << ',' << c.confusion.fn << '\n';
        }
    }
}

void WriteReport(const std::filesystem::path& path,
                 const std::size_t query_count,
                 const std::size_t threads,
                 const BuildResult& build,
                 const Metrics& reference,
                 const Metrics& official,
                 const Metrics& candidate)
{
    std::ofstream stream(path);
    if (!stream) throw std::runtime_error("cannot write " + path.string());
    stream << std::fixed << std::setprecision(6)
           << "# Primitive-leaf BVH Real-scene Comparison\n\n"
           << "- Queries: " << query_count << "\n"
           << "- Threads: " << threads << "\n"
           << "- Truth: original source geometry with the official PQSS Optimized builder\n"
           << "- Official simplified: dataset `obj_simplified`, standard PQSS Optimized builder and query\n"
           << "- Candidate: one BVH per logical model, analytic primitives as leaves\n"
           << "- Internal nodes: PQSS Optimized RSS fitted to descendant responsibility triangles\n"
           << "- Build time is reported only and is not an optimization objective.\n\n"
           << "## Build\n\n"
           << "| Metric | Value |\n|---|---:|\n"
           << "| Original pool build ms |" << build.reference_ms << "|\n"
           << "| Official simplified pool build ms |" << build.official_ms << "|\n"
           << "| Reference import ms |" << build.candidate_import_ms << "|\n"
           << "| Analyzed primitive BVH build ms |" << build.analyzed_bvh_ms << "|\n"
           << "| Primitive BVH cache hits |" << build.cache_stats.cache_hits << "|\n"
           << "| Primitive BVH cache misses |" << build.cache_stats.cache_misses << "|\n"
           << "| Primitive BVH cache load ms |" << build.cache_stats.cache_load_ms << "|\n"
           << "| Primitive BVH uncached build ms |" << build.cache_stats.cache_build_ms << "|\n"
           << "| Primitive BVH cache save ms |" << build.cache_stats.cache_save_ms << "|\n"
           << "| Install into PQSS ModelPool ms |" << build.pqss_install_ms << "|\n"
           << "| Reference pool memory bytes |" << build.reference_memory_bytes << "|\n"
           << "| Official simplified pool memory bytes |" << build.official_memory_bytes << "|\n\n"
           << "## Query Work\n\n"
           << "PQSS's counter omits the model-pair root test, so one root test per query is shown separately.\n\n"
           << "| Metric | Original geometry | Official simplified | Primitive BVH |\n|---|---:|---:|---:|\n"
           << "| Root BV tests |" << reference.root_bv_tests << '|' << official.root_bv_tests << '|'
           << candidate.root_bv_tests << "|\n"
           << "| Internal BV tests |" << reference.internal_bv_tests << '|' << official.internal_bv_tests << '|'
           << candidate.internal_bv_tests << "|\n"
           << "| Total BV tests |" << reference.root_bv_tests + reference.internal_bv_tests << '|'
           << official.root_bv_tests + official.internal_bv_tests << '|'
           << candidate.root_bv_tests + candidate.internal_bv_tests << "|\n"
           << "| Triangle/leaf-pair tests |" << reference.leaf_tests << '|' << official.leaf_tests << '|'
           << candidate.leaf_tests << "|\n"
           << "| Positives |" << reference.positives << '|' << official.positives << '|'
           << candidate.positives << "|\n"
           << "| Parallel query wall ms |" << reference.wall_ms << '|' << official.wall_ms << '|'
           << candidate.wall_ms << "|\n\n"
           << "## Accuracy Against Original Geometry\n\n"
           << "| Candidate | TP | TN | FP | FN | False-positive rate |\n|---|---:|---:|---:|---:|---:|\n"
           << "| Official simplified |" << official.confusion.tp << '|' << official.confusion.tn << '|'
           << official.confusion.fp << '|' << official.confusion.fn << '|'
           << FalsePositiveRate(official.confusion) << "|\n"
           << "| Primitive BVH |" << candidate.confusion.tp << '|' << candidate.confusion.tn << '|'
           << candidate.confusion.fp << '|' << candidate.confusion.fn << '|'
           << FalsePositiveRate(candidate.confusion) << "|\n\n"
           << "A valid conservative result requires `FN = 0`.\n\n"
           << "The candidate uses the same `ModelPool::Query` and `QueryRecurse` implementation as PQSS. "
              "Its leaf-test counter records analytic RSS leaf pairs rather than triangle pairs.\n";
}

} // namespace

int main(int argc, char** argv)
{
    try
    {
        if (argc < 7 || argc > 10)
        {
            std::cerr << "Usage: pqss-primitive-bvh-real-scene-compare <original-base-dir> "
                         "<official-simplified-dir> <workpiece-dir> <primitive-manifest.tsv> "
                         "<queries.txt> <report-dir> [threads] [reference-cache-tag] "
                         "[official-cache-tag]\n";
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
        const std::string reference_cache = argc >= 9 ? argv[8] : "original_91194221491100";
        const std::string official_cache = argc >= 10 ? argv[9] : "official_simplified_v1";
        const auto queries = LoadQueries(argv[5]);
        std::cout << "queries=" << queries.size() << "\nBuilding all three pools...\n" << std::flush;
        const BuildResult build = Build(argv[1], argv[2], argv[3], argv[4],
                                        reference_cache, official_cache);

        std::vector<std::uint8_t> truth(queries.size());
        std::cout << "Running original-geometry queries...\n" << std::flush;
        const Metrics reference = RunQueries(
            queries, threads,
            [&](const Query& query) { return ReferenceSample(*build.reference, query); },
            &truth, nullptr);
        std::cout << "Running official-simplified queries...\n" << std::flush;
        const Metrics official = RunQueries(
            queries, threads,
            [&](const Query& query) { return ReferenceSample(*build.official, query); },
            nullptr, &truth);
        std::cout << "Running primitive-BVH queries...\n" << std::flush;
        const Metrics candidate = RunQueries(
            queries, threads,
            [&](const Query& query) { return CandidateSample(*build.candidate_pqss, query); },
            nullptr, &truth);

        const std::filesystem::path report_dir = argv[6];
        std::filesystem::create_directories(report_dir);
        WriteModels(report_dir / "model_bvh_stats.csv", build);
        WritePairs(report_dir / "query_pair_stats.csv", reference, official, candidate);
        WriteReport(report_dir / "comparison_report.md", queries.size(), threads,
                    build, reference, official, candidate);
        std::cout << "primitive_bvh_build_ms=" << build.analyzed_bvh_ms
                  << "\nprimitive_bvh_cache_hits=" << build.cache_stats.cache_hits
                  << "\nprimitive_bvh_cache_misses=" << build.cache_stats.cache_misses
                  << "\nreference_total_bv_tests="
                  << reference.root_bv_tests + reference.internal_bv_tests
                  << "\nofficial_total_bv_tests="
                  << official.root_bv_tests + official.internal_bv_tests
                  << "\ncandidate_total_bv_tests="
                  << candidate.root_bv_tests + candidate.internal_bv_tests
                  << "\nreference_triangle_tests=" << reference.leaf_tests
                  << "\nofficial_triangle_tests=" << official.leaf_tests
                  << "\ncandidate_leaf_tests=" << candidate.leaf_tests
                  << "\nreference_query_wall_ms=" << reference.wall_ms
                  << "\nofficial_query_wall_ms=" << official.wall_ms
                  << "\ncandidate_query_wall_ms=" << candidate.wall_ms
                  << "\ncandidate_fp=" << candidate.confusion.fp
                  << "\ncandidate_fn=" << candidate.confusion.fn
                  << "\nofficial_fp=" << official.confusion.fp
                  << "\nofficial_fn=" << official.confusion.fn
                  << "\nreport=" << (report_dir / "comparison_report.md").string() << '\n';
        return candidate.confusion.fn == 0 ? 0 : 3;
    }
    catch (const std::exception& error)
    {
        std::cerr << "error: " << error.what() << '\n';
        return 1;
    }
}
