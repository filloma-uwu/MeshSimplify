#include <array>
#include <chrono>
#include <filesystem>
#include <fstream>
#include <string>
#include <stdexcept>

#include "pqss/model_pool.hpp"
#include "pqss/utils/mat_vec.hpp"
#include "pqss_proxy_mesh/primitive_bvh_pool.hpp"

namespace
{

void Check(const bool condition)
{
    if (!condition) throw std::runtime_error("primitive BVH test check failed");
}

void WriteModelOne(const std::filesystem::path& directory)
{
    std::ofstream obj(directory / "1.obj");
    obj << "v -0.1 0 0\nv 0.1 0 0\nv 0 0.1 0\n"
           "v 3.9 0 0\nv 4.1 0 0\nv 4 0.1 0\n"
           "f 1 2 3\nf 4 5 6\n";
    std::ofstream model(directory / "1.tsv");
    model << "PQSS_PRIMITIVE_BVH_MODEL_V1\n"
             "source\t1.obj\nsource_triangles\t2\nprimitive_count\t2\n"
             "primitive\t0\tsphere\t0\t0\t0\t1\t0\t0\t0\t1\t0\t0\t0\t1\t0\t0\t0.25\t0\n"
             "primitive\t1\tcapsule\t3.8\t0\t0\t1\t0\t0\t0\t0\t-1\t0\t1\t0\t0.4\t0\t0.25\t1\n";
}

void WriteModelTwo(const std::filesystem::path& directory)
{
    std::ofstream obj(directory / "2.obj");
    obj << "v -0.1 0 0\nv 0.1 0 0\nv 0 0.1 0\nf 1 2 3\n";
    std::ofstream model(directory / "2.tsv");
    model << "PQSS_PRIMITIVE_BVH_MODEL_V1\n"
             "source\t2.obj\nsource_triangles\t1\nprimitive_count\t1\n"
             "primitive\t0\trss\t-0.1\t-0.1\t0\t1\t0\t0\t0\t1\t0\t0\t0\t1\t0.2\t0.2\t0.05\t0\n";
}

} // namespace

int main()
{
    const auto suffix = std::chrono::steady_clock::now().time_since_epoch().count();
    const std::filesystem::path directory =
        std::filesystem::temp_directory_path() / ("pqss_primitive_bvh_" + std::to_string(suffix));
    std::filesystem::create_directories(directory);
    WriteModelOne(directory);
    WriteModelTwo(directory);
    {
        std::ofstream manifest(directory / "pool.tsv");
        manifest << "PQSS_PRIMITIVE_BVH_POOL_V1\nmodel_id\tmodel_file\n1\t1.tsv\n2\t2.tsv\n";
    }

    pqss_proxy_mesh::PrimitiveBvhPool pool;
    const auto cache_directory = directory / "cache";
    pool.loadAnalyzedPool(directory / "pool.tsv", pqss::BuildStrategy::Optimized, cache_directory);
    Check(pool.cacheStats().cache_hits == 0);
    Check(pool.cacheStats().cache_misses == 2);
    Check(pool.hasModel(1));
    Check(pool.hasModel(2));
    const auto& first_stats = pool.modelStats(1);
    Check(first_stats.analyzed);
    Check(first_stats.internal_containment_certified);
    Check(first_stats.source_triangles == 2);
    Check(first_stats.primitives == 2);
    Check(first_stats.nodes == 3);
    Check(first_stats.primitive_types[0] == 1);
    Check(first_stats.primitive_types[1] == 1);

    const pqss::Mat3 identity_R = pqss::Meident();
    const pqss::Vec3 zero_T = pqss::Veident();
    const auto sphere_rss = pool.query(1, identity_R, zero_T, 2, identity_R, zero_T, 0.0);
    Check(sphere_rss.closer_than_tolerance);
    Check(sphere_rss.leaf_pair_tests == 1);
    Check(sphere_rss.type_pair_tests[0][2] == 1);

    const pqss::Vec3 capsule_position = {4.0, 0.0, 0.0};
    const auto capsule_rss = pool.query(1, identity_R, zero_T, 2, identity_R, capsule_position, 0.0);
    Check(capsule_rss.closer_than_tolerance);
    Check(capsule_rss.type_pair_tests[1][2] == 1);

    const pqss::Vec3 far_position = {20.0, 0.0, 0.0};
    const auto separated = pool.query(1, identity_R, zero_T, 2, identity_R, far_position, 0.0);
    Check(!separated.closer_than_tolerance);
    Check(separated.root_bv_tests == 1);
    Check(separated.internal_bv_tests == 0);
    Check(separated.leaf_pair_tests == 0);

    pqss::ModelPool reference;
    Check(reference.SetSubdivisionEnabled(false) == pqss::Return::Success);
    Check(reference.NewModel(1, pqss::BuildStrategy::Optimized) == pqss::Return::Success);
    Check(reference.AddTri({-0.1, 0.0, 0.0}, {0.1, 0.0, 0.0}, {0.0, 0.1, 0.0}) ==
          pqss::Return::Success);
    Check(reference.EndModel() == pqss::Return::Success);
    Check(reference.EndPool() == pqss::Return::Success);
    pool.importReferencePool(reference, {99});
    Check(pool.hasModel(99));
    Check(!pool.modelStats(99).analyzed);
    const auto imported = pool.query(99, identity_R, zero_T, 2, identity_R, zero_T, 0.0);
    Check(imported.closer_than_tolerance);

    pqss_proxy_mesh::PrimitiveBvhPool cached_pool;
    cached_pool.loadAnalyzedPool(directory / "pool.tsv", pqss::BuildStrategy::Optimized, cache_directory);
    Check(cached_pool.cacheStats().cache_hits == 2);
    Check(cached_pool.cacheStats().cache_misses == 0);
    const auto cached_query = cached_pool.query(1, identity_R, zero_T, 2, identity_R, zero_T, 0.0);
    Check(cached_query.closer_than_tolerance);

    pqss::ModelPool direct_pool;
    Check(direct_pool.SetSubdivisionEnabled(false) == pqss::Return::Success);
    for (int model = 0; model < 2; ++model)
    {
        Check(direct_pool.NewModel(1, pqss::BuildStrategy::Optimized) == pqss::Return::Success);
        Check(direct_pool.AddTri({-0.1, 0.0, 0.0}, {0.1, 0.0, 0.0}, {0.0, 0.1, 0.0}) ==
              pqss::Return::Success);
        Check(direct_pool.EndModel() == pqss::Return::Success);
    }
    Check(direct_pool.EndPool() == pqss::Return::Success);
    cached_pool.installIntoPqssPool(direct_pool, {1, 2});
    pqss::QueryResult direct_collision;
    Check(direct_pool.Query(direct_collision, 0, identity_R, zero_T,
                            1, identity_R, zero_T, 0.0) == pqss::Return::Success);
    Check(direct_collision.CloserThanTolerance());
    Check(direct_collision.NumTriTests() == 1);
    pqss::QueryResult direct_separated;
    Check(direct_pool.Query(direct_separated, 0, identity_R, zero_T,
                            1, identity_R, far_position, 0.0) == pqss::Return::Success);
    Check(!direct_separated.CloserThanTolerance());

    std::filesystem::remove_all(directory);
    return 0;
}
