#include <array>
#include <cassert>
#include <chrono>
#include <filesystem>
#include <fstream>
#include <string>

#include "pqss/model_pool.hpp"
#include "pqss/utils/mat_vec.hpp"
#include "pqss_proxy_mesh/independent_primitive_pool.hpp"

namespace
{

void WriteCube(const std::filesystem::path& path, const double half_extent)
{
    std::ofstream stream(path);
    assert(stream);
    for (const auto& vertex : std::array{
             std::array{-half_extent, -half_extent, -half_extent},
             std::array{ half_extent, -half_extent, -half_extent},
             std::array{ half_extent,  half_extent, -half_extent},
             std::array{-half_extent,  half_extent, -half_extent},
             std::array{-half_extent, -half_extent,  half_extent},
             std::array{ half_extent, -half_extent,  half_extent},
             std::array{ half_extent,  half_extent,  half_extent},
             std::array{-half_extent,  half_extent,  half_extent},
         })
        stream << "v " << vertex[0] << ' ' << vertex[1] << ' ' << vertex[2] << '\n';
    for (const auto& face : std::array{
             std::array{1, 3, 2}, std::array{1, 4, 3},
             std::array{5, 6, 7}, std::array{5, 7, 8},
             std::array{1, 2, 6}, std::array{1, 6, 5},
             std::array{4, 8, 7}, std::array{4, 7, 3},
             std::array{1, 5, 8}, std::array{1, 8, 4},
             std::array{2, 3, 7}, std::array{2, 7, 6},
         })
        stream << "f " << face[0] << ' ' << face[1] << ' ' << face[2] << '\n';
}

} // namespace

int main()
{
    const auto suffix = std::chrono::steady_clock::now().time_since_epoch().count();
    const std::filesystem::path directory =
        std::filesystem::temp_directory_path() / ("pqss_independent_pool_test_" + std::to_string(suffix));
    std::filesystem::create_directories(directory);
    const auto outer = directory / "outer.obj";
    const auto inner = directory / "inner.obj";
    WriteCube(outer, 2.0);
    WriteCube(inner, 0.5);

    const auto manifest = directory / "primitive_pool_manifest.tsv";
    {
        std::ofstream stream(manifest);
        stream << "logical_model_id\tpart_index\trelative_path\tselected\n";
        stream << "1\t0\touter.obj\t1\n";
        stream << "1\t1\tinner.obj\t1\n";
        for (int id = 2; id <= 23; ++id)
            stream << id << "\t0\tinner.obj\t0\n";
    }

    pqss_proxy_mesh::IndependentPrimitivePool pool;
    pool.setBuildCacheTag("independent_primitive_pool_test_" + std::to_string(suffix));
    pool.loadBaseManifest(manifest);
    pool.endBasePool();
    pool.finalizeContainment();
    assert(pool.primitiveCount(0) == 2);
    assert(pool.physicalModelId(0, 0) == 0);
    assert(pool.physicalModelId(0, 1) == 1);
    assert(pool.physicalModelId(1, 0) == 2);
    assert(pool.containmentStats().cross_tree_containment_links > 0);

    const pqss::Mat3 identity_R = pqss::Meident();
    const pqss::Vec3 origin = pqss::Veident();
    const pqss::Vec3 far = {100.0, 0.0, 0.0};

    const auto wrapped = pool.query(1, identity_R, origin, 2, identity_R, origin, 0.0);
    pqss::QueryResult reference;
    const auto status = pool.physicalPool().Query(
        reference,
        pool.physicalModelId(1, 0), identity_R, origin,
        pool.physicalModelId(2, 0), identity_R, origin,
        0.0);
    assert(status == pqss::Return::Success);
    assert(wrapped.closer_than_tolerance == reference.CloserThanTolerance());
    assert(wrapped.internal_bv_tests == reference.NumBvTests());
    assert(wrapped.triangle_tests == reference.NumTriTests());

    const auto unfiltered = pool.query(
        0, identity_R, origin, 1, identity_R, far, 0.0,
        pqss_proxy_mesh::DistanceFilterMode::None);
    const auto containment = pool.query(
        0, identity_R, origin, 1, identity_R, far, 0.0,
        pqss_proxy_mesh::DistanceFilterMode::Containment);
    const auto separating_axis = pool.query(
        0, identity_R, origin, 1, identity_R, far, 0.0,
        pqss_proxy_mesh::DistanceFilterMode::SeparatingAxis);
    const auto hybrid = pool.query(
        0, identity_R, origin, 1, identity_R, far, 0.0,
        pqss_proxy_mesh::DistanceFilterMode::Hybrid);
    assert(!unfiltered.closer_than_tolerance);
    assert(unfiltered.distance_filter_skips == 0);
    assert(unfiltered.primitive_pair_tests == 2);
    assert(containment.closer_than_tolerance == unfiltered.closer_than_tolerance);
    assert(containment.containment_filter_skips > 0);
    assert(containment.primitive_pair_tests < unfiltered.primitive_pair_tests);
    assert(separating_axis.closer_than_tolerance == unfiltered.closer_than_tolerance);
    assert(separating_axis.axis_filter_skips > 0);
    assert(separating_axis.primitive_pair_tests < unfiltered.primitive_pair_tests);
    assert(hybrid.closer_than_tolerance == unfiltered.closer_than_tolerance);
    assert(hybrid.distance_filter_skips > 0);

    const auto early_stop = pool.query(
        0, identity_R, origin, 1, identity_R, origin, 0.0,
        pqss_proxy_mesh::DistanceFilterMode::None);
    assert(early_stop.closer_than_tolerance);
    assert(early_stop.primitive_pair_tests == 1);

    std::filesystem::remove_all(directory);
    return 0;
}
