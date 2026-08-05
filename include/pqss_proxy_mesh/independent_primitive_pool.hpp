#pragma once

#include <cstddef>
#include <filesystem>
#include <memory>
#include <string>

#include "pqss/model_pool.hpp"

namespace pqss_proxy_mesh
{

enum class DistanceFilterMode
{
    None,
    Containment,
    SeparatingAxis,
    Hybrid,
};

[[nodiscard]] const char* distanceFilterModeName(DistanceFilterMode mode);

struct IndependentQueryResult
{
    bool closer_than_tolerance = false;
    std::size_t primitive_pair_tests = 0;
    std::size_t internal_bv_tests = 0;
    std::size_t triangle_tests = 0;
    std::size_t distance_filter_skips = 0;
    std::size_t containment_filter_checks = 0;
    std::size_t containment_filter_skips = 0;
    std::size_t axis_filter_checks = 0;
    std::size_t axis_filter_skips = 0;
    std::size_t axis_generation_checks = 0;
    std::size_t cached_separating_axes = 0;
    double query_time_ms = 0.0;

    [[nodiscard]] std::size_t totalBvTests() const
    {
        return primitive_pair_tests + internal_bv_tests;
    }
};

struct ContainmentStats
{
    std::size_t physical_models = 0;
    std::size_t bvh_nodes = 0;
    std::size_t tree_ancestor_links = 0;
    std::size_t cross_tree_containment_links = 0;
};

class IndependentPrimitivePool
{
public:
    IndependentPrimitivePool();
    ~IndependentPrimitivePool();
    IndependentPrimitivePool(const IndependentPrimitivePool&) = delete;
    IndependentPrimitivePool& operator=(const IndependentPrimitivePool&) = delete;
    IndependentPrimitivePool(IndependentPrimitivePool&&) noexcept;
    IndependentPrimitivePool& operator=(IndependentPrimitivePool&&) noexcept;

    void setBuildCacheTag(const std::string& cache_tag);
    void loadBaseManifest(const std::filesystem::path& manifest,
                          pqss::BuildStrategy strategy = pqss::BuildStrategy::Optimized);
    void endBasePool();
    void appendLogicalModel(std::size_t logical_pool_index,
                            const std::filesystem::path& model,
                            pqss::BuildStrategy strategy = pqss::BuildStrategy::Optimized);
    void finalizeContainment();

    [[nodiscard]] IndependentQueryResult query(std::size_t logical_model_1,
                                               const pqss::Mat3& R1,
                                               const pqss::Vec3& T1,
                                               std::size_t logical_model_2,
                                               const pqss::Mat3& R2,
                                               const pqss::Vec3& T2,
                                               pqss::Real tolerance,
                                               DistanceFilterMode filter_mode = DistanceFilterMode::SeparatingAxis) const;

    [[nodiscard]] const pqss::ModelPool& physicalPool() const;
    [[nodiscard]] const ContainmentStats& containmentStats() const;
    [[nodiscard]] std::size_t logicalModelCount() const;
    [[nodiscard]] std::size_t primitiveCount(std::size_t logical_model) const;
    [[nodiscard]] std::size_t physicalModelId(std::size_t logical_model,
                                              std::size_t primitive_index) const;

private:
    struct Impl;
    std::unique_ptr<Impl> m_impl;
};

} // namespace pqss_proxy_mesh
