#pragma once

#include <array>
#include <cstddef>
#include <filesystem>
#include <memory>
#include <vector>

#include "pqss/model_pool.hpp"

namespace pqss_proxy_mesh
{

enum class AnalyticPrimitiveType : std::size_t
{
    Sphere = 0,
    Capsule = 1,
    Rss = 2,
    Count = 3,
};

struct PrimitiveBvhQueryResult
{
    bool closer_than_tolerance = false;
    std::size_t root_bv_tests = 0;
    std::size_t internal_bv_tests = 0;
    std::size_t leaf_pair_tests = 0;
    std::array<std::array<std::size_t, 3>, 3> type_pair_tests{};
    double query_time_ms = 0.0;

    [[nodiscard]] std::size_t totalBvTests() const
    {
        return root_bv_tests + internal_bv_tests;
    }
};

struct PrimitiveBvhModelStats
{
    std::size_t model_id = 0;
    std::size_t source_triangles = 0;
    std::size_t primitives = 0;
    std::size_t nodes = 0;
    std::size_t max_depth = 0;
    std::array<std::size_t, 3> primitive_types{};
    bool analyzed = false;
    bool internal_containment_certified = false;
};

struct PrimitiveBvhCacheStats
{
    std::size_t eligible_models = 0;
    std::size_t cache_hits = 0;
    std::size_t cache_misses = 0;
    double cache_load_ms = 0.0;
    double cache_build_ms = 0.0;
    double cache_save_ms = 0.0;
};

class PrimitiveBvhPool
{
public:
    PrimitiveBvhPool();
    ~PrimitiveBvhPool();
    PrimitiveBvhPool(const PrimitiveBvhPool&) = delete;
    PrimitiveBvhPool& operator=(const PrimitiveBvhPool&) = delete;
    PrimitiveBvhPool(PrimitiveBvhPool&&) noexcept;
    PrimitiveBvhPool& operator=(PrimitiveBvhPool&&) noexcept;

    void importReferencePool(const pqss::ModelPool& reference_pool,
                             const std::vector<std::size_t>& object_ids);
    void loadAnalyzedPool(const std::filesystem::path& manifest,
                          pqss::BuildStrategy internal_fit_strategy = pqss::BuildStrategy::Optimized,
                          const std::filesystem::path& cache_directory = {});

    [[nodiscard]] PrimitiveBvhQueryResult query(std::size_t model_id_1,
                                                const pqss::Mat3& R1,
                                                const pqss::Vec3& T1,
                                                std::size_t model_id_2,
                                                const pqss::Mat3& R2,
                                                const pqss::Vec3& T2,
                                                pqss::Real tolerance) const;

    [[nodiscard]] bool hasModel(std::size_t model_id) const;
    [[nodiscard]] const PrimitiveBvhModelStats& modelStats(std::size_t model_id) const;
    [[nodiscard]] std::vector<std::size_t> modelIds() const;
    [[nodiscard]] const PrimitiveBvhCacheStats& cacheStats() const;
    void installIntoPqssPool(pqss::ModelPool& pool,
                             const std::vector<std::size_t>& object_ids) const;

private:
    struct Impl;
    std::unique_ptr<Impl> m_impl;
};

[[nodiscard]] const char* analyticPrimitiveTypeName(AnalyticPrimitiveType type);

} // namespace pqss_proxy_mesh
