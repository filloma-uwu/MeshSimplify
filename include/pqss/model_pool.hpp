#pragma once

#include <array>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <numbers>
#include <optional>
#include <string>
#include <string_view>
#include <vector>

#include "pqss/bv.hpp"
#include "pqss/config/config.hpp"
#include "pqss/config/constant.hpp"
#include "pqss/model.hpp"
#include "pqss/query_result.hpp"

namespace pqss
{

struct RssBvhNode
{
    Mat3 R{};
    Vec3 T{};
    std::array<Real, 2> lengths{};
    Real radius = k_zero;
    Real size = k_zero;
    std::size_t left = 0;
    std::size_t right = 0;
    bool leaf = false;
};

class ModelPool
{
public:
    struct BuildCacheStats
    {
        std::size_t eligible_models   = 0;
        std::size_t cache_hits        = 0;
        std::size_t cache_misses      = 0;
        double      area_threshold_ms = 0;
        double      cache_load_ms     = 0;
        double      cache_build_ms    = 0;
        double      cache_save_ms     = 0;
    };

    struct ModelStats
    {
        std::size_t original_tri_count = 0;
        std::size_t built_tri_count    = 0;
        std::size_t subdivision_count  = 0;
        std::size_t bv_count           = 0;
    };

    struct BvhSignature
    {
        std::size_t   original_tri_count = 0;
        std::size_t   built_tri_count    = 0;
        std::size_t   subdivision_count  = 0;
        std::size_t   bv_count           = 0;
        std::size_t   leaf_count         = 0;
        std::size_t   internal_count     = 0;
        std::size_t   max_depth          = 0;
        Real          root_size          = k_zero;
        Real          total_size         = k_zero;
        Real          leaf_size_sum      = k_zero;
        Real          internal_size_sum  = k_zero;
        std::uint64_t topology_signature = 0;
        std::uint64_t geometry_signature = 0;
        std::uint64_t leaf_signature     = 0;
    };

    // ========================================================================
    // Build API
    // ========================================================================

    Return LoadModelFromFile(std::string_view filename, BuildStrategy build_strategy);

    Return Reserve(std::size_t num_models);
    Return ReserveTris(std::size_t num_tris);
    Return NewModel(std::size_t num_tris, BuildStrategy build_strategy = BuildStrategy::Optimized);
    Return AddTri(const Vec3& p1, const Vec3& p2, const Vec3& p3);
    Return ResetBuildStrategy(BuildStrategy build_strategy);
    Return EndModel();
    Return DeleteModel(std::size_t model_id);
    Return EndPool();
    Return Clear();
    Return SetSubdivisionEnabled(bool enabled);
    Return SetAreaThresholdOverride(Real area_threshold);
    Return SetB12vBuildEnabled(bool enabled);
    Return SetBuildCacheTag(std::string_view cache_tag);
    Return ReplaceModelBvhWithRssLeaves(std::size_t model_id, const std::vector<RssBvhNode>& nodes);

    // ========================================================================
    // Query API
    // ========================================================================

    bool                      Processed() const;
    std::size_t               PoolSize() const;
    std::size_t               Info(bool print = true) const;
    Real                      MaxArea() const;
    std::size_t               SubdivisionTimes() const;
    bool                      SubdivisionEnabled() const;
    bool                      B12vBuildEnabled() const;
    std::string_view          BuildCacheTag() const;
    const BuildCacheStats&    GetBuildCacheStats() const;
    std::vector<Real>         ModelAreaThresholds() const;
    std::vector<ModelStats>   GetModelStats() const;
    std::vector<BvhSignature> BvhSignatures() const;

    // Project-local read-only access for the independent-primitive scheduler.
    // The reference single-model Query path remains unchanged.
    const std::vector<build::Model>& Models() const;

    // ========================================================================
    // Collision detection
    // ========================================================================

    Return Query(QueryResult& res,
                 std::size_t  model_id_1,
                 const Mat3&  R1,
                 const Vec3&  T1,
                 std::size_t  model_id_2,
                 const Mat3&  R2,
                 const Vec3&  T2,
                 Real         tolerance) const;

private:
    // ---- Build internal types ----
    enum class PoolBuildState
    {
        Empty,
        Ready,
        Building,
        Processed
    };

    struct ModelSourceInfo
    {
        bool                  from_file = false;
        std::filesystem::path filename;
    };

    // ---- Build helpers ----
    Real   TriArea(const Vec3& p1, const Vec3& p2, const Vec3& p3) const;
    Real   ComputePoolMeanAreaThreshold() const;
    Real   ComputePoolAreaThreshold() const;
    Return RecomputeAreaThresholds();
    Return SetAppendedModelAreaThreshold(std::size_t model_index);
    Real   ModelAreaThreshold(std::size_t model_index) const;
    Return CountTriSubdivisions(const Vec3&  p1,
                                const Vec3&  p2,
                                const Vec3&  p3,
                                Real         area_threshold,
                                std::size_t& subdivision_counter) const;
    Return ComputeModelSubdivisionPlan(const build::Model& model,
                                       Real                area_threshold,
                                       std::size_t&        model_subdivision_times,
                                       std::uint64_t&      subdivision_signature) const;
    void   DiscardTemporaryBuildData(build::Model& model) const;
    Return
    AddTriRecurse(build::Model& model, const Vec3& p1, const Vec3& p2, const Vec3& p3, std::size_t original_tri_id);
    Return AddTriRecurseCounting(build::Model& model,
                                 const Vec3&   p1,
                                 const Vec3&   p2,
                                 const Vec3&   p3,
                                 std::size_t   original_tri_id,
                                 Real          area_threshold,
                                 std::size_t&  subdivision_counter) const;
    Return BuildModelDataWithCurrentThreshold(build::Model& model,
                                              Real          area_threshold,
                                              std::size_t&  model_subdivision_times) const;
    Return BuildModelDataWithoutSubdivision(build::Model& model, std::size_t& model_subdivision_times) const;
    Return BuildModelWithCurrentThreshold(build::Model&          model,
                                          const ModelSourceInfo& source,
                                          Real                   area_threshold,
                                          std::size_t&           model_subdivision_times);
    Return BuildAppendedModel(std::size_t model_index);
    bool   TryLoadBuildCache(build::Model&          model,
                             const ModelSourceInfo& source,
                             std::size_t            model_subdivision_times,
                             std::uint64_t          subdivision_signature);
    bool   SaveBuildCache(const build::Model&    model,
                          const ModelSourceInfo& source,
                          std::size_t            model_subdivision_times,
                          std::uint64_t          subdivision_signature) const;

    // ---- Query helpers ----
    void QueryRecurse(QueryResult&        res,
                      const Mat3&         R,
                      const Vec3&         T,
                      const build::Model& o1,
                      std::size_t         b1,
                      const build::Model& o2,
                      std::size_t         b2,
                      bool                rss_leaf_collision) const;

    // ---- Members ----
    PoolBuildState               m_build_state = PoolBuildState::Empty;
    std::vector<build::Model>    m_models;
    std::vector<ModelSourceInfo> m_model_sources;
    std::vector<Real>            m_model_area_thresholds;
    std::vector<std::uint8_t>    m_rss_leaf_collision_models;
    Real                         m_area_threshold                = k_zero;
    Real                         m_area_threshold_override       = k_zero;
    std::size_t                  m_subdivision_times             = 0;
    bool                         m_append_with_current_threshold = false;
    bool                         m_subdivision_enabled           = true;
    bool                         m_build_b12v_enabled            = true;
    std::string                  m_build_cache_tag               = "bvh";
    BuildCacheStats              m_build_cache_stats;
};

} // namespace pqss
