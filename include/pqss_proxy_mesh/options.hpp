#pragma once

#include <cstddef>
#include <cstdint>
#include <optional>
#include <string>
#include <vector>

namespace pqss_proxy_mesh
{

enum class SimplificationMode
{
    ConservativeOuter
};

struct SimplificationOptions
{
    SimplificationMode mode = SimplificationMode::ConservativeOuter;
    std::uint32_t max_pqss_bvh_depth = 8;
    double relative_alpha = 0.02;
    double relative_offset = 0.005;
    bool allow_component_merge = false;
    bool allow_box_proxy = true;
};

// One source model produces one proxy model at the same pool index. Models
// are never merged across pool boundaries.
struct ModelProxyMetrics
{
    std::string model_name;
    std::size_t source_faces = 0;
    std::size_t proxy_faces = 0;
    std::size_t pqss_built_faces = 0;
    std::size_t pqss_subdivision_count = 0;
    std::size_t bvh_nodes = 0;
    std::uint32_t actual_max_depth = 0;
    double root_rss_size = 0.0;
    double total_internal_rss_size = 0.0;
    double sibling_overlap_score = 0.0;
    bool containment_certified = false;
};

struct ProxyPoolMetrics
{
    std::uint32_t requested_max_depth = 0;
    double pqss_pool_area_threshold = 0.0;
    std::size_t total_source_faces = 0;
    std::size_t total_proxy_faces = 0;
    std::size_t total_pqss_built_faces = 0;
    std::size_t total_pqss_subdivision_count = 0;
    double representative_query_time_ms = 0.0;
    std::vector<ModelProxyMetrics> models;
};

[[nodiscard]] std::optional<std::string> validateOptions(
    const SimplificationOptions& options);

} // namespace pqss_proxy_mesh
