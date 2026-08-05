#pragma once

#include <array>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <string>
#include <vector>

namespace pqss_proxy_mesh
{

struct PrimitiveMeshAnalysisOptions
{
    bool allow_polygon = true;
    bool allow_frustum = true;
    std::uint32_t frustum_segments = 24;
    // Legacy sweep control. The production path uses one structure policy.
    double analysis_strength = 0.5;
    // A non-negative value overrides the threshold derived from analysis_strength.
    double maximum_added_volume_ratio = -1.0;
    double coplanar_relative_tolerance = 1.0e-9;
    double circle_radial_tolerance = 0.01;
    // CAD exports often quantize vertices from one analytic plane independently.
    // A connected surface may therefore miss exact coplanarity by a few parts
    // per million of the model diagonal. This tolerance is used only after a
    // complete connected component passes a best-fit plane certificate.
    double analytic_surface_relative_tolerance = 1.0e-5;
    // Conservatively absorb small boundary details into the surrounding planar
    // proxy. Added area is normalized by the entire input model surface area.
    double tiny_planar_detail_area_ratio = 0.002;
    std::uint32_t tiny_planar_detail_max_vertices = 4;
    // Disconnected coplanar patches may be joined when their added planar area,
    // normalized by the entire input model surface area, is below this ratio.
    double maximum_local_planar_fill_area_ratio = 0.03;
    // Whole-model projection extrusion is diagnostic only. It destroys local
    // depth and must never be selected as the production proxy by default.
    bool enable_volume_evaluated_envelope = false;
    // Per-cavity conservative prism volume, normalized by the model AABB.
    double maximum_cavity_added_volume_ratio = 0.04;
    // Total volume introduced when several disconnected components are replaced
    // by one shared box envelope, normalized by the model AABB.
    double maximum_group_box_added_volume_ratio = 0.05;
    // Added surface of a complete local box envelope, normalized by the whole
    // model surface. A slightly wider budget than a one-sided protrusion is
    // appropriate because the candidate replaces an entire certified volume.
    double maximum_box_envelope_added_area_ratio = 0.04;
    // Maximum volume added by a candidate envelope, normalized by the model AABB.
    double maximum_envelope_added_volume_ratio = 0.40;
    // A fill with substantial volume error must buy a substantial BVH workload reduction.
    double minimum_envelope_primitive_reduction_ratio = 4.0;
    std::uint32_t projection_envelope_resolution = 192;
    double shallow_concavity_closing_ratio = 0.03;
    bool recognize_support_protrusions = true;
    // Added box surface area normalized by the entire input model surface area.
    double protrusion_max_area_excess_ratio = 0.03;
    double protrusion_cluster_gap_relative = 0.01;
    double protrusion_max_inward_relative = 0.02;
    // Collapse a certified thin shell onto its outward support plane while
    // preserving the exact projected polygon union.
    double shallow_parallel_merge_depth_relative = 0.005;
    // After a thin shell is coalesced, absorb minor silhouette steps while
    // retaining its dominant recess. Normalized by the whole-model surface.
    double shallow_shell_silhouette_added_area_ratio = 0.005;
    bool uniform_structure_policy = true;
    // Production pipeline: repair holes, recognize complete surfaces, perform
    // error-bounded fixed-point surface merging, clean the outer boundary, and
    // triangulate. Disable only to exercise the retired unified-frontier path
    // in compatibility tests.
    bool use_staged_surface_pipeline = true;
    // Replacement architecture: all surface and enclosing-volume candidates
    // compete in one hierarchy instead of greedily claiming source faces in
    // geometry-specific passes.
    bool use_unified_candidate_optimizer = false;
    // A complete closed-volume certificate enables internal-overlap removal.
    // Treat full-model certified responsibility as this many saved proxy
    // triangles, but keep the benefit bounded so it can never select a huge
    // fragmented shell merely because its approximation error is zero.
    double unified_enclosure_workload_credit = 64.0;
    // In the uniform policy a high-error whole-body box remains forbidden.
    // Geometry-connected local structures may still merge into a six-polygon
    // box when both globally normalized added-volume and added-area budgets pass.
    bool forbid_main_body_box_approximation = true;
    // Hard one-sided proxy-to-source distance limit in model units. A negative
    // value selects the default of 8% of the model AABB diagonal.
    double maximum_open_error_distance = -1.0;
    // Debug-only bounded export of failed unified leaf candidates. Production
    // analysis keeps this disabled to avoid extra mesh copies and disk output.
    bool write_unified_region_diagnostics = false;
};

struct PrimitiveMeshAnalysisStats
{
    std::size_t source_triangles = 0;
    std::size_t discarded_degenerate_triangles = 0;
    std::size_t primitive_count = 0;
    std::size_t polygon_count = 0;
    std::size_t frustum_count = 0;
    std::size_t disk_count = 0;
    std::size_t annulus_count = 0;
    std::size_t cylindrical_band_count = 0;
    std::size_t conical_band_count = 0;
    std::size_t proxy_triangles = 0;
    std::size_t filled_planar_holes = 0;
    double filled_cavity_volume_ratio = 0.0;
    std::size_t filled_boundary_voids = 0;
    double filled_boundary_void_area = 0.0;
    std::size_t removed_contained_primitives = 0;
    std::size_t removed_sealed_void_wall_primitives = 0;
    std::size_t excluded_sealed_void_wall_triangles = 0;
    std::size_t removed_blind_cavity_primitives = 0;
    std::size_t recognized_closed_box_shells = 0;
    std::size_t recognized_protrusion_box_shells = 0;
    std::size_t merged_local_planar_primitives = 0;
    std::size_t canonicalized_coplanar_groups = 0;
    std::size_t removed_coplanar_redundant_primitives = 0;
    double removed_coplanar_overlap_area = 0.0;
    double minimum_protrusion_candidate_area_excess_ratio = -1.0;
    double selected_envelope_added_volume_ratio = -1.0;
    int selected_envelope_axis = -1;
    std::array<double, 3> envelope_candidate_added_volume_ratios{{-1.0, -1.0, -1.0}};
    std::array<std::size_t, 3> envelope_candidate_primitive_counts{{0, 0, 0}};
    std::array<std::size_t, 3> envelope_candidate_remaining_cavities{{0, 0, 0}};
    std::size_t excluded_redundant_triangles = 0;
    // Unified-candidate optimizer telemetry. Counts and inclusive timings are
    // persisted so expensive offline analyses can be diagnosed without
    // inferring bottlenecks from total wall time.
    std::size_t unified_approximate_planar_regions = 0;
    std::size_t unified_exact_coplanar_regions = 0;
    std::size_t unified_hierarchy_regions = 0;
    std::size_t unified_connected_components = 0;
    std::size_t unified_lateral_sweep_axes = 0;
    std::size_t unified_lateral_sweep_patches = 0;
    std::size_t unified_lateral_sweep_certified = 0;
    std::size_t unified_lateral_sweep_candidates = 0;
    std::size_t unified_lateral_sweep_max_faces = 0;
    std::size_t unified_lateral_significant_candidates = 0;
    double unified_lateral_sweep_max_axial_aspect = 0.0;
    double unified_lateral_sweep_max_area_ratio = 0.0;
    std::size_t unified_lateral_reject_geometry = 0;
    std::size_t unified_lateral_reject_end_ring = 0;
    std::size_t unified_lateral_reject_normal = 0;
    std::size_t unified_lateral_reject_radial = 0;
    std::size_t unified_lateral_best_end_ring_min_points = 0;
    double unified_lateral_best_end_ring_spread_ratio = -1.0;
    double unified_lateral_best_end_ring_gap_ratio = -1.0;
    double unified_lateral_best_end_ring_score = -1.0;
    std::size_t unified_final_frontier_size = 0;
    // True only when the error-constrained optimizer has no admissible
    // candidate and returns the source triangles verbatim.  The export path
    // must keep this representation linear; running the global coplanar union
    // over thousands of one-face primitives is both unnecessary and can have
    // quadratic memory cost.
    bool unified_exact_fallback_selected = false;
    std::size_t unified_final_analytic_choices = 0;
    double unified_selected_workload = -1.0;
    double unified_selected_static_workload = -1.0;
    std::size_t unified_selected_enclosure_overlap_pairs = 0;
    double unified_selected_enclosure_overlap_score = 0.0;
    std::size_t unified_selected_analytic_faces = 0;
    double unified_selected_analytic_area_ratio = 0.0;
    double unified_best_analytic_workload = -1.0;
    std::size_t unified_best_analytic_faces = 0;
    double unified_best_analytic_area_ratio = 0.0;
    double unified_whole_model_box_workload = -1.0;
    bool unified_whole_model_box_passed_local_budget = false;
    std::size_t unified_adaptive_refinements = 0;
    std::size_t unified_skipped_final_audits = 0;
    std::size_t unified_skipped_industrial_fallback_candidates = 0;
    std::vector<double> unified_final_frontier_workloads;
    std::size_t unified_enclosure_extrusions = 0;
    std::size_t unified_occlusion_clipped_primitives = 0;
    std::size_t unified_occlusion_removed_primitives = 0;
    std::size_t unified_occlusion_input_triangles = 0;
    std::size_t unified_occlusion_output_triangles = 0;
    double unified_occlusion_removed_area = 0.0;
    bool unified_occlusion_rolled_back = false;
    double unified_occlusion_seconds = 0.0;
    std::size_t unified_classify_calls = 0;
    std::size_t unified_classify_faces = 0;
    std::size_t unified_frontier_prune_calls = 0;
    std::size_t unified_frontier_prune_choices = 0;
    std::size_t unified_frontier_prune_max_choices = 0;
    std::uint64_t unified_frontier_prune_comparisons = 0;
    double unified_planar_region_seconds = 0.0;
    double unified_connected_component_seconds = 0.0;
    double unified_analytic_recognition_seconds = 0.0;
    double unified_classify_seconds = 0.0;
    double unified_frontier_prune_seconds = 0.0;
    double unified_hierarchy_seconds = 0.0;
    // Final export gate. Local candidates are certified while they are built;
    // these counters verify that post-selection containment removal and planar
    // canonicalization still conservatively cover every source triangle.
    std::size_t coverage_assigned_source_faces = 0;
    std::size_t coverage_enclosure_source_faces = 0;
    std::size_t coverage_planar_source_faces = 0;
    std::size_t coverage_unassigned_source_faces = 0;
    std::size_t coverage_failed_source_faces = 0;
    double coverage_audit_seconds = 0.0;
    // Final, model-independent validation of the exported proxy. Distance is
    // one-sided from proxy to source and excludes certified planar hole fills.
    bool containment_validation_passed = false;
    std::size_t open_error_distance_sample_count = 0;
    double open_mean_distance = 0.0;
    double open_max_distance = 0.0;
    double open_mean_distance_ratio = 0.0;
    double open_max_distance_ratio = 0.0;
    std::array<double, 3> open_max_proxy_point{{0.0, 0.0, 0.0}};
    std::array<double, 3> open_max_source_point{{0.0, 0.0, 0.0}};
    double open_error_audit_seconds = 0.0;
    double maximum_open_error_distance_limit = 0.0;
    double analysis_seconds = 0.0;
};

// Maps [0, 1] to the extended non-negative reals. The formula itself yields
// positive infinity at zero; zero does not select a separate analysis path.
[[nodiscard]] double analysisStrengthToAllowedExcessRatio(double strength);

[[nodiscard]] PrimitiveMeshAnalysisStats analyzePrimitiveMeshObj(
    const std::filesystem::path& input_obj,
    const std::filesystem::path& output_directory,
    const PrimitiveMeshAnalysisOptions& options);

// Builds the expensive merge hierarchy once, then emits strengths 0.00 through
// 1.00 at 0.01 intervals into s000 ... s100.
[[nodiscard]] std::vector<PrimitiveMeshAnalysisStats> analyzePrimitiveMeshObjStrengthSweep(
    const std::filesystem::path& input_obj,
    const std::filesystem::path& output_root,
    const PrimitiveMeshAnalysisOptions& options);

} // namespace pqss_proxy_mesh
