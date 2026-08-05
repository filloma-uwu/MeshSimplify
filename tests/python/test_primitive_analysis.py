from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

import numpy as np
from shapely import Polygon
import trimesh

from pqss_proxy_mesh.obj_io import Mesh, write_obj
from pqss_proxy_mesh.primitive_analysis import (
    Primitive,
    PrimitiveAnalysisOptions,
    PrimitiveRegion,
    _exclude_covered_planar_regions,
    _fill_short_boundary_voids,
    _point_core_distances,
    analyze_directory,
    analyze_mesh,
)
from pqss_proxy_mesh.primitive_bvh_export import export_pool


def _tetrahedron() -> Mesh:
    return Mesh(
        "1.obj",
        np.asarray(
            [
                (0.0, 0.0, 0.0),
                (2.0, 0.0, 0.0),
                (0.0, 1.0, 0.0),
                (0.0, 0.0, 0.5),
            ],
            dtype=np.float64,
        ),
        np.asarray(
            [(0, 2, 1), (0, 1, 3), (0, 3, 2), (1, 2, 3)],
            dtype=np.int64,
        ),
    )


def _bent_strip() -> Mesh:
    return Mesh(
        "strip.obj",
        np.asarray(
            [
                (0.0, 0.0, 0.0),
                (1.0, 0.0, 0.0),
                (0.0, 1.0, 0.0),
                (1.0, 1.0, 0.0),
                (0.0, 1.0, 1.0),
                (1.0, 1.0, 1.0),
            ],
            dtype=np.float64,
        ),
        np.asarray([(0, 1, 2), (1, 3, 2), (2, 3, 4), (3, 5, 4)], dtype=np.int64),
    )


def _subdivided_rectangle() -> Mesh:
    return Mesh(
        "rectangle.obj",
        np.asarray(
            [
                (0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (2.0, 0.0, 0.0),
                (0.0, 1.0, 0.0), (1.0, 1.0, 0.0), (2.0, 1.0, 0.0),
            ],
            dtype=np.float64,
        ),
        np.asarray(
            [(0, 1, 3), (1, 4, 3), (1, 2, 4), (2, 5, 4)],
            dtype=np.int64,
        ),
    )


def _subdivided_l_shape() -> Mesh:
    return Mesh(
        "l_shape.obj",
        np.asarray(
            [
                (0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (2.0, 0.0, 0.0),
                (0.0, 1.0, 0.0), (1.0, 1.0, 0.0), (2.0, 1.0, 0.0),
                (0.0, 2.0, 0.0), (1.0, 2.0, 0.0),
            ],
            dtype=np.float64,
        ),
        np.asarray(
            [
                (0, 1, 3), (1, 4, 3),
                (1, 2, 4), (2, 5, 4),
                (3, 4, 6), (4, 7, 6),
            ],
            dtype=np.int64,
        ),
    )


def _square_ring() -> Mesh:
    vertices = np.asarray(
        [(float(x), float(y), 0.0) for y in range(4) for x in range(4)],
        dtype=np.float64,
    )
    faces: list[tuple[int, int, int]] = []
    for y in range(3):
        for x in range(3):
            if x == 1 and y == 1:
                continue
            lower_left = 4 * y + x
            faces.extend(
                (
                    (lower_left, lower_left + 1, lower_left + 4),
                    (lower_left + 1, lower_left + 5, lower_left + 4),
                )
            )
    return Mesh("ring.obj", vertices, np.asarray(faces, dtype=np.int64))


def _box_frame() -> Mesh:
    parts = []
    for extents, translation in (
        ((4.0, 1.0, 0.2), (0.0, -1.5, 0.0)),
        ((4.0, 1.0, 0.2), (0.0, 1.5, 0.0)),
        ((1.0, 2.0, 0.2), (-1.5, 0.0, 0.0)),
        ((1.0, 2.0, 0.2), (1.5, 0.0, 0.0)),
    ):
        part = trimesh.creation.box(extents=extents)
        part.apply_translation(translation)
        parts.append(part)
    combined = trimesh.util.concatenate(parts)
    return Mesh(
        "frame.obj",
        np.asarray(combined.vertices, dtype=np.float64),
        np.asarray(combined.faces, dtype=np.int64),
    )


class PrimitiveAnalysisTest(unittest.TestCase):
    def test_final_regions_are_an_exact_triangle_partition(self) -> None:
        mesh = _bent_strip()
        result = analyze_mesh(
            mesh,
            PrimitiveAnalysisOptions(max_excess_volume_ratio=0.0),
        )
        assigned = np.concatenate([region.triangle_ids for region in result.regions])
        self.assertEqual(len(assigned), len(mesh.faces))
        np.testing.assert_array_equal(np.sort(assigned), np.arange(len(mesh.faces)))
        self.assertTrue(result.stats["strict_triangle_partition"])

    def test_sharp_planar_faces_do_not_merge_across_right_angle(self) -> None:
        result = analyze_mesh(
            _bent_strip(),
            PrimitiveAnalysisOptions(analysis_strength=0.0, primitive_types=("rss",)),
        )
        self.assertEqual(len(result.regions), 2)

    def test_fully_covered_coplanar_region_is_excluded(self) -> None:
        mesh = Mesh(
            "covered.obj",
            np.asarray(
                [
                    (0.0, 0.0, 0.0), (4.0, 0.0, 0.0),
                    (0.0, 4.0, 0.0), (4.0, 4.0, 0.0),
                    (1.0, 1.0, 0.0), (2.0, 1.0, 0.0),
                    (1.0, 2.0, 0.0), (2.0, 2.0, 0.0),
                ],
                dtype=np.float64,
            ),
            np.asarray(
                [(0, 1, 2), (1, 3, 2), (4, 5, 6), (5, 7, 6)],
                dtype=np.int64,
            ),
        )
        primitive = Primitive(
            "rss", np.zeros(3), np.eye(3), np.asarray((4.0, 4.0)), 1.0e-6, 0.0
        )
        regions = [
            PrimitiveRegion(primitive, np.asarray((0, 1)), np.asarray((0, 1, 2, 3))),
            PrimitiveRegion(primitive, np.asarray((2, 3)), np.asarray((4, 5, 6, 7))),
        ]
        filtered, stats, excluded = _exclude_covered_planar_regions(
            mesh, regions, PrimitiveAnalysisOptions()
        )
        self.assertEqual(len(filtered), 1)
        self.assertEqual(stats["removed_covered_planar_primitives"], 1)
        np.testing.assert_array_equal(excluded, np.asarray((2, 3)))

    def test_each_supported_type_strictly_contains_its_region_vertices(self) -> None:
        mesh = _tetrahedron()
        for kind in ("sphere", "capsule", "rss"):
            with self.subTest(kind=kind):
                result = analyze_mesh(
                    mesh,
                    PrimitiveAnalysisOptions(
                        primitive_types=(kind,),
                        analysis_strength=0.0,
                    ),
                )
                self.assertEqual({region.primitive.kind for region in result.regions}, {kind})
                for region in result.regions:
                    points = mesh.vertices[region.vertex_ids]
                    distances = _point_core_distances(points, region.primitive)
                    self.assertTrue(np.all(distances < region.primitive.radius))

    def test_type_filter_never_emits_a_disabled_type(self) -> None:
        result = analyze_mesh(
            _tetrahedron(),
            PrimitiveAnalysisOptions(
                primitive_types=("sphere", "rss"),
                max_excess_volume_ratio=0.0,
            ),
        )

    def test_degenerate_triangles_are_removed_before_analysis(self) -> None:
        mesh = _tetrahedron()
        vertices = np.vstack((mesh.vertices, np.asarray([(3.0, 3.0, 3.0)] * 3)))
        faces = np.vstack((mesh.faces, np.asarray([(4, 5, 6)], dtype=np.int64)))
        result = analyze_mesh(Mesh("degenerate.obj", vertices, faces), PrimitiveAnalysisOptions())
        self.assertEqual(result.stats["original_source_triangles"], 5)
        self.assertEqual(result.stats["source_triangles"], 4)
        self.assertEqual(result.stats["discarded_degenerate_triangles"], 1)
        assigned = np.concatenate([region.triangle_ids for region in result.regions])
        np.testing.assert_array_equal(np.sort(assigned), np.arange(4))

    def test_near_collinear_faces_are_removed_but_small_well_shaped_faces_are_kept(self) -> None:
        mesh = Mesh(
            "quality.obj",
            np.asarray(
                [
                    (0.0, 0.0, 0.0),
                    (1.0e-20, 0.0, 0.0),
                    (0.0, 1.0e-20, 0.0),
                    (10.0, 0.0, 0.0),
                    (11.0, 0.0, 0.0),
                    (12.0, 1.0e-16, 0.0),
                ],
                dtype=np.float64,
            ),
            np.asarray([(0, 1, 2), (3, 4, 5)], dtype=np.int64),
        )
        result = analyze_mesh(mesh, PrimitiveAnalysisOptions())
        self.assertEqual(result.stats["source_triangles"], 1)
        self.assertEqual(result.stats["discarded_degenerate_triangles"], 1)

    def test_outward_deviation_metrics_are_finite_and_ordered(self) -> None:
        result = analyze_mesh(_tetrahedron(), PrimitiveAnalysisOptions())
        self.assertGreaterEqual(result.stats["mean_sampled_outward_deviation"], 0.0)
        self.assertGreaterEqual(
            result.stats["maximum_sampled_outward_deviation"],
            result.stats["mean_sampled_outward_deviation"],
        )
        for region in result.regions:
            self.assertGreaterEqual(region.primitive.mean_sampled_outward_deviation, 0.0)
            self.assertGreaterEqual(
                region.primitive.maximum_sampled_outward_deviation,
                region.primitive.mean_sampled_outward_deviation,
            )
        self.assertLessEqual(
            {region.primitive.kind for region in result.regions},
            {"sphere", "rss"},
        )
        timings = result.stats["timings_seconds"]
        for phase in (
            "degenerate_filter",
            "face_adjacency",
            "initial_primitive_fitting",
            "candidate_primitive_fitting",
            "merge_loop",
            "outward_deviation_metrics",
            "total",
        ):
            self.assertGreaterEqual(timings[phase], 0.0)
        self.assertGreaterEqual(result.stats["candidate_primitive_fits"], 0)

    def test_loose_concave_planar_region_is_replaced_by_triangle_primitives(self) -> None:
        source = _subdivided_l_shape()
        result = analyze_mesh(
            source,
            PrimitiveAnalysisOptions(analysis_strength=0.0, primitive_types=("rss", "triangle")),
        )
        self.assertEqual(result.stats["triangle_patches"], 1)
        self.assertGreater(result.stats["triangle_primitives"], 0)
        self.assertLessEqual(result.stats["source_triangles"], len(source.faces))
        self.assertEqual({region.primitive.kind for region in result.regions}, {"triangle"})
        self.assertAlmostEqual(result.stats["mean_sampled_outward_deviation"], 0.0)
        self.assertAlmostEqual(result.stats["maximum_sampled_outward_deviation"], 0.0)
        source_area = np.linalg.norm(
            np.cross(
                source.vertices[source.faces][:, 1] - source.vertices[source.faces][:, 0],
                source.vertices[source.faces][:, 2] - source.vertices[source.faces][:, 0],
            ),
            axis=1,
        ).sum() * 0.5
        proxy_triangles = np.asarray(
            [region.primitive.triangle_vertices for region in result.regions],
            dtype=np.float64,
        )
        result_area = np.linalg.norm(
            np.cross(
                proxy_triangles[:, 1] - proxy_triangles[:, 0],
                proxy_triangles[:, 2] - proxy_triangles[:, 0],
            ),
            axis=1,
        ).sum() * 0.5
        self.assertAlmostEqual(result_area, source_area)

    def test_planar_hole_is_filled_without_becoming_outward_error(self) -> None:
        source = _square_ring()
        result = analyze_mesh(
            source,
            PrimitiveAnalysisOptions(analysis_strength=0.0, primitive_types=("rss", "triangle")),
        )
        self.assertEqual(result.stats["filled_planar_holes"], 1)
        self.assertAlmostEqual(result.stats["filled_planar_hole_area"], 1.0)
        self.assertEqual(result.stats["triangle_primitives"], 0)
        self.assertEqual([region.primitive.kind for region in result.regions], ["rss"])
        self.assertEqual(result.stats["source_triangles"], len(source.faces))
        self.assertLess(result.stats["maximum_planar_excess_area_ratio"], 0.02)
        self.assertLess(result.stats["mean_sampled_outward_deviation"], 1.0e-3)

    def test_short_boundary_void_is_conservatively_completed(self) -> None:
        source = Polygon(
            [
                (-4012.38, 3505.75),
                (-3562.38, 3505.75),
                (-3562.38, 3705.75),
                (-3822.38, 4055.75),
                (-4022.38, 4055.75),
                (-4022.38, 3515.75),
            ]
        )
        completed, count, area = _fill_short_boundary_voids(
            source, 1.0e-8, PrimitiveAnalysisOptions()
        )
        self.assertEqual(count, 1)
        self.assertEqual(len(completed.exterior.coords) - 1, 5)
        self.assertTrue(completed.covers(source))
        self.assertAlmostEqual(area, 50.0)

    def test_projection_hole_caps_absorb_inner_wall_responsibility(self) -> None:
        result = analyze_mesh(
            _box_frame(),
            PrimitiveAnalysisOptions(max_excess_volume_ratio=0.0),
        )
        self.assertEqual(result.stats["projection_filled_holes"], 1)
        self.assertEqual(result.stats["projection_fill_primitives"], 4)
        self.assertGreater(result.stats["excluded_void_surface_triangles"], 0)
        self.assertGreater(result.stats["removed_void_surface_primitives"], 0)
        assigned = np.concatenate([region.triangle_ids for region in result.regions])
        combined = np.concatenate((assigned, result.excluded_triangle_ids))
        np.testing.assert_array_equal(np.sort(combined), np.arange(len(result.mesh.faces)))

    def test_complete_planar_rectangle_remains_an_rss(self) -> None:
        result = analyze_mesh(
            _subdivided_rectangle(),
            PrimitiveAnalysisOptions(analysis_strength=0.0, primitive_types=("rss", "triangle")),
        )
        self.assertEqual(result.stats["triangle_primitives"], 0)
        self.assertEqual([region.primitive.kind for region in result.regions], ["rss"])
        self.assertLess(result.stats["maximum_planar_excess_area_ratio"], 0.02)

    def test_higher_strength_never_reduces_primitive_count(self) -> None:
        mesh = _bent_strip()
        counts = [
            len(analyze_mesh(mesh, PrimitiveAnalysisOptions(analysis_strength=value)).regions)
            for value in (0.0, 0.5, 1.0)
        ]
        self.assertEqual(counts, sorted(counts))

    def test_directory_output_contains_viewer_assets_and_metadata(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            source = root / "source"
            output = root / "output"
            source.mkdir()
            write_obj(source / "1.obj", _tetrahedron())
            manifest = analyze_directory(
                source,
                output,
                PrimitiveAnalysisOptions(max_excess_volume_ratio=0.0),
            )
            model_dir = output / "models" / "1"
            self.assertEqual(manifest["model_count"], 1)
            for name in ("source.obj", "regions.obj", "primitives.obj", "model.json"):
                self.assertTrue((model_dir / name).is_file(), name)
            self.assertTrue((output / "viewer_manifest.json").is_file())
            metadata = json.loads((model_dir / "model.json").read_text(encoding="ascii"))
            self.assertEqual(
                sum(record["triangle_count"] for record in metadata["primitive_records"]),
                len(_tetrahedron().faces),
            )
            self.assertTrue(metadata["stats"]["strict_triangle_partition"])
            self.assertTrue(metadata["stats"]["conservative_coverage_by_construction"])
            self.assertIn("mean_sampled_outward_deviation", metadata["stats"])
            self.assertIn("maximum_sampled_outward_deviation", metadata["primitive_records"][0])
            self.assertIn("synthetic_hole_fill", metadata["primitive_records"][0])
            self.assertIn("timings_seconds", metadata["stats"])
            self.assertIn("summary", manifest)
            self.assertEqual(manifest["summary"]["active_source_triangles"], 4)
            self.assertGreaterEqual(manifest["summary"]["mean_sampled_outward_deviation"], 0.0)
            self.assertGreaterEqual(manifest["models"][0]["mesh_analysis_seconds"], 0.0)

    def test_unlimited_strength_writes_standard_json(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            source = root / "source"
            output = root / "output"
            source.mkdir()
            write_obj(source / "1.obj", _tetrahedron())
            analyze_directory(
                source,
                output,
                PrimitiveAnalysisOptions(analysis_strength=0.0),
            )
            manifest_text = (output / "viewer_manifest.json").read_text(encoding="ascii")
            model_text = (output / "models" / "1" / "model.json").read_text(encoding="ascii")
            self.assertNotIn("Infinity", manifest_text)
            self.assertNotIn("Infinity", model_text)
            self.assertIsNone(json.loads(manifest_text)["effective_max_excess_volume_ratio"])

    def test_resume_reuses_completed_models(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            source = root / "source"
            output = root / "output"
            source.mkdir()
            write_obj(source / "1.obj", _tetrahedron())
            first = analyze_directory(source, output, PrimitiveAnalysisOptions())
            resumed = analyze_directory(source, output, PrimitiveAnalysisOptions(), resume=True)
            self.assertTrue(first["complete"])
            self.assertTrue(resumed["complete"])
            self.assertEqual(resumed["model_count"], 1)
            self.assertTrue(resumed["models"][0]["resumed"])

    def test_collision_export_preserves_partition_and_primitive_parameters(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            source = root / "source"
            analysis = root / "analysis"
            collision = root / "collision"
            source.mkdir()
            write_obj(source / "1.obj", _tetrahedron())
            analyze_directory(source, analysis, PrimitiveAnalysisOptions(max_excess_volume_ratio=0.0))
            report = export_pool(analysis, collision)
            pool_lines = (collision / "primitive_bvh_pool.tsv").read_text(encoding="ascii").splitlines()
            model_lines = (collision / "models" / "1.tsv").read_text(encoding="ascii").splitlines()
            primitive_lines = [line for line in model_lines if line.startswith("primitive\t")]
            self.assertEqual(report["model_count"], 1)
            self.assertEqual(report["source_triangles"], 4)
            self.assertEqual(pool_lines[:2], ["PQSS_PRIMITIVE_BVH_POOL_V1", "model_id\tmodel_file"])
            self.assertEqual(len(primitive_lines), report["primitive_count"])
            assigned = []
            for line in primitive_lines:
                fields = line.split("\t")
                self.assertEqual(len(fields), 19)
                assigned.extend(int(value) for value in fields[-1].split(","))
            self.assertEqual(sorted(assigned), list(range(4)))

    def test_collision_export_rejects_triangle_primitive_analysis(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            source = root / "source"
            analysis = root / "analysis"
            collision = root / "collision"
            source.mkdir()
            write_obj(source / "1.obj", _subdivided_l_shape())
            analyze_directory(
                source,
                analysis,
                PrimitiveAnalysisOptions(
                    analysis_strength=0.0,
                    primitive_types=("rss", "triangle"),
                ),
            )
            with self.assertRaisesRegex(ValueError, "triangle primitive distance semantics"):
                export_pool(analysis, collision)


if __name__ == "__main__":
    unittest.main()
