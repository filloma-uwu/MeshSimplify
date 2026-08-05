from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

import numpy as np

from pqss_proxy_mesh.adaptive_outer import GenerationOptions, generate_model, generate_pool
from pqss_proxy_mesh.obj_io import Mesh, read_obj, write_obj
from pqss_proxy_mesh.primitive_decomposition import (
    PrimitiveDecompositionOptions,
    generate_independent_pool,
    generate_model as generate_primitive_model,
    generate_model_parts,
)
from pqss_proxy_mesh.static_bvh import measure_static_bvh


class AdaptiveOuterTest(unittest.TestCase):
    def test_expanded_tetrahedron_is_certified(self) -> None:
        mesh = Mesh(
            "tetra.obj",
            np.asarray(
                [(0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0)],
                dtype=np.float64,
            ),
            np.asarray([(0, 2, 1), (0, 1, 3), (0, 3, 2), (1, 2, 3)], dtype=np.int64),
        )
        proxy, stats = generate_model(mesh, GenerationOptions(exact_passthrough_face_limit=0))
        self.assertTrue(stats["containment_certified"])
        self.assertGreater(stats["minimum_certified_clearance"], 0.0)
        self.assertEqual(stats["methods"], {"local_hull": 1})
        self.assertEqual(len(proxy.faces), 4)

    def test_coplanar_triangle_uses_strict_outer_polytope(self) -> None:
        mesh = Mesh(
            "plane.obj",
            np.asarray([(0.0, 0.0, 0.0), (2.0, 0.0, 0.0), (0.0, 2.0, 0.0)], dtype=np.float64),
            np.asarray([(0, 1, 2)], dtype=np.int64),
        )
        proxy, stats = generate_model(mesh, GenerationOptions(exact_passthrough_face_limit=0))
        self.assertEqual(sum(stats["methods"].values()), 1)
        self.assertTrue(set(stats["methods"]).issubset({"obb", "prism_4", "prism_6", "prism_8"}))
        self.assertGreater(np.ptp(proxy.vertices[:, 2]), 0.0)
        self.assertGreater(stats["minimum_certified_clearance"], 0.0)

    def test_obj_pool_round_trip_and_report(self) -> None:
        mesh = Mesh(
            "1.obj",
            np.asarray(
                [(0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0)],
                dtype=np.float64,
            ),
            np.asarray([(0, 2, 1), (0, 1, 3), (0, 3, 2), (1, 2, 3)], dtype=np.int64),
        )
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            source = root / "source"
            output = root / "output"
            source.mkdir()
            write_obj(source / mesh.name, mesh)
            report = generate_pool(source, output, GenerationOptions())
            loaded = read_obj(output / mesh.name)
            self.assertEqual(report["model_count"], 1)
            self.assertTrue(report["all_models_containment_certified"])
            self.assertEqual(len(loaded.faces), 4)
            self.assertTrue((output / "generation_report.json").is_file())

    def test_primitive_decomposition_keeps_separated_parts(self) -> None:
        first = np.asarray(
            [(0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0)],
            dtype=np.float64,
        )
        vertices = np.vstack((first, first + np.asarray((10.0, 0.0, 0.0))))
        tetra_faces = np.asarray([(0, 2, 1), (0, 1, 3), (0, 3, 2), (1, 2, 3)], dtype=np.int64)
        faces = np.vstack((tetra_faces, tetra_faces + 4))
        mesh = Mesh("two_tetra.obj", vertices, faces)
        proxy, stats = generate_primitive_model(
            mesh,
            PrimitiveDecompositionOptions(
                max_excess_volume_ratio=1.0,
            ),
        )
        self.assertTrue(stats["containment_certified"])
        self.assertEqual(stats["proxy_parts"], 2)
        self.assertGreater(stats["minimum_certified_clearance"], 0.0)
        self.assertGreater(len(proxy.faces), 0)
        bvh = measure_static_bvh(proxy)
        self.assertEqual(bvh.nodes, 2 * len(proxy.faces) - 1)
        self.assertEqual(bvh.internal_nodes, len(proxy.faces) - 1)
        self.assertGreater(bvh.sah_sum, 0.0)

    def test_independent_pool_writes_one_obj_per_primitive(self) -> None:
        first = np.asarray(
            [(0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0)],
            dtype=np.float64,
        )
        vertices = np.vstack((first, first + np.asarray((10.0, 0.0, 0.0))))
        tetra_faces = np.asarray([(0, 2, 1), (0, 1, 3), (0, 3, 2), (1, 2, 3)], dtype=np.int64)
        mesh = Mesh("1.obj", vertices, np.vstack((tetra_faces, tetra_faces + 4)))
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            source = root / "source"
            output = root / "output"
            source.mkdir()
            write_obj(source / mesh.name, mesh)
            report = generate_independent_pool(
                source,
                output,
                PrimitiveDecompositionOptions(
                    max_excess_volume_ratio=1.0,
                ),
                model_ids=(1,),
            )
            manifest = (output / "primitive_pool_manifest.tsv").read_text(encoding="ascii").splitlines()
            self.assertEqual(report["physical_model_count"], 2)
            self.assertEqual(len(manifest), 3)
            self.assertTrue((output / "models" / "1" / "part_0000.obj").is_file())
            self.assertTrue((output / "models" / "1" / "part_0001.obj").is_file())

    def test_analysis_strength_does_not_force_a_part_count(self) -> None:
        mesh = Mesh(
            "tetra.obj",
            np.asarray(
                [(0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0)],
                dtype=np.float64,
            ),
            np.asarray([(0, 2, 1), (0, 1, 3), (0, 3, 2), (1, 2, 3)], dtype=np.int64),
        )
        loose_parts, loose_stats = generate_model_parts(
            mesh,
            PrimitiveDecompositionOptions(analysis_strength=0.0, max_excess_volume_ratio=2.0),
        )
        strict_parts, strict_stats = generate_model_parts(
            mesh,
            PrimitiveDecompositionOptions(analysis_strength=1.0, max_excess_volume_ratio=0.0),
        )
        self.assertEqual(len(loose_parts), 1)
        self.assertGreaterEqual(len(strict_parts), len(loose_parts))
        self.assertEqual(loose_stats["proxy_parts"], 1)
        self.assertEqual(loose_stats["merge_analysis"]["initial_face_primitives"], 4)


if __name__ == "__main__":
    unittest.main()
