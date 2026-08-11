from __future__ import annotations

import os
import tempfile
import unittest
from pathlib import Path

import numpy as np

from pqss_proxy_mesh.topology_fill import (
    BettiNumbers,
    enclosing_topology_fill,
    marching_cubes_enclosure_guard,
    parse_solver_log,
    rasterize_triangle_surface,
    run_fill_only_solver,
    signed_distance_image,
    voxel_betti_numbers,
)


def torus_volume(size: int = 48) -> np.ndarray:
    coordinates = np.arange(size, dtype=np.float64) - (size - 1) * 0.5
    x, y, z = np.meshgrid(coordinates, coordinates, coordinates, indexing="ij")
    major_radius = size * 0.24
    minor_radius = size * 0.10
    return (np.sqrt(x * x + y * y) - major_radius) ** 2 + z * z <= minor_radius ** 2


class TopologyFillTest(unittest.TestCase):
    def test_marching_cubes_guard_surrounds_every_source_cell(self) -> None:
        source = np.zeros((9, 9, 9), dtype=bool)
        source[4, 4, 4] = True
        guarded = marching_cubes_enclosure_guard(source)

        self.assertEqual(int(np.count_nonzero(guarded)), 27)
        self.assertTrue(np.all(guarded[3:6, 3:6, 3:6]))

    def test_enclosing_fill_removes_torus_handle_without_cutting_source(self) -> None:
        source = torus_volume(32)
        output, audit = enclosing_topology_fill(source)

        self.assertEqual(audit.input_betti, BettiNumbers(1, 1, 0))
        self.assertEqual(audit.output_betti, BettiNumbers(1, 0, 0))
        self.assertEqual(audit.removed_voxels, 0)
        self.assertGreater(audit.added_voxels, 0)
        self.assertTrue(np.all(output[source]))

    def test_enclosing_fill_rebuilds_candidates_after_each_axis(self) -> None:
        size = 72
        x, y, z = np.indices((size, size, size), dtype=np.float64)
        major_radius = 7.0
        minor_radius = 2.8
        x_torus = (
            (np.sqrt((y - 36) ** 2 + (z - 36) ** 2) - major_radius) ** 2
            + (x - 16) ** 2
            <= minor_radius ** 2
        )
        y_torus = (
            (np.sqrt((x - 36) ** 2 + (z - 36) ** 2) - major_radius) ** 2
            + (y - 36) ** 2
            <= minor_radius ** 2
        )
        z_torus = (
            (np.sqrt((x - 56) ** 2 + (y - 36) ** 2) - major_radius) ** 2
            + (z - 36) ** 2
            <= minor_radius ** 2
        )
        source = x_torus | y_torus | z_torus
        source[16:57, 33:40, 33:40] = True

        output, audit = enclosing_topology_fill(source)

        self.assertEqual(audit.output_betti, BettiNumbers(1, 0, 0))
        self.assertGreaterEqual(len(audit.steps), 2)
        self.assertGreaterEqual(len({step.axis for step in audit.steps}), 2)
        self.assertTrue(np.all(output[source]))

    def test_triangle_rasterizer_covers_deterministic_surface_samples(self) -> None:
        vertices = np.asarray([
            [0.15, 0.20, 0.25],
            [9.60, 1.10, 7.80],
            [1.25, 8.70, 4.40],
        ])
        grid = rasterize_triangle_surface(
            vertices, np.asarray([[0, 1, 2]]), 1.0, np.zeros(3), (12, 12, 12),
            chunk_cells=7,
        )
        for first in range(21):
            for second in range(21 - first):
                u = first / 20.0
                v = second / 20.0
                point = vertices[0] + u * (vertices[1] - vertices[0]) + v * (vertices[2] - vertices[0])
                cell = tuple(np.floor(point + 0.5).astype(np.int64))
                self.assertTrue(grid[cell], f"uncovered triangle sample at {point}")

    def test_voxel_betti_numbers_distinguish_handle_and_cavity(self) -> None:
        torus = torus_volume()
        self.assertEqual(voxel_betti_numbers(torus), BettiNumbers(1, 1, 0))

        sphere = np.indices((33, 33, 33), dtype=np.float64) - 16.0
        radius = np.sqrt(np.sum(sphere * sphere, axis=0))
        shell = (radius <= 12.0) & (radius >= 6.0)
        self.assertEqual(voxel_betti_numbers(shell), BettiNumbers(1, 0, 1))

    def test_distance_field_preserves_threshold_occupancy(self) -> None:
        occupancy = torus_volume(32)
        encoded = signed_distance_image(occupancy)
        self.assertTrue(np.array_equal(encoded > 127, occupancy))
        self.assertGreater(int(encoded.max()), 127)
        self.assertLessEqual(int(encoded.max()), 254)

    def test_solver_log_parser(self) -> None:
        parsed = parse_solver_log(
            "cuts before: 0 num fills before 5\n"
            "Original Shape Topology: Components: 1 Cavities: 0 Cycles: 1\n"
            "Final topology: Components: 1 Cavities: 0 Cycles: 0\n"
        )
        self.assertEqual(parsed, ((1, 0, 1), (1, 0, 0), 0, 5))

    @unittest.skipUnless(
        os.environ.get("PQSS_TOPO_SIMPLIFIER"),
        "set PQSS_TOPO_SIMPLIFIER to run the paper implementation regression",
    )
    def test_paper_fill_only_torus_is_a_strict_superset(self) -> None:
        executable = Path(os.environ["PQSS_TOPO_SIMPLIFIER"])
        field = signed_distance_image(torus_volume())
        with tempfile.TemporaryDirectory() as temporary:
            output, audit, _ = run_fill_only_solver(
                field, executable, Path(temporary), global_time=30, local_time=10
            )
        self.assertEqual(audit.input_betti, BettiNumbers(1, 1, 0))
        self.assertEqual(audit.output_betti, BettiNumbers(1, 0, 0))
        self.assertEqual(audit.cut_candidates, 0)
        self.assertEqual(audit.removed_voxels, 0)
        self.assertGreater(audit.added_voxels, 0)
        self.assertTrue(np.all(output[field > 127]))


if __name__ == "__main__":
    unittest.main()
