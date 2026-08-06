"""Create a viewer-compatible exact-source fallback for a failed model."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import shutil


def obj_counts(path: Path) -> tuple[int, int]:
    vertices = faces = 0
    with path.open("r", encoding="utf-8", errors="ignore") as stream:
        for line in stream:
            if line.startswith("v "):
                vertices += 1
            elif line.startswith("f "):
                fields = line.split()[1:]
                faces += max(len(fields) - 2, 0)
    return vertices, faces


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--model-id", type=int, required=True)
    parser.add_argument("--reason", required=True)
    args = parser.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)
    _, triangles = obj_counts(args.input)
    assets = (
        "source.obj",
        "phase1_hole_filled.obj",
        "phase2_recognized_surfaces.obj",
        "primitives.obj",
        "proxy.obj",
        "regions.obj",
    )
    for name in assets:
        shutil.copy2(args.input, args.output_dir / name)

    metadata = {
        "stats": {
            "model": args.input.name,
            "source_triangles": triangles,
            "discarded_degenerate_triangles": 0,
            "primitive_count": triangles,
            "primitive_types": {
                "polygon": triangles,
                "disk": 0,
                "annulus": 0,
                "cylindricalband": 0,
                "conicalband": 0,
            },
            "triangulated_proxy_triangles": triangles,
            "proxy_triangles": triangles,
            "containment_validation": {
                "passed": True,
                "method": "exact_source_mesh_fallback",
                "source_triangles": triangles,
                "assigned_source_triangles": triangles,
                "unassigned_source_triangles": 0,
                "failed_source_triangles": 0,
            },
            "simplification_error": {
                "distance_method": "identity_surface",
                "reference": "source.obj",
                "maximum_is_sample_estimate": False,
                "maximum_distance_limit": 0.0,
                "limit_was_default": False,
                "distance_sample_count": 0,
                "mean_distance": 0.0,
                "maximum_distance": 0.0,
                "mean_distance_ratio": 0.0,
                "maximum_distance_ratio": 0.0,
                "maximum_pair": {"proxy": [0, 0, 0], "source": [0, 0, 0]},
                "seconds": 0.0,
            },
            "timings_seconds": {"total": 0.0},
            "fallback": {"used": True, "reason": args.reason},
        },
        "source": "source.obj",
        "phase1_hole_filled": "phase1_hole_filled.obj",
        "phase2_recognized_surfaces": "phase2_recognized_surfaces.obj",
        "phase3_simplified_surfaces": "primitives.obj",
        "phase4_triangulated": "proxy.obj",
        "regions": "regions.obj",
        "primitive_analysis": "primitives.obj",
        "triangulated_proxy": "proxy.obj",
        "proxy": "proxy.obj",
        "open_error_visualization": "open_error.json",
        "proxy_components": [],
    }
    (args.output_dir / "open_error.json").write_text(
        json.dumps({"boundary_points": [], "maximum_pair": None}) + "\n",
        encoding="ascii",
    )
    (args.output_dir / "model.json").write_text(
        json.dumps(metadata, separators=(",", ":")) + "\n", encoding="ascii"
    )
    print(f"fallback_model={args.model_id} triangles={triangles}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
