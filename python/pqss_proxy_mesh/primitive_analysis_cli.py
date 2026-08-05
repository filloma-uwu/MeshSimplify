from __future__ import annotations

import argparse
from pathlib import Path

from .primitive_analysis import PRIMITIVE_TYPES, PrimitiveAnalysisOptions, analyze_directory


def _model_ids(value: str) -> tuple[int, ...]:
    result = tuple(int(field.strip()) for field in value.split(",") if field.strip())
    if not result or any(model_id < 1 for model_id in result):
        raise argparse.ArgumentTypeError("models must be a comma-separated list of positive integers")
    return result


def _primitive_types(value: str) -> tuple[str, ...]:
    result = tuple(field.strip().lower() for field in value.split(",") if field.strip())
    unknown = sorted(set(result) - set(PRIMITIVE_TYPES))
    if not result or unknown:
        raise argparse.ArgumentTypeError(
            "primitive types must be a comma-separated subset of sphere,capsule,rss,triangle"
        )
    return result


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Analyze triangle meshes into a strict partition of conservative primitives"
    )
    parser.add_argument("--input-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--models", type=_model_ids)
    parser.add_argument("--primitive-types", type=_primitive_types, default=PRIMITIVE_TYPES)
    parser.add_argument("--analysis-strength", type=float, default=0.5)
    parser.add_argument("--max-excess-volume-ratio", type=float)
    parser.add_argument("--relative-clearance", type=float, default=1.0e-6)
    parser.add_argument("--capsule-radius-samples", type=int, default=18)
    parser.add_argument("--planar-triangle-max-coverage", type=float, default=0.98)
    parser.add_argument(
        "--model-workers",
        type=int,
        default=0,
        help="models analyzed in parallel; zero uses all available logical CPU cores",
    )
    parser.add_argument(
        "--resume",
        action="store_true",
        help="reuse completed model directories and continue an interrupted analysis",
    )
    args = parser.parse_args()

    options = PrimitiveAnalysisOptions(
        analysis_strength=args.analysis_strength,
        max_excess_volume_ratio=args.max_excess_volume_ratio,
        primitive_types=args.primitive_types,
        relative_clearance=args.relative_clearance,
        capsule_radius_samples=args.capsule_radius_samples,
        model_workers=args.model_workers,
        planar_triangle_max_coverage=args.planar_triangle_max_coverage,
    )
    report = analyze_directory(
        args.input_dir, args.output_dir, options, model_ids=args.models, resume=args.resume
    )
    summary = report["summary"]
    print(f"models={report['model_count']}")
    print(f"original_source_triangles={summary['original_source_triangles']}")
    print(f"active_source_triangles={summary['active_source_triangles']}")
    print(f"discarded_degenerate_triangles={summary['discarded_degenerate_triangles']}")
    print(f"primitives={summary['primitive_count']}")
    print(f"mean_sampled_outward_deviation={summary['mean_sampled_outward_deviation']:.9g}")
    print(f"maximum_sampled_outward_deviation={summary['maximum_sampled_outward_deviation']:.9g}")
    print(f"maximum_planar_excess_area_ratio={summary['maximum_planar_excess_area_ratio']:.9g}")
    print(f"filled_planar_holes={summary['filled_planar_holes']}")
    print(f"filled_planar_hole_area={summary['filled_planar_hole_area']:.9g}")
    print(f"synthetic_fill_primitives={summary['synthetic_fill_primitives']}")
    print(f"filled_boundary_voids={summary['filled_boundary_voids']}")
    print(f"projection_filled_holes={summary['projection_filled_holes']}")
    print(f"excluded_void_surface_triangles={summary['excluded_void_surface_triangles']}")
    print(f"excluded_redundant_triangles={summary['excluded_redundant_triangles']}")
    print(f"removed_void_surface_primitives={summary['removed_void_surface_primitives']}")
    print(f"model_analysis_seconds={summary['model_analysis_seconds']:.6f}")
    print(f"analysis_seconds={report['analysis_seconds']:.6f}")
    print(f"manifest={args.output_dir / 'viewer_manifest.json'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
