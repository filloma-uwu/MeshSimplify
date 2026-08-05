from __future__ import annotations

import argparse
from pathlib import Path

from .primitive_decomposition import (
    PrimitiveDecompositionOptions,
    generate_independent_pool,
    generate_pool,
)


def main() -> int:
    parser = argparse.ArgumentParser(description="Generate PQSS-targeted conservative primitive decomposition")
    parser.add_argument("--input-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument(
        "--analysis-strength",
        type=float,
        default=0.5,
        help="primitive-boundary sensitivity in [0,1]; higher values retain more primitives",
    )
    parser.add_argument(
        "--max-excess-volume-ratio",
        type=float,
        help="override the paper's maximum added volume relative to model AABB volume",
    )
    parser.add_argument(
        "--max-parts",
        type=int,
        default=4096,
        help="resource safety limit, not a target primitive count",
    )
    parser.add_argument(
        "--independent-parts",
        action="store_true",
        help="write every primitive as an independent PQSS model plus a manifest",
    )
    parser.add_argument("--relative-clearance", type=float, default=1.0e-6)
    args = parser.parse_args()
    options = PrimitiveDecompositionOptions(
        analysis_strength=args.analysis_strength,
        max_excess_volume_ratio=args.max_excess_volume_ratio,
        max_parts_per_model=args.max_parts,
        relative_clearance=args.relative_clearance,
    )
    report = (
        generate_independent_pool(args.input_dir, args.output_dir, options)
        if args.independent_parts
        else generate_pool(args.input_dir, args.output_dir, options)
    )
    print(f"models={report['model_count']}")
    print(f"selected_models={report['selected_model_count']}")
    if args.independent_parts:
        print(f"physical_models={report['physical_model_count']}")
    print(f"generation_seconds={report['generation_seconds']:.6f}")
    print(f"report={args.output_dir / 'generation_report.json'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
