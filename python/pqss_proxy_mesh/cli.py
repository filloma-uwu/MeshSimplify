from __future__ import annotations

import argparse
from pathlib import Path

from .adaptive_outer import GenerationOptions, generate_pool


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Generate a conservative PQSS proxy model pool")
    parser.add_argument("--input-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--max-depth", type=int, default=8)
    parser.add_argument("--relative-offset", type=float, default=0.005)
    parser.add_argument("--depth-safety-margin", type=int, default=2)
    parser.add_argument("--max-parts", type=int, default=16)
    parser.add_argument("--min-split-gain", type=float, default=0.02)
    return parser


def main() -> int:
    args = _parser().parse_args()
    options = GenerationOptions(
        max_pqss_bvh_depth=args.max_depth,
        relative_offset=args.relative_offset,
        depth_safety_margin=args.depth_safety_margin,
        max_parts_per_model=args.max_parts,
        min_split_gain=args.min_split_gain,
    )
    report = generate_pool(args.input_dir, args.output_dir, options)
    print(f"models={report['model_count']}")
    print(f"source_faces={report['total_source_faces']}")
    print(f"proxy_faces={report['total_proxy_faces']}")
    print(f"containment_certified={str(report['all_models_containment_certified']).lower()}")
    print(f"report={args.output_dir / 'generation_report.json'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
