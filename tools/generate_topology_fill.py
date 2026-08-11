"""Generate enclosing topology-fill references for arbitrary OBJ soups."""

from __future__ import annotations

import argparse
from pathlib import Path

from pqss_proxy_mesh.topology_fill import generate_batch


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("inputs", nargs="+", type=Path)
    parser.add_argument("--output-root", required=True, type=Path)
    parser.add_argument("--maximum-grid-voxels", type=int, default=3_000_000)
    parser.add_argument("--padding", type=int, default=4)
    args = parser.parse_args()
    manifest = generate_batch(
        args.inputs,
        args.output_root,
        args.maximum_grid_voxels,
        args.padding,
    )
    print(f"manifest={manifest}")


if __name__ == "__main__":
    main()
