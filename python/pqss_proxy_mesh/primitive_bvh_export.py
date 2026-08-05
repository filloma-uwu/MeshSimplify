from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
from typing import Any


POOL_HEADER = "PQSS_PRIMITIVE_BVH_POOL_V1"
MODEL_HEADER = "PQSS_PRIMITIVE_BVH_MODEL_V1"
PRIMITIVE_TYPES = {"sphere", "capsule", "rss"}


def _number_fields(values: Any) -> list[str]:
    if isinstance(values, (int, float)):
        return [format(values, ".17g")]
    result: list[str] = []
    for value in values:
        result.extend(_number_fields(value))
    return result


def export_model(metadata_path: Path, output_path: Path) -> dict[str, Any]:
    metadata = json.loads(metadata_path.read_text(encoding="ascii"))
    records = metadata["primitive_records"]
    if not metadata["stats"].get("collision_export_supported", False):
        raise ValueError(
            f"{metadata_path}: collision export is unavailable because triangle primitive "
            "distance semantics have not been implemented"
        )
    triangle_count = int(metadata["stats"]["source_triangles"])
    assigned = [int(triangle_id) for record in records for triangle_id in record["triangle_ids"]]
    if sorted(assigned) != list(range(triangle_count)):
        raise ValueError(f"{metadata_path}: primitive triangle IDs are not a strict partition")
    source_path = (metadata_path.parent / metadata["source"]).resolve()
    relative_source = Path(os.path.relpath(source_path, output_path.parent.resolve()))
    lines = [
        MODEL_HEADER,
        f"source\t{relative_source.as_posix()}",
        f"source_triangles\t{triangle_count}",
        f"primitive_count\t{len(records)}",
    ]
    type_counts: dict[str, int] = {}
    for expected_id, record in enumerate(records):
        primitive_id = int(record["id"])
        kind = str(record["type"])
        if primitive_id != expected_id:
            raise ValueError(f"{metadata_path}: non-contiguous primitive id {primitive_id}")
        if kind not in PRIMITIVE_TYPES:
            raise ValueError(f"{metadata_path}: unsupported primitive type {kind}")
        type_counts[kind] = type_counts.get(kind, 0) + 1
        fields = ["primitive", str(primitive_id), kind]
        fields.extend(_number_fields(record["origin"]))
        fields.extend(_number_fields(record["axes"]))
        fields.extend(_number_fields(record["lengths"]))
        fields.extend(_number_fields(record["radius"]))
        fields.append(",".join(str(int(value)) for value in record["triangle_ids"]))
        lines.append("\t".join(fields))

    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text("\n".join(lines) + "\n", encoding="ascii")
    return {
        "source_triangles": triangle_count,
        "primitive_count": len(records),
        "primitive_types": type_counts,
    }


def export_pool(
    analysis_directory: Path,
    output_directory: Path,
    overwrite: bool = False,
) -> dict[str, Any]:
    manifest_path = analysis_directory / "viewer_manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="ascii"))
    if not manifest.get("complete"):
        raise ValueError(f"analysis manifest is incomplete: {manifest_path}")
    if output_directory.exists() and any(output_directory.iterdir()) and not overwrite:
        raise FileExistsError(f"output directory is not empty: {output_directory}")
    output_directory.mkdir(parents=True, exist_ok=True)

    rows = [POOL_HEADER, "model_id\tmodel_file"]
    models: list[dict[str, Any]] = []
    for item in manifest["models"]:
        model_id = int(item["id"])
        metadata_path = analysis_directory / item["metadata"]
        model_path = output_directory / "models" / f"{model_id}.tsv"
        stats = export_model(metadata_path, model_path)
        rows.append(f"{model_id}\tmodels/{model_id}.tsv")
        models.append({"id": model_id, **stats})

    pool_path = output_directory / "primitive_bvh_pool.tsv"
    pool_path.write_text("\n".join(rows) + "\n", encoding="ascii")
    report = {
        "analysis_directory": str(analysis_directory.resolve()),
        "output_directory": str(output_directory.resolve()),
        "model_count": len(models),
        "source_triangles": sum(model["source_triangles"] for model in models),
        "primitive_count": sum(model["primitive_count"] for model in models),
        "models": models,
    }
    (output_directory / "export_report.json").write_text(
        json.dumps(report, indent=2, ensure_ascii=True) + "\n", encoding="ascii"
    )
    return report


def main() -> int:
    parser = argparse.ArgumentParser(description="Export primitive analysis for the C++ Primitive BVH")
    parser.add_argument("--analysis-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()
    report = export_pool(args.analysis_dir, args.output_dir, overwrite=args.overwrite)
    print(f"models={report['model_count']}")
    print(f"source_triangles={report['source_triangles']}")
    print(f"primitives={report['primitive_count']}")
    print(f"manifest={args.output_dir / 'primitive_bvh_pool.tsv'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
