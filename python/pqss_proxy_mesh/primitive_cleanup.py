from __future__ import annotations

import argparse
from collections import Counter
import json
from pathlib import Path
from typing import Any

import numpy as np

from .obj_io import Mesh, read_obj, write_obj
from .primitive_analysis import (
    Primitive,
    PrimitiveRegion,
    _apply_outward_deviation_metrics,
    _drop_degenerate_faces,
    _write_grouped_obj,
    _write_region_obj,
    primitive_visualization_mesh,
)


def _primitive(record: dict[str, Any]) -> Primitive:
    return Primitive(
        str(record["type"]),
        np.asarray(record["origin"], dtype=np.float64),
        np.asarray(record["axes"], dtype=np.float64),
        np.asarray(record["lengths"], dtype=np.float64),
        float(record["radius"]),
        float(record["volume"]),
        (
            np.asarray(record["triangle_vertices"], dtype=np.float64)
            if record.get("triangle_vertices") is not None else None
        ),
        float(record.get("mean_sampled_outward_deviation", 0.0)),
        float(record.get("maximum_sampled_outward_deviation", 0.0)),
        (
            float(record["planar_excess_area_ratio"])
            if record.get("planar_excess_area_ratio") is not None else None
        ),
    )


def clean_analysis_model(metadata_path: Path) -> dict[str, Any]:
    metadata = json.loads(metadata_path.read_text(encoding="ascii"))
    source_path = metadata_path.parent / metadata["source"]
    source = read_obj(source_path)
    cleaned, current_source_count, removed_now = _drop_degenerate_faces(source)

    triangles = source.vertices[source.faces]
    diagonal = float(np.linalg.norm(source.vertices.max(axis=0) - source.vertices.min(axis=0)))
    double_areas = np.linalg.norm(
        np.cross(triangles[:, 1] - triangles[:, 0], triangles[:, 2] - triangles[:, 0]),
        axis=1,
    )
    threshold = max(diagonal * diagonal * 2.0e-24, 2.0e-30)
    keep = double_areas > threshold
    remap = np.full(len(source.faces), -1, dtype=np.int64)
    remap[keep] = np.arange(int(np.count_nonzero(keep)), dtype=np.int64)

    regions: list[PrimitiveRegion] = []
    records: list[dict[str, Any]] = []
    for record in metadata["primitive_records"]:
        old_ids = np.asarray(record["triangle_ids"], dtype=np.int64)
        kept_ids = old_ids[keep[old_ids]]
        if len(kept_ids) == 0:
            continue
        triangle_ids = remap[kept_ids]
        vertex_ids = np.unique(cleaned.faces[triangle_ids].reshape(-1))
        primitive = _primitive(record)
        regions.append(PrimitiveRegion(primitive, triangle_ids, vertex_ids))
        updated = dict(record)
        updated["id"] = len(records)
        updated["triangle_count"] = len(triangle_ids)
        updated["triangle_ids"] = triangle_ids.tolist()
        records.append(updated)

    if not regions:
        raise ValueError(f"{metadata_path}: no non-degenerate primitive regions remain")
    assigned = np.concatenate([region.triangle_ids for region in regions])
    if not np.array_equal(np.sort(assigned), np.arange(len(cleaned.faces), dtype=np.int64)):
        raise RuntimeError(f"{metadata_path}: cleaned triangle responsibility is not a partition")

    outward_metrics = _apply_outward_deviation_metrics(cleaned, regions)
    for record, region in zip(records, regions, strict=True):
        record.update(region.primitive.to_dict())

    previous_removed = int(metadata["stats"].get("discarded_degenerate_triangles", 0))
    original_count = int(metadata["stats"].get("original_source_triangles", current_source_count))
    type_counts = dict(Counter(record["type"] for record in records))
    metadata["stats"].update(
        {
            "source_triangles": len(cleaned.faces),
            "original_source_triangles": original_count,
            "discarded_degenerate_triangles": previous_removed + removed_now,
            "primitive_count": len(records),
            "primitive_types": type_counts,
            "total_primitive_volume": sum(float(record["volume"]) for record in records),
            "strict_triangle_partition": True,
            "conservative_coverage_by_construction": True,
            "collision_export_supported": not any(
                record["type"] == "triangle" for record in records
            ),
            **outward_metrics,
        }
    )
    metadata["primitive_records"] = records

    if removed_now:
        write_obj(
            source_path,
            cleaned,
            comments=[f"Removed {previous_removed + removed_now} degenerate source triangles"],
        )
        visual_meshes = [
            (
                f"primitive_{record['id']:05d}_{record['type']}",
                primitive_visualization_mesh(region.primitive),
            )
            for record, region in zip(records, regions, strict=True)
        ]
        _write_grouped_obj(metadata_path.parent / metadata["primitives"], visual_meshes)
        result = type(
            "CleanedResult",
            (),
            {
                "mesh": cleaned,
                "regions": regions,
                "excluded_triangle_ids": np.empty(0, dtype=np.int64),
            },
        )()
        _write_region_obj(metadata_path.parent / metadata["regions"], result)

    metadata_path.write_text(
        json.dumps(metadata, ensure_ascii=True, separators=(",", ":")) + "\n",
        encoding="ascii",
    )
    return {
        "model": metadata["stats"]["model"],
        "removed_now": removed_now,
        "discarded_degenerate_triangles": metadata["stats"]["discarded_degenerate_triangles"],
        "primitive_count": len(records),
        "stats": metadata["stats"],
    }


def clean_analysis_directory(analysis_directory: Path) -> dict[str, Any]:
    manifest_path = analysis_directory / "viewer_manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="ascii"))
    models: list[dict[str, Any]] = []
    for item in manifest["models"]:
        result = clean_analysis_model(analysis_directory / item["metadata"])
        item["stats"] = result["stats"]
        models.append({"id": int(item["id"]), **result})
        print(
            f"model={item['id']} removed={result['removed_now']} "
            f"primitives={result['primitive_count']}",
            flush=True,
        )
    manifest["models"] = manifest["models"]
    manifest["degenerate_cleanup_complete"] = True
    manifest["outward_deviation_metrics_complete"] = True
    manifest_path.write_text(
        json.dumps(manifest, indent=2, ensure_ascii=True) + "\n", encoding="ascii"
    )
    report = {
        "model_count": len(models),
        "removed_now": sum(model["removed_now"] for model in models),
        "discarded_degenerate_triangles": sum(
            model["discarded_degenerate_triangles"] for model in models
        ),
        "primitive_count": sum(model["primitive_count"] for model in models),
        "models": models,
    }
    (analysis_directory / "cleanup_report.json").write_text(
        json.dumps(report, indent=2, ensure_ascii=True) + "\n", encoding="ascii"
    )
    return report


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Remove degenerate source triangles and refresh outward-deviation metrics"
    )
    parser.add_argument("--analysis-dir", type=Path, required=True)
    args = parser.parse_args()
    report = clean_analysis_directory(args.analysis_dir)
    print(f"models={report['model_count']}")
    print(f"removed={report['removed_now']}")
    print(f"primitives={report['primitive_count']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
