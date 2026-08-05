#!/usr/bin/env python3
"""Audit one generated staged-surface viewer manifest and its proxy OBJs."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any


def load_obj(path: Path) -> tuple[list[tuple[float, float, float]], list[tuple[int, ...]]]:
    vertices: list[tuple[float, float, float]] = []
    faces: list[tuple[int, ...]] = []
    for line_number, line in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
        fields = line.split()
        if not fields or fields[0].startswith("#"):
            continue
        if fields[0] == "v":
            if len(fields) < 4:
                raise ValueError(f"{path}:{line_number}: incomplete vertex")
            vertex = tuple(float(value) for value in fields[1:4])
            if not all(math.isfinite(value) for value in vertex):
                raise ValueError(f"{path}:{line_number}: non-finite vertex")
            vertices.append(vertex)
        elif fields[0] == "f":
            indices: list[int] = []
            for token in fields[1:]:
                raw = int(token.split("/", 1)[0])
                index = raw - 1 if raw > 0 else len(vertices) + raw
                if index < 0 or index >= len(vertices):
                    raise ValueError(f"{path}:{line_number}: invalid face index {raw}")
                indices.append(index)
            faces.append(tuple(indices))
    return vertices, faces


def proxy_geometry_audit(path: Path) -> dict[str, int]:
    vertices, faces = load_obj(path)
    if not vertices or not faces:
        raise ValueError(f"{path}: empty proxy geometry")
    lower = [min(vertex[axis] for vertex in vertices) for axis in range(3)]
    upper = [max(vertex[axis] for vertex in vertices) for axis in range(3)]
    diagonal = math.sqrt(sum((upper[axis] - lower[axis]) ** 2 for axis in range(3)))
    area_epsilon_squared = max(diagonal**4 * 1.0e-28, 1.0e-48)
    quantum = max(diagonal * 1.0e-10, 1.0e-12)
    signatures: set[tuple[tuple[int, int, int], ...]] = set()
    duplicate_count = 0
    degenerate_count = 0
    non_triangle_count = 0
    for face in faces:
        if len(face) != 3:
            non_triangle_count += 1
            continue
        first, second, third = (vertices[index] for index in face)
        ab = tuple(second[axis] - first[axis] for axis in range(3))
        ac = tuple(third[axis] - first[axis] for axis in range(3))
        cross = (
            ab[1] * ac[2] - ab[2] * ac[1],
            ab[2] * ac[0] - ab[0] * ac[2],
            ab[0] * ac[1] - ab[1] * ac[0],
        )
        if sum(value * value for value in cross) <= area_epsilon_squared:
            degenerate_count += 1
        signature = tuple(sorted(
            tuple(int(round(value / quantum)) for value in vertices[index])
            for index in face
        ))
        if signature in signatures:
            duplicate_count += 1
        else:
            signatures.add(signature)
    return {
        "vertices": len(vertices),
        "faces": len(faces),
        "non_triangle_faces": non_triangle_count,
        "degenerate_triangles": degenerate_count,
        "duplicate_triangles": duplicate_count,
    }


def read_optional_float(path: Path) -> float | None:
    if not path.is_file():
        return None
    return float(path.read_text(encoding="ascii").strip())


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("manifest", type=Path)
    parser.add_argument("--json-report", type=Path)
    parser.add_argument("--markdown-report", type=Path)
    args = parser.parse_args()

    manifest_path = args.manifest.resolve()
    manifest = json.loads(manifest_path.read_text(encoding="ascii"))
    rows: list[dict[str, Any]] = []
    failures: list[str] = []
    for model in manifest.get("models", []):
        model_id = str(model["id"])
        metadata_path = (manifest_path.parent / model["metadata"]).resolve()
        metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
        stats = metadata["stats"]
        proxy_path = metadata_path.parent / metadata["proxy"]
        geometry = proxy_geometry_audit(proxy_path)
        error = stats["simplification_error"]
        limit = float(error["maximum_distance_limit"])
        maximum = float(error["maximum_distance"])
        limit_tolerance = max(1.0e-8, abs(limit) * 1.0e-10)
        pre_repair_path = metadata_path.parent / "coverage_audit_pre_repair.json"
        pre_repair = json.loads(pre_repair_path.read_text(encoding="utf-8"))
        peak_memory = read_optional_float(metadata_path.parent / "peak_memory_mb.txt")
        peak_private = read_optional_float(
            metadata_path.parent / "peak_private_memory_mb.txt")
        model_failures: list[str] = []
        if not stats["containment_validation"]["passed"]:
            model_failures.append("containment validation failed")
        if stats["coverage_audit"]["failed_source_faces"] != 0:
            model_failures.append("coverage audit has failed source faces")
        if maximum > limit + limit_tolerance:
            model_failures.append(f"maximum error {maximum} exceeds {limit}")
        if geometry["non_triangle_faces"]:
            model_failures.append("proxy contains non-triangle faces")
        if geometry["degenerate_triangles"]:
            model_failures.append("proxy contains degenerate triangles")
        if geometry["duplicate_triangles"]:
            model_failures.append("proxy contains duplicate triangles")
        if geometry["faces"] != stats["proxy_triangles"]:
            model_failures.append("OBJ triangle count differs from model.json")
        failures.extend(f"model {model_id}: {failure}" for failure in model_failures)
        rows.append({
            "model": int(model_id),
            "source_triangles": stats["source_triangles"],
            "primitive_count": stats["primitive_count"],
            "proxy_triangles": stats["proxy_triangles"],
            "maximum_error": maximum,
            "maximum_error_limit": limit,
            "coverage_repairs": pre_repair["repair_face_count"],
            "coverage_repair_output_primitives":
                pre_repair.get("repair_output_primitives", pre_repair["repair_face_count"]),
            "coverage_repair_output_triangles":
                pre_repair.get("repair_output_triangles", pre_repair["repair_face_count"]),
            "analysis_seconds": stats["timings_seconds"]["total"],
            "peak_working_set_mb": peak_memory,
            "peak_private_mb": peak_private,
            **geometry,
            "passed": not model_failures,
        })

    report = {
        "manifest": str(manifest_path),
        "model_count": len(rows),
        "passed": not failures,
        "failures": failures,
        "models": rows,
    }
    json_report = args.json_report or manifest_path.parent / "audit_report.json"
    markdown_report = args.markdown_report or manifest_path.parent / "audit_report.md"
    json_report.write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n",
                           encoding="utf-8")
    header = (
        "| 模型 | 原三角形 | 图元 | 代理三角形 | 最大误差/上限 | 修补面→三角形 | "
        "耗时(s) | 峰值私有内存(MB) | 结果 |\n"
        "|---:|---:|---:|---:|---:|---:|---:|---:|:---:|\n"
    )
    lines = [header]
    for row in rows:
        peak_private_text = (f"{row['peak_private_mb']:.1f}"
                             if row['peak_private_mb'] is not None else "—")
        lines.append(
            f"| {row['model']} | {row['source_triangles']} | {row['primitive_count']} | "
            f"{row['proxy_triangles']} | {row['maximum_error']:.6g}/"
            f"{row['maximum_error_limit']:.6g} | {row['coverage_repairs']}→"
            f"{row['coverage_repair_output_triangles']} | "
            f"{row['analysis_seconds']:.3f} | {peak_private_text} | "
            f"{'通过' if row['passed'] else '失败'} |\n"
        )
    if failures:
        lines.append("\n## 失败项\n\n")
        lines.extend(f"- {failure}\n" for failure in failures)
    markdown_report.write_text("".join(lines), encoding="utf-8")
    print(f"models={len(rows)} passed={not failures}")
    print(f"json_report={json_report}")
    print(f"markdown_report={markdown_report}")
    for failure in failures:
        print(f"FAIL: {failure}")
    return 0 if not failures else 1


if __name__ == "__main__":
    raise SystemExit(main())
