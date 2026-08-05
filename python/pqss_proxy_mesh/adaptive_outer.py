from __future__ import annotations

from dataclasses import asdict, dataclass, replace
import json
import math
from pathlib import Path
from typing import Any

import numpy as np
import numpy.typing as npt
from scipy.spatial import ConvexHull, QhullError

from .obj_io import FloatArray, IndexArray, Mesh, read_obj, write_obj


@dataclass(slots=True)
class GenerationOptions:
    max_pqss_bvh_depth: int = 8
    relative_offset: float = 0.005
    depth_safety_margin: int = 2
    max_parts_per_model: int = 16
    max_hull_faces_per_part: int = 96
    min_triangles_per_part: int = 8
    min_split_gain: float = 0.02
    exact_passthrough_face_limit: int = 12

    def validate(self) -> None:
        if self.max_pqss_bvh_depth < 4:
            raise ValueError("max_pqss_bvh_depth must be at least 4")
        if not math.isfinite(self.relative_offset) or self.relative_offset <= 0.0:
            raise ValueError("relative_offset must be finite and greater than zero")
        if self.depth_safety_margin < 0:
            raise ValueError("depth_safety_margin must not be negative")
        if self.max_parts_per_model < 1:
            raise ValueError("max_parts_per_model must be positive")
        if self.max_hull_faces_per_part < 4:
            raise ValueError("max_hull_faces_per_part must be at least 4")
        if self.min_triangles_per_part < 1:
            raise ValueError("min_triangles_per_part must be positive")
        if not 0.0 <= self.min_split_gain < 1.0:
            raise ValueError("min_split_gain must be in [0, 1)")
        if self.exact_passthrough_face_limit < 0:
            raise ValueError("exact_passthrough_face_limit must not be negative")


@dataclass(slots=True)
class Cluster:
    triangle_ids: IndexArray


@dataclass(slots=True)
class Split:
    left: Cluster
    right: Cluster
    gain: float
    saving: float


@dataclass(slots=True)
class ProxyPart:
    method: str
    vertices: FloatArray
    faces: IndexArray
    volume: float
    minimum_clearance: float
    source_triangle_count: int


def _surface_area(bounds_min: FloatArray, bounds_max: FloatArray) -> float:
    extent = np.maximum(bounds_max - bounds_min, 0.0)
    return float(2.0 * (extent[0] * extent[1] + extent[1] * extent[2] + extent[2] * extent[0]))


def _best_sah_split(mesh: Mesh, triangle_ids: IndexArray, min_count: int) -> Split | None:
    if len(triangle_ids) < 2 * min_count:
        return None

    triangles = mesh.vertices[mesh.faces[triangle_ids]]
    tri_min = triangles.min(axis=1)
    tri_max = triangles.max(axis=1)
    centroids = triangles.mean(axis=1)
    parent_min = tri_min.min(axis=0)
    parent_max = tri_max.max(axis=0)
    parent_cost = _surface_area(parent_min, parent_max) * len(triangle_ids)
    if parent_cost <= 0.0:
        return None

    best: tuple[float, npt.NDArray[np.bool_]] | None = None
    bucket_count = 16
    for axis in range(3):
        low = float(centroids[:, axis].min())
        high = float(centroids[:, axis].max())
        if high <= low:
            continue
        buckets = np.minimum(
            ((centroids[:, axis] - low) * bucket_count / (high - low)).astype(np.int64),
            bucket_count - 1,
        )
        for cut in range(bucket_count - 1):
            left_mask = buckets <= cut
            left_count = int(left_mask.sum())
            right_count = len(triangle_ids) - left_count
            if left_count < min_count or right_count < min_count:
                continue
            left_cost = _surface_area(tri_min[left_mask].min(axis=0), tri_max[left_mask].max(axis=0)) * left_count
            right_mask = ~left_mask
            right_cost = (
                _surface_area(tri_min[right_mask].min(axis=0), tri_max[right_mask].max(axis=0)) * right_count
            )
            cost = left_cost + right_cost
            if best is None or cost < best[0]:
                best = (cost, left_mask)

    if best is None:
        return None
    gain = 1.0 - best[0] / parent_cost
    left_ids = triangle_ids[best[1]]
    right_ids = triangle_ids[~best[1]]
    return Split(Cluster(left_ids), Cluster(right_ids), gain, parent_cost - best[0])


def _partition(mesh: Mesh, options: GenerationOptions, face_budget: int) -> list[Cluster]:
    max_parts = min(options.max_parts_per_model, max(1, face_budget // 12))
    clusters = [Cluster(np.arange(len(mesh.faces), dtype=np.int64))]

    while len(clusters) < max_parts:
        choices: list[tuple[int, Split]] = []
        for index, cluster in enumerate(clusters):
            split = _best_sah_split(mesh, cluster.triangle_ids, options.min_triangles_per_part)
            if split is not None and split.gain >= options.min_split_gain:
                choices.append((index, split))
        if not choices:
            break
        cluster_index, selected = max(choices, key=lambda item: item[1].saving)
        clusters[cluster_index : cluster_index + 1] = [selected.left, selected.right]
    return clusters


_BOX_FACES = np.asarray(
    [
        (0, 2, 1), (0, 3, 2), (4, 5, 6), (4, 6, 7),
        (0, 1, 5), (0, 5, 4), (3, 7, 6), (3, 6, 2),
        (0, 4, 7), (0, 7, 3), (1, 2, 6), (1, 6, 5),
    ],
    dtype=np.int64,
)


def _cluster_points(mesh: Mesh, cluster: Cluster) -> FloatArray:
    return mesh.vertices[mesh.faces[cluster.triangle_ids].reshape(-1)]


def _expanded_obb(mesh: Mesh, cluster: Cluster, margin: float) -> ProxyPart:
    points = _cluster_points(mesh, cluster)
    sample_step = max(1, len(points) // 100_000)
    sample = points[::sample_step]
    centered = sample - sample.mean(axis=0)
    covariance = centered.T @ centered / max(len(centered), 1)
    _, axes = np.linalg.eigh(covariance)
    axes = axes[:, ::-1]
    if np.linalg.det(axes) < 0.0:
        axes[:, 2] *= -1.0

    local = points @ axes
    bounds_min = local.min(axis=0) - margin
    bounds_max = local.max(axis=0) + margin
    corners = np.asarray(
        [
            (bounds_min[0], bounds_min[1], bounds_min[2]),
            (bounds_max[0], bounds_min[1], bounds_min[2]),
            (bounds_max[0], bounds_max[1], bounds_min[2]),
            (bounds_min[0], bounds_max[1], bounds_min[2]),
            (bounds_min[0], bounds_min[1], bounds_max[2]),
            (bounds_max[0], bounds_min[1], bounds_max[2]),
            (bounds_max[0], bounds_max[1], bounds_max[2]),
            (bounds_min[0], bounds_max[1], bounds_max[2]),
        ],
        dtype=np.float64,
    )
    vertices = corners @ axes.T
    lower_clearance = local - bounds_min
    upper_clearance = bounds_max - local
    clearance = float(min(lower_clearance.min(), upper_clearance.min()))
    volume = float(np.prod(bounds_max - bounds_min))
    return ProxyPart("obb", vertices, _BOX_FACES.copy(), volume, clearance, len(cluster.triangle_ids))


def _polygon_from_support_planes(normals: FloatArray, supports: FloatArray) -> FloatArray | None:
    vertices: list[FloatArray] = []
    tolerance = max(float(np.max(np.abs(supports))) * 1.0e-10, 1.0e-10)
    for first in range(len(normals)):
        for second in range(first + 1, len(normals)):
            matrix = np.vstack((normals[first], normals[second]))
            determinant = float(np.linalg.det(matrix))
            if abs(determinant) < 1.0e-12:
                continue
            point = np.linalg.solve(matrix, np.asarray((supports[first], supports[second])))
            if np.all(normals @ point <= supports + tolerance):
                vertices.append(point)
    if len(vertices) < 3:
        return None
    unique = np.unique(np.asarray(vertices), axis=0)
    if len(unique) < 3:
        return None
    hull = ConvexHull(unique)
    return unique[hull.vertices]


def _expanded_prism(mesh: Mesh, cluster: Cluster, margin: float, side_count: int) -> ProxyPart | None:
    points = _cluster_points(mesh, cluster)
    sample_step = max(1, len(points) // 100_000)
    sample = points[::sample_step]
    centered = sample - sample.mean(axis=0)
    covariance = centered.T @ centered / max(len(centered), 1)
    _, eigenvectors = np.linalg.eigh(covariance)

    best: ProxyPart | None = None
    for axial_index in range(3):
        cross_indices = [index for index in range(3) if index != axial_index]
        axes = np.column_stack(
            (eigenvectors[:, axial_index], eigenvectors[:, cross_indices[0]], eigenvectors[:, cross_indices[1]])
        )
        if np.linalg.det(axes) < 0.0:
            axes[:, 2] *= -1.0
        local = points @ axes
        axial_min = float(local[:, 0].min() - margin)
        axial_max = float(local[:, 0].max() + margin)

        for rotation_index in range(12):
            rotation = rotation_index * (2.0 * math.pi / side_count) / 12.0
            angles = rotation + np.arange(side_count) * (2.0 * math.pi / side_count)
            normals = np.column_stack((np.cos(angles), np.sin(angles)))
            supports = np.max(local[:, 1:] @ normals.T, axis=0) + margin
            polygon = _polygon_from_support_planes(normals, supports)
            if polygon is None:
                continue
            cross_hull = ConvexHull(polygon)
            polygon = polygon[cross_hull.vertices]
            count = len(polygon)
            local_vertices = np.vstack(
                (
                    np.column_stack((np.full(count, axial_min), polygon)),
                    np.column_stack((np.full(count, axial_max), polygon)),
                )
            )
            faces: list[tuple[int, int, int]] = []
            for index in range(count):
                following = (index + 1) % count
                faces.append((index, following, count + following))
                faces.append((index, count + following, count + index))
            for index in range(1, count - 1):
                faces.append((0, index + 1, index))
                faces.append((count, count + index, count + index + 1))

            axial_clearance = min(float((local[:, 0] - axial_min).min()), float((axial_max - local[:, 0]).min()))
            radial_clearance = float(np.min(supports - local[:, 1:] @ normals.T))
            candidate = ProxyPart(
                f"prism_{side_count}",
                local_vertices @ axes.T,
                np.asarray(faces, dtype=np.int64),
                float(cross_hull.volume * (axial_max - axial_min)),
                min(axial_clearance, radial_clearance),
                len(cluster.triangle_ids),
            )
            if best is None or candidate.volume < best.volume:
                best = candidate
    return best


def _expanded_hull(
    mesh: Mesh,
    cluster: Cluster,
    relative_offset: float,
    model_diagonal: float,
) -> ProxyPart | None:
    points = np.unique(_cluster_points(mesh, cluster), axis=0)
    if len(points) < 4:
        return None
    centered = points - points.mean(axis=0)
    if np.linalg.matrix_rank(centered, tol=max(model_diagonal, 1.0) * 1.0e-12) < 3:
        return None

    try:
        source_hull = ConvexHull(points, qhull_options="Qt Qx")
    except QhullError:
        return None

    hull_points = points[source_hull.vertices]
    center = hull_points.mean(axis=0)
    required_clearance = max(model_diagonal * 1.0e-10, 1.0e-10)
    expanded_hull: ConvexHull | None = None
    proxy_vertices = hull_points
    clearance = -math.inf

    for attempt in range(10):
        scale = 1.0 + relative_offset * (2**attempt)
        proxy_vertices = center + (hull_points - center) * scale
        try:
            expanded_hull = ConvexHull(proxy_vertices, qhull_options="Qt Qx")
        except QhullError:
            continue
        equations = expanded_hull.equations
        maximum_signed_distance = float(np.max(points @ equations[:, :3].T + equations[:, 3]))
        clearance = -maximum_signed_distance
        if clearance >= required_clearance:
            break
    else:
        return None

    assert expanded_hull is not None
    return ProxyPart(
        "local_hull",
        proxy_vertices,
        np.asarray(expanded_hull.simplices, dtype=np.int64),
        float(expanded_hull.volume),
        clearance,
        len(cluster.triangle_ids),
    )


def _choose_parts(mesh: Mesh, clusters: list[Cluster], options: GenerationOptions, face_budget: int) -> list[ProxyPart]:
    extent = mesh.vertices.max(axis=0) - mesh.vertices.min(axis=0)
    diagonal = float(np.linalg.norm(extent))
    margin = max(diagonal * options.relative_offset, 1.0e-9)

    all_candidates: list[list[ProxyPart]] = []
    for cluster in clusters:
        candidates = [_expanded_obb(mesh, cluster, margin)]
        for side_count in (4, 6, 8):
            prism = _expanded_prism(mesh, cluster, margin, side_count)
            if prism is not None:
                candidates.append(prism)
        hull = _expanded_hull(mesh, cluster, options.relative_offset, diagonal)
        if hull is not None and len(hull.faces) <= options.max_hull_faces_per_part:
            candidates.append(hull)

        best_by_face_count: dict[int, ProxyPart] = {}
        for candidate in candidates:
            face_count = len(candidate.faces)
            previous = best_by_face_count.get(face_count)
            if previous is None or candidate.volume < previous.volume:
                best_by_face_count[face_count] = candidate
        all_candidates.append(list(best_by_face_count.values()))

    states: dict[int, tuple[float, list[ProxyPart]]] = {0: (0.0, [])}
    for candidates in all_candidates:
        next_states: dict[int, tuple[float, list[ProxyPart]]] = {}
        for used_faces, (volume, selected) in states.items():
            for candidate in candidates:
                new_faces = used_faces + len(candidate.faces)
                if new_faces > face_budget:
                    continue
                new_volume = volume + candidate.volume
                previous = next_states.get(new_faces)
                if previous is None or new_volume < previous[0]:
                    next_states[new_faces] = (new_volume, selected + [candidate])
        states = next_states
    if not states:
        raise RuntimeError(f"{mesh.name}: no proxy candidate fits the face budget")
    return min(states.values(), key=lambda state: state[0])[1]


def _merge_parts(name: str, parts: list[ProxyPart]) -> Mesh:
    vertices: list[FloatArray] = []
    faces: list[IndexArray] = []
    vertex_offset = 0
    for part in parts:
        vertices.append(part.vertices)
        faces.append(part.faces + vertex_offset)
        vertex_offset += len(part.vertices)
    return Mesh(name, np.vstack(vertices), np.vstack(faces))


def generate_model(mesh: Mesh, options: GenerationOptions) -> tuple[Mesh, dict[str, Any]]:
    if len(mesh.faces) <= options.exact_passthrough_face_limit:
        stats: dict[str, Any] = {
            "model": mesh.name,
            "source_vertices": len(mesh.vertices),
            "source_faces": len(mesh.faces),
            "proxy_vertices": len(mesh.vertices),
            "proxy_faces": len(mesh.faces),
            "parts": 1,
            "methods": {"exact_passthrough": 1},
            "minimum_certified_clearance": 0.0,
            "containment_certified": True,
            "strict_clearance_certified": False,
        }
        return Mesh(mesh.name, mesh.vertices.copy(), mesh.faces.copy()), stats

    face_budget = max(12, 1 << max(0, options.max_pqss_bvh_depth - options.depth_safety_margin))
    maximum_parts = min(options.max_parts_per_model, max(1, face_budget // 12))
    best_parts: list[ProxyPart] | None = None
    best_volume = math.inf
    for part_count in range(1, maximum_parts + 1):
        clusters = _partition(mesh, replace(options, max_parts_per_model=part_count), face_budget)
        if len(clusters) != part_count:
            break
        parts = _choose_parts(mesh, clusters, options, face_budget)
        volume = sum(part.volume for part in parts)
        if volume < best_volume:
            best_volume = volume
            best_parts = parts
    if best_parts is None:
        raise RuntimeError(f"{mesh.name}: no valid adaptive proxy candidate")
    parts = best_parts
    proxy = _merge_parts(mesh.name, parts)
    minimum_clearance = min(part.minimum_clearance for part in parts)
    if not math.isfinite(minimum_clearance) or minimum_clearance <= 0.0:
        raise RuntimeError(f"{mesh.name}: strict containment certification failed")

    methods: dict[str, int] = {}
    for part in parts:
        methods[part.method] = methods.get(part.method, 0) + 1
    stats: dict[str, Any] = {
        "model": mesh.name,
        "source_vertices": len(mesh.vertices),
        "source_faces": len(mesh.faces),
        "proxy_vertices": len(proxy.vertices),
        "proxy_faces": len(proxy.faces),
        "parts": len(parts),
        "methods": methods,
        "minimum_certified_clearance": minimum_clearance,
        "containment_certified": True,
        "strict_clearance_certified": True,
    }
    return proxy, stats


def _pqss_subdivision_stats(meshes: list[Mesh]) -> tuple[float, list[dict[str, int]]]:
    metrics: list[npt.NDArray[np.float64]] = []
    for mesh in meshes:
        triangles = mesh.vertices[mesh.faces]
        cross = np.cross(triangles[:, 1] - triangles[:, 0], triangles[:, 2] - triangles[:, 0])
        metrics.append(np.einsum("ij,ij->i", cross, cross))
    threshold = float(np.concatenate(metrics).mean())

    result: list[dict[str, int]] = []
    for values in metrics:
        built_faces = 0
        subdivisions = 0
        for value in values:
            level = 0
            reduced = float(value)
            while reduced > threshold:
                level += 1
                reduced /= 16.0
            leaves = 4**level
            built_faces += leaves
            subdivisions += (leaves - 1) // 3
        min_depth = math.ceil(math.log2(built_faces)) if built_faces > 1 else 0
        result.append(
            {
                "pqss_predicted_built_faces": built_faces,
                "pqss_predicted_subdivision_count": subdivisions,
                "binary_depth_lower_bound": min_depth,
            }
        )
    return threshold, result


def _numeric_name(path: Path) -> int:
    try:
        return int(path.stem)
    except ValueError as error:
        raise ValueError(f"model filename must be numeric: {path.name}") from error


def generate_pool(input_dir: Path, output_dir: Path, options: GenerationOptions) -> dict[str, Any]:
    options.validate()
    input_paths = sorted(input_dir.glob("*.obj"), key=_numeric_name)
    if not input_paths:
        raise ValueError(f"no OBJ files in {input_dir}")
    if output_dir.exists() and any(output_dir.iterdir()):
        raise FileExistsError(f"output directory is not empty: {output_dir}")
    output_dir.mkdir(parents=True, exist_ok=True)

    proxies: list[Mesh] = []
    model_stats: list[dict[str, Any]] = []
    for input_path in input_paths:
        source = read_obj(input_path)
        proxy, stats = generate_model(source, options)
        write_obj(
            output_dir / input_path.name,
            proxy,
            [
                f"source={input_path}",
                "mode=ConservativeOuterAdaptive",
                f"parts={stats['parts']}",
                f"minimum_certified_clearance={stats['minimum_certified_clearance']:.17g}",
            ],
        )
        proxies.append(proxy)
        model_stats.append(stats)

    area_threshold, subdivision_stats = _pqss_subdivision_stats(proxies)
    for stats, subdivision in zip(model_stats, subdivision_stats, strict=True):
        stats.update(subdivision)

    report: dict[str, Any] = {
        "algorithm": "ConservativeOuterAdaptive",
        "input_directory": str(input_dir.resolve()),
        "output_directory": str(output_dir.resolve()),
        "options": asdict(options),
        "model_count": len(proxies),
        "total_source_faces": sum(item["source_faces"] for item in model_stats),
        "total_proxy_faces": sum(item["proxy_faces"] for item in model_stats),
        "pqss_pool_area_threshold": area_threshold,
        "all_models_containment_certified": all(item["containment_certified"] for item in model_stats),
        "all_generated_proxies_strictly_contained": all(
            item["strict_clearance_certified"] for item in model_stats
            if "exact_passthrough" not in item["methods"]
        ),
        "models": model_stats,
    }
    (output_dir / "generation_report.json").write_text(
        json.dumps(report, indent=2, ensure_ascii=True), encoding="ascii"
    )
    return report
