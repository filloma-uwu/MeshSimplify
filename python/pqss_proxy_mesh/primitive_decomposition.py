from __future__ import annotations

from dataclasses import asdict, dataclass
import heapq
import json
import math
from pathlib import Path
import shutil
import time
from typing import Any, Iterable

import numpy as np
import trimesh

from .adaptive_outer import (
    Cluster,
    ProxyPart,
    _best_sah_split,
    _expanded_hull,
    _expanded_obb,
    _expanded_prism,
    _merge_parts,
    _surface_area,
)
from .obj_io import FloatArray, IndexArray, Mesh, read_obj, write_obj
from .static_bvh import measure_static_bvh


OFFICIAL_SIMPLIFIED_MODEL_IDS = (2, 3, 4, 5, 12, 13, 14, 15, 16, 17, 18, 19, 20)


@dataclass(slots=True)
class PrimitiveDecompositionOptions:
    analysis_strength: float = 0.5
    max_excess_volume_ratio: float | None = None
    max_parts_per_model: int = 4096
    relative_clearance: float = 1.0e-6
    max_hull_faces_per_part: int = 256
    tangent_quadric_weight: float = 1.0e-2
    union_external_boundary: bool = True

    def validate(self) -> None:
        if not math.isfinite(self.analysis_strength) or not 0.0 <= self.analysis_strength <= 1.0:
            raise ValueError("analysis_strength must be in [0, 1]")
        if self.max_excess_volume_ratio is not None and (
            not math.isfinite(self.max_excess_volume_ratio) or self.max_excess_volume_ratio < 0.0
        ):
            raise ValueError("max_excess_volume_ratio must be finite and non-negative")
        if self.max_parts_per_model < 1:
            raise ValueError("max_parts_per_model must be positive")
        if not math.isfinite(self.relative_clearance) or self.relative_clearance <= 0.0:
            raise ValueError("relative_clearance must be finite and positive")
        if self.max_hull_faces_per_part < 4:
            raise ValueError("max_hull_faces_per_part must be at least 4")
        if not math.isfinite(self.tangent_quadric_weight) or self.tangent_quadric_weight < 0.0:
            raise ValueError("tangent_quadric_weight must be finite and non-negative")

    def excess_volume_ratio(self) -> float:
        if self.max_excess_volume_ratio is not None:
            return self.max_excess_volume_ratio
        # Logarithmic sensitivity matching the paper's useful 1e-2 to 1e-6 range.
        return 10.0 ** (-2.0 - 4.0 * self.analysis_strength)


@dataclass(slots=True)
class _Node:
    node_id: int
    triangle_ids: IndexArray
    vertex_ids: IndexArray
    quadric: FloatArray
    box: ProxyPart
    neighbors: set[int]
    active: bool = True


def _numeric_name(path: Path) -> int:
    try:
        return int(path.stem)
    except ValueError as error:
        raise ValueError(f"model filename must be numeric: {path.name}") from error


def _face_quadrics(mesh: Mesh, tangent_weight: float) -> FloatArray:
    triangles = mesh.vertices[mesh.faces]
    first = triangles[:, 1] - triangles[:, 0]
    second = triangles[:, 2] - triangles[:, 0]
    cross = np.cross(first, second)
    twice_area = np.linalg.norm(cross, axis=1)
    normals = np.zeros_like(cross)
    valid = twice_area > 0.0
    normals[valid] = cross[valid] / twice_area[valid, None]

    edges = np.stack(
        (triangles[:, 1] - triangles[:, 0],
         triangles[:, 2] - triangles[:, 1],
         triangles[:, 0] - triangles[:, 2]),
        axis=1,
    )
    longest = np.argmax(np.einsum("fij,fij->fi", edges, edges), axis=1)
    tangents = edges[np.arange(len(edges)), longest]
    tangent_length = np.linalg.norm(tangents, axis=1)
    tangent_valid = tangent_length > 0.0
    tangents[tangent_valid] /= tangent_length[tangent_valid, None]
    tangents[~tangent_valid] = 0.0

    normal_terms = np.einsum("fi,fj->fij", normals, normals)
    tangent_terms = np.einsum("fi,fj->fij", tangents, tangents)
    return (0.5 * twice_area)[:, None, None] * (normal_terms + tangent_weight * tangent_terms)


def _quadric_axes(quadric: FloatArray) -> FloatArray:
    _, eigenvectors = np.linalg.eigh(quadric)
    axes = eigenvectors[:, ::-1]
    if np.linalg.det(axes) < 0.0:
        axes[:, 2] *= -1.0
    return axes


def _quadric_obb(
    mesh: Mesh,
    triangle_ids: IndexArray,
    vertex_ids: IndexArray,
    quadric: FloatArray,
    margin: float,
) -> ProxyPart:
    axes = _quadric_axes(quadric)
    points = mesh.vertices[vertex_ids]
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
    faces = np.asarray(
        [
            (0, 2, 1), (0, 3, 2), (4, 5, 6), (4, 6, 7),
            (0, 1, 5), (0, 5, 4), (3, 7, 6), (3, 6, 2),
            (0, 4, 7), (0, 7, 3), (1, 2, 6), (1, 6, 5),
        ],
        dtype=np.int64,
    )
    lower_clearance = local - bounds_min
    upper_clearance = bounds_max - local
    return ProxyPart(
        "quadric_obb",
        corners @ axes.T,
        faces,
        float(np.prod(bounds_max - bounds_min)),
        float(min(lower_clearance.min(), upper_clearance.min())),
        len(triangle_ids),
    )


def _oversegment(mesh: Mesh, triangle_limit: int) -> list[Cluster]:
    pending = [Cluster(np.arange(len(mesh.faces), dtype=np.int64))]
    result: list[Cluster] = []
    while pending:
        cluster = pending.pop()
        if len(cluster.triangle_ids) <= triangle_limit:
            result.append(cluster)
            continue
        split = _best_sah_split(mesh, cluster.triangle_ids, 1)
        if split is None:
            result.append(cluster)
            continue
        pending.append(split.right)
        pending.append(split.left)
    return result


def _world_sah_work(part: ProxyPart) -> float:
    return _surface_area(part.vertices.min(axis=0), part.vertices.max(axis=0)) * len(part.faces)


def _make_node(
    node_id: int,
    mesh: Mesh,
    triangle_ids: IndexArray,
    quadric: FloatArray,
    margin: float,
) -> _Node:
    vertex_ids = np.unique(mesh.faces[triangle_ids].reshape(-1))
    box = _quadric_obb(mesh, triangle_ids, vertex_ids, quadric, margin)
    return _Node(node_id, triangle_ids, vertex_ids, quadric, box, set())


def _merged_node(
    node_id: int,
    mesh: Mesh,
    first: _Node,
    second: _Node,
    margin: float,
) -> _Node:
    triangle_ids = np.concatenate((first.triangle_ids, second.triangle_ids))
    vertex_ids = np.union1d(first.vertex_ids, second.vertex_ids)
    quadric = first.quadric + second.quadric
    box = _quadric_obb(mesh, triangle_ids, vertex_ids, quadric, margin)
    return _Node(node_id, triangle_ids, vertex_ids, quadric, box, set())


def _merge_key(first: _Node, second: _Node, merged: _Node, model_volume: float) -> tuple[float, float]:
    excess_volume = merged.box.volume - first.box.volume - second.box.volume
    work_delta = _world_sah_work(merged.box) - _world_sah_work(first.box) - _world_sah_work(second.box)
    return excess_volume / max(model_volume, 1.0e-30), work_delta


def _bounds_overlap(first: ProxyPart, second: ProxyPart) -> tuple[float, float]:
    first_min = first.vertices.min(axis=0)
    first_max = first.vertices.max(axis=0)
    second_min = second.vertices.min(axis=0)
    second_max = second.vertices.max(axis=0)
    extent = np.maximum(np.minimum(first_max, second_max) - np.maximum(first_min, second_min), 0.0)
    area = _surface_area(np.zeros(3), extent)
    return area, float(np.prod(extent))


def _face_adjacency(mesh: Mesh) -> list[set[int]]:
    """Build the paper's face-collapse graph after exact vertex welding."""
    _, welded = np.unique(mesh.vertices, axis=0, return_inverse=True)
    edge_faces: dict[tuple[int, int], list[int]] = {}
    for face_id, face in enumerate(mesh.faces):
        welded_face = welded[face]
        for first, second in (
            (int(welded_face[0]), int(welded_face[1])),
            (int(welded_face[1]), int(welded_face[2])),
            (int(welded_face[2]), int(welded_face[0])),
        ):
            if first == second:
                continue
            edge = (first, second) if first < second else (second, first)
            edge_faces.setdefault(edge, []).append(face_id)
    adjacency = [set() for _ in range(len(mesh.faces))]
    for incident in edge_faces.values():
        for offset, first in enumerate(incident):
            for second in incident[offset + 1:]:
                adjacency[first].add(second)
                adjacency[second].add(first)
    return adjacency


def _coplanar_face_clusters(mesh: Mesh, adjacency: list[set[int]]) -> tuple[list[Cluster], list[set[int]]]:
    """Recover polygonal CAD faces that were triangulated by OBJ export."""
    triangles = mesh.vertices[mesh.faces]
    cross = np.cross(triangles[:, 1] - triangles[:, 0], triangles[:, 2] - triangles[:, 0])
    lengths = np.linalg.norm(cross, axis=1)
    normals = np.zeros_like(cross)
    valid = lengths > 0.0
    normals[valid] = cross[valid] / lengths[valid, None]
    diagonal = float(np.linalg.norm(mesh.vertices.max(axis=0) - mesh.vertices.min(axis=0)))
    plane_tolerance = max(diagonal * 1.0e-10, 1.0e-12)
    normal_tolerance = 1.0e-10

    parent = np.arange(len(mesh.faces), dtype=np.int64)

    def find(index: int) -> int:
        while int(parent[index]) != index:
            parent[index] = parent[int(parent[index])]
            index = int(parent[index])
        return index

    def union(first: int, second: int) -> None:
        first_root = find(first)
        second_root = find(second)
        if first_root != second_root:
            parent[second_root] = first_root

    for first, neighbors in enumerate(adjacency):
        if not valid[first]:
            continue
        for second in neighbors:
            if second <= first or not valid[second]:
                continue
            if float(np.dot(normals[first], normals[second])) < 1.0 - normal_tolerance:
                continue
            offsets = np.abs((triangles[second] - triangles[first, 0]) @ normals[first])
            if float(offsets.max()) <= plane_tolerance:
                union(first, second)

    grouped: dict[int, list[int]] = {}
    for face_id in range(len(mesh.faces)):
        grouped.setdefault(find(face_id), []).append(face_id)
    roots = sorted(grouped, key=lambda root: grouped[root][0])
    root_to_cluster = {root: index for index, root in enumerate(roots)}
    face_to_cluster = np.empty(len(mesh.faces), dtype=np.int64)
    clusters: list[Cluster] = []
    for root in roots:
        face_ids = np.asarray(grouped[root], dtype=np.int64)
        cluster_id = root_to_cluster[root]
        face_to_cluster[face_ids] = cluster_id
        clusters.append(Cluster(face_ids))

    cluster_adjacency = [set() for _ in clusters]
    for face_id, neighbors in enumerate(adjacency):
        first_cluster = int(face_to_cluster[face_id])
        for neighbor in neighbors:
            second_cluster = int(face_to_cluster[neighbor])
            if first_cluster != second_cluster:
                cluster_adjacency[first_cluster].add(second_cluster)
    return clusters, cluster_adjacency


def _merge_faces(
    mesh: Mesh,
    face_quadrics: FloatArray,
    options: PrimitiveDecompositionOptions,
    margin: float,
) -> tuple[list[_Node], dict[str, Any]]:
    face_adjacency = _face_adjacency(mesh)
    clusters, adjacency = _coplanar_face_clusters(mesh, face_adjacency)
    nodes: dict[int, _Node] = {}
    for cluster_id, cluster in enumerate(clusters):
        quadric = face_quadrics[cluster.triangle_ids].sum(axis=0)
        node = _make_node(cluster_id, mesh, cluster.triangle_ids, quadric, margin)
        node.neighbors = adjacency[cluster_id]
        nodes[cluster_id] = node
    next_id = len(nodes)

    bounds_min = mesh.vertices.min(axis=0)
    bounds_max = mesh.vertices.max(axis=0)
    extent = np.maximum(bounds_max - bounds_min, 0.0)
    diagonal = float(np.linalg.norm(extent))
    model_volume = max(float(np.prod(extent)), diagonal ** 3 * 1.0e-12, 1.0e-30)
    threshold_ratio = options.excess_volume_ratio()
    heap: list[tuple[float, float, int, int]] = []

    def push_pair(first_id: int, second_id: int) -> None:
        if first_id == second_id:
            return
        if first_id > second_id:
            first_id, second_id = second_id, first_id
        first = nodes[first_id]
        second = nodes[second_id]
        if not first.active or not second.active:
            return
        trial = _merged_node(-1, mesh, first, second, margin)
        excess_ratio, work_delta = _merge_key(first, second, trial, model_volume)
        heapq.heappush(heap, (excess_ratio, work_delta, first_id, second_id))

    for node in nodes.values():
        for neighbor_id in node.neighbors:
            if node.node_id < neighbor_id:
                push_pair(node.node_id, neighbor_id)

    active_count = len(nodes)
    accepted_merges = 0
    safety_forced_merges = 0
    stopped_excess_ratio: float | None = None
    while heap:
        while heap:
            excess_ratio, _, first_id, second_id = heapq.heappop(heap)
            first = nodes[first_id]
            second = nodes[second_id]
            if not first.active or not second.active or second_id not in first.neighbors:
                continue
            break
        else:
            break

        force_for_safety = active_count > options.max_parts_per_model
        if not force_for_safety and excess_ratio > threshold_ratio:
            stopped_excess_ratio = excess_ratio
            break

        merged = _merged_node(next_id, mesh, first, second, margin)
        next_id += 1
        merged.neighbors = (first.neighbors | second.neighbors) - {first_id, second_id}
        first.active = False
        second.active = False
        for neighbor_id in merged.neighbors:
            neighbor = nodes[neighbor_id]
            neighbor.neighbors.discard(first_id)
            neighbor.neighbors.discard(second_id)
            neighbor.neighbors.add(merged.node_id)
        nodes[merged.node_id] = merged
        for neighbor_id in merged.neighbors:
            push_pair(merged.node_id, neighbor_id)
        active_count -= 1
        accepted_merges += 1
        safety_forced_merges += int(force_for_safety)

    result = [node for node in nodes.values() if node.active]
    return result, {
        "source_triangles": len(mesh.faces),
        "initial_face_primitives": len(clusters),
        "accepted_merges": accepted_merges,
        "max_excess_volume_ratio": threshold_ratio,
        "next_rejected_excess_volume_ratio": stopped_excess_ratio,
        "max_parts_safety_limit": options.max_parts_per_model,
        "safety_forced_merges": safety_forced_merges,
        "topological_components_after_merging": len(result),
    }


def _select_primitive(mesh: Mesh, node: _Node, options: PrimitiveDecompositionOptions, margin: float) -> ProxyPart:
    cluster = Cluster(node.triangle_ids)
    candidates = [node.box, _expanded_obb(mesh, cluster, margin)]
    for side_count in (4, 6, 8):
        candidate = _expanded_prism(mesh, cluster, margin, side_count)
        if candidate is not None:
            candidates.append(candidate)
    diagonal = float(np.linalg.norm(mesh.vertices.max(axis=0) - mesh.vertices.min(axis=0)))
    hull = _expanded_hull(mesh, cluster, options.relative_clearance, diagonal)
    if hull is not None and len(hull.faces) <= options.max_hull_faces_per_part:
        candidates.append(hull)
    def static_key(part: ProxyPart) -> tuple[float, float, float, float, int]:
        candidate_mesh = Mesh(mesh.name, part.vertices, part.faces)
        metrics = measure_static_bvh(candidate_mesh)
        return (
            metrics.sibling_overlap_volume_sum,
            metrics.sibling_overlap_area_sum,
            metrics.sah_sum,
            part.volume,
            len(part.faces),
        )

    return min(candidates, key=static_key)


def _union_external_boundary(name: str, parts: list[ProxyPart]) -> Mesh:
    component_meshes: list[trimesh.Trimesh] = []
    for part in parts:
        component = trimesh.Trimesh(vertices=part.vertices, faces=part.faces, process=True)
        component.fix_normals(multibody=True)
        if component.volume < 0.0:
            component.invert()
        component_meshes.append(component)
    if not all(component.is_watertight for component in component_meshes):
        raise RuntimeError(f"{name}: primitive shell is not closed before union")

    bounds = np.asarray([component.bounds for component in component_meshes])
    neighbors: list[set[int]] = [set() for _ in component_meshes]
    for first in range(len(component_meshes)):
        overlaps = np.all(bounds[first, 0] <= bounds[:, 1], axis=1) & np.all(
            bounds[:, 0] <= bounds[first, 1], axis=1
        )
        for second in np.flatnonzero(overlaps):
            if first != second:
                neighbors[first].add(int(second))

    groups: list[list[int]] = []
    unvisited = set(range(len(component_meshes)))
    while unvisited:
        root = unvisited.pop()
        group = [root]
        pending = [root]
        while pending:
            current = pending.pop()
            discovered = neighbors[current] & unvisited
            unvisited.difference_update(discovered)
            pending.extend(discovered)
            group.extend(discovered)
        groups.append(group)

    external_components: list[trimesh.Trimesh] = []
    for group in groups:
        if len(group) == 1:
            external_components.append(component_meshes[group[0]])
            continue
        union = trimesh.boolean.union(
            [component_meshes[index] for index in group],
            engine="manifold",
            check_volume=False,
        )
        if union is None or len(union.faces) == 0:
            raise RuntimeError(f"{name}: primitive union produced no faces")
        external_components.append(union)
    union = trimesh.util.concatenate(external_components)
    return Mesh(
        name,
        np.asarray(union.vertices, dtype=np.float64),
        np.asarray(union.faces, dtype=np.int64),
    )


def generate_model(
    mesh: Mesh,
    options: PrimitiveDecompositionOptions,
) -> tuple[Mesh, dict[str, Any]]:
    options.validate()
    extent = mesh.vertices.max(axis=0) - mesh.vertices.min(axis=0)
    diagonal = float(np.linalg.norm(extent))
    margin = max(diagonal * options.relative_clearance, 1.0e-9)
    face_quadrics = _face_quadrics(mesh, options.tangent_quadric_weight)
    nodes, merge_analysis = _merge_faces(mesh, face_quadrics, options, margin)
    parts = [_select_primitive(mesh, node, options, margin) for node in nodes]
    minimum_clearance = min(part.minimum_clearance for part in parts)
    if not math.isfinite(minimum_clearance) or minimum_clearance <= 0.0:
        raise RuntimeError(f"{mesh.name}: strict containment certification failed")
    unmerged_proxy = _merge_parts(mesh.name, parts)
    proxy = (
        _union_external_boundary(mesh.name, parts)
        if options.union_external_boundary
        else unmerged_proxy
    )
    static_bvh = measure_static_bvh(proxy)
    methods: dict[str, int] = {}
    for part in parts:
        methods[part.method] = methods.get(part.method, 0) + 1
    stats = {
        "model": mesh.name,
        "source_vertices": len(mesh.vertices),
        "source_faces": len(mesh.faces),
        "seed_parts": merge_analysis["initial_face_primitives"],
        "proxy_parts": len(parts),
        "proxy_vertices": len(proxy.vertices),
        "proxy_faces": len(proxy.faces),
        "pre_union_vertices": len(unmerged_proxy.vertices),
        "pre_union_faces": len(unmerged_proxy.faces),
        "external_boundary_union": options.union_external_boundary,
        "methods": methods,
        "minimum_certified_clearance": minimum_clearance,
        "containment_certified": True,
        "world_sah_work": sum(_world_sah_work(part) for part in parts),
        "proxy_volume_sum": sum(part.volume for part in parts),
        "static_bvh": static_bvh.to_dict(),
        "merge_analysis": merge_analysis,
    }
    return proxy, stats


def generate_model_parts(
    mesh: Mesh,
    options: PrimitiveDecompositionOptions,
) -> tuple[list[Mesh], dict[str, Any]]:
    """Generate separately buildable conservative primitives for one model."""
    options.validate()
    extent = mesh.vertices.max(axis=0) - mesh.vertices.min(axis=0)
    diagonal = float(np.linalg.norm(extent))
    margin = max(diagonal * options.relative_clearance, 1.0e-9)
    face_quadrics = _face_quadrics(mesh, options.tangent_quadric_weight)
    nodes, merge_analysis = _merge_faces(mesh, face_quadrics, options, margin)
    parts = [_select_primitive(mesh, node, options, margin) for node in nodes]
    minimum_clearance = min(part.minimum_clearance for part in parts)
    if not math.isfinite(minimum_clearance) or minimum_clearance <= 0.0:
        raise RuntimeError(f"{mesh.name}: strict containment certification failed")

    methods: dict[str, int] = {}
    part_meshes: list[Mesh] = []
    static_metrics: list[dict[str, Any]] = []
    for index, part in enumerate(parts):
        methods[part.method] = methods.get(part.method, 0) + 1
        part_mesh = Mesh(f"{mesh.name}:part_{index:04d}", part.vertices, part.faces)
        part_meshes.append(part_mesh)
        static_metrics.append(measure_static_bvh(part_mesh).to_dict())

    stats = {
        "model": mesh.name,
        "source_vertices": len(mesh.vertices),
        "source_faces": len(mesh.faces),
        "seed_parts": merge_analysis["initial_face_primitives"],
        "proxy_parts": len(parts),
        "proxy_vertices": sum(len(part.vertices) for part in parts),
        "proxy_faces": sum(len(part.faces) for part in parts),
        "external_boundary_union": False,
        "independent_blas": True,
        "methods": methods,
        "minimum_certified_clearance": minimum_clearance,
        "containment_certified": True,
        "world_sah_work": sum(_world_sah_work(part) for part in parts),
        "proxy_volume_sum": sum(part.volume for part in parts),
        "part_static_bvhs": static_metrics,
        "merge_analysis": merge_analysis,
    }
    return part_meshes, stats


def generate_pool(
    input_dir: Path,
    output_dir: Path,
    options: PrimitiveDecompositionOptions,
    model_ids: Iterable[int] = OFFICIAL_SIMPLIFIED_MODEL_IDS,
) -> dict[str, Any]:
    options.validate()
    selected = frozenset(model_ids)
    input_paths = sorted(input_dir.glob("*.obj"), key=_numeric_name)
    if not input_paths:
        raise ValueError(f"no OBJ files in {input_dir}")
    if output_dir.exists() and any(output_dir.iterdir()):
        raise FileExistsError(f"output directory is not empty: {output_dir}")
    output_dir.mkdir(parents=True, exist_ok=True)

    stats: list[dict[str, Any]] = []
    started = time.perf_counter()
    for input_path in input_paths:
        model_id = _numeric_name(input_path)
        source = read_obj(input_path)
        model_started = time.perf_counter()
        if model_id in selected:
            proxy, model_stats = generate_model(source, options)
            write_obj(output_dir / input_path.name, proxy, ["mode=PQSSPrimitiveDecomposition"])
        else:
            shutil.copy2(input_path, output_dir / input_path.name)
            model_stats = {
                "model": source.name,
                "source_vertices": len(source.vertices),
                "source_faces": len(source.faces),
                "proxy_vertices": len(source.vertices),
                "proxy_faces": len(source.faces),
                "proxy_parts": 1,
                "methods": {"exact_passthrough": 1},
                "containment_certified": True,
            }
        model_stats["selected"] = model_id in selected
        model_stats["generation_seconds"] = time.perf_counter() - model_started
        stats.append(model_stats)

    report: dict[str, Any] = {
        "algorithm": "PQSSWorkloadOrientedPrimitiveDecomposition",
        "input_directory": str(input_dir.resolve()),
        "output_directory": str(output_dir.resolve()),
        "selected_model_ids": sorted(selected),
        "options": asdict(options),
        "model_count": len(stats),
        "selected_model_count": len(selected),
        "generation_seconds": time.perf_counter() - started,
        "models": stats,
    }
    (output_dir / "generation_report.json").write_text(
        json.dumps(report, indent=2, ensure_ascii=True) + "\n", encoding="ascii"
    )
    return report


def generate_independent_pool(
    input_dir: Path,
    output_dir: Path,
    options: PrimitiveDecompositionOptions,
    model_ids: Iterable[int] = OFFICIAL_SIMPLIFIED_MODEL_IDS,
) -> dict[str, Any]:
    """Write one OBJ per primitive plus a logical-model manifest."""
    options.validate()
    selected = frozenset(model_ids)
    input_paths = sorted(input_dir.glob("*.obj"), key=_numeric_name)
    if not input_paths:
        raise ValueError(f"no OBJ files in {input_dir}")
    if output_dir.exists() and any(output_dir.iterdir()):
        raise FileExistsError(f"output directory is not empty: {output_dir}")
    output_dir.mkdir(parents=True, exist_ok=True)

    stats: list[dict[str, Any]] = []
    manifest_rows = ["logical_model_id\tpart_index\trelative_path\tselected"]
    started = time.perf_counter()
    for input_path in input_paths:
        model_id = _numeric_name(input_path)
        source = read_obj(input_path)
        model_started = time.perf_counter()
        model_directory = output_dir / "models" / str(model_id)
        model_directory.mkdir(parents=True, exist_ok=True)
        if model_id in selected:
            part_meshes, model_stats = generate_model_parts(source, options)
            for part_index, part_mesh in enumerate(part_meshes):
                relative_path = Path("models") / str(model_id) / f"part_{part_index:04d}.obj"
                write_obj(
                    output_dir / relative_path,
                    part_mesh,
                    [
                        "mode=PQSSIndependentPrimitive",
                        f"logical_model_id={model_id}",
                        f"part_index={part_index}",
                    ],
                )
                manifest_rows.append(f"{model_id}\t{part_index}\t{relative_path.as_posix()}\t1")
        else:
            relative_path = Path("models") / str(model_id) / "part_0000.obj"
            shutil.copy2(input_path, output_dir / relative_path)
            manifest_rows.append(f"{model_id}\t0\t{relative_path.as_posix()}\t0")
            model_stats = {
                "model": source.name,
                "source_vertices": len(source.vertices),
                "source_faces": len(source.faces),
                "proxy_vertices": len(source.vertices),
                "proxy_faces": len(source.faces),
                "proxy_parts": 1,
                "independent_blas": True,
                "methods": {"exact_passthrough": 1},
                "containment_certified": True,
            }
        model_stats["selected"] = model_id in selected
        model_stats["generation_seconds"] = time.perf_counter() - model_started
        stats.append(model_stats)

    (output_dir / "primitive_pool_manifest.tsv").write_text(
        "\n".join(manifest_rows) + "\n", encoding="ascii"
    )
    report: dict[str, Any] = {
        "algorithm": "PQSSIndependentPrimitiveDecomposition",
        "input_directory": str(input_dir.resolve()),
        "output_directory": str(output_dir.resolve()),
        "manifest": "primitive_pool_manifest.tsv",
        "selected_model_ids": sorted(selected),
        "analysis_strength": options.analysis_strength,
        "max_excess_volume_ratio": options.excess_volume_ratio(),
        "options": asdict(options),
        "model_count": len(stats),
        "physical_model_count": sum(int(model["proxy_parts"]) for model in stats),
        "selected_model_count": len(selected),
        "generation_seconds": time.perf_counter() - started,
        "models": stats,
    }
    (output_dir / "generation_report.json").write_text(
        json.dumps(report, indent=2, ensure_ascii=True) + "\n", encoding="ascii"
    )
    return report
