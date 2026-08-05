from __future__ import annotations

from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import asdict, dataclass
import heapq
import json
import math
import os
from pathlib import Path
import time
from typing import Any, Iterable, Sequence

import numpy as np
from shapely import MultiPoint, Point, Polygon, STRtree, constrained_delaunay_triangles, unary_union
import trimesh

from .obj_io import FloatArray, IndexArray, Mesh, read_obj, write_obj
from .primitive_decomposition import (
    _coplanar_face_clusters,
    _face_adjacency,
    _face_quadrics,
    _quadric_axes,
)


ANALYTIC_PRIMITIVE_TYPES = ("sphere", "capsule", "rss")
PRIMITIVE_TYPES = (*ANALYTIC_PRIMITIVE_TYPES, "triangle")

_DEGENERATE_NORMALIZED_AREA_LIMIT = 128.0 * np.finfo(np.float64).eps


@dataclass(slots=True)
class PrimitiveAnalysisOptions:
    analysis_strength: float = 0.5
    max_excess_volume_ratio: float | None = None
    primitive_types: tuple[str, ...] = PRIMITIVE_TYPES
    relative_clearance: float = 1.0e-6
    tangent_quadric_weight: float = 1.0e-2
    capsule_radius_samples: int = 18
    visualization_sphere_subdivisions: int = 1
    model_workers: int = 0
    planar_triangle_max_coverage: float = 0.98
    planar_relative_tolerance: float = 1.0e-9
    planar_triangle_min_source_faces: int = 4
    planar_boundary_void_max_edge_relative: float = 0.03
    planar_boundary_void_max_area_ratio: float = 0.01
    projection_hole_max_thickness_ratio: float = 0.35
    projection_hole_min_area_ratio: float = 0.01
    planar_sharp_edge_degrees: float = 30.0

    def validate(self) -> None:
        if not math.isfinite(self.analysis_strength) or not 0.0 <= self.analysis_strength <= 1.0:
            raise ValueError("analysis_strength must be in [0, 1]")
        if self.max_excess_volume_ratio is not None and (
            not math.isfinite(self.max_excess_volume_ratio) or self.max_excess_volume_ratio < 0.0
        ):
            raise ValueError("max_excess_volume_ratio must be finite and non-negative")
        if not self.primitive_types:
            raise ValueError("at least one primitive type must be enabled")
        unknown = sorted(set(self.primitive_types) - set(PRIMITIVE_TYPES))
        if unknown:
            raise ValueError(f"unsupported primitive types: {', '.join(unknown)}")
        if len(set(self.primitive_types)) != len(self.primitive_types):
            raise ValueError("primitive types must not be repeated")
        if not set(self.primitive_types).intersection(ANALYTIC_PRIMITIVE_TYPES):
            raise ValueError("at least one sphere, capsule, or rss fallback must be enabled")
        if not math.isfinite(self.relative_clearance) or self.relative_clearance <= 0.0:
            raise ValueError("relative_clearance must be finite and positive")
        if not math.isfinite(self.tangent_quadric_weight) or self.tangent_quadric_weight < 0.0:
            raise ValueError("tangent_quadric_weight must be finite and non-negative")
        if self.capsule_radius_samples < 2:
            raise ValueError("capsule_radius_samples must be at least 2")
        if not 0 <= self.visualization_sphere_subdivisions <= 3:
            raise ValueError("visualization_sphere_subdivisions must be in [0, 3]")
        if self.model_workers < 0:
            raise ValueError("model_workers must be non-negative; zero selects all available CPU cores")
        if not 0.0 <= self.planar_triangle_max_coverage <= 1.0:
            raise ValueError("planar_triangle_max_coverage must be in [0, 1]")
        if not math.isfinite(self.planar_relative_tolerance) or self.planar_relative_tolerance <= 0.0:
            raise ValueError("planar_relative_tolerance must be finite and positive")
        if self.planar_triangle_min_source_faces < 2:
            raise ValueError("planar_triangle_min_source_faces must be at least 2")
        if not 0.0 <= self.planar_boundary_void_max_edge_relative <= 1.0:
            raise ValueError("planar_boundary_void_max_edge_relative must be in [0, 1]")
        if not 0.0 <= self.planar_boundary_void_max_area_ratio <= 1.0:
            raise ValueError("planar_boundary_void_max_area_ratio must be in [0, 1]")
        if not 0.0 <= self.projection_hole_max_thickness_ratio <= 1.0:
            raise ValueError("projection_hole_max_thickness_ratio must be in [0, 1]")
        if not 0.0 <= self.projection_hole_min_area_ratio <= 1.0:
            raise ValueError("projection_hole_min_area_ratio must be in [0, 1]")
        if not 0.0 <= self.planar_sharp_edge_degrees <= 90.0:
            raise ValueError("planar_sharp_edge_degrees must be in [0, 90]")

    def excess_volume_ratio(self) -> float:
        if self.max_excess_volume_ratio is not None:
            return self.max_excess_volume_ratio
        if self.analysis_strength == 0.0:
            return math.inf
        # 0.5 maps to 1e-2, the paper's useful middle example.
        return 10.0 ** (1.0 - 6.0 * self.analysis_strength)


@dataclass(slots=True)
class Primitive:
    kind: str
    origin: FloatArray
    axes: FloatArray
    lengths: FloatArray
    radius: float
    volume: float
    triangle_vertices: FloatArray | None = None
    mean_sampled_outward_deviation: float = 0.0
    maximum_sampled_outward_deviation: float = 0.0
    planar_excess_area_ratio: float | None = None

    def to_dict(self) -> dict[str, Any]:
        return {
            "type": self.kind,
            "origin": self.origin.tolist(),
            "axes": self.axes.tolist(),
            "lengths": self.lengths.tolist(),
            "radius": self.radius,
            "volume": self.volume,
            "triangle_vertices": (
                self.triangle_vertices.tolist() if self.triangle_vertices is not None else None
            ),
            "mean_sampled_outward_deviation": self.mean_sampled_outward_deviation,
            "maximum_sampled_outward_deviation": self.maximum_sampled_outward_deviation,
            "planar_excess_area_ratio": self.planar_excess_area_ratio,
        }


@dataclass(slots=True)
class PrimitiveRegion:
    primitive: Primitive
    triangle_ids: IndexArray
    vertex_ids: IndexArray
    outward_reference_triangles: FloatArray | None = None
    synthetic_fill: bool = False


@dataclass(slots=True)
class PrimitiveAnalysisResult:
    mesh: Mesh
    regions: list[PrimitiveRegion]
    stats: dict[str, Any]
    excluded_triangle_ids: IndexArray


@dataclass(slots=True)
class _Node:
    node_id: int
    triangle_ids: IndexArray
    vertex_ids: IndexArray
    quadric: FloatArray
    primitive: Primitive
    neighbors: set[int]
    active: bool = True


def _frame_from_first_axis(axis: FloatArray) -> FloatArray:
    first = np.asarray(axis, dtype=np.float64)
    first /= max(float(np.linalg.norm(first)), 1.0e-30)
    helper = np.asarray((1.0, 0.0, 0.0)) if abs(float(first[0])) < 0.8 else np.asarray((0.0, 1.0, 0.0))
    second = np.cross(first, helper)
    second /= max(float(np.linalg.norm(second)), 1.0e-30)
    third = np.cross(first, second)
    return np.column_stack((first, second, third))


def _point_core_distances(points: FloatArray, primitive: Primitive) -> FloatArray:
    local = (points - primitive.origin) @ primitive.axes
    closest_x = np.clip(local[:, 0], 0.0, primitive.lengths[0])
    closest_y = np.clip(local[:, 1], 0.0, primitive.lengths[1])
    delta = np.column_stack((local[:, 0] - closest_x, local[:, 1] - closest_y, local[:, 2]))
    return np.linalg.norm(delta, axis=1)


def _nearest_source_distances(points: FloatArray, source_triangles: FloatArray) -> FloatArray:
    source = trimesh.Trimesh(
        vertices=source_triangles.reshape(-1, 3),
        faces=np.arange(source_triangles.size // 3, dtype=np.int64).reshape(-1, 3),
        process=False,
    )
    _, distances, _ = trimesh.proximity.closest_point_naive(source, points)
    return np.asarray(distances, dtype=np.float64)


def _region_outward_deviation_metrics(
    mesh: Mesh,
    region: PrimitiveRegion,
) -> tuple[float, float, float]:
    if region.primitive.kind == "triangle":
        triangle = region.primitive.triangle_vertices
        if triangle is None:
            raise RuntimeError("triangle primitive has no vertices")
        area = 0.5 * float(np.linalg.norm(np.cross(triangle[1] - triangle[0], triangle[2] - triangle[0])))
        return 0.0, 0.0, area

    visual = primitive_visualization_mesh(region.primitive, subdivisions=2)
    proxy_triangles = visual.vertices[visual.faces]
    proxy_areas = 0.5 * np.linalg.norm(
        np.cross(proxy_triangles[:, 1] - proxy_triangles[:, 0],
                 proxy_triangles[:, 2] - proxy_triangles[:, 0]),
        axis=1,
    )
    centroids = proxy_triangles.mean(axis=1)
    source_triangles = (
        region.outward_reference_triangles
        if region.outward_reference_triangles is not None
        else mesh.vertices[mesh.faces[region.triangle_ids]]
    )
    centroid_distances = _nearest_source_distances(centroids, source_triangles)
    area_sum = float(proxy_areas.sum())
    mean_deviation = (
        float(np.dot(proxy_areas, centroid_distances) / area_sum)
        if area_sum > 0.0 else float(centroid_distances.mean())
    )
    maximum_points = np.vstack((visual.vertices, centroids))
    maximum_deviation = float(_nearest_source_distances(maximum_points, source_triangles).max())
    region.primitive.mean_sampled_outward_deviation = mean_deviation
    region.primitive.maximum_sampled_outward_deviation = maximum_deviation
    return mean_deviation, maximum_deviation, area_sum


def _apply_outward_deviation_metrics(mesh: Mesh, regions: Sequence[PrimitiveRegion]) -> dict[str, Any]:
    weighted_deviation = 0.0
    proxy_surface_area = 0.0
    maximum_deviation = 0.0
    maximum_planar_excess = 0.0
    for region in regions:
        mean_deviation, region_maximum, area = _region_outward_deviation_metrics(mesh, region)
        weighted_deviation += mean_deviation * area
        proxy_surface_area += area
        maximum_deviation = max(maximum_deviation, region_maximum)
        if region.primitive.planar_excess_area_ratio is not None:
            maximum_planar_excess = max(
                maximum_planar_excess, region.primitive.planar_excess_area_ratio
            )

    mean_deviation = weighted_deviation / max(proxy_surface_area, 1.0e-30)
    diagonal = float(np.linalg.norm(mesh.vertices.max(axis=0) - mesh.vertices.min(axis=0)))
    return {
        "mean_sampled_outward_deviation": mean_deviation,
        "maximum_sampled_outward_deviation": maximum_deviation,
        "mean_relative_outward_deviation": mean_deviation / max(diagonal, 1.0e-30),
        "maximum_relative_outward_deviation": maximum_deviation / max(diagonal, 1.0e-30),
        "maximum_planar_excess_area_ratio": maximum_planar_excess,
        "outward_deviation_sampling": "proxy vertices and face centroids",
        "proxy_sampled_surface_area": proxy_surface_area,
    }


def _degenerate_face_mask(mesh: Mesh) -> np.ndarray:
    triangles = mesh.vertices[mesh.faces]
    edges = np.stack(
        (
            triangles[:, 1] - triangles[:, 0],
            triangles[:, 2] - triangles[:, 1],
            triangles[:, 0] - triangles[:, 2],
        ),
        axis=1,
    )
    maximum_edge_squared = np.max(np.einsum("tec,tec->te", edges, edges), axis=1)
    double_areas = np.linalg.norm(np.cross(edges[:, 0], -edges[:, 2]), axis=1)
    if not np.all(np.isfinite(double_areas)):
        raise ValueError(f"{mesh.name}: triangle areas contain non-finite values")
    return double_areas <= _DEGENERATE_NORMALIZED_AREA_LIMIT * maximum_edge_squared


def _drop_degenerate_faces(mesh: Mesh) -> tuple[Mesh, int, int]:
    original_count = len(mesh.faces)
    keep = ~_degenerate_face_mask(mesh)
    if not np.any(keep):
        raise ValueError(f"{mesh.name}: all triangles are degenerate")
    cleaned = Mesh(mesh.name, mesh.vertices, mesh.faces[keep])
    return cleaned, original_count, int(np.count_nonzero(~keep))


def _sphere_volume(radius: float) -> float:
    return 4.0 * math.pi * radius ** 3 / 3.0


def _capsule_volume(length: float, radius: float) -> float:
    return math.pi * radius * radius * length + _sphere_volume(radius)


def _rss_volume(lengths: FloatArray, radius: float) -> float:
    return (
        2.0 * radius * float(lengths[0] * lengths[1])
        + math.pi * radius * radius * float(lengths[0] + lengths[1])
        + _sphere_volume(radius)
    )


def _fit_sphere(points: FloatArray, margin: float) -> Primitive:
    bounds_center = 0.5 * (points.min(axis=0) + points.max(axis=0))
    candidates = [points.mean(axis=0), bounds_center]
    if len(points) > 1:
        seed = points[0]
        first = points[int(np.argmax(np.linalg.norm(points - seed, axis=1)))]
        second = points[int(np.argmax(np.linalg.norm(points - first, axis=1)))]
        candidates.append(0.5 * (first + second))
    center = min(candidates, key=lambda value: float(np.linalg.norm(points - value, axis=1).max()))
    radius = float(np.linalg.norm(points - center, axis=1).max()) + margin
    primitive = Primitive(
        "sphere", np.asarray(center), np.eye(3), np.zeros(2), radius,
        _sphere_volume(radius),
    )
    return primitive


def _unique_axes(candidates: Iterable[FloatArray]) -> list[FloatArray]:
    result: list[FloatArray] = []
    for candidate in candidates:
        length = float(np.linalg.norm(candidate))
        if length <= 1.0e-12:
            continue
        axis = np.asarray(candidate, dtype=np.float64) / length
        first_nonzero = next((index for index in range(3) if abs(float(axis[index])) > 1.0e-12), None)
        if first_nonzero is not None and axis[first_nonzero] < 0.0:
            axis = -axis
        if any(abs(float(np.dot(axis, existing))) > 0.9999 for existing in result):
            continue
        result.append(axis)
    return result


def _candidate_axes(points: FloatArray, quadric: FloatArray) -> list[FloatArray]:
    centered = points - points.mean(axis=0)
    covariance = centered.T @ centered
    _, pca = np.linalg.eigh(covariance)
    quadric_axes = _quadric_axes(quadric)
    candidates = [pca[:, index] for index in range(3)]
    candidates.extend(quadric_axes[:, index] for index in range(3))
    candidates.extend(np.eye(3)[:, index] for index in range(3))
    return _unique_axes(candidates)


def _fit_capsule_on_axis(
    points: FloatArray,
    axis: FloatArray,
    margin: float,
    radius_samples: int,
) -> Primitive:
    reference = points.mean(axis=0)
    relative = points - reference
    axial = relative @ axis
    radial_vectors = relative - axial[:, None] * axis
    radial = np.linalg.norm(radial_vectors, axis=1)
    minimum_radius = float(radial.max())
    sphere_radius = float(np.linalg.norm(relative, axis=1).max())
    maximum_radius = max(minimum_radius, sphere_radius)
    if maximum_radius - minimum_radius <= max(maximum_radius, 1.0) * 1.0e-12:
        radii = np.asarray([minimum_radius])
    else:
        radii = np.linspace(minimum_radius, maximum_radius, radius_samples)

    cap = np.sqrt(np.maximum(radii[:, None] ** 2 - radial[None, :] ** 2, 0.0))
    lowers = np.min(axial[None, :] + cap, axis=1)
    uppers = np.max(axial[None, :] - cap, axis=1)
    overlapping = uppers < lowers
    midpoints = 0.5 * (lowers + uppers)
    lowers = np.where(overlapping, midpoints, lowers)
    uppers = np.where(overlapping, midpoints, uppers)
    effective_radii = radii + margin
    lengths = uppers - lowers
    volumes = math.pi * effective_radii ** 2 * lengths + (
        4.0 * math.pi * effective_radii ** 3 / 3.0
    )
    # Preserve the scalar implementation's tuple ordering for exact tie breaks.
    best_index = int(np.lexsort((effective_radii, uppers, lowers, volumes))[0])
    volume = float(volumes[best_index])
    lower = float(lowers[best_index])
    upper = float(uppers[best_index])
    radius = float(effective_radii[best_index])
    frame = _frame_from_first_axis(axis)
    origin = reference + lower * axis
    primitive = Primitive(
        "capsule", origin, frame, np.asarray((upper - lower, 0.0)),
        radius, volume,
    )
    return primitive


def _fit_capsule(
    points: FloatArray,
    axes: Sequence[FloatArray],
    margin: float,
    radius_samples: int,
) -> Primitive:
    return min(
        (_fit_capsule_on_axis(points, axis, margin, radius_samples)
         for axis in axes),
        key=lambda primitive: primitive.volume,
    )


def _fit_rss_on_frame(points: FloatArray, frame: FloatArray, margin: float) -> Primitive:
    local = points @ frame
    bounds_min = local.min(axis=0)
    bounds_max = local.max(axis=0)
    radius = 0.5 * float(bounds_max[2] - bounds_min[2]) + margin
    origin_local = np.asarray((bounds_min[0], bounds_min[1], 0.5 * (bounds_min[2] + bounds_max[2])))
    lengths = np.asarray((bounds_max[0] - bounds_min[0], bounds_max[1] - bounds_min[1]))
    primitive = Primitive(
        "rss", origin_local @ frame.T, frame, lengths, radius,
        _rss_volume(lengths, radius),
    )
    return primitive


def _fit_rss(
    points: FloatArray,
    quadric: FloatArray,
    axes: Sequence[FloatArray],
    margin: float,
) -> Primitive:
    frames: list[FloatArray] = []
    for normal in axes:
        frame = _frame_from_first_axis(normal)
        # RSS uses columns 0/1 as its rectangle and column 2 as normal.
        frames.append(np.column_stack((frame[:, 1], frame[:, 2], frame[:, 0])))
    quadric_frame = _quadric_axes(quadric)
    frames.extend(
        np.column_stack((quadric_frame[:, (normal + 1) % 3],
                         quadric_frame[:, (normal + 2) % 3],
                         quadric_frame[:, normal]))
        for normal in range(3)
    )
    return min((_fit_rss_on_frame(points, frame, margin) for frame in frames),
               key=lambda primitive: primitive.volume)


def _fit_best(
    points: FloatArray,
    quadric: FloatArray,
    margin: float,
    options: PrimitiveAnalysisOptions,
) -> Primitive:
    candidates: list[Primitive] = []
    axes = (
        _candidate_axes(points, quadric)
        if "capsule" in options.primitive_types or "rss" in options.primitive_types
        else []
    )
    for kind in options.primitive_types:
        if kind == "sphere":
            candidates.append(_fit_sphere(points, margin))
        elif kind == "capsule":
            candidates.append(_fit_capsule(points, axes, margin, options.capsule_radius_samples))
        elif kind == "rss":
            candidates.append(_fit_rss(points, quadric, axes, margin))
    order = {kind: index for index, kind in enumerate(options.primitive_types)}
    return min(candidates, key=lambda primitive: (primitive.volume, order[primitive.kind]))


def _make_node(
    node_id: int,
    mesh: Mesh,
    triangle_ids: IndexArray,
    quadric: FloatArray,
    neighbors: set[int],
    margin: float,
    options: PrimitiveAnalysisOptions,
) -> _Node:
    vertex_ids = np.unique(mesh.faces[triangle_ids].reshape(-1))
    primitive = _fit_best(mesh.vertices[vertex_ids], quadric, margin, options)
    return _Node(node_id, triangle_ids, vertex_ids, quadric, primitive, set(neighbors))


def _merged_node(
    node_id: int,
    mesh: Mesh,
    first: _Node,
    second: _Node,
    margin: float,
    options: PrimitiveAnalysisOptions,
) -> _Node:
    triangle_ids = np.concatenate((first.triangle_ids, second.triangle_ids))
    vertex_ids = np.union1d(first.vertex_ids, second.vertex_ids)
    quadric = first.quadric + second.quadric
    primitive = _fit_best(mesh.vertices[vertex_ids], quadric, margin, options)
    return _Node(node_id, triangle_ids, vertex_ids, quadric, primitive, set())


@dataclass(slots=True)
class _PlanarPatch:
    origin: FloatArray
    frame: FloatArray
    geometry: Any
    proxy_coverage: float
    proxy_excess_ratio: float
    filled_hole_count: int
    filled_hole_area: float
    filled_boundary_void_count: int
    filled_boundary_void_area: float
    filled_void_geometry: Any
    tolerance: float
    source_geometry: Any


@dataclass(slots=True)
class _FilledVoid:
    origin: FloatArray
    frame: FloatArray
    geometry: Any
    targets: list[PrimitiveRegion]
    tolerance: float


def _polygons(geometry: Any) -> list[Any]:
    if geometry.geom_type == "Polygon":
        return [geometry]
    if geometry.geom_type == "MultiPolygon":
        return list(geometry.geoms)
    return [part for part in geometry.geoms if part.geom_type == "Polygon"]


def _fill_planar_holes(geometry: Any) -> tuple[Any, int, float]:
    polygons = _polygons(geometry)
    hole_count = sum(len(polygon.interiors) for polygon in polygons)
    if hole_count == 0:
        return geometry, 0, 0.0
    filled = unary_union([Polygon(polygon.exterior) for polygon in polygons])
    return filled, hole_count, max(float(filled.area) - float(geometry.area), 0.0)


def _line_intersection(
    first_start: FloatArray,
    first_end: FloatArray,
    second_start: FloatArray,
    second_end: FloatArray,
) -> FloatArray | None:
    first = first_end - first_start
    second = second_end - second_start
    denominator = float(first[0] * second[1] - first[1] * second[0])
    if abs(denominator) <= 1.0e-14 * max(float(np.linalg.norm(first) * np.linalg.norm(second)), 1.0):
        return None
    delta = second_start - first_start
    parameter = float(delta[0] * second[1] - delta[1] * second[0]) / denominator
    return first_start + parameter * first


def _fill_short_boundary_voids(
    geometry: Any,
    tolerance: float,
    options: PrimitiveAnalysisOptions,
) -> tuple[Any, int, float]:
    output: list[Any] = []
    filled_count = 0
    filled_area = 0.0
    for source_polygon in _polygons(geometry):
        polygon = source_polygon
        while True:
            coordinates = np.asarray(polygon.exterior.coords[:-1], dtype=np.float64)
            if len(coordinates) < 4:
                break
            diagonal = float(np.linalg.norm(np.ptp(coordinates, axis=0)))
            maximum_edge = diagonal * options.planar_boundary_void_max_edge_relative
            maximum_added_area = (
                float(source_polygon.area) * options.planar_boundary_void_max_area_ratio
            )
            lengths = np.linalg.norm(np.roll(coordinates, -1, axis=0) - coordinates, axis=1)
            accepted: tuple[Any, float] | None = None
            for edge_index, edge_length in enumerate(lengths):
                if edge_length > maximum_edge or edge_length <= tolerance:
                    continue
                rotated = np.vstack(
                    (coordinates[edge_index:], coordinates[:edge_index])
                )
                previous, first, second, following = (
                    rotated[-1], rotated[0], rotated[1], rotated[2]
                )
                if min(float(np.linalg.norm(first - previous)),
                       float(np.linalg.norm(following - second))) < 3.0 * edge_length:
                    continue
                intersection = _line_intersection(previous, first, second, following)
                if intersection is None:
                    continue
                if max(float(np.linalg.norm(intersection - first)),
                       float(np.linalg.norm(intersection - second))) > 2.0 * edge_length:
                    continue
                replacement = np.vstack((intersection, rotated[2:]))
                candidate = Polygon(
                    replacement,
                    [np.asarray(ring.coords) for ring in polygon.interiors],
                )
                if not candidate.is_valid:
                    continue
                added_area = float(candidate.area - polygon.area)
                if added_area <= tolerance * tolerance or added_area > maximum_added_area:
                    continue
                if not candidate.buffer(tolerance).covers(polygon):
                    continue
                accepted = candidate, added_area
                break
            if accepted is None:
                break
            polygon, added_area = accepted
            filled_count += 1
            filled_area += added_area
        output.append(polygon)
    return unary_union(output), filled_count, filled_area


def _planar_patch(
    mesh: Mesh,
    region: PrimitiveRegion,
    options: PrimitiveAnalysisOptions,
) -> _PlanarPatch | None:
    triangles = mesh.vertices[mesh.faces[region.triangle_ids]]
    points = np.unique(triangles.reshape(-1, 3), axis=0)
    origin = points.mean(axis=0)
    centered = points - origin
    eigenvalues, eigenvectors = np.linalg.eigh(centered.T @ centered)
    normal = eigenvectors[:, int(np.argmin(eigenvalues))]
    first = eigenvectors[:, int(np.argmax(eigenvalues))]
    second = np.cross(normal, first)
    second /= max(float(np.linalg.norm(second)), 1.0e-30)
    first = np.cross(second, normal)
    frame = np.column_stack((first, second, normal))
    local_triangles = (triangles - origin) @ frame
    diagonal = float(np.linalg.norm(points.max(axis=0) - points.min(axis=0)))
    tolerance = max(diagonal * options.planar_relative_tolerance, 1.0e-10)
    if float(np.abs(local_triangles[:, :, 2]).max()) > tolerance:
        return None

    polygons = [Polygon(triangle[:, :2]) for triangle in local_triangles]
    geometry = unary_union([polygon for polygon in polygons if polygon.area > 0.0])
    if geometry.is_empty or geometry.area <= 0.0:
        return None
    original_area = float(geometry.area)
    geometry = geometry.simplify(tolerance, preserve_topology=True)
    area_tolerance = max(original_area * 1.0e-9, tolerance * tolerance)
    if geometry.is_empty or abs(float(geometry.area) - original_area) > area_tolerance:
        return None

    source_geometry = geometry
    geometry, boundary_void_count, boundary_void_area = _fill_short_boundary_voids(
        geometry, tolerance, options
    )
    geometry, hole_count, filled_hole_area = _fill_planar_holes(geometry)
    filled_void_geometry = geometry.difference(source_geometry)
    proxy_visual = primitive_visualization_mesh(region.primitive, subdivisions=2)
    proxy_local = (proxy_visual.vertices - origin) @ frame
    proxy_area = float(MultiPoint(proxy_local[:, :2]).convex_hull.area)
    if proxy_area <= 0.0:
        return None
    coverage = min(float(geometry.area) / proxy_area, 1.0)
    excess_ratio = max(proxy_area - float(geometry.area), 0.0) / float(geometry.area)
    return _PlanarPatch(
        origin,
        frame,
        geometry,
        coverage,
        excess_ratio,
        hole_count,
        filled_hole_area,
        boundary_void_count,
        boundary_void_area,
        filled_void_geometry,
        tolerance,
        source_geometry,
    )


def _triangulate_planar_patch(patch: _PlanarPatch) -> list[FloatArray]:
    triangulation = constrained_delaunay_triangles(patch.geometry)
    result: list[FloatArray] = []
    area_sum = 0.0
    for candidate in triangulation.geoms:
        if candidate.geom_type != "Polygon" or candidate.area <= 0.0:
            continue
        if not patch.geometry.covers(candidate):
            continue
        coordinates = np.asarray(candidate.exterior.coords[:3], dtype=np.float64)
        local = np.column_stack((coordinates, np.zeros(3)))
        result.append(patch.origin + local @ patch.frame.T)
        area_sum += float(candidate.area)
    area_tolerance = max(float(patch.geometry.area) * 1.0e-9, 1.0e-12)
    if not result or abs(area_sum - float(patch.geometry.area)) > area_tolerance:
        raise RuntimeError(
            f"constrained planar triangulation area mismatch: {area_sum} != {patch.geometry.area}"
        )
    return result


def _triangle_primitive(vertices: FloatArray) -> Primitive:
    first = vertices[1] - vertices[0]
    first /= max(float(np.linalg.norm(first)), 1.0e-30)
    normal = np.cross(vertices[1] - vertices[0], vertices[2] - vertices[0])
    normal /= max(float(np.linalg.norm(normal)), 1.0e-30)
    second = np.cross(normal, first)
    return Primitive(
        "triangle",
        vertices[0].copy(),
        np.column_stack((first, second, normal)),
        np.zeros(2),
        0.0,
        0.0,
        triangle_vertices=vertices.copy(),
        planar_excess_area_ratio=0.0,
    )


def _assign_source_triangles(
    mesh: Mesh,
    region: PrimitiveRegion,
    patch: _PlanarPatch,
    patch_triangles: Sequence[FloatArray],
) -> list[IndexArray]:
    polygons = [
        Polygon(((triangle - patch.origin) @ patch.frame)[:, :2])
        for triangle in patch_triangles
    ]
    tree = STRtree(polygons)
    assignments: list[list[int]] = [[] for _ in patch_triangles]
    centroids = mesh.vertices[mesh.faces[region.triangle_ids]].mean(axis=1)
    local_centroids = (centroids - patch.origin) @ patch.frame
    for triangle_id, local in zip(region.triangle_ids, local_centroids, strict=True):
        point = Point(float(local[0]), float(local[1]))
        candidates = tree.query(point)
        selected = next(
            (int(index) for index in candidates if polygons[int(index)].covers(point)),
            None,
        )
        if selected is None:
            selected = min(range(len(polygons)), key=lambda index: polygons[index].distance(point))
        assignments[selected].append(int(triangle_id))
    return [np.asarray(ids, dtype=np.int64) for ids in assignments]


def _projection_hole_caps(
    mesh: Mesh,
    regions: list[PrimitiveRegion],
    options: PrimitiveAnalysisOptions,
) -> tuple[list[_FilledVoid], dict[str, int | float]]:
    stats: dict[str, int | float] = {
        "projection_filled_holes": 0,
        "projection_fill_primitives": 0,
        "projection_filled_area": 0.0,
    }
    if "triangle" not in options.primitive_types:
        return [], stats
    extent = np.ptp(mesh.vertices, axis=0)
    ordered_axes = np.argsort(extent)
    thin_axis = int(ordered_axes[0])
    reference_extent = float(extent[int(ordered_axes[1])])
    if (
        reference_extent <= 0.0
        or float(extent[thin_axis]) <= float(np.linalg.norm(extent)) * 1.0e-9
        or float(extent[thin_axis]) / reference_extent > (
        options.projection_hole_max_thickness_ratio
        )
    ):
        return [], stats

    plane_axes = [axis for axis in range(3) if axis != thin_axis]
    projected_triangles = mesh.vertices[mesh.faces][:, :, plane_axes]
    polygons = [Polygon(triangle) for triangle in projected_triangles]
    geometry = unary_union([polygon for polygon in polygons if polygon.area > 0.0])
    holes = [
        Polygon(ring)
        for polygon in _polygons(geometry)
        for ring in polygon.interiors
    ]
    projected_bounds = np.ptp(mesh.vertices[:, plane_axes], axis=0)
    minimum_area = (
        float(np.prod(projected_bounds)) * options.projection_hole_min_area_ratio
    )
    holes = [hole for hole in holes if float(hole.area) >= minimum_area]
    if not holes:
        return [], stats

    frame = np.column_stack(
        (
            np.eye(3)[:, plane_axes[0]],
            np.eye(3)[:, plane_axes[1]],
            np.eye(3)[:, thin_axis],
        )
    )
    tolerance = max(float(np.linalg.norm(projected_bounds)) * options.planar_relative_tolerance,
                    1.0e-10)
    descriptors: list[_FilledVoid] = []
    for hole in holes:
        hole = hole.simplify(tolerance, preserve_topology=True)
        targets: list[PrimitiveRegion] = []
        for plane_coordinate in (
            float(mesh.vertices[:, thin_axis].min()),
            float(mesh.vertices[:, thin_axis].max()),
        ):
            origin = np.zeros(3, dtype=np.float64)
            origin[thin_axis] = plane_coordinate
            patch = _PlanarPatch(
                origin,
                frame,
                hole,
                1.0,
                0.0,
                1,
                float(hole.area),
                0,
                0.0,
                hole,
                tolerance,
                hole,
            )
            for triangle in _triangulate_planar_patch(patch):
                target = PrimitiveRegion(
                    _triangle_primitive(triangle),
                    np.empty(0, dtype=np.int64),
                    np.empty(0, dtype=np.int64),
                    synthetic_fill=True,
                )
                regions.append(target)
                targets.append(target)
                stats["projection_fill_primitives"] += 1
        descriptors.append(_FilledVoid(np.zeros(3), frame, hole, targets, tolerance))
        stats["projection_filled_holes"] += 1
        stats["projection_filled_area"] += float(hole.area)
    return descriptors, stats


def _exclude_filled_void_surfaces(
    mesh: Mesh,
    regions: list[PrimitiveRegion],
    descriptors: Sequence[_FilledVoid],
) -> tuple[list[PrimitiveRegion], dict[str, int], IndexArray]:
    if not descriptors:
        return regions, {
            "excluded_void_surface_triangles": 0,
            "removed_void_surface_primitives": 0,
            "paired_void_descriptors": 0,
        }, np.empty(0, dtype=np.int64)
    triangles = mesh.vertices[mesh.faces]
    centroids = triangles.mean(axis=1)
    normals = np.cross(triangles[:, 1] - triangles[:, 0], triangles[:, 2] - triangles[:, 0])
    normals /= np.maximum(np.linalg.norm(normals, axis=1)[:, None], 1.0e-30)
    target_ids = {id(target) for descriptor in descriptors for target in descriptor.targets}
    claimed: set[int] = set()

    def geometry_in_frame(source: _FilledVoid, target: _FilledVoid) -> Any:
        polygons: list[Any] = []
        for polygon in _polygons(source.geometry):
            rings: list[FloatArray] = []
            for ring in (polygon.exterior, *polygon.interiors):
                coordinates = np.asarray(ring.coords, dtype=np.float64)
                local = np.column_stack((coordinates, np.zeros(len(coordinates))))
                world = source.origin + local @ source.frame.T
                rings.append(((world - target.origin) @ target.frame)[:, :2])
            polygons.append(Polygon(rings[0], rings[1:]))
        return unary_union(polygons)

    paired_depths: dict[int, tuple[float, float]] = {}
    for first_index, first in enumerate(descriptors):
        best: tuple[float, float] | None = None
        first_normal = first.frame[:, 2]
        for second_index, second in enumerate(descriptors):
            if first_index == second_index:
                continue
            if abs(float(np.dot(first_normal, second.frame[:, 2]))) < 0.9999:
                continue
            signed_depth = float(np.dot(second.origin - first.origin, first_normal))
            separation = abs(signed_depth)
            if separation <= max(first.tolerance, second.tolerance):
                continue
            transformed = geometry_in_frame(second, first)
            shared_area = float(first.geometry.intersection(transformed).area)
            minimum_area = min(float(first.geometry.area), float(transformed.area))
            if minimum_area <= 0.0 or shared_area / minimum_area < 0.95:
                continue
            if best is None or separation < best[0]:
                best = separation, signed_depth
        if best is not None:
            signed_depth = best[1]
            paired_depths[first_index] = (
                min(0.0, signed_depth), max(0.0, signed_depth)
            )

    for descriptor_index, descriptor in enumerate(descriptors):
        if not descriptor.targets:
            continue
        buffered = descriptor.geometry.buffer(descriptor.tolerance)
        paired_depth = paired_depths.get(descriptor_index)
        bounds = descriptor.geometry.bounds
        fallback_depth = max(
            math.hypot(bounds[2] - bounds[0], bounds[3] - bounds[1]),
            descriptor.tolerance,
        )
        for region in regions:
            if id(region) in target_ids or len(region.triangle_ids) == 0:
                continue
            keep: list[int] = []
            for raw_triangle_id in region.triangle_ids:
                triangle_id = int(raw_triangle_id)
                if triangle_id in claimed:
                    continue
                relative = centroids[triangle_id] - descriptor.origin
                local = relative @ descriptor.frame
                if paired_depth is not None:
                    within_depth = (
                        paired_depth[0] - descriptor.tolerance
                        <= float(local[2])
                        <= paired_depth[1] + descriptor.tolerance
                    )
                else:
                    within_depth = (
                        abs(float(local[2])) <= fallback_depth
                        and abs(float(np.dot(normals[triangle_id], descriptor.frame[:, 2]))) < 0.5
                    )
                if not within_depth:
                    keep.append(triangle_id)
                    continue
                if not buffered.covers(Point(float(local[0]), float(local[1]))):
                    keep.append(triangle_id)
                    continue
                claimed.add(triangle_id)
            region.triangle_ids = np.asarray(keep, dtype=np.int64)
            region.vertex_ids = (
                np.unique(mesh.faces[region.triangle_ids].reshape(-1))
                if len(region.triangle_ids) else np.empty(0, dtype=np.int64)
            )
    removed = sum(
        len(region.triangle_ids) == 0 and not region.synthetic_fill
        for region in regions
    )
    filtered = [
        region for region in regions
        if len(region.triangle_ids) or region.synthetic_fill
    ]
    return filtered, {
        "excluded_void_surface_triangles": len(claimed),
        "removed_void_surface_primitives": removed,
        "paired_void_descriptors": len(paired_depths),
    }, np.asarray(sorted(claimed), dtype=np.int64)


def _project_region_geometry(
    mesh: Mesh,
    region: PrimitiveRegion,
    origin: FloatArray,
    frame: FloatArray,
) -> Any:
    triangles = (mesh.vertices[mesh.faces[region.triangle_ids]] - origin) @ frame
    polygons = [Polygon(triangle[:, :2]) for triangle in triangles]
    return unary_union([polygon for polygon in polygons if polygon.area > 0.0])


def _exclude_covered_planar_regions(
    mesh: Mesh,
    regions: list[PrimitiveRegion],
    options: PrimitiveAnalysisOptions,
) -> tuple[list[PrimitiveRegion], dict[str, int], IndexArray]:
    planar: list[tuple[PrimitiveRegion, _PlanarPatch]] = []
    for region in regions:
        if region.primitive.kind == "triangle" or len(region.triangle_ids) == 0:
            continue
        patch = _planar_patch(mesh, region, options)
        if patch is not None:
            planar.append((region, patch))

    covered_ids: set[int] = set()
    excluded: list[int] = []
    for candidate, candidate_patch in planar:
        if id(candidate) in covered_ids:
            continue
        candidate_points = mesh.vertices[mesh.faces[candidate.triangle_ids]].reshape(-1, 3)
        candidate_area = float(candidate_patch.source_geometry.area)
        for coverer, coverer_patch in planar:
            if candidate is coverer or id(coverer) in covered_ids:
                continue
            normal_dot = abs(float(np.dot(
                candidate_patch.frame[:, 2], coverer_patch.frame[:, 2]
            )))
            if normal_dot < 0.9999:
                continue
            local_points = (candidate_points - coverer_patch.origin) @ coverer_patch.frame
            tolerance = max(candidate_patch.tolerance, coverer_patch.tolerance)
            if float(np.abs(local_points[:, 2]).max()) > tolerance:
                continue
            candidate_in_coverer = _project_region_geometry(
                mesh, candidate, coverer_patch.origin, coverer_patch.frame
            )
            coverer_area = float(coverer_patch.source_geometry.area)
            if coverer_area <= candidate_area * (1.0 + 1.0e-9):
                continue
            if coverer_patch.source_geometry.buffer(tolerance).covers(candidate_in_coverer):
                covered_ids.add(id(candidate))
                excluded.extend(int(value) for value in candidate.triangle_ids)
                break

    filtered = [region for region in regions if id(region) not in covered_ids]
    return filtered, {
        "excluded_covered_planar_triangles": len(excluded),
        "removed_covered_planar_primitives": len(covered_ids),
    }, np.asarray(sorted(excluded), dtype=np.int64)


def _replace_loose_planar_regions(
    mesh: Mesh,
    regions: Sequence[PrimitiveRegion],
    options: PrimitiveAnalysisOptions,
) -> tuple[Mesh, list[PrimitiveRegion], dict[str, int | float], IndexArray]:
    triangles_enabled = "triangle" in options.primitive_types
    output: list[PrimitiveRegion] = []
    converted_patches = 0
    triangle_primitives = 0
    filled_holes = 0
    filled_hole_area = 0.0
    filled_boundary_voids = 0
    filled_boundary_void_area = 0.0
    void_descriptors: list[_FilledVoid] = []

    def record_void(patch: _PlanarPatch, targets: list[PrimitiveRegion]) -> None:
        if not patch.filled_void_geometry.is_empty and patch.filled_void_geometry.area > 0.0:
            for component in _polygons(patch.filled_void_geometry):
                void_descriptors.append(
                    _FilledVoid(
                        patch.origin,
                        patch.frame,
                        component,
                        targets,
                        patch.tolerance,
                    )
                )

    for region in regions:
        patch = (
            _planar_patch(mesh, region, options)
            if len(region.triangle_ids) >= options.planar_triangle_min_source_faces
            else None
        )
        if patch is not None:
            filled_holes += patch.filled_hole_count
            filled_hole_area += patch.filled_hole_area
            filled_boundary_voids += patch.filled_boundary_void_count
            filled_boundary_void_area += patch.filled_boundary_void_area
        if (
            patch is None
            or not triangles_enabled
            or patch.proxy_coverage >= options.planar_triangle_max_coverage
        ):
            if patch is not None:
                region.primitive.planar_excess_area_ratio = patch.proxy_excess_ratio
                if patch.filled_hole_count:
                    region.outward_reference_triangles = np.asarray(
                        _triangulate_planar_patch(patch), dtype=np.float64
                    )
            output.append(region)
            if patch is not None:
                record_void(patch, [region])
            continue

        patch_triangles = _triangulate_planar_patch(patch)
        if len(patch_triangles) > len(region.triangle_ids):
            region.primitive.planar_excess_area_ratio = patch.proxy_excess_ratio
            if patch.filled_hole_count:
                region.outward_reference_triangles = np.asarray(
                    patch_triangles, dtype=np.float64
                )
            output.append(region)
            record_void(patch, [region])
            continue

        converted_patches += 1
        assignments = _assign_source_triangles(mesh, region, patch, patch_triangles)
        patch_targets: list[PrimitiveRegion] = []
        for triangle, triangle_ids in zip(patch_triangles, assignments, strict=True):
            vertex_ids = (
                np.unique(mesh.faces[triangle_ids].reshape(-1))
                if len(triangle_ids) else np.empty(0, dtype=np.int64)
            )
            target = PrimitiveRegion(
                _triangle_primitive(triangle),
                triangle_ids,
                vertex_ids,
                synthetic_fill=len(triangle_ids) == 0,
            )
            output.append(target)
            patch_targets.append(target)
            triangle_primitives += 1
        record_void(patch, patch_targets)

    projection_descriptors, projection_stats = _projection_hole_caps(mesh, output, options)
    void_descriptors.extend(projection_descriptors)
    output, wall_stats, wall_excluded = _exclude_filled_void_surfaces(
        mesh, output, void_descriptors
    )
    output, covered_stats, covered_excluded = _exclude_covered_planar_regions(
        mesh, output, options
    )
    excluded = np.union1d(wall_excluded, covered_excluded)
    triangle_primitives += int(projection_stats["projection_fill_primitives"])
    return mesh, output, {
        "triangle_patches": converted_patches,
        "triangle_primitives": triangle_primitives,
        "filled_planar_holes": filled_holes,
        "filled_planar_hole_area": filled_hole_area,
        "filled_boundary_voids": filled_boundary_voids,
        "filled_boundary_void_area": filled_boundary_void_area,
        **projection_stats,
        **wall_stats,
        **covered_stats,
        "excluded_redundant_triangles": len(excluded),
    }, excluded


def _validate_partition(
    face_count: int,
    nodes: Sequence[_Node],
    excluded_triangle_ids: IndexArray | None = None,
) -> None:
    assigned = np.concatenate([node.triangle_ids for node in nodes])
    excluded = (
        np.empty(0, dtype=np.int64)
        if excluded_triangle_ids is None else np.asarray(excluded_triangle_ids, dtype=np.int64)
    )
    combined = np.concatenate((assigned, excluded))
    if len(combined) != face_count:
        raise RuntimeError(
            f"triangle partition size mismatch: {len(assigned)} active + "
            f"{len(excluded)} excluded != {face_count}"
        )
    if not np.array_equal(np.sort(combined), np.arange(face_count, dtype=np.int64)):
        raise RuntimeError("triangle responsibility is not a strict partition")


def analyze_mesh(mesh: Mesh, options: PrimitiveAnalysisOptions) -> PrimitiveAnalysisResult:
    analysis_started = time.perf_counter()
    timings: dict[str, float] = {}
    options.validate()
    phase_started = time.perf_counter()
    mesh, original_triangle_count, discarded_degenerate = _drop_degenerate_faces(mesh)
    timings["degenerate_filter"] = time.perf_counter() - phase_started
    extent = mesh.vertices.max(axis=0) - mesh.vertices.min(axis=0)
    diagonal = float(np.linalg.norm(extent))
    margin = max(diagonal * options.relative_clearance, 1.0e-9)
    model_volume = max(float(np.prod(np.maximum(extent, 0.0))), diagonal ** 3 * 1.0e-12, 1.0e-30)
    threshold_ratio = options.excess_volume_ratio()

    phase_started = time.perf_counter()
    face_quadrics = _face_quadrics(mesh, options.tangent_quadric_weight)
    timings["face_quadrics"] = time.perf_counter() - phase_started
    phase_started = time.perf_counter()
    face_adjacency = _face_adjacency(mesh)
    timings["face_adjacency"] = time.perf_counter() - phase_started
    phase_started = time.perf_counter()
    clusters, adjacency = _coplanar_face_clusters(mesh, face_adjacency)
    timings["coplanar_clustering"] = time.perf_counter() - phase_started
    phase_started = time.perf_counter()
    nodes: dict[int, _Node] = {}
    for cluster_id, cluster in enumerate(clusters):
        quadric = face_quadrics[cluster.triangle_ids].sum(axis=0)
        nodes[cluster_id] = _make_node(
            cluster_id, mesh, cluster.triangle_ids, quadric,
            adjacency[cluster_id], margin, options,
        )
    timings["initial_primitive_fitting"] = time.perf_counter() - phase_started
    next_id = len(nodes)
    heap: list[tuple[float, float, int, int, Primitive]] = []
    candidate_fit_count = 0
    candidate_fit_seconds = 0.0
    candidate_fit_wall_seconds = 0.0
    heap_pop_count = 0
    stale_heap_pop_count = 0
    peak_heap_entries = 0

    def fit_pair(pair: tuple[int, int]) -> tuple[int, int, Primitive, float]:
        first_id, second_id = pair
        first = nodes[first_id]
        second = nodes[second_id]
        fit_started = time.perf_counter()
        merged = _merged_node(-1, mesh, first, second, margin, options)
        return first_id, second_id, merged.primitive, time.perf_counter() - fit_started

    def crosses_preserved_sharp_edge(first: _Node, second: _Node) -> bool:
        if first.primitive.kind != "rss" or second.primitive.kind != "rss":
            return False
        planar_radius_limit = margin * 1.01
        if (
            first.primitive.radius > planar_radius_limit
            or second.primitive.radius > planar_radius_limit
        ):
            return False
        normal_dot = abs(float(np.dot(
            first.primitive.axes[:, 2], second.primitive.axes[:, 2]
        )))
        threshold = math.cos(math.radians(options.planar_sharp_edge_degrees))
        return normal_dot < threshold

    def push_pairs(pairs: Iterable[tuple[int, int]]) -> None:
        nonlocal candidate_fit_count, candidate_fit_seconds
        nonlocal candidate_fit_wall_seconds, peak_heap_entries
        normalized: list[tuple[int, int]] = []
        for first_id, second_id in pairs:
            if first_id == second_id:
                continue
            if first_id > second_id:
                first_id, second_id = second_id, first_id
            if (
                nodes[first_id].active
                and nodes[second_id].active
                and not crosses_preserved_sharp_edge(nodes[first_id], nodes[second_id])
            ):
                normalized.append((first_id, second_id))
        batch_started = time.perf_counter()
        fitted = map(fit_pair, normalized)
        for first_id, second_id, primitive, fit_seconds in fitted:
            candidate_fit_seconds += fit_seconds
            candidate_fit_count += 1
            first = nodes[first_id]
            second = nodes[second_id]
            added_volume = primitive.volume - first.primitive.volume - second.primitive.volume
            added_ratio = added_volume / model_volume
            heapq.heappush(
                heap,
                (added_ratio, primitive.volume, first_id, second_id, primitive),
            )
            peak_heap_entries = max(peak_heap_entries, len(heap))
        candidate_fit_wall_seconds += time.perf_counter() - batch_started

    phase_started = time.perf_counter()
    initial_pairs: list[tuple[int, int]] = []
    for node in nodes.values():
        for neighbor in node.neighbors:
            if node.node_id < neighbor:
                initial_pairs.append((node.node_id, neighbor))
    push_pairs(initial_pairs)
    timings["initial_merge_candidates"] = time.perf_counter() - phase_started

    accepted_merges = 0
    next_rejected: float | None = None
    phase_started = time.perf_counter()
    while heap:
        while heap:
            added_ratio, _, first_id, second_id, primitive = heapq.heappop(heap)
            heap_pop_count += 1
            first = nodes[first_id]
            second = nodes[second_id]
            if first.active and second.active and second_id in first.neighbors:
                break
            stale_heap_pop_count += 1
        else:
            break
        if added_ratio > threshold_ratio:
            next_rejected = added_ratio
            break

        merged = _Node(
            next_id,
            np.concatenate((first.triangle_ids, second.triangle_ids)),
            np.union1d(first.vertex_ids, second.vertex_ids),
            first.quadric + second.quadric,
            primitive,
            set(),
        )
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
        push_pairs((merged.node_id, neighbor_id) for neighbor_id in merged.neighbors)
        accepted_merges += 1
    timings["merge_loop"] = time.perf_counter() - phase_started

    phase_started = time.perf_counter()
    active = sorted((node for node in nodes.values() if node.active),
                    key=lambda node: int(node.triangle_ids.min()))
    _validate_partition(len(mesh.faces), active)
    regions = [PrimitiveRegion(node.primitive, node.triangle_ids, node.vertex_ids) for node in active]
    timings["partition_validation"] = time.perf_counter() - phase_started
    phase_started = time.perf_counter()
    pre_retriangulation_triangles = len(mesh.faces)
    mesh, regions, planar_stats, excluded_triangle_ids = _replace_loose_planar_regions(
        mesh, regions, options
    )
    _validate_partition(len(mesh.faces), regions, excluded_triangle_ids)
    timings["planar_retriangulation"] = time.perf_counter() - phase_started
    phase_started = time.perf_counter()
    outward_metrics = _apply_outward_deviation_metrics(mesh, regions)
    timings["outward_deviation_metrics"] = time.perf_counter() - phase_started
    timings["candidate_primitive_fitting"] = candidate_fit_seconds
    timings["candidate_primitive_fitting_wall"] = candidate_fit_wall_seconds
    type_counts = {kind: 0 for kind in options.primitive_types}
    for region in regions:
        type_counts[region.primitive.kind] += 1
    type_counts = {kind: count for kind, count in type_counts.items() if count}
    stats = {
        "model": mesh.name,
        "source_vertices": len(mesh.vertices),
        "source_triangles": len(mesh.faces),
        "responsibility_triangles": sum(len(region.triangle_ids) for region in regions),
        "pre_retriangulation_triangles": pre_retriangulation_triangles,
        "original_source_triangles": original_triangle_count,
        "discarded_degenerate_triangles": discarded_degenerate,
        "initial_face_clusters": len(clusters),
        "primitive_count": len(regions),
        "primitive_types": type_counts,
        "accepted_merges": accepted_merges,
        "candidate_primitive_fits": candidate_fit_count,
        "heap_pops": heap_pop_count,
        "stale_heap_pops": stale_heap_pop_count,
        "peak_heap_entries": peak_heap_entries,
        "analysis_strength": options.analysis_strength,
        "max_excess_volume_ratio": threshold_ratio if math.isfinite(threshold_ratio) else None,
        "next_rejected_excess_volume_ratio": next_rejected,
        "total_primitive_volume": sum(region.primitive.volume for region in regions),
        "strict_triangle_partition": True,
        "conservative_coverage_by_construction": True,
        "collision_export_supported": (
            type_counts.get("triangle", 0) == 0 and len(excluded_triangle_ids) == 0
        ),
        "synthetic_fill_primitives": sum(region.synthetic_fill for region in regions),
        "timings_seconds": timings,
        **planar_stats,
        **outward_metrics,
    }
    timings["total"] = time.perf_counter() - analysis_started
    return PrimitiveAnalysisResult(mesh, regions, stats, excluded_triangle_ids)


def primitive_visualization_mesh(primitive: Primitive, subdivisions: int = 1) -> Mesh:
    if primitive.kind == "triangle":
        if primitive.triangle_vertices is None:
            raise RuntimeError("triangle primitive has no vertices")
        return Mesh(
            "triangle",
            primitive.triangle_vertices.copy(),
            np.asarray(((0, 1, 2),), dtype=np.int64),
        )
    directions = trimesh.creation.icosphere(subdivisions=subdivisions, radius=1.0).vertices
    core_corners = []
    for x in (0.0, float(primitive.lengths[0])):
        for y in (0.0, float(primitive.lengths[1])):
            core_corners.append(
                primitive.origin + primitive.axes[:, 0] * x + primitive.axes[:, 1] * y
            )
    core_corners = np.unique(np.asarray(core_corners), axis=0)
    points = np.vstack([corner + primitive.radius * directions for corner in core_corners])
    hull = trimesh.Trimesh(vertices=points, process=False).convex_hull
    return Mesh(
        primitive.kind,
        np.asarray(hull.vertices, dtype=np.float64),
        np.asarray(hull.faces, dtype=np.int64),
    )


def _write_grouped_obj(path: Path, meshes: Sequence[tuple[str, Mesh]]) -> None:
    with path.open("w", encoding="ascii", newline="\n") as stream:
        stream.write("# PQSSProxyMesh primitive analysis visualization\n")
        vertex_offset = 0
        for name, mesh in meshes:
            stream.write(f"g {name}\n")
            for vertex in mesh.vertices:
                stream.write(f"v {vertex[0]:.17g} {vertex[1]:.17g} {vertex[2]:.17g}\n")
            for face in mesh.faces:
                stream.write(
                    f"f {int(face[0]) + 1 + vertex_offset} "
                    f"{int(face[1]) + 1 + vertex_offset} "
                    f"{int(face[2]) + 1 + vertex_offset}\n"
                )
            vertex_offset += len(mesh.vertices)


def _write_region_obj(path: Path, result: PrimitiveAnalysisResult) -> None:
    with path.open("w", encoding="ascii", newline="\n") as stream:
        stream.write("# Source triangle partition by primitive responsibility\n")
        for vertex in result.mesh.vertices:
            stream.write(f"v {vertex[0]:.17g} {vertex[1]:.17g} {vertex[2]:.17g}\n")
        for region_id, region in enumerate(result.regions):
            stream.write(f"g region_{region_id:05d}_{region.primitive.kind}\n")
            for triangle_id in region.triangle_ids:
                face = result.mesh.faces[int(triangle_id)]
                stream.write(f"f {int(face[0]) + 1} {int(face[1]) + 1} {int(face[2]) + 1}\n")
        if len(result.excluded_triangle_ids):
            stream.write("g excluded_redundant_surface\n")
            for triangle_id in result.excluded_triangle_ids:
                face = result.mesh.faces[int(triangle_id)]
                stream.write(f"f {int(face[0]) + 1} {int(face[1]) + 1} {int(face[2]) + 1}\n")


def write_analysis_model(
    source_path: Path,
    model_directory: Path,
    result: PrimitiveAnalysisResult,
    options: PrimitiveAnalysisOptions,
) -> dict[str, Any]:
    model_directory.mkdir(parents=True, exist_ok=True)
    write_obj(
        model_directory / "source.obj",
        result.mesh,
        comments=[f"Source {source_path.name}; degenerate triangles removed"],
    )
    visual_meshes: list[tuple[str, Mesh]] = []
    primitive_records: list[dict[str, Any]] = []
    for primitive_id, region in enumerate(result.regions):
        visual = primitive_visualization_mesh(
            region.primitive, options.visualization_sphere_subdivisions
        )
        visual_meshes.append((f"primitive_{primitive_id:05d}_{region.primitive.kind}", visual))
        record = region.primitive.to_dict()
        record.update({
            "id": primitive_id,
            "triangle_count": len(region.triangle_ids),
            "triangle_ids": region.triangle_ids.tolist(),
            "synthetic_hole_fill": region.synthetic_fill,
        })
        primitive_records.append(record)
    _write_grouped_obj(model_directory / "primitives.obj", visual_meshes)
    _write_region_obj(model_directory / "regions.obj", result)
    metadata = {
        "stats": result.stats,
        "source": "source.obj",
        "regions": "regions.obj",
        "primitives": "primitives.obj",
        "primitive_records": primitive_records,
        "excluded_triangle_ids": result.excluded_triangle_ids.tolist(),
    }
    (model_directory / "model.json").write_text(
        json.dumps(metadata, ensure_ascii=True, separators=(",", ":")) + "\n",
        encoding="ascii",
    )
    return metadata


def _analyze_model_task(
    path: Path,
    model_directory: Path,
    options: PrimitiveAnalysisOptions,
) -> dict[str, Any]:
    model_started = time.perf_counter()
    phase_started = time.perf_counter()
    mesh = read_obj(path)
    read_seconds = time.perf_counter() - phase_started
    phase_started = time.perf_counter()
    result = analyze_mesh(mesh, options)
    mesh_analysis_seconds = time.perf_counter() - phase_started
    phase_started = time.perf_counter()
    write_analysis_model(path, model_directory, result, options)
    visualization_export_seconds = time.perf_counter() - phase_started
    total_seconds = time.perf_counter() - model_started
    return {
        "id": int(path.stem),
        "metadata": f"models/{path.stem}/model.json",
        "stats": result.stats,
        "analysis_seconds": total_seconds,
        "input_read_seconds": read_seconds,
        "mesh_analysis_seconds": mesh_analysis_seconds,
        "visualization_export_seconds": visualization_export_seconds,
    }


def analyze_directory(
    input_directory: Path,
    output_directory: Path,
    options: PrimitiveAnalysisOptions,
    model_ids: Iterable[int] | None = None,
    resume: bool = False,
) -> dict[str, Any]:
    options.validate()
    selected = None if model_ids is None else frozenset(model_ids)
    paths = sorted(input_directory.glob("*.obj"), key=lambda path: int(path.stem))
    if selected is not None:
        paths = [path for path in paths if int(path.stem) in selected]
    if not paths:
        raise ValueError("no selected numeric OBJ models found")
    if output_directory.exists() and any(output_directory.iterdir()) and not resume:
        raise FileExistsError(f"output directory is not empty: {output_directory}")
    output_directory.mkdir(parents=True, exist_ok=True)

    started = time.perf_counter()
    models: list[dict[str, Any]] = []
    worker_count = min(options.model_workers or max(1, os.cpu_count() or 1), len(paths))

    def write_manifest(complete: bool) -> dict[str, Any]:
        effective_ratio = options.excess_volume_ratio()
        ordered_models = sorted(models, key=lambda model: int(model["id"]))
        model_stats = [model["stats"] for model in ordered_models]
        elapsed_seconds = time.perf_counter() - started
        total_surface_area = sum(
            float(stats.get("proxy_sampled_surface_area", 0.0)) for stats in model_stats
        )
        mean_outward_deviation = (
            sum(
                float(stats.get("mean_sampled_outward_deviation", 0.0))
                * float(stats.get("proxy_sampled_surface_area", 0.0))
                for stats in model_stats
            )
            / total_surface_area
            if total_surface_area > 0.0
            else 0.0
        )
        summary = {
            "original_source_triangles": sum(
                int(stats.get("original_source_triangles", stats["source_triangles"]))
                for stats in model_stats
            ),
            "active_source_triangles": sum(
                int(stats.get("responsibility_triangles", stats["source_triangles"]))
                for stats in model_stats
            ),
            "discarded_degenerate_triangles": sum(
                int(stats.get("discarded_degenerate_triangles", 0)) for stats in model_stats
            ),
            "primitive_count": sum(int(stats["primitive_count"]) for stats in model_stats),
            "mean_sampled_outward_deviation": mean_outward_deviation,
            "maximum_sampled_outward_deviation": max(
                (float(stats.get("maximum_sampled_outward_deviation", 0.0)) for stats in model_stats),
                default=0.0,
            ),
            "maximum_planar_excess_area_ratio": max(
                (float(stats.get("maximum_planar_excess_area_ratio", 0.0)) for stats in model_stats),
                default=0.0,
            ),
            "filled_planar_holes": sum(
                int(stats.get("filled_planar_holes", 0)) for stats in model_stats
            ),
            "filled_planar_hole_area": sum(
                float(stats.get("filled_planar_hole_area", 0.0)) for stats in model_stats
            ),
            "synthetic_fill_primitives": sum(
                int(stats.get("synthetic_fill_primitives", 0)) for stats in model_stats
            ),
            "filled_boundary_voids": sum(
                int(stats.get("filled_boundary_voids", 0)) for stats in model_stats
            ),
            "filled_boundary_void_area": sum(
                float(stats.get("filled_boundary_void_area", 0.0)) for stats in model_stats
            ),
            "projection_filled_holes": sum(
                int(stats.get("projection_filled_holes", 0)) for stats in model_stats
            ),
            "projection_filled_area": sum(
                float(stats.get("projection_filled_area", 0.0)) for stats in model_stats
            ),
            "excluded_void_surface_triangles": sum(
                int(stats.get("excluded_void_surface_triangles", 0)) for stats in model_stats
            ),
            "excluded_redundant_triangles": sum(
                int(stats.get("excluded_redundant_triangles", 0)) for stats in model_stats
            ),
            "removed_covered_planar_primitives": sum(
                int(stats.get("removed_covered_planar_primitives", 0)) for stats in model_stats
            ),
            "removed_void_surface_primitives": sum(
                int(stats.get("removed_void_surface_primitives", 0)) for stats in model_stats
            ),
            "paired_void_descriptors": sum(
                int(stats.get("paired_void_descriptors", 0)) for stats in model_stats
            ),
            "model_analysis_seconds": sum(
                float(stats.get("timings_seconds", {}).get("total", 0.0)) for stats in model_stats
            ),
            "wall_seconds_per_completed_model": elapsed_seconds / max(len(models), 1),
        }
        manifest = {
            "algorithm": "HierarchicalConservativePrimitiveAnalysis",
            "input_directory": str(input_directory.resolve()),
            "output_directory": str(output_directory.resolve()),
            "options": asdict(options),
            "effective_max_excess_volume_ratio": (
                effective_ratio if math.isfinite(effective_ratio) else None
            ),
            "requested_model_count": len(paths),
            "model_count": len(models),
            "complete": complete,
            "model_workers": worker_count,
            "analysis_seconds": elapsed_seconds,
            "summary": summary,
            "models": ordered_models,
        }
        (output_directory / "viewer_manifest.json").write_text(
            json.dumps(manifest, indent=2, ensure_ascii=True) + "\n", encoding="ascii"
        )
        return manifest

    pending: list[Path] = []
    for path in paths:
        model_directory = output_directory / "models" / path.stem
        metadata_path = model_directory / "model.json"
        if resume and metadata_path.is_file():
            metadata = json.loads(metadata_path.read_text(encoding="ascii"))
            models.append({
                "id": int(path.stem),
                "metadata": f"models/{path.stem}/model.json",
                "stats": metadata["stats"],
                "analysis_seconds": None,
                "resumed": True,
            })
            print(f"reusing completed model {path.stem}", flush=True)
            write_manifest(False)
            continue
        pending.append(path)

    def record_completed(record: dict[str, Any]) -> None:
        models.append(record)
        stats = record["stats"]
        print(
            f"model={record['id']} primitives={stats['primitive_count']} "
            f"types={stats['primitive_types']} "
            f"discarded_degenerate={stats['discarded_degenerate_triangles']} "
            f"mean_outward_deviation={stats['mean_sampled_outward_deviation']:.9g} "
            f"maximum_outward_deviation={stats['maximum_sampled_outward_deviation']:.9g} "
            f"analysis_seconds={record['mesh_analysis_seconds']:.6f} "
            f"total_seconds={record['analysis_seconds']:.6f}",
            flush=True,
        )
        write_manifest(False)

    if pending:
        active_workers = min(worker_count, len(pending))
        print(f"analyzing {len(pending)} models with {active_workers} worker processes...", flush=True)
        if active_workers == 1:
            for path in pending:
                record_completed(
                    _analyze_model_task(path, output_directory / "models" / path.stem, options)
                )
        else:
            with ProcessPoolExecutor(max_workers=active_workers) as executor:
                futures = {
                    executor.submit(
                        _analyze_model_task,
                        path,
                        output_directory / "models" / path.stem,
                        options,
                    ): path
                    for path in pending
                }
                for future in as_completed(futures):
                    path = futures[future]
                    try:
                        record_completed(future.result())
                    except Exception as error:
                        raise RuntimeError(f"model {path.stem} analysis failed") from error

    return write_manifest(True)
