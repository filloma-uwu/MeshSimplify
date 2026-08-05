from __future__ import annotations

from dataclasses import asdict, dataclass
import math

import numpy as np
import numpy.typing as npt

from .obj_io import FloatArray, IndexArray, Mesh


@dataclass(slots=True)
class StaticBvhMetrics:
    triangles: int
    nodes: int
    internal_nodes: int
    maximum_depth: int
    sah_sum: float
    depth_weighted_area_sum: float
    sibling_overlap_area_sum: float
    sibling_overlap_volume_sum: float

    def to_dict(self) -> dict[str, int | float]:
        return asdict(self)


def _surface_area(bounds_min: FloatArray, bounds_max: FloatArray) -> float:
    extent = np.maximum(bounds_max - bounds_min, 0.0)
    return float(2.0 * (extent[0] * extent[1] + extent[1] * extent[2] + extent[2] * extent[0]))


def _split(
    triangle_ids: IndexArray,
    tri_min: FloatArray,
    tri_max: FloatArray,
    centroids: FloatArray,
) -> tuple[IndexArray, IndexArray]:
    best: tuple[float, npt.NDArray[np.bool_]] | None = None
    bucket_count = 16
    local_centroids = centroids[triangle_ids]
    for axis in range(3):
        low = float(local_centroids[:, axis].min())
        high = float(local_centroids[:, axis].max())
        if high - low < 1.0e-7:
            continue
        buckets = np.minimum(
            ((local_centroids[:, axis] - low) * (bucket_count - 1) / (high - low)).astype(np.int64),
            bucket_count - 1,
        )
        for cut in range(bucket_count - 1):
            left_mask = buckets < cut + 1
            if not left_mask.any() or left_mask.all():
                continue
            left_ids = triangle_ids[left_mask]
            right_ids = triangle_ids[~left_mask]
            left_cost = len(left_ids) * _surface_area(
                tri_min[left_ids].min(axis=0), tri_max[left_ids].max(axis=0)
            )
            right_cost = len(right_ids) * _surface_area(
                tri_min[right_ids].min(axis=0), tri_max[right_ids].max(axis=0)
            )
            cost = left_cost + right_cost
            if best is None or cost < best[0]:
                best = cost, left_mask
    if best is None:
        order = np.argsort(local_centroids[:, 0], kind="stable")
        middle = len(order) // 2
        return triangle_ids[order[:middle]], triangle_ids[order[middle:]]
    return triangle_ids[best[1]], triangle_ids[~best[1]]


def measure_static_bvh(mesh: Mesh) -> StaticBvhMetrics:
    triangles = mesh.vertices[mesh.faces]
    tri_min = triangles.min(axis=1)
    tri_max = triangles.max(axis=1)
    centroids = triangles.mean(axis=1)
    root_ids = np.arange(len(triangles), dtype=np.int64)

    nodes = 0
    internal_nodes = 0
    maximum_depth = 0
    sah_sum = 0.0
    depth_weighted_area_sum = 0.0
    overlap_area_sum = 0.0
    overlap_volume_sum = 0.0
    pending: list[tuple[IndexArray, int]] = [(root_ids, 0)]
    while pending:
        triangle_ids, depth = pending.pop()
        nodes += 1
        maximum_depth = max(maximum_depth, depth)
        bounds_min = tri_min[triangle_ids].min(axis=0)
        bounds_max = tri_max[triangle_ids].max(axis=0)
        area = _surface_area(bounds_min, bounds_max)
        depth_weighted_area_sum += (depth + 1) * area
        if len(triangle_ids) == 1:
            continue
        internal_nodes += 1
        left_ids, right_ids = _split(triangle_ids, tri_min, tri_max, centroids)
        left_min = tri_min[left_ids].min(axis=0)
        left_max = tri_max[left_ids].max(axis=0)
        right_min = tri_min[right_ids].min(axis=0)
        right_max = tri_max[right_ids].max(axis=0)
        sah_sum += len(left_ids) * _surface_area(left_min, left_max)
        sah_sum += len(right_ids) * _surface_area(right_min, right_max)

        overlap_min = np.maximum(left_min, right_min)
        overlap_max = np.minimum(left_max, right_max)
        overlap_extent = np.maximum(overlap_max - overlap_min, 0.0)
        overlap_area_sum += _surface_area(np.zeros(3), overlap_extent)
        overlap_volume_sum += float(np.prod(overlap_extent))
        pending.append((right_ids, depth + 1))
        pending.append((left_ids, depth + 1))

    values = (sah_sum, depth_weighted_area_sum, overlap_area_sum, overlap_volume_sum)
    if not all(math.isfinite(value) for value in values):
        raise RuntimeError(f"{mesh.name}: non-finite static BVH metric")
    return StaticBvhMetrics(
        triangles=len(triangles),
        nodes=nodes,
        internal_nodes=internal_nodes,
        maximum_depth=maximum_depth,
        sah_sum=sah_sum,
        depth_weighted_area_sum=depth_weighted_area_sum,
        sibling_overlap_area_sum=overlap_area_sum,
        sibling_overlap_volume_sum=overlap_volume_sum,
    )
