"""Offline enclosing topological normalization for arbitrary triangle soups."""

from __future__ import annotations

from dataclasses import asdict, dataclass
from fractions import Fraction
import hashlib
import json
import math
import os
from pathlib import Path
import re
import shutil
import subprocess
import time

import numpy as np
from PIL import Image
from scipy import ndimage as ndi
from skimage import measure
import trimesh


SHAPE_THRESHOLD = 127
DISTANCE_INTENSITY_PER_VOXEL = 8.0


@dataclass(frozen=True)
class Grid:
    occupancy: np.ndarray
    pitch: float
    origin: np.ndarray


@dataclass(frozen=True)
class BettiNumbers:
    beta0: int
    beta1: int
    beta2: int


@dataclass(frozen=True)
class SolverAudit:
    input_betti: BettiNumbers
    output_betti: BettiNumbers
    solver_input_topology: tuple[int, int, int]
    solver_output_topology: tuple[int, int, int]
    input_voxels: int
    output_voxels: int
    added_voxels: int
    removed_voxels: int
    cut_candidates: int
    fill_candidates: int
    elapsed_seconds: float


@dataclass(frozen=True)
class EnclosingFillStep:
    axis: int
    radius: int
    added_voxels: int
    topology_before: BettiNumbers
    topology_after: BettiNumbers


@dataclass(frozen=True)
class EnclosingFillAudit:
    input_betti: BettiNumbers
    output_betti: BettiNumbers
    input_voxels: int
    output_voxels: int
    added_voxels: int
    removed_voxels: int
    steps: tuple[EnclosingFillStep, ...]


_TOPOLOGY_RE = re.compile(
    r"(?:Original Shape Topology|Final topology): Components: (\d+) "
    r"Cavities: (\d+) (?:Cycles|Handles): (\d+)"
)
_CANDIDATE_RE = re.compile(r"cuts before: (\d+) num fills before (\d+)")


def choose_pitch(extents: np.ndarray, maximum_grid_voxels: int, padding: int) -> float:
    """Choose one isotropic pitch from a model-independent grid-cell limit."""
    extents = np.asarray(extents, dtype=np.float64)
    if extents.shape != (3,) or np.any(~np.isfinite(extents)) or np.any(extents <= 0):
        raise ValueError(f"invalid mesh extents: {extents}")
    if maximum_grid_voxels <= (2 * padding + 3) ** 3:
        raise ValueError("maximum_grid_voxels is too small for the requested padding")

    def cells(pitch: float) -> int:
        shape = np.ceil(extents / pitch).astype(np.int64) + 1 + 2 * padding
        return int(np.prod(shape, dtype=np.int64))

    low = np.finfo(np.float64).eps * float(np.max(extents))
    high = float(np.max(extents))
    while cells(high) > maximum_grid_voxels:
        high *= 2.0
    for _ in range(80):
        midpoint = (low + high) * 0.5
        if cells(midpoint) > maximum_grid_voxels:
            low = midpoint
        else:
            high = midpoint
    return high


def voxelize_triangle_soup(
    mesh: trimesh.Trimesh,
    maximum_grid_voxels: int,
    padding: int = 4,
) -> Grid:
    """Rasterize source triangles and classify solid cells by exterior flooding."""
    pitch = choose_pitch(mesh.extents, maximum_grid_voxels, padding)
    shape = np.ceil(mesh.extents / pitch).astype(np.int64) + 1 + 2 * padding
    origin = np.asarray(mesh.bounds[0], dtype=np.float64) - padding * pitch
    surface = rasterize_triangle_surface(
        np.asarray(mesh.vertices, dtype=np.float64),
        np.asarray(mesh.faces, dtype=np.int64),
        pitch,
        origin,
        tuple(int(value) for value in shape),
    )
    empty = ~surface
    seeds = np.zeros_like(empty)
    seeds[0, :, :] = empty[0, :, :]
    seeds[-1, :, :] = empty[-1, :, :]
    seeds[:, 0, :] = empty[:, 0, :]
    seeds[:, -1, :] = empty[:, -1, :]
    seeds[:, :, 0] = empty[:, :, 0]
    seeds[:, :, -1] = empty[:, :, -1]
    outside = ndi.binary_propagation(
        seeds,
        structure=ndi.generate_binary_structure(3, 1),
        mask=empty,
    )
    occupancy = ~outside
    return Grid(occupancy=occupancy, pitch=pitch, origin=origin)


def rasterize_triangle_surface(
    vertices: np.ndarray,
    faces: np.ndarray,
    pitch: float,
    origin: np.ndarray,
    shape: tuple[int, int, int],
    chunk_cells: int = 65_536,
) -> np.ndarray:
    """Conservatively mark every grid cube intersected by a source triangle.

    Each triangle is projected along its dominant normal. A 2D separating-axis
    test finds intersecting grid squares, then the triangle plane's full height
    range over each square marks the corresponding cubes. False-positive cells
    are possible at triangle edges, but an intersected source cube is not missed.
    """
    grid = np.zeros(shape, dtype=bool)
    grid_vertices = (np.asarray(vertices, dtype=np.float64) - origin) / pitch
    epsilon = 64.0 * np.finfo(np.float64).eps
    for face in np.asarray(faces, dtype=np.int64):
        triangle = grid_vertices[face]
        normal = np.cross(triangle[1] - triangle[0], triangle[2] - triangle[0])
        dominant = int(np.argmax(np.abs(normal)))
        if abs(normal[dominant]) <= epsilon:
            continue
        projected_axes = [axis for axis in range(3) if axis != dominant]
        first_axis, second_axis = projected_axes
        projected = triangle[:, projected_axes]
        lower = np.maximum(
            np.ceil(np.min(projected, axis=0) - 0.5 - epsilon).astype(np.int64),
            0,
        )
        upper = np.minimum(
            np.floor(np.max(projected, axis=0) + 0.5 + epsilon).astype(np.int64),
            np.asarray([shape[first_axis] - 1, shape[second_axis] - 1]),
        )
        if np.any(lower > upper):
            continue

        first_count = int(upper[0] - lower[0] + 1)
        second_count = int(upper[1] - lower[1] + 1)
        candidate_count = first_count * second_count
        edges = np.roll(projected, -1, axis=0) - projected
        separating_axes = np.column_stack((edges[:, 1], -edges[:, 0]))
        triangle_intervals = projected @ separating_axes.T
        interval_lower = np.min(triangle_intervals, axis=0)
        interval_upper = np.max(triangle_intervals, axis=0)
        square_radius = 0.5 * np.sum(np.abs(separating_axes), axis=1)

        plane_slopes = np.asarray([
            -normal[first_axis] / normal[dominant],
            -normal[second_axis] / normal[dominant],
        ])
        height_radius = 0.5 * float(np.sum(np.abs(plane_slopes)))
        for start in range(0, candidate_count, chunk_cells):
            flat = np.arange(start, min(start + chunk_cells, candidate_count))
            first = lower[0] + flat // second_count
            second = lower[1] + flat % second_count
            centers = np.column_stack((first, second)).astype(np.float64)
            center_intervals = centers @ separating_axes.T
            intersects = np.all(
                (center_intervals + square_radius >= interval_lower - epsilon)
                & (center_intervals - square_radius <= interval_upper + epsilon),
                axis=1,
            )
            if not np.any(intersects):
                continue
            first = first[intersects]
            second = second[intersects]
            centers = centers[intersects]
            center_height = triangle[0, dominant] + (
                centers - triangle[0, projected_axes]
            ) @ plane_slopes
            depth_lower = np.ceil(center_height - height_radius - 0.5 - epsilon).astype(np.int64)
            depth_upper = np.floor(center_height + height_radius + 0.5 + epsilon).astype(np.int64)
            depth_lower = np.maximum(depth_lower, 0)
            depth_upper = np.minimum(depth_upper, shape[dominant] - 1)
            maximum_span = int(np.max(depth_upper - depth_lower, initial=-1))
            for offset in range(maximum_span + 1):
                active = depth_lower + offset <= depth_upper
                indices = [None, None, None]
                indices[first_axis] = first[active]
                indices[second_axis] = second[active]
                indices[dominant] = depth_lower[active] + offset
                grid[tuple(indices)] = True
    return grid


def signed_distance_image(
    occupancy: np.ndarray,
    intensity_per_voxel: float = DISTANCE_INTENSITY_PER_VOXEL,
) -> np.ndarray:
    """Encode a signed voxel distance field while preserving the exact threshold set."""
    occupancy = np.asarray(occupancy, dtype=bool)
    inside = ndi.distance_transform_edt(occupancy)
    outside = ndi.distance_transform_edt(~occupancy)
    signed_distance = inside - outside
    encoded = np.rint(
        SHAPE_THRESHOLD + 0.5 + intensity_per_voxel * signed_distance
    )
    encoded = np.clip(encoded, 1, 254).astype(np.uint8)
    if not np.array_equal(encoded > SHAPE_THRESHOLD, occupancy):
        raise RuntimeError("distance-field encoding changed the input occupancy")
    return encoded


def write_png_volume(volume: np.ndarray, directory: Path) -> None:
    directory.mkdir(parents=True, exist_ok=False)
    for z in range(volume.shape[2]):
        Image.fromarray(volume[:, :, z].T).save(directory / f"{z:04d}.png")


def read_png_occupancy(directory: Path, depth: int) -> np.ndarray:
    paths = sorted(directory.glob("*.png"))
    if len(paths) != depth:
        raise RuntimeError(f"expected {depth} solver slices, found {len(paths)}")
    slices = [np.asarray(Image.open(path).convert("L"), dtype=np.uint8).T for path in paths]
    return np.stack(slices, axis=2) > SHAPE_THRESHOLD


def voxel_betti_numbers(occupancy: np.ndarray) -> BettiNumbers:
    occupancy = np.asarray(occupancy, dtype=bool)
    foreground_structure = ndi.generate_binary_structure(3, 1)
    _, beta0 = ndi.label(occupancy, structure=foreground_structure)
    background_labels, background_count = ndi.label(
        ~occupancy, structure=ndi.generate_binary_structure(3, 3)
    )
    boundary_labels = np.unique(np.concatenate([
        background_labels[0, :, :].ravel(),
        background_labels[-1, :, :].ravel(),
        background_labels[:, 0, :].ravel(),
        background_labels[:, -1, :].ravel(),
        background_labels[:, :, 0].ravel(),
        background_labels[:, :, -1].ravel(),
    ]))
    beta2 = background_count - int(np.count_nonzero(boundary_labels))
    euler = int(measure.euler_number(occupancy, connectivity=1))
    beta1 = int(beta0) + int(beta2) - euler
    return BettiNumbers(int(beta0), int(beta1), int(beta2))


def _fill_cavities(occupancy: np.ndarray) -> np.ndarray:
    """Fill bounded 26-connected background components without removing solids."""
    return ndi.binary_fill_holes(
        np.asarray(occupancy, dtype=bool),
        structure=ndi.generate_binary_structure(3, 3),
    )


def _axis_closing(occupancy: np.ndarray, axis: int, radius: int) -> np.ndarray:
    """Fill grid-line intervals bracketed by occupied cells along one axis."""
    if axis not in (0, 1, 2):
        raise ValueError(f"invalid axis: {axis}")
    if radius <= 0:
        raise ValueError(f"invalid closing radius: {radius}")
    pad = radius + 1
    pad_width = [(0, 0), (0, 0), (0, 0)]
    pad_width[axis] = (pad, pad)
    padded = np.pad(np.asarray(occupancy, dtype=bool), pad_width)
    structure_shape = [1, 1, 1]
    structure_shape[axis] = 3
    closed = ndi.binary_closing(
        padded,
        structure=np.ones(structure_shape, dtype=bool),
        iterations=radius,
        border_value=0,
    )
    crop = [slice(None), slice(None), slice(None)]
    crop[axis] = slice(pad, -pad)
    return np.asarray(occupancy, dtype=bool) | closed[tuple(crop)]


def _topology_defect(topology: BettiNumbers) -> int:
    return abs(topology.beta0 - 1) + topology.beta1 + topology.beta2


def enclosing_topology_fill(
    occupancy: np.ndarray,
    maximum_steps: int = 32,
) -> tuple[np.ndarray, EnclosingFillAudit]:
    """Add bracketed voxel intervals until the solid has topology (1, 0, 0).

    Every candidate is unioned with the current solid, so source occupancy is a
    hard kernel. Axes and scales compete in one search; no model identity or
    coordinate is used. Power-of-two scales bound the search cost, and the first
    improving radius in each scale interval is recovered by binary refinement.
    """
    source = np.asarray(occupancy, dtype=bool)
    if source.ndim != 3 or not np.any(source):
        raise ValueError("occupancy must be a non-empty 3D array")
    current = _fill_cavities(source)
    input_topology = voxel_betti_numbers(source)
    current_topology = voxel_betti_numbers(current)
    steps: list[EnclosingFillStep] = []

    for _ in range(maximum_steps):
        current_defect = _topology_defect(current_topology)
        if current_defect == 0:
            break
        candidates: list[
            tuple[Fraction, int, int, int, int, np.ndarray, BettiNumbers]
        ] = []
        for axis in range(3):
            maximum_radius = int(current.shape[axis])
            radii = []
            radius = 1
            while radius < maximum_radius:
                radii.append(radius)
                radius *= 2
            radii.append(maximum_radius)
            cache: dict[int, tuple[np.ndarray, BettiNumbers, int]] = {}

            def evaluate(candidate_radius: int) -> tuple[np.ndarray, BettiNumbers, int]:
                if candidate_radius not in cache:
                    candidate = _fill_cavities(
                        _axis_closing(current, axis, candidate_radius)
                    )
                    topology = voxel_betti_numbers(candidate)
                    added = int(np.count_nonzero(candidate & ~current))
                    cache[candidate_radius] = candidate, topology, added
                return cache[candidate_radius]

            previous_radius = 0
            previous_defect = current_defect
            candidate_radii = set(radii)
            for coarse_radius in radii:
                _, coarse_topology, _ = evaluate(coarse_radius)
                coarse_defect = _topology_defect(coarse_topology)
                if coarse_defect < current_defect and previous_defect >= current_defect:
                    low = previous_radius + 1
                    high = coarse_radius
                    while low < high:
                        midpoint = (low + high) // 2
                        _, midpoint_topology, _ = evaluate(midpoint)
                        if _topology_defect(midpoint_topology) < current_defect:
                            high = midpoint
                        else:
                            low = midpoint + 1
                    candidate_radii.add(low)
                previous_radius = coarse_radius
                previous_defect = coarse_defect

            for candidate_radius in sorted(candidate_radii):
                candidate, topology, added = evaluate(candidate_radius)
                defect = _topology_defect(topology)
                reduction = current_defect - defect
                if reduction <= 0 or added <= 0:
                    continue
                candidates.append((
                    Fraction(added, reduction),
                    defect,
                    added,
                    axis,
                    candidate_radius,
                    candidate,
                    topology,
                ))

        if not candidates:
            raise RuntimeError(
                f"enclosing topology fill stalled at {current_topology}"
            )
        candidates.sort(key=lambda item: item[:5])
        _, _, added, axis, radius, selected, selected_topology = candidates[0]
        steps.append(EnclosingFillStep(
            axis=axis,
            radius=radius,
            added_voxels=added,
            topology_before=current_topology,
            topology_after=selected_topology,
        ))
        current = selected
        current_topology = selected_topology
    else:
        raise RuntimeError(
            f"enclosing topology fill exceeded {maximum_steps} steps"
        )

    if current_topology != BettiNumbers(1, 0, 0):
        raise RuntimeError(f"topology target was not reached: {current_topology}")
    removed = int(np.count_nonzero(source & ~current))
    if removed != 0:
        raise RuntimeError("enclosing fill removed source occupancy")
    audit = EnclosingFillAudit(
        input_betti=input_topology,
        output_betti=current_topology,
        input_voxels=int(np.count_nonzero(source)),
        output_voxels=int(np.count_nonzero(current)),
        added_voxels=int(np.count_nonzero(current & ~source)),
        removed_voxels=removed,
        steps=tuple(steps),
    )
    return current, audit


def parse_solver_log(log: str) -> tuple[tuple[int, int, int], tuple[int, int, int], int, int]:
    topology = [tuple(int(value) for value in match) for match in _TOPOLOGY_RE.findall(log)]
    candidates = _CANDIDATE_RE.search(log)
    if len(topology) < 2 or candidates is None:
        raise RuntimeError("TopoSimplifier did not report complete topology statistics")
    return topology[0], topology[-1], int(candidates.group(1)), int(candidates.group(2))


def run_fill_only_solver(
    distance_field: np.ndarray,
    executable: Path,
    work_directory: Path,
    global_time: int = 30,
    local_time: int = 10,
) -> tuple[np.ndarray, SolverAudit, str]:
    """Run the paper implementation with finite K and implicit N=-infinity."""
    executable = executable.resolve()
    dictionary = executable.parent / "simpleDictionaryFull.bin"
    if not executable.is_file():
        raise FileNotFoundError(executable)
    if not dictionary.is_file() or dictionary.stat().st_size < 100_000_000:
        raise RuntimeError(f"full simpleDictionaryFull.bin is missing beside {executable}")

    input_directory = work_directory / "solver_input"
    output_directory = work_directory / "solver_output"
    write_png_volume(distance_field, input_directory)
    output_directory.mkdir(parents=True, exist_ok=False)
    command = [
        str(executable),
        "--in", str(input_directory) + os.sep,
        "--out", str(output_directory) + os.sep,
        "--numSlices", str(distance_field.shape[2]),
        "--S", str(SHAPE_THRESHOLD),
        # A finite K equal to S preserves every source component as kernel.
        # Omitting N activates the implementation's fillOnly branch.
        "--K", str(SHAPE_THRESHOLD),
        "--shapeTopo",
        "--finalTopo", "1",
        "--showGeomCost",
        "--globalTime", str(global_time),
        "--localTime", str(local_time),
        "--greedy", "1",
    ]
    started = time.perf_counter()
    result = subprocess.run(
        command,
        cwd=executable.parent,
        capture_output=True,
        text=True,
        errors="replace",
        check=False,
    )
    elapsed = time.perf_counter() - started
    log = result.stdout + result.stderr
    (solver_input, solver_output, cut_candidates, fill_candidates) = parse_solver_log(log)
    if result.returncode != 0:
        raise RuntimeError(f"TopoSimplifier failed with code {result.returncode}\n{log}")

    input_occupancy = distance_field > SHAPE_THRESHOLD
    output_occupancy = read_png_occupancy(output_directory, distance_field.shape[2])
    removed = int(np.count_nonzero(input_occupancy & ~output_occupancy))
    added = int(np.count_nonzero(output_occupancy & ~input_occupancy))
    input_betti = voxel_betti_numbers(input_occupancy)
    output_betti = voxel_betti_numbers(output_occupancy)
    audit = SolverAudit(
        input_betti=input_betti,
        output_betti=output_betti,
        solver_input_topology=solver_input,
        solver_output_topology=solver_output,
        input_voxels=int(np.count_nonzero(input_occupancy)),
        output_voxels=int(np.count_nonzero(output_occupancy)),
        added_voxels=added,
        removed_voxels=removed,
        cut_candidates=cut_candidates,
        fill_candidates=fill_candidates,
        elapsed_seconds=elapsed,
    )
    if cut_candidates != 0 or removed != 0:
        raise RuntimeError(f"fill-only certificate failed: {audit}")
    if output_betti != BettiNumbers(1, 0, 0) or solver_output != (1, 0, 0):
        raise RuntimeError(f"topology target was not reached: {audit}")
    if np.any(output_occupancy[[0, -1], :, :]) or np.any(output_occupancy[:, [0, -1], :]) or np.any(output_occupancy[:, :, [0, -1]]):
        raise RuntimeError("filled occupancy reached the padded grid boundary")
    return output_occupancy, audit, log


def occupancy_mesh(grid: Grid, occupancy: np.ndarray) -> trimesh.Trimesh:
    vertices, faces, _, _ = measure.marching_cubes(
        np.asarray(occupancy, dtype=np.float32),
        level=0.5,
        spacing=(grid.pitch, grid.pitch, grid.pitch),
        allow_degenerate=False,
        method="lewiner",
    )
    vertices += grid.origin
    mesh = trimesh.Trimesh(vertices=vertices, faces=faces, process=True)
    mesh.fix_normals()
    return mesh


def marching_cubes_enclosure_guard(occupancy: np.ndarray) -> np.ndarray:
    """Add one Chebyshev layer so the 0.5 isosurface encloses source voxel cubes."""
    occupancy = np.asarray(occupancy, dtype=bool)
    if occupancy.ndim != 3:
        raise ValueError("occupancy must be a 3D array")
    return ndi.binary_dilation(occupancy, structure=np.ones((3, 3, 3), dtype=bool))


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def generate_filled_model(
    source: Path,
    output_directory: Path,
    maximum_grid_voxels: int,
    padding: int,
) -> dict[str, object]:
    started = time.perf_counter()
    loaded = trimesh.load(source, force="mesh", process=False)
    if not isinstance(loaded, trimesh.Trimesh) or len(loaded.faces) == 0:
        raise RuntimeError(f"no triangle mesh in {source}")
    grid = voxelize_triangle_soup(loaded, maximum_grid_voxels, padding)
    fill_started = time.perf_counter()
    filled_occupancy, audit = enclosing_topology_fill(grid.occupancy)
    fill_elapsed = time.perf_counter() - fill_started
    meshed_occupancy = marching_cubes_enclosure_guard(filled_occupancy)
    meshed_topology = voxel_betti_numbers(meshed_occupancy)
    if meshed_topology != BettiNumbers(1, 0, 0):
        raise RuntimeError(f"enclosure guard changed target topology: {meshed_topology}")
    if (np.any(meshed_occupancy[[0, -1], :, :]) or
            np.any(meshed_occupancy[:, [0, -1], :]) or
            np.any(meshed_occupancy[:, :, [0, -1]])):
        raise RuntimeError("enclosure guard reached the padded grid boundary")
    filled_mesh = occupancy_mesh(grid, meshed_occupancy)
    if not filled_mesh.is_watertight or not filled_mesh.is_winding_consistent:
        raise RuntimeError("marching-cubes boundary is not watertight and consistently wound")

    shutil.copyfile(source, output_directory / "source.obj")
    filled_mesh.export(output_directory / "phase1_hole_filled.obj", file_type="obj")
    elapsed = time.perf_counter() - started
    stats = {
        "model": source.name,
        "source_triangles": int(len(loaded.faces)),
        "primitive_count": 0,
        "primitive_types": {},
        "proxy_triangles": int(len(filled_mesh.faces)),
        "timings_seconds": {"total": elapsed, "topology_fill": fill_elapsed},
        "topology_fill": {
            **asdict(audit),
            "pitch": grid.pitch,
            "grid_shape": list(grid.occupancy.shape),
            "maximum_grid_voxels": maximum_grid_voxels,
            "padding": padding,
            "watertight": bool(filled_mesh.is_watertight),
            "winding_consistent": bool(filled_mesh.is_winding_consistent),
            "mesh_components": int(len(filled_mesh.split(only_watertight=False))),
            "marching_cubes_guard_layers": 1,
            "meshed_occupancy_voxels": int(np.count_nonzero(meshed_occupancy)),
            "meshed_occupancy_betti": asdict(meshed_topology),
            "source_voxel_cubes_enclosed": bool(np.all(meshed_occupancy[grid.occupancy])),
        },
    }
    metadata = {
        "stats": stats,
        "source": "source.obj",
        "phase1_hole_filled": "phase1_hole_filled.obj",
        "phase2_recognized_surfaces": "phase1_hole_filled.obj",
        "phase3_simplified_surfaces": "phase1_hole_filled.obj",
        "phase4_triangulated": "phase1_hole_filled.obj",
        "proxy": "phase1_hole_filled.obj",
        "proxy_components": [],
        "viewer_stages": ["source", "phase1", "split"],
    }
    (output_directory / "model.json").write_text(
        json.dumps(metadata, indent=2) + "\n", encoding="ascii"
    )
    return metadata


def generate_batch(
    sources: list[Path],
    output_root: Path,
    maximum_grid_voxels: int = 3_000_000,
    padding: int = 4,
) -> Path:
    if output_root.exists():
        raise FileExistsError(output_root)
    output_root.mkdir(parents=True)
    models = []
    for source in sources:
        model_id = int(source.stem)
        model_directory = output_root / "models" / str(model_id)
        model_directory.mkdir(parents=True)
        generate_filled_model(
            source,
            model_directory,
            maximum_grid_voxels,
            padding,
        )
        models.append({"id": model_id, "metadata": f"models/{model_id}/model.json"})

    manifest = {
        "algorithm": "EnclosingAxialTopologyFill",
        "topology_reference_doi": "10.1145/3414685.3417854",
        "complete": True,
        "model_count": len(models),
        "models": models,
        "options": {
            "maximum_grid_voxels": maximum_grid_voxels,
            "padding": padding,
            "candidate_axes": [0, 1, 2],
            "scale_schedule": "powers of two with first-improvement binary refinement",
            "selection": "minimum added voxels per reduced topology defect",
            "target_betti": [1, 0, 0],
            "source_occupancy_is_hard_kernel": True,
            "marching_cubes_enclosure_guard_layers": 1,
        },
    }
    manifest_path = output_root / "viewer_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="ascii")
    return manifest_path
