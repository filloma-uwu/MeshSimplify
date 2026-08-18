# MeshSimplify

`MeshSimplify` generates a low-face-count outer approximation of one input OBJ
under a user-specified one-sided Hausdorff error bound.

## Problem

Given a source surface mesh \(M\) and a non-negative error threshold
\(\varepsilon\), produce an output OBJ \(P\) with as few triangles as possible:

\[
\min_P |F(P)|
\]

subject to two hard constraints:

1. **Conservative outer coverage**: the proxy must not cut into or lose source
   geometry. Every replacement must be certified as an outer/conservative
   replacement of its source responsibility.
2. **Directed error bound**:

\[
d_{\rightarrow}(P,M)
=
\sup_{p\in P}\inf_{x\in M}\|p-x\|
\le \varepsilon .
\]

The error threshold is a feasibility constraint, not a weighted objective.
Among feasible proxies, fewer output triangles are always preferred.

The project is geometry-only. PQSS query cost, BVH depth, collision poses,
triangle-test counts, and model-pool behavior are not optimization objectives.

## Current Pipeline

The useful implementation is the staged C++ analyzer in
`src/primitive_mesh_analyzer.cpp`.

- **Stage 0 — source normalization / conservative preparation**  
  Load the OBJ, remove numerical degeneracies, establish geometric
  responsibilities, and prepare the source representation used by later
  conservative certificates. This stage is part of the working pipeline.

- **Stage 1 — conservative outer preprocessing**  
  Apply the currently accepted conservative closure/fill operations and produce
  the stable outer reference used by later approximation-error checks. This
  stage is part of the working pipeline.

- **Stage 2 — surface recognition (exploratory)**  
  Recognize and merge certified surface patches such as planar polygons and
  round analytic surface families. This stage is under active exploration.

- **Stage 3 — Hausdorff-limited simplification (exploratory)**  
  Search for additional conservative surface replacements while respecting the
  user-provided directed error threshold. The long-term objective is minimum
  triangle count among all feasible outputs. This stage is under active
  exploration.

- **Triangulation / export**  
  Convert the selected semantic surface representation to the final OBJ and
  run final containment and directed-distance audits.

The implementation still writes historical intermediate filenames such as
`phase1_hole_filled.obj`, `phase2_recognized_surfaces.obj`, `primitives.obj`,
and `proxy.obj`. Those names are retained temporarily to avoid mixing repository
cleanup with an algorithm/output-format migration.

See [`docs/STAGES.md`](docs/STAGES.md) for the stage contract.


## Project State and Cross-Conversation Continuity

The persistent handoff document is
[`docs/PROJECT_STATE.md`](docs/PROJECT_STATE.md).

Any new conversation or contributor should read that file before changing the
algorithm. Important project decisions, blockers, rejected approaches, stage
status, and workflow changes must be written there rather than left only in
chat history.

## Build

Requirements:

- CMake 3.24+
- a C++20 compiler
- the Clipper2 submodule in `deps/clipper2`

```powershell
git submodule update --init --recursive
cmake -S . -B build
cmake --build build --config Release
ctest --test-dir build -C Release --output-on-failure
```

## Run

```powershell
.\build\Release\mesh-simplify-analyze.exe `
  --input path\to\model.obj `
  --output-dir outputs\model `
  --maximum-open-error-distance 1.0
```

`--maximum-open-error-distance` is expressed in the OBJ's model units.

The current executable still exposes several experimental tuning switches
because Stage 2 and Stage 3 are not finalized. New algorithm work should reduce
those ad-hoc controls rather than make them part of the public problem
definition.

## Visualization

Serve an analyzer output manifest with:

```powershell
python tools\serve_primitive_viewer.py `
  --manifest path\to\viewer_manifest.json
```

The maintained viewer is `viewer/primitive_analysis.html`.

## Repository Scope

This repository intentionally contains only:

- the Stage 0–3 C++ geometry pipeline;
- the final triangulation/validation path;
- model-independent tests and geometry investigation tools;
- the primitive-analysis viewer;
- source OBJ data used for geometric regression.

Historical PQSS/BVH benchmarking, independent primitive-BVH experiments,
official-simplified comparison scopes, and the older Python proxy generators
are not part of the current project.
