# Project State

This file is the cross-conversation source of truth for the current state of
`MeshSimplify`.

A new ChatGPT conversation, developer, or contributor should read this file
before changing the algorithm.

## 1. Core Objective

Input:

- one OBJ triangle mesh `M`;
- one non-negative user-selected error limit `epsilon`.

Output:

- one OBJ proxy `P`.

Optimization problem:

\[
\min_P |F(P)|
\]

subject to:

1. `P` is a conservative outer replacement of the source geometry;
2. the directed proxy-to-source Hausdorff distance obeys

\[
d_{\rightarrow}(P,M)
=
\sup_{p\in P}\inf_{x\in M}\|p-x\|
\le \epsilon.
\]

The directed-error threshold is a hard feasibility constraint, not a weighted
objective. Among feasible proxies, fewer final OBJ triangles is always better.

PQSS query work, BVH depth, primitive count, collision poses, RSS quality,
triangle-test count, and model-pool behavior are not optimization objectives.

## 2. Current Pipeline

### Stage 0 — Source preparation

Status: retained and currently useful.

Purpose:

- parse the input OBJ;
- remove numerical degeneracies;
- normalize the internal geometry representation;
- establish source-triangle responsibility needed by conservative replacement
  and later audits.

### Stage 1 — Conservative outer preprocessing

Status: retained and currently useful.

Purpose:

- perform the accepted conservative closure/fill operations;
- produce the stable outer reference used by later approximation-error checks;
- remove geometry that becomes internal only when conservative coverage remains
  certified.

### Stage 2 — Surface recognition

Status: exploratory.

Purpose:

- recognize larger certified surface patches;
- currently includes planar and round analytic surface families;
- reduce final triangulation complexity without violating Stage 1 coverage.

Stage 2 is not considered a solved/final algorithm.

### Stage 3 — Hausdorff-limited simplification

Status: exploratory.

Purpose:

- generate and compare additional conservative replacements;
- reject candidates that violate conservative coverage or the directed-error
  threshold;
- among feasible candidates, prefer fewer final OBJ triangles.

Stage 3 is not considered a solved/final algorithm.

### Final triangulation and audit

The selected semantic surface representation is triangulated to an OBJ and
validated.

The current implementation still writes historical intermediate names such as:

- `phase1_hole_filled.obj`;
- `phase2_recognized_surfaces.obj`;
- `primitives.obj`;
- `proxy.obj`.

These names are legacy naming only. Renaming them should be a separate
behavior-neutral refactor.

## 3. Authoritative Files

Current core implementation:

- `src/primitive_mesh_analyzer.cpp`
- `include/pqss_proxy_mesh/primitive_mesh_analyzer.hpp`
- `tools/primitive_mesh_analyze.cpp`

Current core regression:

- `tests/primitive_mesh_analyzer_test.cpp`

Current visualization:

- `viewer/primitive_analysis.html`
- `viewer/lib/`
- `tools/serve_primitive_viewer.py`
- viewer-related tests/tools that are still used

Current geometry audit utility:

- `tools/audit_staged_surface_outputs.py`

Current geometric source corpus retained for regression:

- `test_data/real_scene/source_pool/`

Dependency:

- `deps/clipper2/`

## 4. Removed / Obsolete Directions

The following directions are obsolete and should not be reintroduced as project
goals:

- PQSS real-scene query-cost optimization;
- independent primitive model pools;
- primitive RSS/BVH construction experiments;
- BV-test / triangle-test minimization;
- official-simplified MeshLab comparison as an optimization target;
- old Python proxy-generation/decomposition pipeline;
- model-pool subdivision policy as part of simplification;
- scene-pose-driven geometry decisions.

These remain available through Git history if ever needed for archaeology.

## 5. Current Mathematical / Certification Caveat

The project controls the directed proxy-to-source distance.

Any deterministic surface sampling result is only a sampled estimate unless a
conservative continuous upper bound is also established.

Do not describe a sampled maximum as a proof-level continuous Hausdorff
certificate.

Containment/coverage and directed distance are separate properties:

- small directed distance does not prove conservative outer coverage;
- conservative coverage does not prove the directed-error limit.

## 6. Current Research Questions

Stage 2:

- which semantic surface families should remain;
- how to robustly certify recognition under CAD quantization/noise;
- how to compare alternative recognition results by final triangle count;
- how to avoid greedy recognition that blocks a lower-face feasible result.

Stage 3:

- how to construct better conservative replacement candidates;
- how to search globally or near-globally rather than by unrelated heuristics;
- how to obtain a proof-level continuous directed-error upper bound;
- how to minimize final triangulated face count rather than semantic primitive
  count;
- how to preserve conservative coverage while deleting internal/redundant shell
  pieces.

## 7. Cross-Conversation Rule

GitHub is the persistent source of truth.

Do not rely on:

- ChatGPT conversation history;
- ChatGPT temporary execution storage;
- local uncommitted files;
- screenshots;
- temporary benchmark outputs;
- unstated memory.

Any fact needed to continue the project in a new conversation must be written
into the repository.

At minimum, update this file when any of the following changes:

- project objective or mathematical contract;
- stage definition or stage status;
- accepted algorithm invariant;
- rejected algorithm direction and why it was rejected;
- current blocker;
- important experiment result that changes future work;
- build/run/validation workflow;
- authoritative file locations;
- known correctness limitation;
- next research priority;
- ChatGPT/GitHub collaboration workflow.

## 8. New Conversation Bootstrap

At the start of a new conversation about this repository:

1. read `README.md`;
2. read `AGENTS.md`;
3. read `docs/PROJECT_STATE.md`;
4. read `docs/STAGES.md`;
5. only then inspect the relevant source files.

Do not infer project goals from historical filenames such as `pqss_proxy_mesh`
or from old commits unless the user explicitly asks for history.

## 9. ChatGPT Collaboration Workflow

Current ordinary ChatGPT GitHub connector behavior for this project is treated
as read-only for repository contents.

Observed limitation:

- repository reads work;
- repository write operations through the connector have returned
  `403 Resource not accessible by integration`;
- ordinary ChatGPT conversation should therefore not assume it can push commits
  directly.

When direct writes are unavailable, use this workflow:

1. ChatGPT reads the latest GitHub branch before preparing a change.
2. ChatGPT produces a one-click Windows update package/script when practical.
3. The user runs the script locally.
4. The script clones or updates the repository, creates a branch, applies the
   change, commits, and pushes.
5. ChatGPT re-reads the pushed branch from GitHub before continuing.
6. After a major change is accepted, update this file so the next conversation
   can recover the current state without chat history.

The user's local machine is only a transport for authenticated `git push`.
GitHub remains the persistent project state.

## 10. Current Cleanup Status

Repository cleanup is being performed to remove obsolete PQSS/BVH/Python
systems while keeping the current C++ Stage 0–3 geometry pipeline and
visualization.

The cleanup intentionally avoids mixing a large mechanical rename of the core
analyzer with behavioral cleanup. Legacy namespace/output names may remain
temporarily.

After the cleanup branch is pushed, the next task is to inspect
`src/primitive_mesh_analyzer.cpp` for remaining decision logic or comments that
still rank candidates using obsolete PQSS/BVH/workload metrics and remove or
replace those with the current feasibility-plus-final-face-count objective.
