# Pipeline Stages

This document is the authoritative high-level map of the current code.

## Input and Target

Input:

- one OBJ surface mesh \(M\);
- one user-selected non-negative directed-error limit \(\varepsilon\).

Output:

- one OBJ proxy \(P\).

Optimization target:

\[
\min |F(P)|
\]

subject to conservative outer coverage and

\[
d_{\rightarrow}(P,M)\le\varepsilon.
\]

The exact implementation of a proof-level continuous directed Hausdorff
certificate is still research work. Where the current code uses deterministic
surface sampling, the result must be reported as a sampled maximum, not as an
exact continuous maximum.

## Stage 0 — Source Preparation

**Status: retained / working foundation.**

Responsibilities:

- parse the input OBJ;
- remove numerical zero-area input slivers;
- normalize the internal geometric representation;
- establish the source-triangle responsibility information needed for later
  conservative replacement and audit;
- create the source reference used by the rest of the pipeline.

Stage 0 must not simplify according to PQSS/BVH/query metrics.

## Stage 1 — Conservative Outer Preprocessing

**Status: retained / working foundation.**

Responsibilities:

- perform the currently accepted conservative closure/fill operations;
- remove geometry that becomes provably internal only when coverage remains
  certified;
- produce the stable conservative surface against which later outward
  approximation error is measured.

The current implementation emits the historical file
`phase1_hole_filled.obj`.

Stage 1 changes the conservative reference geometry. Its operations therefore
must be justified by outer-coverage semantics, not hidden inside the Stage 2/3
Hausdorff budget.

## Stage 2 — Surface Recognition

**Status: exploratory.**

Goal:

Reduce representation complexity without giving up conservative coverage by
recognizing larger certified surface patches.

Current families include planar polygons and round analytic surface families
such as disks, annuli, cylindrical bands, and conical bands.

Important distinction:

- semantic patch count is not the optimization objective;
- final triangulated OBJ face count is the quantity to minimize.

The current implementation emits the historical file
`phase2_recognized_surfaces.obj`.

Research questions include:

- which surface families are worth keeping;
- how to certify recognition robustly under CAD quantization;
- how to compare alternative recognitions by final triangulation cost;
- how to avoid greedy early recognition blocking a lower-face feasible result.

## Stage 3 — Directed-Error-Limited Simplification

**Status: exploratory.**

Goal:

Starting from the conservative Stage 1 reference and Stage 2 candidates, find
additional outer replacements that minimize final triangle count while obeying
the user threshold \(\varepsilon\).

A candidate must satisfy, independently:

1. conservative coverage/containment;
2. directed proxy-to-source distance at most \(\varepsilon\).

Then it competes by **final triangle count**.

No candidate may be selected because it improves PQSS workload, RSS quality,
BVH depth, primitive count, or a collision benchmark.

The current implementation emits semantic output in `primitives.obj`.

## Final Triangulation and Audit

The selected semantic surfaces are triangulated into `proxy.obj`.

Final acceptance must report at least:

- input triangle count;
- output triangle count;
- containment/coverage audit result;
- configured directed-error limit;
- measured/certified directed error and the method used;
- whether that error is exact or sampling-based.

An \(n\)-vertex simple planar polygon should triangulate to \(n-2\) triangles
without unnecessary Steiner vertices unless a later correctness requirement
forces otherwise.

## Legacy Naming

Several source symbols and output filenames still contain `pqss`, `primitive`,
or `phaseN` terminology. Repository cleanup deliberately removes the standalone
PQSS systems first while keeping the working Stage 0–3 implementation stable.
Renaming core symbols and intermediate files should be done in a separate,
mechanical change after algorithm behavior is covered by tests.
