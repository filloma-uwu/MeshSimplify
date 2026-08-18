# MeshSimplify Development Contract

These rules define the current project. Historical PQSS/BVH goals are obsolete.

1. **Single objective.** Minimize the triangle count of the final output OBJ.
   Runtime, BVH quality, collision-query work, primitive count, and build time
   are not optimization objectives.

2. **Hard feasibility constraints.** Every accepted output must remain a
   conservative outer replacement of the input geometry and must satisfy the
   user-provided directed proxy-to-source Hausdorff threshold.

3. **The Hausdorff limit is not a score.** Do not trade triangle count against
   error with a weighted objective. A candidate is either feasible under the
   threshold or infeasible; among feasible candidates, prefer fewer triangles.

4. **No model-specific behavior.** Never branch on model ID, file name, known
   scene role, hand-picked coordinates, or a screenshot. Fixes must be expressed
   as general geometry rules.

5. **Stages are explicit.**
   - Stage 0 and Stage 1 are the current working foundation.
   - Stage 2 and Stage 3 are exploratory.
   - Do not present an exploratory heuristic as a settled invariant.

6. **Preserve source responsibility.** A source face may disappear from the
   exported surface only when conservative coverage is certified by retained or
   replacement geometry.

7. **No PQSS decision rules.** PQSS workload, BV tests, triangle tests, BVH
   depth, model-pool subdivision, collision poses, and RSS quality must not
   accept, reject, rank, or prune geometry candidates.

8. **Face count means final OBJ triangles.** Semantic primitive count is useful
   for debugging but is not the objective. Candidate comparison must ultimately
   account for the triangulated output.

9. **Distance direction is proxy to source.**
   The controlled quantity is
   `sup_{p in proxy} dist(p, source)`.
   Any sampled audit must be described as an estimate; do not label a sampled
   maximum as a proof-level continuous Hausdorff bound.

10. **Containment and distance are separate checks.** A small directed distance
    does not prove conservative outer coverage. Passing containment does not
    prove the directed error threshold.

11. **Geometry, not OBJ topology, is authoritative.** Coincident or intersecting
    surfaces may use unrelated vertex indices. Do not rely on shared indices as
    a substitute for geometric tests.

12. **Keep planar output canonical.** Geometrically coplanar covered regions
    should be unioned before final triangulation so duplicate/internal coplanar
    faces and unnecessary internal edges do not survive.

13. **Regression tests are model-independent.** Every algorithm change should
    add or update a synthetic/general geometric regression whenever practical.

14. **The viewer is diagnostic, not an oracle.** Visual appearance can reveal a
    bug, but acceptance is determined by geometric certificates, final triangle
    count, and the directed-error audit.

15. **Do not reintroduce removed historical systems.** Independent PQSS model
    pools, RSS/BVH comparison code, official-simplified scene scoring, and the
    old Python proxy-generation pipeline belong in Git history, not the current
    tree.


## Cross-Conversation Persistence

16. GitHub is the project source of truth. Any information required to continue
    development in a fresh conversation must be recorded in
    `docs/PROJECT_STATE.md`; do not rely on chat history or temporary execution
    storage.

17. A new ChatGPT conversation must read `README.md`, `AGENTS.md`,
    `docs/PROJECT_STATE.md`, and `docs/STAGES.md` before changing the algorithm.

18. Update `docs/PROJECT_STATE.md` whenever a major design decision, rejected
    approach, experiment result, blocker, build/run workflow, known correctness
    limitation, stage status, or next-step priority changes.

19. When direct GitHub writes are unavailable, prefer a one-click local update
    script that clones/updates, branches, applies, commits, and pushes. Avoid
    requiring the user to manually reproduce multi-command Git workflows.
