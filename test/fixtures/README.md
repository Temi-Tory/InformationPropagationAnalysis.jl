# Test fixtures

Small networks with known-correct results, used by `../runtests.jl`. All are copied
verbatim from the source project's validation corpus; provenance below.

## `counterexample-n15/`

The canonical bug-#1 case from the diamond-decomposition rewrite. Single source
(node 1), 15 nodes, 23 edges; all node priors and edge probabilities 0.9.

- **Assertion:** `belief(15) ≈ 0.71858899`.
- **Provenance:** exact value established during the rewrite validation campaign,
  cross-checked against an exact CUDD ROBDD (`validation/validate_broad.jl` → 0 wrong
  over 114 random+mutant DAGs). The pre-rewrite hybrid-reuse code gave `0.66916` here.

## `power-network/`

A 23-node, 27-edge synthetic power-distribution DAG. Ships node priors / edge
probabilities in `float/`, `interval/` and `pbox/` form, edge capacities in
`capacity/`, and CPM durations in `cpm/`.

- **Assertions:** float beliefs match `expected-float-beliefs.csv` (golden master,
  regenerate with `../gen_golden.jl`) *and* `v0.1.0-beliefs.json`; `belief ≤ prior`
  elementwise; the interval-valued run brackets the float run.
- **Provenance:** float beliefs on this network agree to ~1e-16 across four independent
  implementations (IPA, path-enumeration + inclusion-exclusion, CUDD, BinaryDecisionDiagrams.jl).
  The golden CSV is this v0.2.0 build's output; `v0.1.0-beliefs.json` is the belief vector
  the v0.1.0 (5-flat-file) release recorded — the two agree to 1.1e-16 across all 23 nodes,
  so the rewrite reproduces the previous release bit-for-bit here.

## `psplib-j301_1/`

PSPLIB single-mode instance `j301_1.sm` (the `j30` set: 30 real activities + 2 dummy
endpoints = 32 nodes, 48 precedence edges). Resources ignored; precedence-only schedule.

- **Assertions:** LONGEST_PATH project value == `38.0` (PSPLIB's own published
  `MPM-Time`); `interval_analyze_split` does not throw on the degenerate start node and
  its per-node slack intervals match `expected-float-bounds.csv`.
- **Provenance:** `expected-float-bounds.csv` is `validation/cpm_v2/psplib_j301_1_float_bounds.csv`,
  the interval-CPM oracle output, itself MC-checked (`validation/cpm_v2/run_mc_check_psplib_j301_1.jl`).
  The interval-split assertion is the regression guard for the `interval_analyze_split`
  degenerate-node fix.

## `genrmf-dag-small/`

A DIMACS `genrmf`-family max-flow instance, converted to a DAG.

- **Assertions:** the three max-flow solvers (Dinic, Edmonds–Karp, push–relabel) agree;
  max-flow == min-cut capacity; `analyze_all` returns a `FlowCapacityResult`.
- **Provenance:** `expected-maxflow.csv` is `validation/flow/dimacs/dimacs_validation_summary.csv`;
  the IPA solvers were validated against `GraphsFlows.jl` there (all PASS).
