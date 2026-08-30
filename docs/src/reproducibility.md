# Reproducibility

## Julia compatibility and `ProbabilityBoundsAnalysis`

The package requires Julia 1.12+. Its dependency `ProbabilityBoundsAnalysis`
(0.2.11, the current release) **does not precompile on Julia 1.12** — it extends
`Distributions.Frechet` without importing it, which is a hard error during
module precompilation on 1.12. The package still loads and runs correctly; it
just falls back to interpreted mode, so the first `using` in a session takes
~40 s and there is no on-disk cache. This is upstream and expected.

To skip the doomed precompile pass (for example in CI), set

```
JULIA_PKG_PRECOMPILE_AUTO=0
```

## Test suite

`Pkg.test("InformationPropagationAnalysis")` runs a curated suite (`Test` is the
only test-only dependency). It has no dependency on CUDD or other BDD packages —
the exact oracle it checks against is a self-contained brute-force state
enumeration (`test/oracle.jl`), tractable because the fixtures are small.

The fixtures and the provenance of each expected value are documented in
`test/fixtures/README.md`. In summary:

| Fixture | Check | Ground truth |
|---|---|---|
| synthetic diamonds (single, nested) | `Probability` belief == exact | brute-force state enumeration |
| `counterexample-n15` | `belief(15) == 0.71858899` | CUDD ROBDD, from the diamond-rewrite validation campaign |
| `power-network` | float beliefs vs golden; interval brackets float | 4-way oracle agreement (IPA / path-enum / CUDD / BDD.jl) |
| `psplib-j301_1` | `LONGEST_PATH` makespan == 38.0 | PSPLIB's published `MPM-Time` |
| `psplib-j301_1` (interval) | `interval_analyze_split` no-throw + slack bounds | recorded interval-CPM oracle, MC-checked |
| `genrmf-dag-small` | 3 solvers agree; max-flow == min-cut | `GraphsFlows.jl` |

## Regenerating golden data

The golden belief vectors are produced by the current build (they are only
meaningful because the algorithm is independently validated — see the table
above). Regenerate with:

```
julia --project=. test/gen_golden.jl
```
